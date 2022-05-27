program genspec
    use constants
    use netcdf_reader
    use netcdf_writer
    use inversion_utils
    use parameters, only : nx, ny, nz, ngrid, vcell
    use field_netcdf, only : field_io_timer, read_netcdf_fields
    use utils, only : setup_domain_and_parameters
    use fields
    use timer
    implicit none

    character(len=512)            :: filename
    integer, allocatable          :: kmag(:, :, :)
    double precision, allocatable :: spec(:)
    integer, allocatable          :: num(:)
    integer                       :: iz, nc, kx, ky, kz, m, kmax
    double precision              :: dk, k, prefactor, snorm
    double precision              :: ens ! enstrophy

    call register_timer('field I/O', field_io_timer)

    call parse_command_line

    ! read domain dimensions
    call setup_domain_and_parameters(trim(filename))

    allocate(kmag(0:nz, 0:ny-1, 0:nx-1))
    allocate(spec(0:max(nz, ny-1, nx-1)))
    allocate(num(0:max(nz, ny-1, nx-1)))

    call field_default

    call read_netcdf_fields(trim(filename))

    print '(a23, i5, a6, i5, a6, i5)', 'Grid dimensions: nx = ', nx, ' ny = ', ny, ' nz = ', nz

    ! use some dummy values for bbdif, nnu and prediss
    call init_inversion(zero, 3, 10.0d0)

    ! (1) compute the 3D spectrum of each vorticity component assuming cosine in z
    do nc = 1, 3
        call fftczp2s(vortg(:, :, :, nc), svortg(:, :, :, nc))
    enddo

    ! (2) sum the squared spectral amplitudes into radial shells in total wavenumber K = sqrt{kx^2 + ky^2 + kz^2}
    do kx = 0, nx-1
        do ky = 0, ny-1
            do kz = 0, nz
                kmag(kz, ky, kx) = nint(dsqrt(rkx(kx+1) ** 2 + rky(ky+1) ** 2 + rkz(kz) ** 2))
            enddo
        enddo
    enddo

    kmax = maxval(kmag)

    ! spacing of the shells
    dk = kmax / dsqrt((f12 * dble(nx)) ** 2 + (f12 * dble(ny)) ** 2 + dble(nz) ** 2)

    ! (3) accumulate spectrum
    spec = zero
    num = 0

    do kx = 0, nx-1
        do ky = 0, ny-1
            do kz = 0, nz
                k = kmag(kz, ky, kx) / dk
                m = int(k)
                spec(m) = svortg(kz, kx, ky, 1) ** 2 + svortg(kz, kx, ky, 2) ** 2 + svortg(kz, kx, ky, 3) ** 2
                num(m) = num(m) + 1
            enddo
        enddo
    enddo

    prefactor = 4.0d0 / 3.0d0 * pi * dK ** 3

    do m = 0, size(spec)
        if (num(m) > 0) then
            spec(m) = spec(m) * prefactor * ((m+1) ** 3 - m ** 3) / num(m)
        else
            print *, "Empty bin!"
            stop
        endif
    enddo

    ! calculate enstrohpy
    ens = f12 * sum(vortg(1:nz-1, :, :, 1) ** 2     &
                  + vortg(1:nz-1, :, :, 2) ** 2     &
                  + vortg(1:nz-1, :, :, 3) ** 2)    &
        + f14 * sum(vortg(0,      :, :, 1) ** 2     &
                  + vortg(0,      :, :, 2) ** 2     &
                  + vortg(0,      :, :, 3) ** 2)    &
        + f14 * sum(vortg(nz,     :, :, 1) ** 2     &
                  + vortg(nz,     :, :, 2) ** 2     &
                  + vortg(nz,     :, :, 3) ** 2)

    ! note: ngrid = nx * ny * (nz+1) and vcell = dx * dy * dz
    ens = ens * ngrid * vcell

    ! calculate spectrum normalisation factor (snorm)
    ! that ensures Parceval's identity, so that the spectrum S(K)
    ! has the property that its integral over K gives the total enstrophy
    snorm = ens / sum(spec * dk)

    ! normalise the spectrum
    spec = spec * snorm


    contains

        subroutine write_spectrum
            logical                   :: exists = .false.
            character(:), allocatable :: fname
            integer                   :: pos, k

            ! 1 October 2021
            ! https://stackoverflow.com/questions/36731707/fortran-how-to-remove-file-extension-from-character
            pos = scan(trim(filename), '.', back=.true.)

            if (pos > 0) then
                fname = filename(1:pos-1) // '_spectrum.asc'
            else
                print *, "Error in reading the filename. File extension not found."
            endif

            inquire(file=fname, exist=exists)
            if (exists) then
                print *, "Error: File '" // trim(fname) // "' already exists."
                stop
            else
                open(unit=1235, file=fname, status='replace')
                write(1235, *) '# The power spectrum.'
                write(1235, *) '# The first column is the wavenumber, the second column the spectrum.'
                write(1235, *) '#         k   P(k)'
            endif

            do k = 0, kmax
                write(1235, *) k * dk, spec(k)
            enddo

            close(1235)

        end subroutine write_spectrum


        ! Get the file name provided via the command line
        subroutine parse_command_line
            integer            :: i, stat
            character(len=512) :: arg
            logical            :: exists

            filename = ''
            i = 0
            do
                call get_command_argument(i, arg)
                if (len_trim(arg) == 0) then
                    exit
                endif

                if (arg == '--filename') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    filename = trim(arg)
                else if (arg == '--help') then
                    print *, 'This program computes the power spectrum and writes it to file.'
                    print *, 'A PS3D field output must be provided.'
                    print *, 'Run code with "genspec --filename [field file]"'
                    stop
                endif
                i = i+1
            enddo

            if (filename == '') then
                print *, 'No file or step provided. Run code with "genspec --help"'
                stop
            endif

            ! check if correct file is passed
            stat = index(trim(filename), '_fields.nc', back=.true.)
            if (stat == 0) then
                print *, "Error: No PS3D field output file provided."
                stop
            endif

            ! check if file exsits
            inquire(file=trim(filename), exist=exists)
            if (.not. exists) then
                print *, "Error: File '" // trim(filename) // "' does not exist."
                stop
            endif
        end subroutine parse_command_line

end program genspec
