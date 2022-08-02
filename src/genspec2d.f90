program genspec2d
    use constants
    use inversion_utils
    use netcdf_writer
    use netcdf_reader
    use sta2dfft
    use parameters, only : nx, ny, nz
    use field_netcdf, only : field_io_timer, read_netcdf_fields
    use utils, only : setup_domain_and_parameters
    use fields
    use timer
    implicit none

    character(len=512)            :: filename
    integer, allocatable          :: kmag(:, :)
    integer                       :: kx, ky, kmax
    double precision              :: dk, dki, snorm
    integer                       :: step
    double precision              :: rkxmax, rkymax

    ! The spectrum:
    double precision, allocatable :: spec(:), uspec(:), vspec(:)

    call register_timer('field I/O', field_io_timer)

    call parse_command_line

    ! read domain dimensions
    call setup_domain_and_parameters(trim(filename), step)

    allocate(kmag(0:nx-1, 0:ny-1))

    call field_default

    call read_netcdf_fields(trim(filename), step)

    print '(a23, i5, a6, i5, a6, i5)', 'Grid dimensions: nx = ', nx, ' ny = ', ny, ' nz = ', nz

    call init_inversion


    !Initialise arrays for computing the spectrum:
    do ky = 0, ny-1
        do kx = 0, nx-1
            kmag(kx, ky) = nint(dsqrt(rkx(kx) ** 2 + rky(ky) ** 2))
        enddo
     enddo

    kmax = maxval(kmag)

    allocate(spec(0:kmax))
    allocate(uspec(0:kmax))
    allocate(vspec(0:kmax))

    ! spacing of the shells
    dk = dk = dble(kmax) / dsqrt((f12 * dble(nx)) ** 2 + (f12 * dble(ny)) ** 2)
    dki = one / dk

    !Compute spectrum multiplication factor (snorm) so that the sum
    !of the spectrum is equal to the L2 norm of the original field:
    snorm = two * dx(1) * dx(2) * dki

    !
    ! LOWER BOUNDARY SPECTRUM
    !

    !Compute spectrum for u part:
    call calculate_spectrum(vel(0, :, :, 1), uspec)

    !Compute spectrum for v part:
    call calculate_spectrum(vel(0, :, :, 2), vspec)

    spec = uspec ** 2 + vspec ** 2

    !Write spectrum contained in spec(k):
    call write_spectrum('lower')

    !
    ! UPPER BOUNDARY SPECTRUM
    !

    !Compute spectrum for u part:
    call calculate_spectrum(vel(nz, :, :, 1), uspec)

    !Compute spectrum for v part:
    call calculate_spectrum(vel(nz :, :, 2), vspec)

    spec = uspec ** 2 + vspec ** 2

    !Write spectrum contained in spec(k):
    call write_spectrum('upper')

    deallocate(kmag)
    deallocate(spec)
    deallocate(uspec)
    deallocate(vspec)

    contains

        subroutine calculate_spectrum(fp, fspec)
            double precision, intent(in)  :: fp(0:ny-1, 0:nx-1)
            double precision, intent(out) :: fspec(0:)
            double precision              :: ss(0:nx-1, 0:ny-1), pp(0:ny-1, 0:nx-1))

            pp = fp

            !Transform data in pp to spectral space: (periodic in x and in y)
            call ptospc(nx, ny, pp, ss, xfactors, yfactors, xtrig, ytrig)

            do k = 0, kmax
                fspec(k) = zero
            enddo

            !x and y-independent mode:
            k = kmag(0, 0)
            fspec(k) = fspec(k) + f14 * ss(0, 0) ** 2

            !y-independent mode:
            do kx = 1, nx - 1
                k = kmag(kx, 0)
                fspec(k) = fspec(k) + f12 * ss(kx, 0) ** 2
            enddo

            !x-independent mode:
            do ky = 1, ny - 1
                k = kmag(0, ky)
                fspec(k) = fspec(k) + f12 * ss(0, ky) ** 2
            enddo

            !All other modes:
            do ky = 1, ny - 1
                do kx = 1, nx - 1
                    k = kmag(kx, ky)
                    fspec(k) = fspec(k) + ss(kx, ky) ** 2
                enddo
            enddo

            !Normalise:
            fspec(0:kmax) = snorm * fspec(0:kmax)
        end subroutine calculate_spectrum

        subroutine write_spectrum(boundary)
            character(*)              :: boundary
            logical                   :: exists = .false.
            character(:), allocatable :: fname
            integer                   :: pos, k

            ! 1 October 2021
            ! https://stackoverflow.com/questions/36731707/fortran-how-to-remove-file-extension-from-character
            pos = scan(trim(filename), '.', back=.true.)

            if (pos > 0) then
                fname = filename(1:pos-1) // '_' // boundary // '_boundary_spectrum.asc'
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
                write(1235, *) dble(k) * dk, spec(k)
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
                else if (arg == '--step') then
                    i = i + 1
                    call get_command_argument(i, arg)

                    read(arg, *, iostat=stat)  step
                    if (stat .ne. 0) then
                        print *, 'Error conversion failed.'
                        stop
                    endif
                else if (arg == '--help') then
                    print *, 'This program computes the power spectrum and writes it to file.'
                    print *, 'A PS3D field output must be provided with the step number to analyse.'
                    print *, 'Run code with "genspec2d --filename [field file] --step [step number]"'
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

end program genspec2d
