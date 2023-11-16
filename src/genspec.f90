program genspec
    use options, only : field_step, filename
    use constants
    use netcdf_reader
    use netcdf_writer
    use inversion_utils
    use sta2dfft, only : dct, dst
    use sta3dfft
    use parameters, only : nx, ny, nz, ncelli
    use field_netcdf, only : field_io_timer, read_netcdf_fields
    use utils, only : setup_domain_and_parameters
    use fields
    use mpi_environment
    use timer
    implicit none

    integer, allocatable          :: kmag(:, :, :)
    double precision, allocatable :: spec(:)
    integer, allocatable          :: num(:)
    integer                       :: nc, kx, ky, kz, m, kmax
    double precision              :: dk, dki, prefactor, snorm
    double precision              :: ke ! kinetic energy

    call mpi_env_initialise

    call register_timer('field I/O', field_io_timer)

    call parse_command_line

    ! read domain dimensions
    call setup_domain_and_parameters

    allocate(kmag(0:nz, box%lo(2):box%hi(2), box%lo(1):box%lo(1)))

    call field_default

    call read_netcdf_fields(trim(filename), field_step)

    if (world%rank == world%root) then
        print '(a23, i5, a6, i5, a6, i5)', 'Grid dimensions: nx = ', nx, ' ny = ', ny, ' nz = ', nz
    endif

    call init_inversion

    ! calculate domain-average kinetic energy
    ke = get_kinetic_energy()

    ! (1) compute the 3D spectrum of each velocity component:
    do nc = 1, 3
        call fftxyp2s(vel(:, :, :, nc), svel(:, :, :, nc))
    enddo
    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            call dct(1, nz, svel(0:nz, ky, kx, 1), ztrig, zfactors) ! u
            call dct(1, nz, svel(0:nz, ky, kx, 2), ztrig, zfactors) ! v
            call dst(1, nz, svel(1:nz, ky, kx, 3), ztrig, zfactors) ! w
        enddo
    enddo

    ! (2) sum the squared spectral amplitudes into radial shells in total wavenumber K = sqrt{kx^2 + ky^2 + kz^2}
    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            do kz = 0, nz
                kmag(kz, ky, kx) = nint(dsqrt(rkx(kx) ** 2 + rky(ky) ** 2 + rkz(kz) ** 2))
            enddo
        enddo
     enddo

    kmax = maxval(kmag)

    call MPI_Allreduce(MPI_IN_PLACE,            &
                       kmax,                    &
                       1,                       &
                       MPI_DOUBLE_PRECISION,    &
                       MPI_MAX,                 &
                       world%comm,              &
                       world%err)

    allocate(spec(0:kmax))
    allocate(num(0:kmax))

    ! spacing of the shells
    dk = dble(kmax) / dsqrt((f12 * dble(nx)) ** 2 + (f12 * dble(ny)) ** 2 + dble(nz) ** 2)
    dki = one / dk

    ! (3) accumulate spectrum
    spec = zero
    num = 0

    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            do kz = 0, nz
                m = int(dble(kmag(kz, ky, kx)) * dki)
                spec(m) = spec(m) &
                        + svel(kz, ky, kx, 1) ** 2 + svel(kz, ky, kx, 2) ** 2 + svel(kz, ky, kx, 3) ** 2
                num(m) = num(m) + 1
            enddo
        enddo
    enddo

    call mpi_blocking_reduce(spec, MPI_SUM, world)
    call mpi_blocking_reduce(num, MPI_SUM, world)

    prefactor = 4.0d0 / 3.0d0 * pi * dK ** 3

    if (world%rank == world%root) then
        do m = 0, kmax
            if (num(m) > 0) then
                spec(m) = spec(m) * prefactor * dble((m+1) ** 3 - m ** 3) / dble(num(m))
            else
                print *, "Bin", m, " is empty!"
            endif
        enddo

        ! calculate spectrum normalisation factor (snorm)
        ! that ensures Parceval's identity, so that the spectrum S(K)
        ! has the property that its integral over K gives the total kinetic energy
        snorm = ke / sum(spec * dk)

        ! normalise the spectrum
        spec = spec * snorm

        call write_spectrum
    endif

    deallocate(kmag)
    deallocate(num)
    deallocate(spec)

    call mpi_env_finalise

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

                    read(arg, *, iostat=stat)  field_step
                    if (stat .ne. 0) then
                        call mpi_stop('Error conversion failed.')
                    endif
                else if (arg == '--help') then
                    if (world%rank == world%root) then
                        print *, 'This program computes the power spectrum and writes it to file.'
                        print *, 'A PS3D field output must be provided with the step number to analyse.'
                        print *, 'Run code with "genspec --filename [field file] --step [step number]"'
                    endif
                endif
                i = i+1
            enddo

            if (filename == '') then
                call mpi_stop('No file or step provided. Run code with "genspec --help"')
            endif

            ! check if correct file is passed
            stat = index(trim(filename), '_fields.nc', back=.true.)
            if (stat == 0) then
                call mpi_stop("Error: No PS3D field output file provided.")
            endif

            ! check if file exsits
            inquire(file=trim(filename), exist=exists)
            if (.not. exists) then
                call mpi_stop("Error: File '" // trim(filename) // "' does not exist.")
            endif
        end subroutine parse_command_line

end program genspec
