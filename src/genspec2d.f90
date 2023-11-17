program genspec2d
    use options, only : field_step, filename
    use constants
    use inversion_utils
    use netcdf_writer
    use netcdf_reader
    use sta2dfft
    use sta3dfft
    use parameters, only : nx, ny, nz
    use field_netcdf, only : field_io_timer, read_netcdf_fields
    use utils, only : setup_domain_and_parameters
    use fields
    use mpi_environment
    use mpi_timer
    implicit none

    integer, allocatable          :: kmag(:, :)
    integer                       :: kx, ky, kmax, k
    double precision              :: dk, dki, snorm, ke

    ! The spectrum:
    double precision, allocatable :: spec(:)

    call mpi_env_initialise

    if (world%size > 1) then
        call mpi_stop("This program does not run in parallel due to call_ptospc!")
    endif

    call register_timer('field I/O', field_io_timer)

    call parse_command_line

    ! read domain dimensions
    call setup_domain_and_parameters

    allocate(kmag(box%lo(2):box%hi(2), box%lo(1):box%lo(1)))

    call field_default

    call read_netcdf_fields(trim(filename), field_step)

    if (world%rank == world%root) then
        print '(a23, i5, a6, i5, a6, i5)', 'Grid dimensions: nx = ', nx, ' ny = ', ny, ' nz = ', nz
    endif

    call init_inversion


    !Initialise arrays for computing the spectrum:
    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            kmag(ky, kx) = nint(dsqrt(rkx(kx) ** 2 + rky(ky) ** 2))
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

    ! spacing of the shells
    dk = dble(kmax) / dsqrt((f12 * dble(nx)) ** 2 + (f12 * dble(ny)) ** 2)
    dki = one / dk

    !
    ! LOWER BOUNDARY SPECTRUM
    !
    call calculate_spectrum(vel(0, :, :, 1), vel(0, :, :, 2))

    if (world%rank == world%root) then
        !Write spectrum contained in spec(k):
        call write_spectrum('lower')
    endif

    !
    ! UPPER BOUNDARY SPECTRUM
    !
    spec = zero

    call calculate_spectrum(vel(nz, :, :, 1), vel(nz, :, :, 2))

    if (world%rank == world%root) then
        !Write spectrum contained in spec(k):
        call write_spectrum('upper')
    endif

    deallocate(kmag)
    deallocate(spec)

    call mpi_env_finalise

    contains

        subroutine calculate_spectrum(u, v)
            double precision, intent(in) :: u(box%lo(2):box%hi(2), box%lo(1):box%lo(1)), &
                                            v(box%lo(2):box%hi(2), box%lo(1):box%lo(1))
            double precision             :: uspec(0:kmax), vspec(0:kmax)

            !Compute spectrum for u part:
            call calculate_spectrum_contribution(u, uspec)

            !Compute spectrum for v part:
            call calculate_spectrum_contribution(v, vspec)

            ! uspec and vspec are already squared
            spec = uspec + vspec

            ! calculate domain-average kinetic energy at surface with u and v part only:
            ke = f12 * sum(u ** 2 + v ** 2)
            ke = ke / dble(nx * ny)

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               ke,                      &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)


            ! calculate spectrum normalisation factor (snorm)
            ! that ensures Parceval's identity, so that the spectrum S(K)
            ! has the property that its integral over K gives the total kinetic energy
            snorm = ke / sum(spec * dk)

            !Normalise:
            spec = snorm * spec

        end subroutine calculate_spectrum

        subroutine calculate_spectrum_contribution(fp, fspec)
            double precision, intent(in)  :: fp(box%lo(2):box%hi(2), box%lo(1):box%lo(1))
            double precision, intent(out) :: fspec(0:kmax)
            double precision              :: ss(box%lo(2):box%hi(2), box%lo(1):box%lo(1)), &
                                             pp(box%lo(2):box%hi(2), box%lo(1):box%lo(1))
            integer                       :: num(0:kmax)

            pp = fp

            !Transform data in pp to spectral space: (periodic in x and in y)
            call call_ptospc(pp, ss)

            do k = 0, kmax
                fspec(k) = zero
                num(k) = 0
            enddo

            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    k = int(dble(kmag(ky, kx)) * dki)
                    fspec(k) = fspec(k) + ss(ky, kx) ** 2
                    num(k) = num(k) + 1
                enddo
            enddo

            call mpi_blocking_reduce(fspec, MPI_SUM, world)
            call mpi_blocking_reduce(num, MPI_SUM, world)

            do k = 0, kmax
                if (num(k) > 0) then
                   fspec(k) = fspec(k) / dble(num(k))
                else
                   print *, "Bin", k, " is empty!"
                endif
            enddo
        end subroutine calculate_spectrum_contribution

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

                    read(arg, *, iostat=stat)  field_step
                    if (stat .ne. 0) then
                        call mpi_stop('Error conversion failed.')
                    endif
                else if (arg == '--help') then
                    if (world%rank == world%root) then
                        print *, 'This program computes the power spectrum and writes it to file.'
                        print *, 'A PS3D field output must be provided with the step number to analyse.'
                        print *, 'Run code with "genspec2d --filename [field file] --step [step number]"'
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

end program genspec2d
