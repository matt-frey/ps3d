program genspec
    use constants, only : pi, twopi, f14, f12, zero, one, four
    use netcdf_reader
    use netcdf_writer
    use inversion_utils
    use parameters, only : nx, ny, nz
    use timer
    implicit none

    character(len=512) :: filename
    integer            :: step
    integer            :: iz

    call register_timer('field I/O', field_io_timer)

    call parse_command_line

    ! read domain dimensions
    call setup_domain_and_parameters(trim(filename))

    call field_default

    call read_netcdf_fields(trim(filename))

    print '(a23, i5, a6, i5, a6, i5)', 'Grid dimensions: nx = ', nx, ' ny = ', ny, ' nz = ', nz

    ! use some dummy values for bbdif, nnu and prediss
    call init_inversion(zero, 3, 10.0d0)


    contains

        subroutine write_spectrum
            logical                   :: exists = .false.
            character(:), allocatable :: fname
            integer                   :: pos

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

            do k = 1, kmax
                write(1235, *) k * delk, spec(k)
            enddo

            close(1235)

        end subroutine write_spectrum


        ! Get the file name provided via the command line
        subroutine parse_command_line
            integer            :: i, stat
            character(len=512) :: arg
            logical            :: exists

            step = -1
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

                    ! 1 October 2021
                    ! https://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90
                    read(arg, *, iostat=stat)  step
                    if (stat .ne. 0) then
                        print *, 'Error conversion failed.'
                        stop
                    endif
                else if (arg == '--help') then
                    print *, 'This program computes the power spectrum and writes it to file.'
                    print *, 'A PS3D field output must be provided with the step number to analyse.'
                    print *, 'Run code with "genspec --filename [field file]" --step [step number]'
                    stop
                endif
                i = i+1
            enddo

            if ((filename == '') .or. (step == -1) ) then
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
