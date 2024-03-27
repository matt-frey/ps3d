! =============================================================================
! This module contains global options that can be set at runtime by the user.
! =============================================================================
module options
    use constants, only : zero, one, two, pi, four, twopi
    use netcdf_writer
    use mpi_utils, only : mpi_stop
    implicit none
    !
    ! global options
    !

    character(len=512) :: field_file = ''

    integer :: field_step = -1

    ! print more info if true
    logical :: verbose = .false.

    ! configuration file
    character(len=512) :: filename = ''

    ! time integrator
    character(len=3) :: stepper = 'IMPL-RK4' ! 'LS-RK4', 'LS-RK3' or 'CN2'
    !
    ! output options
    !
    type info
        double precision                  :: field_freq         = one
        logical                           :: write_fields       = .true.
        character(len=32), dimension(128) :: field_list         = ''
        double precision                  :: field_stats_freq   = one
        logical                           :: write_field_stats  = .true.
        logical                           :: overwrite          = .false.
        character(len=512)                :: basename           = ''
    end type info

    type(info) :: output

    !(Hyper)viscosity parameters:
    type visc_type
        integer :: nnu
        double precision :: prediss
        ! If nnu = 1, this is the molecular viscosity case.  Then, we
        ! choose the viscosity nu = prediss*((b_max-b_min)/k_{x,max}^3)
        ! where k_{x_max} is the maximum x wavenumber.
        ! Note: prediss = 2 is recommended.
        ! ----------------------------------------------------------------
        ! If nnu > 1, this is the hyperviscosity case.  Then, the damping
        ! rate is prediss*zeta_char*(k/k_max)^(2*nnu) on wavenumber k
        ! where k_max is the maximum x or y wavenumber and zeta_char is
        ! a characteristic vorticity (see subroutine adapt of strat.f90).
        ! Note: nnu = 3 and prediss = 10 are recommended.
    end type visc_type

    type(visc_type) :: viscosity

    ! time limit
    type time_info_type
        double precision :: initial     = zero       ! initial time
        double precision :: limit       = zero       ! time limit
        double precision :: alpha       = 0.1d0      ! factor for adaptive time stepping with strain and buoyancy
                                                     ! gradient
                                                     ! Note: EPIC-2D paper recommends alpha = 0.2 for ls-rk4 method
        logical          :: precise_stop = .false.   ! stop at the exact limit
    end type time_info_type

    type(time_info_type) :: time


    contains
        ! parse configuration file
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer :: ios
            integer :: fn = 1
            logical :: exists = .false.

            ! namelist definitions
            namelist /PS3D/ field_file, field_step, stepper, viscosity, output, time

            ! check whether file exists
            inquire(file=filename, exist=exists)

            if (exists .eqv. .false.) then
                call mpi_stop(&
                    'Error: input file "' // trim(filename) // '" does not exist.')
            endif

            ! open and read Namelist file.
            open(action='read', file=filename, iostat=ios, newunit=fn)

            read(nml=PS3D, iostat=ios, unit=fn)

            if (ios /= 0) then
                call mpi_stop('Error: invalid Namelist format.')
            end if

            close(fn)

            ! check whether NetCDF files already exist
            inquire(file=output%basename, exist=exists)

            if (exists) then
                call mpi_stop(&
                    'Error: output file "' // trim(output%basename) // '" already exists.')
            endif

        end subroutine read_config_file

        subroutine write_netcdf_options(ncid)
            integer, intent(in) :: ncid

#ifdef ENABLE_VERBOSE
            call write_netcdf_attribute(ncid, "verbose", verbose)
#endif

            if (viscosity%nnu == 1) then
                call write_netcdf_attribute(ncid, "viscosity", "molecular")
            else
                call write_netcdf_attribute(ncid, "viscosity", "hyperviscosity")
            endif

            call write_netcdf_attribute(ncid, "stepper", stepper)
            call write_netcdf_attribute(ncid, "nnu", viscosity%nnu)
            call write_netcdf_attribute(ncid, "prediss", viscosity%prediss)

            call write_netcdf_attribute(ncid, "field_freq", output%field_freq)
            call write_netcdf_attribute(ncid, "write_fields", output%write_fields)
            call write_netcdf_attribute(ncid, "field_stats_freq", output%field_stats_freq)
            call write_netcdf_attribute(ncid, "write_field_stats", output%write_field_stats)
            call write_netcdf_attribute(ncid, "overwrite", output%overwrite)
            call write_netcdf_attribute(ncid, "basename", trim(output%basename))

            call write_netcdf_attribute(ncid, "limit", time%limit)
            call write_netcdf_attribute(ncid, "initial", time%initial)
            call write_netcdf_attribute(ncid, "precise_stop", time%precise_stop)
            call write_netcdf_attribute(ncid, "alpha", time%alpha)

        end subroutine write_netcdf_options

end module options
