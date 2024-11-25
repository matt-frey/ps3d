module model_manager
    use options, only : output              &
                      , read_config_file
    use model_factory, only : create_model, model_t, model_info_t
    use constants, only : zero
    use mpi_timer
    use fields
    use field_netcdf, only : field_io_timer
    use field_diagnostics_netcdf, only : field_stats_io_timer
    use inversion_mod, only : vor2vel_timer &
                            , vtend_timer   &
                            , vor2vel
    use fields_derived, only : pres_timer   &
                             , delta_timer
    use inversion_utils, only : init_inversion, finalise_inversion
    use advance_mod, only : advance             &
                          , advance_timer
    use utils, only : write_last_step, setup_output_files,   &
                      setup_domain_and_parameters            &
                    , setup_fields
#ifdef ENABLE_BALANCE
    use field_balance, only : initialise_balance, finalise_balance
#endif
    implicit none

    private

    type(model_t) :: model
    integer       :: ps_timer

    public :: pre_run, run, post_run

contains

    subroutine pre_run
        type(model_info_t) :: model_info

        call register_timer('ps', ps_timer)
        call register_timer('field I/O', field_io_timer)
        call register_timer('field diagnostics I/O', field_stats_io_timer)
        call register_timer('vor2vel', vor2vel_timer)
        call register_timer('vorticity tendency', vtend_timer)
        call register_timer('advance', advance_timer)
        call register_timer('pressure calculation', pres_timer)
        call register_timer('horizontal divergence calculation', delta_timer)

        call start_timer(ps_timer)

        ! parse the config file
        call read_config_file

        ! read domain dimensions
        call setup_domain_and_parameters

        call init_inversion

#ifdef ENABLE_BALANCE
        if (output%l_balanced) then
            call initialise_balance
        endif
#endif

        call setup_fields

        call setup_output_files

        model = create_model(model_info)

    end subroutine

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine run
        use options, only : time
        double precision :: t = zero    ! current time

        t = time%initial

        call start_timer(advance_timer)
        do while (t < time%limit)
            call advance(model%stepper, t)
        enddo
        call stop_timer(advance_timer)

        ! write final step (we only write if we really advanced in time)
        if (t > time%initial) then
            call write_last_step(t)
        endif

    end subroutine run

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine post_run
        use options, only : output

        call finalise_inversion

#ifdef ENABLE_BALANCE
        if (output%l_balanced) then
            call finalise_balance
        endif
#endif

        call stop_timer(ps_timer)
        call write_time_to_csv(output%basename)
        call print_timer

    end subroutine

end module model_manager
