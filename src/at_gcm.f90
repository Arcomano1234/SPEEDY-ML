module speedy_main

contains

subroutine agcm_main(input_ndays,start_from_file,internal_state_vector)
    use mod_tsteps, only: ndaysl, ihout, nmonts, sixhrrun
    use mod_utilities, only : state_vector_type

    implicit none

    ! program : agcm_main
    integer, intent(in) :: input_ndays, start_from_file
    type(state_vector_type), intent(inout), optional :: internal_state_vector

    ! experiment identifier
    character(len=3) :: cexp = 'exp'
    integer :: jday, ndays

    ! 1. initialization
    ! ndays = no. of integration days, set by agcm_init
    if(present(internal_state_vector)) then
      print *, 'starting with interal_state_vector'
      call agcm_init(cexp, 0, 0, 0, ndays,start_from_file,internal_state_vector)
    else 
      call agcm_init(cexp, 0, 0, 0, ndays,start_from_file,internal_state_vector)
    endif 

    ! 2. do loop over total no. of integration days
    if(present(internal_state_vector)) then
      ndays = input_ndays
    endif 

    ndays = input_ndays
    print *, 'integration length in days: ', ndays    

    if(present(internal_state_vector)) then
      if(internal_state_vector%is_safe_to_run_speedy) then 
        do jday = 1, ndays
           ! 2.2 run atmospheric model for 1 day
           call agcm_1day(jday, cexp)

           ! 2.1 exchange data with coupler
           call agcm_to_coupler(jday)
           call coupler_to_agcm(jday)
  
        enddo
      endif 
    else
      do jday = 1, ndays
         ! 2.2 run atmospheric model for 1 day
         call agcm_1day(jday, cexp)

         ! 2.1 exchange data with coupler
         call agcm_to_coupler(jday)
         call coupler_to_agcm(jday)

      enddo
    endif 
      
    ! Restart dataset is only written at the end
    !call restart(2)
end subroutine 

subroutine agcm_1day(jday, cexp)
    ! subroutine agcm_1day (jday)
    !
    ! perform atm. model integration for 1 day, 
    ! post-proc. and i/o at selected times 

    use mod_tsteps, only: nsteps, idout, nstout, ihout
    use mod_date, only: iyear, imonth, iday, ndaytot, newdate

    implicit none

    integer, intent(in) :: jday
    character(len=3), intent(in) :: cexp
    integer :: istep

    if (iday == 1) print *, ' start of year/month = ', iyear, imonth

    istep = 1 + (jday - 1) * nsteps

    ! 1. set forcing terms according to date
    call fordate(1)

    ! 2. set daily-average flux arrays to zero
    call dmflux(0)

    ! 3. integrate the atmospheric model for 1 day
    call stloop(istep)

    ! 4. write daily-mean output
    !call dmout(idout)

    ! 5. write time-mean output files and restart file at the end of selected
    ! months
    !if (iday == 1) then
        ! write monthly-mean output for previous month
    !    if (ihout .eqv. .false.) then
    !        if (nstout < 0) call tmout(1)
    !    end if
        
        ! open new output files at the beginning of each year
    !    if (imonth == 1 .and. jday < ndaytot .and. (ihout .eqv. .false.)) call setgrd(1, cexp)
    !endif
end subroutine

end module speedy_main
