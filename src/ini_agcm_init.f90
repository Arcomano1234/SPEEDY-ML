subroutine agcm_init(cexp, inidate, ntimes, irstart, ndays, start_from_file, internal_state_vector)
    !   subroutine agcm_init (cexp,inidate,ntimes,irstart,
    !  &                      ndays)
    !
    !   purpose: initialization of atmos. model and coupling interface 
    !

    use mod_cpl_flags, only: icsea, isstan
    use mod_tsteps
    use mod_date, only: newdate, ndaytot, iyear, imonth, iday, ihour
   
    use mod_utilities, only : state_vector_type 
    implicit none

    ! input (reset by input/include files if inidate = 0):
    character(len=3), intent(inout) :: cexp        ! experiment identifier
    integer, intent(in) :: inidate     ! initial date yyyymm
    integer, intent(in) :: ntimes      ! integr. length in months (< 0) or days (> 0)
    integer, intent(in) :: irstart     ! restart flag: 0 = no, > 0 = yes
    integer, intent(in) :: start_from_file !0 starts from fort.2  if its 1 then uses internal to program type state_vector_type

    type(state_vector_type), intent(inout), optional :: internal_state_vector
    ! output:
    integer, intent(inout) :: ndays       ! total no. of integration days

   

    print *, ' hallo from speedy_agcm'

    ! 1. set run initial time, duration, time-stepping and coupling options

    if(start_from_file == 0) then
      read (2,*) istart !from rest=0, from restartfile=1, from era=2
      read (2,*) era_start !start from grid initial condition=0, Start from grid era_5 re_analysis=1, regrid era=2
      read (2,'(a)') era_file 
      era_file = trim(era_file)

      read (2,*) era_hour !era hour of the month 1 = 00UTC of the first day of the month
      read (2,*) era_hour_plus_one !So I dont have to do calendar stuff in fortran

      ! Read date from fort.2 file
      read (2,*) iyear0
      read (2,*) imont0
      read (2,*) iday
      read (2,*) ihour
    else if(start_from_file == 1) then
      if(present(internal_state_vector)) then
        istart = internal_state_vector%istart
        era_start = internal_state_vector%era_start
        era_file = internal_state_vector%era_file
        era_file = trim(era_file)
  
        iyear0 = internal_state_vector%iyear0
        imont0 = internal_state_vector%imont0
        iday = internal_state_vector%iday
        ihour = internal_state_vector%ihour
      else 
        print *, 'something went horribly wrong check internal_state_vector killing the program' 
        stop
      endif 
    endif 

    iyear = iyear0
    imonth = imont0

    call newdate(0)

    print *, 'start date ', iyear, imonth, iday, ihour

    isst0 = (iyear0 - issty0) * 12 + imont0

    ndays = ndaytot

    ! check consistency of coupling and prescribed SST anomaly flags
    if (icsea >= 4) isstan = 1

    ! 2. initialization of atmospheric model constants and variables 
    call ini_atm(cexp)

    if(present(internal_state_vector)) then
      if(internal_state_vector%is_safe_to_run_speedy) then
         ! 3. initialization of coupled modules (land, sea, ice)
          call ini_coupler(istart)
         
         ! 4. set up the forcing fields for the first time step
         call fordate(0)

         ! 5. do the initial (2nd-order) time step, initialize the semi-impl. scheme
         call stepone
      endif 
     else 
       ! 3. initialization of coupled modules (land, sea, ice)
       call ini_coupler(istart)

       ! 4. set up the forcing fields for the first time step
       call fordate(0)

       ! 5. do the initial (2nd-order) time step, initialize the semi-impl.
       ! scheme
       call stepone
     endif  
end subroutine
