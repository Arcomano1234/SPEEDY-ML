subroutine stloop(istep)
    !   subroutine stloop (istep)
    !
    !   Purpose: Perform a series of time steps calling 
    !            post-processing/output routines at selected steps
    !   Input/output : istep = time step index
    !   Updated common block : lflag2
      
    use mod_lflags, only: lradsw, lrandf
    use mod_tsteps
    use mod_date, only: ihour, newdate, iyear
    use mod_dynvar

    use speedy_res_interface, only : getspeedyvariable
    use mod_reservoir, only : global_time_step

    implicit none

    integer, intent(inout) :: istep
    integer :: iitest = 0, j, jj
    integer :: window_size  
 
    window_size = 24/global_time_step

    ! Break up each day into 24 1-hour windows
    do jj = 1,window_size
        ! Each 1-hour window has nsteps/24 actual timesteps
        do j = 1, nsteps/window_size
            !if (iitest == 0) print*, 'stloop: calling step ', istep

            !Keep track of the timestep and record state vector 
            currentstep = currentstep + 1
            !if(era_start.ne.3) then
            !   call getspeedyvariable()
            !endif 

            print *, 'made it here', currentstep
            ! Set logical flags
            lradsw = (mod(istep,nstrad) == 1)
            lrandf = ((istep <= nstrdf) .or. (nstrdf < 0))
    
            ! Perform one leapfrog time step
            call step(2, 2, delt2, alph, rob, wil)   
    
            ! Do diagnostic, post-processing and I/O tasks 
            call diagns(2, istep)
    
            if (ihout .eqv. .false.) then
                if (mod(istep, nstppr) == 0) call tminc
                if (nstout > 0 .and. mod(istep, nstout) == 0) call tmout(1)
            end if
    
            if(era_start.ne.3) then
               call getspeedyvariable()
            endif
            istep = istep + 1

            if(mod(j,window_size) == 0) then
              ! Increment hour timer (takes values of 0, 6, 12 or 18)
              print *, 'ihour before',ihour

              ihour = mod(ihour + 1 , 24)

              print *, 'ihour after',ihour
              ! If it's a new day...
              if (ihour .eq. 0) then
                 ! Compute new date
                 call newdate(1)
              end if
            endif

        end do

        ! Increment hour timer (takes values of 0, 6, 12 or 18)
        !ihour = mod(ihour + res%timestep , 24)

        ! If it's a new day...
        !if (ihour .eq. 0) then
            ! Compute new date
        !    call newdate(1)
        !end if

        if(onehr_run) then
           !call restart(2)
           call iogrid(69) 
           
           print *,'normal end with 1-hr fcst (yeahhhhhhh!!!!)'
           stop
        endif 
        if(onehr_hybrid) then
         print *, 'exiting speedy via onehr_hybrid'
         call iogrid(31)
         exit 
        endif    
    end do
end
