program main
  use mod_calendar 
  use mod_utilities, only : dp, state_vector_type
  use speedy_res_interface, only : read_era, era_data
  use mpires, only : startmpi, mpi_res, internal_state_vector
  use resdomain, only : initializedomain
  implicit none
  integer :: i
  !type(state_vector_type) :: internal_state_vector

  call startmpi()

  call initialize_calendar(calendar,1981,1,1,0)

  call initializedomain(mpi_res%numprocs,mpi_res%proc_num,0)
  
  call get_current_time_delta_hour(calendar,51000)
  call read_era(calendar%currentyear,calendar%currentyear)

  print *, shape(era_data%era_logp)
  !print *, 'logp',era_data%era_logp(:,:,3)
  where (era_data%eravariables(4,:,:,:,:) < 0.0)
    era_data%eravariables(4,:,:,:,:) = 0.0_dp
  end where
  era_data%eravariables(4,:,:,:,:) = era_data%eravariables(4,:,:,:,:)*1000.0_dp
  do i=1, 100
     call get_current_time_delta_hour(calendar,26305+i)!+i)

     print *, calendar

     if(i == 1) then
       internal_state_vector%variables3d = era_data%eravariables(:,:,:,:,3)
       internal_state_vector%logp = era_data%era_logp(:,:,3)
       print *,'logp',minval(internal_state_vector%logp),maxval(internal_state_vector%logp)
     endif 
      
     print *,'before speedy', minval(internal_state_vector%variables3d(4,:,:,:)), maxval(internal_state_vector%variables3d(4,:,:,:)) 

     internal_state_vector%era_hour = 1 !era hour of the month 1 = 00UTC of the first day of
     internal_state_vector%era_hour_plus_one  = 2!So I dont have to do calendar stuff in

     internal_state_vector%istart = 2
     internal_state_vector%era_start = 3

     internal_state_vector%iyear0 = calendar%currentyear
     internal_state_vector%imont0 = calendar%currentmonth
     internal_state_vector%iday = calendar%currentday
     internal_state_vector%ihour = calendar%currenthour
 
     call agcm_main(1,1,internal_state_vector)
     print *, 'after speedy',minval(internal_state_vector%variables3d(4,:,:,:)), maxval(internal_state_vector%variables3d(4,:,:,:)) 

     rewind(21)
     rewind(22)
     rewind(23)
     rewind(24)
     rewind(26) 
     rewind(30)
     rewind(2)
     rewind(10)
     rewind(11)
     rewind(13)
     rewind(15)
  enddo 
end program main

