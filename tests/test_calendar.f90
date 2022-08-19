program main

 use mod_calendar

 implicit none 

 integer :: discardlength, synclength, traininglength, prediction_length, prediction_num, cycle_length, timestep, i, num_predictions, j

 integer, allocatable :: prediction_markers(:)

 discardlength =  150
 synclength = 24*4
 traininglength = 188280
 prediction_length = 504
 timestep = 6
 num_predictions = 100

 call initialize_calendar(calendar,1990,1,1,0)
 call get_current_time_delta_hour(calendar,discardlength+traininglength+synclength)

 cycle_length = synclength

 allocate(prediction_markers(num_predictions))

 do i=0,num_predictions-1
    prediction_markers(i+1) = cycle_length*i
 enddo

 print *, 'prediction_markers',prediction_markers

 do i=1,num_predictions
    do j=1, prediction_length
       call get_current_time_delta_hour(calendar,traininglength+synclength+prediction_markers(i))
       print *,'prediction file date',calendar%currentmonth,calendar%currentday,calendar%currentyear,calendar%currenthour
    enddo 
 enddo  
end program main
