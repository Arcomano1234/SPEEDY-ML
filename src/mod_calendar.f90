module mod_calendar
   !module to hold calendar stuff for hybrid and parallel reservoir 
   !calculations and this is completely independent from speedy's
   !internal calendar

   use mod_utilities, only : calendar_type

   implicit none 

   type(calendar_type) :: calendar

   contains
  
    subroutine initialize_calendar(datetime,startyear,startmonth,startday,starthour)
       type(calendar_type), intent(inout) :: datetime
       integer, intent(in)                :: startyear,startmonth,startday,starthour 

       datetime%startyear = startyear
       datetime%startmonth = startmonth
       datetime%startday = startday
       datetime%starthour = starthour
    end subroutine 

    subroutine get_current_time_delta_hour(datetime,hours_elapsed)
       !Takes an initialized calendar_type object and updates the current date
       !variables
       type(calendar_type), intent(inout)  :: datetime
       integer, intent(in)                 :: hours_elapsed

       !Local stuff
       integer               :: years_elasped, months_elapsed, days_elapsed, day_of_year
       integer               :: month, day_while_counter, leap_days
       integer               :: i
       integer, parameter    :: hours_in_year=8760 !Average number of hours in a year (includes leap year)

       integer, parameter    :: hours_in_a_day=24

       logical               :: is_leap_year
       !365-day calendar
       integer :: ncal365(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31,&
                                 30, 31 /)

       !Get the new year first
       years_elasped = hours_elapsed/hours_in_year

       datetime%currentyear = years_elasped + datetime%startyear

       !We need to get the leap years between start and now
       leap_days = 0
       do i=0,years_elasped - 1 
          call leap_year_check(datetime%startyear+i,is_leap_year)
          if(is_leap_year) then
            leap_days = leap_days + 1
          endif 
       end do 

       !Lets find what day we are in the year
       call leap_year_check(datetime%currentyear,is_leap_year)
    
       day_of_year = (mod(hours_elapsed,hours_in_year)/hours_in_a_day) - leap_days
       if(is_leap_year) then
          ncal365(2) = 29
       end if

       day_while_counter = day_of_year
       month = 1
       do while(day_while_counter > 0)
          day_while_counter = day_while_counter - ncal365(month)
          month = month + 1
       end do

       month = month - 1

       if(month <= 0) then
         month = 12
         datetime%currentyear = datetime%currentyear - 1
       endif 

       months_elapsed = month 

       !Set the new month
       datetime%currentmonth = months_elapsed

       !Get day of month
       days_elapsed = ncal365(month) + day_while_counter
       datetime%currentday = days_elapsed

       !Get hour of day 
       datetime%currenthour = mod(hours_elapsed,hours_in_a_day)

       return
    end subroutine

    subroutine leap_year_check(year,is_leap_year)
       integer, intent(in)  :: year
       logical, intent(out) :: is_leap_year

       if((mod(year,4) == 0).and.(mod(year,100) /= 0)) then
         is_leap_year = .True.
       else if(mod(year,400) == 0) then
         is_leap_year = .True.
       else
         is_leap_year = .False.
       endif
       return
    end subroutine

    subroutine numof_hours(startyear,endyear,numofhours)
      !Get the number of hours assumes you start of jan 1 of start year and end
      !dec 31 of
      !endyear
      integer, intent(in)   :: startyear, endyear
      integer, intent(out)  :: numofhours

      integer               :: years_elapsed, i 
      integer, parameter    :: hours_in_year=8760 !Number of hours in a 365 day year
      integer, parameter    :: hours_in_year_leap_year = 8784

      logical               :: is_leap_year 

      years_elapsed = endyear-startyear
      numofhours = 0 
      do i=0,years_elapsed
         call leap_year_check(startyear+i,is_leap_year) 
         if(is_leap_year) then
           numofhours = numofhours + hours_in_year_leap_year
         else 
            numofhours = numofhours + hours_in_year
         endif 
      enddo 
    end subroutine

    subroutine numof_hours_into_year(year,month,day,hour,numofhours)
      !Get the number of hours you are into the year assumes you start of jan 1 of year 
      integer, intent(in)   :: year,month,day,hour
      integer, intent(out)  :: numofhours

      integer               :: months_elapsed, i

      logical               :: is_leap_year

      integer :: ncal365(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31,&
                                  30, 31 /)

      integer :: ncal_leap(12) = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31,&
                                    30, 31 /)

      months_elapsed = month 
      
      numofhours = 0
      call leap_year_check(year,is_leap_year)
      !Loop through months
      if (months_elapsed > 1) then
         do i=1,months_elapsed-1
            if(is_leap_year) then 
               numofhours = numofhours + 24*ncal_leap(i)
            else
               numofhours = numofhours + 24*ncal365(i) 
            endif 
         enddo 
      endif 

      !Loop through days 
      if(day > 1) then 
        do i=1,day-1
           numofhours = numofhours + 24
        enddo
      endif
       
      numofhours =  numofhours+hour

      if(numofhours == 0) then
        numofhours = 1
      endif 
    end subroutine

    subroutine time_delta_between_two_dates(start_year,start_month,start_day,start_hour,end_year,end_month,end_day,end_hour,numofhours)
      integer, intent(in)   :: start_year,start_month,start_day,start_hour
      integer, intent(in)   :: end_year,end_month,end_day,end_hour
      integer, intent(out)  :: numofhours 

      !Local variables
      integer :: hours_into_start_year, hours_into_end_year
      integer :: hours_between_years

      if(start_year /= end_year) then 
        call numof_hours(start_year,end_year-1,hours_between_years)

        call numof_hours_into_year(start_year,start_month,start_day,start_hour,hours_into_start_year)

        call numof_hours_into_year(end_year,end_month,end_day,end_hour,hours_into_end_year)

        numofhours = hours_between_years + hours_into_end_year - hours_into_start_year

      else 
        call numof_hours_into_year(start_year,start_month,start_day,start_hour,hours_into_start_year)

        call numof_hours_into_year(end_year,end_month,end_day,end_hour,hours_into_end_year)

        numofhours = hours_into_end_year - hours_into_start_year

      endif  
     
    end subroutine 

    subroutine time_delta_between_two_dates_datetime_type(datatime1,datetime2,timedelta)
      type(calendar_type), intent(inout) :: datatime1,datetime2

      integer, intent(out)               :: timedelta

      call time_delta_between_two_dates(datatime1%currentyear,datatime1%currentmonth,datatime1%currentday,datatime1%currenthour,datetime2%currentyear,datetime2%currentmonth,datetime2%currentday,datetime2%currenthour,timedelta)
      
    end subroutine 
end module mod_calendar
