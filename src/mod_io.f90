module mod_io
   use iso_fortran_env

   use netcdf
   use mod_utilities, only : dp, speedygridnum, xgrid, ygrid, zgrid, &
                             speedylat, model_parameters_type, int_32, &
                             e_constant

   implicit none
   
   integer :: error  
   
   interface write_netcdf
     module procedure write_netcdf_4d
     module procedure write_netcdf_4d_logp
     module procedure write_netcdf_4d_multi_2d
     module procedure write_netcdf_4d_multi_2d_no_sst
   end interface

   contains 
     
     function file_exists(filename) result(result_bool)
       character(len=*), intent(in) :: filename
       logical                      :: result_bool

       !Very basic way to check if a file exists
       inquire(file=trim(filename),exist=result_bool)
       
     end function
    
     !-----------NETCDF SECTION ------------------!
     subroutine write_netcdf_4d(model_parameters,dat,timestep,filename)
       use netcdf

       use stringtype, only : string 
       use mod_utilities, only : main_type

       real(kind=dp), intent(in)    :: dat(:,:,:,:)
       integer, intent(in)          :: timestep
       character(len=*), intent(in) :: filename
       type(model_parameters_type)         :: model_parameters

       integer, parameter :: numdims=4
       integer :: dimsx,dimsy,dimsz    
       
       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       real(kind=dp) :: lon(xgrid), lat(ygrid)

       real(kind=dp), allocatable :: copy(:,:,:,:)

       type(string) :: units(model_parameters%full_predictvars)
       type(string) :: varname(model_parameters%full_predictvars)

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'
       
       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'

       !units(4)%str = 'log(surfacepressure)'
       !varname(4)%str = 'logp'

       !copy data
       allocate(copy,source=dat)

       dimsx = size(dat,2)
       dimsy = size(dat,3)
       dimsz = size(dat,4)

       varcount = [integer:: dimsx, dimsy, dimsz, 1]
       start = [integer:: 1, 1, 1, timestep ]

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = (/-87.159, -83.479, -79.777,  -76.070, -72.362, -68.652, &
               -64.942, -61.232, -57.521, -53.810, -50.099, -46.389, -42.678, -38.967, &
               -35.256, -31.545 , -27.833, -24.122, -20.411, -16.700, -12.989, -9.278, &
               -5.567, -1.856, 1.856, 5.567, 9.278, 12.989, 16.700, 20.411, &
               24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389, 50.099, &
               53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, &
               83.479, 87.159/)

       if(timestep.eq.1) then
           call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER, ncid=file_id)) 

             ! define the dimensions
           call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
           call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))
           call nc_check(nf90_def_dim(file_id, 'Sigma_Level', dimsz, zdim_id))
           call nc_check(nf90_def_dim(file_id, 'Timestep', NF90_UNLIMITED, timedim_id))

           !Assign lat and lon ids and units
           call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
           call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))
 
           call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
           call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
           arrdims = (/ xdim_id, ydim_id, zdim_id, timedim_id /) 
           do i=1, model_parameters%full_predictvars
 
              call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims,array_id))
              ! ...and assign units to them as an attribute

              call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

              call nc_check(nf90_enddef(file_id))

              !Write out the values
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start, count=varcount))

              call nc_check(nf90_put_var(file_id, xvar_id, lon))

              call nc_check(nf90_put_var(file_id, yvar_id, lat))
              
              call nc_check(nf90_redef(file_id))
           enddo 
           ! close; done
           call nc_check(nf90_close(file_id))
        else
           call nc_check(nf90_open(filename,nf90_write,file_id))
           do i=1,model_parameters%full_predictvars
              call nc_check(nf90_inq_varid(file_id,varname(i)%str,array_id))
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start,count=varcount))
           enddo 
           call nc_check(nf90_close(file_id))
        endif
        return
     end subroutine write_netcdf_4d
      
     subroutine write_netcdf_4d_logp(model_parameters,grid4d,grid3d,timestep,filename)
       use, intrinsic :: ieee_arithmetic

       use stringtype, only : string 
       use mod_utilities, only : e_constant

       real(kind=dp), intent(in)    :: grid4d(:,:,:,:)
       real(kind=dp), intent(in)    :: grid3d(:,:) !Yes I know bad name sorry

       integer, intent(in)          :: timestep

       character(len=*), intent(in) :: filename

       type(model_parameters_type), intent(in) :: model_parameters 

       integer, parameter :: numdims4d=4, numdims3d=3
       integer :: dimsx,dimsy,dimsz    
       
       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d),start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims4d(numdims4d),arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(xgrid), lat(ygrid)

       real(kind=dp), allocatable :: copy(:,:,:,:)

       type(string) :: units(5)
       type(string) :: varname(5)

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'
       
       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'

       units(4)%str = 'g/kg'
       varname(4)%str = 'Specific-Humidity'
      
       units(5)%str = 'log(surfacepressure)'
       varname(5)%str = 'logp'

       !copy data
       allocate(copy,source=grid4d)

      
       !if(res%specific_humidity_log_bool) then    
       !  copy(4,:,:,:) = e_constant**copy(4,:,:,:) - res%specific_humidity_epsilon
       !endif 

       dimsx = size(grid4d,2)
       dimsy = size(grid4d,3)
       dimsz = size(grid4d,4)

       varcount4d =  [integer:: dimsx, dimsy, dimsz, 1]
       start4d =  [integer:: 1, 1, 1, timestep]
    
       varcount3d =  [integer:: dimsx, dimsy, 1]
       start3d =  [integer:: 1, 1, timestep]

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = (/-87.159, -83.479, -79.777,  -76.070, -72.362, -68.652, &
               -64.942, -61.232, -57.521, -53.810, -50.099, -46.389, -42.678, -38.967, &
               -35.256, -31.545 , -27.833, -24.122, -20.411, -16.700, -12.989, -9.278, &
               -5.567, -1.856, 1.856, 5.567, 9.278, 12.989, 16.700, 20.411, &
               24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389, 50.099, &
               53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, &
               83.479, 87.159/)

       if(timestep.eq.1) then
           call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER, ncid=file_id)) 

             ! define the dimensions
           call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
           call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))
           call nc_check(nf90_def_dim(file_id, 'Sigma_Level', dimsz, zdim_id))
           call nc_check(nf90_def_dim(file_id, 'Timestep', NF90_UNLIMITED, timedim_id))

           !Assign lat and lon ids and units
           call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
           call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))
 
           call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
           call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
           arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /) 
           arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

           do i=1, model_parameters%full_predictvars
              call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims4d,array_id))
              ! ...and assign units to them as an attribute

              call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

              call nc_check(nf90_enddef(file_id))

              !Write out the values
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d, count=varcount4d))

              call nc_check(nf90_put_var(file_id, xvar_id, lon))

              call nc_check(nf90_put_var(file_id, yvar_id, lat))
              
              call nc_check(nf90_redef(file_id))
           enddo 
           !Lets do logp
           call nc_check(nf90_def_var(file_id,varname(5)%str,NF90_REAL,arrdims3d,array_id))
           
           call nc_check(nf90_put_att(file_id, array_id, "units", units(5)%str))

           call nc_check(nf90_enddef(file_id))

           !Write out the values
           call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d, count=varcount3d))

           call nc_check(nf90_put_var(file_id, xvar_id, lon))

           call nc_check(nf90_put_var(file_id, yvar_id, lat))
           ! close; done
           call nc_check(nf90_close(file_id))
        else
           call nc_check(nf90_open(filename,nf90_write,file_id))
           do i=1,model_parameters%full_predictvars
              call nc_check(nf90_inq_varid(file_id,varname(i)%str,array_id))
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d,count=varcount4d))
           enddo 

           call nc_check(nf90_inq_varid(file_id,varname(5)%str,array_id))
           call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d,count=varcount3d))

           call nc_check(nf90_close(file_id))
        endif
        return
     end subroutine write_netcdf_4d_logp

     subroutine write_netcdf_4d_multi_2d(model_parameters,grid4d,grid3d,timestep,filename,ocean_model)
       use, intrinsic :: ieee_arithmetic

       use stringtype, only : string
       use mod_utilities, only : unstandardize_data, reservoir_type, e_constant

       real(kind=dp), intent(in)    :: grid4d(:,:,:,:)
       real(kind=dp), intent(in)    :: grid3d(:,:,:) !Yes I know bad name sorry

       integer, intent(in)          :: timestep

       character(len=*), intent(in) :: filename

       logical :: ocean_model

       type(model_parameters_type), intent(in) :: model_parameters

       integer, parameter :: numdims4d=4, numdims3d=3
       integer :: dimsx,dimsy,dimsz

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d),start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims4d(numdims4d),arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(xgrid), lat(ygrid)

       real(kind=dp), allocatable :: copy(:,:,:,:), copy3d(:,:,:)

       type(string) :: units(7)
       type(string) :: varname(7)

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'

       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'

       units(4)%str = 'g/kg'
       varname(4)%str = 'Specific-Humidity'
       
       units(5)%str = 'log(surfacepressure)'
       varname(5)%str = 'logp'

       units(6)%str = 'Kelvin'
       varname(6)%str = 'SST'
 
       units(7)%str = 'in/6hr'
       varname(7)%str = 'p6hr'

       !copy data
       allocate(copy,source=grid4d)
       allocate(copy3d, source=grid3d)

       !unstandarize the data
       !call unstandardize_data(res,copy,res%mean,res%std)

       !if(res%specific_humidity_log_bool) then
       !  copy(4,:,:,:) = e_constant**copy(4,:,:,:) -
       !  res%specific_humidity_epsilon
       !endif

       dimsx = size(grid4d,2)
       dimsy = size(grid4d,3)
       dimsz = size(grid4d,4)

       varcount4d = (/ dimsx, dimsy, dimsz, 1 /)
       start4d = (/ 1, 1, 1, timestep /)

       varcount3d = (/ dimsx, dimsy, 1 /)
       start3d = (/ 1, 1, timestep /)

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = speedylat

       if(timestep.eq.1) then
           call nc_check(nf90_create(path=filename,cmode=ior(nf90_clobber,nf90_64bit_offset),ncid=file_id))

             ! define the dimensions
           call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
           call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))
           call nc_check(nf90_def_dim(file_id, 'Sigma_Level', dimsz, zdim_id))
           call nc_check(nf90_def_dim(file_id, 'Timestep', NF90_UNLIMITED, timedim_id))

           !Assign lat and lon ids and units
           call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
           call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

           call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
           call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
           arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /)
           arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

           do i=1, model_parameters%full_predictvars
              call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims4d,array_id))
              ! ...and assign units to them as an attribute

              call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

              call nc_check(nf90_enddef(file_id))

              !Write out the values
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d, count=varcount4d))

              call nc_check(nf90_put_var(file_id, xvar_id, lon))

              call nc_check(nf90_put_var(file_id, yvar_id, lat))

              call nc_check(nf90_redef(file_id))
           enddo
           !Lets do logp

           call nc_check(nf90_def_var(file_id,varname(5)%str,NF90_REAL,arrdims3d,array_id))

           call nc_check(nf90_put_att(file_id, array_id, "units", units(5)%str))

           call nc_check(nf90_enddef(file_id))

           !Write out the values
           call nc_check(nf90_put_var(file_id, array_id, grid3d(1,:,:),start=start3d, count=varcount3d))

           call nc_check(nf90_put_var(file_id, xvar_id, lon))

           call nc_check(nf90_put_var(file_id, yvar_id, lat))

           if(model_parameters%slab_ocean_model_bool) then
             call nc_check(nf90_redef(file_id)) 

             !Lets do SST
             call nc_check(nf90_def_var(file_id,varname(6)%str,NF90_REAL,arrdims3d,array_id))

             call nc_check(nf90_put_att(file_id, array_id, "units", units(6)%str))

             call nc_check(nf90_enddef(file_id))

             !Write out the values
             call nc_check(nf90_put_var(file_id, array_id, grid3d(2,:,:),start=start3d, count=varcount3d))

             call nc_check(nf90_put_var(file_id, xvar_id, lon))

             call nc_check(nf90_put_var(file_id, yvar_id, lat))
           endif

          
           if(model_parameters%precip_bool) then
             call nc_check(nf90_redef(file_id)) 

             !Lets do 6 hour precip
             call nc_check(nf90_def_var(file_id,varname(7)%str,NF90_REAL,arrdims3d,array_id))

             call nc_check(nf90_put_att(file_id, array_id, "units", units(7)%str))

             call nc_check(nf90_enddef(file_id))

             !Write out the values
             copy3d(3,:,:) = model_parameters%precip_epsilon * ( e_constant**grid3d(3,:,:) - 1)
             copy3d(3,:,:) = copy3d(3,:,:) * 25.4 * 1000.0!* 39.3701

             call nc_check(nf90_put_var(file_id, array_id, copy3d(3,:,:), start=start3d, count=varcount3d))

             call nc_check(nf90_put_var(file_id, xvar_id, lon))

             call nc_check(nf90_put_var(file_id, yvar_id, lat))
           endif

           ! close; done
           call nc_check(nf90_close(file_id))
        else
           call nc_check(nf90_open(filename,nf90_write,file_id))
           do i=1, model_parameters%full_predictvars
              call nc_check(nf90_inq_varid(file_id,varname(i)%str,array_id))
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d,count=varcount4d))
           enddo

           call nc_check(nf90_inq_varid(file_id,varname(5)%str,array_id))
           call nc_check(nf90_put_var(file_id, array_id, grid3d(1,:,:),start=start3d,count=varcount3d))

           if(ocean_model) then
             call nc_check(nf90_inq_varid(file_id,varname(6)%str,array_id)) 
             call nc_check(nf90_put_var(file_id, array_id, grid3d(2,:,:), start=start3d, count=varcount3d))
           endif

           if(model_parameters%precip_bool) then
             copy3d(3,:,:) = model_parameters%precip_epsilon * (e_constant**grid3d(3,:,:) - 1) 
             copy3d(3,:,:) = copy3d(3,:,:) * 39.3701

             call nc_check(nf90_inq_varid(file_id,varname(7)%str,array_id))
             call nc_check(nf90_put_var(file_id, array_id, copy3d(3,:,:), start=start3d, count=varcount3d))
           endif

           call nc_check(nf90_close(file_id))
        endif
        return
     end subroutine

     subroutine write_netcdf_4d_multi_2d_no_sst(model_parameters,grid4d,grid3d,timestep,filename)
       use, intrinsic :: ieee_arithmetic

       use stringtype, only : string
       use mod_utilities, only : unstandardize_data, reservoir_type, e_constant

       real(kind=dp), intent(in)    :: grid4d(:,:,:,:)
       real(kind=dp), intent(in)    :: grid3d(:,:,:) !Yes I know bad name sorry

       integer, intent(in)          :: timestep

       character(len=*), intent(in) :: filename

       logical :: ocean_model

       type(model_parameters_type), intent(in) :: model_parameters

       integer, parameter :: numdims4d=4, numdims3d=3
       integer :: dimsx,dimsy,dimsz

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d),start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims4d(numdims4d),arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(xgrid), lat(ygrid)

       real(kind=dp), allocatable :: copy(:,:,:,:), copy3d(:,:,:)

       type(string) :: units(6)
       type(string) :: varname(6)

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'

       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'

       units(4)%str = 'g/kg'
       varname(4)%str = 'Specific-Humidity'
       
       units(5)%str = 'log(surfacepressure)'
       varname(5)%str = 'logp'
 
       units(6)%str = 'in/6hr'
       varname(6)%str = 'p6hr'

       !copy data
       allocate(copy,source=grid4d)
       allocate(copy3d,source=grid3d)

       !unstandarize the data
       !call unstandardize_data(res,copy,res%mean,res%std)

       !if(res%specific_humidity_log_bool) then
       !  copy(4,:,:,:) = e_constant**copy(4,:,:,:) -
       !  res%specific_humidity_epsilon
       !endif

       dimsx = size(grid4d,2)
       dimsy = size(grid4d,3)
       dimsz = size(grid4d,4)

       varcount4d = (/ dimsx, dimsy, dimsz, 1 /)
       start4d = (/ 1, 1, 1, timestep /)

       varcount3d = (/ dimsx, dimsy, 1 /)
       start3d = (/ 1, 1, timestep /)

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = speedylat

       if(timestep.eq.1) then
           call nc_check(nf90_create(path=filename,cmode=ior(nf90_clobber,nf90_64bit_offset),ncid=file_id))

             ! define the dimensions
           call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
           call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))
           call nc_check(nf90_def_dim(file_id, 'Sigma_Level', dimsz, zdim_id))
           call nc_check(nf90_def_dim(file_id, 'Timestep', NF90_UNLIMITED, timedim_id))

           !Assign lat and lon ids and units
           call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
           call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

           call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
           call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
           arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /)
           arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

           do i=1, model_parameters%full_predictvars
              call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims4d,array_id))
              ! ...and assign units to them as an attribute

              call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

              call nc_check(nf90_enddef(file_id))

              !Write out the values
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d, count=varcount4d))

              call nc_check(nf90_put_var(file_id, xvar_id, lon))

              call nc_check(nf90_put_var(file_id, yvar_id, lat))

              call nc_check(nf90_redef(file_id))
           enddo
           !Lets do logp
           call nc_check(nf90_def_var(file_id,varname(5)%str,NF90_REAL,arrdims3d,array_id))

           call nc_check(nf90_put_att(file_id, array_id, "units", units(5)%str))

           call nc_check(nf90_enddef(file_id))

           !Write out the values
           call nc_check(nf90_put_var(file_id, array_id, grid3d(1,:,:),start=start3d, count=varcount3d))

           call nc_check(nf90_put_var(file_id, xvar_id, lon))

           call nc_check(nf90_put_var(file_id, yvar_id, lat))

           call nc_check(nf90_redef(file_id))

           if(model_parameters%precip_bool) then
             !Lets do 6 hour precip
             call nc_check(nf90_def_var(file_id,varname(6)%str,NF90_REAL,arrdims3d,array_id))

             call nc_check(nf90_put_att(file_id, array_id, "units", units(6)%str))

             call nc_check(nf90_enddef(file_id))

             !Write out the values
             copy3d(2,:,:) = model_parameters%precip_epsilon * (e_constant**grid3d(2,:,:) - 1)
             copy3d(2,:,:) = copy3d(2,:,:) * 39.3701

             call nc_check(nf90_put_var(file_id, array_id, copy3d(2,:,:),start=start3d, count=varcount3d))

             call nc_check(nf90_put_var(file_id, xvar_id, lon))

             call nc_check(nf90_put_var(file_id, yvar_id, lat))
           endif

           ! close; done
           call nc_check(nf90_close(file_id))
        else
           call nc_check(nf90_open(filename,nf90_write,file_id))
           do i=1, model_parameters%full_predictvars
              call nc_check(nf90_inq_varid(file_id,varname(i)%str,array_id))
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d,count=varcount4d))
           enddo

           call nc_check(nf90_inq_varid(file_id,varname(5)%str,array_id))
           call nc_check(nf90_put_var(file_id, array_id, grid3d(1,:,:),start=start3d,count=varcount3d))

           if(model_parameters%precip_bool) then
             copy3d(2,:,:) = model_parameters%precip_epsilon * (e_constant**grid3d(2,:,:) - 1) 
             copy3d(2,:,:) = copy3d(2,:,:) * 39.3701

             call nc_check(nf90_inq_varid(file_id,varname(6)%str,array_id))
             call nc_check(nf90_put_var(file_id, array_id, copy3d(2,:,:),start=start3d, count=varcount3d))
           endif

           call nc_check(nf90_close(file_id))
        endif
        return
     end subroutine

     subroutine write_netcdf_4d_multi_2d_sst_only(model_parameters,grid4d,grid3d,timestep,filename)
       use, intrinsic :: ieee_arithmetic

       use stringtype, only : string
       use mod_utilities, only : unstandardize_data, reservoir_type, e_constant

       real(kind=dp), intent(in)    :: grid4d(:,:,:,:)
       real(kind=dp), intent(in)    :: grid3d(:,:,:) !Yes I know bad name sorry

       integer, intent(in)          :: timestep

       character(len=*), intent(in) :: filename

       logical :: ocean_model

       type(model_parameters_type), intent(in) :: model_parameters

       integer, parameter :: numdims4d=4, numdims3d=3
       integer :: dimsx,dimsy,dimsz

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d),start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims4d(numdims4d),arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(xgrid), lat(ygrid)

       real(kind=dp), allocatable :: copy(:,:,:,:), copy3d(:,:,:)

       type(string) :: units(6)
       type(string) :: varname(6)

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'

       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'

       units(4)%str = 'g/kg'
       varname(4)%str = 'Specific-Humidity'

       units(5)%str = 'log(surfacepressure)'
       varname(5)%str = 'logp'

       units(6)%str = 'Kelvin'
       varname(6)%str = 'sst'

       !copy data
       allocate(copy,source=grid4d)
       allocate(copy3d,source=grid3d)

       !unstandarize the data
       !call unstandardize_data(res,copy,res%mean,res%std)

       !if(res%specific_humidity_log_bool) then
       !  copy(4,:,:,:) = e_constant**copy(4,:,:,:) -
       !  res%specific_humidity_epsilon
       !endif

       dimsx = size(grid4d,2)
       dimsy = size(grid4d,3)
       dimsz = size(grid4d,4)

       varcount4d = (/ dimsx, dimsy, dimsz, 1 /)
       start4d = (/ 1, 1, 1, timestep /)

       varcount3d = (/ dimsx, dimsy, 1 /)
       start3d = (/ 1, 1, timestep /)

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = speedylat

       if(timestep.eq.1) then
           call nc_check(nf90_create(path=filename,cmode=ior(nf90_clobber,nf90_64bit_offset),ncid=file_id))

             ! define the dimensions
           call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
           call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))
           call nc_check(nf90_def_dim(file_id, 'Sigma_Level', dimsz, zdim_id))
           call nc_check(nf90_def_dim(file_id, 'Timestep', NF90_UNLIMITED, timedim_id))

           !Assign lat and lon ids and units
           call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
           call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

           call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
           call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
           arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /)
           arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

           do i=1, model_parameters%full_predictvars
              call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims4d,array_id))
              ! ...and assign units to them as an attribute

              call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

              call nc_check(nf90_enddef(file_id))

              !Write out the values
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d, count=varcount4d))

              call nc_check(nf90_put_var(file_id, xvar_id, lon))

              call nc_check(nf90_put_var(file_id, yvar_id, lat))

              call nc_check(nf90_redef(file_id))
           enddo
           !Lets do logp
           call nc_check(nf90_def_var(file_id,varname(5)%str,NF90_REAL,arrdims3d,array_id))

           call nc_check(nf90_put_att(file_id, array_id, "units", units(5)%str))

           call nc_check(nf90_enddef(file_id))

           !Write out the values
           call nc_check(nf90_put_var(file_id, array_id, grid3d(1,:,:),start=start3d, count=varcount3d))

           call nc_check(nf90_put_var(file_id, xvar_id, lon))

           call nc_check(nf90_put_var(file_id, yvar_id, lat))

           call nc_check(nf90_redef(file_id))

           if(model_parameters%slab_ocean_model_bool) then
             !Lets do sst
             call nc_check(nf90_def_var(file_id,varname(6)%str,NF90_REAL,arrdims3d,array_id))

             call nc_check(nf90_put_att(file_id, array_id, "units", units(6)%str))

             call nc_check(nf90_enddef(file_id))

             !Write out the values
             call nc_check(nf90_put_var(file_id, array_id, copy3d(2,:,:),start=start3d, count=varcount3d))
 
             call nc_check(nf90_put_var(file_id, xvar_id, lon))

             call nc_check(nf90_put_var(file_id, yvar_id, lat))
           endif

           ! close; done
           call nc_check(nf90_close(file_id))
        else
           call nc_check(nf90_open(filename,nf90_write,file_id))
           do i=1, model_parameters%full_predictvars
              call nc_check(nf90_inq_varid(file_id,varname(i)%str,array_id))
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d,count=varcount4d))
           enddo

           call nc_check(nf90_inq_varid(file_id,varname(5)%str,array_id))
           call nc_check(nf90_put_var(file_id, array_id, grid3d(1,:,:),start=start3d,count=varcount3d))

           if(model_parameters%slab_ocean_model_bool) then
             call nc_check(nf90_inq_varid(file_id,varname(6)%str,array_id))
             call nc_check(nf90_put_var(file_id, array_id, copy3d(2,:,:),start=start3d, count=varcount3d))
           endif

           call nc_check(nf90_close(file_id))
        endif
        return
     end subroutine

     subroutine write_netcdf_speedy_full_mpi(timestep,model_parameters,filename,mpi_res,grid4d,grid3d)
       use mpi
   
       use mod_utilities, only : mpi_type, model_parameters_type
       use stringtype, only : string

       real(kind=dp), intent(in), optional :: grid4d(:,:,:,:)
       real(kind=dp), intent(in), optional :: grid3d(:,:) !Yes I know bad name sorry

       integer, intent(in)          :: timestep

       character(len=*), intent(in) :: filename

       type(model_parameters_type), intent(in) :: model_parameters

       type(mpi_type), intent(in)   :: mpi_res

       !local netcdf stuff
       integer, parameter :: numdims4d=4, numdims3d=3, numspeedyvars=5
       integer :: dimsx,dimsy,dimsz

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d),start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims4d(numdims4d),arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(xgrid), lat(ygrid)

       real(kind=dp), allocatable :: copy(:,:,:,:)

       type(string) :: units(numspeedyvars)
       type(string) :: varname(numspeedyvars)

       logical :: makefile

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'

       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'

       units(4)%str = 'g/kg'
       varname(4)%str = 'Specific_Humidity'

       units(5)%str = 'log(surfacepressure)'
       varname(5)%str = 'logp'

       dimsx = xgrid!size(grid4d,2)
       dimsy = ygrid!size(grid4d,3)
       dimsz = zgrid!size(grid4d,4)

       varcount4d = [integer:: dimsx, dimsy, dimsz, 1]
       start4d =  [integer:: 1, 1, 1, timestep]

       varcount3d =  [integer:: dimsx, dimsy, 1 ]
       start3d =  [integer:: 1, 1, timestep ]

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = (/-87.159, -83.479, -79.777,  -76.070, -72.362, -68.652, &
               -64.942, -61.232, -57.521, -53.810, -50.099, -46.389, -42.678,-38.967, &
               -35.256, -31.545 , -27.833, -24.122, -20.411, -16.700, -12.989, -9.278, &
               -5.567, -1.856, 1.856, 5.567, 9.278, 12.989, 16.700, 20.411, &
               24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389, 50.099, &
               53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, &
               83.479, 87.159/)

       if(timestep == 1) then 
         makefile = .True.
       else
         makefile = .False.
       endif 

       if(makefile) then
           !call nc_check(nf90_create(filename,IOR(NF90_NETCDF4,NF90_MPIIO),file_id,comm=mpi_res%mpi_world%MPI_VAL, info=MPI_INFO_NULL))
           call nc_check(nf90_create(filename,IOR(NF90_NETCDF4,NF90_MPIIO),file_id,comm=mpi_res%mpi_world, info=MPI_INFO_NULL))

             ! define the dimensions
           call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
           call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))
           call nc_check(nf90_def_dim(file_id, 'Sigma_Level', dimsz, zdim_id))
           call nc_check(nf90_def_dim(file_id, 'Timestep', int(model_parameters%predictionlength/model_parameters%timestep+1,kind=int32), timedim_id))

           !Assign lat and lon ids and units
           call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
           call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

           call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
           call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
           arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /)
           arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

           do i=1, numspeedyvars-1
              call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims4d,array_id))
              ! ...and assign units to them as an attribute

              call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

              call nc_check(nf90_enddef(file_id))
 
              call nc_check(nf90_var_par_access(file_id, array_id, nf90_collective))

              !Write out the values

              call nc_check(nf90_put_var(file_id, xvar_id, lon))

              call nc_check(nf90_put_var(file_id, yvar_id, lat))

              call nc_check(nf90_redef(file_id))

           enddo
           !Lets do logp
           call nc_check(nf90_def_var(file_id,varname(5)%str,NF90_REAL,arrdims3d,array_id))

           call nc_check(nf90_put_att(file_id, array_id, "units", units(5)%str))

           call nc_check(nf90_enddef(file_id))

           call nc_check(nf90_var_par_access(file_id, array_id, nf90_collective))

           !Write out the values

           call nc_check(nf90_put_var(file_id, xvar_id, lon))

           call nc_check(nf90_put_var(file_id, yvar_id, lat)) 
           call nc_check(nf90_close(file_id))
        else
           !call nc_check(nf90_open(filename,IOR(NF90_WRITE,NF90_MPIIO),file_id,comm=mpi_res%mpi_world%MPI_VAL,info=MPI_INFO_NULL))
           call nc_check(nf90_open(filename,IOR(NF90_WRITE,NF90_MPIIO),file_id,comm=mpi_res%mpi_world,info=MPI_INFO_NULL))

           if(mpi_res%is_root) then 
             !copy data
             allocate(copy,source=grid4d)
             do i=1, numspeedyvars-1
                call nc_check(nf90_inq_varid(file_id,varname(i)%str,array_id))
                call nc_check(nf90_var_par_access(file_id,array_id,nf90_independent))
                call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d,count=varcount4d))
             enddo

             call nc_check(nf90_inq_varid(file_id,varname(5)%str,array_id))

             call nc_check(nf90_var_par_access(file_id,array_id,nf90_independent))
           
             call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d,count=varcount3d))
             deallocate(copy)
           endif 

           call nc_check(nf90_close(file_id))
        endif

        return
     end subroutine

     subroutine write_netcdf_speedy_full(grid4d,grid3d,timestep,filename)
       use stringtype, only : string

       real(kind=dp), intent(in)    :: grid4d(:,:,:,:)
       real(kind=dp), intent(in)    :: grid3d(:,:) !Yes I know bad name sorry
       integer, intent(in)          :: timestep
       character(len=*), intent(in) :: filename

       integer, parameter :: numdims4d=4, numdims3d=3, numspeedyvars=5
       integer :: dimsx,dimsy,dimsz

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d),start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims4d(numdims4d),arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(xgrid), lat(ygrid)

       real(kind=dp), allocatable :: copy(:,:,:,:)

       type(string) :: units(numspeedyvars)
       type(string) :: varname(numspeedyvars)

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'

       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'
       
       units(4)%str = 'g/kg'
       varname(4)%str = 'Specific_Humidity'

       units(5)%str = 'log(surfacepressure)'
       varname(5)%str = 'logp'

       !copy data
       allocate(copy,source=grid4d)

       dimsx = size(grid4d,2)
       dimsy = size(grid4d,3)
       dimsz = size(grid4d,4)

       varcount4d =  [integer:: dimsx, dimsy, dimsz, 1]
       start4d =  [integer:: 1, 1, 1, timestep ]

       varcount3d =  [integer:: dimsx, dimsy, 1]
       start3d =  [integer:: 1, 1, timestep ]

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = (/-87.159, -83.479, -79.777,  -76.070, -72.362, -68.652, &
               -64.942, -61.232, -57.521, -53.810, -50.099, -46.389, -42.678, -38.967, &
               -35.256, -31.545 , -27.833, -24.122, -20.411, -16.700, -12.989, -9.278, &
               -5.567, -1.856, 1.856, 5.567, 9.278, 12.989, 16.700, 20.411, &
               24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389, 50.099, &
               53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, &
               83.479, 87.159/)

       if(timestep.eq.1) then
           call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER, ncid=file_id))

             ! define the dimensions
           call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
           call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))
           call nc_check(nf90_def_dim(file_id, 'Sigma_Level', dimsz, zdim_id))
           call nc_check(nf90_def_dim(file_id, 'Timestep', NF90_UNLIMITED, timedim_id))

           !Assign lat and lon ids and units
           call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
           call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

           call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
           call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
           arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /)
           arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

           do i=1, numspeedyvars-1
              call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims4d,array_id))
              ! ...and assign units to them as an attribute

              call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

              call nc_check(nf90_enddef(file_id))
  
              !Write out the values
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d, count=varcount4d))

              call nc_check(nf90_put_var(file_id, xvar_id, lon))

              call nc_check(nf90_put_var(file_id, yvar_id, lat))

              call nc_check(nf90_redef(file_id))
           enddo
           !Lets do logp
           call nc_check(nf90_def_var(file_id,varname(5)%str,NF90_REAL,arrdims3d,array_id))

           call nc_check(nf90_put_att(file_id, array_id, "units", units(5)%str))

           call nc_check(nf90_enddef(file_id))

           !Write out the values
           call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d, count=varcount3d))

           call nc_check(nf90_put_var(file_id, xvar_id, lon))

           call nc_check(nf90_put_var(file_id, yvar_id, lat))
           ! close; done
           call nc_check(nf90_close(file_id))
        else
           call nc_check(nf90_open(filename,nf90_write,file_id))
           do i=1, numspeedyvars-1
              call nc_check(nf90_inq_varid(file_id,varname(i)%str,array_id))
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d,count=varcount4d))
           enddo

           call nc_check(nf90_inq_varid(file_id,varname(5)%str,array_id))
           call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d,count=varcount3d))

           call nc_check(nf90_close(file_id))
        endif
        
        deallocate(copy)
        return
     end subroutine 

     subroutine read_netcdf_4d(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename
      
        real(kind=dp), allocatable, intent(out) :: var(:,:,:,:) 
        
        !Parmeter
        integer, parameter :: numofdims = 4 !We can assume its a 4d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid
 
        integer :: dimids(numofdims), dim_length(numofdims)
    
        integer :: i
 
        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))
        
        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo 
      
        allocate(var(dim_length(1),dim_length(2),dim_length(3),dim_length(4)))

        call nc_check(nf90_get_var(ncid,varid,var))

        return 
     end subroutine 

     subroutine read_netcdf_3d(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        real(kind=dp), allocatable, intent(out) :: var(:,:,:)

        !Parmeter
        integer, parameter :: numofdims = 3 !We can assume its a 3d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1),dim_length(2),dim_length(3)))

        call nc_check(nf90_get_var(ncid,varid,var))

        return
     end subroutine

     subroutine write_netcdf_2d(var,varname,filename,units)
       use netcdf

       real(kind=dp), intent(in)    :: var(:,:)
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units
       

       integer, parameter :: numdims=2
       integer :: dimsx,dimsy

       integer :: file_id, xdim_id, ydim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       real(kind=dp) :: lon(96), lat(48)

       dimsx = size(var,1)
       dimsy = size(var,2)

       varcount =  [integer:: dimsx, dimsy ]
       start =  [integer:: 1, 1]

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = (/-87.159, -83.479, -79.777,  -76.070, -72.362, -68.652, &
               -64.942, -61.232, -57.521, -53.810, -50.099, -46.389, -42.678, -38.967, &
               -35.256, -31.545 , -27.833, -24.122, -20.411, -16.700, -12.989,-9.278, &
               -5.567, -1.856, 1.856, 5.567, 9.278, 12.989, 16.700, 20.411, &
               24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389, 50.099, &
               53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, &
               83.479, 87.159/)

      call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id))

             ! define the dimensions
      call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
      call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))

           !Assign lat and lon ids and units
      call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
      call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

      call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
      call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
      arrdims = (/ xdim_id, ydim_id/)

      call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims,array_id))

      call nc_check(nf90_put_att(file_id, array_id, "units", units))

      call nc_check(nf90_enddef(file_id))

      !Write out the values
      call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount))

      call nc_check(nf90_put_var(file_id, xvar_id, lon))

      call nc_check(nf90_put_var(file_id, yvar_id, lat))
      ! close; done
      call nc_check(nf90_close(file_id))
     end subroutine 

     subroutine write_netcdf_2d_non_met_data(var,varname,filename,units)
       real(kind=dp), intent(in)    :: var(:,:)
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units

       integer, parameter :: numdims=2
       integer :: dimsx,dimsy

       integer :: file_id, xdim_id, ydim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       dimsx = size(var,1)
       dimsy = size(var,2)

       varcount =  [integer:: dimsx, dimsy ]
       start =  [integer:: 1, 1]

       call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id))

       call nc_check(nf90_def_dim(file_id, 'X', dimsx, xdim_id))
       call nc_check(nf90_def_dim(file_id, 'Y', dimsy, ydim_id))

       arrdims = (/ xdim_id, ydim_id/)

       call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims,array_id))

       call nc_check(nf90_put_att(file_id, array_id, "units", units))

       call nc_check(nf90_enddef(file_id))

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount))

       call nc_check(nf90_close(file_id))
     end subroutine

     subroutine write_netcdf_2d_non_met_data_new(var,varname,filename,units,x_dim,y_dim)
       real(kind=dp), intent(in)    :: var(:,:)
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units
       character(len=*), intent(in) :: x_dim, y_dim

       integer, parameter :: numdims=2
       integer :: dimsx,dimsy

       integer :: file_id, xdim_id, ydim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       dimsx = size(var,1)
       dimsy = size(var,2)

       varcount =  [integer:: dimsx, dimsy ]
       start =  [integer:: 1, 1]

       if(.not.file_exists(filename)) then
          call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id),message='write_netcdf_2d_non_met_data create')
       else
          call nc_check(nf90_open(filename,nf90_write,file_id),message='write_netcdf_2d_non_met_data nf90_open')
          call nc_check(nf90_redef(file_id))
       endif

       call nc_check(nf90_def_dim(file_id, x_dim, dimsx, xdim_id),message='write_netcdf_2d_non_met_data nf90_def_dimx'//filename//'_'//x_dim)
       call nc_check(nf90_def_dim(file_id, y_dim, dimsy, ydim_id),message='write_netcdf_2d_non_met_data nf90_def_dimy'//filename)

       arrdims = (/ xdim_id, ydim_id/)

       call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims,array_id),message='write_netcdf_2d_non_met_data nf90_def_var'//filename)

       call nc_check(nf90_put_att(file_id, array_id, "units", units),message='write_netcdf_2d_non_met_data nf90_put_att'//filename)

       call nc_check(nf90_enddef(file_id),message='write_netcdf_2d_non_met_data nf90_enddef'//filename)

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount),message='write_netcdf_2d_non_met_data nf90_put_var'//filename)

       call nc_check(nf90_close(file_id),message='write_netcdf_2d_non_met_data nf90_close')
     end subroutine

     subroutine write_netcdf_2d_non_met_data_new_double(var,varname,filename,units,x_dim,y_dim)
       real(kind=dp), intent(in)    :: var(:,:)
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units
       character(len=*), intent(in) :: x_dim, y_dim

       integer, parameter :: numdims=2
       integer :: dimsx,dimsy

       integer :: file_id, xdim_id, ydim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       dimsx = size(var,1)
       dimsy = size(var,2)

       varcount =  [integer:: dimsx, dimsy ]
       start =  [integer:: 1, 1]

       if(.not.file_exists(filename)) then
          call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id),message='write_netcdf_2d_non_met_data create')
       else
          call nc_check(nf90_open(filename,nf90_write,file_id),message='write_netcdf_2d_non_met_data nf90_open')
          call nc_check(nf90_redef(file_id))
       endif

       call nc_check(nf90_def_dim(file_id, x_dim, dimsx, xdim_id),message='write_netcdf_2d_non_met_data nf90_def_dimx'//filename//'_'//x_dim)
       call nc_check(nf90_def_dim(file_id, y_dim, dimsy, ydim_id),message='write_netcdf_2d_non_met_data nf90_def_dimy'//filename)

       arrdims = (/ xdim_id, ydim_id/)

       call nc_check(nf90_def_var(file_id,varname,NF90_DOUBLE,arrdims,array_id),message='write_netcdf_2d_non_met_data nf90_def_var'//filename)

       call nc_check(nf90_put_att(file_id, array_id, "units", units),message='write_netcdf_2d_non_met_data nf90_put_att'//filename)

       call nc_check(nf90_enddef(file_id),message='write_netcdf_2d_non_met_data nf90_enddef'//filename)

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount),message='write_netcdf_2d_non_met_data nf90_put_var'//filename)

       call nc_check(nf90_close(file_id),message='write_netcdf_2d_non_met_data nf90_close')
     end subroutine

     subroutine read_netcdf_2d_dp(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        real(kind=dp), allocatable, intent(out) :: var(:,:)

        !Parmeter
        integer, parameter :: numofdims = 2 !We can assume its a 2d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1),dim_length(2)))

        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))
        return
     end subroutine

     subroutine write_netcdf_1d_non_met_data_int(var,varname,filename,units)
       integer,       intent(in)    :: var(:)
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units

       integer, parameter :: numdims=1
       integer :: dimsx,dimsy

       integer :: file_id, xdim_id, ydim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       dimsx = size(var,1)

       varcount = [dimsx]
       start =  [integer:: 1]

       call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id))

       call nc_check(nf90_def_dim(file_id, 'X', dimsx, xdim_id))

       arrdims = (/ xdim_id/)

       call nc_check(nf90_def_var(file_id,varname,NF90_INT,arrdims,array_id))

       call nc_check(nf90_put_att(file_id, array_id, "units", units))

       call nc_check(nf90_enddef(file_id))

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount))

       call nc_check(nf90_close(file_id))
     end subroutine

     subroutine write_netcdf_1d_non_met_data_int_new(var,varname,filename,units,x_dim)
       integer,       intent(in)    :: var(:)
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units
       character(len=*), intent(in) :: x_dim

       integer, parameter :: numdims=1
       integer :: dimsx,dimsy

       integer :: file_id, xdim_id, ydim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       dimsx = size(var,1)

       varcount = [dimsx]
       start =  [integer:: 1]


       if(.not.file_exists(filename)) then
          call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id),message='write_netcdf_1d_non_met_data_int nf90_create')
       else
          call nc_check(nf90_open(filename,nf90_write,file_id),message='write_netcdf_1d_non_met_data_int nf90_open')
          call nc_check(nf90_redef(file_id),message='write_netcdf_1d_non_met_data_int nf_redef')
       endif

       call nc_check(nf90_def_dim(file_id, x_dim, dimsx, xdim_id),message='write_netcdf_1d_non_met_data_int nf90_def_dim')

       arrdims = (/ xdim_id/)

       call nc_check(nf90_def_var(file_id,varname,NF90_INT,arrdims,array_id),message='write_netcdf_1d_non_met_data_int nf90_def_var')

       call nc_check(nf90_put_att(file_id, array_id, "units", units),message='write_netcdf_1d_non_met_data_int nf90_put_att')

       call nc_check(nf90_enddef(file_id),message='write_netcdf_1d_non_met_data_int nf90_enddef')

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount),message='write_netcdf_1d_non_met_data_int nf90_put_var')

       call nc_check(nf90_close(file_id))
     end subroutine

     subroutine write_netcdf_1d_non_met_data_real(var,varname,filename,units)
       real(kind=dp), intent(in)    :: var(:)
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units

       integer, parameter :: numdims=1
       integer :: dimsx,dimsy

       integer :: file_id, xdim_id, ydim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       dimsx = size(var,1)

       varcount = [dimsx]
       start =  [integer:: 1]

       call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id))

       call nc_check(nf90_def_dim(file_id, 'X', dimsx, xdim_id))

       arrdims = (/ xdim_id/)

       call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims,array_id))

       call nc_check(nf90_put_att(file_id, array_id, "units", units))

       call nc_check(nf90_enddef(file_id))

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount))

       call nc_check(nf90_close(file_id))
     end subroutine

     subroutine write_netcdf_1d_non_met_data_real_new(var,varname,filename,units,x_dim)
       real(kind=dp), intent(in)    :: var(:)
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units
       character(len=*), intent(in) :: x_dim

       integer, parameter :: numdims=1
       integer :: dimsx,dimsy

       integer :: file_id, xdim_id, ydim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       dimsx = size(var,1)

       varcount = [dimsx]
       start =  [integer:: 1]

       if(.not.file_exists(filename)) then
          call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id),message='write_netcdf_1d_non_met_data_real nf90_create')
       else
          call nc_check(nf90_open(filename,nf90_write,file_id),message='write_netcdf_1d_non_met_data_real nf90_open')
          call nc_check(nf90_redef(file_id),message='write_netcdf_1d_non_met_data_real nf_redef')
       endif

       call nc_check(nf90_def_dim(file_id, x_dim, dimsx, xdim_id),message='write_netcdf_1d_non_met_data_real nf90_def_dim')

       arrdims = (/ xdim_id/)

       call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims,array_id),message='write_netcdf_1d_non_met_data_real nf90_def_var')

       call nc_check(nf90_put_att(file_id, array_id, "units", units),message='write_netcdf_1d_non_met_data_real nf90_put_att')

       call nc_check(nf90_enddef(file_id),message='write_netcdf_1d_non_met_data_real nf90_enddef')

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount),message='write_netcdf_1d_non_met_data_real nf90_put_var')

       call nc_check(nf90_close(file_id),message='write_netcdf_1d_non_met_data_real nf90_close')
     end subroutine
 
     subroutine write_netcdf_1d_non_met_data_real_new_double(var,varname,filename,units,x_dim)
       real(kind=dp), intent(in)    :: var(:)
       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units
       character(len=*), intent(in) :: x_dim

       integer, parameter :: numdims=1
       integer :: dimsx,dimsy

       integer :: file_id, xdim_id, ydim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       dimsx = size(var,1)

       varcount = [dimsx]
       start =  [integer:: 1]

       if(.not.file_exists(filename)) then
          call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id),message='write_netcdf_1d_non_met_data_real nf90_create')
       else
          call nc_check(nf90_open(filename,nf90_write,file_id),message='write_netcdf_1d_non_met_data_real nf90_open')
          call nc_check(nf90_redef(file_id),message='write_netcdf_1d_non_met_data_real nf_redef')
       endif

       call nc_check(nf90_def_dim(file_id, x_dim, dimsx, xdim_id),message='write_netcdf_1d_non_met_data_real nf90_def_dim')

       arrdims = (/ xdim_id/)

       call nc_check(nf90_def_var(file_id,varname,NF90_DOUBLE,arrdims,array_id),message='write_netcdf_1d_non_met_data_real nf90_def_var')

       call nc_check(nf90_put_att(file_id, array_id, "units", units),message='write_netcdf_1d_non_met_data_real nf90_put_att')

       call nc_check(nf90_enddef(file_id),message='write_netcdf_1d_non_met_data_real nf90_enddef')

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount),message='write_netcdf_1d_non_met_data_real nf90_put_var')

       call nc_check(nf90_close(file_id),message='write_netcdf_1d_non_met_data_real nf90_close')
     end subroutine

     subroutine write_netcdf_2d_reservoir_matrices_mpi(mpi_res,var,varname,filename,units)
       use mpi
       use netcdf

       use mod_utilities, only : mpi_type, reservoir_type

       type(mpi_type), intent(in)   :: mpi_res

       real(kind=dp), intent(in)    :: var(:,:)

       character(len=*), intent(in) :: varname
       character(len=*), intent(in) :: filename
       character(len=*), intent(in) :: units

       integer, parameter :: numdims=3
       integer :: dimsx,dimsy,dimsworker

       integer :: file_id, xdim_id, ydim_id, workerdim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start(numdims),varcount(numdims)
       integer, dimension(numdims) :: arrdims

       integer :: i, counter

       dimsx = size(var,1)
       dimsy = size(var,2)
       dimsworker = mpi_res%numprocs

       varcount = [integer:: 1, dimsx, dimsy ]
       start = [integer:: mpi_res%proc_num, 1, 1]

       if(.not.file_exists(filename)) then
          !call nc_check(nf90_create(path=filename,cmode=IOR(NF90_NETCDF4, NF90_MPIIO),ncid=file_id,comm=mpi_res%mpi_world%MPI_VAL, info=MPI_INFO_NULL))
          call nc_check(nf90_create(path=filename,cmode=IOR(NF90_NETCDF4, NF90_MPIIO),ncid=file_id,comm=mpi_res%mpi_world, info=MPI_INFO_NULL))

          call nc_check(nf90_def_dim(file_id, 'worker', dimsworker, workerdim_id))
          call nc_check(nf90_def_dim(file_id, 'X', dimsx, xdim_id))
          call nc_check(nf90_def_dim(file_id, 'Y', dimsy, ydim_id))

          arrdims = (/ workerdim_id,xdim_id, ydim_id/)

          call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims,array_id))

          call nc_check(nf90_put_att(file_id, array_id, "units", units))

          call nc_check(nf90_enddef(file_id))
       endif
      
       
       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var,start=start, count=varcount))

       call nc_check(nf90_close(file_id))
     end subroutine

     subroutine read_netcdf_1d_dp(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        real(kind=dp), allocatable, intent(out) :: var(:)

        !Parmeter
        integer, parameter :: numofdims = 1!We can assume its a 1d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1)))

        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))
        return
     end subroutine

     subroutine read_netcdf_1d_int(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        integer, allocatable, intent(out) :: var(:)

        !Parmeter
        integer, parameter :: numofdims = 1!We can assume its a 1d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1)))

        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))
        return
     end subroutine

     subroutine nc_check(status,worker,message)
         use netcdf
         
         integer, intent (in) :: status
         integer, intent(in), optional :: worker
         character(len=*), intent(in), optional :: message

         if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            if(present(worker)) print *, 'processor',worker
            if(present(message)) print *, message
            stop "Stopped"
         end if
     end subroutine nc_check 

    !-----------Parallel netcdf Routines----------!

    subroutine read_era_data_parallel(filename,model_parameters,mpi_res,grid,era_data,start_time_arg,stride_arg)
      use mpi
      use netcdf

      use mod_utilities, only : era_data_type, mpi_type, grid_type, model_parameters_type
      use stringtype, only    : string

      type(mpi_type), intent(in)              :: mpi_res
      type(grid_type), intent(in)             :: grid
      type(model_parameters_type), intent(in) :: model_parameters
       
      character(len=*), intent(in)            :: filename

      type(era_data_type), intent(inout)        :: era_data

      integer, intent(in), optional           :: start_time_arg, stride_arg

      !Local netcdf stuff
      integer :: start4d(4), start3d(3)
      integer :: count4d(4), count3d(3)
      integer :: stride4d(4), stride3d(3)
      integer :: ncid, varid
      integer :: dimids(4), dim_length(4)
      integer :: i,j
      integer :: start_time, stride

      integer, parameter :: maxnum_of_dims=4
      integer, parameter :: num_of_era_vars=5

      type(string) :: era_var_names(num_of_era_vars)

      if(present(start_time_arg)) then
        start_time = start_time_arg
      else 
        start_time = 1
      endif 

      if(present(stride_arg)) then
        stride = stride_arg
      else
        stride = 1
      endif

      era_var_names(1)%str = 'Temperature'
      era_var_names(2)%str = 'U-wind'
      era_var_names(3)%str = 'V-wind'
      era_var_names(4)%str = 'Specific_Humidity'
      era_var_names(5)%str = 'logp'

      !call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world%MPI_VAL, info=MPI_INFO_NULL)) 
      call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world, info=MPI_INFO_NULL))
    
      do i=1, num_of_era_vars-1
         call nc_check(nf90_inq_varid(ncid,era_var_names(i)%str,varid))
    
         !call nc_check(nf90_var_par_access(ncid, varid, nf90_collective))
 
         call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids)) 

    
         if(i == 1) then 
           do j=1,maxnum_of_dims
              call nc_check(nf90_inquire_dimension(ncid,dimids(j),len=dim_length(j)))
           enddo
           !lets allocate grid4
           !allocate(era_data%eravariables(num_of_era_vars-1,grid%inputxchunk,grid%inputychunk,grid%inputzchunk,dim_length(4)/stride))
           allocate(era_data%eravariables(num_of_era_vars-1,grid%inputxchunk,grid%inputychunk,zgrid,dim_length(4)/stride))
         endif 

         !check if its a periodicboundary or not 
         if(grid%periodicboundary) then 
           !Do stuff
           !start4d =  [integer:: grid%input_xstart,grid%input_ystart,grid%input_zstart,start_time]
           !count4d =  [integer:: xgrid-grid%input_xstart+1,grid%inputychunk,grid%inputzchunk,dim_length(4)/stride]
           start4d =  [integer:: grid%input_xstart,grid%input_ystart,1,start_time] 
           count4d =  [integer:: xgrid-grid%input_xstart+1,grid%inputychunk,zgrid,dim_length(4)/stride]
           stride4d = [integer:: 1,1,1,stride]

           print *, 'start4d_1',start4d,mpi_res%proc_num
           print *, 'count4d_1',count4d,mpi_res%proc_num
       
           call nc_check(nf90_get_var(ncid,varid,era_data%eravariables(i,1:(xgrid-(grid%input_xstart-1)),:,:,:),start=start4d,count=count4d))!,stride=stride4d),mpi_res%proc_num,'first')

           !start4d = [integer:: 1,grid%input_ystart,grid%input_zstart,start_time]
           !count4d = [integer:: grid%input_xend,grid%inputychunk,grid%inputzchunk,dim_length(4)/stride]
           start4d = [integer:: 1,grid%input_ystart,1,start_time]
           count4d = [integer:: grid%input_xend,grid%inputychunk,zgrid,dim_length(4)/stride]

           print *, 'start4d_2',start4d,mpi_res%proc_num
           print *, 'count4d_2',count4d,mpi_res%proc_num
           call nc_check(nf90_get_var(ncid,varid,era_data%eravariables(i,(xgrid-(grid%input_xstart-1))+1:grid%inputxchunk,:,:,:),start=start4d,count=count4d))!,stride=stride4d),mpi_res%proc_num,'second')    
         else
           !start4d = [integer:: grid%input_xstart,grid%input_ystart,grid%input_zstart,start_time]
           !count4d = [integer:: grid%inputxchunk,grid%inputychunk,grid%inputzchunk,dim_length(4)/stride] 
           start4d = [integer:: grid%input_xstart,grid%input_ystart,1,start_time]
           count4d = [integer:: grid%inputxchunk,grid%inputychunk,zgrid,dim_length(4)/stride]
           stride4d = [integer:: 1,1,1,stride]
 
           call nc_check(nf90_get_var(ncid,varid,era_data%eravariables(i,:,:,:,:),start=start4d,count=count4d))!,stride=stride4d))      
        endif 
      enddo 
     
      !Time to do logp
      call nc_check(nf90_inq_varid(ncid,era_var_names(num_of_era_vars)%str,varid)) 
     
      !call nc_check(nf90_var_par_access(ncid, varid, nf90_collective))
 
      allocate(era_data%era_logp(grid%inputxchunk,grid%inputychunk,dim_length(4)/stride))
 
      if(grid%periodicboundary) then
 
        start3d = [integer:: grid%input_xstart,grid%input_ystart,start_time]
        count3d = [integer:: xgrid-grid%input_xstart+1,grid%inputychunk,dim_length(4)/stride]
        stride3d = [integer:: 1,1,stride]

        call nc_check(nf90_get_var(ncid,varid,era_data%era_logp(1:(xgrid-(grid%input_xstart-1)),:,:),start=start3d,count=count3d))!,stride=stride3d))

        start3d = [integer:: 1,grid%input_ystart,start_time]
        count3d = [integer:: grid%input_xend,grid%inputychunk,dim_length(4)/stride]


        call nc_check(nf90_get_var(ncid,varid,era_data%era_logp((xgrid-(grid%input_xstart-1))+1:grid%inputxchunk,:,:),start=start3d,count=count3d))!,stride=stride3d))

      else
        start3d = [integer:: grid%input_xstart,grid%input_ystart,start_time]
        count3d = [integer:: grid%inputxchunk,grid%inputychunk,dim_length(4)/stride] 
        stride3d = [integer:: 1,1,stride]

        call nc_check(nf90_get_var(ncid,varid,era_data%era_logp,start=start3d,count=count3d))!,stride=stride3d))
      endif 

      call nc_check(nf90_close(ncid))
    end subroutine  
   
    subroutine read_era_data_parallel_old(filename,mpi_res,grid,era_data)
      use mpi
      use netcdf

      use mod_utilities, only : era_data_type, mpi_type, grid_type, reservoir_type
      use stringtype, only    : string

      type(mpi_type), intent(in)          :: mpi_res
      type(grid_type), intent(in)         :: grid

      character(len=*), intent(in)        :: filename

      type(era_data_type), intent(out) :: era_data

      !Local netcdf stuff
      integer :: start4d(4), start3d(3)
      integer :: count4d(4), count3d(3)
      integer :: ncid, varid
      integer :: dimids(4), dim_length(4)
      integer :: i,j

      integer, parameter :: maxnum_of_dims=4
      integer, parameter :: num_of_era_vars=5

      type(string) :: era_var_names(num_of_era_vars)

      era_var_names(1)%str = 'Temperature'
      era_var_names(2)%str = 'U-wind'
      era_var_names(3)%str = 'V-wind'
      era_var_names(4)%str = 'Specific_Humidity'
      era_var_names(5)%str = 'logp'

      call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO),ncid,comm=mpi_res%mpi_world, info=MPI_INFO_NULL))

      do i=1, num_of_era_vars-1
         call nc_check(nf90_inq_varid(ncid,era_var_names(i)%str,varid))

         call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

         do j=1,maxnum_of_dims
            call nc_check(nf90_inquire_dimension(ncid,dimids(j),len=dim_length(j)))
         enddo 
   
         if(i == 1) then
           !lets allocate era data
           allocate(era_data%eravariables(num_of_era_vars-1,grid%inputxchunk,grid%inputychunk,dim_length(3),dim_length(4)))
         endif

         !check if its a periodicboundary or not
         if(grid%periodicboundary) then
           !Do stuff
           start4d = [grid%input_xstart,grid%input_ystart,1,1]
           count4d = [xgrid-grid%input_xstart+1,grid%inputychunk,dim_length(3),dim_length(4)]

           call nc_check(nf90_get_var(ncid,varid,era_data%eravariables(i,1:(xgrid-(grid%input_xstart-1)),:,:,:),start=start4d,count=count4d),mpi_res%proc_num,'first')

           start4d = [1,grid%input_ystart,1,1]
           count4d = [grid%input_xend,grid%inputychunk,dim_length(3),dim_length(4)]

           call nc_check(nf90_get_var(ncid,varid,era_data%eravariables(i,(xgrid-(grid%input_xstart-1))+1:grid%inputxchunk,:,:,:),start=start4d,count=count4d),mpi_res%proc_num,'second')
         else
           start4d = [grid%input_xstart,grid%input_ystart,1,1]
           count4d = [grid%inputxchunk,grid%inputychunk,dim_length(3),dim_length(4)]

           call nc_check(nf90_get_var(ncid,varid,era_data%eravariables(i,:,:,:,:),start=start4d,count=count4d))
         endif
      enddo

      !Time to do logp
      call nc_check(nf90_inq_varid(ncid,era_var_names(num_of_era_vars)%str,varid))

      allocate(era_data%era_logp(grid%inputxchunk,grid%inputychunk,dim_length(4)))

      if(grid%periodicboundary) then

        start3d = [grid%input_xstart,grid%input_ystart,1]
        count3d = [xgrid-grid%input_xstart+1,grid%inputychunk,dim_length(4)]

        call nc_check(nf90_get_var(ncid,varid,era_data%era_logp(1:(xgrid-(grid%input_xstart-1)),:,:),start=start3d,count=count3d))

        start3d = [1,grid%input_ystart,1]
        count3d = [grid%input_xend,grid%inputychunk,dim_length(4)]

        call nc_check(nf90_get_var(ncid,varid,era_data%era_logp((xgrid-(grid%input_xstart-1))+1:grid%inputxchunk,:,:),start=start3d,count=count3d))
       endif

       call nc_check(nf90_close(ncid))
    end subroutine

    subroutine read_speedy_data_parallel_old(filename,mpi_res,grid,speedy_data)
      use mpi
      use netcdf

      use mod_utilities, only : speedy_data_type, mpi_type, grid_type
      use stringtype, only    : string

      type(mpi_type), intent(in)          :: mpi_res
      type(grid_type), intent(in)         :: grid

      character(len=*), intent(in)        :: filename

      type(speedy_data_type), intent(out) :: speedy_data

      !Local netcdf stuff
      integer :: start4d(4), start3d(3)
      integer :: count4d(4), count3d(3)
      integer :: ncid, varid
      integer :: dimids(4), dim_length(4)
      integer :: i,j

      integer, parameter :: maxnum_of_dims=4
      integer, parameter :: num_of_speedy_vars=5

      type(string) :: speedy_var_names(num_of_speedy_vars)

      speedy_var_names(1)%str = 'Temperature'
      speedy_var_names(2)%str = 'U-wind'
      speedy_var_names(3)%str = 'V-wind'
      speedy_var_names(4)%str = 'Specific_Humidity'
      speedy_var_names(5)%str = 'logp'

      call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO),ncid,comm=mpi_res%mpi_world, info=MPI_INFO_NULL))

      do i=1, num_of_speedy_vars-1
         call nc_check(nf90_inq_varid(ncid,speedy_var_names(i)%str,varid))

         call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

         if(i == 1) then
           do j=1,maxnum_of_dims
              call nc_check(nf90_inquire_dimension(ncid,dimids(j),len=dim_length(j)))
           enddo
           
           !lets allocate speedy data
           allocate(speedy_data%speedyvariables(num_of_speedy_vars-1,grid%resxchunk,grid%resychunk,dim_length(3),dim_length(4)))
         endif

         start4d = [grid%res_xstart,grid%res_ystart,1,1]
         count4d = [grid%resxchunk,grid%resychunk,dim_length(3),dim_length(4)]

         call nc_check(nf90_get_var(ncid,varid,speedy_data%speedyvariables(i,:,:,:,:),start=start4d,count=count4d))
      enddo

      !Time to do logp
      call nc_check(nf90_inq_varid(ncid,speedy_var_names(num_of_speedy_vars)%str,varid))

      allocate(speedy_data%speedy_logp(grid%resxchunk,grid%resychunk,dim_length(4)))

      start3d = [grid%res_xstart,grid%res_ystart,1]
      count3d = [grid%resxchunk,grid%resychunk,dim_length(4)]

      call nc_check(nf90_get_var(ncid,varid,speedy_data%speedy_logp,start=start3d,count=count3d))

      call nc_check(nf90_close(ncid))
    end subroutine

    subroutine read_speedy_data_parallel(filename,mpi_res,grid,speedy_data,start_time_arg,stride_arg)
      use mpi
      use netcdf

      use mod_utilities, only : speedy_data_type, mpi_type, grid_type
      use stringtype, only    : string

      type(mpi_type), intent(in)          :: mpi_res
      type(grid_type), intent(in)         :: grid
       
      character(len=*), intent(in)        :: filename

      type(speedy_data_type), intent(out) :: speedy_data
 
      integer, intent(in), optional           :: start_time_arg, stride_arg
 
      !Local netcdf stuff
      integer :: start4d(4), start3d(3)
      integer :: count4d(4), count3d(3)
      integer :: stride4d(4), stride3d(3)
      integer :: ncid, varid
      integer :: dimids(4), dim_length(4)
      integer :: i,j
      integer :: start_time, stride

      integer, parameter :: maxnum_of_dims=4
      integer, parameter :: num_of_speedy_vars=5

      type(string) :: speedy_var_names(num_of_speedy_vars)

      if(present(start_time_arg)) then
        start_time = start_time_arg
      else
        start_time = 1
      endif

      if(present(stride_arg)) then
        stride = stride_arg
      else
        stride = 1
      endif

      speedy_var_names(1)%str = 'Temperature'
      speedy_var_names(2)%str = 'U-wind'
      speedy_var_names(3)%str = 'V-wind'
      speedy_var_names(4)%str = 'Specific_Humidity'
      speedy_var_names(5)%str = 'logp'

      !call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world%MPI_VAL, info=MPI_INFO_NULL)) 
      call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world, info=MPI_INFO_NULL))
    
      do i=1, num_of_speedy_vars-1
         call nc_check(nf90_inq_varid(ncid,speedy_var_names(i)%str,varid))
    
         !call nc_check(nf90_var_par_access(ncid, varid, nf90_collective))
 
         call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids)) 

         if(i == 1) then
           do j=1,maxnum_of_dims
              call nc_check(nf90_inquire_dimension(ncid,dimids(j),len=dim_length(j)))
           enddo

           !lets allocate speedy data
           allocate(speedy_data%speedyvariables(num_of_speedy_vars-1,grid%resxchunk,grid%resychunk,zgrid,dim_length(4)/stride))
         endif 
         
         !start4d = [integer:: grid%res_xstart,grid%res_ystart,grid%res_zstart,start_time]
         !count4d = [integer:: grid%resxchunk,grid%resychunk,grid%reszchunk,dim_length(4)/stride]
         start4d = [integer:: grid%res_xstart,grid%res_ystart,1,start_time]
         count4d = [integer:: grid%resxchunk,grid%resychunk,zgrid,dim_length(4)/stride]
         stride4d = [integer:: 1,1,1,stride]
 
         call nc_check(nf90_get_var(ncid,varid,speedy_data%speedyvariables(i,:,:,:,:),start=start4d,count=count4d,stride=stride4d))      
      enddo 
     
      !Time to do logp
      call nc_check(nf90_inq_varid(ncid,speedy_var_names(num_of_speedy_vars)%str,varid)) 

      !call nc_check(nf90_var_par_access(ncid, varid, nf90_collective))
 
      allocate(speedy_data%speedy_logp(grid%resxchunk,grid%resychunk,dim_length(4)/stride))
 
      start3d = [integer:: grid%res_xstart,grid%res_ystart,start_time]
      count3d = [integer:: grid%resxchunk,grid%resychunk,dim_length(4)/stride] 
      stride3d = [integer:: 1,1,stride]

      call nc_check(nf90_get_var(ncid,varid,speedy_data%speedy_logp,start=start3d,count=count3d,stride=stride3d))

      call nc_check(nf90_close(ncid))
    end subroutine

    subroutine write_prediction_local_region_vert_level_mpi(grid,model_parameters,mpi_res,grid4d,grid3d,timestep,filename,make_file)
       use mpi

       use mod_utilities, only : mpi_type, grid_type, model_parameters_type
       use stringtype, only : string

       type(grid_type), intent(in)             :: grid
       type(model_parameters_type), intent(in) :: model_parameters 
       type(mpi_type), intent(in)              :: mpi_res
   
       real(kind=dp), intent(inout)            :: grid4d(:,:,:,:), grid3d(:,:)
       integer, intent(in)                     :: timestep
     
       character(len=*), intent(in)            :: filename 

       logical, intent(in)                     :: make_file

       integer, parameter :: numdims4d=4, numdims3d=3, numspeedyvars=5
       integer :: dimsx,dimsy,dimsz

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d),start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims4d(numdims4d),arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(xgrid), lat(ygrid)

       real(kind=dp), allocatable :: copy(:,:,:,:)

       type(string) :: units(numspeedyvars)
       type(string) :: varname(numspeedyvars)

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'

       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'

       units(4)%str = 'g/kg'
       varname(4)%str = 'Specific_Humidity'

       units(5)%str = 'log(surfacepressure)'
       varname(5)%str = 'logp'

       !copy data
       allocate(copy,source=grid4d)

       dimsx = size(grid4d,2)
       dimsy = size(grid4d,3)
       dimsz = size(grid4d,4)

       varcount4d = [integer:: dimsx, dimsy, dimsz, 1]
       start4d = [integer:: grid%res_xstart, grid%res_ystart, grid%res_zstart, timestep]

       varcount3d = [integer:: dimsx, dimsy, 1]
       start3d = [integer:: grid%res_xstart, grid%res_ystart, timestep]

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = (/-87.159, -83.479, -79.777,  -76.070, -72.362, -68.652, &
               -64.942, -61.232, -57.521, -53.810, -50.099, -46.389, -42.678,-38.967, &
               -35.256, -31.545 , -27.833, -24.122, -20.411, -16.700, -12.989, -9.278, &
               -5.567, -1.856, 1.856, 5.567, 9.278, 12.989, 16.700, 20.411, &
               24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389, 50.099, &
               53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, &
               83.479, 87.159/)

       if(make_file) then
           !call nc_check(nf90_create(filename,IOR(NF90_NETCDF4,NF90_MPIIO),file_id,comm=mpi_res%mpi_world%MPI_VAL, info=MPI_INFO_NULL))
           call nc_check(nf90_create(filename,IOR(NF90_NETCDF4,NF90_MPIIO),file_id,comm=mpi_res%mpi_world, info=MPI_INFO_NULL))

             ! define the dimensions
           call nc_check(nf90_def_dim(file_id, 'Lon', xgrid, xdim_id))
           call nc_check(nf90_def_dim(file_id, 'Lat', ygrid, ydim_id))
           call nc_check(nf90_def_dim(file_id, 'Sigma_Level', zgrid, zdim_id))
           call nc_check(nf90_def_dim(file_id, 'Timestep', int(model_parameters%predictionlength/model_parameters%timestep,kind=int32), timedim_id))

           !Assign lat and lon ids and units
           call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
           call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

           call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
           call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
           arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /)
           arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

           do i=1, numspeedyvars-1
              call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims4d,array_id))
              ! ...and assign units to them as an attribute

              call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

              call nc_check(nf90_enddef(file_id))
              
              !parallel io is more complicated
              ! Unlimited dimensions require collective writes
              call nc_check(nf90_var_par_access(file_id, array_id, nf90_collective))

              !Write out the values
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d, count=varcount4d))

              call nc_check(nf90_put_var(file_id, xvar_id, lon))

              call nc_check(nf90_put_var(file_id, yvar_id, lat))

              call nc_check(nf90_redef(file_id))
           enddo
           !Lets do logp
           call nc_check(nf90_def_var(file_id,varname(5)%str,NF90_REAL,arrdims3d,array_id))

           call nc_check(nf90_put_att(file_id, array_id, "units", units(5)%str))

           call nc_check(nf90_enddef(file_id))

           call nc_check(nf90_var_par_access(file_id,array_id,nf90_independent))

           if(grid%logp_bool) then 
             !Write out the values
             call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d, count=varcount3d))
           endif 
            
           ! close; done
           call nc_check(nf90_close(file_id))
        else
           !call nc_check(nf90_open(filename,IOR(NF90_WRITE,NF90_MPIIO),file_id,comm=mpi_res%mpi_world%MPI_VAL,info=MPI_INFO_NULL),message='nf90_open')
           call nc_check(nf90_open(filename,IOR(NF90_WRITE,NF90_MPIIO),file_id,comm=mpi_res%mpi_world,info=MPI_INFO_NULL),message='nf90_open')
 
           do i=1, numspeedyvars-1
              call nc_check(nf90_inq_varid(file_id,varname(i)%str,array_id))
              call nc_check(nf90_var_par_access(file_id,array_id,nf90_independent))
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d,count=varcount4d))
           enddo

           call nc_check(nf90_inq_varid(file_id,varname(5)%str,array_id))

           call nc_check(nf90_var_par_access(file_id,array_id,nf90_independent))
 
           if(grid%logp_bool) then
             call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d,count=varcount3d))
           endif 

           call nc_check(nf90_close(file_id))
        endif

        deallocate(copy)
        return
   end subroutine write_prediction_local_region_vert_level_mpi
   
   subroutine read_prediction_local_region_vert_level_mpi(grid,model_parameters,mpi_res,grid4d,grid2d,timestep,filename)
      use mpi
      use netcdf

      use mod_utilities, only : era_data_type, mpi_type, grid_type, model_parameters_type
      use stringtype, only    : string

      type(mpi_type), intent(in)                :: mpi_res
      type(grid_type), intent(in)               :: grid
      type(model_parameters_type), intent(in)   :: model_parameters

      integer, intent(in)                       :: timestep

      real(kind=dp), allocatable, intent(out)   :: grid4d(:,:,:,:), grid2d(:,:)

      character(len=*), intent(in)              :: filename


      !Local netcdf stuff
      integer :: start4d(4), start3d(3)
      integer :: count4d(4), count3d(3)
      integer :: ncid, varid
      integer :: dimids(4), dim_length(4)
      integer :: i,j

      integer, parameter :: maxnum_of_dims=4
      integer, parameter :: num_of_era_vars=5

      type(string) :: era_var_names(num_of_era_vars)

      era_var_names(1)%str = 'Temperature'
      era_var_names(2)%str = 'U-wind'
      era_var_names(3)%str = 'V-wind'
      era_var_names(4)%str = 'Specific_Humidity'
      era_var_names(5)%str = 'logp'

      !call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world%MPI_VAL, info=MPI_INFO_NULL))
      call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world, info=MPI_INFO_NULL))
   
      if(allocated(grid4d)) deallocate(grid4d)
      if(allocated(grid2d)) deallocate(grid2d)

      !lets allocate grid4
      if(allocated(grid4d)) print *, 'grid4d in mod io',shape(grid4d)
      allocate(grid4d(num_of_era_vars-1,grid%inputxchunk,grid%inputychunk,grid%inputzchunk))
      do i=1, num_of_era_vars-1
         call nc_check(nf90_inq_varid(ncid,era_var_names(i)%str,varid))

         call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

         do j=1,maxnum_of_dims
            call nc_check(nf90_inquire_dimension(ncid,dimids(j),len=dim_length(j)))
         enddo

         !check if its a periodicboundary or not
         if(grid%periodicboundary) then
           !Do stuff
           start4d = [integer:: grid%input_xstart,grid%input_ystart,grid%input_zstart,timestep]
           count4d = [integer:: xgrid-grid%input_xstart+1,grid%inputychunk,grid%inputzchunk,1]

           call nc_check(nf90_get_var(ncid,varid,grid4d(i,1:(xgrid-(grid%input_xstart-1)),:,:),start=start4d,count=count4d),mpi_res%proc_num,'first')

           start4d = [integer:: 1,grid%input_ystart,grid%input_zstart,timestep]
           count4d = [integer:: grid%input_xend,grid%inputychunk,grid%inputychunk,1]

           call nc_check(nf90_get_var(ncid,varid,grid4d(i,(xgrid-(grid%input_xstart-1))+1:grid%inputxchunk,:,:),start=start4d,count=count4d),mpi_res%proc_num,'second')
         else
           start4d = [integer:: grid%input_xstart,grid%input_ystart,grid%input_zstart,timestep]
           count4d = [integer:: grid%inputxchunk,grid%inputychunk,grid%inputzchunk,1]

           call nc_check(nf90_get_var(ncid,varid,grid4d(i,:,:,:),start=start4d,count=count4d))
        endif
      enddo

      !Time to do logp
      call nc_check(nf90_inq_varid(ncid,era_var_names(num_of_era_vars)%str,varid))

      allocate(grid2d(grid%inputxchunk,grid%inputychunk))

      if(grid%periodicboundary) then

        start3d = [integer:: grid%input_xstart,grid%input_ystart,timestep]
        count3d = [integer:: xgrid-grid%input_xstart+1,grid%inputychunk,1]

        call nc_check(nf90_get_var(ncid,varid,grid2d(1:(xgrid-(grid%input_xstart-1)),:),start=start3d,count=count3d))

        start3d = [integer:: 1,grid%input_ystart,timestep]
        count3d = [integer:: grid%input_xend,grid%inputychunk,1]

        call nc_check(nf90_get_var(ncid,varid,grid2d((xgrid-(grid%input_xstart-1))+1:grid%inputxchunk,:),start=start3d,count=count3d))

      else
        start3d = [integer:: grid%input_xstart,grid%input_ystart,timestep]
        count3d = [integer:: grid%inputxchunk,grid%inputychunk,1]

        call nc_check(nf90_get_var(ncid,varid,grid2d,start=start3d,count=count3d))
      endif
      call nc_check(nf90_close(ncid))
   end subroutine 
   
   subroutine write_truth_local_region_vert_level_mpi(grid,model_parameters,mpi_res,grid4d,grid3d,timestep,filename)
       use mpi

       use mod_utilities, only : mpi_type, grid_type, model_parameters_type
       use stringtype, only : string

       type(grid_type), intent(in)             :: grid
       type(model_parameters_type), intent(in) :: model_parameters
       type(mpi_type), intent(in)              :: mpi_res

       real(kind=dp), intent(inout)            :: grid4d(:,:,:,:), grid3d(:,:)
       integer, intent(in)                     :: timestep

       character(len=*), intent(in)            :: filename

       integer, parameter :: numdims4d=4, numdims3d=3, numspeedyvars=5
       integer :: dimsx,dimsy,dimsz

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d),start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims4d(numdims4d),arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(xgrid), lat(ygrid)

       real(kind=dp), allocatable :: copy(:,:,:,:)

       type(string) :: units(numspeedyvars)
       type(string) :: varname(numspeedyvars)

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'

       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'

       units(4)%str = 'g/kg'
       varname(4)%str = 'Specific_Humidity'

       units(5)%str = 'log(surfacepressure)'
       varname(5)%str = 'logp'
    
       !copy data
       allocate(copy,source=grid4d)

       dimsx = size(grid4d,2)
       dimsy = size(grid4d,3)
       dimsz = size(grid4d,4)

       varcount4d = [integer:: dimsx, dimsy, dimsz, 1]
       start4d =  [integer:: grid%res_xstart, grid%res_ystart, grid%res_zstart, timestep ]

       varcount3d =  [integer:: dimsx, dimsy, 1]
       start3d =  [integer:: grid%res_xstart, grid%res_ystart, timestep ]

       lon = (/(real(counter)*3.75,counter=0,95)/)
       lat = (/-87.159, -83.479, -79.777,  -76.070, -72.362, -68.652, &
               -64.942, -61.232, -57.521, -53.810, -50.099, -46.389, -42.678,-38.967, &
               -35.256, -31.545 , -27.833, -24.122, -20.411, -16.700, -12.989, -9.278, &
               -5.567, -1.856, 1.856, 5.567, 9.278, 12.989, 16.700, 20.411, &
               24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389, 50.099, &
               53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, &
               83.479, 87.159/)

       if(timestep.eq.1) then
           !call nc_check(nf90_create(filename,IOR(NF90_NETCDF4,NF90_MPIIO),file_id,comm=mpi_res%mpi_world%MPI_VAL,info=MPI_INFO_NULL))
           call nc_check(nf90_create(filename,IOR(NF90_NETCDF4,NF90_MPIIO),file_id,comm=mpi_res%mpi_world,info=MPI_INFO_NULL))

             ! define the dimensions
           call nc_check(nf90_def_dim(file_id, 'Lon', xgrid, xdim_id))
           call nc_check(nf90_def_dim(file_id, 'Lat', ygrid, ydim_id))
           call nc_check(nf90_def_dim(file_id, 'Sigma_Level', zgrid, zdim_id))
           call nc_check(nf90_def_dim(file_id, 'Timestep', int(model_parameters%predictionlength/model_parameters%timestep + 1,kind=int32), timedim_id))

           !Assign lat and lon ids and units
           call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
           call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

           call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
           call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
           ! now that the dimensions are defined, we can define variables on
           ! them,...
           arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /)
           arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

           do i=1, numspeedyvars-1
              call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims4d,array_id))
              ! ...and assign units to them as an attribute

              call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

              call nc_check(nf90_enddef(file_id))
              
              !parallel io is more complicated
              ! Unlimited dimensions require collective writes
              call nc_check(nf90_var_par_access(file_id, array_id, nf90_collective))

              !Write out the values
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d, count=varcount4d))

              call nc_check(nf90_put_var(file_id, xvar_id, lon))

              call nc_check(nf90_put_var(file_id, yvar_id, lat))

              call nc_check(nf90_redef(file_id))
           enddo
           !Lets do logp
           call nc_check(nf90_def_var(file_id,varname(5)%str,NF90_REAL,arrdims3d,array_id))

           call nc_check(nf90_put_att(file_id, array_id, "units", units(5)%str))

           call nc_check(nf90_enddef(file_id))

           call nc_check(nf90_var_par_access(file_id,array_id,nf90_collective))
           !Write out the values

           if(grid%logp_bool) then
             call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d, count=varcount3d))
           endif 
           ! close; done
           call nc_check(nf90_close(file_id))
        else
           !call nc_check(nf90_open(filename,IOR(NF90_WRITE,NF90_MPIIO),file_id,comm=mpi_res%mpi_world%MPI_VAL,info=MPI_INFO_NULL))
           call nc_check(nf90_open(filename,IOR(NF90_WRITE,NF90_MPIIO),file_id,comm=mpi_res%mpi_world,info=MPI_INFO_NULL))

           do i=1, numspeedyvars-1
              call nc_check(nf90_inq_varid(file_id,varname(i)%str,array_id))
              call nc_check(nf90_var_par_access(file_id,array_id,nf90_collective))
              call nc_check(nf90_put_var(file_id, array_id, copy(i,:,:,:),start=start4d,count=varcount4d))
           enddo

           call nc_check(nf90_inq_varid(file_id,varname(5)%str,array_id))

           call nc_check(nf90_var_par_access(file_id,array_id,nf90_collective))

           call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d,count=varcount3d))

           call nc_check(nf90_close(file_id))
        endif

        deallocate(copy)
        return
   end subroutine write_truth_local_region_vert_level_mpi

   subroutine read_full_file_4d(grid4d,grid2d,timestep,filename)
      use stringtype, only    : string

      real(kind=dp), intent(inout)            :: grid4d(:,:,:,:), grid2d(:,:)
      integer, intent(in)                     :: timestep

      character(len=*), intent(in)            :: filename
    
      !Local netcdf stuff
      integer :: start4d(4), start3d(3)
      integer :: count4d(4), count3d(3)
      integer :: ncid, varid
      integer :: dimids(4), dim_length(4)
      integer :: i,j

      integer, parameter :: maxnum_of_dims=4
      integer, parameter :: num_of_era_vars=5 
      integer, parameter :: num_of_height_levels=8

      type(string) :: era_var_names(num_of_era_vars)

      era_var_names(1)%str = 'Temperature'
      era_var_names(2)%str = 'U-wind'
      era_var_names(3)%str = 'V-wind'
      era_var_names(4)%str = 'Specific_Humidity'
      era_var_names(5)%str = 'logp' 

      call nc_check(nf90_open(filename, nf90_nowrite, ncid))

      start4d = [integer:: 1,1,1,timestep]
      count4d = [integer:: xgrid,ygrid,num_of_height_levels,1]

      do i=1, num_of_era_vars - 1
         call nc_check(nf90_inq_varid(ncid,era_var_names(i)%str,varid))

         call nc_check(nf90_get_var(ncid,varid,grid4d(i,:,:,:)))
      enddo 

      call nc_check(nf90_inq_varid(ncid,era_var_names(num_of_era_vars)%str,varid))

      call nc_check(nf90_get_var(ncid,varid,grid2d(:,:)))

      call nc_check(nf90_close(ncid)) 
      return
    end subroutine

    subroutine read_3d_file_parallel(filename,varname,mpi_res,grid,var3d,start_time_arg,stride_arg,time_length)
      use mpi
      use netcdf

      use mod_utilities, only : speedy_data_type, mpi_type, grid_type
      use stringtype, only    : string

      type(mpi_type), intent(in)          :: mpi_res
      type(grid_type), intent(in)         :: grid

      character(len=*), intent(in)        :: filename
      character(len=*), intent(in)        :: varname

      real(kind=dp), allocatable, intent(inout) :: var3d(:,:,:)

      integer, intent(in), optional             :: start_time_arg, stride_arg, time_length

      !Local netcdf stuff
      integer :: start3d(3)
      integer :: count3d(3) 
      integer :: stride3d(3)
      integer :: ncid, varid
      integer :: dimids(3), dim_length(3)
      integer :: i,j,maxnum_of_dims
      integer :: start_time, stride 
 
      if(present(start_time_arg)) then
        start_time = start_time_arg
      else
        start_time = 1
      endif

      if(present(stride_arg)) then
        stride = stride_arg
      else
        stride = 1
      endif

      maxnum_of_dims = 3

      !call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world%MPI_VAL, info=MPI_INFO_NULL),message='nf90_open')
      call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world, info=MPI_INFO_NULL),message='nf90_open')

      call nc_check(nf90_inq_varid(ncid,varname,varid))

      call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

      do j=1,maxnum_of_dims
         call nc_check(nf90_inquire_dimension(ncid,dimids(j),len=dim_length(j)),message='nf90_inquire_dimension')
      enddo


      if(present(time_length)) then
        dim_length(3) = time_length
      endif

      call nc_check(nf90_inq_varid(ncid,varname,varid),message='nf90_inq_varid')

      allocate(var3d(grid%inputxchunk,grid%inputychunk,dim_length(3)/stride))

      if(grid%periodicboundary) then

        start3d = [integer:: grid%input_xstart,grid%input_ystart,start_time]
        count3d = [integer:: xgrid-grid%input_xstart+1,grid%inputychunk,dim_length(3)/stride]
        stride3d = [integer:: 1,1,stride]

        call nc_check(nf90_get_var(ncid,varid,var3d(1:(xgrid-(grid%input_xstart-1)),:,:),start=start3d,count=count3d,stride=stride3d))

        start3d = [integer:: 1,grid%input_ystart,start_time]
        count3d = [integer:: grid%input_xend,grid%inputychunk,dim_length(3)/stride]
        stride3d = [integer:: 1,1,stride]

        call nc_check(nf90_get_var(ncid,varid,var3d((xgrid-(grid%input_xstart-1))+1:grid%inputxchunk,:,:),start=start3d,count=count3d,stride=stride3d))

      else
        start3d = [integer:: grid%input_xstart,grid%input_ystart,start_time]
        count3d = [integer:: grid%inputxchunk,grid%inputychunk,dim_length(3)/stride]
        stride3d = [integer:: 1,1,stride]

        call nc_check(nf90_get_var(ncid,varid,var3d,start=start3d,count=count3d))
      endif
      call nc_check(nf90_close(ncid))
    end subroutine read_3d_file_parallel

    subroutine read_3d_file_parallel_res(filename,varname,mpi_res,grid,var3d,time_index)
      !Parallel IO routine to read a file containing a 3d variable (x,y,t)

      !This routine needs to be called by all processors even if some of the
      !processors dont use that data. (e.g. for sst variables some workers dont
      !need it but in order for parallel io to work in netcdf all processors
      !need to open the file)
      use mpi
      use netcdf

      use mod_utilities, only : speedy_data_type, mpi_type, grid_type
      use stringtype, only    : string

      type(mpi_type), intent(in)          :: mpi_res
      type(grid_type), intent(in)         :: grid

      character(len=*), intent(in)        :: filename
      character(len=*), intent(in)        :: varname

      real(kind=dp), allocatable, intent(inout) :: var3d(:,:,:)

      integer, intent(in), optional       :: time_index

      !Local netcdf stuff
      integer :: start3d(3)
      integer :: count3d(3)
      integer :: ncid, varid
      integer :: dimids(3), dim_length(3)
      integer :: i,j,maxnum_of_dims
      integer :: start_time_index

      maxnum_of_dims = 3

      call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world, info=MPI_INFO_NULL),message='nf90_open')

      call nc_check(nf90_inq_varid(ncid,varname,varid))

      call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

      do j=1,maxnum_of_dims
         call nc_check(nf90_inquire_dimension(ncid,dimids(j),len=dim_length(j)),message='nf90_inquire_dimension')
      enddo

      if(present(time_index)) then
        start_time_index = time_index
        dim_length(3) = 1
      else
        start_time_index = 1
      endif

      call nc_check(nf90_inq_varid(ncid,varname,varid),message='nf90_inq_varid')

      allocate(var3d(grid%resxchunk,grid%resychunk,dim_length(3)))

      start3d = [grid%res_xstart,grid%res_ystart,start_time_index]
      count3d = [grid%resxchunk,grid%resychunk,dim_length(3)]

      call nc_check(nf90_get_var(ncid,varid,var3d,start=start3d,count=count3d))

      call nc_check(nf90_close(ncid))
    end subroutine read_3d_file_parallel_res

    subroutine read_prediction_local_model_vert_level_mpi(grid,model_parameters,mpi_res,grid4d,grid2d,timestep,filename)
      use mpi
      use netcdf

      use mod_utilities, only : era_data_type, mpi_type, grid_type, model_parameters_type
      use stringtype, only    : string

      type(mpi_type), intent(in)                :: mpi_res
      type(grid_type), intent(in)               :: grid
      type(model_parameters_type), intent(in)   :: model_parameters

      integer, intent(in)           :: timestep

      real(kind=dp), allocatable, intent(out) :: grid4d(:,:,:,:), grid2d(:,:)

      character(len=*), intent(in)              :: filename


      !Local netcdf stuff
      integer :: start4d(4), start3d(3)
      integer :: count4d(4), count3d(3)
      integer :: ncid, varid
      integer :: dimids(4), dim_length(4)
      integer :: i,j

      integer, parameter :: maxnum_of_dims=4
      integer, parameter :: num_of_speedy_vars=5

      type(string) :: speedy_var_names(num_of_speedy_vars)

      speedy_var_names(1)%str = 'Temperature'
      speedy_var_names(2)%str = 'U-wind'
      speedy_var_names(3)%str = 'V-wind'
      speedy_var_names(4)%str = 'Specific_Humidity'
      speedy_var_names(5)%str = 'logp'

      !call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world%MPI_VAL, info=MPI_INFO_NULL))
      call nc_check(nf90_open(filename, IOR(NF90_NOWRITE, NF90_MPIIO), ncid,comm=mpi_res%mpi_world, info=MPI_INFO_NULL))
      
      allocate(grid4d(num_of_speedy_vars-1,grid%resxchunk,grid%resychunk,grid%reszchunk))
      allocate(grid2d(grid%resxchunk,grid%resychunk))

      do i=1, num_of_speedy_vars-1
         call nc_check(nf90_inq_varid(ncid,speedy_var_names(i)%str,varid))

         call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

         start4d = [integer:: grid%res_xstart,grid%res_ystart,grid%res_zstart,timestep]
         count4d = [integer:: grid%resxchunk,grid%resychunk,grid%reszchunk,1]

         call nc_check(nf90_get_var(ncid,varid,grid4d(i,:,:,:),start=start4d,count=count4d))
      enddo

      !Time to do logp
      call nc_check(nf90_inq_varid(ncid,speedy_var_names(num_of_speedy_vars)%str,varid))

      start3d = [integer:: grid%res_xstart,grid%res_ystart,timestep]
      count3d = [integer:: grid%resxchunk,grid%resychunk,1]

      call nc_check(nf90_get_var(ncid,varid,grid2d,start=start3d,count=count3d))

      call nc_check(nf90_close(ncid))
    end subroutine

    subroutine read_trained_res(reservoir,model_parameters,grid)
      use mod_utilities, only : reservoir_type,model_parameters_type,grid_type

      type(reservoir_type), intent(inout) :: reservoir
      type(model_parameters_type), intent(in) :: model_parameters
      type(grid_type), intent(inout)             :: grid

      character(len=:), allocatable :: file_path
      character(len=4) :: worker_char
      character(len=1) :: height_char
      character(len=:), allocatable :: filename

      integer :: ncid

      file_path = '/scratch/user/troyarcomano/ML_SPEEDY_WEIGHTS/'

      write(worker_char,'(i0.4)') reservoir%assigned_region
      print *, 'reservoir%assigned_region',reservoir%assigned_region
      write(height_char,'(i0.1)') grid%level_index

      filename = file_path//'worker_'//worker_char//'_level_'//height_char//'_'//trim(model_parameters%trial_name)//'.nc'

      print *, 'reservoir%assigned_region',reservoir%assigned_region, filename
      call nc_check(nf90_open(filename, nf90_nowrite, ncid))

      print *, 'reading win'
      call read_netcdf_2d_dp_opened('win',ncid,reservoir%win)

      print *, 'reading wout'
      call read_netcdf_2d_dp_opened('wout',ncid,reservoir%wout)

      print *, 'reading rows'
      call read_netcdf_1d_int_opened('rows',ncid,reservoir%rows)

      print *, 'reading cols'
      call read_netcdf_1d_int_opened('cols',ncid,reservoir%cols)

      call read_netcdf_1d_dp_opened('vals',ncid,reservoir%vals)

      call read_netcdf_1d_dp_opened('mean',ncid,grid%mean)

      call read_netcdf_1d_dp_opened('std',ncid,grid%std)

      call nc_check(nf90_close(ncid))

    end subroutine

    subroutine read_trained_ocean_res(reservoir,model_parameters,grid)
      use mod_utilities, only : reservoir_type,model_parameters_type,grid_type

      type(reservoir_type), intent(inout)     :: reservoir
      type(model_parameters_type), intent(in) :: model_parameters
      type(grid_type), intent(inout)          :: grid

      character(len=:), allocatable :: file_path
      character(len=4) :: worker_char
      character(len=1) :: height_char
      character(len=:), allocatable :: filename

      integer :: ncid

      file_path = '/scratch/user/troyarcomano/ML_SPEEDY_WEIGHTS/'

      write(worker_char,'(i0.4)') reservoir%assigned_region
      print *, 'reservoir%assigned_region',reservoir%assigned_region

      filename = file_path//'worker_'//worker_char//'_ocean_'//trim(model_parameters%trial_name)//'.nc'
      print *, 'filename',filename

      if(file_exists(filename)) then
        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        print *, 'reading ocean win'
        call read_netcdf_2d_dp_opened('win',ncid,reservoir%win)

        print *, 'reading ocean wout'
        call read_netcdf_2d_dp_opened('wout',ncid,reservoir%wout)

        print *, 'reading ocean rows'
        call read_netcdf_1d_int_opened('rows',ncid,reservoir%rows)

        print *, 'reading ocean cols'
        call read_netcdf_1d_int_opened('cols',ncid,reservoir%cols)

        call read_netcdf_1d_dp_opened('vals',ncid,reservoir%vals)

        call read_netcdf_1d_dp_opened('mean',ncid,grid%mean)

        call read_netcdf_1d_dp_opened('std',ncid,grid%std)

        call nc_check(nf90_close(ncid))

        reservoir%sst_bool_input = .True.
        reservoir%sst_bool_prediction = .True.
      else
        reservoir%sst_bool_input = .False.
        reservoir%sst_bool_prediction = .False.
      endif
    end subroutine

    subroutine read_netcdf_2d_dp_opened(varname,ncid,var)
        character(len=*), intent(in) :: varname
        integer, intent(in)          :: ncid

        real(kind=dp), allocatable, intent(out) :: var(:,:)

        !Parmeter
        integer, parameter :: numofdims = 2 !We can assume its a 2d variable

        !Local netcdf variables
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1),dim_length(2)))

        call nc_check(nf90_get_var(ncid,varid,var))
        return
     end subroutine

    subroutine read_netcdf_1d_dp_opened(varname,ncid,var)
        character(len=*), intent(in) :: varname
        integer, intent(in)          :: ncid

        real(kind=dp), allocatable, intent(out) :: var(:)

        !Parmeter
        integer, parameter :: numofdims = 1!We can assume its a 1d variable

        !Local netcdf variables
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1)))

        call nc_check(nf90_get_var(ncid,varid,var))
        return
     end subroutine

     subroutine read_netcdf_1d_int_opened(varname,ncid,var)
        character(len=*), intent(in) :: varname
        integer, intent(in)          :: ncid

        integer, allocatable, intent(out) :: var(:)

        !Parmeter
        integer, parameter :: numofdims = 1!We can assume its a 1d variable

        !Local netcdf variables
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1)))

        call nc_check(nf90_get_var(ncid,varid,var))
        return
     end subroutine
end module mod_io
