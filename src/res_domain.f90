module resdomain
   use mod_utilities, only : dp, grid_type, main_type, &
                             speedygridnum, xgrid, ygrid, zgrid, &
                             speedylat, reservoir_type, &
                             model_parameters_type
   
   implicit none
 
   integer, parameter :: one=1
   
   !Interface to the tiler function
   interface tileoverlapgrid
     module procedure tileoverlapgrid5d
     module procedure tileoverlapgrid4d
     module procedure tileoverlapgrid3d
     module procedure tileoverlapgrid2d
   end interface tileoverlapgrid
  
   interface tile_full_input_to_target_data
     module procedure tile_full_input_to_target_data2d
     module procedure tile_full_input_to_target_data1d
   end interface tile_full_input_to_target_data

   interface tile_full_input_to_target_data_ocean_model
     module procedure tile_full_input_to_target_data2d_ocean_model
     module procedure tile_full_input_to_target_data1d_ocean_model
   end interface tile_full_input_to_target_data_ocean_model

   contains

      subroutine processor_decomposition(model_parameters)
         type(model_parameters_type), intent(inout) :: model_parameters

         integer :: num_regions_per_processor, left_over
         integer :: i 

         num_regions_per_processor = model_parameters%number_of_regions/model_parameters%numprocs

         left_over = mod(model_parameters%number_of_regions,model_parameters%numprocs)

         if((model_parameters%irank >= left_over+1).and.(model_parameters%irank > 0)) then 
            allocate(model_parameters%region_indices(num_regions_per_processor))
            model_parameters%num_of_regions_on_proc = num_regions_per_processor
            do i=1,num_regions_per_processor
               model_parameters%region_indices(i) = num_regions_per_processor*model_parameters%irank + i - 1
            enddo 
         elseif(model_parameters%irank == 0) then 
            allocate(model_parameters%region_indices(num_regions_per_processor))
            model_parameters%num_of_regions_on_proc = num_regions_per_processor
            do i=1,num_regions_per_processor
               model_parameters%region_indices(i) = i - 1
            enddo 
         else 
            allocate(model_parameters%region_indices(num_regions_per_processor+1))
            model_parameters%num_of_regions_on_proc = num_regions_per_processor + 1
            do i=1,num_regions_per_processor
               model_parameters%region_indices(i) = num_regions_per_processor*model_parameters%irank + i - 1
            enddo
            model_parameters%region_indices(i) = model_parameters%number_of_regions - left_over + model_parameters%irank - 1
         endif 
            
      end subroutine       

      subroutine processor_decomposition_manual(proc_number,numprocs,number_of_regions,region_indices)
         integer, intent(in) :: proc_number, numprocs, number_of_regions
 
         integer, allocatable, intent(out) :: region_indices(:)

         integer :: num_regions_per_processor, left_over
         integer :: i

         num_regions_per_processor = number_of_regions/numprocs

         left_over = mod(number_of_regions,numprocs)

         if((proc_number >= left_over+1).and.(proc_number > 0)) then
            allocate(region_indices(num_regions_per_processor))
            do i=1,num_regions_per_processor
               region_indices(i) = num_regions_per_processor*proc_number + i - 1
            enddo
         elseif(proc_number == 0) then
            allocate(region_indices(num_regions_per_processor))
            do i=1,num_regions_per_processor
               region_indices(i) = i - 1
            enddo
         else
            allocate(region_indices(num_regions_per_processor+1))
            do i=1,num_regions_per_processor
               region_indices(i) = num_regions_per_processor*proc_number + i - 1
            enddo
            region_indices(i) = number_of_regions - left_over + proc_number - 1
         endif

      end subroutine

      subroutine initializedomain(num_regions,region_num,overlap,num_vert_levels,vert_level,vert_overlap,grid)
         !Initializes the object grid with the processors domain info
         type(grid_type), intent(inout) :: grid 

         integer, intent(in) :: num_regions,region_num,overlap
         integer, intent(in) :: num_vert_levels,vert_level,vert_overlap

         logical, parameter :: setpole=.True.
        
         call getxyresextent(num_regions,region_num,grid%res_xstart,grid%res_xend, grid%res_ystart, grid%res_yend, grid%resxchunk,grid%resychunk)
       
         call get_z_res_extent(num_vert_levels,vert_level,grid%res_zstart,grid%res_zend,grid%reszchunk)

         call getoverlapindices(num_regions,region_num,overlap,grid%input_xstart,grid%input_xend,grid%input_ystart,grid%input_yend,grid%inputxchunk,grid%inputychunk,grid%pole,grid%periodicboundary,setpole)
         call getoverlapindices_vert(num_vert_levels,vert_level,vert_overlap,grid%input_zstart,grid%input_zend,grid%inputzchunk,grid%top,grid%bottom,setpole)

         call get_trainingdataindices(num_regions,region_num,overlap,grid%tdata_xstart,grid%tdata_xend,grid%tdata_ystart,grid%tdata_yend)

         call get_trainingdataindices_vert(num_vert_levels,vert_level,vert_overlap,grid%tdata_zstart,grid%tdata_zend)
         
         grid%overlap = overlap
         grid%num_vert_levels = num_vert_levels
         grid%vert_overlap = vert_overlap
         grid%number_of_regions = num_regions
        
      end subroutine 

      subroutine getxyresextent(num_regions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk)
          !Get the ij indices on the global for the local processor of just the
          !data that is fitted by the reservoir 

          integer, intent(in) :: num_regions,region_num
          integer, intent(out) :: localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk
          integer :: cornerx,cornery
        
          call domaindecomposition(num_regions,localxchunk,localychunk)
          
          call getworkerlower_leftcorner(region_num,localychunk,cornerx,cornery)

          localres_xstart = cornerx*localxchunk + 1
          localres_xend = (cornerx+1)*localxchunk

          localres_ystart = cornery*localychunk + 1
          localres_yend = (cornery+1)*localychunk
          return
      end subroutine
     
      subroutine get_z_res_extent(num_vert_levels,vert_level,local_res_zstart,local_res_zend,local_reszchunk) 
          integer, intent(in) :: num_vert_levels,vert_level
          integer, intent(out) :: local_res_zstart,local_res_zend,local_reszchunk
          
          local_reszchunk = zgrid/num_vert_levels

          local_res_zstart = (vert_level - 1)*local_reszchunk + 1
          local_res_zend = (vert_level)*local_reszchunk

          return
      end subroutine
        
      subroutine getoverlapindices(numregions,region_num,overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinputxchunk,localinputychunk,localpole,localperiodicboundary,setflag)
        !Gets the ij indices in 2d on the global for the overlapping (input)
        !region of each sub-domain

        integer, intent(in) :: region_num, numregions, overlap
        integer, intent(out) :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinputxchunk,localinputychunk
        logical, intent(out) :: localpole, localperiodicboundary
        integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk
        logical, optional :: setflag
          
        call getxyresextent(numregions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk)
        localinputxchunk = localxchunk + 2*overlap !There will always be a set number from this because there is always a left and right worker 
        localinputychunk = localychunk + 2*overlap !Start with the max and work your way down if the conditions are meet


        localperiodicboundary = .False.
        localpole = .False.

        if((localres_xstart - overlap).lt.one) then
          localinput_xstart = xgrid - overlap + 1
          localperiodicboundary = .True.
        else
           localinput_xstart = localres_xstart - overlap
        endif 
        
        if((localres_xend + overlap).gt.xgrid) then 
           localinput_xend = overlap
           localperiodicboundary = .True.
        else 
           
           localinput_xend = overlap + localres_xend
        endif 

        if((localres_ystart - overlap).lt.one) then
          localinput_ystart = one !localres_ystart
          localinputychunk = localychunk + overlap + (localres_ystart - 1)!localinputychunk - overlap
          localpole = .True.
        else
           localinput_ystart = localres_ystart - overlap 
        endif

        if((localres_yend + overlap).gt.ygrid) then
           localinput_yend = ygrid !localres_yend
           localinputychunk = localychunk + overlap + (ygrid - localres_yend)  !localinputychunk - overlap
           localpole = .True.
        else
           localinput_yend = overlap + localres_yend
        endif
        return
      end subroutine 
     
      subroutine getoverlapindices_vert(num_vert_levels,vert_level,vert_overlap,input_zstart,input_zend,inputzchunk,top,bottom,setflag) 
        integer, intent(in)  :: num_vert_levels, vert_level, vert_overlap
        integer, intent(out) :: input_zstart,input_zend,inputzchunk
        logical, intent(out) :: top, bottom
        logical, optional :: setflag
         
        integer :: local_res_zstart,local_res_zend,local_reszchunk
 

        call get_z_res_extent(num_vert_levels,vert_level,local_res_zstart,local_res_zend,local_reszchunk)

        top = .false.
        bottom = .false.

     
        if(local_res_zstart == one) then 
          top = .true.
        endif 

        if(local_res_zend == zgrid) then
          bottom = .true.  
        endif 
 
        if((local_res_zstart - vert_overlap >= one).and.(local_res_zend + vert_overlap <= zgrid)) then
          input_zstart = local_res_zstart - vert_overlap
          input_zend = local_res_zend + vert_overlap
    
          inputzchunk = local_reszchunk + 2*vert_overlap
        elseif(local_res_zstart - vert_overlap < one) then
          input_zstart = one
          input_zend = local_res_zend + vert_overlap
  
          inputzchunk = local_reszchunk + vert_overlap + (local_res_zstart - 1) 

        elseif(local_res_zend + vert_overlap > zgrid) then
          input_zstart = local_res_zstart - vert_overlap
          input_zend = zgrid

          inputzchunk = local_reszchunk + vert_overlap + (zgrid - local_res_zend) !(local_res_zend - zgrid)

        else
          print *, 'something is wrong'
          print *, 'vert_level,local_res_zstart,local_res_zend,vert_overlap',vert_level,local_res_zstart,local_res_zend,vert_overlap
        endif 

        !if(inputzchunk > zgrid .or. input_zstart < 1 .or. input_zend > zgrid) then
        !   inputzchunk = min(inputzchunk,zgrid) 
        !   input_zstart = max(input_zstart,one)
        !   input_zend = min(input_zend,zgrid)
        !endif 
      end subroutine 
          
      subroutine domaindecomposition(numregions,factorx,factory)
        !Breaks the global in 2d into each grid point rectangles based 
        !on the number of regions (sub-domains)
        integer, intent(in) :: numregions
        integer, intent(out) :: factorx,factory
        integer :: i,n,factorMax,check

        n = speedygridnum/numregions
        factorMax = floor(SQRT(real(n)))

        do i=factorMax,0, -1
           if(MOD(ygrid,i).eq.0) then
              factory = i
              if(mod(n,factory).eq.0) then
                 factorx = n/factory
                 if(mod(xgrid,factorx).eq.0) then
                   exit
                 endif
              endif
           endif
        enddo
        return  
      end subroutine        
         
      subroutine getworkerlower_leftcorner(region_num,factory,row,col)
         !Gets lower left corner ij index for a specific region 

         integer, intent(in) :: region_num, factory
         integer, intent(out) :: row,col

         col = mod(region_num,ygrid/factory)
         row = floor(real(region_num)/(real(ygrid)/real(factory)))

         return 
      end subroutine 

      subroutine tileoverlapgrid5d(grid,numregions,region_num,overlap,localgrid)
         !tiler for the full grid and time dimension
         real(kind=dp), intent(in)               :: grid(:,:,:,:,:) 
         integer, intent(in)                     :: numregions,region_num,overlap
         real(kind=dp), intent(out), allocatable :: localgrid(:,:,:,:,:)

         integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk
         integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk
         integer :: i,counter
         integer :: gridshape(5)
         logical :: localpole, localperiodicboundary

         gridshape = shape(grid)
         call getxyresextent(numregions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk)
         call getoverlapindices(numregions,region_num,overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk,localpole,localperiodicboundary,.FALSE.)

         !Lets get the reservoir extent
         if(.NOT.localpole) then
           if(.NOT.localperiodicboundary) then
             allocate(localgrid,MOLD=grid(:,localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:,:))
             localgrid = grid(:,localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:,:)
           endif 
           if(localperiodicboundary) then 
              allocate(localgrid(gridshape(1),localinput_xchunk,localinput_ychunk,gridshape(4),gridshape(5)))
              !Now lets determine if localinput_xstart>localinput_xend because
              !if so we need to tile localgrid from localinput_xstart:xgrid
              !first
              if(localres_xend.gt.localinput_xend.or.localinput_xstart.gt.localres_xstart) then
                localgrid(:,1:(xgrid-(localinput_xstart-1)),:,:,:) = grid(:,localinput_xstart:xgrid,localinput_ystart:localinput_yend,:,:)
                localgrid(:,(xgrid-(localinput_xstart-1))+1:localinput_xchunk,:,:,:) = grid(:,1:localinput_xend,localinput_ystart:localinput_yend,:,:)
              else
                 localgrid = grid(:,localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:,:) 
              endif 
           endif
         elseif(localpole) then
            if(.NOT.localperiodicboundary) then
              allocate(localgrid(gridshape(1),localinput_xchunk,localinput_ychunk,gridshape(4),gridshape(5))) 
              localgrid = grid(:,localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:,:)
            elseif(localperiodicboundary) then
              allocate(localgrid(gridshape(1),localinput_xchunk,localinput_ychunk,gridshape(4),gridshape(5)))
              !Now lets determine if localinput_xstart>localinput_xend because
              !if so we need to tile localgrid from localinput_xstart:xgrid
              !first
              if(localres_xend.gt.localinput_xend.or.localinput_xstart.gt.localres_xstart) then
                 localgrid(:,1:(xgrid-(localinput_xstart-1)),:,:,:) = grid(:,localinput_xstart:xgrid,localinput_ystart:localinput_yend,:,:)
                 localgrid(:,(xgrid-(localinput_xstart-1))+1:localinput_xchunk,:,:,:) = grid(:,1:localinput_xend,localinput_ystart:localinput_yend,:,:)
              else
                 !This means something is wrong or there is no overlap
                 localgrid = grid(:,localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:,:)
              endif
           endif
         endif  
      end subroutine

      subroutine tileoverlapgrid4d(grid,numregions,region_num,overlap,num_vert_levels,vert_level,vert_overlap,localgrid)
         !tiler that gets the full grid and gives back the local grid when there
         !is overlap 
         !NOTE only for one time step 
         !NOTE Im almost 100% this works 

         real(kind=dp), intent(in) :: grid(:,:,:,:)
         integer, intent(in) :: numregions,region_num,overlap,num_vert_levels,vert_level,vert_overlap
         real(kind=dp), intent(out),allocatable :: localgrid(:,:,:,:)

         integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk
         integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk
 
         integer :: localinput_zstart,localinput_zend,localinput_zchunk
         integer :: localres_zstart,localres_zend,localreszchunk

         integer :: i,counter
         integer :: gridshape(4)
         logical :: localpole,localperiodicboundary
         logical :: top, bottom

         gridshape = shape(grid)

         call getxyresextent(numregions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk)

         call get_z_res_extent(num_vert_levels,vert_level,localres_zstart,localres_zend,localreszchunk)

         call getoverlapindices(numregions,region_num,overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk,localpole,localperiodicboundary,.FALSE.)

         call getoverlapindices_vert(num_vert_levels,vert_level,vert_overlap,localinput_zstart,localinput_zend,localinput_zchunk,top,bottom,.False.)


         !Lets get the reservoir extent
         if(.NOT.localpole) then
           if(.NOT.localperiodicboundary) then
             allocate(localgrid,source=grid(:,localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,localinput_zstart:localinput_zend))
           endif
           if(localperiodicboundary) then
              allocate(localgrid(gridshape(1),localinput_xchunk,localinput_ychunk,localinput_zchunk))
 
              !Now lets determine if localinput_xstart>localinput_xend because
              !if so we need to tile localgrid from localinput_xstart:xgrid
              !first
              if(localres_xend.gt.localinput_xend.or.localinput_xstart.gt.localres_xstart) then
                 localgrid(:,1:(xgrid-(localinput_xstart-1)),:,:) = grid(:,localinput_xstart:xgrid,localinput_ystart:localinput_yend,localinput_zstart:localinput_zend)
                 localgrid(:,(xgrid-(localinput_xstart-1))+1:localinput_xchunk,:,:) = grid(:,1:localinput_xend,localinput_ystart:localinput_yend,localinput_zstart:localinput_zend)
              else
                 !This means something is wrong or there is no overlap
                 print *,'This means something is wrong not pole',region_num
                 localgrid = grid(:,localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,localinput_zstart:localinput_zend)
              endif
           endif
         elseif(localpole) then
            if(.NOT.localperiodicboundary) then
              allocate(localgrid,source=grid(:,localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,localinput_zstart:localinput_zend))
            elseif(localperiodicboundary) then
              allocate(localgrid(gridshape(1),localinput_xchunk,localinput_ychunk,localinput_zchunk))
              
              !Now lets determine if localinput_xstart>localinput_xend because
              !if so we need to tile localgrid from localinput_xstart:xgrid
              !first
              if(localres_xend.gt.localinput_xend.or.localinput_xstart.gt.localres_xstart) then
                 localgrid(:,1:(xgrid-(localinput_xstart-1)),:,:) = grid(:,localinput_xstart:xgrid,localinput_ystart:localinput_yend,localinput_zstart:localinput_zend)
                 localgrid(:,(xgrid-(localinput_xstart-1))+1:localinput_xchunk,:,:) = grid(:,1:localinput_xend,localinput_ystart:localinput_yend,localinput_zstart:localinput_zend)
              else
                 !This means something is wrong or there is no overlap
                 print *,'This means something is wrong pole',region_num
                 localgrid = grid(:,localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,localinput_zstart:localinput_zend)
              endif
           endif
         endif
         return
      end subroutine
       
      subroutine tileoverlapgrid3d(speedygrid,numprocs,worker,overlap,localgrid)
         !tiler that gets the full grid and gives back the local grid when there
         !is overlap
         !NOTE 3d grid x y time

         real(kind=dp), intent(in) :: speedygrid(:,:,:)
         integer, intent(in) :: numprocs,worker,overlap
         real(kind=dp), intent(out),allocatable :: localgrid(:,:,:)

         integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk
         integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk
         integer :: i,counter
         integer :: speedygridshape(3)
         logical :: localpole,localperiodicboundary

         speedygridshape = shape(speedygrid)
         call getxyresextent(numprocs,worker,localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk)
         call getoverlapindices(numprocs,worker,overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk,localpole,localperiodicboundary,.FALSE.)
         !Lets get the reservoir extent
         if(.NOT.localpole) then
           if(.NOT.localperiodicboundary) then
             allocate(localgrid,MOLD=speedygrid(localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:))
             localgrid = speedygrid(localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:)
           endif
           if(localperiodicboundary) then
              allocate(localgrid(localinput_xchunk,localinput_ychunk,speedygridshape(3)))

              !Now lets determine if localinput_xstart>localinput_xend because
              !if so we need to tile localgrid from localinput_xstart xgrid
              !first
              if(localres_xend.gt.localinput_xend.or.localinput_xstart.gt.localres_xstart) then
                 localgrid(1:(xgrid-(localinput_xstart-1)),:,:) = speedygrid(localinput_xstart:xgrid,localinput_ystart:localinput_yend,:)
                 localgrid((xgrid-(localinput_xstart-1))+1:localinput_xchunk,:,:) = speedygrid(1:localinput_xend,localinput_ystart:localinput_yend,:)
              else
                 !This means something is wrong or there is no overlap
                 print *,'This means something is wrong not pole',worker
                 localgrid = speedygrid(localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:)
              endif
           endif
         elseif(localpole) then
            if(.NOT.localperiodicboundary) then
              allocate(localgrid(localinput_xchunk,localinput_ychunk,speedygridshape(3)))
              localgrid = speedygrid(localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:)
            elseif(localperiodicboundary) then
              allocate(localgrid(localinput_xchunk,localinput_ychunk,speedygridshape(3)))

              !Now lets determine if localinput_xstart>localinput_xend because
              !if so we need to tile localgrid from localinput_xstart:xgrid
              !first
              if(localres_xend.gt.localinput_xend.or.localinput_xstart.gt.localres_xstart) then
                 localgrid(1:(xgrid-(localinput_xstart-1)),:,:) = speedygrid(localinput_xstart:xgrid,localinput_ystart:localinput_yend,:)
                 localgrid((xgrid-(localinput_xstart-1))+1:localinput_xchunk,:,:) = speedygrid(1:localinput_xend,localinput_ystart:localinput_yend,:)
              else
                 !This means something is wrong or there is no overlap
                 print *,'This means something is wrong pole',worker
                 localgrid = speedygrid(localinput_xstart:localinput_xend,localinput_ystart:localinput_yend,:)
              endif
           endif
         endif
         return
      end subroutine

      subroutine tileoverlapgrid2d(speedygrid,numregions,region_num,overlap,localgrid)
         !tiler that gets the full grid and gives back the local grid when there
         !is overlap
         !NOTE 2d grid x y time

         real(kind=dp), intent(in) :: speedygrid(:,:)
         integer, intent(in) :: numregions,region_num,overlap
         real(kind=dp), intent(out),allocatable :: localgrid(:,:)

         integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk
         integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk
         integer :: i,counter
         integer :: speedygridshape(2)
         logical :: localpole,localperiodicboundary

         speedygridshape = shape(speedygrid)

         call getxyresextent(numregions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk)
         call getoverlapindices(numregions,region_num,overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk,localpole,localperiodicboundary,.FALSE.)

         !Lets get the reservoir extent
         if(.NOT.localpole) then
           if(.NOT.localperiodicboundary) then
             allocate(localgrid,source=speedygrid(localinput_xstart:localinput_xend,localinput_ystart:localinput_yend))
           endif
           if(localperiodicboundary) then
              allocate(localgrid(localinput_xchunk,localinput_ychunk))

              !Now lets determine if localinput_xstart>localinput_xend because
              !if so we need to tile localgrid from localinput_xstart:xgrid
              !first
              if(localres_xend.gt.localinput_xend.or.localinput_xstart.gt.localres_xstart) then
                 localgrid(1:(xgrid-(localinput_xstart-1)),:) = speedygrid(localinput_xstart:xgrid,localinput_ystart:localinput_yend)
                 localgrid((xgrid-(localinput_xstart-1))+1:localinput_xchunk,:) = speedygrid(1:localinput_xend,localinput_ystart:localinput_yend)
              else
                 !This means something is wrong or there is no overlap
                 print *,'This means something is wrong not pole',region_num
                 localgrid = speedygrid(localinput_xstart:localinput_xend,localinput_ystart:localinput_yend)
              endif
           endif
         elseif(localpole) then
            if(.NOT.localperiodicboundary) then
              allocate(localgrid(localinput_xchunk,localinput_ychunk))
              localgrid = speedygrid(localinput_xstart:localinput_xend,localinput_ystart:localinput_yend)
            elseif(localperiodicboundary) then
              allocate(localgrid(localinput_xchunk,localinput_ychunk))

              !Now lets determine if localinput_xstart>localinput_xend because
              !if so we need to tile localgrid from localinput_xstart:xgrid
              !first
              if(localres_xend.gt.localinput_xend.or.localinput_xstart.gt.localres_xstart) then
                 localgrid(1:(xgrid-(localinput_xstart-1)),:) = speedygrid(localinput_xstart:xgrid,localinput_ystart:localinput_yend)
                 localgrid((xgrid-(localinput_xstart-1))+1:localinput_xchunk,:) = speedygrid(1:localinput_xend,localinput_ystart:localinput_yend)
              else
                 !This means something is wrong or there is no overlap
                 print *,'This means something is wrong pole',region_num
                 localgrid = speedygrid(localinput_xstart:localinput_xend,localinput_ystart:localinput_yend)
              endif
           endif
         endif
         return
      end subroutine 

      subroutine get_trainingdataindices(num_regions,region_num,overlap,xstart,xend,ystart,yend)
         !Get the fitting data indices in memory space

         integer, intent(in) :: region_num, num_regions, overlap
         integer, intent(out) :: xstart,xend,ystart,yend 
         integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,local_resxchunk,local_resychunk
         integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk
         logical :: localpole,localperiodicboundary

         call getxyresextent(num_regions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,local_resxchunk,local_resychunk)
         call getoverlapindices(num_regions,region_num,overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk,localpole,localperiodicboundary,.FALSE.)       

         xstart = 1 + overlap
         xend =  localinput_xchunk - overlap 
         
         !print *,worker,xstart,xend,localinput_xchunk 
         if((localres_ystart - overlap).lt.one) then
          ystart = 1 + (localres_ystart - 1) 
          yend = localinput_ychunk - overlap
         elseif((localres_yend + overlap).gt.ygrid) then
          ystart = 1 + overlap
          yend = localinput_ychunk - (ygrid - localres_yend)
         else 
          ystart= 1 + overlap
          yend = localinput_ychunk - overlap
         endif
         return
      end subroutine 
     
      subroutine get_trainingdataindices_vert(num_vert_levels,vert_level,vert_overlap,zstart,zend)
         !Get the fitting data indices in memory space

         integer, intent(in) :: num_vert_levels,vert_level,vert_overlap
         integer, intent(out) :: zstart,zend
         integer :: localres_zstart,localres_zend,local_reszchunk
         integer :: localinput_zstart,localinput_zend,localinput_zchunk
         logical :: top,bottom

         call get_z_res_extent(num_vert_levels,vert_level,localres_zstart,localres_zend,local_reszchunk)

         call getoverlapindices_vert(num_vert_levels,vert_level,vert_overlap,localinput_zstart,localinput_zend,localinput_zchunk,top,bottom,.false.) 

         if((localres_zstart - vert_overlap).lt.one) then
          zstart = 1 + (localres_zstart - 1)
          zend = localinput_zchunk - vert_overlap
         elseif((localres_zend + vert_overlap).gt.zgrid) then
          zstart = 1 + vert_overlap
          zend = localinput_zchunk - (zgrid - localres_zend)
         else
          zstart= 1 + vert_overlap
          zend = localinput_zchunk - vert_overlap
         endif
         return
      end subroutine 

      subroutine tile_full_input_to_target_data2d(reservoir,grid,statevec,tiledstatevec)
         !Takes 2d input array and tiles it to a 2d target (res) array 
         !Second dimension is time. Dimension of statevec(reservoir_numinputs,time)
         real(kind=dp), intent(in)  :: statevec(:,:)
         type(reservoir_type)       :: reservoir
         type(grid_type)            :: grid

         real(kind=dp), allocatable, intent(out) :: tiledstatevec(:,:)
        
         real(kind=dp), allocatable :: temp5d(:,:,:,:,:), temp3d(:,:,:)
         integer :: i, j
          
         i = size(statevec,1)
         j = size(statevec,2)
          
         if(reservoir%logp_bool) then
           allocate(temp5d(reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,reservoir%local_heightlevels_input,j))
           allocate(temp3d(grid%inputxchunk,grid%inputychunk,j))
           allocate(tiledstatevec(reservoir%chunk_size_prediction,j))
  
           temp5d = reshape(statevec(1:reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*reservoir%local_heightlevels_input,:),(/reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,reservoir%local_heightlevels_input,j/))

           temp3d = reshape(statevec(grid%logp_start:grid%logp_end,:),(/grid%inputxchunk,grid%inputychunk,j/))

           tiledstatevec(1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res,:) = reshape(temp5d(:,grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,grid%tdata_zstart:grid%tdata_zend,:),(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res,j/))
           tiledstatevec(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk,:) = reshape(temp3d(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:),(/grid%resxchunk*grid%resychunk,j/))

           deallocate(temp5d)
           deallocate(temp3d)
         else 
           allocate(temp5d(reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,reservoir%local_heightlevels_input,j))
           allocate(tiledstatevec(reservoir%chunk_size_prediction,j))

           temp5d = reshape(statevec,(/reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,reservoir%local_heightlevels_input,j/))
           tiledstatevec = reshape(temp5d(:,grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,grid%tdata_zstart:grid%tdata_zend,:),(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res,j/))   
     
           deallocate(temp5d)
         endif 
        
         if(reservoir%precip_bool) then  
           allocate(temp3d(grid%inputxchunk,grid%inputychunk,j)) 

           temp3d = reshape(statevec(grid%precip_start:grid%precip_end,:),(/grid%inputxchunk,grid%inputychunk,j/))

           tiledstatevec(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk+1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk*2,:)  = reshape(temp3d(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:),(/grid%resxchunk*grid%resychunk,j/))
         endif 
            
         return 

      end subroutine 

      subroutine tile_full_input_to_target_data1d(reservoir,grid,statevec,tiledstatevec)
         !Takes 1d input array and tiles it to a 1d target (res) array 
         real(kind=dp), intent(in)  :: statevec(:)

         type(reservoir_type)       :: reservoir
         type(grid_type)            :: grid

         real(kind=dp), allocatable, intent(out) :: tiledstatevec(:)

         real(kind=dp), allocatable :: temp4d(:,:,:,:), temp2d(:,:)
         integer :: i, j

         i = size(statevec,1)

         if(reservoir%logp_bool) then
           allocate(temp4d(reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,reservoir%local_heightlevels_input))
           allocate(temp2d(grid%inputxchunk,grid%inputychunk))
           allocate(tiledstatevec(reservoir%chunk_size_prediction))

           temp4d = reshape(statevec(1:reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*reservoir%local_heightlevels_input),(/reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,reservoir%local_heightlevels_input/))

           temp2d = reshape(statevec(reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*reservoir%local_heightlevels_input+1:reservoir%reservoir_numinputs-reservoir%tisr_size_input),(/grid%inputxchunk,grid%inputychunk/))

           tiledstatevec(1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res) = reshape(temp4d(:,grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,grid%tdata_zstart:grid%tdata_zend),(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res/))
           tiledstatevec(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+1:reservoir%chunk_size_prediction) = reshape(temp2d(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend),(/grid%resxchunk*grid%resychunk/))

         else
           allocate(temp4d(reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,reservoir%local_heightlevels_input))
           allocate(tiledstatevec(reservoir%chunk_size_prediction))

           temp4d = reshape(statevec,(/reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,reservoir%local_heightlevels_input/))
           tiledstatevec = reshape(temp4d(1:reservoir%local_predictvars,grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,grid%tdata_zstart:grid%tdata_zend),(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res/))
         endif

         return

      end subroutine
     
      subroutine tile_full_input_to_target_data2d_ocean_model(reservoir,grid,statevec,tiledstatevec)
         !Takes 2d input array and tiles it to a 2d target (res) array
         !Second dimension is time. Dimension of
         !statevec(reservoir_numinputs,time)
         real(kind=dp), intent(in)  :: statevec(:,:)
         type(reservoir_type)       :: reservoir
         type(grid_type)            :: grid

         real(kind=dp), allocatable, intent(out) :: tiledstatevec(:,:)

         real(kind=dp), allocatable :: temp3d(:,:,:)
         integer :: i, j

         i = size(statevec,1)
         j = size(statevec,2)

         allocate(temp3d(grid%inputxchunk,grid%inputychunk,j))
         allocate(tiledstatevec(reservoir%chunk_size_prediction,j))

         temp3d = reshape(statevec(grid%sst_start:grid%sst_end,:),(/grid%inputxchunk,grid%inputychunk,j/))

         tiledstatevec = reshape(temp3d(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:),(/grid%resxchunk*grid%resychunk,j/))

         return

      end subroutine

      subroutine tile_full_input_to_target_data1d_ocean_model(reservoir,grid,statevec,tiledstatevec)
         !Takes 1d input array and tiles it to a 1d target (res) array
         real(kind=dp), intent(in)  :: statevec(:)

         type(reservoir_type)       :: reservoir
         type(grid_type)            :: grid

         real(kind=dp), allocatable, intent(out) :: tiledstatevec(:)

         real(kind=dp), allocatable :: temp4d(:,:,:,:), temp2d(:,:)
         integer :: i

         i = size(statevec,1)

         allocate(temp2d(grid%inputxchunk,grid%inputychunk))
         allocate(tiledstatevec(reservoir%chunk_size_prediction))

         temp2d = reshape(statevec(grid%sst_start:grid%sst_end),(/grid%inputxchunk,grid%inputychunk/))

         tiledstatevec = reshape(temp2d(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend),(/grid%resxchunk*grid%resychunk/))

         return
      end subroutine 
 
      subroutine tile_4d_and_logp_state_vec_res1d(reservoir,num_regions,statevec,worker,grid4d,grid2d)
         !Takes the 1d res vector and converts it to the 4d and 2d grid target
         !(res) grid 
         real(kind=dp), intent(in)        :: statevec(:)
         integer, intent(in)              :: worker,num_regions
         type(reservoir_type), intent(in) :: reservoir
         
         real(kind=dp), intent(inout) :: grid4d(:,:,:,:), grid2d(:,:)
         
         integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk
         integer :: length
    
         length = size(statevec,1)

         call getxyresextent(num_regions,worker,localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk)

         if(reservoir%logp_bool) then
           grid4d = reshape(statevec(1:reservoir%local_predictvars*localres_xchunk*localres_ychunk*reservoir%local_heightlevels_res),(/reservoir%local_predictvars,localres_xchunk,localres_ychunk,reservoir%local_heightlevels_res/))
           grid2d = reshape(statevec(reservoir%local_predictvars*localres_xchunk*localres_ychunk*reservoir%local_heightlevels_res+1:length),(/localres_xchunk,localres_ychunk/))
         else  
           grid4d = reshape(statevec,(/reservoir%local_predictvars,localres_xchunk,localres_ychunk,reservoir%local_heightlevels_res/))
           grid2d = 0 
         endif 
        
         return
      end subroutine 
     
      subroutine tile_full_grid_with_local_state_vec_res1d(model_parameters,region_num,vert_level,statevec,wholegrid4d,wholegrid2d,wholegrid_precip)
         !Takes the 1d res vector and converts it to the 4d and 2d grid target
         !(res) grid
         real(kind=dp), intent(in)        :: statevec(:)
         integer, intent(in)              :: region_num, vert_level
         type(model_parameters_type), intent(in) :: model_parameters

         real(kind=dp), intent(inout) :: wholegrid4d(:,:,:,:), wholegrid2d(:,:), wholegrid_precip(:,:)

         integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk
         integer :: localres_zstart,localres_zend,localres_zchunk
         integer :: atmo3d_length
         integer :: length

         length = size(statevec,1)

         call getxyresextent(model_parameters%number_of_regions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk)         
 
         call get_z_res_extent(model_parameters%num_vert_levels,vert_level,localres_zstart,localres_zend,localres_zchunk)

         if(localres_zend == zgrid) then
           wholegrid4d(:,localres_xstart:localres_xend,localres_ystart:localres_yend,localres_zstart:localres_zend) = reshape(statevec(1:model_parameters%full_predictvars*localres_xchunk*localres_ychunk*localres_zchunk),(/model_parameters%full_predictvars,localres_xchunk,localres_ychunk,localres_zchunk/))
           wholegrid2d(localres_xstart:localres_xend,localres_ystart:localres_yend) = reshape(statevec(model_parameters%full_predictvars*localres_xchunk*localres_ychunk*localres_zchunk+1:model_parameters%full_predictvars*localres_xchunk*localres_ychunk*localres_zchunk+localres_xchunk*localres_ychunk),(/localres_xchunk,localres_ychunk/))
           if(model_parameters%precip_bool) then
             atmo3d_length = model_parameters%full_predictvars*localres_xchunk*localres_ychunk*localres_zchunk
             wholegrid_precip(localres_xstart:localres_xend,localres_ystart:localres_yend) = reshape(statevec(atmo3d_length+localres_xchunk*localres_ychunk+1:length),(/localres_xchunk,localres_ychunk/))
           endif 
         else
           wholegrid4d(:,localres_xstart:localres_xend,localres_ystart:localres_yend,localres_zstart:localres_zend) = reshape(statevec,(/model_parameters%full_predictvars,localres_xchunk,localres_ychunk,localres_zchunk/))
           wholegrid2d(localres_xstart:localres_xend,localres_ystart:localres_yend) = 0
         endif

        

         return
      end subroutine

      subroutine tile_full_2d_grid_with_local_res(model_parameters,region_num,statevec,wholegrid2d)
         !Takes the 1d res vector and converts it to the 4d and 2d grid target
         !(res) grid
         real(kind=dp), intent(in)        :: statevec(:)
         integer, intent(in)              :: region_num
         type(model_parameters_type), intent(in) :: model_parameters

         real(kind=dp), intent(inout) :: wholegrid2d(:,:)

         integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk
         integer :: localres_zstart,localres_zend,localres_zchunk
         integer :: length

         length = size(statevec,1)

         call getxyresextent(model_parameters%number_of_regions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk)

         !print *, 'localres_xstart:localres_xend,localres_ystart:localres_yend',localres_xstart,localres_xend,localres_ystart,localres_yend
         !print *, 'statevec(1:length)',statevec(1:length)
         wholegrid2d(localres_xstart:localres_xend,localres_ystart:localres_yend) = reshape(statevec(1:length),(/localres_xchunk,localres_ychunk/))

         return
      end subroutine
 
      subroutine tile_4d_and_logp_state_vec_input1d(reservoir,grid,statevec,grid4d,grid2d)
         !Tiler that takes a 1d local state vector and returns the 4d grid and
         !1d input grid
         real(kind=dp), intent(in)        :: statevec(:)

         type(reservoir_type), intent(in) :: reservoir
         type(grid_type), intent(in)      :: grid

         real(kind=dp), allocatable, intent(inout) :: grid4d(:,:,:,:), grid2d(:,:)
         
         integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk
         integer :: localinput_zstart,localinput_zend,localinput_zchunk

         integer :: length

         logical :: localpole,localperiodicboundary,top,bottom
    
         length = size(statevec,1)

         call getoverlapindices_vert(grid%num_vert_levels,grid%level_index,grid%vert_overlap,localinput_zstart,localinput_zend,localinput_zchunk,top,bottom)
         call getoverlapindices(grid%number_of_regions,reservoir%assigned_region,grid%overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk,localpole,localperiodicboundary)
 
         if(.not.allocated(grid4d)) then
           allocate(grid4d(reservoir%local_predictvars,localinput_xchunk,localinput_ychunk,localinput_zchunk))
         endif 

         if(.not.allocated(grid2d)) then
           allocate(grid2d(localinput_xchunk,localinput_ychunk))
         endif

         if(reservoir%logp_bool) then
           grid4d = reshape(statevec(grid%atmo3d_start:grid%atmo3d_end),(/reservoir%local_predictvars,localinput_xchunk,localinput_ychunk,localinput_zchunk/))
           grid2d = reshape(statevec(grid%logp_start:grid%logp_end),(/localinput_xchunk,localinput_ychunk/))
         else
           grid4d = reshape(statevec,(/reservoir%local_predictvars,localinput_xchunk,localinput_ychunk,localinput_zchunk/))
           grid2d = 0
         endif
         return
      end subroutine

      subroutine tile_4d_and_logp_state_vec_input1d_precip(reservoir,grid,statevec,grid4d,grid2d,grid_precip)
         !Tiler that takes a 1d local state vector and returns the 4d grid and
         !1d input grid
         real(kind=dp), intent(in)        :: statevec(:)

         type(reservoir_type), intent(in) :: reservoir
         type(grid_type), intent(in)      :: grid

         real(kind=dp), allocatable, intent(inout) :: grid4d(:,:,:,:), grid2d(:,:), grid_precip(:,:)

         integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk
         integer :: localinput_zstart,localinput_zend,localinput_zchunk

         integer :: length

         logical :: localpole,localperiodicboundary,top,bottom

         length = size(statevec,1)

         call getoverlapindices_vert(grid%num_vert_levels,grid%level_index,grid%vert_overlap,localinput_zstart,localinput_zend,localinput_zchunk,top,bottom)
         call getoverlapindices(grid%number_of_regions,reservoir%assigned_region,grid%overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk,localpole,localperiodicboundary)

         if(.not.allocated(grid4d)) then
           allocate(grid4d(reservoir%local_predictvars,localinput_xchunk,localinput_ychunk,localinput_zchunk))
         endif

         if(.not.allocated(grid2d)) then
           allocate(grid2d(localinput_xchunk,localinput_ychunk))
         endif

         if(.not.allocated(grid2d)) then
           allocate(grid_precip(localinput_xchunk,localinput_ychunk))
         endif

         if(reservoir%logp_bool) then
           grid4d = reshape(statevec(grid%atmo3d_start:grid%atmo3d_end),(/reservoir%local_predictvars,localinput_xchunk,localinput_ychunk,localinput_zchunk/))
           grid2d = reshape(statevec(grid%logp_start:grid%logp_end),(/localinput_xchunk,localinput_ychunk/))
         else
           grid4d = reshape(statevec,(/reservoir%local_predictvars,localinput_xchunk,localinput_ychunk,localinput_zchunk/))
           grid2d = 0
         endif

         if(reservoir%precip_bool) then 
           grid_precip = reshape(statevec(grid%precip_start:grid%precip_end),(/localinput_xchunk,localinput_ychunk/))
         endif 
         return
      end subroutine

      subroutine tile_4d_and_logp_state_vec_input1d_global(model_parameters,region_num,vert_level,statevec,grid4d,grid2d)
         !Tiler that takes a 1d local state vector and returns the 4d grid and
         !1d input grid
         !For a specific reservoir it can by called by any processor just takes
         !a little longer than using tile_4d_and_logp_state_vec_input1d
         real(kind=dp), intent(in)        :: statevec(:)

         integer, intent(in)              :: region_num,vert_level

         type(model_parameters_type), intent(in)      :: model_parameters

         real(kind=dp), allocatable, intent(inout) :: grid4d(:,:,:,:), grid2d(:,:)

         integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk
         integer :: localinput_zstart,localinput_zend,localinput_zchunk

         integer :: length

         logical :: localpole,localperiodicboundary,top,bottom

         length = size(statevec,1)

         call getoverlapindices_vert(model_parameters%num_vert_levels,vert_level,model_parameters%vert_loc_overlap,localinput_zstart,localinput_zend,localinput_zchunk,top,bottom)
         call getoverlapindices(model_parameters%number_of_regions,region_num,model_parameters%overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinput_xchunk,localinput_ychunk,localpole,localperiodicboundary)

         if(.not.allocated(grid4d)) then
           allocate(grid4d(model_parameters%full_predictvars,localinput_xchunk,localinput_ychunk,localinput_zchunk))
         endif

         if(.not.allocated(grid2d)) then
           allocate(grid2d(localinput_xchunk,localinput_ychunk))
         endif

         if(bottom) then
           grid4d = reshape(statevec(1:model_parameters%full_predictvars*localinput_xchunk*localinput_ychunk*localinput_zchunk),(/model_parameters%full_predictvars,localinput_xchunk,localinput_ychunk,localinput_zchunk/))
           grid2d = reshape(statevec(model_parameters%full_predictvars*localinput_xchunk*localinput_ychunk*localinput_zchunk+1:length),(/localinput_xchunk,localinput_ychunk/))
         else
           grid4d = reshape(statevec,(/model_parameters%full_predictvars,localinput_xchunk,localinput_ychunk,localinput_zchunk/))
           grid2d = 0
         endif
         return
      end subroutine

      subroutine tile_4d_and_logp_res_state_vec_res1d(reservoir,numregions,region_num,statevec,grid4d,grid2d)
         !Tiler that takes a 1d local (res) state vector and returns the 4d grid and
         !2d res grid
         real(kind=dp), intent(in)        :: statevec(:)

         integer, intent(in)              :: numregions,region_num

         type(reservoir_type), intent(in)        :: reservoir

         real(kind=dp), allocatable, intent(inout) :: grid4d(:,:,:,:), grid2d(:,:)

         !local variables
         integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk

         integer :: length

         logical :: localpole,localperiodicboundary

         call getxyresextent(numregions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk)
         length = size(statevec,1)

         if(.not.allocated(grid4d)) then
           allocate(grid4d(reservoir%local_predictvars,localres_xchunk,localres_ychunk,reservoir%local_heightlevels_res))
         endif

         if(.not.allocated(grid2d)) then
           allocate(grid2d(localres_xchunk,localres_ychunk))
         endif

         if(reservoir%logp_bool) then
           grid4d = reshape(statevec(1:reservoir%local_predictvars*localres_xchunk*localres_ychunk*reservoir%local_heightlevels_res),(/reservoir%local_predictvars,localres_xchunk,localres_ychunk,reservoir%local_heightlevels_res/))
           grid2d = reshape(statevec(reservoir%local_predictvars*localres_xchunk*localres_ychunk*reservoir%local_heightlevels_res+1:length),(/localres_xchunk,localres_ychunk/))
         else
           grid4d = reshape(statevec,(/reservoir%local_predictvars,localres_xchunk,localres_ychunk,reservoir%local_heightlevels_res/))
           grid2d = 0
         endif
         return
      end subroutine

      subroutine tile_4d_and_logp_full_grid_to_local_res_vec(model_parameters,region_num,vert_level,grid4d,grid2d,statevec)
        !Tiler that takes the global 4d and 2d grid to local state res vec
        real(kind=dp), intent(in)        :: grid4d(:,:,:,:), grid2d(:,:)
        integer, intent(in)              :: vert_level, region_num
        type(model_parameters_type), intent(in) :: model_parameters

        real(kind=dp), intent(out)    :: statevec(:)

        !local variables
        integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk 
        integer :: localres_zstart,localres_zend,localres_zchunk
        integer :: numvars
        
        real(kind=dp), allocatable :: temp4d(:,:,:,:), temp2d(:,:)

        numvars = size(grid4d,1)

        call getxyresextent(model_parameters%number_of_regions,region_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localres_xchunk,localres_ychunk)
        call get_z_res_extent(model_parameters%num_vert_levels,vert_level,localres_zstart,localres_zend,localres_zchunk)

        allocate(temp4d(numvars,localres_xchunk,localres_ychunk,localres_zchunk),temp2d(localres_xchunk,localres_ychunk))

        temp4d = grid4d(:,localres_xstart:localres_xend,localres_ystart:localres_yend,localres_zstart:localres_zend)
        temp2d = grid2d(localres_xstart:localres_xend,localres_ystart:localres_yend)

        statevec(1:numvars*localres_xchunk*localres_ychunk*localres_zchunk) = reshape(temp4d,(/numvars*localres_xchunk*localres_ychunk*localres_zchunk/))

        if(localres_zend == zgrid) then 
          statevec(numvars*localres_xchunk*localres_ychunk*localres_zchunk+1:numvars*localres_xchunk*localres_ychunk*localres_zchunk+localres_xchunk*localres_ychunk) = reshape(temp2d,(/localres_xchunk*localres_ychunk/))
        endif 
        
      end subroutine 

      subroutine tile_4d_and_logp_state_vec_input_to_local_grids(model_parameters,statevec,region_num,vert_level,grid4d,grid2d)
         !Tiler that takes a 1d local state vector and returns the 4d grid and
         !2d grid but only for the data with out any overlap
         real(kind=dp), intent(in)        :: statevec(:)
         integer, intent(in)              :: region_num, vert_level
         type(model_parameters_type), intent(in) :: model_parameters

         real(kind=dp), intent(inout) :: grid4d(:,:,:,:), grid2d(:,:)
    
         real(kind=dp), allocatable   :: temp4d(:,:,:,:), temp2d(:,:)
         
         integer :: local_tdata_xstart,local_tdata_xend,local_tdata_ystart,local_tdata_yend
         integer :: local_tdata_zstart,local_tdata_zend

         call tile_4d_and_logp_state_vec_input1d_global(model_parameters,region_num,vert_level,statevec,temp4d,temp2d)

         call get_trainingdataindices(model_parameters%number_of_regions,region_num,model_parameters%overlap,local_tdata_xstart,local_tdata_xend,local_tdata_ystart,local_tdata_yend)

         call get_trainingdataindices_vert(model_parameters%num_vert_levels,vert_level,model_parameters%vert_loc_overlap,local_tdata_zstart,local_tdata_zend)

         grid4d = temp4d(:,local_tdata_xstart:local_tdata_xend,local_tdata_ystart:local_tdata_yend,local_tdata_zstart:local_tdata_zend)
         grid2d = temp2d(local_tdata_xstart:local_tdata_xend,local_tdata_ystart:local_tdata_yend) 
         
         return 
      end subroutine 
         
      subroutine tile_4d_and_logp_to_local_state_input(model_parameters,region_num,vert_level,grid4d,grid2d,precip_grid,inputvec)
         !Takes the 4d and 2d grids and makes an input vector
         type(model_parameters_type), intent(in) :: model_parameters

         real(kind=dp), intent(in)        :: grid4d(:,:,:,:), grid2d(:,:), precip_grid(:,:)

         integer, intent(in)              :: region_num, vert_level

         real(kind=dp) , intent(inout)    :: inputvec(:)

         real(kind=dp), allocatable :: localgrid4d(:,:,:,:)
         real(kind=dp), allocatable :: localgrid2d(:,:)
         real(kind=dp), allocatable :: localprecipgrid(:,:)

         integer :: numvars, x, y, z, veclength
         integer :: localres_zstart,localres_zend,localreszchunk

         call get_z_res_extent(model_parameters%num_vert_levels,vert_level,localres_zstart,localres_zend,localreszchunk)

         call tileoverlapgrid(grid4d,model_parameters%number_of_regions,region_num,model_parameters%overlap,model_parameters%num_vert_levels,vert_level,model_parameters%vert_loc_overlap,localgrid4d)

         if(localres_zend == zgrid) then
           call tileoverlapgrid(grid2d,model_parameters%number_of_regions,region_num,model_parameters%overlap,localgrid2d)
           if(model_parameters%precip_bool) then
              call tileoverlapgrid(precip_grid,model_parameters%number_of_regions,region_num,model_parameters%overlap,localprecipgrid)
           endif 
         endif 
          
         numvars = size(localgrid4d,1)
         x = size(localgrid4d,2)
         y = size(localgrid4d,3)
         z = size(localgrid4d,4)

         veclength = size(inputvec)

         inputvec(1:numvars*x*y*z) = reshape(localgrid4d,(/numvars*x*y*z/))
         if(localres_zend == zgrid) then
            inputvec(numvars*x*y*z+1:numvars*x*y*z+x*y) = reshape(localgrid2d,(/x*y/))
            if(model_parameters%precip_bool) then
              inputvec(numvars*x*y*z+x*y+1:numvars*x*y*z+x*y*2) = reshape(localprecipgrid,(/x*y/))
            endif 
         endif 
         
         return 
      end subroutine 

      subroutine tile_4d_and_logp_to_local_state_input_slab(model_parameters,region_num,grid2d,inputvec)
        !Takes the 4d and 2d grids and makes an input vector
         type(model_parameters_type), intent(in) :: model_parameters

         real(kind=dp), intent(in)        :: grid2d(:,:)

         integer, intent(in)              :: region_num

         real(kind=dp) , intent(inout)    :: inputvec(:)

         real(kind=dp), allocatable :: localgrid2d(:,:)

         integer :: numvars, x, y, z, veclength
         integer :: localres_zstart,localres_zend,localreszchunk

         call tileoverlapgrid(grid2d,model_parameters%number_of_regions,region_num,model_parameters%overlap,localgrid2d)

         x = size(localgrid2d,1)
         y = size(localgrid2d,2)

         veclength = size(inputvec)

         inputvec(1:x*y) = reshape(localgrid2d,(/x*y/))

         return
      end subroutine

      subroutine tile_4d_to_local_state_input(model_parameters,region_num,vert_level,grid4d,grid2d,inputvec)
         !Takes the 4d and 2d grids and makes an input vector
         type(model_parameters_type), intent(in) :: model_parameters

         real(kind=dp), intent(in)        :: grid4d(:,:,:,:), grid2d(:,:)

         integer, intent(in)              :: region_num, vert_level

         real(kind=dp) , intent(inout)    :: inputvec(:)

         real(kind=dp), allocatable :: localgrid4d(:,:,:,:)
         real(kind=dp), allocatable :: localgrid2d(:,:)

         integer :: numvars, x, y, z, veclength

         call tileoverlapgrid(grid4d,model_parameters%number_of_regions,region_num,model_parameters%overlap,model_parameters%num_vert_levels,vert_level,model_parameters%vert_loc_overlap,localgrid4d)


         numvars = size(localgrid4d,1)
         x = size(localgrid4d,2)
         y = size(localgrid4d,3)
         z = size(localgrid4d,4)

         veclength = size(inputvec)

         inputvec(1:numvars*x*y*z) = reshape(localgrid4d,(/numvars*x*y*z/))

         return
      end subroutine

      subroutine tile_4d_and_logp_to_local_state_res(worker,grid4d,grid2d,inputvec)
         !Takes the 4d and 2d grids and makes a res vector
         real(kind=dp), intent(in)        :: grid4d(:,:,:,:), grid2d(:,:)
         integer, intent(in)              :: worker

         real(kind=dp) , intent(inout)    :: inputvec(:)

         real(kind=dp), allocatable :: localgrid4d(:,:,:,:)
         real(kind=dp), allocatable :: localgrid2d(:,:)

         integer :: numvars, x, y, z, veclength


         numvars = size(localgrid4d,1)
         x = size(localgrid4d,2)
         y = size(localgrid4d,3)
         z = size(localgrid4d,4)

         veclength = size(inputvec)

         inputvec(1:numvars*x*y*z) = reshape(localgrid4d,(/numvars*x*y*z/))
         inputvec(numvars*x*y*z+1:veclength) = reshape(localgrid2d,(/x*y/))

         return
      end subroutine

      
      subroutine standardize_state_vec_input(reservoir,grid,state_vec)
        !Subroutine to take state vector input and standardize the data
        !returns input state_vec standardized upon completetion

        type(reservoir_type)                    :: reservoir
        type(grid_type)                         :: grid

        real(kind=dp), intent(inout)     :: state_vec(:)

        !Local stuff
        real(kind=dp), allocatable :: temp4d(:,:,:,:), temp2d(:,:), temp_precip(:,:)

        integer :: num_of_standardized_vars, num_of_heights
        integer :: i, j, l
 
        call tile_4d_and_logp_state_vec_input1d(reservoir,grid,state_vec,temp4d,temp2d)

        call input_grid_to_input_statevec_and_standardization(reservoir,grid,temp4d,temp2d,state_vec)
      end subroutine

      subroutine input_grid_to_input_statevec_and_standardization(reservoir,grid,temp4d,temp2d,state_vec)
        !Subroutine to take state vector input and standardize the data
        !returns input state_vec standardized upon completetion
        use mod_utilities, only : standardize_data_given_pars3d, standardize_data_given_pars2d

        type(reservoir_type)                    :: reservoir
        type(grid_type)                         :: grid

        real(kind=dp), intent(inout)            :: state_vec(:)

        real(kind=dp), intent(inout)            :: temp4d(:,:,:,:), temp2d(:,:)

        !Local stuff
        integer :: num_of_standardized_vars, num_of_heights
        integer :: i, j, l

        num_of_standardized_vars = size(temp4d,1)
        num_of_heights = size(temp4d,4)

        l = 1
        do i=1, num_of_standardized_vars
           do j=1, num_of_heights
              call standardize_data_given_pars2d(temp4d(i,:,:,j),grid%mean(l),grid%std(l))
              l = l + 1
           enddo
        enddo

        if(reservoir%logp_bool) then
          call standardize_data_given_pars2d(temp2d,grid%mean(l),grid%std(l))
        endif 

        if(reservoir%logp_bool) then
           state_vec(grid%atmo3d_start:grid%atmo3d_end) = reshape(temp4d,(/reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk/))
           state_vec(grid%logp_start:grid%logp_end) = reshape(temp2d,(/grid%inputxchunk*grid%inputychunk/))
         else
           state_vec = reshape(temp4d,(/reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk/))
         endif
      end subroutine

      subroutine standardize_state_vec_res(reservoir,grid,state_vec)
        !Subroutine to take state vector input and standardize the data
        !returns input state_vec standardized upon completetion
        use mod_utilities, only : standardize_data_given_pars3d, standardize_data_given_pars2d

        type(reservoir_type), intent(in) :: reservoir
        type(grid_type), intent(in)      :: grid

        real(kind=dp), intent(inout)     :: state_vec(:)

        !Local stuff
        real(kind=dp), allocatable :: temp4d(:,:,:,:), temp2d(:,:)

        integer :: num_of_standardized_vars, num_of_heights
        integer :: i, j, l, data_height, length, height

        call tile_4d_and_logp_res_state_vec_res1d(reservoir,grid%number_of_regions,reservoir%assigned_region,state_vec,temp4d,temp2d)

        num_of_standardized_vars = size(temp4d,1)
        num_of_heights = size(temp4d,4)

        length = reservoir%local_predictvars
        height = reservoir%local_heightlevels_input
        l = 1
        do i=1, length
           data_height = 1
           do j=1, height
              if((j >= grid%tdata_zstart).and.(j <= grid%tdata_zend)) then
                 call standardize_data_given_pars2d(temp4d(i,:,:,data_height),grid%mean(l),grid%std(l))
                 data_height = data_height + 1
              endif
              l = l + 1
           enddo
        end do

        if(reservoir%logp_bool) then
          call standardize_data_given_pars2d(temp2d,grid%mean(l),grid%std(l))
        endif

        if(reservoir%logp_bool) then
           state_vec(1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res) = reshape(temp4d,(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res/))
           state_vec(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res + grid%resxchunk*grid%resychunk) = reshape(temp2d,(/grid%resxchunk*grid%resychunk/))
         else
           state_vec = reshape(temp4d,(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res/))
         endif
      end subroutine

      subroutine standardize_speedy_data(reservoir,grid,speedy_data)
        use mod_utilities, only : standardize_data_given_pars3d, speedy_data_type

        type(reservoir_type), intent(inout)   :: reservoir 
        type(speedy_data_type), intent(inout) :: speedy_data 
        type(grid_type), intent(inout)        :: grid

        integer       :: length, height, speedy_height
        integer       :: i, j, l

        length = reservoir%local_predictvars
        height = reservoir%local_heightlevels_input
        l = 1
        do i=1, length 
           speedy_height = 1 
           do j=1, height
              if((j >= grid%tdata_zstart).and.(j <= grid%tdata_zend)) then
                 call standardize_data_given_pars3d(speedy_data%speedyvariables(i,:,:,speedy_height,:),grid%mean(l),grid%std(l))
                 speedy_height = speedy_height + 1
              endif 
              l = l + 1
           enddo
        end do
 
        if(reservoir%logp_bool) then
          call standardize_data_given_pars3d(speedy_data%speedy_logp,grid%mean(l),grid%std(l)) 
        endif 
        return
      end subroutine

      subroutine standardize_speedy_data_gp_by_gp(reservoir,grid,speedy_data)
        use mod_utilities, only : standardize_data_given_pars3d, speedy_data_type

        type(reservoir_type), intent(inout)   :: reservoir
        type(speedy_data_type), intent(inout) :: speedy_data
        type(grid_type), intent(inout)        :: grid

        integer       :: vars, width, length, height, speedy_height, speedy_x, speedy_y
        integer       :: i, j, k, l, counter

        vars = reservoir%local_predictvars
        width = grid%inputxchunk
        length = grid%inputychunk
        height = reservoir%local_heightlevels_input
        counter = 1

        !TODO 
        do l=1, vars 
           speedy_height = 1   
           speedy_x = 1
           speedy_y = 1
           do i=1, width
                do j=1, length
                      do k=1, height
                         if((k >= grid%tdata_zstart).and.(k <= grid%tdata_zend).and.(i >= grid%tdata_xstart).and.(i <= grid%tdata_xend).and.(j >= grid%tdata_ystart).and.(j <= grid%tdata_yend)) then !TODO) then
                            call standardize_data_given_pars3d(speedy_data%speedyvariables(i,:,:,speedy_height,:),grid%mean(l),grid%std(l))
                            speedy_height = speedy_height + 1
                         endif
                        counter = counter + 1
                      enddo 
                enddo
            enddo 
        end do

        if(reservoir%logp_bool) then
          call standardize_data_given_pars3d(speedy_data%speedy_logp,grid%mean(l),grid%std(l))
        endif
        return
      end subroutine
 
      subroutine standardize_grid_res_tile_statevec(reservoir,grid,temp4d,temp2d,state_vec)
        !Subroutine to take state vector input and standardize the data
        !returns input state_vec standardized upon completetion
        use mod_utilities, only : standardize_data_given_pars3d, standardize_data_given_pars2d

        type(reservoir_type), intent(in) :: reservoir
        type(grid_type), intent(in)      :: grid

        real(kind=dp), intent(inout)     :: state_vec(:)

        real(kind=dp), intent(inout)     :: temp4d(:,:,:,:), temp2d(:,:)

        !local vars
        integer :: num_of_standardized_vars, num_of_heights
        integer :: i, j, l

        num_of_standardized_vars = size(temp4d,1)
        num_of_heights = size(temp4d,4)

        l = 1
        do i=1, num_of_standardized_vars
           do j=1, num_of_heights
              call standardize_data_given_pars2d(temp4d(i,:,:,j),grid%mean(l),grid%std(l))
              l = l + 1
           enddo
        enddo

        call standardize_data_given_pars2d(temp2d,grid%mean(l),grid%std(l))

        if(reservoir%logp_bool) then
           state_vec(1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res) = reshape(temp4d,(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res/))
           state_vec(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+1:reservoir%chunk_size_prediction) = reshape(temp2d,(/grid%resxchunk*grid%resychunk/))
         else
           state_vec = reshape(temp4d,(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res/))
         endif
      end subroutine

      subroutine unstandardize_state_vec_res(reservoir,grid,state_vec)
        !Subroutine to take state vector and unstandardize the data
        !returns state_vec unstandardized upon completetion
        use mod_utilities, only : unstandardize_data, unstandardize_data_2d, unstandardize_data_1d

        type(reservoir_type), intent(in) :: reservoir
        type(grid_type), intent(in)      :: grid

        real(kind=dp), intent(inout)     :: state_vec(:)

        !Local stuff
        real(kind=dp), allocatable :: temp4d(:,:,:,:), temp2d(:,:)

        integer :: num_of_unstandardized_vars, i, j, data_height, l, length, height

        allocate(temp4d(reservoir%local_predictvars,grid%resxchunk,grid%resychunk,reservoir%local_heightlevels_res))
        allocate(temp2d(grid%resxchunk,grid%resychunk))
     
        call tile_4d_and_logp_state_vec_res1d(reservoir,grid%number_of_regions,state_vec,reservoir%assigned_region,temp4d,temp2d)

        length = reservoir%local_predictvars
        height = reservoir%local_heightlevels_input
        l = 1
        do i=1, length
           data_height = 1
           do j=1, height
              if((j >= grid%tdata_zstart).and.(j <= grid%tdata_zend)) then
                 call unstandardize_data_2d(temp4d(i,:,:,data_height),grid%mean(l),grid%std(l))
                 data_height = data_height + 1
              endif
              l = l + 1
           enddo
        end do

        if(reservoir%logp_bool) then
          call unstandardize_data_2d(temp2d,grid%mean(grid%logp_mean_std_idx),grid%std(grid%logp_mean_std_idx))
        endif

        if(reservoir%logp_bool) then
           state_vec(1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res) = reshape(temp4d,(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res/))
           state_vec(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk) = reshape(temp2d,(/grid%resxchunk*grid%resychunk/))
        else
           state_vec = reshape(temp4d,(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res/))
        endif

        if(reservoir%precip_bool) then
          call unstandardize_data_1d(state_vec(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk+1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk*2),grid%mean(grid%precip_mean_std_idx),grid%std(grid%precip_mean_std_idx))          

          !state_vec(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk+1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk*2) = reshape(temp2d,(/grid%resxchunk*grid%resychunk/)) 
        endif
 
      end subroutine


      subroutine unstandardize_state_vec_res_and_tile_grids(reservoir,grid,state_vec,grid4d,grid2d)
        !Subroutine to take state vector and unstandardize the data
        !returns state_vec unstandardized upon completetion
        use mod_utilities, only : unstandardize_data

        type(reservoir_type), intent(in) :: reservoir
        type(grid_type), intent(in)      :: grid

        real(kind=dp), intent(inout)     :: state_vec(:)

        real(kind=dp), intent(out), allocatable :: grid4d(:,:,:,:), grid2d(:,:)
  
        allocate(grid4d(reservoir%local_predictvars,grid%resxchunk,grid%resychunk,reservoir%local_heightlevels_res))
        allocate(grid2d(grid%resxchunk,grid%resychunk))

        call tile_4d_and_logp_state_vec_res1d(reservoir,grid%number_of_regions,state_vec,reservoir%assigned_region,grid4d,grid2d)

        call unstandardize_data(reservoir,grid4d,grid2d,grid%mean,grid%std)

      end subroutine

      subroutine unstandardize_state_vec_input(reservoir,grid,state_vec)
        !Subroutine to take state vector and standardize the data
        !returns state_vec standardized upon completetion
        use mod_utilities, only : unstandardize_data

        type(reservoir_type), intent(in) :: reservoir
        type(grid_type), intent(in)      :: grid

        real(kind=dp), intent(inout)     :: state_vec(:)

        !Local stuff
        real(kind=dp), allocatable :: temp4d(:,:,:,:), temp2d(:,:)

        integer :: num_of_unstandardized_vars, i

        allocate(temp4d(reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,grid%inputzchunk))
        allocate(temp2d(grid%inputxchunk,grid%inputychunk))

        call tile_4d_and_logp_state_vec_input1d(reservoir,grid,state_vec,temp4d,temp2d)

        if(reservoir%logp_bool) then
           call unstandardize_data(reservoir,temp4d,temp2d,grid%mean,grid%std)
        else
           call unstandardize_data(reservoir,temp4d,grid%mean,grid%std)
        endif 

        if(reservoir%logp_bool) then
           state_vec(grid%atmo3d_start:grid%atmo3d_end) = reshape(temp4d,(/reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk/))
           state_vec(grid%logp_start:grid%logp_end) = reshape(temp2d,(/grid%inputxchunk*grid%inputychunk/))
        else
           state_vec = reshape(temp4d,(/reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk/))
         endif
      end subroutine

      subroutine unstandardize_state_vec_input_to_grid(reservoir,grid,state_vec,grid4d,grid2d)
        !Subroutine to take state vector and standardize the data
        !returns state_vec standardized upon completetion
        use mod_utilities, only : unstandardize_data

        type(reservoir_type), intent(in) :: reservoir
        type(grid_type), intent(in)      :: grid

        real(kind=dp), intent(inout)     :: state_vec(:)

        real(kind=dp), intent(out), allocatable :: grid4d(:,:,:,:), grid2d(:,:)

        !Local stuff
        real(kind=dp), allocatable :: temp4d(:,:,:,:), temp2d(:,:)

        integer :: num_of_unstandardized_vars, i

        allocate(temp4d(reservoir%local_predictvars,grid%inputxchunk,grid%inputychunk,grid%inputzchunk))
        allocate(temp2d(grid%inputxchunk,grid%inputychunk))

        call tile_4d_and_logp_state_vec_input1d(reservoir,grid,state_vec,temp4d,temp2d)

        call unstandardize_data(reservoir,temp4d,temp2d,grid%mean,grid%std)

        allocate(grid4d(reservoir%local_predictvars,grid%resxchunk,grid%resychunk,reservoir%local_heightlevels_res))
        allocate(grid2d(grid%resxchunk,grid%resychunk))
        
        grid4d = temp4d(:,grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,grid%tdata_zstart:grid%tdata_zend)
        grid2d = temp2d(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend)
      end subroutine

      subroutine set_region(grid)
        !Determines which region this worker is (tropics, extratropics, or
        !polar region)
        !This routine changes the local copy of grid so this can only be called
        !by the processor with is grid object and not from any other processor

        !We are calling tropics (25S to 25N)
        !Extratropics (25N to 60N and 25S to 60S)
        !Polar region (60N to 90N and 60S to 90S)
   
        type(grid_type), intent(inout)  :: grid
 
        real(kind=dp), parameter :: pole_lat_N = 60, pole_lat_S = -60
        real(kind=dp), parameter :: extra_lat_N = 30, extra_lat_S = -30

        real(kind=dp)            :: starting_lat, ending_lat

        !Lets determine whats the smallest and largest latitude that is inputted
        !into the reservoir
        starting_lat = speedylat(grid%res_ystart)
        ending_lat = speedylat(grid%res_yend)

        if((starting_lat <= pole_lat_S).or.(ending_lat >= pole_lat_N)) then
          grid%region_char = 'polar'
        elseif((starting_lat <= extra_lat_S).and.(starting_lat > pole_lat_S)) then
          grid%region_char = 'extratropic'
        elseif((starting_lat >= extra_lat_N).and.(starting_lat < pole_lat_N)) then
          grid%region_char = 'extratropic'
        else
          grid%region_char = 'tropic'
        endif 
        
      end subroutine
         
      subroutine set_reservoir_by_region(reservoir,grid)
        !Set reservoir parameters such as noise, spectral radius, etc by
        !geographic regions : tropics, extratropics, and polar regions
   
        type(reservoir_type), intent(inout) :: reservoir
        type(grid_type), intent(inout)      :: grid 
       
        call set_region(grid) 

        if(grid%region_char == 'tropic') then
          !res%specific_humidity_log_bool = .True.  
          reservoir%noisemag = 0.20!0
        elseif(grid%region_char == 'extratropic') then
          !res%specific_humidity_log_bool = .True.
          reservoir%noisemag = 0.20
        elseif(grid%region_char == 'polar') then
          !res%specific_humidity_log_bool = .False.
          reservoir%noisemag = 0.20
        else
          print *, 'something is wrong worker',reservoir%assigned_region,'doesnt have a region'
        endif 
  
        reservoir%radius = get_radius_by_lat(speedylat(grid%res_ystart),speedylat(grid%res_yend))
      end subroutine
  
      function get_radius_by_lat(startlat,endlat) result(radius)
        !First attempt at continuously increasing spectral radius from tropics
        !to extratropics. We keep the radius for areas above the highest_lat as
        !the poles and extra tropics seem to like spectral radius around 0.9 and
        !tropics like 0.3ish 
        real(kind=dp), intent(in) :: startlat, endlat
        real(kind=dp)             :: radius

        !descriptors to define how the function looks
        !Currently looks like this:
        !        r |    ___      ____  max_radius
        !        a |       \    /
        !        d |        \  /
        !        i |         \/
        !        u |      min_radius
        !        s ----------------------
        !            -90     0     90
        !                  lat
        real(kind=dp), parameter  :: highest_lat = 45.0_dp !Latitude where you want
                                                      !the spectral radius to be constant for any reservoir above this latitude  
        real(kind=dp), parameter  :: max_radius = 0.7_dp!0.7_dp !Max spectral radius 
        real(kind=dp), parameter  :: min_radius = 0.3_dp !Min spectral radius
   
        !Local stuff 
        real(kind=dp) :: smallest_lat
        real(kind=dp) :: largest_lat
   
        smallest_lat = abs(min(startlat,endlat))
        largest_lat  = abs(max(startlat,endlat))

        if(smallest_lat >= highest_lat) then
          radius = max_radius
        else
          radius = (max_radius - min_radius)/highest_lat + min_radius 
        endif 
        
        return
      end function get_radius_by_lat  

end module resdomain
