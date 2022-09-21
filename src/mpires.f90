module mpires
   use iso_fortran_env

   !use mpi_f08

   use mpi

   use mod_utilities, only : dp, main_type, mpi_type, state_vector_type, &
                             speedygridnum, xgrid, ygrid, &
                             zgrid, speedylat, model_parameters_type
   use resdomain, only : getxyresextent,getoverlapindices,tileoverlapgrid4d
   use mod_io, only : write_netcdf
   use mod_calendar, only : calendar, get_current_time_delta_hour
   implicit none
   
   type(mpi_type) :: mpi_res
   type(state_vector_type) :: internal_state_vector


   contains
     subroutine startmpi()
       !Your basic mpi starting routine for fortran

       call mpi_init(mpi_res%ierr)
       call mpi_comm_size(MPI_COMM_WORLD,mpi_res%numprocs,mpi_res%ierr)
       call mpi_comm_rank(MPI_COMM_WORLD, mpi_res%proc_num,mpi_res%ierr)
   
       mpi_res%mpi_world = MPI_COMM_WORLD
       
       if(mpi_res%proc_num == 0) then
          mpi_res%is_root = .true.
       endif 
   
       if(mpi_res%numprocs == 1) then
         mpi_res%is_serial = .True.
       endif  
     end subroutine

     subroutine killmpi()
       if (mpi_res%is_root) print *, 'cleaning up mpi'

       !call MPI_Abort(mpi_res%mpi_world,0, mpi_res%ierr)
       call mpi_finalize(mpi_res%ierr)
       if(mpi_res%is_root) then
         print *, 'killing program'
       endif 
       stop
     end subroutine 

     subroutine stop_mpi_safe()
       if(mpi_res%is_root) print *, 'safely cleaning up mpi'
       call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)

       call mpi_finalize(mpi_res%ierr)
       stop
     end subroutine

     subroutine predictionmpicontroller(res,timestep)
        !main routine to do all of the reading and writing for a single time
        !step prediction. Also runs the SPEEDY componement of the hybrid
      
        use mod_io, only : write_truth_local_region_vert_level_mpi, &
                           write_prediction_local_region_vert_level_mpi, &
                           read_prediction_local_region_vert_level_mpi, & 
                           read_prediction_local_model_vert_level_mpi, &
                           write_netcdf_speedy_full_mpi 
        use resdomain, only : unstandardize_state_vec_res_and_tile_grids, &
                              input_grid_to_input_statevec_and_standardization, &
                              unstandardize_state_vec_input_to_grid, &
                              standardize_grid_res_tile_statevec
        use mod_calendar 

        integer, intent(in)              :: timestep
        type(main_type), intent(inout)   :: res
 
        !local variables
        integer :: i, j

        real(kind=dp), allocatable :: grid4d(:,:,:,:), grid2d(:,:)

        !
        character(len=21) :: hybrid_out_root
        character(len=3)  :: file_end
        character(len=6)  :: trial_word
        character(len=2)  :: month
        character(len=4)  :: year
        character(len=2)  :: day
        character(len=2)  :: hour
        character(len=:), allocatable :: date_file
        character(len=:), allocatable :: hybrid_out_file_name

        character(len=9)  :: truth_out_root
        character(len=:), allocatable :: truth_out_file_name

        character(len=:), allocatable :: speedy_file
        
        character(len=:), allocatable :: file_path

        logical :: make_file

        file_path = '/scratch/user/troyarcomano/Predictions/Hybrid/'
        speedy_file = file_path//'hybrid_speedy_out.nc'
        hybrid_out_root='hybrid_prediction_era'
        truth_out_root = 'era_truth'
        trial_word = 'trial_'
        file_end = '.nc'


      
        call get_current_time_delta_hour(calendar,res%model_parameters%traininglength+res%model_parameters%synclength+res%model_parameters%prediction_markers(res%model_parameters%current_trial_number))!*res%model_parameters%timestep+timestep*res%model_parameters%timestep)

        write(year,'(I4.4)') calendar%currentyear
        write(month,'(I2.2)') calendar%currentmonth
        write(day,'(I2.2)') calendar%currentday
        write(hour,'(I2.2)') calendar%currenthour

        date_file = month//'_'//day//'_'//year//'_'//hour
        hybrid_out_file_name = file_path//hybrid_out_root//res%model_parameters%trial_name//trial_word//date_file//file_end

        truth_out_file_name = file_path//truth_out_root//res%model_parameters%trial_name//trial_word//date_file//file_end
   
        do i=1, res%model_parameters%num_of_regions_on_proc
           do j=1, res%model_parameters%num_vert_levels
              !print *, 'calling write_prediction_local_region_vert_level_mpi'
              call unstandardize_state_vec_res_and_tile_grids(res%reservoir(i,j),res%grid(i,j),res%reservoir(i,j)%outvec,grid4d,grid2d)
              !print *, 'grid4d',grid4d(:,1,1,1),'region',res%model_parameters%region_indices(i)
              if((timestep == 1).and.(i == 1).and.(j == 1)) then
                !print *, 'timestep',timestep,'i',i,'j',j
                make_file = .True.
              else
                make_file = .False.
              endif
              call write_prediction_local_region_vert_level_mpi(res%grid(i,j),res%model_parameters,mpi_res,grid4d,grid2d,timestep,hybrid_out_file_name,make_file)

              deallocate(grid4d)
              deallocate(grid2d)
           enddo 
        enddo  

        if(timestep == 1) then
          make_file = .True.
        else
          make_file = .False.
        endif 

        !if(make_file) then
          !call write_netcdf_speedy_full_mpi(timestep,res%model_parameters,speedy_file,mpi_res,make_file) 
        !endif 

        !if(mpi_res%is_root) then 
          !call run_model(res%model_parameters,timestep,hybrid_out_file_name,grid4d,grid2d)
        !endif 

        call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)

        !if(mpi_res%is_root) then
          !call write_netcdf_speedy_full_mpi(timestep,res%model_parameters,speedy_file,mpi_res,.False.,grid4d=grid4d,grid3d=grid2d)
        !  deallocate(grid4d)  
        !  deallocate(grid2d)
        !else
          !call write_netcdf_speedy_full_mpi(timestep,res%model_parameters,speedy_file,mpi_res,.False.)
        !endif 

        call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)

        do i=1, res%model_parameters%num_of_regions_on_proc
           do j=1, res%model_parameters%num_vert_levels
              call read_prediction_local_region_vert_level_mpi(res%grid(i,j),res%model_parameters,mpi_res,grid4d,grid2d,timestep,hybrid_out_file_name)
              call input_grid_to_input_statevec_and_standardization(res%reservoir(i,j),res%grid(i,j),grid4d,grid2d,res%reservoir(i,j)%feedback) 
              if(res%reservoir(i,j)%tisr_input_bool) then 
                call get_tisr_by_date(res%reservoir(i,j),res%grid(i,j),res%model_parameters,timestep-1,res%reservoir(i,j)%feedback(res%reservoir(i,j)%reservoir_numinputs-res%reservoir(i,j)%tisr_size_input+1:res%reservoir(i,j)%reservoir_numinputs))
              endif 
              deallocate(grid4d)
              deallocate(grid2d)
           enddo
        enddo

        call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)

        print *, 'read_prediction_local_model_vert_level_mpi'
        do i=1, res%model_parameters%num_of_regions_on_proc
           do j=1, res%model_parameters%num_vert_levels
              call read_prediction_local_model_vert_level_mpi(res%grid(i,j),res%model_parameters,mpi_res,grid4d,grid2d,timestep,speedy_file)
              call standardize_grid_res_tile_statevec(res%reservoir(i,j),res%grid(i,j),grid4d,grid2d,res%reservoir(i,j)%local_model)
              deallocate(grid4d)
              deallocate(grid2d)
           enddo
        enddo

        call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)

        print *, 'write_truth_local_region_vert_level_mpi'
        if(timestep == 1) then
          do i=1, res%model_parameters%num_of_regions_on_proc
             do j=1, res%model_parameters%num_vert_levels
                !Unstandardize here
                call unstandardize_state_vec_input_to_grid(res%reservoir(i,j),res%grid(i,j),res%reservoir(i,j)%predictiondata(:,res%model_parameters%synclength/res%model_parameters%timestep),grid4d,grid2d)
                call write_truth_local_region_vert_level_mpi(res%grid(i,j),res%model_parameters,mpi_res,grid4d,grid2d,timestep,truth_out_file_name)
             enddo
          enddo

          do i=1, res%model_parameters%num_of_regions_on_proc
             do j=1, res%model_parameters%num_vert_levels
                !Unstandardize here
                call unstandardize_state_vec_input_to_grid(res%reservoir(i,j),res%grid(i,j),res%reservoir(i,j)%predictiondata(:,res%model_parameters%synclength/res%model_parameters%timestep+timestep),grid4d,grid2d)
                call write_truth_local_region_vert_level_mpi(res%grid(i,j),res%model_parameters,mpi_res,grid4d,grid2d,timestep+1,truth_out_file_name)
             enddo
          enddo
        else!elseif(timestep == -99 ) then
          do i=1, res%model_parameters%num_of_regions_on_proc
             do j=1, res%model_parameters%num_vert_levels
                call unstandardize_state_vec_input_to_grid(res%reservoir(i,j),res%grid(i,j),res%reservoir(i,j)%predictiondata(:,res%model_parameters%synclength/res%model_parameters%timestep+timestep),grid4d,grid2d)
                call write_truth_local_region_vert_level_mpi(res%grid(i,j),res%model_parameters,mpi_res,grid4d,grid2d,timestep+1,truth_out_file_name)
             enddo
          enddo
        endif

     end subroutine 

     subroutine sendrecievegrid(res,timestep,ocean_model)
        use resdomain, only : tile_4d_and_logp_to_local_state_input, tile_4d_and_logp_state_vec_res1d, unstandardize_state_vec_input, tile_4d_and_logp_full_grid_to_local_res_vec, &
                              processor_decomposition_manual, getxyresextent, tile_full_grid_with_local_state_vec_res1d, standardize_state_vec_input, standardize_state_vec_res, &
                              tile_4d_to_local_state_input, tile_full_2d_grid_with_local_res, tile_4d_and_logp_to_local_state_input_slab

        use mod_utilities, only : unstandardize_data, standardize_data_given_pars3d, standardize_data_given_pars1d

        type(main_type), intent(inout) :: res

        integer, intent(in) :: timestep

        logical, intent(in) :: ocean_model 

        real(kind=dp), allocatable :: wholegrid4d(:,:,:,:), wholegrid2d(:,:)
        real(kind=dp), allocatable :: wholegrid_sst(:,:), wholegrid_precip(:,:)
        
        real(kind=dp), allocatable :: forecast_4d(:,:,:,:), forecast_2d(:,:)
        real(kind=dp), allocatable :: sendreceivedata(:), temp4d(:,:,:,:), temp2d(:,:), temp3d(:,:,:), temp1d(:)

        integer, parameter :: root=0
        integer :: i,j, recieverequest, sendrequest, local_domain_size, receive_size
        
        integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localresxchunk,localresychunk
        integer :: localres_zstart,localres_zend,localreszchunk
        integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinputxchunk,localinputychunk
        integer :: to, from, tag, tag2
        integer :: status(MPI_STATUS_SIZE)
        integer :: counter
        integer :: proc_num, number_of_regions
        integer :: full_grid_num_2ds
        
        integer, allocatable :: region_indices(:)

        logical, parameter :: setflag=.False.
        logical            :: localpole,localperiodicboundary
        logical, allocatable :: local_sst_flag(:)

        character(len=3)  :: file_end
        character(len=6)  :: trial_word
        character(len=2)  :: month
        character(len=4)  :: year
        character(len=2)  :: day
        character(len=2)  :: hour
        character(len=:), allocatable :: date_file
        character(len=:), allocatable :: hybrid_out_file_name
        character(len=:), allocatable :: file_path
        character(len=:), allocatable :: hybrid_out_root

        !The receiving part of the routine
        !Gets all of the outvecs from each worker
        !and gives it to the master node (worker == 0)
        !Master node reconstructs the whole global set vector and
        !then writes it out to the disk

        !print *, "starting receivesend", res%model_parameters%irank
        if(mpi_res%is_root) then
           allocate(wholegrid4d(res%model_parameters%full_predictvars,xgrid,ygrid,zgrid))
           allocate(wholegrid2d(xgrid,ygrid))

           allocate(forecast_4d(res%model_parameters%full_predictvars,xgrid,ygrid,zgrid))
           allocate(forecast_2d(xgrid,ygrid))

           if(ocean_model) then
             allocate(wholegrid_sst(xgrid,ygrid))  
             wholegrid_sst = res%model_parameters%base_sst_grid
             print *, 'allocated local_sst_flag'
             allocate(local_sst_flag(mpi_res%numprocs))
           endif

           if(res%model_parameters%precip_bool) then 
             allocate(wholegrid_precip(xgrid,ygrid))
             wholegrid_precip = 0.0_dp
           endif 

           !print *, 'root allocated full grids',res%model_parameters%irank
           wholegrid4d = 0
           wholegrid2d = 0

           do i=1, res%model_parameters%num_of_regions_on_proc
              do j=1, res%model_parameters%num_vert_levels

                 !print *, 'i,j,root,outvec',i,j,res%reservoir(i,j)%outvec
  
                 call tile_full_grid_with_local_state_vec_res1d(res%model_parameters,res%model_parameters%region_indices(i),j,res%reservoir(i,j)%outvec,wholegrid4d,wholegrid2d,wholegrid_precip)

              enddo 
              if(ocean_model) then
                if(res%reservoir_special(i,1)%sst_bool_prediction) then
                  call tile_full_2d_grid_with_local_res(res%model_parameters,res%model_parameters%region_indices(i),res%reservoir_special(i,1)%outvec,wholegrid_sst)
                endif 
              endif 
           enddo
        endif

        tag = 11
        tag2 = 12

        if(ocean_model) then
          if(.not. allocated(local_sst_flag)) allocate(local_sst_flag(1))

          !print *, 'shape(res%reservoir_special(:,1)%sst_bool_prediction)',shape(res%reservoir_special(:,1)%sst_bool_prediction)
          !print *, 'res%reservoir_special(:,1)%sst_bool_prediction',res%reservoir_special(:,1)%sst_bool_prediction
          call MPI_Gather(res%reservoir_special(1,1)%sst_bool_prediction, 1, MPI_LOGICAL,local_sst_flag, 1, MPI_LOGICAL, 0 ,mpi_res%mpi_world,mpi_res%ierr)
        endif 

        call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)

        counter = 1
        if(.not.(mpi_res%is_root)) then
           !print *, 'sending outvec',res%model_parameters%irank
           do i=1, res%model_parameters%num_of_regions_on_proc
              do j=1, res%model_parameters%num_vert_levels
                 local_domain_size = size(res%reservoir(i,j)%outvec)

                 allocate(sendreceivedata(local_domain_size))

                 sendreceivedata = res%reservoir(i,j)%outvec
                 to = root

                 tag = counter

                 call MPI_SEND(sendreceivedata,local_domain_size,MPI_DOUBLE_PRECISION,to,tag,mpi_res%mpi_world,mpi_res%ierr)

                 deallocate(sendreceivedata)

                 counter = counter + 1
              enddo 

              if(ocean_model) then 

                 if(res%reservoir_special(i,1)%sst_bool_prediction) then 
                    local_domain_size = size(res%reservoir_special(i,1)%outvec)

                    allocate(sendreceivedata(local_domain_size))

                    sendreceivedata = res%reservoir_special(i,1)%outvec
                 else 
                    local_domain_size = res%grid_special(i,1)%resxchunk * res%grid_special(i,1)%resychunk 

                    allocate(sendreceivedata(local_domain_size))

                    sendreceivedata = 272.0_dp
                 endif 
                 to = root

                 tag = counter
                 
                 !print *, 'worker',mpi_res%proc_num,'tag',tag,'sending sst_data',sendreceivedata
                 call MPI_SEND(sendreceivedata,local_domain_size,MPI_DOUBLE_PRECISION,to,tag,mpi_res%mpi_world,mpi_res%ierr)

                 deallocate(sendreceivedata)

                 counter = counter + 1
              endif 
           enddo  
        endif

        call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)
        
        if((mpi_res%is_root)) then
          do proc_num=1,mpi_res%numprocs-1

             !print *, 'root receiving from data from processor',proc_num
             call processor_decomposition_manual(proc_num,mpi_res%numprocs,res%model_parameters%number_of_regions,region_indices)
          
             number_of_regions = size(region_indices)

             counter = 1
             do i=1, number_of_regions
                do j=1, res%model_parameters%num_vert_levels

                    !print *, 'root i,j',i,j
                    call getsend_receive_size_res(res%model_parameters,region_indices(i),j,receive_size)

                    allocate(sendreceivedata(receive_size))

                    from = proc_num

                    tag = counter

                    !print *, 'from',from,'tag',tag
                    call MPI_RECV(sendreceivedata,receive_size,MPI_DOUBLE_PRECISION,from,tag,mpi_res%mpi_world,MPI_STATUS_IGNORE,mpi_res%ierr)

                    !print *, 'successfully received from',from,'tag',tag
                    call tile_full_grid_with_local_state_vec_res1d(res%model_parameters,region_indices(i),j,sendreceivedata,wholegrid4d,wholegrid2d,wholegrid_precip)

                    deallocate(sendreceivedata)

                    counter = counter + 1
                 enddo
                 if(ocean_model) then 
                    !print *, 'root i,j',i,j
                    call getsend_receive_size_res_slab(res%model_parameters,region_indices(i),j,receive_size)

                    allocate(sendreceivedata(receive_size))

                    from = proc_num

                    tag = counter

                    call MPI_RECV(sendreceivedata,receive_size,MPI_DOUBLE_PRECISION,from,tag,mpi_res%mpi_world,MPI_STATUS_IGNORE,mpi_res%ierr)

                    !print *, 'from',from,'tag',tag,'sst data',sendreceivedata
                    !print *, 'region_indices(i)',region_indices(i),i 
                    call tile_full_2d_grid_with_local_res(res%model_parameters,region_indices(i),sendreceivedata,wholegrid_sst)

                    deallocate(sendreceivedata)

                    counter = counter + 1
                 endif 
              enddo
           enddo
        endif

        if(mpi_res%is_root) then
           !Make sure we dont wander too far away from realistic values of
           !specific humidity

           where(wholegrid4d(4,:,:,:) < 0.000001)
                 wholegrid4d(4,:,:,:) = 0.000001
           endwhere

           !where(wholegrid4d(4,:,:,:) > 22)
           !      wholegrid4d(4,:,:,:) = 22
           !endwhere

           !print *, 'root starting writing hybrid prediction'

           if(ocean_model) then
             do i=1, xgrid
                do j=1, ygrid
                   if(res%model_parameters%sea_mask(i,j) > 0.0) then
                      wholegrid_sst(i,j) = res%model_parameters%base_sst_grid(i,j)
                   endif
                enddo
             enddo 
           endif 

           !NOTE TODO change back 
           if(ocean_model .and. .not. res%model_parameters%train_on_sst_anomalies) then
              where(wholegrid_sst < 272.0)
                 wholegrid_sst = 272.0_dp
              endwhere
           endif 

           if(res%model_parameters%precip_bool) then
             where(wholegrid_precip < 0.00001 )
                 wholegrid_precip = 0.0_dp
             endwhere 
           endif 

           if(.not. res%model_parameters%ml_only) then
             hybrid_out_root='hybrid_prediction_era'
           else
             hybrid_out_root='ml_prediction_era'
           endif 

           trial_word = 'trial_'
           file_end = '.nc'

           call get_current_time_delta_hour(calendar,res%model_parameters%traininglength+res%model_parameters%synclength+res%model_parameters%prediction_markers(res%model_parameters%current_trial_number))
           write(year,'(I4.4)') calendar%currentyear
           write(month,'(I2.2)') calendar%currentmonth
           write(day,'(I2.2)') calendar%currentday
           write(hour,'(I2.2)') calendar%currenthour

           file_path = '/scratch/user/troyarcomano/Predictions/Hybrid/'
           date_file = month//'_'//day//'_'//year//'_'//hour
           hybrid_out_file_name = file_path//hybrid_out_root//res%model_parameters%trial_name//trial_word//date_file//file_end

           if(timestep == 1) then 
             res%model_parameters%prediction_file = hybrid_out_file_name
           endif 
           if(timestep > 1) then
             hybrid_out_file_name = res%model_parameters%prediction_file
           endif  
           !Write the hybrid prediction out
           print *, 'writing hybrid to', hybrid_out_file_name
           print *, 'res%model_parameters%traininglength+res%model_parameters%synclength+res%model_parameters%prediction_markers(res%model_parameters%current_trial_number)',res%model_parameters%traininglength,res%model_parameters%synclength,res%model_parameters%prediction_markers(res%model_parameters%current_trial_number),res%model_parameters%current_trial_number
 
           full_grid_num_2ds = 1
           if(ocean_model) then
             full_grid_num_2ds = full_grid_num_2ds + 1
           endif

           if(res%model_parameters%precip_bool) then
             full_grid_num_2ds = full_grid_num_2ds + 1 
           endif 

           full_grid_num_2ds = 3 

           allocate(temp3d(full_grid_num_2ds,xgrid,ygrid))
         
           temp3d = 0.0_dp

           temp3d(1,:,:) = wholegrid2d

           if(ocean_model) then
             temp3d(2,:,:) = wholegrid_sst
           endif 

           if(res%model_parameters%precip_bool) then 
             temp3d(3,:,:) = wholegrid_precip
           endif 

           !print *, 'wholegrid_sst',wholegrid_sst
           if(ocean_model .or. res%model_parameters%precip_bool) then
              call write_netcdf(res%model_parameters,wholegrid4d,temp3d,timestep,hybrid_out_file_name,ocean_model)
           !elseif(res%model_parameters%precip_bool .and. .not. ocean_model) then 
           !   call write_netcdf(res%model_parameters,wholegrid4d,temp3d,timestep,hybrid_out_file_name)
           !elseif( .not. res%model_parameters%precip_bool .and. ocean_model) then 
           !   call write_netcdf_4d_multi_2d_sst_only(res%model_parameters,wholegrid4d,temp3d,timestep,hybrid_out_file_name)
           else 
              call write_netcdf(res%model_parameters,wholegrid4d,wholegrid2d,timestep,hybrid_out_file_name)
           endif 
           deallocate(temp3d)

           !Run SPEEDY in hybrid configuration 
           if(.not. res%model_parameters%ml_only) then
             print *, 'root running model'
             call run_model(res%model_parameters,timestep,wholegrid4d,wholegrid2d,wholegrid_sst,forecast_4d,forecast_2d)
           endif 
        endif

        if(timestep == 1) then
          !print *, 'calling write_truth_data', res%model_parameters%irank
          call write_truth_data(res,timestep)
          call write_truth_data(res,timestep+1)
        elseif(timestep == -99) then !else
          call write_truth_data(res,timestep+1)
        endif

      
        if(mpi_res%is_root) then
           do i=1, res%model_parameters%num_of_regions_on_proc
              do j=1, res%model_parameters%num_vert_levels

                 !print *, 'tile_4d_and_logp_to_local_state_input root'
                 
                 call tile_4d_and_logp_to_local_state_input(res%model_parameters,res%model_parameters%region_indices(i),j,wholegrid4d,wholegrid2d,wholegrid_precip,res%reservoir(i,j)%feedback(1:res%reservoir(i,j)%reservoir_numinputs-res%reservoir(i,j)%tisr_size_input-res%reservoir(i,j)%sst_size_input))

                 if(.not. res%model_parameters%ml_only) then
                   call tile_4d_and_logp_full_grid_to_local_res_vec(res%model_parameters,res%model_parameters%region_indices(i),j,forecast_4d,forecast_2d,res%reservoir(i,j)%local_model)

                   call standardize_state_vec_res(res%reservoir(i,j),res%grid(i,j),res%reservoir(i,j)%local_model)
                 endif 
              enddo
           enddo 

           do proc_num=1,mpi_res%numprocs-1 
              
              call processor_decomposition_manual(proc_num,mpi_res%numprocs,res%model_parameters%number_of_regions,region_indices)

              number_of_regions = size(region_indices) 

              counter = 1
              do i=1, number_of_regions
                 do j=1, res%model_parameters%num_vert_levels

                    call getsend_receive_size_input(res%model_parameters,region_indices(i),j,receive_size)
                    !print *, 'sending from root to',proc_num,'region num',region_indices(i),receive_size

                    allocate(sendreceivedata(receive_size))

                    call tile_4d_and_logp_to_local_state_input(res%model_parameters,region_indices(i),j,wholegrid4d,wholegrid2d,wholegrid_precip,sendreceivedata)

                    to = proc_num

                    tag = counter

                    call MPI_SEND(sendreceivedata,size(sendreceivedata),MPI_DOUBLE_PRECISION,to,tag,mpi_res%mpi_world,mpi_res%ierr)

                    counter = counter + 1
                    
                    deallocate(sendreceivedata)

                    !Only send speedy data if this is hybrid 
                    if(.not. res%model_parameters%ml_only) then
                      call getsend_receive_size_speedy(res%model_parameters,region_indices(i),j,receive_size)

                      allocate(sendreceivedata(receive_size))

                      call tile_4d_and_logp_full_grid_to_local_res_vec(res%model_parameters,region_indices(i),j,forecast_4d,forecast_2d,sendreceivedata)

                      tag = counter

                      call MPI_SEND(sendreceivedata,size(sendreceivedata),MPI_DOUBLE_PRECISION,to,tag,mpi_res%mpi_world,mpi_res%ierr)

                      deallocate(sendreceivedata)

                      counter = counter + 1
                    endif 

                 enddo
                 if(ocean_model) then
                    !print *, 'sending slab from root to',proc_num,'region num',region_indices(i)
                    call getsend_receive_size_input_slab(res%model_parameters,region_indices(i),receive_size)

                    allocate(sendreceivedata(receive_size))

                    call tile_4d_and_logp_to_local_state_input_slab(res%model_parameters,region_indices(i),wholegrid_sst,sendreceivedata)

                    to = proc_num

                    tag = counter

                    call MPI_SEND(sendreceivedata,size(sendreceivedata),MPI_DOUBLE_PRECISION,to,tag,mpi_res%mpi_world,mpi_res%ierr)

                    counter = counter + 1

                    deallocate(sendreceivedata)
                  endif 
               enddo 
           enddo 
        endif 
        
        call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)

        if(.not.(mpi_res%is_root)) then
          counter = 1
          do i=1, res%model_parameters%num_of_regions_on_proc
             do j=1, res%model_parameters%num_vert_levels
                call getsend_receive_size_input(res%model_parameters,res%model_parameters%region_indices(i),j,receive_size)

                allocate(sendreceivedata(receive_size))

                from = root

                tag = counter

                call MPI_RECV(sendreceivedata,receive_size,MPI_DOUBLE_PRECISION,from,tag,mpi_res%mpi_world,MPI_STATUS_IGNORE,mpi_res%ierr)

                res%reservoir(i,j)%feedback(1:res%reservoir(i,j)%reservoir_numinputs-res%reservoir(i,j)%tisr_size_input-res%reservoir(i,j)%sst_size_input) = sendreceivedata
                !print *, res%model_parameters%irank,'receiving feedback from root for region, level',res%model_parameters%region_indices(i),j,res%reservoir(i,j)%feedback(res%grid(i,j)%precip_start:res%grid(i,j)%precip_end)

                deallocate(sendreceivedata)

                counter = counter + 1 

                if(.not. res%model_parameters%ml_only) then
                  allocate(sendreceivedata(res%reservoir(i,j)%chunk_size_speedy))!prediction))

                  tag = counter

                  call MPI_RECV(sendreceivedata,res%reservoir(i,j)%chunk_size_speedy,MPI_DOUBLE_PRECISION,from,tag,mpi_res%mpi_world,MPI_STATUS_IGNORE,mpi_res%ierr) !prediction,MPI_DOUBLE_PRECISION,from,tag,mpi_res%mpi_world,MPI_STATUS_IGNORE,mpi_res%ierr)

                  res%reservoir(i,j)%local_model = sendreceivedata
                  deallocate(sendreceivedata)
    
                  call standardize_state_vec_res(res%reservoir(i,j),res%grid(i,j),res%reservoir(i,j)%local_model) 
                  counter = counter + 1
                endif 
              enddo

              if(ocean_model) then
                call getsend_receive_size_input_slab(res%model_parameters,res%model_parameters%region_indices(i),receive_size)

                allocate(sendreceivedata(receive_size))
                allocate(temp1d(receive_size)) !Holds input sst data for this
                                               !process from here until the end of the routine  

                from = root

                tag = counter

                !print *, res%model_parameters%irank,'receiving slab feedback from root for region, level',res%model_parameters%region_indices(i),j
                call MPI_RECV(sendreceivedata,receive_size,MPI_DOUBLE_PRECISION,from,tag,mpi_res%mpi_world,MPI_STATUS_IGNORE,mpi_res%ierr)

                if(res%reservoir_special(i,1)%sst_bool_prediction) then
                   call standardize_data_given_pars1d(sendreceivedata,res%grid_special(i,1)%mean(res%grid_special(i,1)%sst_mean_std_idx),res%grid_special(i,1)%std(res%grid_special(i,1)%sst_mean_std_idx))

                   res%reservoir_special(i,1)%feedback(res%grid_special(i,1)%sst_start:res%grid_special(i,1)%sst_end) = sendreceivedata
                endif 

                deallocate(sendreceivedata)
                deallocate(temp1d)

                counter = counter + 1
              endif 
           enddo
        endif

        !Distribute the logical flag to let all workers know if we are going to
        !try to make the next prediction

        call MPI_Bcast(res%model_parameters%run_speedy,1,MPI_LOGICAl,0,mpi_res%mpi_world,mpi_res%ierr)

        !Stupid parallel netcdf requires all works so we cannot do it during the
        !sending and receiving loops because the root will not be able to
        !call get_tisr_by_date
        do i=1, res%model_parameters%num_of_regions_on_proc
           do j=1, res%model_parameters%num_vert_levels
              if(res%reservoir(i,j)%tisr_input_bool) then
                 call get_tisr_by_date(res%reservoir(i,j),res%grid(i,j),res%model_parameters,timestep-1,res%reservoir(i,j)%feedback(res%grid(i,j)%tisr_start:res%grid(i,j)%tisr_end))
                 if(res%reservoir(i,j)%assigned_region == 36) print *, 'tisr in prediction,level',res%grid(i,j)%level_index,res%reservoir(i,j)%feedback(res%grid(i,j)%tisr_start)
              endif
 
              if(ocean_model) then 
                 if(res%reservoir(i,j)%sst_bool_input) then
                    res%reservoir(i,j)%feedback(res%grid(i,j)%sst_start:res%grid(i,j)%sst_end) = res%reservoir_special(i,1)%feedback(res%grid_special(i,1)%sst_start:res%grid_special(i,1)%sst_end) 
                 endif 
              endif 

              if(res%reservoir(i,j)%logp_bool) then
                call standardize_state_vec_input(res%reservoir(i,j),res%grid(i,j),res%reservoir(i,j)%feedback(1:res%grid(i,j)%logp_end))
              else 
                call standardize_state_vec_input(res%reservoir(i,j),res%grid(i,j),res%reservoir(i,j)%feedback(1:res%grid(i,j)%atmo3d_end))
              endif 

              if(res%reservoir(i,j)%precip_bool) then
                call standardize_data_given_pars1d(res%reservoir(i,j)%feedback(res%grid(i,j)%precip_start:res%grid(i,j)%precip_end),res%grid(i,j)%mean(res%grid(i,j)%precip_mean_std_idx),res%grid(i,j)%std(res%grid(i,j)%precip_mean_std_idx))
              endif 
                
           enddo
           if(ocean_model) then
              if(res%reservoir_special(i,1)%sst_bool_prediction) then
                 res%reservoir_special(i,1)%averaged_atmo_input_vec(:,mod(timestep-1,res%model_parameters%timestep_slab/res%model_parameters%timestep-1)+1) = res%reservoir(i,j-1)%feedback(res%reservoir_special(i,1)%atmo_training_data_idx)
                 if(res%reservoir_special(i,1)%assigned_region == 10) print *,'averaged_atmo_input_vec(1,:)',res%reservoir_special(i,1)%averaged_atmo_input_vec(1,:)
                 res%reservoir_special(i,1)%feedback = sum(res%reservoir_special(i,1)%averaged_atmo_input_vec,dim=2)/(res%model_parameters%timestep_slab/res%model_parameters%timestep-1) !res%reservoir(i,j-1)%feedback(res%reservoir_special(i,1)%atmo_training_data_idx)
                 if(res%reservoir_special(i,1)%assigned_region == 10) print *,'res%reservoir_special(i,1)%feedback(1)',res%reservoir_special(i,1)%feedback(1)
              endif
           endif 
        enddo 
        return
     end subroutine 
           
     subroutine getsend_receive_size_res(model_parameters,region,vert_level,arraysize)
       use resdomain, only : get_z_res_extent, getxyresextent

       type(model_parameters_type), intent(in) :: model_parameters

       integer, intent(in)  :: region, vert_level
       integer, intent(out) :: arraysize

       !local stuff
       integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localresxchunk,localresychunk
       integer :: localres_zstart,localres_zend,localreszchunk

       call getxyresextent(model_parameters%number_of_regions,region,localres_xstart,localres_xend,localres_ystart,localres_yend,localresxchunk,localresychunk)

       call get_z_res_extent(model_parameters%num_vert_levels,vert_level,localres_zstart,localres_zend,localreszchunk)

       arraysize = localresxchunk*localresychunk*localreszchunk*model_parameters%full_predictvars

       if(localres_zend == zgrid) then
         arraysize = arraysize + localresxchunk*localresychunk
         if(model_parameters%precip_bool) then
           arraysize = arraysize + localresxchunk*localresychunk
         endif
       endif 
     end subroutine

     subroutine getsend_receive_size_speedy(model_parameters,region,vert_level,arraysize)
       use resdomain, only : get_z_res_extent, getxyresextent

       type(model_parameters_type), intent(in) :: model_parameters

       integer, intent(in)  :: region, vert_level
       integer, intent(out) :: arraysize

       !local stuff
       integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localresxchunk,localresychunk
       integer :: localres_zstart,localres_zend,localreszchunk

       call getxyresextent(model_parameters%number_of_regions,region,localres_xstart,localres_xend,localres_ystart,localres_yend,localresxchunk,localresychunk)

       call get_z_res_extent(model_parameters%num_vert_levels,vert_level,localres_zstart,localres_zend,localreszchunk)

       arraysize = localresxchunk*localresychunk*localreszchunk*model_parameters%full_predictvars

       if(localres_zend == zgrid) then
         arraysize = arraysize + localresxchunk*localresychunk
       endif
     end subroutine

     subroutine getsend_receive_size_res_slab(model_parameters,region,vert_level,arraysize)
       use resdomain, only : get_z_res_extent, getxyresextent

       type(model_parameters_type), intent(in) :: model_parameters

       integer, intent(in)  :: region, vert_level
       integer, intent(out) :: arraysize

       !local stuff
       integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localresxchunk,localresychunk
       integer :: localres_zstart,localres_zend,localreszchunk

       call getxyresextent(model_parameters%number_of_regions,region,localres_xstart,localres_xend,localres_ystart,localres_yend,localresxchunk,localresychunk)

       arraysize = localresxchunk*localresychunk

     end subroutine 

     subroutine getsend_receive_size_input(model_parameters,region,vert_level,arraysize)
       use resdomain, only : getoverlapindices_vert, getoverlapindices

       type(model_parameters_type), intent(in) :: model_parameters
       integer, intent(in)  :: region, vert_level
       integer, intent(out) :: arraysize

       !local stuff
       integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinputxchunk,localinputychunk
       integer :: localinput_zstart,localinput_zend,localinputzchunk
       logical :: localpole,localperiodicboundary,localtop,localbottom

       call getoverlapindices(model_parameters%number_of_regions,region,model_parameters%overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinputxchunk,localinputychunk,localpole,localperiodicboundary,.False.)

       call getoverlapindices_vert(model_parameters%num_vert_levels,vert_level,model_parameters%vert_loc_overlap,localinput_zstart,localinput_zend,localinputzchunk,localtop,localbottom,.False.)

       arraysize = localinputxchunk*localinputychunk*localinputzchunk*model_parameters%full_predictvars 

       if(localbottom) then
         arraysize = arraysize + localinputxchunk*localinputychunk
         if(model_parameters%precip_bool) then
           arraysize = arraysize + localinputxchunk*localinputychunk
         endif 
       endif 
       !print *, 'arraysize',arraysize
     end subroutine

     subroutine getsend_receive_size_input_slab(model_parameters,region,arraysize)
       use resdomain, only : getoverlapindices_vert, getoverlapindices

       type(model_parameters_type), intent(in) :: model_parameters
       integer, intent(in)  :: region
       integer, intent(out) :: arraysize

       !local stuff
       integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinputxchunk,localinputychunk
       integer :: localinput_zstart,localinput_zend,localinputzchunk
       logical :: localpole,localperiodicboundary,localtop,localbottom

       call getoverlapindices(model_parameters%number_of_regions,region,model_parameters%overlap,localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinputxchunk,localinputychunk,localpole,localperiodicboundary,.False.)


       arraysize = localinputxchunk*localinputychunk

     end subroutine

               
     subroutine distribute_prediction_marker(model_parameters)
       !Routine to get random prediction marker and distribute to the network
       use mod_utilities, only : init_random_seed, shuffle, model_parameters_type

       type(model_parameters_type), intent(inout) :: model_parameters

       integer :: i, markertag, status, root
       integer :: cycle_length !cylc forecast ever x number of hours
     
       markertag = 2
       root = 0
       cycle_length = model_parameters%synclength 

       allocate(model_parameters%prediction_markers(model_parameters%num_predictions))

       do i=0,model_parameters%num_predictions-1
          model_parameters%prediction_markers(i+1) = cycle_length*i
       enddo
     
       return 
     end subroutine 

     subroutine write_truth_data(res,timestep)
       !subroutine to write the true data to a file

       use resdomain, only : tile_4d_and_logp_state_vec_res1d,tile_4d_and_logp_state_vec_input_to_local_grids,get_trainingdataindices, &
                             processor_decomposition_manual, get_z_res_extent, unstandardize_state_vec_input

       !input variables
       integer, intent(in)            :: timestep

       type(main_type), intent(in)    :: res

       !local variables
       real(kind=dp), allocatable :: wholegrid4d(:,:,:,:), wholegrid2d(:,:)
       real(kind=dp), allocatable :: sendreceivedata(:), temp4d(:,:,:,:), temp2d(:,:)
       real(kind=dp), allocatable :: copy_truth_vec(:)

       !mpi and grid local variables
       integer, parameter :: root=0
       integer :: i, recieverequest, sendrequest
       integer(kind=int32) :: from, to, tag, tag2
       integer :: local_domain_size, receive_size
       integer :: localres_xstart,localres_xend,localres_ystart,localres_yend,localresxchunk,localresychunk
       integer :: localinput_xstart,localinput_xend,localinput_ystart,localinput_yend,localinputxchunk,localinputychunk
       integer :: local_res_zstart,local_res_zend,local_reszchunk
       integer :: status(MPI_STATUS_SIZE)

       integer :: counter, proc_num, j, number_of_regions

       integer, allocatable :: region_indices(:)

       logical, parameter :: setflag=.False.
       logical            :: localpole,localperiodicboundary

       !Filename stuff
       character(len=9)  :: truth_out_root
       character(len=3)  :: file_end
       character(len=6)  :: trial_word
       character(len=2)  :: month
       character(len=4)  :: year
       character(len=2)  :: day
       character(len=2)  :: hour
       character(len=:), allocatable :: date_file
       character(len=:), allocatable :: truth_out_file_name
       character(len=:), allocatable :: file_path


        !The receiving part of the routine
        !Gets all of the outvecs from each worker
        !and gives it to the master node (worker == 0)
        !Master node reconstructs the whole global set vector and
        !then writes it out to the disk
        if(mpi_res%is_root) then
           allocate(wholegrid4d(res%model_parameters%full_predictvars,xgrid,ygrid,zgrid))
           allocate(wholegrid2d(xgrid,ygrid))

           wholegrid4d = 0
           wholegrid2d = 0

           do i=1, res%model_parameters%num_of_regions_on_proc
              do j=1, res%model_parameters%num_vert_levels

                 allocate(temp4d(res%reservoir(i,j)%local_predictvars,res%grid(i,j)%resxchunk,res%grid(i,j)%resychunk,res%reservoir(i,j)%local_heightlevels_res))
                 allocate(temp2d(res%grid(i,j)%resxchunk,res%grid(i,j)%resychunk))

                 allocate(copy_truth_vec(res%reservoir(i,j)%reservoir_numinputs-res%reservoir(i,j)%tisr_size_input))

                 !print *, 'shape(res%reservoir(i,j)%predictiondata)',shape(res%reservoir(i,j)%predictiondata)
                 !print *, 'res%model_parameters%synclength/res%model_parameters%timestep+timestep',res%model_parameters%synclength/res%model_parameters%timestep+timestep
                 copy_truth_vec = res%reservoir(i,j)%predictiondata(1:res%reservoir(i,j)%reservoir_numinputs-res%reservoir(i,j)%sst_size_input-res%reservoir(i,j)%tisr_size_input,res%model_parameters%synclength/res%model_parameters%timestep+timestep)
                 call unstandardize_state_vec_input(res%reservoir(i,j),res%grid(i,j),copy_truth_vec)

                 call tile_4d_and_logp_state_vec_input_to_local_grids(res%model_parameters,copy_truth_vec,res%model_parameters%region_indices(i),j,temp4d,temp2d)

                 wholegrid4d(:,res%grid(i,j)%res_xstart:res%grid(i,j)%res_xend,res%grid(i,j)%res_ystart:res%grid(i,j)%res_yend,res%grid(i,j)%res_zstart:res%grid(i,j)%res_zend) = temp4d
                 if(res%reservoir(i,j)%logp_bool) then
                   wholegrid2d(res%grid(i,j)%res_xstart:res%grid(i,j)%res_xend,res%grid(i,j)%res_ystart:res%grid(i,j)%res_yend) = temp2d
                 endif 
                 deallocate(temp4d)
                 deallocate(temp2d)
                 deallocate(copy_truth_vec)
              enddo
           enddo
        endif

        tag = 11
        tag2 = 12

        call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)
        !TODO instead of doing the do loop then the if statement try
        !if(mpi_res%is_root) loop over all workers
        !else send data

        counter = 1
        if(.not.(mpi_res%is_root)) then
           do i=1, res%model_parameters%num_of_regions_on_proc
              do j=1, res%model_parameters%num_vert_levels
                 local_domain_size = size(res%reservoir(i,j)%predictiondata(1:res%reservoir(i,j)%reservoir_numinputs-res%reservoir(i,j)%tisr_size_input-res%reservoir(i,j)%sst_size_input,res%model_parameters%synclength/res%model_parameters%timestep+timestep))

                 allocate(copy_truth_vec(res%reservoir(i,j)%reservoir_numinputs-res%reservoir(i,j)%tisr_size_input-res%reservoir(i,j)%sst_size_input))

                 copy_truth_vec = res%reservoir(i,j)%predictiondata(1:res%reservoir(i,j)%reservoir_numinputs-res%reservoir(i,j)%tisr_size_input-res%reservoir(i,j)%sst_size_input,res%model_parameters%synclength/res%model_parameters%timestep+timestep)
                 call unstandardize_state_vec_input(res%reservoir(i,j),res%grid(i,j),copy_truth_vec) 

                 allocate(sendreceivedata(local_domain_size))
                  
                 sendreceivedata = copy_truth_vec

                 deallocate(copy_truth_vec)

                 !if(res%model_parameters%irank == 1) then
                 !print *,'proc_num,i,j,local_domain_size,counter',res%model_parameters%irank,i,j,local_domain_size,counter
                 !print *,'sendreceivedata(1:10)',sendreceivedata(1:10)
                 !endif  
                 to = root

                 tag = counter

                 call MPI_SEND(sendreceivedata,local_domain_size,MPI_DOUBLE_PRECISION,to,tag,mpi_res%mpi_world,mpi_res%ierr)

                 deallocate(sendreceivedata)

                 counter = counter + 1
              enddo
           enddo
        endif

        call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)

        if((mpi_res%is_root)) then
          do proc_num=1,mpi_res%numprocs-1
             call processor_decomposition_manual(proc_num,mpi_res%numprocs,res%model_parameters%number_of_regions,region_indices)

             number_of_regions = size(region_indices)

             counter = 1
             do i=1, number_of_regions
                do j=1, res%model_parameters%num_vert_levels
                    call getxyresextent(res%model_parameters%number_of_regions,region_indices(i),localres_xstart,localres_xend,localres_ystart,localres_yend,localresxchunk,localresychunk)

                    call get_z_res_extent(res%model_parameters%num_vert_levels,j,local_res_zstart,local_res_zend,local_reszchunk)

                    call getsend_receive_size_input(res%model_parameters,region_indices(i),j,receive_size)

                    allocate(sendreceivedata(receive_size))
                    allocate(temp4d(res%model_parameters%full_predictvars,localresxchunk,localresychunk,local_reszchunk))
                    allocate(temp2d(localresxchunk,localresychunk))

                    from = proc_num

                    tag = counter

                    !print *, 'root receiving from',region_indices(i),'j',j,'receive_size',receive_size

                    call MPI_RECV(sendreceivedata,receive_size,MPI_DOUBLE_PRECISION,from,tag,mpi_res%mpi_world,MPI_STATUS_IGNORE,mpi_res%ierr)

                    call tile_4d_and_logp_state_vec_input_to_local_grids(res%model_parameters,sendreceivedata,region_indices(i),j,temp4d,temp2d)

                    wholegrid4d(:,localres_xstart:localres_xend,localres_ystart:localres_yend,local_res_zstart:local_res_zend) = temp4d

                    if(local_res_zend == zgrid) then
                      wholegrid2d(localres_xstart:localres_xend,localres_ystart:localres_yend) = temp2d
                    endif 

                    deallocate(sendreceivedata)
                    deallocate(temp4d)
                    deallocate(temp2d)

                    counter = counter + 1
                 enddo
              enddo
           enddo
        endif

       if(mpi_res%is_root) then
           print *, 'root writing truth'
           truth_out_root = 'era_truth'
           file_end = '.nc'
           trial_word = 'trial_'
           file_path = '/scratch/user/troyarcomano/Predictions/Hybrid/'

           call get_current_time_delta_hour(calendar,res%model_parameters%traininglength+res%model_parameters%synclength+res%model_parameters%prediction_markers(res%model_parameters%current_trial_number))
           write(year,'(I4.4)') calendar%currentyear
           write(month,'(I2.2)') calendar%currentmonth
           write(day,'(I2.2)') calendar%currentday
           write(hour,'(I2.2)') calendar%currenthour

           date_file = month//'_'//day//'_'//year//'_'//hour
           truth_out_file_name = file_path//truth_out_root//res%model_parameters%trial_name//trial_word//date_file//file_end

           print *, 'writing truth to',truth_out_file_name
           call write_netcdf(res%model_parameters,wholegrid4d,wholegrid2d,timestep,truth_out_file_name)
       endif

       return
     end subroutine

     subroutine run_model(model_parameters,timestep,grid4d,grid2d,sst_grid,speedy_grid4d,speedy_grid2d)
       use mod_utilities, only : unstandardize_data, standardize_data_given_pars3d, e_constant, model_parameters_type
       use mod_io, only : write_netcdf_speedy_full_mpi, read_full_file_4d
       use speedy_main

       integer, intent(in)                        :: timestep

       type(model_parameters_type), intent(inout) :: model_parameters

       real(kind=dp), intent(in)        :: grid4d(:,:,:,:)
       real(kind=dp), intent(in)        :: grid2d(:,:)
       real(kind=dp), intent(in), optional        :: sst_grid(:,:)
       
       real(kind=dp), allocatable, intent(out)    :: speedy_grid4d(:,:,:,:), speedy_grid2d(:,:)

       !local stuff
       real(kind=dp), allocatable :: copy(:,:,:,:)

       integer :: num_of_vars, size_x, size_y, size_z, num_of_standardized_vars
       integer :: i

       character(len=:), allocatable :: file_path 
       character(len=:), allocatable :: speedy_file

       logical :: make_file

       file_path = '/scratch/user/troyarcomano/Predictions/Hybrid/'
       speedy_file = file_path//'hybrid_speedy_out.nc'

       call get_current_time_delta_hour(calendar,model_parameters%traininglength+model_parameters%prediction_markers(model_parameters%current_trial_number)+model_parameters%synclength+timestep*model_parameters%timestep)

       print *,'before speedy hour forecast specific humidity',minval(grid4d(4,:,:,:)),maxval(grid4d(4,:,:,:))

       allocate(copy,source=grid4d) 

       where(copy(4,:,:,:) < 0.000001)
         copy(4,:,:,:) = 0.000001
       endwhere

       !where(grid4d(4,:,:,:) > 25.0)
       !  grid4d(4,:,:,:) = 25.0
       !endwhere

       if(.not. allocated(internal_state_vector%variables3d)) then
         num_of_vars = size(grid4d,1)
         size_x = size(grid4d,2)
         size_y = size(grid4d,3)
         size_z = size(grid4d,4)

         allocate(internal_state_vector%variables3d(num_of_vars,size_x,size_y,size_z))
        endif 
        if(.not. allocated(internal_state_vector%logp)) then
          size_x = size(grid4d,2)
          size_y = size(grid4d,3)
          allocate(internal_state_vector%logp(size_x,size_y))
        endif
   
        internal_state_vector%is_safe_to_run_speedy = .True.
        model_parameters%run_speedy = .True. 

        internal_state_vector%hybrid_slab = model_parameters%slab_ocean_model_bool !Doesnt do anything yet

        if(model_parameters%slab_ocean_model_bool) then
          if(.not. allocated(internal_state_vector%sst_hybrid)) then 
            allocate(internal_state_vector%sst_hybrid(size_x,size_y))
          endif 
          internal_state_vector%sst_hybrid = sst_grid
          print *, 'mean hybrid sst',sum(sst_grid)/size(sst_grid)
        endif 

        internal_state_vector%variables3d = copy !grid4d
        internal_state_vector%logp = grid2d

        internal_state_vector%era_hour = 1 !era hour of the month 1 = 00UTC of the first day of
        internal_state_vector%era_hour_plus_one  = 2!So I dont have to do calendar stuff in

        internal_state_vector%istart = 2
        internal_state_vector%era_start = 3
     
        internal_state_vector%iyear0 = calendar%currentyear
        internal_state_vector%imont0 = calendar%currentmonth
        internal_state_vector%iday = calendar%currentday
        internal_state_vector%ihour = calendar%currenthour

        !Slowly adding sst_bias over the whole run 
        internal_state_vector%sst_bias = 0.0_dp!(2.0/(model_parameters%predictionlength/3))*(model_parameters%timestep*timestep)
        print *, 'timestep, internal_state_vector%sst_bias',timestep, internal_state_vector%sst_bias

        call agcm_main(1,1,internal_state_vector)  

        call clean_up_speedy()

        print *, 'after speedy specific humidity',minval(internal_state_vector%variables3d(4,:,:,:)),maxval(internal_state_vector%variables3d(4,:,:,:))
        print *, 'speedy_file',speedy_file

        !call write_netcdf(model_parameters,internal_state_vector%variables3d,internal_state_vector%logp,timestep,speedy_file)
        !call write_netcdf_speedy_full_mpi(timestep,model_parameters,speedy_file,mpi_res,internal_state_vector%variables3d,internal_state_vector%logp)

        allocate(speedy_grid4d, source=internal_state_vector%variables3d)
        allocate(speedy_grid2d, source=internal_state_vector%logp)

        where(internal_state_vector%variables3d(4,:,:,:) < 0.000001) 
             internal_state_vector%variables3d(4,:,:,:) = 0.000001_dp
        endwhere

        if(internal_state_vector%is_safe_to_run_speedy .eqv. .False.) then
          model_parameters%run_speedy = .False.
        else
          model_parameters%run_speedy = .True.
        endif 
     end subroutine
  
     subroutine clean_up_speedy()
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
     end subroutine 

  subroutine get_tisr_by_date(reservoir,grid,model_parameters,timestep,var1d)
    use mod_utilities, only : standardize_data_given_pars1d, reservoir_type, grid_type, model_parameters_type

    use mod_calendar, only : numof_hours_into_year,get_current_time_delta_hour

    use mod_io, only : read_3d_file_parallel 
   
    type(reservoir_type), intent(inout)     :: reservoir
    type(grid_type), intent(inout)          :: grid
    type(model_parameters_type), intent(in) :: model_parameters

    integer, intent(in)                 :: timestep

    real(kind=dp), intent(inout)        :: var1d(:)

    !local
    integer          :: date_into_year_index

    call get_current_time_delta_hour(calendar,(model_parameters%traininglength+model_parameters%prediction_markers(model_parameters%current_trial_number)+model_parameters%synclength+timestep*model_parameters%timestep))
    call numof_hours_into_year(calendar%currentyear,calendar%currentmonth,calendar%currentday,calendar%currenthour,date_into_year_index)

    if(reservoir%assigned_region == 0) print *,'current date plus timestep',calendar%currentyear,calendar%currentmonth,calendar%currentday,calendar%currenthour,'date_into_year_index',date_into_year_index

    var1d = reshape(reservoir%full_tisr(:,:,date_into_year_index),(/grid%inputxchunk*grid%inputychunk/))

    !call standardize_data_given_pars1d(var1d,grid%mean(grid%tisr_mean_std_idx),grid%std(grid%tisr_mean_std_idx))
  end subroutine
end module mpires

