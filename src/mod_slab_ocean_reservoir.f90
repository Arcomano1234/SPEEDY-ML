module mod_slab_ocean_reservoir
  USE, INTRINSIC :: IEEE_ARITHMETIC
  use MKL_SPBLAS
  use mod_utilities, only : dp, main_type, reservoir_type, grid_type, model_parameters_type

  implicit none
contains

subroutine initialize_slab_ocean_model(reservoir,grid,model_parameters)
  !Routine to allocate the special slab ocean model reservoir and 
  !variables needed for that

  use resdomain, only : set_reservoir_by_region

  type(reservoir_type), intent(inout)     :: reservoir
  type(grid_type), intent(inout)          :: grid
  type(model_parameters_type), intent(inout) :: model_parameters

  integer :: nodes_per_input

  reservoir%local_predictvars = model_parameters%full_predictvars
  reservoir%local_heightlevels_input = grid%inputzchunk

  reservoir%local_heightlevels_res = grid%reszchunk

  model_parameters%ml_only_ocean = .True.

  reservoir%m = 4000!6000

  reservoir%deg = 6
  reservoir%radius = 0.9
  reservoir%beta_res = 0.0001_dp
  reservoir%beta_model = 1.0_dp
  reservoir%sigma = 0.6_dp !0.5_dp

  reservoir%prior_val = 0.0_dp

  reservoir%density = reservoir%deg/reservoir%m

  reservoir%noisemag = 0.10

  reservoir%leakage = 1.0_dp!/14.0_dp !/4.0_dp !12.0_dp

  !call set_reservoir_by_region(reservoir,grid)

  reservoir%sst_bool = .True.
  reservoir%sst_bool_prediction = .True.
  reservoir%sst_bool_input = .True.
  
  reservoir%tisr_input_bool = .True.
  reservoir%atmo_to_ocean_coupled = .True. 
  
  reservoir%ohtc_input = .False.
  reservoir%ohtc_prediction = .False.

  reservoir%num_atmo_levels = 1!grid%inputzchunk

  if(reservoir%sst_bool_input) then
    reservoir%sst_size_res = grid%resxchunk*grid%resychunk
  else
    reservoir%sst_size_res = 0
  endif

  if(reservoir%sst_bool_input) then
    reservoir%sst_size_input = grid%inputxchunk*grid%inputychunk
  else
    reservoir%sst_size_input = 0
  endif

  if(reservoir%sst_climo_input) then
     reservoir%sst_climo_res = grid%inputxchunk*grid%inputychunk
  else
     reservoir%sst_climo_res = 0
  endif

  if(reservoir%tisr_input_bool) then
    reservoir%tisr_size_res = grid%resxchunk*grid%resychunk
  else
    reservoir%tisr_size_res = 0
  endif

  if(reservoir%tisr_input_bool) then
    reservoir%tisr_size_input = grid%inputxchunk*grid%inputychunk
  else
    reservoir%tisr_size_input = 0
  endif

  if(reservoir%atmo_to_ocean_coupled) then
    reservoir%atmo_size_input = grid%inputxchunk*grid%inputychunk*reservoir%local_predictvars + grid%inputxchunk*grid%inputychunk!*reservoir%num_atmo_levels + grid%inputxchunk*grid%inputychunk 
    if(reservoir%precip_input_bool) then
      print *, 'reservoir%precip_input_bool slab',reservoir%precip_input_bool
      reservoir%atmo_size_input = reservoir%atmo_size_input + grid%inputxchunk*grid%inputychunk
    endif 
    !reservoir%atmo_size_input = grid%resxchunk*grid%resychunk*reservoir%local_predictvars*reservoir%num_atmo_levels + grid%resxchunk*grid%resychunk
  else
    reservoir%atmo_size_input = 0
  endif

  if(reservoir%ohtc_input) then
    reservoir%ohtc_input_size = grid%inputxchunk*grid%inputychunk
  else
    reservoir%ohtc_input_size = 0
  endif

  if(reservoir%ohtc_prediction) then
    reservoir%ohtc_res_size = grid%resxchunk*grid%resychunk
  else
    reservoir%ohtc_res_size = 0
  endif 


  if(model_parameters%ml_only_ocean) then
    reservoir%chunk_size_speedy = 0
  endif

  reservoir%chunk_size = reservoir%sst_size_res + reservoir%ohtc_res_size

  reservoir%chunk_size_prediction = reservoir%sst_size_res + reservoir%ohtc_res_size

  reservoir%locality = 0

  reservoir%locality = reservoir%atmo_size_input + reservoir%sst_size_input + reservoir%sst_climo_res + reservoir%tisr_size_input - reservoir%chunk_size

  nodes_per_input = NINT(dble(reservoir%m)/(dble(reservoir%chunk_size)+dble(reservoir%locality)))
  reservoir%n = nodes_per_input*(reservoir%chunk_size+reservoir%locality)
  reservoir%k = reservoir%density*reservoir%n*reservoir%n
  reservoir%reservoir_numinputs = reservoir%chunk_size+reservoir%locality

  if(.not. allocated(reservoir%vals)) allocate(reservoir%vals(reservoir%k))
  if(.not. allocated(reservoir%win))  allocate(reservoir%win(reservoir%n,reservoir%reservoir_numinputs))
  if(.not. allocated(reservoir%wout)) allocate(reservoir%wout(reservoir%chunk_size_prediction,reservoir%n+reservoir%chunk_size_speedy))
  if(.not. allocated(reservoir%rows)) allocate(reservoir%rows(reservoir%k))
  if(.not. allocated(reservoir%cols)) allocate(reservoir%cols(reservoir%k))
end subroutine


subroutine gen_res(reservoir)
   use mod_linalg, only : mklsparse, sparse_eigen, makesparse

   type(reservoir_type)       :: reservoir

   real(kind=dp)              :: eigs,average
   real(kind=dp), allocatable :: newvals(:)

  
   print *,'makesparse'
   call makesparse(reservoir) 

   print *,'sparse_eigen'
   call sparse_eigen(reservoir,reservoir%n*10,6,eigs)
 
   allocate(newvals(reservoir%k)) 
   newvals = (reservoir%vals/eigs)*reservoir%radius
   reservoir%vals = newvals
   
   print *,'remake'
   call mklsparse(reservoir)
   
   if(reservoir%assigned_region == 0) then
     print *, 'region num', reservoir%assigned_region 
     print *, 'radius', reservoir%radius
     print *, 'degree', reservoir%deg
     print *, 'max res', maxval(reservoir%vals)
     print *, 'min res', minval(reservoir%vals)
     print *, 'average', sum(reservoir%vals)/size(reservoir%vals)
     print *, 'k',reservoir%k
     print *, 'eig',eigs
   endif 
   
   return 
end subroutine   

subroutine train_slab_ocean_model(reservoir,grid,model_parameters)
   use mod_utilities, only : init_random_seed
   use resdomain, only : tile_full_input_to_target_data_ocean_model

   type(reservoir_type), intent(inout)        :: reservoir
   type(model_parameters_type), intent(inout) :: model_parameters
   type(grid_type), intent(inout)             :: grid
 
   integer :: q,i,num_inputs,j,k
   integer :: un_noisy_sync
   integer :: betas_res, betas_model,priors
   integer :: vert_loop
    
   real(kind=dp), allocatable :: ip(:),rand(:),average
   real(kind=dp), allocatable :: targetdata_1d(:)
   real(kind=dp), allocatable :: targetdata_2d(:,:)

   character(len=:), allocatable :: base_trial_name
   character(len=50) :: beta_res_char,beta_model_char,prior_char
     
   call gen_res(reservoir)

   q = reservoir%n/reservoir%reservoir_numinputs

   if(reservoir%assigned_region == 0) print *, 'q',q,'n',reservoir%n,'num inputs',reservoir%reservoir_numinputs

   allocate(ip(q))
   allocate(rand(q))

   reservoir%win = 0.0_dp

   do i=1,reservoir%reservoir_numinputs

      call random_number(rand)

      ip = (-1d0 + 2*rand) 
      
      reservoir%win((i-1)*q+1:i*q,i) = reservoir%sigma*ip
   enddo
   
   deallocate(rand)
   deallocate(ip) 

   print *,'starting reservoir_layer'   
  
   call initialize_chunk_training(reservoir,model_parameters) 
 
   if(.not. model_parameters%ml_only_ocean) then 
     allocate(reservoir%imperfect_model_states(reservoir%chunk_size_prediction,size(reservoir%trainingdata,2)))
   
     call tile_full_input_to_target_data_ocean_model(reservoir,grid,reservoir%trainingdata(:,1:model_parameters%timestep_slab),targetdata_2d)
     reservoir%imperfect_model_states(:,1:model_parameters%timestep_slab) = targetdata_2d
   
     deallocate(targetdata_2d)

     call tile_full_input_to_target_data_ocean_model(reservoir,grid,reservoir%trainingdata(:,1:size(reservoir%trainingdata,2)-model_parameters%timestep_slab),targetdata_2d)  
     reservoir%imperfect_model_states(:,model_parameters%timestep_slab+1:size(reservoir%trainingdata,2)) = targetdata_2d

     deallocate(targetdata_2d)
   endif 

   do i=1,model_parameters%timestep_slab
      print *, 'loop number',i
      !if(reservoir%assigned_region == 954) print *, 'reservoir%trainingdata(eservoir_numinputs,1:40)',reservoir%trainingdata(reservoir%reservoir_numinputs,1:40)
      !if(reservoir%assigned_region == 690) print *, 'reservoir%imperfect_model_states(:,i) slab',reservoir%imperfect_model_states(:,i)
      if(reservoir%assigned_region == 690) print *, 'reservoir%trainingdata(:,i) slab', reservoir%trainingdata(:,i)
      if(reservoir%assigned_region == 690) print *, 'shape(reservoir%trainingdata) slab',shape(reservoir%trainingdata)
      if(reservoir%assigned_region == 690) print *, 'shape(reservoir%trainingdata(:,i:model_parameters%traininglength:model_parameters%timestep_slab))',shape(reservoir%trainingdata(:,i:model_parameters%traininglength:model_parameters%timestep_slab))

      if(reservoir%assigned_region == 690) print *, 'reservoir%trainingdata(grid_atmo%sst_start,i:model_parameters%traininglength:model_parameters%timestep_slab)',reservoir%trainingdata(grid%sst_start,i:model_parameters%traininglength:model_parameters%timestep_slab)

      if(model_parameters%ml_only_ocean) then 
        call reservoir_layer_chunking_ml(reservoir,model_parameters,grid,reservoir%trainingdata(:,i:model_parameters%traininglength:model_parameters%timestep_slab))
      else
        call reservoir_layer_chunking_hybrid(reservoir,model_parameters,grid,reservoir%trainingdata(:,i:model_parameters%traininglength:model_parameters%timestep_slab),reservoir%imperfect_model_states(:,i:model_parameters%traininglength:model_parameters%timestep_slab))
      endif

   enddo

   
   !TODO NOTE need to change this !deallocate(reservoir%trainingdata)
    
   if(.not. model_parameters%ml_only_ocean) then
     deallocate(reservoir%imperfect_model_states)
   endif 
 
   print *, 'fitting slab',reservoir%assigned_region

   if(model_parameters%ml_only_ocean) then
     call fit_chunk_ml(reservoir,model_parameters,grid)
   else
     call fit_chunk_hybrid(reservoir,model_parameters,grid)
   endif

   print *, 'cleaning up', reservoir%assigned_region
   call clean_batch(reservoir)
 
end subroutine 

subroutine get_training_data_from_atmo(reservoir,model_parameters,grid,reservoir_atmo,grid_atmo)
   use mod_utilities, only : rolling_average_over_a_period_2d
   use mod_calendar
   use speedy_res_interface, only : read_era, read_model_states
   use resdomain, only : standardize_speedy_data

   type(reservoir_type), intent(inout)        :: reservoir, reservoir_atmo
   type(model_parameters_type), intent(inout) :: model_parameters
   type(grid_type), intent(inout)             :: grid, grid_atmo
  
   integer :: sst_res_input_size, i, counter 
   real(kind=dp), allocatable :: ohtc_var(:,:,:)
    
   allocate(grid%mean,source=grid_atmo%mean)
   allocate(grid%std, source=grid_atmo%std)
 
   print *, 'slab grid%mean',grid%mean
   print *, 'slab grid%std', grid%std
   if(reservoir_atmo%sst_bool_input) then
     print *, 'atmo has sst'
     reservoir%sst_bool_input = .True.
     reservoir%sst_bool_prediction = .True.
   else 
     reservoir%sst_bool_input = .False.
     reservoir%sst_bool_prediction = .False.
   endif 
    
   if(reservoir%sst_bool_prediction) then  

     if((reservoir_atmo%sst_climo_bool).and.(reservoir_atmo%sst_bool_prediction)) then
        reservoir%sst_climo_input = .False.!.True.
     else
        reservoir%sst_climo_input = .False.
     endif 

     if(reservoir_atmo%precip_bool) then
         reservoir%precip_input_bool = .False.
     else 
         reservoir%precip_input_bool = .False.
     endif 

     sst_res_input_size = grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars + grid_atmo%inputxchunk*grid_atmo%inputychunk + grid_atmo%inputxchunk*grid_atmo%inputychunk + grid_atmo%inputxchunk*grid_atmo%inputychunk
   
     print *, 'grid_atmo%sst_start,grid_atmo%sst_end',grid_atmo%sst_start,grid_atmo%sst_end

     grid%atmo3d_start = grid_atmo%atmo3d_start
     grid%atmo3d_end = grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars!grid_atmo%atmo3d_end

     grid%logp_start = grid%atmo3d_end + 1!grid_atmo%logp_start
     grid%logp_end = grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars + grid_atmo%inputxchunk*grid_atmo%inputychunk!grid_atmo%logp_end

     grid%sst_start = grid%logp_end + 1!grid_atmo%sst_start
     grid%sst_end = grid%sst_start + grid_atmo%inputxchunk*grid_atmo%inputychunk - 1!grid_atmo%sst_end

     grid%tisr_start = grid%sst_end + 1
     grid%tisr_end = grid%tisr_start + grid_atmo%inputxchunk*grid_atmo%inputychunk - 1

     grid%sst_mean_std_idx = grid_atmo%sst_mean_std_idx

     print *, 'slab grid%sst_mean_std_idx',grid%sst_mean_std_idx
     print *, 'slab grid%mean(grid%sst_mean_std_idx)',grid%mean(grid%sst_mean_std_idx)
     print *, 'slab grid%std(grid%sst_mean_std_idx)', grid%std(grid%sst_mean_std_idx)

     allocate(reservoir%trainingdata(sst_res_input_size,size(reservoir_atmo%trainingdata,2)))

     allocate(reservoir%atmo_training_data_idx(sst_res_input_size))

     counter = 0
     do i=grid_atmo%atmo3d_end-grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars+1,grid_atmo%logp_end
        counter = counter + 1
        reservoir%atmo_training_data_idx(counter) = i
     enddo 
     do i=grid_atmo%sst_start,grid_atmo%sst_end
        counter = counter + 1
        reservoir%atmo_training_data_idx(counter) = i
     enddo
     do i=grid_atmo%tisr_start,grid_atmo%tisr_end
        counter = counter + 1
        reservoir%atmo_training_data_idx(counter) = i
     enddo

     print *, 'counter',counter 
     print *, 'shape(reservoir%trainingdata)',shape(reservoir%trainingdata)
     print *, 'grid%atmo3d_start,grid%logp_end',grid%atmo3d_start,grid%logp_end
     print *, 'tisr_end',grid%tisr_end

     reservoir%ohtc_input = .False.
     if(reservoir%ohtc_input) then 
       call read_ohtc_parallel_training(reservoir,model_parameters,grid,ohtc_var)
     endif 

     reservoir%trainingdata(grid%atmo3d_start:grid%logp_end,:) = reservoir_atmo%trainingdata(grid_atmo%atmo3d_end - grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars+1:grid_atmo%logp_end,:)
     print *, 'grid_atmo%sst_start:grid_atmo%sst_end,',grid_atmo%sst_start,grid_atmo%sst_end
     print *, '(grid%sst_start:grid%sst_end',grid%sst_start,grid%sst_end
     reservoir%trainingdata(grid%sst_start:grid%sst_end,:) = reservoir_atmo%trainingdata(grid_atmo%sst_start:grid_atmo%sst_end,:)
     reservoir%trainingdata(grid%tisr_start:grid%tisr_end,:) = reservoir_atmo%trainingdata(grid_atmo%tisr_start:grid_atmo%tisr_end,:)
 
     if(reservoir%assigned_region == 10) print *, 'before reservoir%trainingdata(grid%sst_start,1:100)',reservoir%trainingdata(grid%sst_start,1:100)
     call rolling_average_over_a_period_2d(reservoir%trainingdata,model_parameters%timestep_slab) 
     if(reservoir%assigned_region == 10) print *, 'after reservoir%trainingdata(grid%sst_start,1:100)',reservoir%trainingdata(grid%sst_start,1:100)
     !print *, 'better slab training data',reservoir%trainingdata(:,1000)
  endif 

  deallocate(reservoir_atmo%trainingdata)
end subroutine     

subroutine get_prediction_data_from_atmo(reservoir,model_parameters,grid,reservoir_atmo,grid_atmo,delete_atmo_data)
   use mod_utilities, only : rolling_average_over_a_period_2d
   use mod_calendar
   use speedy_res_interface, only : read_era, read_model_states
   use resdomain, only : standardize_speedy_data

   type(reservoir_type), intent(inout)        :: reservoir, reservoir_atmo
   type(model_parameters_type), intent(inout) :: model_parameters
   type(grid_type), intent(inout)             :: grid, grid_atmo

   logical, optional :: delete_atmo_data

   integer :: atmo_ocean_tstep_ratio
   integer :: num_syncs 

   real(kind=dp), allocatable :: temp(:,:)

   atmo_ocean_tstep_ratio = model_parameters%timestep_slab/model_parameters%timestep

   num_syncs = size(reservoir_atmo%predictiondata(:,1:size(reservoir_atmo%predictiondata,2):atmo_ocean_tstep_ratio),2)
   allocate(reservoir%predictiondata(reservoir%reservoir_numinputs,num_syncs))

   print *, 'num_syncs',num_syncs
   print *, 'size(reservoir_atmo%predictiondata,2)',size(reservoir_atmo%predictiondata,2)
   print *, 'atmo_ocean_tstep_ratio',atmo_ocean_tstep_ratio
   print *, 'model_parameters%timestep_slab,model_parameters%timestep',model_parameters%timestep_slab,model_parameters%timestep
   
   print *, 'shape(reservoir%predictiondate)',shape(reservoir%predictiondata)
   print *, 'shape(reservoir_atmo%predictiondata(,1:size(reservoir_atmo%predictiondata,2):atmo_ocean_tstep_ratio))',shape(reservoir_atmo%predictiondata(grid_atmo%atmo3d_end - grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars+1:grid_atmo%logp_end,1:size(reservoir_atmo%predictiondata,2):atmo_ocean_tstep_ratio))

   !Remember reservoir_atmo%predictiondata has a time resolution of
   !model_parameters%timestep so its not hourly (probably)
    allocate(temp,source = reservoir_atmo%predictiondata)
    call rolling_average_over_a_period_2d(temp,atmo_ocean_tstep_ratio)

   reservoir%predictiondata(grid%atmo3d_start:grid%logp_end,:) = temp(grid_atmo%atmo3d_end - grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars+1:grid_atmo%logp_end,1:size(reservoir_atmo%predictiondata,2):atmo_ocean_tstep_ratio)
   reservoir%predictiondata(grid%sst_start:grid%sst_end,:) = temp(grid_atmo%sst_start:grid_atmo%sst_end,1:size(reservoir_atmo%predictiondata,2):atmo_ocean_tstep_ratio)
   reservoir%predictiondata(grid%tisr_start:grid%tisr_end,:) = temp(grid_atmo%tisr_start:grid_atmo%tisr_end,1:size(reservoir_atmo%predictiondata,2):atmo_ocean_tstep_ratio)

   deallocate(temp)

   print *, 'shape(reservoir%predictiondata) slab',shape(reservoir%predictiondata)
   if(.not.(present(delete_atmo_data))) then
     deallocate(reservoir_atmo%predictiondata)
   endif 
end subroutine 

subroutine get_training_data(reservoir,model_parameters,grid,loop_index)
   use mod_utilities, only : era_data_type, speedy_data_type, &
                             standardize_data_given_pars_5d_logp_tisr, &
                             standardize_data_given_pars_5d_logp, &
                             standardize_data_given_pars5d, &
                             standardize_data
   use mod_calendar
   use speedy_res_interface, only : read_era, read_model_states
   use resdomain, only : standardize_speedy_data

   type(reservoir_type), intent(inout)        :: reservoir
   type(model_parameters_type), intent(inout) :: model_parameters
   type(grid_type), intent(inout)             :: grid

   integer, intent(in) :: loop_index

   type(era_data_type)    :: era_data
   type(speedy_data_type) :: speedy_data

   call initialize_calendar(calendar,1990,1,1,0)
   call get_current_time_delta_hour(calendar,model_parameters%discardlength+model_parameters%traininglength+model_parameters%synclength)

   print *, 'reading era states' 

   !Read data in stride and whats only needed for this loop of training
   call read_era(reservoir,grid,model_parameters,1990,calendar%currentyear,era_data)

   !Match units for specific humidity
   era_data%eravariables(4,:,:,:,:) = era_data%eravariables(4,:,:,:,:)*1000.0_dp
   where (era_data%eravariables(4,:,:,:,:) < 0.0)
    era_data%eravariables(4,:,:,:,:) = 0.0_dp
   end where
   !print *, 'shape(era_data%eravariables)',shape(era_data%eravariables)
   !print *, 'era_data',era_data%eravariables(1,:,:,:,10:11)
   !Make sure tisr doesnt have zeroes 
   if(reservoir%tisr_input_bool) then
    where(era_data%era_tisr < 0.0_dp)
      era_data%era_tisr = 0.0_dp
    end where
  endif

  !if(reservoir%assigned_region == 954) print *, 'era_data%eravariables(4,2,2,:,1)', era_data%eravariables(4,2,2,:,1)
  !if(reservoir%assigned_region == 954) print *, ' era_data%era_tisr(:,1:6)',era_data%era_tisr(4,4,1:6)

  if(reservoir%assigned_region == 690) then
    print *, 'era max min temp before',maxval(era_data%eravariables(1,:,:,:,:)),minval(era_data%eravariables(1,:,:,:,:))
    print *, 'era max min u-wind before',maxval(era_data%eravariables(2,:,:,:,:)),minval(era_data%eravariables(2,:,:,:,:))
    print *, 'era max min v-wind before',maxval(era_data%eravariables(3,:,:,:,:)),minval(era_data%eravariables(3,:,:,:,:))
    print *, 'era max min sp before',maxval(era_data%eravariables(4,:,:,:,:)),minval(era_data%eravariables(4,:,:,:,:))
    if(reservoir%logp_bool) print *, 'era max min logp before',maxval(era_data%era_logp),minval(era_data%era_logp)

    if(reservoir%tisr_input_bool) print *, 'era max min tisr before',maxval(era_data%era_tisr),minval(era_data%era_tisr)
  endif
   !Get mean and standard deviation for the first stride of data and use those
   !values for the rest of the program
   if(loop_index == 1) then   
     !Standardize each variable using local std and mean and save the std and
     !mean
     if((reservoir%tisr_input_bool).and.(reservoir%logp_bool)) then
        allocate(grid%mean(reservoir%local_predictvars*reservoir%local_heightlevels_input+2),grid%std(reservoir%local_predictvars*reservoir%local_heightlevels_input+2))

        grid%tisr_mean_std_idx = reservoir%local_predictvars*reservoir%local_heightlevels_input+2 
        grid%logp_mean_std_idx = reservoir%local_predictvars*reservoir%local_heightlevels_input+1

        call standardize_data(reservoir,era_data%eravariables,era_data%era_logp,era_data%era_tisr,grid%mean,grid%std)
     elseif(reservoir%logp_bool) then
        allocate(grid%mean(reservoir%local_predictvars*reservoir%local_heightlevels_input+1),grid%std(reservoir%local_predictvars*reservoir%local_heightlevels_input+1))
        grid%logp_mean_std_idx = reservoir%local_predictvars*reservoir%local_heightlevels_input+1

        call standardize_data(reservoir,era_data%eravariables,era_data%era_logp,grid%mean,grid%std)
     elseif(reservoir%tisr_input_bool) then
        allocate(grid%mean(reservoir%local_predictvars*reservoir%local_heightlevels_input+1),grid%std(reservoir%local_predictvars*reservoir%local_heightlevels_input+1))
 
        grid%tisr_mean_std_idx = reservoir%local_predictvars*reservoir%local_heightlevels_input+1

        call standardize_data(reservoir,era_data%eravariables,era_data%era_tisr,grid%mean,grid%std)
     else
        allocate(grid%mean(reservoir%local_predictvars*reservoir%local_heightlevels_input),grid%std(reservoir%local_predictvars*reservoir%local_heightlevels_input))
        call standardize_data(reservoir,era_data%eravariables,grid%mean,grid%std)
     endif 
   else 
     !Standardize the data from the first stride's std and mean 
     if((reservoir%tisr_input_bool).and.(reservoir%logp_bool)) then
        call standardize_data_given_pars_5d_logp_tisr(grid%mean,grid%std,era_data%eravariables,era_data%era_logp,era_data%era_tisr)
     elseif(reservoir%logp_bool) then
        call standardize_data_given_pars_5d_logp(grid%mean,grid%std,era_data%eravariables,era_data%era_logp)
     elseif(reservoir%tisr_input_bool) then
        call standardize_data_given_pars_5d_logp(grid%mean,grid%std,era_data%eravariables,era_data%era_tisr)
     else
        call standardize_data_given_pars5d(grid%mean,grid%std,era_data%eravariables)
     endif
   endif 

   if(reservoir%assigned_region == 954) then
    print *, 'era max min temp after',maxval(era_data%eravariables(1,:,:,:,:)),minval(era_data%eravariables(1,:,:,:,:))
    print *, 'era max min u-wind after',maxval(era_data%eravariables(2,:,:,:,:)),minval(era_data%eravariables(2,:,:,:,:))
    print *, 'era max min v-wind after',maxval(era_data%eravariables(3,:,:,:,:)),minval(era_data%eravariables(3,:,:,:,:))
    print *, 'era max min sp after',maxval(era_data%eravariables(4,:,:,:,:)),minval(era_data%eravariables(4,:,:,:,:))
    if(reservoir%logp_bool) print *, 'era max min logp after',maxval(era_data%era_logp),minval(era_data%era_logp)

    if(reservoir%tisr_input_bool) print *, 'era max min tisr after',maxval(era_data%era_tisr),minval(era_data%era_tisr)
    print *, 'res%mean,res%std',grid%mean,grid%std
  endif
   !Lets get some training data
   allocate(reservoir%trainingdata(reservoir%reservoir_numinputs,size(era_data%eravariables,5)))

   print *, 'reservoir%reservoir_numinputs',reservoir%assigned_region,grid%level_index
   reservoir%trainingdata(1:reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk,:) = reshape(era_data%eravariables,(/reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk,size(era_data%eravariables,5)/))
   
   if(reservoir%logp_bool) then
     reservoir%trainingdata(reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk+1:reservoir%reservoir_numinputs-reservoir%tisr_size_input,:) = reshape(era_data%era_logp,(/grid%inputxchunk*grid%inputychunk,size(era_data%eravariables,5)/))
   endif 
   
   if(reservoir%tisr_input_bool) then
     reservoir%trainingdata(reservoir%reservoir_numinputs-reservoir%tisr_size_input+1:reservoir%reservoir_numinputs,:) = reshape(era_data%era_tisr,(/grid%inputxchunk*grid%inputychunk,size(era_data%eravariables,5)/))
   endif

   print *, 'reservoir%trainingdata(:,500)',reservoir%trainingdata(:,500)

   
   deallocate(era_data%eravariables)
   deallocate(era_data%era_logp)

   if(allocated(era_data%era_tisr)) then
     deallocate(era_data%era_tisr)
   endif

   !Portion of the routine for getting speedy (imperfect model) data
   print *, 'reading model states'
   call read_model_states(reservoir,grid,model_parameters,1990,calendar%currentyear,speedy_data)
   
   !Lets get imperfect model states
   where(speedy_data%speedyvariables(4,:,:,:,:) < 0.0)
      speedy_data%speedyvariables(4,:,:,:,:) = 0.0_dp
   end where

   if(reservoir%assigned_region == 954) then
     print *, 'speedy max min temp',maxval(speedy_data%speedyvariables(1,:,:,:,:)),minval(speedy_data%speedyvariables(1,:,:,:,:))
     print *, 'speedy max min u-wind',maxval(speedy_data%speedyvariables(2,:,:,:,:)),minval(speedy_data%speedyvariables(2,:,:,:,:))
     print *, 'speedy max min v-wind',maxval(speedy_data%speedyvariables(3,:,:,:,:)),minval(speedy_data%speedyvariables(3,:,:,:,:))
     print *, 'speedy max min sp',maxval(speedy_data%speedyvariables(4,:,:,:,:)),minval(speedy_data%speedyvariables(4,:,:,:,:))
     if(reservoir%logp_bool) print *, 'speedy max min logp',maxval(speedy_data%speedy_logp),minval(speedy_data%speedy_logp)

     print *, 'res%mean,res%std',grid%mean, grid%std
   endif

   if(reservoir%assigned_region == 954) print *, 'speedy_data%speedyvariables(:,1,1,1,1)',speedy_data%speedyvariables(:,1,1,1,1)

   call standardize_speedy_data(reservoir,grid,speedy_data)

   if(reservoir%assigned_region == 954) print *, 'speedy_data%speedyvariables(:,1,1,1,1) after',speedy_data%speedyvariables(:,1,1,1,1)
   allocate(reservoir%imperfect_model_states(reservoir%chunk_size_prediction,size(speedy_data%speedyvariables,5)))
   reservoir%imperfect_model_states = 0.0_dp

   reservoir%imperfect_model_states(1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk,:) = reshape(speedy_data%speedyvariables,(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk,size(speedy_data%speedyvariables,5)/))

   if(reservoir%logp_bool) then
      reservoir%imperfect_model_states(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk+1:reservoir%chunk_size_prediction,:) = reshape(speedy_data%speedy_logp,(/grid%resxchunk*grid%resychunk,size(speedy_data%speedyvariables,5)/))
   endif 

   deallocate(speedy_data%speedyvariables)
   deallocate(speedy_data%speedy_logp)
end subroutine 

subroutine get_prediction_data(reservoir,model_parameters,grid,start_index,length)
   use mod_utilities, only : era_data_type, speedy_data_type, &
                             standardize_data_given_pars_5d_logp_tisr, &
                             standardize_data_given_pars_5d_logp, &
                             standardize_data_given_pars5d, &
                             standardize_data
   use mod_calendar
   use speedy_res_interface, only : read_era, read_model_states 
   use resdomain, only : standardize_speedy_data

   type(reservoir_type), intent(inout)        :: reservoir
   type(model_parameters_type), intent(inout) :: model_parameters
   type(grid_type), intent(inout)             :: grid

   integer, intent(in) :: start_index,length
   
   integer                :: hours_into_first_year, start_year
   integer                :: start_time_memory_index, end_time_memory_index

   type(era_data_type)    :: era_data
   type(speedy_data_type) :: speedy_data  

   call get_current_time_delta_hour(calendar,start_index)

   call numof_hours_into_year(calendar%currentyear,calendar%currentmonth,calendar%currentday,calendar%currenthour,hours_into_first_year)
 
   start_year = calendar%currentyear

   call get_current_time_delta_hour(calendar,start_index+length) 

   !Read data in stride and whats only needed for this loop of training
   call read_era(reservoir,grid,model_parameters,start_year,calendar%currentyear,era_data,1)

   start_time_memory_index = hours_into_first_year   
   end_time_memory_index = start_time_memory_index + length

   !Match units for specific humidity
   era_data%eravariables(4,:,:,:,:) = era_data%eravariables(4,:,:,:,:)*1000.0_dp

   where (era_data%eravariables(4,:,:,:,:) < 0.0)
    era_data%eravariables(4,:,:,:,:) = 0.0_dp
   end where

   !Make sure tisr doesnt have zeroes 
   if(reservoir%tisr_input_bool) then
    where(era_data%era_tisr < 0.0_dp)
      era_data%era_tisr = 0.0_dp
    end where
  endif

  !Standardize the data from mean and std of training data
  if((reservoir%tisr_input_bool).and.(reservoir%logp_bool)) then
     call standardize_data_given_pars_5d_logp_tisr(grid%mean,grid%std,era_data%eravariables,era_data%era_logp,era_data%era_tisr)
  elseif(reservoir%logp_bool) then
     call standardize_data_given_pars_5d_logp(grid%mean,grid%std,era_data%eravariables,era_data%era_logp)
  elseif(reservoir%tisr_input_bool) then
     call standardize_data_given_pars_5d_logp(grid%mean,grid%std,era_data%eravariables,era_data%era_tisr)
  else
     call standardize_data_given_pars5d(grid%mean,grid%std,era_data%eravariables)
  endif

 
   if(allocated(reservoir%predictiondata)) then
     deallocate(reservoir%predictiondata)
   endif 

   !Lets get some prediction data
   allocate(reservoir%predictiondata(reservoir%reservoir_numinputs,length/model_parameters%timestep_slab))

   reservoir%predictiondata(1:reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk,:) = reshape(era_data%eravariables(:,:,:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep_slab),[reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk,length/model_parameters%timestep_slab])
   
   if(reservoir%logp_bool) then
     reservoir%predictiondata(reservoir%local_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk+1:reservoir%reservoir_numinputs-reservoir%tisr_size_input,:) = reshape(era_data%era_logp(:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep_slab),[grid%inputxchunk*grid%inputychunk,length/model_parameters%timestep_slab])
   endif 
   
   if(reservoir%tisr_input_bool) then
     reservoir%predictiondata(reservoir%reservoir_numinputs-reservoir%tisr_size_input+1:reservoir%reservoir_numinputs,:) = reshape(era_data%era_tisr(:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep),[grid%inputxchunk*grid%inputychunk,length/model_parameters%timestep])
   endif

   deallocate(era_data%eravariables)
   deallocate(era_data%era_logp)

   if(allocated(era_data%era_tisr)) then
     deallocate(era_data%era_tisr)
   endif

   !Portion of the routine for getting speedy (imperfect model) data
   print *, 'reading model states'
   call read_model_states(reservoir,grid,model_parameters,start_year,calendar%currentyear,speedy_data,1)
   
   !Lets get imperfect model states
   where(speedy_data%speedyvariables(4,:,:,:,:) < 0.0)
      speedy_data%speedyvariables(4,:,:,:,:) = 0.0_dp
   end where

   call standardize_speedy_data(reservoir,grid,speedy_data) 

   if(allocated(reservoir%imperfect_model_states)) then
     deallocate(reservoir%imperfect_model_states)
   endif 

   allocate(reservoir%imperfect_model_states(reservoir%chunk_size_prediction,length/model_parameters%timestep_slab))

   reservoir%imperfect_model_states = 0.0_dp

   reservoir%imperfect_model_states(1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk,:) = reshape(speedy_data%speedyvariables(:,:,:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep_slab),[reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk,length/model_parameters%timestep_slab])

   if(reservoir%logp_bool) then
      reservoir%imperfect_model_states(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk+1:reservoir%chunk_size_prediction,:) = reshape(speedy_data%speedy_logp(:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep_slab),[grid%resxchunk*grid%resychunk,length/model_parameters%timestep_slab])
   endif 
   deallocate(speedy_data%speedyvariables)
   deallocate(speedy_data%speedy_logp)
end subroutine 

subroutine initialize_prediction_slab(reservoir,model_parameters,grid,atmo_reservoir,atmo_grid)
   use mod_calendar
   
   type(reservoir_type), intent(inout)        :: reservoir
   type(grid_type), intent(inout)             :: grid 
   type(reservoir_type), intent(inout)        :: atmo_reservoir
   type(grid_type), intent(inout)             :: atmo_grid

   type(model_parameters_type), intent(inout) :: model_parameters

   integer :: q,i,num_inputs,j,k
   integer :: un_noisy_sync
   integer :: betas_res, betas_model,priors
   integer :: vert_loop

   real(kind=dp), allocatable :: ip(:),rand(:),average
   real(kind=dp), allocatable :: test_beta_res(:), test_beta_model(:), test_priors(:)
   real(kind=dp), allocatable :: states_x_states_original_copy(:,:)

   character(len=:), allocatable :: base_trial_name
   character(len=50) :: beta_res_char,beta_model_char,prior_char

   !Try syncing on un-noisy data
   if(.not.(allocated(reservoir%saved_state))) allocate(reservoir%saved_state(reservoir%n))
   reservoir%saved_state = 0
   un_noisy_sync = 2160!700

   !From this point on reservoir%trainingdata and reservoir%imperfect_model have a temporal resolution model_parameters%timestep_slab
   !instead of 1 hour resolution and atmo_reservoir%trainingdata has a temporal
   !resolution model_parameters%timestep
   call get_prediction_data_from_atmo(reservoir,model_parameters,grid,atmo_reservoir,atmo_grid)

   call synchronize(reservoir,reservoir%predictiondata,reservoir%saved_state,un_noisy_sync/(model_parameters%timestep_slab)-1)

   if(reservoir%tisr_input_bool) then
     allocate(reservoir%full_tisr,source=atmo_reservoir%full_tisr)
   endif 

   deallocate(reservoir%predictiondata)

   allocate(reservoir%local_model(reservoir%chunk_size_prediction))
   allocate(reservoir%outvec(reservoir%chunk_size_prediction))
   allocate(reservoir%feedback(reservoir%reservoir_numinputs))
   allocate(reservoir%averaged_atmo_input_vec(reservoir%reservoir_numinputs,model_parameters%timestep_slab/model_parameters%timestep-1))
   reservoir%averaged_atmo_input_vec = 0.0_dp
end subroutine 

subroutine get_full_tisr(reservoir,model_parameters,grid)
   use mpires, only : mpi_res
   use mod_io, only : read_3d_file_parallel

   type(reservoir_type), intent(inout)        :: reservoir
   type(grid_type), intent(inout)             :: grid
   type(model_parameters_type), intent(inout) :: model_parameters

   character(len=:), allocatable :: file_path
   character(len=:), allocatable :: tisr_file
 
   file_path = '/scratch/user/troyarcomano/ERA_5/2012/'
   tisr_file = file_path//'toa_incident_solar_radiation_2012_regridded_classic4.nc'

   call read_3d_file_parallel(tisr_file,'tisr',mpi_res,grid,reservoir%full_tisr,1,1)

end subroutine  

subroutine start_prediction_slab(reservoir,model_parameters,grid,atmo_reservoir,atmo_grid,prediction_number)
   use resdomain, only : tile_full_input_to_target_data_ocean_model

   type(reservoir_type), intent(inout)        :: reservoir, atmo_reservoir
   type(grid_type), intent(inout)             :: grid, atmo_grid
   type(model_parameters_type), intent(inout) :: model_parameters

   integer, intent(in)                        :: prediction_number

   model_parameters%current_trial_number = prediction_number

   call get_prediction_data_from_atmo(reservoir,model_parameters,grid,atmo_reservoir,atmo_grid,.False.)

   call synchronize(reservoir,reservoir%predictiondata(:,1:model_parameters%synclength/model_parameters%timestep_slab),reservoir%saved_state,model_parameters%synclength/model_parameters%timestep_slab)

   print *, 'model_parameters%synclength/model_parameters%timestep_slab',model_parameters%synclength/model_parameters%timestep_slab
   print *, 'shape(reservoir%predictiondata)',shape(reservoir%predictiondata)
   reservoir%feedback = reservoir%predictiondata(:,model_parameters%synclength/model_parameters%timestep_slab)
  
   !This is a trick so that we can store the last sst era5 data for plotting
   !when the hybrid model intergration step is less than the slab_timestep
   !e.g. the first x days of whole model sst are those from the era5 data where
   !x == model_parameters%timestep_slab
   call tile_full_input_to_target_data_ocean_model(reservoir,grid,reservoir%predictiondata(:,model_parameters%synclength/model_parameters%timestep_slab),reservoir%outvec)
 
   reservoir%outvec = reservoir%outvec*grid%std(grid%sst_mean_std_idx) + grid%mean(grid%sst_mean_std_idx)

   if(.not. model_parameters%ml_only_ocean) then
     call tile_full_input_to_target_data_ocean_model(reservoir,grid,reservoir%predictiondata(:,model_parameters%synclength/model_parameters%timestep_slab),reservoir%local_model)
     print *, 'start_prediction_slab local_model',reservoir%local_model
   endif 
end subroutine 

subroutine reservoir_layer_chunking_ml(reservoir,model_parameters,grid,trainingdata)
   use mpires
   use mod_utilities, only : gaussian_noise_1d_function


   type(reservoir_type), intent(inout)      :: reservoir
   type(model_parameters_type) , intent(in) :: model_parameters
   type(grid_type) , intent(in)            :: grid

   real(kind=dp), intent(in) :: trainingdata(:,:)

   integer :: i,info
   integer :: training_length, batch_number

   real(kind=dp), allocatable :: temp(:),x(:),y(:),x_(:)
   real(kind=dp), parameter   :: alpha=1.0,beta=0.0
   real(kind=dp), allocatable :: gaussian_noise

   allocate(temp(reservoir%n),x(reservoir%n),x_(reservoir%n),y(reservoir%n))

   x = 0
   y = 0
   do i=1, model_parameters%discardlength/model_parameters%timestep_slab
      info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,x,beta,y)
      print *, 'shape(reservoir%win)',shape(reservoir%win),'shape(trainingdata(:,i))',shape(trainingdata(:,i)), 'worker',reservoir%assigned_region
      temp = matmul(reservoir%win,gaussian_noise_1d_function(trainingdata(:,i),reservoir%noisemag))
      
      x_ = tanh(y+temp)

      x = (1_dp-reservoir%leakage)*x + reservoir%leakage*x_

      y = 0
   enddo

   !call initialize_chunk_training()
 
   reservoir%states(:,1) = x
   batch_number = 0

   training_length = size(trainingdata,2) - model_parameters%discardlength/model_parameters%timestep_slab

   do i=1, training_length-1
      if(reservoir%assigned_region == 1) print *, 'i',i,'batch size',reservoir%batch_size, 'training_length-1',training_length-1,'ize(trainingdata,2)',size(trainingdata,2)
      if(mod(i+1,reservoir%batch_size).eq.0) then
        print *,'chunking slab',i,'region',reservoir%assigned_region
        batch_number = batch_number + 1

        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,reservoir%states(:,mod(i,reservoir%batch_size)),beta,y)
        temp = matmul(reservoir%win,gaussian_noise_1d_function(trainingdata(:,model_parameters%discardlength/model_parameters%timestep_slab+i),reservoir%noisemag))

        x_ = tanh(y+temp)
        x = (1-reservoir%leakage)*x + reservoir%leakage*x_

        reservoir%states(:,reservoir%batch_size) = x

        reservoir%saved_state = reservoir%states(:,reservoir%batch_size)

        reservoir%states(2:reservoir%n:2,:) = reservoir%states(2:reservoir%n:2,:)**2

        call chunking_matmul_ml(reservoir,model_parameters,grid,batch_number,trainingdata)

      elseif (mod(i,reservoir%batch_size).eq.0) then
        print *,'new state',i, 'region',reservoir%assigned_region

        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,reservoir%states(:,reservoir%batch_size),beta,y)
        temp = matmul(reservoir%win,gaussian_noise_1d_function(trainingdata(:,model_parameters%discardlength/model_parameters%timestep_slab+i),reservoir%noisemag)) 
         
        x_ = tanh(y+temp)
        x = (1-reservoir%leakage)*x + reservoir%leakage*x_

        reservoir%states(:,1) = x
      else 
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,reservoir%states(:,mod(i,reservoir%batch_size)),beta,y)
        temp = matmul(reservoir%win,gaussian_noise_1d_function(trainingdata(:,model_parameters%discardlength/model_parameters%timestep_slab+i),reservoir%noisemag))

        x_ = tanh(y+temp)
        x = (1-reservoir%leakage)*x + reservoir%leakage*x_

        reservoir%states(:,mod(i+1,reservoir%batch_size)) = x
      endif 

      y = 0
   enddo
   
   return
end subroutine

subroutine reservoir_layer_chunking_hybrid(reservoir,model_parameters,grid,trainingdata,imperfect_model)
   use mpires
   use mod_utilities, only : gaussian_noise_1d_function,gaussian_noise_1d_function_precip


   type(reservoir_type), intent(inout)      :: reservoir
   type(model_parameters_type) , intent(in) :: model_parameters
   type(grid_type) , intent(in)            :: grid

   real(kind=dp), intent(in) :: trainingdata(:,:)
   real(kind=dp), intent(in) :: imperfect_model(:,:)

   integer :: i,info
   integer :: training_length, batch_number

   real(kind=dp), allocatable :: temp(:), x(:), x_(:), y(:)
   real(kind=dp), parameter   :: alpha=1.0,beta=0.0
   real(kind=dp), allocatable :: gaussian_noise

   allocate(temp(reservoir%n),x(reservoir%n),x_(reservoir%n),y(reservoir%n))

   x = 0
   y = 0
   do i=1, model_parameters%discardlength/model_parameters%timestep
      info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,x,beta,y)

      if(model_parameters%precip_bool) then
        temp = matmul(reservoir%win,gaussian_noise_1d_function_precip(trainingdata(:,i),reservoir%noisemag,grid,model_parameters))
      else
        temp = matmul(reservoir%win,gaussian_noise_1d_function(trainingdata(:,i),reservoir%noisemag))
      endif

      x_ = tanh(y+temp)
      x = (1_dp-reservoir%leakage)*x + reservoir%leakage*x_

      y = 0
   enddo

   !call initialize_chunk_training()

   reservoir%states(:,1) = x
   batch_number = 0

   training_length = size(trainingdata,2) - model_parameters%discardlength/model_parameters%timestep

   do i=1, training_length-1
      if(mod(i+1,reservoir%batch_size).eq.0) then
        print *,'chunking',i, 'region',reservoir%assigned_region
        batch_number = batch_number + 1

        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,reservoir%states(:,mod(i,reservoir%batch_size)),beta,y)

        if(model_parameters%precip_bool) then
          temp = matmul(reservoir%win,gaussian_noise_1d_function_precip(trainingdata(:,model_parameters%discardlength/model_parameters%timestep+i),reservoir%noisemag,grid,model_parameters))
        else
          temp = matmul(reservoir%win,gaussian_noise_1d_function(trainingdata(:,model_parameters%discardlength/model_parameters%timestep+i),reservoir%noisemag))
        endif

        x_ = tanh(y+temp)
        x = (1-reservoir%leakage)*x + reservoir%leakage*x_
 
        reservoir%states(:,reservoir%batch_size) = x

        reservoir%saved_state = reservoir%states(:,reservoir%batch_size)

        reservoir%states(2:reservoir%n:2,:) = reservoir%states(2:reservoir%n:2,:)**2

        call chunking_matmul_hybrid(reservoir,model_parameters,grid,batch_number,trainingdata,imperfect_model)

      elseif (mod(i,reservoir%batch_size).eq.0) then
        print *,'new state',i, 'region',reservoir%assigned_region

        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,reservoir%saved_state,beta,y)

        if(model_parameters%precip_bool) then
          temp = matmul(reservoir%win,gaussian_noise_1d_function_precip(trainingdata(:,model_parameters%discardlength/model_parameters%timestep+i),reservoir%noisemag,grid,model_parameters))
        else
          temp = matmul(reservoir%win,gaussian_noise_1d_function(trainingdata(:,model_parameters%discardlength/model_parameters%timestep+i),reservoir%noisemag))
        endif

        x_ = tanh(y+temp)
        x = (1-reservoir%leakage)*x + reservoir%leakage*x_

        reservoir%states(:,1) = x
      else
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,reservoir%states(:,mod(i,reservoir%batch_size)),beta,y)

        if(model_parameters%precip_bool) then
          temp = matmul(reservoir%win,gaussian_noise_1d_function_precip(trainingdata(:,model_parameters%discardlength/model_parameters%timestep+i),reservoir%noisemag,grid,model_parameters))
        else
          temp = matmul(reservoir%win,gaussian_noise_1d_function(trainingdata(:,model_parameters%discardlength/model_parameters%timestep+i),reservoir%noisemag))
        endif

        x_ = tanh(y+temp)
        x = (1-reservoir%leakage)*x + reservoir%leakage*x_

        reservoir%states(:,mod(i+1,reservoir%batch_size)) = x
      endif

      y = 0
   enddo

   return
end subroutine

subroutine fit_chunk_ml(reservoir,model_parameters,grid)
    !This solves for Wout using least squared solver for the ml only version
    !This should be called only if you are chunking the training
    !There is an option to train using a Prior
    !The prior would try to force the weights of Wout
    !for the numerical model to be near reservoir%prior_val

    use mod_linalg, only : pinv_svd, mldivide
    use mpires
    use mod_io, only : write_netcdf_2d_non_met_data

    type(reservoir_type), intent(inout)        :: reservoir
    type(model_parameters_type), intent(inout) :: model_parameters
    type(grid_type), intent(inout)             :: grid

    integer :: i

    real(kind=dp), allocatable  :: a_trans(:,:), b_trans(:,:), invstates(:,:)
    real(kind=dp), allocatable  :: prior(:,:), temp_beta(:,:)

    real(kind=dp), parameter    :: alpha=1.0, beta=0.0

    character(len=2) :: level_char

    !Do regularization

    do i=1, reservoir%n
         reservoir%states_x_states_aug(i,i) = reservoir%states_x_states_aug(i,i) + reservoir%beta_res
    enddo

    !NOTE moving to mldivide not using pinv anymore
    print *, 'trying mldivide'
    allocate(a_trans(size(reservoir%states_x_states_aug,2),size(reservoir%states_x_states_aug,1)))
    allocate(b_trans(size(reservoir%states_x_trainingdata_aug,2),size(reservoir%states_x_trainingdata_aug,1)))
    a_trans = transpose(reservoir%states_x_states_aug)
    b_trans = transpose(reservoir%states_x_trainingdata_aug)

    if(reservoir%assigned_region == 954)  print *, 'a_trans(1:20,1:20)',a_trans(1:4,1:4)
    if(reservoir%assigned_region == 954)  print *, 'b_trans(1:20,1:20)',b_trans(1:4,1:4)

    call mldivide(a_trans,b_trans)
    reservoir%wout = transpose(b_trans)

    if(reservoir%assigned_region == 954)  print *, 'worker',reservoir%assigned_region,'wout(1,1:20)',reservoir%wout(1,1:20)

    deallocate(a_trans)
    deallocate(b_trans)

    write(level_char,'(i0.2)') grid%level_index
    !if(reservoir%assigned_region == 954) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_954_ocean_'//level_char//'wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')
    !if(reservoir%assigned_region == 217) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_217_ocean_'//level_char//'wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')
    !if(reservoir%assigned_region == 218) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_218_ocean_'//level_char//'wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')

    call write_trained_res(reservoir,model_parameters,grid)

    print *, 'finish fit'
end subroutine

subroutine fit_chunk_hybrid(reservoir,model_parameters,grid)
    !This solves for Wout using least squared solver 
    !This should be called only if you are chunking the training 
    !There is an option to train using a Prior
    !The prior would try to force the weights of Wout 
    !for the numerical model to be near reservoir%prior_val

    use mod_linalg, only : pinv_svd, mldivide
    use mpires
    use mod_io, only : write_netcdf_2d_non_met_data

    type(reservoir_type), intent(inout)        :: reservoir
    type(model_parameters_type), intent(inout) :: model_parameters
    type(grid_type), intent(inout)             :: grid

    integer :: i

    real(kind=dp), allocatable  :: a_trans(:,:), b_trans(:,:), invstates(:,:)
    real(kind=dp), allocatable  :: prior(:,:), temp_beta(:,:)

    real(kind=dp), parameter    :: alpha=1.0, beta=0.0
  
    character(len=2) :: level_char

   

    !If we have a prior we need to make the prior matrix
    if(model_parameters%using_prior) then
      allocate(prior(size(reservoir%states_x_trainingdata_aug,1),size(reservoir%states_x_trainingdata_aug,2)))

      prior = 0.0_dp
       
      do i=1, reservoir%chunk_size_prediction
            prior(i,i) =  reservoir%prior_val*reservoir%beta_model**2.0_dp
      enddo 

    endif 

    !Do regularization

    !If we are doing a prior we need beta^2 
    if(model_parameters%using_prior) then 
      do i=1, reservoir%n+reservoir%chunk_size_prediction
         if(i <= reservoir%chunk_size_prediction) then
            reservoir%states_x_states_aug(i,i) = reservoir%states_x_states_aug(i,i) + reservoir%beta_model**2.0_dp
         else
            reservoir%states_x_states_aug(i,i) = reservoir%states_x_states_aug(i,i) + reservoir%beta_res**2.0_dp
         endif
      enddo
    else
      do i=1, reservoir%n+reservoir%chunk_size_prediction
         if(i <= reservoir%chunk_size_prediction) then 
            reservoir%states_x_states_aug(i,i) = reservoir%states_x_states_aug(i,i) + reservoir%beta_model
         else 
            reservoir%states_x_states_aug(i,i) = reservoir%states_x_states_aug(i,i) + reservoir%beta_res
         endif 
      enddo
    endif 

    !NOTE moving to mldivide not using pinv anymore
    print *, 'trying mldivide'
    allocate(a_trans(size(reservoir%states_x_states_aug,2),size(reservoir%states_x_states_aug,1)))
    allocate(b_trans(size(reservoir%states_x_trainingdata_aug,2),size(reservoir%states_x_trainingdata_aug,1)))
    a_trans = transpose(reservoir%states_x_states_aug)
    b_trans = transpose(reservoir%states_x_trainingdata_aug)

    if(reservoir%assigned_region == 690)  print *, 'slab a_trans(1:4,1:4)',a_trans(1:4,1:4)
    if(reservoir%assigned_region == 690)  print *, 'slab b_trans(1:4,1:4)',b_trans(1:4,1:4)
    if(any(IEEE_IS_NAN(reservoir%states_x_states_aug))) print *, 'reservoir%states_x_states_aug nan', reservoir%assigned_region
    if(any(IEEE_IS_NAN(reservoir%states_x_trainingdata_aug))) print *, 'reservoir%states_x_states_aug nan', reservoir%assigned_region
    if(any(IEEE_IS_NAN(a_trans))) print *, 'a_trans has nan',reservoir%assigned_region
    if(any(IEEE_IS_NAN(b_trans))) print *, 'b_trans has nan',reservoir%assigned_region
    
    !If we are trying a prior then we need to add it to b_trans
    if(model_parameters%using_prior) then
      b_trans = b_trans + transpose(prior)
    endif  

    call mldivide(a_trans,b_trans)
    reservoir%wout = transpose(b_trans)

    if(any(IEEE_IS_NAN(reservoir%wout))) print *, 'wout has nan', reservoir%assigned_region
    if(IEEE_IS_NAN(reservoir%wout(1,1))) print *, 'wout element 1 has nan', reservoir%assigned_region


    if(reservoir%assigned_region == 690)  print *, 'worker',reservoir%assigned_region,'slab wout(1,1:4)',reservoir%wout(1,1:4)

    deallocate(a_trans)
    deallocate(b_trans)

    write(level_char,'(i0.2)') grid%level_index
    if(reservoir%assigned_region == 690) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_690_slab_ocean_wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')
    if(reservoir%assigned_region == 217) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_217_slab_ocean_wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')
    if(reservoir%assigned_region == 218) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_218_slab_ocean_wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')

    !call write_trained_res(reservoir,model_parameters,grid)

    print *, 'finish fit'
end subroutine 

subroutine predictcontroller(reservoir,model_parameters,grid,imperfect_model_in)
    type(reservoir_type), intent(inout)     :: reservoir
    type(model_parameters_type), intent(in) :: model_parameters
    type(grid_type), intent(in)            :: grid

    real(kind=dp), intent(inout) :: imperfect_model_in(:)

    real(kind=dp), allocatable :: x(:)

    allocate(x(reservoir%n))

    x = reservoir%saved_state
 
    call synchronize(reservoir,reservoir%predictiondata(:,1:model_parameters%synclength/model_parameters%timestep_slab),x,model_parameters%synclength/model_parameters%timestep_slab)

    call predict_slab(reservoir,model_parameters,grid,x,imperfect_model_in)    
end subroutine 

subroutine synchronize(reservoir,input,x,length)
    type(reservoir_type), intent(inout) :: reservoir
    
    real(kind=dp), intent(in)     :: input(:,:)
    real(kind=dp), intent(inout)  :: x(:)

    integer, intent(in)           :: length

    real(kind=dp), allocatable    :: y(:), temp(:), x_(:)
    real(kind=dp), parameter      :: alpha=1.0,beta=0.0

    integer :: info,i

    allocate(y(reservoir%n))
    allocate(temp(reservoir%n))
    allocate(x_(reservoir%n))

    y=0
    do i=1, length
       info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,x,beta,y)

       temp = matmul(reservoir%win,input(:,i))

       x_ = tanh(y+temp)
       x = (1-reservoir%leakage)*x + reservoir%leakage*x_

    enddo 
    
    return 
end subroutine  

subroutine predict_slab(reservoir,model_parameters,grid,x,local_model_in)
    use mpires, only : predictionmpicontroller
    use resdomain, only : unstandardize_state_vec_res

    type(reservoir_type), intent(inout)     :: reservoir
    type(model_parameters_type), intent(in) :: model_parameters
    type(grid_type), intent(in)             :: grid

    real(kind=dp), intent(inout) :: x(:)
    real(kind=dp), intent(inout) :: local_model_in(:)

    real(kind=dp), allocatable :: y(:), temp(:)
    real(kind=dp), allocatable :: local_model_temp(:)
    real(kind=dp), allocatable :: x_temp(:),x_augment(:)

    real(kind=dp), parameter :: alpha=1.0,beta=0.0

    integer :: info,i,j

    allocate(y(reservoir%n),temp(reservoir%n))
    allocate(x_augment(reservoir%n+reservoir%chunk_size_prediction))

    y = 0

    info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,x,beta,y)
    temp = matmul(reservoir%win,reservoir%feedback)

    x = tanh(y + temp)

    x_temp = x
    x_temp(2:reservoir%n:2) = x_temp(2:reservoir%n:2)**2

    x_augment(1:reservoir%chunk_size_prediction) = reservoir%local_model
    x_augment(reservoir%chunk_size_prediction+1:reservoir%chunk_size_prediction+reservoir%n) = x_temp

    reservoir%outvec = matmul(reservoir%wout,x_augment)
    reservoir%local_model = reservoir%outvec
    !call unstandardize_state_vec_res(reservoir,grid,reservoir%outvec)
    reservoir%outvec = reservoir%outvec*grid%std(grid%sst_mean_std_idx) + grid%mean(grid%sst_mean_std_idx)

    if((reservoir%assigned_region == 954).and.(mod(i,14*24) == 0)) then
      print *, '*******'
      print *, 'reservoir%predictiondata(grid%sst_start:grid%sst_end,model_parameters%synclength/model_parameters%timestep_slab+10)',reservoir%predictiondata(grid%sst_start:grid%sst_end,model_parameters%synclength/model_parameters%timestep_slab+10)
      print *, 'local_model slab region', reservoir%assigned_region, reservoir%local_model
      print *, 'outvec slab region', reservoir%assigned_region, reservoir%outvec
      print *, '*******'
      print *, 'slab feedback',reservoir%feedback
    endif
end subroutine

subroutine predict_slab_ml(reservoir,model_parameters,grid,x)
    use mpires, only : predictionmpicontroller
    use resdomain, only : unstandardize_state_vec_res
    use mod_utilities, only : e_constant

    type(reservoir_type), intent(inout)     :: reservoir
    type(model_parameters_type), intent(in) :: model_parameters
    type(grid_type), intent(in)             :: grid

    real(kind=dp), intent(inout) :: x(:)

    real(kind=dp), allocatable :: y(:), temp(:), x_(:)
    real(kind=dp), allocatable :: x_temp(:),x_augment(:)

    real(kind=dp), parameter :: alpha=1.0,beta=0.0
  
    integer :: info,i,j

    allocate(y(reservoir%n),temp(reservoir%n),x_(reservoir%n))
    allocate(x_augment(reservoir%n))!reservoir%chunk_size_prediction))

    y = 0

    info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,x,beta,y)
    temp = matmul(reservoir%win,reservoir%feedback)

    x_ = tanh(y + temp)
    x = (1-reservoir%leakage)*x + reservoir%leakage*x_

    x_temp = x
    x_temp(2:reservoir%n:2) = x_temp(2:reservoir%n:2)**2

    x_augment(1:reservoir%n) = x_temp

    reservoir%outvec = matmul(reservoir%wout,x_augment)

    reservoir%outvec = reservoir%outvec*grid%std(grid%sst_mean_std_idx) + grid%mean(grid%sst_mean_std_idx)

    if((reservoir%assigned_region == 954).and.(mod(i,14*24) == 0)) then
      print *, '*******'
      print *, 'reservoir%predictiondata(grid%sst_start:grid%sst_end,model_parameters%synclength/model_parameters%timestep_slab+10)',reservoir%predictiondata(grid%sst_start:grid%sst_end,model_parameters%synclength/model_parameters%timestep_slab+10)
      print *, 'outvec slab region', reservoir%assigned_region, reservoir%outvec
      print *, '*******'
      print *, 'slab feedback',reservoir%feedback
    endif
end subroutine

subroutine clean_sparse(reservoir)
   type(reservoir_type), intent(inout) :: reservoir

   deallocate(reservoir%vals)
   deallocate(reservoir%rows)
   deallocate(reservoir%cols)
end subroutine

subroutine clean_batch(reservoir)
   type(reservoir_type), intent(inout) :: reservoir

   deallocate(reservoir%states_x_trainingdata_aug)
   deallocate(reservoir%states_x_states_aug)
end subroutine

subroutine clean_prediction(reservoir)
   type(reservoir_type), intent(inout) :: reservoir
  
   deallocate(reservoir%local_model)
   deallocate(reservoir%outvec)
   deallocate(reservoir%feedback)

end subroutine 
  
subroutine initialize_chunk_training(reservoir,model_parameters)
   use mod_utilities, only : find_closest_divisor

   type(reservoir_type), intent(inout)     :: reservoir
   type(model_parameters_type), intent(in) :: model_parameters

   integer :: num_of_batches !the number of chunks we want
   integer :: approx_batch_size !approximate size of the batch

   num_of_batches = 1!*2!2 !20*6
   approx_batch_size = (model_parameters%traininglength - model_parameters%discardlength)/(num_of_batches*model_parameters%timestep_slab)

   !routine to get the closest reservoir%batch_size to num_of_batches that
   !divides into reservoir%traininglength
   call find_closest_divisor(approx_batch_size,(model_parameters%traininglength - model_parameters%discardlength)/model_parameters%timestep_slab,reservoir%batch_size)
   print *, 'num_of_batches,reservoir%traininglength,reservoir%batch_size slab',num_of_batches,model_parameters%traininglength,reservoir%batch_size

   !Should be reservoir%n+ reservoir%chunk_size
   allocate(reservoir%states_x_trainingdata_aug(reservoir%chunk_size_prediction,reservoir%n+reservoir%chunk_size_speedy))
   allocate(reservoir%states_x_states_aug(reservoir%n+reservoir%chunk_size_speedy,reservoir%n+reservoir%chunk_size_speedy))
   allocate(reservoir%states(reservoir%n,reservoir%batch_size))
   allocate(reservoir%augmented_states(reservoir%n+reservoir%chunk_size_speedy,reservoir%batch_size))
   allocate(reservoir%saved_state(reservoir%n))

   reservoir%states_x_trainingdata_aug = 0.0_dp
   reservoir%states_x_states_aug = 0.0_dp
   reservoir%states = 0.0_dp
   reservoir%augmented_states = 0.0_dp

end subroutine

subroutine chunking_matmul_ml(reservoir,model_parameters,grid,batch_number,trainingdata)
   use mod_utilities, only : gaussian_noise
   use resdomain, only : tile_full_input_to_target_data_ocean_model
   use mpires

   type(reservoir_type), intent(inout)     :: reservoir
   type(model_parameters_type), intent(in) :: model_parameters
   type(grid_type), intent(in)             :: grid

   integer, intent(in)          :: batch_number

   real(kind=dp), intent(in)    :: trainingdata(:,:)

   real(kind=dp), allocatable   :: temp(:,:), targetdata(:,:)
   real(kind=dp), parameter     :: alpha=1.0, beta=0.0

   integer                      :: n, m, l

   n = size(reservoir%augmented_states,1)
   m = size(reservoir%augmented_states,2)

   reservoir%augmented_states(1:reservoir%n,:) = reservoir%states

   call tile_full_input_to_target_data_ocean_model(reservoir,grid,trainingdata(:,model_parameters%discardlength/model_parameters%timestep_slab+(batch_number-1)*m+1:batch_number*m+model_parameters%discardlength/model_parameters%timestep_slab),targetdata)

   allocate(temp(reservoir%chunk_size_prediction,n))
   temp = 0.0_dp

   temp = matmul(targetdata,transpose(reservoir%augmented_states))
   !TODO make this matmul DGEMM

   reservoir%states_x_trainingdata_aug = reservoir%states_x_trainingdata_aug + temp

   deallocate(temp)
   deallocate(targetdata)

   allocate(temp(n,n))

   call DGEMM('N','N',n,n,m,alpha,reservoir%augmented_states,n,transpose(reservoir%augmented_states),m,beta,temp,n)
   reservoir%states_x_states_aug = reservoir%states_x_states_aug + temp

   deallocate(temp)

   return
end subroutine

subroutine chunking_matmul_hybrid(reservoir,model_parameters,grid,batch_number,trainingdata,imperfect_model)
   use mod_utilities, only : gaussian_noise
   use resdomain, only : tile_full_input_to_target_data_ocean_model
   use mpires

   type(reservoir_type), intent(inout)     :: reservoir
   type(model_parameters_type), intent(in) :: model_parameters
   type(grid_type), intent(in)             :: grid

   integer, intent(in)          :: batch_number

   real(kind=dp), intent(in)    :: trainingdata(:,:)
   real(kind=dp), intent(in)    :: imperfect_model(:,:)

   real(kind=dp), allocatable   :: temp(:,:), targetdata(:,:)
   real(kind=dp), parameter     :: alpha=1.0, beta=0.0

   integer                      :: n, m, l

   n = size(reservoir%augmented_states,1)
   m = size(reservoir%augmented_states,2)

   reservoir%augmented_states(1:reservoir%chunk_size_prediction,:) = imperfect_model(:,model_parameters%discardlength/model_parameters%timestep_slab+(batch_number-1)*m+1:batch_number*m+model_parameters%discardlength/model_parameters%timestep_slab)
   reservoir%augmented_states(reservoir%chunk_size_prediction+1:reservoir%n+reservoir%chunk_size_prediction,:) = reservoir%states

   call tile_full_input_to_target_data_ocean_model(reservoir,grid,trainingdata(:,model_parameters%discardlength/model_parameters%timestep_slab+(batch_number-1)*m+1:batch_number*m+model_parameters%discardlength/model_parameters%timestep_slab),targetdata)

   allocate(temp(reservoir%chunk_size_prediction,n))
   temp = 0.0_dp
   
   temp = matmul(targetdata,transpose(reservoir%augmented_states))
   !TODO make this matmul DGEMM
    
   reservoir%states_x_trainingdata_aug = reservoir%states_x_trainingdata_aug + temp
   deallocate(temp)
   deallocate(targetdata)
 
   allocate(temp(n,n))
   call DGEMM('N','N',n,n,m,alpha,reservoir%augmented_states,n,transpose(reservoir%augmented_states),m,beta,temp,n)
   !call DGEMM('N','T',n,n,m,alpha,reservoir%augmented_states,n,reservoir%augmented_states,m,beta,temp,n)
   !print *, 'slab shape(reservoir%states_x_states_aug),shape(temp)',shape(reservoir%states_x_states_aug),shape(temp)
   reservoir%states_x_states_aug = reservoir%states_x_states_aug + temp
   deallocate(temp)
   
   return 
end subroutine  

subroutine write_controller_file(model_parameters)
   type(model_parameters_type), intent(in) :: model_parameters

   character(len=:), allocatable :: file_path

   file_path = '/scratch/user/troyarcomano/ML_SPEEDY_WEIGHTS/'//trim(model_parameters%trial_name)//'_controller_file.txt'
 
   open (10, file=file_path, status='unknown')

   ! write to file
   write(10,*)"-----------------------------------------------------------"
   write(10,*)
   write(10,*)"num_vert_levels:",model_parameters%num_vert_levels
   write(10,*)"-----------------------------------------------------------"

   ! close file
   close(10) 

end subroutine 

subroutine write_trained_res(reservoir,model_parameters,grid)
  use mod_io, only : write_netcdf_2d_non_met_data, write_netcdf_1d_non_met_data_int, write_netcdf_1d_non_met_data_real

  type(reservoir_type), intent(in) :: reservoir
  type(model_parameters_type), intent(in) :: model_parameters
  type(grid_type), intent(in)             :: grid

  character(len=:), allocatable :: file_path
  character(len=4) :: worker_char

  file_path = '/scratch/user/troyarcomano/ML_SPEEDY_WEIGHTS/'

  write(worker_char,'(i0.4)') reservoir%assigned_region

  print *, 'writing ML ocean',reservoir%assigned_region
  call write_netcdf_2d_non_met_data(reservoir%win,'win',file_path//'worker_'//worker_char//'_ocean_'//trim(model_parameters%trial_name)//'.nc','unitless','win_x','win_y')
  call write_netcdf_2d_non_met_data(reservoir%wout,'wout',file_path//'worker_'//worker_char//'_ocean_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')

  call write_netcdf_1d_non_met_data_int(reservoir%rows,'rows',file_path//'worker_'//worker_char//'_ocean_'//trim(model_parameters%trial_name)//'.nc','unitless','rows_x')
  call write_netcdf_1d_non_met_data_int(reservoir%cols,'cols',file_path//'worker_'//worker_char//'_ocean_'//trim(model_parameters%trial_name)//'.nc','unitless','cols_x')

  call write_netcdf_1d_non_met_data_real(reservoir%vals,'vals',file_path//'worker_'//worker_char//'_ocean_'//trim(model_parameters%trial_name)//'.nc','unitless','vals_x')

  call write_netcdf_1d_non_met_data_real(grid%mean,'mean',file_path//'worker_'//worker_char//'_ocean_'//trim(model_parameters%trial_name)//'.nc','unitless','mean_x')
  call write_netcdf_1d_non_met_data_real(grid%std,'std',file_path//'worker_'//worker_char//'_ocean_'//trim(model_parameters%trial_name)//'.nc','unitless','std_x')

end subroutine

subroutine trained_ocean_reservoir_prediction(reservoir,model_parameters,grid,reservoir_atmo,grid_atmo)
  use mod_linalg, only : mklsparse
  use mod_io, only : read_trained_ocean_res

  type(reservoir_type), intent(inout)     :: reservoir, reservoir_atmo
  type(model_parameters_type), intent(inout) :: model_parameters
  type(grid_type), intent(inout)          :: grid, grid_atmo

  integer :: mean_std_length

  integer :: sst_res_input_size, i, counter

  call read_trained_ocean_res(reservoir,model_parameters,grid)

  if(reservoir%sst_bool_input) then
     print *, 'has sst',reservoir%assigned_region
     reservoir%sst_bool_input = .True.
     reservoir%sst_bool_prediction = .True.
   else
     print *, 'doesnt have sst',reservoir%assigned_region
     reservoir%sst_bool_input = .False.
     reservoir%sst_bool_prediction = .False.
   endif

   if(reservoir%sst_bool_prediction) then

     if((reservoir_atmo%sst_climo_bool).and.(reservoir_atmo%sst_bool_prediction)) then
        reservoir%sst_climo_input = .False.!.True.
     else
        reservoir%sst_climo_input = .False.
     endif

     if(reservoir_atmo%precip_bool) then
         reservoir%precip_input_bool = .False.
     else
         reservoir%precip_input_bool = .False.
     endif

     sst_res_input_size = grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars + grid_atmo%inputxchunk*grid_atmo%inputychunk + grid_atmo%inputxchunk*grid_atmo%inputychunk + grid_atmo%inputxchunk*grid_atmo%inputychunk

     print *, 'grid_atmo%sst_start,grid_atmo%sst_end',grid_atmo%sst_start,grid_atmo%sst_end

     grid%atmo3d_start = grid_atmo%atmo3d_start
     grid%atmo3d_end = grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars!grid_atmo%atmo3d_end

     grid%logp_start = grid%atmo3d_end + 1!grid_atmo%logp_start
     grid%logp_end = grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars + grid_atmo%inputxchunk*grid_atmo%inputychunk!grid_atmo%logp_end

     grid%sst_start = grid%logp_end + 1!grid_atmo%sst_start
     grid%sst_end = grid%sst_start + grid_atmo%inputxchunk*grid_atmo%inputychunk - 1!grid_atmo%sst_end

     grid%tisr_start = grid%sst_end + 1
     grid%tisr_end = grid%tisr_start + grid_atmo%inputxchunk*grid_atmo%inputychunk - 1

     grid%sst_mean_std_idx = grid_atmo%sst_mean_std_idx

     allocate(reservoir%atmo_training_data_idx(sst_res_input_size))
     counter = 0
     do i=grid_atmo%atmo3d_end-grid_atmo%inputxchunk*grid_atmo%inputychunk*reservoir_atmo%local_predictvars+1,grid_atmo%logp_end
        counter = counter + 1
        reservoir%atmo_training_data_idx(counter) = i
     enddo
     do i=grid_atmo%sst_start,grid_atmo%sst_end
        counter = counter + 1
        reservoir%atmo_training_data_idx(counter) = i
     enddo
     do i=grid_atmo%tisr_start,grid_atmo%tisr_end
        counter = counter + 1
        reservoir%atmo_training_data_idx(counter) = i
     enddo

     call initialize_slab_ocean_model(reservoir,grid,model_parameters)

     call mklsparse(reservoir)
  endif
end subroutine

subroutine read_ohtc_parallel_training(reservoir,model_parameters,grid,ohtc_var)
   use mpires, only : mpi_res
   use mod_io, only : read_3d_file_parallel
   use mod_calendar 

   type(reservoir_type), intent(inout)     :: reservoir
   type(model_parameters_type), intent(in) :: model_parameters
   type(grid_type), intent(inout)          :: grid
 
   real(kind=dp), allocatable, intent(out) :: ohtc_var(:,:,:)

   !Local
   type(calendar_type) :: ohtc_calendar


   integer :: start_index, end_index

   character :: ohtc_file 

   ohtc_file = '/scratch/user/troyarcomano/ORAS5/regridded_sohtc300_control_monthly_highres_2D_CONS_v0.1_hourly.nc'

   !Starting date of the ohtc data
   call initialize_calendar(ohtc_calendar,1979,1,16,0)

   call get_current_time_delta_hour(ohtc_calendar,0)

   call get_current_time_delta_hour(calendar,model_parameters%traininglength+model_parameters%synclength)

   call time_delta_between_two_dates_datetime_type(ohtc_calendar,calendar,start_index) 

   allocate(ohtc_var(grid%inputxchunk,grid%inputychunk,model_parameters%traininglength+model_parameters%synclength+100))

   call read_3d_file_parallel(ohtc_file,'sohtc300',mpi_res,grid,ohtc_var,start_index,1,model_parameters%traininglength+model_parameters%synclength+100) 
end subroutine 

subroutine read_ohtc_parallel_prediction(reservoir,model_parameters,grid,ohtc_var)
   use mpires, only : mpi_res
   use mod_io, only : read_3d_file_parallel
   use mod_calendar

   type(reservoir_type), intent(inout)     :: reservoir
   type(model_parameters_type), intent(in) :: model_parameters
   type(grid_type), intent(inout)          :: grid

   real(kind=dp), allocatable, intent(out) :: ohtc_var(:,:,:)

   !Local
   type(calendar_type) :: ohtc_calendar


   integer :: start_index, end_index

   character :: ohtc_file

   ohtc_file = '/scratch/user/troyarcomano/ORAS5/regridded_sohtc300_control_monthly_highres_2D_CONS_v0.1_hourly.nc'

   !Starting date of the ohtc data
   call initialize_calendar(ohtc_calendar,1979,1,16,0)

   call get_current_time_delta_hour(ohtc_calendar,0)

   call get_current_time_delta_hour(calendar,model_parameters%traininglength+model_parameters%prediction_markers(model_parameters%current_trial_number))

   call time_delta_between_two_dates_datetime_type(ohtc_calendar,calendar,start_index)

   allocate(ohtc_var(grid%inputxchunk,grid%inputychunk,model_parameters%synclength+100))

   call read_3d_file_parallel(ohtc_file,'sohtc300',mpi_res,grid,ohtc_var,start_index,1,model_parameters%synclength+100)

end subroutine
 
end module 
