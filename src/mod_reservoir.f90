module mod_reservoir 
  USE, INTRINSIC :: IEEE_ARITHMETIC
  use MKL_SPBLAS
  use mod_utilities, only : dp, main_type, reservoir_type, grid_type, model_parameters_type

  implicit none

  integer :: global_time_step
  logical :: write_training_weights
contains

subroutine initialize_model_parameters(model_parameters,processor,num_of_procs)
   use mpires, only : distribute_prediction_marker

   !bunch of reservoir parameters
   type(model_parameters_type)     :: model_parameters
   integer, intent(in) :: processor, num_of_procs

   model_parameters%number_of_regions = 1152!*4 

   model_parameters%ml_only = .False.
   
   write_training_weights = .True.

   model_parameters%num_predictions = 40
   model_parameters%trial_name = '6000_20_20_20_sigma0.5_beta_res0.001_beta_model_1.0_prior_0.0_overlap1_vertlevel_1_precip_epsilon0.001_ohtc_multiple_leakage_test_oceantimestep_24hr_'!2kbias_10_year_then_platue_speedy_atmo_only' !14d_0.9rho_10noise_beta0.001_20years'  
   !model_parameters%trial_name = '6000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_4_vertlap_6_slab_ocean_model_true_precip_true'
   !'4000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_4_vertlap_2_full_timestep_1'
   !model_parameters%trial_name = '4000_20_20_20_beta_res0.01_beta_model_1.0_prior_0.0_overlap1_vertlevels_4_vertlap_2_full_test_climate_all_tisr_longer'
   model_parameters%trial_name_extra_end = ''!'climo_2kbias_10_year_then_platue_speedy_bc_atmo_no_ice_2k_sst_mean_20std_increase_'

   model_parameters%discardlength = 24*10!7
   model_parameters%traininglength = 227760 - 24*10 !- 40*24!166440 - 24*10  !87600*2+24*10!3+24*10!188280 !254040 !81600!188280!0!0!0!166600!81600 !00!58000!67000!77000
   model_parameters%predictionlength = 8760*2!*70!8760*70!8760*3!1 + 24*5!8760*30 + 24*5!504!8760*11 + 24*5 !504!0
   model_parameters%synclength = 24*14*2*3!*4 + 3*24!24*14*2 !+ 180*24
   model_parameters%timestep = 6!1 !6
   model_parameters%timestep_slab = 24!*14!*2!*7

   global_time_step = model_parameters%timestep

   model_parameters%slab_ocean_model_bool = .True.
   model_parameters%train_on_sst_anomalies = .False.

   model_parameters%ohtc_bool_input = .True.

   model_parameters%non_stationary_ocn_climo = .False.
   model_parameters%final_sst_bias = 2.0

   model_parameters%precip_bool = .True. !.True.  
   model_parameters%precip_epsilon = 0.001!0.0005

   model_parameters%timeofday_bool = .False.
 
   model_parameters%full_predictvars = 4
   model_parameters%full_heightlevels = 8

   model_parameters%num_vert_levels = 1
   model_parameters%vert_loc_overlap = 0!6!1!2

   model_parameters%overlap = 1

   model_parameters%irank = processor
   model_parameters%numprocs = num_of_procs

   model_parameters%noisy = .True.

   model_parameters%regional_vary = .True.  

   model_parameters%using_prior = .True.

   model_parameters%model_noise = 0.0_dp
  
   model_parameters%outvec_component_contribs = .False.

   call distribute_prediction_marker(model_parameters)

end subroutine 
  

subroutine allocate_res_new(reservoir,grid,model_parameters)
  !Routine to allocate all of the reservoir arrays

  use resdomain, only : set_reservoir_by_region

  type(reservoir_type), intent(inout)     :: reservoir
  type(grid_type), intent(inout)          :: grid 
  type(model_parameters_type), intent(in) :: model_parameters

  integer :: nodes_per_input

  reservoir%m = 6000!6000

  reservoir%deg = 6
  reservoir%radius = 0.9
  reservoir%beta_res = 0.001_dp
  reservoir%beta_model = 1.0_dp !0.01_dp!1.0_dp
  reservoir%sigma = 0.5_dp!/6.0_dp

  reservoir%leakage = 1.0_dp!/3.0_dp!1.0!1.0_dp/12.0_dp!6.0_dp
   
  reservoir%prior_val = 0.0_dp
 
  reservoir%density = reservoir%deg/real(reservoir%m,kind=dp)

  call set_reservoir_by_region(reservoir,grid)

  if(reservoir%logp_bool) then
    reservoir%logp_size_input = grid%inputxchunk*grid%inputychunk !((grid%resxchunk+grid%overlap*2)*(grid%resychunk+1*grid%overlap))
  else  
    reservoir%logp_size_input = 0
  endif 

  if(reservoir%sst_bool_input) then
    reservoir%sst_size_input = grid%inputxchunk*grid%inputychunk
  else
    reservoir%sst_size_input = 0
  endif
   
  if(reservoir%logp_bool) then
    reservoir%logp_size_res = grid%resxchunk*grid%resychunk
  else
    reservoir%logp_size_res = 0
  endif

  if(reservoir%precip_input_bool) then
    reservoir%precip_size_res = grid%resxchunk*grid%resychunk
  else
    reservoir%precip_size_res = 0
  endif

  if(reservoir%precip_input_bool) then
    reservoir%precip_size_input = grid%inputxchunk*grid%inputychunk
  else
    reservoir%precip_size_input = 0
  endif

  if(reservoir%sst_bool_input) then
    reservoir%sst_size_res = grid%resxchunk*grid%resychunk
  else
    reservoir%sst_size_res = 0
  endif

  if(reservoir%tisr_input_bool) then
    reservoir%tisr_size_res = grid%resxchunk*grid%resychunk
  else 
    reservoir%tisr_size_res = 0
  endif 

  if(reservoir%tisr_input_bool) then
    reservoir%tisr_size_input = grid%inputxchunk*grid%inputychunk!((grid%resxchunk+grid%overlap*2)*(grid%resychunk+1*grid%overlap))
  else
    reservoir%tisr_size_input = 0
  endif 

  reservoir%chunk_size = grid%resxchunk*grid%resychunk*reservoir%local_predictvars*grid%reszchunk + reservoir%logp_size_res + reservoir%precip_size_res
  reservoir%chunk_size_prediction = grid%resxchunk*grid%resychunk*reservoir%local_predictvars*grid%reszchunk + reservoir%logp_size_res + reservoir%precip_size_res
  reservoir%chunk_size_speedy = grid%resxchunk*grid%resychunk*reservoir%local_predictvars*grid%reszchunk + reservoir%logp_size_res

  if(model_parameters%ml_only) then 
    reservoir%chunk_size_speedy = 0
  endif 
  
  reservoir%locality = 0

  reservoir%locality = grid%inputxchunk*grid%inputychunk*grid%inputzchunk*reservoir%local_predictvars + reservoir%logp_size_input + reservoir%precip_size_input + reservoir%tisr_size_input + reservoir%sst_size_input - reservoir%chunk_size

  print *, 'grid%inputxchunk,grid%inputychunk,grid%inputzchunk,reservoir%local_predictvars',grid%inputxchunk,grid%inputychunk,grid%inputzchunk,reservoir%local_predictvars
  print *,'grid%inputxchunk*grid%inputychunk*grid%inputzchunk*reservoir%local_predictvars,reservoir%logp_size_input, reservoir%precip_size_input,reservoir%tisr_size_input,reservoir%sst_size_input',grid%inputxchunk*grid%inputychunk*grid%inputzchunk*reservoir%local_predictvars,reservoir%logp_size_input, reservoir%precip_size_input,reservoir%tisr_size_input,reservoir%sst_size_input
  nodes_per_input = NINT(dble(reservoir%m)/(dble(reservoir%chunk_size)+dble(reservoir%locality)))
  reservoir%n = nodes_per_input*(reservoir%chunk_size+reservoir%locality)
  print *, 'reservoir%chunk_size,reservoir%locality,nodes_per_input,n',reservoir%chunk_size,reservoir%locality, nodes_per_input,reservoir%n
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
   
   call makesparse(reservoir) 

   call sparse_eigen(reservoir,reservoir%n*10,6,eigs)
 
   allocate(newvals(reservoir%k)) 
   newvals = (reservoir%vals/eigs)*reservoir%radius
   reservoir%vals = newvals
   
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

subroutine train_reservoir(reservoir,grid,model_parameters)
   use mod_utilities, only : init_random_seed

   type(reservoir_type), intent(inout)        :: reservoir
   type(model_parameters_type), intent(inout) :: model_parameters
   type(grid_type), intent(inout)             :: grid

   integer :: q,i,num_inputs,j,k
   integer :: un_noisy_sync
   integer :: betas_res, betas_model,priors
   integer :: vert_loop
    
   real(kind=dp), allocatable :: ip(:),rand(:),average
   real(kind=dp), allocatable :: test_beta_res(:), test_beta_model(:), test_priors(:)
   real(kind=dp), allocatable :: states_x_states_original_copy(:,:)

   character(len=:), allocatable :: base_trial_name
   character(len=50) :: beta_res_char,beta_model_char,prior_char
  
   if(grid%bottom) then
    reservoir%logp_bool = .True.
    reservoir%tisr_input_bool = .True.
    grid%logp_bool = .True.

    reservoir%sst_bool = model_parameters%slab_ocean_model_bool
    reservoir%sst_climo_bool = .False.!.True. !.False.

    reservoir%precip_input_bool = model_parameters%precip_bool
    reservoir%precip_bool = model_parameters%precip_bool

    reservoir%m = 6000
    
  else
    reservoir%sst_climo_bool = .False.
    reservoir%logp_bool = .False.
    reservoir%tisr_input_bool = .True.
    reservoir%sst_bool = .False.
    reservoir%precip_input_bool = .False.
    reservoir%precip_bool = .False.
  endif
 
   call get_training_data(reservoir,model_parameters,grid,1)
 
   !NOTE moving this to get_training_data 
   !call allocate_res_new(reservoir,grid,model_parameters)
 
   call gen_res(reservoir)

   q = reservoir%n/reservoir%reservoir_numinputs
   !q = reservoir%reservoir_numinputs

   if(reservoir%assigned_region == 0) print *, 'q',q,'n',reservoir%n,'num inputs',reservoir%reservoir_numinputs

   allocate(ip(q))
   allocate(rand(q))

   reservoir%win = 0.0_dp

   do i=1,reservoir%reservoir_numinputs

      call random_number(rand)

      ip = (-1d0 + 2*rand) 
      
      reservoir%win((i-1)*q+1:i*q,i) = reservoir%sigma*ip
      !reservoir%win(:,i) = reservoir%sigma*ip
   enddo
   
   deallocate(rand)
   deallocate(ip) 

   print *,'starting reservoir_layer'   
  
   call initialize_chunk_training(reservoir,model_parameters) 
   
   do i=1,model_parameters%timestep
      print *, 'loop number',i
      !if(reservoir%assigned_region == 954) print *, 'reservoir%trainingdata(eservoir_numinputs,1:40)',reservoir%trainingdata(reservoir%reservoir_numinputs,1:40)
      if(reservoir%assigned_region == 954 .and. .not. model_parameters%ml_only) print *, 'reservoir%imperfect_model_states(:,i)',reservoir%imperfect_model_states(:,i)
      if(reservoir%assigned_region == 954) print *, 'reservoir%trainingdata(:,i)', reservoir%trainingdata(:,i)

      if(model_parameters%ml_only) then
        call reservoir_layer_chunking_ml(reservoir,model_parameters,grid,reservoir%trainingdata(:,i:model_parameters%traininglength:model_parameters%timestep))
      else
        call reservoir_layer_chunking_hybrid(reservoir,model_parameters,grid,reservoir%trainingdata(:,i:model_parameters%traininglength:model_parameters%timestep),reservoir%imperfect_model_states(:,i:model_parameters%traininglength:model_parameters%timestep))
      endif 

   enddo

   print *, 'fitting',reservoir%assigned_region

   if(model_parameters%ml_only) then
     call fit_chunk_ml(reservoir,model_parameters,grid)
   else
     call fit_chunk_hybrid(reservoir,model_parameters,grid)
   endif

   !Get training error of last batch
   call training_error(reservoir,model_parameters,grid,reservoir%trainingdata(:,model_parameters%timestep:model_parameters%traininglength:model_parameters%timestep),reservoir%imperfect_model_states(:,model_parameters%timestep:model_parameters%traininglength:model_parameters%timestep))
 
   if((model_parameters%slab_ocean_model_bool).and.(grid%bottom)) then
     if(allocated(reservoir%imperfect_model_states)) deallocate(reservoir%imperfect_model_states)
   else
     deallocate(reservoir%trainingdata)
     if(allocated(reservoir%imperfect_model_states)) deallocate(reservoir%imperfect_model_states)
   endif

 
   print *, 'cleaning up', reservoir%assigned_region
   call clean_batch(reservoir)
 
end subroutine 

subroutine training_error(reservoir,model_parameters,grid,trainingdata,imperfect_model_states)
   use mod_utilities, only: unstandardize_data_4d !, unstandardize_data_2d
   use resdomain, only : tile_full_input_to_target_data, tile_4d_and_logp_state_vec_res1d
   use mpires, only : killmpi, mpi_res
   use mod_io, only : write_netcdf_parallel_mpi

   type(reservoir_type), intent(inout)     :: reservoir
   type(model_parameters_type), intent(in) :: model_parameters
   type(grid_type), intent(in)             :: grid
   real(kind=dp), intent(in)               :: trainingdata(:,:)
   real(kind=dp), intent(in)               :: imperfect_model_states(:,:)

   real(kind=dp), allocatable :: x_augment(:)
   real(kind=dp), allocatable :: prediction(:,:)
   real(kind=dp), allocatable :: targetdata(:,:)
   real(kind=dp), allocatable :: grid5d_truth(:,:,:,:,:), grid3d_truth(:,:,:)
   real(kind=dp), allocatable :: grid5d_prediction(:,:,:,:,:), grid3d_prediction(:,:,:)
   real(kind=dp), allocatable :: grid5d_difference(:,:,:,:,:), grid3d_difference(:,:,:)

   integer :: time_length
   integer :: i
   integer :: batch_number

   time_length = size(reservoir%states,2)
   batch_number = (size(trainingdata,2)-model_parameters%discardlength/model_parameters%timestep)/reservoir%batch_size !model_parameters%traininglength/model_parameters%timestep/reservoir%batch_size
   print *, 'Calculated number of batches in training_error routine: ', batch_number
   allocate(prediction(reservoir%chunk_size_prediction,time_length))
   allocate(x_augment(reservoir%n+reservoir%chunk_size_speedy))

   do i=1, time_length
      x_augment(1:reservoir%chunk_size_speedy) = imperfect_model_states(:,model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*time_length+i)
      x_augment(reservoir%chunk_size_speedy+1:reservoir%chunk_size_speedy+reservoir%n) = reservoir%states(:,i)
      prediction(:,i) = matmul(reservoir%wout,x_augment)
   enddo

   call tile_full_input_to_target_data(reservoir,grid,trainingdata(1:grid%predict_end,model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*time_length+1:batch_number*time_length+model_parameters%discardlength/model_parameters%timestep),targetdata)

   allocate(grid5d_truth(reservoir%local_predictvars,grid%resxchunk,grid%resychunk,grid%reszchunk,time_length))
   allocate(grid3d_truth(grid%resxchunk,grid%resychunk,time_length))

   allocate(grid5d_prediction(reservoir%local_predictvars,grid%resxchunk,grid%resychunk,grid%reszchunk,time_length))
   allocate(grid3d_prediction(grid%resxchunk,grid%resychunk,time_length))

   allocate(grid5d_difference(reservoir%local_predictvars,grid%resxchunk,grid%resychunk,grid%reszchunk,time_length))
   allocate(grid3d_difference(grid%resxchunk,grid%resychunk,time_length))

   do i=1,time_length
      call tile_4d_and_logp_state_vec_res1d(reservoir,grid%number_of_regions,targetdata(:,i),reservoir%assigned_region,grid5d_truth(:,:,:,:,i),grid3d_truth(:,:,i))
      call tile_4d_and_logp_state_vec_res1d(reservoir,grid%number_of_regions,prediction(:,i),reservoir%assigned_region,grid5d_prediction(:,:,:,:,i),grid3d_prediction(:,:,i))
   enddo

   if(reservoir%assigned_region==1) then
    print *,'target data ',targetdata
    print *,'prediction data ',prediction
   endif
   do i=1,time_length
      call unstandardize_data_4d(reservoir,grid5d_truth(:,:,:,:,i),grid%mean,grid%std)
      call unstandardize_data_4d(reservoir,grid5d_prediction(:,:,:,:,i),grid%mean,grid%std)
      !call unstandardize_data_2d(grid3d_truth(:,:,i),grid%mean(grid%logp_mean_std_idx),grid%std(grid%logp_mean_std_idx))
      !call unstandardize_data_2d(grid3d_prediction(:,:,i),grid%mean(grid%logp_mean_std_idx),grid%std(grid%logp_mean_std_idx))
   enddo
   if(reservoir%assigned_region==1) then
    print *,'target data unstandardized ',targetdata
    print *,'prediction data unstandardized ',prediction
   endif
   grid5d_difference = 0.0_dp
   grid3d_difference = 0.0_dp

   grid5d_difference = grid5d_prediction - grid5d_truth
   grid3d_difference = grid3d_prediction - grid3d_truth

   call write_netcdf_parallel_mpi(grid%res_xstart,grid%res_ystart,grid5d_difference,grid3d_difference,1,'training_error_'//model_parameters%trial_name//'.nc',mpi_res,.True.)
   call write_netcdf_parallel_mpi(grid%res_xstart,grid%res_ystart,grid5d_prediction,grid3d_prediction,1,'training_prediction_'//model_parameters%trial_name//'.nc',mpi_res,.True.)
   call write_netcdf_parallel_mpi(grid%res_xstart,grid%res_ystart,grid5d_truth,grid3d_truth,1,'training_truth_'//model_parameters%trial_name//'.nc',mpi_res,.True.)
   deallocate(reservoir%states)
end subroutine


subroutine get_training_data(reservoir,model_parameters,grid,loop_index)
   use mod_utilities, only : era_data_type, speedy_data_type, &
                             standardize_sst_data_3d, & 
                             standardize_data_given_pars_5d_logp_tisr, &
                             standardize_data_given_pars_5d_logp, &
                             standardize_data_given_pars5d, &
                             standardize_data, &
                             total_precip_over_a_period, &
                             standardize_data_3d, &
                             unstandardize_data_2d, &
                             e_constant, &
                             rolling_average_over_a_period

   use mod_calendar
   use speedy_res_interface, only : read_era, read_model_states
   use resdomain, only : standardize_speedy_data

   type(reservoir_type), intent(inout)        :: reservoir
   type(model_parameters_type), intent(inout) :: model_parameters
   type(grid_type), intent(inout)             :: grid

   integer, intent(in) :: loop_index

   type(era_data_type)    :: era_data
   type(speedy_data_type) :: speedy_data

   integer :: mean_std_length

   reservoir%local_predictvars = model_parameters%full_predictvars
   reservoir%local_heightlevels_input = grid%inputzchunk

   reservoir%local_heightlevels_res = grid%reszchunk

   call initialize_calendar(calendar,1981,1,1,0)

   call get_current_time_delta_hour(calendar,model_parameters%discardlength+model_parameters%traininglength+model_parameters%synclength)

   print *, 'reading era states' 

   !Read data in stride and whats only needed for this loop of training
   call read_era(reservoir,grid,model_parameters,calendar%startyear,calendar%currentyear,era_data)

   !Match units for specific humidity
   era_data%eravariables(4,:,:,:,:) = era_data%eravariables(4,:,:,:,:)*1000.0_dp
   where (era_data%eravariables(4,:,:,:,:) < 0.000001)
    era_data%eravariables(4,:,:,:,:) = 0.000001_dp
   end where
   !print *, 'shape(era_data%eravariables)',shape(era_data%eravariables)
   !print *, 'era_data',era_data%eravariables(1,:,:,:,10:11)
   !Make sure tisr doesnt have zeroes 
   if(reservoir%tisr_input_bool) then
    where(era_data%era_tisr < 0.0_dp)
      era_data%era_tisr = 0.0_dp
    end where
  endif

  if(reservoir%precip_bool) then
    !era_data%era_precip = era_data%era_precip * 39.3701
    where(era_data%era_precip < 0.0_dp)
      era_data%era_precip = 0.0_dp
    end where
    if(reservoir%assigned_region == 954) print *, 'era_data%era_precip(1,1,100) before',era_data%era_precip(1,1,100:150)
    call total_precip_over_a_period(era_data%era_precip,model_parameters%timestep)
    if(reservoir%assigned_region == 954) print *, 'era_data%era_precip(1,1,100) after average',era_data%era_precip(1,1,100:150)
    if(reservoir%assigned_region == 954) print *, 'era_data%era_precip/model_parameters%precip_epsilon',era_data%era_precip(1,1,100:150)/model_parameters%precip_epsilon
    era_data%era_precip =  log(1 + era_data%era_precip/model_parameters%precip_epsilon)
    if(reservoir%assigned_region == 954) print *, 'era_data%era_precip(1,1,100:150) after',era_data%era_precip(1,1,100:150) 
  endif

  if(reservoir%sst_bool .and. .not. model_parameters%train_on_sst_anomalies) then   
    where(era_data%era_sst < 272.0_dp)
      era_data%era_sst = 272.0_dp
    end where
  endif

  if(reservoir%assigned_region == 954) then
    print *, 'era max min temp before',maxval(era_data%eravariables(1,:,:,:,:)),minval(era_data%eravariables(1,:,:,:,:))
    print *, 'era max min u-wind before',maxval(era_data%eravariables(2,:,:,:,:)),minval(era_data%eravariables(2,:,:,:,:))
    print *, 'era max min v-wind before',maxval(era_data%eravariables(3,:,:,:,:)),minval(era_data%eravariables(3,:,:,:,:))
    print *, 'era max min sp before',maxval(era_data%eravariables(4,:,:,:,:)),minval(era_data%eravariables(4,:,:,:,:))
    if(reservoir%logp_bool) print *, 'era max min logp before',maxval(era_data%era_logp),minval(era_data%era_logp)

    if(reservoir%tisr_input_bool) print *, 'era max min tisr before',maxval(era_data%era_tisr),minval(era_data%era_tisr)
    if(reservoir%sst_bool)  print *, 'era max min sst before',maxval(era_data%era_sst),minval(era_data%era_sst)
    if(reservoir%precip_bool) print *, 'era max min precip rate before',maxval(era_data%era_precip), minval(era_data%era_precip)
  endif
  !Get mean and standard deviation for the first stride of data and use those
  !values for the rest of the program
  if(loop_index == 1) then   
     !Standardize each variable using local std and mean and save the std and
     !mean

     !Get number of height levels * vars + 2d variables 
     mean_std_length = model_parameters%full_predictvars*grid%inputzchunk
     if(reservoir%logp_bool) then
       mean_std_length = mean_std_length + 1
       grid%logp_mean_std_idx = mean_std_length
     endif 

     if(reservoir%tisr_input_bool) then
        mean_std_length = mean_std_length + 1
        grid%tisr_mean_std_idx = mean_std_length
     endif

     if(reservoir%precip_bool)  then
       mean_std_length = mean_std_length + 1
       grid%precip_mean_std_idx = mean_std_length
     endif 

     if(reservoir%sst_bool) then
        mean_std_length = mean_std_length + 1 
        grid%sst_mean_std_idx = mean_std_length 
     endif 

     allocate(grid%mean(mean_std_length),grid%std(mean_std_length))

     if((reservoir%tisr_input_bool).and.(reservoir%logp_bool)) then
        !grid%tisr_mean_std_idx = reservoir%local_predictvars*reservoir%local_heightlevels_input+2 
        !grid%logp_mean_std_idx = reservoir%local_predictvars*reservoir%local_heightlevels_input+1

        print *, 'mean_std_length',mean_std_length
        call standardize_data(reservoir,era_data%eravariables,era_data%era_logp,era_data%era_tisr,grid%mean,grid%std)
        !call standardize_data(reservoir,grid,era_data%eravariables,era_data%era_logp,era_data%era_tisr,grid%mean,grid%std)
     elseif(reservoir%logp_bool) then
        !grid%logp_mean_std_idx = reservoir%local_predictvars*reservoir%local_heightlevels_input+1

        call standardize_data(reservoir,era_data%eravariables,era_data%era_logp,grid%mean,grid%std)
     elseif(reservoir%tisr_input_bool) then
        !grid%tisr_mean_std_idx = reservoir%local_predictvars*reservoir%local_heightlevels_input+1

        call standardize_data(reservoir,era_data%eravariables,era_data%era_tisr,grid%mean,grid%std)
     else
        call standardize_data(reservoir,era_data%eravariables,grid%mean,grid%std)
     endif 
   
     if(reservoir%sst_bool) then
        !grid%sst_mean_std_idx = mean_std_length
        call standardize_sst_data_3d(era_data%era_sst,grid%mean(grid%sst_mean_std_idx),grid%std(grid%sst_mean_std_idx),reservoir%sst_bool_input)
     else 
        reservoir%sst_bool_input = .False.
     endif
  
     if(reservoir%precip_bool) then
       ! if(reservoir%assigned_region == 954) print *, 'era_data%era_precip before',era_data%era_precip(1,1,1:100)
        !call total_precip_over_a_period(era_data%era_precip,model_parameters%timestep)
        !if(reservoir%assigned_region == 954) print *, 'era max min precip rate after summing',maxval(era_data%era_precip), minval(era_data%era_precip)
        !if(reservoir%assigned_region == 954) print *, 'era_data%era_precip after',era_data%era_precip(1,1,1:100)
        call standardize_data_3d(era_data%era_precip,grid%mean(grid%precip_mean_std_idx),grid%std(grid%precip_mean_std_idx))
        if(reservoir%assigned_region == 954) print *, 'precip mean and std',grid%mean(grid%precip_mean_std_idx),grid%std(grid%precip_mean_std_idx)
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
     if(reservoir%sst_bool)  print *, 'era max min sst after',maxval(era_data%era_sst),minval(era_data%era_sst)
     if(reservoir%precip_bool) print *, 'era max min precip rate after',maxval(era_data%era_precip), minval(era_data%era_precip)
     print *, 'res%mean,res%std',grid%mean,grid%std
   endif

   !NOTE moving this here to we can get the training data 
   call allocate_res_new(reservoir,grid,model_parameters)

   !Lets get some training data
   allocate(reservoir%trainingdata(reservoir%reservoir_numinputs,size(era_data%eravariables,5)))

   print *, 'reservoir%reservoir_numinputs',reservoir%reservoir_numinputs,reservoir%assigned_region,grid%level_index

   grid%logp_start = 0
   grid%logp_end = 0
   grid%sst_start = 0
   grid%sst_end = 0
  
   grid%atmo3d_start = 1
   grid%atmo3d_end = model_parameters%full_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk

   grid%predict_start = 1
   grid%predict_end = grid%atmo3d_end
   print *, 'grid%atmo3d_start,grid%atmo3d_end',grid%atmo3d_start,grid%atmo3d_end
   print *, 'shape(reservoir%trainingdata(grid%atmo3d_start:grid%atmo3d_end,:))',shape(reservoir%trainingdata(grid%atmo3d_start:grid%atmo3d_end,:)) 
   print *, 'shape(reshape(era_data%eravariables,(/grid%atmo3d_end,size(era_data%eravariables,5)/)))',shape(reshape(era_data%eravariables,(/grid%atmo3d_end,size(era_data%eravariables,5)/)))
   reservoir%trainingdata(grid%atmo3d_start:grid%atmo3d_end,:) = reshape(era_data%eravariables,(/grid%atmo3d_end,size(era_data%eravariables,5)/))
  
   if(reservoir%logp_bool) then 
     grid%logp_start = grid%atmo3d_end + 1
     grid%logp_end = grid%atmo3d_end + reservoir%logp_size_input 
     reservoir%trainingdata(grid%logp_start:grid%logp_end,:) = reshape(era_data%era_logp,(/reservoir%logp_size_input,size(era_data%eravariables,5)/))

     grid%predict_end = grid%logp_end
   endif 

   if(reservoir%precip_bool) then
     grid%precip_start = grid%atmo3d_end + reservoir%logp_size_input + 1
     grid%precip_end = grid%precip_start + reservoir%precip_size_input - 1
     reservoir%trainingdata(grid%precip_start:grid%precip_end,:) = reshape(era_data%era_precip,(/reservoir%precip_size_input,size(era_data%eravariables,5)/))

     grid%predict_end = grid%precip_end
   endif
    
   if(reservoir%sst_bool_input) then
     grid%sst_start = grid%atmo3d_end + reservoir%logp_size_input + reservoir%precip_size_input + 1
     grid%sst_end =  grid%sst_start + reservoir%sst_size_input - 1
     reservoir%trainingdata(grid%sst_start:grid%sst_end,:) = reshape(era_data%era_sst,(/reservoir%sst_size_input,size(era_data%eravariables,5)/))
   endif 

   if(reservoir%tisr_input_bool) then
     grid%tisr_start = grid%atmo3d_end + reservoir%logp_size_input + reservoir%precip_size_input + reservoir%sst_size_input + 1
     grid%tisr_end = grid%tisr_start + reservoir%tisr_size_input - 1
     reservoir%trainingdata(grid%tisr_start:grid%tisr_end,:) = reshape(era_data%era_tisr,(/reservoir%tisr_size_input,size(era_data%eravariables,5)/))
   endif

   
   if(reservoir%assigned_region == 954) print *, 'reservoir%trainingdata(:,1000)',reservoir%trainingdata(:,1000)
   if(reservoir%assigned_region == 954) print *, 'reservoir%trainingdata(grid%tisr_start,1000:1100)',reservoir%trainingdata(grid%tisr_start,1000:1100)
   deallocate(era_data%eravariables)
   deallocate(era_data%era_logp)

   if(allocated(era_data%era_tisr)) then
     deallocate(era_data%era_tisr)
   endif

   if(allocated(era_data%era_sst)) then
     deallocate(era_data%era_sst)
   endif


   !Portion of the routine for getting speedy (imperfect model) data
   
   if(.not. model_parameters%ml_only) then
     print *, 'reading model states'
     call read_model_states(reservoir,grid,model_parameters,calendar%startyear,calendar%currentyear,speedy_data)
   
     !Lets get imperfect model states
     where(speedy_data%speedyvariables(4,:,:,:,:) < 0.000001)
        speedy_data%speedyvariables(4,:,:,:,:) = 0.000001_dp
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

     if(reservoir%assigned_region == 36)  print *, 'reservoir%trainingdata(grid%tisr_start:grid%tisr_end,1000:1100)', reservoir%trainingdata(grid%tisr_start,1000:1100)
     call standardize_speedy_data(reservoir,grid,speedy_data)

     if(reservoir%assigned_region == 954) print *, 'speedy_data%speedyvariables(:,1,1,1,1) after',speedy_data%speedyvariables(:,1,1,1,1)
     allocate(reservoir%imperfect_model_states(reservoir%chunk_size_speedy,size(speedy_data%speedyvariables,5)))
     reservoir%imperfect_model_states = 0.0_dp

     reservoir%imperfect_model_states(1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk,:) = reshape(speedy_data%speedyvariables,(/reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk,size(speedy_data%speedyvariables,5)/))

     if(reservoir%logp_bool) then
        reservoir%imperfect_model_states(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk+1:reservoir%chunk_size_speedy,:) = reshape(speedy_data%speedy_logp,(/grid%resxchunk*grid%resychunk,size(speedy_data%speedyvariables,5)/))
     endif 

     deallocate(speedy_data%speedyvariables)
     deallocate(speedy_data%speedy_logp)
   endif 
end subroutine 

subroutine get_prediction_data(reservoir,model_parameters,grid,start_index,length)
   use mod_utilities, only : era_data_type, speedy_data_type, &
                             standardize_data_given_pars_5d_logp_tisr, &
                             standardize_data_given_pars_5d_logp, &
                             standardize_data_given_pars5d, &
                             standardize_data, &
                             standardize_data_given_pars3d, &
                             total_precip_over_a_period, &
                             opened_netcdf_type

   use mod_calendar
   use speedy_res_interface, only : read_era, read_model_states, &
                                    read_era_netcdf_opened
   use resdomain, only : standardize_speedy_data

   type(reservoir_type), intent(inout)        :: reservoir
   type(model_parameters_type), intent(inout) :: model_parameters
   type(grid_type), intent(inout)             :: grid

   integer, intent(in) :: start_index,length
   
   integer                :: hours_into_first_year, start_year
   integer                :: start_time_memory_index, end_time_memory_index

   integer                :: max_netcdf_files

   type(era_data_type)    :: era_data
   type(speedy_data_type) :: speedy_data  

   call get_current_time_delta_hour(calendar,start_index)

   call numof_hours_into_year(calendar%currentyear,calendar%currentmonth,calendar%currentday,calendar%currenthour,hours_into_first_year)
 
   start_year = calendar%currentyear

   call get_current_time_delta_hour(calendar,start_index+length) 

   if(model_parameters%num_of_regions_on_proc > 1153) then
     if(.not. allocated(model_parameters%opened_netcdf_files)) then
        max_netcdf_files = (calendar%currentyear - start_year + 1) * 6 
        allocate(model_parameters%opened_netcdf_files(max_netcdf_files))
        model_parameters%opened_netcdf_files(:)%is_opened = .False.
     endif 
     call read_era_netcdf_opened(reservoir,grid,model_parameters,start_year,calendar%currentyear,era_data,model_parameters%opened_netcdf_files,1)
   else 
     !Read data in stride and whats only needed for this loop of training
     !This is the normal routine when we just have 1 region per processor
     !If we have more then we can get a big speed up by keep files opened 
     !but it requires over head and more complication so its not worth doing for
     !1 region per processor
     call read_era(reservoir,grid,model_parameters,start_year,calendar%currentyear,era_data,1)
   endif 
   start_time_memory_index = hours_into_first_year   
   end_time_memory_index = start_time_memory_index + length

   !Match units for specific humidity
   era_data%eravariables(4,:,:,:,:) = era_data%eravariables(4,:,:,:,:)*1000.0_dp

   where (era_data%eravariables(4,:,:,:,:) < 0.000001)
    era_data%eravariables(4,:,:,:,:) = 0.000001_dp
   end where

   !Make sure tisr doesnt have zeroes 
   if(reservoir%tisr_input_bool) then
    where(era_data%era_tisr < 0.0_dp)
      era_data%era_tisr = 0.0_dp
    end where
   endif

   if(reservoir%sst_bool .and. .not. model_parameters%train_on_sst_anomalies) then 
     !NOTE Trying to give a distinction between land and sea ice
     !where(era_data%era_sst < -1) 
     !  era_data%era_sst = 271.0_dp
     !end where 
     !Note maybe we should make it 271.5 
     where(era_data%era_sst < 272.0_dp)
       era_data%era_sst = 272.0_dp
     end where
     era_data%era_sst = era_data%era_sst !+ 4.0_dp
   endif

  if(reservoir%precip_bool) then
    !era_data%era_precip = era_data%era_precip * 39.3701
    where(era_data%era_precip < 0.0_dp)
      era_data%era_precip = 0.0_dp
    end where
    call total_precip_over_a_period(era_data%era_precip,model_parameters%timestep)
    if(reservoir%assigned_region == 954) print *, 'era_data%era_precip(1,1,100) before',era_data%era_precip(1,1,100:150)
    if(reservoir%assigned_region == 954) print *, 'era_data%era_precip/model_parameters%precip_epsilon',era_data%era_precip(1,1,100:150)/model_parameters%precip_epsilon
    era_data%era_precip =  log(1 + era_data%era_precip/model_parameters%precip_epsilon)
    if(reservoir%assigned_region == 954) print *, 'era_data%era_precip(1,1,100:150) after',era_data%era_precip(1,1,100:150)
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

  if(reservoir%sst_bool_input) then
     call standardize_data_given_pars3d(era_data%era_sst,grid%mean(grid%sst_mean_std_idx),grid%std(grid%sst_mean_std_idx))
  endif 

  if(reservoir%precip_bool) then
     call standardize_data_given_pars3d(era_data%era_precip,grid%mean(grid%precip_mean_std_idx),grid%std(grid%precip_mean_std_idx))
  endif
 
  if(allocated(reservoir%predictiondata)) then
     deallocate(reservoir%predictiondata)
  endif 

   !Lets get some prediction data
   allocate(reservoir%predictiondata(reservoir%reservoir_numinputs,length/model_parameters%timestep))

   reservoir%predictiondata(grid%atmo3d_start:grid%atmo3d_end,:) = reshape(era_data%eravariables(:,:,:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep),[grid%atmo3d_end,length/model_parameters%timestep])

   if(reservoir%logp_bool) then 
     reservoir%predictiondata(grid%logp_start:grid%logp_end,:) = reshape(era_data%era_logp(:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep),[reservoir%logp_size_input,length/model_parameters%timestep])
   endif
  
   if(reservoir%precip_bool) then
     reservoir%predictiondata(grid%precip_start:grid%precip_end,:) = reshape(era_data%era_precip(:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep),[reservoir%precip_size_input,length/model_parameters%timestep])
   endif

   if(reservoir%sst_bool_input) then
     reservoir%predictiondata(grid%sst_start:grid%sst_end,:) = reshape(era_data%era_sst(:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep),[reservoir%sst_size_input,length/model_parameters%timestep])
   endif

   if(reservoir%tisr_input_bool) then
     reservoir%predictiondata(grid%tisr_start:grid%tisr_end,:) = reshape(era_data%era_tisr(:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep),[reservoir%tisr_size_input,length/model_parameters%timestep])
   endif

   deallocate(era_data%eravariables)
   deallocate(era_data%era_logp)

   if(allocated(era_data%era_tisr)) then
     deallocate(era_data%era_tisr)
   endif

   if(allocated(era_data%era_sst)) then
     deallocate(era_data%era_sst)
   endif 

   if(allocated(era_data%era_precip)) then
     deallocate(era_data%era_precip)
   endif
 
   print *, 'shape(reservoir%predictiondata) mod_res',shape(reservoir%predictiondata)

   if(.not. model_parameters%ml_only) then
     !Portion of the routine for getting speedy (imperfect model) data
     print *, 'reading model states'
     call read_model_states(reservoir,grid,model_parameters,start_year,calendar%currentyear,speedy_data,1)
   
     !Lets get imperfect model states
     where(speedy_data%speedyvariables(4,:,:,:,:) < 0.000001)
        speedy_data%speedyvariables(4,:,:,:,:) = 0.000001_dp
     end where

     call standardize_speedy_data(reservoir,grid,speedy_data) 

     if(allocated(reservoir%imperfect_model_states)) then
       deallocate(reservoir%imperfect_model_states)
     endif 

     allocate(reservoir%imperfect_model_states(reservoir%chunk_size_speedy,length/model_parameters%timestep))

     reservoir%imperfect_model_states = 0.0_dp

     reservoir%imperfect_model_states(1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk,:) = reshape(speedy_data%speedyvariables(:,:,:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep),[reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk,length/model_parameters%timestep])

     if(reservoir%logp_bool) then
        reservoir%imperfect_model_states(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*grid%reszchunk+1:reservoir%chunk_size_speedy,:) = reshape(speedy_data%speedy_logp(:,:,start_time_memory_index:end_time_memory_index:model_parameters%timestep),[grid%resxchunk*grid%resychunk,length/model_parameters%timestep])
     endif 
     deallocate(speedy_data%speedyvariables)
     deallocate(speedy_data%speedy_logp)
   endif   

end subroutine 

subroutine initialize_prediction(reservoir,model_parameters,grid)
   use mod_calendar
   use mpires, only        : mpi_res
   use mod_io, only        : read_netcdf_3d
   use mod_utilities, only : xgrid, ygrid
   
   type(reservoir_type), intent(inout)        :: reservoir
   type(grid_type), intent(inout)             :: grid 
   type(model_parameters_type), intent(inout) :: model_parameters

   integer :: q,i,num_inputs,j,k
   integer :: un_noisy_sync
   integer :: betas_res, betas_model,priors
   integer :: vert_loop

   real(kind=dp), allocatable :: ip(:),rand(:),average
   real(kind=dp), allocatable :: test_beta_res(:), test_beta_model(:), test_priors(:)
   real(kind=dp), allocatable :: states_x_states_original_copy(:,:)
   real(kind=dp), allocatable :: temp3d(:,:,:), temp3d_2(:,:,:)

   character(len=:), allocatable :: base_trial_name, file_path, sst_file
   character(len=50) :: beta_res_char,beta_model_char,prior_char
   character(len=4)  :: year

   !Try syncing on un-noisy data
   if(.not.(allocated(reservoir%saved_state))) allocate(reservoir%saved_state(reservoir%n))
   reservoir%saved_state = 0
   un_noisy_sync = 24*364 !98!2160!700

   !From this point on reservoir%trainingdata and reservoir%imperfect_model have a temporal resolution of 1
   !hour instead of a resolution of model_parameters%timestep
   call get_prediction_data(reservoir,model_parameters,grid,model_parameters%traininglength-un_noisy_sync,un_noisy_sync)

   call synchronize(reservoir,reservoir%predictiondata,reservoir%saved_state,un_noisy_sync/model_parameters%timestep-1)

   if(reservoir%tisr_input_bool) then
      call get_full_tisr(reservoir,model_parameters,grid)
   endif 

   if(model_parameters%non_stationary_ocn_climo) then 
      call get_full_sst_climo(reservoir,model_parameters,grid)
   endif 

   if(.not.((model_parameters%slab_ocean_model_bool).and.(grid%bottom))) then
      deallocate(reservoir%predictiondata)
   endif 

   allocate(reservoir%local_model(reservoir%chunk_size_speedy))
   allocate(reservoir%outvec(reservoir%chunk_size_prediction))
   allocate(reservoir%feedback(reservoir%reservoir_numinputs))

   if(model_parameters%outvec_component_contribs) then
     allocate(reservoir%v_ml(reservoir%chunk_size_prediction))
     allocate(reservoir%v_p(reservoir%chunk_size_prediction))
   endif 

   if(model_parameters%slab_ocean_model_bool .and. mpi_res%is_root .and. grid%bottom .and. reservoir%assigned_region == 0) then
     print *, 'doing sea mask and default values'
     allocate(model_parameters%base_sst_grid(xgrid, ygrid))
     allocate(model_parameters%sea_mask(xgrid, ygrid))

     write(year,'(I4)') calendar%currentyear

     file_path = '/scratch/user/troyarcomano/ERA_5/'//year//'/'
     sst_file = file_path//'era_5_y'//year//'_sst_regridded_fixed_var_gcc.nc' !'_sst_regridded_mpi_fixed_var_gcc.nc'

     call read_netcdf_3d('sst',sst_file,temp3d)

     if(model_parameters%train_on_sst_anomalies) then
       call read_netcdf_3d('sst','/scratch/user/troyarcomano/ERA_5/regridded_era_sst_climatology1981_1999_gcc.nc',temp3d_2)
     endif 
 
     if(.not. model_parameters%train_on_sst_anomalies) then
     where(temp3d < 272.0_dp)
       temp3d = 272.0_dp
     end where
     endif 

     model_parameters%base_sst_grid = temp3d(:,:,1) !this is a major bug should be the current forecast start time 

     if(model_parameters%train_on_sst_anomalies) then
        model_parameters%base_sst_grid =  model_parameters%base_sst_grid - temp3d_2(:,:,1)
     endif 

     print *, 'shape(temp3d)',shape(temp3d)
     do i=1,xgrid
        do j=1, ygrid  
           if(all(temp3d(i,j,:) < 273.1)) then
             print *, 'land/permanent ice at',i,j,model_parameters%base_sst_grid(i,j)
             model_parameters%sea_mask(i,j) = 1.0
           else
             model_parameters%sea_mask(i,j) = 0.0
           endif
         enddo
     enddo 
     deallocate(temp3d)
   endif
end subroutine 

subroutine get_full_tisr(reservoir,model_parameters,grid)
   use mpires, only        : mpi_res
   use mod_io, only        : read_3d_file_parallel 
   use mod_utilities, only : standardize_data_given_pars3d

   type(reservoir_type), intent(inout)        :: reservoir
   type(grid_type), intent(inout)             :: grid
   type(model_parameters_type), intent(inout) :: model_parameters

   character(len=:), allocatable :: file_path
   character(len=:), allocatable :: tisr_file
 
   file_path = '/scratch/user/troyarcomano/ERA_5/2007/' !2012/'
   tisr_file = file_path//'toa_incident_solar_radiation_2007_regridded_classic4.nc'

   call read_3d_file_parallel(tisr_file,'tisr',mpi_res,grid,reservoir%full_tisr,1,1)
   print *, 'tisr_mean_std_idx,grid%mean(grid%tisr_mean_std_idx),grid%std(grid%tisr_mean_std_idx)',grid%tisr_mean_std_idx,grid%mean(grid%tisr_mean_std_idx),grid%std(grid%tisr_mean_std_idx)
   call standardize_data_given_pars3d(reservoir%full_tisr,grid%mean(grid%tisr_mean_std_idx),grid%std(grid%tisr_mean_std_idx))

end subroutine  

subroutine get_full_sst_climo(reservoir,model_parameters,grid)
   use mpires, only        : mpi_res
   use mod_io, only        : read_3d_file_parallel
   use mod_utilities, only : standardize_data_given_pars3d

   type(reservoir_type), intent(inout)        :: reservoir
   type(grid_type), intent(inout)             :: grid
   type(model_parameters_type), intent(inout) :: model_parameters

   character(len=:), allocatable :: file_path
   character(len=:), allocatable :: tisr_file

   file_path = '/scratch/user/troyarcomano/ERA_5/'
   tisr_file = file_path//'sst_regridded_climo1981_2021_gcc.nc'

   print *, 'sst_file',tisr_file
   call read_3d_file_parallel(tisr_file,'sst',mpi_res,grid,reservoir%full_sst,1,1)

   print *, 'tisr_mean_std_idx,grid%mean(grid%tisr_mean_std_idx),grid%std(grid%tisr_mean_std_idx)',grid%sst_mean_std_idx,grid%mean(grid%sst_mean_std_idx),grid%std(grid%sst_mean_std_idx)
 
   where(reservoir%full_sst < 272.0_dp)
       reservoir%full_sst = 272.0_dp
   end where
   if(reservoir%sst_bool_input) then
      call standardize_data_given_pars3d(reservoir%full_sst,grid%mean(grid%sst_mean_std_idx),grid%std(grid%sst_mean_std_idx))
   endif 

end subroutine

subroutine start_prediction(reservoir,model_parameters,grid,prediction_number)
   type(reservoir_type), intent(inout)        :: reservoir
   type(grid_type), intent(inout)             :: grid
   type(model_parameters_type), intent(inout) :: model_parameters

   integer, intent(in)                        :: prediction_number

   model_parameters%current_trial_number = prediction_number

   call get_prediction_data(reservoir,model_parameters,grid,model_parameters%traininglength+model_parameters%prediction_markers(prediction_number),model_parameters%synclength+100)
  
   call synchronize_print(reservoir,grid,reservoir%predictiondata(:,1:model_parameters%synclength/model_parameters%timestep-1),reservoir%saved_state,model_parameters%synclength/model_parameters%timestep-1)
  
   if(reservoir%assigned_region == 36)  print *, 'reservoir%predictiondata(:,model_parameters%synclength/model_parameters%timestep)',reservoir%predictiondata(:,model_parameters%synclength/model_parameters%timestep)
   if(reservoir%assigned_region == 36 .and. .not. model_parameters%ml_only)  print *, 'reservoir%imperfect_model_states(:,model_parameters%synclength/model_parameters%timestep-1)',reservoir%imperfect_model_states(:,model_parameters%synclength/model_parameters%timestep-1)

   reservoir%feedback = reservoir%predictiondata(:,model_parameters%synclength/model_parameters%timestep)

   if(.not. model_parameters%ml_only) then
     reservoir%local_model = reservoir%imperfect_model_states(:,model_parameters%synclength/model_parameters%timestep+1)
   endif 
end subroutine 

subroutine reservoir_layer_chunking_ml(reservoir,model_parameters,grid,trainingdata)
   use mpires
   use mod_utilities, only : gaussian_noise_1d_function, gaussian_noise_1d_function_precip

   type(reservoir_type), intent(inout)      :: reservoir
   type(model_parameters_type) , intent(in) :: model_parameters
   type(grid_type) , intent(in)            :: grid

   real(kind=dp), intent(in) :: trainingdata(:,:)

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
          !if(reservoir%assigned_region == 956) 
        else
          temp = matmul(reservoir%win,gaussian_noise_1d_function(trainingdata(:,model_parameters%discardlength/model_parameters%timestep+i),reservoir%noisemag))
        endif

        x_ = tanh(y+temp)
        x = (1-reservoir%leakage)*x + reservoir%leakage*x_

        reservoir%states(:,reservoir%batch_size) = x

        reservoir%saved_state = reservoir%states(:,reservoir%batch_size)

        reservoir%states(2:reservoir%n:2,:) = reservoir%states(2:reservoir%n:2,:)**2

        call chunking_matmul_ml(reservoir,model_parameters,grid,batch_number,trainingdata)

      elseif (mod(i,reservoir%batch_size).eq.0) then
        print *,'new state',i, 'region',reservoir%assigned_region

        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,reservoir%states(:,reservoir%batch_size),beta,y)

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

subroutine reservoir_layer_chunking_hybrid(reservoir,model_parameters,grid,trainingdata,imperfect_model)
   use mpires
   use mod_utilities, only : gaussian_noise_1d_function,gaussian_noise_1d_function_precip
   use mod_io, only : write_netcdf_2d_non_met_data_timeseries


   type(reservoir_type), intent(inout)      :: reservoir
   type(model_parameters_type) , intent(in) :: model_parameters
   type(grid_type) , intent(in)            :: grid

   real(kind=dp), intent(in) :: trainingdata(:,:)
   real(kind=dp), intent(in) :: imperfect_model(:,:)

   integer :: i,info
   integer :: training_length, batch_number

   real(kind=dp), allocatable :: temp(:), x(:), x_(:), y(:), temp2(:)
   real(kind=dp), parameter   :: alpha=1.0,beta=0.0
   real(kind=dp), allocatable :: gaussian_noise

   allocate(temp(reservoir%n),x(reservoir%n),x_(reservoir%n),y(reservoir%n))

   if(reservoir%assigned_region == 1094) allocate(temp2(reservoir%reservoir_numinputs))

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

        call chunking_matmul(reservoir,model_parameters,grid,batch_number,trainingdata,imperfect_model)

      elseif (mod(i,reservoir%batch_size).eq.0) then
        print *,'new state',i, 'region',reservoir%assigned_region

        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,reservoir%saved_state,beta,y) !reservoir%states(:,reservoir%batch_size),beta,y)

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
      if(reservoir%assigned_region == 1094) temp2 = gaussian_noise_1d_function(trainingdata(:,model_parameters%discardlength/model_parameters%timestep+i),reservoir%noisemag)
      if(reservoir%assigned_region == 1094) call write_netcdf_2d_non_met_data_timeseries(temp2,'input','inputvec_noisy.nc',i,'unitless','x','t')
      if(reservoir%assigned_region == 1094) call write_netcdf_2d_non_met_data_timeseries(trainingdata(:,model_parameters%discardlength/model_parameters%timestep+i),'input','inputvec.nc',i,'unitless','x','t')
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

    if(reservoir%assigned_region == 954)  print *, 'a_trans(1:20,1:20)',a_trans(1:20,1:20)
    if(reservoir%assigned_region == 954)  print *, 'b_trans(1:20,1:20)',b_trans(1:20,1:20)

    call mldivide(a_trans,b_trans)
    reservoir%wout = transpose(b_trans)

    if(reservoir%assigned_region == 954)  print *, 'worker',reservoir%assigned_region,'wout(1,1:20)',reservoir%wout(1,1:20)

    deallocate(a_trans)
    deallocate(b_trans)

    write(level_char,'(i0.2)') grid%level_index
    if(reservoir%assigned_region == 954) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_954_level_'//level_char//'wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')
    if(reservoir%assigned_region == 217) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_217_level_'//level_char//'wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')
    if(reservoir%assigned_region == 218) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_218_level_'//level_char//'wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')

    !call write_trained_res(reservoir,model_parameters,grid)

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
       
      do i=1, reservoir%chunk_size_speedy!reservoir%chunk_size_prediction
            prior(i,i) =  reservoir%prior_val*reservoir%beta_model**2.0_dp
      enddo 

    endif 

    !Do regularization

    !If we are doing a prior we need beta^2 
    if(model_parameters%using_prior) then 
      do i=1, reservoir%n+reservoir%chunk_size_speedy!prediction
         if(i <= reservoir%chunk_size_speedy) then!_prediction) then
            reservoir%states_x_states_aug(i,i) = reservoir%states_x_states_aug(i,i) + reservoir%beta_model**2.0_dp
         else
            reservoir%states_x_states_aug(i,i) = reservoir%states_x_states_aug(i,i) + reservoir%beta_res**2.0_dp
         endif
      enddo
    else
      do i=1, reservoir%n+reservoir%chunk_size_speedy!prediction
         if(i <= reservoir%chunk_size_speedy) then!prediction) then 
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

    if(reservoir%assigned_region == 954)  print *, 'a_trans(1:20,1:20)',a_trans(1:20,1:20)
    if(reservoir%assigned_region == 954)  print *, 'b_trans(1:20,1:20)',b_trans(1:20,1:20)
    !if(any(IEEE_IS_NAN(reservoir%states_x_states_aug))) print *, 'reservoir%states_x_states_aug nan', reservoir%assigned_region
    !if(any(IEEE_IS_NAN(reservoir%states_x_trainingdata_aug))) print *, 'reservoir%states_x_states_aug nan', reservoir%assigned_region
    !if(any(IEEE_IS_NAN(a_trans))) print *, 'a_trans has nan',reservoir%assigned_region
    !if(any(IEEE_IS_NAN(b_trans))) print *, 'b_trans has nan',reservoir%assigned_region
    
    !If we are trying a prior then we need to add it to b_trans
    if(model_parameters%using_prior) then
      b_trans = b_trans + transpose(prior)
    endif  

    call mldivide(a_trans,b_trans)
    reservoir%wout = transpose(b_trans)

    !if(any(IEEE_IS_NAN(reservoir%wout))) print *, 'wout has nan', reservoir%assigned_region
    !if(IEEE_IS_NAN(reservoir%wout(1,1))) print *, 'wout element 1 has nan', reservoir%assigned_region


    if(reservoir%assigned_region == 954)  print *, 'worker',reservoir%assigned_region,'wout(1,1:20)',reservoir%wout(1,1:20)

    deallocate(a_trans)
    deallocate(b_trans)

    write(level_char,'(i0.2)') grid%level_index
    !if(reservoir%assigned_region == 954) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_954_level_'//level_char//'wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')
    !if(reservoir%assigned_region == 217) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_217_level_'//level_char//'wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')
    !if(reservoir%assigned_region == 218) call write_netcdf_2d_non_met_data(reservoir%wout,'wout','region_218_level_'//level_char//'wout_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')

    if(write_training_weights) then 
       call write_trained_res(reservoir,model_parameters,grid)
    endif 

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
 
    call synchronize(reservoir,reservoir%predictiondata(:,1:model_parameters%synclength/model_parameters%timestep),x,model_parameters%synclength/model_parameters%timestep)

    call predict(reservoir,model_parameters,grid,x,imperfect_model_in)    
end subroutine 

subroutine synchronize(reservoir,input,x,length)
    type(reservoir_type), intent(inout) :: reservoir
    
    real(kind=dp), intent(in)     :: input(:,:)
    real(kind=dp), intent(inout)  :: x(:)

    integer, intent(in)           :: length

    real(kind=dp), allocatable    :: y(:), temp(:),  x_(:)
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
       x = (1.0_dp-reservoir%leakage)*x + reservoir%leakage*x_
    enddo 
    
    return 
end subroutine  

subroutine synchronize_print(reservoir,grid,input,x,length)
    type(reservoir_type), intent(inout) :: reservoir
    type(grid_type), intent(inout)      :: grid

    real(kind=dp), intent(in)     :: input(:,:)
    real(kind=dp), intent(inout)  :: x(:)

    integer, intent(in)           :: length

    real(kind=dp), allocatable    :: y(:), temp(:),  x_(:)
    real(kind=dp), parameter      :: alpha=1.0,beta=0.0

    integer :: info,i

    allocate(y(reservoir%n))
    allocate(temp(reservoir%n))
    allocate(x_(reservoir%n))

    y=0

    print *, 'shape(input) sync',reservoir%assigned_region,shape(input)
    print *, 'length sync',reservoir%assigned_region,length
    do i=1, length
       info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,x,beta,y)
       if(reservoir%assigned_region == 36)  print *, 'i',i
       temp = matmul(reservoir%win,input(:,i))
 
       x_ = tanh(y+temp)
       x = (1-reservoir%leakage)*x + reservoir%leakage*x_

    enddo

    return
end subroutine

subroutine predict(reservoir,model_parameters,grid,x,local_model_in)
    use mpires, only : predictionmpicontroller
    use resdomain, only : unstandardize_state_vec_res
    use mod_utilities, only : e_constant

    type(reservoir_type), intent(inout)     :: reservoir
    type(model_parameters_type), intent(in) :: model_parameters
    type(grid_type), intent(in)             :: grid

    real(kind=dp), intent(inout) :: x(:)
    real(kind=dp), intent(inout) :: local_model_in(:)

    real(kind=dp), allocatable :: y(:), temp(:), x_(:)
    real(kind=dp), allocatable :: local_model_temp(:)
    real(kind=dp), allocatable :: x_temp(:),x_augment(:)
    real(kind=dp), allocatable :: v_ml(:), v_p(:), v_h(:)

    real(kind=dp), parameter :: alpha=1.0,beta=0.0

    integer :: info,i,j

    allocate(y(reservoir%n),temp(reservoir%n),x_(reservoir%n))
    allocate(x_augment(reservoir%n+reservoir%chunk_size_speedy))!reservoir%chunk_size_prediction))

    y = 0

    info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,x,beta,y)
    temp = matmul(reservoir%win,reservoir%feedback)

    x_ = tanh(y + temp)
    x = (1-reservoir%leakage)*x + reservoir%leakage*x_

    x_temp = x
    x_temp(2:reservoir%n:2) = x_temp(2:reservoir%n:2)**2

    x_augment(1:reservoir%chunk_size_speedy) = reservoir%local_model !prediction) = reservoir%local_model
    x_augment(reservoir%chunk_size_speedy+1:reservoir%chunk_size_speedy+reservoir%n) = x_temp !prediction+1:reservoir%chunk_size_prediction+reservoir%n) = x_temp

    reservoir%outvec = matmul(reservoir%wout,x_augment)

    if(model_parameters%outvec_component_contribs) then
       reservoir%v_p = matmul(reservoir%wout(:,1:reservoir%chunk_size_speedy),reservoir%local_model)
       reservoir%v_ml = matmul(reservoir%wout(:,reservoir%chunk_size_speedy+1:reservoir%chunk_size_speedy+reservoir%n),x_temp)
    endif 
      
    if(reservoir%assigned_region == 954 .and. model_parameters%outvec_component_contribs) then
      print *,'reservoir%v_p',reservoir%v_p
      print *,'reservoir%v_ml',reservoir%v_ml
      print *,'reservoir%v_ml + reservoir%v_p',reservoir%v_ml + reservoir%v_p
      print *,'reservoir%outvec',reservoir%outvec
      print *,'diff',reservoir%outvec - (reservoir%v_ml + reservoir%v_p)
    endif 

    call unstandardize_state_vec_res(reservoir,grid,reservoir%outvec)

    !if(model_parameters%outvec_component_contribs) then 
    !  call unstandardize_state_vec_res(reservoir,grid,reservoir%v_p)
    !  call unstandardize_state_vec_res(reservoir,grid,reservoir%v_ml)
    !endif 

    if((reservoir%assigned_region == 954).and.(mod(i,24) == 0)) then
      print *, '*******' 
      print *, 'feedback',reservoir%feedback
      print *, '*******'
      print *, 'local_model',reservoir%local_model
      print *, '*******'
      print *, 'outvec atmo',reservoir%outvec
      print *, '*******'
      !print *, 'precip',model_parameters%precip_epsilon * (e_constant**reservoir%outvec(reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk+1:reservoir%local_predictvars*grid%resxchunk*grid%resychunk*reservoir%local_heightlevels_res+grid%resxchunk*grid%resychunk*2) - 1)
      !print *, 'feedback',reservoir%feedback
    endif
end subroutine

subroutine predict_ml(reservoir,model_parameters,grid,x)
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
    allocate(x_augment(reservoir%n+reservoir%chunk_size_speedy))!reservoir%chunk_size_prediction))

    y = 0

    info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,x,beta,y)
    temp = matmul(reservoir%win,reservoir%feedback)

    x_ = tanh(y + temp)
    x = (1-reservoir%leakage)*x + reservoir%leakage*x_

    x_temp = x
    x_temp(2:reservoir%n:2) = x_temp(2:reservoir%n:2)**2

    x_augment(1:reservoir%chunk_size_speedy+reservoir%n) = x_temp !prediction+1:reservoir%chunk_size_prediction+reservoir%n) = x_temp

    reservoir%outvec = matmul(reservoir%wout,x_augment)

    call unstandardize_state_vec_res(reservoir,grid,reservoir%outvec)

    if((reservoir%assigned_region == 954).and.(mod(i,1) == 0)) then
      print *, '*******'
      print *, 'feedback',reservoir%feedback
      print *, '*******'
      print *, 'outvec atmo',reservoir%outvec
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

   num_of_batches = 20!10*6!2!10!*5!0!10*6 !20*6
   approx_batch_size = (model_parameters%traininglength - model_parameters%discardlength)/(num_of_batches*model_parameters%timestep)

   !routine to get the closest reservoir%batch_size to num_of_batches that
   !divides into reservoir%traininglength
   call find_closest_divisor(approx_batch_size,(model_parameters%traininglength - model_parameters%discardlength)/model_parameters%timestep,reservoir%batch_size)
   print *, 'num_of_batches,approx_batch_size,reservoir%traininglength,reservoir%batch_size',num_of_batches,approx_batch_size,model_parameters%traininglength-model_parameters%discardlength,reservoir%batch_size

   !Should be reservoir%n+ reservoir%chunk_size
   allocate(reservoir%states_x_trainingdata_aug(reservoir%chunk_size_prediction,reservoir%n+reservoir%chunk_size_speedy))!prediction))
   allocate(reservoir%states_x_states_aug(reservoir%n+reservoir%chunk_size_speedy,reservoir%n+reservoir%chunk_size_speedy))!prediction,reservoir%n+reservoir%chunk_size_prediction))
   allocate(reservoir%states(reservoir%n,reservoir%batch_size))
   allocate(reservoir%augmented_states(reservoir%n+reservoir%chunk_size_speedy,reservoir%batch_size))!prediction,reservoir%batch_size))
   allocate(reservoir%saved_state(reservoir%n))
   allocate(reservoir%current_state(reservoir%n))

   reservoir%states_x_trainingdata_aug = 0.0_dp
   reservoir%states_x_states_aug = 0.0_dp
   reservoir%states = 0.0_dp
   reservoir%current_state = 0.0_dp
   reservoir%augmented_states = 0.0_dp

end subroutine

subroutine chunking_matmul_ml(reservoir,model_parameters,grid,batch_number,trainingdata)
   use mod_utilities, only : gaussian_noise
   use resdomain, only : tile_full_input_to_target_data
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

   reservoir%augmented_states(reservoir%chunk_size_speedy+1:reservoir%n+reservoir%chunk_size_speedy,:) = reservoir%states

   print *, 'grid%predict_end',grid%predict_end,model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*m+1,batch_number*m+model_parameters%discardlength/model_parameters%timestep

   call tile_full_input_to_target_data(reservoir,grid,trainingdata(1:grid%predict_end,model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*m+1:batch_number*m+model_parameters%discardlength/model_parameters%timestep),targetdata)

   if(reservoir%assigned_region == 954)  print *, 'shape(trainingdata)',shape(trainingdata)
   if(reservoir%assigned_region == 954)  print *, 'shape(targetdata)',shape(targetdata)

   allocate(temp(reservoir%chunk_size_prediction,n))
   temp = 0.0_dp

   temp = matmul(targetdata,transpose(reservoir%augmented_states))
   !TODO make this matmul DGEMM

   reservoir%states_x_trainingdata_aug = reservoir%states_x_trainingdata_aug + temp
   deallocate(temp)
   deallocate(targetdata)

   allocate(temp(n,n))
   call DGEMM('N','N',n,n,m,alpha,reservoir%augmented_states,n,transpose(reservoir%augmented_states),m,beta,temp,n)
   !call
   !DGEMM('N','T',n,n,m,alpha,reservoir%augmented_states,n,reservoir%augmented_states,m,beta,temp,n)
   reservoir%states_x_states_aug = reservoir%states_x_states_aug + temp
   deallocate(temp)

   return
end subroutine


subroutine chunking_matmul(reservoir,model_parameters,grid,batch_number,trainingdata,imperfect_model)
   use mod_utilities, only : gaussian_noise
   use resdomain, only : tile_full_input_to_target_data
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

   !reservoir%augmented_states(1:reservoir%chunk_size_prediction,:) = imperfect_model(:,model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*m+1:batch_number*m+model_parameters%discardlength/model_parameters%timestep)
   reservoir%augmented_states(1:reservoir%chunk_size_speedy,:) = imperfect_model(:,model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*m+1:batch_number*m+model_parameters%discardlength/model_parameters%timestep)
   !if(any(IEEE_IS_NAN(imperfect_model))) print *, 'imperfect_model has nan',reservoir%assigned_region,batch_number
   
   !reservoir%augmented_states(reservoir%chunk_size_prediction+1:reservoir%n+reservoir%chunk_size_prediction,:) = reservoir%states
   reservoir%augmented_states(reservoir%chunk_size_speedy+1:reservoir%n+reservoir%chunk_size_speedy,:) = reservoir%states
   !if(any(IEEE_IS_NAN(reservoir%states))) print *, 'reservoir%states has nan',reservoir%assigned_region,batch_number

   print *, 'grid%predict_end',grid%predict_end,model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*m+1,batch_number*m+model_parameters%discardlength/model_parameters%timestep

   call tile_full_input_to_target_data(reservoir,grid,trainingdata(1:grid%predict_end,model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*m+1:batch_number*m+model_parameters%discardlength/model_parameters%timestep),targetdata)
   if(any(IEEE_IS_NAN(targetdata))) print *, 'targetdata has nan',reservoir%assigned_region,batch_number
   !if(reservoir%assigned_region == 954)  print *, 'model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*m+1:batch_number*m+model_parameters%discardlength/model_parameters%timestep',model_parameters%discardlength/model_parameters%timestep+(batch_number-1)*m+1,batch_number*m+model_parameters%discardlength/model_parameters%timestep
   !if(reservoir%assigned_region == 954)  print *, 'model_parameters%discardlength',model_parameters%discardlength,'batch_number',batch_number,'m',m
   if(reservoir%assigned_region == 954)  print *, 'shape(trainingdata)',shape(trainingdata)
   if(reservoir%assigned_region == 954)  print *, 'shape(targetdata)',shape(targetdata)

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
   reservoir%states_x_states_aug = reservoir%states_x_states_aug + temp
   deallocate(temp)
   
   return 
end subroutine  

subroutine write_trained_res(reservoir,model_parameters,grid)
  use mod_io, only : write_netcdf_2d_non_met_data, write_netcdf_1d_non_met_data_int, write_netcdf_1d_non_met_data_real, write_netcdf_1d_non_met_data_real_shared_dim

  type(reservoir_type), intent(in) :: reservoir
  type(model_parameters_type), intent(in) :: model_parameters
  type(grid_type), intent(in)             :: grid

  character(len=:), allocatable :: file_path
  character(len=4) :: worker_char
  character(len=1) :: height_char
  character(len=2) :: ens_member

  integer :: i 

  file_path = '/scratch/user/dpp94/ML_SPEEDY_WEIGHTS/'

  write(worker_char,'(i0.4)') reservoir%assigned_region
  write(height_char,'(i0.1)') grid%level_index

  if((reservoir%assigned_region == 0).and.(grid%level_index == 1)) then
    call write_controller_file(model_parameters,reservoir) 
  endif 

  print *, 'write_trained_res',reservoir%assigned_region
  call write_netcdf_2d_non_met_data(reservoir%win,'win',file_path//'worker_'//worker_char//'_level_'//height_char//'_'//trim(model_parameters%trial_name)//'.nc','unitless','win_x','win_y')
  call write_netcdf_2d_non_met_data(reservoir%wout,'wout',file_path//'worker_'//worker_char//'_level_'//height_char//'_'//trim(model_parameters%trial_name)//'.nc','unitless','wout_x','wout_y')

  call write_netcdf_1d_non_met_data_int(reservoir%rows,'rows',file_path//'worker_'//worker_char//'_level_'//height_char//'_'//trim(model_parameters%trial_name)//'.nc','unitless','rows_x')
  call write_netcdf_1d_non_met_data_int(reservoir%cols,'cols',file_path//'worker_'//worker_char//'_level_'//height_char//'_'//trim(model_parameters%trial_name)//'.nc','unitless','cols_x')

  call write_netcdf_1d_non_met_data_real(reservoir%vals,'vals',file_path//'worker_'//worker_char//'_level_'//height_char//'_'//trim(model_parameters%trial_name)//'.nc','unitless','vals_x')

  call write_netcdf_1d_non_met_data_real(grid%mean,'mean',file_path//'worker_'//worker_char//'_level_'//height_char//'_'//trim(model_parameters%trial_name)//'.nc','unitless','mean_x')
  call write_netcdf_1d_non_met_data_real(grid%std,'std',file_path//'worker_'//worker_char//'_level_'//height_char//'_'//trim(model_parameters%trial_name)//'.nc','unitless','std_x')

end subroutine

subroutine write_controller_file(model_parameters,reservoir)
   type(model_parameters_type), intent(in) :: model_parameters
   type(reservoir_type), intent(in)        :: reservoir

   character(len=:), allocatable :: file_path

   file_path = '/scratch/user/dpp94/ML_SPEEDY_WEIGHTS/'//trim(model_parameters%trial_name)//'_controller_file.txt'
 
   open (10, file=file_path, status='unknown')

   ! write to file
   write(10,*)"-----------------------------------------------------------"
   write(10,*) "num_hor_regions:",model_parameters%number_of_regions
   write(10,*)"ml_only:",model_parameters%ml_only
   write(10,*)"num_vert_levels:",model_parameters%num_vert_levels
   write(10,*)"atmo_timestep:",model_parameters%timestep
   write(10,*)"ocean_timestep:",model_parameters%timestep_slab
   write(10,*)"ocean_model_bool:",model_parameters%slab_ocean_model_bool 
   write(10,*)"train_on_sst_anomalies:",model_parameters%train_on_sst_anomalies
   write(10,*)"precip_bool:",model_parameters%precip_bool
   write(10,*)"precip_epsilon:",model_parameters%precip_epsilon
   write(10,*)"full_predictvars:",model_parameters%full_predictvars
   write(10,*)"full_heightlevels:",model_parameters%full_heightlevels
   write(10,*)"num_vert_levels:",model_parameters%num_vert_levels
   write(10,*)"vert_loc_overlap:",model_parameters%vert_loc_overlap
   write(10,*)"overlap:",model_parameters%overlap
   write(10,*)"regional_vary:",model_parameters%regional_vary
   write(10,*)"using_prior:",model_parameters%using_prior
   write(10,*)"reservoir_nodes:",reservoir%m
   write(10,*)"deg:",reservoir%deg
   write(10,*)"radius:",reservoir%radius
   write(10,*)"beta_res:",reservoir%beta_res
   write(10,*)"beta_model:",reservoir%beta_model
   write(10,*)"sigma:",reservoir%sigma
   write(10,*)"leakage:",reservoir%leakage
   write(10,*)"prior_val:",reservoir%prior_val
   write(10,*)"-----------------------------------------------------------"

   ! close file
   close(10) 

end subroutine 

subroutine trained_reservoir_prediction(reservoir,model_parameters,grid)
  use mod_linalg, only : mklsparse 
  use mod_io, only : read_trained_res

  type(reservoir_type), intent(inout)     :: reservoir
  type(model_parameters_type), intent(in) :: model_parameters
  type(grid_type), intent(inout)          :: grid

  integer :: mean_std_length

  if(grid%bottom) then
    reservoir%logp_bool = .True.
    reservoir%tisr_input_bool = .True.
    grid%logp_bool = .True.

    reservoir%sst_bool = model_parameters%slab_ocean_model_bool
    reservoir%sst_climo_bool = .False. !.False.

    reservoir%precip_input_bool = model_parameters%precip_bool
    reservoir%precip_bool = model_parameters%precip_bool

    reservoir%m = 6000

  else
    reservoir%sst_climo_bool = .False.
    reservoir%logp_bool = .False.
    reservoir%tisr_input_bool = .True.
    reservoir%sst_bool = .False.
    reservoir%precip_input_bool = .False.
    reservoir%precip_bool = .False.
  endif
 
  reservoir%local_predictvars = model_parameters%full_predictvars
  reservoir%local_heightlevels_input = grid%inputzchunk

  reservoir%local_heightlevels_res = grid%reszchunk

  call read_trained_res(reservoir,model_parameters,grid)

  !Get number of height levels * vars + 2d variables
  mean_std_length = model_parameters%full_predictvars*grid%inputzchunk
  if(reservoir%logp_bool) then
    mean_std_length = mean_std_length + 1
    grid%logp_mean_std_idx = mean_std_length
  endif

  if(reservoir%tisr_input_bool) then
    mean_std_length = mean_std_length + 1
    grid%tisr_mean_std_idx = mean_std_length
  endif

  if(reservoir%precip_bool)  then
    mean_std_length = mean_std_length + 1
    grid%precip_mean_std_idx = mean_std_length
  endif

  if(reservoir%sst_bool) then
    mean_std_length = mean_std_length + 1
    grid%sst_mean_std_idx = mean_std_length
    print *, 'grid%std(grid%sst_mean_std_idx)',grid%std(grid%sst_mean_std_idx)
    if(grid%std(grid%sst_mean_std_idx) > 0.2) then
      reservoir%sst_bool_input = .True.
    else    
      reservoir%sst_bool_input = .False.  
    endif 
  endif

  call allocate_res_new(reservoir,grid,model_parameters)

  call mklsparse(reservoir)
   
  grid%logp_start = 0
  grid%logp_end = 0
  grid%sst_start = 0
  grid%sst_end = 0

  grid%atmo3d_start = 1
  grid%atmo3d_end = model_parameters%full_predictvars*grid%inputxchunk*grid%inputychunk*grid%inputzchunk

  grid%predict_start = 1
  grid%predict_end = grid%atmo3d_end

   if(reservoir%logp_bool) then
     grid%logp_start = grid%atmo3d_end + 1
     grid%logp_end = grid%atmo3d_end + reservoir%logp_size_input
     grid%predict_end = grid%logp_end
   endif

   if(reservoir%precip_bool) then
     grid%precip_start = grid%atmo3d_end + reservoir%logp_size_input + 1
     grid%precip_end = grid%precip_start + reservoir%precip_size_input - 1
     grid%predict_end = grid%precip_end
   endif

   if(reservoir%sst_bool_input) then
     grid%sst_start = grid%atmo3d_end + reservoir%logp_size_input + reservoir%precip_size_input + 1
     grid%sst_end =  grid%sst_start + reservoir%sst_size_input - 1
   endif

   if(reservoir%tisr_input_bool) then
     grid%tisr_start = grid%atmo3d_end + reservoir%logp_size_input + reservoir%precip_size_input + reservoir%sst_size_input + 1
     grid%tisr_end = grid%tisr_start + reservoir%tisr_size_input - 1
   endif
end subroutine

end module 
