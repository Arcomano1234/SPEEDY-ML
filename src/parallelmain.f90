program main
  !use mpi_f08
  use mpi
  use, intrinsic :: ieee_arithmetic

  use mpires, only : mpi_res, startmpi, distribute_prediction_marker, killmpi, predictionmpicontroller, sendrecievegrid, send_outvec_ml_contrib, send_outvec_speedy_contrib
  use mod_reservoir, only : initialize_model_parameters, allocate_res_new, train_reservoir, start_prediction, initialize_prediction, predict, trained_reservoir_prediction, predict_ml
  use mod_slab_ocean_reservoir, only : initialize_slab_ocean_model, train_slab_ocean_model, get_training_data_from_atmo, initialize_prediction_slab, start_prediction_slab, predict_slab, predict_slab_ml, trained_ocean_reservoir_prediction
  use speedy_res_interface, only : startspeedy
  use resdomain, only : processor_decomposition, initializedomain, set_reservoir_by_region
  use mod_utilities, only : main_type, init_random_seed, dp, gaussian_noise, standardize_data_given_pars4d, standardize_data_given_pars3d, standardize_data, init_random_marker
  use mod_calendar
  !use mod_unit_tests, only : test_linalg, test_res_domain #TODO not working yet

  implicit none 

  integer :: standardizing_vars, i, j, k , t, prediction_num

  logical :: runspeedy = .False.
  logical :: trained_model = .True.
  logical :: slab_model
 
  real(kind=dp), allocatable :: test_state(:), test_feedback(:)

  type(main_type) :: res

  !Fortran has command line augs TODO 

  !Starts the MPI stuff and initializes mpi_res
  call startmpi()
 
  !mpi_res%numprocs = 1152 

  !Makes the object called res and declares all of the main parameters 
  call initialize_model_parameters(res%model_parameters,mpi_res%proc_num,mpi_res%numprocs)
  
  !Do domain decomposition based off processors and do vertical localization of
  !reservoir
  call processor_decomposition(res%model_parameters)

  !Need this for each worker gets a new random seed
  call init_random_marker(33)

  !Allocate atmo3d reservoirs and any special ones 
  allocate(res%reservoir(res%model_parameters%num_of_regions_on_proc,res%model_parameters%num_vert_levels))
  allocate(res%grid(res%model_parameters%num_of_regions_on_proc,res%model_parameters%num_vert_levels))

  if(res%model_parameters%slab_ocean_model_bool) then
    res%model_parameters%special_reservoirs = .True.
    res%model_parameters%num_special_reservoirs = 1
  endif 

  !NOTE one day may make precip its own special reservoir 
  !if(res%model_parameters%precip_bool) then
  !  res%model_parameters%num_special_reservoirs = res%model_parameters%num_special_reservoirs + 1
  !endif

  if(res%model_parameters%special_reservoirs) then
    allocate(res%reservoir_special(res%model_parameters%num_of_regions_on_proc,res%model_parameters%num_special_reservoirs))
    allocate(res%grid_special(res%model_parameters%num_of_regions_on_proc,res%model_parameters%num_special_reservoirs))
  endif  


  !---This is for debugging----!
  !You can run the code with a small number of processors and look at a few
  !regions of the globe 
  !if(res%model_parameters%irank == 4) res%model_parameters%region_indices(1) = 954
  !if(res%model_parameters%irank == 2) res%model_parameters%region_indices(1) = 552
  !if(res%model_parameters%irank == 3)  res%model_parameters%region_indices(1) = 36

  if(.not.(trained_model)) then
  !-----------Main Training Loop-------------!
  !Main training loop we loop through each region and each level at that region
  !for each processor. Every processor has its specific regions
  !
  !First we initialize the domain by calling initializedomain which populates
  !grid(i,j) with the necessary information about the grid
  !
  !Second we initialize the derived type called reservoir by call
  !allocate_res_new and makes reservoir(i,j)
  !
  !Finally we train the reservoir in the subroutine train_reservoir

  !Loop 1: Loop over all sub domains (regions) on each processor
  do i=1,res%model_parameters%num_of_regions_on_proc

     !Loop 2: Loop over each vertical level for a particular sub domain 
     do j=1,res%model_parameters%num_vert_levels

        call initializedomain(res%model_parameters%number_of_regions,res%model_parameters%region_indices(i), & 
                              res%model_parameters%overlap,res%model_parameters%num_vert_levels,j,res%model_parameters%vert_loc_overlap, &
                              res%grid(i,j))


        res%grid(i,j)%level_index = j

        print *,'region,level,input_zstart,nput_zend,inputzchunk,b,t',res%model_parameters%region_indices(i),j,res%grid(i,j)%input_zstart,res%grid(i,j)%input_zend,res%grid(i,j)%inputzchunk,res%grid(i,j)%bottom,res%grid(i,j)%top
        print *,'region,level,resxchunk,resychunk,reszchunk',res%model_parameters%region_indices(i),j,res%grid(i,j)%resxchunk,res%grid(i,j)%resychunk,res%grid(i,j)%reszchunk
        print *,'region,level,res_zstart,res_zend',res%model_parameters%region_indices(i),j,res%grid(i,j)%res_zstart,res%grid(i,j)%res_zend
        print *,'region,level,input_ystart,input_yend,inputychunk',res%model_parameters%region_indices(i),j,res%grid(i,j)%input_ystart,res%grid(i,j)%input_yend,res%grid(i,j)%inputychunk
        print *, 'region,level,tdata_yend',res%grid(i,j)%tdata_yend

        res%reservoir(i,j)%assigned_region = res%model_parameters%region_indices(i)

        !print *, 'res%reservoir(i,j):reservoir_numinputs,chunk_size_prediction',res%reservoir(i,j)%reservoir_numinputs, res%reservoir(i,j)%chunk_size_prediction
        call train_reservoir(res%reservoir(i,j),res%grid(i,j),res%model_parameters)
     enddo 
    
     !Lets do all of the training for these special reservoirs for each sub
     !region
     if(res%model_parameters%slab_ocean_model_bool) then
       call initializedomain(res%model_parameters%number_of_regions,res%model_parameters%region_indices(i), &
                             res%model_parameters%overlap,res%model_parameters%num_vert_levels,j-1,res%model_parameters%vert_loc_overlap, &
                             res%grid_special(i,1))


       res%grid_special(i,1)%level_index = j-1
     
       res%reservoir_special(i,1)%assigned_region = res%model_parameters%region_indices(i)

       print *, 'i,j',i,j,res%grid(i,j-1)%level_index,res%grid(i,j-1)%bottom,res%reservoir(i,j-1)%sst_bool_input
       print *, 'shape(res%reservoir(i,j-1)%trainingdata)',shape(res%reservoir(i,j-1)%trainingdata)

       call get_training_data_from_atmo(res%reservoir_special(i,1),res%model_parameters,res%grid_special(i,1),res%reservoir(i,j-1),res%grid(i,j-1))
       !We only want to train regions with oceans and or lakes
       if(res%reservoir_special(i,1)%sst_bool_prediction) then 
         call initialize_slab_ocean_model(res%reservoir_special(i,1),res%grid_special(i,1),res%model_parameters)    
         call train_slab_ocean_model(res%reservoir_special(i,1),res%grid_special(i,1),res%model_parameters) 
         
         test_feedback = res%reservoir_special(i,1)%trainingdata(:,res%model_parameters%traininglength)
         test_state = res%reservoir_special(i,1)%saved_state
         deallocate(res%reservoir_special(i,1)%trainingdata) 
       else 
         print *, 'i,j not training slab ocean',i,j
       endif 
     endif 

  enddo 
  endif 

  !If we already trained and are just reading in files then we go here 
  if(trained_model) then
    !Loop 1: Loop over all sub domains (regions) on each processor
    print *, 'res%model_parameters%num_of_regions_on_proc',res%model_parameters%num_of_regions_on_proc
    do i=1, res%model_parameters%num_of_regions_on_proc
       print *, 'i', i
       !Loop 2: Loop over each vertical level for a particular sub domain
        do j=1,res%model_parameters%num_vert_levels
           print *, 'j', j
           print *, 'doing initializedomain'
           call initializedomain(res%model_parameters%number_of_regions,res%model_parameters%region_indices(i), &
                                 res%model_parameters%overlap,res%model_parameters%num_vert_levels,j,res%model_parameters%vert_loc_overlap, &
                                 res%grid(i,j))


           res%reservoir(i,j)%assigned_region = res%model_parameters%region_indices(i)
           res%grid(i,j)%level_index = j

           print *, 'doing trained_reservoir_prediction'

           call initialize_calendar(calendar,1981,1,1,0)

           call trained_reservoir_prediction(res%reservoir(i,j),res%model_parameters,res%grid(i,j))
            
        enddo
  
        !Lets read in special reservoir 
        if(res%model_parameters%slab_ocean_model_bool) then 
          call initializedomain(res%model_parameters%number_of_regions,res%model_parameters%region_indices(i), &
                             res%model_parameters%overlap,res%model_parameters%num_vert_levels,j-1,res%model_parameters%vert_loc_overlap, &
                             res%grid_special(i,1))


          res%grid_special(i,1)%level_index = j-1

          res%reservoir_special(i,1)%assigned_region = res%model_parameters%region_indices(i)
 
          call trained_ocean_reservoir_prediction(res%reservoir_special(i,1),res%model_parameters,res%grid_special(i,1),res%reservoir(i,j-1),res%grid(i,j-1))
        endif 
     enddo
     print *, 'done reading trained model'
  endif  
    
  !Initialize Prediction 
  !Loop through all of the regions and vertical levels 
  do i=1,res%model_parameters%num_of_regions_on_proc
     do j=1,res%model_parameters%num_vert_levels
        print *,'initialize prediction region,level',res%reservoir(i,j)%assigned_region,res%grid(i,j)%level_index
        call initialize_prediction(res%reservoir(i,j),res%model_parameters,res%grid(i,j))  
     enddo 

     if(res%model_parameters%slab_ocean_model_bool) then
        !if(res%reservoir_special(i,1)%sst_bool_prediction) then
          print *,'ocean model initialize prediction region,i',res%reservoir_special(i,1)%assigned_region,i
          print *, 'shape(res%reservoir_special)',shape(res%reservoir_special)
          call initialize_prediction_slab(res%reservoir_special(i,1),res%model_parameters,res%grid_special(i,1),res%reservoir(i,j-1),res%grid(i,j-1))
        !endif
     endif 
  enddo 

  !Main prediction loop. 
  !Loop 1 through the user specified number of predictions
  !Loop 2 through time for a specific prediction
  !Loop 3/4 over the number of regions on the processor and all of the vertical
  !levels for a region
  do prediction_num=1, res%model_parameters%num_predictions
     do t=1, res%model_parameters%predictionlength/res%model_parameters%timestep
        if(t == 1) then 
          do i=1, res%model_parameters%num_of_regions_on_proc
             do j=1,res%model_parameters%num_vert_levels
                if(res%reservoir(i,j)%assigned_region == 954) print *, 'starting start_prediction region',res%model_parameters%region_indices(i),'prediction_num prediction_num',prediction_num
                call start_prediction(res%reservoir(i,j),res%model_parameters,res%grid(i,j),prediction_num)
           
                res%reservoir(i,j)%current_state = res%reservoir(i,j)%saved_state 
             enddo
             if(res%model_parameters%slab_ocean_model_bool) then
               !if(res%reservoir_special(i,1)%sst_bool_prediction) then 
                 call start_prediction_slab(res%reservoir_special(i,1),res%model_parameters,res%grid_special(i,1),res%reservoir(i,j-1),res%grid(i,j-1),prediction_num) 

                 if(res%reservoir_special(i,1)%sst_bool_prediction) then
                   res%reservoir_special(i,1)%current_state = res%reservoir_special(i,1)%saved_state
                 endif 
             endif 
          enddo 
        endif
        do i=1, res%model_parameters%num_of_regions_on_proc
           do j=1, res%model_parameters%num_vert_levels
              if(res%reservoir(i,j)%assigned_region == 954) print *, 'calling predict'
              if(res%model_parameters%ml_only) then
                call predict_ml(res%reservoir(i,j),res%model_parameters,res%grid(i,j),res%reservoir(i,j)%current_state)
                res%model_parameters%run_speedy = .True.
              else
                call predict(res%reservoir(i,j),res%model_parameters,res%grid(i,j),res%reservoir(i,j)%current_state,res%reservoir(i,j)%local_model)
              endif 
           enddo
           !print *, 'mod((t-1)*res%model_parameters%timestep,res%model_parameters%timestep_slab)',mod((t-1)*res%model_parameters%timestep,res%model_parameters%timestep_slab)
           if(res%model_parameters%slab_ocean_model_bool) then !NOTE TODO ML ocean doesnt get called until t = mod(res%model_parameters%timestep,res%model_parameters%timestep_slab) == 0
             if(mod((t)*res%model_parameters%timestep,res%model_parameters%timestep_slab) == 0 .and. res%reservoir_special(i,1)%sst_bool_prediction .and. .not. res%model_parameters%non_stationary_ocn_climo ) then
                if(res%reservoir_special(i,1)%assigned_region == 954) print *, 'calling predict slab'
                !TODO rolling_average_over_a_period(grid,period)
                !if( t > 28) then
                !  res%reservoir_special(i,1)%local_model = res%reservoir_special(i,1)%outvec 
                !endif 
                if(res%model_parameters%ml_only_ocean) then
                  call predict_slab_ml(res%reservoir_special(i,1),res%model_parameters,res%grid_special(i,1),res%reservoir_special(i,1)%current_state)
                else
                  call predict_slab(res%reservoir_special(i,1),res%model_parameters,res%grid_special(i,1),res%reservoir_special(i,1)%current_state,res%reservoir_special(i,1)%local_model)
                endif 
              endif 
           endif 
        enddo

  
        if(res%model_parameters%slab_ocean_model_bool) then !if(mod(t*res%model_parameters%timestep,res%model_parameters%timestep_slab) == 0) then
           slab_model = .True.
        else
           slab_model = .False.
        endif

        if(mpi_res%is_root) print *, 'sending data and writing predictions','prediction_num prediction_num',prediction_num,'time',t

        call sendrecievegrid(res,t,slab_model)
     
        if(res%model_parameters%outvec_component_contribs) then
           call send_outvec_ml_contrib(res,t)
           call send_outvec_speedy_contrib(res,t)
        endif 

        if(res%model_parameters%run_speedy .eqv. .False.) then
          exit
        endif
      enddo
  enddo 

  call MPI_Barrier(mpi_res%mpi_world, mpi_res%ierr)

  call mpi_finalize(mpi_res%ierr)

  if(res%model_parameters%irank == 0) then
     print *, 'program finished correctly'
  endif   

end program

