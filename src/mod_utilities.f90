module mod_utilities
  use iso_fortran_env
  use MKL_SPBLAS
  !use mpi_f08
  use mpi

  implicit none
  
  integer, parameter :: dp=selected_real_kind(14)  !Double precision (64 bit)
  integer, parameter :: sp = selected_real_kind(6, 37) !Single precison (32 bit)
  integer, parameter :: int_32 = int32

  real(kind=dp), parameter :: e_constant = 2.7182818284590452353602874_dp

  !every processor needs this and its static
  !Layers out the model grid
  integer :: speedygridnum=96*48
  integer(kind=int32) :: xgrid=96
  integer(kind=int32) :: ygrid=48
  integer(kind=int32) :: zgrid=8

  !SPEEDY lat and lons
  real(kind=dp) :: speedylat(48) = (/-87.159, -83.479, -79.777,  -76.070, -72.362, -68.652, &
                      -64.942, -61.232, -57.521, -53.810, -50.099, -46.389, -42.678, -38.967, &
                      -35.256, -31.545 , -27.833, -24.122, -20.411, -16.700, -12.989, -9.278, &
                      -5.567, -1.856, 1.856, 5.567, 9.278, 12.989, 16.700, 20.411, &
                      24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389, 50.099, &
                      53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, &
                      83.479, 87.159/)


  type grid_type
     !Grid derived type for a specific reservoir 

     !i j indices in 2d space of the data that the reservoir is fitted for the processor (wholegrid)
     integer :: res_xstart
     integer :: res_xend
     integer :: res_ystart
     integer :: res_yend
     integer :: res_zstart
     integer :: res_zend
     integer :: resxchunk
     integer :: resychunk
     integer :: reszchunk

     !Memory space 2d indices for the fitting data (not the wholegrid just local
     !grid)
     integer :: tdata_xstart
     integer :: tdata_xend
     integer :: tdata_ystart
     integer :: tdata_yend
     integer :: tdata_zstart
     integer :: tdata_zend
     
     !ij indices for the input data for the processor (wholegrid indices)
     integer :: input_xstart
     integer :: input_xend
     integer :: input_ystart
     integer :: input_yend
     integer :: input_zstart
     integer :: input_zend
     integer :: inputxchunk
     integer :: inputychunk
     integer :: inputzchunk

     !helpful logical flags
     logical :: pole
     logical :: periodicboundary
     logical :: top_vert_level
     logical :: bottom_vert_level


     !overlap
     integer :: overlap
    
     !Vertical localization 
     integer :: num_vert_levels 
     integer :: vert_overlap

     !Descriptor of geographical region (options for now are tropics,
     !extratropics, and polar)
     character(len=:), allocatable :: region_char

     !scale factors when we standardize
     real(kind=dp), allocatable :: mean(:)  !has length n where n=number of
     real(kind=dp), allocatable :: std(:)   !physical components in the state
                                            !vector

     !Index of tisr in the mean and std arrays
     integer :: tisr_mean_std_idx  

     !Index of logp in the mean and std arrays
     integer :: logp_mean_std_idx

     !Index of sst in the mean and std arrays
     integer :: sst_mean_std_idx

     !Index of precip in the mean and std arrays 
     integer :: precip_mean_std_idx

     !Index of oceanic heat content in the mean and std arrays
     integer :: ohtc_mean_std_idx

     !Number of regions 
     integer :: number_of_regions

     !Logical flag to indicate if the grid is the top or bottom of the
     !atmosphere
     logical :: top
     logical :: bottom

     !Integer to tell the grid which level it is could be 1 ... n n=number of
     !vertical levels
     integer :: level_index

     !Boolean whether this reservoir associated with this grid type has the logp
     !variable
     logical :: logp_bool

     !Where the 3-d atmosphere start and stop  indices in u(t) where u looks something like  u(atmos...logp...sst...sst_climo...tisr)
     !start and end will depend if the variable is from the input region u(t) or
     !just
     !the prediction region v(t)
     integer :: atmo3d_start
     integer :: atmo3d_end

     !SST indices in u(t) where u(atmos...logp...sst...sst_climo...tisr) 
     !start and end will depend if the variable is from the input region u(t) or just
     !the prediction region v(t)
     integer :: sst_start
     integer :: sst_end
     
     !logp indices in u(t) where u(atmos...logp...sst...sst_climo...tisr)
     !start and end will depend if the variable is from the input region u(t) or just
     !the prediction region v(t)
     integer :: logp_start
     integer :: logp_end

     !precipitation indices in u(t) where u(atmos...logp...precip...sst...sst_climo...tisr)
     !start and end will depend if the variable is from the input region u(t) or
     !just
     !the prediction region v(t)
     integer :: precip_start
     integer :: precip_end

     !tisr indices in u(t) where u(atmos...logp...sst...sst_climo...tisr)
     !start and end will depend if the variable is from the input region u(t) or
     !just the prediction region v(t)
     integer :: tisr_start
     integer :: tisr_end

     !variables that need predicting indices in u(t) where u(atmos...logp...sst...sst_climo...tisr)
     !so if its a region that only predicts the atmosphere then it will be
     !predict_start = 1, predict_end = atmo3d_end
     !if the reservoir predicts the atmos, logp and sst then predict_start = 1,
     !predict_end = sst_end
     !start and end will depend if the variable is from the input region u(t) or
     !just the prediction region v(t)
     integer :: predict_start
     integer :: predict_end

     !ocean heat content (ohtc) indices in u(t) where
     !u(atmos...logp...sst...sst_climo...tisr...othc)
     !This will only be for ocean reservoirs
     integer :: ohtc_start
     integer :: ohtc_end
  end type grid_type

  type reservoir_type
    !Reservoir derived type that holds all of the variables needed to run a
    !reservoir routine including the training data, reservoir hyper-parameters,
    !and the matrices for fitting the reservoir 

    !Integer that is used to index which region this specific reservoir is
    !assigned. For example 512 processors and 1152 regions and this was
    !processor 1 and region 2 on that processor that assigned_region == 4
    integer                    :: assigned_region

    !Holds the indices of the sigma levels for this specific reservoir 
    integer, allocatable       :: vert_indices_res(:)
    integer, allocatable       :: vert_indices_input(:)


    !Training data for this specific reservior
    real(kind=dp), allocatable :: trainingdata(:,:)

    !reservoir stuff relating to the adjancency matrix
    type(SPARSE_MATRIX_T) cooA
    type(MATRIX_DESCR) descrA

   
    real(kind=dp)              :: deg
    real(kind=dp)              :: radius
    real(kind=dp)              :: beta_res !reservoir ridge regression parameter
    real(kind=dp)              :: beta_model !model ridge regression parameter
    real(kind=dp)              :: density

    !Sigma for input into reservoir
    real(kind=dp)              :: sigma

    !Leakage
    real(kind=dp)              :: leakage 
    real(kind=dp), allocatable :: leakage_slab(:)
    real(kind=dp)              :: leakage_slab_upper
    real(kind=dp)              :: leakage_slab_lower

    !COO sparse matrix holds row and col indexs
    integer, allocatable       :: rows(:)
    integer, allocatable       :: cols(:)
    real(kind=dp), allocatable :: vals(:)

    integer                    :: k
 
    integer                    :: reservoir_numinputs

    integer                    :: locality

    !reservoir node size m=what you want n=closest possible num of nodes to m
    integer                    :: m
    integer                    :: n

    !Reservoir components Win, Wout, etc
    real(kind=dp), allocatable :: win(:,:)
    real(kind=dp), allocatable :: wout(:,:)
    real(kind=dp), allocatable :: states(:,:)
    real(kind=dp), allocatable :: augmented_states(:,:)

    !Batch stuff
    integer                    :: batch_size
    real(kind=dp), allocatable :: states_x_states(:,:)
    real(kind=dp), allocatable :: states_x_trainingdata(:,:)

    !batch hybrid stuff
    real(kind=dp), allocatable :: states_x_states_aug(:,:)
    real(kind=dp), allocatable :: states_x_trainingdata_aug(:,:)


    !Local num of height levels and variables predicted per level
    !may change on level and processor.
    !NOTE full_heightlevels /= local_heightlevels possible
    integer :: local_heightlevels_res
    integer :: local_heightlevels_input
    integer :: local_predictvars

    !Becase logp is a special 2d variable its is complicated
    !and needs its own parameters
    integer :: logp_size_res
    integer :: logp_size_input
  
    !LOGP bool gotta set this to false for vertical reservoirs away from the
    !surface
    logical :: logp_bool

    !Save reservoir state for later use
    real(kind=dp), allocatable :: saved_state(:)

    !Current reservoir state for later use
    real(kind=dp), allocatable :: current_state(:)

    !TISR relievent 
    logical :: tisr_input_bool !This is the actual variable that determines
                               !which reservoir gets tisr. Just because
                               !model_parameter%toa_isr_bool is true doesnt
                               !mean tisr_input_bool is true. It will depend on
                               !what part of the atmosphere the reservoir
                               !predicts

    integer :: tisr_size_input !Saves the size of the tisr input length
    integer :: tisr_size_res !size of non overlap tisr length



    !Variables related to precipitation
    logical :: precip_bool


    logical :: precip_input_bool

    integer :: precip_size_res
    integer :: precip_size_input



    !Becase sst is a special 2d variable its is complicated
    !and needs its own parameters
    logical :: sst_bool !Whether this reservoir should even have SST as a
                        !possible input

    logical :: sst_bool_input !Whether this reservoir has sst data as an input
                              !If there isnt any variance in the domain
                              !this is set to false even if sst_bool is
                              !true

    logical :: sst_bool_prediction !The slab ocean model reservoir should be the
                                   !only one predicting sst data

    integer :: sst_size_res
    integer :: sst_size_input


    !Input SST climatology
    logical :: sst_climo_bool

    logical :: sst_climo_input
    integer :: sst_climo_res

    !Atmosphere to ocean coupling parameters
    logical :: atmo_to_ocean_coupled
   
    integer :: atmo_size_input
    integer :: num_atmo_levels
   
    !Variables related to inputting the 0-300 meter oceanic heat content
    logical :: ohtc_prediction

    integer :: ohtc_res_size
    integer :: ohtc_input_size      

    integer, allocatable :: atmo_training_data_idx(:)

    !A 2d var to hold the atmospheric inputs to the ocean model during the
    !prediction phase that then gets averaged over 
    !This will be held by the ocean reservoirs not the atmo ones 
    real(kind=dp), allocatable :: averaged_atmo_input_vec(:,:) !shape = (num_inputs of atmo variables into ml ocean model,timestep_slab)

    !Standard variables needed for the reservoir and hybrid related to length of
    !the predicted vector 
    integer                    :: chunk_size
    integer                    :: chunk_size_prediction
    integer                    :: chunk_size_speedy !When prognostic variables
                                                    !of the hybrid model arent
                                                    !in speedy

    !Chunk size and chunk size prediction are the same unless toa_isr_bool is
    !TruA

   
    !Contains the imperfect model states
    !e.g. feeding the training data through speedy and saving the output
    real(kind=dp), allocatable :: imperfect_model_states(:,:)

    !Holds predictions by reservoir 
    real(kind=dp), allocatable :: predictiondata(:,:)

    !stuff related to noise
    real(kind=dp) :: noisemag !magnitude of gaussian noise

    !The value of the prior [0,1]
    real(kind=dp) :: prior_val
   
    !Current local speedy state
    real(kind=dp), allocatable :: local_model(:)

    !Current outvec (prediction) for this reservoir
    real(kind=dp), allocatable :: outvec(:)

    !Current ml component contribution of the outvec for this reservoir
    real(kind=dp), allocatable :: v_ml(:)

    !Current SPEEDY contribution of the outvec for this reservoir
    real(kind=dp), allocatable :: v_p(:)

    !Current feedback for this reservoir 
    real(kind=dp), allocatable :: feedback(:)
    
    !TISR for the full year to save time reading and writing during prediction
    real(kind=dp), allocatable :: full_tisr(:,:,:) 

    real(kind=dp), allocatable :: full_sst(:,:,:)

    integer :: predictvars2d
  end type reservoir_type

  type model_parameters_type 
    !Holds the static variables for a specific processor 

    !ML-only Atmo Model Setup 
    logical :: ml_only

    !If this is a ocean reservoir make it ml only if true 
    logical :: ml_only_ocean
    
    !Vars for vertical localization and domain decomposition     
    integer                           :: num_vert_levels
    integer                           :: vert_loc_overlap
    integer                           :: number_of_regions
    integer                           :: num_of_regions_on_proc

    !array that holds the indices of which regions this processor has
    !e.g. 1152 regions and 576 processors then res%irank 0 will have regions 0
    !and 1
    integer, allocatable              :: region_indices(:)

    !Various parameters global variables of the global model state
    !e.g. when set to the main_type these are the max number and are static
    !however each reservoir in the vertical may have a different number than
    !these
    integer :: full_heightlevels
    integer :: full_predictvars

    !data for reservoir
    integer                    :: traininglength

    integer                    :: discardlength
    integer                    :: synclength

    integer                    :: predictionlength

    integer                    :: overlap
    integer, allocatable       :: prediction_markers(:)
    integer                    :: num_predictions
    integer                    :: current_trial_number

    !mpi stuff
    integer :: irank
    integer :: numprocs

    !prediction
    real(kind=dp), allocatable :: prediction(:,:)

    logical                    :: specific_humidity_log_bool
    real(kind=dp)              :: specific_humidity_epsilon=0.3_dp !we are
    !saving specific as log(sp + epsilon)

    !Pole model forecast only
    logical :: pole_only

    !File writing strings
    character(len=3)  :: trial_number
    character(len=10) :: trial_date
    character(len=:), allocatable :: trial_name
    character(len=:), allocatable :: trial_name_extra_end

    !Make sure SPEEDY code doesnt break ours
    logical :: run_speedy

    !Input time of day into reservoir
    logical :: timeofday_bool

    !Vary the hyper-parameters by region in hopes of getting better results
    logical :: regional_vary

    !Number of reading and training chunks
    !This is independent of the batch training stuff
    !The memory limitation will not allow us to read in
    !multiple decades of era and speedy data
    !instead we will chunk the reading and training into parts
    !while also batching the training

    !---Section for variables related to using a prior for the LS solver---!

    !Boolean saying if we are using a prior or not
    logical :: using_prior

    !Noise to be added to the imperfect model states
    real(kind=dp) :: model_noise

    !Time step of the reservoir
    integer :: timestep

    !Time step of the ocean model
    integer :: timestep_slab

    !Input top of atmosphere incident solar radiation into reservoir
    logical :: toa_isr_bool

    !Input and predict precipitation 
    logical :: precip_bool

    !Transforming Precip into log space see Pathek et al 2022
    real :: precip_epsilon

    !stuff related to noise
    logical :: noisy !Flag to add noise

    character(len=:), allocatable :: prediction_file

    !We can assume "special" reservoirs that are different than our main reservoirs
    !which predict the 3d prognostic variables (T,U,V,R) at each height level
    !and grid point
    !
    !These reservoirs and the associated grid object (if present) also wont be defined on the usual xyz grid
    !Things like precip and SST and their grids are significantly different
    !enough to have there own designation of these special reservoirs
    !
    !They can also have there own completely different hyper-parameters such as
    !time step etc
    logical :: special_reservoirs

    !Num of these special reservoirs (e.g. if using precip + SST then
    !num_special_reservoirs == 2) 
    integer :: num_special_reservoirs

    logical :: slab_ocean_model_bool

    logical :: ohtc_bool_input
 
    logical :: train_on_sst_anomalies

    real(kind=dp), allocatable :: base_sst_grid(:,:)
    real(kind=dp), allocatable :: sea_mask(:,:)

    type(opened_netcdf_type), allocatable :: opened_netcdf_files(:)

    logical :: non_stationary_ocn_climo
  
    real(kind=dp) :: final_sst_bias

    real(kind=dp) :: current_sst_bias

    logical :: outvec_component_contribs
  end type model_parameters_type
 
  type main_type
    !The main derive type of the program. This holds basically everything needed
    !to run the program including the model parameters, grid type for each
    !reservoir, and each reservoir it self

    !The shape of grid and reservoir are equal and grid(i,j) is the grid type for
    !reservoir(i,j). Shape(reservoir) is (num_regions_per_processor,num_vertical_levels)
     
 
    !Allocatable grid type for multiple sub-domains per processor
    type(grid_type), allocatable      :: grid(:,:)
      
    !Allocatable reservoir type for vert loc
    type(reservoir_type), allocatable :: reservoir(:,:)

    !Allocatable grid type for multiple sub-domains per processor
    type(grid_type), allocatable      :: grid_special(:,:)

    !Allocatable reservoir type for vert loc
    type(reservoir_type), allocatable :: reservoir_special(:,:)

    !A simple derived type to hold the static parameters of the model 
    !needed to passing to subroutines were reservoir and or grid are also passed
    type(model_parameters_type)       :: model_parameters 

  end type main_type

  type speedy_data_type
    !Holds all of the 4d speedy variables in one 5d array
    real(kind=dp), allocatable :: speedyvariables(:,:,:,:,:)
     
    !LogP is a special speedy variable
    real(kind=dp), allocatable :: speedy_logp(:,:,:)
  end type speedy_data_type

  type era_data_type
    !Holds all of the 4d era variables in one 5d array
     real(kind=dp), allocatable :: eravariables(:,:,:,:,:)

     !LogP is a special era variable
     real(kind=dp), allocatable :: era_logp(:,:,:)

     !TISR
     real(kind=dp), allocatable :: era_tisr(:,:,:)
 
     !SST
     real(kind=dp), allocatable :: era_sst(:,:,:)

     !Climatological SST (hourly average iof ERA5 from 1981 - 2000)
     real(kind=dp), allocatable :: era_sst_climo(:,:,:)

     !hourly Precip
     real(kind=dp), allocatable :: era_precip(:,:,:)

     !Oceanic HeaT Content
     real(kind=dp), allocatable :: era_ohtc(:,:,:)
  end type era_data_type

  type state_vector_type
     !Holds all of the 3d variables in one 4d array
     real(kind=dp), allocatable :: variables3d(:,:,:,:)

     !LogP is a special variable
     real(kind=dp), allocatable :: logp(:,:)
  
     !Lets also keep some meta data in this type to make SPEEDY run
     integer :: istart !from rest=0, from restartfile=1, from era=2
     integer :: era_start !start from grid initial condition=0, Start from grid=1
     character(len=100) :: era_file

     integer :: era_hour !era hour of the month 1 = 00UTC of the first day of
     integer :: era_hour_plus_one !So I dont have to do calendar stuff in

      !Speedy internal calenadr stuff
     integer :: iyear0
     integer :: imont0
     integer :: iday
     integer :: ihour 

     !A logical flag to see if it is safe to run speedy anymore
     logical :: is_safe_to_run_speedy

     logical :: hybrid_slab

     real(kind=dp), allocatable :: sst_hybrid(:,:)

     real(kind=dp) :: sst_bias = 0.0_dp
  end type state_vector_type

  type mpi_type
     !Holds the necessary mpi stuff one per processor
     integer(kind=int32) :: ierr
     integer(kind=int32) :: numprocs
     integer(kind=int32) :: proc_num
     integer :: mpi_world
     !type(MPI_Comm) :: mpi_world

     logical :: is_root = .False.

     logical :: is_serial  = .False.
  end type mpi_type

  type calendar_type
     !Calendar derived type to do calendar calculations and keep track of
     !current and starting times 

     !Initialize calendar_type with start date
     integer :: startyear
     integer :: startmonth
     integer :: startday
     integer :: starthour
  
     !Current date 
     integer :: currentyear
     integer :: currentmonth
     integer :: currentday
     integer :: currenthour 
 
  end type calendar_type 

  type opened_netcdf_type
     logical :: is_opened
     logical :: is_closed

     integer :: ncid
     
     character(len=:), allocatable :: filename
  end type opened_netcdf_type

  !Overload the standardize_data routine 
  interface standardize_data
   module procedure standardize_data_1d
   module procedure standardize_data_2d
   module procedure standardize_data_3d
   module procedure standardize_data_4d
   module procedure standardize_data_5d
   module procedure standardize_data_5d_logp
   module procedure standardize_data_5d_logp_tisr
   module procedure standardize_data_5d_logp_tisr_local_region
  end interface 

  !Overload unstandardize_data routine  
  interface unstandardize_data
    module procedure unstandardize_data_4d
    module procedure unstandardize_data_4d_logp
    module procedure unstandardize_data_4d_multi_2d
  end interface

  !Overloads gaussian_noise routine 
  interface gaussian_noise
    module procedure gaussian_noise_2d
    module procedure gaussian_noise_1d
  end interface 

  contains
    
    subroutine unstandardize_data_4d(reservoir,inputdata,mean,std)
      !Unstandardize inputdata
      !mean and std needs to be the same length as the first dimension of
      !inputdata
      type(reservoir_type), intent(in) :: reservoir
      real(kind=dp), intent(inout)     :: inputdata(:,:,:,:)
      real(kind=dp), intent(in)        :: mean(:), std(:)

      integer                          :: i, j, k
      integer                          :: n, m, l

      n = size(inputdata,1)

      k = size(inputdata,4) !Number of height variables

      !If we are doing log(specific_humidity) no need to unstandardize because
      !we didnt standardize it
      !if(res%specific_humidity_log_bool) then
        !For plotting and regional variety lets make sure we are out of log
        !space
        !inputdata(4,:,:,:) = e_constant**inputdata(4,:,:,:) - res%specific_humidity_epsilon
      !endif

      l = 1
      do i=1, n
         do j=1, k
            inputdata(i,:,:,j) = inputdata(i,:,:,j)*std(l)
            inputdata(i,:,:,j) = inputdata(i,:,:,j) + mean(l)
            l = l + 1
         enddo
      end do

      return
    end subroutine

    subroutine unstandardize_data_4d_logp(reservoir,inputdata,logp,mean,std)
      !Unstandardize inputdata
      !mean and std needs to be the same length as the first dimension of
      !inputdata
      type(reservoir_type), intent(in)       :: reservoir

      real(kind=dp), intent(inout)           :: inputdata(:,:,:,:)
      real(kind=dp), intent(inout)           :: logp(:,:)
      real(kind=dp), intent(in)              :: mean(:), std(:)

      integer                          :: i, j, k
      integer                          :: n, m, l

      n = size(inputdata,1)

      k = size(inputdata,4) !Number of height variables

      l = 1
      do i=1, n
         do j=1, k
            inputdata(i,:,:,j) = inputdata(i,:,:,j)*std(l)
            inputdata(i,:,:,j) = inputdata(i,:,:,j) + mean(l)
            l = l + 1
         enddo
      end do

      logp = logp*std(l)
      logp = logp + mean(l)
      return
    end subroutine

    subroutine unstandardize_data_4d_logp_2d(reservoir,inputdata,logp,mean,std)
      !Unstandardize inputdata
      !mean and std needs to be the same length as the first dimension of
      !inputdata
      type(reservoir_type), intent(in)       :: reservoir

      real(kind=dp), intent(inout)           :: inputdata(:,:,:,:,:)
      real(kind=dp), intent(inout)           :: logp(:,:,:)
      real(kind=dp), intent(in)              :: mean(:), std(:)

      integer                          :: i, j, k
      integer                          :: n, m, l

      n = size(inputdata,1)

      k = size(inputdata,4) !Number of height variables

      l = 1
      do i=1, n
         do j=1, k
            inputdata(i,:,:,j,:) = inputdata(i,:,:,j,:)*std(l)
            inputdata(i,:,:,j,:) = inputdata(i,:,:,j,:) + mean(l)
            l = l + 1
         enddo
      end do

      logp = logp*std(l)
      logp = logp + mean(l)
      return
    end subroutine

    subroutine unstandardize_data_4d_multi_2d(reservoir,inputdata,inputdata2d,mean,std)
      !Unstandardize inputdata
      !mean and std needs to be the same length as the first dimension of
      !inputdata
      type(reservoir_type), intent(in)       :: reservoir
      real(kind=dp), intent(inout)           :: inputdata(:,:,:,:)
      real(kind=dp), intent(inout)           :: inputdata2d(:,:,:)
      real(kind=dp), intent(in)              :: mean(:), std(:)

      integer                          :: i, j, k
      integer                          :: n, m, l

      n = size(inputdata,1)

      k = size(inputdata,4) !Number of height variables


      l = 1
      do i=1, n
         do j=1, k
            inputdata(i,:,:,j) = inputdata(i,:,:,j)*std(l)
            inputdata(i,:,:,j) = inputdata(i,:,:,j) + mean(l)
            l = l + 1
         enddo
      end do

      do i=1, reservoir%predictvars2d
         !print *, 'i',i,'l',l,'mean,std',std(l),mean(l)
         inputdata2d(i,:,:) = inputdata2d(i,:,:)*std(l)
         inputdata2d(i,:,:) = inputdata2d(i,:,:) +  mean(l)
         l = l + 2
      enddo
      return
    end subroutine

    subroutine unstandardize_data_2d(inputdata,mean,std)
      !Unstandardize 2-dimensional inputdata
      !mean and std needs to be the same length as the first dimension of
      !inputdata

      real(kind=dp), intent(inout)           :: inputdata(:,:)
      real(kind=dp), intent(in)              :: mean, std

      integer                          :: i, j, k
      integer                          :: n, m, l

      inputdata = inputdata*std
      inputdata = inputdata + mean

      return
    end subroutine

    subroutine unstandardize_data_1d(inputdata,mean,std)
      !Unstandardize 2-dimensional inputdata
      !mean and std needs to be the same length as the first dimension of
      !inputdata

      real(kind=dp), intent(inout)           :: inputdata(:)
      real(kind=dp), intent(in)              :: mean, std

      integer                          :: i, j, k
      integer                          :: n, m, l

      inputdata = inputdata*std
      inputdata = inputdata + mean

      return
    end subroutine

    subroutine standardize_data_2d(inputdata,mean,std)
      !Standardizes the input data by subtracting the mean out 
      !and then dividing by the std 
      real(kind=dp), intent(inout) :: inputdata(:,:) !2d var with data being 1st and 2nd being time
      real(kind=dp), intent(out)   :: mean, std

      real(kind=dp)                :: mean_, std_

      mean_ = sum(inputdata)/size(inputdata)
      std_ = sqrt((sum(inputdata**2)-sum(inputdata)**2/size(inputdata))/size(inputdata))  

      !Subtracts the mean out and then divides by the std
      inputdata = inputdata - mean_
      inputdata = inputdata/std_
      
      mean = mean_
      std = std_
      return 
    end subroutine

    subroutine standardize_sst_data_3d(inputdata,mean,std,any_change)
      !Standardizes the input data by subtracting the mean out
      !and then dividing by the std
      real(kind=dp), intent(inout) :: inputdata(:,:,:) !3d var with data being 1st and 3nd being time
      real(kind=dp), intent(out)   :: mean, std

      logical, intent(out)         :: any_change

      real(kind=dp)                :: mean_, std_


      if((sum(inputdata**2)-sum(inputdata)**2/size(inputdata)) > 0) then
        mean_ = sum(inputdata)/size(inputdata)
 
        mean_ = mean_ !+ 2.0
        !sqrt(sum((inputdata(i,:,:,j,:) - mean_)**2)/size(inputdata(i,:,:,j,:)))
        std_ = sqrt(sum((inputdata-mean_)**2)/size(inputdata))
 
        std_ = std_ !* 1.2
        print *, 'std_',std_
        if(std_ > 0.2) then
          !Subtracts the mean out and then divides by the std
          inputdata = inputdata - mean_
          inputdata = inputdata/std_

          mean = mean_
          std = std_
          any_change = .True.
        else
          mean = 0
          std = 0
          any_change = .False.
         endif
      else
        mean = 0
        std = 0
        any_change = .False.
      endif
      return
    end subroutine
 
    subroutine standardize_data_3d(inputdata,mean,std)
      !Standardizes the input data by subtracting the mean out
      !and then dividing by the std
      real(kind=dp), intent(inout) :: inputdata(:,:,:) !3d var with data being 1st and 3nd being time
      real(kind=dp), intent(out)   :: mean, std

      real(kind=dp)                :: mean_, std_

      mean_ = sum(inputdata)/size(inputdata)
      std_ = sqrt((sum(inputdata**2)-sum(inputdata)**2/size(inputdata))/size(inputdata))

      !Subtracts the mean out and then divides by the std
      inputdata = inputdata - mean_
      inputdata = inputdata/std_

      mean = mean_
      std = std_
      return
    end subroutine

    subroutine standardize_data_4d(inputdata,mean,std)
      !Standardizes the input data by subtracting the mean out
      !and then dividing by the std
      real(kind=dp), intent(inout) :: inputdata(:,:,:,:) !4d var with data being 1-3 and 4th being time
      real(kind=dp), intent(out)   :: mean, std

      real(kind=dp)                :: mean_, std_

      mean_ = sum(inputdata)/size(inputdata)
      std_ = sqrt((sum(inputdata**2)-sum(inputdata)**2/size(inputdata))/size(inputdata))

      !Subtracts the mean out and then divides by the std
      inputdata = inputdata - mean_
      inputdata = inputdata/std_

      mean = mean_
      std = std_
      return
    end subroutine 

    subroutine standardize_data_5d(reservoir,inputdata,mean,std)
      type(reservoir_type), intent(in) :: reservoir

      real(kind=dp), intent(inout)     :: inputdata(:,:,:,:,:) !5d var with data being 1 being a variable type 2-4 xyz and 5th being time
      real(kind=dp), intent(out)       :: mean(:), std(:)

      real(kind=dp)                :: mean_, std_

      integer       :: length, height
      integer       :: i, j, l

      length = size(inputdata,1)
      height = size(inputdata,4)

      l = 1
      do i=1, length
         do j=1, height
            mean_ = sum(inputdata(i,:,:,j,:))/size(inputdata(i,:,:,j,:))
            std_ = sqrt(sum((inputdata(i,:,:,j,:) - mean_)**2)/size(inputdata(i,:,:,j,:)))

            call standardize_data_given_pars3d(inputdata(i,:,:,j,:),mean_,std_)

            mean(l) = mean_
            std(l) = std_
            l = l + 1
         enddo
      end do

      return
    end subroutine

    subroutine standardize_data_1d(inputdata,mean,std)
      !Standardizes the input data by subtracting the mean out
      !and then dividing by std
      real(kind=dp), intent(inout) :: inputdata(:) !1d var data var
      real(kind=dp), intent(out)   :: mean, std

      real(kind=dp)                :: mean_, std_

      mean_ = sum(inputdata)/size(inputdata)
      std_ = sqrt((sum(inputdata**2)-sum(inputdata)**2/size(inputdata))/size(inputdata))

      !Subtracts the mean out and then divides by the std
      inputdata = inputdata - mean_
      inputdata = inputdata/std_

      mean = mean_
      std = std_
      return
    end subroutine 

    subroutine standardize_data_5d_logp(reservoir,inputdata,logp,mean,std)
      type(reservoir_type), intent(in) :: reservoir

      real(kind=dp), intent(inout)     :: inputdata(:,:,:,:,:) !5d var with data being 1 being a variable type 2-4 xyz and 5th being time
      real(kind=dp), intent(inout)     :: logp(:,:,:)
      real(kind=dp), intent(out)       :: mean(:), std(:)

      real(kind=dp)                :: mean_, std_

      integer       :: length, height
      integer       :: i, j, l

      length = size(inputdata,1)
      height = size(inputdata,4)

      l = 1
      do i=1, length
         do j=1, height
            mean_ = sum(inputdata(i,:,:,j,:))/size(inputdata(i,:,:,j,:))
            std_ = sqrt(sum((inputdata(i,:,:,j,:) - mean_)**2)/size(inputdata(i,:,:,j,:)))

            call standardize_data_given_pars3d(inputdata(i,:,:,j,:),mean_,std_)

            mean(l) = mean_
            std(l) = std_
            l = l + 1
         enddo
      end do

      mean_ = sum(logp)/size(logp)
      std_ = sqrt(sum((logp - mean_)**2)/size(logp))

      call standardize_data_given_pars3d(logp,mean_,std_)

      mean(l) = mean_
      std(l) = std_

      return
    end subroutine

    subroutine standardize_data_5d_logp_tisr_local_region(reservoir,grid,inputdata,logp,tisr,mean,std)
      type(reservoir_type), intent(in) :: reservoir 
      type(grid_type), intent(in)      :: grid

      real(kind=dp), intent(inout)     :: inputdata(:,:,:,:,:) !5d var with data being 1 being a variable type 2-4 xyz and 5th being time
      real(kind=dp), intent(inout)     :: logp(:,:,:), tisr(:,:,:)
      real(kind=dp), intent(out)       :: mean(:), std(:)

      real(kind=dp)                :: mean_, std_

      integer       :: length, height
      integer       :: i, j, l

      length = size(inputdata,1)
      height = size(inputdata,4)

      l = 1
      do i=1, length
         do j=1, height
            mean_ = sum(inputdata(i,grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,j,:))/size(inputdata(i,grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,j,:))
            std_ = sqrt(sum((inputdata(i,grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,j,:) - mean_)**2)/size(inputdata(i,grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,j,:)))

            call standardize_data_given_pars3d(inputdata(i,:,:,j,:),mean_,std_)

            mean(l) = mean_
            std(l) = std_
            l = l + 1
         enddo
      end do

      !logp
      mean_ =sum(logp(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:))/size(logp(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:))
      std_ = sqrt(sum((logp(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:) - mean_)**2)/size(logp(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:)))

      call standardize_data_given_pars3d(logp,mean_,std_)

      mean(l) = mean_
      std(l) = std_

      !Tisr
      l = l + 1

      mean_ = sum(tisr(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:))/size(tisr(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:))
      std_ = sqrt(sum((tisr(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:) - mean_)**2)/size(tisr(grid%tdata_xstart:grid%tdata_xend,grid%tdata_ystart:grid%tdata_yend,:)))

      call standardize_data_given_pars3d(tisr,mean_,std_)

      mean(l) = mean_
      std(l) = std_
      return
    end subroutine 
   
    subroutine standardize_data_5d_logp_tisr_gp_by_gp(reservoir,grid,inputdata,logp,tisr,mean,std)
      type(reservoir_type), intent(in) :: reservoir
      type(grid_type), intent(in)      :: grid

      real(kind=dp), intent(inout)     :: inputdata(:,:,:,:,:) !5d var with data being 1 being a variable type 2-4 xyz and 5th being time
      real(kind=dp), intent(inout)     :: logp(:,:,:), tisr(:,:,:)
      real(kind=dp), intent(out)       :: mean(:), std(:)

      real(kind=dp)                :: mean_, std_

      integer       :: vars, height, width, length
      integer       :: i, j, k, l, counter

      vars = size(inputdata,1)
      width = size(inputdata,2)
      length = size(inputdata,3)
      height = size(inputdata,4)
      

      counter = 1
      do l=1, vars
         do i=1, width
            do j=1, length
               do k=1, height
                  mean_ = sum(inputdata(l,i,j,k,:))/size(inputdata(l,i,j,k,:))
                  std_ = sqrt(sum((inputdata(l,i,j,k,:) - mean_)**2)/size(inputdata(l,i,j,k,:)))

                  call standardize_data_given_pars1d(inputdata(l,i,j,k,:),mean_,std_)

                  mean(counter) = mean_
                  std(counter) = std_
                  counter = counter + 1
               enddo 
            enddo 
         enddo
      end do

      !logp
      do i=1, width
         do j=1, length
            mean_ =sum(logp(i,j,:))/size(logp(i,j,:))
            std_ = sqrt(sum((logp(i,j,:) - mean_)**2)/size(logp(i,j,:)))

            call standardize_data_given_pars1d(logp(i,j,:),mean_,std_)

            mean(counter) = mean_
            std(counter) = std_ 
  
            counter = counter + 1
         enddo 
      enddo 

      !Tisr
      do i=1, width
         do j=1, length
            mean_ =sum(tisr(i,j,:))/size(tisr(i,j,:))
            std_ = sqrt(sum((tisr(i,j,:) - mean_)**2)/size(tisr(i,j,:)))

            call standardize_data_given_pars1d(tisr(i,j,:),mean_,std_)

            mean(counter) = mean_
            std(counter) = std_

            counter = counter + 1
         enddo
      enddo
    end subroutine 
    subroutine standardize_data_5d_logp_tisr(reservoir,inputdata,logp,tisr,mean,std)
      type(reservoir_type), intent(in) :: reservoir

      real(kind=dp), intent(inout)     :: inputdata(:,:,:,:,:) !5d var with data being 1 being a variable type 2-4 xyz and 5th being time
      real(kind=dp), intent(inout)     :: logp(:,:,:), tisr(:,:,:)
      real(kind=dp), intent(out)       :: mean(:), std(:)

      real(kind=dp)                :: mean_, std_

      integer       :: length, height
      integer       :: i, j, l

      length = size(inputdata,1)
      height = size(inputdata,4)

      l = 1
      do i=1, length
         do j=1, height
            mean_ = sum(inputdata(i,:,:,j,:))/size(inputdata(i,:,:,j,:))
            std_ = sqrt(sum((inputdata(i,:,:,j,:) - mean_)**2)/size(inputdata(i,:,:,j,:)))

            call standardize_data_given_pars3d(inputdata(i,:,:,j,:),mean_,std_)

            mean(l) = mean_
            std(l) = std_
            l = l + 1
         enddo
      end do

      !logp
      mean_ = sum(logp)/size(logp)
      std_ = sqrt(sum((logp - mean_)**2)/size(logp))

      call standardize_data_given_pars3d(logp,mean_,std_)
  
      mean(l) = mean_
      std(l) = std_

      !Tisr
      l = l + 1

      mean_ = sum(tisr)/size(tisr)
      std_ = sqrt(sum((tisr - mean_)**2)/size(tisr))

      call standardize_data_given_pars3d(tisr,mean_,std_)

      mean(l) = mean_
      std(l) = std_
      return
    end subroutine
    
    subroutine standardize_data_given_pars_5d_logp_tisr(mean,std,input_data,input_logp,input_tisr)
      real(kind=dp), intent(in)        :: mean(:), std(:)

      real(kind=dp), intent(inout)     :: input_data(:,:,:,:,:) !5d var with data being 1 being a variable type 2-4 xyz and 5th being time
      real(kind=dp), intent(inout)     :: input_logp(:,:,:), input_tisr(:,:,:)

      integer       :: length, height
      integer       :: i, j, l

      length = size(input_data,1)
      height = size(input_data,4)

      l = 1
      do i=1, length
         do j=1, height

            call standardize_data_given_pars3d(input_data(i,:,:,j,:),mean(l),std(l))

            l = l + 1
         enddo
      end do

      !Logp
      call standardize_data_given_pars3d(input_logp,mean(l),std(l))
  
      !Tisr
      l = l + 1
  
      call standardize_data_given_pars3d(input_tisr,mean(l),std(l))

      return
    end subroutine
 
    subroutine standardize_data_given_pars_5d_logp(mean,std,input_data,input_logp,mean_logp,std_logp)
      real(kind=dp), intent(in)           :: mean(:), std(:)
      real(kind=dp), intent(in), optional :: mean_logp, std_logp
      
      real(kind=dp), intent(inout)     :: input_data(:,:,:,:,:) !5d var with data being 1 being a variable type 2-4 xyz and 5th being time
      real(kind=dp), intent(inout)     :: input_logp(:,:,:)

      integer       :: length, height
      integer       :: i, j, l

      length = size(input_data,1)
      height = size(input_data,4)

      l = 1
      do i=1, length
         do j=1, height

            call standardize_data_given_pars3d(input_data(i,:,:,j,:),mean(l),std(l))

            l = l + 1
         enddo
      end do

      if((present(mean_logp).and.present(std_logp))) then
        call standardize_data_given_pars3d(input_logp,mean_logp,std_logp)
      else
        call standardize_data_given_pars3d(input_logp,mean(l),std(l))
      endif
      return
    end subroutine

    subroutine standardize_data_given_pars5d(mean,std,input_data)
      real(kind=dp), intent(in)        :: mean(:), std(:)

      real(kind=dp), intent(inout)     :: input_data(:,:,:,:,:) !5d var with data being 1 being a variable type 2-4 xyz and 5th being time

      integer       :: length, height
      integer       :: i, j, l

      length = size(input_data,1)
      height = size(input_data,4)

      l = 1
      do i=1, length
         do j=1, height

            call standardize_data_given_pars3d(input_data(i,:,:,j,:),mean(l),std(l))

            l = l + 1
         enddo
      end do

      return
    end subroutine

    subroutine standardize_data_given_pars4d(inputdata,mean,std)
      !Standardizes input data by the given std and mean
      real(kind=dp), intent(inout) :: inputdata(:,:,:,:)
      real(kind=dp), intent(in)    :: mean, std

      !Subtracts the mean out and then divides by the std
      inputdata = inputdata - mean
      inputdata = inputdata/std
   
      return 
    end subroutine 

    subroutine standardize_data_given_pars3d(inputdata,mean,std)
      !Standardizes input data by the given std and mean
      real(kind=dp), intent(inout) :: inputdata(:,:,:)
      real(kind=dp), intent(in)    :: mean, std

      !Subtracts the mean out and then divides by the std
      inputdata = inputdata - mean
      inputdata = inputdata/std

      return
    end subroutine

    subroutine standardize_data_given_pars2d(inputdata,mean,std)
      !Standardizes input data by the given std and mean
      real(kind=dp), intent(inout) :: inputdata(:,:)
      real(kind=dp), intent(in)    :: mean, std

      !Subtracts the mean out and then divides by the std
      inputdata = inputdata - mean
      inputdata = inputdata/std

      return
    end subroutine

    subroutine standardize_data_given_pars1d(inputdata,mean,std)
      !Standardizes input data by the given std and mean
      real(kind=dp), intent(inout) :: inputdata(:)
      real(kind=dp), intent(in)    :: mean, std

      !Subtracts the mean out and then divides by the std
      inputdata = inputdata - mean
      inputdata = inputdata/std

      return
    end subroutine

    subroutine normalize_data(inputdata,min_data,max_data)
      !normalize the input data (inputdata-min)/(max-min)
      real(kind=dp), intent(inout) :: inputdata(:,:)
      real(kind=dp), intent(out)   :: min_data,max_data
      
      min_data = minval(inputdata)
      max_data = maxval(inputdata)

      inputdata = (inputdata - min_data)/(max_data-min_data)
     
      return
    end subroutine 
     
           
    subroutine gaussian_noise_2d(inputdata,noisemag)
      !Adds gaussian noise to the input data
      real(kind=dp), intent(inout) :: inputdata(:,:)
      real(kind=dp), intent(in)    :: noisemag

      real(kind=dp), allocatable   :: gaussnoise(:,:)
      real(kind=dp), parameter     :: sigma=1.0, mean=0.0

      allocate(gaussnoise(size(inputdata,1),size(inputdata,2)))

      !call init_random_seed(33)
     
      call random_gaussian_gen_2d(gaussnoise,sigma,mean)

      inputdata = inputdata+gaussnoise*noisemag*inputdata

      deallocate(gaussnoise)

      return
    end subroutine 

    subroutine gaussian_noise_1d(inputdata,noisemag)
      !Adds gaussian noise to the input data
      real(kind=dp), intent(inout) :: inputdata(:)
      real(kind=dp), intent(in)    :: noisemag

      real(kind=dp), allocatable   :: gaussnoise(:)
      real(kind=dp), parameter     :: sigma=1.0, mean=0.0

      allocate(gaussnoise(size(inputdata,1)))

      !call init_random_seed(33)

      call random_gaussian_gen_1d(gaussnoise,sigma,mean)

      inputdata = inputdata+gaussnoise*noisemag*inputdata

      deallocate(gaussnoise)

      return
    end subroutine
   
    function gaussian_noise_1d_function(inputdata,noisemag) result(noisy_data)
      !Adds gaussian noise to the input data
      real(kind=dp), intent(in)  :: inputdata(:)
      real(kind=dp), intent(in)  :: noisemag
      real(kind=dp), allocatable  :: noisy_data(:)

      real(kind=dp), allocatable   :: gaussnoise(:)
      real(kind=dp), parameter     :: sigma=1.0, mean=0.0

      allocate(gaussnoise(size(inputdata,1)))
      allocate(noisy_data(size(inputdata,1)))

      !call init_random_seed(33)

      call random_gaussian_gen_1d(gaussnoise,sigma,mean)

      !TODO
      noisy_data = inputdata+gaussnoise*noisemag*inputdata

      deallocate(gaussnoise)

      return
    end function
  
    function gaussian_noise_1d_function_precip(inputdata,noisemag,grid,model_parameters) result(noisy_data)
      !Adds gaussian noise to the input data and makes sure the noise to
      !precipitation is not done in log space
      real(kind=dp), intent(in)  :: inputdata(:)
      real(kind=dp), intent(in)  :: noisemag

      type(grid_type), intent(in) :: grid
      type(model_parameters_type), intent(in) :: model_parameters

      real(kind=dp), allocatable  :: noisy_data(:)

      real(kind=dp), allocatable   :: gaussnoise(:)
      real(kind=dp), allocatable   :: temp(:)

      real(kind=dp), parameter     :: sigma=1.0, mean=0.0

      allocate(gaussnoise(size(inputdata,1)))
      allocate(noisy_data(size(inputdata,1)))

      call random_gaussian_gen_1d(gaussnoise,sigma,mean)

      !!!!TODO NOTE 
      noisy_data(1:grid%precip_start-1) = inputdata(1:grid%precip_start-1)+gaussnoise(1:grid%precip_start-1)*noisemag*inputdata(1:grid%precip_start-1)

      !Precip stuff

      allocate(temp(grid%precip_end - grid%precip_start))

      temp = inputdata(grid%precip_start:grid%precip_end)

      temp = temp*grid%std(grid%precip_mean_std_idx) + grid%mean(grid%precip_mean_std_idx)

      temp = model_parameters%precip_epsilon * (e_constant**temp - 1)

      temp = temp + gaussnoise(grid%precip_start:grid%precip_end)*noisemag*temp

      temp = abs(temp) !NOTE make sure we dont get any negative numbers

      temp = log(1 + temp/model_parameters%precip_epsilon)

      temp = temp - grid%mean(grid%precip_mean_std_idx)

      temp = temp/grid%std(grid%precip_mean_std_idx)

      noisy_data(grid%precip_start:grid%precip_end) = temp

 
      !Rest of stuff
      noisy_data(grid%precip_end + 1:size(inputdata,1)) = inputdata(grid%precip_end + 1:size(inputdata,1))+gaussnoise(grid%precip_end + 1:size(inputdata,1))*noisemag*inputdata(grid%precip_end + 1:size(inputdata,1))

      deallocate(gaussnoise) 
     
      return
    end function
 
    subroutine random_gaussian_gen_2d(array,sigma,mean)
       !Returns a 2d array of normally distributed 
       !random numbers and gives you control over sigma
       !mean
       real(kind=dp), intent(inout) :: array(:,:)
       real(kind=dp), intent(in)    :: sigma, mean

       integer                      :: n, m, i, j
       real(kind=dp)                :: noise
       
       
       !dims of array
       n = size(array,1)
       m = size(array,2)
       
       do i=1,n
          do j=1,m 

             noise = gaussian_noise_maker(mean,sigma)
              
             !make sure its not greater than 2 sigma away
             !if((noise > 2.0).or.(noise < -2.0)) then
             !  noise = gaussian_noise_maker(mean,sigma)   
             !endif 
             
             array(i,j) = noise
           enddo 
       enddo
       return 
    end subroutine

    subroutine random_gaussian_gen_1d(array,sigma,mean)
       !Returns a 2d array of normally distributed
       !random numbers and gives you control over sigma
       !mean
       real(kind=dp), intent(inout) :: array(:)
       real(kind=dp), intent(in)    :: sigma, mean

       integer                      :: n, m, i, j
       real(kind=dp)                :: noise


       !dims of array
       n = size(array,1)

       do i=1,n
          noise = gaussian_noise_maker(mean,sigma)

          array(i) = noise
       enddo
       return
    end subroutine

    function gaussian_noise_maker(mean,sigma) result(noise)
      !Box-mueller method to get gaussian distributed noise

      real(kind=dp), intent(in)    :: sigma, mean 
      real(kind=dp)                :: noise

      real(kind=dp)                :: u1, u2
      real(kind=dp), parameter     :: twopi=8.0_dp*atan(1.0_dp)
     
      call random_number(u1)
      call random_number(u2)

      noise = mean + sigma*sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
      return 
    end function 

    subroutine init_random_seed(worker)
      !Get random seed not thread safe!!
      integer :: i, n, clock, worker
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = clock + (18+worker*12) * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)

      deallocate(seed)

      return 
    end

    subroutine init_random_marker(input)
      !Get random seed thats thread safe 
      integer :: i, n, clock, input
      integer, dimension(:), allocatable :: seed
 
      call random_seed(size = n)
      allocate(seed(n))
 
      seed = (3+input*2) * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
 
      deallocate(seed)

      return
    end

    subroutine shuffle(n,returnsize,shufflereturn)
      !K-shuffle used to get random choice 
      integer, intent(in) :: n, returnsize

      integer, dimension(1:n) :: choices, choiceshuffle
      integer, dimension(returnsize), intent(out):: shufflereturn
      real :: a
 
      integer :: n_chosen
      integer :: this
      integer :: tmp
      integer :: i

      choices = (/ ( i, i = 1, n ) /)
      n_chosen = 0
      do i = 1, n
         call random_number( a )
         this = a * ( n - n_chosen ) + 1
         tmp = choices( this )
         choiceshuffle(i) = tmp
         choices( this ) = choices( n - n_chosen )
         choices( n - n_chosen ) = tmp
         n_chosen = n_chosen + 1
      enddo
      shufflereturn = choiceshuffle(1:returnsize)

      return 
    end subroutine
    
    subroutine find_closest_divisor(target_,number_,divisor) 
      !Find closest divisor to target_ for number_
      !Example target_ = 7 number_= 64 then divisor = 8
      
      integer, intent(in)   :: target_,number_
      integer, intent(out)  :: divisor
 
      integer               :: i, radius
     
      logical               :: not_found
      
      if(modulo(number_,target_) == 0) then
        !Horray number_ is divisible by target_
        divisor = target_

      else 
        !While loop until divisor is found
        !We loop through potential divisors 
        !until we find one

        not_found = .True.    
        radius = 2

        do while(not_found)
     
           do i=target_-radius,target_+radius
              if(mod(number_,i) == 0) then
                 divisor = i
                 not_found = .False.
                 exit 
              endif 
           enddo 
           radius = radius + 1
        enddo     
     
      endif  

      return 
    end subroutine   
    
    subroutine lorenz63(stepcount,dt,xp,yp,zp)
    integer, intent(in) :: stepcount
    integer :: i 
    real(kind=dp), intent(in) :: dt
    real(kind=dp), parameter :: xi=0.0_dp, yi=1.0_dp, zi=1.05_dp, alpha=0.0_dp
    real(kind=dp) :: xdot, ydot, zdot
    real(kind=dp), intent(out) ::  xp(stepcount+1), yp(stepcount+1),zp(stepcount+1)
    
    xp(1) = xi
    yp(1) = yi
    zp(1) = zi

    do i=1,stepcount
       call advancelorenz(xp(i),yp(i),zp(i),xdot,ydot,zdot,alpha)
       
       xp(i+1) = xp(i) + (xdot * dt)
       yp(i+1) = yp(i) + (ydot * dt)
       zp(i+1) = zp(i) + (zdot * dt)
    enddo 
    return 
    
    end subroutine  

    subroutine advancelorenz(x,y,z,xdot,ydot,zdot,alpha)
      real(kind=dp), parameter :: s=10.0_dp, r=28.0_dp, b=2.66667_dp
      real(kind=dp), intent(in) :: x, y, z, alpha
      real(kind=dp), intent(out) :: xdot, ydot, zdot 
      xdot = s*(y - x)
      ydot = (r+alpha)*x - y - x*z
      zdot = x*y - b*z
      return

    end subroutine     
   
    subroutine tick(t)
      integer, intent(OUT) :: t

      call system_clock(t)
    end subroutine tick

    ! returns time in seconds from now to time described by t
    real function tock(t)
      integer, intent(in) :: t
      integer :: now, clock_rate

      call system_clock(now,clock_rate)

      tock = real(now - t)/real(clock_rate)
    end function tock 

    subroutine total_precip_over_a_period(precip_grid,period)
      !This routine takes hourly precip and will compute the total precip
      !between i-period to i. Will assume the first period elements of precip
      !will be padded with 0
      !
      !Period should equal timestep of the hybrid model
      !
      !Example precip_grid = [1,0,0,0,0,0,0,1,1,0,0,0,0,0,0] and period = 6
      !hours
      !This routine will output precip_grid = [1,1,1,1,1,1,1,0,2,2,2,2,2,2,1,0]
      real(kind=dp), intent(inout) :: precip_grid(:,:,:) !Hourly precip data (x,y,t) that
                                                         !at exit will be total
                                                         !precip over a period

      integer, intent(in) :: period

      !local vars
      real(kind=dp), allocatable :: copy(:,:,:)

      integer :: i, j, t
      integer :: x_len, y_len, t_len

      x_len = size(precip_grid,1)
      y_len = size(precip_grid,2)
      t_len = size(precip_grid,3)

      allocate(copy,source=precip_grid)

      print *, 'shape(precip_grid)',shape(precip_grid)
      do i=1, x_len
         do j=1, y_len
            do t=1, t_len
               if(t-period < 1) then
                 precip_grid(i,j,t) = sum(copy(i,j,1:t))
               else
                 precip_grid(i,j,t) = sum(copy(i,j,t-period:t))
               endif 
            enddo
         enddo
       enddo
       deallocate(copy)
     end subroutine

     subroutine rolling_average_over_a_period(grid,period)
      !This routine takes hourly 2d variable and will compute the running
      !average. Will assume the first period elements will be 
      !will be padded with 0
      !
      !Period should equal timestep of the hybrid model
      !
      !Example precip_grid = [1,0,0,0,0,0,0,1,1,0,0,0,0,0,0] and period = 6
      !hours
      !This routine will output grid = [1,0.5,0.33,0.25,0.2,0.03125,0,0.166,0.333,0.333,0.333,0.333,0.166,0,0]
      real(kind=dp), intent(inout) :: grid(:,:,:) !Hourly data (x,y,t) that
                                                         !at exit will be 
                                                         !averag over a period

      integer, intent(in) :: period

      !local vars
      real(kind=dp), allocatable :: copy(:,:,:)

      integer :: i, j, t
      integer :: x_len, y_len, t_len

      x_len = size(grid,1)
      y_len = size(grid,2)
      t_len = size(grid,3)

      allocate(copy,source=grid)

      do i=1, x_len
         do j=1, y_len
            do t=1, t_len
               if(t-period < 1) then
                 grid(i,j,t) = sum(copy(i,j,1:t))/t
               else
                 grid(i,j,t) = sum(copy(i,j,t-period:t))/period
               endif
            enddo
         enddo
       enddo
       deallocate(copy)
     end subroutine

     subroutine rolling_average_over_a_period_2d(grid,period)
      !This routine takes hourly 2d variable and will compute the running
      !average. Will assume the first period elements will be
      !will be padded with 0
      !
      !Period should equal timestep of the hybrid model
      !
      !Example precip_grid = [1,0,0,0,0,0,0,1,1,0,0,0,0,0,0] and period = 6
      !hours
      !This routine will output grid =
      ![1,0.5,0.33,0.25,0.2,0.03125,0,0.166,0.333,0.333,0.333,0.333,0.166,0,0]
      real(kind=dp), intent(inout) :: grid(:,:) !Hourly data (x,y,t) that
                                                         !at exit will be
                                                         !averag over a period

      integer, intent(in) :: period

      !local vars
      real(kind=dp), allocatable :: copy(:,:)

      integer :: i, j, t
      integer :: x_len, y_len, t_len

      x_len = size(grid,1)
      t_len = size(grid,2)

      allocate(copy,source=grid)

      do i=1, x_len
         do t=1, t_len
            if(t-period < 1) then
              grid(i,t) = sum(copy(i,1:t))/t
            else
              if(abs(sum(copy(i,t-period:t))) > 0.0000001) then 
                 grid(i,t) = sum(copy(i,t-period:t))/period
              else
                 grid(i,t) = grid(i,t)
              endif 
            endif
         enddo
       enddo
       deallocate(copy)
     end subroutine

     function linear_increase_then_platue(start_val,final_val,current_time,platue_time) result(val)
       !Function that should be used to make a linearly increasing sst bias that
       !starts at start_val and linearly increases to final_val at platue_time
       real(kind=dp), intent(in)  :: start_val
       real(kind=dp), intent(in)  :: final_val
       real(kind=dp), intent(in)  :: current_time
       real(kind=dp), intent(in)  :: platue_time

       real(kind=dp)  :: val  

       if(current_time > platue_time) then 
        val = final_val
       else
        val = (final_val/(platue_time))*(current_time) + start_val
       endif 

       return 
    end function 
end module mod_utilities
