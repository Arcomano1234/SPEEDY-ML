module mod_unit_tests
  !We need something to test varies stuff

  use mod_utilities, only : dp

  contains
  !------------Linear Algebra Tests-------------!
  subroutine test_linalg()
     !Main routine to test various linear algebra functions
     logical :: passer
     
     call test_pinv(passer) 

  end subroutine
   
  subroutine test_pinv(passer)
    !Just tests pinv 
    use mod_linalg, only : pinv_svd 

    real(kind=dp) :: A(10,10)
    real(kind=dp) :: Ainv(10,10)
    real(kind=dp) :: realinv(10,10)
    real(kind=dp) :: threshold = 1e-10

    integer       :: i 

    logical, intent(out) :: passer 

    passer = .false.

    A = 0
    Ainv = 0
    realinv = 0
 
    do i=1, 10 
       A(i,i) = i
       realinv(i,i) = 1_dp/i       
    enddo 
    
    call pinv_svd(A,n,n,Ainv)
    
    if(sum(Ainv-realinv) < threshold) then
       passer = .true.
    endif
 
    return        
  end subroutine
   
  !-----------End of Linear Algebra Tests ----!



  !---------Grid and Tiling Tests ------------!
  
  subroutine test_res_domain()
    !Main routine to test various
    logical :: passer

    call test_getxyresextent(passer)

  end subroutine
  
  subroutine test_getxyresextent(passer)
    integer :: numprocs,proc_num
    integer :: localres_xstart,localres_xend,localres_ystart,localres_yend
    integer :: localxchunk,localychunk 
   
    integer :: truth_xstart, truth_xend, truth_ystart, truth_yend
    integer :: truth_xchunk, truth_ychunk

    logical :: passer
   
    passer = .false.

    truth_xchunk = 4
    truth_ychunk = 4
    
    truth_xstart =  49
    truth_xend = 52
    
    truth_ystart = 9
    truth_yend = 12
    
    numprocs = 288
    proc_num = 145
   
    call getxyresextent(numprocs,proc_num,localres_xstart,localres_xend,localres_ystart,localres_yend,localxchunk,localychunk)
    
    if((truth_xchunk==localxchunk).and.(truth_ychunk==localychunk).and.(truth_xstart==localres_xstart).and.(truth_xend==localres_xend).and.(truth_ystart==localres_ystart).and.(truth_yend==localres_yend)) then
      passer = .true.
    else 
      print *, truth_xchunk,localxchunk,truth_ychunk,localychunk,truth_xstart,localres_xstart,truth_xend,localres_xend,truth_ystart,localres_ystart,truth_yend,localres_yend
    endif
    
    return  
  end subroutine
  
  !---------End Grid and Tiling Tests ------------!

  !---------MPI Tests  ------------!
  
  subroutine test_mpi(dat)
  !Main routine to test mpi stuff
  !A higher level test needs real data to work as well
  !as the things like mpi_type, grid type, and res type
  !Call this routine after type initialization and reading in the data

  !Input dat which is a single time step of speedy full grid

    real(kind=dp), intent(in) :: dat(:,:,:,:,:)
  
    logical                   :: passer
  end subroutine

  subroutine test_predictionmpicontroller(dat,passer)
    use mod_reservoir, only : res

    real(kind=dp), intent(in)  :: dat(:,:,:,:,:)
    real(kind=dp), allocatable :: outvec(:,:), feedbackvec(:,:)

    integer                    :: timestep

    logical, intent(out)       :: passer  

    timestep = 1
   
    call predictionmpicontroller(res,timestep,outvec,feedbackvec)
   
  end subroutine 
end module 
