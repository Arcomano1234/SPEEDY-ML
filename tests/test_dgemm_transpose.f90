program main
 use mod_utilities, only : dp
 use blas95 
 implicit none 

 real(kind=dp), allocatable :: states(:,:)
 real(kind=dp), allocatable :: temp(:,:)

 real(kind=dp), parameter :: beta=0.0_dp, alpha = 1.0_dp

 integer :: n,m

 n = 7000+136
 m = n
 allocate(states(n,m))
 allocate(temp(n,m))
 call random_number(states)
 call DGEMM('N','T',n,n,m,alpha,states,n,states,m,beta,temp,n) 

end program main
