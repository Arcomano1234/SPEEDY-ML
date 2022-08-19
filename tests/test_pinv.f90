program main
  use mod_linalg
  use mod_utilities, only : dp

  implicit none

  real(kind=dp), allocatable :: a(:,:), b(:,:), a_trans(:,:), b_trans(:,:)

  real(kind=dp), allocatable :: invstates(:,:)
  integer :: n,m,workernuminputs

  n=5000
  m=5000

  workernuminputs = 640

  allocate(a(n+workernuminputs,m+workernuminputs))
  allocate(b(400,n+workernuminputs))
 
  allocate(a_trans(n+workernuminputs,m+workernuminputs))
  allocate(b_trans(n+workernuminputs,400))
  call random_number(a)
  call random_number(b)

  a_trans = transpose(a)
  b_trans = transpose(b)

  call mldivide(a_trans,b_trans)
  !call pinv_svd(a,n+workernuminputs,n+workernuminputs,invstates) 
  print *, shape(b_trans), shape(transpose(b_trans))
  !print *, invstates(300:310,800:810)
end program 
