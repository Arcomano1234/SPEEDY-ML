module mod_linalg

  use MKL_SPBLAS
  use mod_utilities, only : dp, main_type, reservoir_type
  
  implicit none

  contains
  
      subroutine mklsparse(reservoir)
         !routine to make a MKL sparse array

         type(reservoir_type), intent(inout) :: reservoir

         integer :: stat

         stat = mkl_sparse_d_create_coo(reservoir%cooA, SPARSE_INDEX_BASE_ONE, reservoir%n, reservoir%n, reservoir%k, reservoir%rows, reservoir%cols, reservoir%vals)
         if(stat.ne.0) then
            print *, 'MKL sparse creation failed because of stat error',stat,'exiting'
            print *, 'res%n,res%k',reservoir%n,reservoir%k
            stop
         endif

         reservoir%descrA%TYPE = SPARSE_MATRIX_TYPE_GENERAL
     end subroutine
    
     subroutine pinv_svd(A,m,n,Ainv)
       !A Moore-Penrose pseudo-inverse using single value decomposition
 
       external DLANGE
       real(kind=dp) :: DLANGE
 
       integer,intent(in) :: m, n
       real(kind=dp),intent(in) :: A(m,n)
       real(kind=dp),intent(out) :: Ainv(n,m)
       real(kind=dp),dimension(m,n):: A1, A2, sigma, temp
       integer :: i, j, k, l, lwork, info
       real(kind=dp), allocatable, dimension(:,:) :: U
       real(kind=dp), allocatable, dimension(:,:) :: VT
       real(kind=dp), dimension(N,N) :: BUFF
       real(kind=dp), allocatable, dimension(:) :: WORK
       real(kind=dp), allocatable, dimension(:) :: RWORK
       real(kind=dp), allocatable, dimension(:) :: S
       real(kind=dp),parameter :: thres=1e-2 !TODO maybe change this 
       real(kind=dp) :: normA, normAPA, normPAP, tester
       logical, parameter :: verbose= .False.
 
       k=min(m,n)
       l=max(m,n)
       lwork = max(1,3*k+l,5*k)*4
 
       allocate(U(m,k))
       allocate(VT(k,n))
       allocate(work(lwork))
       allocate(S(k))
 
       !This is needed because LAPACK likes destroying input matrices
       A1 = A
       A2 = A
 
       ! Compute the SVD of A1
       call DGESVD('S', 'S', M, N, A1, M, S, U, M, VT, K, WORK, LWORK,INFO) !Major problem was fixed we didnt need RWORK and led to memory problems 
       if(INFO /= 0) then 
         print *, 'PINV DGESVD has a problem error is info=',INFO
         
       endif 

       !  Compute PINV = VT**T * SIGMA * U**T in two steps
       do j = 1, K
          tester = s(j)
          !Need this because if tester is smaller than threshold makes S(j) zero
          !Because when we take the reciprical we dont want to get too large of a
          !number. See numerical linear algebra book on single value decomposition
          !and svd to get an pseudo-inverse

  
         if(tester.gt.thres) then
            call DSCAL( M, real(1 / S( j ),kind=dp), U( 1, j ), 1 )
         else
            call DSCAL( M, 0.0_dp, U( 1, j ), 1 )
         endif
      end do

      call DGEMM( 'C', 'C', N, M, K, real(1.0,kind=dp), VT, K, U, M, real(0.0,kind=dp), Ainv, N)

      !If you really don't believe me that the above code works you can
      !check by making verbose true
      if(verbose) then
          !check if the diagonals of temp are one
          temp = matmul(A,Ainv)
          print *,'max of what should be an identity maxtrix if not close to one then problem', maxval(temp)
          !  check the result
          normA = DLANGE( 'F', M, N, A2, M, lwork )
          call DGEMM( 'N', 'N', N, N, M, real(1.0), Ainv, N, A2, M, real(0.0), BUFF, N )

          call DGEMM( 'N', 'N', M, N, N, real(-1.0), A2, M, BUFF, N, real(1.0), A2, M )
          normAPA = DLANGE( 'F', M, N, A2, M, lwork )

          call DGEMM( 'N', 'N', N, M, N, real(-1.0), BUFF, N, Ainv, N, real(1.0), Ainv, N );
          normPAP = DLANGE( 'F', N, M, Ainv, N, lwork )

          write(*,"(A, e11.4)") '|| A - A*P*A || = ', normAPA/normA
          write(*,"(A, e11.4)") '|| P - P*A*P || = ', normPAP/normA
       endif

       return
     end subroutine

     subroutine mldivide(A,B)
       !What should be the same as matlab mldivide also known as \
       !Solves A*x = B given A and B
       !A is n-by-n and B 
       !B becomes X at the exit if info == 0 
  
        real(kind=dp), intent(inout)  :: A(:,:), B(:,:)


        !Local parameters 
        real(kind=dp), allocatable :: temp2d(:,:)

        integer                    :: n, m, k, l

        !lapack stuff
        integer              :: nrhs, lda, ldb, info
        integer, allocatable :: ipiv(:)
 

        !check if inputs are the right shape
        n = size(A,1)
        m = size(A,2) 
        l = size(B,1)
        k = size(B,2)
  
        if(n /= l) then
           print *, 'Column of A is not the same size of column of B. Cant compute solution returning A and B unchanged'
           return
        endif 

        nrhs = k
        lda = max(1,n)
        ldb = max(1,n)

        allocate(ipiv(n))

        call dgesv(n, nrhs, A, lda, ipiv, B, ldb, info ) 
       
        if(info /= 0) then
          print *, 'something went wrong with dgesv info = ',info
          print *, 'B is not the solution'
        endif  
     end subroutine

     subroutine eigval(A,x,y,maxeig)
       !Outdated routine to get eigenvalues of a 
       !dense array
       integer, intent(in) :: x,y
       real, intent(in),dimension(x,y) :: A
       integer :: ldVL, ldVR, lworker, ierr, info, i, ldA
       real, allocatable :: work(:), VR(:,:), VL(:,:)
       real, dimension(x) :: eig, wr, wi
       character        :: jobVL, jobVR
       real :: maxeig
       external SGEEV

       print *,'test'
       jobVL = 'N' ! The left eigenvector u(j) of A satisfies: u(j)**H * A = lambda(j) * u(j)**H. 'N' to not compute.
       jobVR = 'N' ! The right eigenvector v(j) of A satisfies: A * v(j) = lambda(j) * v(j). 'V' to compute.
       ldA = x; ldVL = 1; ldVR = 1
       lworker = max(1,5*x)
       allocate(work(lworker))
       call SGEEV(jobVL,jobVR,x,A,ldA,WR,WI,VL,ldVL,VR,ldVR,work,lworker,info) !TODO this is not working
  
       do i=1,x
          eig(i) = (wr(i)**2+wi(i)**2)**0.5
       enddo
       maxeig = maxval(eig)
       return
     end subroutine

     subroutine makesparse(reservoir)
       !This subroutine 100% works for making a random
       !sparse matrix
       !Hand checked this
       use mod_utilities, only : shuffle 

       type(reservoir_type), intent(inout) :: reservoir

       integer ::  counter, leftover, i

       !Get random vals
       call RANDOM_NUMBER(reservoir%vals)

       !This block makes a random choice with no repeat (for small M >100 or
       !very dense sparse matrix technically there could be a repeat
       !row/column pair but its low chance) and doesnt affect anything.
       !The random choice uses a kshuffle to make the random choice
       if(reservoir%k.gt.reservoir%n) then
          counter = floor(real(reservoir%k/reservoir%n))
          leftover = mod(reservoir%k,reservoir%n)
          do i = 1,counter
             call shuffle(reservoir%n,reservoir%n,reservoir%rows((i-1)*reservoir%n+1:i*reservoir%n))
             call shuffle(reservoir%n,reservoir%n,reservoir%cols((i-1)*reservoir%n+1:i*reservoir%n))
          enddo
     
         if(leftover.ne.0) then
            call shuffle(reservoir%n,leftover,reservoir%rows((i-1)*reservoir%n+1:reservoir%k))
            call shuffle(reservoir%n,leftover,reservoir%cols((i-1)*reservoir%n+1:reservoir%k))
         endif
       else
         call shuffle(reservoir%n,reservoir%k,reservoir%rows)

         call shuffle(reservoir%n,reservoir%k,reservoir%cols)
       endif

       call mklsparse(reservoir)
       return

     end subroutine
     
     subroutine sparse_eigen(reservoir,maxn,k,eigs)
        !Horrible to follow routine to get the max eigenvalue of a sparse array 
        !Do not touch this 
        !---------------------------------------------------------------------------
        !
        !     %-----------------------------%
        !     | Define maximum dimensions   |
        !     | for all arrays.             |
        !     | MAXN:   Maximum dimension   |
        !     |         of the A allowed.   |
        !     | MAXNEV: Maximum NEV allowed |
        !     | MAXNCV: Maximum NCV allowed |
        !     %-----------------------------%
        !
        type(reservoir_type) :: reservoir
        integer, intent(in) :: k, maxn

        real(kind=dp), intent(out) :: eigs
        !     %--------------%
        !     | Local Arrays |
        !     %--------------%
        !
        integer           maxnev, maxncv, ldv
        parameter (maxncv=30)
        integer           iparam(11), ipntr(14)
        logical           select(maxncv)
        real(kind=dp) ::  ax(maxn), d(maxncv,3), resid(maxn),v(maxn,maxncv), workd(3*maxn),workev(3*maxncv),workl(3*maxncv*maxncv+6*maxncv)
        !
        !     %---------------%
        !     | Local Scalars |
        !     %---------------%
        !
        character         bmat*1, which*2
        integer           ido, n, nx, nev, ncv, lworkl, info, j, ierr, nconv, maxitr, ishfts, mode
        real(kind=dp) ::  tol, sigmar, sigmai
        logical           first, rvec
        !
        !     %------------%
        !     | Parameters |
        !     %------------%
        !
        real(kind=dp) ::  zero
        parameter         (zero = 0.0D+0)
        !
        !     %-----------------------------%
        !     | BLAS & LAPACK routines used |
        !     %-----------------------------%
        real(kind=dp) ::  dlapy2, dnrm2
        external          dlapy2, dnrm2, daxpy
        !
        !     %--------------------%
        !     | Intrinsic function |
        !     %--------------------%
        !
        intrinsic         abs
        !
        !     %-----------------------%
        !     | Executable Statements |
        !     %-----------------------%
        !     |                   N <= MAXN                      |
        !     |                 NEV <= MAXNEV                    |
        !     |           NEV + 2 <= NCV <= MAXNCV               |
        !     %--------------------------------------------------%
        !
      
        maxnev=k
        ldv=maxn
        nx    = reservoir%n
        n     = nx
        nev   = 4
        ncv   = 20
        if ( n .gt. maxn ) then
           print *, ' ERROR with _NDRV1: N is greater than MAXN '
           go to 9000
        else if ( nev .gt. maxnev ) then
           print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
           go to 9000
        else if ( ncv .gt. maxncv ) then
           print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
           go to 9000
        end if
        bmat  = 'I'
        which = 'LM'
        !
        !     %-----------------------------------------------------%
        !     | The work array WORKL is used in DNAUPD as           |
        !     | workspace.  Its dimension LWORKL is set as          |
        !     | illustrated below.  The parameter TOL determines    |
        !     | the stopping criterion. If TOL<=0, machine          |
        !     | precision is used.  The variable IDO is used for    |
        !     | reverse communication, and is initially set to 0.   |
        !     | Setting INFO=0 indicates that a random vector is    |
        !     | generated in DNAUPD to start the Arnoldi iteration. |
        !     %-----------------------------------------------------%
        !
        lworkl  = 3*ncv**2+6*ncv
        tol    = zero
        ido    = 0
        info   = 0
        !
        !     %---------------------------------------------------%
        !     | This program uses exact shifts with respect to    |
        !     | the current Hessenberg matrix (IPARAM(1) = 1).    |
        !     | IPARAM(3) specifies the maximum number of Arnoldi |
        !     | iterations allowed.  Mode 1 of DNAUPD is used     |
        !     | (IPARAM(7) = 1). All these options can be changed |
        !     | by the user. For details see the documentation in |
        !     | DNAUPD.                                           |
        !     %---------------------------------------------------%
        !
        ishfts = 1
        maxitr = 300
        mode   = 1
        !
        iparam(1) = ishfts
        iparam(3) = maxitr
        iparam(7) = mode
        !
        !     %-------------------------------------------%
        !     | M A I N   L O O P (Reverse communication) |
        !     %-------------------------------------------%
        !
        10   continue
        !
        !        %---------------------------------------------%
        !        | Repeatedly call the routine DNAUPD and take |
        !        | actions indicated by parameter IDO until    |
        !        | either convergence is indicated or maxitr   |
        !        | has been exceeded.                          |
        !        %---------------------------------------------%
        !
             call dnaupd ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
             !
             if (ido .eq. -1 .or. ido .eq. 1) then
             !
             !           %-------------------------------------------%
             !           | Perform matrix vector multiplication      |
             !           |                y <--- OP*x                |
             !           | The user should supply his/her own        |
             !           | matrix vector multiplication routine here |
             !           | that takes workd(ipntr(1)) as the input   |
             !           | vector, and return the matrix vector      |
             !           | product to workd(ipntr(2)).               |
             !           %-------------------------------------------%
             !
               call smatrix_vector(reservoir, nx, workd(ipntr(1)), workd(ipntr(2)))
             !           %-----------------------------------------%
             !           | L O O P   B A C K to call DNAUPD again. |
             !           %-----------------------------------------%
             !
               go to 10
             !
              end if
             !
             !     %----------------------------------------%
             !     | Either we have convergence or there is |
             !     | an error.                              |
             !     %----------------------------------------%
             !
             if ( info .lt. 0 ) then
             !
             !        %--------------------------%
             !        | Error message, check the |
             !        | documentation in DNAUPD. |
             !        %--------------------------%
             !
               print *, ' '
               print *, ' Error with _naupd, info = ', info
               print *, ' Check the documentation of _naupd'
               print *, ' '
            !
            else
            !
            !        %-------------------------------------------%
            !        | No fatal errors occurred.                 |
            !        | Post-Process using DNEUPD.                |
            !        |                                           |
            !        | Computed eigenvalues may be extracted.    |
            !        |                                           |
            !        | Eigenvectors may also be computed now if  |
            !        | desired.  (indicated by rvec = .true.)    |
            !        %-------------------------------------------%
            !
               rvec = .true.
            !
               call dneupd(rvec, 'A', select, d, d(1,2), v, ldv, sigmar, sigmai, workev, bmat, n, which, nev,tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr )
            !
            !        %-----------------------------------------------%
            !        | The real part of the eigenvalue i returned   |
            !        | in the first column of the two dimensional   |
            !        | array D, and the imaginary part is returned  |
            !        | in the second column of D. The corresponding |
            !        | eigenvectors are returned in the first NEV   |
            !        | columns of the two dimensional array V if    |
            !        | requested.  Otherwise, an orthogonal  basis  |
            !        | for the invariant subspace corresponding to  |
            !        | the eigenvalues in D is returned in V.       |
            !        %-----------------------------------------------%
            !
               if ( ierr .ne. 0) then
                        
            !           %------------------------------------%
            !           | Error condition:                  |
            !           | Check the documentation of DNEUPD.|
            !           %------------------------------------%
            !
                   print *, ' '
                   print *, ' Error with _neupd, info = ', ierr
                   print *, ' Check the documentation of _neupd. '
                   print *, ' '
            !
               else
!
                  first  = .true.
                  nconv  = iparam(5)
                  do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                     if (d(j,2) .eq. zero)  then
!
!                  %--------------------%
!                  | Ritz value is real |
!                  %--------------------%
!
                       call smatrix_vector(reservoir, nx, v(1,j), ax)
                       call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                       d(j,3) = dnrm2(n, ax, 1)
                       d(j,3) = d(j,3) / abs(d(j,1))
!
                     else if (first) then
!
!                  %------------------------%
!                  | Ritz value is complex. |
!                  | Residual of one Ritz   |
!                  | value of the conjugate |
!                  | pair is computed.      |
!                  %------------------------%
                     !
                       call smatrix_vector(reservoir, nx, v(1,j), ax)
                       call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                       call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                       d(j,3) = dnrm2(n, ax, 1)
                       call smatrix_vector(reservoir, nx, v(1,j+1), ax)
                       call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                       call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                       d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                       d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                       d(j+1,3) = d(j,3)
                       first = .false.
                     else
                       first = .true.
                     end if
!
 20          continue
!
             end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
            if ( info .eq. 1) then
                print *, ' '
                print *, ' Maximum number of iterations reached.'
                print *, ' '
            else if ( info .eq. 3) then
                print *, ' '
                print *, ' No shifts could be applied during implicit',' Arnoldi update, try increasing NCV.'
                print *, ' '
            end if
!
        end if
!
!     %---------------------------%
!     | Done with program dndrv1. |
!     %---------------------------%
!
 9000 continue
!
       eigs = maxval(d)
       return
     
       end subroutine

       subroutine smatrix_vector(reservoir, col, x, y)

         type(reservoir_type), intent(in) :: reservoir

         integer, intent(in)              :: col
         real(kind=dp), intent(inout)     ::  x(col), y(col)

         real(kind=dp) :: alpha, beta
         integer       :: info
         alpha = 1.0
         beta = 0.0

         info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,reservoir%cooA,reservoir%descrA,x,beta,y)
         
         return
       end subroutine
end module 
