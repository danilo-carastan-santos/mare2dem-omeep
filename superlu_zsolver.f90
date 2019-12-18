! SuperLU_zsolver
! A Fortran module for calling the sparse LU solver library SuperLU.
! Separate subroutines are available for factoring and solving a system
! so that the factors can be reused for multiple right-hand-sides.
!
!  Kerry Key 
!  Scripps Institution of Oceanography
 
!

! DGM 2/23/2012 For Windows compilation
#IF DEFINED(_WIN32) .OR. DEFINED(_WIN64)
#DEFINE c_fortran_zgssv c_fortran_zgssv_
#ENDIF

module superlu_zsolver
 
    implicit none
!
! Derived type to store superlu variables for a given matrix factorization
!
    type superlu_z 
        integer, public		:: error        
        integer(8), private	:: factors  ! ptr to the factors
    end type superlu_z
 
    public :: superlu_zfactor  ! Step 1. Factor a complex sparse matrix: A = LU.
    public :: superlu_zsolve   ! Step 2. Solve the linear system: A x = LU x = b.
    public :: superlu_zfree    ! Step 3. Deallocate the memory for the factors LU.
    
    contains

!------------------------------------------------------------------------- 
!----------------------------------------------------------superlu_zfactor 
!------------------------------------------------------------------------- 
    subroutine superlu_zfactor(this,Az,Ai,Ap)
!
! Factors a complex matrix A, with pointers to solution being stored in 
! superlu derived type variable "this".
!

!
! Arguments:
!
    type(superlu_z), intent(inout)        :: this
    complex(8), dimension(:), intent(in)  :: Az    
    integer, dimension(:), intent(inout)  :: Ai, Ap     ! modified for c indexing here, corrected before returing

!
! Local variables:
!
    integer    :: nz, info
    integer(4) :: n

!
! Get matrix size
!
    n = ubound(Ap,1) - 1
    nz = Ap(n+1) - 1

!
! Call superlu c-fortran interface to do factorization:
!
    call c_fortran_zgssv( 1, n, nz,  , Az, Ai, Ap,   ,  , this%factors, info )
    
    if (info .eq. 0) then
        this%error = 0
    else
       this%error = 1
    endif
 
    end subroutine superlu_zfactor

!------------------------------------------------------------------------- 
!-----------------------------------------------------------superlu_zsolve
!------------------------------------------------------------------------- 
    subroutine superlu_zsolve (this,b)
    implicit none
    
!
! Arguments:
!
    type(superlu_z), intent(inout)        :: this
    complex(8), dimension(:,:), intent(inout)  :: b   

!
! Local variables:
!
    integer :: n, nrhs, info
  
!
! Solve the linear system Ax=b:
!
    n = ubound(b,1)
    
    nrhs = ubound(b,2)
 
    call c_fortran_zgssv( 2, n, n , nrhs,  ,  ,  ,  b, n, this%factors, info )
    
    if (info .eq. 0) then
        this%error = 0
    else
       this%error = 1
    endif
    
    end subroutine superlu_zsolve    
    
!------------------------------------------------------------------------- 
!------------------------------------------------------------superlu_zfree
!------------------------------------------------------------------------- 
    subroutine superlu_zfree(this)
    
    implicit none
        
    type(superlu_z), intent(inout)       :: this
    
    integer :: info
!
! Free numeric factorization memory:
!
    call c_fortran_zgssv( 3,  ,   ,  ,  ,  ,  ,   ,  , this%factors, info )
    
    if (info .eq. 0) then
        this%error = 0
    else
       this%error = 1
    endif

    end subroutine superlu_zfree


end module superlu_zsolver
 
