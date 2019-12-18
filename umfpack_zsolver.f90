! UMFPack_zsolver
! A Fortran module for calling the sparse LU solver library
! UMFPack.  Separate subroutines are available for factoring and solving a system
! so that the factors can be reused for multiple right-hand-sides.
!
!  Kerry Key 
!  Scripps Institution of Oceanography
!
! Uses Timothy A. Davis 's umf4_f77zwrapper.c for wrappers to the C functions in UMFpack.
! 
!

module umfpack_zsolver
 
    implicit none
!
! Derived type to store UMFPack variables for a given matrix factorization
!
    type umfpack_z 
        integer, public        :: error
        integer(8), private    :: numeric  ! the factors
        real(8), dimension(20) :: control
        real(8), dimension(90) :: info 
    end type umfpack_z
 
    public :: umfpack_zfactor  ! Step 1. Factor a complex sparse matrix: A = LU.
    public :: umfpack_zsolve   ! Step 2. Solve the linear system: A x = LU x = b.
    public :: umfpack_zfree    ! Step 3. Deallocate the memory for the factors LU.
    
    contains

!------------------------------------------------------------------------- 
!----------------------------------------------------------umfpack_zfactor 
!------------------------------------------------------------------------- 
    subroutine umfpack_zfactor(this,Az,Ai,Ap)
!
! Factors a complex matrix A, with pointers to solution being stored in 
! umfpack derived type variable "this".
!

!
! Arguments:
!
    type(umfpack_z), intent(inout)        :: this
    complex(8), dimension(:), intent(in)  :: Az    
    integer, dimension(:), intent(inout)  :: Ai, Ap     ! modified for c indexing here, corrected before returing

!
! Local variables:
!
    integer    :: nz,n
    integer(8) :: symbolic 
  
!
! Get matrix size
!
    n = ubound(Ap,1) - 1
    nz = Ap(n+1) - 1
!
! Convert to C's zero based indices:
!
    Ai = Ai - 1
    Ap = Ap - 1
!
! Start calling the UMFPack routines:
!
    call umf4zdef (this%control)
    
!
! Tweak the controls for forcing symmetric matrix:
!
! from umfpack.h:
!#define UMFPACK_STRATEGY 5		/* auto, symmetric, or unsym. */
!#define UMFPACK_ORDERING 10             /* ordering method to use */
!/* Control [UMFPACK_STRATEGY] is one of the following: */
!#define UMFPACK_STRATEGY_AUTO 0		/* use sym. or unsym. strategy */
!#define UMFPACK_STRATEGY_UNSYMMETRIC 1	/* COLAMD(A), coletree postorder,
!					   not prefer diag*/
!#define UMFPACK_STRATEGY_OBSOLETE 2     /* 2-by-2 is no longer available */
!#define UMFPACK_STRATEGY_SYMMETRIC 3	/* AMD(A+A'), no coletree postorder,
!					   prefer diagonal */
!#define UMFPACK_ORDERING_CHOLMOD 0      /* use CHOLMOD (AMD/COLAMD then METIS)*/
!#define UMFPACK_ORDERING_AMD 1          /* use AMD/COLAMD */
!#define UMFPACK_ORDERING_GIVEN 2        /* user-provided Qinit */
!#define UMFPACK_ORDERING_METIS 3        /* use METIS */
!#define UMFPACK_ORDERING_BEST 4         /* try many orderings, pick best */
!#define UMFPACK_ORDERING_NONE 5         /* natural ordering */
!#define UMFPACK_ORDERING_USER 6         /* user-provided function */					   
!    this%control (5+1)  =  3 ! SYMMETRIC strategy
!    this%control (10+1) =  1 ! AMD/COLAMD fast 
!
! Print UMFPack controls:
!
     this%control (1) = 1  
!
! Symbolic analysis:
!
    call umf4zsym (n, n, Ap, Ai, Az, , symbolic, this%control, this%info)

    if (this%info (1) < 0) then
       this%error = this%info (1) 
       call umf4zfsym (symbolic)
       return
    end if
!
! Numeric factorization:
!
    call umf4znum (Ap, Ai, Az, , symbolic, this%numeric, this%control,this%info)

    if (this%info (1) < 0) then
       this%error = this%info (1) 
       call umf4zfsym (symbolic)
       return
    end if
    
!    call umf4zpinf (this%control, this%info)
!    call umf4zpcon (this%control)

!
! Convert back to Fortran's one based indices:
!
    Ai = Ai + 1
    Ap = Ap + 1
!
! Free symbolic analysis memory:
!
    call umf4zfsym (symbolic)

    end subroutine umfpack_zfactor

!------------------------------------------------------------------------- 
!-----------------------------------------------------------umfpack_zsolve
!------------------------------------------------------------------------- 
    subroutine umfpack_zsolve (this,b,x)
    implicit none
    
!
! Arguments:
!
    type(umfpack_z), intent(inout)        :: this
    complex(8), dimension(:), intent(in)  :: b   
    complex(8), dimension(:), intent(out) :: x 
!
! Local variables:
!
    integer(4) :: sys
!
! Solve the linear system Ax=b:
!
    sys = 0 ! Ax = b

   call umf4zsol (sys, x,  , b,   , this%numeric, this%control, this%info)
   ! Note the calling sequence uses the 'packed' complex format described in the 
   ! UMFPack UserGuide.pdf

    if (this%info (1) < 0) then
       this%error = 1
       return
    end if
    
    end subroutine umfpack_zsolve    
    
!------------------------------------------------------------------------- 
!------------------------------------------------------------umfpack_zfree
!------------------------------------------------------------------------- 
    subroutine umfpack_zfree(this)
    
    implicit none
        
    type(umfpack_z), intent(inout)       :: this
!
! Free numeric factorization memory:
!
    call umf4zfnum (this%numeric)

    end subroutine umfpack_zfree


end module umfpack_zsolver
 
