!-----------------------------------------------------------------------
!
!    Copyright 2008-2013
!    Kerry Key
!    Scripps Institution of Oceanography
!    kkey@ucsd.edu
!
!    This file is part of MARE2DEM.
!
!    MARE2DEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MARE2DEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MARE2DEM.  If not, see <http://www.gnu.org/licenses/>.
!
!-----------------------------------------------------------------------
!
! This file contains the 2.5D EM routines that compute the responses
! at a given wavenumber (kx) and optionally a finite element error estimate.
!
! Note that this is a serial routine that does not have any MPI constructs
!
! Fall 2013 Updates:
!   - tweaks to error estimator routines
!   - added option for fast diagonal approximation of bump system
!   - local apriori mesh refinement options
!
! Fall 2012 Updates:
!   - replaced all gen_lhs routines with versions that make better use of mesh topology information. 
! 
! Fall 2011 Updates:
!   - removed QMR routines since we never use them anymore
!   - added linearSolver as required input parameter ('umfpack' or 'superlu' supported)
!   - added lprintDebug_em2dkx as optional input parameter
!   - added option for creating the sparse unstructured LHS matrices using a tree of pointers or the older fixed matrix approach
!   - cleaned up the main routine for improved readability
!   - minimum inner angle constraint for Triangle.c is now a required input parameter
!
! Spring 2011 Updates:
!  - added triaxial anisotropy
!  - now supports UMFPACK sparse solver in addition to SuperLU
!  - solvers now called via external f90 modules.
!  - cleaned up stiffness matrix comps
!  - cleaned up many of the gen_rhs routines
!  - gen_1dmt subroutine was generalized 
!
!==================================================================================================================================! 
!======================================================================================================================== em2dkx_mod
!==================================================================================================================================!     
    module em2dkx_mod 
!
! Note that you need to allocate the arrays below before calling em2dkx 
! and have to deallocate them after calling.  Also note that variables listed 
! below have global scope withing this module.
!
!

    use EM_constants
        
    use fem2d_utilities  ! for generic FE operations
    use MT1D             ! for the 1D boundary conditions for 2D MT computations
    use kdtree2_module   ! for fast searches
    use triangle_mesh    ! for trimesh structure 
    
    use superlu_zsolver  ! for superlu sparse solver
    use umfpack_zsolver  ! for umfpack sparse solver
    
    use sort             ! for quicksort routines
    
    implicit none        ! I declare that you must declare everything
    
    private  ! The default is that all variables, subroutines and functions are private
    

!-----------------------------------------------------------------------
! Public interface:
!-----------------------------------------------------------------------


!--------  
! Inputs:
!--------
    integer,public :: myID  

!
! Conductivity parameters:  The code uses a lookup table for both fixed and free conductivity parameters.
!
! For cAnisotropy:  
!
!  'isotropic':      sigParams(i)         is sig  for mesh%attr = i
!  'triaxial':       sigParams(3*(i-1)+1) is sigx for mesh%attr = i
!                    sigParams(3*(i-1)+2) is sigy for mesh%attr = i
!                    sigParams(3*(i-1)+3) is sigz for mesh%attr = i   
! 'tix','tiy','tiz': 
!                    sigParams(2*(i-1)+1) is sig in symmetry axis (sigx for tix, sigy for ti sigz for tiz) for mesh%attr = i, 
!                    sigParams(2*(i-1)+2) is sig in transverse plane (sigyz for tix,sigyz for tiy, sigxy for tiz) for mesh%attr = i
!
!            
    character(9), public                        :: cAnisotropy  = 'isotropic'  !  ('isotropic','triaxial','tix','tiy','tiz') or '' blank 
    integer, public                             :: nSigParams = 0    ! number of conductivities in sigParams (=nRegions*nRhoPerRegion)
    real(8), dimension(:), allocatable, public  :: sigParams         ! array of conductivities for for fixed and free parameters    
    integer, dimension(:), allocatable, public  :: iFreeParam        ! Integer of free parameter # for each value in sigParams

    
!
! Transmitter parameters:
!
    real(8), public :: kx   ! 2.5D Fourier transformation wavenumber, ignored for MT computations
    real(8), public :: w    ! Angular frequency 
    integer, public :: nTx  
    logical, public :: lMT  ! True if MT computation, else CSEM
                                        
    real(8), dimension(:), allocatable, public :: xTx, yTx, zTx ! x,y,z location of the transmitter, ignored for MT computations

    real(8), dimension(:), allocatable, public :: jx,jy,jz  ! Electric transmitter  dipole moments for each nTx.  
    ! You can easily mess this up, so be careful!  If jx is non-zero, then jy and jz have to be zero. Likewise,
    ! if jy and jz are non-zero, then jx has to be zero. This has to do with the symmetry used in the 2.5D FT.  
    ! kwk debug: put in routine to check that jx and jy,jz are sensible...
    
    real(8), dimension(:), allocatable, public :: mx,my,mz  ! Magnetic transmitter dipole moments for each nTx.
    ! these need to be compatible with jx,jy,jz. ie, if j/=0 then m = 0 and vice versa. You can only have either electric or 
    ! magnetic in a given line, not both. Also, if mx is non-zero, then my and mz have to be zero. Likewise,
    ! if my and mz are non-zero, then mx has to be zero. This has to do with the symmetry used in the 2.5D FT.  

!
! Receiver parameters:
!
    integer, public :: nRx                                     ! Number of receivers       
    real(8), dimension(:), allocatable, public :: xRx,yRx,zRx  ! Receiver locations

!
! Data masks. Set to true if EM fields to be computed at particular Rx-Tx pair.
!
    integer, dimension(:,:), allocatable, public :: iDataMask ! 0 = no data, 1 = no data but include for refinement, 2=data present

!
! Mesh parameters are stored in the derived type 'trimesh' for ease of use
!
    type (trimesh), public :: mesh, newmesh    
    character(256), public :: fileroot         ! base name for mesh files 
    integer, public        :: meshnumber       ! keeps track of mesh refinement number (1 = starting mesh)
    
!
! Local a priori refinement parameters:
!
    logical, public :: lLocalRefine  = .true.        ! Set to true to refine grid on input using skin depth rules
    real(8), public :: nTxEleArea = 1d0, nRxEleArea = 1d2, nSkinDepthsRegional = 1d0/2d0 ! Tolerances for a priori mesh refinement around Tx's and Rx's
                                                                                                              
!
!  Adjustable parameters controlling the adaptive refinement and error estimation:
!
    real(8), public :: errortolerance    = 1.0       ! tolerance for relative error
    real(8), public :: minRangeProtector = 10        ! (m) protects agains over refinement close to the Tx
    real(8), public :: rxCloud           = 0         ! 
    real(8), public :: minArea           = 1d-2      ! don't let refinement get carried away to make really small elements  

    integer, public :: max_nsubrefine    = 3         ! max # of projection subrefinements for a given error estimator    
    real(8), public :: pct_refine        = 5.        ! percent of worst elements to refine per step (each of max_nsubrefine)    
    real(8), public :: ecutoff           = 1.d-18    ! E(kx) error floor in V/m/Am  
    real(8), public :: hcutoff           = 1.d-12    ! H(kx) error floor A/m/Am, note  B = 1.26x10^-6 * H
    integer, public :: maxMeshNodes      = 75000     ! max# nodes in mesh, stop refinement if mesh grows this big. This will keep 
                                                     ! the code from running out of memory, but repone accuracy may be compromised                        
!
! Model response derivatives:
!    
    logical, public :: lCompDerivs ! set =.true. to output the field derivatives with respect to conductivity
    ! Sensitivities, only output if lCompDerivs = .true. and the component flags are true.
    ! set these to false on input to skip their computations (good for saving memory for large inversion data sets).
    ! lComps(1:6) = [lex ley lez lhx lhy lhz]
    logical, dimension(6), public  :: lComps = .true.           

!
! Other adjustable parameters:
!
    real(8), public         :: minqangle          = 25.    ! quality angle for triangle.c (minimun inner angle for triangulation)
    character(32), public   :: linearSolver       = 'umfpack'    ! linear solver for FE systems: 'superlu' or 'umfpack'
    logical, public         :: lprintDebug_em2dkx = .false.   
    logical, public         :: lUseBumpFields     = .true.
    logical, public         :: lSaveLinearSystem  = .false.  ! write out linear system Ax=b to files A = *.lhs, b=*.rhs, x=*.exhx
    logical, public         :: lSaveMeshFiles     = .false. ! set to .true. to write out the mesh to .poly, .ele and .node files   
    logical, public         :: lComputeErrorEst             ! Set to true to compute error estimates for each element
    integer, public         :: idual_func         = 2       ! 0,1,2 for G0,G1,G2 error functionals from our paper
    logical, public         :: lBumpDiagApprox    = .true.  ! only use diagonal of bump matrix, so no factorization needed...

!-----------
! Outputs:
! ---------
 
! 
! EM fields in kx domain:
!
! *** Note these are single precision since there's no point in outputting double precision since the FE results will never
!     be that accurate (but yet we NEED to use DP within the FE method). SP saves memory and reduces the MPI load for passing these
!     back to the master node.  
!
    
    complex(8), dimension(:,:), allocatable, public  ::  ex_kx,ey_kx,ez_kx,hx_kx,hy_kx,hz_kx    ! nRx x nTx size arrays                                          
    
    type :: kxdsig
         complex(4),dimension(:), allocatable   :: dsig      ! (iparam)  ! this is single precision to reduce memory footprint
    end type kxdsig   
    type(kxdsig),dimension(:,:), allocatable, public :: dex_kx,dey_kx,dez_kx,dhx_kx,dhy_kx,dhz_kx ! (isite,iTx)    
    
!
! Element error indicator:
!
    real(8), public :: maxerr ! maximum error at any site       
    logical, public :: lmadeNewMesh     ! true if a new mesh has been made
        
!
! Public subroutines:
!        
    public :: em2dkx, localRefinement, deallocate_em2dkx

   
!-----------------------------------------------------------------------
! Private variables:
!-----------------------------------------------------------------------
!
! These are allocated and handled only inside em2dkx and the private attribute 
! means you can't access these from outside the em2dkx module.
!
!
! A few adjustable parameters (don't change these unless you know what you are doing): 
!
    logical, private            :: lusedualerror  = .true.    ! for debugging purposes only, don't change this

!
! Variables and arrays used across multiple subroutines:
!
    real(8), private                                 :: kx2
    complex(8), private                              :: ikx, ommu, iomeps 
    integer, dimension(:), allocatable, private      :: eRx, eTx
    logical, dimension(:), allocatable, private      :: lCloud
    integer, private                                 :: nedges 
    integer, dimension(:,:), allocatable, private    :: edges
    integer, dimension(:), allocatable, private      :: node2tri

!
! Sparse LHS and RHS arrays for linear FE's
!
    complex(8), dimension(:), allocatable, private   :: val
    integer,    dimension(:), allocatable, private   :: col, irw
    integer, private                                 :: nnz, info
    complex(8), dimension(:,:), allocatable, private :: rhs
    
    type(superlu_z) :: superlu_p  ! Derived type that holds pointers to superlu factorization  for primary/dual system  
    type(umfpack_z) :: umfpack_p  ! Derived type that holds pointers to umfpack factorization  for primary/dual system         
!
! Arrays for uncoupled E and H MT solver:
!
    complex(8), dimension(:), allocatable, private   :: val_lE
    integer,    dimension(:), allocatable, private   :: col_lE, irw_lE
    integer, private                                 :: nnz_lE, info_lE 
    
    complex(8), dimension(:), allocatable, private   :: val_lH
    integer,    dimension(:), allocatable, private   :: col_lH, irw_lH 
    integer, private                                 :: nnz_lH, info_lH 
    
    complex(8), dimension(:,:), allocatable, private :: rhs_lE, rhs_lH
    
    type(superlu_z) :: superlu_p_E, superlu_p_H  ! Derived type that holds pointers to superlu factorization for primary/dual system
    type(umfpack_z) :: umfpack_p_E, umfpack_p_H  ! Derived type that holds pointers to umfpack factorization for primary/dual system
         
!
! Sparse LHS and RHS arrays for quadratic bump FE's
!
! Coupled anisotropic:
    complex(8), dimension(:), allocatable, private   :: val_q 
    integer,    dimension(:), allocatable, private   :: col_q, irw_q 
    integer, private                                 :: nnz_q, info_q

    complex(8), dimension(:,:), allocatable, private :: rhs_q 
    complex(8), dimension(:,:), allocatable, private :: rhs_wh 
    

    type(superlu_z) :: superlu_q  ! Derived type that holds pointers to superlu factorization for bump system  
    type(umfpack_z) :: umfpack_q  ! Derived type that holds pointers to umfpack factorization for bump system       
    

    integer, dimension(:), allocatable, private      :: nodlab,bcnod   
    real(8), dimension(:), allocatable, private      :: errnrm          
    character(32) :: retricommand

!
! Kd-tree pointer for fast element searches:
!   
    type(kdtree2), pointer, private  :: tree    
    
!
! Timers (only needed for debugging and performance optimizations)
!
    logical, private, parameter :: lprinttimers = .false. 
    real(8), private :: t0,t1,t2,t3,t4,t5,t6,t7,t7b,t8,t9
    real(8), private :: t10,t11,t12,t13,t14,t15,t16,t17,t18,t19, t20    
 
    
    contains

!==================================================================================================================================! 
!================================================================================================================= deallocate_em2dkx
!==================================================================================================================================!
    subroutine deallocate_em2dkx

!
! Deallocate public i/o variables in this module
!    
    if ( allocated( sigParams ) )       deallocate( sigParams )   
    if ( allocated( iFreeParam ) )      deallocate( iFreeParam )      
    if ( allocated(xTx) )               deallocate ( xTx,yTx,zTx ) 
    if ( allocated(jx) )                deallocate ( jx,jy,jz )  
    if ( allocated(mx) )                deallocate ( mx,my,mz )
    if ( allocated(xRx) )               deallocate ( xRx,yRx,zRx )
    if ( allocated(iDataMask) )         deallocate ( iDataMask )
    if ( allocated(ex_kx) )             deallocate ( ex_kx,ey_kx,ez_kx,hx_kx,hy_kx,hz_kx )                                                     
    if ( allocated(dex_kx) )            deallocate ( dex_kx)
    if ( allocated(dey_kx) )            deallocate ( dey_kx)
    if ( allocated(dez_kx) )            deallocate ( dez_kx)
    if ( allocated(dhx_kx) )            deallocate ( dhx_kx)
    if ( allocated(dhy_kx) )            deallocate ( dhy_kx)
    if ( allocated(dhz_kx) )            deallocate ( dhz_kx)
        
    if ( allocated( mesh%attr ) )       call deallocate_trimesh(mesh,.false.)
    if ( allocated( newmesh%attr ) )    call deallocate_trimesh(newmesh,.false.)
    
    end subroutine deallocate_em2dkx   
 
!==================================================================================================================================! 
!============================================================================================================================ em2dkx
!==================================================================================================================================!
    subroutine em2dkx   
!
! The main subroutine for computing the 2.5D EM fields in the kx domain
! for the input transmitter, receivers, frequency, wavenumber and mesh.
! Also computes an error estimate for each element if 
! lComputeErrorEst == .true.
!
! This routine does both 2D MT and 2.5D CSEM computations.
!
! Kerry Key 
! Scripps Institution of Oceanography
!

    implicit none
    
    if (lprintDebug_em2dkx) write(*,*) myID,': entered em2dkx...'     

!
! Initializations:
!
    call cpu_time(t0)
    
    call em2dkx_initialize
 
    call cpu_time(t1)

!
! Get element based lists (for node2tri, Rx and Tx locations):
!        
    call get_element_lists
    
!
! Generate the node labels and boundary condition flags for nodes on mesh boundaries
!
    call cpu_time(t4)
    
    call gen_nodlab  

    call gen_bcnod     

!
! Generate LHS matrix, rhs source vectors and solve the linear system(s):
!
    call cpu_time(t5)
    
    call gen_lhs 
    
    call cpu_time(t6)
    
    call gen_rhs
    
    call cpu_time(t7)
     
    call solve_primal
    
!
! If requested, compute the adjoint partial derivatives:
!    
    call cpu_time(t7b)
    if ( lCompDerivs )  then
     
        call comp_adj_derivs     
    
    endif  

!    
! Compute error estimator and refine mesh if needed; also computes accurate gradients 
! and fills in the solution arrays for all field components:
!
    call cpu_time(t8)

    call errorindicator 
  
!
! Deallocate arrays:
! 
    call cpu_time(t19)
    
    call em2dkx_deallocate  

!
! Display timing info:
!
    call cpu_time(t20)

    if (lprinttimers)  call em2dkx_printTimers   
    
!
! All done, goodbye
!    
    if (lprintDebug_em2dkx) write(*,*) myID,': leaving em2dkx...'
    
    end subroutine em2dkx
    
!==================================================================================================================================! 
!==================================================================================================================== errorindicator
!==================================================================================================================================! 
    subroutine errorindicator()
! 
! The dual-weighted residual (DWR) error estimation subroutine.
!
! Kerry Key                                               
! Scripps Institution of Oceanography         
!
! Versions:
!  Jan-April 2009
!  Summer  2008
!
    implicit none 
 
!
! Local variables:
!   
    integer         :: e, irefine, nlast
    character(256)  :: filename, cadapt
    real(8) :: area

    if (lprintDebug_em2dkx) write(*,*) myID,': errorindicator...'
    
!
!  Allocate local arrays:
!

    allocate (rhs_q(2*nedges,nTx))  
    rhs_q = 0d0
!
! Generate RHS residual vector in bump space:
!    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_B(0) ...'
    call gen_rhs_bump_B(0)  
    call cpu_time(t9)
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_F ...'
    call gen_rhs_bump_F 
    call cpu_time(t10)    
!
! Generate LHS matrix in bump space:
!
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_lhs_bump  ...'
 
    if (lBumpDiagApprox) then
        call gen_lhs_bump_diag
    else
        call gen_lhs_bump
    endif
    
    call cpu_time(t11)
!
! Solve the linear system using either SuperLU or UMFPack
!

    if (lBumpDiagApprox) then
        
        call solve_varepsilon_diag
        
    else
    
        select case (trim(linearSolver))
    
        case ('superlu')
           call factor_lhs_bump_superlu
           call solve_varepsilon_superlu
      
        case ('umfpack')
           call factor_lhs_bump_umfpack
           call solve_varepsilon_umfpack

        case default
            write(*,*) 'Error, no linear system solver parameter specified, stopping!'
            stop
              
        end select    
    
    endif
    
    call cpu_time(t12)

!
! Fill in field vectors for receivers:
!   
    call fill_solmtx  
    
    call cpu_time(t13)
    
!
! Now we do the goal oriented error estimator if a refinement iteration:
!
    if ( lusedualerror .and. lComputeErrorEst) then
!
! Compute relative errors to use as sourcing function in dual problem:
!           
     !
     ! 1. Solve dual problem in linear basis 
     !
        allocate(rhs_wh(2*mesh%nnod,nTx))
        ! DGM 2/23/2012 without the write() below, the next call fails because
        ! rhs_wh is not allocated ... but only with UMFPack and not superLU.
        ! VERY weird. This is true, even if the DEBUG is off, i.e. just the
        ! presence of this line in the compile causes something to change...
        ! Most likely, the optimizer is delaying the allocation until it sees
        ! that the variable is used. "Bad optimizer. Bad."
        if ( lprintDebug_em2dkx ) write(*,*) 'alloc: (', 2*mesh%nnod, ',', nTx, ')'
        
        call gen_rhs_dual  
        call cpu_time(t14)
        
        select case (trim(linearSolver))
        
        case ('superlu')
            
            if (lMT) then
                call solve_dual_superlu_mt 
            else
                call solve_dual_superlu
            endif
 
            
        case ('umfpack')
     
            if (lMT) then
                call solve_dual_umfpack_mt 
            else
                call solve_dual_umfpack 
            endif    

        end select    
        
        call cpu_time(t15)
     !
     ! 2.  Approximate dual error in bump basis
     !

        call gen_rhs_bump_G  

        call gen_rhs_bump_B(1)     
   
        call cpu_time(t16)   
           
     !      
     ! 3. Solve linear system for \delta_h:
     !
        if (lBumpDiagApprox) then
        
            call solve_varepsilon_diag
        
        else     
            select case (trim(linearSolver))
        
            case ('superlu')
               call solve_varepsilon_superlu
            case ('umfpack')
               call solve_varepsilon_umfpack                          
            end select    
        
        endif
                
        call cpu_time(t17)  
       
     !
     ! 4. Compute local error estimate
     !
        call copy_trimesh( mesh, newmesh )
        
        do irefine = 1, max_nsubrefine 
            
            deallocate(errnrm)
            allocate(errnrm(newmesh%nele))
            
            call comp_drw_error_proj()  
            
            if ( (irefine == 1 ).and.lSaveMeshFiles) then
            
                write(cadapt,'(i6)') meshnumber  
                cadapt = adjustl(cadapt)    
                filename  = trim(fileroot)//'.'//trim(cadapt)//'.error'
                open(16,file=trim(filename),status='REPLACE') 
                
                do e = 1,mesh%nele
                    write(16,*) errnrm(e)
                enddo
                close(16)
                
                filename  = trim(fileroot)//'.'//trim(cadapt)//'.rhs_ex'
                open(16,file=trim(filename),status='REPLACE') 
                
                do e = 1,mesh%nnod
                    write(16,*) abs(rhs(2*e-1,1))
                enddo
                close(16)
                
                filename  = trim(fileroot)//'.'//trim(cadapt)//'.rhs_hx'
                open(16,file=trim(filename),status='REPLACE') 
                
                do e = 1,mesh%nnod
                    write(16,*) abs(rhs(2*e  ,1))
                enddo
                close(16)
            
            endif
            maxerr = maxval(errnrm) 
 
            e = maxloc(errnrm,1)
            call getArea(e,area,1)
 
            
            if ( maxerr < errortolerance/1d2) exit
            
            nlast = newmesh%nnod
            
            call refineMesh()
            if (newmesh%nnod > nlast)    then   
                lmadeNewMesh= .true.
            else
                exit
            endif 
            
            if ( newmesh%nnod > maxMeshNodes) exit
        
        enddo
        call cpu_time(t18)
 
          
     !
     ! Deallocate dual array
     !
        deallocate (rhs_wh)
    
    endif ! lusedualerror
    
    
!
! Smooth the mesh:
!     
    if (lmadeNewMesh)  call smooth_TriMesh(newmesh,50,.false.)
 
    !
    ! Free up Solver memory:
    !
    if (.not.lBumpDiagApprox) then
        select case (trim(linearSolver))
        case ('superlu')
            call superlu_zfree (superlu_q)        
        case ('umfpack')
            call umfpack_zfree (umfpack_q)     
        end select 
    endif   
!        
! Deallocate arrays:
!
    if (allocated(rhs_q))    deallocate (rhs_q)
    if (allocated(val_q))    deallocate (val_q)
    
    end subroutine errorindicator    
    
!==================================================================================================================================! 
!=================================================================================================================== comp_adj_derivs
!==================================================================================================================================! 
    subroutine comp_adj_derivs   
!
! Computes the model parameter sensitivity matrix using the adjoint method
!
! Kerry Key                                               
! Scripps Institution of Oceanography         
!
    implicit none   
                               
    real(8), dimension(6) :: jax = [1,0,0,0,0,0]  ! adjoint source terms (1 for each component)
    real(8), dimension(6) :: jay = [0,1,0,0,0,0]
    real(8), dimension(6) :: jaz = [0,0,1,0,0,0]
    real(8), dimension(6) :: max = [0,0,0,1,0,0]
    real(8), dimension(6) :: may = [0,0,0,0,1,0]
    real(8), dimension(6) :: maz = [0,0,0,0,0,1]
    
    integer         ::  isite, esite, e, n(3), iTx, isrc, iparamnum(3), iRx 
    
    real(8)         :: ye(3),ze(3),a(3),b(3),c(3),area, l(3),dldy(3), dldz(3), yp,zp 
    complex(8)      :: fex(3), fhx(3), fey, fez, dexdz, dexdy, dhxdy, dhxdz,ee(3),hh(3)
    complex(8)      :: fex_adj(3), fhx_adj(3) , fey_adj, fez_adj 
    
    complex(8), dimension(6)    :: sxx, syy, szz 
    real(8)                     :: sigx,sigy,sigz 
    complex(8)                  :: gammay2, gammaz2      
    real(8), dimension(3,3)     :: M = [ 2d0, 1d0, 1d0,   1d0, 2d0, 1d0,   1d0, 1d0, 2d0] ! mass matrix, column major order

 
    integer                     :: i,ncomps, ict
    
   ! real(8) :: t0,t1
    
    if (lprintDebug_em2dkx) write(*,*) myID,': comp_adj_derivs...'     
 
    do iRx = 1,nRx
        do iTx = 1,nTx
            if (iDataMask(iRx,iTx)==2) then
                if (lComps(1)) dex_kx(iRx,iTx)%dsig = 0
                if (lComps(2)) dey_kx(iRx,iTx)%dsig = 0
                if (lComps(3)) dez_kx(iRx,iTx)%dsig = 0
                if (lComps(4)) dhx_kx(iRx,iTx)%dsig = 0
                if (lComps(5)) dhy_kx(iRx,iTx)%dsig = 0
                if (lComps(6)) dhz_kx(iRx,iTx)%dsig = 0        
            endif
        enddo
    enddo

  
    ncomps = 0
    do i = 1,6
        if (lComps(i)) ncomps = ncomps+1
    enddo
    
    allocate(rhs_wh(2*mesh%nnod,ncomps)) ! 3 Jx dipoles and 3 Mx dipoles, using dual matrix since it already has solvers coded
    
!
! Loop over each receiver and compute the adjoint sensitivity:
!
    do isite=1,nRx
        
        if ( .not.(any(iDataMask(isite,:)==2))) cycle ! no data, no need for adjoint
        
        yp = yRx(isite)
        zp = zRx(isite)
        
        esite = eRx(isite)
        call getSigs(esite,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
        gammay2 = kx2-ommu*(sigy - iomeps)
        gammaz2 = kx2-ommu*(sigz - iomeps)              

    !
    ! Set up the basis functions and evaluate at the receiver location:
    !
        n   = mesh%emap(1:3,esite)
        ye  = mesh%y(n)
        ze  = mesh%z(n)
        call get_linear_basis(ye,ze,yp,zp,l,dldy,dldz)    
                   
    !
    ! Insert the 6 sources (1 for each component):
    !
        rhs_wh = 0 ! Reset to zero 
        
        ict = 0   
        do isrc = 1,6
            
            if (lComps(isrc)) then
                ict = ict + 1
                ee =    -l*jax(isrc) - ikx  * ( dldy*jay(isrc)/gammay2      + dldz*jaz(isrc)/gammaz2 )  & 
                    &                - ommu * ( dldy*maz(isrc)*(sigy - iomeps)/gammay2 - dldz*may(isrc)*(sigz - iomeps)/gammaz2 )  
                hh = -(                          - ommu * ( dldy*jaz(isrc)/gammaz2      - dldz*jay(isrc)/gammay2 )  &
                    &   -l*max(isrc)*ommu -  ikx * ommu * ( dldy*may(isrc)/gammaz2      + dldz*maz(isrc)/gammay2 )  )  
          
                rhs_wh(2*n-1, ict) = ee  ! Ex  
                rhs_wh(2*n  , ict) = hh  ! Hx  
            endif
        enddo

    !
    ! Solve the linear system:
    !
       ! call cpu_time(t0)
        
        select case (trim(linearSolver))         
        case ('superlu')
            if (lMT) then
                call solve_dual_superlu_mt
            else
                call solve_dual_superlu
            endif
        case ('umfpack')
            if (lMT) then
                call solve_dual_umfpack_mt
            else
                call solve_dual_umfpack 
            endif
                      
        end select     
        
      !  call cpu_time(t1)
        
          
    
    !
    ! Loop over all elements and integrate them:
    !    
        do e = 1, mesh%nele
            call getParamNums(e,iparamnum)
 
            if ( all(iparamnum == 0) ) cycle  ! only include free parameters 
                 
            n    = mesh%emap(1:3,e)
            ye   = mesh%y(n)
            ze   = mesh%z(n)
      
            call get_abc_coeffs(ye,ze,a,b,c,area)   
            dldy = b/2d0/area
            dldz = c/2d0/area           
                                        
            call getSigs(e,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
            gammay2 = kx2-ommu*(sigy - iomeps)
            gammaz2 = kx2-ommu*(sigz - iomeps)    
 
            !
            ! Loop over all transmitters:
            !
            
            do iTx = 1,nTx
            
                if ( .not.(iDataMask(isite,iTx)==2)) cycle ! no data, no need for adjoint
        
            
            ! Get primal fields from source at transmitter:
                fex  = rhs(2*n-1,iTx) 
                fhx  = rhs(2*n  ,iTx)   
                
                dexdy =  sum(fex*dldy)
                dexdz =  sum(fex*dldz)
                dhxdy =  sum(fhx*dldy)
                dhxdz =  sum(fhx*dldz)
                                
                fey  = (-ommu*dhxdz - ikx*dexdy )/gammay2   ! constant across element since basis is linear
                fez  = ( ommu*dhxdy - ikx*dexdz )/gammaz2  
                    
                
                ! Loop through components and extract the adjoint fields:
                ict = 0
                do isrc = 1,6
                
                    if (lComps(isrc)) then
                    
                        ict = ict + 1
                     ! Get adjoint fields from source at receiver:
                        fex_adj  = rhs_wh(2*n-1, ict) 
                        fhx_adj  = rhs_wh(2*n  , ict)   

                                        
                        ! Sxx sensitivity can be computed using a mass matrix:
                        sxx(isrc) = sum(fex_adj*matmul(M,fex))*area/12d0
                    
                        ! Syy and Szz sensitivities:
                        dexdy =  sum(fex_adj*dldy)
                        dexdz =  sum(fex_adj*dldz)
                        dhxdy =  sum(fhx_adj*dldy)
                        dhxdz =  sum(fhx_adj*dldz)
  
                        fey_adj  = (-ommu*dhxdz - ikx*dexdy )/gammay2    
                        fez_adj  = ( ommu*dhxdy - ikx*dexdz )/gammaz2  
  
                        syy(isrc) = fey*fey_adj*area
                        szz(isrc) = fez*fez_adj*area
                                                
                        ! Apply correct sign depending on field symmetry about x/kx axis. This arise from the -ikx in the integral:
                        ! integral ( E_adj(-ikx) dot E(kx)) dkx
                        if (.not.lMT) then ! not MT
                            select case (isrc)
                        
                            case(1,5,6)  ! Jx, My and Mz have same symmetry
                        
                                syy(isrc) = -syy(isrc)
                                szz(isrc) = -szz(isrc)
                            
                            case(2,3,4) ! Jy,Jz, Mx have same symmetry
                            
                                sxx(isrc) = -sxx(isrc)
                        
                            end select  
                        endif                   
                        
                    else
                       sxx(isrc) = 0
                       syy(isrc) = 0
                       szz(isrc) = 0
                    endif
                enddo ! isrc = 1,6

!                
                ! Scale magnetic terms by ommu since  M_s^a = ommu*H_s^a            
                sxx(4:6) = sxx(4:6) / ommu
                syy(4:6) = syy(4:6) / ommu
                szz(4:6) = szz(4:6) / ommu      
                ! At this point we have 6 components each with 3 terms for the sensitivity (sxx+syy+szz). 
                ! Now we need to dole them out into the correct parameters, depending on anisotropy...
                ! hang in there...we're almost done!
                ! *** note conversion to single precision:   
                select case (trim(cAnisotropy))

                case ('isotropic')
            if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(1)) = dex_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(1) + syy(1) + szz(1) )
            if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(1)) = dey_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(2) + syy(2) + szz(2) )
            if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(1)) = dez_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(3) + syy(3) + szz(3) )
            if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(1)) = dhx_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(4) + syy(4) + szz(4) )
            if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(1)) = dhy_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(5) + syy(5) + szz(5) )
            if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(1)) = dhz_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(6) + syy(6) + szz(6) ) 
                                    
                case ('triaxial')
                    if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(1))  = dex_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(1) )
                    if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(2))  = dex_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(1) ) 
                    if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(3))  = dex_kx(isite,iTx)%dsig(iparamnum(3)) + cmplx( szz(1) ) 
                      
                    if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(1))  = dey_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(2) )
                    if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(2))  = dey_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(2) )
                    if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(3))  = dey_kx(isite,iTx)%dsig(iparamnum(3)) + cmplx( szz(2) ) 
                        
                    if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(1))  = dez_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(3) )
                    if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(2))  = dez_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(3) )
                    if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(3))  = dez_kx(isite,iTx)%dsig(iparamnum(3)) + cmplx( szz(3) )
                      
                    if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(1))  = dhx_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(4) )
                    if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(2))  = dhx_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(4) ) 
                    if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(3))  = dhx_kx(isite,iTx)%dsig(iparamnum(3)) + cmplx( szz(4) ) 
                      
                    if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(1))  = dhy_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(5) )
                    if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(2))  = dhy_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(5) )
                    if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(3))  = dhy_kx(isite,iTx)%dsig(iparamnum(3)) + cmplx( szz(5) )
                    
                    if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(1))  = dhz_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(6) ) 
                    if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(2))  = dhz_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(6) ) 
                    if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(3))  = dhz_kx(isite,iTx)%dsig(iparamnum(3)) + cmplx( szz(6) )      
                                  
                                  
                case ('tix') ! For transversely isotropic:
                     ! First param is sig along symmetry axis; second param is along transverse plane
                    if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(1))  = dex_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(1) ) 
                    if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(2))  = dex_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(1) + szz(1) ) 
 
                    if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(1))  = dey_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(2) ) 
                    if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(2))  = dey_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(2) + szz(2) )
                        
                    if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(1))  = dez_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(3) ) 
                    if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(2))  = dez_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(3) + szz(3) )
                      
                    if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(1))  = dhx_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(4) ) 
                    if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(2))  = dhx_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(4) + szz(4) )
                      
                    if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(1))  = dhy_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(5) ) 
                    if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(2))  = dhy_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(5) + szz(5) )
         
                    if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(1))  = dhz_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( sxx(6) ) 
                    if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(2))  = dhz_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(6) + szz(6) )
                                                     
    
                case ('tiy')
                
                    if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(1))  = dex_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( syy(1) )
                    if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(2))  = dex_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( sxx(1) + szz(1) )  
 
                    if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(1))  = dey_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( syy(2) ) 
                    if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(2))  = dey_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( sxx(2) + szz(2) )
                        
                    if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(1))  = dez_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( syy(3) ) 
                    if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(2))  = dez_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( sxx(3) + szz(3) )
                      
                    if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(1))  = dhx_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( syy(4) ) 
                    if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(2))  = dhx_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( sxx(4) + szz(4) )
                      
                    if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(1))  = dhy_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( syy(5) ) 
                    if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(2))  = dhy_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( sxx(5) + szz(5) )
         
                    if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(1))  = dhz_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( syy(6) ) 
                    if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(2))  = dhz_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( sxx(6) + szz(6) )
                                                            
                case ('tiz')
                
                    if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(1))  = dex_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( szz(1) ) 
                    if (lComps(1)) dex_kx(isite,iTx)%dsig(iparamnum(2))  = dex_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(1) + sxx(1) ) 
 
                    if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(1))  = dey_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( szz(2) ) 
                    if (lComps(2)) dey_kx(isite,iTx)%dsig(iparamnum(2))  = dey_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(2) + sxx(2) )
                        
                    if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(1))  = dez_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( szz(3) ) 
                    if (lComps(3)) dez_kx(isite,iTx)%dsig(iparamnum(2))  = dez_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(3) + sxx(3) )
                      
                    if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(1))  = dhx_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( szz(4) ) 
                    if (lComps(4)) dhx_kx(isite,iTx)%dsig(iparamnum(2))  = dhx_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(4) + sxx(4) )
                      
                    if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(1))  = dhy_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( szz(5) ) 
                    if (lComps(5)) dhy_kx(isite,iTx)%dsig(iparamnum(2))  = dhy_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(5) + sxx(5) )
         
                    if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(1))  = dhz_kx(isite,iTx)%dsig(iparamnum(1)) + cmplx( szz(6) ) 
                    if (lComps(6)) dhz_kx(isite,iTx)%dsig(iparamnum(2))  = dhz_kx(isite,iTx)%dsig(iparamnum(2)) + cmplx( syy(6) + sxx(6) )
                                                      
            
                end select              
 
 
            enddo ! itx = 1,nTx 
            
            
            
        enddo ! e = 1 nele      
        
        !call cpu_time(t2)
        
       ! write(*,'(i5,a,g16.6,1x,a,g16.6)') isite,'solve time: ',t1-t0,'integrate time: ', t2-t1                      
    enddo  !  isite=1,nRx 
 
     
    deallocate(rhs_wh)
    
    end subroutine comp_adj_derivs

!==================================================================================================================================! 
!================================================================================================================= em2dkx_initialize
!==================================================================================================================================! 
    subroutine em2dkx_initialize   
!
! Initializes variables, checks input options and allocates some arrays
!
! Kerry Key                                               
! Scripps Institution of Oceanography         
!

    implicit none       

    character(32)  :: cend   
    

! 
! Check kx:
!
    if(lMT)  kx = 0             
 
! 
! Check input linear solver parameter:
!
    select case (trim(linearSolver))
        case ('superlu','umfpack')
        case default
            write(*,*) ' Error in em2dkx, unknown linear solver: ',trim(linearSolver)
            write(*,*) ' Stopping !'
            stop
    end select
 
  
! 
! Initialize a few things:
!        
    kx2    = kx**2
    ikx    = ic*kx
    ommu   = ic*w*mu0
    iomeps = ic*w*EPSILON
    
    lmadeNewMesh = .false.
    
    ! Triangle command:
    write(cend,fmt='(F4.0)') minQangle
    retricommand = 'q'//trim(adjustl(cend))//'rpnajQ'//CHAR(0) 
 
!
! Allocate private arrays used only in em2dkx:
!
    if (lprintDebug_em2dkx) write(*,*) myID,': allocating in em2dkx...'
    
    allocate ( errnrm(mesh%nele) )
    allocate ( nodlab(mesh%nnod), bcnod(mesh%nnod)  )
    allocate ( rhs(2*mesh%nnod,nTx) )
    allocate ( eRx(nRx) )
    allocate ( lCloud(mesh%nele) )
    allocate ( edges(3,mesh%nele) )
    allocate ( eTx(nTx) )


    end subroutine em2dkx_initialize
        
!==================================================================================================================================! 
!================================================================================================================= em2dkx_deallocate
!==================================================================================================================================! 
    subroutine em2dkx_deallocate   
!
! Deallocates all remaining internal arrays at end of em2dkx
!
! Kerry Key                                               
! Scripps Institution of Oceanography         
!    
    implicit none   
    
!
! Deallocate all superlu memory and associated arrays:
! 
    if (lprintDebug_em2dkx) write(*,*) myID,': deallocating sparse matrices...'

    if (.not.lMT) then

        select case (trim(linearSolver))
         
        case ('superlu')
            call superlu_zfree (superlu_p)   
            
        case ('umfpack')
            call umfpack_zfree (umfpack_p)   
            
        end select  
        
        if ( allocated(val) )       deallocate(val)
        if ( allocated(col) )       deallocate(col)
        if ( allocated(irw) )       deallocate(irw)    
 
    
    else  ! MT    
        
        select case (trim(linearSolver))
        
        case ('superlu')
            call superlu_zfree (superlu_p_E)   
            call superlu_zfree (superlu_p_H)  
            
        case ('umfpack')
            call umfpack_zfree (umfpack_p_E)   
            call umfpack_zfree (umfpack_p_H)   
            
        end select 
        
        if ( allocated(rhs_lE) )    deallocate(rhs_lE) 
        if ( allocated(val_lE) )    deallocate(val_lE)
        if ( allocated(col_lE) )    deallocate(col_lE)
        if ( allocated(irw_lE) )    deallocate(irw_lE) 
        
        if ( allocated(rhs_lH) )    deallocate(rhs_lH)                
        if ( allocated(val_lH) )    deallocate(val_lH)
        if ( allocated(col_lH) )    deallocate(col_lH)
        if ( allocated(irw_lH) )    deallocate(irw_lH)                                           
 
    end if 



    if (lprintDebug_em2dkx) write(*,*) myID,': deallocating other arrays...'
    deallocate ( errnrm ) 
    deallocate ( nodlab, bcnod ) 
    deallocate ( rhs )  
    deallocate ( eRx )
    deallocate( lCloud)
    deallocate ( edges )
    deallocate ( eTx )
    call kdtree2_destroy(tree)  
    deallocate( node2tri )
    
    end subroutine em2dkx_deallocate  
    
!==================================================================================================================================! 
!================================================================================================================ em2dkx_printTimers
!==================================================================================================================================! 
    subroutine em2dkx_printTimers 
!
! Kerry Key                                               
! Scripps Institution of Oceanography         
!
    implicit none   
      
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx 1st allocations:  ',t1-t0
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx computeEdges:     ',t2-t1
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx kdtree2:          ',t3-t2
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx findElement:      ',t4-t3
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx gen_nodlab, bcnod:',t5-t4        
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx gen_lhs:          ',t6-t5
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx gen_rhs:          ',t7-t6
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx solve_primal:     ',t7b-t7
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx comp_adj_derivs : ',t8-t7b
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx gen_rhs_bump_F:   ',t9-t8
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx gen_rhs_bump_B:   ',t10-t9
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx gen_lhs_bump:     ',t11-t10
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx solve_varepsilon: ',t12-t11
    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx fill_solmtx:      ',t13-t12
    
    if (lComputeErrorEst) then
        write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx gen_rhs_dual:     ',t14-t13
        write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx solve_dual:       ',t15-t14
        write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx gen_rhs_bump_GB:  ',t16-t15        
        write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx solve_varepsilon: ',t17-t16    
        write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx comp_drw_error:   ',t18-t17  
        write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx deallocations:    ',t20-t18    
    else
         write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx deallocations:    ',t20-t19    
    endif   

    write(*,'(i4,1x,a24,1x,f8.4)') myID,' em2dkx total:            ',t20-t0

    end subroutine em2dkx_printTimers

!==================================================================================================================================! 
!=================================================================================================================== localRefinement
!==================================================================================================================================! 
    subroutine localRefinement
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
! 
    implicit none
    
    integer                 :: iTx, iRx, i,e, n(3), inele(1)
    real(8)                 :: area, skinDepth,r, rmin, scaling, rminRx,rminTx 
    real(8), dimension(3)   :: ye,ze 
    character(32)           :: cend  
    logical                 :: lMakeNewMesh 
    real(kdkind)            :: qv(2) 
    type(kdtree2_result)    :: results(1)   
    type(kdtree2), pointer  :: treeTx, treeRx, treeNodes      
    real(8), dimension(:), allocatable  :: sigTemp
    real(8)                             :: sigx,sigy,sigz
            
    ! Triangle command:
    write(cend,fmt='(F4.0)') minQangle
    retricommand = 'q'//trim(adjustl(cend))//'rpnajQ'//CHAR(0) 
    
    allocate ( eRx(nRx), eTx(nTx) )  

!
! Create some kdtrees of the Tx's and Rx's for fast closest point searches:
! 
    if (nRx > 1) then
        treeRx => kdtree2_create( transpose(reshape( [yRx, zRx]  , [nRx,2]) ), sort=.true.,rearrange=.true.)
        nullify(treeRx%the_data)  
    endif    
    if (nTx > 1) then
        treeTx => kdtree2_create( transpose(reshape( [yTx, zTx]  , [nTx,2]) ), sort=.true.,rearrange=.true.)  
        nullify(treeTx%the_data)  
    endif    

!
! Loop on a priori refinement until all elements are smaller than requested tolerances:
!    
    do  
    
        lMakeNewMesh = .false.

        allocate ( mesh%area(mesh%nele) )
        
        mesh%area  = -1d0
        
        if (lprintDebug_em2dkx) write(*,*) 'local refinement loop, mesh%nnod: ',mesh%nnod
        
      
        !
        ! Refine elements containing receivers:
        !
        allocate (sigTemp(mesh%nele))
                    
        do e = 1,mesh%nele
            call getSigs(e,sigx,sigy,sigz)
            sigTemp(e) = (sigx+sigy+sigz)/3d0
        enddo
        
        allocate (node2tri(mesh%nnod))
        call getNode2Tri( mesh%nele, mesh%emap,mesh%nnod, node2tri)                
 
        treeNodes => kdtree2_create( transpose(reshape( [mesh%y, mesh%z]  , [size(mesh%y),2]) ), sort=.true.,rearrange=.true.) 
        nullify(treeNodes%the_data) 
                   
        do iRx = 1,nRx
  
            call findElement(treeNodes, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
            & sigTemp, 1, yRx(iRx), zRx(iRx), inele)
 
            call getArea(inele(1),area,0) 
 
            if ( area >  nRxEleArea ) then  
            
                mesh%area(inele) = area/2.0
                lMakeNewMesh= .true.  
                
            endif
            
        enddo      
        
        deallocate (sigTemp,node2tri)
        call kdtree2_destroy(treeNodes)
                   
          
        do e = 1,mesh%nele
  
            call getArea(e,area,0) 
            call getSkinDepth(e,skinDepth,0,.false.)
 
       
            rmin = 1d100
            
            rminRx = 1d100
            rminTx = 1d100
            
            n    = mesh%emap(1:3,e)
            ye   = mesh%y(n)
            ze   = mesh%z(n)    

            do i = 1,3

                qv(1) = ye(i)
                qv(2) = ze(i)
                            
                !    
                ! Find closest Rx:
                !
                if (nRx > 1) then
                    call kdtree2_n_nearest(tp=treeRx,qv=qv,nn=1,results=results) 
                else
                    results(1)%idx = 1
                endif
                
                r = sqrt( (yRx(results(1)%idx) - ye(i))**2 + (zRx(results(1)%idx) - ze(i))**2)
            
                if (r < rminRx) rminRx = r  
            
                !
                ! Find closest Tx:
                !
                if (.not.lMT) then
                    if (nTx > 1) then
                        call kdtree2_n_nearest(tp=treeTx,qv=qv,nn=1,results=results) 
                    else
                        results(1)%idx = 1
                    endif                
                    r = sqrt( (yTx(results(1)%idx) - ye(i))**2 + (zTx(results(1)%idx) - ze(i))**2)
            
                    if (r < rminTx) rminTx = r  
                
                 endif                  
                 
            enddo
            
            !
            ! Refinement in a cloud surrounding the Rx's and Tx's:
            ! 
            rmin = min(rminRx,rminTx)
            
            if  ( (rmin/(2.*skinDepth)) > 50.)  then  ! catch to avoid overflow
                scaling = exp(50.)
            else
                scaling = exp(rmin/(2.*skinDepth))
            endif
 
            if ( sqrt(area) > (  skinDepth*nSkinDepthsRegional*scaling) ) then
                mesh%area(e) = area/2d0
                lMakeNewMesh= .true.  
            endif  
                
    
              
        enddo ! Loop over elements

        if (.not.lMT) then
        
            treeNodes => kdtree2_create( transpose(reshape( [mesh%y, mesh%z]  , [size(mesh%y),2]) ), sort=.true.,rearrange=.true.) 
            nullify(treeNodes%the_data) 
            
            allocate (sigTemp(mesh%nele))
             
            do e = 1,mesh%nele
                call getSigs(e,sigx,sigy,sigz)
                sigTemp(e) = (sigx+sigy+sigz)/3d0
            enddo
    
            allocate (node2tri(mesh%nnod))
    
            call getNode2Tri( mesh%nele, mesh%emap,mesh%nnod, node2tri)
    
            do iTx = 1,nTx
 
                call findElement(treeNodes, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
                & sigTemp, 1, yTx(iTx), zTx(iTx), inele)
                
                call getArea(inele(1),area,0) 
                if ( area > nTxEleArea) then
                
                    mesh%area(inele) = area/2.0
                    lMakeNewMesh= .true.  
                    
                endif
                
            enddo
            
            deallocate(node2tri,sigTemp)           
            call kdtree2_destroy(treeNodes)
 
        
        endif
             
        ! 
        ! Refine the mesh:
        !
        if (lMakeNewMesh ) then        
            call call_triangle(retricommand, mesh )
            if (allocated(mesh%area))   deallocate(mesh%area)     
            
            ! Smooth the new mesh:
            call smooth_TriMesh(mesh,50,.false.)

        endif
 
        if (.not.lMakeNewMesh) then
            if (lprintDebug_em2dkx) write(*,*) 'a priori refinement: ',mesh%nnod, mesh%nele
            exit
        endif
                
    enddo
    
    if (allocated(mesh%area))   deallocate(mesh%area) 
    if (allocated(eRx))         deallocate(eRx) 
    if (allocated(eTx))         deallocate(eTx) 
 
    if (nRx > 1) call kdtree2_destroy(treeRx)  
    if (nTx > 1) call kdtree2_destroy(treeTx)
      
    end subroutine localRefinement
    
!==================================================================================================================================! 
!====================================================================================================================== getSkinDepth
!==================================================================================================================================! 
    subroutine getSkinDepth(e,skinDepth,iopt,lTx)
    
    integer, intent(in)   :: e
    integer, intent(in)   :: iopt
    logical, intent(in), optional :: lTx
    real(8), intent(out)  :: skinDepth 
    
    real(8)               :: sig,sigx,sigy,sigz
         
    !
    ! Get conductivity of element e:
    !    
    call getSigs(e,sigx,sigy,sigz,iopt)
    sig = maxval([sigx,sigy,sigz])
       
    skinDepth = sqrt(2d0 / (sig*MU0*w) )  
    
    end subroutine getSkinDepth

!==================================================================================================================================! 
!=========================================================================================================================== getArea
!==================================================================================================================================! 
    subroutine getArea(e,area,iopt)
    
    integer, intent(in)   :: e
    real(8), intent(out)  :: area
    integer, intent(in)    :: iopt
        
    integer               :: n(3)  
    real(8), dimension(3) :: ye,ze,a,b,c
    
    !
    ! Area of element e:        
    !
    if (iopt == 1) then
        n    = newmesh%emap(1:3,e)
        ye   = newmesh%y(n)
        ze   = newmesh%z(n)              
    else
        n    = mesh%emap(1:3,e)
        ye   = mesh%y(n)
        ze   = mesh%z(n)   
    endif  
    call get_abc_coeffs(ye,ze,a,b,c,area)    
    
    
    end subroutine getArea
           
!==================================================================================================================================! 
!================================================================================================================= get_element_lists
!==================================================================================================================================! 
    subroutine get_element_lists
!
! Computes the edge list, site element list and the transmitter element lists
!    
! Kerry Key                                               
! Scripps Institution of Oceanography         
!    
    implicit none
    
!
! Local variables:
!
    integer                             :: e
    real(8), dimension(:), allocatable  :: sigTemp
    real(8)                             :: sigx,sigy,sigz
     
!
! Create element edge list:
!
    if (lprintDebug_em2dkx) write(*,*) myID,': computeEdges...'
    
    
    call computeEdges( mesh%nele, mesh%emap, mesh%neighborlist, edges, nedges )      
    
    call cpu_time(t2)
        
!
! Get site element lists
!
    if (lprintDebug_em2dkx) write(*,*) myID,': getting site element lists...'
    
    !
    ! First make the kdtree:
    !  
    
    if (lprintDebug_em2dkx) write(*,*) myID,': getting site element lists...kdtree2',size(mesh%y),size(mesh%z)
    
    tree => kdtree2_create( transpose(reshape( [mesh%y, mesh%z]  , [size(mesh%y),2]) ), sort=.true.,rearrange=.true.) 
    ! DGM 2/23/2012: the_data is *STACK* based from above line and immediately disappears after the call.
    !       Don't keep a pointer around to some ficticious place in the stack. Bad juju.
    nullify(tree%the_data)
    
    call cpu_time(t3)
    
    ! Get node2tri array:
    
    if (lprintDebug_em2dkx) write(*,*) myID,': getting site element lists...getNode2Tri'
    
    allocate (node2tri(mesh%nnod))
    
    call getNode2Tri( mesh%nele, mesh%emap,mesh%nnod, node2tri)
 
    !
    ! Now call the new jump-and-walk findElement routine:
    !
    
    if (lprintDebug_em2dkx) write(*,*) myID,': getting site element lists...findElement'
    
    allocate (sigTemp(mesh%nele))
    
    ! We will use the average conductivity to define the element for the searches
    ! since if the search finds the point is incident on a node or edge, it chooses the most conductive element.
    do e = 1,mesh%nele
        call getSigs(e,sigx,sigy,sigz)
        sigTemp(e) = (sigx+sigy+sigz)/3d0
    enddo
    
    call findElement(tree, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
    & sigTemp, nRx, yRx, zRx, eRx)
    
    deallocate(sigTemp)    
    
!
! If CSEM computation, get transmitter element list eTx:
!    
    if (.not.lMT) then 
 
        allocate (sigTemp(mesh%nele)) 
    
        do e = 1,mesh%nele
          call getSigs(e,sigx,sigy,sigz)
          sigTemp(e) = (sigx+sigy+sigz)/3d0
        enddo
        
        call findElement(tree,mesh%nnod,mesh%y,mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
        & sigTemp, nTx,yTx,zTx,eTx)
    
        deallocate(sigTemp)
    
    endif
    
    end subroutine get_element_lists
           
!==================================================================================================================================! 
!======================================================================================================================== gen_nodlab
!==================================================================================================================================! 
    subroutine gen_nodlab
!
! Subroutine to label the mesh boundaries.  Assumes that the mesh
! has an exact rectangular outer boundary aligned along the coordinate 
! axes!
! 
!  nodlab array flags:       
!   0       Mesh interior
!   1       Top edge, including left and right most corners
!   2       Bottom edge, including left and right most corners
!   3       Left edge 
!   4       Right edge
!
! y is positive right, z is positive down
!
!
    implicit none
    
    real(8)  :: ymin,ymax,zmin,zmax 

    if (lprintDebug_em2dkx) write(*,*) myID,': gen_nodlab...' 
      
    nodlab = 0 ! default is interior node
    
    ymin = minval(mesh%y)
    ymax = maxval(mesh%y)
    zmin = minval(mesh%z)
    zmax = maxval(mesh%z)
    
    ! First set left and right side flags:
    where (mesh%y == ymax) nodlab = 4  ! right
    where (mesh%y == ymin) nodlab = 3  ! left
    where (mesh%z == zmax) nodlab = 2  ! bottom 
    where (mesh%z == zmin) nodlab = 1  ! top
   
    end subroutine gen_nodlab      


!==================================================================================================================================! 
!========================================================================================================================= gen_bcnod
!==================================================================================================================================! 
    subroutine gen_bcnod
!                                                                      
!  Routine to remap the vector of node labels in 'nodlab' to a 
!  new vector of labels bcnod that is appropriate for particular boundary
!  conditions. 
!
! Mapping:
!
!  nodlab   bcnod  
!   0           0   Mesh interior
!   1           1   Top edge, including left and right most corners
!   2           1   Bottom edge, including left and right most corners
!   3           1   Left edge 
!   4           1   Right edge
!           
! 
!  Values within bcnod indicate the fate of the node they represent:   
!                                                                      
!    0  ...  interior node where FE solution is obtained.              
!    1  ...  boundary node where Dirichlet condition is applied.       
!    
!
    implicit none

    if (lprintDebug_em2dkx) write(*,*) myID,': gen_bcnod...'  
      
    bcnod = 0
!
!  Remap the node labels 
!
    where ( nodlab > 0) bcnod = 1
    
    end subroutine gen_bcnod
    
!==================================================================================================================================! 
!=========================================================================================================================== getSigs
!==================================================================================================================================!  

    subroutine getSigs(e,sigx,sigy,sigz,iopt) 
!    
! Gets anisotropic conductivity parameters for element e
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Nov 2011   - modified for new positive only element attributes
! April 2011 - implemented
!                                                               
! 

    implicit none
 
    integer, intent(in) :: e
    integer, intent(in), optional :: iopt
    real(8), intent(out)   :: sigx,sigy,sigz
    
    ! Local variables:
    integer :: iparamnum
      
 !
 ! Get the anisotropy type and its associated conductivities
 !
    if (present(iopt)) then
        if (iopt==1) then 
            iparamnum =  nint(abs(newmesh%attr(e)))
        else
             iparamnum =  nint(abs(mesh%attr(e)))
        endif
    else
        iparamnum =  nint(abs(mesh%attr(e)))
    endif
   
    select case (trim(cAnisotropy))
    
    case ('isotropic')
        sigx = sigParams(iparamnum)
        sigy = sigParams(iparamnum) 
        sigz = sigParams(iparamnum)
    case ('triaxial')
        sigx = sigParams(3*(iparamnum-1)+1)
        sigy = sigParams(3*(iparamnum-1)+2)
        sigz = sigParams(3*(iparamnum-1)+3)
    case ('tix')                               ! For transversely isotropic:
        sigx = sigParams(2*(iparamnum-1)+1)     ! First param is sig along symmetry axis
        sigy = sigParams(2*(iparamnum-1)+2)     ! Second param is along transverse plane
        sigz = sigParams(2*(iparamnum-1)+2)         
    case ('tiy')
        sigx = sigParams(2*(iparamnum-1)+2)
        sigy = sigParams(2*(iparamnum-1)+1)
        sigz = sigParams(2*(iparamnum-1)+2)         
    case ('tiz')
        sigx = sigParams(2*(iparamnum-1)+2)
        sigy = sigParams(2*(iparamnum-1)+2)
        sigz = sigParams(2*(iparamnum-1)+1)         
    case default
        write(*,*) 'Error decoding anisotropy in getSigs in em2dkx.f90'
        write(*,*) 'Unknown or unsupported anisotropy code: ',trim(cAnisotropy)
        stop

    end select
    
     !write(*,*) iparamnum, 1./real(sigx),1./real(sigy),1./real(sigz)  
     
    end subroutine getSigs
!==================================================================================================================================! 
!====================================================================================================================== getParamNums
!==================================================================================================================================!  

    subroutine getParamNums(e,iparamnum) 
!
! Gets parameter indices for elment e
! 

    implicit none
 
    integer, intent(in)  :: e
    integer, intent(out) :: iparamnum(3)
    
    integer :: iregion
    iparamnum = 0
 !
 ! Get the anisotropy type and its associated conductivities
 !
    iregion =  nint(abs(mesh%attr(e)))
       
    select case (trim(cAnisotropy))
    
    case ('isotropic')
        iparamnum(1) = iFreeParam(iregion)
    case ('triaxial')
        iparamnum(1) = iFreeParam(3*(iregion-1) + 1)
        iparamnum(2) = iFreeParam(3*(iregion-1) + 2)
        iparamnum(3) = iFreeParam(3*(iregion-1) + 3)
    case ('tix','tiy','tiz')                               ! For transversely isotropic:
        iparamnum(1) = iFreeParam(2*(iregion-1) + 1)
        iparamnum(2) = iFreeParam(2*(iregion-1) + 2)    
        
    case default
        write(*,*) 'Error decoding anisotropy in getSigs in em2dkx.f90'
        write(*,*) 'Unknown or unsupported anisotropy code: ',trim(cAnisotropy)
        stop

    end select
    
     !write(*,*) iparamnum, real(sigx),real(sigy),real(sigz)  
     
    end subroutine getParamNums
        
!==================================================================================================================================! 
!========================================================================================================================= gen_1d_mt 
!==================================================================================================================================! 
      subroutine gen_1d_mt(model1D,cSide,cSigDirection)
!                                                                      
! Routine to generate 1D layered models based on the conductivity     
! profile of the right side of the mesh (y-max) and another along the
! left side (y-min).  Now also computes the 1D MT layer coefficients 
! at the current angular frequency "w".  
!
! This routine originated from one Chester
! Weiss and I wrote for MARE2DMT, but has been updated here using new
! Fortran 2003 intrinsics.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! September 2009: Created
! April 22, 2011: Modified for anisotropic models. TE uses sigx, TM uses sigy
!                 Now has model derived type as IO object.
!                                                               
! 

    implicit none
  
    character(5), intent(in)  :: cSide         ! 'left' or 'right'
    character(1), intent(in)  :: cSigDirection ! 'x', 'y' for TE or TM anisotropic 1D conductivity
    type(PW1D), intent(inout) :: model1D     
    

    real(8)                  :: yposition  ! position where 1D model to be generated (needs to be leftmost or rightmost side)    
    integer                  :: nSideNodes,n(3)
    real(8), dimension(:), allocatable :: sigSide, zSide 
    
    integer  :: i,e, eSide
    real(8)  :: d1,d2, sigx,sigy,sigz

    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_1d_mt ...'  
!
! First count the number of nodes on the requested side
!
    select case (trim(cSide))
    case ('left')
        yposition =  minval(mesh%y)
    case('right')
        yposition =  maxval(mesh%y)  
    case default
        write(*,*) 'Error in gen_1d_mt, undefined cSide: ', trim(cSide) 
        write(*,*) 'stopping!'
        stop
    end select
    
    nSideNodes  = count(mesh%y == yposition)
    
    if (nSideNodes <= 0) then
        write(*,*) 'Error in gen_1d_mt, number of side nodes: ', nSideNodes, ' y position for profile: ', yposition
        write(*,*) 'stopping!'
        stop
    endif
 
!
! Allocate buffer arrays:
!
    allocate(sigSide(nSideNodes), zSide(nSideNodes) )
    

!
!  Extract an unsorted list of element conductivities along the requested side of the mesh:
!
    eSide = 0
    do e=1,mesh%nele  
        n = mesh%emap(1:3,e)
        
        do i = 1,3
        
            d1 = abs(mesh%y(n(eps(i)))   - yposition) 
            d2 = abs(mesh%y(n(eps(i+1))) - yposition) 
            
            if ( (d1 < 1d-6) .and. (d2 < 1d-6) ) then
                eSide = eSide + 1
                
                call getSigs(e,sigx,sigy,sigz)
                
                select case (trim(cSigDirection))
                case('x')
                    sigSide(eSide) = sigx
                case('y')
                    sigSide(eSide) = sigy
                case default
                    write(*,*) 'Error in gen_1d_mt, undefined cSide: ', trim(cSigDirection) 
                    write(*,*) 'stopping!'
                    stop                   
                end select 
                
                zSide(eSide)   = min(mesh%z(n(eps(i))),mesh%z(n(eps(i+1))))            ! Stores top depth
            
            endif
       
        enddo
    enddo
!
! Sort the lists on z:
!
    call quick_sort_d(zSide(1:eSide),sigSide(1:eSide))
  
!
! Find layers of constant conductivity:
!
     nSideNodes = 1
     do i=2,eSide
        if ( sigSide(i) /= sigSide(nSideNodes)  ) then
            nSideNodes = nSideNodes + 1 
            ! Record new layer top:
            zSide(nSideNodes) = zSide(i)
            sigSide(nSideNodes) = sigSide(i)
        endif
     enddo
  
!    write(*,*) ' 1D model for side:',trim(cSide)
!    do i=1,nSideNodes
!       write(*,*) i,zSide(i),(sigSide(i))
!    enddo    

!
! Now insert these into the PW1D derived type: 
!
    model1D%nlayer = nSideNodes
    allocate (model1D%zlay(nSideNodes),model1D%sig(nSideNodes))
    model1D%zlay  = zSide(1:nSideNodes)
    model1D%sig   = sigSide(1:nSideNodes)
    model1D%omega = w
 
!
! Deallocate local arrays:
!
    deallocate( zSide, sigSide )
      
    end subroutine gen_1d_mt
!==================================================================================================================================! 
!====================================================================================================================== getStiff_lin 
!==================================================================================================================================! 
    subroutine getStiff_lin(e,mtx)
!
! Generates the 6 x 6 nodal linear FE stiffness matrix for the 2.5D EM problem.
! 
! Kerry Key
! Scripps Institution of Oceanography
!
! April 21, 2011: Modified to handle triaxial anisotropy, woo-hoo!
!
      
    implicit none
    
    integer, intent(in) :: e
    complex(8),dimension(21), intent(out) :: mtx   

! Local variables:  
    integer, dimension(3)   :: n
    real(8), dimension(3)   :: ye,ze
    real(8)                 :: sigx,sigy,sigz
    complex(8)              :: gammay2,gammaz2
    real(8), dimension(3)   :: a,b,c
    real(8)                 :: ar
    complex(8)              :: aey,aez,bex,cey,cez,ahy,ahz,bhx,chy,chz ! pde coefficients for e and h terms
  
!
! Get triaxial conductivity for element e
!
    call getSigs(e,sigx,sigy,sigz) 
    gammay2 = kx2-ommu*(sigy - iomeps)
    gammaz2 = kx2-ommu*(sigz - iomeps)
        
!
! Get Linear basis Coefficients:
!
    n    = mesh%emap(1:3,e)
    ye   = mesh%y(n)
    ze   = mesh%z(n)   
    call get_abc_coeffs(ye,ze,a,b,c,ar)
 
!
! 2.5D EM  coefficients for triaxial conductivity:
!
    aey = (sigy - iomeps) / (gammay2*4d0*ar)
    aez = (sigz - iomeps) / (gammaz2*4d0*ar)
    bex = (sigx - iomeps) * ar / 6d0
    cey =  ikx / (gammay2*4d0*ar)
    cez =  ikx / (gammaz2*4d0*ar)

    ahz = -ommu / (gammay2*4d0*ar)    ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    ahy = -ommu / (gammaz2*4d0*ar)    ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    bhx = -ommu * ar / 6d0
    chy = -cez 
    chz = -cey 

!
! Elements of the stiffness matrix (upper triangular part only since its symmetric):
!
    ! E1
    mtx(1)  = aey*b(1)*b(1) + aez*c(1)*c(1) + bex      ! e1 e1
    mtx(2)  = cez*b(1)*c(1) - cey*c(1)*b(1)            ! e1 h1 not zero when anisotropic!
    mtx(3)  = aey*b(2)*b(1) + aez*c(2)*c(1) + bex/2d0  ! e1 e2
    mtx(4)  = cez*b(2)*c(1) - cey*c(2)*b(1)            ! e1 h2 
    mtx(5)  = aey*b(3)*b(1) + aez*c(3)*c(1) + bex/2d0  ! e1 e3
    mtx(6)  = cez*b(3)*c(1) - cey*c(3)*b(1)            ! e1 h3
    ! H1      
    mtx(7)  = ahy*b(1)*b(1) + ahz*c(1)*c(1) + bhx      ! h1 h1
    mtx(8)  = chz*b(2)*c(1) - chy*c(2)*b(1)            ! h1 e2 
    mtx(9)  = ahy*b(2)*b(1) + ahz*c(2)*c(1) + bhx/2d0  ! h1 h2
    mtx(10) = chz*b(3)*c(1) - chy*c(3)*b(1)            ! h1 e3
    mtx(11) = ahy*b(3)*b(1) + ahz*c(3)*c(1) + bhx/2d0  ! h1 h3
    ! E2
    mtx(12) = aey*b(2)*b(2) + aez*c(2)*c(2) + bex      ! e2 e2
    mtx(13) = cez*b(2)*c(2) - cey*c(2)*b(2)            ! e2 h2 not zero when anisotropic!
    mtx(14) = aey*b(3)*b(2) + aez*c(3)*c(2) + bex/2d0  ! e2 e3 
    mtx(15) = cez*b(3)*c(2) - cey*c(3)*b(2)            ! e2 h3
    !H2
    mtx(16) = ahy*b(2)*b(2) + ahz*c(2)*c(2) + bhx      ! h2 h2
    mtx(17) = chz*b(3)*c(2) - chy*c(3)*b(2)            ! h2 e3
    mtx(18) = ahy*b(3)*b(2) + ahz*c(3)*c(2) + bhx/2d0  ! h2 h3
    !E3
    mtx(19) = aey*b(3)*b(3) + aez*c(3)*c(3) + bex      ! e3 e3
    mtx(20) = cez*b(3)*c(3) - cey*c(3)*b(3)            ! e3 h3 not zero when anisotropic!
    !H3 
    mtx(21) = ahy*b(3)*b(3) + ahz*c(3)*c(3) + bhx      ! h3 h3

    end subroutine getStiff_lin

!==================================================================================================================================! 
!===================================================================================================================== getStiff_bump 
!==================================================================================================================================! 
    subroutine getStiff_bump(e,mtx)
!
! Routine to compute the element stiffness matrix for the coupled 2.5D
! CSEM system in the space of quadratic edge bump functions.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! June 2008:        Created
! April 21, 2011:   Modified to handle triaxial anisotropy. Now has cross
!                   coupling terms between E-H that vanish if isotropic.
       
    implicit none
    
    integer, intent(in)     :: e
    complex(8),intent(out)  :: mtx(21)   

! Local variables:  
    integer, dimension(3)   :: n
    real(8), dimension(3)   :: ye,ze
    real(8)                 :: sigx,sigy,sigz
    complex(8)              :: gammay2,gammaz2
    real(8), dimension(3)   :: a,b,c
    real(8)                 :: ar, ar1o3,ar2o3
    complex(8)              :: aey,aez,bex,cey,cez,ahy,ahz,bhx,chy,chz ! pde coefficients for e and h terms

!
! Get triaxial conductivity for element e
!
    call getSigs(e,sigx,sigy,sigz) 
    gammay2 = kx2-ommu*(sigy - iomeps)
    gammaz2 = kx2-ommu*(sigz - iomeps)
        
!
! Get Linear basis Coefficients:
!
    n    = mesh%emap(1:3,e)
    ye   = mesh%y(n)
    ze   = mesh%z(n)   
    call get_abc_coeffs(ye,ze,a,b,c,ar)
 
!
! 2.5D EM  coefficients for triaxial conductivity:
!
    aey = (sigy - iomeps) / gammay2
    aez = (sigz - iomeps) / gammaz2
    bex = (sigx - iomeps) * ar / 45d0
    cey =  ikx / gammay2  
    cez =  ikx / gammaz2

    ahz = -ommu / gammay2     ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    ahy = -ommu / gammaz2     ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    bhx = -ommu * ar / 45d0
    chy = -cez 
    chz = -cey 
   
    ar2o3 = 2d0/(3d0*ar)
    ar1o3 = 1d0/(3d0*ar)
    
!
! Elements of the stiffness matrix (upper triangular part only since its symmetric):
!
    ! E1
    mtx(1)  = ar2o3*(aey*(b(1)*b(1) - b(2)*b(3)) + aez*(c(1)*c(1) - c(2)*c(3)) ) + bex*8d0 ! e1 e1
    mtx(2)  = ar1o3*(cez-cey)*(b(1)*c(1) + b(2)*c(2) + b(3)*c(3) )                         ! e1 h1, not zero when anisotropic!
    mtx(3)  = ar2o3*(aey*b(2)*b(1) + aez*c(2)*c(1) )                             + bex*4d0 ! e1 e2
    mtx(4)  = ar1o3*(cez-cey)*(b(1)*c(2) + b(2)*c(1))                                      ! e1 h2 
    mtx(5)  = ar2o3*(aey*b(1)*b(3) + aez*c(1)*c(3) )                             + bex*4d0 ! e1 e3
    mtx(6)  = ar1o3*(cez-cey)*(b(1)*c(3) + b(3)*c(1))                                      ! e1 h3
    ! H1     
    mtx(7)  = ar2o3*(ahy*(b(1)*b(1) - b(2)*b(3)) + ahz*(c(1)*c(1) - c(2)*c(3)) ) + bhx*8d0 ! h1 h1
    mtx(8)  = ar1o3*(chz-chy)*(b(2)*c(1) + b(1)*c(2))                                      ! h1 e2 
    mtx(9)  = ar2o3*(ahy*b(2)*b(1) + ahz*c(2)*c(1) )                             + bhx*4d0 ! h1 h2
    mtx(10) = ar1o3*(chz-chy)*(b(3)*c(1) + c(3)*b(1))                                      ! h1 e3
    mtx(11) = ar2o3*(ahy*b(3)*b(1) + ahz*c(3)*c(1) )                             + bhx*4d0 ! h1 h3
    ! E2
    mtx(12) = ar2o3*(aey*(b(2)*b(2) - b(1)*b(3)) + aez*(c(2)*c(2) - c(1)*c(3)) ) + bex*8d0 ! e2 e2
    mtx(13) = ar1o3*(cez-cey)*(b(1)*c(1) + b(2)*c(2) + b(3)*c(3) )                         ! e2 h2, not zero when anisotropic! 
    mtx(14) = ar2o3*(aey*b(3)*b(2) + aez*c(3)*c(2))                              + bex*4d0 ! e2 e3 
    mtx(15) = ar1o3*(cez-cey)*(b(3)*c(2) + c(3)*b(2))                                      ! e2 h3
    !H2
    mtx(16) = ar2o3*(ahy*(b(2)*b(2) - b(1)*b(3)) + ahz*(c(2)*c(2) - c(1)*c(3)) ) + bhx*8d0 ! h2 h2
    mtx(17) = ar1o3*(chz-chy)*(b(3)*c(2) + c(3)*b(2) )                                     ! h2 e3
    mtx(18) = ar2o3*(ahy*b(3)*b(2) + ahz*c(3)*c(2) )                             + bhx*4d0 ! h2 h3
    !E3
    mtx(19) = ar2o3*(aey*(b(3)*b(3) - b(1)*b(2)) + aez*(c(3)*c(3) - c(1)*c(2)) ) + bex*8d0 ! e3 e3
    mtx(20) = ar1o3*(cez-cey)*(b(1)*c(1) + b(2)*c(2) + b(3)*c(3) )                         ! e3 h3, not zero when anisotropic! 
    !H3 
    mtx(21) = ar2o3*(ahy*(b(3)*b(3) - b(1)*b(2)) + ahz*(c(3)*c(3) - c(1)*c(2)) ) + bhx*8d0 ! h3 h3
        
 
    end subroutine getStiff_bump

!==================================================================================================================================! 
!================================================================================================================  getStiff_lin_bump
!==================================================================================================================================! 
    subroutine getStiff_lin_bump(e,mtx)
!
! Routine to compute the element stiffness matrix for the coupled 2.5D
! CSEM system for B(u,v) where u is linear and v is bump space.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! April 2011:   Created
       
    implicit none
     
    integer, intent(in)    :: e
    complex(8),intent(out) :: mtx(21)   

! Local variables:  
    integer, dimension(3)   :: n
    real(8), dimension(3)   :: ye,ze,a,b,c
    real(8)                 :: sigx,sigy,sigz, ar, aro3 
    complex(8)              :: gammay2,gammaz2
    complex(8)              :: aey,aez,bex,cey,cez,ahy,ahz,bhx,chy,chz ! pde coefficients for e and h terms

!
! Get triaxial conductivity for element e
!
    call getSigs(e,sigx,sigy,sigz) 
    gammay2 = kx2-ommu*(sigy - iomeps)
    gammaz2 = kx2-ommu*(sigz - iomeps)
        
!
! Get Linear basis Coefficients:
!
    n    = mesh%emap(1:3,e)
    ye   = mesh%y(n)
    ze   = mesh%z(n)   
    call get_abc_coeffs(ye,ze,a,b,c,ar)
 
!
! 2.5D EM  coefficients for triaxial conductivity:
!
    aey = (sigy - iomeps) / gammay2
    aez = (sigz - iomeps) / gammaz2
    bex = (sigx - iomeps) * ar / 15d0
    cey =  ikx / gammay2  
    cez =  ikx / gammaz2

    ahz = -ommu / gammay2     ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    ahy = -ommu / gammaz2     ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    bhx = -ommu * ar / 15d0
    chy = -cez 
    chz = -cey 
   
    aro3 = 1d0/(3d0*ar)
    
!
! Elements of the stiffness matrix (upper triangular part only since its symmetric):
!
    ! E1
    mtx(1)  = -aro3*(aey*b(1)*b(1) + aez*c(1)*c(1))  + bex      ! e1 e1
    mtx(2)  =  aro3*(cey*b(1)*c(1) - cez*b(1)*c(1))             ! e1 h1, not zero when anisotropic!
    mtx(3)  = -aro3*(aey*b(1)*b(2) + aez*c(2)*c(1))  + bex*2d0  ! e1 e2
    mtx(4)  =  aro3*(cey*b(1)*c(2) - cez*b(2)*c(1))             ! e1 h2 
    mtx(5)  = -aro3*(aey*b(1)*b(3) + aez*c(3)*c(1))  + bex*2d0  ! e1 e3
    mtx(6)  =  aro3*(cey*b(1)*c(3) - cez*b(3)*c(1))             ! e1 h3
    ! H1     
    mtx(7)  = -aro3*(ahy*b(1)*b(1) + ahz*c(1)*c(1))  + bhx      ! h1 h1
    mtx(8)  =  aro3*(chy*b(1)*c(2) - chz*b(2)*c(1))             ! h1 e2 
    mtx(9)  = -aro3*(ahy*b(1)*b(2) + ahz*c(2)*c(1))  + bhx*2d0  ! h1 h2
    mtx(10) =  aro3*(chy*b(1)*c(3) - chz*b(3)*c(1))             ! h1 e3
    mtx(11) = -aro3*(ahy*b(1)*b(3) + ahz*c(3)*c(1))  + bhx*2d0  ! h1 h3
    ! E2
    mtx(12) = -aro3*(aey*b(2)*b(2) + aez*c(2)*c(2))  + bex      ! e2 e2
    mtx(13) =  aro3*(cey*b(2)*c(2) - cez*b(2)*c(2))             ! e2 h2, not zero when anisotropic! 
    mtx(14) = -aro3*(aey*b(2)*b(3) + aez*c(3)*c(2))  + bex*2d0  ! e2 e3 
    mtx(15) =  aro3*(cey*b(2)*c(3) - cez*b(3)*c(2))             ! e2 h3
    !H2
    mtx(16) = -aro3*(ahy*b(2)*b(2) + ahz*c(2)*c(2))  + bhx      ! h2 h2
    mtx(17) =  aro3*(chy*b(2)*c(3) - chz*b(3)*c(2))             ! h2 e3
    mtx(18) = -aro3*(ahy*b(2)*b(3) + ahz*c(3)*c(2))  + bhx*2d0  ! h2 h3
    !E3
    mtx(19) = -aro3*(aey*b(3)*b(3) + aez*c(3)*c(3))  + bex      ! e3 e3
    mtx(20) =  aro3*(cey*b(3)*c(3) - cez*b(3)*c(3))             ! e3 h3, not zero when anisotropic! 
    !H3 
    mtx(21) = -aro3*(ahy*b(3)*b(3) + ahz*c(3)*c(3))  + bhx      ! h3 h3
        
    end subroutine getStiff_lin_bump
     
 
!==================================================================================================================================! 
!=========================================================================================================================== gen_rhs 
!==================================================================================================================================! 
    subroutine gen_rhs 
!
! Generates the rhs vectors for CSEM and MT problems using optional 
! methods

    if (lMT) then ! MT using decoupled linear systems:           
        call gen_rhs_mt  
    else
        call gen_rhs_dipole
    endif    
 

         
    end subroutine gen_rhs          
!==================================================================================================================================! 
!==================================================================================================================== gen_rhs_dipole 
!==================================================================================================================================!     
    subroutine  gen_rhs_dipole
!                                                                        
! Routine to compute the RHS of the FEM system of equations using 
! a total field formulation for an electric dipole.
!
! Kerry Key
! Scripps Institution of Oceanography
!
    implicit none 
    
    integer                 :: i,n(3),e, iTx
    real(8)                 :: sigx,sigy,sigz
    complex(8)              :: gammay2,gammaz2        
    real(8), dimension(3)   :: l,dldy,dldz, ye,ze   
    complex(8),dimension(3) :: ee,hh  

    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_dipole ...'  

!
!  Initialize the RHS 
!
    rhs = (0d0,0d0)     
   
!
! Loop over transmitters:
!
    do iTx = 1,nTx
    
    !
    ! Get a single element that has the transmitter:
    !
        e = eTx(iTx)
        
        call getSigs(e,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
        gammay2 = kx2-ommu*(sigy - iomeps)
        gammaz2 = kx2-ommu*(sigz - iomeps)              
    
    !
    ! Set up the basis functions:
    !
        n   = mesh%emap(1:3,e)
        ye  = mesh%y(n)
        ze  = mesh%z(n)

    !
    ! Evaluate the basis functions and their derivatives at the source location:
    !
        call get_linear_basis(ye,ze,yTx(iTx),zTx(iTx),l,dldy,dldz)    
        
    !write(*,*)'gen_rhs_dipole:', iTx,jx(iTx),jy(iTx),jz(iTx),mx(iTx),my(iTx),mz(iTx)
    !
    ! Apply coefficients:
    !
        ee =    -l*jx(iTx) - ikx  * ( dldy*jy(iTx)/gammay2      + dldz*jz(iTx)/gammaz2 )  & 
           &               - ommu * ( dldy*mz(iTx)*(sigy - iomeps)/gammay2 - dldz*my(iTx)*(sigz - iomeps)/gammaz2 )  ! ommu 4 H src
        hh = -(                        - ommu * ( dldy*jz(iTx)/gammaz2      - dldz*jy(iTx)/gammay2 )  &
           &    -l*mx(iTx)*ommu -  ikx * ommu * ( dldy*my(iTx)/gammaz2      + dldz*mz(iTx)/gammay2 )  )  
        ! Note minus sign applied to hh in order to make the matrix symmetric

    !
    ! Insert these into the correct rhs location  
    !
        rhs(2*n-1, iTx) = ee  ! Ex  
        rhs(2*n  , iTx) = hh  ! Hx  
    
    enddo

! Lastly, insert the dirichlet boundary condition into the entries 
! of the RHS whose nodes lie on the border.  
!
! In the total field formulation this isn't necessary unless somebody 
! mistakenly places a Tx on an exterior node boundary, but then the solution
! is going to be incorrect anyway...
!
    do i=1,mesh%nnod
        if(bcnod(i).ne.0) then
          rhs(2*i-1,1:nTx) = (0d0, 0d0)
          rhs(2*i  ,1:nTx) = (0d0, 0d0)
        endif
    enddo    
    
    end subroutine gen_rhs_dipole
!==================================================================================================================================! 
!======================================================================================================================== gen_rhs_mt 
!==================================================================================================================================!     
    subroutine  gen_rhs_mt
!                                                                       
! Routine to compute the RHS of the FEM system of equations using 
! a total field formulation for a magnetotelluric plane wave.
! Here I insert the RHS for both TE and TM modes into the same vector.
! Note that later on the TE and TM modes are split apart and the 
! smaller uncoupled systems solved separately for better performance.
! See subroutine solve_primal_superlu_mt
!
! Original routine written for MARE2DMT by Chet Weiss and me. Heavily
! modified since.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! April 2011:       Added support for triaxial anisotropy (4 separate TE/TM and left/right models)  
! Sept 2009:     
!
    implicit none 
    
    integer                     :: i,e, n(3)
    integer,dimension(3,2)      :: imtxE, imtxH ! maps for BCs
    real(8)                     :: ymin, ymax, ediff, eta      
    complex(8)                  :: ehlTE(2), ehrTE(2),ehlTM(2), ehrTM(2), comboTE(2), comboTM(2)
    complex(8),dimension(21)    :: mtx   

    type(PW1D)                  :: Left1DTE,Left1DTM, Right1DTE,Right1DTM 


    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_mt ...'
    
!
! Set up mapping arrays for pulling out parts of 6 x 6 symmetric stiffness matrix stored as 21 upper triangular coeffs:
!
    imtxE(1,1) = 3   ! node 1 on boundary, select inner product with nodes 2 and 3
    imtxE(1,2) = 5
    imtxE(2,1) = 14  ! node 2 on boundary, select inner product with nodes 3 and 1
    imtxE(2,2) = 3 
    imtxE(3,1) = 5   ! node 3 on boundary, select inner product with nodes 1 and 2
    imtxE(3,2) = 14    
    imtxH(1,1) = 9
    imtxH(1,2) = 11
    imtxH(2,1) = 18
    imtxH(2,2) = 9 
    imtxH(3,1) = 11
    imtxH(3,2) = 18       
  
!
!  Initialize the RHS
!
    rhs = (0d0,0d0)     
! 
!  Compute the min and max x-values.
!
    ymin = minval(mesh%y)
    ymax = maxval(mesh%y)       
     
!
! Setup 1D MT models for BC's and get layer propagation coefficients
!
    call gen_1d_mt(Left1DTE , 'left ' , 'x')   ! setup 1D models for BC's
    call gen_1d_mt(Right1DTE, 'right', 'x')           
    call gen_1d_mt(Left1DTM , 'left ' , 'y')      
    call gen_1d_mt(Right1DTM, 'right', 'y')           

!
! Now call the PW1D initialize function to pre-compute the 1D recursions so that we can 
! get speedy solutions later on in gen_rhs_MT:
!
    call getPW1Dcoeffs(Left1DTE)  
    call getPW1Dcoeffs(Right1DTE) 
    call getPW1Dcoeffs(Left1DTM)  
    call getPW1Dcoeffs(Right1DTM)     
    
!
!  Loop over elements.
!
    do e=1,mesh%nele
        n = mesh%emap(1:3,e)
        
        !
        ! Loop over nodes and apply BC if node is on boundary:
        !
        do i=1,3
        
            if (bcnod(n(i)) == 1) then 
            
                ! Get stiffness matrix:
                call getStiff_lin(e, mtx)

                ! Get the 1D MT fields:
                call getPW1Dfields(Left1DTE , mesh%z(n(i)), ehlTE) 
                call getPW1Dfields(Right1DTE, mesh%z(n(i)), ehrTE) 
                call getPW1Dfields(Left1DTM , mesh%z(n(i)), ehlTM) 
                call getPW1Dfields(Right1DTM, mesh%z(n(i)), ehrTM) 
                
                ! Cosine taper between left and right sides:          
                eta = (mesh%y(n(i))-ymin)/(ymax-ymin)            
                ediff = (1d0-dcos(PI*eta))/2d0    
                comboTE = (ehrTE - ehlTE)*ediff + ehlTE
                comboTM = (ehrTM - ehlTM)*ediff + ehlTM
                
                ! Now the hard part is subtracting this off the RHS appropriately:
                ! remembering that there is no Ex-Hx coupling for 2D MT...
                ! TE Ex part:
                rhs(2*n(eps(i+1))-1,1) = rhs(2*n(eps(i+1))-1,1) - mtx(imtxE(i,1))*comboTE(1)
                rhs(2*n(eps(i+2))-1,1) = rhs(2*n(eps(i+2))-1,1) - mtx(imtxE(i,2))*comboTE(1)

                ! TM Hx part:
                rhs(2*n(eps(i+1))  ,1) = rhs(2*n(eps(i+1)) ,1 ) - mtx(imtxH(i,1))*comboTM(2)
                rhs(2*n(eps(i+2))  ,1) = rhs(2*n(eps(i+2)) ,1 ) - mtx(imtxH(i,2))*comboTM(2)              
    
            endif
        enddo  ! i = 1,3
    enddo  ! e = 1,nele

!
! Lastly, we need to set the boundary node values.
!
!  Out of sheer laziness, the nodes which lie on the boundary are still
!  present in the linear system.  However, the corresponding rows
!  of the matrix contain only '1' on the diagonal.  Thus, the RHS
!  for these nodes specifies the exact Dirichlet condition.  This
!  scheme has the advantage that when the FE system is solved, all
!  function values are present in the solution, not just those in
!  the mesh interior.
!
    do i = 1,mesh%nnod
        if (bcnod(i) == 1) then
            
            ! Get the 1D MT fields:
            call getPW1Dfields(Left1DTE , mesh%z(i), ehlTE) 
            call getPW1Dfields(Right1DTE, mesh%z(i), ehrTE)  
            call getPW1Dfields(Left1DTM , mesh%z(i), ehlTM) 
            call getPW1Dfields(Right1DTM, mesh%z(i), ehrTM)  
            
            ! Cosine taper between left and right sides:          
            eta = (mesh%y(i)-ymin)/(ymax-ymin)            
            ediff = (1d0-dcos(PI*eta))/2d0        
            comboTE = (ehrTE - ehlTE)*ediff + ehlTE   
            comboTM = (ehrTM - ehlTM)*ediff + ehlTM
            rhs(2*i - 1 ,1) = comboTE(1)
            rhs(2*i     ,1) = comboTM(2)
                
        endif
    enddo
    
    ! Deallocate 1D arrays:
    call PW1Dfinalize(Left1DTE)   
    call PW1Dfinalize(Right1DTE)
    call PW1Dfinalize(Left1DTM)   
    call PW1Dfinalize(Right1DTM)
    
    end subroutine gen_rhs_mt

!==================================================================================================================================! 
!==================================================================================================================== gen_rhs_bump_f 
!==================================================================================================================================!     
    subroutine  gen_rhs_bump_f
!
! Routine to compute F(v)  for all v in W_h (bump space)
!
! Kerry Key
! Scripps Institution of Oceanography
!
! April 2011    Added support for triaxial anisotropy    
! Sept 2009     Added support for 2DMT total field rhs, super difficult to code, see below...
! June 2008
!  

    implicit none

!
! Local variables:
!
    integer               :: n(3),e,iTx,j, iedge
    real(8), dimension(3) :: q,dqdy,dqdz,ye,ze  
    real(8)               :: sigx,sigy,sigz
    complex(8)            :: gammay2,gammaz2    
 
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_f ...'  
 
    do iTx = 1,nTx
    
        if (.not.lMT) then ! CSEM E dipole, x,y,z or yz 
        
        
        !
        ! Compute F(v) terms using a delta function for the source:
        !
         
        !
        ! Get a single element that has the transmitter:
        !           
            e = eTx(iTx)
            
            call getSigs(e,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
            gammay2 = kx2-ommu*(sigy - iomeps)
            gammaz2 = kx2-ommu*(sigz - iomeps)
 
        !
        ! Set up the basis functions:
        ! 
            n   = mesh%emap(1:3,e)
            ye  = mesh%y(n)
            ze  = mesh%z(n)
            
        !
        ! Evaluate bump basis at transmitter location:
        !
            call get_bump_basis(ye,ze,yTx(iTx),zTx(iTx),q,dqdy,dqdz)  
            
        !
        ! Apply coefficients:
        !
        
!            rhs_qE(edges(1:3,e),iTx) =  -q*jx(iTx) - ikx/gamma2 * ( dqdy*jy(iTx) + dqdz*jz(iTx) )
!            rhs_qH(edges(1:3,e),iTx) =             -ommu/gamma2 * ( dqdy*jz(iTx) - dqdz*jy(iTx) )  
!            rhs_qE(edges(1:3,e),iTx) =  -q*jx(iTx) - ikx * ( dqdy*jy(iTx)/gammay2 + dqdz*jz(iTx)/gammaz2 )
!            rhs_qH(edges(1:3,e),iTx) =             -ommu * ( dqdy*jz(iTx)/gammaz2 - dqdz*jy(iTx)/gammay2 )  

             rhs_q(2*edges(1:3,e)-1,iTx) =  rhs_q(2*edges(1:3,e)-1,iTx) &
                                         &    -q*jx(iTx) - ikx  * ( dqdy*jy(iTx)/gammay2      + dqdz*jz(iTx)/gammaz2 )  & 
                                         &               - ommu * ( dqdy*mz(iTx)*(sigy - iomeps)/gammay2 - dqdz*my(iTx)* &
                                         &                          (sigz - iomeps)/gammaz2 ) 
             rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx)  &
                                         &            -(            - ommu * ( dqdy*jz(iTx)/gammaz2      - dqdz*jy(iTx)/gammay2 ) &
                                         &    -q*mx(iTx)*ommu -  ikx * ommu *( dqdy*my(iTx)/gammaz2      + dqdz*mz(iTx)/gammay2 ) )
            
                
        else! MT  ! this is essentially the same as gen_rhs_mt, except that now we have a bump basis function
        ! Since for total field MT F = 0, all we have to do here is account for the boundary conditions.
        ! Tried using homogeneous Dirichlet and it works okay, but now we account for small curvature error along edges...
            call gen_rhs_bump_f_MT

        endif
        
    enddo

!
! Set 0 Dirichlet condition for CSEM:
!    
    if (.not.lMT) then
     do e = 1,mesh%nele
         do j = 1,3
            !
            ! ex and hx on boundaries are fixed to 0 so boundaries are fixed to 0 for error terms as well:
             if ( (bcnod( mesh%emap(eps(j+1),e) ) /= 0) .and. (bcnod( mesh%emap(eps(j+2),e) ) /= 0) ) then
                iedge = edges(j,e) ! Global edge number:
                rhs_q(2*iedge-1,1:nTx) = 0d0
                rhs_q(2*iedge  ,1:nTx) = 0d0
            endif
        enddo
    enddo
    endif
    
    
    end subroutine gen_rhs_bump_f
!==================================================================================================================================! 
!================================================================================================================= gen_rhs_bump_f_MT 
!==================================================================================================================================!     
    subroutine  gen_rhs_bump_f_MT 
    
    implicit none 
    
    integer                     :: i,e, n(3)
    integer,dimension(3,2)      :: imtxE, imtxH ! maps for BCs
    real(8)                     :: ymin, ymax, z1, z2, z0, y2, y1  
    complex(8)                  :: TE0(2),TE1(2),TE2(2), TM0(2),TM1(2),TM2(2), errTE, errTM
    complex(8),dimension(21)    :: mtx   
    type(PW1D)                  :: Left1DTE,Left1DTM, Right1DTE,Right1DTM 
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_f_MT ...'
    
!
! Set up mapping arrays for pulling out parts of 6 x 6 symmetric stiffness matrix stored as 21 upper triangular coeffs:
!
    imtxE(1,1) = 3   ! node 1 on boundary, select inner product with nodes 2 and 3
    imtxE(1,2) = 5
    imtxE(2,1) = 14  ! node 2 on boundary, select inner product with nodes 3 and 1
    imtxE(2,2) = 3 
    imtxE(3,1) = 5   ! node 3 on boundary, select inner product with nodes 1 and 2
    imtxE(3,2) = 14    
    imtxH(1,1) = 9
    imtxH(1,2) = 11
    imtxH(2,1) = 18
    imtxH(2,2) = 9 
    imtxH(3,1) = 11
    imtxH(3,2) = 18       
  
!
!  Initialize the RHS
!
   ! rhs = (0d0,0d0)     
! 
!  Compute the min and max x-values.
!
    ymin = minval(mesh%y)
    ymax = maxval(mesh%y)       
     
!
! Setup 1D MT models for BC's and get layer propagation coefficients
!
    call gen_1d_mt(Left1DTE , 'left ', 'x')   ! setup 1D models for BC's
    call gen_1d_mt(Right1DTE, 'right', 'x')           
    call gen_1d_mt(Left1DTM , 'left ', 'y')      
    call gen_1d_mt(Right1DTM, 'right', 'y')           

!
! Now call the PW1D initialize function to pre-compute the 1D recursions so that we can 
! get speedy solutions later on in gen_rhs_MT:
!
    call getPW1Dcoeffs(Left1DTE)  
    call getPW1Dcoeffs(Right1DTE) 
    call getPW1Dcoeffs(Left1DTM)  
    call getPW1Dcoeffs(Right1DTM)     
            
    do e = 1,mesh%nele

        n    = mesh%emap(1:3,e)
        
        do i = 1,3
    
        ! Apply boundary conditions for i-j pairs (i-i terms are set to 1 afterwards):
        if ( ( bcnod(n(eps(i+1))) /= 0 ) .and. ( bcnod(n(eps(i+2))) /= 0 )  ) then ! we got a boundary:
 
            ! Ge bounding node depths and the midpoint depth:
            z1 = mesh%z(n(eps(i+1)))
            z2 = mesh%z(n(eps(i+2)))

            
            z0 = (z1 + z2)/2d0
            
            y1 = mesh%y(n(eps(i+1)))
            y2 = mesh%y(n(eps(i+2)))
               
            if (y1 /= y2) cycle !write(*,*) 'error,y1/=y1!!!, element: ',e,i, y1, y2
            if (z1 == z2) cycle ! we only care about vertical boundaries
                     
             ! Get the 1D MT fields:   
            if (y1 == ymin) then ! left side:
                call getPW1Dfields(Left1DTE , z0, TE0)   
                call getPW1Dfields(Left1DTE , z1, TE1)    
                call getPW1Dfields(Left1DTE , z2, TE2)   
                call getPW1Dfields(Left1DTM , z0, TM0)   
                call getPW1Dfields(Left1DTM , z1, TM1)    
                call getPW1Dfields(Left1DTM , z2, TM2)       
            else ! right side
                call getPW1Dfields(Right1DTE , z0, TE0)   
                call getPW1Dfields(Right1DTE , z1, TE1)    
                call getPW1Dfields(Right1DTE , z2, TE2)   
                call getPW1Dfields(Right1DTM , z0, TM0)   
                call getPW1Dfields(Right1DTM , z1, TM1)    
                call getPW1Dfields(Right1DTM , z2, TM2)                 
            endif
            
            ! The error in the bump basis is then:
            errTE = TE0(1) - (TE1(1) + TE2(1) ) / 2d0               
            errTM = TM0(2) - (TM1(2) + TM2(2) ) / 2d0      
           ! write(*,*) abs(errTE) / abs(TE0(1)) , abs(errTM) / abs(TM0(2)) 
            
            ! Set RHS equal to this:
            rhs_q(2*edges(i,e)-1,1) = errTE
            rhs_q(2*edges(i,e)  ,1) = errTM
            
            ! Subtract off RHS of other edges:
            call getStiff_bump(e, mtx)
            
            ! Note that E-H coupling is 0 with anisotropic MT since kx=0, so only use decoupled terms for TE and TM modes: 
            rhs_q(2*edges(eps(i+1),e)-1,1) = rhs_q(2*edges(eps(i+1),e)-1,1) - mtx(imtxE(i,1))*errTE
            rhs_q(2*edges(eps(i+2),e)-1,1) = rhs_q(2*edges(eps(i+2),e)-1,1) - mtx(imtxE(i,2))*errTE
            
            rhs_q(2*edges(eps(i+1),e)  ,1) = rhs_q(2*edges(eps(i+1),e)  ,1) - mtx(imtxH(i,1))*errTM   
            rhs_q(2*edges(eps(i+2),e)  ,1) = rhs_q(2*edges(eps(i+2),e)  ,1) - mtx(imtxH(i,2))*errTM
            
        endif
       
       enddo !i                  
                               
    enddo ! e
 
    ! Deallocate 1D arrays:
    call PW1Dfinalize(Left1DTE)   
    call PW1Dfinalize(Right1DTE)
    call PW1Dfinalize(Left1DTM)   
    call PW1Dfinalize(Right1DTM)
        
    end subroutine gen_rhs_bump_f_MT
!==================================================================================================================================! 
!==================================================================================================================== gen_rhs_bump_b 
!==================================================================================================================================!     
    subroutine  gen_rhs_bump_b(imode)
!
! Routine to compute F(v) - B(u_h,v) for all v in W_h (bump space)
! F(v) contribution is handled in another subroutine.
! Also computes G(v) - B(v,wh) for dual error.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! April 2011    Modified for triaxial anisotropy, major clean up.
!               Now uses call to getStiff_lin_bump.
! June 2008     Created 
!  
    implicit none

    integer, intent(in) :: imode

!
! Local variables:
!
    integer         :: i, j, k, e, iedge, iTx, n(3), ict, ij(6,6) 
    complex(8)      :: ul(6), mtx(21)
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_b ...'   
!
! Create 6 x 6 lookup table for indices to symmetric upper triangular stiffness matrix:
!
    ict = 0
    do i = 1,6    
        do k = i,6
            ict = ict + 1
            ij(i , k ) = ict
            ij(k , i ) = ict
       enddo
    enddo
  
! 
! Create int( -B(u_h,v) ) where u_h is in V_h and v is in W_h   
!    
    
!
! Loop over elements:
!
    do e = 1,mesh%nele
    
        !
        ! Get B(u,v) lin-bump stiffness matrix for this element:
        !
        call getStiff_lin_bump(e,mtx)
 
        !  
        ! Compute B(u,v) for this element, subtract from RHS vector (rhs = F - B(u,v)):
        !
        ! B(u,v) =  MTX*u
        !
        n  = mesh%emap(1:3,e)        
        do iTx = 1,nTx
        
            do i = 1,3 ! nodes
                if (imode == 0 ) then ! use rhs from primal problem:
                    ul(2*i-1) = (rhs( 2*n(i)-1, iTx ))
                    ul(2*i)   = (rhs( 2*n(i)  , iTx ))
                else
                    ul(2*i-1) = (rhs_wh( 2*n(i)-1, iTx ))
                    ul(2*i)   = (rhs_wh( 2*n(i)  , iTx ))
                endif  
            enddo
            
            do j = 1,3 ! edges
                iedge = edges(j,e) ! Global edge number
                rhs_q(2*iedge-1, iTx) = rhs_q(2*iedge-1, iTx) - sum(mtx(ij(2*j-1,1:6))*ul)
                rhs_q(2*iedge  , iTx) = rhs_q(2*iedge  , iTx) - sum(mtx(ij(2*j  ,1:6))*ul)
            enddo              
            
        enddo ! iTx
 
    enddo ! e = 1,nele
 
    do e = 1,mesh%nele
        do j = 1,3
            !
            ! ex and hx on boundaries are fixed to 0, so same for error terms on boundaries:
             if ( (bcnod( mesh%emap(eps(j+1),e) ) /= 0) .and. (bcnod( mesh%emap(eps(j+2),e) ) /= 0) ) then
                iedge = edges(j,e) ! Global edge number:
                rhs_q(2*iedge-1,1:nTx) = 0d0
                rhs_q(2*iedge  ,1:nTx) = 0d0
            endif
        enddo
    enddo

    
    end subroutine gen_rhs_bump_b         
    
!==================================================================================================================================! 
!==================================================================================================================== gen_rhs_bump_g 
!==================================================================================================================================!  
    subroutine  gen_rhs_bump_g
!
! Computes G(v)  for dual error.
! Since E and H are decoupled on LHS (in bump basis), we generate
! separate RHS vectors so that the linear systems are smaller.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Created June 2008
!  

    implicit none

!
! Local variables:
!
    integer                  :: n(3),e, iTx, iRx,i
    real(8)                  :: smtx(3,3), mmtx(3,3),denomEx,denomHx,denomEg,denomHg
    real(8)                  :: wtEx,wtHx,wtEg,wtHg, ye(3),ze(3),a(3),b(3),c(3),area,yp,zp
    complex(8), dimension(3) :: fel,fhl,feq,fhq
    complex(8), dimension(3) :: exrn,etrn
    complex(8), dimension(3) :: hxrn,htrn
    complex(8)               :: exrd, etrd, hxrd, htrd
 
    complex(8) , dimension(:,:), allocatable:: rhs_q0
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_g ...'  
        
 
    allocate(rhs_q0(2*nedges,1:nTx)   )
    rhs_q0 = rhs_q
    
! Reset to zero:
    rhs_q = 0
     
!
! Loop over transmitters and receivers:

    do iTx = 1,nTx
        
        do iRx = 1,nRx
        
            call markCloud(iRx,iTx)
 
            denomEx = 0.
            denomHx = 0.
            denomEg = 0.
            denomHg = 0. 
        
        !
        ! Integrate denominator:
        ! 
            do e = 1,mesh%nele

                if (.not.lCloud(e)) cycle
            
                n    = mesh%emap(1:3,e)
                ye   = mesh%y(n)
                ze   = mesh%z(n)
                call get_abc_coeffs(ye,ze,a,b,c,area)      
 
                fel(1:3) = rhs(2*n-1,iTx)
                fhl(1:3) = rhs(2*n  ,iTx)         
            
                ! l_i l_j
                mmtx      =  area/12d0
                mmtx(1,1) =  area/6d0
                mmtx(2,2) =  area/6d0
                mmtx(3,3) =  area/6d0                  
              
                exrd = (sum(matmul(mmtx,fel)*conjg(fel)))
                hxrd = (sum(matmul(mmtx,fhl)*conjg(fhl)))
                       
            
                ! grad l_i dot grad l_j
        
                do i = 1,3
                    smtx(i,:) = (b*b(i)+c*c(i))/(4d0*area)
                enddo
 
                etrd = dble(sum(matmul(smtx,fel)*conjg(fel)))
                htrd = dble(sum(matmul(smtx,fhl)*conjg(fhl)))

                denomEx = denomEx + exrd + abs(ecutoff)**2*area
                denomHx = denomHx + hxrd + abs(hcutoff)**2*area

                denomEg = denomEg + etrd + abs(ecutoff)**2*area
                denomHg = denomHg + htrd + abs(hcutoff)**2*area
            
            enddo

    !
    ! Integrate numerator:
    !    
                            
            do e = 1,mesh%nele
        
                if (.not.lCloud(e)) cycle
   
                n    = mesh%emap(1:3,e)   
                ye   = mesh%y(n)
                ze   = mesh%z(n)
                call get_abc_coeffs(ye,ze,a,b,c,area)      
                yp = sum(ye)/3.
                zp = sum(ze)/3.

                fel(1:3) = rhs(2*n-1,iTx)
                fhl(1:3) = rhs(2*n  ,iTx)    
                feq(1:3) = rhs_q0( 2*edges(1:3,e)-1,iTx )
                fhq(1:3) = rhs_q0( 2*edges(1:3,e)  ,iTx )
 
 
            ! Compute weights: 
                call getWeights(iRx,iTx,yp,zp,area,ye,ze,fel,feq,fhl,fhq,wtEx,wtHx,wtEg,wtHg)
       
            ! Numerator:
                ! q_i * q_j:    
                mmtx      = area*4d0/45d0
                mmtx(1,1) = area*8d0/45d0
                mmtx(2,2) = area*8d0/45d0
                mmtx(3,3) = area*8d0/45d0      
                exrn = matmul(mmtx,(feq))
                hxrn = matmul(mmtx,(fhq))
        
                ! grad_q_i dot grad_q_j:
                do i = 1,3
                    smtx(i,:) = (b*b(i)+c(i)*c)*2d0/(3d0*area)
                enddo           
                smtx(1,1) = (b(1)*b(1)-b(2)*b(3) + c(1)*c(1) - c(2)*c(3) )*2d0/(3d0*area)
                smtx(2,2) = (b(2)*b(2)-b(3)*b(1) + c(2)*c(2) - c(3)*c(1) )*2d0/(3d0*area)
                smtx(3,3) = (b(3)*b(3)-b(1)*b(2) + c(3)*c(3) - c(1)*c(2) )*2d0/(3d0*area)
                etrn = matmul(smtx,feq)
                htrn = matmul(smtx,fhq)     
         
                exrn = abs(exrn)
                hxrn = abs(hxrn)
                etrn = abs(etrn)
                htrn = abs(htrn)
 
                !
                ! User selected functional:
                !           
                select case (idual_func)
            
                ! G0:  residual error only:
                case (0)
                         
                    rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + wtEx*( exrn/denomEx )
                    rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + wtHx*( hxrn/denomHx )
 
                ! G1:  gradient of residual error only:
                case(1)
                
                    rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + wtEg*( etrn/denomEg )
                    rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + wtHg*( htrn/denomHg )
                
                ! G2: residual error plus gradient of residual error:
                case(2)
   
                    rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + wtEx*( exrn/denomEx )
                    rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + wtHx*( hxrn/denomHx )
 
                    rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + wtEg*( etrn/denomEg )
                    rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + wtHg*( htrn/denomHg )
                                                       
                end select
  
            enddo
        enddo  
    enddo  
 
    deallocate (  rhs_q0) 
        
    end subroutine gen_rhs_bump_g
 
!==================================================================================================================================! 
!====================================================================================================================== gen_rhs_dual 
!==================================================================================================================================! 
    subroutine gen_rhs_dual
!
!  Computes RHS vector for dual system used for DRW
!                                                                    
! Kerry Key                                             
! Scripps Institution of Oceanography                   
! 
    
    implicit none

    integer                  :: i,iRx,e,n(3),iTx
    real(8)                  :: smtx(3,3), mmtx(3,3), denomEx, denomHx, denomEg, denomHg
    real(8)                  :: wtEx,wtHx,wtEg,wtHg,ye(3),ze(3), a(3),b(3),c(3),area,yp,zp 
    complex(8), dimension(3) :: fel,fhl,feq,fhq 
    complex(8), dimension(3) :: exrn,etrn
    complex(8), dimension(3) :: hxrn,htrn
    complex(8)               :: exrd, etrd, hxrd, htrd
    
 
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_dual ...'     
    
!
!  Initialize the RHS
!
    rhs_wh = 0d0
     
!
! Loop over transmitters and receivers:      
!    
    do iTx = 1,nTx
    
        do iRx = 1,nRx

        call markCloud(iRx,iTx)
        
        denomEx = 0.
        denomHx = 0.
        denomEg = 0.
        denomHg = 0. 
        
    !
    ! Integrate denominator:
    ! 
        do e = 1,mesh%nele

            if (.not.lCloud(e)) cycle
            
            n    = mesh%emap(1:3,e)
            ye   = mesh%y(n)
            ze   = mesh%z(n)
            call get_abc_coeffs(ye,ze,a,b,c,area)      
 
            fel(1:3) = rhs(2*n-1,iTx)  
            fhl(1:3) = rhs(2*n  ,iTx)    
            
            ! l_i l_j
            mmtx      =  area/12d0
            mmtx(1,1) =  area/6d0
            mmtx(2,2) =  area/6d0
            mmtx(3,3) =  area/6d0                  
              
            exrd = (sum(matmul(mmtx,fel)*conjg(fel)))
            hxrd = (sum(matmul(mmtx,fhl)*conjg(fhl)))                     
            
            ! grad l_i dot grad l_j
            do i = 1,3
                smtx(i,:) = (b*b(i)+c*c(i))/(4d0*area)
            enddo
 
            etrd = dble(sum(matmul(smtx,fel)*conjg(fel)))
            htrd = dble(sum(matmul(smtx,fhl)*conjg(fhl)))

            denomEx = denomEx + exrd + abs(ecutoff)**2*area
            denomHx = denomHx + hxrd + abs(hcutoff)**2*area


            denomEg = denomEg + etrd + abs(ecutoff)**2*area
            denomHg = denomHg + htrd + abs(hcutoff)**2*area

            
        enddo
        
    !
    ! Integrate numerator:
    !       
       do e = 1,mesh%nele
       
            if (.not.lCloud(e)) cycle
            
            !
            ! Set up the basis functions:
            ! 
            n    = mesh%emap(1:3,e)
            ye   = mesh%y(n)
            ze   = mesh%z(n)
            call get_abc_coeffs(ye,ze,a,b,c,area)      
 
            yp = sum(ye)/3.
            zp = sum(ze)/3.
 
            fel(1:3) = rhs(2*n-1,iTx)
            fhl(1:3) = rhs(2*n  ,iTx)
            feq(1:3) = rhs_q(2*edges(1:3,e)-1,iTx)
            fhq(1:3) = rhs_q(2*edges(1:3,e)  ,iTx)   
            
      ! Compute weights: 
            call getWeights(iRx,iTx,yp,zp,area,ye,ze,fel,feq,fhl,fhq,wtEx,wtHx,wtEg,wtHg)

        ! Numerator:   
       
            ! q_i * l_j:    
            mmtx      =  area*2d0/15d0
            mmtx(1,1) =  area/15d0   
            mmtx(2,2) =  area/15d0
            mmtx(3,3) =  area/15d0      
            
            exrn = matmul(mmtx,(feq))
            hxrn = matmul(mmtx,(fhq))
            
            ! grad_q_i dot grad_l_i:
            do i = 1,3
                smtx(i,:) = -(b*b(i)+c(i)*c)/(3d0*area)
            enddo            
            etrn = matmul(smtx,feq)
            htrn = matmul(smtx,fhq)
 
            exrn = abs(exrn)
            hxrn = abs(hxrn)
            etrn = abs(etrn)
            htrn = abs(htrn)
 
            !    
            ! User selected functional:
            !            
            select case (idual_func)
            
            ! G0:  residual error only:
            case (0)
                
                rhs_wh(2*n-1,iTx) = rhs_wh(2*n-1,iTx) + wtEx*( exrn/denomEx )
                rhs_wh(2*n  ,iTx) = rhs_wh(2*n  ,iTx) + wtHx*( hxrn/denomHx )
                       
            ! G1:  gradient of residual error only:
            case(1)

                rhs_wh(2*n-1,iTx) = rhs_wh(2*n-1,iTx) + wtEg*( etrn/denomEg ) 
                rhs_wh(2*n  ,iTx) = rhs_wh(2*n  ,iTx) + wtHg*( htrn/denomHg )              
                
            ! G2: residual error plus gradient of residual error:
            case(2)
 
                rhs_wh(2*n-1,iTx) = rhs_wh(2*n-1,iTx) + wtEx*( exrn/denomEx )
                rhs_wh(2*n  ,iTx) = rhs_wh(2*n  ,iTx) + wtHx*( hxrn/denomHx )
 
                rhs_wh(2*n-1,iTx) = rhs_wh(2*n-1,iTx) + wtEg*( etrn/denomEg ) 
                rhs_wh(2*n  ,iTx) = rhs_wh(2*n  ,iTx) + wtHg*( htrn/denomHg )   
      
            end select
 
        enddo    ! loop over elements
        
        enddo ! loop over receivers          
    enddo ! loop over transmitters
 
!
!  Set sourcing function to 0 on outer boundaries
!
    do i=1,mesh%nnod
        if (bcnod(i).ge.1) then
            rhs_wh(2*i-1,1:nTx) = 0d0
            rhs_wh(2*i  ,1:nTx) = 0d0
        endif
    enddo
               

    end subroutine gen_rhs_dual      
       
!==================================================================================================================================! 
!========================================================================================================================= markCloud 
!==================================================================================================================================! 
 
 subroutine markCloud(iRx,iTx)

    integer :: e,i,n(3), iTx, iRx
    real(8) :: r,ye(3),ze(3)
 
      
    lCloud = .false.
!
! Mark all elements that are close to the receiver:
! 
    do e = 1,mesh%nele
 
        n    = mesh%emap(1:3,e)
        ye   = mesh%y(n)
        ze   = mesh%z(n)    

        do i = 1,3
 
            r = sqrt( (yRx(iRx) - ye(i))**2 + (zRx(iRx) - ze(i))**2)
    
            if (r < rxCloud)  lCloud(e) = .true.
        
        enddo
    enddo            
 
 !
 ! Also mark all the element containing the receiver:
 !       
 
     if ( iDataMask(iRx,iTx) > 0 ) then
        lCloud(eRx(iRx)) = .true.
     endif   
 
 
 end subroutine markCloud
                
                
!==================================================================================================================================! 
!======================================================================================================================== getWeights 
!==================================================================================================================================! 
    subroutine getWeights(iRx,iTx,yp,zp,area,ye,ze,fel,feq,fhl,fhq,wtEx,wtHx,wtEg,wtHg)

!
! This subroutine sets some damping weights for the dual source functions used to guide adaptive refinement.
!
    
    integer, intent(in)                     :: iTx, iRx
    real(8), intent(in)                     :: area,yp,zp 
    real(8), dimension(3), intent(in)       :: ye,ze
    complex(8), dimension(3), intent(in)    :: fel,feq,fhl,fhq 
    real(8), intent(out)                    :: wtEx,wtHx,wtEg,wtHg 
  
    real(8)                     :: rTx,wrTx
    complex(8), dimension(3)    :: fiex,fihx,fieq,fihq 
       
 !
 ! First get a range dependent weight, to avoid over-refinement where receivers are located close to the source:
 !   
    wrTx = 1d0
        
    if (.not.lMT) then
        rTx  = sqrt( ( yp - yTx(iTx))**2 + ( zp - zTx(iTx))**2 )
        
        if (rTx < minRangeProtector)   wrTx = (rTx/minRangeProtector)**5
   
    endif           

    ! Data mask weights. If isite,iTx combo not needed, set weight to zero         
    if ( iDataMask(iRx,iTx) == 0 ) wrTx = 0d0  
    
    if ( area < minArea ) wrTx = 0d0
 
 
    ! map to component weights:
    wtEx = wrTx
    wtHx = wrTx
    wtEg = wrTx
    wtHg = wrTx

    !       
    ! Interpolate field to yp,zp location:        
    !
    call fget_interp_lin(yp,zp,fel,ye,ze,fiex)
    call fget_interp_lin(yp,zp,fhl,ye,ze,fihx)
    call fget_interp_bump(yp,zp,feq,ye,ze,fieq)
    call fget_interp_bump(yp,zp,fhq,ye,ze,fihq)        
            
    ! Set weights to zero if fields or errors are below cutoff thresholds:
    if ( maxval(abs(fel)) < ecutoff ) then
         wtEx = 0d0
         wtEg = 0d0
    endif
    if ( maxval(abs(fhl)) < hcutoff )  then
        wtHx = 0d0
        wtHg = 0d0
    endif
    if ( maxval(abs(feq)) < ecutoff )  then
        wtEx = 0d0
        wtEg = 0d0
    endif
    if ( maxval(abs(fhq)) < hcutoff )  then
        wtHx = 0d0
        wtHg = 0d0
    endif    
    

    if (lMT) then
  
      ! Set weights to zero if ratio of error/field is small:
        if ( abs(fieq(1))/( abs(fiex(1)) ) < 1d-6 )  then
            wtEx = 0d0
            wtEg = 0d0
        endif
        if ( abs(fihq(1))/( abs(fihx(1)) ) < 1d-6 )  then
            wtHx = 0d0
            wtHg = 0d0
        endif
            
    else ! CSEM data:
    
      ! Set weights to zero if ratio of error/field is small:  
        if ( abs(fieq(1))/( abs(fiex(1)) ) < 1d-6 )  then
            wtEx = 0d0
            wtEg = 0d0
        endif
        if ( abs(fihq(1))/( abs(fihx(1)) ) < 1d-6 )  then
            wtHx = 0d0
            wtHg = 0d0
        endif
             
        
    ! Catch for degenerate cases where E or H is unstable, usually for low f, low sigma and broadside transmitters:
        if ( (abs(jx(iTx))>0).or.(abs(jy(iTx))>0).or.(abs(jz(iTx))>0)) then
            if ( abs(fihx(1)) / abs(fiex(1)) < 1d0) then
                 wtHx = 0d0
                 wtHg = 0d0    
             endif            
        else
            if ( abs(fiex(1)) / abs(fihx(1)) < 1d0) then
                 wtEx = 0d0   
                 wtEg = 0d0
             endif
        endif
  
     
    endif 
          
    end subroutine getWeights
     
!==================================================================================================================================! 
!=========================================================================================================================== gen_lhs 
!==================================================================================================================================! 
    subroutine gen_lhs 
!
! Generates the LHS matrix for CSEM and MT problems using optional 
! methods
!
! Kerry Key
! Scripps Institution of Oceanography
!
 
    if (lMT) then ! MT using decoupled linear systems:
    
        kx = 0d0 ! Make sure kx=0    
        call gen_lhs_mt
    else
        call gen_lhs_csem
    endif
    
   

    end subroutine gen_lhs        

!==================================================================================================================================! 
!====================================================================================================================== gen_lhs_csem 
!==================================================================================================================================! 
      subroutine gen_lhs_csem 
! 
! Generates the lhs matrix for the coupled system in compressed sparse row format. 
!
! Kerry Key
! Scripps Institution of Oceanography
 
    
!
! Local parameters:
!
    integer                             :: i,j,e, n(3),v1,v2, indi, indj, indm, ind, i0,i1,ii,jj
    complex(8), dimension(21)           :: mtx
    integer,dimension(:), allocatable   :: nAdj
        
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_lhs ...'
 
!
! Allocate the sparse matrix arrays:
!    
    nnz = 4*mesh%nnod + 8*nedges        
    allocate ( val(nnz), col(nnz), irw(2*mesh%nnod+1) )
    allocate (   nAdj(mesh%nnod) )  ! temporary for CSR matrix
    
!
! Initialize a few things
!
    val  = (0d0,0d0)     ! Initialize the lhs matrix in compressed sparse row format:
    col = 0
    irw = 0
 
!
! Loop through all elements and count node adjacencies:
!     
    nAdj = 1  ! for the self adjacency
    do e = 1,mesh%nele
        n = mesh%emap(1:3,e)   
        do j = 1,3
            v1 = n(eps(j))
            v2 = n(eps(j+1))       
            if (v1 > v2) then ! use direction of edge to ensure we count each connection only once
                nAdj(v1) = nAdj(v1) + 1
                nAdj(v2) = nAdj(v2) + 1              
             elseif (mesh%neighborlist(eps(j+2),e) == -1) then ! this is a boundary edge with reverse order, add it:
                nAdj(v1) = nAdj(v1) + 1
                nAdj(v2) = nAdj(v2) + 1     
            endif
        enddo
    enddo   

!    
! Make irw row pointer array:
!
    irw(1) =  1
   
    do e =  1,mesh%nnod
        irw(2*e      )  = irw(2*e-1) + 2*nAdj(e)  ! Hx that goes with last Ex
        irw(2*(e+1)-1)  = irw(2*e  ) + 2*nAdj(e)  ! Ex for next row
    enddo
   
!
! Loop over all elements and insert entries into the sparse matrix:
!
    do e = 1,mesh%nele
    
        ! Get stiffness matrix for this element:
        call getStiff_lin(e, mtx)
            
        n = mesh%emap(1:3,e)    

        do i = 1,3
        
            indi = n(i)
                
            do j = 1,3
                
                indj = n(j)                
                
                indm = getUpperTriIndex(6,2*i-1,2*j-1)  ! Ex-Ex
                
                do ind = irw(2*indi-1),irw(2*indi)-1                 
                    if ( (col(ind) == 0) .or. (col(ind) == 2*indj-1) ) then  
                        val(ind) = val(ind) + mtx(indm)    
                        col(ind) = 2*indj-1
                        exit
                    endif 
                enddo
                
                indm = getUpperTriIndex(6,2*i-1,2*j  )  ! Ex-Hx
                
                do ind = irw(2*indi-1),irw(2*indi)-1
                    if ( (col(ind) == 0) .or. (col(ind) == 2*indj) ) then  
                        val(ind) = val(ind) + mtx(indm)    
                        col(ind) = 2*indj
                        exit
                    endif 
                enddo
 
                indm = getUpperTriIndex(6,2*i  ,2*j  )  ! Hx-Hx
                
                do ind = irw(2*indi),irw(2*indi+1)-1                    
                    if ( (col(ind) == 0) .or. (col(ind) == 2*indj  ) ) then  
                        val(ind) = val(ind) + mtx(indm)    
                        col(ind) = 2*indj  
                        exit
                    endif 
                enddo
                
                indm = getUpperTriIndex(6,2*i  ,2*j-1 )  ! Hx-Ex
                
                do ind = irw(2*indi),irw(2*indi+1)-1                  
                    if ( (col(ind) == 0) .or. (col(ind) == 2*indj-1  ) ) then  
                        val(ind) = val(ind) + mtx(indm)    
                        col(ind) = 2*indj-1  
                        exit
                    endif 
                enddo
            
 
            enddo ! j
            
        enddo  ! i  
        
   enddo  ! e      

!
! Apply boundary conditions:
! 
    do i = 1,2*mesh%nnod

        do i1 = irw(i),irw(i+1)-1
        
             j  = col(i1) 
             ii = floor(dble(i+1.)/2.)
             jj = floor(dble(j+1.)/2.) 
             
             if (   (  (bcnod(ii)>0).or.(bcnod(jj)>0) ) .and. (i /= j ) ) then
                val(i1) = 0d0    
             elseif (  (bcnod(ii)>0 ) .and. (i == j) ) then
                val(i1) = 1d0   
                !write(*,*)  i,j,ii,jj,i1,val(i1)
             endif 
 
       enddo
    enddo
 
   !   
   ! Finish up by sorting each row by column number:
   !     
    do i = 1,2*mesh%nnod
        i0 = irw(i)
        i1 = irw(i+1)-1
        call quick_sort(col(i0:i1),zlist1=val(i0:i1)) 
    enddo
   
    deallocate (nAdj) 
 
         
    end subroutine gen_lhs_csem 
!==================================================================================================================================! 
!======================================================================================================================== gen_lhs_mt
!==================================================================================================================================! 
    subroutine gen_lhs_mt
!  
! Generates the lhs matrix in compressed sparse row format for the 
! uncoupled MT E and H systems. 
!
! Kerry Key
! Scripps Institution of Oceanography
! 
    
    implicit none
    
    integer                             :: nnz, e, indm, n(3), i, j, ind, i0, i1, indi, indj, v1, v2
    complex(8), dimension(21)           :: mtx
    integer,dimension(:), allocatable   :: nAdj 
       
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_lhs_mt...'

!
! Initialize a few things
!
    nnz =  mesh%nnod + 2*nedges  
    
    allocate ( rhs_lE(mesh%nnod,nTx), rhs_lH(mesh%nnod,nTx)  ) 
    allocate ( val_lE(nnz), col_lE(nnz), irw_lE(mesh%nnod+1) )
    allocate ( val_lH(nnz), col_lH(nnz), irw_lH(mesh%nnod+1) )
    
    allocate (   nAdj(mesh%nnod) )  ! temporary for CSR matrix
    
         
    val_lE = (0d0,0d0)     ! Initialize the lhs matrices in compressed sparse row format:
    col_lE = 0
    irw_lE = 0
    nnz_lE = nnz
    
    val_lH = (0d0,0d0)     ! Initialize the lhs matrices in compressed sparse row format:
    col_lH = 0
    irw_lH = 0
    nnz_lH = nnz  
 
!
! Loop through all elements and count node adjacencies:
!    
 
    nAdj = 1  ! for the self adjacency
    do e = 1,mesh%nele
        n = mesh%emap(1:3,e)   
        do j = 1,3
            v1 = n(eps(j))
            v2 = n(eps(j+1))       
            if (v1 > v2) then ! use direction of edge to ensure we count each connection only once
                nAdj(v1) = nAdj(v1) + 1
                nAdj(v2) = nAdj(v2) + 1              
             elseif (mesh%neighborlist(eps(j+2),e) == -1) then ! this is a boundary edge with reverse order, add it:
                nAdj(v1) = nAdj(v1) + 1
                nAdj(v2) = nAdj(v2) + 1     
            endif
        enddo
    enddo   
 

!    
! Make irw row pointer arrays:
!
    irw_lE(1) =  1
    irw_lH(1) =  1 
    do e =  1,mesh%nnod
         irw_lE(e+1)  = irw_lE(e) + nAdj(e)
         irw_lH(e+1)  = irw_lH(e) + nAdj(e)
    enddo
   
!
! Loop over all elements and insert entries into the sparse matrix:
!
    do e = 1,mesh%nele
    
        ! Get stiffness matrix for this element:
        call getStiff_lin(e, mtx)

        n = mesh%emap(1:3,e)    

        do i = 1,3
        
            indi = n(i)
                
            do j = 1,3
                
                indj = n(j)                
                
                ! Ex:
                indm = getUpperTriIndex(6,2*i-1,2*j-1) 
                do ind = irw_lE(indi),irw_lE(indi+1)-1
                    if ( (col_lE(ind) < 1) .or. (col_lE(ind) == indj) ) then  
                        val_lE(ind) = val_lE(ind) + mtx(indm)    
                        col_lE(ind) = indj
                        exit
                    endif 
                enddo
 
                
                ! Hx:
                indm = getUpperTriIndex(6,2*i,2*j) 
                do ind = irw_lH(indi),irw_lH(indi+1)-1
                    if ( (col_lH(ind) < 1) .or. (col_lH(ind) == indj) ) then  
                        val_lH(ind) = val_lH(ind) + mtx(indm)    
                        col_lH(ind) = indj
                        exit
                    endif 
                enddo
 
            enddo ! j
            
        enddo  ! i  
        
   enddo  ! e      

  
    !
    ! Apply boundary conditions:
    ! 
    do i = 1,mesh%nnod

        do i1 = irw_lE(i),irw_lE(i+1)-1
        
             j = col_lE(i1)

             if (   (  (bcnod(i)>0).or.(bcnod(j)>0) ) .and. (i /= j ) ) then
                val_lE(i1) = 0d0  
                val_lH(i1) = 0d0      
             elseif (  (bcnod(i)>0 ) .and. (i == j) ) then
                val_lE(i1) = 1d0  
                val_lH(i1) = 1d0     
             endif   
 
       enddo
    enddo    

   !
   ! Finish up by sorting each row by column number:
   !     
    do i = 1,mesh%nnod
        i0 = irw_lE(i)
        i1 = irw_lE(i+1)-1
        call quick_sort(col_lE(i0:i1),zlist1=val_lE(i0:i1))
        call quick_sort(col_lH(i0:i1),zlist1=val_lH(i0:i1)) 
    enddo
    
    deallocate (  nAdj)  
        
    end subroutine gen_lhs_mt

!==================================================================================================================================! 
!================================================================================================================== getUpperTriIndex
!==================================================================================================================================!  
    integer function getUpperTriIndex(n,i,j) result(ind)
    
    integer :: n,i,j
    
    if (i<=j) then
        ind  = ((i-1)*(2*n-i))/2 + j
    else
        ind  = ((j-1)*(2*n-j))/2 + i
    endif
    
    end function getUpperTriIndex

!==================================================================================================================================! 
!====================================================================================================================== gen_lhs_bump
!==================================================================================================================================! 
      subroutine gen_lhs_bump
! 
! Generates the lhs matrix for the bump space in compressed sparse row
! format.  Coupled anisotropic matrix is 2nedges x 2nedges.  
!
! Kerry Key
! Scripps Institution of Oceanography
 
!
! Local parameters for element stiffness comps:
!
    integer                             :: i,j,e, nnz, indm, indi,indj, n(3),  i1, i0, ind 
    integer, dimension(:),allocatable   :: bcedge, eAdj
    complex(8), dimension(21)           :: mtx
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_lhs_bump ...'   
 
!
! Allocate the sparse matrix arrays:
!    
    nnz = 4*nedges + 24*mesh%nele          
    allocate (  bcedge(2*nedges)   )  ! temporary for CSR matrix
    allocate ( col_q(nnz), val_q(nnz),irw_q(2*nedges+1)  )
    
!
! Initialize a few things
!
    val_q  = (0d0,0d0)     ! Initialize the lhs matrix in compressed sparse row format:
    col_q = 0
    irw_q = 0
    nnz_q = nnz
 
    bcedge = 0

  ! create flag for edges on boundary:
    do e = 1,mesh%nele
        n = mesh%emap(1:3,e)       
        do i = 1,3
            if ( (bcnod(n(eps(i+1))) > 0) .and. (bcnod(n(eps(i+2))) > 0) ) then
                bcedge(2*edges(i,e)-1) = 1
                bcedge(2*edges(i,e)  ) = 1
            endif
        enddo    
    enddo  
    
    
!
! Loop through all elements and count edge adjacencies:
!    
    allocate(eAdj(2*nedges))
    eAdj = 2  ! for the self connections
    do e = 1,mesh%nele
        do i = 1,3
            indi = edges(i,e)            
            eAdj(2*indi-1) = eAdj(2*indi-1) + 4  ! Ex edges
            eAdj(2*indi  ) = eAdj(2*indi  ) + 4  ! Hx edges
        enddo
     
    enddo

!    
! Make irw_q row pointer array:
!
    irw_q(1) =  1 
    do e = 1,2*nedges
         irw_q(e+1)  = irw_q(e) + eAdj(e)
    enddo
 
!    
! Loop over all elements and insert stiffness matrix components:
!
        
    do e = 1,mesh%nele
    
        ! Get stiffness matrix for this element:
        call getStiff_bump(e,mtx)     
       

        ! Loop over rows and columns and insert 2 x 2 entries:
            
        do i = 1,3  ! row
            
            indi = edges(i,e)
            
            do j = 1,3  ! column
                
                indj = edges(j,e)
                             
                indm = getUpperTriIndex(6,2*i-1,2*j-1)  ! Ex-Ex
                
                do ind = irw_q(2*indi-1),irw_q(2*indi)-1
                    
                    if ( (col_q(ind) == 0) .or. (col_q(ind) == 2*indj-1) ) then  
                        val_q(ind) = val_q(ind) + mtx(indm)    
                        col_q(ind) = 2*indj-1
                        exit
                    endif 
                enddo
                
                indm = getUpperTriIndex(6,2*i-1,2*j  )  ! Ex-Hx
                
                do ind = irw_q(2*indi-1),irw_q(2*indi)-1
                    
                    if ( (col_q(ind) == 0) .or. (col_q(ind) == 2*indj) ) then  
                        val_q(ind) = val_q(ind) + mtx(indm)    
                        col_q(ind) = 2*indj
                        exit
                    endif 
                enddo
 
                indm = getUpperTriIndex(6,2*i  ,2*j  )  ! Hx-Hx
                
                do ind = irw_q(2*indi),irw_q(2*indi+1)-1
                    
                    if ( (col_q(ind) == 0) .or. (col_q(ind) == 2*indj  ) ) then  
                        val_q(ind) = val_q(ind) + mtx(indm)    
                        col_q(ind) = 2*indj  
                        exit
                    endif 
                enddo
                
                indm = getUpperTriIndex(6,2*i  ,2*j-1 )  ! Hx-Ex
                
                do ind = irw_q(2*indi),irw_q(2*indi+1)-1
                    
                    if ( (col_q(ind) == 0) .or. (col_q(ind) == 2*indj-1  ) ) then  
                        val_q(ind) = val_q(ind) + mtx(indm)    
                        col_q(ind) = 2*indj-1  
                        exit
                    endif 
                enddo
                                   
            enddo
            
        enddo
    enddo

    
    !
    ! Apply boundary conditions:
    !
    do i = 1,2*nedges

        do i1 = irw_q(i),irw_q(i+1)-1
        
             j = col_q(i1)
             
             if ( ( (bcedge(i)>0) .or. (bcedge(j)>0) ) .and. (i /= j ) ) then
                val_q(i1) = 0.    
             elseif ( (bcedge(i)>0 ) .and. (i == j) ) then
                val_q(i1) = 1.   
             endif   
 
       enddo
    enddo
 
 
   !   
   ! Finish up by sorting each row by column number:
   !     
 
    do i = 1,2*nedges
        i0 = irw_q(i)
        i1 = irw_q(i+1)-1
        call quick_sort(col_q(i0:i1),zlist1=val_q(i0:i1)) 
    enddo
 
    deallocate (bcedge,eAdj)
 
         
    end subroutine gen_lhs_bump
!==================================================================================================================================! 
!================================================================================================================= gen_lhs_bump_diag
!==================================================================================================================================! 
      subroutine gen_lhs_bump_diag
! 
! *** testing version that only stores diagonal of matrix***
!
! Generates the lhs matrix for the bump space in compressed sparse row
! format.  Coupled anisotropic matrix is 2nedges x 2nedges.  
!
! Kerry Key
! Scripps Institution of Oceanography
 
!
! Local parameters for element stiffness comps:
!
    integer                             :: i,e, indm, indi, n(3) 
    complex(8), dimension(21)           :: mtx
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_lhs_bump_diag ...'   
 
!
! Allocate array to store diagonal:
!
    allocate ( val_q(2*nedges)  )
    
!
! Initialize 
!
    val_q  = 0d0     ! Initialize the lhs matrix in compressed sparse row format: 
 
!    
! Loop over all elements and insert stiffness matrix components:
!
        
    do e = 1,mesh%nele
    
        ! Get stiffness matrix for this element:
        call getStiff_bump(e,mtx)     
       
        ! Loop over rows and columns and insert 2 x 2 entries:
            
        do i = 1,3  ! row
            
            indi = edges(i,e)
            
            indm = getUpperTriIndex(6,2*i-1,2*i-1)              ! Ex-Ex
            val_q(2*indi-1) = val_q(2*indi-1) + mtx(indm)       ! Ex-Ex 
            
            indm = getUpperTriIndex(6,2*i  ,2*i  )              ! Hx-Hx
            val_q(2*indi  ) = val_q(2*indi  ) + mtx(indm)       ! Hx-Hx
             
        enddo
 
    enddo
 
    !
    ! Apply boundary conditions:
    !
    do e = 1,mesh%nele
        n = mesh%emap(1:3,e)       
        do i = 1,3
            if ( (bcnod(n(eps(i+1))) > 0) .and. (bcnod(n(eps(i+2))) > 0) ) then
                val_q(2*edges(i,e)-1) = 1.   
                val_q(2*edges(i,e)  ) = 1.   
            endif
        enddo    
    enddo  
   
    end subroutine gen_lhs_bump_diag

!==================================================================================================================================! 
!======================================================================================================================= fill_solmtx
!==================================================================================================================================! 
    subroutine fill_solmtx
!
! Computes the kx domain fields at the receivers for the current wavenumber
! using a total field formulation.
!
! Kerry Key
! Scripps Institution of Oceanography
!

    implicit none
    
  ! local variables  
    real(8)         :: yp,zp
 
    integer         :: isite,e, iTx
    complex(8)      :: exx,hxx,dexdy,dexdz,dhxdy,dhxdz   
    
    real(8)               :: sigx,sigy,sigz
    complex(8)            :: gammay2,gammaz2    
    
  ! FE Nodal quantities:
    integer,   dimension(3) :: n 
    complex(8),dimension(3) :: fex,fhx,fiex,fihx ,ge,gh,fige,figh
    real(8),   dimension(3) :: ye,ze  

    if (lprintDebug_em2dkx) write(*,*) myID,': fill_solmtx ...'  
        
!
! Loop over the sites and fill in the kx domain solutions:
!
    do isite=1,nRx

        yp = yRx(isite)
        zp = zRx(isite)
    
!
! Get the site element index
!
        e = eRx(isite)
        
        call getSigs(e,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
        gammay2 = kx2-ommu*(sigy - iomeps)
        gammaz2 = kx2-ommu*(sigz - iomeps)         
       
!
! Extract solution and gradient vectors:
!
       ! nodes
        n    = mesh%emap(1:3,e)
        ye   = mesh%y(n)
        ze   = mesh%z(n)     
        
        do iTx = 1,nTx
        
            ! Get linear fields:
            fex = rhs(2*n-1,iTx) 
            fhx = rhs(2*n  ,iTx)     
            
            ! Get bump fields: 
            ge = rhs_q(2*edges(1:3,e)-1,iTx)
            gh = rhs_q(2*edges(1:3,e)  ,iTx)
               
    !       
    ! Interpolate field to site location:        
    !
            call fget_interp_lin(yp,zp,fex,ye,ze,fiex)
            call fget_interp_lin(yp,zp,fhx,ye,ze,fihx)
    
    !
    ! Set the fields and gradients:
    !    
            if (lUseBumpFields) then
                ! Interpolate error bump terms to site location: 
                call fget_interp_bump(yp,zp,ge,ye,ze,fige)
                call fget_interp_bump(yp,zp,gh,ye,ze,figh)
            else
                fige = 0d0
                figh = 0d0
            endif
 
         ! 
         ! Linear + bump interpolated fields:
         !       
            exx   = fiex(1) + fige(1)
            dexdy = fiex(2) + fige(2)
            dexdz = fiex(3) + fige(3)
            hxx   = fihx(1) + figh(1)
            dhxdy = fihx(2) + figh(2)
            dhxdz = fihx(3) + figh(3)  

    !
    ! Compute the total fields at current wavenumber, assuming Rx is not located at location of point sources, so those terms 
    ! ignored here
    !
    ! Note the conversion to single precision here:
          ! Ex:
            ex_kx(isite,iTx) = cmplx( exx )
          
          ! Ey
            ey_kx(isite,iTx) = cmplx( (-ommu*dhxdz - ikx*dexdy )/gammay2 )                  
        
          ! Ez         
            ez_kx(isite,iTx) = cmplx( ( ommu*dhxdy - ikx*dexdz )/gammaz2 ) 
                        
          ! Hx
            hx_kx(isite,iTx) = cmplx( hxx )              
    
          ! Hy                                  
            hy_kx(isite,iTx) = cmplx( (-ikx*dhxdy - (sigz - iomeps)*dexdz )/gammaz2 )                   
                                                        
          ! Hz    
            hz_kx(isite,iTx) = cmplx( (-ikx*dhxdz + (sigy - iomeps)*dexdy )/gammay2 )      
                
               
        enddo ! iTx

    enddo  ! loop over nRx    
    
    end subroutine fill_solmtx
 
!==================================================================================================================================! 
!=============================================================================================================== comp_drw_error_proj 
!==================================================================================================================================!     
    subroutine  comp_drw_error_proj( )
!
! Routine to compute  | F(delta_h) - B(u_h,delta_h) |_{\tau}
! where u_h is linear FE solution and delta_h is dual error in quadratic 
! bump space.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Updates:
! Dec 2009  Now uses projection method and quadrature so that 
!           multiple refinement iterations can be done with a
!           single error estimate computation
! Fall 2008
! June 2008
!  
    implicit none
 
    
    integer                 :: nele, e, n(3), iTx,k
    real(8)                 :: sigx,sigy,sigz
    complex(8)              :: gammay2,gammaz2        
    complex(8)              :: aey,aez,bex,cey,cez,ahy,ahz,bhx,chy,chz ! pde coefficients for e and h terms
    real(8)                 :: areanew, erre,errh

    real(8), dimension(3)                :: ye,ze,yeq,zeq,a,b,c 
    complex(8), dimension(3)             :: ex,hx,deltaE,deltaH,feli,fhli,feqi,fhqi 
    
    real(8), dimension(:), allocatable    :: yp,zp
    integer, dimension(:), allocatable    :: inme
    complex(8), dimension(:), allocatable :: ee,hh 

    real(8), dimension(:), allocatable :: sigTemp
    
    complex(8), dimension(4) :: quade,quadh 
    real(8), dimension(3)    :: bary
  
    integer :: nqp 
    real(8), dimension(4,3) :: qp
    real(8), dimension(4)   :: qwt 
        
    real(8) :: yq,zq
 
    nqp = 4
    qp(1,:) = [ 1d0/3d0, 1d0/3d0, 1d0/3d0 ]
    qp(2,:) = [ 3d0/5d0, 1d0/5d0, 1d0/5d0 ]
    qp(3,:) = [ 1d0/5d0, 3d0/5d0, 1d0/5d0 ]
    qp(4,:) = [ 1d0/5d0, 1d0/5d0, 3d0/5d0 ]  
    qwt = [-0.5625d0, 0.520833333333333d0, 0.520833333333333d0, 0.520833333333333d0]

!    nqp = 1
!    qp(1,:) = [ 1d0/3d0, 1d0/3d0, 1d0/3d0 ]
!    qwt = [1d0,0d0,0d0,0d0]
            
    if (lprintDebug_em2dkx) write(*,*) myID,': comp_drw_error_proj ...'  
       
    !
    ! Get list of elements in old mesh that contain centroids of new elements:
    !   
    nele = newmesh%nele
    
    allocate( yp(nele), zp(nele), inme(nele))
 
    yp = ( newmesh%y(newmesh%emap(1,:)) +  newmesh%y(newmesh%emap(2,:)) +  newmesh%y(newmesh%emap(3,:)) ) /3d0
    zp = ( newmesh%z(newmesh%emap(1,:)) +  newmesh%z(newmesh%emap(2,:)) +  newmesh%z(newmesh%emap(3,:)) ) /3d0
    
    allocate (sigTemp(mesh%nele))
    do e = 1,mesh%nele
        call getSigs(e,sigx,sigy,sigz)
        sigTemp(e) = (sigx+sigy+sigz)/3d0
    enddo
    
    call findElement(tree, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri,  &
    & sigTemp, newmesh%nele,yp, zp,inme)
    
    deallocate(sigTemp)
   
    allocate( ee(nTx),hh(nTx)  )
 
!  
! Loop over number of elements in newmesh:
!   
    do e = 1,nele 
            
        !
        ! Generate the element coefficients from the old mesh:
        !        
        call getSigs(inme(e),sigx,sigy,sigz) ! gets anisotropic conductivity components for element e
        gammay2 = kx2-ommu*(sigy - iomeps)
        gammaz2 = kx2-ommu*(sigz - iomeps) 
        
        !write(*,*) yp(inme(e)),zp(inme(e)),sigx,sigy,sigz
        
        aey = (sigy - iomeps) / gammay2
        aez = (sigz - iomeps) / gammaz2
        bex = (sigx - iomeps)  
        cey =  ikx / gammay2  
        cez =  ikx / gammaz2
        
        ahy =  ommu / gammay2  
        ahz =  ommu / gammaz2  
        bhx =  ommu  
        chy =  ikx / gammay2
        chz =  ikx / gammaz2         
                             
        n  = newmesh%emap(1:3,e)
        ye = newmesh%y(n)        
        ze = newmesh%z(n)                   
        call get_abc_coeffs(ye,ze,a,b,c,areanew)  
        
        n  = mesh%emap(1:3,inme(e))
        yeq = mesh%y(n)
        zeq = mesh%z(n)                       

        do iTx = 1,nTx
 
            ex     = (rhs( 2*n-1 ,iTx)) ! linear
            hx     = (rhs( 2*n   ,iTx)) 
                                 
            deltaE = (rhs_q(2*edges(1:3,inme(e))-1,iTx))   ! bump 
            deltaH = (rhs_q(2*edges(1:3,inme(e))  ,iTx))     
 
            do k = 1,nqp
   
              ! Get quadrature point (in element of the new mesh):
                bary = qp(k,1:3)
                yq   = sum(bary*ye)
                zq   = sum(bary*ze)
           
                ! Get linear fields:             
                call fget_interp_lin(yq,zq,ex,yeq,zeq,feli)   
                call fget_interp_lin(yq,zq,hx,yeq,zeq,fhli)
                    
                ! Get bump fields:
                call fget_interp_bump(yq,zq,deltaE,yeq,zeq,feqi)
                call fget_interp_bump(yq,zq,deltaH,yeq,zeq,fhqi)                                     
      
                ! Evaluate form:
                quade(k) = aey*feli(2)*feqi(2) + aez*feli(3)*feqi(3) + bex*feli(1)*feqi(1) + cez*fhli(2)*feqi(3)-cey*fhli(3)*feqi(2)  
                quadh(k) = ahz*fhli(2)*fhqi(2) + ahy*fhli(3)*fhqi(3) + bhx*fhli(1)*fhqi(1) + chy*feli(2)*fhqi(3)-chz*feli(3)*fhqi(2)   
                                     
                quade(k) =  (quade(k))*qwt(k)
                quadh(k) =  (quadh(k))*qwt(k)
         
            enddo
           
            ee(iTx) = sum((quade))*areanew  
            hh(iTx) = sum((quadh))*areanew
         
        enddo ! iTx = 1,nTx
 
        erre =  maxval(sqrt(abs(ee)))    ! ignoring F(delta) term since point sources only contribute as delta functions
        errh =  maxval(sqrt(abs(hh))) 
 
        errnrm(e) = max(erre,errh) !  sqrt(erre**2 + errh**2)  
              
    enddo ! do e = 1,nele   in newmesh
 
 
    deallocate(ee,hh,yp,zp,inme)
    
    end subroutine comp_drw_error_proj
   
!==================================================================================================================================!
!======================================================================================================================== refineMesh
!==================================================================================================================================! 
    subroutine refineMesh()  
!
! Refines mesh associated with iTx, updates edges and site element lists.
!
!
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
! 

    implicit none
    
    integer               :: i,e
    real(8)               :: myarea

 
    integer, dimension(:),allocatable  ::   indx    

    if (lprintDebug_em2dkx) write(*,*) myID,': refineMesh...'   
    
!
! Compute area of each element and set to 1/2 if element needs refinement
!
    allocate(indx(newmesh%nele))
    call indexx(newmesh%nele, errnrm,indx)    
    
    allocate (newmesh%area(newmesh%nele))
    newmesh%area = -1d0 
 
     do i = newmesh%nele,nint(newmesh%nele*(1d0 - pct_refine/100d0)),-1
        e = indx(i)
    
        if (errnrm(e) > errortolerance/1d2) then  

            call getArea(e,myarea,1)
 
            if ( myarea > minArea ) then  ! in case refinement goes bonkers, limit smallest elements possible
                newmesh%area(e) = myarea/2d0
            endif

        endif
    enddo

!    
! Refine the mesh:
!
    call call_triangle(retricommand, newmesh )

    deallocate(indx)

    end subroutine refineMesh     
    
!==================================================================================================================================! 
!=========================================================================================================== factor_lhs_bump_umfpack 
!==================================================================================================================================!   
    subroutine factor_lhs_bump_umfpack()
!    
!  Written by Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
    
    implicit none    

    if (lprintDebug_em2dkx) write(*,*) myID,': factor_lhs_bump_umfpack ...' 
    
    call umfpack_zfactor(umfpack_q,val_q,col_q,irw_q)
    
    if (umfpack_q%error /= 0) then
        write(*,*) ' Error in umfpack_zfactor from calling in subr factor_lhs_bump_umfpack' 
        write(*,*) ' Error code = ',umfpack_q%error
        write(*,*) ' Stopping !!!'
        stop
    endif
    
    deallocate (col_q,val_q,irw_q)
    
    end subroutine factor_lhs_bump_umfpack
    
!==================================================================================================================================! 
!=========================================================================================================== factor_lhs_bump_superlu 
!==================================================================================================================================!   
    subroutine factor_lhs_bump_superlu()
!    
!  Written by Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
    
    implicit none    

    if (lprintDebug_em2dkx) write(*,*) myID,': factor_lhs_bump_superlu ...' 
    
    call superlu_zfactor(superlu_q,val_q,col_q,irw_q)
      
    if (superlu_q%error /= 0) then
        write(*,*) ' Error in superlu_zfactor from calling in subr factor_lhs_bump_superlu' 
        write(*,*) ' Error code = ',superlu_q%error
        write(*,*) ' Stopping !!!'
        stop
    endif
    
    deallocate (col_q,val_q,irw_q)
    
    end subroutine factor_lhs_bump_superlu
    
!==================================================================================================================================! 
!========================================================================================================================== write_Ab 
!==================================================================================================================================!    
    subroutine write_Ab()
!
! Writes A and b from Ax=b out to file. 
! Matrix A is saved as:  irow,jcol, real, imag 
!
    implicit none
    
    integer         :: i,j
    
    character(256)  :: filename, cadapt 
        
    write(cadapt,'(i6)') meshnumber  
    cadapt = adjustl(cadapt)       
        
        
    if (.not.lMT) then
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.lhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) nnz
        do i = 1,2*mesh%nnod
            do j = irw(i),irw(i+1)-1
                write(16,*) i, col(j), dble(val(j)), aimag(val(j)) 
            enddo
        enddo
        close(16)       
        
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.rhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) 2*mesh%nnod, nTx
        do j = 1,nTx
            do i = 1,2*mesh%nnod
                write(16,*)   dble(rhs(i,j)), aimag(rhs(i,j)) 
            enddo
        enddo
        close(16)                      
        
    else  ! MT using decoupled linear systems:
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.hx_lhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*)nnz_lH
        do i = 1,mesh%nnod
            do j = irw_lH(i),irw_lH(i+1)-1
                write(16,*) i, col_lH(j), dble(val_lH(j)), aimag(val_lH(j)) 
            enddo
        enddo
        close(16)           
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.hx_rhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) mesh%nnod
        do i = 1,mesh%nnod
            write(16,*) dble(rhs(2*i,1)), aimag(rhs(2*i,1)) 
        enddo
        close(16)    
            
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.ex_lhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) nnz_lE
        do i = 1,mesh%nnod
            do j = irw_lE(i),irw_lE(i+1)-1
                write(16,*) i, col_lE(j), dble(val_lE(j)), aimag(val_lE(j)) 
            enddo
        enddo
        close(16)           
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.ex_rhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) mesh%nnod
        do i = 1,mesh%nnod
            write(16,*) dble(rhs(2*i-1,1)), aimag(rhs(2*i-1,1)) 
        enddo
        close(16)    
                               
    end if
               
    end subroutine write_Ab
!==================================================================================================================================! 
!=========================================================================================================================== write_x
!==================================================================================================================================!    
    subroutine write_x()
!
! Writes x from Ax=b out to file. 
! Matrix A is saved as:  irow,jcol, real, imag 
!
    implicit none
    
    integer         :: i,j
    
    character(256)  :: filename, cadapt 
        
    write(cadapt,'(i6)') meshnumber  
    cadapt = adjustl(cadapt)       
        
        
    if (.not.lMT) then  ! CSEM dipole
    
        
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.hxex'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) 2*mesh%nnod, nTx
        do j = 1,nTx
            do i = 1,2*mesh%nnod
                write(16,*)   dble(rhs(i,j)), aimag(rhs(i,j)) 
            enddo
        enddo
        close(16)                      
        
    else ! MT using decoupled linear systems:
    
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.hx'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) mesh%nnod
        do i = 1,mesh%nnod
            write(16,*) dble(rhs(2*i,1)), aimag(rhs(2*i,1)) 
        enddo
        close(16)    
            
        
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.ex'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) mesh%nnod
        do i = 1,mesh%nnod
            write(16,*) dble(rhs(2*i-1,1)), aimag(rhs(2*i-1,1)) 
        enddo
        close(16)    
                               
    end if         
               
    end subroutine write_x
                    
!==================================================================================================================================! 
!====================================================================================================================== solve_primal 
!==================================================================================================================================! 
    subroutine solve_primal 
!
! Solves the linear system AX=B for CSEM and MT problems using optional 
! methods
!
! Kerry Key
! Scripps Institution of Oceanography
!
    if (lSaveLinearSystem) call write_Ab
             

    if (.not.lMT) then ! CSEM dipole
    
 
        select case (trim(linearSolver))
        case ('superlu')
            call solve_primal_superlu
        case ('umfpack')
            call solve_primal_umfpack      
        end select
         
        ! Matrix has been factored, so we can deallocate the CSR arrays:
        deallocate(val,col,irw)       
        
    else ! MT using decoupled linear systems:
    
        select case (trim(linearSolver)) 
        case ('superlu')
            call solve_primal_superlu_mt
        case ('umfpack')
            call solve_primal_umfpack_mt
        end select  
        
        ! Matrices have been factored, so we can deallocate the CSR arrays:
        deallocate(val_lE,col_lE,irw_lE)  
        deallocate(val_lH,col_lH,irw_lH)  
                      
    end if

    if (lSaveLinearSystem) call write_x
                 
    end subroutine solve_primal   
                    
!==================================================================================================================================! 
!============================================================================================================== solve_primal_umfpack 
!==================================================================================================================================!    
    subroutine solve_primal_umfpack()
!
! Solves AX=B for all vectors in matrix B=rhs=[rhs1, rhs2,...].  
! The lhs A is passed in in compressed sparse row format.
! umfpack takes the lhs and factors it, then  solves for all vectors 
! in X corresponding to vectors in the rhs. X is returned in rhs to
! conserve memory.  I leave the umfpack factors in memory since the error
! estimator will use these for solving the dual problem later on.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! May 2011: Updated to use my new umfpack_solver module
!

    implicit none
    
    complex(8) , dimension(:), allocatable :: x
    integer :: j
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_primal_umfpack...'
              
!
! Factor the matrix first:
!
    call umfpack_zfactor(umfpack_p,val,col,irw)
    
    if (umfpack_p%error /= 0) then
        write(*,*) ' Error in umfpack_zfactor from calling in subr  solve_primal_umfpack()' 
        write(*,*) ' Error code = ',umfpack_p%error
        write(*,*) ' Stopping !!!'
        stop
    endif
    
!
! Now solve  AX=B for all column vectors in B:
!
    allocate ( x(size(rhs,1)) )
    
    do j = 1,size(rhs,2)   
        
        call umfpack_zsolve (umfpack_p,rhs(:,j),x )
        
        if (umfpack_p%error /= 0) then
            write(*,*) ' Error in umfpack_zsolve from calling in subr solve_primal_umfpack()' 
            write(*,*) ' Error code = ',umfpack_p%error
            write(*,*) ' Stopping !!!'
            stop
        endif
        
        rhs(:,j) = x
    enddo
    deallocate(x)
    
    end subroutine solve_primal_umfpack
    
    
!==================================================================================================================================! 
!============================================================================================================== solve_primal_superlu 
!==================================================================================================================================!    
    subroutine solve_primal_superlu()
!
! Solves AX=B for all vectors in matrix B=rhs=[rhs1, rhs2,...].  
! The lhs A is passed in in compressed sparse row format.
! SuperLU takes the lhs and factors it, then  solves for all vectors 
! in X corresponding to vectors in the rhs. X is returned in rhs to
! conserve memory.  I leave the SuperLU factors in memory since the error
! estimator will use these for solving the dual problem later on.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! May 2011: Updated to use my new superlu_solver module
!

    implicit none

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_primal_superlu ...'
          
!
! Factor the matrix first:
!
    call superlu_zfactor(superlu_p,val,col,irw)
    
    if (superlu_p%error /= 0) then
        write(*,*) ' Error in superlu_zfactor from calling in subr  solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p%error
        write(*,*) ' Stopping !!!'
        stop
    endif
!
! Now solve  AX=B for all column vectors in B:
!
    call superlu_zsolve (superlu_p,rhs )
    
    if (superlu_p%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p%error
        write(*,*) ' Stopping !!!'
        stop
    endif
 
    end subroutine solve_primal_superlu
    
!==================================================================================================================================! 
!=========================================================================================================== solve_primal_umfpack_mt 
!==================================================================================================================================!    
    subroutine solve_primal_umfpack_mt()
!
! Solves AX=B for all vectors in matrix B=rhs=[rhs1, rhs2,...].  
! The lhs A is passed in in compressed sparse row format.
! umfpack takes the lhs and factors it, then  solves for all vectors 
! in X corresponding to vectors in the rhs. X is returned in rhs to
! conserve memory.  I leave the umfpack factors in memory since the error
! estimator will use these for solving the dual problem later on.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! May 2011: Updated to use my new umfpack_solver module
!

    implicit none
    
    complex(8) , dimension(:), allocatable :: x

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_primal_umfpack_mt ...'
     
!
! Factor the matrix first:
!
    call umfpack_zfactor(umfpack_p_E,val_lE,col_lE,irw_lE)
    
    if (umfpack_p_E%error /= 0) then
        write(*,*) ' Error in umfpack_zfactor from calling in subr  solve_primal_umfpack_mt for E' 
        write(*,*) ' Error code = ',umfpack_p_E%error
        write(*,*) ' Stopping !!!'
        stop
    endif
    
    call umfpack_zfactor(umfpack_p_H,val_lH,col_lH,irw_lH)
    
    if (umfpack_p_H%error /= 0) then
        write(*,*) ' Error in umfpack_zfactor from calling in subr solve_primal_umfpack_mt for H' 
        write(*,*) ' Error code = ',umfpack_p_H%error
        write(*,*) ' Stopping !!!'
        stop
    endif    
!
! Now solve  AX=B for all column vectors in B:
!
    allocate ( x(size(rhs_lE,1)) )
    
    rhs_lE = rhs(1:2*mesh%nnod:2,1:1) 
            
    call umfpack_zsolve (umfpack_p_E,rhs_lE(:,1),x )
    
    if (umfpack_p_E%error /= 0) then
        write(*,*) ' Error in umfpack_zsolve trisolv from calling in subr solve_primal_umfpack_mt for E' 
        write(*,*) ' Error code = ',umfpack_p_E%error
        write(*,*) ' Stopping !!!'
        stop
    endif
    
    rhs_lE(:,1) = x

    rhs_lH = rhs(2:2*mesh%nnod:2,1:1) 
            
    call umfpack_zsolve (umfpack_p_H,rhs_lH(:,1),x )
    
    if (umfpack_p_H%error /= 0) then
        write(*,*) ' Error in umfpack_zsolve trisolv from calling in subr solve_primal_umfpack_mt for H' 
        write(*,*) ' Error code = ',umfpack_p_H%error
        write(*,*) ' Stopping !!!'
        stop
    endif
    
    rhs_lH(:,1) = x
     
    deallocate(x)
 
    
    ! Push these back into the full RHS vector: 
    rhs(1:2*mesh%nnod:2,1) = rhs_lE(:,1)
    rhs(2:2*mesh%nnod:2,1) = rhs_lH(:,1) 
    
       
    end subroutine solve_primal_umfpack_mt
        
!==================================================================================================================================! 
!=========================================================================================================== solve_primal_superlu_mt 
!==================================================================================================================================!    
    subroutine solve_primal_superlu_mt()
!
!  Solves AX=B for all vectors in matrix B=rhs=[rhs1, rhs2,...].  
! The lhs A is passed in in compressed sparse row format.
! SuperLU takes the lhs and factors it, then  solves for all vectors 
! in X corresponding to vectors in the rhs. X is returned in rhs to
! conserve memory.  I leave the SuperLU factors in memory since the error
! estimator will use these for solving the dual problem later on.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! May 2011: Updated to use my new superlu_solver module
!

    implicit none

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_primal_superlu_mt ...'
        
!
! Decoupled Matrix for MT problem:
!
! Factor  E matrix
!
    call superlu_zfactor(superlu_p_E,val_lE,col_lE,irw_lE)
    
    if (superlu_p_E%error /= 0) then
        write(*,*) ' Error in superlu_zfactor from calling in subr  solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p_E%error
        write(*,*) ' Stopping !!!'
        stop
    endif
!
! Factor  H matrix
!
    call superlu_zfactor(superlu_p_H,val_lH,col_lH,irw_lH)
    
    if (superlu_p_H%error /= 0) then
        write(*,*) ' Error in superlu_zfactor from calling in subr  solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p_H%error
        write(*,*) ' Stopping !!!'
        stop
    endif    
!
! Now solve  AX=B for all column vectors in B:
!
! E system:
!
    rhs_lE = rhs(1:2*mesh%nnod:2,1:1) 
    
    call superlu_zsolve (superlu_p_E,rhs_lE )
    
    if (superlu_p_E%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p_E%error
        write(*,*) ' Stopping !!!'
        stop
    endif
!
! H system:
!  
    rhs_lH = rhs(2:2*mesh%nnod:2,1:1) 
    
    call superlu_zsolve (superlu_p_H,rhs_lH )  
    if (superlu_p_H%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p_H%error
        write(*,*) ' Stopping !!!'
        stop
    endif 
   
   ! Push these back into the full RHS vector: 
    rhs(1:2*mesh%nnod:2,1) = rhs_lE(1:mesh%nnod,1)
    rhs(2:2*mesh%nnod:2,1) = rhs_lH(1:mesh%nnod,1) 
    
    
    end subroutine solve_primal_superlu_mt
    
!==================================================================================================================================! 
!========================================================================================================== solve_varepsilon_umfpack 
!==================================================================================================================================!   
    subroutine solve_varepsilon_umfpack
!
! Solves AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in umfpack_q derived type.
!   
!  Written by Kerry Key
!  Scripps Institution of Oceanography
!
! May 2011:         Updated to use my new umfpack_solver module  
! Summer 2008:      Implemented
 

    implicit none
    
    integer :: j 
    complex(8) , dimension(:), allocatable :: x

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_varepsilon_umfpack ...'

!
! Solves AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in umfpack_q derived type.
!
    allocate ( x(size(rhs_q,1)) )   
    
    do j = 1,size(rhs,2)   
        call umfpack_zsolve (umfpack_q,rhs_q(:,j),x)
    
        if (umfpack_q%error /= 0) then
            write(*,*) ' Error in umfpack_zsolve from calling in subr solve_varepsilon_umfpack()' 
            write(*,*) ' Error code = ',umfpack_q%error
            write(*,*) ' Stopping !!!'
            stop
        endif
        
        rhs_q(:,j) = x
    
    enddo
    
    deallocate(x)
    
    end subroutine solve_varepsilon_umfpack
!==================================================================================================================================! 
!========================================================================================================== solve_varepsilon_superlu
!==================================================================================================================================!   
    subroutine solve_varepsilon_superlu
!    
!  Written by Kerry Key
!  Scripps Institution of Oceanography
!
! May 2011:         Updated to use my new superlu_solver module  
! Summer 2008:      Implemented
 

    implicit none

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_varepsilon_superlu ...'
    
!
! Solves AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in superlu_q derived type.
!

    call superlu_zsolve (superlu_q,rhs_q)
    
    if (superlu_q%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_varepsilon_superlu()' 
        write(*,*) ' Error code = ',superlu_q%error
        write(*,*) ' Stopping !!!'
        stop
    endif

    end subroutine solve_varepsilon_superlu

!==================================================================================================================================! 
!============================================================================================================= solve_varepsilon_diag 
!==================================================================================================================================!   
    subroutine solve_varepsilon_diag
!
! *** A is approximated with a diagonal matrix ***
!
! Solves AX=B for all column vectors in B:
!
 
 
    implicit none
    
    integer :: j 
    complex(8) , dimension(:), allocatable :: x

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_varepsilon_umfpack ...'

!
! Solves AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in umfpack_q derived type.
!
    allocate ( x(size(rhs_q,1)) )   
    
    do j = 1,size(rhs,2)   
 
        rhs_q(:,j) =   rhs_q(:,j) / val_q
    
    enddo
    
    deallocate(x)
    
    end subroutine solve_varepsilon_diag
!==================================================================================================================================! 
!================================================================================================================ solve_dual_umfpack
!==================================================================================================================================!     
    subroutine solve_dual_umfpack
!
!  Solve dual problem for error weighting with umfpack.
!                                                                            
! Kerry Key                                             
! Scripps Institution of Oceanography                                     
!
    
    implicit none
    
    integer :: j
    complex(8) , dimension(:), allocatable :: x 

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_dual_umfpack ...'
    
!
! Solve AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in umfpack_p derived type.
!
    allocate ( x(size(rhs_wh,1)) )
    
    do j = 1,size(rhs_wh,2)   
         
        call umfpack_zsolve (umfpack_p,rhs_wh(:,j),x)
        
        if (umfpack_p%error /= 0) then
            write(*,*) ' Error in umfpack_zsolve from calling in subr solve_dual_umfpack()' 
            write(*,*) ' Error code = ',umfpack_p%error
            write(*,*) ' Stopping !!!'
            stop
        endif
        
        rhs_wh(:,j) = x
        
    enddo
    
    deallocate(x)
    
    end subroutine solve_dual_umfpack
    
!==================================================================================================================================! 
!================================================================================================================ solve_dual_superlu
!==================================================================================================================================!     
    subroutine solve_dual_superlu
!
!  Solve dual problem for error weighting with superlu.
!                                                                            
! Kerry Key                                             
! Scripps Institution of Oceanography                                     
!
    
    implicit none
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_dual_superlu ...'
     
!
! Solve AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in superlu_p derived type.
!
 
    call superlu_zsolve (superlu_p,rhs_wh)
    
    if (superlu_p%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_dual_superlu()' 
        write(*,*) ' Error code = ',superlu_p%error
        write(*,*) ' Stopping !!!'
        stop
    endif

    end subroutine solve_dual_superlu
!==================================================================================================================================! 
!============================================================================================================= solve_dual_superlu_mt
!==================================================================================================================================!     
    subroutine solve_dual_superlu_mt
!
! Solve dual problem for error weighting with superlu.
!                                                                            
! Kerry Key                                             
! Scripps Institution of Oceanography                                     
!
    
    implicit none
    integer :: j  

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_dual_superlu_mt ...'
  
!
! Solve AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in superlu_p derived type.
!
    do j = 1,size(rhs_wh,2)  
    
        rhs_lE(:,1) = rhs_wh(1:2*mesh%nnod:2,j) 
        call superlu_zsolve (superlu_p_E,rhs_lE)
        
        if (superlu_p_E%error /= 0) then
            write(*,*) ' Error in superlu_zsolve from calling in subr solve_dual_superlu()' 
            write(*,*) ' Error code = ',superlu_p_E%error
            write(*,*) ' Stopping !!!'
            stop
        endif
    
        rhs_lH(:,1) = rhs_wh(2:2*mesh%nnod:2,j) 
        call superlu_zsolve (superlu_p_H,rhs_lH)
        
        if (superlu_p_H%error /= 0) then
            write(*,*) ' Error in superlu_zsolve from calling in subr solve_dual_superlu()' 
            write(*,*) ' Error code = ',superlu_p_H%error
            write(*,*) ' Stopping !!!'
            stop
        endif
    
     !
     ! Fill back into rhs_wh array
     !
        rhs_wh(1:2*mesh%nnod:2,j) = rhs_lE(:,1)
        rhs_wh(2:2*mesh%nnod:2,j) = rhs_lH(:,1)  
    
    enddo 
    
    end subroutine solve_dual_superlu_mt
!==================================================================================================================================! 
!============================================================================================================= solve_dual_umfpack_mt
!==================================================================================================================================!    
    subroutine solve_dual_umfpack_mt()
 

    implicit none
    
    complex(8) , dimension(:), allocatable :: x
    
    integer :: j
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_dual_umfpack_mt ...'
     
    allocate ( x(size(rhs_lE,1)) )
    
    do j = 1,size(rhs_wh,2)  
    
        rhs_lE(:,1) = rhs_wh(1:2*mesh%nnod:2,j) 
                
        call umfpack_zsolve (umfpack_p_E,rhs_lE(:,1),x )
        
        if (umfpack_p_E%error /= 0) then
            write(*,*) ' Error in umfpack_zsolve from calling in subr solve_primal_umfpack()' 
            write(*,*) ' Error code = ',umfpack_p_E%error
            write(*,*) ' Stopping !!!'
            stop
        endif
        
        rhs_lE(:,1) = x
    
        rhs_lH(:,1) = rhs_wh(2:2*mesh%nnod:2,j) 
                
        call umfpack_zsolve (umfpack_p_H,rhs_lH(:,1),x )
        
        if (umfpack_p_H%error /= 0) then
            write(*,*) ' Error in umfpack_zsolve from calling in subr solve_primal_umfpack()' 
            write(*,*) ' Error code = ',umfpack_p_H%error
            write(*,*) ' Stopping !!!'
            stop
        endif
        
        rhs_lH(:,1) = x

     
     !
     ! Fill back into rhs_wh array
     !
        rhs_wh(1:2*mesh%nnod:2,j) = rhs_lE(:,1)
        rhs_wh(2:2*mesh%nnod:2,j) = rhs_lH(:,1) 
       
    enddo
            
    deallocate(x)
    
    end subroutine solve_dual_umfpack_mt


    end module em2dkx_mod
    

