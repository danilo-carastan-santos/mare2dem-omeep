!
! This version of Occam has been rewritten from the ground up in the spring of 2012. It is based on the original Occam 
! algorithm written by Steven Constable, but includes many newer features added by David Myer and Kerry Key 
! for the Occam1DCSEM and Occam2DMT codes, as well many new features created specifically for MARE2DEM.
!
! Highlights:
! 
! - This rewrite has focused on creating an easily readable source code.
! - Only modern Fortran commands are used (for example, there are no GOTO statements).
! - The minimization algorithms now have escape hatches that result in vastly fewer forward calls than the original Occam algorithm. 
! - The logfile and terminal output texts are now identical.
! - Model parameters are bound using non-linear transformations.
! - Multithreaded LAPACK and BLAS libraries are used for the linear algebra operations.
! - The starting model response is now output to a file.
! - Parameters can have hard bounds implemented with non-linear transforms
! - Parameter ratios can be prejudiced. This is actually a prejudice against parameter differences but for log10(resistivity)
!   parameters this is equivalent to the resistivity ratios. This is useful for having an anisotropy ratio prejudice.
!
!                 This work was supported by: 
!
!            The Seafloor Electromagnetic Methods Consortium   
!               at Scripps Institution of Oceanography  
!
!         See this URL for a list of current SEMC funding sponsors:
!
!               http://marineemlab.ucsd.edu/semc.html
!
!    Copyright (C) 2012-2013
!    Kerry Key
!    Scripps Institution of Oceanography
!    University of California, San Diego
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
!==================================================================================================================================! 
!===================================================================================================================== Occam Module
!==================================================================================================================================! 
    module Occam

!#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
!    implicit none   ! since this is up top in the module, all contained subroutines are also implicit none by default.
!    include 'mpif.h'
!#else
    use mpi
    implicit none   ! since this is up top in the module, all contained subroutines are also implicit none by default.
!#endif 
    
    private         ! default is to make everything private to this module
 
! 
! But there are a few public subroutines that can be called externally:
! 
    public :: computeOccamIteration 
    public :: get_time_offset,openOccamLogFile, closeOccamLogFile, printOccamLog
    public :: deallocateOccam, transformToBound, transformToUnbound
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
    public :: mpi_worker_ScaLapackMM, mpi_worker_ScaLapackMS 
#endif
           
!-----------------------------------------------------------------------------------------------------------------------------------
! Public data:
!-----------------------------------------------------------------------------------------------------------------------------------
!
! These can be accessed outside Occam by placing a "use Occam" statement in your subroutine.
!    
    integer :: nMKL_numThreads = 8 ! adjust this for the number of cores on a given cluster node
                                   ! Used for Intel Math Kernel Library multithreaded LAPACK and BLAS
    
    integer, parameter, public :: RealPrec = kind(0d0)   ! Double precision:  RealPrec = kind(0d0) 
                                                         ! Single precision:  RealPrec = kind(0e0)  
    ! ***Leave RealPrec set to double precision for MARE2DEM or bad things will happen***
    
    integer, public        :: currentIteration = 0      ! Current Occam iteration number 
    integer, public        :: maxOccamIterations = 100  ! Just what is says.    
    integer,public         :: occamPrintLevel  = 1      ! 0: print nothing, 1: also print logfile statements to the terminal window.     
    real(RealPrec), public :: modelMu          = 5      ! input and output model mu, large default    
    real(RealPrec), public :: targetRMS                 ! target RMS misfit  
    real(RealPrec), public :: modelRoughness            ! output model roughness 
    real(RealPrec), public :: modelRMS                  ! output model RMS misfit  

    
    integer, public        :: convergenceFlag
    ! Definition:  
    !  0 - Successful iteration: Model RMS was lowered but target RMS has not been obtained
    !  1 - Successful iteration: Target RMS has been obtained and Occam smoothing is in progress
    !  2 - Normal convergence:   Target RMS has been obtained, model has been smoothed until the step size is below a threshold 
    !  3 - Unusual convergence:  A perfectly smooth model has been found. Occam is done.
    !  4 - Convergence problem:  Occam could not find a model with an RMS less than the last iteration
    !  5 - Convergence problem:  Occam found the target misfit, but the model roughness couldn't be decreased.
    !  6 - Convergence problem:  Maximum number of iterations achieved.

    ! Data and model parameter variables:
    integer, public                                     :: nd = 0     ! number of data
    real(RealPrec), dimension(:), allocatable, public   :: d          ! data array d(nd) 
    real(RealPrec), dimension(:), allocatable, public   :: sd         ! standard error array sd(nd)    
    integer, dimension(:,:), allocatable, public        :: dp         ! data parameters (nd x nDataParams) lookup table
                                                                      ! of parameter indices (frequencies, positions, etc)    
    real(RealPrec), dimension(:), allocatable, public   :: d_wt       ! Data weights applied with sd but not to misfit,
                                                                      ! used to balance the influence of different data types
                                                                      ! for example a few MT data with LOTS of csem, use data 
                                                                      ! weights to upweight the MT influence. 
                                                                     
    integer, public                                     :: nParams = 0  ! number of model parameters
    real(RealPrec), dimension(:), allocatable, public   :: pm         ! model parameter array pm(nParams)
    real(RealPrec), dimension(:), allocatable, public   :: dm         ! model response vector dm(nd)
    real(RealPrec), dimension(:,:), allocatable, public :: wj         ! model response Jacobian matrix
    real(RealPrec), dimension(:), allocatable, public   :: premod     ! model preference values premod(nParams)
    real(RealPrec), dimension(:), allocatable, public   :: prewts     ! weights for model preference prewts(nParams). Set to 0 for 
                                                                      ! no preference  
                                                                      
    integer, public                                     :: nDiff = 0  ! number of parameter difference preference values                                                               
    real(RealPrec), dimension(:), allocatable, public   :: preDiff    ! model parameter difference preference values
    real(RealPrec), dimension(:), allocatable, public   :: preDiffwts ! weights for model parameter difference. Set to 0 for 
                                                                      ! no preference.  preDiffwts(nParams)                                             
    integer,        dimension(:,:), allocatable, public :: ijDiff     ! nDiff x 2 array of i,j indices of model parameters to 
                                                                      ! difference. ijDiff(nParams,2). Set to 0 for no difference                                                                 
                                                                                                                                            
    integer, public                                     :: npenalty   ! number of rows in the penalty matrix                                                            
    real(RealPrec), dimension(:), allocatable, public   :: penaltywts ! vector of penalty weights to apply  
    integer,        dimension(:,:), allocatable, public :: ijpenalty  ! npenalty x 2 array of i,j indices of model parameters to 
                                                                      ! difference
                                                                      
    ! June 1 2012: new penalty matrix stored as CSR matrix:
    type PenaltyMatrix
        integer, public                                     :: nnz, nrows
        integer, dimension(:), allocatable, public          :: colind
        integer, dimension(:), allocatable, public          :: rowptr
        real(RealPrec), dimension(:), allocatable, public   :: val 
    end type PenaltyMatrix
    
    type(PenaltyMatrix), public                             :: pen          
                                                             
    logical, public                                         :: lMGS = .false.   ! option for minimum gradient support regularization   
    real(RealPrec), public                                  :: beta_mgs = 1d-2  ! not yet fully supported, so leave it .false.                                                               
    
    integer, public                                     :: npm_assoc = 0 ! Number of associate (pass through) parameters
    real(RealPrec), dimension(:), allocatable, public   :: pm_assoc      ! Associated parameter arrays pm_assoc(npm_assoc)
 
 
    ! David Myer's autoconverge feature (my implementation is still a work in progress):
    logical                                             :: lAutoConverge = .true. 
    real(RealPrec)                                      :: nAutoConvergeFactor = 1.5
 
   ! David Myer's large RMS decrease escape hatch. Minimization stops if RMS decreases significantly before the minimum is found, 
   ! so that the next iteration can then proceed with a new Jacobian. Iteration ends when rms < rmsThreshold*startingRMS 
    real(RealPrec), public                              :: rmsThreshold = 0  ! ( 0 <= rmsThreshold < 1), 0.7-0.9 works well
                                                                                 !              
    ! Option to slowly obtain the target misfit setting the target misfit to be max(targetRMS,startingRMS*rmsThreshold)
    ! so that Occam finds the smoothest model at each step. 
    logical, public                                     :: lConvergeSlowly  = .false. 

    
    ! Non-linear transformation for constraining parameters with upper and lower bounds:
    real(RealPrec), parameter, private                  :: bandPassfactor = 15.             ! factor in exp(bandPassFactor/(a-b)*x)
    character(80), public                               :: cBoundsTransform ='bandpass'     ! 'exponential' or 'bandpass'.  
    !
    ! Can be modified by an optional parameter in the iteration file. 
    ! 'exponential' - from Habashy and Abubakar (2004), Commer and Newman (2008)
    ! 'bandpass' - designed by K. Key in March 2011
    !
    real(RealPrec), dimension(:), allocatable, public   :: lowerBound,upperBound ! public so that it can be overridden by readModel
    logical, dimension(:), allocatable, public          :: lBoundMe

    integer, parameter, public                          :: ioUnitOccamLogFile = 22  
     
!-----------------------------------------------------------------------------------------------------------------------------------
! Private data:
!
! You can't (and don't need to) externally access anything below here
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Parameters controlling Occam's behavior:
!             
    real(RealPrec), parameter                   :: rmsTol = 0.01       ! Tolerance for hitting target RMS 
                                                                       ! iteration complete if: abs(RMS-targetRMS) < rmsTol    
    integer, parameter                          :: maxNumStepCuts = 8  ! don't cut the step-size in half more than this    
    real(RealPrec), parameter                   :: stepSizeTol    = 2.**(1-maxNumStepCuts) ! Stop if targetRMS obtained and stepSize is smaller than this    
    integer                                     :: numStepCuts    = 0
    
    integer                                     :: nParallelProc
    integer, parameter                          :: MinProcScaLAPACK = 12 
    
!
! Working arrays and variables:
! 
    logical                                     :: lTargetHit  ! true if input model is at target RMS
    real(RealPrec)                              :: startingRMS, startingRoughness                  
    real(RealPrec)                              :: stepCut ! current reduction in model update size
    real(RealPrec)                              :: stepSize      
    real(RealPrec), dimension(:,:), allocatable :: wjtwj   ! product of weighted jacobian matrix and its transpose
    real(RealPrec), dimension(:),   allocatable :: wjtwd   ! product of weighted jacobian matrix and weighted translated data
    real(RealPrec), dimension(:),   allocatable :: resid   ! model fit residual, normalized by sd
    real(RealPrec), dimension(:,:), allocatable :: amat    ! matrix to invert for model update
 
    real(RealPrec), dimension(:),   allocatable :: pm_b, dm_b, pm_assoc_b  ! *_b arrays keep track of optimal model and response
    real(RealPrec)                              :: mu_b, rms_b, rough_b 

    real(RealPrec), dimension(:),   allocatable :: pm_a, dm_a, pm_assoc_a  ! used during intercept search
    real(RealPrec)                              :: mu_a, rms_a, rough_a 
    real(RealPrec), dimension(:),   allocatable :: pm_c, dm_c, pm_assoc_c  
    real(RealPrec)                              :: mu_c, rms_c, rough_c 
    
    real(RealPrec), dimension(:),   allocatable :: pm_test                 ! used in misfitOfMu
    real(RealPrec)                              :: rough_test
 
    real(RealPrec)                              :: tracker1_rms, tracker1_mu
    real(RealPrec), parameter                   :: tracker_Dummy = 1d15    ! dummy initialization flag
    logical                                     :: lTrackerOn
               
    integer                                     :: numForwardCalls, numForwardCallsCumulative
    real(8)                                     :: timeIter0, timeIterEnd, timeOccamStart, timeOccamEnd    

    
    contains   ! Contained functions and subroutines are given below                                                               

 
!==================================================================================================================================! 
!=========================================================================================================== allocateOccamIteration
!==================================================================================================================================! 
    subroutine allocateOccamIteration
!
! Allocates arrays during an Occam iteration:
!     
    integer :: ierr
    character(128)  :: cStr   

    allocate ( amat(nParams,nParams), pm_b(nParams), pm_c(nParams), pm_test(nParams), dm_b(nd),  dm_c(nd), stat=ierr)    
    if (ierr .ne. 0) then
        write(cStr,'(a,i8)')  'Out of memory in Occam.  Too many free parameters: ', nParams 
        call printOccamLog(cStr)  
        stop 
    endif      
 

    if (npm_assoc > 0)  allocate (pm_assoc_b(npm_assoc), pm_assoc_c(npm_assoc), stat=ierr)         
    if (ierr .ne. 0) then
        write(cStr,'(a,i8)')  'Out of memory in Occam.  Too many associated parameters: ', npm_assoc
        call printOccamLog(cStr)          
        stop 
    endif           

 
    
    end subroutine allocateOccamIteration
    
!==================================================================================================================================! 
!========================================================================================================= deallocateOccamIteration
!==================================================================================================================================! 
    subroutine deallocateOccamIteration
!
! Deallocates arrays used during an Occam iteration
!    
    if (allocated(wjtwj))       deallocate(wjtwj)
    if (allocated(wjtwd))       deallocate(wjtwd)
    if (allocated(resid))       deallocate(resid)
    if (allocated(amat))        deallocate(amat)
    if (allocated(pm_b))        deallocate(pm_b)
    if (allocated(dm_b))        deallocate(dm_b)
    if (allocated(pm_test))     deallocate(pm_test)    
    if (allocated(pm_assoc_b))  deallocate(pm_assoc_b) 
    if (allocated(pm_c))        deallocate(pm_c)
    if (allocated(dm_c))        deallocate(dm_c)
    if (allocated(pm_assoc_c))  deallocate(pm_assoc_c) 
 
        
    end subroutine deallocateOccamIteration    
    
!==================================================================================================================================! 
!========================================================================================================= deallocateOccamIteration
!==================================================================================================================================! 
    subroutine deallocateOccam
!
! Deallocates arrays used during the entire Occam inversion
! 
    if (allocated(pm))          deallocate(pm)
    if (allocated(prewts))      deallocate(prewts)
    if (allocated(premod))      deallocate(premod)
    if (allocated(preDiff))     deallocate(preDiff)
    if (allocated(preDiffwts))  deallocate(preDiffwts)    
    if (allocated(ijDiff))      deallocate(ijDiff)       
    if (allocated(d))           deallocate(d)
    if (allocated(sd))          deallocate(sd)
    if (allocated(dp))          deallocate(dp)
    if (allocated(dm))          deallocate(dm)    
    if (allocated(pm_assoc))    deallocate(pm_assoc)    
    if (allocated(penaltywts))  deallocate(penaltywts) 
    if (allocated(ijpenalty))   deallocate(ijpenalty)   
    if (allocated(pen%colind))  deallocate(pen%colind)   
    if (allocated(pen%val))     deallocate(pen%val)   
    if (allocated(pen%rowptr))  deallocate(pen%rowptr)   
    if (allocated(lowerBound))  deallocate(lowerBound)
    if (allocated(upperBound))  deallocate(upperBound)
    
    end subroutine deallocateOccam
  
!==================================================================================================================================! 
!============================================================================================================ computeOccamIteration
!==================================================================================================================================! 
    subroutine computeOccamIteration 
!
! Computes an Occam iteration. 
!    

!
! Local variables:
! 
    logical         :: lGetIntercept   ! true if target RMS obtained during this iteration
    character(128)  :: cStr   
    real(RealPrec)  :: targetRMS_input 
 
!
! Say hello and initialize a few quantities
!        
    call printOccamIntroMessage  
     
!
! Compute forward response and Jacobian matrix, also allocate arrays
! used in the model update equations during Occam's lagrange multiplier search:
!
    call computeJacobian

! 
! Save the starting model response if this is the first call:
! 
    if (currentIteration == 1) call writeResponse(0)       
    
    
!
! If using Slow Occam, reset the target misfit:
!
    if (lConvergeSlowly) then
        
        targetRMS_input = targetRMS
        targetRMS       = max(startingRMS*rmsThreshold,targetRMS)
        
        if (targetRMS /= targetRMS_input) then  
            write(cStr,'(a,g16.4)')  ' Occam set to smoothly converge, for this iteration the target RMS is: ',targetRMS
            call printOccamLog(cStr) 
        endif
        
    endif



!
! Is the starting model at the target RMS already? If so, let the smoothing begin.
!    
    if ( startingRMS  <  targetRMS + rmsTol ) lTargetHit = .true.  
        
!
! Allocate some working arrays used throughout the rest of the iteration
! This is done after the Jacobian computation in order to keep the memory usage overlap as low as possible
!
    call allocateOccamIteration                

!
! Get minimum RMS as a function of mu:
!
    call getMinimum(modelMu)       ! optimal model variables are returned in _b arrays
           
    if (convergenceFlag == 4) then ! too many step size cuts, give up         
        call deallocateOccamIteration
        return
    endif  
    

!
! Minimum has been obtained, now analyze it:
!

    ! Is the minimum less than the target RMS? If so, get the intercept with the target RMS       
    if ( rms_b  <= targetRMS + rmsTol  )  then 
        
        convergenceFlag =  1
        
        lGetIntercept = .true.
        
        do while (lGetIntercept)     
    
            call getIntercept ! On return the _b arrays have the intercept model and responses
            
            if ( lTargetHit )  then ! Occam is in the smoothing phase, see if the new optimal model is smoother:
                
                if (rough_b == 0 ) then ! perfectly smooth model found
                    
                    convergenceFlag =  3   
                    exit  ! the do loop  
                    
                elseif ( rough_b < startingRoughness )  then ! smoother model found
                
                    convergenceFlag =  1                 
                    exit  ! the do loop  
                     
                else ! ( rough_b >= startingRoughness ) ! Roughness has increased, try cutting the stepsize:
                        
                    stepCut     = stepCut/2.
                    numStepCuts = numStepCuts + 1
                    call initializeTrackers
                    
                    if (numStepCuts < maxNumStepCuts ) then 
                        
                        ! Otherwise plow on
                        write(cStr,'(a,g16.4)')  ' Cutting step size due to model roughness, new step size is: ',stepCut
                        call printOccamLog(cStr)
                        
                        ! Get a new minimum             
                        call getMinimum(mu_b) 
                        
                         ! kwk debug: what if this throws a convergenceFlag = 4 ?
                         
                        cycle ! move on to the next lGetIntercept while loop iteration  
                        
                    else ! We've already cut the step size several times to no avail, give up:        
                        
                        convergenceFlag = 5
                
                        cStr =  ' Roughness not decreased by using small steps, giving up.'
                        call printOccamLog(cStr)
                        
                        exit ! the do loop       
                               
                    endif    
                                
                endif
                
            else ! This is the first time the target RMS has been obtained, just return a happy camper:
                
                convergenceFlag = 1       
                exit ! the do loop    
            
            endif
            
        enddo ! while (lGetIntercept)     
    
    endif   
 
!
! Compute step size (how much did the model change): 
!
    pm       = pm_b - pm
    stepSize = sqrt(dot_product(pm,pm)/nParams)  
 
    if ( (convergenceFlag == 1).and.(lTargetHit) ) then ! Check model step size
        if ( stepSize <  stepSizeTol ) convergenceFlag = 2
    endif

!
! Move optimal model to output:
!
    pm              = pm_b
    dm              = dm_b
    modelRoughness  = rough_b
    modelMu         = mu_b    
    modelRMS        = rms_b
    if (npm_assoc > 0) pm_assoc = pm_assoc_b
    
!
! Are we at the last allowable iteration?
!    
    if (currentIteration >= maxOccamIterations)  convergenceFlag = 6

!
! Shall we apply David Myer's auto-convergence operation:
!           
!    if ( lAutoConverge .and. ((convergenceFlag == 4 ).or.(convergenceFlag == 5)) ) then
!        
!        targetRMS       = modelRMS*nAutoConvergeFactor
!        convergenceFlag = 0
!        
!        write(cStr,'(a,g16.4)')  ' Auto-convergence trigged, backing off to new target misfit:  ' , targetRMS
!        call printOccamLog(cStr)
!                        
!        if( currentIteration > maxOccamIterations - 10 ) maxOccamIterations = currentIteration + 10
!        
!    endif   
    
    if (lConvergeSlowly) then
        targetRMS = targetRMS_input 
        if ( modelRMS > targetRMS + rmsTol ) convergenceFlag = 0
    endif

    
!
! Print out some summary information:
! 
    call printOccamExitMessage

!
! Deallocate working arrays:
!        
    call deallocateOccamIteration    
    

    end subroutine computeOccamIteration

!==================================================================================================================================! 
!================================================================================================================== computeJacobian
!==================================================================================================================================! 
    subroutine computeJacobian 
!
! Computes the forward response and Jacobian matrix for the input model at the start of an Occam iteration.
! Inserts these into various arrays required by the Occam model update equations.
!
      
!
! Local variables:
!
    real(RealPrec), allocatable :: dhat(:)
    integer                     :: i, ierr
    real(RealPrec)              :: beta, t0, t1, t2, minv, maxv
    character(128)              :: cStr

           
!
! Allocate arrays:
!      
    allocate( wj(nd,nParams), dhat(nd), resid(nd), stat=ierr )
    if (ierr .ne. 0) then
        write(cStr,'(a,1x,i8,1x,i8)') 'Out of memory in Occam.  Too many free parameters and data: ', nParams, nd
        call printOccamLog(cStr)
        stop 
    endif
    wj = 0
 
!
! Compute the forward response and the Jacobian matrix: results returned in arrays dm and wj
!     
    call get_time_offset(0.d0,t0)
    
    call computeFwd( .true., transformToBound(pm,lowerBound,upperBound,lBoundMe)) ! don't forget to convert back to bound params:  
    
    call get_time_offset(t0,t1)
    
    numForwardCalls = numForwardCalls + 1     
    
    call get_time_offset(0d0,t0)
    
!
! Scale Jacobian for nonlinear transformation to bound model parameters:
!
    call transformWJ()
 
!
! Compute residual and RMS misfit:
!  
    resid = (d - dm) / sd   
      
    startingRMS       = sqrt(dot_product(resid,resid)/nd)
    startingRoughness = getRoughness(pm)
    
!           
! Weight the Jacobian matrix by the data uncertainty, also apply optional data weights:
!
    do i = 1,nd
        wj(i,:) = (wj(i,:)* d_wt(i) / sd(i) )   ! Array math - possible alloc problem
    enddo
    
!
! Form arrays needed for the Occam model update equations:
! 
    
    cStr = ' ... Forming matrix products... '
    call printOccamLog(cStr)
    
    allocate(wjtwj(nParams,nParams) , wjtwd(nParams), stat = ierr )   
    if (ierr .ne. 0) then
        write(cStr,'(a,1x,i8)') 'Out of memory in Occam.  Too many free parameters: ', nParams 
        call printOccamLog(cStr)
        stop 
    endif    
    wjtwj = 0
    
    !
    ! Form (WJ)^T WJ  
    !
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
    if (nParallelProc < MinProcScaLAPACK) then
#endif
        call multiplyATA(wj,wjtwj)  
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
    else
        call multiplyATA_scalapack 
    endif
#endif
    
    ! Form (WJ)^T W \hat(d)
    beta = 1.0
    dhat = resid*d_wt 
    call multiplyAx(wj,pm,dhat,'n',beta)     ! dhat = wj*pm + beta*dhat
    beta = 0.0
    call multiplyAx(wj,dhat,wjtwd,'t',beta)  ! wjtwd = wj*dhat 
    deallocate( wj, dhat )
    
!
! Display the input model information:
!    
    call get_time_offset(t0,t2)
    minv              = minval(transformToBound(pm,lowerBound,upperBound,lBoundMe))
    maxv              = maxval(transformToBound(pm,lowerBound,upperBound,lBoundMe))

        
    cStr = ' ... Starting model and matrix assembly:'
    call printOccamLog(cStr)
    write(cStr,'(7(a16,1x))') '     Misfit     ','   Roughness    ','   Log10(mu)    ', &
                            & ' Min log10(rho) ',' Max log10(rho) ', '  Jacobian (s)  ', ' Matrix Ops. (s)'
    call printOccamLog(cStr)
    write(cStr,'(7(g16.5,1x))')  startingRMS , startingRoughness, modelMu,  minv, maxv, t1, t2 
    call printOccamLog(cStr)
                
    end subroutine computeJacobian   

!==================================================================================================================================! 
!======================================================================================================================= misfitOfMu
!==================================================================================================================================! 
    real(RealPrec) function misfitOfMu(alogmu)   
!
! Solves the Occam model update equations for a given log10(mu)
!    
    
!
! Input arguments:
!  
    real(RealPrec), intent(in)  :: alogmu
    
!
! Local variables:
!    
    real(RealPrec)  :: minv, maxv, mu, t0, t1, t2, w_i, w_j
    integer         :: istat,i,j,ipen, irow, i1,i2, istart, iend
    real(RealPrec)  :: w2, chi2, w_mgs
    character(128)  :: cStr

    call get_time_offset(0.d0,t0)
     
! 
! Check for very large or small mu, if so set large misfit and return:
!
    if (abs(alogmu) > 20.) then
    
        cStr = ' Skipping forward call due to extremely large or small mu in misfitOfMu ...'
        call printOccamLog(cStr)
        misfitOfMu = 1000. + abs(alogmu)  ! return large misfit that grows with small or large mu 
        ! don't print misfit to logfile since this number is artificial
        return  
        
    endif
    
    mu = 10d0**(alogmu)
    
!
! Construct rhs vector:
!            
    pm_test  = wjtwd + mu*(prewts**2)*premod     
    
    ! Add on parameter difference preferences, if any given:
    do irow = 1,nDiff
        i  = ijDiff(irow,1) ! +1
        j  = ijDiff(irow,2) ! -1
        w2 = preDiffwts(irow)**2
         
        pm_test(i) = pm_test(i) + mu*w2*preDiff(irow)
        pm_test(j) = pm_test(j) - mu*w2*preDiff(irow)
    enddo
    
    
!
! Construct matrix A:
!
    amat  = 0.0
    
  ! First insert penalty terms for del^T del:   
  
    if (allocated(pen%val)) then
  
        do irow = 1,pen%nrows 
            
 
            istart       = pen%rowptr(irow)
            iend         = pen%rowptr(irow+1)-1
            
            ! Re-weight if using minimum gradient support regularization:
            if (lMGS) then
                w_mgs = 1d0/( dot_product( pen%val(istart:iend), pm(pen%colind(istart:iend)) )**2 + beta_mgs**2  )  
            else
                w_mgs = 1d0
            endif
                        
            do i1 = istart,iend
                do i2 = istart,iend
                    
                    i = pen%colind(i1)
                    j = pen%colind(i2)
                    
                    w_i =  pen%val(i1)
                    w_j =  pen%val(i2) 
                    
                    amat(i,j) = amat(i,j) + w_i*w_j*w_mgs 
 
                    
                enddo
            enddo
                 
 
        enddo  
  
    else ! my previous simple sparse matrix of differences:  
  
        do ipen = 1,npenalty
            i = ijpenalty(ipen,1)
            j = ijpenalty(ipen,2)
            w2 = penaltywts(ipen)**2
            
            ! kwk mgs test:
            !w2 =  w2 / ( (pm(i) - pm(j) )**2 + 0.01**2 )
            
            amat(i,j) = amat(i,j)  - w2
            amat(j,i) = amat(j,i)  - w2
            amat(i,i) = amat(i,i)  + w2
            amat(j,j) = amat(j,j)  + w2
        enddo
        
    endif
         
    ! Add on the penalty terms for p^T p, where p is diagonal:
    forall (i = 1:nParams, prewts(i) .ne. 0)   ! prewts is usually entirely zero
        amat(i,i) = amat(i,i) + prewts(i)**2
    end forall                
 
    ! Add on parameter difference preference terms, if any given:
    do irow = 1,nDiff
        i  = ijDiff(irow,1) ! +1
        j  = ijDiff(irow,2) ! -1
        w2 = preDiffwts(irow)**2
        
        amat(i,j) = amat(i,j)  - w2
        amat(j,i) = amat(j,i)  - w2
        amat(i,i) = amat(i,i)  + w2
        amat(j,j) = amat(j,j)  + w2
    enddo    
 
    ! Add on the WJ^T WJ   
    amat = mu*amat + wjtwj
   
 
!
! Solve the linear system: 
!
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
    if (nParallelProc < MinProcScaLAPACK) then
#endif
        call solveLinearSystem(nParams, amat, pm_test ,istat)
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
    else
        call solveLinearSystem_scalapack(istat) 
    endif
#endif


    if (istat > 0) then
        cStr = ' Lapack Cholesky factorization failed for this mu, returning large misfit...'
        call printOccamLog(cStr)
        misfitOfMu = 1000. + abs(log10(mu))  ! return large misfit that grows with small or large mu 
        return  
    endif    
 
!    
! Cut step size if necessary
! 
    if (stepCut < 1.0)  pm_test = (1.-stepCut)*pm + stepCut*pm_test  

    call get_time_offset(t0,t1)
    
!
! Compute forward response:
!    
    call get_time_offset(0d0,t0)  
      
    call computefwd( .false.,transformToBound(pm_test,lowerBound,upperBound,lBoundMe)) ! dont forget to convert back to bounded
    
    call get_time_offset(t0,t2)
        
    numForwardCalls = numForwardCalls + 1  

        
!
! Compute RMS misfit:
!
    resid       = (d - dm)/sd 
    chi2        = dot_product(resid,resid)
    misfitOfMu  = sqrt(chi2 / dble(nd))
        
    ! If the misfit is nan then exit early.  Inversion has failed.
    if (misfitOfMu >= huge(misfitOfMu)) then
        cStr = 'Error: nan misfit.  Process aborted.'
        call printOccamLog(cStr)
        stop
    endif
                
!
! Compute test model roughness:
!
    rough_test = getRoughness(pm_test)

!
! Update the (mu,rms) tracker 
!    
    call updateTrackers(alogmu,misfitOfMu)
    
!
! Report to terminal and logfile:
!         
    minv        = minval(transformToBound(pm_test,lowerBound,upperBound,lBoundMe))
    maxv        = maxval(transformToBound(pm_test,lowerBound,upperBound,lBoundMe))
               
    write(cStr,'(7(g16.5,1x))')  misfitOfMu, rough_test, alogmu,  minv, maxv, t2, t1  
    call printOccamLog(cStr)
            

    end function misfitOfMu   

!==================================================================================================================================! 
!================================================================================================================ initializeTrackers
!==================================================================================================================================! 
    subroutine initializeTrackers
!
! The trackers keep track of the best fitting mu and the mu immediately to the right, so that if Occam meets the target RMS,
! it will already have bounding mu's for the right side intercept.
!    
!  The two trackers need to be reinitialized from a few places, so they are placed in this subroutine to hide the distracting
! details from the rest of the code.
!
    tracker1_mu  = -tracker_Dummy  ! left and right bounds on mu
    mu_c         =  tracker_Dummy
    tracker1_rms =  tracker_Dummy  ! large positive misfits for dummy values
    rms_c        =  tracker_Dummy    
    lTrackerOn   = .true.
    
    end subroutine initializeTrackers

!==================================================================================================================================! 
!=================================================================================================================== updateTrackers
!==================================================================================================================================! 
    subroutine updateTrackers(alogmu,rms)
    
    real(RealPrec), intent(in)  :: alogmu,rms
    real(RealPrec)              :: mOld

    if (.not.lTrackerOn) return
    
    if (rms < tracker1_rms) then ! tracker (1) has lowest RMS 
       
        mOld = tracker1_mu
        
        tracker1_mu  = alogmu  ! model will go to pm_b, dm_b so no need to save it here
        tracker1_rms = rms
        
        if ( mOld > tracker1_mu  ) then
        
            mu_c    = mu_b     ! here we need to save the model parameters so they don't conflict with the _a, _b and _c arrays
            rms_c   = rms_b     
            rough_c = rough_b      
            dm_c    = dm_b
            pm_c    = pm_b
           if (npm_assoc > 0) pm_assoc_c = pm_assoc_b
 
        endif 
        
    elseif  ( (alogmu > tracker1_mu) .and. (alogmu < mu_c ) ) then
        mu_c    = alogmu
        rms_c   = rms    
        rough_c = rough_test
        dm_c    = dm       
        pm_c    = pm_test
        if (npm_assoc > 0) pm_assoc_c = pm_assoc 
  
    endif      
  
 
    
    end subroutine updateTrackers       
!==================================================================================================================================! 
!===================================================================================================================== getRoughness
!==================================================================================================================================! 
    logical function lEscapeFromMinimization(rms)  
!
! Test various conditions for which the Occam minimization should be terminated early.
!    
    real(RealPrec), intent(in) :: rms
    
    lEscapeFromMinimization = .false.
    
    !
    ! Condition 1: the rms is less than or equal to the target
    !
    if (rms <= targetRMS + rmsTol) then
        
        lEscapeFromMinimization = .true.
     
    !    
    ! Condition 2: a significant decrease in RMS has occurred
    !
    elseif    ( rms <= rmsThreshold*startingRMS) then
        
         lEscapeFromMinimization = .true.
    
    endif   
 
    end function lEscapeFromMinimization 
    
!==================================================================================================================================! 
!===================================================================================================================== getRoughness
!==================================================================================================================================! 
    real(RealPrec) function getRoughness(model)  
!
! Computes the roughness for the input model
!    
    real(RealPrec), dimension(nParams), intent(in) :: model
    
    integer         :: ipen, i, j, istart, iend, irow
    real(RealPrec)  :: rtemp, w_mgs
    
    getRoughness = 0. 
      
    !
    ! KWK: June 1 2012 modification to now use CSR format for input penalty matrix (more flexible than previous):
    !  
    if (allocated(pen%val)) then
    
        do i = 1,pen%nrows 
        
                istart       = pen%rowptr(i)
                iend         = pen%rowptr(i+1)-1
                
                rtemp        = dot_product( pen%val(istart:iend), model(pen%colind(istart:iend)) )
                
                ! Re-weight if using minimum gradient support regularization:
                if (lMGS) then
                    w_mgs = 1d0/sqrt(  dot_product( pen%val(istart:iend), pm(pen%colind(istart:iend)) )**2 + beta_mgs**2   )
                    rtemp = rtemp*w_mgs  
                endif
                
                getRoughness = getRoughness + rtemp**2
                
        enddo  
        
        
    else ! my previous simple sparse matrix of differences:  
        do ipen = 1,npenalty
           i = ijpenalty(ipen,1)
           j = ijpenalty(ipen,2)            
           
           ! kwk mgs test:
           !getRoughness=     getRoughness + (penaltywts(ipen)/sqrt( (pm(i) - pm(j) )**2 + 0.01**2 )*(model(i) - model(j)))**2 
     
           getRoughness = getRoughness + (penaltywts(ipen)*(model(i) - model(j)))**2 
        enddo
    endif
         
  ! Now add on the norm for the preference model (which is 0 if not used since prewts=0 then):
    getRoughness = getRoughness + sum( (prewts*(premod - model))**2 )
 
  ! Add on norm for model parameter difference preferences, if any given:
    do irow = 1,nDiff
        i  = ijDiff(irow,1) ! +1
        j  = ijDiff(irow,2) ! -1
        getRoughness = getRoughness + (preDiffwts(irow)*(model(i) - model(j)))**2 
    enddo    
    
 
    end function getRoughness
    
!==================================================================================================================================! 
!======================================================================================================================= getMinimum
!==================================================================================================================================! 
    subroutine getMinimum(lastMu) 
!    
! Sweeps through mu to find optimal RMS
!

!
!   Input arguments:
!
    real(RealPrec), intent(inout)  :: lastMu
 
!
! Local variables:
!
    real(RealPrec)              :: mu1,mu2,mu3, rms2
    logical                     :: lGetMinimum    
   
       
    character(128) :: cStr 
    
    cStr = ' Searching Lagrange multiplier...'  
    call printOccamLog(cStr)

     
    lGetMinimum = .true.    
 
    do while (lGetMinimum)
                        
!
! Bracket the minimum using parabolic interpolation or golden section search:
!
        mu1 = lastMu- 1.0
        mu2 = lastMu 
        
        call bracketMinimum(mu1,mu2,mu3,misfitOfMu)  ! best model is returned in _b arrays
        
!    
! Now find the minimum if the RMS is still big:
!
        if  (rms_b  <= targetRMS + rmsTol)  then ! RMS is good enough already 
            ! nothing to do...
        
        elseif ( rms_b <= rmsThreshold*startingRMS )  then    
            ! nothing to do...
            
        elseif ( rms_b > targetRMS + rmsTol ) then
        
            ! Find the minimum using Brent's method:
            rms2 = rms_b
            call findMinimum(mu1,mu2,mu3,rms2,misfitOfMu)            
        
        endif

!
! Analyze minimum or whatever rms is returned if shortcuts were taken:
!
        if  ( (rms_b >= startingRMS .and.(.not.lTargetHit) ) .or. ( (lTargetHit) .and. (rms_b > targetRMS + rmsTol ) ))then  
           ! Minimum isn't as good as starting model, cut step size:
            
            stepCut         = stepCut/2.0
            numStepCuts     = numStepCuts + 1
            lastMu          = mu_b  ! save last minimum mu for next call to bracketMinimum
            call initializeTrackers
            
            ! Have we cut the step size too many times to no avail?
            if (numStepCuts > maxNumStepCuts) then 
                 
                convergenceFlag = 4
                lGetMinimum     = .false.
           
            else 
                
                write(cStr,'(a,g16.4)')  ' Cutting the size of the model update due to divergence, new fractional size: ',stepCut
                call printOccamLog(cStr)    
           
            endif
            
        elseif (rms_b <= targetRMS + rmsTol ) then 
        
            lGetMinimum = .false.
            
        elseif ( rms_b <= rmsThreshold*startingRMS )  then

            write(cStr,'(a,f5.2)')  ' Large misfit decrease detected, ending minimization search. Percent decrease in misfit: ', &
                                    & ( startingRMS-rms_b )/startingRMS*100
            call printOccamLog(cStr)
            
            lGetMinimum = .false.
            
        else    ! no special case, the minimum rms is less than or equal to the starting model so exit the while loop   
        
            lGetMinimum = .false.
   
        endif
        
    enddo ! do while (lGetMinimum) 
    
    end subroutine getMinimum 
    
!==================================================================================================================================! 
!===================================================================================================================== getIntercept
!==================================================================================================================================! 
    subroutine getIntercept 
!
! Given input mu_b with rms_b < targetRMS, this routine finds a larger mu above the intercept and then calls
! the NR root finding routine to find the intercept to within a specified tolerance.  
!

!    
! Local variables:
!    
    real(RealPrec)  :: interceptMu, fac,  mu_b_start
    character(180)  :: cStr
    logical         :: lSkipIntercept

!
! Check to see if pm_c or pm_b are close enough to the targetRMS:
!
    if ( abs(rms_c - targetRMS) < rmsTol) then ! _c is from larger mu, so give it precedence
    
        rms_b   = rms_c
        mu_b    = mu_c
        pm_b    = pm_c   
        dm_b    = dm_c
        rough_b = rough_c
        if (npm_assoc > 0) pm_assoc_b = pm_assoc_c       
        
        return  
           
    elseif   ( abs(rms_b - targetRMS) < rmsTol)  then 
    
        return
        
    endif
    
!
! Allocate arrays for carrying around models associated with a (b uses _b arrays from module)
!
    allocate( pm_a(nParams), dm_a(nd) )
    if (npm_assoc > 0) allocate( pm_assoc_a(npm_assoc) )

    cStr = ' ... Finding intercept:'
    call printOccamLog(cStr)
    write(cStr,'(7(a16,1x))') '     Misfit     ','   Roughness    ','   Log10(mu)    ', &
                            & ' Min log10(rho) ',' Max log10(rho) ', '  Fwd Call (s) ', ' Matrix Ops. (s)'
    call printOccamLog(cStr)      
    
    lTrackerOn     = .false.
    lSkipIntercept = .false.    
    
!
! If rms_c has been set, use it to bound intercept, otherwise find a mu with an rms above the intercept
!   
 
    if ( (rms_c < tracker_Dummy) .and. (rms_c > targetRMS) ) then
        
        rms_a   = rms_b
        mu_a    = mu_b
        pm_a    = pm_b   
        dm_a    = dm_b
        rough_a = rough_b
        if (npm_assoc > 0) pm_assoc_a = pm_assoc_b
           
        rms_b   = rms_c
        mu_b    = mu_c
        pm_b    = pm_c   
        dm_b    = dm_c
        rough_b = rough_c
        if (npm_assoc > 0) pm_assoc_b = pm_assoc_c    
            
    else  ! find mu_c with rms_c above the intercept                                                                            
        
        fac =  0.0     
        mu_b_start = mu_b
                     
        do while (rms_b < targetRMS)   ! Keep on increasing mu until we cross the targetRMS
            
            ! Update arrays keep track of models near intercept (_b goes to _a before updating _b):
            rms_a   = rms_b
            mu_a    = mu_b
            pm_a    = pm_b   
            dm_a    = dm_b
            rough_a = rough_b
            if (npm_assoc > 0) pm_assoc_a = pm_assoc_b
            
            ! Test new mu_b:
            fac     = fac+1.0
            mu_b    = mu_b + fac*0.699d0  ! Move up by 1/2 decade, then full decade .... 
            
            if (mu_b  < 12.0 ) then 
               
                rms_b   = misfitOfMu(mu_b)
                
                ! Save results in _b arrays:
                call updateB(mu_b,rms_b)      
                
            else
            
                lSkipIntercept = .true.
                cStr = ' Intercept not found by increasing mu, giving up...'
                call printOccamLog(cStr)  
                mu_b = mu_b_start  ! reset mu to reasonable value
                
                exit ! while loop since we've gone far enough
                
            endif 
            
        enddo
    
    endif
        
!
! mu_a and mu_b bound the intercept, now find the intercept more precisely:
!
    if (.not.lSkipIntercept) then

        call findIntercept( misfitOfMu,interceptMu)  

        ! Returning from findIntercept(), the intercept model is in pm_b
        write(cStr,'(a,g16.4)')  ' Intercept is at mu: ',interceptMu
        call printOccamLog(cStr)    
              
    endif    
!
! There's nothing more to do except deallocate: 
!
    deallocate( pm_a, dm_a )
    if (npm_assoc > 0) deallocate( pm_assoc_a )
       

                
    end subroutine getIntercept

!==================================================================================================================================! 
!=========================================================================================================== Linear Algebra Routines
!==================================================================================================================================! 
!
! Various helper subroutines used for the linear algebra aspects of Occam.
!
! All of these call LAPACK and BLAS routines, and are set to use the multithreaded Intel Math Kernel Library if available.
!
! multiplyATA:          Computes A^t A, uses calls to BLAS routines for speed (faster when nparams is large)
! multiplyAx:           Computes A*x or A^t*x
! solveLinearSystem:    Solves system Ax=b using LAPACK's Cholesky factorization
!     
!-----------------------------------------------------------------------------------------------------------------------------------    
    subroutine multiplyAx(A,x,c,trans,beta)
!-----------------------------------------------------------------------------------------------------------------------------------  

    real(RealPrec), intent(in),  dimension(:,:)  :: A
    real(RealPrec), intent(in),  dimension(:)    :: x
    real(RealPrec), intent(in)                   :: beta
    character(1),   intent(in)                   :: trans 
    real(RealPrec), intent(out), dimension(:)    :: c
    
    
    
    integer        :: inc, m, n
    real(RealPrec) :: alpha
    
 
    alpha  = 1.
    inc    = 1
    m      = size(A,dim=1);
    n      = size(A,dim=2);
    
#if defined(__INTEL_COMPILER)  
    call  mkl_set_num_threads(nMKL_numThreads)
#endif

    if (RealPrec == kind(0d0) ) then  ! Double
        call dgemv(trans,m,n,alpha,A,m,x,inc,beta,c,inc)
    else ! (RealPrec == kind(0e0) ) Single
        call sgemv(trans,m,n,alpha,A,m,x,inc,beta,c,inc)
    endif

#if defined(__INTEL_COMPILER)  
    call mkl_set_num_threads ( 1 ) 
#endif
     
     
    end subroutine multiplyAx

!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine multiplyATA(A,C)
!----------------------------------------------------------------------------------------------------------------------------------- 
!
! Arguments:
!
    real(RealPrec), intent(in),  dimension(:,:)  :: A
    real(RealPrec), intent(out), dimension(:,:)  :: C
   
!
! Local variables:
!
    integer        :: i,m,n
    real(RealPrec) :: alpha, beta
 
    c      = 0.
    alpha  = 1.
    beta   = 0.
    m      = size(A,dim=1);
    n      = size(A,dim=2);
 
    
#if defined(__INTEL_COMPILER)  
    call  mkl_set_num_threads(nMKL_numThreads)
#endif
    
    if (RealPrec == kind(0d0) ) then  ! Double
       call dsyrk('u','t',n,m,alpha,a,m,beta,c,n)
    else ! (RealPrec == kind(0e0) ) Single
       call ssyrk('u','t',n,m,alpha,a,m,beta,c,n)
    endif

#if defined(__INTEL_COMPILER)  
    call mkl_set_num_threads ( 1 ) 
#endif
     
    do i=2,n
        C(i,1:i-1) = C(1:i-1,i)
    enddo  
 
           
    end subroutine multiplyATA
    
!-----------------------------------------------------------------------------------------------------------------------------------    
    subroutine solveLinearSystem(m,A,b,istat)
!-----------------------------------------------------------------------------------------------------------------------------------    
!
! Solves Ax=b using LAPACK's Cholesky factorization
!
! The solution x is written to the right hand side vector b on output.
!
!
! Input arguments:
!
    integer, intent(in)           :: m
    integer, intent(out)          :: istat
    real(RealPrec), intent(in)    :: A(m,m)
    real(RealPrec), intent(inout) :: b(m) 

#if defined(__INTEL_COMPILER)  
    call  mkl_set_num_threads(nMKL_numThreads)
#endif

    if (RealPrec == kind(0d0) ) then  ! Double
    
        call dpotrf( 'u', m, A, m, istat )
        if (istat > 0) return
        
        call dpotrs( 'u', m, 1, A, m, b, m, istat )
        if (istat > 0) return
              
    else    ! Single
    
        call spotrf( 'u', m, A, m, istat )
        if (istat > 0) return  
       
        call spotrs( 'u', m, 1, A, m, b, m, istat )     
        if (istat > 0) return  
        
    endif 
   
#if defined(__INTEL_COMPILER)  
    call  mkl_set_num_threads(1)
#endif
    
       
    end subroutine solveLinearSystem
    

#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mpi_worker_ScaLapackMM()
!-----------------------------------------------------------------------------------------------------------------------------------
    
!
! Matrix Matrix Multiply
!       
 
    
    integer :: i, istat, numroc, info, ictxt
    integer :: myid, ierr, nprocs 
    
    integer   :: mycol, myrow, nb
    integer   :: npcol, nprow
    integer,parameter :: descriptor_len=9
    
    
    integer   :: l_nrows_wj,l_ncols_wj    
    integer   :: l_nrows_wjtwj,l_ncols_wjtwj
    
    integer   :: desc_wj( descriptor_len ),desc_wjg( descriptor_len )
    integer   :: desc_wjtwj( descriptor_len )
    integer   :: desc_wjtwjg( descriptor_len )

    integer   :: rootNodeContext, sizeM, sizeN
    
    real(RealPrec), dimension(:,:),allocatable :: wj_loc,wjtwj_loc
 
 
!
! Which processor am I and how many others are there?
!
    call blacs_pinfo( myid, nprocs )
        
!
! Broadcast matrix size to all processors:
! 
    if (myid == 0)  sizeM = size(wj,dim=1)
    if (myid == 0)  sizeN = size(wj,dim=2)
    call mpi_bcast(sizeM, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(sizeN, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

!    
! Set the dimension of the 2d processors grid:
!
    call gridsetup(nprocs,nprow,npcol)
    ! The small processor grid, short, fat matrix tweak:
    ! This is for SOME cpu's (such as Triton) where pdsyrk stalls for fat matrices and small grids
    if ( (nprocs <= 16).and.(sizeN/sizeM > 2) ) then
        nprow = 1
        npcol = nprocs
    endif

!
! Initialize a single blacs context:
!
    call blacs_get(-1,0,ictxt) 
    call blacs_gridinit(ictxt,'r',nprow,npcol) 
    call blacs_gridinfo(ictxt, nprow,npcol,myrow,mycol) 

!   
! Calculate the blocking factor for the matrix:
!
    call blockset( nb, 64, sizeM, sizeN, nprow, npcol) 
    
!      
! Distributed matrices: get num. local rows/cols. Create description:
!
    l_nrows_wj = numroc(sizeM,nb,myrow,0,nprow)
    l_ncols_wj = numroc(sizeN,nb,mycol,0,npcol)
    call descinit( desc_wj, sizeM, sizeN, nb, nb, 0, 0, ictxt, l_nrows_wj, info )  

    l_nrows_wjtwj = numroc(sizeN,nb,myrow,0,nprow)
    l_ncols_wjtwj = numroc(sizeN,nb,mycol,0,npcol)
    call descinit( desc_wjtwj, sizeN, sizeN, nb, nb, 0, 0, ictxt, l_nrows_wjtwj, info ) 

    allocate (wj_loc(l_nrows_wj, l_ncols_wj ), stat=istat)
    if (istat /= 0) then
        write(*,*) 'Error allocating wj_loc! l_nrows_wj,l_ncols_wj =',l_nrows_wj,l_ncols_wj
        stop
    endif
   
    allocate (wjtwj_loc(l_nrows_wjtwj, l_ncols_wjtwj ), stat=istat) 
    if (istat /= 0) then
        write(*,*) 'Error allocating wjtwj_loc! l_nrows_wjtwj,l_ncols_wjtwj =',l_nrows_wjtwj, l_ncols_wjtwj
        stop
    endif
    
! 
! Distribute sub-blocks of global wj from master to all process:
!
    call blacs_get( -1, 0, rootNodeContext )
    call blacs_gridinit( rootNodeContext, 'row-major', 1, 1 )

     if (myid == 0) then
        call descinit( desc_wjg, sizeM, sizeN, sizeM, sizeN, 0, 0, rootNodeContext, sizeM, info )
        call descinit( desc_wjtwjg, sizeN, sizeN, sizeN, sizeN, 0, 0, rootNodeContext, sizeN, info )
    else 
        desc_wjg        = 0
        desc_wjg(2)     = -1
        desc_wjtwjg     = 0
        desc_wjtwjg(2)  = -1
    endif
 
    call pdgemr2d( sizeM, sizeN, wj, 1, 1, desc_wjg, wj_loc, 1, 1, desc_wj, desc_wj( 2 ) )
    
!
! Parallel matrix multiply:
!           
    call pdsyrk('u','t',sizeN,sizeM,1d0,wj_loc,1,1,desc_wj,0d0,wjtwj_loc,1,1,desc_wjtwj) 
        
!
! Send solution in wjtwj_loc back to matrix wjtwj:
!
     
    call pdgemr2d( sizeN, sizeN, wjtwj_loc, 1, 1, desc_wjtwj, wjtwj, 1, 1, desc_wjtwjg, ictxt )
    if (myid == 0) then  
        ! fill in lower triangle
        do i=2,sizeN
            wjtwj(i,1:i-1) = wjtwj(1:i-1,i)
        enddo    
    endif
    if (myid == 0)  call blacs_gridexit( rootNodeContext ) 

!
! Deallocate local arrays...
!
    if (allocated(wj_loc))    deallocate(wj_loc)
    if (allocated(wjtwj_loc)) deallocate(wjtwj_loc)
        
!
! Shut down blacs:
!    
    call blacs_gridexit(ictxt) 

    end subroutine mpi_worker_ScaLapackMM

 
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mpi_worker_ScaLapackMS(istat)
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Matrix Factorization and Solve
!    
! This is called from solveLinearSystem_scalapack() in mare2dem_mpi.f90
!

  
    integer :: numroc, info,ictxt
    integer :: myid, ierr, nprocs, istat
    
    integer   :: mycol, myrow, nb
    integer   :: npcol, nprow
    integer,parameter :: descriptor_len=9

    integer   :: l_nrows_amat,l_ncols_amat    
    integer   :: l_nrows_pm_test,l_ncols_pm_test
    
    integer   :: desc_amat( descriptor_len )
    integer   :: desc_pm_test( descriptor_len )
    integer   :: desc_amatg( descriptor_len )
    integer   :: desc_pm_testg( descriptor_len )    

    integer   ::  rootNodeContext,   sizeN
     
    real(RealPrec), dimension(:,:),allocatable :: amat_loc
    real(RealPrec), dimension(:,:)  ,allocatable :: pm_test_loc
 
!
! Which processor am I and how many others are there?
!
    call blacs_pinfo( myid, nprocs )
    
!
! Broadcast matrix N x N size to all nodes:
! 
    if (myid == 0)  sizeN = size(amat,dim=1)
    call mpi_bcast(sizeN, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        
!    
! Set the dimension of the 2d processors grid:
!
    call gridsetup(nprocs,nprow,npcol)

!
! Initialize a single blacs context:
!
    call blacs_get(-1,0,ictxt) 
    call blacs_gridinit(ictxt,'r',nprow,npcol) 
    call blacs_gridinfo(ictxt, nprow,npcol,myrow,mycol) 
 
!   
! Calculate the blocking factor for the matrix:
!
    call blockset( nb, 64, sizeN, sizeN, nprow, npcol) !
  
!      
! Distributed matrices: get num. local rows/cols. Create description:
!
    l_nrows_amat = numroc(sizeN,nb,myrow,0,nprow)
    l_ncols_amat = numroc(sizeN,nb,mycol,0,npcol)
    call descinit( desc_amat, sizeN, sizeN, nb, nb, 0, 0, ictxt,l_nrows_amat, info ) 

    l_nrows_pm_test = numroc(sizeN,nb,myrow,0,nprow)
    l_ncols_pm_test = 1  
    call descinit( desc_pm_test, sizeN, 1, nb, 1, 0, 0, ictxt, l_nrows_pm_test, info )  
 
    allocate (amat_loc(l_nrows_amat, l_ncols_amat ), stat=istat)
    if (istat /= 0) then
        write(*,*) 'Error allocating amat_loc! l_nrows_amat,l_ncols_amat =',l_nrows_amat,l_ncols_amat
        stop
    endif 
    
    allocate (pm_test_loc(l_nrows_pm_test, l_ncols_pm_test ), stat=istat) 
    if (istat /= 0) then
        write(*,*) 'Error allocating pm_test_loc! l_nrows_pm_test,l_ncols_pm_test =',l_nrows_pm_test,l_ncols_pm_test
        stop
    endif     
! 
! Distribute sub-blocks of global amat and pm_test from master to all process:
!
    call blacs_get( -1, 0, rootNodeContext )
    call blacs_gridinit( rootNodeContext, 'row-major', 1, 1 )
    if (myid == 0) then
        call descinit( desc_amatg, sizeN, sizeN, sizeN, sizeN, 0, 0, rootNodeContext, sizeN, info )
        call descinit( desc_pm_testg, sizeN, 1, sizeN, 1, 0, 0, rootNodeContext, sizeN, info )
    else 
        desc_amatg        = 0
        desc_amatg(2)     = -1
        desc_pm_testg     = 0
        desc_pm_testg(2)  = -1
    endif
    
    call pdgemr2d( sizeN, sizeN, amat, 1, 1, desc_amatg, amat_loc, 1, 1, desc_amat, desc_amat( 2 ) )

    call pdgemr2d( sizeN, 1, pm_test, 1, 1, desc_pm_testg, pm_test_loc, 1, 1, desc_pm_test, desc_pm_test( 2 ) )
    

!
! Parallel matrix solver:
!
    call pdpotrf( 'u',  sizeN, amat_loc, 1, 1, desc_amat, info )
    if (info /= 0) then
        write(*,*) 'Error calling pdpotrf! info: ',info
        stop
    endif      
    
    call pdpotrs( 'u',  sizeN, 1, amat_loc, 1, 1, desc_amat, pm_test_loc, 1, 1, desc_pm_test,  info )
    if (info /= 0) then
        write(*,*) 'Error calling pdpotrs! info: ',info
        stop
    endif  
   
   
!
! Send solution in vector back to master:
!   
    call pdgemr2d( sizeN, 1, pm_test_loc, 1, 1, desc_pm_test, pm_test, 1, 1, desc_pm_testg, ictxt )

    if (myid == 0) call blacs_gridexit( rootNodeContext ) 
    
!
! Deallocate local arrays...
!
    if (allocated(amat_loc))     deallocate(amat_loc)
    if (allocated(pm_test_loc))  deallocate(pm_test_loc)
        
!
! Shut down blacs:
!    
    call blacs_gridexit(ictxt) 
     
    
    end subroutine mpi_worker_ScaLapackMS
#endif


!
!-----------------------------------------------------------------------
!    
          subroutine gridsetup(nproc,nprow,npcol)
!
! This subroutine factorizes the number of processors (nproc)
! into nprow and npcol,  that are the sizes of the 2d processors mesh.
!
! Written by Carlo Cavazzoni
!
      integer nproc,nprow,npcol
      integer sqrtnp,i
 
      sqrtnp = int( sqrt( dble(nproc) ) + 1 )
      do i=1,sqrtnp
        if(mod(nproc,i).eq.0) nprow = i
      end do
      npcol = nproc/nprow
 
      return
      end  subroutine gridsetup
!
!-----------------------------------------------------------------------
!
      subroutine blockset( nb, nbuser, m,n, nprow, npcol)
!
!     This subroutine try to choose an optimal block size
!     for the distributd matrix.
!
!     Written by Carlo Cavazzoni, CINECA
!
      integer :: m, n
      integer nb, nprow, npcol, nbuser
 
      nb = min ( m/nprow, n/npcol )
      if(nbuser.gt.0) then
        nb = min ( nb, nbuser )
      endif
      nb = max(nb,1)
 
      return
      end subroutine blockset  
!==================================================================================================================================! 
!============================================================================================ Non-linear Model Parameter Transforms
!==================================================================================================================================! 
!
! Nonlinear transforms to bound model parameters:
!
! Computes x(m) and m(x) where b < m < a
!
! Kerry Key
! Scripps Institution of Oceanography
!
! March 2011    Implemented with exponential and bandpass options.
!
! Usage:  
!
! m = transformToBound(x,a,b)
! x = transformToUnbound(m)
! call transformWJ() - converts Jacobian matrix WJ using current x
!
! Uses module variables upperBound, lowerBound, cBoundsTransform
!
! cBoundsTransform must be:
!
! 'exponential' is from Habashy and Abubakar(2004); Commer and Newman (2008)
!
! 'bandpass' is one I designed in March 2011, which has a flat pass 
!            band for the sensitivity transform rather than the peak 
!            of the exponential.
!  
!---------------------------------------------------------------------
      elemental real(RealPrec) function transformToBound(xin,bin,ain,lbind)
!---------------------------------------------------------------------
!
! Converts the unbound parameter x to the bound parameter b < m < a    
!        
      implicit none
      real(RealPrec), intent(in) :: xin,ain,bin   
      logical, intent(in)        :: lbind 
      
      real(8) :: a, b, c, p, q, x ! using double no matter what here to avoid exp issues
   
        
      ! Convert input to double variables:
      
        a = ain
        b = bin
        x = xin   
        
        if (lbind) then
            select case (cBoundsTransform)
            
            case ('exponential')  
             
             if ( x.le.0) then
                 x =  ( a*exp( x) + b )  / ( exp( x) + 1.0 )
             else
                 x =  ( b*exp(-x) + a )  / ( exp(-x) + 1.0 )
             endif
             
            case ('bandpass')
             
             c = bandPassFactor / (a-b)
             p = c*(1.0-exp(c*(b-a)))
             
             if (x.le.0) then
                q =  log( (exp(c*x)+exp(c*b)) / (exp(c*x)+exp(c*a)) )
             else
                q =  log( (1.0+exp(c*(b-x))) / (1.0+exp(c*(a-x))) )
             endif       
             x = (a*c + q) / p
             
            end select
        endif        
        transformToBound = x
        
    end function
    
!---------------------------------------------------------------------
    elemental real(RealPrec) function transformToUnbound(mmin,bin,ain,lbind)
!---------------------------------------------------------------------    
! Converts the bound parameter b < m < a to the unbound parameter x    
!
        implicit none
        real(RealPrec), intent(in) :: mmin, ain, bin   
        logical, intent(in)        :: lbind           
        real(8) :: a, b, c, p, m

      ! Convert input to double variables:
       
        a = ain
        b = bin
        m = mmin   

        if (lbind) then
            
            select case (cBoundsTransform)
            
            case ('exponential')  
            
                m = log(m-b) - log(a-m)
            
            case ('bandpass')
             
                 c = bandPassFactor / (a-b)
                 p = m*c*(1.0-exp(c*(b-a)))
                 
                 m =  (log(exp(a*c)*(exp(p)-exp(b*c)))-log(exp(a*c)-exp(p))) / c
             
            end select
        endif
        
        transformToUnbound = m
        
    end function
    
 !---------------------------------------------------------------------   
    subroutine transformWJ()
 !---------------------------------------------------------------------   
 !
 ! Converts WJ from bound parameter m bound parameter x: dF/dx = dm/dx *dF/dm
 !  
 
    implicit none
    
    integer        :: j
    real(8) :: a, b, c, p, q, dmdx, pdp 
 
    
        do j = 1,nParams
            
            if (lBoundMe(j)) then
            
                a = upperBound(j)
                b = lowerBound(j)
                
                pdp = pm(j)
                
                select case (cBoundsTransform)
                
                case ('exponential')  
                    if (pdp.le.0) p = exp( pdp)
                    if (pdp.gt.0) p = exp(-pdp)
                    dmdx = (a-b)*p / (1.0 + p)**2
                    
                    wj(:,j) = wj(:,j) * dmdx
                
                case ('bandpass')
                 
                     c = bandPassFactor / (a-b)
                    
                     if (pdp.le.0)  then
                        p = exp( c*(pdp-b))
                        q = exp( c*(pdp-a))
                        dmdx = p / ( (1.0+p)*(1.0+q) )
                     elseif (pdp.gt.0) then
                        p = exp(-c*(pdp-b)) 
                        q = exp(-c*(pdp-a)) 
                        dmdx = q / ( (1.0+p)*(1.0+q) )
                     endif
                      wj(:,j) = wj(:,j) * dmdx
                  
            
                end select
            endif 
        enddo ! do j = 1,nParams
  
             
    end subroutine transformWJ

!==================================================================================================================================! 
!=========================================================================================================== printOccamIntroMessage 
!==================================================================================================================================! 
    subroutine printOccamIntroMessage  
!
! Messages written at the start of an Occam iteration
!       
    character(128)  :: cStr       
    

    
!
! Say hello:
! 
    cStr ='-----------------------------------------------------------------------------------------------------------------------'
    call printOccamLog(cStr)                
    write(cStr,'(a,i5,a)')  '** Iteration ',currentIteration,' **'
    call printOccamLog(cStr)
    cStr = ' '
    call printOccamLog(cStr)    
    cStr =  ' Constructing derivative dependent matrices...'
    call printOccamLog(cStr)
   
   
!
! Start the iteration timers:
!
    call get_time_offset(0d0,timeIter0)
       
!
! Initialize a few things:
!
    numStepCuts     = 0
    stepCut         = 1. 
    numForwardCalls = 0   
    convergenceFlag = 0
    call initializeTrackers   
    
    end subroutine printOccamIntroMessage 
    
!==================================================================================================================================! 
!=========================================================================================================== printOccamExitMessage 
!==================================================================================================================================! 
    subroutine printOccamExitMessage  
!
! Messages written at the end of an Occam iteration
!    
    character(128)  :: cStr       
    real(RealPrec)  :: timeOccamCumulative
    
    !
    ! Stop the timer:
    !
    call get_time_offset(timeIter0,timeIterEnd)
    call get_time_offset(timeOccamStart,timeOccamCumulative)
    numForwardCallsCumulative = numForwardCallsCumulative + numForwardCalls
    
    !
    ! Print out the summary information about this iteration:
    !  
    
    cStr = ' '
    call printOccamLog(cStr)    
    
    if ( convergenceFlag /= 4 .and. convergenceFlag /= 5 ) then    
        cStr = ' Occam iteration completed successfully: '
        call printOccamLog(cStr)    
        cStr = ' '
        call printOccamLog(cStr)      
        write(cStr,'(a32,g16.4)') 'Target Misfit: ',targetRMS
        call printOccamLog(cStr)         
        write(cStr,'(a32,g16.4)') 'Model Misfit: ',modelRMS
        call printOccamLog(cStr)      
        write(cStr,'(a32,g16.4)')  'Roughness: ',modelRoughness
        call printOccamLog(cStr)     
        write(cStr,'(a32,g16.4)') 'Optimal Mu: ',modelMu
        call printOccamLog(cStr)           
        write(cStr,'(a32,g16.4)') 'Stepsize: ',stepSize
        call printOccamLog(cStr)
        write(cStr,'(a32,7x,i1)') 'Convergence Status: ',convergenceFlag
        call printOccamLog(cStr)      
        write(cStr,'(a32,4x,2(i4,1x))') ' # Forward Calls, Cumulative: ', numForwardCalls,numForwardCallsCumulative
        call printOccamLog(cStr)    
        write(cStr,'(a32,2(g16.4,1x))') 'Iteration Time, Cumulative (s): ',timeIterEnd,timeOccamCumulative
        call printOccamLog(cStr)
        cStr = ' '
        call printOccamLog(cStr) 
    endif    
    !
    ! Special messages based on convence status:
    !
    
    select case (convergenceFlag)
    case (0) ! Do  nothing since target misfit not yet obtained
    case (1) ! Do  nothing since Occam is still smoothing the model
    case (2)
        cStr = ' Stopping on normal convergence.'
        call printOccamLog(cStr)
    case (3)    
        cStr = ' Perfectly smooth model found, stopping.'
        call printOccamLog(cStr)         
    case (4)
        cStr =  ' RMS misfit can not be decreased further, stopping.'
        call printOccamLog(cStr)
        cStr =  ' May have reached minimum possible RMS misfit level!'
        call printOccamLog(cStr)        
    case (5)
        cStr =  ' Model roughness can not be decreased further, stopping. '
        call printOccamLog(cStr)  
    case(6) 
        cStr =  ' Maximum number of iterations achieved, stopping. '
        call printOccamLog(cStr)   
    end select   


                
    
    end subroutine printOccamExitMessage 

!==================================================================================================================================! 
!================================================================================================================= openOccamLogFile
!==================================================================================================================================! 
    subroutine openOccamLogFile(cFileRoot)
!
! Opens the Occam log file, appending to an existing file if the inversion has been restarted from an intermediate result.
!
    integer        :: ierr 
    logical        :: lFileExists
    character(256) :: cFileName, cFileRoot
    
    character(80)  :: dateStamp
    character(180) :: cStr
          
    
    
    ! The log file name:
    cFileName = trim(cFileRoot)//'.logfile'
    
    ! See if it exists:
    inquire(file=cFileName,exist=lFileExists)
     
     
    ! If the file exists then append the new log to it: 
    if ( (currentIteration > 0).and.(lFileExists) ) then
    
        open (ioUnitOccamLogFile, file=cFileName, position='append',iostat=ierr)
      
      ! Drop a line letting the user know this is a restart:
     cStr='-----------------------------------------------------------------------------------------------------------------------'
   
        write(ioUnitOccamLogFile,'(a)') cStr
        write(ioUnitOccamLogFile,'(a)') cStr
        write(ioUnitOccamLogFile,'(a)') cStr
        write(ioUnitOccamLogFile,*) ' *** Occam restarted here from iteration: ',currentIteration
        write(ioUnitOccamLogFile,*) ' '
        
    else  ! Open a brand spanking new file, woohoo!  
        open (ioUnitOccamLogFile, file=cFileName, iostat=ierr)
        write(ioUnitOccamLogFile,'(a)') 'Format:     OccamLog.2012.0'
        
    endif
    if (ierr .ne. 0) then
        write(cStr,'(a,a)') ' Error opening log file:', trim(cFileName)
        call printOccamLog(cStr) 
        stop 
    end if

    !
    ! Insert a time stamp:
    !
    call getDateStamp(dateStamp)  
    write(cStr,'(a12,a)')  'Start time: ',trim(adjustl(dateStamp))
    call printOccamLog(cStr) 
    cStr = ' ' 
    call printOccamLog(cStr)        
    
    !
    ! Print number of processors:
    !
    call mpi_comm_size( MPI_COMM_WORLD, nParallelProc, ierr )   
    write(cStr,'(a32,i8)') ' Number of parallel processors:',  nParallelProc
    call printOccamLog(cStr) 
    
    !
    ! Write number of data and model parameters:
    !
    write(cStr,'(a32,i8)') ' Number of input data:',  nd
    call printOccamLog(cStr) 
    write(cStr,'(a32,i8)') ' Number of model parameters:', nParams
    call printOccamLog(cStr)  
    !
    ! Also initialize a few things for the Occam run:
    !
    
    !
    ! Start the timer for the entire program:
    !
    call get_time_offset(0d0,timeOccamStart) 
    
    numForwardCallsCumulative = 0
    
    end subroutine openOccamLogFile
    
!==================================================================================================================================! 
!================================================================================================================ closeOccamLogFile
!==================================================================================================================================!         
    subroutine closeOccamLogFile 
    
    character(80)  :: dateStamp
    character(180) :: cStr
    
      
    call getDateStamp(dateStamp)  
    
    cStr = ' ' 
    call printOccamLog(cStr)
    
    write(cStr,'(a32,a)')  ' End time: ',trim(adjustl(dateStamp)) 
    call printOccamLog(cStr)
    
    call get_time_offset(timeOccamStart,TimeOccamEnd) 
    write(cStr,'(a32,g16.4)')  ' Total time for all iterations: ',TimeOccamEnd
    call printOccamLog(cStr)

           
    close(ioUnitOccamLogFile)    
    
    end subroutine closeOccamLogFile 
    
!==================================================================================================================================! 
!===================================================================================================================== getDateStamp
!==================================================================================================================================!         
    subroutine getDateStamp(dateStamp)    
    
    character(80), intent(out)  :: dateStamp
     
    character(10)  :: cdate,ctime
    
    !
    ! Standard intrinsic Fortran90 data_and_time call:
    !
    call date_and_time(cdate,ctime) !this ouputs only values, ignoring other optional arguments
    
    dateStamp = cdate(5:6)//'/'//cdate(7:8)//'/'//cdate(1:4)//' '// ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:) 
    
    end subroutine getDateStamp        
!==================================================================================================================================! 
!==================================================================================================================== printOccamLog
!==================================================================================================================================! 
    subroutine printOccamLog(cStr)
!
! Writes messages to the Occam log file and also optionally to the terminal.
!    
    character(128), intent(in) :: cStr
    
    logical                    :: lopen
    
    inquire(ioUnitOccamLogFile,opened=lopen)    
     
    if (lopen) write(ioUnitOccamLogFile,'(a)') trim(cStr)  ! since sometimes this is called on forward runs when no log file exists
    
    if ( occamPrintLevel > 0)  write(*,'(a)') trim(cStr)
    
    
    end subroutine printOccamLog
!==================================================================================================================================! 
!================================================================================================================== get_time_offset
!==================================================================================================================================! 
    subroutine  get_time_offset(timein,timeout)

!    
! timein is the clock start time or 0.
! timeout is the time in seconds since the input time
!
! Version 2.0  February 25, 2008  Now uses date_and_time Fortran intrinsic
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    
    integer, dimension(8) :: values
    integer               :: i,j,k,mjd
    
    real(8) :: timein, timeout, fracday
    
    !
    ! Fortran90 Time function:
    !
    call date_and_time(values=values) !this ouputs only values, ignoring other optional arguments
    
    ! Convert year, month day to modified julian day:
      
    i = values(1)
    j = values(2)
    k = values(3)
    mjd = -678927 + k + 1461*(i+(j-14)/12)/4 + 367*(j-2-12 * &
          & ((j-14)/12))/12 + (24002-12*i-j)/1200
    
               
    ! Add on fractional day:
                ! hour            ! minute          ! sec       ! millisec
    fracday = ((values(5)*60.d0 + values(6))*60.d0 + values(7) + values(8)/1000.d0 )/86400.d0
    
    timeout = mjd + fracday
    
    ! Finally, convert timeout to time difference between it and timein:  
    timeout =  timeout*86400.d0  - timein               
           
    end subroutine  get_time_offset

!==================================================================================================================================! 
!============================================================================================== Modified Numerical Recipes Routines
!==================================================================================================================================! 
! The routines below are modified from Numerical Recipes in Fortran 90 (Press et al.)
!
!-----------------------------------------------------------------------------------------------------------------------------------    
    subroutine bracketMinimum(ax,bx,cx,func)
!-----------------------------------------------------------------------------------------------------------------------------------    
! Modified from NR's mnbrak() routine
!
! Given a function func, and given distinct initial points ax and bx, this routine searches in the downhill direction (defined by 
! the function as evaluated at the initial points) and returns new points ax, bx, cx that bracket a minimum of the function. Also 
! returned are the function values at the three points, fa, fb, and fc.
! Parameters: GOLD is the default ratio by which successive intervals are magnified; GLIMIT
! is the maximum magnification allowed for a parabolic-fit step.    
    
 
    real(RealPrec), intent(inout)   :: ax,bx
    real(RealPrec), intent(out)     :: cx
    real(RealPrec), external        :: func
    
    real(RealPrec)                  :: fa,fb,fc
    real(RealPrec), parameter       :: gold=1.618034, glimit=100.0, tiny=1.0e-20
    real(RealPrec)                  :: fu,q,r,u,ulim
    character(128)                  :: cStr    
 
 
    cStr = ' ... Bracketing minimum:'
    call printOccamLog(cStr)
    write(cStr,'(7(a16,1x))') '     Misfit     ','   Roughness    ','   Log10(mu)    ', &
                            & ' Min log10(rho) ',' Max log10(rho) ', '  Fwd Call (s) ', ' Matrix Ops. (s)'
    call printOccamLog(cStr)
        
    
    fb=func(bx)
    call updateB(bx,fb)
    if (lEscapeFromMinimization(fb)) return
            
    fa=func(ax)
    
    if ( fa < fb) then   ! Switch roles of a and b so that we can go downhill in the direction from a to b.
        call updateB(ax,fa)
        call swap(ax,bx)
        call swap(fa,fb)
        if (lEscapeFromMinimization(fb)) return
    end if
    
    cx=bx+gold*(bx-ax)  ! First guess for c.
    fc=func(cx)
    
    do                  ! Do-while-loop: Keep returning here until we bracket.
        
        if (fb < fc) return  ! ie fa > fb < fc, so we've bracketed a minimum
        
        call updateB(cx,fc)  ! cx/fc has lower misfit, so save it
        if (lEscapeFromMinimization(fc)) then
            bx=cx
            fb=fc
            return
        endif
        
        ! Compute u by parabolic extrapolation from a, b, c. TINY is used to prevent any possible division by zero.
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*sign(max(abs(q-r),tiny),q-r))
        ulim=bx+glimit*(cx-bx)
        ! We won’t go farther than this. Test various possibilities:
        if ((bx-u)*(u-cx) > 0.0) then               ! Parabolic u is between b and c: try it.
            fu=func(u)
            if (fu < fc) then                       ! Got a minimum between b and c.
                ax=bx
                fa=fb
                bx=u
                fb=fu
                call updateB(bx,fb) ! u/fu model saved
                return
            else if (fu > fb) then                  ! Got a minimum between a and u.
                cx=u
                fc=fu
                return
            end if
            u=cx+gold*(cx-bx)                       ! Parabolic fit was no use. Use default magnification.
            fu=func(u)
        else if ((cx-u)*(u-ulim) > 0.0) then        ! Parabolic fit is between c and its allowed limit.
            fu=func(u)
            if (fu < fc) then
                call updateB(u,fu) ! u/fu model saved
                bx=cx
                cx=u
                u=cx+gold*(cx-bx)   
                call shft(fb,fc,fu,func(u))  
            end if
        else if ((u-ulim)*(ulim-cx) >= 0.0) then    ! Limit parabolic u to maximum allowed value.
            u=ulim
            fu=func(u)
        else                                        ! Reject parabolic u, use default magnification.
            u=cx+gold*(cx-bx)
            fu=func(u)
        end if
        call shft(ax,bx,cx,u)                       ! Eliminate oldest point and continue.
        call shft(fa,fb,fc,fu)
    end do
 
    end subroutine bracketMinimum
    
!-----------------------------------------------------------------------------------------------------------------------------------    
    subroutine findMinimum(ax,bx,cx,fbx,func)
!-----------------------------------------------------------------------------------------------------------------------------------    
! Modified from NR's brent() routine
!
    real(RealPrec), intent(in)  :: ax,bx,cx,fbx
    !real(RealPrec), intent(out) :: xmin,rmsmin
    real(RealPrec), external    :: func
        
    real(RealPrec)              :: tol = 0.1  ! isolates minimum to about this tolerance     
    integer, parameter          :: itmax=100
    real(RealPrec) , parameter  :: cgold=0.3819660,zeps=1.0e-3*epsilon(ax)
    integer                     :: iter
    real(RealPrec)              :: a,b,dd,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    character(128)              :: cStr

    cStr = ' ... Finding minimum:'
    call printOccamLog(cStr)
    write(cStr,'(7(a16,1x))') '     Misfit     ','   Roughness    ','   Log10(mu)    ', &
                            & ' Min log10(rho) ',' Max log10(rho) ', '  Fwd Call (s) ', ' Matrix Ops. (s)'
    call printOccamLog(cStr)        

    a  = min(ax,cx)
    b  = max(ax,cx)
    v  = bx
    w  = v
    x  = v
    e  = 0.0
    fx = fbx
    fv = fx
    fw = fx
    do iter=1,itmax
        xm   = 0.5*(a+b)
        tol1 = tol*abs(x)+zeps
        tol2 = 2.0*tol1
        if (abs(x-xm) <= (tol2-0.5*(b-a))) then 
!           xmin=x
!           rmsmin=fx
            return
        end if
        if (abs(e) > tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.0*(q-r)
            if (q > 0.0) p = -p
            q     = abs(q)
            etemp = e
            e     = dd
            if (abs(p) >= abs(0.5*q*etemp) .or. &
                p <= q*(a-x) .or. p >= q*(b-x)) then
                e  = merge(a-x,b-x, x >= xm )
                dd = cgold*e
            else
                dd = p/q
                u  = x+dd 
                if (u-a < tol2 .or. b-u < tol2) dd = sign(tol1,xm-x)
            end if
        else
            e = merge(a-x,b-x, x >= xm )
            dd = cgold*e
        end if
        u=merge(x + dd, x+sign(tol1,dd), abs(dd) >= tol1 )
        fu=func(u)    ! func call 
        if (fu <= fx) then
            call updateB(u,fu)
            if (lEscapeFromMinimization(fu)) then
                return ! return if current misfit is good enough
            endif
            if (u >= x) then
                a=x
            else
                b=x
            end if
            call shft(v,w,x,u)
            call shft(fv,fw,fx,fu)
        else
            if (u < x) then
                a=u
            else
                b=u
            end if
            if (fu <= fw .or. w == x) then
                v=w
                fv=fw
                w=u
                fw=fu
            else if (fu <= fv .or. v == x .or. v == w) then
                v=u
                fv=fu
            end if
        end if
        ! David Myer's edit that ends the minimization based on misfit flattening, rather than by mu changing below the tolerance:
        if  (abs(fw-fx) < rmsTol ) return

    end do
        
    end subroutine findMinimum      
    
    
!-----------------------------------------------------------------------------------------------------------------------------------    
    subroutine findIntercept(func,xintercept )
!-----------------------------------------------------------------------------------------------------------------------------------    
! Modified from NR's zbrent() routine
! 
 
    implicit none
    
    real(RealPrec), intent(out) :: xintercept 
 
    real(RealPrec), external    :: func
    
    real(RealPrec)              :: fa,fb,rb  
    
    integer, parameter          :: itmax=100
    real(RealPrec),parameter    :: tol = 0.001 
    real(RealPrec),parameter    :: tolRMS = 0.005 ! RMS tolerance for intercept, stop if at least this close to target
    real(RealPrec), parameter   :: eps=epsilon(tolRMS)
    
    integer                     :: iter
    real(RealPrec)              :: a,b,c,dd,e,fc,p,q,r,s,tol1,xm
    
    character(180)              :: cStr

    
    a  = mu_a !x1
    b  = mu_b !x2
    fa = rms_a - targetRMS ! func(a)
    fb = rms_b - targetRMS ! func(b) 
    
    if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
        write(cStr,'(a)')  'Error in subroutine findIntercept:  root must be bracketed on entry!'
        call printOccamLog(cStr)   
    endif   
    c  = b
    fc = fb
    
    ! Update arrays keep track of models near intercept:
    rms_c   = rms_b
    mu_c    = mu_b
    pm_c    = pm_b   
    dm_c    = dm_b
    rough_c = rough_b
    if (npm_assoc > 0) pm_assoc_c = pm_assoc_b  
    
    do iter=1,itmax
    
        if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
            c  = a
            fc = fa
            dd = b-a
            e  = dd
            
            ! a was moved to c: 
            ! Update arrays keep track of models near intercept:
            rms_c   = rms_a
            mu_c    = mu_a
            pm_c    = pm_a   
            dm_c    = dm_a
            rough_c = rough_a
            if (npm_assoc > 0) pm_assoc_c = pm_assoc_a  

        end if
        
        if (abs(fc) < abs(fb)) then
            a  = b
            b  = c
            c  = a
            fa = fb
            fb = fc
            fc = fa
                
            ! b moved to a, c moved to b, a moved to c (so a=c...)
            ! Update arrays keep track of models near intercept:
            rms_a   = rms_b
            mu_a    = mu_b
            pm_a    = pm_b   
            dm_a    = dm_b
            rough_a = rough_b
            if (npm_assoc > 0) pm_assoc_a = pm_assoc_b  
            
            rms_b   = rms_c
            mu_b    = mu_c
            pm_b    = pm_c   
            dm_b    = dm_c
            rough_b = rough_c
            if (npm_assoc > 0) pm_assoc_b = pm_assoc_c
            
            rms_c   = rms_a
            mu_c    = mu_a
            pm_c    = pm_a   
            dm_c    = dm_a
            rough_c = rough_a
            if (npm_assoc > 0) pm_assoc_c = pm_assoc_a
                            
        end if
        tol1=2.0*eps*abs(b)+0.5*tol
        xm=0.5*(c-b)
        if (abs(xm) <= tol1 .or. abs(fb) <= tolRMS) then ! fb == 0.0) then
            xintercept = b          
            return
        end if
        if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
            s=fb/fa
            if (a == c) then
                p=2.0*xm*s
                q=1.0-s
            else
                q=fa/fc
                r=fb/fc
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
                q=(q-1.0)*(r-1.0)*(s-1.0)
            end if
            if (p > 0.0) q=-q
            p=abs(p)
            if (2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
                e  = dd
                dd = p/q
            else
                dd = xm
                e  = dd
            end if
        else
            dd = xm
            e = dd
        end if
        a  = b
        fa = fb     
        
        ! b moved to a 
        ! Update arrays keep track of models near intercept:
        rms_a   = rms_b
        mu_a    = mu_b
        pm_a    = pm_b   
        dm_a    = dm_b
        rough_a = rough_b
        if (npm_assoc > 0) pm_assoc_a = pm_assoc_b  
                    
        b = b+merge(dd,sign(tol1,xm), abs(dd) > tol1 )
        rb = func(b)   ! func call
        call updateB(b,rb)  
        fb = rb - targetRMS 
    
    end do
    
    cStr = 'Error in subroutine findIntercept: exceeded maximum iterations trying to find intercept'
    call printOccamLog(cStr)    
    xintercept=b
   
    end subroutine findIntercept
    
!-----------------------------------------------------------------------------------------------------------------------------------    
    subroutine shft(a,b,c,d)
!-----------------------------------------------------------------------------------------------------------------------------------    
        real(RealPrec), intent(out)   :: a
        real(RealPrec), intent(inout) :: b,c
        real(RealPrec), intent(in)    :: d
        a=b
        b=c
        c=d
    end subroutine shft

!-----------------------------------------------------------------------------------------------------------------------------------        
    subroutine swap(a,b)
!-----------------------------------------------------------------------------------------------------------------------------------    
    real(RealPrec), intent(inout) :: a,b 
    real(RealPrec)                :: dum
        dum=a
        a=b
        b=dum
    end subroutine swap
        
!-----------------------------------------------------------------------------------------------------------------------------------        
    subroutine updateB(mu,rms)
!-----------------------------------------------------------------------------------------------------------------------------------    
! Helper routine for updating the *_b arrays that keep track of the optimal model and its response.
!      
    real(RealPrec), intent(in) :: mu, rms
 
    pm_b    = pm_test
    dm_b    = dm
    rms_b   = rms
    mu_b    = mu 
    rough_b = rough_test
    if (npm_assoc > 0)  pm_assoc_b = pm_assoc 
            
    end subroutine updateB       
    
!-----------------------------------------------------------------------------------------------------------------------------------
end module Occam
