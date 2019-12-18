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

 
!==================================================================================================================================! 
!=================================================================================================================== mare2dem_global
!==================================================================================================================================! 
    module mare2dem_global
 
    use EM_constants
 
    use triangle_mesh  ! for trimesh structure

!
!  Variables read in from the .resistivity file that are stored here  
!
    character(256)  :: resistivityFile  = ''     
    character(256)  :: outputFileRoot   = ''     
    character(256)  :: modelFileName    = ''      
    character(256)  :: penaltyFileName  = ''          
    character(256)  :: dataFileName     = ''
    character(256)  :: settingsFileName = ''                 
    real(8)         :: lowerBoundGlobal, upperBoundGlobal  ! There are used to initialize the bound arrays from module Occam
    character(256)  :: cParamHeader     = ''    ! This is not used by MARE2D, just passes thru to the output files.


    integer                               :: nRhoPerRegion = 1  ! 1, 2 or 3 depending on anisotropy   
    integer                               :: nRegions               
    integer                               :: nFree = 0         ! Number of conductivities in sigParams that are free parameters     
    real(8), dimension(:,:), allocatable  :: boundsTemp, PrejTemp, ratioTemp  ! input arrays, keep them around for easy output 
    
! 
! Variables that can be changed using optional inputs in the settings file:
!
    logical         :: lprintDecomposition     = .false.  ! set to true to display the parallel decomposition settings
    logical         :: lDisplayRefinementStats = .false.  ! set to true to show adaptive refinement iteration status                                              
    logical         :: lprintMPItimers         = .false.   ! set to true to display timing info on MPI send and recv commands
    logical         :: lPrintDebug             = .false.   
    logical         :: lPrintSetup             = .false.   
    
    integer         :: maxnadapt  = 30         ! maximum number of refinements        
    integer         :: nwave      = 30         ! # of wavenumbers for 2.5D Fourier transformation. 
    real(8)         :: loglower   = -5         ! log10 lower limit of wavenumbers
    real(8)         :: logupper   = -1         ! log10 upper limit of wavenumbers
    
    character(256)  :: scratchFolder  = '/tmp'    
   
    !
    ! Group parameters for data decomposition parallelization:
    !
    integer         :: nTxPerGroup       = 10     ! number of transmitters to group together for mesh refinement 
    integer         :: nRxPerGroupCSEM   = 40     ! number of receivers to group together for mesh refinement
    integer         :: nRxPerGroupMT     = 40     ! number of receivers to group together for mesh refinement        
    integer         :: nKxPerGroup       = 5      ! number of wavenumbers to group together for mesh refinement
    integer         :: nFreqPerGroupCSEM = 1      ! number of frequencies to group together for mesh refinement
    integer         :: nFreqPerGroupMT   = 1      ! number of frequencies to group together for mesh refinement
    
    !
    ! For finite Tx and Rx dipoles we need to specify the quadrature order to use (can be overridden by mare2dem.settings file):
    !   
    ! * Use an odd number if in addition to E wires, you are modeling point dipole magnetics since the 
    !   the code uses the central value for B point dipole... *
    !
    integer         :: nQuadTxCSEM = 3, nQuadRxCSEM = 3,  nQuadRxMT = 3          

!
! Variables set from command line arguments:
!
    logical         :: lFwdFields  = .false.  ! output all possible field components to a file. use MARE2DEM -FF to make true.  
     
 
!-----------------------------------------------------------------------------------------------------------------------------------
! Everything below here is an internal variable only
!    

!
! MPI worker status array ( true means ready for next task)
!   
    integer                                 :: nworkers
    logical,dimension(:), allocatable       :: lworker_status 
     
    type(trimesh)                           :: inputmodel, inputmesh ! see call_triangle.f90 for the definition of derived type


    integer                                 :: nMT, nCSEM            ! number of MT and CSEM data    
   
    logical, dimension(:,:,:), allocatable  :: lDataMaskCSEM ! nRxCSEM x nFreqCSEM x nTx. true where data file has data    
    logical, dimension(:,:),   allocatable  :: lDataMaskMT   ! nRxMT   x nFreqMT. true where data file has data    
       
         

    logical                                 :: linversion  = .false.  
    logical, dimension(:,:),allocatable     :: lCompCSEM, lCompMT  ! (6,nRx), request fields to compute (ex ey ez hx hy hz)
    ! Note this corresponds to input data code but the rx could be rotated/tilted, so in reality we need to compute say ex,ey,ez

    real(8), dimension(:),allocatable       :: wavenum    ! array of wavenumbers 
                       
!
! Arrays for RxTx and Frequency groups:
!
    type  :: RxTxG   ! Derived type to store each checkboard subgrid of Rx and Tx combinations:
        integer                             :: nTx      ! number of transmitters in this group 
                                                       ! Note this can be larget the nTxPerGroup to handle modal transmitters
        integer, dimension(:), allocatable  :: iTx      ! index of global transmitter number, used to get x,y,z of tx
        logical                             :: lMT      ! true if this is an MT RxTxG  s
        integer                             :: nRx
        integer, dimension(:), allocatable  :: iRx      ! index of global receiver numbers           
    end type
 
    integer                                 :: nRxTxGroups
    type(RxTxG), dimension(:), allocatable  :: RxTxGroups   
    
    type  :: FqG   ! Derived type to store each checkboard subgrid of Rx and Tx combinations:
        integer                             :: nFq      ! number of frequencies in this group   
        integer, dimension(:), allocatable  :: iFq      ! index of global frequencies for this group
        integer                             :: iFqRefine  ! local index to iFq for refinement frequency for this group         
        logical                             :: lMT      ! true if this is an MT FqG       
    end type
 
    integer                                 :: nFqGroups 
    type(FqG), dimension(:), allocatable    :: FqGroups   

    type  :: RefG   ! Derived type to store adaptively refined mesh and settings for each refinement group 
        integer                             :: iRxTxGroup    ! index to the group of Rx and Tx
        integer                             :: iFqGroup      ! index to the group of Kx and Fq
        logical                             :: lMT           ! true if this is an MT RefG    
        logical                             :: lReturned = .false. ! true when returned from worker                        
    end type
    
    integer                                 :: nRefinementGroups   
    integer                                 :: iPtr_refGroups   ! points to next refinement group to send out
    type(RefG), dimension(:), allocatable   :: refinementGroups ! nRefinementGroups rows of  iRxTxGroup , iFqGroup

   
!
!  A few variables used by the worker nodes:
!
     integer     :: iRefinementGrp

    
    contains
    
    subroutine deallocate_mare2dem_global  

    if ( allocated( boundsTemp ) )          deallocate( boundsTemp )
    if ( allocated( PrejTemp ) )            deallocate( PrejTemp )
    if ( allocated( ratioTemp ) )           deallocate( ratioTemp ) 
    if ( allocated( lworker_status ) )      deallocate (lworker_status)
    if ( allocated( inputmodel%attr ) )     call deallocate_trimesh(inputmodel,.false.)
    if ( allocated( inputmesh%attr ) )      call deallocate_trimesh(inputmesh,.false.)
    if ( allocated( lDataMaskCSEM ) )       deallocate( lDataMaskCSEM )
    if ( allocated( lDataMaskMT ) )         deallocate( lDataMaskMT ) 
    if ( allocated( wavenum ) )             deallocate( wavenum ) 
    if ( allocated( RxTxGroups ) )          deallocate( RxTxGroups )
    if ( allocated( FqGroups ) )            deallocate( FqGroups ) 
    if ( allocated( refinementGroups ) )    deallocate( refinementGroups ) 
    if ( allocated( lCompCSEM ) )           deallocate( lCompCSEM )
    if ( allocated( lCompMT ) )             deallocate( lCompMT )
        
    end subroutine deallocate_mare2dem_global    
        
    end module mare2dem_global   
 
!==================================================================================================================================! 
!======================================================================================================== mare2dem_input_data_params
!==================================================================================================================================! 
    module mare2dem_input_data_params
    
! While the data parameter array dp is passed in/out of Occam, it is only used
! in the forward routines.  So we set cnDataParams in the forward routines
! and then allocate the dp array in readData.
!        
    integer         :: cnDataParams
 
    character(256)  :: cUTMLine        = ''         ! the UTM line from the data file (not used by MARE2D, just passes thru)

    character(4)    :: phaseConvention = 'lag'      ! 'lag' or 'lead'
    character(4)    :: reciprocityUsed = ' '        ! 'yes' or 'no' (not used, just passes through for later plotting ease)
    
!
! Data dependent Tx-Rx variables read from data file in subroutine readData:
!
    ! CSEM sites:
    integer                                   :: nTxCSEM=0, nRxCSEM=0, nFreqCSEM=0
    real(8), dimension(:), allocatable        :: azimuthTxCSEM,dipTxCSEM,lengthTxCSEM
    real(8), dimension(:), allocatable        :: xTxCSEM,yTxCSEM,zTxCSEM,xRxCSEM,yRxCSEM,zRxCSEM        
    real(8), dimension(:), allocatable        :: ThetaRxCSEM,AlphaRxCSEM, BetaRxCSEM  
    real(8), dimension(:), allocatable        :: lengthRxCSEM      ! nRx  for x,y,z dipoles all same length. Use new Rx if otherwise 
    real(8), dimension(:), allocatable        :: fTxCSEM  
    

    character(8),   dimension(:), allocatable :: cSourceType  
    character(128), dimension(:), allocatable :: cRxNamesCSEM
    character(128), dimension(:), allocatable :: cTxNamesCSEM  
    ! MT sites:
    integer                                   :: nRxMT=0, nFreqMT=0
    real(8), dimension(:), allocatable        :: xRxMT,yRxMT,zRxMT   
    real(8), dimension(:), allocatable        :: ThetaRxMT,AlphaRxMT, BetaRxMT, lengthRxMT      
    real(8), dimension(:), allocatable        :: fTxMT  
    integer, dimension(:), allocatable        :: iSolveStatic  ! 1 = TE + TM, 2 = TE only, 3 = TM only. otherwise ignored
   
    character(128), dimension(:), allocatable :: cRxNamesMT

    contains
    
    subroutine deallocate_mare2dem_input_data_params
    
    if (allocated (fTxCSEM))        deallocate( fTxCSEM ) 
    if (allocated (xTxCSEM) )       deallocate( xTxCSEM, yTxCSEM, zTxCSEM, azimuthTxCSEM, dipTxCSEM, cSourceType, lengthTxCSEM )  
    if (allocated (cRxNamesCSEM) )  deallocate( cRxNamesCSEM )
    if (allocated (cTxNamesCSEM) )  deallocate( cTxNamesCSEM )

    if (allocated (fTxMT))          deallocate( fTxMT )    
    if (allocated (xRxCSEM))        deallocate( xRxCSEM, yRxCSEM, zRxCSEM, ThetaRxCSEM, AlphaRxCSEM, BetaRxCSEM, lengthRxCSEM )   
    if (allocated (xRxMT))          deallocate( xRxMT, yRxMT, zRxMT, ThetaRxMT, AlphaRxMT, BetaRxMT, lengthRxMT )         
    if (allocated (cRxNamesMT) )    deallocate( cRxNamesMT )
    if (allocated (iSolveStatic) )  deallocate( iSolveStatic )
        
    
    end subroutine deallocate_mare2dem_input_data_params
        
    end module mare2dem_input_data_params


!==================================================================================================================================! 
!=================================================================================================================== mare2dem_output   
!==================================================================================================================================! 
    module mare2dem_output     
    
    integer, parameter, public :: ZPrec = 4  ! note that this is single precision to reduce the scratchfile disk footprint
!
! CSEM fields at (isite,iTx,ifreq):
!
    complex(ZPrec),dimension(:,:,:), allocatable    :: ex,ey,ez,hx,hy,hz
!
! Partial derivatives with respect to sigma:
!
 
    type :: dsigma
         complex(ZPrec),dimension(:), allocatable   :: dsig      ! (iparam)
    end type dsigma   
    type(dsigma),dimension(:,:,:), allocatable      :: dex,dey,dez,dhx,dhy,dhz ! (isite,iTx,ifreq)
 

    complex(ZPrec)                                  :: zte,ztm, htipper   
    complex(ZPrec), dimension(:), allocatable       :: dzte_dsig,dztm_dsig,dhtipper_dsig  ! (isite,ifreq,iparam)
!   
! MT fields at (isite,ifreq):
!
! ex, hy and hz are for the TE mode
! hx, ey and ez are for the TM mode 
!
    complex(ZPrec), dimension(:,:), allocatable     :: ex_mt,ey_mt,ez_mt,hx_mt,hy_mt,hz_mt ! (isite,ifreq)
    complex(ZPrec),dimension(:,:,:), allocatable    :: dex_mt_dsig,dey_mt_dsig,dez_mt_dsig ! (isite,ifreq,iparam)
    complex(ZPrec),dimension(:,:,:), allocatable    :: dhx_mt_dsig,dhy_mt_dsig,dhz_mt_dsig 
 
    
    contains
    
    subroutine deallocate_mare2dem_output
    
    if ( allocated(ex) )                deallocate (ex,ey,ez,hx,hy,hz)
    if ( allocated(dex) )               deallocate (dex)
    if ( allocated(dey) )               deallocate (dey)
    if ( allocated(dez) )               deallocate (dez)
    if ( allocated(dhx) )               deallocate (dhx)
    if ( allocated(dhy) )               deallocate (dhy) 
    if ( allocated(dhz) )               deallocate (dhz)       
        
    if ( allocated(dzte_dsig) )         deallocate (dzte_dsig)  
    if ( allocated(dztm_dsig) )         deallocate (dztm_dsig)  
    if ( allocated(dhtipper_dsig) )     deallocate (dhtipper_dsig)  
    
    if ( allocated(ex_mt) )             deallocate (ex_mt,ey_mt,ez_mt,hx_mt,hy_mt,hz_mt)    
    if ( allocated(dex_mt_dsig) )       deallocate (dex_mt_dsig)
    if ( allocated(dey_mt_dsig) )       deallocate (dey_mt_dsig)
    if ( allocated(dez_mt_dsig) )       deallocate (dez_mt_dsig)
    if ( allocated(dhx_mt_dsig) )       deallocate (dhx_mt_dsig)
    if ( allocated(dhy_mt_dsig) )       deallocate (dhy_mt_dsig)
    if ( allocated(dhz_mt_dsig) )       deallocate (dhz_mt_dsig)
                        
    
    end subroutine deallocate_mare2dem_output
        
    end module mare2dem_output     
 
!==================================================================================================================================! 
!=============================================================================================================== mare2dem_deallocate
!==================================================================================================================================!     
    subroutine mare2dem_deallocate
!
! Deallocates all I/O variables for MARE2DEM and a few others. This routine is called at the end of Occam.
!        
    use em2dkx_mod 
    use mare2dem_output
    use mare2dem_global
    use mare2dem_input_data_params
    
     
    call deallocate_em2dkx 
    
    call deallocate_mare2dem_output ! should already be deleted at this point, but just in case....
     
    call deallocate_mare2dem_global
    
    call deallocate_mare2dem_input_data_params

        
    end subroutine mare2dem_deallocate
        

!==================================================================================================================================! 
!=========================================================================================================== display_MARE2DEM_Params
!==================================================================================================================================!     
    subroutine display_MARE2DEM_Params
 
    use em2dkx_mod
    use mare2dem_global
    
    implicit none
    
    character(80)   :: ctemp

    write(6,*) '========== Parameters settings for MARE2DEM =============================='
    write(*,*) ' '     
 
 
       
    write(*,*) ' '
    write(*,*) ' ----------------------------'
    write(*,*) '  Adaptive refinement params:'
    write(*,*) ' ----------------------------'    
    write(ctemp,'(i5)') maxnadapt     
    write(6,fmt='(a32,a)') 'Max # refinements:  ',trim(adjustl(ctemp))         
    
    write(ctemp,'(f6.3)') errortolerance
    write(*,fmt='(a32,a)') 'Tolerance (%):  ',trim(adjustl(ctemp))     

    
    write(ctemp,'(f5.1)') minQangle 
    write(6,fmt='(a32,a)') 'Mesh quality angle:  ',trim(adjustl(ctemp))      
   
    
    write(ctemp,'(i5)') idual_func      
    write(*,fmt='(a32,a)') 'Dual function:  ', trim(adjustl(ctemp))         
 
    if (lUseBumpFields) then
        write(6,fmt='(a32,a3)') 'Use bump fields:  ','yes'  
    else
        write(6,fmt='(a32,a2)') 'Use bump fields:  ','no' 
    endif  
 
    if (lSaveMeshFiles) then
        write(6,fmt='(a32,a3)') 'Save meshes:  ','yes'  
    else
        write(6,fmt='(a32,a2)') 'Save meshes:  ','no' 
    endif 
    if (lprintDebug) then
        write(6,fmt='(a32,a3)') 'Print debug:  ','yes'  
    else
        write(6,fmt='(a32,a2)') 'Print debug:  ','no' 
    endif     
    if (lDisplayRefinementStats) then
        write(6,fmt='(a32,a3)') 'Print adaptive:  ','yes'  
    else
        write(6,fmt='(a32,a2)') 'Print adaptive:  ','no' 
    endif  

    write(*,*) ' '    
    write(*,*) ' ----------------------------'
    write(*,*) '          Numerical settings:'
    write(*,*) ' ----------------------------'     
    write(*,fmt='(a32,f6.1,1a,f6.1,1a,i6,32a)') 'Wavenumbers:  ', loglower,',', logupper,',', nwave, '   ! log10 lower, upper, # '  
!    write(6,fmt='(a32,a)') ' (Co)Sine Transform Filters:  ', trim(FCTfilter)  
    write(6,fmt='(a32,a)') ' Linear Solver:  ', trim(linearSolver)      
        
 
    
    write(*,*) ' '    
    write(*,*) ' ----------------------------'
    write(*,*) '    Error Estimator settings:'
    write(*,*) ' ----------------------------'     
    write(ctemp,'(f8.2)') minRangeProtector
    write(6,fmt='(a32,a)') 'Minimum error range:  ', trim(adjustl(ctemp)) 

    write(ctemp,'(i5)') max_nsubrefine
    write(6,fmt='(a32,a)') 'Max # subrefinements:  ', trim(adjustl(ctemp))

    write(ctemp,'(f8.2)') pct_refine
    write(6,fmt='(a32,a)') 'Percent refine:  ', trim(adjustl(ctemp))

    write(ctemp,'(e11.4)') minArea
    write(6,fmt='(a32,a)') 'Minimum area:  ', trim(adjustl(ctemp)) 

    write(ctemp,'(e11.4)') ecutoff
    write(6,fmt='(a32,a)') 'E noise floor:  ', trim(adjustl(ctemp)) 

    write(ctemp,'(e11.4)') hcutoff
    write(6,fmt='(a32,a)') 'H noise floor:  ', trim(adjustl(ctemp))   


    write(*,*) ' '    
    write(*,*) ' ----------------------------'
    write(*,*) '      Finite dipole settings:'
    write(*,*) ' ----------------------------'   
    
    write(ctemp,'(i5)') nQuadTxCSEM
    write(6,fmt='(a32,a)') ' Transmitter quadrature order:  ', trim(adjustl(ctemp))
     
    write(ctemp,'(i5)') nQuadRxCSEM
    write(6,fmt='(a32,a)') 'CSEM receiver quadrature order: ', trim(adjustl(ctemp)) 
    
    write(ctemp,'(i5)') nQuadRxMT
    write(6,fmt='(a32,a)') ' MT receiver quadrature order:  ', trim(adjustl(ctemp))    
    
    write(*,*) ' ' 
    

     
    end subroutine display_MARE2DEM_Params
    
!==================================================================================================================================! 
!======================================================================================================================== checkModel
!==================================================================================================================================!       
    subroutine checkModel 
    
    use em2dkx_mod
    use mare2dem_global    
    use triangle_mesh
    
    implicit none
 
    logical :: lHasSlivers
         
    character(32) :: cend
    !
    write(*,*) '======== Inspecting Input Values ========================================='
    write(*,*) ' '     
    !
    ! First check to make sure that if model has negative conductivities (i.e., indices to params), then the rho file has
    ! been defined and read in:
    !
 
    if ( nSigParams == 0) then
 
        write(*,*) ' '
        write(*,*) ' Error, no resistivity parameters defined! '
        write(*,*) ' Please include a resistivity file!'
        write(*,*) ' Stopping!'
        stop        
  
    endif
    
    !
    ! Check to make sure the mesh attributes agree with the input number of fixed and free parameters:
    !
!    maxParam = 0
!    do i = 1,inputmodel%numberofregions
!        if ( inputmodel%regionlist(4*i-1) > 0 ) then   
!            maxParam = max(maxParam,nint(inputmodel%regionlist(4*i-1)))
!        endif
!    enddo
!    
!    if (maxParam /= nSigParams) then
!        write(*,*) ' '
!        write(*,*) ' Error, number of model parameters disagrees with the model file (.poly): '
!        write(*,*) ' Model file number of parameters:       ', maxParam
!        write(*,*) ' Resistivity file number of parameters: ', nSigParams
!        write(*,*) ' Stopping!'
!        stop    
!    endif
    !
    ! Check the model for slivers:
    !
 
    write(*,*) ' ... Checking the input model for slivers ...'

    call check_model(inputmodel,minQangle, lHasSlivers)
 
    if (lHasSlivers) then
    
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '        
        write(cend,fmt='(F4.0)') minQangle
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) '!!!     Warning Model has slivers... stopping        !!!'
        write(*,*) '!!! Please fix the mesh at the locations noted above !!!'
        write(*,*) '!!! so that segment intersections are > ',trim(adjustl(cend)),' degrees  !!!'  
        write(*,*) '!!!    >>> Better luck next time, bucko <<<          !!!' ! creds to CJW
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '                
        stop
    else
        write(*,*) ' '
        write(*,*) ' A-okay buddy:'
        write(*,*) ' '
        write(*,*) ' Your model has been checked for slivers and  '
        write(*,*) ' none were found.  '
        write(*,*) ' '   
    endif

    end subroutine checkModel       
         
!==================================================================================================================================! 
!==================================================================================================================== getWavenumbers
!==================================================================================================================================!       
    subroutine getWavenumbers
 
    use mare2dem_global
    use mare2dem_input_data_params ! for nFreqCSEM
    
    implicit none
 
    integer :: i 
    real(8) ::dlg
    
    if (nfreqCSEM < 1) return
    
    if (lprintSetup) then
        write(6,*) '========== Setting up wavenumbers ==========='
        write(*,*) ' '     
    endif

    allocate (wavenum(nwave))
    dlg = abs(logupper - loglower) / (nwave - 1)  
    
    do i = 1,nwave
        wavenum(i) = 10**(min(loglower,logupper) +(i-1)*dlg)
    enddo  
    
    if (lprintSetup) then  
        write(*,'(a32)') ' Wavenumbers:  '
        do i = 1,nwave            
             write(*,'(i32,1x,es12.3)') i,wavenum(i)
        enddo    
        write(*,*) ' '  
    endif    
    
    end subroutine getWavenumbers

!==================================================================================================================================! 
!===================================================================================================================== getRxTxGroups
!==================================================================================================================================!       
    subroutine getRxTxGroups 
!
! Sets up the Rx-Tx groups:    
!
    use mare2dem_global
    use mare2dem_input_data_params
    
    implicit none
 
 
    integer ::  nTxGroups, nRxGroupsCSEM, nRxGroupsMT, iGroup
    integer ::  itxG, ntxg, nt0, irxG, nrx0, nrxg
 

!
! Get total number of RxTxGroups:
!
    nTxGroups = 0
    if (nTxCSEM > 0 ) nTxGroups = ceiling( dble(nTxCSEM) / dble(nTxPerGroup) )   
 
    nRxGroupsCSEM     = ceiling( dble(nRxCSEM) / dble(nRxPerGroupCSEM) )   
    nRxGroupsMT       = ceiling( dble(nRxMT)   / dble(nRxPerGroupMT) )   
 
    nRxTxGroups       = nRxGroupsCSEM*nTxGroups + nRxGroupsMT
 
!
! Allocate storage for RxTxGroups
!   
    allocate ( RxTxGroups(nRxTxGroups) )
 
!
! CSEM RxTxGroups:
!
    iGroup = 0  
    do itxG = 1,nTxGroups
    
        nt0     = (itxG-1)*nTxPerGroup              ! last input transmitter number from previous group
        ntxg    = min(nTxPerGroup, nTxCSEM - nt0)    ! number of transmitters in this group
 
        do irxG = 1,nRxGroupsCSEM
    
            nrx0    =  (irxG-1)*nRxPerGroupCSEM                 ! last input receiver number from previous group or 0 if first group
            nrxg    =  min(nRxPerGroupCSEM, nRxCSEM - nrx0 )  ! number of receivers in this group
            
            !
            ! Insert arrays into RxTxGroups(iGroup):
            !
            
            iGroup = iGroup + 1
            
            !
            ! Tx indices:
            !
            RxTxGroups(iGroup)%lMT  = .false.
            
            RxTxGroups(iGroup)%nTx  = ntxg
            
            allocate( RxTxGroups(iGroup)%iTx(ntxg) )
 
            RxTxGroups(iGroup)%iTx  = [1:ntxg] + nt0    
             
            
            !
            ! Rx indices:
            !            
            RxTxGroups(iGroup)%nRx  = nrxg
            
            allocate( RxTxGroups(iGroup)%iRx(nrxg) )
            
            RxTxGroups(iGroup)%iRx  = [1:nrxg] + nrx0    
 
        enddo
    enddo
    
!
! MT RxTxGroups:
!
    if (nFreqMT > 0 ) then
    
        do irxG = 1,nRxGroupsMT
            
            nrx0  =  (irxG-1)*nRxPerGroupMT             ! last input receiver number from previous group or 0 if first group
            nrxg =  min(nRxPerGroupMT, nRxMT - nrx0 )   ! number of receivers in this group    
        
            iGroup = iGroup + 1
    
            !
            ! Tx stuff:
            !                           
            RxTxGroups(iGroup)%nTx      = 1         
                
            allocate ( RxTxGroups(iGroup)%iTx(1)  )
            RxTxGroups(iGroup)%iTx(1)   = 0   
            RxTxGroups(iGroup)%lMT      = .true.
         
            !
            ! Rx stuff:
            !
            RxTxGroups(iGroup)%nRx      = nrxg
            allocate( RxTxGroups(iGroup)%iRx(nrxg) )
                
            RxTxGroups(iGroup)%iRx      =  [1:nrxg] + nrx0  
                      
                        
        enddo
    endif    
        
    end subroutine getRxTxGroups
    
!==================================================================================================================================! 
!======================================================================================================================= getFqGroups
!==================================================================================================================================!        
    subroutine getFqGroups 
!
! Sets up the Fq groups:    
!
    use mare2dem_global
    use mare2dem_input_data_params
       
    implicit none
 
    integer :: nFreqGroupsCSEM          ! Total number of frequency groups    
    integer :: nFreqGroupsMT            ! Total number of frequency groups    
    integer :: iGroup, i 
    integer :: ifreqG, nf0, nfg
 
    
!
! Get total number of FqGroups:
!    
 
    nFreqGroupsCSEM = 0
    if (nTxCSEM > 0 )  nFreqGroupsCSEM = ceiling( dble(nfreqCSEM)  / dble(nFreqPerGroupCSEM)   )    
 
    nFreqGroupsMT = 0         
    if (nfreqMT > 0 ) nFreqGroupsMT = ceiling( dble(nfreqMT)    / dble(nFreqPerGroupMT) ) 
 
    nFqGroups = nFreqGroupsCSEM + nFreqGroupsMT
  
    allocate( FqGroups(nFqGroups) ) 
    
 !
 ! CSEM:
 !   
    iGroup = 0 
    
    do ifreqG = 1,nFreqGroupsCSEM       
            
        nf0 =  (ifreqG-1)*nFreqPerGroupCSEM          ! last input frequency number from previous group
        nfg = min(nFreqPerGroupCSEM, nfreqCSEM - nf0)  ! number of freqs in this group
     
        iGroup = iGroup + 1
         
        FqGroups(iGroup)%nFq = nfg
        
        allocate ( FqGroups(iGroup)%iFq(nfg) )
        
        do i = 1,nfg
            
            FqGroups(iGroup)%iFq(i) = nf0 + i    
            
            if ( mod(i+ floor(dble(nfg)/dble(2) ),nfg) == 0 ) FqGroups(iGroup)%iFqRefine = i
              
        enddo
                        
        FqGroups(iGroup)%lMT = .false.       
                            
    enddo
    
 !
 ! MT:
 !   
    do ifreqG = 1,nFreqGroupsMT   
    
        nf0 = (ifreqG-1)*nFreqPerGroupMT          ! last input frequency number from previous group
        nfg = min(nFreqPerGroupMT, nfreqMT - nf0) 

        iGroup = iGroup + 1

        FqGroups(iGroup)%nFq = nfg
                
        allocate ( FqGroups(iGroup)%iFq(nfg) )
 
        do i = 1,nfg
        
            FqGroups(iGroup)%iFq(i) = nf0 + i    
            if ( mod(i+ floor(dble(nfg)/dble(2) ),nfg) == 0 ) FqGroups(iGroup)%iFqRefine = i
          
        enddo     
        
        FqGroups(iGroup)%lMT = .true.
             
    enddo
     
    end subroutine getFqGroups     
    
!==================================================================================================================================! 
!=============================================================================================================== getRefinementGroups
!==================================================================================================================================!      
    subroutine getRefinementGroups
!
! Sets up the refinement groups by looking at RxTx and Fq combos that
! have input data. Also creates the data mask arrays. 
!
    use mare2dem_global
    use mare2dem_input_data_params
    
    implicit none  
    
    integer :: iRxTxG, iFqG, rx0,rx1,tx0,tx1,fq0,fq1, iPass
    integer :: nFq, nRx, nTx
    
 !
 ! Use a first pass to count how many groups to allocate, then we use a second pass to allocate and insert:
 !   
 
    do iPass = 1,2
    
        if (iPass == 2) then
            allocate( refinementGroups(nRefinementGroups) )
        endif    
        
        nRefinementGroups = 0 
        
        do iRxTxG = 1,nRxTxGroups
        
            nRx = RxTxGroups(iRxTxG)%nRx
            nTx = RxTxGroups(iRxTxG)%nTx
            
            rx0 = RxTxGroups(iRxTxG)%iRx(1)
            rx1 = RxTxGroups(iRxTxG)%iRx(nRx)
            tx0 = RxTxGroups(iRxTxG)%iTx(1)    ! note for MT, iTx(1) = 0
            tx1 = RxTxGroups(iRxTxG)%iTx(nTx)

            do iFqG = 1,nFqGroups
        
                nFq = FqGroups(iFqG)%nFq
                fq0 = FqGroups(iFqG)%iFq(1)
                fq1 = FqGroups(iFqG)%iFq(nFq)

                !
                ! MT refinement groups
                !
                if ( ( RxTxGroups(iRxTxG)%lMT ) .and. ( FqGroups(iFqG)%lMT ) ) then
                
                    if ( any( lDataMaskMT(rx0:rx1,fq0:fq1)) )  then
                        
                        nRefinementGroups = nRefinementGroups + 1  
                        
                        if (iPass == 2) then
                        
                            refinementGroups(nRefinementGroups)%iRxTxGroup = iRxTxG
                            refinementGroups(nRefinementGroups)%iFqGroup   = iFqG
                            refinementGroups(nRefinementGroups)%lMT        = .true.
      
                        endif    
                
                    endif  
            
                !
                ! CSEM refinement groups
                !
                elseif ( ( .not.RxTxGroups(iRxTxG)%lMT ) .and. ( .not.FqGroups(iFqG)%lMT ) ) then
               
                    if ( any( lDataMaskCSEM(rx0:rx1,fq0:fq1,tx0:tx1)) )  then

                        nRefinementGroups = nRefinementGroups + 1  
                        
                        if (iPass == 2) then
                        
                            refinementGroups(nRefinementGroups)%iRxTxGroup = iRxTxG
                            refinementGroups(nRefinementGroups)%iFqGroup   = iFqG
                            refinementGroups(nRefinementGroups)%lMT        = .false.
            
                        endif    
                     
                    endif
                
                endif
                
            enddo ! iFq
                       
        enddo ! iRxTxG
 
    enddo ! iPass   

!
! Point to first group:
!
    iPtr_refGroups = 1
    
    end subroutine getRefinementGroups
    
!==================================================================================================================================! 
!================================================================================================================ printDecomposition
!==================================================================================================================================!       
    subroutine printDecomposition

    use mare2dem_global
    use mare2dem_input_data_params
    
    implicit none  

    integer :: nTxGroups, nRxGroupsCSEM, nRxGroupsMT
    integer :: nFreqGroupsCSEM              ! Total number of frequency groups    
    integer :: nFreqGroupsMT            ! Total number of frequency groups    

    integer :: ict,i,j,k
    
    character (24) :: ctemp, str1, str2
    
         
    nTxGroups = 0
    if (nTxCSEM > 0 ) then
        nTxGroups     = ceiling( dble(nTxCSEM)     / dble(nTxPerGroup)     )   
    endif 
    
    
    nRxGroupsCSEM     = ceiling( dble(nRxCSEM)     / dble(nRxPerGroupCSEM)     )   
    nRxGroupsMT       = ceiling( dble(nRxMT)     / dble(nRxPerGroupMT)     )   

    nRxTxGroups       = nRxGroupsCSEM*nTxGroups + nRxGroupsMT

    nFreqGroupsCSEM = 0
    if (nTxCSEM > 0 )  nFreqGroupsCSEM   = ceiling( dble(nfreqCSEM)  / dble(nFreqPerGroupCSEM)   )    
    
    nFreqGroupsMT = 0         
    if (nfreqMT > 0 )  nFreqGroupsMT = ceiling( dble(nfreqMT)    / dble(nFreqPerGroupMT) ) 
    
    nFqGroups = nFreqGroupsCSEM + nFreqGroupsMT
 
    if (lprintDecomposition) then
        write(*,*) ' '    
        write(6,*) '============== Parallel decomposition settings: =========================='
        write(*,*) ' '    
        if (nTxCSEM > 0 ) write(ctemp,'(i5)') nTxPerGroup
        if (nTxCSEM > 0 ) write(6,fmt='(a32,a)') 'Transmitters per group:  ', trim(adjustl(ctemp))
        
        write(ctemp,'(i5)') nRxPerGroupMT
        if (nfreqMT > 0 ) write(6,fmt='(a32,a)') 'MT Receivers per group:  ', trim(adjustl(ctemp))

        write(ctemp,'(i5)') nRxPerGroupCSEM
        if (nTxCSEM > 0 )  write(6,fmt='(a32,a)') 'CSEM Receivers per group:  ', trim(adjustl(ctemp))
                
        if (nTxCSEM > 0 ) write(ctemp,'(i5)') nFreqPerGroupCSEM
        if (nTxCSEM > 0 ) write(6,fmt='(a32,a)') 'CSEM frequencies per group:  ', trim(adjustl(ctemp))
        
        
        if (nfreqMT > 0 )  write(ctemp,'(i5)') nFreqPerGroupMT  
        if (nfreqMT > 0 ) write(6,fmt='(a32,a)') 'MT frequencies per group:  ', trim(adjustl(ctemp))
    
        write(*,*) ' '
        write(*,*) ' Maximum possible groups if all Tx-Rx-Frequency '
        write(*,*) ' combinations are present in the input data: '
        write(*,*) ' '
        if (nTxCSEM > 0 ) write(ctemp,'(i5)') nTxGroups
        if (nTxCSEM > 0 ) write(*,fmt='(a32,a)') '# transmitter groups:  ', trim(adjustl(ctemp))

        
        write(ctemp,'(i5)') nRxGroupsCSEM
        write(*,fmt='(a32,a)') '# CSEM receiver groups:  ', trim(adjustl(ctemp))
 
        write(ctemp,'(i5)') nRxGroupsMT
        write(*,fmt='(a32,a)') '# MT receiver groups:  ', trim(adjustl(ctemp))
               
        if (nTxCSEM > 0 ) write(ctemp,'(i5)') nFreqGroupsCSEM
        if (nTxCSEM > 0 ) write(*,fmt='(a32,a)') '# CSEM frequency groups:  ', trim(adjustl(ctemp))
        
        if (nfreqMT > 0 ) write(ctemp,'(i5)') nFreqGroupsMT
        if (nfreqMT > 0 ) write(*,fmt='(a32,a)') '# MT frequency groups:  ', trim(adjustl(ctemp))
        
        write(*,fmt='(a32,a)') ' -','--'
        write(ctemp,'(i7)') nTxGroups*nRxGroupsCSEM*nFreqGroupsCSEM + nRxGroupsMT*nFreqGroupsMT
        write(*,fmt='(a32,a)') ' # refinement groups possible:  ', trim(adjustl(ctemp))    
        write(*,*) ' '
        
        
        write(*,*) ' After inspecting the input data, I have found: '
        write(*,*) ' ' 
        write(ctemp,'(i7)') nRefinementGroups
        write(*,fmt='(a32,a)') ' # of refinement groups:  ', trim(adjustl(ctemp))      
        write(*,*) ' '  
        write(*,*) ' '  
        
        if (nTxCSEM > 0 ) then
    
            ! compute fill-in of CSEM data:
            ict = 0
            do i=1,nRxCSEM
                do j = 1,nFreqCSEM
                    do k = 1,nTxCSEM
                        if ( lDataMaskCSEM(i,j,k)) ict = ict + 1
                    enddo
                enddo
            enddo
 
            write(str1,*) ict
            write(str2,*) nRxCSEM*nFreqCSEM*nTxCSEM
            write(*,'(a32,a,1x,a6,1x,a,1x,a1,f6.1,a4)') ' CSEM Data fill in: ',trim(adjustl(str1)), 'out of', trim(adjustl(str2)),&
            &'(', real(ict)/real(nRxCSEM*nFreqCSEM*nTxCSEM)*100.,' % )'         
     
        endif
        
        if ( nFreqMT > 0 ) then
            
            ! compute fill-in of CSEM data:
            ict = 0
            do i=1,nRxMT
                do j = 1,nFreqMT
                    if ( lDataMaskMT(i,j)) ict = ict + 1
                enddo
            enddo
 
            write(str1,*) ict
            write(str2,*) nRxMT*nFreqMT
            write(*,'(a32,a,1x,a6,1x,a,1x,a1,f6.1,a4)') ' MT Data fill in: ', trim(adjustl(str1)), 'out of', trim(adjustl(str2)), &
                 & '(', real(ict)/real(nRxMT*nFreqMT)*100.,' % )'
                            
        endif   
        write(*,*) '   '
    endif

    end subroutine printDecomposition
    
!==================================================================================================================================! 
!=================================================================================================================== getStartingMesh
!==================================================================================================================================!     
    subroutine getStartingMesh
!
 
    use triangle_mesh
    use mare2dem_global
    use em2dkx_mod
 
    implicit none
 
    character (24) :: cend,tricommand 
    
!
! Create the starting mesh:
!
    call copy_trimesh(inputmodel,inputmesh)

    if (inputmesh%nele == 0 ) then 
        if (lprintSetup) then
            write(*,*)  ' '
            write(*,*) ' Generating mesh from input polygon model by calling Triangle:'
            write(*,*) ' '
        endif
        
        write(cend,fmt='(F4.0)') minQangle
        tricommand = 'q'//trim(adjustl(cend))//'QpanjA'//CHAR(0)
        call call_triangle(tricommand,inputmesh)  
    else
        write(*,*) 'Error: input mesh already has been triangulated, stopping'
        write(*,*) ' '
        stop
    endif
    if (lprintSetup) then
        write(*,'(a38,1x,i6,a2,i6)') ' Starting mesh: # nodes, # elements: ',inputmesh%nnod,', ',inputmesh%nele
   
        write(*,'(a38,1x,g11.4,1x,g11.4)') ' Min/max resistivities: ',1./maxval(sigParams),1./minval(sigParams)
        write(*,*) ' '  
    endif
    
    !
    ! Check to make sure regions don't have 0 sig (usually regions that you've forgotten to assigned during mesh design)
    !
    if ( any( inputmesh%attr == 0) )then
        write(*,*) ''
        write(*,*) ' !!! Error reading mesh !!!'
        Write(*,*) ' Elements have 0 attributes!'
        ! write(*,*) ' Please fix the finite element mesh to have nonzero integers for fixed ( < 0 ) or free ( > 0 ) parameters'
        ! kwk debug: put in loop here that shows where the problem element(s) are located
        write(*,*) ' Stopping. '
        stop    
    endif   
 
    if (lprintSetup) then
        write(*,*) '  '
        write(*,*) ' Set up done, starting 2.5D EM computations...'
        write(*,*) '  '
    endif
    
    
    end subroutine getStartingMesh



!==================================================================================================================================! 
!=================================================================================================================== get_time_offset
!==================================================================================================================================!  
    subroutine  get_time_offset(timein,timeout)
!    
! timein is the clock start time or 0.
! timeout is the time in seconds since the input time
!
! Version 2.0  February 25, 2008  Now uses date_and_time Fortran intrinsic
!
! Written by
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
! ---------------------------------------------------------------------
   implicit none
   
   integer, dimension(8) :: values
   integer               :: i,j,k,mjd
   
   real(8) :: timein, timeout, fracday
    

 !
 ! New standard Fortran90 Time function:
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
                         
    end  subroutine  get_time_offset
    
!==================================================================================================================================! 
!============================================================================================================== applyPhaseConvention
!==================================================================================================================================!     
    subroutine applyPhaseConvention

    use mare2dem_input_data_params        
    use mare2dem_output
    use mare2dem_global   

    implicit none 

    integer :: iTx, iRx, ifreq    

    select case (trim(phaseConvention))

    case ('lag')
            ! do nothing since em2dkx uses lag by default
    case ('lead')
        
        ex = conjg(ex)
        ey = conjg(ey)
        ez = conjg(ez)
        hx = conjg(hx)
        hy = conjg(hy)
        hz = conjg(hz)
    
        if  (linversion) then
            
            do iTx = 1,nTxCSEM         

                do ifreq = 1,nfreqCSEM
       
                    do iRx=1,nRxCSEM
                         if (lCompCSEM(1,iRx)) dex(iRx,iTx,ifreq)%dsig = conjg( dex(iRx,iTx,ifreq)%dsig )
                         if (lCompCSEM(2,iRx)) dey(iRx,iTx,ifreq)%dsig = conjg( dey(iRx,iTx,ifreq)%dsig )
                         if (lCompCSEM(3,iRx)) dez(iRx,iTx,ifreq)%dsig = conjg( dez(iRx,iTx,ifreq)%dsig )
                         if (lCompCSEM(4,iRx)) dhx(iRx,iTx,ifreq)%dsig = conjg( dhx(iRx,iTx,ifreq)%dsig )
                         if (lCompCSEM(5,iRx)) dhy(iRx,iTx,ifreq)%dsig = conjg( dhy(iRx,iTx,ifreq)%dsig )
                         if (lCompCSEM(6,iRx)) dhz(iRx,iTx,ifreq)%dsig = conjg( dhz(iRx,iTx,ifreq)%dsig )
                    enddo !   
         
                enddo  ! ifreq
        
            enddo ! iTx
    
        endif 
        
    end select   
                
    end subroutine applyPhaseConvention

           
!==================================================================================================================================! 
!============================================================================================================================ indexx
!==================================================================================================================================!
! From N.R.
       
    SUBROUTINE indexx(n,arr,indx) 
! Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that 
! arr(indx(j)) is in ascending order for j=1,2,...,N. 
! The input quantities n and arr are not changed.
      implicit none
      INTEGER n,indx(n),M,NSTACK
      REAL(8) arr(n)
      PARAMETER (M=7,NSTACK=1000)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL(8) a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

