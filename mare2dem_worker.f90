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
!=================================================================================================================== mare2dem_worker
!==================================================================================================================================! 
    module mare2dem_worker   
    
!
! Variables stored here are *ONLY* used at the local worker cpu level.
!
    use em2dkx_mod                      ! for settings
    use mare2dem_global                 ! for settings 
    use mare2dem_input_data_params      ! for data input 
    use mare2dem_output                 ! for data output 
    use spline_integ_kx
    
    implicit none   

    private

!
! Local transmitter and receiver arrays for handling finite dipoles
!
    integer                             :: ntrans = 0              ! # of transmitter locations and Tx parameters
    real(8), dimension(:), allocatable  :: azimuthtrans,diptrans,xtrans,ytrans,ztrans 
    character(32), allocatable          :: sourceType(:)  !'edipole' or 'bdipole',
        
    integer                             :: nsiteCSEM = 0
    real(8), dimension(:), allocatable  :: xsiteCSEM, ysiteCSEM, zsiteCSEM
    real(8), dimension(:), allocatable  :: thetaCSEM, alphaCSEM, betaCSEM  ! for Rx rotations for point and finite dipoles    

    integer                             :: nsiteMT = 0
    real(8), dimension(:), allocatable  :: xsiteMT, ysiteMT, zsiteMT
    real(8), dimension(:), allocatable  :: thetaMT, alphaMT, betaMT  ! for Rx rotations for point and finite dipoles  
            
!  
! y-z modes are combined since field symmetry is the same:
!    
    integer, dimension(7)               :: txcount = (/1,1,2,1,2,1,2/) ! number of x,y,z dipole computations need for each mode  
    integer, dimension(:), allocatable  :: TxMode  ! for each nTx, gives 1-7 index to txcount array
        !
        ! Here are 7 possible cases, using abs(azimuth) and abs(dip)
        !
        !       Azimuth     Dip         kxmode(s)   Dipoles     #Modes
        ! 1     0  | 180    0           1           x           1
        ! 2     90 | 270    0           2           y           1
        ! 3     any         0           1,2         x,y         2
        ! 4     any         90          3           z           1
        ! 5     0  | 180    >0          1,3         x,z         2
        ! 6     90 | 270    >0          4           yz          1
        ! 7     any         >0          1,4         x,yz        2
        
 
!
! Quadrature weights used for finite dipole integrations:
!
    real(8), dimension(:), allocatable  :: quad_weights_Tx,     quad_xint_Tx     ! integration weights and local position 
    real(8), dimension(:), allocatable  :: quad_weights_RxCSEM, quad_xint_RxCSEM ! integration weights and local position
    real(8), dimension(:), allocatable  :: quad_weights_RxMT,   quad_xint_RxMT   ! integration weights and local position
 
    integer, dimension(:), allocatable  :: indexTrans,indexSite ! ntrans arrays storing input iTx and iRx
    integer, dimension(:), allocatable  :: indexTrans_em2dkx    ! nTx array storing input iTx index
     
    type  :: KxFqG   ! Derived type to store subgrid of Kx and Fq combinations:
        integer                             :: nKxFq    ! number of non-refinement wavenumber,frequency pairs in this group   
        integer, dimension(:), allocatable  :: iFq      ! index of non-refinement frequencies 
        integer, dimension(:), allocatable  :: iKx      ! index of wavenumbers for this group
        integer                             :: iFqRefine  ! local index to freq for refinement frequency for this group         
        integer                             :: iKxRefine  ! local index to wavenumber for refinement for this group  
    end type
 
    integer                                 :: nKxFqGroups
    type(KxFqG), dimension(:), allocatable  :: KxFqGroups      
 
    type  :: ScrG   ! Derived type to store Rx,Tx combinations written to scratch files (for sensitivity derivatives):
        integer                             :: nRxTx    ! number of Rx,Tx pairs in this group
        integer, dimension(:), allocatable  :: iRx      ! index of receivers  (for em2dkx arrays)
        integer, dimension(:), allocatable  :: iTx      ! index of transmitters (for em2dkx arrays)
    end type
 
    integer                                 :: nScratchGroups  
    type(ScrG), dimension(:), allocatable   :: ScratchGroups        
    
    logical     :: lRefine

    integer, parameter ::  nRxTxPerIFT = 10
    
!
! Local storage for kx domain arrays:
!    
    complex(4), dimension(:,:,:,:), allocatable     ::  ex_kx_s,ey_kx_s,ez_kx_s,hx_kx_s,hy_kx_s,hz_kx_s    ! nRx x nTx x nkx x nFreq      
    
    complex(4), dimension(:,:), allocatable         ::  readBuffer   

!
! Temporary storage for wavenumber domain arrays:
!    
    type(dsigma),dimension(:,:,:), allocatable      :: dexds,deyds,dezds,dhxds,dhyds,dhzds! (iRx,iTx,iFq), defined in mare2dem_common
    
    logical, dimension(:,:), allocatable :: lCSEMworker   ! 6 by nRxCSEM for logical flags on which components needed to get output
       
    public :: worker_EM2D
    
    contains
!==================================================================================================================================! 
!================================================================================================================= worker_deallocate
!==================================================================================================================================!       

    subroutine worker_deallocate
    
    call deallocate_em2dkx  ! deallocate i/o arrays for em2dkx
 
    if ( allocated(quad_weights_Tx) )       deallocate( quad_weights_Tx, quad_xint_Tx )       
    if ( allocated(quad_weights_RxCSEM) )   deallocate( quad_weights_RxCSEM, quad_xint_RxCSEM )       
    if ( allocated(quad_weights_RxMT) )     deallocate( quad_weights_RxMT, quad_xint_RxMT )    
 
    if ( allocated(indexTrans) )            deallocate( indexTrans )        
    if ( allocated(indexSite) )             deallocate( indexSite )  
    if ( allocated(indexTrans_em2dkx) )     deallocate( indexTrans_em2dkx )           
        
    if ( allocated(azimuthtrans) )          deallocate( azimuthtrans, diptrans, xtrans, ytrans, ztrans, sourcetype )  
    if ( allocated(xsiteCSEM) )             deallocate( xsiteCSEM, ysiteCSEM, zsiteCSEM )  
    if ( allocated(thetaCSEM) )             deallocate (thetaCSEM, alphaCSEM, betaCSEM)
    
    if ( allocated(xsiteMT) )               deallocate( xsiteMT, ysiteMT, zsiteMT )  
    if ( allocated(thetaMT) )               deallocate (thetaMT, alphaMT, betaMT)    

    if ( allocated(TxMode) )                deallocate( TxMode ) 
    if ( allocated(KxFqGroups) )            deallocate( KxFqGroups )     
    if ( allocated(ScratchGroups) )         deallocate( ScratchGroups )      
      
    if ( allocated(ex_kx_s) )               deallocate( ex_kx_s )    
    if ( allocated(ey_kx_s) )               deallocate( ey_kx_s )    
    if ( allocated(ez_kx_s) )               deallocate( ez_kx_s )  
    if ( allocated(hx_kx_s) )               deallocate( hx_kx_s )    
    if ( allocated(hy_kx_s) )               deallocate( hy_kx_s )    
    if ( allocated(hz_kx_s) )               deallocate( hz_kx_s )  
    
    if ( allocated(readBuffer) )            deallocate( readBuffer )                    
    if ( allocated(dexds) )                 deallocate( dexds )    
    if ( allocated(deyds) )                 deallocate( deyds )  
    if ( allocated(dezds) )                 deallocate( dezds )  
    if ( allocated(dhxds) )                 deallocate( dhxds )    
    if ( allocated(dhyds) )                 deallocate( dhyds )  
    if ( allocated(dhzds) )                 deallocate( dhzds )     
                           
    if ( allocated(lCSEMworker) )           deallocate( lCSEMworker )                                 
  
!
! Zero out a few local counters so they are ready for the next call:
!  
    ntrans    = 0 
    nsiteCSEM = 0
    nsiteMT   = 0
    
    end subroutine worker_deallocate    
            
!==================================================================================================================================! 
!======================================================================================================================= worker_EM2D
!==================================================================================================================================!       

    subroutine worker_EM2D()

!
! Set up some logical arrays about which components are needed:
!
    if (lPrintDebug) write(*,*) 'call setup_component_flags ', myID
    call setup_component_flags 
    
!
! Set up Rx and Tx arrays to deal with any finite Rx and Tx dipoles:
!
    if (lPrintDebug) write(*,*) 'call setup_finite_rx_tx ', myID
    call setup_finite_rx_tx   
 
!
! Set up the modal transmitters needed by em2dkx:
!
    if (lPrintDebug) write(*,*) 'call setup_em2dkx_rxtx ', myID
    call setup_em2dkx_rxtx ! this decomposes arbitrary Tx's into x and yz pointing transmitters
 
!
! Create the KxGroups structure:
!
    if (lPrintDebug) write(*,*) 'call setup_KxFqGroups ', myID
    call setup_KxFqGroups  
    
!
! Create the ScratchGroups structure:
!
    if (lPrintDebug) write(*,*) 'call setup_ScratchGroups ', myID
     call setup_ScratchGroups  
        
! 
! Compute the 2.5D wavenumber domain responses for all wavenumbers and frequencies:
!
    if (lPrintDebug) write(*,*) 'call computeAllKxFq ', myID
    call computeAllKxFq
    
!
! Inverse Fourier transform from the (kx,y,z) domain to  the (x,y,z) domain
!    
    if (lPrintDebug) write(*,*) 'call inverseFourierTransform ', myID
    call inverseFourierTransform

    if (linversion) then
        if (lPrintDebug) write(*,*) 'call inverseFourierTransform_Derivs ', myID
        call inverseFourierTransform_Derivs
    endif

!
! Rotate to requested receiver orientations:
!
    if (lPrintDebug) write(*,*) 'call rotateToRxOrientation ', myID
    call rotateToRxOrientation

!
! Allocate solution arrays: Note that these are for the point dipoles, which are latter integrated for finite length dipoles
!
    if (lPrintDebug) write(*,*) 'call setup_solutionArrays ', myID
    call setup_solutionArrays   
           
!
! Integrate the finite length Tx and Rx antennas:
!    
    if (lPrintDebug) write(*,*) 'call integrate_dipoles ', myID
    call integrate_dipoles
    

!
! Deallocate:
!
    call worker_deallocate
    
    end subroutine worker_EM2D

!==================================================================================================================================! 
!============================================================================================================= setup_component_flags
!==================================================================================================================================!       
    subroutine setup_component_flags
    
    integer :: iRx
    
    if (.not.lMT) then
    
        allocate( lCSEMworker(6,nRxCSEM) )
        lCSEMworker = .false.
        
        do iRx = 1,nRxCSEM
 
    
            if ( any(lCompCSEM(1:3,iRx)) ) then  ! any electrics requested
        
                if ( lCompCSEM(1,iRx) )          lCSEMworker(1,iRx) = .true.
                if ( lCompCSEM(2,iRx) )          lCSEMworker(2,iRx) = .true.
                if ( lCompCSEM(3,iRx) )          lCSEMworker(3,iRx) = .true.
            
                if ( thetaRxCSEM(iRx) /= 0 )     lCSEMworker(1,iRx) = .true.
                if ( thetaRxCSEM(iRx) /= 0 )     lCSEMworker(2,iRx) = .true.
    
                if ( AlphaRxCSEM(iRx) /= 0 )     lCSEMworker(1,iRx) = .true.
                if ( AlphaRxCSEM(iRx) /= 0 )     lCSEMworker(3,iRx) = .true.
                
                if ( BetaRxCSEM(iRx)  /= 0 )     lCSEMworker(2,iRx) = .true.
                if ( BetaRxCSEM(iRx)  /= 0 )     lCSEMworker(3,iRx) = .true.
                
            endif
        
            if ( any(lCompCSEM(4:6,iRx)) ) then  ! any electrics requested
        
                if ( lCompCSEM(4,iRx) )          lCSEMworker(4,iRx) = .true.
                if ( lCompCSEM(5,iRx) )          lCSEMworker(5,iRx) = .true.
                if ( lCompCSEM(6,iRx) )          lCSEMworker(6,iRx) = .true.
            
                if ( thetaRxCSEM(iRx) /= 0 )     lCSEMworker(4,iRx) = .true.
                if ( thetaRxCSEM(iRx) /= 0 )     lCSEMworker(5,iRx) = .true.
    
                if ( AlphaRxCSEM(iRx) /= 0 )     lCSEMworker(4,iRx) = .true.
                if ( AlphaRxCSEM(iRx) /= 0 )     lCSEMworker(6,iRx) = .true.
                
                if ( BetaRxCSEM(iRx)  /= 0 )     lCSEMworker(5,iRx) = .true.
                if ( BetaRxCSEM(iRx)  /= 0 )     lCSEMworker(6,iRx) = .true.
                
            endif
                
        enddo
        
    endif
    
    end subroutine setup_component_flags       
!==================================================================================================================================! 
!================================================================================================================ setup_finite_rx_tx
!==================================================================================================================================!       
    subroutine setup_finite_rx_tx

!
! Allocates arrays for dealing with finite length Tx's and Rx's
!
  
    integer                     :: i,j,  icnt,  icomp
    real(8)                     :: xx, azm, dip
    real(8)                     :: theta, alpha, beta           ! site rotation angles
    real(8) , dimension(3,3)    :: RotR 
 
 
if (lMT) then
        
!-----------------
! MT Rx dipoles:
!-----------------
 
    icnt = 0    
    do i = 1,nRxMT
 
        if ( lengthRxMT(i) > 0 ) then ! finite dipole    
        
            do j = 1,3
                if (lCompMT(j,i)) then ! this receiver & electric dipole component has data, add to counter:
                    icnt = icnt + nQuadRxMT
                endif
            enddo        
        
        else
        
            icnt = icnt + 1 ! single point dipole
        
        endif
    enddo
    nsiteMT = icnt

    allocate ( xsiteMT(nsiteMT), ysiteMT(nsiteMT), zsiteMT(nsiteMT) ) 
    allocate ( thetaMT(nsiteMT), alphaMT(nsiteMT), betaMT(nsiteMT) ) 
    allocate ( indexSite(nsiteMT) )   
    
    !
    ! Get Rx quadrature points:
    !
    allocate(quad_weights_RxMT(nQuadRxMT), quad_xint_RxMT(nQuadRxMT) )   
    
    call legendre_compute_dr(nQuadRxMT,quad_xint_RxMT,quad_weights_RxMT) 

    !
    ! Set Rx points:
    !
    icnt = 0
    do i = 1,nRxMT
        
        if ( lengthRxMT(i) > 0 ) then ! finite dipole              
            
            theta = ThetaRxMT(i)
            alpha = AlphaRxMT(i)
            beta  = BetaRxMT(i)        
 
            call getRotationMatrix(theta,alpha,beta, RotR)
             
            ! Row 1 corresponds to mapping of x dipole into rotated space, row 2 is for y dipole and row z is for z dipole.           
          
            ! So loop over x,y,z dipole components and create virtual receiver point dipoles:

            do icomp = 1,3
            
                if (lCompMT(icomp,i)) then
                
                    do j = 1,nQuadRxMT

                        icnt = icnt + 1

                        xx = quad_xint_RxMT(j)*lengthRxMT(i)/2d0 ! xint goes from -1 to 1, so scale to position point dipole array from -l/2 to +l/2

                        indexSite(icnt) = i
                        
                        xsiteMT(icnt) = xRxMT(i) + xx*RotR(icomp,1)
                        ysiteMT(icnt) = yRxMT(i) + xx*RotR(icomp,2)
                        zsiteMT(icnt) = zRxMT(i) + xx*RotR(icomp,3)

                        thetaMT(icnt) = ThetaRxMT(i)
                        alphaMT(icnt) = AlphaRxMT(i)
                        betaMT(icnt)  = BetaRxMT(i)  

                        !if (lPrintDebug) write(*,'(a,1x,i6,i6,3(f9.2,1x))') 'MT Rx fd: ',i, icnt,  xsiteMT(icnt),ysiteMT(icnt),zsiteMT(icnt) 

                    enddo

                endif
                
            enddo
                           
        else  ! Point dipole

            icnt = icnt + 1

            indexSite(icnt) = i
            
            xsiteMT(icnt) = xRxMT(i)
            ysiteMT(icnt) = yRxMT(i)
            zsiteMT(icnt) = zRxMT(i)            
            thetaMT(icnt) = ThetaRxMT(i)
            alphaMT(icnt) = AlphaRxMT(i)
            betaMT(icnt)  = BetaRxMT(i)    
         
            !if (lPrintDebug) write(*,'(a,1x,i6,i6,3(f9.2,1x))') 'MT Rx pd: ',i,icnt,  xsiteMT(icnt),ysiteMT(icnt),zsiteMT(icnt)   
       
        endif            
           
    enddo 
 
 else 
!-----------------
! CSEM Tx dipoles:
!-----------------
 
    !
    ! Get total number of point dipoles needed for point and finite length dipole calculations:
    !
    icnt = 0
    do i = 1,nTxCSEM
        
        if ( any(lDataMaskCSEM(:,:,i)) ) then ! only add dipoles needed by data
        
            if (lengthTxCSEM(i) > 0) then   
                icnt = icnt + nQuadTxCSEM
            else
                icnt = icnt + 1
            endif
        
        endif
    enddo
    ntrans = icnt  

    !
    ! Allocate Tx arrays:
    !
    allocate ( azimuthtrans(ntrans),diptrans(ntrans),xtrans(ntrans),ytrans(ntrans),ztrans(ntrans),sourceType(ntrans) )     
    allocate ( indexTrans(ntrans) )   
    
    !
    ! Get Tx quadrature points:
    !
    allocate( quad_weights_Tx(nQuadTxCSEM), quad_xint_Tx(nQuadTxCSEM) )   
    
    call legendre_compute_dr( nQuadTxCSEM,quad_xint_Tx,quad_weights_Tx )

    !
    ! Set Tx points:
    !
    icnt = 0
    do i = 1,nTxCSEM
        
        if ( any(lDataMaskCSEM(:,:,i)) ) then ! only add dipoles needed by data
        
            if (lengthTxCSEM(i) > 0) then ! Finite length dipole using point dipoles at Gauss-Legendre points

                dip = dipTxCSEM(i) * deg2rad
                azm = azimuthTxCSEM(i) * deg2rad

                do j = 1,nQuadTxCSEM
                    icnt = icnt + 1
                    
                    indexTrans(icnt)   = i
                    azimuthtrans(icnt) = azimuthTxCSEM(i)
                    diptrans(icnt)     = dipTxCSEM(i) 
                    sourceType(icnt)   = cSourceType(i)

                    xx = quad_xint_Tx(j)*lengthTxCSEM(i)/2d0 ! xint goes from -1 to 1, so scale to position point dipole array from -l/2 to +l/2

                    xtrans(icnt)       = xTxCSEM(i) + (xx*cos(dip) * cos(azm))   
                    ytrans(icnt)       = yTxCSEM(i) + (xx*cos(dip) * sin(azm))
                    ztrans(icnt)       = zTxCSEM(i) + (xx*sin(dip))
                
                   !if (lPrintDebug) write(*,'(a,1x,i4,i4,3(f9.2,1x))') 'Tx fd: ',i,  icnt,   xtrans(icnt),ytrans(icnt),ztrans(icnt) 
                enddo

            else  ! Point dipole
        
                icnt = icnt + 1
                indexTrans(icnt)   = i
                azimuthtrans(icnt) = azimuthTxCSEM(i)
                diptrans(icnt)     = dipTxCSEM(i)
                xtrans(icnt)       = xTxCSEM(i)
                ytrans(icnt)       = yTxCSEM(i)
                ztrans (icnt)      = zTxCSEM(i)
                sourceType(icnt)   = cSourceType(i)
               ! if (lPrintDebug) write(*,'(a,1x,i4,i4,3(f9.2,1x))') 'Tx pd: ',i, icnt,  xtrans(icnt),ytrans(icnt),ztrans(icnt)  
            
            endif
        endif
    enddo 
               
!-----------------
! CSEM Rx dipoles:
!-----------------

    icnt = 0    
    do i = 1,nRxCSEM

        if ( any(lDataMaskCSEM(i,:,:)) ) then ! only add dipoles needed by data
        
            if ( lengthRxCSEM(i) > 0 ) then ! finite dipole    
                    
                do j = 1,3
                    if (lCompCSEM(j,i)) then ! this receiver & electric dipole component has data, add to counter:
                        icnt = icnt + nQuadRxCSEM 
                    endif
                enddo        
        
            else
        
                icnt = icnt + 1 ! single point dipole
        
            endif
            
        endif
    enddo
    
    nsiteCSEM = icnt

    allocate ( xsiteCSEM(nsiteCSEM), ysiteCSEM(nsiteCSEM), zsiteCSEM(nsiteCSEM) ) 
    allocate ( thetaCSEM(nsiteCSEM), alphaCSEM(nsiteCSEM), betaCSEM(nsiteCSEM) ) 
    allocate ( indexSite(nsiteCSEM) )   
    !
    ! Get Rx quadrature points:
    !
    allocate(quad_weights_RxCSEM(nQuadRxCSEM), quad_xint_RxCSEM(nQuadRxCSEM) )   
    
    call legendre_compute_dr(nQuadRxCSEM,quad_xint_RxCSEM,quad_weights_RxCSEM) 

    !
    ! Set Rx points:
    !
    icnt = 0
    do i = 1,nRxCSEM
    
        if ( any(lDataMaskCSEM(i,:,:)) ) then ! only add dipoles needed by data

            if ( lengthRxCSEM(i) > 0 ) then ! finite dipole              
            
                theta = ThetaRxCSEM(i)
                alpha = AlphaRxCSEM(i)
                beta  = BetaRxCSEM(i)        
 
                call getRotationMatrix(theta,alpha,beta, RotR)
             
                ! Row 1 corresponds to mapping of x dipole into rotated space, row 2 is for y dipole and row z is for z dipole.           
          
                ! So loop over x,y,z dipole components and create virtual receiver point dipoles:

                do icomp = 1,3
            
                    if (lCompCSEM(icomp,i)) then
                
                        do j = 1,nQuadRxCSEM

                            icnt = icnt + 1

                            xx = quad_xint_RxCSEM(j)*lengthRxCSEM(i)/2d0 ! xint goes from -1 to 1, so scale to position point dipole array from -l/2 to +l/2
                            
                            indexSite(icnt) = i
                            xsiteCSEM(icnt) = xRxCSEM(i) + xx*RotR(icomp,1)
                            ysiteCSEM(icnt) = yRxCSEM(i) + xx*RotR(icomp,2)
                            zsiteCSEM(icnt) = zRxCSEM(i) + xx*RotR(icomp,3)

                            thetaCSEM(icnt) = ThetaRxCSEM(i)
                            alphaCSEM(icnt) = AlphaRxCSEM(i)
                            betaCSEM(icnt)  = BetaRxCSEM(i)  

                        ! if (lPrintDebug) write(*,'(a,1x,i6,i6,3(f9.2,1x))') 'Rx fd: ',i,icnt,  xsiteCSEM(icnt),ysiteCSEM(icnt),zsiteCSEM(icnt)   

                        enddo

                    endif
                
                enddo
                           
            else  ! Point dipole

                icnt = icnt + 1
                
                indexSite(icnt) = i
                xsiteCSEM(icnt) = xRxCSEM(i)
                ysiteCSEM(icnt) = yRxCSEM(i)
                zsiteCSEM(icnt) = zRxCSEM(i)            
                thetaCSEM(icnt) = ThetaRxCSEM(i)
                alphaCSEM(icnt) = AlphaRxCSEM(i)
                betaCSEM(icnt)  = BetaRxCSEM(i)    
         
               ! if (lPrintDebug) write(*,'(a,1x,i6,i6,3(f9.2,1x))') 'Rx pd: ',i,icnt,  xsiteCSEM(icnt),ysiteCSEM(icnt),zsiteCSEM(icnt)   
       
            endif            
              
        endif 
              
    enddo 

endif ! lMT    
    

    end subroutine setup_finite_rx_tx
    
!==================================================================================================================================! 
!================================================================================================================= getRotationMatrix
!==================================================================================================================================!       
    subroutine getRotationMatrix(theta,alpha,beta,RotR)
!
! Outputs the rotation matrix is for left multiplication: F_rotated = RotR * F_xyz
!    
    real(8), intent(in)                     :: theta, alpha, beta   ! site rotation angles
    real(8), dimension(3,3), intent(out)    :: RotR                 ! 3x3 rotation matrices
 
    real(8)                                 :: cct, sst, cca, ssa, ccb, ssb ! cos and sin values for rotation angles
    real(8) , dimension(3,3)                :: Rotx,Roty,Rotz          ! 3x3 rotation matrices
   
    
    cct = cos(deg2rad* theta)
    sst = sin(deg2rad* theta)
    cca = cos(deg2rad* alpha)
    ssa = sin(deg2rad* alpha)
    ccb = cos(deg2rad* beta)
    ssb = sin(deg2rad* beta)                

    Rotx = 0d0
    Rotx(1,1) = 1.d0
    Rotx(2,2) =  ccb
    Rotx(2,3) =  ssb
    Rotx(3,2) = -ssb
    Rotx(3,3) =  ccb

    Roty = 0d0
    Roty(1,1) =  cca  
    Roty(1,3) =  ssa
    Roty(2,2) =  1.d0 
    Roty(3,1) = -ssa
    Roty(3,3) =  cca


    Rotz = 0d0
    Rotz(1,1) =  cct
    Rotz(1,2) =  sst
    Rotz(2,1) = -sst
    Rotz(2,2) =  cct
    Rotz(3,3) =  1.d0    

    RotR = matmul(Roty,Rotz)  
    RotR = matmul(Rotx,RotR)       
    
    
    end subroutine getRotationMatrix 
    
!==================================================================================================================================! 
!============================================================================================================== setup_solutionArrays
!==================================================================================================================================!       
    subroutine setup_solutionArrays
 

!
! The arrays allocated here are from module mare2dem_output
!     
     
!
! Allocate arrays for CSEM solutions at receivers: Note that these are the rotated spatial domain fields & are passed to manager 
!
    allocate( ex(nRxCSEM,nTxCSEM,nfreqCSEM),  &
              ey(nRxCSEM,nTxCSEM,nfreqCSEM),  &
              ez(nRxCSEM,nTxCSEM,nfreqCSEM),  &
              hx(nRxCSEM,nTxCSEM,nfreqCSEM),  &
              hy(nRxCSEM,nTxCSEM,nfreqCSEM),  &
              hz(nRxCSEM,nTxCSEM,nfreqCSEM)  )
 
    !
    ! If requested, allocate space for the sensitivity derivatives:
    !        
    if (linversion) then 
       
       allocate( dex(nRxCSEM,nTxCSEM,nfreqCSEM) )  
       allocate( dey(nRxCSEM,nTxCSEM,nfreqCSEM) ) 
       allocate( dez(nRxCSEM,nTxCSEM,nfreqCSEM) ) 
       allocate( dhx(nRxCSEM,nTxCSEM,nfreqCSEM) ) 
       allocate( dhy(nRxCSEM,nTxCSEM,nfreqCSEM) ) 
       allocate( dhz(nRxCSEM,nTxCSEM,nfreqCSEM) )   
       
          
    endif
                 
      
!
! MT arrays:
!    
    allocate( ex_mt(nRxMT,nfreqMT),ey_mt(nRxMT,nfreqMT), ez_mt(nRxMT,nfreqMT) )
    allocate( hx_mt(nRxMT,nfreqMT),hy_mt(nRxMT,nfreqMT), hz_mt(nRxMT,nfreqMT) )
    
    if (linversion) then         
         allocate( dex_mt_dsig(nRxMT,nfreqMT,nFree) )
         allocate( dhy_mt_dsig(nRxMT,nfreqMT,nFree) )
         allocate( dhz_mt_dsig(nRxMT,nfreqMT,nFree) )
 
         allocate( dhx_mt_dsig(nRxMT,nfreqMT,nFree) )
         allocate( dey_mt_dsig(nRxMT,nfreqMT,nFree) )
         allocate( dez_mt_dsig(nRxMT,nfreqMT,nFree) )
            
    endif

    end subroutine setup_solutionArrays    

!==================================================================================================================================! 
!================================================================================================================= setup_em2dkx_rxtx
!==================================================================================================================================!       
    subroutine setup_em2dkx_rxtx
!
! Allocates and specifies the Tx and Rx arrays needed by em2dkx
!

    integer                 :: j,iTx, icnt,nf,nk, iRx,jRx,jTx
    integer, dimension(7,3) :: txkxtab    
    real(8)                 :: cosdip, sinazi,sindip 

!
! Define a lookup table for setting the kxmodes based on the TxMode for each transmitter
!
    txkxtab(1,1) = 1 ! x
    txkxtab(2,1) = 2 ! y
    txkxtab(3,1) = 1 ! x,y
    txkxtab(3,2) = 2 
    txkxtab(4,1) = 3 ! z
    txkxtab(5,1) = 1 ! x,z
    txkxtab(5,2) = 3    
    txkxtab(6,1) = 4 ! yz
    txkxtab(7,1) = 1 ! x,yz
    txkxtab(7,2) = 4 


!
! Get modes needed for the point dipoles in xtrans,ytrans,ztrans arrays:
!    
    call getTxModes

!
! Create Rx arrays needed by em2dkx:
!    
    if (lMT) then
        
        nRx = nsiteMT
        
        allocate( xRx(nRx),yRx(nRx),zRx(nRx) )
        
        xRx = xsiteMT
        yRx = ysiteMT
        zRx = zsiteMT
        
    else
        
        nRx = nsiteCSEM
        
        allocate( xRx(nRx),yRx(nRx),zRx(nRx) )
        
        xRx = xsiteCSEM
        yRx = ysiteCSEM
        zRx = zsiteCSEM
        
    endif
    
!
! Create Tx arrays needed by em2dkx:
!    
    nTx = 0
    do iTx = 1,ntrans
        nTx = nTx + txcount(TxMode(iTx))        
    enddo
    
    allocate( jx(nTx),jy(nTx),jz(nTx),mx(nTx),my(nTx),mz(nTx) )
    allocate( xTx(nTx),yTx(nTx),zTx(nTx) )
    allocate( indexTrans_em2dkx(nTx) ) 
        
    icnt = 1
    do iTx = 1,ntrans
        
        do j = 1,txcount(TxMode(iTx))      
           
            
            !
            ! Store index to original input Tx:
            !
            indexTrans_em2dkx(icnt) = indexTrans(iTx)
            
            !
            ! Assign dipole moments:
            !
            jx(icnt) = 0d0; jy(icnt) = 0d0; jz(icnt) = 0d0
            mx(icnt) = 0d0; my(icnt) = 0d0; mz(icnt) = 0d0
            
            select case (trim(sourceType(iTx)))
            case ('edipole')
                select case ( txkxtab(TxMode(iTx),j) )
                case(1)
                    jx(icnt) = 1.d0
                case(2)
                    jy(icnt) = 1.d0
                case(3)
                    jz(icnt) = 1.d0
                case(4) ! yz transmitter, adjust jy and jz for dip
                    sinazi = sin(pi/180d0*azimuthtrans(iTx) )
                    cosdip = cos(pi/180d0*diptrans(iTx))
                    sindip = sin(pi/180d0*diptrans(iTx))
                  
                    jy(icnt) = sinazi*cosdip
                    jz(icnt) = sindip             
                end select
            case ('bdipole')
                select case ( txkxtab(TxMode(iTx),j) )
                case(1)
                    mx(icnt) = 1.d0
                case(2)
                    my(icnt) = 1.d0
                case(3)
                    mz(icnt) = 1.d0
                case(4) ! yz transmitter, adjust jy and jz for dip
                    sinazi = sin(pi/180d0*azimuthtrans(iTx))
                    cosdip = cos(pi/180d0*diptrans(iTx))
                    sindip = sin(pi/180d0*diptrans(iTx))
                  
                    my(icnt) = sinazi*cosdip
                    mz(icnt) = sindip             
                end select                            
            end select 
        
            !
            ! Assign dipole position:
            !
            xTx(icnt) = xtrans(iTx)
            yTx(icnt) = ytrans(iTx)
            zTx(icnt) = ztrans(iTx)
            
            icnt = icnt + 1
            
        enddo
          
    enddo    
    
    if (lMT) nTx = 1
    
    fileroot = outputFileroot
    
!
! Allocate kx domain solution arrays:  
!
    allocate( ex_kx(nRx,nTx),ey_kx(nRx,nTx),ez_kx(nRx,nTx),hx_kx(nRx,nTx),hy_kx(nRx,nTx),hz_kx(nRx,nTx)  ) 

 
    if (lMT) then 
        nf = nfreqMT
        nk = 1
    else
        nf = nFreqCSEM
        nk = nwave
    endif           
    
    allocate( ex_kx_s(nRx,nTx,nk,nf),ey_kx_s(nRx,nTx,nk,nf),ez_kx_s(nRx,nTx,nk,nf), &
            & hx_kx_s(nRx,nTx,nk,nf),hy_kx_s(nRx,nTx,nk,nf),hz_kx_s(nRx,nTx,nk,nf)  ) 
 
  
!
! Set lComps flags for em2dkx based on input data:
!  
    lComps = .false.
    
    if (lMT) then   
    
        lComps = .true.  ! compute them all since this is low bandwidth
    
    else
        
        if ( any(lCSEMworker(1:3,:)) ) then  ! any electrics requested
        
            if ( any(lCSEMworker(1,:)) )        lComps(1) = .true.
            if ( any(lCSEMworker(2,:)) )        lComps(2) = .true.
            if ( any(lCSEMworker(3,:)) )        lComps(3) = .true.
 
        endif
        
        if ( any(lCSEMworker(4:6,:)) ) then  ! any magnetics requested
        
            if ( any(lCSEMworker(4,:)) )        lComps(4) = .true.
            if ( any(lCSEMworker(5,:)) )        lComps(5) = .true.
            if ( any(lCSEMworker(6,:)) )        lComps(6) = .true.  
 
                               
        endif
                    
    endif

    
    if (linversion) then
        allocate( dex_kx(nRx,nTx) )   ! % dsig fields will be allocated elsewhere
        allocate( dey_kx(nRx,nTx) )
        allocate( dez_kx(nRx,nTx) )
        allocate( dhx_kx(nRx,nTx) ) 
        allocate( dhy_kx(nRx,nTx) )
        allocate( dhz_kx(nRx,nTx) ) 
        
        
        do iRx = 1,nRx
            jRx = indexSite(iRx)
            
            do iTx = 1,nTx
                
                if (lMT) then
                    if ( any(lDataMaskMT(jRx,:)) ) then
                        if (lComps(1)) allocate( dex_kx(iRx,iTx)%dsig(nFree))
                        if (lComps(2)) allocate( dey_kx(iRx,iTx)%dsig(nFree))
                        if (lComps(3)) allocate( dez_kx(iRx,iTx)%dsig(nFree))
                        if (lComps(4)) allocate( dhx_kx(iRx,iTx)%dsig(nFree))
                        if (lComps(5)) allocate( dhy_kx(iRx,iTx)%dsig(nFree))
                        if (lComps(6)) allocate( dhz_kx(iRx,iTx)%dsig(nFree))
                    endif
                else
                    jTx = indexTrans_em2dkx(iTx)
                    if ( any(lDataMaskCSEM(jRx,:,jTx)) ) then
                        if (lComps(1)) allocate( dex_kx(iRx,iTx)%dsig(nFree))  ! need to remember to zero these out in em2dkx...
                        if (lComps(2)) allocate( dey_kx(iRx,iTx)%dsig(nFree))
                        if (lComps(3)) allocate( dez_kx(iRx,iTx)%dsig(nFree))
                        if (lComps(4)) allocate( dhx_kx(iRx,iTx)%dsig(nFree))
                        if (lComps(5)) allocate( dhy_kx(iRx,iTx)%dsig(nFree))
                        if (lComps(6)) allocate( dhz_kx(iRx,iTx)%dsig(nFree))                        
                    endif     
                endif
                
            enddo
        enddo 
    endif 
        
    end subroutine setup_em2dkx_rxtx 
    
!==================================================================================================================================! 
!================================================================================================================== setup_KxFqGroups
!==================================================================================================================================!       
    subroutine setup_KxFqGroups

    integer :: i, icnt, ikxG, nk0, nkg, ikx, iFq
    
!
! Count how many KxFq groups we need:
!
    if (lMT) then
        nKxFqGroups = 1
    else
        nKxFqGroups = ceiling( dble(nwave)      / dble(nKxPerGroup)  )   
    endif
!
! Allocate the structure:
!
    allocate ( KxFqGroups(nKxFqGroups) )
    
!
! Fill in the values:
!     
    if (lMT) then
    
        KxFqGroups(1)%nKxFq = nFreqMT - 1
        
        allocate( KxFqGroups(1)%iFq(KxFqGroups(1)%nKxFq),KxFqGroups(1)%iKx(KxFqGroups(1)%nKxFq) )
        
        icnt = 0
        do i = 1,nFreqMT
            if ( mod(i+floor(dble(nFreqMT)/dble(2)),nFreqMT) ==0 ) then
                KxFqGroups(1)%iFqRefine = i
                KxFqGroups(1)%iKxRefine = 1
            else
                icnt = icnt + 1
                KxFqGroups(1)%iFq(icnt) = i
                KxFqGroups(1)%iKx(icnt) = 1
                if (lPrintDebug) write(*,'(a,2(i4,1x))') 'MT KxFqG refine icnt,iFq  : ',icnt,i
            endif
            
        enddo    
    
    else ! CSEM:
    
        do ikxG = 1,nKxFqGroups   
        
            nk0 =  (ikxG-1)*nKxPerGroup
            nkg = min(nKxPerGroup, nwave - nk0)   ! number of wavenumbers in this group 
            
            KxFqGroups(ikxG)%nKxFq = nkg*nFreqCSEM - 1 ! -1 for the refinement pair
            
            allocate( KxFqGroups(ikxG)%iFq(KxFqGroups(ikxG)%nKxFq),KxFqGroups(ikxG)%iKx(KxFqGroups(ikxG)%nKxFq) )
            
            icnt = 0
            do iFq = 1,nfreqCSEM
            
                do i = 1,nkg
                    
                    ikx = nk0 + i
                    
                    if ( ( mod(iFq + floor(dble(nFreqCSEM)/dble(2)),nFreqCSEM) ==0 ) .and.   &
                    &    ( mod(i + floor(dble(nkg)/dble(2)),nkg) ==0 ) ) then
                        KxFqGroups(ikxG)%iFqRefine = iFq
                        KxFqGroups(ikxG)%iKxRefine = ikx
                        if (lPrintDebug) write(*,'(a,3(i4,1x))') 'CSEM KxFqG refine ikxG,iFq,ikx  : ',ikxG,iFq,ikx  
                    else
                        icnt = icnt + 1
                        KxFqGroups(ikxG)%iFq(icnt) = iFq
                        KxFqGroups(ikxG)%iKx(icnt) = ikx
                        if (lPrintDebug) write(*,'(a,4(i4,1x))') 'CSEM KxFqG compute ikxG,icnt,iFq,ikx  : ',ikxG,icnt,iFq,ikx  
                    endif         
                    
                   
                enddo
            enddo
        enddo
  
    endif
 
    end subroutine setup_KxFqGroups 

!==================================================================================================================================! 
!=============================================================================================================== setup_ScratchGroups
!==================================================================================================================================!       
    subroutine setup_ScratchGroups
 
    integer     :: i,iTx,iRx,jTx,jRx, nRxTx, iPass
    logical     :: lIncrement
    
    if (lMT) return
!
! Count how many ScratchGroups groups we need:
!
    do iPass = 1,2 ! First pass counts nScratchGroup and nRxTx for each group, second pass fills in the indices
    
        nRxTx           = 0
        nScratchGroups  = 0
    
        do iTx = 1, nTx
        
            jTx = indexTrans_em2dkx(iTx) ! kwk debug this could be moved into .not.LMT below
        
            do iRx = 1,nRx
            
                jRx = indexSite(iRx)
                                
                lIncrement = .false.
                
                if ( any(lDataMaskCSEM(jRx,:,jTx)) ) lIncrement = .true.
    
                
                if (lIncrement) then
                    nRxTx = nRxTx + 1
                    if (nRxTx == nRxTxPerIFT + 1 ) nRxTx = 1
                    if (nRxTx == 1 ) nScratchGroups = nScratchGroups + 1 
                    
                    if (iPass == 2) then
                        ScratchGroups(nScratchGroups)%iRx(nRxTx) = iRx
                        ScratchGroups(nScratchGroups)%iTx(nRxTx) = iTx    
                    endif
                endif
    
            
            enddo ! iTx
        enddo  ! iRx
        
         if (iPass == 1) then
            allocate ( ScratchGroups(nScratchGroups) )       
            
            do i = 1,nScratchGroups-1
                allocate( ScratchGroups(i)%iRx(nRxTxPerIFT), ScratchGroups(i)%iTx(nRxTxPerIFT) )
                ScratchGroups(i)%nRxTx = nRxTxPerIFT
            enddo
            allocate( ScratchGroups(nScratchGroups)%iRx(nRxTx), ScratchGroups(nScratchGroups)%iTx(nRxTx) )
            ScratchGroups(nScratchGroups)%nRxTx = nRxTx     
                  
         endif   
               
    enddo ! iPass
 
 
    end subroutine setup_ScratchGroups     
!==================================================================================================================================! 
!======================================================================================================================= getTxModes
!==================================================================================================================================!       
    subroutine getTxModes
!
!  Sets the mode flag for each input CSEM transmitter:
! 
 
    integer :: i 
 
    real(8) :: cosdip, cosazi 
 
    allocate (TxMode(ntrans))
    
    do i = 1,ntrans
        cosdip = abs(cos(pi/180.*diptrans(i)))
        cosazi = abs(cos(pi/180.*azimuthtrans(i)))
        if (abs(cosdip-1.).le.1d-7) then       ! dip is n*pi, HED only
            if (abs(cosazi-1.).le.1d-7)   then ! azimuth is n*pi, +- x dipole only
                Txmode(i) = 1
            elseif  (abs(cosazi).le.1d-7)   then  ! y dipole only
                Txmode(i) = 2
            else                ! x and y dipole needed
                Txmode(i) = 3
            endif       
        elseif (abs(cosdip).le.1d-7) then  ! VED only
            Txmode(i) = 4
        else                        ! arbitrary dip, we need z and x and/or y dipoles 
            if (abs(cosazi-1.).le.1d-7)  then ! x and z needed
                Txmode(i) = 5
            elseif (abs(cosazi).le.1d-7) then  ! y and z needed
                Txmode(i) = 6
            else                ! x,y,z  needed 
                Txmode(i) = 7
            endif   
            
        endif     
                   
    enddo
    
    end subroutine getTxModes     
   
!==================================================================================================================================! 
!=================================================================================================================== computeAllKxFq
!==================================================================================================================================!       
    subroutine computeAllKxFq
 
    integer     :: ikx,ifq,ikxG,iKxFq, isubset
 
    isubset = 0
    
    do ikxG = 1,nKxFqGroups
        
        ikx  = KxFqGroups(ikxG)%iKxRefine
        ifq =  KxFqGroups(ikxG)%iFqRefine
        
        if (lMT) then
            w  = 2*pi*fTxMT(ifq)
            kx = 0
        else
            w  = 2*pi*fTxCSEM(ifq)
            kx = wavenum(ikx)        
        endif
        
        lRefine    = .true.
        meshnumber = 1
        
        if (lPrintDebug) write(*,*) 'call make_iDataMask(ifq): ',myid
        call make_iDataMask(ifq)
      
        isubset = isubset + 1
        
        if (lPrintDebug) write(*,*) 'call refinementKernel(isubset): ',myid
        call refinementKernel(isubset,ikx,ifq)
        
        do iKxFq = 1,KxFqGroups(ikxG)%nKxFq
    
            ikx  = KxFqGroups(ikxG)%iKx(iKxFq)
            ifq =  KxFqGroups(ikxG)%iFq(iKxFq)

            if (lMT) then
                w  = 2*pi*fTxMT(ifq)
                kx = 0
            else
                w  = 2*pi*fTxCSEM(ifq)
                kx = wavenum(ikx)        
            endif
            
            lRefine = .false.
            
            if (lPrintDebug) write(*,*) 'call make_iDataMask(ifq): ',myid
            call make_iDataMask(ifq)
            
            isubset = isubset + 1
            
            if (lPrintDebug) write(*,*) 'call refinementKernel(isubset): ',myid
            call refinementKernel(isubset,ikx,ifq)

        enddo
         
    enddo
 
    end subroutine computeAllKxFq

!==================================================================================================================================! 
!==================================================================================================================== make_iDataMask
!==================================================================================================================================!       
    subroutine make_iDataMask(iFq)
    
    integer, intent(in) :: iFq
    
    integer :: iTx,iRx,jRx,jTx
    
    if (.not.allocated(iDataMask)) allocate(iDataMask(nRx,nTx))
    
    iDataMask = 0
    
    if (lMT) then
    
        !
        ! MT Rx:
        !
        do iRx = 1,nRx
        
            jRx = indexSite(iRx)
                        
            if ( lDataMaskMT(jRx,iFq) ) then
                
                iDataMask(iRx,1) = 2
                
            elseif ( any(lDataMaskMT(jRx,:) ) ) then
                
                iDataMask(iRx,1) = 1
            
            endif 
            
        enddo     
         
    else
    
        !
        ! CSEM Rx and Tx:
        !
        do iTx = 1,nTx
            
            jTx = indexTrans_em2dkx(iTx)
            
            do iRx = 1,nRx
                
                jRx = indexSite(iRx)
                                          
                if ( lDataMaskCSEM(jRx,iFq,jTx) ) then
            
                    iDataMask(iRx,iTx) = 2
                
                elseif ( any(lDataMaskCSEM(jRx,:,jTx) ) ) then
                
                    iDataMask(iRx,iTx) = 1
            
                endif 
            enddo
        enddo       
        
 
    endif
    
    end subroutine make_iDataMask    
    
!==================================================================================================================================! 
!================================================================================================================== refinementKernel
!==================================================================================================================================!       
    subroutine refinementKernel(isubset,ikx,ifq)

    integer, intent(in) :: isubset, ikx,ifq
    integer             :: oldnodes, newnodes
    logical             :: lhaveNewMesh, lconv, lincreased
    character(256)      :: cadapt, filename, cGroup, cSubset     
    real(8)             :: lasterr, tstart, tend
    
    integer :: i,j
    real(8) :: wt,rTx
    complex(8), dimension(:,:), allocatable  ::  exlast,eylast,ezlast,hxlast,hylast,hzlast   
 
            
!
! Now we will compute the 2.5 response for this kx, and possibly do some mesh refinement:
!

    lComputeErrorEst    = lRefine     
    lconv               = .false.
    lasterr             = 1d100
    lincreased          = .false.
      
    write(cGroup,'(i6)') myID
    write(cSubset,'(i6)') isubset      
    fileroot = trim(outputFileroot)//'.'//trim(adjustl(cGroup))//'.'//trim(adjustl(cSubset))   
    
!
! Local refinement:
!
    
    if (meshnumber == 1) then
        
        call cpu_time(tstart)
        
        !
        ! Copy inputmesh to mesh:
        !
        if ( allocated( mesh%attr ) )       call deallocate_trimesh(mesh,.false.)
        if ( allocated( newmesh%attr ) )    call deallocate_trimesh(newmesh,.false.) 
        
        if (lPrintDebug) write(*,*) 'call copy_trimesh: ',myid
        call copy_trimesh(inputmesh,mesh)
        
        oldnodes = mesh%nnod  
       
        if ( lLocalRefine ) then
        
            if (lPrintDebug) write(*,*) 'call localRefinement: ',myid
            call localRefinement()  
            newnodes = mesh%nnod  
       
            if (mesh%nnod > oldnodes) then
            
                lhaveNewMesh = .true.  ! used to send new mesh back to manager later on    
                meshnumber   = meshnumber + 1   
             
            endif   
        
            call cpu_time(tend)
 
            if (lDisplayRefinementStats) call displayRefinementStats(isubset,meshnumber,lconv,oldnodes,newnodes,-2d0, tend-tstart)
        
        endif
        
        if (lSaveMeshFiles)   then  
            write(cadapt,'(i6)') meshnumber  
 
            cadapt = adjustl(cadapt)    
            filename  = trim(fileroot)//'.'//trim(cadapt)
            
            call write_triangle(mesh,filename)
        endif        
           
    endif

!
! Adaptive refinement:
!
    lCompDerivs  = .false. ! don't compute derivatives until we have a good mesh    
    lmadeNewMesh = .false. 
 
    if (meshnumber > maxnadapt ) lComputeErrorEst = .false.
    
    if (lComputeErrorEst) then
                
        allocate( exlast(nRx,nTx),eylast(nRx,nTx),ezlast(nRx,nTx),hxlast(nRx,nTx),hylast(nRx,nTx),hzlast(nRx,nTx)  ) 
        exlast = 0; eylast = 0; ezlast = 0; hxlast= 0; hylast= 0 ; hzlast = 0
                 
        do while (.not.lconv)
            
            !
            ! Compute the EM response:
            !
            call cpu_time(tstart)
            
            call em2dkx
 
            !
            ! As a backup to the error estimator, compute the numerical difference in the computed fields for each
            ! subsequent adaptive refinement iteration.
            !
            do i = 1,nTx
                do j = 1,nRx
                    
                    wt = 1.
                    if (.not.lMT) then
                        rTx  = sqrt( ( yRx(j)- yTx(i))**2 + ( zRx(j)- zTx(i))**2 )
                        if (rTx < minRangeProtector)  wt = 0.
                    endif
                    exlast(j,i) = wt*abs(ex_kx(j,i) - exlast(j,i)) / ( abs(ex_kx(j,i)) + ecutoff) 
                    eylast(j,i) = wt*abs(ey_kx(j,i) - eylast(j,i)) / ( abs(ey_kx(j,i)) + ecutoff) 
                    ezlast(j,i) = wt*abs(ez_kx(j,i) - ezlast(j,i)) / ( abs(ez_kx(j,i)) + ecutoff) 
                    hxlast(j,i) = wt*abs(hx_kx(j,i) - hxlast(j,i)) / ( abs(hx_kx(j,i)) + hcutoff) 
                    hylast(j,i) = wt*abs(hy_kx(j,i) - hylast(j,i)) / ( abs(hy_kx(j,i)) + hcutoff) 
                    hzlast(j,i) = wt*abs(hz_kx(j,i) - hzlast(j,i)) / ( abs(hz_kx(j,i)) + hcutoff) 
                enddo
            enddo
            
            !
            ! Test to see if the iterations are changing the responses significantly.
            ! if so the adaptive refinement has not converged to an accurate solution yet:
            !
            lconv = .true. 
            if ( lComps(1) .and. (maxval(abs(exlast)) > errortolerance/100. ) )  lconv = .false. 
            if ( lComps(2) .and. (maxval(abs(eylast)) > errortolerance/100. ) )  lconv = .false.
            if ( lComps(3) .and. (maxval(abs(ezlast)) > errortolerance/100. ) )  lconv = .false.
            if ( lComps(4) .and. (maxval(abs(hxlast)) > errortolerance/100. ) )  lconv = .false.
            if ( lComps(5) .and. (maxval(abs(hylast)) > errortolerance/100. ) )  lconv = .false.
            if ( lComps(6) .and. (maxval(abs(hzlast)) > errortolerance/100. ) )  lconv = .false.
             
            exlast = ex_kx; eylast = ey_kx; ezlast = ez_kx; hxlast= hx_kx; hylast= hy_kx ; hzlast = hz_kx
                       
            !
            ! If new mesh was made:
            !
            newnodes = 0
            oldnodes = mesh%nnod  
            if (lmadeNewMesh) then
                
                meshnumber   = meshnumber + 1  
                
                lhaveNewMesh = .true.  ! used to send new mesh back to manager later on    
                
                newnodes = newmesh%nnod
                oldnodes = mesh%nnod   
                
                ! Tolerance threshold, give up if only a "few" elements refined and its not the first mesh:
                if ( (newnodes - oldnodes < 5).and.(meshnumber>1) )   lconv = .true.
                
                ! Give up if we've already refined the mesh many times:
                if (meshnumber > maxnadapt ) lconv = .true.
                
                ! Stop refinement if the mesh has grown too large:
                if ( newnodes > maxMeshNodes) lconv = .true.
                               
                !
                ! Move newmesh to mesh:
                !
                call deallocate_trimesh( mesh, .false. )
                call copy_trimesh( newmesh,mesh )
                call deallocate_trimesh( newmesh, .false. )
                newmesh%nnod = 0
                if (lSaveMeshFiles)   then  
                    write(cadapt,'(i6)') meshnumber  
                    write(cGroup,'(i6)') myID
                    write(cSubset,'(i6)') isubset   
                    cadapt = adjustl(cadapt)    
                    filename  = trim(fileroot)//'.'//trim(cadapt)
            
                    call write_triangle(mesh,filename)
                endif   
                
            else ! no new mesh made, must have converged
                
                lconv  = .true.                                                     
 
            endif      
            
            !
            ! Check for convergence:
            !
            if ( maxerr < errortolerance/1d2)  lconv  = .true.          
            
            !
            ! Catch for runaway refinement:
            !
            if ( (maxerr > lasterr) .and. (lincreased) ) then
                lconv = .true. ! error has increased twice, abort adaptive refinement...
            elseif (maxerr > lasterr) then
                lincreased = .true.
            else
                lincreased = .false.
            endif
            
            lasterr = maxerr
         
            call cpu_time(tend)
                   
            if (lDisplayRefinementStats) call displayRefinementStats(iSubset,meshnumber,lconv,oldnodes,newnodes,maxerr,tend-tstart)

            
        enddo  
        
        deallocate( exlast,eylast,ezlast,hxlast,hylast,hzlast  ) 
        
        !
        ! At this point lconv = .true. If a new mesh was made we need to compute its solution, without requiring the error estimates
        !
        if ( lmadeNewMesh ) then
           
            if (linversion) lCompDerivs = .true. 
            
            lComputeErrorEst = .false.
            maxerr           = -1d-2      
            
            call cpu_time(tstart)
            
            call em2dkx
            
            oldnodes = mesh%nnod    

            call cpu_time(tend)
                    
            if (lDisplayRefinementStats) call displayRefinementStats(iSubset,meshnumber,lconv,oldnodes,0,maxerr,tend-tstart)

         
                    
        elseif (linversion) then ! The while loop didn't make a new mesh,  but we still need to compute the inversion derivatives:
            
            lCompDerivs      = .true.     
            lComputeErrorEst = .false.
            maxerr           = -1d-2  
            
            call cpu_time(tstart)
            
            call em2dkx           
            
            oldnodes = mesh%nnod   

            call cpu_time(tend)

            if (lDisplayRefinementStats) call displayRefinementStats(iSubset,meshnumber,lconv,oldnodes,0,maxerr,tend-tstart)

                
        endif
        
    else  ! (.not. lComputeErrorEst) then  ! no mesh refinement, just a straightforward computation:  
         
        if ( linversion )  lCompDerivs = .true. 
        maxerr = -1d0  ! special flag for displayRef
         
        call cpu_time(tstart)
            
        call em2dkx            
        
        lconv = .true.
        oldnodes = mesh%nnod
        
        call cpu_time(tend)

        if (lDisplayRefinementStats) call displayRefinementStats(iSubset,meshnumber,lconv,oldnodes,0,maxerr, tend-tstart)
  
    endif
 
        
!
! Move kx domain fields into holding arrays for temporary storage:
!
    ex_kx_s(:,:,ikx,ifq) = ex_kx
    ey_kx_s(:,:,ikx,ifq) = ey_kx
    ez_kx_s(:,:,ikx,ifq) = ez_kx    
    hx_kx_s(:,:,ikx,ifq) = hx_kx
    hy_kx_s(:,:,ikx,ifq) = hy_kx
    hz_kx_s(:,:,ikx,ifq) = hz_kx     
    
!
! Save kx domain sensitivity derivatives to scratch space
!
    if (linversion) then
    
    !
    ! If MT, just copy the derivatives across dsig arrays:
    !
        if (lMT) then
        
            call moveDerivsMT(iFq,isubset)
            
        else   
        
        !
        ! For CSEM data, save to scratch space:
        !
        
            call cpu_time(tstart)

            call writeScratchFile(iFq,iKx)  

            call cpu_time(tend)

            if (lprintMPItimers)  write(*, '((a5,2x,i6,2x,a,f12.5))')  'Proc:',myID,  ' Worker writescratch:',  tend - tstart

        endif
             
     endif       
     
     end subroutine refinementKernel


!==================================================================================================================================! 
!=============================================================================================================== legendre_compute_dr
!==================================================================================================================================!   
    subroutine legendre_compute_dr ( order, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_DR computes a Gauss-Legendre rule, Davis-Rabinowitz method.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be greater than 0.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
!
    implicit none
    
    integer ( kind = 4 ) order
    
    real    ( kind = 8 ) d1
    real    ( kind = 8 ) d2pn
    real    ( kind = 8 ) d3pn
    real    ( kind = 8 ) d4pn
    real    ( kind = 8 ) dp
    real    ( kind = 8 ) dpn
    real    ( kind = 8 ) e1
    real    ( kind = 8 ) fx
    real    ( kind = 8 ) h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) iback
    integer ( kind = 4 ) k
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mp1mi
    integer ( kind = 4 ) ncopy
    integer ( kind = 4 ) nmove
    real    ( kind = 8 ) p
    real    ( kind = 8 ) :: pi = 3.141592653589793D+00
    real    ( kind = 8 ) pk
    real    ( kind = 8 ) pkm1
    real    ( kind = 8 ) pkp1
    real    ( kind = 8 ) t
    real    ( kind = 8 ) u
    real    ( kind = 8 ) v
    real    ( kind = 8 ) x0
    real    ( kind = 8 ) xtab(order)
    real    ( kind = 8 ) xtemp
    real    ( kind = 8 ) weight(order)
    
    if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
    end if
    
    e1 = real ( order * ( order + 1 ), kind = 8 )
    
    m = ( order + 1 ) / 2
    
    do i = 1, m
    
    mp1mi = m + 1 - i
    
    t = real ( 4 * i - 1, kind = 8 ) * pi &
      / real ( 4 * order + 2, kind = 8 )
    
    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( order, kind = 8 ) ) &
      / real ( 8 * order * order, kind = 8 ) )
    
    pkm1 = 1.0D+00
    pk = x0
    
    do k = 2, order
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do
    
    d1 = real ( order, kind = 8 ) * ( pkm1 - x0 * pk )
    
    dpn = d1 / ( 1.0D+00 - x0 * x0 )
    
    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 * x0 )
    
    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 * x0 )
    
    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) / &
      ( 1.0D+00 - x0 * x0 )
    
    u = pk / dpn
    v = d2pn / dpn
    !
    !  Initial approximation H:
    !
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
    !
    !  Refine H using one step of Newton's method:
    !
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )
    
    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )
    
    h = h - p / dp
    
    xtemp = x0 + h
    
    xtab(mp1mi) = xtemp
    
    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )
    
    weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp * xtemp ) / ( fx * fx )
    
    end do
    
    if ( mod ( order, 2 ) == 1 ) then
    xtab(1) = 0.0D+00
    end if
    !
    !  Shift the data up.
    !
    nmove = ( order + 1 ) / 2
    ncopy = order - nmove
    
    do i = 1, nmove
    iback = order + 1 - i
    xtab(iback) = xtab(iback-ncopy)
    weight(iback) = weight(iback-ncopy)
    end do
    !
    !  Reflect values for the negative abscissas.
    !
    do i = 1, order - nmove
    xtab(i) = - xtab(order+1-i)
    weight(i) = weight(order+1-i)
    end do
    
    end subroutine legendre_compute_dr

!==================================================================================================================================!
!============================================================================================================ displayRefinementStats
!==================================================================================================================================!
    subroutine displayRefinementStats( isubset, meshnumber,lconv, oldnodes,newnodes,maxerr,time)
    
    integer ::  meshnumber, oldnodes, newnodes,isubset
    real(8) :: maxerr, time
    logical :: lconv
    
    character(256) :: sFmt
    
    character(32) :: smaxerr, snewnodes, status, stime
    
    status = 'adaptive'
     
    if (maxerr > 0) then
        write(smaxerr,'(a11,2x,g9.2,1x,a1)') 'max(Error):',maxerr*1d2,'%'
    elseif (maxerr == -2d0) then
        status = 'local' ! special flag for local refinement mode  
        smaxerr = ' '
     elseif (maxerr == -1d0) then
        status = 'shared' ! special flag for local refinement mode  
        smaxerr = ' '   
        lconv = .false.
    else
        smaxerr = ' '
        
    endif
    if (newnodes > 0) then
        write(snewnodes,'(a3,2x,i6)') '-->',newnodes
    else
        snewnodes   = ' '
    endif       
   
    if (lconv) then
        status = 'converged'
    endif       
    
    write(stime,'(a8,f9.2,a2)') ' Timer: ',time, ' s'
    
    
    sFmt = '(a5,2x,i6,2x,a6,2x,i6,2x,a7,2x,i3,2x,a7,2x,i2,2x,a7,2x,a9,2x,a8,2x,i6,2x,a11,2x,a23,2x,a19)'     
    write(*,sFmt)  'Proc:',myID,'Group:',iRefinementGrp,'Subset:',isubset, 'Mesh #:',meshnumber,'Status:', trim(status), &
                & '# Nodes:',oldnodes,  trim(snewnodes), trim(smaxerr), trim(stime)
    
    end subroutine displayRefinementStats 


!==================================================================================================================================! 
!================================================================================================================== writeScratchFile
!==================================================================================================================================! 
    subroutine writeScratchFile(iFq,iKx)
!
! Writes scratch files by peeling off nRxTxPerIFT Rx,Tx pairs from the RxTxGroup iGrp  
!    
    use em2dkx_mod
    
    implicit none 
    
    integer, intent(in)     :: iKx, iFq 

    integer                 :: iScratch, ios, posStart, iRxTx, iRx,iTx,jRx, jTx
    character(256)          :: cFilename     
 
        
    do iScratch = 1,nScratchGroups
        
        !
        ! Get a name for the file:    
        !
        call getScratchFilename(iRefinementGrp,iScratch,iFq,iKx,cFilename)
         
        !  
        ! Open the file:            
        !
        open(unit=21,file=trim(cFilename),form="unformatted",status="replace",access="stream",iostat=ios)
        if (ios /= 0) then
            write(*,*) 'error, creating scratch file: ',trim(cFilename)
            stop
        endif
        inquire(21,pos=posStart)
        write(21,iostat=ios) 0      ! dummy value for now     

        !
        ! Write out the data:
        !
        do iRxTx = 1,ScratchGroups(iScratch)%nRxTx
        
            iTx = ScratchGroups(iScratch)%iTx(iRxTx)
            iRx = ScratchGroups(iScratch)%iRx(iRxTx)
            
            jRx = indexSite(iRx)
            jTx = indexTrans_em2dkx(iTx)
            
            if (lDataMaskCSEM(jRx,iFq,jTx)) then
        
                !
                ! Write data from this Rx,Tx to the file:
                !
                if (lCSEMworker(1,jRx)) write(21,iostat=ios)  dex_kx(iRx,iTx)%dsig(1:nFree)
                if (lCSEMworker(2,jRx)) write(21,iostat=ios)  dey_kx(iRx,iTx)%dsig(1:nFree)
                if (lCSEMworker(3,jRx)) write(21,iostat=ios)  dez_kx(iRx,iTx)%dsig(1:nFree)
                if (lCSEMworker(4,jRx)) write(21,iostat=ios)  dhx_kx(iRx,iTx)%dsig(1:nFree)
                if (lCSEMworker(5,jRx)) write(21,iostat=ios)  dhy_kx(iRx,iTx)%dsig(1:nFree)
                if (lCSEMworker(6,jRx)) write(21,iostat=ios)  dhz_kx(iRx,iTx)%dsig(1:nFree)   
            
                if (ios /= 0) then
                    write(*,*) 'error, writing scratch file: ',trim(cFilename), 'ios: ',ios
                    stop
                endif
                
            endif

        enddo


        ! Rewind and write status indicator:
        write(21,iostat=ios, pos=posStart) nFree  ! insert this at the end of i/o so it can be used as a file i/o status indicator
    
        ! Close the file:
        close(21,iostat=ios)   
        if (ios /= 0) then
            write(*,*) 'error, closing scratch file: ',trim(cFilename)
            stop
        endif
     
    enddo
  
    end subroutine writeScratchFile
    
!==================================================================================================================================! 
!================================================================================================================ getScratchFilename
!==================================================================================================================================! 
    subroutine getScratchFilename(iGrp,iFileNum,iFq,iKx,cFilename)  
    
    implicit none 
    
    integer, intent(in)         :: iGrp,iKx,iFq,iFileNum
    character(256), intent(out) :: cFilename
    
    character(256) ::  cGrp, cFileNum, cKx, cFq
 
    write(cGrp,'(i12)') iGrp  
    cGrp = adjustl(cGrp)   
     
    write(cFileNum,'(i12)') iFileNum  
    cFileNum = adjustl(cFileNum)    
    
     
    write(cFq,'(i12)') iFq  
    cFq = adjustl(cFq)     
    
    write(cKx,'(i12)') iKx  
    cKx = adjustl(cKx)    
    
    cFilename  = trim(scratchFolder)//'mare2dem_tmp.'//trim(cGrp)//'.'//trim(cFileNum)//'.'//trim(cFq)//'.'//trim(cKx) 
    
    end subroutine getScratchFilename    

!==================================================================================================================================! 
!=========================================================================================================== inverseFourierTransform
!==================================================================================================================================! 
    subroutine inverseFourierTransform
 
    use SinCosFilters
 
    integer         :: iTx,iRx, iFq, i0,i1,sym1,sym2, jTx,jRx, kTx, nt, nift
    real(8)         :: xr, tStart,tEnd    
  
    character(256)  :: sFmt
    character(32)   :: stime
    
    if (lMT) return
    
    call cpu_time(tStart) 
    
 !
 ! Initialize digital filters
 !
    call spline_integ_kx_setup()
 
    nift = 0
    
    do iTx = 1,nTx
        
        jTx = indexTrans_em2dkx(iTx)
        
        !
        ! Get dipole type:
        !
        if ( (abs(jx(iTx))>0) .or.  (abs(jy(iTx))>0) .or.  (abs(jz(iTx))>0) ) then ! electric dipole
            i0 = 0  ! 0 is for even cosine transform
            i1 = 1  ! 1 is for odd sine transform
        else  ! magnetic dipole
            i0 = 1
            i1 = 0
        end if       
        
        !
        ! Set symmetry variables based on dipole type (e or b) and direction of dipole (x or yz)
        !
        if ( (abs(jx(iTx))>0) .or. (abs(mx(iTx))>0) ) then  
            sym1 = i0   
            sym2 = i1  
        else ! yz plane
            sym1 = i1
            sym2 = i0            
        endif    

        do iRx = 1,nRx
            
            jRx = indexSite(iRx)
        

            !
            ! Compute the inverse Fourier transform using spline itegration with digital filters:
            !     
            xr  = xRx(iRx) - xTx(iTx)

            do iFq = 1,nFreqCSEM

                if (lDataMaskCSEM( jRx, iFq, jTx) ) then

                    if (lCSEMworker(1,jRx)) call spline_integ_fkx(xr, ex_kx_s(iRx,iTx,:,iFq), sym1, ex_kx_s(iRx,iTx,1,iFq))   
                    if (lCSEMworker(2,jRx)) call spline_integ_fkx(xr, ey_kx_s(iRx,iTx,:,iFq), sym2, ey_kx_s(iRx,iTx,1,iFq))
                    if (lCSEMworker(3,jRx)) call spline_integ_fkx(xr, ez_kx_s(iRx,iTx,:,iFq), sym2, ez_kx_s(iRx,iTx,1,iFq))          
                    if (lCSEMworker(4,jRx)) call spline_integ_fkx(xr, hx_kx_s(iRx,iTx,:,iFq), sym2, hx_kx_s(iRx,iTx,1,iFq))
                    if (lCSEMworker(5,jRx)) call spline_integ_fkx(xr, hy_kx_s(iRx,iTx,:,iFq), sym1, hy_kx_s(iRx,iTx,1,iFq))
                    if (lCSEMworker(6,jRx)) call spline_integ_fkx(xr, hz_kx_s(iRx,iTx,:,iFq), sym1, hz_kx_s(iRx,iTx,1,iFq))
        
                    nift = nift + count(lCSEMworker(:,jRx))

                endif

            enddo ! iFq

        enddo ! iRx = 1,nRx 

    enddo ! iTx
 

    call spline_integ_kx_deallocate
 
! 
! Now combine the modal transmitters to get the responses of arbitrarily oriented transmitters:
! 
    
    kTx = 1
    
    do iTx = 1,ntrans
        
        jTx = indexTrans(iTx)
           
        select case (Txmode(iTx))
        case(1,2,4,6)
            nt = 1
        case(3,5,7)
            nt = 2
        end select
       
         
        do iRx = 1,nRx
         
            jRx = indexSite(iRx)
            
            do iFq = 1,nFreqCSEM

                if (lDataMaskCSEM( jRx, iFq, jTx) ) then

                    call combineTxModes(ex_kx_s(iRx,kTx:kTx+nt-1,1,ifq),azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                    ex_kx_s(iRx,iTx,1,ifq) = ex_kx_s(iRx,kTx,1,ifq) ! slide result over to left
        
                    call combineTxModes(ey_kx_s(iRx,kTx:kTx+nt-1,1,ifq),azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                    ey_kx_s(iRx,iTx,1,ifq) = ey_kx_s(iRx,kTx,1,ifq) ! slide result over to left
                                                
                    call combineTxModes(ez_kx_s(iRx,kTx:kTx+nt-1,1,ifq),azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                    ez_kx_s(iRx,iTx,1,ifq) = ez_kx_s(iRx,kTx,1,ifq) ! slide result over to left                                

                    call combineTxModes(hx_kx_s(iRx,kTx:kTx+nt-1,1,ifq),azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                    hx_kx_s(iRx,iTx,1,ifq) = hx_kx_s(iRx,kTx,1,ifq) ! slide result over to left
        
                    call combineTxModes(hy_kx_s(iRx,kTx:kTx+nt-1,1,ifq),azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                    hy_kx_s(iRx,iTx,1,ifq) = hy_kx_s(iRx,kTx,1,ifq) ! slide result over to left
                                                
                    call combineTxModes(hz_kx_s(iRx,kTx:kTx+nt-1,1,ifq),azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                    hz_kx_s(iRx,iTx,1,ifq) = hz_kx_s(iRx,kTx,1,ifq) ! slide result over to left   
                                    
                endif
            
            enddo ! iFq
        enddo ! iRx
        
        kTx = kTx + nt
        
    enddo
 
    
    call cpu_time(tEnd) 
 
    write(stime,'(a8,f9.2,a2)') ' Timer: ',tEnd - tStart, ' s'
  
    sFmt = '(a5,2x,i6,2x,        a6,2x,i6,2x,  5x,  9x,  46x, a32,2x,i7,  2x,a19)'     
    if (lDisplayRefinementStats) write(*,sFmt)  'Proc:',myID,'Group:',iRefinementGrp, '# Field Transforms:',nift,   trim(stime)
                      
    end subroutine inverseFourierTransform
    
!==================================================================================================================================! 
!==================================================================================================== inverseFourierTransform_Derivs
!==================================================================================================================================! 
    subroutine inverseFourierTransform_Derivs

    integer     :: iScratch, nf, iRx,jRx,iTx,jTx,iFq
    real(8)     :: tstart, tend

    if (lMT) return
!
! Allocate storage for large working arrays: 
! 
    allocate( readBuffer(count(lComps)*nFree*min(nRxTxPerIFT,nRx*nTx),nwave) )   
    
 
    nf = nfreqCSEM
    
        
    call cpu_time(tStart)     

    
    if (lComps(1)) allocate( dexds(nRx,nTx,nf) )   
    if (lComps(2)) allocate( deyds(nRx,nTx,nf) )   
    if (lComps(3)) allocate( dezds(nRx,nTx,nf) )   
    if (lComps(4)) allocate( dhxds(nRx,nTx,nf) )   
    if (lComps(5)) allocate( dhyds(nRx,nTx,nf) )   
    if (lComps(6)) allocate( dhzds(nRx,nTx,nf) )   
    
    
    do iTx = 1,nTx
        jTx = indexTrans_em2dkx(iTx)
        
        do iRx = 1,nRx
            jRx = indexSite(iRx)
        
            do iFq = 1,nf
                         
                if ( lDataMaskCSEM(jRx,iFq,jTx) ) then
                    if (lCSEMworker(1,jRx)) allocate( dexds(iRx,iTx,iFq)%dsig(nFree) )   
                    if (lCSEMworker(2,jRx)) allocate( deyds(iRx,iTx,iFq)%dsig(nFree) )   
                    if (lCSEMworker(3,jRx)) allocate( dezds(iRx,iTx,iFq)%dsig(nFree) )   
                    if (lCSEMworker(4,jRx)) allocate( dhxds(iRx,iTx,iFq)%dsig(nFree) )   
                    if (lCSEMworker(5,jRx)) allocate( dhyds(iRx,iTx,iFq)%dsig(nFree) )   
                    if (lCSEMworker(6,jRx)) allocate( dhzds(iRx,iTx,iFq)%dsig(nFree) )   
                endif             
          
                        
            enddo
        enddo
    enddo
    
    call cpu_time(tend)    
    if (lPrintDebug) write(*,*) 'time to allocate derivative arrays, myID,  : ', myID  , tend- tstart
   
!
! Loop over frequencies:
!     
    do iFq = 1,nf
    
        !
        ! Read in all scratch files:
        !
        do iScratch = 1,nScratchGroups
            
            ! Load data:
            if (lPrintDebug) write(*,*) 'call loadScratch myID, iScratch: ', myID, iScratch
            call cpu_time(tStart)  
            
            call loadScratch(iScratch,iFq)
            
            call cpu_time(tend)    
            if (lPrintDebug) write(*,*) 'time to load scratch files, myID,  : ', myID  , tend- tstart
            
            ! Move from readBuffer into deriv arrays (dexds etc), performing inverse FT if CSEM data:
            if (lPrintDebug) write(*,*) 'call moveIntoDerivs myID, iScratch: ', myID, iScratch
            
            call cpu_time(tStart)  
            
            call moveIntoDerivs(iScratch,iFq)              
            
            call cpu_time(tend)    
            if (lPrintDebug) write(*,*) 'time to inverse FT deriv arrays, myID,  : ', myID  , tend- tstart
            
        enddo
        
        !
        ! Combine the modal Tx's to get the derivative responses for arbitrarily oriented transmitters:
        !
        if (lPrintDebug) write(*,*) 'call combineDerivModes myID, iFq: ', myID, iFq
        call combineDerivModes(iFq)  
        
    enddo
    
    deallocate(readBuffer)

    if (lPrintDebug) write(*,*) 'done with inverseFourierTransform_Derivs ', myID 
    
    end subroutine inverseFourierTransform_Derivs    

!==================================================================================================================================! 
!==================================================================================================================== combineTxModes
!==================================================================================================================================!          
   subroutine combineTxModes(f,aziTx,dipTx,modeTx)
!   
! Combines the fundamental kxmode results into the fields from arbitrarily 
! oriented transmitters. 
!
! Results returned in f(1)
 
    integer,intent(in)       :: modeTx
    complex(4),intent(inout) :: f(2)
    real(8),intent(in)       :: aziTx,dipTx 
 
 
    real(4)     :: ca,sa,cd,sd 
 
    
    ! get Tx angles:
    ca = cos( aziTx*pi/180. )   
    sa = sin( aziTx*pi/180. )  
    cd = cos( dipTx*pi/180. )
    sd = sin( dipTx*pi/180. )
    
 
    select case ( modeTx )

    case (1) ! x dipole

        ! ca gives correct phase for 0 and +-180 degree azimuths
        f(1) = ca*f(1)
                
    case (2) ! y dipole

        ! sa gives correct phase for +-90,270 degree azimuths              
        f(1) = sa*f(1)
            
    case (3)  ! we need to combine x and y dipoles  

        f(1) = ca*f(1) + sa*f(2)
            
    case (4) ! z dipole
        
        f(1) = sd*f(1)

    case (5)  ! we need to combine x and z dipoles

        f(1) = cd*f(1) + sd*f(2)
  
    case (6) ! dipping in y-z plane, em2dkx handles yz dip already
        !  Full sign of azimuth and dip accounted for already
        !f(1) = f(1)

    case (7) ! completely an arbitrary angle

         ! Combine x and y-z transmitter fields:
         ! Note yz fields already have sa*cd applied.  
         f(1) = ca*cd*f(1) + f(2)    

    end select
 
     
    end subroutine combineTxModes     
    
!==================================================================================================================================! 
!================================================================================================================= integrate_dipoles
!==================================================================================================================================! 
    subroutine integrate_dipoles 
  
    integer :: iTx, iRx, iFq, jTx, jRx, nt, nr, jRx0, kTx, l, i, nm,j, nskip 
    
    complex(8), dimension(:), allocatable :: temparray
    
    if (linversion) allocate( temparray(nFree) )
       
! 
! Integrate finite length transmitter dipoles:
!    

    if (lPrintDebug) write(*,*) 'integrating finite Tx dipoles...'
    jTx = 1
    kTx = 1 ! for derivs
    do iTx = 1,nTxCSEM ! Loop over input transmitters
        
        if ( any(lDataMaskCSEM(:,:,iTx)) ) then ! only add dipoles needed by data
            
            select case (Txmode(jTx))
            case(1,2,4,6)
                nm = 1
            case(3,5,7)
                nm = 2
            end select
            
            if (lengthTxCSEM(iTx) > 0) then  ! only integrate if finite dipole present
                nt = nQuadTxCSEM 
                
                do iRx = 1,nsiteCSEM
                    jRx = indexSite(iRx)
                
                    do iFq = 1,nfreqCSEM    
                
                        if (lDataMaskCSEM(jRx,iFq,iTx)) then
  
                            if (lCSEMworker(1,jRx))  ex_kx_s(iRx,jTx,1,iFq) = sum( ex_kx_s(iRx,jTx:jTx+nt-1,1,iFq)*quad_weights_Tx ) / 2.d0  
                            if (lCSEMworker(2,jRx))  ey_kx_s(iRx,jTx,1,iFq) = sum( ey_kx_s(iRx,jTx:jTx+nt-1,1,iFq)*quad_weights_Tx ) / 2.d0  
                            if (lCSEMworker(3,jRx))  ez_kx_s(iRx,jTx,1,iFq) = sum( ez_kx_s(iRx,jTx:jTx+nt-1,1,iFq)*quad_weights_Tx ) / 2.d0  
                        
                            if (lCSEMworker(4,jRx))  hx_kx_s(iRx,jTx,1,iFq) = sum( hx_kx_s(iRx,jTx:jTx+nt-1,1,iFq)*quad_weights_Tx ) / 2.d0  
                            if (lCSEMworker(5,jRx))  hy_kx_s(iRx,jTx,1,iFq) = sum( hy_kx_s(iRx,jTx:jTx+nt-1,1,iFq)*quad_weights_Tx ) / 2.d0  
                            if (lCSEMworker(6,jRx))  hz_kx_s(iRx,jTx,1,iFq) = sum( hz_kx_s(iRx,jTx:jTx+nt-1,1,iFq)*quad_weights_Tx ) / 2.d0                           
                             
                            if (linversion) then
                            
                                if (lCSEMworker(1,jRx)) then
                                    temparray = 0
                                    do l = 1,nQuadTxCSEM
                                        temparray =  temparray + dexds(iRx,kTx+nm*(l-1),iFq)%dsig*quad_weights_Tx(l)      
                                    enddo
                                    dexds(iRx,kTx,iFq)%dsig = temparray/2.d0
                                endif
                                
                                if (lCSEMworker(2,jRx)) then
                                    temparray = 0
                                    do l = 1,nQuadTxCSEM
                                        temparray =  temparray + deyds(iRx,kTx+nm*(l-1),iFq)%dsig*quad_weights_Tx(l)      
                                    enddo
                                    deyds(iRx,kTx,iFq)%dsig = temparray/2.d0
                                endif

                                if (lCSEMworker(3,jRx)) then
                                    temparray = 0
                                    do l = 1,nQuadTxCSEM
                                        temparray =  temparray + dezds(iRx,kTx+nm*(l-1),iFq)%dsig*quad_weights_Tx(l)      
                                    enddo
                                    dezds(iRx,kTx,iFq)%dsig = temparray/2.d0
                                endif

                                if (lCSEMworker(4,jRx)) then
                                    temparray = 0
                                    do l = 1,nQuadTxCSEM
                                        temparray =  temparray + dhxds(iRx,kTx+nm*(l-1),iFq)%dsig*quad_weights_Tx(l)      
                                    enddo
                                    dhxds(iRx,kTx,iFq)%dsig = temparray/2.d0
                                endif
                                
                                if (lCSEMworker(5,jRx)) then
                                    temparray = 0
                                    do l = 1,nQuadTxCSEM
                                        temparray =  temparray + dhyds(iRx,kTx+nm*(l-1),iFq)%dsig*quad_weights_Tx(l)      
                                    enddo
                                    dhyds(iRx,kTx,iFq)%dsig = temparray/2.d0
                                endif

                                if (lCSEMworker(6,jRx)) then
                                    temparray = 0
                                    do l = 1,nQuadTxCSEM
                                        temparray =  temparray + dhzds(iRx,kTx+nm*(l-1),iFq)%dsig*quad_weights_Tx(l)      
                                    enddo
                                    dhzds(iRx,kTx,iFq)%dsig = temparray/2.d0
                                endif
                                                                                                
                            endif   

                        endif   ! if (lDataMaskCSEM(jRx,iFq,iTx)) then                    
                                    
                    enddo ! iFq
                enddo ! iRx
                
            else
                nt = 1
            endif
            
            kTx = kTx + nt*nm
            jTx = jTx + nt
        
        endif
  
    enddo

!
! Integrate finite length CSEM receiver dipoles:
!  
    if (lPrintDebug) write(*,*) 'integrating finite Rx dipoles...'               
 
    jRx = 1       
    do iRx = 1,nRxCSEM
                  
        if ( any(lDataMaskCSEM(iRx,:,:)) ) then ! only add dipoles needed by data
        
            if ( lengthRxCSEM(iRx) > 0 ) then ! finite dipole    
                 
                jRx0 = jRx ! all integrated fields will be in first jRx for this receiver
                
                jTx = 1
                kTx = 1  ! for derivs
                
                do iTx = 1,nTxCSEM
                    
                    if ( any(lDataMaskCSEM(:,:,iTx)) ) then ! only increment if Tx has data present 
 
                        select case (Txmode(jTx))
                        case(1,2,4,6)
                            nm = 1
                        case(3,5,7)
                            nm = 2
                        end select
                    
                        do iFq = 1,nfreqCSEM
        
                            ! For magnetic fields, take the middle quadrature point (assumes odd order!) for point dipole coil
                            if (lComps(4))  hx_kx_s(jRx,jTx,1,iFq) = hx_kx_s(jRx+(nQuadRxCSEM-1)/2,jTx,1,iFq)
                            if (lComps(5))  hy_kx_s(jRx,jTx,1,iFq) = hy_kx_s(jRx+(nQuadRxCSEM-1)/2,jTx,1,iFq)
                            if (lComps(6))  hz_kx_s(jRx,jTx,1,iFq) = hz_kx_s(jRx+(nQuadRxCSEM-1)/2,jTx,1,iFq)      
                        
                            jRx = jRx0
                        
                            ! For electrics, integrate along 3 possible dipole lines                   
                            do j = 1,3

                            
                                if (lCompCSEM(j,iRx)) then ! this receiver & electric dipole component has data, add to counter:
                                              
                                    if (lCSEMworker(1,iRx)) then ! ex finite dipole:
                                        ex_kx_s(jRx,jTx,1,iFq) = sum( ex_kx_s(jRx:jRx+nQuadRxCSEM-1,jTx,1,iFq)*quad_weights_RxCSEM ) / 2.d0  
                                    endif
                                    if (lCSEMworker(2,iRx)) then ! ey finite dipole:
                                        ey_kx_s(jRx,jTx,1,iFq) = sum( ey_kx_s(jRx:jRx+nQuadRxCSEM-1,jTx,1,iFq)*quad_weights_RxCSEM ) / 2.d0  
                                    endif                    
                                    if (lCSEMworker(3,iRx)) then ! ez finite dipole:
                                        ez_kx_s(jRx,jTx,1,iFq) = sum( ez_kx_s(jRx:jRx+nQuadRxCSEM-1,jTx,1,iFq)*quad_weights_RxCSEM ) / 2.d0  
                                    endif                               
                                
                                    if ( lDataMaskCSEM(iRx,iFq,iTx) ) then  ! since these buggers are only allocated if data present:
                                
                                    if (linversion) then
                             
                                        ! For magnetic fields, take the middle quadrature point (assumes odd order!) for point dipole coil
                                        if (lCSEMworker(4,iRx))  then
                                            do i = 1,nFree
                                                dhxds(jRx,kTx,iFq)%dsig(i)  = dhxds(jRx+(nQuadRxCSEM-1)/2,kTx,iFq)%dsig(i)
                                            enddo
                                        endif
                                        if (lCSEMworker(5,iRx))  then
                                            do i = 1,nFree
                                                dhyds(jRx,kTx,iFq)%dsig(i)  = dhyds(jRx+(nQuadRxCSEM-1)/2,kTx,iFq)%dsig(i)
                                            enddo
                                        endif
                                        if (lCSEMworker(6,iRx))  then
                                            do i = 1,nFree
                                                dhzds(jRx,kTx,iFq)%dsig(i)  = dhzds(jRx+(nQuadRxCSEM-1)/2,kTx,iFq)%dsig(i)
                                            enddo
                                        endif                  
                                      
                                                  
                                        ! Electric fields:
                                        if (lCSEMworker(1,iRx)) then
                                            temparray = 0
                                            do l = 1,nQuadRxCSEM
                                                temparray =  temparray + dexds(jRx+l-1,kTx,iFq)%dsig*quad_weights_RxCSEM(l)      
                                            enddo
                                            dexds(jRx,kTx,iFq)%dsig = temparray/2.d0
                                        endif  
                                   
                                        if (lCSEMworker(2,iRx)) then
                                            temparray = 0
                                            do l = 1,nQuadRxCSEM
                                                temparray =  temparray + deyds(jRx+l-1,kTx,iFq)%dsig*quad_weights_RxCSEM(l)      
                                            enddo                                    
                                            deyds(jRx,kTx,iFq)%dsig = temparray/2.d0
                                        endif    
 
                                        if (lCSEMworker(3,iRx)) then
                                            temparray = 0
                                            do l = 1,nQuadRxCSEM
                                                temparray =  temparray + dezds(jRx+l-1,kTx,iFq)%dsig*quad_weights_RxCSEM(l)      
                                            enddo
                                            dezds(jRx,kTx,iFq)%dsig = temparray/2.d0
                                        endif              
                                                                
                                    endif ! linversion
                                
                                    endif
                                
                                     jRx = jRx + nQuadRxCSEM   
                            
                                endif
                        
                            enddo
   
                        enddo  
                        if (lengthTxCSEM(iTx) > 0) then  
                            nt = nQuadTxCSEM            
                        else 
                            nt = 1
                        endif
                    
                        kTx = kTx + nt*nm
                    
                        jTx = jTx + nt     
                                 
                    endif  !   if ( any(lDataMaskCSEM(:,:,iTx)) ) then 
                    
                enddo ! iTx
                
                ! increment counter by the number of receiver point dipoles used for up to 3 components:
                jRx = jRx0
                do j = 1,3
                    if (lCompCSEM(j,iRx)) then
                        jRx = jRx + nQuadRxCSEM
                    endif
                enddo
            else
        
               jRx = jRx + 1 ! single point dipole
        
            endif
            
        endif  ! if ( any(lDataMaskCSEM(iRx,:,:)) ) then ! only add dipoles needed by data
        
    enddo        

!
! Transfer CSEM results into output solution arrays
!  
    if (lPrintDebug) write(*,*) 'Transferring CSEM results into output solution arrays...'
    
    jTx = 1
    kTx = 1
    do iTx = 1,nTxCSEM ! Loop over input transmitters

        
        if ( any(lDataMaskCSEM(:,:,iTx)) ) then ! only add dipoles needed by data
                           
            if (lengthTxCSEM(iTx) > 0) then  ! only integrate if finite dipole present
                nt = nQuadTxCSEM 
            else
                nt = 1
            endif
            
            select case (Txmode(jTx))
            case(1,2,4,6)
                nm = 1
            case(3,5,7)
                nm = 2
            end select
                                
            jRx = 1  
            
            do iRx = 1,nRxCSEM
                
                if ( any(lDataMaskCSEM(iRx,:,:)) ) then ! only add dipoles needed by data
                
                    if ( lengthRxCSEM(iRx) > 0 ) then ! finite dipole    
                        nr = 0 
                        if (lCompCSEM(1,iRx))   nr = nr + nQuadRxCSEM    
                        if (lCompCSEM(2,iRx))   nr = nr + nQuadRxCSEM    
                        if (lCompCSEM(3,iRx))   nr = nr + nQuadRxCSEM 
                    else
                        nr = 1
                    endif     
                          
                      
                    do iFq = 1,nfreqCSEM
                       
                        
                        if (lDataMaskCSEM(iRx,iFq,iTx)) then
                            
                            nskip = 0 
                            
                            ex(iRx,iTx,iFq) = ex_kx_s(jRx,jTx,1,iFq)  
                            if ( (lCompCSEM(1,iRx)) .and. (lengthRxCSEM(iRx) > 0) ) nskip = nskip + nQuadRxCSEM
                            
                            ey(iRx,iTx,iFq) = ey_kx_s(jRx+nskip,jTx,1,iFq)  
                            if ( (lCompCSEM(2,iRx)) .and. (lengthRxCSEM(iRx) > 0) ) nskip = nskip + nQuadRxCSEM
                            
                            ez(iRx,iTx,iFq) = ez_kx_s(jRx+nskip,jTx,1,iFq)  
                            
                            hx(iRx,iTx,iFq) = hx_kx_s(jRx,jTx,1,iFq)  
                            hy(iRx,iTx,iFq) = hy_kx_s(jRx,jTx,1,iFq)  
                            hz(iRx,iTx,iFq) = hz_kx_s(jRx,jTx,1,iFq)  
                    
                    
                            if (linversion) then
                                nskip = 0 
                                    
                                if (lCompCSEM(1,iRx)) then
                                    allocate (dex(iRx,iTx,iFq)%dsig(nFree) )
                                    dex(iRx,iTx,iFq)%dsig = dexds(jRx,kTx,iFq)%dsig   
                                    deallocate(dexds(jRx,kTx,iFq)%dsig)
                                    if ( (lCompCSEM(1,iRx)) .and. (lengthRxCSEM(iRx) > 0) ) nskip = nskip + nQuadRxCSEM
                                endif
                                if (lCompCSEM(2,iRx)) then
                                    allocate ( dey(iRx,iTx,iFq)%dsig(nFree) )
                                    dey(iRx,iTx,iFq)%dsig = deyds(jRx+nskip,kTx,iFq)%dsig
                                    deallocate( deyds(jRx,kTx,iFq)%dsig )
                                    if ( (lCompCSEM(2,iRx)) .and. (lengthRxCSEM(iRx) > 0) ) nskip = nskip + nQuadRxCSEM
                                endif  
                                if (lCompCSEM(3,iRx)) then
                                    allocate ( dez(iRx,iTx,iFq)%dsig(nFree) )
                                    dez(iRx,iTx,iFq)%dsig = dezds(jRx+nskip,kTx,iFq)%dsig
                                    deallocate( dezds(jRx,kTx,iFq)%dsig )
                                endif
                                if (lCompCSEM(4,iRx)) then
                                    allocate ( dhx(iRx,iTx,iFq)%dsig(nFree) )
                                    dhx(iRx,iTx,iFq)%dsig = dhxds(jRx,kTx,iFq)%dsig  
                                    deallocate ( dhxds(jRx,kTx,iFq)%dsig ) 
                                endif
                                if (lCompCSEM(5,iRx)) then
                                    allocate ( dhy(iRx,iTx,iFq)%dsig(nFree) )
                                    dhy(iRx,iTx,iFq)%dsig = dhyds(jRx,kTx,iFq)%dsig
                                    deallocate( dhyds(jRx,kTx,iFq)%dsig )
                                endif
                                if (lCompCSEM(6,iRx)) then
                                    allocate ( dhz(iRx,iTx,iFq)%dsig(nFree) )
                                    dhz(iRx,iTx,iFq)%dsig = dhzds(jRx,kTx,iFq)%dsig 
                                    deallocate( dhzds(jRx,kTx,iFq)%dsig )  
                                endif             
                            endif
                    
                        endif
                
                    enddo ! iFq
                    
                    jRx = jRx + nr
                
                endif !  any(lDataMaskCSEM(iRx,:,:))
                            
            enddo ! iRx
         
            jTx = jTx + nt
            kTx = kTx + nt*nm
            
        endif ! any(lDataMaskCSEM(:,:,iTx))
        
    enddo ! iTx
    
    
!
! Integrate finite length MT receiver dipoles:
!   
    if (lPrintDebug)  write(*,*) 'Integrating finite length MT receiver dipoles...' ! kwk debug: need to fix this for finite MT Rx dipoles! 
    jRx = 1   
    do iRx = 1,nRxMT
                  
        if ( lengthRxMT(iRx) > 0 ) then ! finite dipole    
            
            jRx0 = jRx
            
            do iFq = 1,nfreqMT
         
                ! For magnetic fields, take the middle quadrature point (assumes odd order!) for point dipole coil
                if (lComps(4))  hx_kx_s(jRx,1,1,iFq) = hx_kx_s(jRx+(nQuadRxMT-1)/2,1,1,iFq)
                if (lComps(5))  hy_kx_s(jRx,1,1,iFq) = hy_kx_s(jRx+(nQuadRxMT-1)/2,1,1,iFq)
                if (lComps(6))  hz_kx_s(jRx,1,1,iFq) = hz_kx_s(jRx+(nQuadRxMT-1)/2,1,1,iFq)           
                 
                jRx = jRx0
                 
                ! For electrics, integrate along 3 possible dipole lines                   
                do j = 1,3
                         
                    if (lCompMT(1,iRx)) then ! this receiver & electric dipole component has data 
                        ex_kx_s(jRx,1,1,iFq) = sum( ex_kx_s(jRx:jRx+nQuadRxMT-1,1,1,iFq)*quad_weights_RxMT ) / 2.d0 
                    endif
                    if (lCompMT(2,iRx)) then ! this receiver & electric dipole component has data 
                        ey_kx_s(jRx,1,1,iFq) = sum( ey_kx_s(jRx:jRx+nQuadRxMT-1,1,1,iFq)*quad_weights_RxMT ) / 2.d0 
                    endif
                    if (lCompMT(3,iRx)) then ! this receiver & electric dipole component has data 
                        ez_kx_s(jRx,1,1,iFq) = sum( ez_kx_s(jRx:jRx+nQuadRxMT-1,1,1,iFq)*quad_weights_RxMT ) / 2.d0 
                    endif
      
                    if (linversion) then
                
                    
                        ! For magnetic fields, take the middle quadrature point (assumes odd order!) for point dipole coil
                        if (lCompMT(4,iRx))  then
                            do i = 1,nFree
                                dhxds(jRx,1,iFq)%dsig(i)  = dhxds(jRx+(nQuadRxMT-1)/2,1,iFq)%dsig(i)
                            enddo
                        endif
                        if (lCompMT(5,iRx))  then
                            do i = 1,nFree
                                dhyds(jRx,1,iFq)%dsig(i)  = dhyds(jRx+(nQuadRxMT-1)/2,1,iFq)%dsig(i)
                            enddo
                        endif
                        if (lCompMT(6,iRx))  then
                            do i = 1,nFree
                                dhzds(jRx,1,iFq)%dsig(i)  = dhzds(jRx+(nQuadRxMT-1)/2,1,iFq)%dsig(i)
                            enddo
                        endif                  
                                                                                    
                        ! Electric fields:
                        if (lCompMT(1,iRx)) then
                            temparray = 0
                            do l = 1,nQuadRxMT
                                temparray =  temparray + dexds(jRx+l-1,1,iFq)%dsig*quad_weights_RxMT(l)      
                            enddo
                            dexds(jRx,1,iFq)%dsig = temparray/2.d0
                        endif  
                                          
                        if (lCompMT(2,iRx)) then
                            temparray = 0
                            do l = 1,nQuadRxMT
                                temparray =  temparray + deyds(jRx+l-1,1,iFq)%dsig*quad_weights_RxMT(l)      
                            enddo
                            deyds(jRx,1,iFq)%dsig = temparray/2.d0
                        endif    
                    
                        if (lCompMT(3,iRx)) then
                            temparray = 0
                            do l = 1,nQuadRxMT
                                temparray =  temparray + dezds(jRx+l-1,1,iFq)%dsig*quad_weights_RxMT(l)      
                            enddo
                            dezds(jRx,1,iFq)%dsig = temparray/2.d0
                        endif              
                                                                  
                    endif ! linversion    
                    
                    jRx = jRx + nQuadRxMT              
 
                enddo ! j=1,3
            enddo
            
        else
    
            jRx = jRx + 1 ! single point dipole
    
        endif
  
    enddo                  
    
!
! Transfer MT results into output solution arrays
!           
    if (lPrintDebug)  write(*,*) 'transferring MT results into solution arrays...'       
  
    jRx = 1  
    do iRx = 1,nRxMT
         
       ! if ( any(lDataMaskMT(iRx,:)) ) then ! only add dipoles needed by data    
                
            if ( lengthRxMT(iRx) > 0 ) then ! finite dipole    
                nr = 0 
                if (lCompMT(1,iRx))   nr = nr + nQuadRxMT     
                if (lCompMT(2,iRx))   nr = nr + nQuadRxMT    
                if (lCompMT(3,iRx))   nr = nr + nQuadRxMT 
            else
                nr = 1    
            endif
 
            do iFq = 1,nfreqMT
                
                nskip = 0
                if (lDataMaskMT(iRx,iFq)) then
                    
                    ex_mt(iRx,iFq) = conjg(ex_kx_s(jRx,1,1,iFq))  
                    if ( (lCompMT(1,iRx)) .and. (lengthRxMT(iRx) > 0) ) nskip = nskip + nQuadRxMT
                    ey_mt(iRx,iFq) = conjg(ey_kx_s(jRx+nskip,1,1,iFq))  
                    if ( (lCompMT(2,iRx)) .and. (lengthRxMT(iRx) > 0) ) nskip = nskip + nQuadRxMT
                    ez_mt(iRx,iFq) = conjg(ez_kx_s(jRx+nskip,1,1,iFq))  
                    
                    hx_mt(iRx,iFq) = conjg(hx_kx_s(jRx,1,1,iFq))  
                    hy_mt(iRx,iFq) = conjg(hy_kx_s(jRx,1,1,iFq))  
                    hz_mt(iRx,iFq) = conjg(hz_kx_s(jRx,1,1,iFq))  

                    if (linversion) then
                        nskip = 0 
                         
                        if (lCompMT(1,iRx)) then
                            dex_mt_dsig(iRx,iFq,:) = conjg(dexds(jRx,1,iFq)%dsig) 
                            if (lengthRxMT(iRx) > 0) nskip = nskip + nQuadRxMT 
                        endif

                        if (lCompMT(2,iRx)) then
                            dey_mt_dsig(iRx,iFq,:) = conjg(deyds(jRx+nskip,1,iFq)%dsig) 
                            if (lengthRxMT(iRx) > 0) nskip = nskip + nQuadRxMT  
                        endif    

                        if (lCompMT(3,iRx)) then
                            dez_mt_dsig(iRx,iFq,:) = conjg(dezds(jRx+nskip,1,iFq)%dsig)   
                        endif
                        if (lCompMT(4,iRx)) dhx_mt_dsig(iRx,iFq,:) = conjg(dhxds(jRx,1,iFq)%dsig)   
                        if (lCompMT(5,iRx)) dhy_mt_dsig(iRx,iFq,:) = conjg(dhyds(jRx,1,iFq)%dsig)   
                        if (lCompMT(6,iRx)) dhz_mt_dsig(iRx,iFq,:) = conjg(dhzds(jRx,1,iFq)%dsig)                 
                
                    endif
                                
                endif
        
            enddo ! iFq
            jRx = jRx + nr
      !  endif        
    enddo ! iRx
    
    
    if ( allocated( temparray ) )   deallocate (temparray)
  
       
    end subroutine integrate_dipoles
    
!==================================================================================================================================! 
!======================================================================================================================= loadScratch
!==================================================================================================================================! 
    subroutine loadScratch(iScratch,iFq)
 
    integer, intent(in) :: iScratch, iFq   
    
    integer             :: iKx, nData, nk, iRx, jRx, iTx, jTx, i
    
    character(256)      :: cFilename  
     
  
    nk = nwave  
    
    nData = 0
    do i = 1,ScratchGroups(iScratch)%nRxTx
    
        iRx = ScratchGroups(iScratch)%iRx(i)
        jRx = indexSite(iRx)
        iTx = ScratchGroups(iScratch)%iTx(i)
        jTx = indexTrans_em2dkx(iTx)
        
        if (lDataMaskCSEM(jRx,iFq,jTx)) then
           if (lCSEMworker(1,jRx) ) nData = nData + nFree
           if (lCSEMworker(2,jRx) ) nData = nData + nFree
           if (lCSEMworker(3,jRx) ) nData = nData + nFree
           if (lCSEMworker(4,jRx) ) nData = nData + nFree
           if (lCSEMworker(5,jRx) ) nData = nData + nFree
           if (lCSEMworker(6,jRx) ) nData = nData + nFree 
        endif
    
    enddo
 
    
    do iKx = 1,nk
    
        ! Get a name for the file:     
        call getScratchFilename(iRefinementGrp,iScratch,iFq,iKx,cFilename)
 
        ! Read in the scratch file when it is ready to be read:
        call readScratchFileWhenReady(cFilename, nData, iKx)
        
    enddo
       
    end subroutine loadScratch

!==================================================================================================================================! 
!========================================================================================================== readScratchFileWhenReady
!==================================================================================================================================! 
    subroutine readScratchFileWhenReady(cFilename, nData, iKx)

#if defined(__INTEL_COMPILER)
    use IFPORT          ! DGM 8/23/2013 Windows Intel compiler requirement for the sleepqq routine
#endif
 
    integer, intent(in)         :: iKx, nData
    character(256), intent(in)  :: cFilename  
 
    logical :: lReady
    integer :: nFreeTest, ios 

!
! First inquire to see if the file exists:
! 
    lReady = .false. 
    
    do while(.not.lReady)
    
        inquire(file=cFileName,exist=lReady)    

        if (lReady) exit 

        !
        ! Else take a chill pill:
        !   
#if defined(__INTEL_COMPILER)  
        call sleepqq(10)     ! sleepqq time is in milliseconds.
#endif    
    
    
    enddo

!
! Next open file and see if marker integer has been written, indication file i/o is done
!        
    do  
    
        open(unit=21,file=trim(cFilename),status='old',access="stream",form="unformatted",action="read",iostat=ios)
        if (ios /= 0) then
            write(*,*) 'error, opening scratch file: ',trim(cFilename),  ' iostat:' ,ios
            stop
        endif

        read(unit=21,iostat=ios) nFreeTest
         if (ios /= 0) then
            write(*,*) 'error, reading nFreetest from scratch file: ',trim(cFilename),  ' iostat:' ,ios
            stop
        endif        
         
        if (nFreeTest > 0) exit

        !
        ! Else take a chill pill:
        !   
#if defined(__INTEL_COMPILER)  
        call sleepqq(10)     ! sleepqq time is in milliseconds.
#endif    
    
    
    enddo    
 
!
! Now read in the file's data:
!    
    read(21,iostat=ios)  readBuffer(1:nData,iKx)   

    if (ios /= 0) then
        write(*,*) 'error, reading from scratch file on unit: ',21
        stop
    endif     
    
!
! Finally, close and delete the file:
!       

    close(21, status='delete')
    
    end subroutine readScratchFileWhenReady
!==================================================================================================================================! 
!====================================================================================================================== moveDerivsMT
!==================================================================================================================================! 
    subroutine moveDerivsMT(iFq,isubset)

    integer, intent(in) ::  iFq, isubset
    
    
    
    integer         :: iRx,jRx
 
    if (isubset == 1) then 
        if (lComps(1)) allocate( dexds(nRx,1,nfreqMT) )   
        if (lComps(2)) allocate( deyds(nRx,1,nfreqMT) )   
        if (lComps(3)) allocate( dezds(nRx,1,nfreqMT) )   
        if (lComps(4)) allocate( dhxds(nRx,1,nfreqMT) )   
        if (lComps(5)) allocate( dhyds(nRx,1,nfreqMT) )   
        if (lComps(6)) allocate( dhzds(nRx,1,nfreqMT) )    
    endif 
    
!
! Loop over all Rx,Tx pairs and work it:
!        
    do iRx = 1,nRx
        jRx = indexSite(iRx)

        if (lDataMaskMT(jRx,iFq)) then

            if (lComps(1)) then 
                allocate(dexds(iRx,1,iFq)%dsig(nFree))
                dexds(iRx,1,iFq)%dsig =  dex_kx(iRx,1)%dsig
            endif
          
            if (lComps(2)) then 
                allocate(deyds(iRx,1,iFq)%dsig(nFree))
                deyds(iRx,1,iFq)%dsig =  dey_kx(iRx,1)%dsig
            endif
          
            if (lComps(3)) then 
                allocate(dezds(iRx,1,iFq)%dsig(nFree))
                dezds(iRx,1,iFq)%dsig =  dez_kx(iRx,1)%dsig
            endif
 
            if (lComps(4)) then 
                allocate(dhxds(iRx,1,iFq)%dsig(nFree))
                dhxds(iRx,1,iFq)%dsig =  dhx_kx(iRx,1)%dsig
            endif
          
            if (lComps(5)) then 
                allocate(dhyds(iRx,1,iFq)%dsig(nFree))
                dhyds(iRx,1,iFq)%dsig =  dhy_kx(iRx,1)%dsig
            endif
          
            if (lComps(6)) then 
                allocate(dhzds(iRx,1,iFq)%dsig(nFree))
                dhzds(iRx,1,iFq)%dsig =  dhz_kx(iRx,1)%dsig
            endif
      
      endif                             
    enddo
         
    end subroutine moveDerivsMT
!==================================================================================================================================! 
!==================================================================================================================== moveIntoDerivs
!==================================================================================================================================! 
    subroutine moveIntoDerivs(iScratch,iFq)

    integer, intent(in) :: iScratch, iFq
    
    integer         :: iRx,iTx,jRx,jTx,iRxTx, icnt,i0,i1,sym1,sym2, i
    real(8)         :: xr, tStart,tEnd    
    complex(4)      :: temp
  
    character(256)  :: sFmt
    character(32)   :: stime
    

    call cpu_time(tStart)  
        
    icnt = 0 

!
! Initialize digital filters
!
    call spline_integ_kx_setup()
 
!
! Loop over all Rx,Tx pairs and work it:
!        
    do iRxTx = 1,ScratchGroups(iScratch)%nRxTx
        
        iRx = ScratchGroups(iScratch)%iRx(iRxTx)
        iTx = ScratchGroups(iScratch)%iTx(iRxTx)
        
        jRx = indexSite(iRx)
        jTx = indexTrans_em2dkx(iTx)
        
        if (lDataMaskCSEM(jRx,iFq,jTx)) then
            
            xr  = xRx(iRx) - xTx(iTx)
            
            !
            ! Get dipole type:
            !
            if ( (abs(jx(iTx))>0) .or.  (abs(jy(iTx))>0) .or.  (abs(jz(iTx))>0) ) then ! electric dipole
                i0 = 0  ! 0 is for even cosine transform
                i1 = 1  ! 1 is for odd sine transform
            else  ! magnetic dipole
                i0 = 1
                i1 = 0
            end if       
        
            !
            ! Set symmetry variables based on dipole type (e or b) and direction of dipole (x or yz)
            !
            if ( (abs(jx(iTx))>0) .or. (abs(mx(iTx))>0) ) then  
                sym1 = i0   
                sym2 = i1  
            else ! yz plane
                sym1 = i1
                sym2 = i0            
            endif    
    
            if (lCSEMworker(1,jRx)) then
                do i = 1,nFree
                    call spline_integ_fkx(xr, readBuffer(icnt+1,1:nwave), sym1, temp )
                    dexds(iRx,iTx,iFq)%dsig(i) = temp 
                    icnt = icnt + 1
                enddo  
            endif
            if (lCSEMworker(2,jRx)) then
                do i = 1,nFree
                    call spline_integ_fkx(xr, readBuffer(icnt+1,1:nwave), sym2, temp ) 
                    deyds(iRx,iTx,iFq)%dsig(i) = temp                    
                    icnt = icnt + 1
                enddo  
            endif
             if (lCSEMworker(3,jRx)) then
                do i = 1,nFree
                    call spline_integ_fkx(xr, readBuffer(icnt+1,1:nwave), sym2, temp ) 
                    dezds(iRx,iTx,iFq)%dsig(i) = temp
                    icnt = icnt + 1
                enddo  
            endif                    
            if (lCSEMworker(4,jRx)) then
                do i = 1,nFree
                    call spline_integ_fkx(xr, readBuffer(icnt+1,1:nwave), sym2, temp )
                     dhxds(iRx,iTx,iFq)%dsig(i) = temp 
                    icnt = icnt + 1
                enddo  
            endif
            if (lCSEMworker(5,jRx)) then
                do i = 1,nFree
                    call spline_integ_fkx(xr, readBuffer(icnt+1,1:nwave), sym1, temp )
                    dhyds(iRx,iTx,iFq)%dsig(i) = temp 
                    icnt = icnt + 1
                enddo  
            endif
             if (lCSEMworker(6,jRx)) then
                do i = 1,nFree
                    call spline_integ_fkx(xr, readBuffer(icnt+1,1:nwave), sym1, temp )
                    dhzds(iRx,iTx,iFq)%dsig(i) = temp 
                    icnt = icnt + 1
                enddo  
            endif                    
     
        endif        

    enddo

    call spline_integ_kx_deallocate
    
    call cpu_time(tEnd) 
             
    if (.not.lMT) then
    
        write(stime,'(a8,f9.2,a2)') ' Timer: ',tEnd - tStart, ' s'
  
        sFmt = '(a5,2x,i6,2x,        a6,2x,i6,2x,  5x, 46x, 9x, a32,2x,i7,    2x,a19)'     
        if (lDisplayRefinementStats) write(*,sFmt) 'Proc:',myID,'Group:',iRefinementGrp, '# Derivative Transforms:',icnt,trim(stime)
    
    endif
         
    end subroutine moveIntoDerivs

!==================================================================================================================================! 
!================================================================================================================= combineDerivModes
!==================================================================================================================================! 
    subroutine combineDerivModes(iFq)
    
    integer, intent(in) :: iFq
    
    integer             :: iTx, jTx, iRx, jRx, nt, kTx, i, j
    
    complex(4) :: field(2)
    
    kTx = 1
    do iTx = 1,ntrans
        
        jTx = indexTrans(iTx)
           
        select case (Txmode(iTx))
        case(1,2,4,6)
            nt = 1
        case(3,5,7)
            nt = 2
        end select
       
         
        do iRx = 1,nRx
         
            jRx = indexSite(iRx)   

            if (lDataMaskCSEM( jRx, iFq, jTx) ) then
                
                if (lCSEMworker(1,jRx)) then
                    do i = 1,nFree                            
                        do j = 1,nt
                            field(j) = dexds(iRx,kTx+j-1,iFq)%dsig(i)
                        enddo
                        call combineTxModes(field,azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                        dexds(iRx,kTx,iFq)%dsig(i) = field(1)   
                    enddo
                endif
                
                if (lCSEMworker(2,jRx)) then
                    do i = 1,nFree                            
                        do j = 1,nt
                            field(j) = deyds(iRx,kTx+j-1,iFq)%dsig(i)
                        enddo
                        call combineTxModes(field,azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                        deyds(iRx,kTx,iFq)%dsig(i) = field(1)    
                    enddo
                endif

                if (lCSEMworker(3,jRx)) then
                    do i = 1,nFree                            
                        do j = 1,nt
                            field(j) = dezds(iRx,kTx+j-1,iFq)%dsig(i)
                        enddo
                        call combineTxModes(field,azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                        dezds(iRx,kTx,iFq)%dsig(i) = field(1)    
                    enddo
                endif                                                
                
                if (lCSEMworker(4,jRx)) then
                    do i = 1,nFree                            
                        do j = 1,nt
                            field(j) = dhxds(iRx,kTx+j-1,iFq)%dsig(i)
                        enddo
                        call combineTxModes(field,azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                        dhxds(iRx,kTx,iFq)%dsig(i) = field(1)   
                    enddo
                endif
                
                if (lCSEMworker(5,jRx)) then
                    do i = 1,nFree                            
                        do j = 1,nt
                            field(j) = dhyds(iRx,kTx+j-1,iFq)%dsig(i)
                        enddo
                        call combineTxModes(field,azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                        dhyds(iRx,kTx,iFq)%dsig(i) = field(1)   
                    enddo
                endif

                if (lCSEMworker(6,jRx)) then
                    do i = 1,nFree                            
                        do j = 1,nt
                            field(j) = dhzds(iRx,kTx+j-1,iFq)%dsig(i)
                        enddo
                        call combineTxModes(field,azimuthtrans(iTx),diptrans(iTx), TxMode(iTx))
                        dhzds(iRx,kTx,iFq)%dsig(i) = field(1)   
                    enddo
                endif              
            endif
            
   
        enddo ! iRx
        
        kTx = kTx + nt
        
    enddo   
    
    end subroutine combineDerivModes       
         
!==================================================================================================================================! 
!============================================================================================================= rotateToRxOrientation
!==================================================================================================================================! 
    subroutine rotateToRxOrientation

    integer                     :: iRx, jRx, iTx, jTx, nf, iFq  
    real(8)                     :: theta, alpha, beta           ! site rotation angles
    real(8) , dimension(3,3)    :: RotR 
    complex(8), dimension(3)    :: tempvec                      ! temporary vector of ex,ey,ez for rotations   
    complex(8), dimension(:,:), allocatable :: temparray        ! temporary array rotations  
    
    logical   :: lRotate, lCC(6)
    
    allocate(temparray(nFree,3))      
     
    if (lMT) then
        nf = nfreqMT
    else 
        nf = nFreqCSEM
    endif     

    do iRx = 1,nRx
        
        jRx = indexSite(iRx)
        
        ! Get site rotation parameters:
        if (lMT) then
            theta = ThetaRxMT(jRx)
            alpha = AlphaRxMT(jRx)
            beta  = BetaRxMT(jRx)          
        else        
            theta = ThetaRxCSEM(jRx)
            alpha = AlphaRxCSEM(jRx)
            beta  = BetaRxCSEM(jRx)        
        endif
        
        if ( (theta /= 0) .or. (alpha /= 0) .or. ( beta /= 0) ) then

            call getRotationMatrix(theta,alpha,beta, RotR)
            RotR = transpose(RotR) ! since i'm using right multiplication below, take transpose to get correct rotation matrix
            
            ! Rotate Responses:
            do iTx = 1,nTx
                if (.not.lMT) jTx = indexTrans_em2dkx(iTx)
                
                do iFq = 1,nf  
     
                    lRotate = .false. 
                    if (lMT) then
                        if (lDataMaskMT(jRx,iFq))  lRotate = .true.
                    else
                        if (lDataMaskCSEM(jRx,iFq,jTx)) lRotate = .true.
                    endif
                    
                    if (lRotate) then 

                        tempvec = 0d0 

                        tempvec(1) = ex_kx_s(iRx,iTx,1,iFq)
                        tempvec(2) = ey_kx_s(iRx,iTx,1,iFq)
                        tempvec(3) = ez_kx_s(iRx,iTx,1,iFq)

                        tempvec = matmul(tempvec,RotR)

                        ex_kx_s(iRx,iTx,1,iFq) = tempvec(1) 
                        ey_kx_s(iRx,iTx,1,iFq) = tempvec(2) 
                        ez_kx_s(iRx,iTx,1,iFq) = tempvec(3) 
 
                        tempvec = 0d0 

                        tempvec(1) = hx_kx_s(iRx,iTx,1,iFq)
                        tempvec(2) = hy_kx_s(iRx,iTx,1,iFq)
                        tempvec(3) = hz_kx_s(iRx,iTx,1,iFq)

                        tempvec = matmul(tempvec,RotR)

                        hx_kx_s(iRx,iTx,1,iFq) = tempvec(1) 
                        hy_kx_s(iRx,iTx,1,iFq) = tempvec(2) 
                        hz_kx_s(iRx,iTx,1,iFq) = tempvec(3)             
            
                      
                        if (linversion) then   ! Rotate sensitivities as well:
                                
                                lCC = .false.
                                
                                if (lMT) then
                                    lCC = lComps
                                else
                                    lCC = lCSEMworker(1:6,jRx)
                                endif   
                                                                 
                                temparray = 0
                                              
                                if (lCC(1)) temparray(:,1) = dexds(iRx,iTx,iFq)%dsig
                                if (lCC(2)) temparray(:,2) = deyds(iRx,iTx,iFq)%dsig
                                if (lCC(3)) temparray(:,3) = dezds(iRx,iTx,iFq)%dsig
                            
                                temparray  = matmul(temparray,RotR)
                 
                                if (lCC(1))  dexds(iRx,iTx,iFq)%dsig = temparray(:,1) 
                                if (lCC(2))  deyds(iRx,iTx,iFq)%dsig = temparray(:,2) 
                                if (lCC(3))  dezds(iRx,iTx,iFq)%dsig = temparray(:,3) 
                            
                                temparray = 0
                           
                                if (lCC(4)) temparray(:,1) = dhxds(iRx,iTx,iFq)%dsig
                                if (lCC(5)) temparray(:,2) = dhyds(iRx,iTx,iFq)%dsig
                                if (lCC(6)) temparray(:,3) = dhzds(iRx,iTx,iFq)%dsig
                            
                                temparray  = matmul(temparray,RotR)
                 
                                if (lCC(4))  dhxds(iRx,iTx,iFq)%dsig = temparray(:,1) 
                                if (lCC(5))  dhyds(iRx,iTx,iFq)%dsig = temparray(:,2) 
                                if (lCC(6))  dhzds(iRx,iTx,iFq)%dsig = temparray(:,3) 
                       
                         
                        endif ! lInversion  
                    
                    endif !lDataMask
                    
                enddo ! k
            enddo ! j
        endif ! rotated
    enddo ! i 
    
    
    if ( allocated(temparray)) deallocate(temparray)
                         
    
    end subroutine rotateToRxOrientation 
            
 
    end module mare2dem_worker   

