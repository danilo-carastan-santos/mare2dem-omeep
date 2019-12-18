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
!=============================================================================================================== dataTransformations
!==================================================================================================================================!
!
! Module for various data transformations
!
    module dataTransformations
    
    use Occam, only: RealPrec
    use EM_constants
    use mare2dem_output, only : Zprec
    implicit none
    
    
    contains 

!
! MT Apparent resistivity
!
        pure real(RealPrec) function getAppRes(z,freq)
        
          implicit none
          complex(Zprec), intent(in)   :: z    ! impedance Z, not C
          real(4), intent(in)      :: freq ! linear frequency
        
          
          getAppRes = abs(z)**2  / (2d0*pi*freq*mu0)
          
        end function

        elemental real(RealPrec) function getAppResDeriv(z,dz,freq)
        
          implicit none
          complex(Zprec), intent(in)   :: z    ! impedance Z, not C
          complex(Zprec), intent(in)   :: dz   ! derivative of impedance Z 
          real(4), intent(in)      :: freq ! linear frequency
          
          getAppResDeriv = 2d0 /  (2d0*pi*freq*mu0) * (dble(z)*dble(dz) + aimag(z)*aimag(dz))

        end function

!
! log10 MT Apparent Resistivity
!
        elemental real(RealPrec) function getlog10AppResDeriv(z,dz)
        
          implicit none
          complex(Zprec), intent(in)   :: z    ! impedance Z, not C
          complex(Zprec), intent(in)   :: dz   ! derivative of impedance Z 
 
          getlog10AppResDeriv = ( 2d0* log10(exp(1d0)) / (dble(z)*dble(z) + aimag(z)*aimag(z)) )* &
                                & (dble(z)*dble(dz) + aimag(z)*aimag(dz))
        
        end function
        
!
! Phase with unwrapping:
!  
        pure real(RealPrec) function getPhase(z,dataphi)
        
        ! JH hacky method of getting around lack of modulo subtraction in Occam.
        
          implicit none
          complex(Zprec), intent(in)        :: z
          real(RealPrec), intent(in)    :: dataphi
          real(4)                       :: tmp 
    
          tmp = atan2( aimag(z), real(z) )*rad2deg  ! phase in degrees
          
          getPhase = tmp
          if( abs(tmp + 360d0 - dataphi) .le. abs(tmp - dataphi) ) getPhase = (tmp + 360d0)
          if( abs(tmp - 360d0 - dataphi) .le. abs(tmp - dataphi) ) getPhase = (tmp - 360d0)
          
    
        end function
        
        elemental real(RealPrec) function getPhaseDeriv(z,dz) 
        
          implicit none
          complex(Zprec),intent(in) :: z, dz
          getPhaseDeriv = rad2deg * (dble(z)*aimag(dz) - aimag(z)*dble(dz) ) / (dble(z)*dble(z) + aimag(z)*aimag(z) )    
          
        end function          
         
!
! Absolute value:
!
        elemental real(RealPrec) function absDeriv(phi,dphi) 
        
          implicit none
          complex(Zprec),intent(in) :: phi, dphi
    
          absDeriv = ( dble(phi)*dble(dphi) + aimag(phi)*aimag(dphi))  / sqrt( dble(phi)*dble(phi) + aimag(phi)*aimag(phi)) 
          
        end function

!
! Log10 Absolute value:
!
        elemental real(RealPrec) function log10absDeriv(phi,dphi) 
        
          implicit none
          complex(Zprec),intent(in) :: phi, dphi
    
          log10absDeriv = (dble(phi)*dble(dphi) + aimag(phi)*aimag(dphi))  / &
                        & ( log(10d0)*(dble(phi)*dble(phi) + aimag(phi)*aimag(phi)) )
          ! log10absDeriv = absDeriv(phi,dphi)  / (log(10.)*abs(phi))
 
 
        end function
          
!
! Polarization ellipse parameters:
!

        pure real(RealPrec) function getPE(e1,e2,comp)
        
      ! Taken from Steve's old Occam routines - does not match  Smith and Ward 
      ! but probably works OK.
        implicit none
        complex(Zprec), intent(in)   :: e1    ! E/B-field channel 1
        complex(Zprec), intent(in)   :: e2    ! E/B-field channel 2
        character(4), intent(in) :: comp  ! flag  'pmax' or 'pmin' for output 
        
      ! Local:
        real(4)                  :: x1, x2, y1, y2
        real(4)                  :: a, b, phi, s, c, p1, p2, pmin, pmax


      ! bust up the complex components of the electric field:
        x1 = real(e1)
        y1 = aimag(e1)
        x2 = real(e2)
        y2 = aimag(e2)
      
      ! find the critical angle
        a   = -(x1*y1 + x2*y2)
        b   = 0.5*(y1*y1 + y2*y2 -x1*x1 - x2*x2)
        phi = atan2(a,b)/2.
      
      ! find the two critical axes of the pe 
        s  = sin(phi)
        c  = cos(phi)
        p1 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p1 = sqrt(abs(p1))
        
        s  = sin(phi+pi/2.d0)
        c  = cos(phi+pi/2.d0)
        p2 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p2 = sqrt(abs(p2))

        pmax = max(p1,p2)
        pmin = min(p1,p2)
        
        select case (comp)
        
          case('pmin')
            getPE = pmin
            
          case('pmax')
            getPE = pmax
            
        end select    

        end function  
        
        elemental real(RealPrec) function getPEDeriv(e1,e2, de1, de2,comp)
      ! Taken from Steve's old Occam routines - does not match Smith and Ward 
      ! but probably works OK
        implicit none
        complex(Zprec), intent(in)   :: e1    ! E/B-field channel 1
        complex(Zprec), intent(in)   :: e2    ! E/B-field channel 2
        complex(Zprec), intent(in)   :: de1   ! E/B-field deriv channel 1
        complex(Zprec), intent(in)   :: de2   ! E/B-field deriv channel 2
        character(4), intent(in) :: comp  ! flag  'pmax' or 'pmin' for output 
        
        ! Local:         
        real(4)                  :: x1, x2, y1, y2, dx1, dx2, dy1, dy2
        real(4)                  :: a, b, phi, s, c, p1, d1, p2, d2
        
      ! bust up the complex components of the electric field:
        x1  = real(e1)
        y1  = aimag(e1)
        x2  = real(e2)
        y2  = aimag(e2)
        dx1 = real(de1)
        dy1 = aimag(de1)
        dx2 = real(de2)
        dy2 = aimag(de2)
      
      ! find the critical angle
        a = -(x1*y1 + x2*y2)
        b = 0.5*(y1*y1 + y2*y2 -x1*x1 - x2*x2)
        phi = atan2(a,b)/2.
      
      ! find the two critical axes of the pe and derivatives
        s = sin(phi)
        c = cos(phi)
        p1 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p1 = sqrt(abs(p1))
        d1 = c*c*(x1*dx1+x2*dx2) + s*s*(y1*dy1+y2*dy2) + s*c*(x1*dy1+dx1*y1+x2*dy2+dx2*y2)
        d1 = d1/p1
        
        s = sin(phi+pi/2.d0)
        c = cos(phi+pi/2.d0)
        p2 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p2 = sqrt(abs(p2))
        d2 = c*c*(x1*dx1+x2*dx2) + s*s*(y1*dy1+y2*dy2) + s*c*(x1*dy1+dx1*y1+x2*dy2+dx2*y2)
        d2 = d2/p2


        if (p1 >= p2) then
        
            select case (comp)
            
                case('pmin')
                getPEDeriv = d2
                
                case('pmax')
                getPEDeriv = d1
            
            end select
        
        else ! ( p1 < p2)
        
            select case (comp)
            
                case('pmin')
                getPEDeriv = d1
                
                case('pmax')
                getPEDeriv = d2
            
            end select
        
        end if     
        
        end function  

        
    end module dataTransformations

!==================================================================================================================================! 
!======================================================================================================================== computeFwd
!==================================================================================================================================!  
    subroutine computeFwd( bDoPartials, currentMod )
!
! Routine to compute the 2D forward response and model Jacobian matrix
! This is called from Occam.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    use Occam        ! Passes out dm and wj here
    use em2dkx_mod      ! for sigparams 
    use mare2dem_global
    
    implicit none

!
! Arguments
!
    logical, intent(in)        :: bDoPartials           ! logical flag for computing Jacobians
    real(RealPrec), intent(in) :: currentMod(nFree)   ! input model parameters

!
! Local variables:
!
    integer :: i,j, ict
 
!
! Copy the model parameters (log10(resistivity)) in currentMod into the local arrays used by the fwd code:
!
    ict = 0
    do i=1,nRegions
        do j = 1,nRhoPerRegion    
            if (iFreeParam((i-1)*nRhoPerRegion + j) > 0 ) then
                ict = ict+1
                sigParams((i-1)*nRhoPerRegion + j) = 1d0/(10d0**currentMod(ict))
            endif 
        enddo        
    enddo        
 
! 
! Compute responses and sensitivities for 2D CSEM and MT data:
! 
    call computeFwd_MARE2DEM( bDoPartials,currentMod )
     
!
! If inversion update data weights for joint inversion (d_wt(:) = 1 if only CSEM or MT data present)
!
    if (bDoPartials) call compute_DataWeights
 
    
    end subroutine computeFwd

!==================================================================================================================================! 
!============================================================================================================== computeFwd_MARE2DEM
!==================================================================================================================================!   
    subroutine computeFwd_MARE2DEM( bDoPartials, currentMod )
!
! Routine to compute the 2D MT and CSEM forward response and model Jacobian matrix
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
! 
!------------------------------------------------------------------------------

    use Occam
 
    use mare2dem_global
    use dataTransformations
    use mare2dem_io
    
    implicit none

!
! Arguments
!
   logical, intent(in)        :: bDoPartials         ! logical flag for computing Jacobians  
   real(RealPrec), intent(in) :: currentMod(nFree)   ! input model parameters
 
!
! Initialize:
!  
    linversion = bDoPartials
 
!
!  Copy across data parameters to MARE2D parameters: 
!     
    if (lPrintDebug) write(*,*) 'call allocate_computeFwd_MARE2DEM...'
    call allocate_computeFwd_MARE2DEM(bDoPartials)  

!
! Call the MARE2DEM forward kernel:
!
    if (lPrintDebug) write(*,*) 'call mpi_mare2dem...'
    call mpi_mare2dem
    
!
! If computing Jacobians, convert df/dsigma to df/dlog10(rho):
!   
    if ( (bDoPartials).and.(lPrintDebug) ) write(*,*) 'call scale_derivative...'
    if (bDoPartials) call scale_derivatives(currentMod)
 
!
! Convert H to B:
!
    if (lPrintDebug) write(*,*) 'convert_H_to_B...'
    call convert_H_to_B(bDoPartials)
 
!
! Finally, extract the CSEM and MT responses in the correct data formats:
!
    if (lPrintDebug) write(*,*) 'convert_to_data_format...'
    call convert_to_data_format(bDoPartials)
 
!
! Save CSEM and MT fields to files if lFwdFields:
!
    if (lFwdFields) call writeFields_CSEM(currentIteration)
    if (lFwdFields) call writeFields_MT(currentIteration)    
          
!
! Deallocate solution arrays:
!   
    call deallocate_computeFwd_MARE2DEM
    
!
! Compute the MT static shifts, if requested:
!
    call compute_MT_staticShifts
       
      
    end subroutine computeFwd_MARE2DEM
!==================================================================================================================================! 
!=============================================================================================================== compute_DataWeights
!==================================================================================================================================!      
 
    subroutine compute_DataWeights 
    
    use Occam
    use mare2dem_global
    use mare2dem_input_data_params
    use mare2dem_io
    
    implicit none
     
    integer         :: i
    character(32)   :: cNum
    real(RealPrec)  :: weightMT, weightCSEM   , chiMT, chiCS
 
 
 !
 ! Modify data weights if joint CSEM/MT inversion:
 !   
 
    d_wt = 1.0
    
    if ( (nFreqCSEM > 0 ) .and. (nFreqMT > 0 ) ) then  ! joint inversion, compute weights:       

!
! First get the Chi^2 misfit:
!           
        chiMT  = 0
        chiCS = 0
        do i = 1,nd
            if (dp(i,1) < 100)     then       
                chiCS = chiCS + ( (d(i) - dm(i)) /sd(i))**2
            else
                chiMT = chiMT + ( (d(i) - dm(i)) /sd(i))**2
            endif 
        enddo

! Unity weights:
!        weightMT   =  1d0 
!        weightCSEM =  1d0   
!        write(*,*) ' '
!        write(*,*) 'Joint inversion: not using any CSEM nor MT data weighting'
!        write(*,*) ' '
!  

! Normalized Chi^2:      
!        weightMT   =  1d0/sqrt(real(nMT))    
!        weightCSEM =  1d0 /sqrt(real(nCSEM))  
!        write(*,*) ' '
!        write(*,*) 'Using joint inversion data weights based on normalized Chi^2 '
!        write(*,*) ' '
! 
        
! Normalized Chi^2 with misfit balancing:         
        weightMT   =  1d0/sqrt(real(nMT))*sqrt(chiMT/nMT) 
        weightCSEM =  1d0/sqrt(real(nCSEM))*sqrt(chiCS/nCSEM)               
        write(*,*) ' '
        write(*,'(a)') 'Using joint inversion data weights based on normalized Chi^2 and misfit balancing'
        write(*,*) ' '
 
    
        !
        ! Insert appropriate weights into d_wt vector:
        !  
        write(cnum,'(g12.3)') weightCSEM
        call print2col32(' CSEM Relative Weights: ',cnum,6)
        write(cnum,'(g12.3)') weightMT
        call print2col32(' MT   Relative Weights: ',cnum,6)
        write(*,*) ' '
        do i = 1,nd
            select case (  dp(i,1) ) ! what data type is this?
        
            case (1:99) 
        
                d_wt(i) = weightCSEM
        
            case(100:199)
        
                d_wt(i) = weightMT              
        
            end select

        enddo    
    
    
    endif
    
    end subroutine compute_DataWeights 
    
!==================================================================================================================================! 
!================================================================================================================ scale_derivatives
!==================================================================================================================================!    

    subroutine scale_derivatives(currentMod)

    use Occam, only : RealPrec
    
    use mare2dem_global
    use mare2dem_input_data_params
    use mare2dem_output
    
    implicit none
     
    real(RealPrec), intent(in) :: currentMod(nFree)   ! input model parameters
   
    integer                     :: iTx,iFq,iRx 

    real(8), dimension(:), allocatable :: temparray
            
 
    allocate(temparray(nFree))
    temparray = (-1.d0*dlog(10.d0)/(10**currentMod) ) 
 
    !
    ! CSEM:
    !
    do iRx = 1,nRxCSEM            
        do iTx = 1,nTxCSEM
            do iFq = 1,nfreqCSEM
 
                 if ( lDataMaskCSEM(iRx,iFq,iTx)) then 
 
                    ! Convert 
                    if (lCompCSEM(1,iRx)) dex(iRx,iTx,iFq)%dsig = dex(iRx,iTx,iFq)%dsig*temparray
                    if (lCompCSEM(2,iRx)) dey(iRx,iTx,iFq)%dsig = dey(iRx,iTx,iFq)%dsig*temparray
                    if (lCompCSEM(3,iRx)) dez(iRx,iTx,iFq)%dsig = dez(iRx,iTx,iFq)%dsig*temparray 
               
                    ! Convert B's
                    if (lCompCSEM(4,iRx)) dhx(iRx,iTx,iFq)%dsig = dhx(iRx,iTx,iFq)%dsig*temparray
                    if (lCompCSEM(5,iRx)) dhy(iRx,iTx,iFq)%dsig = dhy(iRx,iTx,iFq)%dsig*temparray 
                    if (lCompCSEM(6,iRx)) dhz(iRx,iTx,iFq)%dsig = dhz(iRx,iTx,iFq)%dsig*temparray
                    
                endif
            enddo
        enddo
    enddo 
    
    !
    ! MT: 
    !
    do iRx = 1,nRxMT            
        do iFq = 1,nfreqMT
            if ( lDataMaskMT(iRx,iFq)) then 
                if (allocated(dex_mt_dsig)) dex_mt_dsig(iRx,iFq,:) = dex_mt_dsig(iRx,iFq,:)*temparray
                if (allocated(dey_mt_dsig)) dey_mt_dsig(iRx,iFq,:) = dey_mt_dsig(iRx,iFq,:)*temparray
                if (allocated(dez_mt_dsig)) dez_mt_dsig(iRx,iFq,:) = dez_mt_dsig(iRx,iFq,:)*temparray
                if (allocated(dhx_mt_dsig)) dhx_mt_dsig(iRx,iFq,:) = dhx_mt_dsig(iRx,iFq,:)*temparray 
                if (allocated(dhy_mt_dsig)) dhy_mt_dsig(iRx,iFq,:) = dhy_mt_dsig(iRx,iFq,:)*temparray
                if (allocated(dhz_mt_dsig)) dhz_mt_dsig(iRx,iFq,:) = dhz_mt_dsig(iRx,iFq,:)*temparray    
            endif
        enddo
    enddo
 
    deallocate (temparray)    
    end subroutine scale_derivatives 
    
!==================================================================================================================================! 
!==================================================================================================================== convert_H_to_B
!==================================================================================================================================!
    subroutine convert_H_to_B( bDoPartials )    


    use EM_constants
    use mare2dem_global
    use mare2dem_input_data_params
    use mare2dem_output
    
    implicit none
     
    logical, intent(in)         :: bDoPartials         ! logical flag for computing Jacobians   
 
    integer                     :: iTx,iFq,iRx 
            
    hx = hx*mu0
    hy = hy*mu0
    hz = hz*mu0
    if (bDoPartials) then
        do iRx = 1,nRxCSEM            
            do iTx = 1,nTxCSEM
                do iFq = 1,nfreqCSEM
                    if (lCompCSEM(4,iRx)) then
                        if (allocated(dhx(iRx,iTx,iFq)%dsig)) dhx(iRx,iTx,iFq)%dsig = dhx(iRx,iTx,iFq)%dsig*mu0
                    endif
                    if (lCompCSEM(5,iRx)) then
                        if (allocated(dhy(iRx,iTx,iFq)%dsig)) dhy(iRx,iTx,iFq)%dsig = dhy(iRx,iTx,iFq)%dsig*mu0
                    endif
                    if (lCompCSEM(6,iRx)) then
                        if (allocated(dhz(iRx,iTx,iFq)%dsig)) dhz(iRx,iTx,iFq)%dsig = dhz(iRx,iTx,iFq)%dsig*mu0
                    endif
                enddo
            enddo
        enddo 
    endif
    
    end subroutine convert_H_to_B
    
!==================================================================================================================================! 
!===================================================================================================== allocate_computeFwd_MARE2DEM
!==================================================================================================================================!    
    subroutine allocate_computeFwd_MARE2DEM(bDoPartials )     
!
! Allocates and sets the I/O arrays needed by the forward kernel, and also handles some of the
! finite length dipole quadrature point setup
!
    
    use Occam, only : nd, dp
    use mare2dem_global  
    use mare2dem_output
    use mare2dem_input_data_params
    
    implicit none
     
    logical, intent(in)         :: bDoPartials         ! logical flag for computing Jacobians   
 
    integer                     :: i,iTx,iFq,iRx 
   
     
    !
    ! Allocate arrays for CSEM solutions at receivers:
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
    if (bDoPartials) then 
               

        allocate( dex(nRxCSEM,nTxCSEM,nfreqCSEM) )   
        allocate( dey(nRxCSEM,nTxCSEM,nfreqCSEM) ) 
        allocate( dez(nRxCSEM,nTxCSEM,nfreqCSEM) ) 
        allocate( dhx(nRxCSEM,nTxCSEM,nfreqCSEM) ) 
        allocate( dhy(nRxCSEM,nTxCSEM,nfreqCSEM) ) 
        allocate( dhz(nRxCSEM,nTxCSEM,nfreqCSEM) )   
        
                
        do i = 1,nd
        
            iFq =  dp(i,2)
            iTx =  dp(i,3)
            iRx =  dp(i,4)
                  
            ! Allocate Frechet derivative arrays only where the data is present:                 
            if (dp(i,1) < 100) then
                            
                if (lCompCSEM(1,iRx)) then
                    if (.not.allocated(dex(iRx,iTx,iFq)%dsig)) then
                        allocate( dex(iRx,iTx,iFq)%dsig(nFree) ) 
                        dex(iRx,iTx,iFq)%dsig = 0
                    endif
                endif
                if (lCompCSEM(2,iRx)) then 
                    if (.not.allocated(dey(iRx,iTx,iFq)%dsig))  then
                        allocate( dey(iRx,iTx,iFq)%dsig(nFree) ) 
                        dey(iRx,iTx,iFq)%dsig = 0 
                    endif
                endif  
                if (lCompCSEM(3,iRx)) then
                    if (.not.allocated(dez(iRx,iTx,iFq)%dsig))  then
                        allocate( dez(iRx,iTx,iFq)%dsig(nFree) ) 
                        dez(iRx,iTx,iFq)%dsig = 0 
                    endif 
                endif
                if (lCompCSEM(4,iRx)) then
                    if (.not.allocated(dhx(iRx,iTx,iFq)%dsig))  then
                        allocate( dhx(iRx,iTx,iFq)%dsig(nFree) ) 
                        dhx(iRx,iTx,iFq)%dsig = 0
                    endif 
                endif
                if (lCompCSEM(5,iRx)) then
                    if (.not.allocated(dhy(iRx,iTx,iFq)%dsig))  then
                        allocate( dhy(iRx,iTx,iFq)%dsig(nFree) ) 
                        dhy(iRx,iTx,iFq)%dsig = 0 
                    endif
                endif
                if (lCompCSEM(6,iRx)) then
                    if (.not.allocated(dhz(iRx,iTx,iFq)%dsig))  then
                        allocate( dhz(iRx,iTx,iFq)%dsig(nFree) )  
                        dhz(iRx,iTx,iFq)%dsig = 0
                    endif 
                endif
                 
            endif 
           
        enddo
               
    endif

    !
    ! Specify CSEM data mask array:
    !
    if (nTxCSEM > 0 ) then
        
        allocate(lDataMaskCSEM(nRxCSEM,nFreqCSEM,nTxCSEM))  ! kwk debug: this could be generated at bottom of readData
        
        lDataMaskCSEM = .false.
        
        do i = 1,nd
            iRx = dp(i,4)  
            iTx = dp(i,3) 
            iFq = dp(i,2)            
            if ( dp(i,1) < 100 ) lDataMaskCSEM(iRx,iFq,iTx) = .true. 
        enddo

        if (lFwdFields) lDataMaskCSEM = .true.    
        
    endif        
    
    
    
    !
    ! MT:
    ! 
 
    allocate( ex_mt(nRxMT,nfreqMT),ey_mt(nRxMT,nfreqMT), ez_mt(nRxMT,nfreqMT) )
    allocate( hx_mt(nRxMT,nfreqMT),hy_mt(nRxMT,nfreqMT), hz_mt(nRxMT,nfreqMT) )
    if (bDoPartials) then           
        if (any(lCompMT(1,:))) allocate( dex_mt_dsig(nRxMT,nfreqMT,nFree) )
        if (any(lCompMT(5,:))) allocate( dhy_mt_dsig(nRxMT,nfreqMT,nFree) )   
        if (any(lCompMT(6,:))) allocate( dhz_mt_dsig(nRxMT,nfreqMT,nFree) )
        if (any(lCompMT(1,:))) allocate(   dzte_dsig(nFree) )
        if (any(lCompMT(6,:))) allocate(   dhtipper_dsig(nFree) )
        
        if (any(lCompMT(4,:))) allocate( dhx_mt_dsig(nRxMT,nfreqMT,nFree) )
        if (any(lCompMT(2,:))) allocate( dey_mt_dsig(nRxMT,nfreqMT,nFree) )
        if (any(lCompMT(3,:))) allocate( dez_mt_dsig(nRxMT,nfreqMT,nFree) )
        if (any(lCompMT(2,:))) allocate(   dztm_dsig(nFree) )            
    endif

    
    !
    ! MT data masks:
    !
    if (nFreqMT > 0 ) then
    
        allocate(lDataMaskMT(nRxMT,nFreqMT))  
        
        !         Type        Freq#        Tx#           Rx#   
        lDataMaskMT = .false.
        do i = 1,nd
            iRx = dp(i,4)
            iTx = dp(i,3)
            iFq = dp(i,2)
            
            if ( dp(i,1) > 100 ) then ! only consider MT data here:
                 lDataMaskMT(iRx,iFq) = .true. 
                 if (iTx > 0 ) lDataMaskMT(iTx,iFq) = .true.  ! for magnetics at receiver iTx
            endif
        enddo
        
        if (lFwdFields) lDataMaskMT = .true.
        
    endif
    
    end subroutine allocate_computeFwd_MARE2DEM
  
!==================================================================================================================================! 
!============================================================================================================ convert_to_data_format
!==================================================================================================================================! 
    subroutine convert_to_data_format( bDoPartials )
    
  
    use Occam
    use mare2dem_global
    use mare2dem_output
    use mare2dem_input_data_params    
    use dataTransformations

    implicit none
     
    include 'EM_parameters.inc'
           
    logical, intent(in)         :: bDoPartials         ! logical flag for computing Jacobians   
 
    integer                     :: i,iTx,iFq,iRx

    
!
! This is a giant block that handles many different data types:
! 
    do i = 1,nd
    
        iFq =  dp(i,2)
        iTx =  dp(i,3)
        iRx =  dp(i,4)         
            

!
! Special for MT data:
!
!       magnetic fields come from iTx. Normally iTx = iRx for standard MT responses. If <= 0 then iTx := iRx
!       electric fields come from iRx. Use iTx /=iRx for hydrid MT or differential tilts between E and H 
!
        if ( dp(i,1) > 100 ) then
    
            if (iTx <= 0) iTx = iRx  ! if iTx not specified, make it the current receiver 
 
            zte     = ex_mt(iRx,iFq) / hy_mt(iTx,iFq)
            ztm     = ey_mt(iRx,iFq) / hx_mt(iTx,iFq)
            htipper = hz_mt(iTx,iFq) / hy_mt(iTx,iFq)
            
                       
        endif
    
        select case (  dp(i,1) ) ! what data type is this?

!
!  CSEM data:                                     
!                              
        ! Electric: Real and Imaginary Data:
                
                case (indRealEx)  
                  dm(i) = dble(ex(iRx,iTx,iFq))
                case (indImagEx)  
                  dm(i) = aimag(ex(iRx,iTx,iFq))
                case (indRealEy)  
                  dm(i) = dble(ey(iRx,iTx,iFq))
                case (indImagEy)  
                  dm(i) = aimag(ey(iRx,iTx,iFq))
                case (indRealEz)  
                  dm(i) = dble(ez(iRx,iTx,iFq))  
                case (indImagEz)  
                  dm(i) = aimag(ez(iRx,iTx,iFq)) 
                    
        ! Electric: Amplitude and Phase Data:
                           
                case (indAmpEx)   
                  dm(i) = abs(ex(iRx,iTx,iFq))    !     |Ex|
                case (indAmpEy)   
                  dm(i) = abs(ey(iRx,iTx,iFq))    !     |Ey|
                case (indAmpEz)   
                  dm(i) = abs(ez(iRx,iTx,iFq))    !     |Ez|
                case (indPhsEx) 
                  dm(i) = getPhase(ex(iRx,iTx,iFq),d(i))
                case (indPhsEy) 
                  dm(i) = getPhase(ey(iRx,iTx,iFq),d(i))
                case (indPhsEz) 
                  dm(i) = getPhase(ez(iRx,iTx,iFq),d(i))

        ! Electric: Log10 Amplitude Data:
                           
                case (indLog10AmpEx)   
                  dm(i) = log10(abs(ex(iRx,iTx,iFq)))    ! log10|Ex|
                case (indLog10AmpEy)   
                  dm(i) = log10(abs(ey(iRx,iTx,iFq)))    ! log10|Ey|
                case (indLog10AmpEz)   
                  dm(i) = log10(abs(ez(iRx,iTx,iFq)))    ! log10|Ez|                      
                  
        ! Magnetic: Real and Imaginary Data:

                case (indRealBx) 
                  dm(i) = dble(hx(iRx,iTx,iFq))   !        real(Bx)
                case (indImagBx) 
                  dm(i) = aimag(hx(iRx,iTx,iFq))  !        imag(Bx)
                case (indRealBy) 
                  dm(i) = dble(hy(iRx,iTx,iFq))   !        real(By)
                case (indImagBy) 
                  dm(i) = aimag(hy(iRx,iTx,iFq))  !        imag(By)
                case (indRealBz) 
                  dm(i) = dble(hz(iRx,iTx,iFq))   !        real(Bz)
                case (indImagBz) 
                  dm(i) = aimag(hz(iRx,iTx,iFq))  !        imag(Bz)
                    
        ! Magnetic: Amplitude and Phase Data:
                    
                case (indAmpBx) 
                  dm(i) = abs(hx(iRx,iTx,iFq))   !     |Bx|
                case (indAmpBy) 
                  dm(i) = abs(hy(iRx,iTx,iFq))   !     |By|
                case (indAmpBz) 
                  dm(i) = abs(hz(iRx,iTx,iFq))   !     |Bz|
                case (indPhsBx) 
                  dm(i) = getPhase(hx(iRx,iTx,iFq),(d(i)))
                case (indPhsBy) 
                  dm(i) = getPhase(hy(iRx,iTx,iFq),(d(i)))
                case (indPhsBz) 
                  dm(i) = getPhase(hz(iRx,iTx,iFq),(d(i)))
                  
        ! Magnetic: Log10 Amplitude Data:
                           
                case (indLog10AmpBx)   
                  dm(i) = log10(abs(hx(iRx,iTx,iFq)))    ! log10|Bx|
                case (indLog10AmpBy)   
                  dm(i) = log10(abs(hy(iRx,iTx,iFq)))    ! log10|By|
                case (indLog10AmpBz)   
                  dm(i) = log10(abs(hz(iRx,iTx,iFq)))    ! log10|Bz|                      
                  
                                 
         ! Polarization Ellipse Parameters:
         
                case (iPEmax)
                  dm(i) = getPE(ex(iRx,iTx,iFq),ey(iRx,iTx,iFq),'pmax') ! E-field max
                case (iPEmin)
                  dm(i) = getPE(ex(iRx,iTx,iFq),ey(iRx,iTx,iFq),'pmin') ! E-field min
                case (iPBmax)
                  dm(i) = getPE(hx(iRx,iTx,iFq),hy(iRx,iTx,iFq),'pmax') ! B-field max
                case (iPBmin)
                  dm(i) = getPE(hx(iRx,iTx,iFq),hy(iRx,iTx,iFq),'pmin') ! B-field min
                  
!
!  MT data:                                     
!                     
                  
              ! TE:            
                case (indRhoZXY)                                        ! TE Apparent resistivity
                    dm(i) = getAppRes(zte,real(fTxMT(iFq)))    
                case (indlog10RhoZXY)                                   ! log10 TE Apparent resistivity
                    dm(i) = log10(getAppRes(zte,real(fTxMT(iFq))))                                                         
                case (indPhsZXY)                                        ! TE Phase
                    dm(i) = getPhase(zte,(d(i)))                  
                case (indRealZXY)                                       ! real(Zxy)
                    dm(i) = dble( zte ) 
                case (indImagZXY)                                       ! imag(Zxy)
                    dm(i) = aimag(zte )      
                                        
                !TM:    
                case (indRhoZYX)                                        ! TM Apparent resistivity
                    dm(i) = getAppRes( -ztm,real(fTxMT(iFq)))                         
                case (indlog10RhoZYX)                                   ! log10 TM Apparent resistivity
                    dm(i) = log10(getAppRes( -ztm,real(fTxMT(iFq))))     
                case (indPhsZYX)                                        ! TM Phase (moved to first quadrant from third)   
                    dm(i) = getPhase( -ztm,(d(i)))                  
                case (indRealZYX)                                       ! real(Zyx)
                    dm(i) = dble( ztm )             
                case (indImagZYX)                                       ! imag(Zyx)
                    dm(i) = aimag(ztm ) 
                
                                             
                ! Magnetic Tipper:    
                case (indRealMZY)                                       ! real(Mzy) ! TE mode Hz = Mzy*Hy
                    dm(i) = dble(  htipper )        
                case (indImagMZY)                                       ! imag(Mzy)
                    dm(i) = aimag( htipper )           
                                      
             
        end select  ! case dp(i,1)
      
      enddo  ! loop over nd

      if( bDoPartials) then  ! calculate the partial derivatives as well.
        
        do i=1,nd 
    
            iFq =  dp(i,2)
            iTx =  dp(i,3)
            iRx =  dp(i,4)        
 
          
!
! Special tweaks for MT data:
!
            if ( dp(i,1) > 100 ) then
                if (iTx <= 0) iTx = iRx  ! if iTx not specified, make it the current receiver 
        
                zte     = ex_mt(iRx,iFq) / hy_mt(iTx,iFq)
                ztm     = ey_mt(iRx,iFq) / hx_mt(iTx,iFq)
                htipper = hz_mt(iTx,iFq) / hy_mt(iTx,iFq)
        
                if (bDoPartials) then

                    if (lCompMT(1,iRx)) dzte_dsig(:)     = (dex_mt_dsig(iRx,iFq,:) - zte    *dhy_mt_dsig(iTx,iFq,:))/hy_mt(iTx,iFq)
                    if (lCompMT(6,iRx)) dhtipper_dsig(:) = (dhz_mt_dsig(iTx,iFq,:) - htipper*dhy_mt_dsig(iTx,iFq,:))/hy_mt(iTx,iFq)
                    if (lCompMT(2,iRx)) dztm_dsig(:)     = (dey_mt_dsig(iRx,iFq,:) - ztm    *dhx_mt_dsig(iTx,iFq,:))/hx_mt(iTx,iFq)

                endif
            endif
            
           select case (  dp(i,1) )  ! what data type is this?
 
!
! CSEM data:
!              
           ! Electric: Real and Imaginary Data:
                      
              case (indRealEx)  
                  wj(i, :) = dble(dex(iRx,iTx,iFq)%dsig(:))    !        real(Ex)
              case (indImagEx)  
                  wj(i, :) = aimag(dex(iRx,iTx,iFq)%dsig(:))   !        imag(Ex)
              case (indRealEy)  
                  wj(i, :) = dble(dey(iRx,iTx,iFq)%dsig(:))    !        real(Ey)
              case (indImagEy)  
                  wj(i, :) = aimag(dey(iRx,iTx,iFq)%dsig(:))   !        imag(Ey)
              case (indRealEz)  
                  wj(i, :) = dble( dez(iRx,iTx,iFq)%dsig(:))   !        real(Ez)
              case (indImagEz)  
                  wj(i, :) = aimag( dez(iRx,iTx,iFq)%dsig(:))  !        imag(Ez)
                              
           ! Electric: Amplitude and Phase Data:
                                 
              case (indAmpEx)   
                  wj(i, :) = absDeriv( ex(iRx,iTx,iFq), dex(iRx,iTx,iFq)%dsig(:) )   !     |Ex|
              case (indAmpEy)   
                  wj(i, :) = absDeriv( ey(iRx,iTx,iFq), dey(iRx,iTx,iFq)%dsig(:) )   !     |Ey|
              case (indAmpEz)   
                  wj(i, :) = absDeriv( ez(iRx,iTx,iFq), dez(iRx,iTx,iFq)%dsig(:) )   !     |Ez|
              case (indPhsEx) 
                  wj(i, :) = getPhaseDeriv( ex(iRx,iTx,iFq), dex(iRx,iTx,iFq)%dsig(:) )   !     Ex phase 
              case (indPhsEy)
                  wj(i, :) = getPhaseDeriv( ey(iRx,iTx,iFq), dey(iRx,iTx,iFq)%dsig(:) )   !     Ey phase 
              case (indPhsEz) 
                  wj(i, :) = getPhaseDeriv( ez(iRx,iTx,iFq), dez(iRx,iTx,iFq)%dsig(:) )   !     Ez phase 

        ! Electric: Log10 Amplitude Data:
                        
              case (indLog10AmpEx)   
                  wj(i, :) = log10absDeriv(ex(iRx,iTx,iFq), dex(iRx,iTx,iFq)%dsig(:))  ! log10|Ex|
              case (indLog10AmpEy)   
                  wj(i, :) = log10absDeriv(ey(iRx,iTx,iFq), dey(iRx,iTx,iFq)%dsig(:))  ! log10|Ey|
              case (indLog10AmpEz)   
                  wj(i, :) = log10absDeriv(ez(iRx,iTx,iFq), dez(iRx,iTx,iFq)%dsig(:))  ! log10|Ez|    

           ! Magnetic: Real and Imaginary Data:

              case (indRealBx) 
                  wj(i, :) = dble(dhx(iRx,iTx,iFq)%dsig(:))   !        real(Bx)
              case (indImagBx) 
                  wj(i, :) = aimag(dhx(iRx,iTx,iFq)%dsig(:))  !        imag(Bx)
              case (indRealBy) 
                  wj(i, :) = dble(dhy(iRx,iTx,iFq)%dsig(:))   !        real(By)
              case (indImagBy) 
                  wj(i, :) = aimag(dhy(iRx,iTx,iFq)%dsig(:))  !        imag(By)
              case (indRealBz) 
                  wj(i, :) = dble(dhz(iRx,iTx,iFq)%dsig(:))   !        real(Bz)
              case (indImagBz) 
                  wj(i, :) = aimag(dhz(iRx,iTx,iFq)%dsig(:))  !        imag(Bz)
                              
           ! Magnetic: Amplitude and Phase Data:
                          
              case (indAmpBx) 
                  wj(i, :) = absDeriv( hx(iRx,iTx,iFq), dhx(iRx,iTx,iFq)%dsig(:) )   !     |Bx|
              case (indAmpBy) 
                  wj(i, :) = absDeriv( hy(iRx,iTx,iFq), dhy(iRx,iTx,iFq)%dsig(:) )   !     |By|
              case (indAmpBz) 
                  wj(i, :) = absDeriv( hz(iRx,iTx,iFq), dhz(iRx,iTx,iFq)%dsig(:) )   !     |Bz|
              case (indPhsBx) 
                  wj(i, :) = getPhaseDeriv( hx(iRx,iTx,iFq), dhx(iRx,iTx,iFq)%dsig(:) )   !     Bx phase 
              case (indPhsBy) 
                  wj(i, :) = getPhaseDeriv( hy(iRx,iTx,iFq), dhy(iRx,iTx,iFq)%dsig(:) )   !     By phase 
              case (indPhsBz) 
                  wj(i, :) = getPhaseDeriv( hz(iRx,iTx,iFq), dhz(iRx,iTx,iFq)%dsig(:) )   !     Bz phase 

          ! Magnetic: Log10 Amplitude Data:
                        
              case (indLog10AmpBx)   
                  wj(i, :) = log10absDeriv(hx(iRx,iTx,iFq), dhx(iRx,iTx,iFq)%dsig(:))  ! log10|Bx|
              case (indLog10AmpBy)   
                  wj(i, :) = log10absDeriv(hy(iRx,iTx,iFq), dhy(iRx,iTx,iFq)%dsig(:))  ! log10|By|
              case (indLog10AmpBz)   
                  wj(i, :) = log10absDeriv(hz(iRx,iTx,iFq), dhz(iRx,iTx,iFq)%dsig(:))  ! log10|Bz|    

       
           ! Polarization Ellipse Parameters:   
           
              case (iPEmax)
                  wj(i, :) = getPEDeriv(ex(iRx,iTx,iFq),ey(iRx,iTx,iFq), &
                                        dex(iRx,iTx,iFq)%dsig(:),dey(iRx,iTx,iFq)%dsig(:),'pmax') ! E-field max
              case (iPEmin)
                  wj(i, :) = getPEDeriv(ex(iRx,iTx,iFq),ey(iRx,iTx,iFq), &
                                        dex(iRx,iTx,iFq)%dsig(:),dey(iRx,iTx,iFq)%dsig(:),'pmin') ! E-field min
              case (iPBmax)
                  wj(i, :) = getPEDeriv(hx(iRx,iTx,iFq),hy(iRx,iTx,iFq), &
                                        dhx(iRx,iTx,iFq)%dsig(:),dhy(iRx,iTx,iFq)%dsig(:),'pmax') ! B-field max
              case (iPBmin)
                  wj(i, :) = getPEDeriv(hx(iRx,iTx,iFq),hy(iRx,iTx,iFq), &
                                        dhx(iRx,iTx,iFq)%dsig(:),dhy(iRx,iTx,iFq)%dsig(:),'pmin') ! B-field min

!
! MT data:
!
            ! TE mode:
                case (indRhoZXY)                                        ! TE Apparent resistivity
                    wj(i, :) = getAppResDeriv(zte,dzte_dsig(:),real(fTxMT(iFq)) )          
                case (indlog10RhoZXY)                                    
                    wj(i, :) = getlog10AppResDeriv(zte,dzte_dsig(:))  
                case (indPhsZXY)                                        ! TE Phase
                    wj(i, :) = getPhaseDeriv( zte,dzte_dsig(:) )        
                case (indRealZXY)                                       ! real(Zxy)
                   wj(i, :) = real( (dzte_dsig(:)) )
                case (indImagZXY)                                       ! imag(Zxy)
                    wj(i, :) = aimag( (dzte_dsig(:)) )                           

            ! TM mode:                                       
                case (indRhoZYX)                                        ! TM Apparent resistivity
                    wj(i, :) = getAppResDeriv( -ztm,-dztm_dsig(:),real(fTxMT(iFq)) )    
                case (indlog10RhoZYX)                                        ! TM Apparent resistivity
                    wj(i, :) = getlog10AppResDeriv( -ztm,-dztm_dsig(:))                      
                case (indPhsZYX)                                        ! TM Phase    
                    wj(i, :) = getPhaseDeriv(  -ztm,-dztm_dsig(:) )    
                case (indRealZYX)                                       ! real(Zyx)
                    wj(i, :) = real( dztm_dsig(:) )     
                case (indImagZYX)                                       ! imag(Zyx)
                   wj(i, :) = aimag( dztm_dsig(:) )    
   
            ! Magnetic Tipper:    
                case (indRealMZY)                                       ! real(Mzy) ! TE mode Hz = Mzy Hy
                   wj(i, :) = real(  dhtipper_dsig(:) )         
                case (indImagMZY)                                       ! imag(Mzy)
                   wj(i, :) = aimag( dhtipper_dsig(:) )         



           end select  ! case dp(i,1)

         end do

      endif !bDoPartials
    
    end subroutine convert_to_data_format      

!==================================================================================================================================! 
!==================================================================================================== deallocate_computeFwd_MARE2DEM
!==================================================================================================================================!     
    subroutine deallocate_computeFwd_MARE2DEM
    
    use mare2dem_global
    use mare2dem_output             
    !
    ! Deallocate solution arrays:
    !
    call deallocate_mare2dem_output 
    
    if (allocated(lDataMaskCSEM)) deallocate(lDataMaskCSEM) 
    if (allocated(lDataMaskMT))   deallocate(lDataMaskMT)  
        
    
     end subroutine deallocate_computeFwd_MARE2DEM    
        
!==================================================================================================================================! 
!=========================================================================================================== compute_MT_staticShifts
!==================================================================================================================================! 
    subroutine compute_MT_staticShifts
!
! A simplified static shift estimator based solely on the size of the mean residual. 
!
! ****** USE WITH CAUTION! *****
! 
    use Occam
    use mare2dem_input_data_params
    
    implicit none

    include 'EM_parameters.inc'
    
    integer         :: iRx, i, nte, ntm
    
    real(RealPrec)  :: tesum, tmsum
    
!
! Solve for the static shift parameters at each site, if requested:
!
    do iRx = 1,nRxMT
        
        tesum     = 0
        tmsum     = 0
        pm_assoc(2*iRx-1:2*iRx) = 0
        
        if (iSolveStatic(iRx) /= 0) then
        
            nte = 0
            ntm = 0
            
            do i=1,nd

                if ( dp(i,4)  == iRx ) then
            
                    select case (  dp(i,1) )  ! what data type is this?   
                    
                    case (indRhoZXY)    ! TE Apparent resistivity
                        tesum = tesum + log10(d(i)) - log10(dm(i))         
                        nte   = nte + 1
                    case (indlog10RhoZXY)                                    
                        tesum = tesum + d(i) - dm(i)
                        nte   = nte + 1     
                    case (indRhoZYX)   ! TM Apparent resistivity
                        tmsum = tmsum + log10(d(i)) - log10(dm(i)) 
                        ntm   = ntm + 1
                    case (indlog10RhoZYX) 
                        tmsum = tmsum + d(i) - dm(i)
                        ntm   = ntm + 1                    
                    end select
                endif     
            enddo
            if (nte > 0 ) then
                tesum = tesum / nte
            endif
            if (ntm > 0 ) then
                tmsum = tmsum / ntm
            endif            
            
            ! If only TE or TM statics requested, zero the other static estimate:
            if (iSolveStatic(iRx) == 2) then ! TE only
                tmsum = 0
            elseif (iSolveStatic(iRx) == 3) then ! TM only
                tesum = 0
            endif 
            
            !write(*,*) 'tesum, tmsum: ',iRx, tesum, tmsum
             
            
            
            ! Store the static shifts:
            pm_assoc(2*iRx-1) = tesum
            pm_assoc(2*iRx  ) = tmsum    
            
            
            
            ! Finally do a second pass to apply the static shifts to the forward responses:
            do i=1,nd

                if ( dp(i,4)  == iRx ) then
            
                    select case (  dp(i,1) )  ! what data type is this?   
                    
                    case (indRhoZXY)    ! TE Apparent resistivity
                        dm(i) = dm(i)*10**tesum       
                    case (indlog10RhoZXY)                                    
                        dm(i) = dm(i) + tesum
                    case (indRhoZYX)   ! TM Apparent resistivity
                        dm(i) = dm(i)*10**tmsum 
                    case (indlog10RhoZYX) 
                        dm(i) = dm(i) + tmsum    
                                   
                    end select
                endif     
            enddo 
        endif
    enddo
 
    end subroutine compute_MT_staticShifts
    
 