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
!=================================================================================================================== spline_integ_kx
!==================================================================================================================================!         
    module spline_integ_kx
!
! Stores temporary arrays used for the spline integration along kx
!
    use SinCosFilters
     
    use mare2dem_global  ! use nwave and wavenum 
  
    implicit none

    private
        
    complex(8),dimension(:)             :: fkxbuff,coeff
    allocatable                         :: fkxbuff,coeff 
    
    integer                             :: nfiltcoef 
    real(8)                             :: filter(1001),base(1001),cosfc(1001),sinfc(1001)
    
    real(8), dimension(:), allocatable  :: ureal,uimag, yreal, yimag, y2imag, y2real, x2   
      

    
    public :: spline_integ_kx_setup
    public :: spline_integ_kx_deallocate
    public :: spline_integ_fkx
    
    contains
    
    
!==================================================================================================================================! 
!============================================================================================================= spline_integ_kx_setup
!==================================================================================================================================! 
    subroutine spline_integ_kx_setup
    
!
! Allocate arrays passed in module spline_integ
!
    allocate(coeff(nwave),fkxbuff(nwave) )    
    

!
! Always use 601 point filter here to make sure final integration is robust
!
        call init_SinCosFilters ! this is done in mpi_slave_receive_invFT

        
        nfiltcoef       = 601               ! these filters seem to be the most robust.
        base(1:601)     = kk_ct_base_601
        cosfc(1:601)    = kk_ct_cos_601
        sinfc(1:601)    = kk_ct_sin_601     

!        nfiltcoef       = 201
!        base(1:201)     = kk_ct_base_201
!        cosfc(1:201)    = kk_ct_cos_201
!        sinfc(1:201)    = kk_ct_sin_201     

!        nfiltcoef       = 901
!        base(1:901)     = kk_ct_base_901
!        cosfc(1:901)    = kk_ct_cos_901
!        sinfc(1:901)    = kk_ct_sin_901   
        
!        nfiltcoef       = 241
!        base(1:241)     = kk_ct_base_241
!        cosfc(1:241)    = kk_ct_cos_241
!        sinfc(1:241)    = kk_ct_sin_241     
        
!        nfiltcoef       = 81
!        base(1:81)     = kk_ct_base_81
!        cosfc(1:81)    = kk_ct_cos_81
!        sinfc(1:81)    = kk_ct_sin_81     

    
    allocate (x2(nwave), ureal(nwave),uimag(nwave),yreal(nwave), yimag(nwave), y2imag(nwave), y2real(nwave))
      
    x2 = log10(wavenum)
               
    end subroutine spline_integ_kx_setup
    
!==================================================================================================================================! 
!======================================================================================================== spline_integ_kx_deallocate
!==================================================================================================================================! 
    subroutine spline_integ_kx_deallocate
    
 
    deallocate(coeff,fkxbuff)    
    deallocate (x2, ureal,uimag,yreal, yimag, y2imag, y2real)
          
    end subroutine spline_integ_kx_deallocate
          
          
!==================================================================================================================================! 
!================================================================================================================== spline_integ_fkx
!==================================================================================================================================!              
    subroutine spline_integ_fkx(xr, fkx, nord, result)
! 
!
! Calculates the integral of constant*f(kx)*w(x*kx) from
! kx = 0 to infinity where w() is cos() if nord = 0 or 
! sin() if nord = 1.  constant is 1/pi or sqrt(-1)/pi
! for cosine and sine transforms.  Uses spline interpolation
! to fill in f(kx) at values between discrete solutions.
!   
! input fkx can be exkx, hxkx, hykx, etc...
!
! Uses the new sine and cosine digital filters I made
! using a method similar to the F.N. Kong GP (2007) paper.
    
 
!    
! Input:
!
    real(8), intent(in)                      :: xr
    integer, intent(in)                      :: nord
    complex(4), dimension(nwave), intent(in) :: fkx
    complex(4), intent(out)                  :: result
    
!
! Local:
!
    integer                       :: i, xsgn
    real(8)                       :: r, kx
    complex(8)                    :: c1, f, sumf
 
!
! call spline_integ_setup() before calling this routine
! 
 

!
! Set up the constant:
!   
    xsgn = 1 
    if (nord.eq.0) then  ! even function, use cos integral
        c1 = cmplx(1.0d0/pi,0.0d0)
        xsgn = 1    ! no sign reversal for even (cos) integral
        filter(1:nfiltcoef) = cosfc(1:nfiltcoef)
    elseif (nord.eq.1) then ! odd function, use sine integral
        c1 = cmplx(0.d0,1.0d0/pi)
        filter(1:nfiltcoef) = sinfc(1:nfiltcoef)
        if (xr.lt.0) xsgn = -1 ! sign reversal for odd integral
    endif

           
!    
! Copy exkx to a double precision buffer:
!
    fkxbuff(1:nwave) = fkx(1:nwave)
!
! Get some spline coefficients for the exkx points
!
    call spline_i(fkxbuff,nwave,coeff)
 
!
! Make sure the xr location of the receiver isn't too small
!
    r = abs(xr)
    if(r.lt.1d0) then  ! limit small ranges
        r = 1d0
    endif

! 
! Perform the cosine (or sine) transform using the digital filter method
!    
    sumf = 0d0
    
    do i = 1,nfiltcoef
 
!
! Get kx value from base
!
        kx = base(i)/r
!
! Get interpolated field in wavenumber domain: 
! 
        f = splint_i(log10(kx))
!
! Apply digital filter:
!
        sumf = sumf + f*filter(i) ! Note that this loop could be more sophisticated, for the case of small fields, 
                                      ! but it works for now. 
          
    enddo
!
! Normalize by r
!
    sumf = sumf/r 
  
   
!
! Apply coefficients:
!          
    result=xsgn*c1*sumf   
 
    end subroutine spline_integ_fkx
       
!==================================================================================================================================! 
!========================================================================================================================== spline_i
!==================================================================================================================================!          
     subroutine spline_i(y,n,y2)
!
! Subroutine to generate second derivative of spline interpolating function
! Modified from Numerical Recipes section 3.3.  Use this before calling splint.
!
! x (real, length n), y (complex, length n) are the precomputed values of 
! the function to be splined.  y2 (complex, length n) is the second derivative 
! of the interpolating function
!
!
! For EM transforms, use x = log10(kx)  
! 
    
   
    integer, intent(in)     :: n
    complex(8), intent(in)  :: y(n)  
    complex(8), intent(out) :: y2(n)    
    
    ! Local:
    
    integer i,k
    real(8) sig, dx1,dx2,dx3,dyr1,dyr2,dyi1,dyi2  
    real(8) preal,pimag  
    
 
!
! Compute real parts:
!
    y2real = 0d0    
    yreal  = dble(y)
    ureal  = 0d0
    
    y2imag = 0d0    
    yimag  = aimag(y)   
    uimag  = 0d0

       
    do i=2,n-1
    
        dx1 = x2(i+1) - x2(i) 
        dx2 = x2(i)   - x2(i-1) 
        dx3 = x2(i+1) - x2(i-1)
        sig = dx2 / dx3
        
        preal = sig*y2real(i-1) + 2d0
        pimag = sig*y2imag(i-1) + 2d0
        
        y2real(i) = (sig-1.0d0)/preal
        y2imag(i) = (sig-1.0d0)/pimag
         
        dyr1 = yreal(i+1) - yreal(i)
        dyr2 = yreal(i)   - yreal(i-1)
        dyi1 = yimag(i+1) - yimag(i)
        dyi2 = yimag(i)   - yimag(i-1) 

        ureal(i) = ( 6d0*( dyr1/dx1 - dyr2/dx2 )/ dx3 - sig*ureal(i-1) )/preal
        uimag(i) = ( 6d0*( dyi1/dx1 - dyi2/dx2 )/ dx3 - sig*uimag(i-1) )/pimag
        
    enddo
  
    do  k=n-1,1,-1
        y2real(k) = y2real(k)*y2real(k+1) + ureal(k)
        y2imag(k) = y2imag(k)*y2imag(k+1) + uimag(k)
    enddo
 
     y2 = cmplx(y2real,y2imag)
 
    end subroutine spline_i
    
!==================================================================================================================================! 
!========================================================================================================================== splint_i
!==================================================================================================================================!
    complex(8) function splint_i(x)
 
! Performs spline interpolation on set of
! fkxbuff at wavenumbers 
! 
   
    integer     :: klo, khi, k
    real(8)     :: x, a, b, h2
    complex(8)  :: dydx

 
! locate which interval we are in
! x gt maximum tabulated value
    if( x >= x2(nwave) ) then ! use zero for higher wavenumbers
        splint_i = 0d0
    else
! x in tabulated range
    if( x >= x2(1) ) then
        klo = 1
        khi = nwave
 
 1      if ( (khi-klo) > 1) then
            k = (khi+klo)/2
            if ( x2(k) > x) then
                khi = k
            else
                klo = k
            endif
            goto 1
        endif 
    
        h2 = x2(khi) - x2(klo)        
        a  = ( x2(khi) - x      ) / h2
        b  = ( x      - x2(klo) ) / h2
        splint_i = a*fkxbuff(klo) + b*fkxbuff(khi) + ((a**3-a)*coeff(klo) + (b**3-b)*coeff(khi))*(h2**2)/6.0
    
    else
!    
! x is less than the minimum tabulated value:
!
       !  splint_i = 0d0   
!  Modified to use trapezoidal extrapolation, this helps a lot for LIN CSEM responses at long offsets
         dydx     = (fkxbuff(2)-fkxbuff(1) ) / (wavenum(2)-wavenum(1)) 
         splint_i = fkxbuff(1) + dydx*(10.**x - wavenum(1)) 

        endif
      endif

      return
      end  function splint_i          

    end module  spline_integ_kx  
 
    