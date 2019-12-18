!-----------------------------------------------------------------------
!
!    Copyright 2008-2012
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
!======================================================================!
!=======================================================! module PW1D
!======================================================================!
! Fortran module for 1D EM plane-wave computations
!
! This is a work in progress towards making a Fortran OO code for various
! 1D MT related computations.  Some of the object oriented functionality 
! of Fortran 2003 isn't yet implemented in today's compilers, so this is 
! only a partial implementation of the MT1D object. 
! 
! 
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 1.0, September 2009
!
!
    module MT1D
    
    implicit none
    
    type :: PW1D 

        ! Public: 
        ! users need to allocate and set these:
        integer                                :: nlayer
        real(8), dimension(:), allocatable     :: zlay, sig  ! top depth and conductivity   
        real(8)                                :: omega 
        
        ! Private (but compiler doesn't support private attr yet for derived types...)
        complex(8), dimension(:), allocatable  :: a,b,c,d  ! layer coefficients
       
       ! contains
       ! <put procedure pointers here>...when the compilers finally support this F2003 standard method!
       ! then you could do this:
       ! type(pw1d) :: mymodel
       ! call mymodel%initialze
       ! call mymodel%getEH(z)
       
    end type PW1D

    ! Constants:
    real(8), parameter, private     :: pi = 3.141592653589793d0
    real(8), parameter, private     :: mu0 = 4d-7*pi    
    complex(8), parameter, private  :: ii = (0d0,1d0)       
    
    contains 
        

!======================================================================!
!=======================================================! getMT1Dcoeffs
!======================================================================!

    subroutine getPW1Dcoeffs(model)
!
! Computes the layer propagation coefficients for 
! plane-wave magnetotelluric electric and magnetic fields as
! a function of depth in a 1D model.  Uses a stable recursions similar to the 
! ones I used in Dipole1D, but modified for the simpler case of a 
! plane wave. 
!
!
! The routine finds coefficients a,b,c,d for each layer i, so that fields 
! at depth z have the form:
!
!   E_i(z) = a_i * exp( +k_i * (z - z_i+1) ) + b_i * exp( -k_i * (z - z_i) )
!   H_i(z) = c_i * exp( +k_i * (z - z_i+1) ) + d_i * exp( -k_i * (z - z_i) )
!
!   where k_i = sqrt( -sqrt(-1)*omega*mu0*sig_i )
!
! Fields are normalized so that E and H at the top of the model have unit magnitude.
!
! After calling this routine, you can get the fields by calling 
! subroutine getMT1Dfields().
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 1.1, March 2010       Fixed error with top layer by adding a 0th layer.
! Version 1.0, September 2009
!
    implicit none
    
    ! I/O:
    type (PW1D), intent(inout) :: model
 

    ! Local variables:
    integer    :: i
    complex(8) :: gmogp, rjexp, Atop 
    complex(8), dimension(:), allocatable :: kk, Rp, expmgh
    
!
! Compute recursions:
!
    allocate( kk(0:model%nlayer), Rp(0:model%nlayer), expmgh(0:model%nlayer) )
        
    kk(1:model%nlayer)  = sqrt(-ii*model%omega*mu0*model%sig(1:model%nlayer))
    kk(0) = sqrt(-ii*model%omega*mu0*model%sig(1))
    
    expmgh = 0d0    
    do i = 1,model%nlayer-1
        expmgh(i) = exp( -kk(i)*( model%zlay(i+1) - model%zlay(i) ) )
    enddo
   
    Rp = 0d0
    
    do i = model%nlayer-1,1,-1
       gmogp   = (-ii*model%omega*mu0*(model%sig(i)-model%sig(i+1) )) / (kk(i)+kk(i+1))**2
       rjexp   = Rp(i+1) * expmgh(i+1)
       Rp(i)   =  (gmogp + rjexp)   / ( 1.d0 + gmogp*rjexp)   * expmgh(i)
    enddo
    
!
! Back propagate after setting top layer to boundary conditions (a=0,b=1), (c=0,d=1):
!   
    if ( (.not.allocated(model%a) ) .or. ( size(model%a) /= model%nlayer ) ) then
        allocate( model%a(0:model%nlayer), model%b(0:model%nlayer), model%c(0:model%nlayer), model%d(0:model%nlayer))
    endif
    model%a = 0; model%b = 0; model%c = 0; model%d = 0;   
    model%a(0) = 1; model%c(0) = 1;   
    
    ! This propagates the solution across each interface via field continuity conditions:
    do i = 1,model%nlayer
    
         ! E coefficients
        Atop = (model%a(i-1) + model%b(i-1)*expmgh(i-1) )   
        model%b(i)  = Atop / (1.d0 + Rp(i)*expmgh(i) )
        model%a(i)  = model%b(i)*Rp(i)
    
        ! H coefficients:
        Atop = (model%c(i-1) + model%d(i-1)*expmgh(i-1) )    ! H uses same recursion coeff, but with reversed sign (-):
        model%d(i)  = Atop / (1.d0 - Rp(i)*expmgh(i) )
        model%c(i)  = -model%d(i)*Rp(i)       
        
    enddo
!
! Deallocate:
!
    deallocate ( kk, Rp, expmgh )
    
    end subroutine getPW1Dcoeffs
    



!======================================================================!
!=======================================================! getMT1Dfields
!======================================================================!

    subroutine getPW1Dfields(model,z,eh) 
!
! Call this after calling getMT1Dcoeffs(model)
!
! Outputs the E and H fields at absolute depth z using:
!
!   E_i(z) = a_i * exp( +k_i * (z - z_i+1) ) + b_i * exp( -k_i * (z - z_i) )
!   H_i(z) = c_i * exp( +k_i * (z - z_i+1) ) + d_i * exp( -k_i * (z - z_i) )
!
!   where k_i = sqrt( -sqrt(-1)*omega*mu0*sig_i )
!   
! Fields are normalized so that E and H at top of the model (i.e., z = z_1)
! have unit magnitude.
! 
! Note that auxiliary fields for MT impedances can be computed using 
! Ey = 1/sigma dz Hx and Hy = dz Ex, but I haven't coded impedances yet.
! Do not use eh(1) / eh(2) since both E and H are normalized here to unity
! at the top of the model, which is useful for creating 1D fields as 
! boundary conditions in 2D finite element models.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Version 1.0, September 2009
!
    implicit none
    
    ! I/O:
    type (PW1D), intent(in)              :: model   
    real(8), intent(in)                  :: z        ! depth to output E and H at
    complex(8),dimension(2), intent(out) :: eh

    ! Local variables:
    integer    :: i,ilay
    complex(8) :: kk, expp, expm
    
!
! First find the layer:
!
    ilay = 1
    do i = 1,model%nlayer
        if (z >= model%zlay(i)) then
            ilay = i
        endif
    enddo

!
! Compute the fields at depth z:
!

    kk = sqrt(-ii*model%omega*mu0*model%sig(ilay))
    
    expp = 0
    expm = 0
    
    expm = exp( -kk * (z - model%zlay(ilay) ) )
    if (ilay /= model%nlayer) expp = exp( kk * (z - model%zlay(ilay+1) ))
    
    eh(1) = model%a(ilay) * expp  + model%b(ilay) * expm  ! E
    eh(2) = model%c(ilay) * expp  + model%d(ilay) * expm  ! H

    ! Zero out low values:
    if (abs(eh(1)) < 1d-80) eh(1) = 0d0
    if (abs(eh(2)) < 1d-80) eh(2) = 0d0
    
    end subroutine getPW1Dfields
    
!======================================================================!
!=======================================================! PW1Dfinalize
!======================================================================!   
    subroutine PW1Dfinalize(model)
 !
 ! Deallocates allocatable parts of derived type model
 !
    
    implicit none
    
    ! I/O:
    type (PW1D), intent(inout)              :: model   
    
    deallocate(model%a,model%b,model%c,model%d,model%sig,model%zlay)
    
    end subroutine PW1Dfinalize

    end module MT1D