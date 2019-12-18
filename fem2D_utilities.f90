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
!
! This module contains some subroutines for various generic 2D finite element
! computations.
! 
!

    module fem2D_utilities

    implicit none

!
! Subroutines:
!
    public :: fget_interp_bump   ! Computes the value of bump interpolant at a point.
    public :: fget_interp_lin    ! Computes the value of linear interpolant at a point.
    public :: get_abc_coeffs     ! Gets the linear basis function coefficients for an element.
    public :: get_linear_basis   ! Evaluates linear basis functions and its derivatives at a point.
    public :: get_bump_basis     ! Evaluates bump basis functions and its derivatives at a point.
    public :: findElement        ! Find which element a point is located in.
    public :: computeEdges       ! Given an element vertex and neighbor lists, this returns global edge list with dim = (3, nele)
    public :: getNode2Tri        ! Lists one triangle that each node is incident on

    private :: loopwalk, inElement, direction, onsegment

!
! Variables:
!
    integer, parameter, private :: eps(6) = (/ 1,2,3,1,2,3 /)      
    real(8), parameter, private :: edgetolerance = 1d-5  ! distance tolerance for considering whether a point is located on an edge

    contains


!======================================================================!
!================================================================! FGET
!======================================================================!
      subroutine fget_interp_lin(yp,zp,fe,ye,ze,fi)
!  
! Function to get the interpolated value of the field and its gradient
! in a linear basis element.
!
!======================================================================!
 
     
    implicit none
    
    real(8), intent(in)                  :: yp,zp
    complex(8), dimension(3), intent(in) :: fe
    complex(8), dimension(3),intent(out) :: fi             
    
    real(8), dimension(3)    :: ye,ze
    real(8), dimension(3)    :: l,dldy,dldz
        
!
! Evaluate linear basis functions at yp,zp:
!
    call get_linear_basis(ye,ze,yp,zp,l,dldy,dldz)
    
! 
!  Compute the value of the interpolant and the y and z derivatives
!

    fi(1) = sum(l*fe) 
    fi(2) = sum(dldy*fe) 
    fi(3) = sum(dldz*fe) 
    
    end subroutine fget_interp_lin
    
    
!!======================================================================!
!!====================================================! fget_interp_bump
!!======================================================================!  
   subroutine fget_interp_bump(yp,zp,fe,ye,ze,fi)
!  
! Function to get the interpolated value of the field and its gradient
! in a quadratic bump space element.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!  
! Created, June 2008.
!
!======================================================================!
     
    implicit none

    real(8), intent(in)                  :: yp,zp
    complex(8), dimension(3), intent(in) :: fe
    real(8),    dimension(3), intent(in) :: ye,ze    
    complex(8), dimension(3),intent(out) :: fi             
    
    real(8), dimension(3)    :: q,dqdy,dqdz

    integer,parameter        :: eps(6) = (/ 1,2,3,1,2,3 /)

!
!  Evaluate bump basis functions at yp,zp:
!
    call get_bump_basis(ye,ze,yp,zp,q,dqdy,dqdz)      
      
! 
!  Compute the value of the interpolant and the y and z derivatives
!  
    fi(1) = sum(q*fe)  
    fi(2) = sum(dqdy*fe)  
    fi(3) = sum(dqdz*fe)  
    
    end subroutine fget_interp_bump

    
!======================================================================!
!====================================================! get_abc_coeffs
!======================================================================!    
    subroutine get_abc_coeffs(ye,ze,a,b,c,area)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!       
    implicit none
    real(8), dimension(3)    :: ye,ze,a,b,c
    real(8)                  :: area
    
    a    = ye(eps(2:4))*ze(eps(3:5)) - ye(eps(3:5))*ze(eps(2:4))
    b    = ze(eps(2:4)) - ze(eps(3:5)) 
    c    = ye(eps(3:5)) - ye(eps(2:4))      
    
    area = dabs(a(1) + b(1)*ye(1) +c(1)*ze(1))/2.d0 
    
    end subroutine get_abc_coeffs
    
!======================================================================!
!====================================================! get_linear_basis
!======================================================================!    
    subroutine get_linear_basis(ye,ze,yp,zp,l,dldy,dldz)    
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!       
    implicit none 
    real(8), dimension(3)   :: ye,ze,a,b,c
    real(8)                 :: area, yp,zp
    real(8), dimension(3)   :: l,dldy,dldz    
    
    call get_abc_coeffs(ye,ze,a,b,c,area)
        
    l    = ( a + b*yp + c*zp) / 2.d0/area
    dldy = b/2.d0/area
    dldz = c/2.d0/area
    
    end subroutine get_linear_basis         
    
!======================================================================!
!====================================================! get_bump_basis
!======================================================================!    
    subroutine get_bump_basis(ye,ze,yp,zp,q,dqdy,dqdz)  
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!       
    implicit none 
    real(8), dimension(3)   :: ye,ze,a,b,c
    real(8)                 :: area, yp,zp
    real(8), dimension(3)   :: l  
    real(8), dimension(3)   :: q,dqdy,dqdz      

    !
    ! First get the linear basis:
    !
    call get_abc_coeffs(ye,ze,a,b,c,area)       
    
    l    = ( a + b*yp + c*zp ) / 2.d0/area
    
    !
    ! Now compute bump basis
    !
    q    = 4.d0*( l(eps(2:4)) * l(eps(3:5)) )
    dqdy = 4.d0*( l(eps(2:4)) * b(eps(3:5)) + l(eps(3:5)) * b(eps(2:4))) / 2.d0/area
    dqdz = 4.d0*( l(eps(2:4)) * c(eps(3:5)) + l(eps(3:5)) * c(eps(2:4))) / 2.d0/area    
      
    end subroutine get_bump_basis       
    

!======================================================================! 
!========================================================! computeEdges
!======================================================================! 
    subroutine computeEdges( nele, emap, neigh, edges, nedges )
!
! Computes global edge number list for a triangular grid.
! Edges 1,2,3 are opposite nodes 1,2,3.
!
! Revisions
! Nov. 18, 2009     New efficient version that uses neigh array 
! 2008              Original version
!
!  Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
!======================================================================! 
 
    
    implicit none 
   
    integer, intent(in)  :: nele 
    integer, intent(in)  :: emap(3,nele), neigh(3,nele)
    integer, intent(out) :: edges(3,nele), nedges
 
    ! Local:
    integer    :: e, i, j, n, v1, v2
 
 
    ! Initialize:
    edges = 0   ! no indices yet
    nedges = 0  
    
    ! Loop over elements and add the edges when v1 > v2, also add the index for the neighbor
    
    do e = 1,nele
        
        do i = 1,3   ! loop over three edges
        
            v1 = emap(eps(i+1),e)
            v2 = emap(eps(i+2),e)
            
            if ( (v2.gt.v1) ) then
            
                nedges = nedges + 1
                edges(i,e) = nedges
                
                n = neigh(i,e)

                if (n > 0) then  ! only if we have a neighbor... n = -1 means no neighbor
                    
                    do j = 1,3                      
                        if ( (emap(j,n) /= v1).and.(emap(j,n) /= v2) ) then  ! this is v3
                            edges(j,n) = nedges
                            exit
                        endif
                    enddo
                    
                endif
            elseif (neigh(i,e) == -1) then ! this is a boundary edge with reverse order, add it:
                nedges = nedges + 1
                edges(i,e) = nedges         
            endif  ! ( (v2.gt.v1) ) 
        
        enddo  ! i = 1,3
        
    enddo ! e = 1,nele
    
    return

    
    end subroutine computeEdges    

!======================================================================!
!==========================================================! findElement
!======================================================================!       
    subroutine findElement(tree, n,x,y, ne, emap, neigh, node2tri, attr, np, xp,yp, inelements) 
!
! Routine to quickly find the elements containing points in xp,yp arrays.
! Uses the jump-and-walk algorithm with the jump performed using a kd-tree 
! search to locate the nearest node to the query point.  If point is 
! on an edge ( dist < edgetolerance ) then the most conductive element
! is used. This works well for EM computations since field gradients 
! more robustly computed in conductors, and seafloor EM sites are located
! in conductive seawater, land EM receivers located in conductive earth not
! air.
!
! This routine scales much better than brute force "left-of-all-edges" 
! testing over every element. Assumes kdtree has been precomputed and 
! stored in kdtree2 type pointer "tree".
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Nov. 26, 2009     Added loopWalk test for point on a node case.
! Nov. 18, 2009     Created and testing.
!  
    
    use kdtree2_module
    use, intrinsic :: ieee_arithmetic
    
    implicit none
    
    type(kdtree2), pointer                :: tree     
    integer, intent(in)                   :: n, ne
    real(8), dimension(n), intent(in)     :: x, y
    integer, dimension(3,ne), intent(in)  :: emap, neigh
    integer,dimension(n), intent(in)      :: node2tri    
    real(8), dimension(ne), intent(in)    :: attr
    integer, intent(in)                   :: np
    real(8), dimension(np), intent(in)    :: xp,yp    
    integer, dimension(np), intent(out)   :: inelements

    !integer :: myid
        
    integer                               :: i, etest, inme, e0, n0  
    real(kdkind)                          :: qv(2) 
    type(kdtree2_result)                  :: results(1)   
  
   inelements = -1
    
 
    !
    ! Loop over query points (sites):
    !
    do i = 1,np
       
        if ( ieee_is_nan(xp(i)) .or. ieee_is_nan(yp(i)) )  then
            write(*,*) 'error, findElement has a nan: ',i,xp(i), yp(i)
        endif
        qv(1) = xp(i)
        qv(2) = yp(i)
 
         
        !
        ! Jump by finding nearest node:
        !
        call kdtree2_n_nearest(tp=tree,qv=qv,nn=1,results=results) 
        
        !
        ! Get an element this node is incident on:
        !
        etest =  node2tri(results(1)%idx)
 
        !
        ! Now do a careful walk to find the element containing this point:
        !
        if ( sqrt(results(1)%dis) <= edgetolerance ) then !  query point is on the nearest node
            
            ! We need to find the most conductive element containing this node, so lets do 
            ! a loop walk around all elements containing this node:
            
            ! Run around all the elements that use this node and pick most conductive: 
            e0 = etest
            n0 = results(1)%idx
            inme = etest
            
            call loopWalk(etest,e0,n0, ne, emap, neigh, attr, inme )
               
            inelements(i) = inme
              
        else ! Find the nearby element containing the query point
            
            call inElement(etest, xp(i),yp(i), n,x,y, ne, emap, neigh, attr, inelements(i) )
            
        endif
        
        if (inelements(i) == -1 ) then
            write(*,*) 'Error in findElement! could not find element for point: ',xp(i),yp(i)
            stop
        endif
    enddo
  
    
    end subroutine findElement

!======================================================================!
!===========================================================! loopWalk
!======================================================================!     
  recursive subroutine loopWalk(etest,e0,n0, ne, emap, neigh, attr, inme )
!
! Recursive subroutine that walks through all elements containing node n0
! and returns the element with the largest attr value (ie. conductivity).
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Nov. 26, 2009     Created and testing.
!  
          
    implicit none
    
    integer, intent(in)                   :: e0,n0,ne,etest
    integer, dimension(3,ne), intent(in)  :: emap, neigh
    real(8), dimension(ne), intent(in)    :: attr 
    integer, intent(inout)                :: inme 
    
    integer :: i, next
 
  
    !
    ! Compate attribute of test element and an earlier pick:
    !
    if ( attr(etest) > attr(inme) ) inme = etest
    
    !
    ! Select the next left element containing node n0:
    !
    do i = 1,3
        if ( emap(i,etest) == n0 ) then 
            next = neigh(eps(i+1),etest)
            if ( (next == e0) .or. (next == -1) ) then ! Made it back to the starting element, just return.
                return                                 ! Also if neighbor doesn't exist (-1) just return 
                                                       ! since this is boundary element and we don't really care 
                                                       ! about solution here anyway...
            else ! walk to the neighbor element:
            
                call loopWalk(next,e0,n0, ne, emap, neigh, attr, inme )
                exit ! do loop
            endif
        endif
    enddo
  
    end subroutine loopWalk
  
!======================================================================!
!===========================================================! inElement
!======================================================================!     
  recursive subroutine inElement(etest,xp,yp,nn,x,y,ne,emap,neigh,attr,inme)
!
! Recursive subroutine to test if element etest contains the query
! point xp,yp.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Nov. 18, 2009     Created and testing.
!  
     ! use mare2d_mpi_definitions  
      
    implicit none
    
    integer, intent(in)                   :: nn, ne
    integer, intent(inout)                :: etest
    real(8), dimension(nn), intent(in)    :: x, y
    integer, dimension(3,ne), intent(in)  :: emap, neigh
    real(8), dimension(ne), intent(in)    :: attr
    real(8), intent(in)                   :: xp,yp    
    integer, intent(out)                  :: inme 
    
    integer :: n(3),i
    real(8) :: dirc(3), dist, d1, d2
    logical :: lonsegment(3)

    !integer  myid 
!        call mpi_comm_rank( MPI_COMM_WORLD, myID, i )
        
    !
    ! Get direction from each element edge to the query point:
    !
    n = emap(1:3,etest)  ! vertices of element tested   
    do i = 1,3
        dirc(i) =  direction(x(n(eps(i+1))),y(n(eps(i+1))),x(n(eps(i+2))),y(n(eps(i+2))),xp,yp)  
    enddo

    if ( all(dirc >= edgetolerance) ) then
    
        ! Point is to the left of all segments so it is inside element etest
        inme = etest
        return
        
    else  ! See if its on an edge:
        do i = 1,3
            lonsegment(i) =  onsegment( x(n(eps(i+1))),y(n(eps(i+1))),x(n(eps(i+2))),y(n(eps(i+2))),xp,yp )
        enddo
         
!         if (myid == 4) then
!       
!           write(*,*) etest, real(dirc),lonsegment, neigh(1:3,etest)
!       endif
        
        do i = 1,3
            
            if   ( ( abs(dirc(i)) <= edgetolerance ) .and. ( lonsegment(i) ) ) then
                
                if ( neigh(i,etest) > 0 ) then
                  
                    if ( attr(etest) > attr( neigh(i,etest)) ) then ! choose most conductive
                        inme = etest
                    else
                        inme =  neigh(i,etest)
                    endif
 
                    return
                else
                    inme = etest ! its on an outer boundary edge...
                    return
                endif
                
            endif
            
        enddo

        
        ! Its not on an edge, so jump to the neighbor across the edge that the point is to the right-of.  
        do i = 1,3                
            if ( (dirc(i) < 0)  .and. (abs(dirc(i)) >= edgetolerance ) ) then ! point is to the right of edges
               
                if ( neigh(i,etest) > 0 ) then
                     etest = neigh(i,etest)
                     call inElement(etest, xp,yp, nn,x,y, ne, emap, neigh, attr, inme)
                     return
!                else
!                    inme = -1 ! just give up...
!                    return
                endif
            endif
        enddo
        
        do i = 1,3
              d1 = sqrt(  (x(n(eps(i+1))) - xp)**2 +  (y(n(eps(i+1))) - yp)**2 )
              d2 = sqrt(  (x(n(eps(i+2))) - xp)**2 +  (y(n(eps(i+2))) - yp)**2 )
              dist = min(d1,d2)
             if ( (dirc(i) < 0 ) .and. (dist > edgetolerance ) ) then
                 if ( neigh(i,etest) > 0 ) then
                     etest = neigh(i,etest)
                     call inElement(etest, xp,yp, nn,x,y, ne, emap, neigh, attr, inme)
                     return
!                else
!                    inme = -1 ! just give up...
!                    return
                endif               
             endif
        enddo
  
        write(*,*) 'inElement: etest, neight(1:3): ', etest,  neigh(1:3,etest)
        write(*,*) 'inElement: dirc(1:3): ',   real(dirc)
        write(*,*) 'inElement: xp, xe(1:3): ', real(xp),real(x(emap(:,etest))) 
        write(*,*) 'inElement: yp, ye(1:3): ', real(yp),real(y(emap(:,etest)))
        write(*,*) 'inElement: lonsegment: '  , lonsegment       
        write(*,*) ' '
        !
        ! kwk debug: need a check for when point is outside entire mesh!
        ! 
    endif
    
    
    end subroutine inElement
    
!======================================================================!
!===========================================================! direction
!======================================================================!      
     real(8) function direction(p1x,p1y,p2x,p2y,p3x,p3y) 
!
! Function that computes the signed orthongal intersection distance using
! the cross product (p2-p1) x (p3-p1) normalized by p3-p1.
!
! Result:
!          positive:  p3 is left of line extending along p2-p1
!          negative:  p3 is right of line extending along p2-p1
!          zero:      p3 lies on line extending along p2-p1
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Nov. 18, 2009     Created and testing.
!  
    implicit none
    
    real(8) :: p1x,p1y,p2x,p2y,p3x,p3y
    
    direction = ((p2x - p1x)*(p3y - p1y) - (p3x - p1x)*(p2y - p1y) ) / &
              &     sqrt( (p3y - p1y)**2 + (p3x - p1x)**2)
    
    
    end function direction

!======================================================================!
!===========================================================! onsegment
!======================================================================!      
     logical function onsegment(p1x,p1y,p2x,p2y,p3x,p3y)    
!
! Function that determines if a colinear point p3 is between end points
! p1 and p2.  Uses edgetolerance to allow for inexact point data.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! Nov. 18, 2009     Created and testing.
!  
    implicit none
    
    real(8) :: p1x,p1y,p2x,p2y,p3x,p3y

    if ( (min(p1x,p2x) - edgetolerance <= p3x ) .and. (max(p1x,p2x) + edgetolerance >= p3x ) .and. &
       & (min(p1y,p2y) - edgetolerance <= p3y ) .and. (max(p1y,p2y) + edgetolerance >= p3y ) ) then
        onsegment = .true.
        
    else
        onsegment = .false.
    endif
    
    end function onsegment    
!======================================================================! 
!=========================================================! getNode2Tri
!======================================================================! 
    subroutine getNode2Tri( nele, emap, nnod, node2tri)
!
! Computes a nnod length list of one triangle incident with each vertex.
!
! Nov. 18, 2009     Created and testing.
!
!  Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
!======================================================================! 
 
    
    implicit none 
   
    integer, intent(in)  :: nele,nnod
    integer, intent(in)  :: emap(3,nele)
    integer, intent(out) :: node2tri(nnod)
    
    integer :: e
    
    !
    ! Loop over all elements and just fill in node2tri, overwritting 
    ! any existing entries.  This is a fast and easy way to create this 
    ! index array.
    !
    
    do e=1,nele
        node2tri(emap(1:3,e)) = e
    enddo
 
    end subroutine getNode2Tri
    
    
    
    end module fem2D_utilities