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
! This file contains subroutines for interfacing with Triangle.c
!
!
!  Uses external routines: c_fortran_triangluate_input, c_fortran_triangluate_output
!
! DGM 2/23/2012 - Windows compilation support
#IF DEFINED(_WIN32) .OR. DEFINED(_WIN64)
#DEFINE c_fortran_triangluate_input  c_fortran_triangluate_input_
#DEFINE c_fortran_triangluate_output c_fortran_triangluate_output_
#ENDIF

    module triangle_mesh
    
    implicit none 
    
    type trimesh
    
        integer                                 :: nnod, nele  
        real(8), dimension(:),   allocatable    :: y,z             ! nnod node positions for this mesh
        integer, dimension(:,:), allocatable    :: emap            ! 3 by nele elments mapping and 3 by edges edge mapping     
        real(8), dimension(:),   allocatable    :: attr             ! nele element conductivity
 
        real(8), dimension(:), allocatable      :: pointattributelist
        integer, dimension(:), allocatable      :: pointmarkerlist
        integer                                 :: numberofpointattributes
 
        real(8), dimension(:), allocatable      :: area ! only use in refinement comp, input only. NOT allocated in this module
        integer, dimension(:,:), allocatable    :: neighborlist
        integer                                 :: numberofcorners
        integer                                 :: numberoftriangleattributes
        integer                                 :: numberofsegments
        integer, dimension(:), allocatable      :: segmentlist
        integer, dimension(:), allocatable      :: segmentmarkerlist
        real(8), dimension(:), allocatable      :: holelist
        integer                                 :: numberofholes
        real(8), dimension(:), allocatable      :: regionlist
        integer                                 :: numberofregions      
    
    end type trimesh
    
    public  :: call_triangle, copy_trimesh, deallocate_trimesh, smooth_TriMesh 
    public  :: read_triangle, write_triangle 
    public  :: check_model   !---> checks all input segment intersections for slivers, good safety check before refinement! 
    private :: read_mesh
    
    contains
    
!======================================================================! 
!========================================================! call_triangle
!======================================================================!    
    subroutine call_triangle(tricommand, in, out)       
!
! Subroutine to pass mesh structure to Triangle.c so that
! a mesh can be generated or refined. Kludgey, yes.  
!
! "out" trimesh structure is optional.  If not present, then "in" trimesh 
! structure is overwritten.  If it is present, then its arrays are allocated
! here.
!
!  Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
!  Versions:
!  2.0  Nov. 20, 2009.  Modified for new trimesh structure
!  1.0  April 2009
!
    
    implicit none

    character(24), intent(in)               :: tricommand     
    type(trimesh), intent(inout)            :: in
    type(trimesh), intent(inout), optional  :: out
    
    ! local variables:
    
    integer :: nnod,  numberofpointattributes 
    integer :: nele,  numberofcorners, numberoftriangleattributes, numberofsegments 
    
    
    
!
! Step 1: Pass mesh to triangle and receive back the dimensions of the new arrays and a pointer to them:
!         This seems to be the only way to deal with a C structure that is allocated in Triangle.c.
!
    if ( .not.allocated( in%area ) )  then   ! allocate dummy array so Intel compiler doesn't complain when doing -traceback -check
        allocate( in%area(in%nele)  )        ! However this code works without this dummy allocation.
        in%area = -1
    endif

    call c_fortran_triangluate_input( tricommand,  &
        &  in%y,in%z, in%pointattributelist, in%pointmarkerlist, &
        &  in%nnod,  in%numberofpointattributes,  &
        &  in%emap, in%attr, in%area, &
        &  in%nele,  in%numberofcorners, in%numberoftriangleattributes, &
        &  in%segmentlist, in%segmentmarkerlist, in%numberofsegments, &
        &  in%holelist,  in%numberofholes, &
        &  in%regionlist, in%numberofregions, &       
        &  nnod,  numberofpointattributes,                         &
        &  nele,  numberofcorners, numberoftriangleattributes, &
        &  numberofsegments  )

  !
  ! Now we need to either allocate output arrays, or else deallocate and then overwrite input arrays
  ! 
    if (present(out)) then  ! arrays will go into the out structure:
 
        !
        ! Allocate output Fortran arrays for the new grid:
        !
        ! kwk debug: should put some allocation error checks here...
        
        out%nnod = nnod
        out%numberofpointattributes  = numberofpointattributes
        out%nele = nele
        out%numberofcorners = in%numberofcorners
        out%numberoftriangleattributes = in%numberoftriangleattributes
        out%numberofsegments = numberofsegments
 
        
        allocate ( out%y( out%nnod ) , out%z( out%nnod )  )
        allocate ( out%pointattributelist( out%numberofpointattributes*out%nnod )  )    
        allocate ( out%pointmarkerlist( out%nnod )  ) 
        allocate ( out%emap(3, out%nele  ) )
        allocate ( out%attr(out%nele) )
        allocate ( out%neighborlist(3, out%nele) )
        allocate ( out%segmentlist( 2*out%numberofsegments)  )    
        allocate ( out%segmentmarkerlist( out%numberofsegments)  )          
        
        !
        ! Step 2: Receive the new mesh arrays from c_fortran_triangle.c
        !
        call c_fortran_triangluate_output( &
        &   out%y,out%z,  out%pointattributelist, out%pointmarkerlist, &
        &   out%emap,   out%attr, &
        &   out%neighborlist, &
        &   out%segmentlist, &
        &   out%segmentmarkerlist)
                
                
        ! copy across region and hole lists:
        out%numberofregions = in%numberofregions
        allocate (  out%regionlist(4*out%numberofregions ) )
        out%regionlist = in%regionlist
        
        out%numberofholes = in%numberofholes
        allocate (  out%holelist(2*out%numberofholes ) ) 
        out%holelist = in%holelist
 
        
        
    else ! overwrite the in structure
        !
        ! Deallocate old  mesh arrays:
        !
 
        if ( allocated( in%y ) )                      deallocate ( in%y )
        if ( allocated( in%z ) )                      deallocate ( in%z )    
        if ( allocated( in%pointattributelist ) )     deallocate ( in%pointattributelist   )    
        if ( allocated( in%pointmarkerlist ) )        deallocate ( in%pointmarkerlist   ) 
        if ( allocated( in%emap ) )                   deallocate ( in%emap  )
        if ( allocated( in%attr ) )                   deallocate ( in%attr  )
        if ( allocated( in%area ) )                   deallocate ( in%area  )
        if ( allocated( in%neighborlist ) )           deallocate ( in%neighborlist  )
        if ( allocated( in%segmentlist  ) )           deallocate ( in%segmentlist   )    
        if ( allocated( in%segmentmarkerlist ) )      deallocate ( in%segmentmarkerlist  )          
        ! Note, regions and holes stay the same, so no deallocation of them

        !
        ! Allocate output Fortran arrays for the new grid:
        !
        in%nnod = nnod
        in%numberofpointattributes  = numberofpointattributes
        in%nele = nele
       ! in%numberofcorners = in%numberofcorners
       ! in%numberoftriangleattributes = in%numberoftriangleattributes
        in%numberofsegments = numberofsegments
        

        allocate ( in%y( in%nnod ) , in%z( in%nnod )  )
        allocate ( in%pointattributelist( in%numberofpointattributes*in%nnod )  )    
        allocate ( in%pointmarkerlist( in%nnod )  ) 
        allocate ( in%emap(3, in%nele  ) )
        allocate ( in%attr(in%nele) )
        !allocate ( in%area(in%nele) )
        allocate ( in%neighborlist(3, in%nele) )
        allocate ( in%segmentlist( 2*in%numberofsegments)  )    
        allocate ( in%segmentmarkerlist( in%numberofsegments)  )        
        

        !
        ! Step 2: Receive the new mesh arrays from c_fortran_triangle.c
        !
        call c_fortran_triangluate_output( &
        &   in%y,in%z,  in%pointattributelist, in%pointmarkerlist, &
        &   in%emap,   in%attr, &
        &   in%neighborlist, &
        &   in%segmentlist, &
        &   in%segmentmarkerlist)
                

        
    endif
 
    end subroutine call_triangle
    
    
!======================================================================! 
!========================================================! copy_trimesh
!======================================================================!    
    subroutine copy_trimesh( in,out )
        
!
! Copies trimesh from one structure to another with the allocations done here.
!
!  Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
!  Nov. 20, 2009   Created and tested.
!
    
    implicit none
   
    type(trimesh), intent(in)      :: in
    type(trimesh), intent(inout)  :: out

    out%nnod = in%nnod
    
    allocate (out%y(out%nnod))
    out%y = in%y
    allocate (out%z(out%nnod))  
    out%z = in%z    
    
    out%nele = in%nele
    allocate( out%attr(in%nele), out%emap(3,in%nele) ) !
    
    out%attr = in%attr
    out%emap = in%emap
    
    allocate (out%neighborlist(3,in%nele) )

    out%neighborlist= in%neighborlist

    out%numberofpointattributes = in%numberofpointattributes
    allocate(out%pointattributelist(size(in%pointattributelist)))
    out%pointattributelist = in%pointattributelist

    allocate(out%pointmarkerlist(size(in%pointmarkerlist)))
    out%pointmarkerlist = in%pointmarkerlist 
    
    
!    allocate(out%area(out%nele))
!    out%area = in%area

    out%numberofcorners = in%numberofcorners
    
    out%numberoftriangleattributes = in%numberoftriangleattributes
    
    out%numberofsegments = in%numberofsegments
    

    allocate( out%segmentlist (size(in%segmentlist)) )
    out%segmentlist = in%segmentlist
    allocate( out%segmentmarkerlist (size(in%segmentmarkerlist)) )
    out%segmentmarkerlist = in%segmentmarkerlist    
 

    out%numberofholes = in%numberofholes
    allocate(out%holelist(size(in%holelist)))
    out%holelist = in%holelist
     

    out%numberofregions = in%numberofregions
    allocate(out%regionlist(size(in%regionlist)))
    out%regionlist = in%regionlist

    end subroutine copy_trimesh
!======================================================================! 
!==================================================! deallocate_trimesh 
!======================================================================!    
    subroutine deallocate_trimesh( model , lsave)
!        
! Deallocates all the allocatable fields in structure trimesh
!
!
!  Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
!  Versions:
!  1.0 Nov. 20, 2009.  Created and tested.

 
    
    implicit none
   
    type(trimesh), intent(inout)  :: model
    logical, intent(in), optional :: lsave

    if ( allocated ( model%y) )                    deallocate ( model%y )
    if ( allocated ( model%z) )                    deallocate ( model%z )   
    if ( allocated ( model%attr ) )                 deallocate ( model%attr )
    if ( allocated ( model%emap ) )                deallocate ( model%emap )
    if ( allocated ( model%neighborlist ) )        deallocate ( model%neighborlist )  
    if ( allocated ( model%pointattributelist )  ) deallocate ( model%pointattributelist )
    if ( allocated ( model%pointmarkerlist )  )    deallocate ( model%pointmarkerlist )
    if ( allocated ( model%area)  )                deallocate ( model%area )
    if ( allocated ( model%segmentlist) )          deallocate ( model%segmentlist )
    if ( allocated ( model%segmentmarkerlist )  )  deallocate ( model%segmentmarkerlist )
    if ( present(lsave) ) then
        if (lsave) return
    endif
    if ( allocated ( model%holelist )  )           deallocate ( model%holelist )
    if ( allocated ( model%regionlist )  )         deallocate ( model%regionlist )


    end subroutine deallocate_trimesh
    
!======================================================================! 
!======================================================! write_triangle
!======================================================================!    
    subroutine write_triangle(model, fileroot)  
!
! Writes out the .poly, .node and .ele files for the mesh in  structure 'in'
!
!  Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
!  Versions:
!  2.0 Nov. 20, 2009        Modified for new trimesh structure
!  1.0  April 2009
!   
 
    implicit none
    
    type(trimesh), intent(in)  :: model
    character(256),intent(in)  :: fileroot
  
    character(80) :: nfmt
    character(80)  :: cnpa
    
    integer :: i,j
    
 !
 ! .poly file:
 !
    open(16,file=trim(fileroot) //'.poly',status='REPLACE') 
    write(16,*) 0, 2, 0, 1  ! no nodes for for the intermediate .poly file, nodes are in .node file
    write(16,*) model%numberofsegments,1  ! segment heading   
    do i = 1,model%numberofsegments
        write(16,*) i,  model%segmentlist( 2*i-1),model%segmentlist( 2*i),model%segmentmarkerlist(i)
    enddo
    ! holes
    write (16,*) model%numberofholes
    do i = 1,model%numberofholes
        !Following lines: <hole #> <x> <y>
        write(16,'(i0,1x,2(E22.14,1x))') i,model%holelist(2*i-1), model%holelist(2*i)
    enddo       
    ! regions:
    write (16,*) model%numberofregions
    do i = 1,model%numberofregions
     
        !Following lines: <region #> <x> <y> attribute area_constraint
        write(16,'(i0,1x,4(E22.14,1x))') i, model%regionlist(4*i - 3 : 4*i)
    enddo               

    close(16)
    
 !
 ! .node file:
 !   
    open(16,file=trim(fileroot) //'.node',status='REPLACE') 
    write(16,*) model%nnod, 2, model%numberofpointattributes, 1
    
    
    write(cnpa,*) model%numberofpointattributes
    nfmt = '(1(i0,1x),2(E22.14,1x),'//trim(adjustl(cnpa))//'(E22.14,1x),'// 'E22.14)'
    !write(*,*) 'nfmt: ',trim(nfmt)
    do i = 1,model%nnod
         write(16,fmt=nfmt) i, model%y(i),model%z(i), &
                        & ( model%pointattributelist( (i-1)*model%numberofpointattributes + j), &
                        &  j=1,model%numberofpointattributes) , &
                        &  model%pointmarkerlist(i) 
    enddo
    close(16)
    
!
! .ele file:
!
 
    open(16,file=trim(fileroot) //'.ele',status='REPLACE') 
    write(16,*) model%nele,3,model%numberoftriangleattributes 
    do i = 1,model%nele
        write(16,fmt='(4(i0,1x),(E22.14,1x))') i, model%emap(1:3,i),model%attr(i)
    enddo 
    close(16)
    
! 
! KWK debug: write out neighbors?    
!    

    end subroutine write_triangle
    
 

!======================================================================! 
!=======================================================! read_triangle
!======================================================================! 
    subroutine read_triangle(fileroot,model)
!
! Reads the .poly file for the model and optionally any .node and .ele 
! files. Mesh is stored in trimesh structure "out".
!
!  Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
   
    implicit none
  
    type(trimesh), intent(inout) :: model
    character(256),intent(in)    :: fileroot
  
    integer             :: ierr, iskip, nbndmarkers, i, j, npa
    character(256)      :: filename
    
    !
    ! Open the file:
    !
    filename = trim(fileroot)//'.poly'
     
    write(*,fmt='(a,a)') ' Opening Poly file:              ',trim(filename)
    open(unit=10,file=trim(filename),status='old',iostat = ierr)
    if (ierr .ne. 0) then
        write(*,*) ' Error opening POLY file: ',filename
        stop
    end if    
    read (10,*) model%nnod, iskip, model%numberofpointattributes, nbndmarkers
  !  write(*,*) 'model%nnod, iskip, model%numberofpointattributes, nbndmarkers:',model%nnod, iskip, &
  !             & model%numberofpointattributes, nbndmarkers
    
    
    if (model%nnod == 0) then
        call read_mesh(fileroot,model)  ! Read in nodes and elements from <fileroot>.node and <fileroot>.ele files
        ! allocates ypoint and zpoint arrays and fills them in
        ! allocates attr and emap and fills them in.
    else
        model%nele = 0
        model%numberofcorners  = 3
        model%numberoftriangleattributes = 1
        
        allocate (model%emap(3,model%nele), model%attr(model%nele) )   
        allocate (model%neighborlist(3,model%nele) )
         
        !
        ! Allocate point sized arrays:
        !
        allocate ( model%y(model%nnod), model%z(model%nnod) )
        allocate ( model%pointattributelist( model%numberofpointattributes*model%nnod )  )
        allocate ( model%pointmarkerlist( model%nnod )  )
        !
        ! Read in points:
        !
        npa = model%numberofpointattributes
        do i = 1,model%nnod
            ! : <vertex #> <x> <y> [attributes] [boundary marker]
               read(10,*) iskip, model%y(i),model%z(i), &
                        & ( model%pointattributelist((i-1)*npa+j) , j=1,npa) , &
                        & ( model%pointmarkerlist(i), j = 1,nbndmarkers)  ! nbndmarkers = 0|1
            
        enddo
    endif 
    
    allocate (model%area(model%nele) )

    !
    ! Read in segments:
    !
    read (10,*) model%numberofsegments, nbndmarkers
    !write(*,*) 'model%numberofsegments, nbndmarkers:',model%numberofsegments, nbndmarkers
    ! allocate segment sized arrays:
    allocate ( model%segmentlist(model%numberofsegments*2))
    allocate (model%segmentmarkerlist(model%numberofsegments) ) 
    if (nbndmarkers > 0 ) then
        do i = 1,model%numberofsegments
           !Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
            read(10,*) iskip,  model%segmentlist(2*i-1), model%segmentlist(2*i) , model%segmentmarkerlist(i)
           !  write(*,*) 'segs:', i,  model%segmentlist(2*i-1), model%segmentlist(2*i) 
        enddo
    else
        do i = 1,model%numberofsegments
           !Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
            read(10,*) iskip,  model%segmentlist(2*i-1), model%segmentlist(2*i) 
        enddo
    endif           
        
    !
    ! Read in the holes:
    !
    read (10,*) model%numberofholes
    
    allocate (model%holelist(2*model%numberofholes) )
    do i = 1,model%numberofholes
        !Following lines: <hole #> <x> <y>
        read(10,*) iskip,model%holelist(2*i-1), model%holelist(2*i)
    enddo       
    
    !
    ! Read in regions:
    !
    read (10,*) model%numberofregions
    allocate (  model%regionlist(model%numberofregions*4 ) )
    do i = 1,model%numberofregions
     
        !Following lines: <region #> <x> <y> attribute area_constraint
        read(10,*) iskip, model%regionlist(4*i - 3 : 4*i)
    enddo           
 
    !
    ! Close the file:
    !
    close(10)
 
    
    end subroutine read_triangle
    
!======================================================================! 
!===========================================================! read_mesh
!======================================================================! 
    subroutine read_mesh(fileroot,model) 
    
    implicit none
  
    type(trimesh), intent(inout) :: model
    character(256),intent(in)  :: fileroot
        
  ! Local:
    integer        ::  ierr, e,i,j, iskip, nbndmarkers, npa
    character(256) :: file2read

!
!  Read in NODE file
!
    file2read = trim(fileroot)//'.node'
    write(*,*) 'Opening Node file: ',trim(file2read)
    open (unit=11,file=trim(file2read),status='old',iostat=ierr)
    if (ierr .ne. 0) then
        write(*,*) ' Error opening NODE file: '//trim(file2read) 
        stop
    end if
    read (11,*) model%nnod, iskip, model%numberofpointattributes,nbndmarkers
   ! write (*,*) model%nnod, iskip, model%numberofpointattributes,nbndmarkers
    allocate ( model%y(model%nnod), model%z(model%nnod),  model%pointmarkerlist(model%nnod))    
    allocate ( model%pointattributelist( model%numberofpointattributes*model%nnod )  )
    !
    ! Read in points:
    !
    npa = model%numberofpointattributes
    do i = 1,model%nnod
   ! write(*,*) i, npa,nbndmarkers
        ! : <vertex #> <x> <y> [attributes] [boundary marker]
           read(11,*) iskip, model%y(i),model%z(i), &
                    & ( model%pointattributelist((i-1)*npa+j) , j=1,npa) , &
                    & ( model%pointmarkerlist(i), j = 1,nbndmarkers)  ! nbndmarkers = 0|1
 !                   write(*,*) iskip,model%y(i),model%z(i), &
!                    & ( model%pointattributelist((i-1)*npa+j) , j=1,npa) , &
!                    & ( model%pointmarkerlist(i), j = 1,nbndmarkers)
        
    enddo
    
    close(11)
  
!
!  Read in ELEMENT file
!  Read in the element map where emap(i,e) contains the index of the
!  i'th node (i=1,3) which defines the e'th element (e=1,nele).
!
    file2read = trim(fileroot)//'.ele'
    write(*,*) 'Opening Element file: ',trim(file2read)
    open (unit=11,file=trim(file2read),status='old',iostat=ierr)
    if (ierr .ne. 0) then
        write(*,*) ' Error opening element file: '//trim(file2read)
        stop
    end if
    
    read (11,*) model%nele,model%numberofcorners,model%numberoftriangleattributes


    allocate ( model%emap(3,model%nele),model%attr(model%nele) ) 
       
    do e=1,model%nele
        read(11,*) iskip, (model%emap(i,e),i=1,3),model%attr(e)
        !model%attr(e) = 1/model%attr(e) ! convert from input resistivity to conductivity
    enddo


    close(11)
    
!
! Read in NEIGHBOR file
! This file lists the 3 neighboring triangles for each triangle and is 
! useful for traversing the mesh.
!

    file2read = trim(fileroot)//'.neigh'
    write(*,*) 'Opening Neighbor file: ',trim(file2read)
    open (unit=11,file=trim(file2read) ,status='old',iostat=ierr)
    if (ierr .ne. 0) then
        write(*,*)  ' Error opening element neighbor file: '//trim(file2read) 
         write(*,*) ' Did you remember to include the ''n'' flag when calling Triangle?'
         write(*,*) ' Stopping! '
        stop
    end if
    
    read (11,*) model%nele,iskip
    
    
    allocate ( model%neighborlist(3,model%nele) ) 
    
    read(11,*) (iskip, (model%neighborlist(i,e),i=1,3) , e = 1,model%nele)
    
    close(11)
    
   
    
    end subroutine read_mesh        
    
    
!======================================================================! 
!======================================================! smooth_TriMesh
!======================================================================!        
    subroutine smooth_TriMesh(this,nsmooth, bTalk)
    
    implicit none
    type(trimesh),  intent(inout)   :: this   
    integer, intent(in)             :: nsmooth
    logical, intent(in)             :: bTalk
    
    real(8),dimension(:,:), allocatable :: icenters 
    
    integer                      :: i,j,k,e,n(3), ering(50), estart, nring,ismooth
    real(8)                      :: x(2), xpart(2), alpha, xtest(2), areatest(50)
    real(8)                      :: dx, maxdx,verts(2,3),lastmaxdx
    
    integer, dimension(:), allocatable   :: lfixedNodes,nadj,node2tri
    integer, dimension(:,:), allocatable ::  nodeAdj

    integer, parameter          :: epss(6) = (/ 1,2,3,1,2,3 /)  
    
    allocate(icenters(2,this%nele),lfixedNodes(this%nnod))
    
    allocate(nodeAdj(40,this%nnod),nadj(this%nnod))
 
    allocate(node2tri(this%nnod))
    
 
    lastmaxdx = 1d100
    
    !
    ! Loop over all elements and just fill in node2tri, overwritting 
    ! any existing entries.  This is a fast and easy way to create this 
    ! index array.
    !
    
    do e=1,this%nele
        node2tri(this%emap(1:3,e)) = e
    enddo
 
    !
    ! Flag fixed nodes:
    !    
    lfixedNodes = 0
    lfixedNodes(this%segmentlist) = 1
 
 
    do i = 1,this%nele
        n     = this%emap(:,i)
        verts(1,:) = this%y(n) 
        verts(2,:) = this%z(n)        
        icenters(1:2,i) = sum(verts,2)/3d0
    enddo 
 
     
    nodeAdj = 0
    nadj = 0  
    
    do i = 1,this%nnod  
    
        if ( lfixedNodes(i) /= 1) then
        
            estart = node2tri(i)
            e = estart
 
            do 
                nadj(i) =  nadj(i)  + 1
 
                nodeAdj(nadj(i),i) = e
                do j=1,3
 
                    if (i == this%emap(j,e))  then
            
                        e = this%neighborlist(epss(j+1),e)
                        exit
                    endif
                enddo
            
                if (e ==estart) exit
        
            enddo     
       
        endif
    enddo   
 
     
    do ismooth = 1,nsmooth    

        maxdx = 0.
 
        do i = 1,this%nnod
 
            ! if vertex is free
             if ( lfixedNodes(i) /= 1) then
             
                ! Get ring of elements   
                nring = nadj(i)
                ering(1:nring) = nodeAdj(1:nring,i)
 
                !x  = this%vertices(1:2,i)
                x(1)  = this%y(i)
                x(2)  = this%z(i)
                
                xpart(1) = sum( icenters(1,ering(1:nring))  ) / dble(nring)
                xpart(2) = sum( icenters(2,ering(1:nring))  ) / dble(nring)
             
                alpha = 1.d0
                
                do  ! alpha
                
                    xtest = x*(1.d0-alpha) + alpha*xpart
                    
                    do j = 1,nring
                        
                        n     = this%emap(:,ering(j))
!                        verts = this%vertices(1:2,n)
                        verts(1,:) = this%y(n) 
                        verts(2,:) = this%z(n)   
                        do k = 1,3
                            if (n(k) == i) then
                                verts(:,k) = xtest
                                !write(*,*) 'i,j,k:',i,j,k
                                exit
                            endif
                        enddo
                        
                        areatest(j) = getTriArea(verts)  
                       
                    enddo  
                             
                    if (any(areatest(1:nring)<=0)) then
                        alpha = alpha/2d0
!                        write(*,*) 'splitting alpha,i:',i,alpha
                        if (alpha < 1d-3) exit ! all hope is lost, give up trying to make this one better...
                        
                    else

!                        this%vertices(1:2,i) = xtest
                        this%y(i) = xtest(1)
                        this%z(i) = xtest(2)         
                                        
                        dx = sqrt(sum((xtest-x)**2))
                        
                        if (dx > maxdx) maxdx = dx
                        
                        do j = 1,nring
                            
                            e     = ering(j)
                            n     = this%emap(:,e)
!                            verts = this%vertices(1:2,n)
                            verts(1,:) = this%y(n) 
                            verts(2,:) = this%z(n)                 
                            icenters(1:2,e) = sum(verts,2)/3d0
  
                        enddo           
                            
                        
                        exit
                    
                    endif
                    
                enddo  ! alpha
                
        
             endif ! if ( lfixedNodes(i) /= 1) then
 
                
        enddo ! i = 1,this%nnod
        
        if (bTalk) write(*,*) 'ismooth, max dx: ',ismooth,maxdx

        if (abs(lastmaxdx- maxdx)/maxdx < 1d-6) then 
            exit
        else
            lastmaxdx = maxdx
        endif
    enddo ! ismooth = 1,nsmooth    
    
    if (bTalk) write(*,*) '# smoothing iterations, max dx: ',ismooth,maxdx
 
    deallocate (icenters,lfixedNodes,nodeAdj,nadj,node2tri)
 
!    
    end subroutine smooth_TriMesh
!-----------------------------------------------------------------------------------------------------------------------------------        
!--------------------------------------------------------------------------------------------------------------------------- getArea       
!-----------------------------------------------------------------------------------------------------------------------------------     
    real(8) function getTriArea(verts) result(area)

    real(8), intent(in) :: verts(2,3)
    
    real(8) :: acx,bcx,acy,bcy
    integer, parameter          :: eps(6) = (/ 1,2,3,1,2,3 /) 
    acx = verts(1,eps(1)) - verts(1,eps(3))
    bcx = verts(1,eps(2)) - verts(1,eps(3))
    acy = verts(2,eps(1)) - verts(2,eps(3))
    bcy = verts(2,eps(2)) - verts(2,eps(3))
    
    area = (acx * bcy - acy * bcx)/2.d0

    end function getTriArea
    
!======================================================================! 
!=========================================================! check_model
!======================================================================! 
    subroutine check_model(model,minAngle,lHasSlivers)
!
! Checks a trimesh model to see if any segments intersect with 
! small angles (< minAngle). 
!
!  Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
! 
!
!
    use binaryTree
       
    implicit none
    
    ! Arguments:   
    type(trimesh), intent(inout)            :: model
    real(8), intent(in)                     :: minAngle      
    logical, intent(out)                    :: lHasSlivers        
    
    ! Local:
 
    type  :: tree
        type(node), pointer :: start =>null()
    end type tree
    
    type(tree),dimension(:),allocatable :: mytree
    
    
    integer                                 :: i,j, v1,v2, list_length, nSmallAngles
    real(8)                                 :: dy, dz, ang
    character(32)                           :: sortby
    real(8), dimension(:), allocatable      :: values
    integer, dimension(:), allocatable      :: indices
    real(8), parameter                      :: pi = 3.141592653589793d0

    
    lHasSlivers = .false.
    
!---------------------------------------------------------------------- 
! Step 1: Form an adjacency list using an array of linked list 
! (since adj matrix non-zeros are jagged):
!---------------------------------------------------------------------- 

    ! Allocate storage:
    allocate(mytree(model%nnod))  

    !
    ! Loop over all segments and fill in the adjacency list rows:
    !
    do i = 1,model%numberofsegments
    
        v1 = model%segmentlist(2*i-1)
        v2 = model%segmentlist(2*i  )
        
        ! Compute polar angle of segment v1-v2:
        dy = model%y(v2) - model%y(v1)
        dz = model%z(v2) - model%z(v1)
        ang = 180/pi*atan2(dz,dy)

        sortby  = 'value'
        
        ! Add v1 --> v2
        call insert_by_value(mytree(v1)%start,v2, ang)

        ! Add v2 --> v1
        ang = ang+180
        if (ang>=180.d0)  then  ! keep -180 <= ang <= 180
            ang = ang-360
        endif

        call insert_by_value(mytree(v2)%start,v1, ang)            
    
    enddo
!
! Print check:
!    
!    do i = 1,model%nnod
!        write(*,*) 'Node ',i,' is connected to:'
!        call print_tree(  mytree(i)%start)
!        write(*,*) ' '
!    enddo
    
!---------------------------------------------------------------------- 
! Step 2: Loop through each list (i.e. adjacency matrix row) and test 
! polar angles for small inner angles:
!---------------------------------------------------------------------- 
    nSmallAngles = 0
    do i = 1,model%nnod
     
        ! Get values by first getting length, allocating and then retrieving 
        ! them from the linked list:
        list_length =  get_length(mytree(i)%start)
       
        if (list_length==1) cycle  ! ignore unconnected segments
        
        allocate(values(list_length), indices(list_length))
        
        call get_data(mytree(i)%start, indices,values)
        
        ! Now test the values for small angles:
        do j = 2,list_length+1
            if (j == list_length+1) then
                ang = 360- abs(values(1) - values(j-1))
            else
                ang = abs(values(j) - values(j-1))  
            endif
            if ( (ang < minAngle) ) then
                lHasSlivers = .true.
                nSmallAngles = nSmallAngles + 1
                if (nSmallAngles == 1) then
                    write(*,*) '    !!! Warning Small Angles Found !!!  '
                    write(*,*) ' '
                    write(*,fmt='(a6,4x,a6,1x,a12,1x,a12,1x)') '#','Angle', 'x location','y location'
                endif
                write(*,'(i6,4x,f6.1,1x,f12.1,1x,f12.1)') nSmallAngles,ang, model%y(i),model%z(i)
            
            endif
        enddo
        
        deallocate(values,indices)
    
    enddo

!---------------------------------------------------------------------- 
! Deallocate:
!---------------------------------------------------------------------- 

    do i = 1,model%nnod
        call delete_tree(mytree(i)%start)
    enddo
 
    deallocate(mytree)
    
    
    end subroutine check_model
    



    end module  triangle_mesh