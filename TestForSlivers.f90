!
! A program to test out the Fortran-C interface for Triangle.c.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!

program TestForSlivers
	
	use triangle_mesh
	
	implicit none
	
	type(trimesh) :: in 
	
	character (24) :: tristr,cend	
	character (256) :: fileroot 
	
	logical :: lHasSlivers
	real(8) :: minQangle 
	
	!
	! Get name of .poly file to read:
	!
	write(*,*) 'Enter base PSLG file name (without .poly):'
	read(*,*)  fileroot

	!
	! Get quality angle
	!
	write(*,*) 'Enter triangulation quality angle (minimum inner angle):'
	read(*,*)  minQangle
	minQangle = abs(minQangle)
	
	!
	! Read in the poly file:
	!
	write(*,*) ' '
	call read_triangle(fileroot,in)  ! read poly file into structure in
    write(*,*) ' ' 

    !
    ! Check the model for slivers:
    !
	call check_model(in,(abs(minQangle)),lHasSlivers)
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
	    write(*,*) '!!!    >>> Better luck next time, bucko <<<          !!!'
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
        write(*,*) '  A-okay buddy:'
        write(*,*) ' '
	    write(*,*) ' Your model has been checked for slivers and  '
	    write(*,*) ' none were found.  '
    endif
    write(*,*) ' '
    

end program TestForSlivers

 
!    
	