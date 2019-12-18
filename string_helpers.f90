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
!==================================================================================================================== string_helpers
!==================================================================================================================================!         
    module string_helpers

!
! This has a bunch of useful functions written by David Myer and myself for decomposing input strings into components.
! 
  
    implicit none
    
    private
    
    public  :: parseCode
    public  :: parseFields
    public  :: parseLine
    public  :: fileParts
    public  :: lower
    public  :: convertStrToInteger
    
    contains
    
!==================================================================================================================================! 
!========================================================================================================================= ParseCode
!==================================================================================================================================! 
    subroutine ParseCode(   sLine, sCode, sValue, bComment )

! David Myer IGPP/SIO La Jolla CA 92093-0225
! Subroutine Revision 3.0, November 2006
! DGM Nov 2006 - parse a line read from a file into a code & value.
! Force the code to be all lowercase with no ending colon.  Terminate
! the line at a '%' or '!' sign (these allow for user comments!)
!
    implicit none
    
    ! Args
  
    character(len=*)            :: sLine
    character(len=*) , intent(out) :: sCode, sValue
    logical, intent(out)        :: bComment
    
    ! Local vars
    integer :: iFrom, iTo, nLen
    
    ! Init returns
    bComment = .false.
    sCode = ' '
    sValue = ' '
    
    nLen = len_trim(sLine)
    
    ! Convert all tab characters to spaces
    forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
    
    ! Skip any beginning blanks
    do iFrom = 1,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    
    ! If the first char is a comment char, then the whole line is a comment.
    if( sLine(iFrom:iFrom) == '%' .or. sLine(iFrom:iFrom) == '!' ) then
        bComment = .true.
        return
    endif
    
    ! Pull off the code value.  Cvt to lowercase as we go.
    iTo = index(sLine,':') - 1
    if (iTo < iFrom) then
       ! write(*,*) 'Parsing Error: missing colon in line below:'
        !write(*,*) sLine
        bComment = .true. ! KWK Feb 2008 treat blank line like a comment 
        return
    endif
    sCode = sLine(iFrom:iTo)
    call Lower(sCode)
    
    ! Skip spaces after the colon
    do iFrom = iTo+2,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    
    ! Get the rest, up to any comment
    sValue = sLine(iFrom:)
    iTo = len_trim(sValue)
    iFrom = index(sValue,'%')
    if (iFrom > 0 .and. iFrom < iTo) then
        sValue(iFrom:iTo) = ' '
    endif
    iFrom = index(sValue,'!')
    if (iFrom > 0 .and. iFrom < iTo) then
        sValue(iFrom:iTo) = ' '
    endif
    
    !call Lower(sValue)     ! No: Some values are filenames which are case-sensitive on UNIX!
    
    end subroutine ParseCode

!==================================================================================================================================! 
!============================================================================================================================= Lower
!==================================================================================================================================! 
    subroutine Lower( s )

! David Myer IGPP/SIO La Jolla CA 92093-0225
! DGM Nov 2006 - convert string to lower case
    implicit none
    character(*), intent(out)  :: s
    integer i

    do  i=1,len_trim(s)
      if  ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) then
        s(i:i) = char(ichar(s(i:i)) + 32)
      endif
    enddo
    
    end subroutine Lower  
    

!==================================================================================================================================! 
!======================================================================================================================= parseFields
!==================================================================================================================================! 
subroutine parseFields(  sLine, nFields, sFields)
!
! Routine to parse out mixed numeric and character fields from a line
! This is useful for mixed format input tables, which are awkward for Fortran
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! Version 1.0.   November 21, 2008. 
!
    implicit none
    
    ! Args

    character(len=*)        :: sLine
    integer, intent(in)     :: nFields
    character(len=*) , intent(out) :: sFields(nFields)
    
    
   
    
    ! Local vars
    integer :: iFrom, iTo, i, nLen
    
    nLen = len_trim(sLine)
    
    ! Convert all tab characters to spaces
    forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
    
    iFrom = 1
    
    ! Loop through the line and get the nFields:
    do i=1,nFields
    
        ! Skip any beginning blanks
        do iFrom = iFrom,nLen
            if (sLine(iFrom:iFrom) .ne. ' ') exit
        enddo
        
        ! Pull out nonblank character string:
        do iTo = iFrom,nLen
             if (sLine(iTo:iTo) .eq. ' ') exit
        enddo
        sFields(i) = trim(sLine(iFrom:iTo-1))
        iFrom = iTo
    enddo  
    
    end subroutine parseFields     
    
!==================================================================================================================================! 
!========================================================================================================================= fileParts
!==================================================================================================================================! 
    subroutine fileParts(   sStr, sRoot, sExt1, sExt2 )
!
! Like the MATLAB function fileparts.
!    
    implicit none
    
    character(len=*),intent(in)  :: sStr 
    character(len=*),intent(out) :: sRoot
    character(len=*),intent(out), optional :: sExt1 , sExt2 
    
     
    integer :: idot
     
    sRoot = adjustl(sStr)
  
    idot = index(sRoot,'.',.true.)
    
    sExt1 = ' '
    sExt2 = ' '

    if (idot > 0) then 
        
            sExt1 = sRoot(idot+1:)
            sExt1 = adjustl(sExt1)
            sRoot = sRoot(1:idot-1)
            idot = index(sRoot,'.',.true.)
            if (idot == 0) then
                sExt2 = ''
            else
                sExt2 = sRoot(idot+1:)
                sExt2 = adjustl(sExt2)
                sRoot = sRoot(1:idot-1)
            endif 
    
  
    endif 
    
    end  subroutine fileParts 
    
!==================================================================================================================================! 
!========================================================================================================================= parseLine
!==================================================================================================================================! 
    subroutine parseLine( nLen, sLine, bComment )
    !   
    ! Subroutine to check if the sLine is blank or a comment line, and if it isn't
    ! then any comment at the end of the line is blanked. 
    !
    ! This is styled after D. Myer's parseCode function.
    !
    ! This is useful for reading in data tables that have user comment lines 
    ! or comments at the end of the line, denoted by the ! and % symbols.
    ! 
    ! If the entire line is a comment, bComment = .true.
    !
    ! Kerry Key
    ! Scripps Institution of Oceanography
    ! kkey@ucsd.edu
    !
    ! Version 1.0.   April, 2008.
    !
 
    
    ! Args
    integer, intent(in)     :: nLen
    character(nLen)         :: sLine
    logical, intent(out)    :: bComment
    
    ! Local vars
    integer :: iFrom, iTo
    
    ! Init returns
    bComment = .false.
    
    
    ! Convert all tab characters to spaces
    forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
    
    ! Skip any beginning blanks
    do iFrom = 1,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    
    ! If the first char is a comment char, then the whole line is a comment.
    ! DGM April 2008 Also, if the line is blank, consider it a comment.
    if( iFrom > nLen .or. sLine(iFrom:iFrom) == '%' &
        .or. sLine(iFrom:iFrom) == '!' ) then
        bComment = .true.
        return
    endif
    
    ! Now trim off any comments at the end of the line  
    iTo = len_trim(sLine)
    iFrom = index(sLine,'%')
    if (iFrom > 0 .and. iFrom < iTo) then
        sLine(iFrom:iTo) = ' '
    endif
    iFrom = index(sLine,'!')
    if (iFrom > 0 .and. iFrom < iTo) then
        sLine(iFrom:iTo) = ' '
    endif
    
    end subroutine parseLine    
    

!==================================================================================================================================! 
!=============================================================================================================== convertStrToInteger
!==================================================================================================================================! 
    subroutine convertStrToInteger(string,iVal,lokay)
    
    implicit none 
    
    character(len=*)            :: string
    character(len=len(string))  :: first
    character(len=20)           :: form
    integer                     :: ierr
    integer                     :: iVal
    logical                     :: lokay

    lokay = .false.
    read(string, *, iostat = ierr ) first

    if ( ierr == 0 ) then
        write( form, '(a,i0,a)' ) '(i', len(string), ')'
        read( first, form, iostat = ierr ) iVal
    endif

    lokay = ierr == 0   
    
    end subroutine convertStrToInteger
        
    end module string_helpers