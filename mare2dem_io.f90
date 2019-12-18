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

    module mare2dem_io
    
    use string_helpers
    use em2dkx_mod
    use mare2dem_global
    use mare2dem_input_data_params
    use Occam 
         
    implicit none
    
    private
    
    public :: readModel
    public :: readData
    !public :: writeResponse
    public :: writeJacobian
    public :: writeResistivityFile
    public :: writeFields_MT
    public :: writeFields_CSEM
    public :: print2col32
    
    contains
    
!==================================================================================================================================! 
!=============================================================================================================== readResistivityFile
!==================================================================================================================================!       
    subroutine readResistivityFile
! 
! Reads in the resistivity file that specifies all resistivity regions
! for all fixed and free parameters.
! 
 
 
    integer             :: ierr, iLine, i, j, iAllocErr, iostat_err, iskip, ict, iFmt, iparam1, iparam2
    character(256)      :: filename, vers
    character(256)      :: sLine, sCode, sValue 
    character(32)       :: cnum
    logical             :: bComment
    real(RealPrec)      :: pmtol, rtemp 
    real(RealPrec)      :: inputMisfit, inputRoughness ! only locally read in and displayed here since Occam will recompute them...
 
    
    !
    ! Open the file:
    !
    filename = trim(resistivityFile)//'.resistivity' 
    
    write(6,*) '=============== Reading in Resistivity File =============================='
    write(*,*) ' '            
    call print2col32('Opening file: ',filename,6)
    open(unit=10,file=trim(filename),status='old',iostat = ierr)
    if (ierr .ne. 0) then
        write(*,*) ' Error opening Resistivity File: ',trim(filename)
        call exitMARE2DEM
    end if  
    write(*,*) ' '  
    
    nRhoPerRegion  = 0
    lowerBoundGlobal = 0
    upperBoundGlobal = 0
    
    ! Read in the file a line at a time and decode the sCode/sValue pairs that are separated by a semicolon
    
    iLine = 0
    do  ! infinite while loop
        
        ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        
        read( 10, '(A)', iostat = ierr ) sLine
        
        if (ierr /= 0) exit  ! end of file read, escape from while loop and 
                        ! proceed to checking that required inputs defined
        
        iLine = iLine + 1   ! advance line counter
        
        call ParseCode(  sLine, sCode, sValue, bComment )
        if( bComment ) cycle
        
        !
        ! What do we have?
        !
        select case (trim(sCode))
        
        case ('version','format')
        
            vers = sValue(1:len(vers))    
            call Lower(vers)
            select case (trim(vers))
            
            case ('mare2dem_1.0')
                iFmt = 1
            case ('mare2dem_1.1')   
                iFmt = 2
            case default
           
                write(*,*) 'Error, unsupported resistivity file version: ',trim(vers)
                write(*,*) ' Try using mare2dem_1.1 format.   Stopping!'
                call exitMARE2DEM
            end select
                         
        case ('model file')   
        
            modelFileName = sValue(1:len(modelFileName))         
            call print2col32('Model File: ', modelFileName,6)
        
        case ('data file')   
        
            dataFileName = sValue(1:len(dataFileName))         
            call print2col32('Data File: ', dataFileName,6)
        
        case ('settings file')   
        
            settingsFileName = sValue(1:len(settingsFileName))         
            call print2col32('Settings File: ', settingsFileName,6)
                    
                    
        case ('penalty file')   
        
            penaltyFileName = sValue(1:len(penaltyFileName))         
            call print2col32('Penalty File: ', penaltyFileName,6)
        
        case ('maximum iterations')   
            read(sValue,*) maxOccamIterations
        
        case ('bounds transform')
            cBoundsTransform = trim(sValue)              
            
            ! Check to make sure input is valid:
            !
            select case (cBoundsTransform)
            case ('exponential','bandpass')
                call print2col32('Model bounds transform: ', cBoundsTransform,6) 
            case default 
                write(*,*) '  Error, unrecognized model bounds transform:  ', trim(cBoundsTransform)
                call exitMARE2DEM
            end select
        
        case ('global bounds')    
            ! Limit value string must have two values: min, max
            ! Comma is required.
            i = index( sValue, ',' )
            
            if (i == 0) then
                write(*,*) 'Invalid "global Bounds" entry in resistivity file.'
                write(*,*) 'The format for this line is'
                write(*,*) '   Global Bounds: <min>,<max>'
                write(*,*) 'Do not include angle brackets.  Comma required.'
                write(*,*) 'Min must be less than max.'

                call exitMARE2DEM
            endif
             
            read(sValue(:i-1),*) lowerBoundGlobal
            
            ! Skip space padding
            do i=i+1,len_trim(sValue)
                if (sValue(i:i) .ne. ' ') exit
            enddo
            read(sValue(i:),*) upperBoundGlobal
            
            ! Check for mistake:
            if (lowerBoundGlobal > upperBoundGlobal) then
                rtemp = lowerBoundGlobal
                lowerBoundGlobal = upperBoundGlobal
                upperBoundGlobal = rtemp 
            endif
            
            write(cnum,'(es11.4,1x,es11.4)') lowerBoundGlobal, upperBoundGlobal
            call print2col32(' Global bounds: ',cnum,6)
        
            
        case ('debug level', 'print level')
            read(sValue,*) occamPrintLevel
            write(cnum,'(i3)') occamPrintLevel
            call print2col32('Print Level: ',cnum,6)
                
        case ('target misfit')
            read(sValue,*) targetRMS
            write(cnum,'(es11.4)') targetRMS
            call print2col32('Target Misfit: ',cnum,6)       
            
        case ('iteration')    
        ! kwk: iteration number is now in file name, ignore any input numbers
 
        case  ('misfit decrease threshold')
            read(sValue,*) rmsThreshold
            write(cnum,'(g11.4)') rmsThreshold     
            call print2col32('Misfit Decrease Threshold: ', cnum,6)           
 
        case  ('converge slowly')
            call Lower(sValue)
            select case (trim(sValue))
            case ('yes')
                  lConvergeSlowly = .true.
            end select  
            call print2col32('Converge Slowly: ', trim(sValue),6)
            

        case ('lagrange value')  
            if (len_trim(sValue) > 0 ) then
                read(sValue,*) modelMu  
                write(cnum,'(es11.4)') modelMu
                call print2col32('Lagrange Value: ',cnum,6)
            endif
            
        case ('model roughness') 
            if (len_trim(sValue) > 0 ) then
                read(sValue,*) inputRoughness 
                write(cnum,'(es11.4)') inputRoughness
                call print2col32('Model Roughness: ',cnum,6)
            endif
        
        case ('model misfit')
            if (len_trim(sValue) > 0 ) then
                read(sValue,*) inputMisfit
                write(cnum,'(es11.4)') inputMisfit
                call print2col32('Model Misfit: ',cnum,6)
            endif            
        case  ('date/time')
            ! skip
            
!        case  ('inversion method')
!            call Lower(sValue)
!            sInversionMethod = sValue(1:len(sInversionMethod))         
!            call print2col32('Inversion Method: ', sInversionMethod,6)
            
!        case  ('fixed mu cut')
!            read(sValue,*) rMuCut
!            write(cnum,'(es11.4)') rMuCut     
!            call print2col32('Fixed Mu Cut: ', cnum,6)                        
!        
        case  ('anisotropy') 
            cAnisotropy = sValue(1:len(cAnisotropy))  
            ! Check to make sure this makes sense:    
            select case (trim(cAnisotropy ))
                    
            case ('isotropic')
                nRhoPerRegion = 1
            case ('triaxial')
                nRhoPerRegion = 3
            case ('tix','tiy','tiz')
                nRhoPerRegion = 2
            
            case default
                write(*,*) ' Error decoding anisotropy in resistivity file'
                write(*,*) ' Unknown or unsupported anisotropy code:', trim(cAnisotropy )
                write(*,*) ' Stopping!'
                call exitMARE2DEM
            end select    
            
            call print2col32('Anisotropy: ',cAnisotropy,6)
          
!    
        case  ('number of regions')     
            read(sValue,*) nRegions
            write(cnum,'(i8)') nRegions
            call print2col32('Number of regions: ',cnum,6)
                            
            nSigParams =   nRegions*nRhoPerRegion             
            if (nSigParams < 1)  then
                write(*,*) ' Error in resistivity file, nRegions, nRho per region: ',nRegions, nRhoPerRegion
                write(*,*) ' Stopping!'
                call exitMARE2DEM    
            endif
            
          !
          ! Now read in the regions:
          !
        
          !
          ! First allocate storage:
          !
            allocate(sigParams(nSigParams), iFreeParam(nSigParams), stat=iAllocErr)
            if (iAllocErr .ne. 0) then
                write(*,*) ' Out of memory.  Too many resistivity parameters (', nSigParams, ')'
                call exitMARE2DEM 
            endif      
            
            allocate (boundsTemp(nRegions,2*nRhoPerRegion), &
                      PrejTemp(nRegions,2*nRhoPerRegion), stat=iAllocErr)   
            if (iAllocErr .ne. 0) then
                write(*,*) ' Out of memory.  Too many resistivity parameters (', nRegions*nRhoPerRegion, ')'
                call exitMARE2DEM 
            endif
            
            if ( nRhoPerRegion > 1 ) then 
                allocate (ratioTemp(nRegions,nRhoPerRegion*(nRhoPerRegion-1)), stat=iAllocErr)  ! 2 or 6 columns
                if (iAllocErr .ne. 0) then
                    write(*,*) ' Out of memory.  Too many resistivity parameters (', nRegions*nRhoPerRegion, ')'
                    call exitMARE2DEM 
                endif            
                ratioTemp = 0 
            endif
                                   
            !
            ! Read in each line:
            !  
            ! Skip comment line:
            read( 10, '(A)', iostat = ierr ) cParamHeader
            
            if ( (iFmt == 1)  .or. (nRhoPerRegion == 1) )then ! Format: ', ' MARE2DEM_1.0
                do i = 1,nRegions
                    read(10,*,iostat = iostat_err) iskip, ( sigParams((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  & ! sig is linear
                    &                                      (iFreeParam((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),    &
                    &                                      (boundsTemp(i,j),j=1,2*nRhoPerRegion),  &
                    &                                      (  prejTemp(i,j),j=1,2*nRhoPerRegion) 
 
                    if (iostat_err /= 0) then
                        write(*,*) ' Error reading parameters, stopping!'
                        call exitMARE2DEM
                    endif        
                enddo
            
            else  ! nRhoPerRegion > 1 and new format with anisotropy ratio support:
            
                do i = 1,nRegions
                    read(10,*,iostat = iostat_err) iskip, ( sigParams((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  & ! sig is linear
                    &                                      (iFreeParam((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),    &
                    &                                      (boundsTemp(i,j),j=1,2*nRhoPerRegion),  &
                    &                                      (  prejTemp(i,j),j=1,2*nRhoPerRegion), &
                    &                                      (  ratioTemp(i,j),j=1,nRhoPerRegion*(nRhoPerRegion-1))     
 
                    if (iostat_err /= 0) then
                        write(*,*) ' Error reading parameters, stopping!'
                        call exitMARE2DEM
                    endif        
                enddo
                            
            endif
                        
        case default
            write(*,*) 'Error reading resistivity file!'
            write(*,*) ' On line :', iLine            
            write(*,*) 'Unknown or unsupported code:'
            write(*,*) sCode
        
            call exitMARE2DEM
        
        end select 
        
    enddo ! read while loop
    
    write(*,*) ' ' 

!
! Close resistivity file
!
    close(10)
    
 !
 ! Perform a few odds and ends here:
 !    
    
    nFree = count(iFreeParam > 0)
    nParams = nFree
    
    write(*,*) 'Done reading resistivity file, here is what I found: '
    write(*,*) ' '
    write(cnum,'(i8)') nFree
    call print2col32('# Free parameters: ',cnum,6)
    write(cnum,'(i8)') nSigParams-nFree
    call print2col32('# Fixed parameters: ',cnum,6)
    write(*,*) ' '
    
    
    !
    ! Now fill in the free parameter arrays
    !
    allocate( pm(nFree), lowerBound(nFree),upperBound(nFree),lBoundMe(nFree),premod(nFree),prewts(nFree), stat=iAllocErr)
    if (iAllocErr .ne. 0) then
        write(*,*) ' Out of memory.  Too many free parameters (', nFree, ')'
        call exitMARE2DEM 
    endif      
    
    ! Initialize these to the globals, then override if local value is nonzero:
    ! Also don't forget that everything needs to be converted to log10 here:
    
    lowerBound = 0 
    upperBound = 0 
    lBoundMe = .false.
    
    if ( (lowerBoundGlobal > 0) .and. (upperBoundGlobal > 0) ) then
        lowerBound = log10(lowerBoundGlobal)
        upperBound = log10(upperBoundGlobal)
        lBoundMe         = .true. 
    endif 
    
    premod     = 0
    prewts     = 0           
    
    ict = 0
    do i = 1,nRegions
        do j = 1,nRhoPerRegion
            
            if (iFreeParam((i-1)*nRhoPerRegion+j)>0) then
                
                ict = ict + 1
                
                if ( ( boundsTemp(i,2*j-1) > 0 ) .and. ( boundsTemp(i,2*j  ) > 0 ) ) then
                    lowerBound(ict) = log10(boundsTemp(i,2*j-1)) ! (i,1 3 5)
                    upperBound(ict) = log10(boundsTemp(i,2*j  )) ! (i,2 4 6)    
                    lBoundMe(ict)   = .true.
                endif      
                
                pm(ict)     = log10(sigParams((i-1)*nRhoPerRegion+j)) ! inversion for log10(rho)
                
                prewts(ict) =     abs(  prejTemp(i,2*j  ) ) ! (i,2 4 6) ! weights are linear      
                if (   prewts(ict) > 0) then
                    premod(ict) = log10(prejTemp(i,2*j-1)) ! (i,1 3 5)
                endif                
                !write(*,*) ict, pm(ict),premod(ict),prewts(ict), lowerBound(ict),upperBound(ict),lBoundMe(ict)       
            endif
        enddo
    enddo
    
    if (ict /= nFree) write(*,*) 'Error in readResistivityFile, ict /=nFree: ', ict, nFree
    ! Convert input rho to sigma since that's what em2dkx uses
    sigParams = 1./sigParams 
    
    
    !
    ! If model is anisotropic, create parameter difference preference array (log of ratio), if any
    !
    nDiff = 0
    if ( nRhoPerRegion > 1 ) then 
        
        allocate(preDiff(nFree),preDiffwts(nFree),ijDiff(nFree,2)) ! over-allocating, only up to nDiff values will be set
        preDiffwts = 0
        preDiff    = 0
        
        ict = 0
        do i = 1,nRegions
            
            do j = 1,nRhoPerRegion

                if (iFreeParam((i-1)*nRhoPerRegion+j)>0) then
    
                    ict = ict + 1 ! free parameter counter
                    
                    if ( (iFreeParam((i-1)*nRhoPerRegion+j)>0).and.(iFreeParam((i-1)*nRhoPerRegion+eps(j+1))>0) ) then
                    ! if this parameter and next one in this region are free:
                        
                        if ( (nRhoPerRegion==2) .and. (j>1) )  then
                        
                         ! do nothing
                         
                         elseif  (  (ratioTemp(i,2*j) > 0) .and. (ratioTemp(i,2*j-1) > 0 ) ) then ! add a ratio penalty for two anisotropic components:
                         
                        
                            nDiff = nDiff + 1
                    
                            iparam1 = ict
                            iparam2 = ict+1
                            if (j==3) iparam2 = ict - 2  ! x/y, y/z, and z/x, so skip back to 1 if j==3
                    
                            preDiffwts(nDiff) = ratioTemp(i,2*j  ) 
                            preDiff(nDiff)    = log10(ratioTemp(i,2*j-1))
                            ijDiff(nDiff,1)   = iparam1
                            ijDiff(nDiff,2)   = iparam2 
                        
                           ! write(*,'(4(i8,1x),2(g12.4,1x))') j,nDiff,iparam1,iparam2,preDiff(nDiff),preDiffwts(nDiff) 
                        
                        endif                        
                        
                    endif    
                endif  
            enddo ! j
        enddo ! i       
        
    endif
    
!
! Double check the input pm to make sure its within the bounds, then transform it to the unbound parameters
!      
    do i = 1,nFree
       
        if (lBoundMe(i)) then
        
            pmtol = (upperBound(i) - lowerBound(i)) / 232.
            if (pm(i) <= lowerBound(i) + pmtol) pm(i) = lowerBound(i) + pmtol
            if (pm(i) >= upperBound(i) - pmtol) pm(i) = upperBound(i) - pmtol
!       
            ! now transform it:
            pm(i)  =  transformToUnbound(pm(i) ,lowerBound(i) ,upperBound(i),lBoundMe(i) ) 
       
            
        endif   
      
    enddo 
 
    end subroutine readResistivityFile           
!==================================================================================================================================! 
!============================================================================================================== writeResistivityFile
!==================================================================================================================================!       
    subroutine writeResistivityFile(nCurrentIter,misfit,misfitRequested,lagrange,roughness)
 
! 
! Reads in the resistivity file that specifies all resistivity regions
! for all fixed and free parameters.
! 
 
    integer, intent(in)        :: nCurrentIter
    real(RealPrec), intent(in) :: misfit,lagrange,roughness,misfitRequested 
    
    real(RealPrec) :: logrho 
    character(256) :: cNum, sFmt 
    integer        :: lerr, ioUnit = 21, i, ict, j
    character(80)  :: dateAndTime

!    
! Open resistivity file
!
    write (cNum,*) nCurrentIter
    open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.resistivity', iostat=lerr)
  
    call print2col32('Format: ', ' MARE2DEM_1.1',ioUnit)
    call print2col32('Model File: ', modelFileName,ioUnit)
    call print2col32('Data File: ', dataFileName,ioUnit)
    call print2col32('Settings File: ', settingsFileName,ioUnit)
    call print2col32('Penalty File: ', penaltyFileName,ioUnit)
    write (cNum,'(i6)') maxOccamIterations
    call print2col32('Maximum Iterations: ',cnum ,ioUnit)
    call print2col32('Bounds Transform: ', cBoundsTransform,ioUnit)
    write(cnum,'(es12.4,1a,1x,es12.4)') lowerBoundGlobal,',', upperBoundGlobal
    call print2col32('Global Bounds: ', cnum,ioUnit)
    write (cNum,'(i6)') occamPrintLevel
    call print2col32('Print Level: ', cnum,ioUnit)
    write(cnum,'(g12.4)') misfitRequested !
    call print2col32('Target Misfit: ',cnum,ioUnit) 
    
 
    write(cnum,'(g12.4)') rmsThreshold !
    call print2col32('Misfit Decrease Threshold: ',cnum,ioUnit)        
    if (lConvergeSlowly) then
        cnum = 'yes'
    else
        cnum = 'no'
    endif
    call print2col32('Converge Slowly: ',cnum,ioUnit)        
     
    write(cnum,'(e12.4)') lagrange
    call print2col32('Lagrange Value: ',cnum,ioUnit)
    write(cnum,'(g12.4)') roughness
    call print2col32('Model Roughness: ',cnum,ioUnit)    
    write(cnum,'(f11.4)') misfit
    call print2col32('Model Misfit: ',cnum,ioUnit)
    call datime(dateAndTime) 
    call print2col32('Date/Time: ',dateAndTime,ioUnit)
    
    ! kwk debug: this is currently a test feature and might not ever be a supported feature, so don't output 'fixed' flag unless
    !            it is being used:
!    if (trim(sInversionMethod) =='fixed') then
!        call print2col32('Inversion Method: ', sInversionMethod,ioUnit)
!        write(cnum,'(g12.4)')  rMuCut
!        call print2col32('Fixed Mu Cut: ', cnum,ioUnit)
!    endif
    
    call print2col32('Anisotropy: ',cAnisotropy,ioUnit) 
    write(cnum,'(i8)') nRegions  
    call print2col32('Number of regions: ',cnum,ioUnit)  
    
    ! First we need to copy pm into sigParams:
    ict = 0
    do i=1,nRegions
        do j = 1,nRhoPerRegion  
            if (iFreeParam((i-1)*nRhoPerRegion + j) > 0 ) then
                ict = ict+1
                logrho  =  transformToBound(pm(ict),lowerBound(ict),upperBound(ict),lBoundMe(ict) )   
                sigParams((i-1)*nRhoPerRegion + j) =  1./(10.d0**logrho)
           endif
        enddo    
    enddo 
      
    ! Now loop over the regions and write out a line for each:  
     sFmt = '(i8,1x'
     do j=1,nRhoPerRegion ! rho
        sFmt = trim(sFmt)//',es12.4,1x'
     enddo
     do j=1,nRhoPerRegion ! free
        sFmt = trim(sFmt)//',i8,1x'
     enddo   
     do j=1,2*2*nRhoPerRegion ! bounds, prej
        sFmt = trim(sFmt)//',es12.4,1x'
     enddo
     if (nRhoPerRegion > 1) then
         do j=1,nRhoPerRegion*(nRhoPerRegion-1) ! ratio prejudice
            sFmt = trim(sFmt)//',es12.4,1x'
         enddo   
     endif
     sFmt = trim(sFmt)//')'   
      
!     write(*,*) sFmt
     write(ioUnit,'(a)') trim(adjustl(cParamHeader))
     do i = 1,nRegions
    
        if (nRhoPerRegion == 1) then
            write(ioUnit,sFmt)  i, (1./sigParams((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &  
                    &              (  iFreeParam((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &
                    &              (   boundsTemp(i,j),j=1,2*nRhoPerRegion),  &
                    &              (     prejTemp(i,j),j=1,2*nRhoPerRegion) 
        else
            write(ioUnit,sFmt)  i, (1./sigParams((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &  
                &                  (  iFreeParam((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &
                &                  (   boundsTemp(i,j),j=1,2*nRhoPerRegion),  &
                &                  (     prejTemp(i,j),j=1,2*nRhoPerRegion),  &                
                &                  (    ratioTemp(i,j),j=1,nRhoPerRegion*(nRhoPerRegion-1) ) 
                
        endif            
     enddo  
      
! Close the file:    
    close(21)
    
    end subroutine writeResistivityFile    
    
!==================================================================================================================================! 
!============================================================================================================================ datime
!==================================================================================================================================!       
    subroutine datime(datetm)
    
    character(80) datetm
    character(10) cdate,ctime
    
    !
    ! Standard intrinsic Fortran90 data_and_time call:
    !
    call date_and_time(cdate,ctime) !this ouputs only values, ignoring other optional arguments
    
    datetm = cdate(5:6)//'/'//cdate(7:8)//'/'//cdate(1:4)//' '// ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:)               
    
    end subroutine datime  
    

  
!==================================================================================================================================! 
!========================================================================================================================== readPoly
!==================================================================================================================================!       
    subroutine readPoly
! 
! Reads in the .poly file used to describe the boundary model.
! October, 2011: retooled to use routines in call_triangle.f90 for most
! file reading, just does some minor adjustments to input model here
! 
 

    use triangle_mesh   
    use mare2dem_global, only: modelFileName,inputmodel
    
    implicit none
 
    character(256)      :: filename
       
    integer             :: idot
!
! Read in the .poly file (and optionally a .nod and .ele file if they exist):
!
 
    write(*,*) '======== Reading in Poly file ============================================'
    write(*,*) ' '     
 
    modelFileName = adjustl(modelFileName)
    idot = index(modelFileName,'.',.true.)
    filename = modelFileName(1:idot-1)
 
    call read_triangle(filename,inputmodel)
    write(*,*) 'Done reading Poly file'
    write(*,*) ' '
              
    end subroutine readPoly   
    
!==================================================================================================================================! 
!======================================================================================================================= readPenalty
!==================================================================================================================================!        
    subroutine readPenalty
 
    
    integer             :: ierr, i
    character(8)        :: cnum
    character(256)      :: sLine, sCode, sValue 
    logical             :: bComment 
    
    write(*,*) '======== Reading in Penalty file ========================================='
    write(*,*) ' '     

    call print2col32(' Opening Penalty file: ',penaltyFileName,6)
    open(unit=10,file=trim(penaltyFileName),status='old',iostat = ierr)
    
    read( 10, '(A)' ) sLine
    call ParseCode(  sLine, sCode, sValue, bComment )
    select case (trim(sCode))
    
    case ( 'format') ! only one 'format' for now, so no need to check the sValue
        
        call print2col32(' reading new CSR format... ','',6)
                
        read( 10, * ) pen%nnz, pen%nrows
       ! write( *, * ) pen%nnz, pen%nrows
        
        write(cnum,'(i8)') pen%nnz
        call print2col32('Number of penalties: ',cnum,6)
        
        allocate(pen%colind(pen%nnz), pen%val(pen%nnz), pen%rowptr(pen%nrows+1) ) 
         
        do i = 1,pen%nnz
            read( 10, * ) pen%colind(i), pen%val(i)
            !write(*,*) pen%colind(i), pen%val(i)
        enddo
        
        do i = 1,pen%nrows+1
            read( 10, * ) pen%rowptr(i)
            !write(*,*) pen%rowptr(i)
        enddo        
    
    case default ! old style penalty file:  npen \n [i,j,weight]
    
        read (sLine,*) npenalty
         
        write(cnum,'(i8)') npenalty
        call print2col32('Number of penalties: ',cnum,6)
        
        allocate(ijpenalty(npenalty,2),penaltywts(npenalty)) 
        do i = 1,npenalty
            read(10,*) ijpenalty(i,1:2),penaltywts(i)
            !write(*,*) i,ijpenalty(i,1:2),penaltywts(i)
        enddo
   
        
        if (any(ijpenalty > nFree)) then
            write(*,*) ' Error, penalty matrix points to parameter numbers larger than the input free parameters!'
            write(*,*) ' Largest index: ',maxval(ijpenalty),', Number of free paramters: ', nFree
            write(*,*) ' Stopping!'
            call exitMARE2DEM
        endif
        
    end select
    
    close(10)
   
    write(*,*) ' '
    write(*,*) 'Done reading penalty file '
    write(*,*) ' '
    
    end subroutine readPenalty
    
!==================================================================================================================================! 
!====================================================================================================================== readSettings
!==================================================================================================================================!      
    subroutine readSettings  
 
    
    integer             :: err, iLine = 0
    character(256)      :: sLine, sCode, sValue 
    logical             :: bComment
    
    integer             :: nRxPerGroup 
    
 
    write(6,*) '========== Reading in Parallel and Adaptive Refinement Settings =========='
    write(*,*) ' '     
    call print2col32('Reading settings file:  ',settingsFileName,6)
    write(*,*) ' '
        
    open (unit=10,file=settingsFileName,status='old',iostat=err)
    if (err .ne. 0) then
        write(*,*) ' Error opening settings file'
        call exitMARE2DEM
    end if
    
! Read in the file a line at a time and decode the sCode/sValue pairs that are separated by a semicolon
! Certain sCodes are followed by lists such as model parameters, site locations, etc.
! Parsecode is D Myer's useful code tidebit for separating sCode and sValue pairs
 
    do  ! infinite while loop
      
      ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        
        read( 10, '(A)', iostat = err ) sLine
        
        if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
        iLine = iLine + 1   ! advance line counter
        
        call ParseCode(  sLine, sCode, sValue, bComment )
        if( bComment ) cycle
          
        !
        ! What do we have?
        !
        select case (trim(sCode))
       
        case ('version')
        
        case ('scratch folder','scratchfolder')    
            
            scratchFolder = trim(adjustl(sValue))
                                  
        case ('wavenumbers','csem wavenumbers')    
            
            read(sValue,*)  loglower, logupper, nwave 
   
        case ('max # refinement','max # refinements') 
        
            read(sValue,*) maxnadapt
        
        case ('tolerance (%)','tolerance')
        
            read(sValue,*) errortolerance        

        case ('mesh quality angle')
        
            read(sValue,*) minQangle
            
            if (minQangle > 33) then
                write(6,fmt='(a32,f6.3)') 'Mesh quality angle:  ',minQangle  
                write(*,*) ' '
                write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write(*,*) ' !!!! Error: quality angle should be 10 < minQangle < 33 degrees !!!!'
                write(*,*) ' minimum angles > 33 degrees are nearly impossible to mesh ! '
                write(*,*) ' while angles < 10 degrees lead to numerical problems with the FE method '
                write(*,*) ' Sorry buddy, but try specifying a better angle!'
                write(*,*) ' ' 
                write(*,*) ' STOPPING!  '
                write(*,*) ' ' 
                write(*,*) ' ' 
                write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                call exitMARE2DEM
            endif           
 
        case ('save meshes')
        select case (trim(sValue))
            case ('yes')
              lSaveMeshFiles = .true.
            end select  

        case ('dual function')

        read(sValue,*) idual_func        

        !
        ! Grouping parameters for parallelization:
        !
        case ('transmitters per group')   
        read(sValue,*) nTxPerGroup 

        case ('receivers per group')   !  here for backwards compatibility       
        read(sValue,*) nRxPerGroup 
        nRxPerGroupMT   = nRxPerGroup
        nRxPerGroupCSEM = nRxPerGroup

        case ('mt receivers per group')   
        read(sValue,*) nRxPerGroupMT

        case ('csem receivers per group')    
        read(sValue,*) nRxPerGroupCSEM
    
        case ('csem frequencies per group')         
        read(sValue,*) nFreqPerGroupCSEM 

        case ('wavenumbers per group')         
        read(sValue,*) nKxPerGroup 

        case ('mt frequencies per group')         
        read(sValue,*) nFreqPerGroupMT 


        case (' (Co)Sine Transform Filters','ct filters','ct filter','ctfilters','ctfilter')

        !           FCTfilter = sValue(1:len(FCTfilter))  ! sorry, no longer supporting this. Filters are fixed to be 601 pts...  

        case ('solver','linear solver','linearsolver')
        linearSolver =  sValue(1:len(linearSolver))  

        case ('print debug','printdebug','debug')  
        select case (trim(sValue))
        case ('yes')
            lprintDebug = .true.
            lprintDebug_em2dkx = .true.
        end select 
        case ('print setup','printsetup')  
        select case (trim(sValue))
        case ('yes')
            lprintSetup = .true.
        end select 
        case ('print decomposition')  
        select case (trim(sValue))
        case ('yes')
            lprintDecomposition = .true.
        end select                                   
        case ('use bump fields')  
        select case (trim(sValue))
        case ('yes')
            lUseBumpFields = .true.
        case('no')
            lUseBumpFields = .false.                
        end select  

        case ('print adaptive')  
        select case (trim(sValue))
        case ('yes')
            lDisplayRefinementStats = .true.
        case('no')
            lDisplayRefinementStats = .false.
        end select             
        
        !
        ! Items below are for fine tuning the error estimator and adaptive refinement: advanced users only!
        !    
        case ('minimum error range')  
            read(sValue,*) minRangeProtector  
 
        case ('minimum area')
            read(sValue,*) minArea
             
        case ('max # subrefinements')
            read(sValue,*) max_nsubrefine 
            
         case ('percent refine')
            read(sValue,*) pct_refine 
            
        case ('e noise floor')  
            read(sValue,*) ecutoff 
               
        case ('h noise floor')  
            read(sValue,*) hcutoff
                          
        case ('transmitter quadrature order')
             read(sValue,*) nQuadTxCSEM  
             
        case ('mt quadrature order','mt receiver quadrature order') 
            read(sValue,*) nQuadRxMT  
            
        case ('csem quadrature order','csem receiver quadrature order') 
            read(sValue,*) nQuadRxCSEM  

 
                                 
        case default
            write(*,*) 'Error reading RUNFILE file!'
            write(*,*) ' On line :', iLine            
            write(*,*) 'Unknown or unsupported code:'
            write(*,*) sCode

            call exitMARE2DEM
        
        end select 
        
    enddo ! read while loop
!
! Close settings file:
!
    close(10)
    
    write(*,*) 'Done reading settings file '
    write(*,*) ' '
 
    end subroutine readSettings
    
!==================================================================================================================================! 
!========================================================================================================================= readModel
!==================================================================================================================================!       
    subroutine readModel
!
! Subroutine to read in an Occam1DCSEM model file and allocate model dependent
! arrays.
 


! Step 1: read the .resistivity file: 
    call readResistivityFile
 
! Step 2: read in the Model Boundary .poly file:
    call readPoly
     
! Step 3: read in the Model Penalty .penalty file:
    call readPenalty
    
! Step 4: read the MARE2DEM settings file:
    call readSettings
! 
!  Check boundary model for slivers
!
    call checkModel
           
    end subroutine readModel 
    
!==================================================================================================================================! 
!========================================================================================================================== readData
!==================================================================================================================================!  
    subroutine readData
!
! Reads in the EMData_2.0 format data file used in the MARE2DEM inversion
!
 
    
    implicit none

    include 'EM_parameters.inc'
          
!
! Local variables:
! 
    integer         :: i, err, iTxRead, iRxRead, iFreqRead, iDataRead, nAllocErr, iLen, iRx, iTx
    character(180)  :: sLine, sCode, sValue ! These are for reading the lines of the file
    character(180)  :: sFields(9),  sFmt
    logical         :: bComment, lerror 
    character(32)   :: cNum
    character(180)  :: sDataFormat 
    
    integer, parameter :: iof = 15   ! I/O File identifier  
    


    write(*,*) ' '    
    write(6,*) '=============== Reading in Data File ====================================='
    write(*,*) ' '    
    
!
! Open the data file:
!
    open (iof, file=trim(dataFileName), status='old', iostat=err)
    if (err /= 0) then
        write(*,*) ' Error opening data file:', trim(dataFileName)
        call exitMARE2DEM
    else
        call print2col32('Reading data from file: ', dataFileName,6)
    end if  
    
!
! Loop through the model file looking for the Format and Number of Layers flags
! and skipping over any comments
!
    do while (.true.) 
    
        ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        
        read( iof, '(A180)', iostat = err ) sLine
       
        if (err /= 0) exit  ! end of file read, escape from while(.true.) loop and 
                            ! proceed to checking that required inputs are defined
        call ParseCode( sLine, sCode, sValue, bComment )

        if( bComment ) cycle
       
!
! What do we have?
!
        select case (trim(sCode))
       
        case ('format')         ! Set cnDataParams for the input data format
        
            call lower(sValue) 
            sDataFormat = trim(sValue)
            
            select case (trim(sDataFormat) )
            
            case ('emdata_2.0','emdata_2.1','emdata_2.2')  
                cnDataParams = 4    
                
 
            case default
            
                write(*,*) ' Error: data format unsupported: ', trim(sDataFormat)
                write(*,*) ' Try using EMDATA_2.2'
                close(iof)
                call exitMARE2DEM
            end select
            
            call print2col32('Format: ',trim(sValue),6) 
     
            
        case ('utm of x,y origin (utm zone, n, e, 2d strike)','utm','origin')
            ! These aren't used in MARE2DEM, but save the line so we can include it in the .resp files:
            cUTMLine = trim(adjustl(sValue))
            call print2col32('UTM of origin (zone,N,E,theta): ',trim(sValue),6) 
            
        case ('phase','phase convention')  
            call lower(sValue)
            
            select case (trim(sValue))
            
            case ('lag')
                    phaseConvention = 'lag'
            case ('lead')
                    phaseConvention = 'lead'
            case default
                write(*,*) ' Error: unrecognized phase convenction: ', trim(sValue)
                write(*,*) ' Try using lag or lead'
                close(iof)
                call exitMARE2DEM          
            end select
            call print2col32('Phase convention: ', phaseConvention,6)

        case ('reciprocity','reciprocity used')  
            call lower(sValue)
            reciprocityUsed = trim(sValue) ! this is just passed through MARE2DEM to the output files
!            select case (trim(sValue))
!            case ('yes')
!                    reciprocityUsed = 'yes'
!            case ('no')
!                    reciprocityUsed = 'no'
!            case default
                !write(*,*) ' Error: unrecognized reciprocity used value: ', trim(sValue)
                !write(*,*) ' Try using yes or no'
                !close(iof)
                !call exitMARE2DEM          
            !end select
            call print2col32('Reciprocity used: ', reciprocityUsed,6)      
                        
                               
        case ('# transmitters')
        
            read(sValue,*) nTxCSEM
            
            write(*,*) ' '
            write(cnum,'(i3)') nTxCSEM
            call print2col32('# Transmitters: ',cnum,6)
            !write(*,*) ''    
            !write(*,*) '  Transmitter Positions: '

            sFmt = '(a1,1x,7(a12,2x),a12)' 
            write(*,sFmt) ' ','X','Y','Z','Azimuth','Dip','Length','Type','Name'

            
             
            !
            !  Allocate arrays dependent on nTx
            !
            allocate ( azimuthTxCSEM(nTxCSEM),dipTxCSEM(nTxCSEM), &
                     & xTxCSEM(nTxCSEM),yTxCSEM(nTxCSEM),zTxCSEM(nTxCSEM), lengthTxCSEM(nTxCSEM), &
                     & cSourceType(nTxCSEM), cTxNamesCSEM(nTxCSEM),   stat=nAllocErr )   
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many transmitters (', nTxCSEM, ')'
                call exitMARE2DEM 
            endif        
            cTxNamesCSEM = ' '
            lengthTxCSEM = 0
                
            !
            ! Now read in block of transmitters, skipping any comment lines:
            !
            
            iTxRead = 0
            do while (iTxRead < nTxCSEM) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iTxRead  = iTxRead + 1    
                
                select case (trim(sDataFormat) )
                
                case ('emdata_2.0')  
                    
                    read(sLine,*) xTxCSEM(iTxRead),yTxCSEM(iTxRead),zTxCSEM(iTxRead), &
                                & azimuthTxCSEM(iTxRead),dipTxCSEM(iTxRead),cSourceType(iTxRead)                     
 
                    
                case ('emdata_2.1') 
                
                    call parseFields( sLine, 7, sFields)
                     
                    read(sFields(1),*) xTxCSEM(iTxRead)
                    read(sFields(2),*) yTxCSEM(iTxRead) 
                    read(sFields(3),*) zTxCSEM(iTxRead) 
                    
                    read(sFields(4),*) azimuthTxCSEM(iTxRead)
                    read(sFields(5),*) dipTxCSEM(iTxRead)
                    read(sFields(6),*) cSourceType(iTxRead)

                    if (len_trim( sFields(7) ) > 0 ) then
                        read(sFields(7),*) cTxNamesCSEM(iTxRead)    
                    else
                        cTxNamesCSEM(iTxRead)     = ' '
                    endif
 
 
                case ('emdata_2.2') 
                
                    call parseFields( sLine, 8, sFields)
                     
                    read(sFields(1),*) xTxCSEM(iTxRead)
                    read(sFields(2),*) yTxCSEM(iTxRead) 
                    read(sFields(3),*) zTxCSEM(iTxRead) 
                    
                    read(sFields(4),*) azimuthTxCSEM(iTxRead)
                    read(sFields(5),*) dipTxCSEM(iTxRead)
                    read(sFields(6),*) lengthTxCSEM(iTxRead)
                    read(sFields(7),*) cSourceType(iTxRead)

                    if (len_trim( sFields(8) ) > 0 ) then
                        read(sFields(8),*) cTxNamesCSEM(iTxRead)    
                    else
                        cTxNamesCSEM(iTxRead)     = ' '
                    endif
                

                end select                
                              
                 
                sFmt = '(2x,6(f12.1,2x),a12,2x,a12)'            


                if ( iTxRead < 101 ) then
                    write(*,sFmt) xTxCSEM(iTxRead),yTxCSEM(iTxRead),zTxCSEM(iTxRead), &
                            & azimuthTxCSEM(iTxRead),dipTxCSEM(iTxRead), lengthTxCSEM(iTxRead), &
                            & trim(adjustl(cSourceType(iTxRead))), trim(adjustl(cTxNamesCSEM(iTxRead)))
                elseif( iTxRead == 101 ) then
                    write(*,*) 'More than 101 CSEM transmitters. Not writing the rest.'  
                endif
                                                         
            enddo
            !
            ! Check and make sure all the transmitters were read in:
            !
            if (iTxRead/=nTxCSEM) then
                write(*,*) ' Error reading transmitters: iTxRead/=nTxCSEM:',iTxRead, nTxCSEM
                close(iof)
                call exitMARE2DEM            
            
            endif       
            
        
                
        case ('# csem frequencies')
            
            read(sValue,*) nFreqCSEM
            write(cnum,'(i3)') nFreqCSEM
            write(*,*) ' '
            call print2col32('# CSEM Frequencies: ',cnum,6) 
            allocate (fTxCSEM(nFreqCSEM), stat=nAllocErr )   

            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many CSEM frequencies (', nFreqCSEM, ')'
                call exitMARE2DEM 
            endif   
            
            !
            ! Now read in block of frequencies, skipping any comment lines:
            !
            iFreqRead = 0
            do while (iFreqRead < nFreqCSEM) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
               
               ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iFreqRead  = iFreqRead + 1                  
                read(sLine,*)  fTxCSEM(iFreqRead)
                write(*,'(g12.3)') fTxCSEM(iFreqRead)
               
        
            enddo
            !
            ! Check and make sure all the frequencies were read in:
            !
            if (iFreqRead/=nFreqCSEM) then
                write(*,*) ' Error reading CSEM frequencies: iFreqRead/=nFreqCSEM:',iFreqRead, nFreqCSEM
                close(iof)
                call exitMARE2DEM                        
            endif          
             
        case ('# csem receivers')
         
            read(sValue,*) nRxCSEM
            write(cnum,'(i6)') nRxCSEM
            write(*,*) ' '
            call print2col32(' # CSEM Receivers: ',cnum,6) 
            
 
            write(*,'(a1,1x,7(a12,2x),a12)') ' ','X','Y', 'Z','Theta','Alpha','Beta','Length', 'Name' 
               
        
            !
            ! Allocate the field arrays: 
            !
            allocate(xRxCSEM(nRxCSEM), yRxCSEM(nRxCSEM), zRxCSEM(nRxCSEM), lengthRxCSEM(nRxCSEM),  &
                      ThetaRxCSEM(nRxCSEM),AlphaRxCSEM(nRxCSEM),BetaRxCSEM(nRxCSEM), cRxNamesCSEM(nRxCSEM),  &
                      stat=nAllocErr )   ! rxnames go here
                   
            cRxNamesCSEM = ' '
            lengthRxCSEM = 0
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many CSEM receivers (', nRxCSEM, ')'
                call exitMARE2DEM 
            endif         

            
 
         
            !
            ! Now read in block of receiver locations, skipping any comment lines:
            !
            iRxRead = 0
 
            do while (iRxRead < nRxCSEM) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1                        
     
                select case (trim(sDataFormat) )
            
                case ('emdata_2.0')  
                
     
                    call parseFields( sLine, 6, sFields)
                    read(sFields(1),*) xRxCSEM(iRxRead)
                    read(sFields(2),*) yRxCSEM(iRxRead)
                    read(sFields(3),*) zRxCSEM(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxCSEM(iRxRead)
                    read(sFields(5),*) AlphaRxCSEM(iRxRead)
                    read(sFields(6),*) BetaRxCSEM(iRxRead)

  
                case ('emdata_2.1')  
                
                    call parseFields( sLine, 7, sFields)
                     
                    read(sFields(1),*) xRxCSEM(iRxRead)
                    read(sFields(2),*) yRxCSEM(iRxRead)
                    read(sFields(3),*) zRxCSEM(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxCSEM(iRxRead)
                    read(sFields(5),*) AlphaRxCSEM(iRxRead)
                    read(sFields(6),*) BetaRxCSEM(iRxRead)

                    if (len_trim( sFields(7) ) > 0 ) then
                        read(sFields(7),*) cRxNamesCSEM(iRxRead)
                    else
                        cRxNamesCSEM(iRxRead) = ' '
                    endif
  
                  case ('emdata_2.2')  
                
                    call parseFields( sLine, 8, sFields)
                     
                    read(sFields(1),*) xRxCSEM(iRxRead)
                    read(sFields(2),*) yRxCSEM(iRxRead)
                    read(sFields(3),*) zRxCSEM(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxCSEM(iRxRead)
                    read(sFields(5),*) AlphaRxCSEM(iRxRead)
                    read(sFields(6),*) BetaRxCSEM(iRxRead)
                    read(sFields(7),*) lengthRxCSEM(iRxRead) 
                    
                    if (len_trim( sFields(8) ) > 0 ) then
                        read(sFields(8),*) cRxNamesCSEM(iRxRead)
                    else
                        cRxNamesCSEM(iRxRead) = ' '
                    endif
                      
                end select
  
                if ( iRxRead < 101 ) then
                    write(*,'(2x,7(f12.1,2x),a12)')  xRxCSEM(iRxRead),yRxCSEM(iRxRead),zRxCSEM(iRxRead), & 
                         & ThetaRxCSEM(iRxRead),AlphaRxCSEM(iRxRead),BetaRxCSEM(iRxRead),lengthRxCSEM(iRxRead), & 
                         & trim(adjustl(cRxNamesCSEM(iRxRead)))
                elseif( iRxRead == 101 ) then
                    write(*,*) 'More than 101 CSEM receivers. Not writing the rest.'  
                endif
                                          
            enddo
            !
            ! Check and make sure all the rx were read in:
            !
            if (iRxRead/=nRxCSEM) then
                write(*,*) ' Error reading receivers: iRxRead/=nRxCSEM:',iRxRead,nRxCSEM
                close(iof)
                call exitMARE2DEM                        
            endif       
            
        case ('csem receiver names') 
           
       
            cRxNamesCSEM = ''   
            iRxRead = 0
            write(*,*) ' '
            call print2col32('CSEM Receiver Names: ', ' ',6)
            do while (iRxRead < nRxCSEM) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1                        
     
                 cRxNamesCSEM(iRxRead) = trim(adjustl(sLine))
                 write(*,'(1x,a)') trim(cRxNamesCSEM(iRxRead) )
  
                
            enddo  
                  
        case ('# mt frequencies')
            
            read(sValue,*) nFreqMT
            write(cnum,'(i3)') nFreqMT
            write(*,*) ' '
            call print2col32('# MT Frequencies: ',cnum,6) 
            allocate (fTxMT(nFreqMT), stat=nAllocErr )   

            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many MT frequencies (', fTxMT, ')'
                call exitMARE2DEM 
            endif   
            
            !
            ! Now read in block of frequencies, skipping any comment lines:
            !
            iFreqRead = 0
            do while (iFreqRead < nFreqMT) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
               
               ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iFreqRead  = iFreqRead + 1                  
                read(sLine,*)  fTxMT(iFreqRead)
                write(*,'(g12.3)') fTxMT(iFreqRead)
               
        
            enddo
            !
            ! Check and make sure all the frequencies were read in:
            !
            if (iFreqRead/=nFreqMT) then
                write(*,*) ' Error reading MT frequencies: iFreqRead/=nFreqMT:',iFreqRead, nFreqMT
                close(iof)
                call exitMARE2DEM                        
            endif          
             
      
  
        case ('# mt receivers')
         
            read(sValue,*) nRxMT
            write(cnum,'(i6)') nRxMT
            write(*,*) ' '
            call print2col32('# MT Receivers: ',cnum,6) 
            
 
            write(*,'(a1,1x,9(a12,2x))') ' ','x','y', 'z','Theta','Alpha','Beta', 'Length', 'SolveStatic', 'Name'

                
            !
            ! Allocate the field arrays: 
            !
            allocate( xRxMT(nRxMT), yRxMT(nRxMT), zRxMT(nRxMT),ThetaRxMT(nRxMT), AlphaRxMT(nRxMT),BetaRxMT(nRxMT), &
                    &  cRxNamesMT(nRxMT), iSolveStatic(nRxMT),  lengthRxMT(nRxMT), stat=nAllocErr )
            iSolveStatic = 0 
            
            cRxNamesMT = ' ' 
            lengthRxMT = 0
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many MT receivers (', nRxMT, ')'
                call exitMARE2DEM 
            endif         
 
!         
            !
            ! Now read in block of receiver locations, skipping any comment lines:
            !
            iRxRead = 0
 
            do while (iRxRead < nRxMT) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1        
         
                 
                select case (trim(sDataFormat) )
            
                case ('emdata_2.0')  
                
                    call parseFields( sLine, 6, sFields)
                    read(sFields(1),*) xRxMT(iRxRead)
                    read(sFields(2),*) yRxMT(iRxRead)
                    read(sFields(3),*) zRxMT(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxMT(iRxRead)
                    read(sFields(5),*) AlphaRxMT(iRxRead)
                    read(sFields(6),*) BetaRxMT(iRxRead)         
                               
       
           
                case ('emdata_2.1')  
                
                    call parseFields( sLine, 8, sFields)
                     
                    read(sFields(1),*) xRxMT(iRxRead)
                    read(sFields(2),*) yRxMT(iRxRead)
                    read(sFields(3),*) zRxMT(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxMT(iRxRead)
                    read(sFields(5),*) AlphaRxMT(iRxRead)
                    read(sFields(6),*) BetaRxMT(iRxRead)
                     
                    read(sFields(7),*) iSolveStatic(iRxRead)
                     
                    if (len_trim( sFields(8) ) > 0 ) then
                        read(sFields(8),*) cRxNamesMT(iRxRead)
                    else
                        cRxNamesMT(iRxRead) = ' '
                    endif

                case ('emdata_2.2')  
                
                    call parseFields( sLine, 9, sFields)
                     
                    read(sFields(1),*) xRxMT(iRxRead)
                    read(sFields(2),*) yRxMT(iRxRead)
                    read(sFields(3),*) zRxMT(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxMT(iRxRead)
                    read(sFields(5),*) AlphaRxMT(iRxRead)
                    read(sFields(6),*) BetaRxMT(iRxRead)
                    read(sFields(7),*) lengthRxMT(iRxRead)
                    read(sFields(8),*) iSolveStatic(iRxRead)
                     
                    if (len_trim( sFields(9) ) > 0 ) then
                        read(sFields(9),*) cRxNamesMT(iRxRead)
                    else
                        cRxNamesMT(iRxRead) = ' '
                    endif
               
                end select
            
                if ( iRxRead < 101 ) then
                        write(*,'(2x,7(f12.1,2x),i12,2x,a12)')  xRxMT(iRxRead),yRxMT(iRxRead),zRxMT(iRxRead), & 
                             & ThetaRxMT(iRxRead),AlphaRxMT(iRxRead),BetaRxMT(iRxRead), lengthRxMT(iRxRead),iSolveStatic(iRxRead), &
                                trim(adjustl(cRxNamesMT(iRxRead)))  
                elseif( iRxRead == 101 ) then
                    write(*,*) 'More than 101 MT receivers. Not writing the rest.'  
                endif
                                          
              
            enddo
            !
            ! Check and make sure all the rx were read in:
            !
            if (iRxRead/=nRxMT) then
                write(*,*) ' Error reading MT receivers: iRxRead/=nRxMT:',iRxRead,nRxMT
                close(iof)
                call exitMARE2DEM                        
            endif                    
            
        case ('mt receiver names')             
               
            cRxNamesMT = ' '   
            iRxRead = 0
            write(*,*) ' '
            call print2col32('MT Receiver Names: ', ' ',6)
            do while (iRxRead < nRxMT) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1                        
     
                 cRxNamesMT(iRxRead) = trim(adjustl(sLine))
                 write(*,'(1x,a)') trim(cRxNamesMT(iRxRead)) 
 
         enddo  
   
    !
    ! Lastly we need to read in the DATA block.
    !
         case ('# data','#data','ndata')
         
            read(sValue,*) nd
            write(*,*) ' ' 
            write(cnum,'(i6)') nd
            call print2col32('Number of Data: ',cnum,6)
                      
            !
            ! Allocate the nd data arrays: 
            !
            allocate(dp(nd,cnDataParams), d(nd), d_wt(nd), sd(nd), dm(nd), stat=nAllocErr ) 
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many data points (', nd, ')'
                call exitMARE2DEM 
            endif        
            d_wt = 1. ! standard data weight, modified later on if joint MT + CSEM inversion
             
            iDataRead = 0
            write(*,'(4(a12,1x),a12,2x,a12,2x)') 'Type','Time/Freq #','Tx #','Rx #', 'Data','Std Error'    
            
            do while (iDataRead < nd) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
               
               ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iDataRead  = iDataRead + 1     
                
                ! Read in 6 fields as character strings,
                ! this is clunky but works for mixed character/numeric tables:
                call parseFields( sLine, 6, sFields)    

                ! Test first field to see if it is character data type, 
                ! if so convert it to numeric code
                iLen = len_trim(sFields(1))
                
                if (iLen < 5) then  ! Field is integer code for data type
                    read(sFields(1),*) dp(iDataRead,1)  
                    
                else ! Field is ascii character string describing data:
                    call getDataCode(sFields(1),dp(iDataRead,1))
                
                endif

                ! Parse the remaining fields for the data parameters and data 
                read(sFields(2),*) dp(iDataRead,2)
                read(sFields(3),*) dp(iDataRead,3)
                read(sFields(4),*) dp(iDataRead,4)
                read(sFields(5),*) d(iDataRead)
                read(sFields(6),*) sd(iDataRead)

                
                ! DGM 10/2010 Don't spew everything to the screen - only need to
                ! check the first so many.
                if( iDataRead < 101 ) then
                    write(*,fmt= '(4(i12,1x),ES12.3,2x,ES12.3,2x)') (dp(iDataRead,i),i=1,4),d(iDataRead),sd(iDataRead)
                elseif( iDataRead == 101 ) then
                    write(*,*) 'More than 101 data. Not writing the rest.'  
                endif
            enddo    
            write(*,*) ' '    
            !
            ! Check and make sure all the data was read in:
            !
            if (iDataRead/=nd) then
                write(*,*) ' Error reading data: iDataRead/=nd:',iDataRead,nd
                close(iof)
                call exitMARE2DEM                        
            endif      
            !

            
            
        case default
                write(*,*) 'Error, unknown code in Data file: ',trim(sCode)
                write(*,*) 'skipping this, hopefully its not crucial!'
                !call exitMARE2DEM
            
        end select ! case (trim(sCode))
        
        
    enddo ! while (.true.) ! The main data reading loop
!
! Close the data file
!
    close(iof)


!
! Check that parameters have been defined:
!
    lerror = .false.
    if ((nTxCSEM < 1).and.(nRxMT<1)) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of transmitters: nTxCSEM = ',nTxCSEM
    endif
    if ( (nRxCSEM < 1) .and. (nRxMT<1) ) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of receivers: nRxCSEM, nRxMT = ',nRxCSEM, nRxMT
    endif
    if ((nFreqCSEM < 1) .and. (nFreqMT < 1) ) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of frequencies: nFreqCSEM,nFreqMT = ',nFreqCSEM,nFreqMT
    endif
    if ( (nd < 1).and.(.not.lFwdFields) ) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of data: # Data = ',nd
    endif   
    
    if (lerror) then
        write(*,*)  ' Stopping ! '
        call exitMARE2DEM
    endif
 
!
! Check to see which data components are present and set logical flags to false for the unused field components in order to 
! save memory & time during the inversion calls to mare2dem
!
    
    allocate( lCompCSEM(6, nRxCSEM ), lCompMT(6,nRxMT) )
    lCompCSEM   = .false.
    lCompMT     = .false.
     
    nMT   = 0
    nCSEM = 0 
    
    do i = 1,nd
        
        iRx = dp(i,4)
        iTx = dp(i,3) 
        if (iTx <= 0) iTx = iRx  ! if iTx not specified, make it the current receiver 
            
        select case (  dp(i,1) ) ! what data type is this?
            

                                            
            ! CSEM types:                                  
            case (indRealEx,indImagEx,indAmpEx,indPhsEx,indLog10AmpEx)  
                lCompCSEM(1,iRx)    = .true. 
                nCSEM               = nCSEM + 1 
            case (indRealEy,indImagEy,indAmpEy,indPhsEy,indLog10AmpEy)  
                lCompCSEM(2,iRx)    = .true. 
                nCSEM               = nCSEM + 1 
            case (indRealEz,indImagEz,indAmpEz,indPhsEz,indLog10AmpEz)  
                lCompCSEM(3,iRx)    = .true. 
                nCSEM               = nCSEM + 1 
            case (indRealBx,indImagBx,indAmpBx,indPhsBx,indLog10AmpBx)  
                lCompCSEM(4,iRx)    = .true. 
                nCSEM               = nCSEM + 1 
            case (indRealBy,indImagBy,indAmpBy,indPhsBy,indLog10AmpBy)  
                lCompCSEM(5,iRx)    = .true. 
                nCSEM               = nCSEM + 1 
            case (indRealBz,indImagBz,indAmpBz,indPhsBz,indLog10AmpBz)                      
                lCompCSEM(6,iRx)    = .true. 
                nCSEM               = nCSEM + 1    
            case (iPEmax,iPEmin)
                lCompCSEM(1,iRx)    = .true.
                lCompCSEM(2,iRx)    = .true. 
                nCSEM               = nCSEM + 1 
            case (iPBmax,iPBmin)
                lCompCSEM(4,iRx)    = .true.
                lCompCSEM(5,iRx)    = .true. 
                nCSEM               = nCSEM + 1 
                
             ! MT types:   
            case (indRhoZXY,indPhsZXY,indRealZXY,indImagZXY,indlog10RhoZXY, indRealMZY, indImagMZY )    
                lCompMT(1,iRx)  = .true.
                lCompMT(5,iRx)  = .true.
                lCompMT(6,iRx)  = .true.
                lCompMT(1,iTx)  = .true. ! for hybrid MT stations
                lCompMT(5,iTx)  = .true. ! for hybrid MT stations
                lCompMT(6,iTx)  = .true. ! for hybrid MT stations                
                nMT             = nMT + 1

            case (indRhoZYX,indPhsZYX,indRealZYX,indImagZYX,indlog10RhoZYX)    
                lCompMT(2,iRx)  = .true.
                lCompMT(3,iRx)  = .true.
                lCompMT(4,iRx)  = .true.
                lCompMT(2,iTx)  = .true. ! for hybrid MT stations
                lCompMT(3,iTx)  = .true. ! for hybrid MT stations
                lCompMT(4,iTx)  = .true. ! for hybrid MT stations               
                nMT             = nMT + 1
            
            case (indRhoZXX,indPhsZXX,indRealZXX,indImagZXX,indlog10RhoZXX,indRhoZYY,indPhsZYY,indRealZYY,indImagZYY,indlog10RhoZYY)
                write(*,*) ' Error in input data file, diagonal impedance elements not support in 2D. Stopping!'
                call exitMARE2DEM
                 
        end select  ! case dp(i,1)
    enddo  ! loop over nd

!
! If lFwdFields then user wants all field components computed regardless of input data array:
!
    if (lFwdFields) then
        if ( nFreqCSEM > 0 ) then
            lCompCSEM   = .true.
            nCSEM       = nFreqCSEM*nRxCSEM*nTxCSEM*6
        endif
      
        if ( nFreqMT > 0 ) then
            nMT         = nFreqMT*nRxMT*6
            lCompMT     = .true.
        endif
                
    
    endif
 
    do i = 1,nRxMT
        if (ThetaRxMT(i) /= 0 ) then
            write(*,*) 'Ignoring theta rotation angle for MT Receiver #',i
        endif
        if (AlphaRxMT(i) /= 0 ) then
            write(*,*) 'Ignoring alpha rotation angle for MT Receiver #',i
        endif           
        if (BetaRxMT(i) /= 0 )  then
             !this is allowed since its in the y-z plane
        endif
    enddo    
    
     write(*,*) 'Done reading the data file, here is what I found:'
     write(*,*)  ' '
     
 
    write(cnum,'(i6)') nCSEM
    call print2col32(' # CSEM Data: ',cnum,6)
    write(cnum,'(i6)') nMT
    call print2col32(' # MT Data: ',cnum,6)

    
 

!
! Allocate associated parameters that are carried through Occam for the best fitting model. Currently these are only used
! for the MT static shift parameter estimates:
!
    if (nRxMT > 0 ) then
        npm_assoc = nRxMT*2
        allocate( pm_assoc(npm_assoc) )
    endif 
    
!
! Finally, run a consistency check on the data parameter indices
!
! 'Type','Freq #','Tx #','Rx #', 'Data','Std Error' 
    do i = 1,nd
        select case (  dp(i,1) ) ! what data type is this?
                                          
        case (1:99)
            if (dp(i,2) > nFreqCSEM) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Frequency index exceeds # of CSEM frequencies'
                call exitMARE2DEM
            elseif (dp(i,3) > nTxCSEM) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Transmitter index exceeds # of CSEM transmitters'
                call exitMARE2DEM
            elseif (dp(i,4) > nRxCSEM) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Receiver index exceeds # of CSEM receivers'
                call exitMARE2DEM            
            endif
            
        case (101:199 ) 
                 
            if (dp(i,2) > nFreqMT) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Frequency index exceeds # of MT frequencies'
                call exitMARE2DEM
           
            elseif (dp(i,4) > nRxMT) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Receiver index exceeds # of MT receivers'
                call exitMARE2DEM            
            endif
        end select
    enddo
    write(*,*) ' '
    
    end subroutine readData       
  

!==================================================================================================================================! 
!==================================================================================================================== writeFields_MT
!==================================================================================================================================! 
    subroutine writeFields_MT( nCurrentIter)
!
! Writes out all 6 MT field components to a file 
! 

    use Occam 
    use mare2dem_input_data_params
    use mare2dem_global 
    use mare2dem_output 
    
    implicit none
    
    integer, intent(in)     :: nCurrentIter
    character(50)           :: cNum
    
    character(128)  :: sfmt
    integer         :: lerr, i, iFreq, iRx
 
    
    if (nFreqMT < 1) return
    
!    
! Open response file
!
    write (cNum,*) nCurrentIter
    open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.fieldsMT', iostat=lerr)

!
! Catch error opening response file:
!
    if (lerr /= 0) then
        write(*,*) ' Error opening MT fields file, stopping!'
        call exitMARE2DEM
    end if
    
!
! No error, write out the responses:
!    
 
    write(21,*) 'Format:    MTfields_1.0'
    
    ! UTM section:  
    if (len_trim(cUTMLine) > 0 ) then
        write(21,'(a,a)') 'UTM of x,y origin (UTM zone, N, E, 2D strike): ',trim(cUTMLine)
    endif       
    
    
    if (nRxMT>0) then    
        write(cnum,'(i6)') nFreqMT
        write(21,'(a,a)') '# MT Frequencies: ', adjustl(cnum)    
        do i = 1,nFreqMT
            write(21,'(g13.5)') fTxMT(i)
        enddo        
        write(cnum,'(i6)') nRxMT
        write(21,'(a,a)') '# MT Receivers: ', adjustl(cnum)   
        write(21,'(a1,8(a12,1x))') '!','X','Y','Z','Theta','Alpha','Beta','SolveStatic','Name'
        do i = 1,nRxMT
            write(21,'(1x,6(f12.1,1x),i12,1x,a12)') xRxMT(i),yRxMT(i),zRxMT(i), ThetaRxMT(i),AlphaRxMT(i),BetaRxMT(i), &
                                       & iSolveStatic(i), trim(adjustl(cRxNamesMT(i)))
        enddo   
 
        
    endif    
    
    sfmt = '(12(E22.14E3,1x))'
    do iFreq = 1,nFreqMT
        write(cnum,'(i6)') iFreq
        write(21,'(a,a)') 'Frequency # ', adjustl(cnum)    
        do iRx = 1,nRxMT
               
            write(21,sfmt)  real(ex_mt(iRx,iFreq)),aimag(ex_mt(iRx,iFreq)), &
                            real(ey_mt(iRx,iFreq)),aimag(ey_mt(iRx,iFreq)), &
                            real(ez_mt(iRx,iFreq)),aimag(ez_mt(iRx,iFreq)), &
                            real(hx_mt(iRx,iFreq)),aimag(hx_mt(iRx,iFreq)), &
                            real(hy_mt(iRx,iFreq)),aimag(hy_mt(iRx,iFreq)), &
                            real(hz_mt(iRx,iFreq)),aimag(hz_mt(iRx,iFreq)) 
 
        enddo
        
    enddo
    
    close(21)
 
    end subroutine writeFields_MT
    
!==================================================================================================================================! 
!================================================================================================================== writeFields_CSEM
!==================================================================================================================================! 
    subroutine writeFields_CSEM( nCurrentIter)
!
! Writes out all 6 CSEM field components to a file 
! 

    use Occam 
    use mare2dem_input_data_params
    use mare2dem_global 
    use mare2dem_output
 
    
    implicit none
    
    integer, intent(in)     :: nCurrentIter
    character(50)           :: cNum
    
    character(128)  :: sfmt
    integer         :: lerr, i, iFreq, iRx, iTx
 
    
    if (nFreqCSEM < 1) return  
!    
! Open response file
!
    write (cNum,*) nCurrentIter
    open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.fieldsCSEM', iostat=lerr)

!
! Catch error opening response file:
!
    if (lerr /= 0) then
        write(*,*) ' Error opening CSEM fields file, stopping!'
        call exitMARE2DEM
    end if
    
!
! No error, write out the responses:
!    
 
    write(21,*) 'Format:    CSEMfields_1.0'
    
    
    ! UTM section:  
    if (len_trim(cUTMLine) > 0 ) then
        write(21,'(a,a)') 'UTM of x,y origin (UTM zone, N, E, 2D strike): ',trim(cUTMLine)
    endif       

    write(21,'(a,a)') 'Phase convention: ',trim(adjustl(phaseConvention))
    write(cnum,'(i6)') nTxCSEM
    write(21,'(a,a)') '# Transmitters: ', adjustl(cnum)
    write(21,'(a1,7(a12,1x))') '!','X','Y','Z','Azimuth','Dip','Type', 'Name'   
    do i = 1,nTxCSEM
       write(21,'(1x,5(f12.1,1x),a12,1x,a12)') xTxCSEM(i),yTxCSEM(i),zTxCSEM(i),azimuthTxCSEM(i),dipTxCSEM(i),cSourceType(i), &
                                            & trim(adjustl(cTxNamesCSEM(i)))
    enddo          
    write(cnum,'(i6)') nFreqCSEM
    write(21,'(a,a)') '# CSEM Frequencies: ', adjustl(cnum)    
    do i = 1,nFreqCSEM
        write(21,'(g13.5)') fTxCSEM(i)
    enddo        
    write(cnum,'(i6)') nRxCSEM
    write(21,'(a,a)') '# CSEM Receivers: ', adjustl(cnum)   
    write(21,'(a1,7(a12,1x))') '!','X','Y','Z','Theta','Alpha','Beta', 'Name'  
    do i = 1,nRxCSEM
        write(21,'(1x,6(f12.1,1x),a12)') xRxCSEM(i),yRxCSEM(i),zRxCSEM(i), ThetaRxCSEM(i),AlphaRxCSEM(i),BetaRxCSEM(i), &
                                 &  trim(adjustl(cRxNamesCSEM(i)))
    enddo   

    sfmt = '(12(E22.14E3,1x))'
    
    do iFreq = 1,nFreqCSEM
        write(cnum,'(i6)') iFreq
        write(21,'(a,a)') 'Frequency # ', adjustl(cnum)    
        
        do iTx = 1,nTxCSEM
            write(cnum,'(i6)') iTx
            write(21,'(a,a)') 'Transmitter # ', adjustl(cnum)   
            
            do iRx = 1,nRxCSEM
               
                write(21,sfmt)  real(ex(iRx,iTx,iFreq)),aimag(ex(iRx,iTx,iFreq)), &
                                real(ey(iRx,iTx,iFreq)),aimag(ey(iRx,iTx,iFreq)), &
                                real(ez(iRx,iTx,iFreq)),aimag(ez(iRx,iTx,iFreq)), &
                                real(hx(iRx,iTx,iFreq)),aimag(hx(iRx,iTx,iFreq)), &
                                real(hy(iRx,iTx,iFreq)),aimag(hy(iRx,iTx,iFreq)), &
                                real(hz(iRx,iTx,iFreq)),aimag(hz(iRx,iTx,iFreq)) 
 
            enddo
        enddo        
    enddo
    
    close(21)
 
    end subroutine writeFields_CSEM
  
!==================================================================================================================================! 
!===================================================================================================================== writeJacobian
!==================================================================================================================================! 
    subroutine writeJacobian( nCurrentIter)
!     
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! 2011-2012     Updated for MARE2DEM inversion code
!
    use Occam 
    use mare2dem_input_data_params
    use mare2dem_global 
 
    
    implicit none
    
    integer         :: nCurrentIter
    character(50)   :: cNum
    integer         :: lerr, i, j
 

 
!----------------------------
!    
! Open Jacobian file
!
    write (cNum,*) nCurrentIter
    open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.jacobian', iostat=lerr)

!
! Catch error opening file:
!
    if (lerr /= 0) then
        write(*,*) ' Error opening Jacobian file, stopping!'
        call exitMARE2DEM
    end if
!
! A-okay, write out the Jacobian
! 
!
! Note that this routine was called from RunMARE2DEM and so WJ is just J (i.e. it has not been normalized by the data error bars):
!
    write(21,*)  nd, nParams
 
    do j = 1,nParams    
        do i = 1,nd
            write(21,'(g12.4)') wj(i,j)  
        enddo
    enddo
     
    close(21)
      
    end subroutine writeJacobian

!==================================================================================================================================! 
!======================================================================================================================= getDataCode
!==================================================================================================================================! 
    subroutine getDataCode( sField, iCode)
!
! Management routine to get data format magic number code from ascii text label.
!
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
!
! Revision      Date            Author          Notes
! 1.0           November 2008   Kerry Key       Created
!
    
    implicit none
 
    include 'EM_parameters.inc'
        
    character(180), intent(inout) :: sField
    integer, intent(out) :: iCode
    
    !
    ! Get lower case string:
    !
    call lower(sField) 
    
    !
    ! Decode string:
    !
    select case (trim(sField))
    
    !
    ! CSEM Data:
    !
        ! Electric Real/Imag Data:
        case ('realex','reex')
            iCode = indRealEx
        case ('imagex','imex','imaginaryex')    
            iCode = indImagEx
        case ('realey','reey')
            iCode = indRealEy
        case ('imagey','imey','imaginaryey')    
            iCode = indImagEy   
        case ('realez','reez')
            iCode = indRealEz   
        case ('imagez','imez','imaginaryez')    
            iCode = indImagEz
        
        ! Magnetic Real/Imag Data:
        case ('realbx','rebx')
            iCode = indRealBx
        case ('imagbx','imbx','imaginarybx')    
            iCode = indImagBx
        case ('realby','reby')
            iCode = indRealBy
        case ('imagby','imby','imaginaryby')    
            iCode = indImagBy       
        case ('realbz','rebz')
            iCode = indRealBz
        case ('imagbz','imbz','imaginarybz')    
            iCode = indImagBz   
        
        ! Electric Amplitude and Phase Data:
        case ('ampex','amplitudeex')
            iCode = indAmpEx
        case ('phsex','phaseex')
            iCode = indPhsEx
        case ('ampey','amplitudeey')
            iCode = indAmpEy
        case ('phsey','phaseey')
            iCode = indPhsEy
        case ('ampez','amplitudeez')
            iCode = indAmpEz                
        case ('phsez','phaseez')
            iCode = indPhsEz    
        
        ! Log10 Electric Amplitude Data:
        case ('log10ampex','log10amplitudeex')        
            iCode = indAmpEx           
        case ('log10ampey','log10amplitudeey')        
            iCode = indAmpEy           
        case ('log10ampez','log10amplitudeez')        
            iCode = indAmpEz        

        ! Log10 Magnetic Amplitude Data:
        case ('log10ampbx','log10amplitudebx')        
            iCode = indAmpBx           
        case ('log10ampby','log10amplitudeby')        
            iCode = indAmpBy           
        case ('log10ampbz','log10amplitudebz')        
            iCode = indAmpBz        


        ! Magnetic Amplitude and Phase Data:
        case ('ampbx','amplitudebx')
            iCode = indAmpBx        
        case ('phsbx','phasebx')
            iCode = indPhsBx
        case ('ampby','amplitudeby')
            iCode = indAmpBy
        case ('phsby','phaseby')
            iCode = indPhsBy
        case ('ampbz','amplitudebz')
            iCode = indAmpBz        
        case ('phsbz','phasebz')
            iCode = indPhsBz    
            
        ! CSEM polarization ellipse parameters:
        case('pemax')
            icode = iPEmax
        case('pemin')
            icode = iPEmin
        case('pbmax')
            icode = iPBmax
        case('pbmin')
            icode = iPBmin
    !   
    ! Magnetotelluric Data:
    !
        ! Real/Imag
        case ('realzxx')        ! ignored in 2D
            iCode = indRealZXX  
        case ('imagzxx')        ! ignored in 2D
            iCode = indImagZXX      
        case ('realzxy')
            iCode = indRealZXY  
        case ('imagzxy')
            iCode = indImagZXY      
        case ('realzyx')
            iCode = indRealZYX  
        case ('imagzyx')
            iCode = indImagZYX          
        case ('realzyy')        ! ignored in 2D
            iCode = indRealZYY  
        case ('imagzyy')        ! ignored in 2D
            iCode = indImagZYY                  
        
        ! Apparent resistivity and phase:
        case ('rhozxx')         ! ignored in 2D
            iCode = indRhoZXX   
        case ('phszxx')         ! ignored in 1D
            iCode = indPhsZXX       
        case ('rhozxy')         
            iCode = indRhoZXY   
        case ('phszxy')         
            iCode = indPhsZXY       
        case ('rhozyx')         
            iCode = indRhoZYX   
        case ('phszyx')         
            iCode = indPhsZYX       
        case ('rhozyy')         ! ignored in 2D
            iCode = indRhoZYY   
        case ('phszyy')         ! ignored in 2D
            iCode = indPhsZYY       
                
        case ('log10rhozxx')          
            iCode = indlog10RhoZXX   
        case ('log10rhozxy')         
            iCode = indlog10RhoZXY   
        case ('log10rhozyx')         
            iCode = indlog10RhoZYX   
        case ('log10rhozyy')          
            iCode = indlog10RhoZYY   
        
        ! Magnetic Tipper:
        case ('realmzy')
            iCode = indRealMZY 
        case ('imagmzy')
            iCode = indImagMZY      
          
        ! Error if unknown type:
        case default
            write(*,*) ' Error decoding data type: ', sField
            write(*,*) ' Stopping!'
            call exitMARE2DEM
        
    
    end select

    
    end subroutine getDataCode       
          
!==================================================================================================================================! 
!======================================================================================================================= print2col32
!==================================================================================================================================! 
    subroutine print2col32(str1,str2,ioUnit) 
    
    implicit none
    
    character(len=*),intent(in)  :: str1, str2
    integer, intent(in)          :: ioUnit
    
 
    character(32) :: str1b
    
    str1b = str1 
    write(ioUnit,'(1x,a,a)') (adjustl(str1b)), trim(adjustl(str2))
     
    end subroutine print2col32
    
!-----------------------------------------------------------------------------------------------------------------------------------

    end module mare2dem_io
    

! This routine has to exist outside the module so that Occam can find it:    
!==================================================================================================================================! 
!===================================================================================================================== writeResponse
!==================================================================================================================================! 
    subroutine writeResponse( nCurrentIter)
!     
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! 2011-2012     Updated for MARE2DEM inversion code
!
    use Occam 
    use mare2dem_input_data_params
    use mare2dem_global
 
    
    implicit none
    
    integer         :: nCurrentIter
    character(50)   :: cNum
    character(128)  :: cStr
    integer         :: lerr, i, j, ncolumns
    
    real(RealPrec) :: residMT, residCSEM, resid
    

!----------------------------------------------
! Display MT and CSEM misfits:
!
    residMT  = 0
    residCSEM = 0
    do i = 1,nd
        if (sd(i) == 0 ) cycle ! skip if data uncertainty is 0
        
        if ( dp(i,1) < 100 ) then ! CSEM data:
               
             residCSEM = residCSEM + ( (d(i) - dm(i)) /sd(i))**2
             
        else ! MT data: dp(i,1) > 100 
 
              residMT = residMT + ( (d(i) - dm(i)) /sd(i))**2
              
        endif

    enddo
    
    ! If joint CSEM/MT inversion, display misfit for each type separately:
    if ((nCSEM> 0).and.(nMT> 0)) then
        write(cStr,'(a32,g16.4)') 'CSEM Misfit: ', sqrt(residCSEM/nCSEM)
        call printOccamLog(cStr)    
        write(cStr,'(a32,g16.4)') 'MT Misfit: ', sqrt(residMT/nMT)
        call printOccamLog(cStr)    
         
    endif
     
     ! Display any MT static shift solutions:
     
     if ( nRxMT > 0 ) then
     
         if (any(iSolveStatic > 0)) then
             
             write(cStr,'(a)') ' MT static shift estimates: '
             call printOccamLog(cStr)   
             write(cStr,'(4(a12,2x),a34)') ' Site #' ,'Name', 'TE Factor' , 'TM Factor', ' where linear ApRes = ApRes*Factor'
             call printOccamLog(cStr)   
             do i = 1,nRxMT
                 if (iSolveStatic(i) > 0) then
                    write(cStr,'(i12,2x,a12,2x,f12.3,2x,f12.3)') i, trim(cRxNamesMT(i)), 10**pm_assoc(2*i-1), 10**pm_assoc(2*i)
                    call printOccamLog(cStr) 
                 endif
             
             enddo
         endif
     
     endif
         
        
!----------------------------
!    
! Open response file
!
    write (cNum,*) nCurrentIter
    open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.resp', iostat=lerr)

!
! Catch error opening response file:
!
    if (lerr /= 0) then
        write(*,*) ' Error opening response file, stopping!'
        call exitMARE2DEM
    end if
!
! No error, write out the responses:
!    
 
    write(21,'(a)') 'Format:    EMResp_2.2'
     
    ! UTM section:  
    if (len_trim(cUTMLine) > 0 ) then
        write(21,'(a,a)') 'UTM of x,y origin (UTM zone, N, E, 2D strike): ',trim(cUTMLine)
    endif   
    
    ! CSEM section:
    if (nTxCSEM>0) then
        write(21,'(a,a)') 'Phase convention: ',trim(adjustl(phaseConvention))
        write(21,'(a,a)') 'Reciprocity used: ',trim(adjustl(reciprocityUsed))
        write(cnum,'(i6)') nTxCSEM
        write(21,'(a,a)') '# Transmitters: ', adjustl(cnum)
        write(21,'(a1,8(a12,1x))') '!','X','Y','Z','Azimuth','Dip','Length','Type', 'Name'   
        do i = 1,nTxCSEM
           write(21,'(1x,6(f12.1,1x),a12,1x,a12)') xTxCSEM(i),yTxCSEM(i),zTxCSEM(i),azimuthTxCSEM(i),dipTxCSEM(i),lengthTxCSEM(i),&
                                                & cSourceType(i), trim(adjustl(cTxNamesCSEM(i)))
        enddo          
        write(cnum,'(i6)') nFreqCSEM
        write(21,'(a,a)') '# CSEM Frequencies: ', adjustl(cnum)    
        do i = 1,nFreqCSEM
            write(21,'(g13.5)') fTxCSEM(i)
        enddo        
        write(cnum,'(i6)') nRxCSEM
        write(21,'(a,a)') '# CSEM Receivers: ', adjustl(cnum)   
        write(21,'(a1,8(a12,1x))') '!','X','Y','Z','Theta','Alpha','Beta','Length', 'Name'  
        do i = 1,nRxCSEM
            write(21,'(1x,7(f12.1,1x),a12)') xRxCSEM(i),yRxCSEM(i),zRxCSEM(i), ThetaRxCSEM(i),AlphaRxCSEM(i),BetaRxCSEM(i), &
                                           & lengthRxCSEM(i), trim(adjustl(cRxNamesCSEM(i)))  
        enddo   
 
    endif
    
    ! MT section:
    if (nRxMT>0) then    
        write(cnum,'(i6)') nFreqMT
        write(21,'(a,a)') '# MT Frequencies: ', adjustl(cnum)    
        do i = 1,nFreqMT
            write(21,'(g13.5)') fTxMT(i)
        enddo        
        write(cnum,'(i6)') nRxMT
        write(21,'(a,a)') '# MT Receivers: ', adjustl(cnum)   
        write(21,'(a1,9(a12,1x))') '!','X','Y','Z','Theta','Alpha','Beta','Length','SolveStatic','Name'
        do i = 1,nRxMT
            write(21,'(1x,7(f12.1,1x),i12,1x,a12)') xRxMT(i),yRxMT(i),zRxMT(i), ThetaRxMT(i),AlphaRxMT(i),BetaRxMT(i), &
                                       & lengthRxMT(i),iSolveStatic(i), trim(adjustl(cRxNamesMT(i)))
        enddo   
 
        
    endif    
    
 
   ! Data and Model response section:
    
    write(cnum,'(i6)') nd      
    write(21,'(a,a)')'# Data: ', adjustl(cnum)   
    write(21,'(a1,4(a12,1x),4(a13,1x))') '!','Type','Freq#','Tx#','Rx#','Data','StdError','Response','Residual' 
    ncolumns = size(dp,2)
    do i = 1,nd
       if (sd(i) == 0) then  ! catch for dummy data file with 0 errors
            resid = 0
       else
            resid = (d(i)-dm(i))/sd(i)
       endif
       write(21,'(1x,4(i12,1x),3(g15.7,1x),g13.2)') (dp(i,j),j=1,ncolumns), d(i),sd(i), dm(i), resid
    enddo   
    
    close(21)
    
    end subroutine writeResponse
    