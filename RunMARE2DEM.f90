!=======================================================================
! MARE2DEM: Modeling with Adaptively Refined Elements for 2.5D EM
!
! Version: 10 March 2014
!
! A parallel goal-oriented adaptive finite element forward and inverse
! modeling code for electromagnetic fields from electric dipoles, magnetic
! dipoles and magnetotelluric sources in triaxially anisotropic conducting
! media. Iterative adaptive mesh refinement is accomplished using the
! goal-oriented error estimation method described in Key and Ovall (2011).
! Inversion is accomplished with Occam's method (Constable et al., 1987).    
!
!                 This work was supported by: 
!
!            The Seafloor Electromagnetic Methods Consortium   
!               at Scripps Institution of Oceanography  
!
!         See this URL for a list of current SEMC funding sponsors:
!
!               http://marineemlab.ucsd.edu/semc.html
!
!    Copyright (C) 2008-2014 
!    Kerry Key
!    Scripps Institution of Oceanography
!    University of California, San Diego
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
    program runMARE2DEM
 
    use Occam  
    use mare2d_mpi_definitions    ! for MPI communicator tags etc
    use mare2dem_global           ! for setting defaults and storing worker status array
    use mare2dem_io
    use mare2dem_input_data_params
    
    implicit none

    logical         :: lFwdOnly = .false.       ! use command line argument -F to set this to true. 
    logical         :: lCompJacobian = .false.  ! use command line argument -FJ to set this to true
    logical         :: lrunMARE2DEM
    real(8)         :: timeOffset, time0 ! Timer temp
    integer         :: ierr, nproc, myID !  Local MPI variables
    character(12)   :: ctemp   
    
!
! Initialize MPI:
!
    call mpi_init( ierr )
    
!
! Start the timer:
! 
    call get_time_offset(0d0, time0)
    
!
! Get each processor's rank:
!
    call mpi_comm_rank( mcomm, myID, ierr )
!
! For Intel math kernel library, set #threads to 1. Occam matrix routines will override this since they are done 
! after parallel MARE2DEM calls have finished: 
#if defined(__INTEL_COMPILER)  
    call mkl_set_num_threads ( 1 )  ! Intel MKL specific  
! 25 Feb 2013: note that threading is only used on small machines, MARE2DEM now calls ScaLAPACK
! if a large number of processors are being used. This is handled in module Occam...
#endif
   

!
! Get total number of processes and number of workers:
!
    call mpi_comm_size( MPI_COMM_WORLD, nProc, ierr )   
    if (nProc == 1) then
        write(*,*) ' '
        write(*,*) '!!!!! Error: MARE2DEM needs to be run using at least 2 processors !!!'
        write(*,*) '             You specified only one, sorry, stopping!'
        write(*,*) ' '
        call exitMARE2DEM()
    endif
    nWorkers = nproc - 1   

!
! Launch the worker controllers:
!
    if (myID /= manager) call mpi_worker_controller

!
! The manager node does everything here and the workers are off listening for forward modeling task assignments...
!
 
!
! Launch the Manager controller:
!
    if (myID == manager) then 

    call displayBanner()    

    !
    ! Allocate the worker status array:
    !
    allocate (lWorker_status(nworkers))      
    lWorker_status  = .true.  
    
    !
    ! Tell the user about the cluster:
    !
    write(*,*) ' ' 
    write(ctemp,'(i8)') nworkers
    write(*,fmt='(a,a,a)') ' MARE2DEM is using one manager node and ',trim(adjustl(ctemp)), ' compute nodes'
    write(*,*) ' '     

    
    !
    ! Get any command line arguments:
    !
    call getMARE2DEMCommandLineArguments(lFwdOnly,lCompJacobian,currentIteration)    
    
    !
    ! Read in the model files (.resistivity, .poly,  .settings and .penalty):
    !
    call readModel 
    
    if (nParams == 0) then
        lFwdOnly = .true. ! No free parameters to invert for...
        write(*,*) ' ' 
        write(*,*) ' No free parameters found in the input model, therefore I will only '
        write(*,*) ' compute the forward response of the input model. '
        write(*,*) ' ' 
    endif
    write(*,fmt='(a32,a)') 'Scratch folder:  ',trim(adjustl(scratchFolder))    
    write(*,*) ' '    
    ! 
    ! Read in the data file:
    ! 
    call readData
        
    !
    ! Forward Call:
    !
    if (lFwdOnly) then
    
        if (lCompJacobian) then
            allocate( wj(nd,nParams) , stat = ierr)
            if (ierr .ne. 0) then
                write(*,*) ' Out of memory.  Too many data and free parameters (', nd,nParams, ')'
                call exitMARE2DEM()
            endif
            wj = 0
            
            !
            ! Initialize the scratch folders
            ! 
            if (nTxCSEM > 0) call mpi_manager_initializeScratch   
        
        
        endif
        
        !
        ! Compute forward response:
        ! 
        call computeFwd( lCompJacobian, transformToBound(pm,lowerBound,upperBound,lBoundMe)) 
        
        ! Show the misfit  
        write(*,*) ' '
        if (any(abs(d)>0)) then ! Only show misfit if input data is sensible:
            write(ctemp,'(f12.2)')   sqrt(sum( ((d-dm)/sd)**2 )/nd )  
            call print2col32('Forward Model Misfit:', ctemp,6)
        endif 
        call get_time_offset( time0,timeOffset)
        write(ctemp,'(f8.2)')  timeOffset
        call print2col32(' Total Time (s):', ctemp,6)
       
            
        ! Write out the model response:
        call writeResponse(currentIteration) 
        write(*,*) ' '
       
        ! Write out the model Jacobian matrix:
        if (lCompJacobian) then
            call writeJacobian(currentIteration) 
            deallocate (wj)
        endif     
        
    !
    ! Inversion Call:
    !    
    else 
     
        !
        ! Initialize the scratch folders
        ! 
        if (nTxCSEM > 0) call mpi_manager_initializeScratch   
        
        !
        ! Open a logfile:
        !
        call openOccamLogFile(outputFileRoot)
        
        !
        ! Initialize the inversion settings:
        !
        lrunMARE2DEM       = .true.
        
        !
        ! Run the main loop
        !
        do while (lrunMARE2DEM)
            
            currentIteration = currentIteration + 1 
            
            !
            ! Compute an Occam iteration:
            !
             call computeOccamIteration 
             
            !
            ! Write out the results:
            ! 
              if ( convergenceFlag /= 4 .and. convergenceFlag /= 5 ) then
            
                ! output an iteration file:
                call writeResistivityFile(currentIteration,modelRMS,targetRMS,modelMu,modelRoughness)
                
                ! output the model response:
                call writeResponse(currentIteration)            

            endif

            !
            ! Check the convergence flag:
            ! 
            ! Compute a single Occam iteration for profiling purposes, the below line is the
            ! orginal if statement
            !if (convergenceFlag > 1) then
            if (.true.) then
                exit
            end if

        enddo
        
        ! Close the Occam log file:
        call closeOccamLogFile


    endif
    
    !
    ! Deallocate memory
    !
        write(*,*) 'Deallocating memory...'
        call deallocateOccam
        call mare2dem_deallocate
            
    !
    ! Shut down nicely:
    !
        call exitMARE2DEM()
        
    !
    ! End of Manager node section  
    !
    endif ! (myid == 0)
   
    call mpi_finalize(ierr)
    
end program runMARE2DEM



!==================================================================================================================================! 
!====================================================================================================================== exitMARE2DEM
!==================================================================================================================================!    
    subroutine exitMARE2DEM()
!    
! Handles a graceful and elegant exit from MARE2EM rather than letting 
! chaos ensue 
!    
    integer :: ierr
    
!
! Shut down all workers:
!
    call mpi_shutDownWorkers
    
!
! Finalize MPI
!    
    call mpi_finalize(ierr)
   
!
! Stop...hammer time:
! 
    stop
  
! ..can't touch this...
    
    end subroutine exitMARE2DEM

!==================================================================================================================================! 
!=================================================================================================== getMARE2DEMCommandLineArguments
!==================================================================================================================================!      
    subroutine getMARE2DEMCommandLineArguments(lFwdOnly,lCompJacobian,nCurrentIter)
    !
    ! Subroutine to get command line arguments for MARE2DEM.
    ! 
    use mare2dem_global, only : resistivityFile, outputFileRoot, lFwdFields
    use string_helpers
 
    logical, intent(out)   :: lFwdOnly, lCompJacobian
    integer, intent(out)   :: nCurrentIter
    
    
    character(256) :: arg, sExt1, rfile='', ofile='', rbase
    integer        :: n, iExt1, idot
    logical        :: lokay1, lBad=.false.
        
    n = command_argument_count()
    
    if ( (n == 0) .or. (n > 3) ) then
        write(*,*) ' '
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
        write(*,*) 'MARE2DEM error, no command line arguments given!' 
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
        write(*,*) ' '
        call mare2dem_help
    endif
    
    ! 
    ! n > 0:
    !
    call get_command_argument(1, arg)
    
    select case (arg)
    
        case ('?') ! Help info requested:
            call mare2dem_help
            
        case ('-f','-F','/F','/f',char(92)//'F',char(92)//'f') ! Forward run only
            lFwdOnly = .true.
            
            if (n == 1) then
                write(*,*) ' '
                write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
                write(*,*) 'MARE2DEM error, no command line arguments given!' 
                write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
                write(*,*) ' '
                call mare2dem_help
            endif
            !write(*,*) ' '
            write(*,*) ' -F argument detected, MARE2DEM forward only call'
            !write(*,*) ' '
            if ( n  == 2 ) then
                call get_command_argument(2, rfile)
            elseif ( n  == 3 ) then
                call get_command_argument(2, rfile)
                call get_command_argument(3, ofile)
            endif
        case ('-fj','-FJ','-J','/FJ','/fj',char(92)//'FJ',char(92)//'fj') ! Forward run + Jacobian computation only
            lFwdOnly = .true.
            lCompJacobian = .true.
            if (n == 1) then
                write(*,*) ' '
                write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
                write(*,*) 'MARE2DEM error, no command line arguments given!' 
                write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
                write(*,*) ' '
                call mare2dem_help
            endif
            !write(*,*) ' '
            write(*,*) ' -FJ argument detected, MARE2DEM forward + Jacobian matrix only call'
            !write(*,*) ' '
            if ( n  == 2 ) then
                call get_command_argument(2, rfile)
            elseif ( n  == 3 ) then
                call get_command_argument(2, rfile)
                call get_command_argument(3, ofile)
            endif        
        case ('-ff','-FF','/FF','/ff',char(92)//'FF',char(92)//'ff') ! Forward run + all field components computed
            lFwdOnly   = .true.
            lFwdFields = .true. 
            if (n == 1) then
                write(*,*) ' '
                write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
                write(*,*) 'MARE2DEM error, no command line arguments given!' 
                write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
                write(*,*) ' '
                call mare2dem_help
            endif
            !write(*,*) ' '
            write(*,*) ' -FF argument detected, MARE2DEM forward only call with all fields computed'
            write(*,*) '    - MT fields will be output to <'
            !write(*,*) ' '
            if ( n  == 2 ) then
                call get_command_argument(2, rfile)
            elseif ( n  == 3 ) then
                call get_command_argument(2, rfile)
                call get_command_argument(3, ofile)
            endif               
        case default
            rfile = arg
            if ( n == 2 ) then 
                call get_command_argument(2, ofile)
            endif 
            
    end select   
    
    !
    ! Now parse the input names to remove any .resistivity extensions:
    !
    rfile = adjustl(rfile)
    ofile = adjustl(ofile)
    
    ! Remove .resistivity from filename:        
    resistivityFile = trim(rfile)
    idot = index(resistivityFile,'.',.true.)
    if (idot > 0) then 
        sExt1 = resistivityFile(idot+1:)
        if (trim(sExt1) =='resistivity')  resistivityFile = resistivityFile(1:idot-1)
    endif
    
    ! At this point there is no resistivity extension. Now get the iteration number:
    idot = index(resistivityFile,'.',.true.)
    if (idot > 0) then 
        sExt1 = resistivityFile(idot+1:)
        call convertStrToInteger(sExt1,iExt1,lokay1)
        if (lokay1) then ! extension is an integer, omit it
            nCurrentIter = iExt1
            rbase = resistivityFile(1:idot-1)
        else
            lBad = .true.
        endif
    else
        lBad = .true.
    endif                 
    if (lBad) then
        write(*,*) ' ' 
        write(*,*) 'Error with command line arguments, no iteration number in resistivity file: '
        write(*,*) trim(adjustl(rfile))
        write(*,*) ' ' 
        write(*,*) 'Try again!'
        call exitMARE2DEM()
    endif
 
    if (len_trim(ofile) == 0)  then
        outputFileRoot = rbase
    else
        outputFileRoot = ofile
    endif
    
    
        
 
    write(*,*) ' '
    write(*,*) ' Input resistivity file name root: ', trim(resistivityFile)
    write(*,*) 'Output resistivity file name root: ', trim(outputFileRoot)
    write(*,*) ' '

    end subroutine  getMARE2DEMCommandLineArguments

!==================================================================================================================================! 
!===================================================================================================================== mare2dem_help
!==================================================================================================================================! 
    subroutine mare2dem_help
    
    write(*,*) ' '
    write(*,*) 'Usage:  MARE2DEM <-F, -FJ or -FF> <resistivity file> <output root name>'
    write(*,*) ' '
    write(*,*) 'MARE2DEM has three optional flags (*normal inversion uses no flags):'
    write(*,*) ' '
    write(*,*) '      -F    Computes the forward response of the input model only.'
    write(*,*) '            The forward response is output to <resistivity file>.resp'
    write(*,*) ' '
    write(*,*) '      -FJ   Computes the forward response and the model Jacobian matrix.' 
    write(*,*) '            The forward response is output to <resistivity file>.resp'
    write(*,*) '            The Jacaobian matrix is output to <resistivity file>.jacobian'      
    write(*,*) ' '
    write(*,*) '      -FF   Computes the forward response of the input model and '
    write(*,*) '            outputs all 3 components of E and H for all combinations of '
    write(*,*) '            transmitters, receivers and frequencies in the data file.  '
    write(*,*) '            Ignores the data array section of the data file. The fields  '    
    write(*,*) '            The fields are output to <resistivity file>.fieldsMT and '      
    write(*,*) '            and <resistivity file>.fieldsCSEM '    
    write(*,*) ' '
    write(*,*) 'MARE2DEM has one required parameter:'
    write(*,*) ' '    
    write(*,*) '      <resistivity file root> - This is the name of the input resistivity  '
    write(*,*) '      file. By convention, this file should have the extension .resistivity. '
    write(*,*) '      For example inputModel.0.resistivity.  The model found by each'
    write(*,*) '      inversion iteration is then output to a new resistivity file with the '
    write(*,*) '      iteration number incremented.  For example: inputModel.1.resistivity, '
    write(*,*) '      inputModel.2.resistivity, ... The corresponding model responses are '
    write(*,*) '      written to inputModel.1.resp, inputModel.2.resp,... '
    write(*,*) ' ' 
    write(*,*) 'MARE2DEM has the optional parameter:'  
    write(*,*) ' '              
    write(*,*) '      <output file root> - With this option, the output files are '
    write(*,*) '      named <outputfileroot>.1.resistivity, <outputfileroot>.1.resp, ' 
    write(*,*) '      named <outputfileroot>.2.resistivity, <outputfileroot>.2.resp,... ' 
    write(*,*) ' ' 
    
    call exitMARE2DEM()
     
    end subroutine mare2dem_help
    
!==================================================================================================================================! 
!===================================================================================================================== displayBanner
!==================================================================================================================================!       
    subroutine displayBanner
 
    implicit none
    
    write(*,*) ' ' 
    write(*,*) '============================= MARE2DEM ==================================='
    write(*,*) ' '  
    write(*,*) ' MARE2DEM: Modeling with Adaptively Refined Elements for 2.5D EM'
    write(*,*) ''
    write(*,*) ' Version: 10 March 2014'
    write(*,*) ''
    write(*,*) ' A parallel goal-oriented adaptive finite element forward and inverse'
    write(*,*) ' modeling code for electromagnetic fields from electric dipoles, magnetic'
    write(*,*) ' dipoles and magnetotelluric sources in triaxially anisotropic conducting'
    write(*,*) ' media. Iterative adaptive mesh refinement is accomplished using the'
    write(*,*) ' goal-oriented error estimation method described in Key and Ovall (2011).'
    write(*,*) ' Inversion is accomplished by the Occam method (Constable et al., 1987).  '
    write(*,*) ''
    write(*,*) ' This work was supported by: '
    write(*,*) ''
    write(*,*) ' The Seafloor Electromagnetic Methods Consortium   '
    write(*,*) ' at Scripps Institution of Oceanography  '
    write(*,*) ''
    write(*,*) ' See this URL for a list of current SEMC funding sponsors:'
    write(*,*) ''
    write(*,*) ' http://marineemlab.ucsd.edu/semc.html'
    write(*,*) ''
    write(*,*) ' Copyright (C) 2008-2014 '
    write(*,*) ' Kerry Key'
    write(*,*) ' Scripps Institution of Oceanography'
    write(*,*) ' University of California, San Diego'
    write(*,*) ''
    write(*,*) ' This file is part of MARE2DEM.'
    write(*,*) ''
    write(*,*) ' MARE2DEM is free software: you can redistribute it and/or modify'
    write(*,*) ' it under the terms of the GNU General Public License as published by'
    write(*,*) ' the Free Software Foundation, either version 3 of the License, or'
    write(*,*) ' (at your option) any later version.'
    write(*,*) ''
    write(*,*) ' MARE2DEM is distributed in the hope that it will be useful,'
    write(*,*) ' but WITHOUT ANY WARRANTY; without even the implied warranty of'
    write(*,*) ' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
    write(*,*) ' GNU General Public License for more details.'
    write(*,*) ''
    write(*,*) ' You should have received a copy of the GNU General Public License'
    write(*,*) ' along with MARE2DEM.  If not, see <http://www.gnu.org/licenses/>. '
    write(*,*) ' ' 
    write(*,*) '==========================================================================' 
    !write(*,*) ' '
    end subroutine displayBanner
    