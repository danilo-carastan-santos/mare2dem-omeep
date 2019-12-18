!-----------------------------------------------------------------------
!
!    Copyright 2008-2014
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
!============================================================================================================ mare2d_mpi_definitions
!==================================================================================================================================!
    module mare2d_mpi_definitions

#if defined(__INTEL_COMPILER)
    use IFPORT          ! DGM 2/23/2012 for the sleepqq routine
#endif
    
!#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
!    include 'mpif.h'
!#else
    use mpi
!#endif

!
! MPI definitions common to manager and Workers:
!
    
    integer, parameter :: manager = 0        ! manager is 0, all others are Workers
    integer, parameter :: mcomm = MPI_COMM_WORLD
    integer            :: status(MPI_STATUS_SIZE)   
              
    integer, parameter :: tag_M_SendingData              = 1     ! tags for manager sending a flag data    
    integer, parameter :: tag_W_SendingData              = 2
            
    ! iflags are like tags, but i put them in the data field of the mpi command
    integer, parameter :: iflag_TellWorkerToQuit            = 0
    integer, parameter :: iflag_TellWorkerToRunSubset       = 1   
    integer, parameter :: iflag_TellWorkerToScaLapackMM     = 2      
    integer, parameter :: iflag_TellWorkerToScaLapackMS     = 3
    integer, parameter :: iflag_TellWorkerToScratch         = 4
       
    end module mare2d_mpi_definitions

 
!==================================================================================================================================! 
!====================================================================================================================== mpi_mare2dem
!==================================================================================================================================!
    subroutine mpi_mare2dem()
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu 
!
    use mare2d_mpi_definitions
    use mare2dem_global    
    use mare2dem_input_data_params
        
    implicit none
 
    logical :: lkeepgoing, lflag
    
    integer :: ierr,  iWorker, nProc

    call mpi_comm_size( MPI_COMM_WORLD, nProc, ierr )  
         
!
! Set up the wavenumbers, starting mesh and the refinement group arrays:
!
    call getStartingMesh   
    call getWavenumbers    
    call getRxTxGroups
    call getFqGroups   
    call getRefinementGroups
    
    call printDecomposition 
    
    if (lPrintSetup)   call display_MARE2DEM_Params
    
!
! Launch the manager controller:
!
    do    
                
        status = 0
        
        !
        ! Check for a message from any Worker:
        !
        call  mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mcomm, lflag, status, ierr)
        
        !
        ! If message from Worker, then decode it:
        !
        if (lflag) then 
        
            iWorker = status(MPI_SOURCE)

            if  (status(MPI_TAG) == tag_W_SendingData) then
                
                !
                ! Receive data from Worker (results for this Worker's subset)
                !
                call mpi_manager_receive_results(iWorker)  
        
                                     
            endif 
            
        endif ! lflag 
                
        !
        ! Try sending out a job to the first available worker in lWorker_status
        !
        call mpi_manager_send_job()
        
        !
        ! David Myer's tweak for idling the manager, which can be useful when running MARE2DEM on a small
        ! multicore desktop or laptop. This keeps the manager from consuming 100% cpu while waiting for
        ! the workers to finish:
        !
#if defined(__INTEL_COMPILER)  
        if (nProc <= 12) call sleepqq(1)
#endif
   
        !
        ! Check to see if we're done with all the models.  If so then terminate the worker controllers.
        !
        lkeepgoing = .true.
        if ( all( lWorker_status ) ) then
             lkeepgoing = .false.
             exit
        endif
     
        
    enddo ! manager controller loop
        
!
! Done with the adaptive mesh refinements, now let's extract the responses:
!
 
!
! CSEM:
!
    if ( nfreqCSEM > 0 ) call applyPhaseConvention 
  
!
! Deallocate all internal arrays:
!
    deallocate( RxTxGroups, FqGroups, refinementGroups )
    ! this also deallocates the allocatable fields within the derived types
    
    if (allocated(wavenum))         deallocate (wavenum)

    if (allocated(inputmesh%attr) ) call deallocate_trimesh(inputmesh,.false.)
    
    if (lprintDebug) write(*,*) 'leaving mpi_mare2dem '
  
    
    end subroutine mpi_mare2dem

!==================================================================================================================================! 
!==================================================================================================== mpi_manager_initializeScratch
!==================================================================================================================================!
    subroutine mpi_manager_initializeScratch

 
    use mare2d_mpi_definitions
    use mare2dem_global    

    implicit none
    
    integer :: ierr, iWorker 

    
    do iWorker = 1,nWorkers
        call mpi_send( iflag_TellWorkerToScratch,   1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm, ierr )
        call mpi_send( scratchFolder,             256, MPI_CHARACTER, iWorker, tag_M_SendingData, mcomm,ierr )        
        
    enddo
    
    end subroutine mpi_manager_initializeScratch


!==================================================================================================================================! 
!============================================================================================================== mpi_manager_send_job
!==================================================================================================================================!
    subroutine mpi_manager_send_job

 
    use mare2d_mpi_definitions
    use mare2dem_global    

    implicit none
    
    integer :: i, iWorker
    
    !
    ! Send out work until there are no more workers are available:
    !
    
    do 
        
        !
        ! Find first available worker:
        !
        iWorker = 0
        do i = 1,nWorkers
            if (lWorker_status(i)) then
                iWorker = i
                exit
            endif
        enddo

        if (iWorker == 0) return ! no workers available

        
        !
        ! Now find a job for that worker:
        !
        if ( iPtr_refGroups <= nRefinementGroups ) then

            call mpi_manager_send_subset(iWorker, iPtr_refGroups)            
            iPtr_refGroups = iPtr_refGroups + 1
            cycle
      
        endif    
     
        !
        ! There are no tasks to do, return to main loop
        !
        return
        
    enddo
    
    end subroutine mpi_manager_send_job

!==================================================================================================================================! 
!=========================================================================================================== mpi_manager_send_subset
!==================================================================================================================================!
   subroutine mpi_manager_send_subset(iWorker,iGrp)
   
    use mare2d_mpi_definitions
    use em2dkx_mod 
    use mare2dem_global
    use mare2dem_input_data_params
    
    implicit none
    
    integer, intent(in) :: iWorker, iGrp
    
    !
    ! Local:
    !
    integer         :: nRR, nTT, ierr ,nnod, nele, nn, iRxTxG, iFqG, nFq, rx0,rx1,tx0,tx1,fq0,fq1
    real(8)         :: tStart, tEnd
    character(256)  :: sFmt
    character(32)   :: stime
       
    call cpu_time(tStart) 
    
    lWorker_status(iWorker) = .false. 
    
    if (lprintDebug) write(*,'(a,2(i6,1x))') 'manager sending job, iGrp, iWorker: ', iGrp, iWorker
    
    !
    ! Tell the worker that your going to send it a model to run:
    !
    call mpi_send( iflag_TellWorkerToRunSubset, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm, ierr )
 
    !
    ! Tell the worker which refinement group is being sent:
    !
    call mpi_send( iGrp,     1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr ) 
 
    iRxTxG  = refinementGroups(iGrp)%iRxTxGroup 
    
    nRR     = RxTxGroups(iRxTxG)%nRx          
    rx0     = RxTxGroups(iRxTxG)%iRx(1)
    rx1     = RxTxGroups(iRxTxG)%iRx(nRR)
    
    nTT     = RxTxGroups(iRxTxG)%nTx  
    tx0     = RxTxGroups(iRxTxG)%iTx(1)    ! note for MT, iTx(1) = 0
    tx1     = RxTxGroups(iRxTxG)%iTx(nTT)

    iFqG    = refinementGroups(iGrp)%iFqGroup
        
    nFq     = FqGroups(iFqG)%nFq
    fq0     = FqGroups(iFqG)%iFq(1)
    fq1     = FqGroups(iFqG)%iFq(nFq)
   
    lMT     = refinementGroups(iGrp)%lMT
    
    call mpi_send( lMT,      1, MPI_LOGICAL, iWorker, tag_M_SendingData, mcomm,ierr )
              
    if ( .not. lMT ) then ! send CSEM info:

        !
        ! Send source positions and directions:
        !
        call mpi_send( nTT,                                     1, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( xTxCSEM(RxTxGroups(iRxTxG)%iTx),       nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( yTxCSEM(RxTxGroups(iRxTxG)%iTx),       nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( zTxCSEM(RxTxGroups(iRxTxG)%iTx),       nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( azimuthTxCSEM(RxTxGroups(iRxTxG)%iTx), nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( dipTxCSEM(RxTxGroups(iRxTxG)%iTx),     nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( lengthTxCSEM(RxTxGroups(iRxTxG)%iTx),  nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( cSourceType(RxTxGroups(iRxTxG)%iTx), 8*nTT, MPI_CHARACTER,       iWorker,tag_M_SendingData,mcomm,ierr )

        !
        ! Send receivers
        !
        call mpi_send( nRR,                                    1, MPI_INTEGER,          iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( xRxCSEM(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( yRxCSEM(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( zRxCSEM(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( ThetaRxCSEM(RxTxGroups(iRxTxG)%iRx),  nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( AlphaRxCSEM(RxTxGroups(iRxTxG)%iRx),  nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( BetaRxCSEM(RxTxGroups(iRxTxG)%iRx),   nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( lengthRxCSEM(RxTxGroups(iRxTxG)%iRx), nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        !
        ! Send wavenumbers:
        ! 
        call mpi_send(   nwave,      1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )  
        call mpi_send( wavenum,  nwave, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        
        !
        ! Send frequencies:
        !
        nFq = FqGroups(iFqG)%nFq
        call mpi_send(   nFq,                              1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr ) 
        call mpi_send( fTxCSEM(FqGroups(iFqG)%iFq),  nFq, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )

        call mpi_send( nKxPerGroup,       1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( nQuadRxCSEM,       1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( nQuadTxCSEM,       1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
 
        !
        ! Send the data masks:
        !
        call mpi_send( lDataMaskCSEM(rx0:rx1,fq0:fq1,tx0:tx1), nRR*nTT*nFq, MPI_LOGICAL, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( lCompCSEM(1:6,rx0:rx1),                       6*nRR, MPI_LOGICAL, iWorker, tag_M_SendingData,mcomm,ierr)
 
    else  ! MT  
        
        !
        ! Send receivers
        !
        call mpi_send( nRR,                                  1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( xRxMT(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( yRxMT(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( zRxMT(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( ThetaRxMT(RxTxGroups(iRxTxG)%iRx),  nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( AlphaRxMT(RxTxGroups(iRxTxG)%iRx),  nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( BetaRxMT(RxTxGroups(iRxTxG)%iRx),   nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( lengthRxMT(RxTxGroups(iRxTxG)%iRx), nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )

        !
        ! Send frequencies:
        !        
        nFq = FqGroups(iFqG)%nFq
        call mpi_send(   nFq,                              1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )  
        call mpi_send( fTxMT(FqGroups(iFqG)%iFq),    nFq, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        
        call mpi_send( nQuadRxMT,         1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )    

        !
        ! Send the data masks:
        ! 
        call mpi_send( lDataMaskMT(rx0:rx1,fq0:fq1), nRR*1*nFq, MPI_LOGICAL, iWorker, tag_M_SendingData, mcomm, ierr )      
        call mpi_send( lCompMT(1:6,rx0:rx1)        ,     6*nRR, MPI_LOGICAL, iWorker, tag_M_SendingData, mcomm, ierr ) 
                           
    endif
 
    !
    ! Send the mesh: kwk debug: could send only model and have slaves compute initial triangulation...
    !    
    nnod   = inputmesh%nnod
    nele   = inputmesh%nele
    call mpi_send( nnod, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( nele, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
          
    call mpi_send( inputmesh%y,              nnod, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmesh%z,              nnod, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmesh%attr,           nele, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmesh%emap,         3*nele, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmesh%neighborlist, 3*nele, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
        
    call mpi_send( inputmesh%numberofpointattributes,    1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmesh%numberofcorners,            1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmesh%numberoftriangleattributes, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmesh%numberofsegments,           1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmesh%numberofholes,              1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmesh%numberofregions,            1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )


    nn =  inputmesh%numberofpointattributes*inputmesh%nnod
    if (nn > 0 ) call mpi_send( inputmesh%pointattributelist, nn, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )

    call mpi_send( inputmesh%pointmarkerlist,   nnod,                         MPI_INTEGER, iWorker, tag_M_SendingData,mcomm,ierr)
    call mpi_send( inputmesh%segmentlist,       2*inputmesh%numberofsegments, MPI_INTEGER, iWorker, tag_M_SendingData,mcomm,ierr)
    call mpi_send( inputmesh%segmentmarkerlist,   inputmesh%numberofsegments, MPI_INTEGER, iWorker, tag_M_SendingData,mcomm,ierr)

    if ( inputmesh%numberofholes > 0 ) then
        call mpi_send( inputmesh%holelist, 2*inputmesh%numberofholes, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr)
    endif

    call mpi_send( inputmesh%regionlist, 4*inputmesh%numberofregions, MPI_DOUBLE_PRECISION, iWorker,  &
          & tag_M_SendingData, mcomm,ierr )

    ! send fixed and free resistivity values:
    call mpi_send( cAnisotropy,         9, MPI_CHARACTER,        iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( nSigParams ,         1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( sigParams , nSigParams, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( iFreeParam, nSigParams, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )  
    call mpi_send( nFree ,              1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    
    call mpi_send( maxnadapt,                            1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( errortolerance,  1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( idual_func,      1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( minQangle,       1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )

    call mpi_send( outputFileRoot,         256, MPI_CHARACTER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send(lSaveMeshFiles,            1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )  
    call mpi_send(linversion,                1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )  
    call mpi_send(lDisplayRefinementStats ,  1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )

    
    call mpi_send( linearSolver,  32, MPI_CHARACTER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( lprintDebug,    1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )  
    call mpi_send( lUseBumpFields, 1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )      

    call mpi_send( minRangeProtector, 1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( minArea,           1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( ecutoff,           1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( hcutoff,           1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( max_nsubrefine,    1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( pct_refine,        1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
      
 
    !
    ! Worker will then compute the 2D EM response and send its results back to the manager...
    !           
    call cpu_time(tEnd)     
    write(stime,'(a8,f9.4,a2)') ' Timer: ',tEnd - tStart, ' s'
    sFmt = '(a5,2x,i6,2x,       a6,2x,i6,2x,    69x , a32,  2x,a19)'    
    if (lprintMPItimers) write(*,sFmt)  'Proc:',0,'Group:',iGrp,   'mpi_manager_send_subset:',   trim(stime)

    end subroutine mpi_manager_send_subset
!==================================================================================================================================!
!=========================================================================================================== mpi_manager_send_subset
!==================================================================================================================================!
    subroutine mpi_manager_receive_results(iWorker)
 
    use em2dkx_mod 
    use mare2d_mpi_definitions
    use mare2dem_global
    use mare2dem_output   
    
    implicit none
    
    integer, intent(in) :: iWorker
    
    integer :: iGrp, iRxTxG, iFqG, nTT, nRR, nFq, nrec, iRx,iTx,iFq, ierr, rx0,rx1,tx0,tx1,fq0,fq1
    
    real(8)         :: tStart, tEnd
    character(256)  :: sFmt
    character(32)   :: stime
    
    call cpu_time(tStart) 
!
! First receive the refinement group index:
!
    call mpi_recv( iGrp,   1, MPI_INTEGER, iWorker, tag_W_SendingData , mcomm, status, ierr )
 
    iRxTxG  = refinementGroups(iGrp)%iRxTxGroup 
    
    nRR     = RxTxGroups(iRxTxG)%nRx          
    rx0     = RxTxGroups(iRxTxG)%iRx(1)
    rx1     = RxTxGroups(iRxTxG)%iRx(nRR)
    
    nTT     = RxTxGroups(iRxTxG)%nTx  
    tx0     = RxTxGroups(iRxTxG)%iTx(1)    ! note for MT, iTx(1) = 0
    tx1     = RxTxGroups(iRxTxG)%iTx(nTT)

    iFqG    = refinementGroups(iGrp)%iFqGroup
        
    nFq     = FqGroups(iFqG)%nFq
    fq0     = FqGroups(iFqG)%iFq(1)
    fq1     = FqGroups(iFqG)%iFq(nFq)
   
    lMT     = refinementGroups(iGrp)%lMT
    
    
    if (.not.lMT) then
      
      nrec    = nTT*nFq
        
      do iRx = rx0,rx1
        
          if (lCompCSEM(1,iRx)) call mpi_recv(ex(iRx,tx0:tx1,fq0:fq1),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
          if (lCompCSEM(2,iRx)) call mpi_recv(ey(iRx,tx0:tx1,fq0:fq1),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
          if (lCompCSEM(3,iRx)) call mpi_recv(ez(iRx,tx0:tx1,fq0:fq1),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
          if (lCompCSEM(4,iRx)) call mpi_recv(hx(iRx,tx0:tx1,fq0:fq1),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
          if (lCompCSEM(5,iRx)) call mpi_recv(hy(iRx,tx0:tx1,fq0:fq1),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
          if (lCompCSEM(6,iRx)) call mpi_recv(hz(iRx,tx0:tx1,fq0:fq1),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)  
      
      enddo 
    
      if (linversion) then
      
       do iFq = fq0,fq1  
         do iTx = tx0,tx1
           do iRx = rx0,rx1
                
                if ( lDataMaskCSEM(iRx,iFq,iTx) ) then
                
                if ( lCompCSEM(1,iRx) ) then
                    call mpi_recv( dex(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX,  iWorker, tag_W_SendingData , mcomm, status,ierr )
                endif
                if ( lCompCSEM(2,iRx) ) then
                    call mpi_recv( dey(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX,  iWorker, tag_W_SendingData , mcomm, status,ierr )
                endif
                if ( lCompCSEM(3,iRx) ) then
                    call mpi_recv( dez(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX,  iWorker, tag_W_SendingData , mcomm, status,ierr )
                endif      
                if ( lCompCSEM(4,iRx) ) then
                    call mpi_recv( dhx(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX,  iWorker, tag_W_SendingData , mcomm, status,ierr )
                endif
                if ( lCompCSEM(5,iRx) ) then
                    call mpi_recv( dhy(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX,  iWorker, tag_W_SendingData , mcomm, status,ierr )
                endif
                if ( lCompCSEM(6,iRx) ) then
                    call mpi_recv( dhz(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX,  iWorker, tag_W_SendingData , mcomm, status,ierr )
                endif                 
                
                endif
                
            enddo
          enddo
        enddo
 
      endif
      
    else
      
 ! MT:
    
    do iRx = rx0,rx1
          
      if (lCompMT(1,iRx)) call mpi_recv(ex_mt(iRx,fq0:fq1),nFq,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
      if (lCompMT(2,iRx)) call mpi_recv(ey_mt(iRx,fq0:fq1),nFq,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
      if (lCompMT(3,iRx)) call mpi_recv(ez_mt(iRx,fq0:fq1),nFq,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
      if (lCompMT(4,iRx)) call mpi_recv(hx_mt(iRx,fq0:fq1),nFq,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
      if (lCompMT(5,iRx)) call mpi_recv(hy_mt(iRx,fq0:fq1),nFq,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
      if (lCompMT(6,iRx)) call mpi_recv(hz_mt(iRx,fq0:fq1),nFq,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)  
    
    enddo 
    
      if (linversion) then
          
       nrec  = nFq*nFree
       
       do iRx = rx0,rx1
      
        if (lCompMT(1,iRx)) call mpi_recv(dex_mt_dsig(iRx,fq0:fq1,:),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
        if (lCompMT(2,iRx)) call mpi_recv(dey_mt_dsig(iRx,fq0:fq1,:),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
        if (lCompMT(3,iRx)) call mpi_recv(dez_mt_dsig(iRx,fq0:fq1,:),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
        if (lCompMT(4,iRx)) call mpi_recv(dhx_mt_dsig(iRx,fq0:fq1,:),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
        if (lCompMT(5,iRx)) call mpi_recv(dhy_mt_dsig(iRx,fq0:fq1,:),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)
        if (lCompMT(6,iRx)) call mpi_recv(dhz_mt_dsig(iRx,fq0:fq1,:),nrec,MPI_COMPLEX,iWorker,tag_W_SendingData,mcomm,status,ierr)  
      
       enddo
            
      endif    
    
    endif
 
     
    call cpu_time(tEnd)     
    write(stime,'(a8,f9.4,a2)') ' Timer: ',tEnd - tStart, ' s'
    sFmt = '(a5,2x,i6,2x,       a6,2x,i6,2x,    a24,2x,i7,  14x,  54x,   2x,a19)'     
    if (lprintMPItimers) write(*,sFmt)  'Proc:',0,'Group:',iGrp,   'Manager recv results:',0,   trim(stime )
     
!
! Set the availability status for this worker:
!
    lWorker_status(iWorker) = .true.   
        
    end subroutine mpi_manager_receive_results 
        
!==================================================================================================================================! 
!============================================================================================================= mpi_worker_controller
!==================================================================================================================================!
    subroutine mpi_worker_controller

    use mare2d_mpi_definitions    
    use occam 
 
    
    implicit none
    
    logical :: lflag
    
    integer :: myID, ierr, iflag
    
    
    !
    ! Determine which worker I am:
    !
    call mpi_comm_rank( mcomm, myID, ierr )
    
    
    !
    ! Create an infinite loop for this worker controller:
    ! 
    do 

        !
        ! David Myer's tweak for idling the worker, which can be useful when running MARE2DEM on a small
        ! multicore desktop or laptop. This keeps the worker from consuming 100% cpu while waiting for
        ! a new task from the manager:
        !
#if defined(__INTEL_COMPILER)            
        !
        ! Check for a message from the manager:
        !
        call  mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mcomm, lflag, status, ierr)
        
        !
        ! If no message from the manager, take a valium and try later:
        !
        if (.not.lflag) then 
            call sleepqq(1) ! kwk debug: this is specific to the intel compiler
            cycle
        endif
#endif        
        !
        ! Get the command from the manager:
        !
        call mpi_recv( iflag, 1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
        
        !
        ! What is the manager telling me?
        !
        select case (iflag)
        
        case ( iflag_TellWorkerToQuit ) 
            write(*,'(a12,1x,i5,1x,a32)') ' Worker: ', myID, ' is done with mare2dem...'     
            exit    ! Worker now leaves this subroutine
           
        case ( iflag_TellWorkerToRunSubset )     
            
            call mpi_worker_run_subset()  
            
         case ( iflag_TellWorkerToScratch )     
            
            call mpi_worker_intializeScratch()                      
            
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
        case ( iflag_TellWorkerToScaLapackMM )     
            
            call mpi_worker_ScaLapackMM()

        case ( iflag_TellWorkerToScaLapackMS )     
            
            call mpi_worker_ScaLapackMS(ierr)
#endif
                           
        case default
            write(*,*) 'Error in Worker controller, received unknown command from the manager !'
            write(*,*) 'iflag is: ',iflag
            write(*,*) ' quit mpi and debug the code ...'
            stop
        end select
        
    enddo   ! the infinite while loop
    
    end subroutine mpi_worker_controller
    
!==================================================================================================================================! 
!======================================================================================================= mpi_worker_intializeScratch
!==================================================================================================================================!  
    subroutine mpi_worker_intializeScratch


    use mare2d_mpi_definitions    
    use mare2dem_global
   
    implicit none
 
   
    integer         :: myID, ierr
    
    character(256)  :: cFilename, strID
    integer         :: ios
    
    logical         :: lstatus
    !
    ! Determine which worker I am:
    !
    call mpi_comm_rank( mcomm, myID, ierr )

    !
    ! Get scratch folder base name:
    ! 
    call mpi_recv( scratchFolder,         256, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr)
             
             

    scratchFolder = trim(scratchFolder)//'/'  
 
    
    ! 
    ! Test for scratch folder existence by trying to write a dummy file:
    !
    write(strID,*) myID       
    cFilename  = trim(scratchFolder)//'test'//trim(adjustl(strID)) 
 
    open(unit=21,file=trim(cFilename),form="unformatted",status="replace",access="stream",iostat=ios)
              
    if (ios /= 0) then
    
        !
        ! No scratch folder, so create one:
        !

#if defined(__INTEL_COMPILER)  
 
        if( .not. makedirqq(trim(scratchFolder)) ) then
        
            write(*,'(a)') ' ------------------------------------------------------------------------------------------------------'
            write(*,'(a)') ' Error on processor: ',myID
            write(*,'(a)') ' MARE2DEM cannot create a new scratch folder. Stopping! '
            write(*,'(a,a)') ' The scratch folder MARE2DEM tried to use is: ',trim(scratchFolder)
            write(*,'(a)') ' '
            write(*,'(a)') ' What has gone wrong? '
            write(*,'(a)') '  (1) The default scratch folder does not exist on your cluster. '
            write(*,'(a)') '  or '
            write(*,'(a)') '  (2) You specified a scratch folder in the settings file, but it does not exist on your system.  '
            write(*,'(a)') '     Solution: create this folder before running MARE2DEM. '
            write(*,*)    ' '
            write(*,'(a)') ' Note that on large clusters, an existing local folder for each node should be used (for example /tmp) '
            write(*,'(a)') ' so that the scratch data is not transmitted across the network and instead only uses the local disk.   '
            write(*,'(a)') ' If you are running MARE2DEM on a single computer (laptop or desktop), then you can use a local '
            write(*,'(a)') ' directory for scratch space. '
            write(*,'(a)') ' Or if you are running MARE2DEM on cluster with a superfast parallel filesystem, you might be able to use'
            write(*,'(a)') ' a networked directory for scratch space, but a local directory is still probably more efficient. '            
            write(*,*)    ' '
            write(*,*)    ' So long and thanks for all the fish! '
            write(*,*)    ' '
            stop
        endif
 
!#elif defined(__GNUC__) 
!        
!        call system("mkdir scratch")
                  
#endif
    
    
    else
    
        !
        ! A-okay, folder exists, delete the dummy file:
        !
        close(21,status="delete")
        
#if defined(__INTEL_COMPILER)  
        
        ! Erase its contents:
        lstatus = delfilesqq(trim(scratchFolder)//'mare2dem_tmp*')
 
 
!#elif defined(__GNUC__) 
!        
!        call system("rmdir -fR scratch") ! this should work but i've not tested in on gfortran yet...
!        call system("mkdir scratch")
      
                  
#endif   
     
    endif 
   
 
    
    end subroutine mpi_worker_intializeScratch  

!==================================================================================================================================!
!============================================================================================================= mpi_worker_run_subset
!==================================================================================================================================!
    subroutine mpi_worker_run_subset

    use mare2d_mpi_definitions 
    use mare2dem_worker
    use mare2dem_global
    
    implicit none
 

!
! Receive data from the manager:
!
    call mpi_worker_receive_subset 
    
!
! Compute the 2D EM response:
!
    call worker_EM2D

!
! Send results back to the manager
!
    call mpi_worker_send_results 
    
!
! Deallocate:
!  
    call mare2dem_deallocate  ! this routine deallocates everything from the worker
    
    
    end subroutine mpi_worker_run_subset
    
!==================================================================================================================================! 
!========================================================================================================= mpi_worker_receive_subset
!==================================================================================================================================!
    subroutine mpi_worker_receive_subset

    use mare2d_mpi_definitions    
    use em2dkx_mod 
    use mare2dem_global
    use mare2dem_input_data_params
 
    implicit none
 
    integer        :: ierr, nele, nnod, nn 
 
    real(8)         :: tStart, tEnd
    character(256)  :: sFmt
    character(32)   :: stime
    
    call cpu_time(tStart) 
       
    call mpi_comm_rank( mcomm, myID, ierr )
    
    !
    ! Receive data for em2dkx from the manager:
    !
    call mpi_recv( iRefinementGrp,    1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( lMT,               1, MPI_LOGICAL, manager, tag_M_SendingData, mcomm, status, ierr )

    nRxCSEM   = 0
    nTxCSEM   = 0
    nRxMT     = 0
    nfreqCSEM = 0 
    nfreqMT   = 0
    
    if (.not.lMT) then ! receive CSEM arrays:
        
        !
        ! Receive the transmitters:
        !
        call mpi_recv( nTxCSEM,  1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  
        
        allocate ( azimuthTxCSEM(nTxCSEM),dipTxCSEM(nTxCSEM), &
                     & xTxCSEM(nTxCSEM),yTxCSEM(nTxCSEM),zTxCSEM(nTxCSEM), lengthTxCSEM(nTxCSEM), &
                     & cSourceType(nTxCSEM) )    
                     
        call mpi_recv( xTxCSEM ,        nTxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )  
        call mpi_recv( yTxCSEM ,        nTxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( zTxCSEM ,        nTxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( azimuthTxCSEM ,  nTxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( dipTxCSEM ,      nTxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( lengthTxCSEM ,   nTxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( cSourceType ,  8*nTxCSEM, MPI_CHARACTER,        manager, tag_M_SendingData, mcomm, status, ierr )      

        
        !
        ! Receive the receivers:
        !
        call mpi_recv( nRxCSEM,  1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  

        allocate( xRxCSEM(nRxCSEM), yRxCSEM(nRxCSEM), zRxCSEM(nRxCSEM), lengthRxCSEM(nRxCSEM),  &
                      ThetaRxCSEM(nRxCSEM),AlphaRxCSEM(nRxCSEM),BetaRxCSEM(nRxCSEM) )   
                              
        call mpi_recv( xRxCSEM ,      nRxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )  
        call mpi_recv( yRxCSEM ,      nRxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( zRxCSEM ,      nRxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( ThetaRxCSEM ,  nRxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( AlphaRxCSEM ,  nRxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( BetaRxCSEM ,   nRxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )    
        call mpi_recv( lengthRxCSEM , nRxCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )      
    
        !
        ! Receive the wavenumbers:
        !
        call mpi_recv( nwave,        1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
        allocate( wavenum(nwave) )
        call mpi_recv( wavenum,  nwave, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        
        ! 
        ! Receive the frequencies:
        !
        call mpi_recv( nFreqCSEM,        1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
        allocate( fTxCSEM(nFreqCSEM) )
        call mpi_recv( fTxCSEM,  nFreqCSEM, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        
        
        call mpi_recv( nKxPerGroup,  1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  
        call mpi_recv( nQuadRxCSEM,  1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  
        call mpi_recv( nQuadTxCSEM,  1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  
        
        !
        ! Receive the data masks:
        !
        allocate ( lDataMaskCSEM(nRxCSEM,nFreqCSEM,nTxCSEM), lCompCSEM(6,nRxCSEM) )
        
        call mpi_recv( lDataMaskCSEM, nRxCSEM*nFreqCSEM*nTxCSEM, MPI_LOGICAL, manager, tag_M_SendingData, mcomm, status, ierr )
        call mpi_recv(     lCompCSEM,                 6*nRxCSEM, MPI_LOGICAL, manager, tag_M_SendingData, mcomm, status, ierr )
 
 
    else ! receive MT arrays:
    
    
        call mpi_recv( nRxMT,  1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  

        allocate( xRxMT(nRxMT), yRxMT(nRxMT), zRxMT(nRxMT),ThetaRxMT(nRxMT), AlphaRxMT(nRxMT),BetaRxMT(nRxMT), &
                    &   lengthRxMT(nRxMT)  )
                      
        call mpi_recv( xRxMT ,      nRxMT, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )  
        call mpi_recv( yRxMT ,      nRxMT, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( zRxMT ,      nRxMT, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( ThetaRxMT ,  nRxMT, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( AlphaRxMT ,  nRxMT, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( BetaRxMT ,   nRxMT, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )    
        call mpi_recv( lengthRxMT , nRxMT, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )    

        ! 
        ! Receive the frequencies:
        !   
        call mpi_recv( nFreqMT,        1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr ) 
         
        allocate( fTxMT(nFreqMT) )
        call mpi_recv( fTxMT,    nFreqMT, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 

        call mpi_recv( nQuadRxMT,      1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  
        
        !
        ! Receive the data masks:
        !
        allocate ( lDataMaskMT(nRxMT,nFreqMT), lCompMT(6,nRxMT) )
        
        call mpi_recv( lDataMaskMT, nRxMT*nFreqMT, MPI_LOGICAL, manager, tag_M_SendingData, mcomm, status, ierr )
        call mpi_recv( lCompMT ,          6*nRxMT, MPI_LOGICAL, manager, tag_M_SendingData, mcomm, status, ierr )
        
        nKxPerGroup = 1 ! since there is only one
        
    endif
    
    
    !
    ! Receive the mesh:
    !
    call mpi_recv( inputmesh%nnod, 1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmesh%nele, 1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )

    nnod = inputmesh%nnod
    nele = inputmesh%nele

    allocate ( inputmesh%y(nnod),inputmesh%z(nnod) )
    allocate ( inputmesh%emap(3,nele), inputmesh%attr(nele) )
    allocate ( inputmesh%neighborlist(3,nele) ) 
    
    call mpi_recv( inputmesh%y,              nnod, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmesh%z,              nnod, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmesh%attr,           nele, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmesh%emap,         3*nele, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmesh%neighborlist, 3*nele, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )
 
   call mpi_recv( inputmesh%numberofpointattributes,    1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
   call mpi_recv( inputmesh%numberofcorners,            1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
   call mpi_recv( inputmesh%numberoftriangleattributes, 1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
   call mpi_recv( inputmesh%numberofsegments,           1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
   call mpi_recv( inputmesh%numberofholes,              1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
   call mpi_recv( inputmesh%numberofregions,            1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )

   nn = inputmesh%numberofpointattributes*inputmesh%nnod
   allocate (  inputmesh%pointattributelist(nn), inputmesh%pointmarkerlist(inputmesh%nnod) ) 
   inputmesh%pointattributelist = 0
 
   if (nn > 0 ) call mpi_recv( inputmesh%pointattributelist, nn, MPI_DOUBLE_PRECISION, manager, & 
                             & tag_M_SendingData, mcomm, status,ierr )

   call mpi_recv( inputmesh%pointmarkerlist, nnod, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
   
   allocate (inputmesh%segmentlist(2*inputmesh%numberofsegments), inputmesh%segmentmarkerlist(inputmesh%numberofsegments)) 
   
   call mpi_recv( inputmesh%segmentlist,    2*inputmesh%numberofsegments,MPI_INTEGER,manager,tag_M_SendingData,mcomm,status,ierr)
   call mpi_recv( inputmesh%segmentmarkerlist,inputmesh%numberofsegments,MPI_INTEGER,manager,tag_M_SendingData,mcomm,status,ierr)


   allocate (inputmesh%holelist (2*inputmesh%numberofholes) )
   if (inputmesh%numberofholes > 0 ) then
        call mpi_recv( inputmesh%holelist, 2*inputmesh%numberofholes,MPI_DOUBLE_PRECISION, & 
                       & manager,tag_M_SendingData,mcomm,status,ierr)
   endif

 
   allocate (inputmesh%regionlist(4*inputmesh%numberofregions) ) 
   call mpi_recv( inputmesh%regionlist, 4*inputmesh%numberofregions, MPI_DOUBLE_PRECISION, manager,  &
          & tag_M_SendingData, mcomm, status, ierr )
 
    call mpi_recv( cAnisotropy, 9, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr)
    call mpi_recv( nSigParams,  1, MPI_INTEGER,   manager, tag_M_SendingData, mcomm, status, ierr ) 
    
    allocate ( sigParams(nSigParams) ,iFreeParam(nSigParams))
    
    call mpi_recv( sigParams ,  nSigParams, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( iFreeParam , nSigParams, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( nFree,                1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
  
   call mpi_recv( maxnadapt,                            1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
   call mpi_recv(errortolerance, 1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
   call mpi_recv(idual_func ,    1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr ) 
   call mpi_recv(minQangle,      1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
 
    
    call mpi_recv( outputFileRoot,        256, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr)
    call mpi_recv( lSaveMeshFiles,          1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( linversion,              1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( lDisplayRefinementStats, 1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr ) 

    
    call mpi_recv( linearSolver,  32, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr)
    call mpi_recv( lprintDebug ,   1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( lUseBumpFields, 1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr )
    
    if (lPrintDebug) then
        lprintDebug_em2dkx = .true.
    else
        lprintDebug_em2dkx = .false.
    endif        
    
    call mpi_recv(minRangeProtector, 1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(minArea,           1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(ecutoff,           1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(hcutoff,           1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(max_nsubrefine,    1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(pct_refine,        1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
  
    call cpu_time(tEnd)     
    
    write(stime,'(a8,f9.4,a2)') ' Timer: ',tEnd - tStart, ' s'
    sFmt = '(a5,2x,i6,2x,       a6,2x,i6,2x,   69x, a32,     2x,a19)'     
    if (lprintMPItimers) write(*,sFmt)  'Proc:',myID,'Group:',iRefinementGrp,   'mpi_worker_receive_subset:',   trim(stime )
          
    end subroutine mpi_worker_receive_subset

!==================================================================================================================================! 
!=========================================================================================================== mpi_worker_send_results
!==================================================================================================================================!
    subroutine mpi_worker_send_results

    use mare2d_mpi_definitions    
    use em2dkx_mod 
    use mare2dem_global
    use mare2dem_input_data_params
    use mare2dem_output

    implicit none
    
    integer         ::  nrec, ierr, iRx, iTx, iFq

    real(8)         :: tStart, tEnd
    character(256)  :: sFmt
    character(32)   :: stime
    
    call cpu_time(tStart) 

    call mpi_comm_rank( mcomm, myID, ierr )
        
    call mpi_send( iRefinementGrp, 1, MPI_INTEGER,  manager, tag_W_SendingData, mcomm, ierr )  

    if (.not.lMT) then

        nrec =  nTxCSEM*nFreqCSEM

        do iRx = 1,nRxCSEM
                    
            if (lCompCSEM(1,iRx)) call mpi_send(ex(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompCSEM(2,iRx)) call mpi_send(ey(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompCSEM(3,iRx)) call mpi_send(ez(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompCSEM(4,iRx)) call mpi_send(hx(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompCSEM(5,iRx)) call mpi_send(hy(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompCSEM(6,iRx)) call mpi_send(hz(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)    
        
        enddo
        
        if (linversion) then
 
            do iFq = 1,nFreqCSEM
              do iTx = 1,nTxCSEM
                do iRx = 1,nRxCSEM
                    
                    if ( lDataMaskCSEM(iRx,iFq,iTx) ) then
                    
                    if ( lCompCSEM(1,iRx) ) then
                        call mpi_send(dex(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                    endif
                    if ( lCompCSEM(2,iRx) ) then
                        call mpi_send(dey(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                    endif
                    if ( lCompCSEM(3,iRx) ) then
                        call mpi_send(dez(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                    endif      
                    if ( lCompCSEM(4,iRx) ) then
                        call mpi_send(dhx(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                    endif
                    if ( lCompCSEM(5,iRx) ) then
                        call mpi_send(dhy(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                    endif
                    if ( lCompCSEM(6,iRx) ) then
                        call mpi_send(dhz(iRx,iTx,iFq)%dsig, nFree, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                    endif   
            
                    endif
                    
                enddo
              enddo
            enddo
        
        endif
 
    else ! MT:
    
        nrec =  nFreqMT
        
        do iRx = 1,nRxMT            
        
            if (lCompMT(1,iRx)) call mpi_send(ex_mt(iRx,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompMT(2,iRx)) call mpi_send(ey_mt(iRx,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompMT(3,iRx)) call mpi_send(ez_mt(iRx,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompMT(4,iRx)) call mpi_send(hx_mt(iRx,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompMT(5,iRx)) call mpi_send(hy_mt(iRx,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
            if (lCompMT(6,iRx)) call mpi_send(hz_mt(iRx,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)    
        
        enddo

        if (linversion) then
            
            nrec =  nFreqMT*nFree
            do iRx = 1,nRxMT 
             
                if (lCompMT(1,iRx)) call mpi_send(dex_mt_dsig(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                if (lCompMT(2,iRx)) call mpi_send(dey_mt_dsig(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                if (lCompMT(3,iRx)) call mpi_send(dez_mt_dsig(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                if (lCompMT(4,iRx)) call mpi_send(dhx_mt_dsig(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                if (lCompMT(5,iRx)) call mpi_send(dhy_mt_dsig(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr)
                if (lCompMT(6,iRx)) call mpi_send(dhz_mt_dsig(iRx,:,:), nrec, MPI_COMPLEX, manager, tag_W_SendingData, mcomm, ierr) 
            
            enddo
                      
        endif
            
    endif
   
    call cpu_time(tEnd)     
    
    write(stime,'(a8,f9.4,a2)') ' Timer: ',tEnd - tStart, ' s'
    sFmt = '(a5,2x,i6,2x,       a6,2x,i6,2x,   69x, a32,     2x,a19)'     
    if (lprintMPItimers) write(*,sFmt)  'Proc:',myID,'Group:',iRefinementGrp,   'mpi_worker_send_results:',   trim(stime )
     
      
    end subroutine mpi_worker_send_results
    
!==================================================================================================================================! 
!=============================================================================================================== mpi_shutDownWorkers
!==================================================================================================================================!  
    subroutine mpi_shutDownWorkers
    
    use mare2d_mpi_definitions
    use mare2dem_global
    
    implicit none
    
    integer :: ierr, iworker
    
    write(*,*) ' '
    write(*,*) ' Shutting down the worker processors:'
    write(*,*) ' '
    do iworker= 1,nworkers
        call mpi_send( iflag_TellWorkerToQuit, 1, MPI_INTEGER, iworker, tag_M_SendingData, mcomm, ierr )
    enddo    
    write(*,*) ' '  

    end subroutine mpi_shutDownWorkers    

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! 
! Wrapper routines for Occam matrix-matrix and linear solves using Scalapack:
!
! multiplyATA_scalapack()
! solveLinearSystem_scalapack()
!
!-----------------------------------------------------------------------------------------------------------------------------------
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
    subroutine multiplyATA_scalapack()
!----------------------------------------------------------------------------------------------------------------------------------- 
    
    use mare2d_mpi_definitions
    use mare2dem_global
    use occam 
    
    implicit none 

!
! Local variables:
!
    integer :: iWorker, ierr
          
    !
    ! Launch Scalapack Matrix Multiply subroutine on Workers
    !
    do iWorker= 1,nWorkers
        call mpi_send( iflag_TellWorkerToScaLapackMM, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm, ierr )
    enddo 
    
    !
    ! Launch Scalapack subroutine on manager as well
    !    
    call mpi_Worker_ScaLapackMM( )

 
    end subroutine multiplyATA_scalapack

!-----------------------------------------------------------------------------------------------------------------------------------     
    subroutine solveLinearSystem_scalapack(istat)
!----------------------------------------------------------------------------------------------------------------------------------- 
    
    use mare2d_mpi_definitions  
    use mare2dem_global   
    use occam  
    
    implicit none
    
    integer :: iWorker, istat, ierr
       
    !
    ! Launch Scalapack Matrix Solver subroutine on Workers
    !
    do iWorker= 1,nWorkers
        call mpi_send( iflag_TellWorkerToScaLapackMS, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm, ierr )
    enddo 

    !
    ! Launch Scalapack subroutine on manager as well
    !    
    call mpi_Worker_ScaLapackMS(istat)
        
          
    end subroutine solveLinearSystem_scalapack
    
#endif
