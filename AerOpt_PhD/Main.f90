program AerOpt
   
    ! ****Initializing Parameters and Implement Modules**** !
    use CreateInitialNests
    use GenerateInitialMeshes
    use Toolbox
    use Optimization
    use ReadData
    use InputData
    
    implicit none
    integer :: ii
    
    ! ****User Input****** !
    call SubInputData(IV)
    
    ! Check xmax, ymax, zmax & NoDim Input
    if (IV%NoDim == 1) then
        if (IV%ymax /= 0 .or. IV%zmax /= 0) then
            print *, 'ymax and/or zmax remain unconsidered in 1 Dimension'
            print *, 'Input any and click enter to continue'
            read(*, *)
        end if
    elseif (IV%NoDim == 2) then
        if (IV%zmax /= 0) then
            print *, 'zmax remains unconsidered in 2 Dimensions'
            print *, 'Input any and click enter to continue'
            read(*, *)
        end if
    end if
    
    ! Automatically generates a random initial number based on time and date
    call RANDOM_SEED
    
    ! Get Time and Date for File and Folder Name creation
    call DATE_AND_TIME(date, time)
    newdir = '2DEngInletSnapshots_'//IV%version//'_'//date(3:8)//'_'//time(1:4)
    
    
    ! ****Read Input Data(Fine Mesh, Coarse Mesh, CP Coordinates, Influence Box/Rectangle (IB)**** !
    print *, 'Start Read Data'
    call SubReadData()
    ! Output: Boundf, Coord, Connecf, Coord_CP
    
    
    ! ****Sub-Section: Create Initial Nests for the CFD Solver****** ! 
    ! ***********included in CreateInitialNests module************** !
    print *, 'Start LHS Sampling - Create Initial Nests'
    call SubCreateInitialNests()                !Sampling of initial points/nests via LHC    
    ! Output: InitialNests - Sampling Points for initial Nests
    
    
!!!!!! IMPLEMENT double-check, wether Dimension of file and Input are compliant OR error check while Reading files
    
    ! ****Generate initial Meshes/Snapshots**** !
    allocate(coord_temp(np,IV%NoDim))
    allocate(boundff(nbf,(IV%NoDim+1)))
    boundff(:,1:2) = boundf
    do ii = 1, IV%NoNests
        print *, "Generating Mesh", ii, "/", IV%NoNests
        coord_temp = coord
        call SubGenerateInitialMeshes(coord_temp, connecf, boundf, coarse, connecc, Coord_CP,Rect, InitialNests(ii,:))
        ! Output: New Coordinates - 30 Snapshots with moved boundaries based on initial nests
        
        call IdentifyBoundaryFlags()
        ! Output: Boundary Matrix incluing flags of adiabatic viscous wall, far field & engine inlet (boundff)
        
!!!!! IMPLEMENT Mesh Quality Test

        ! Determine correct String      
        call DetermineStrLen(istr, ii) 
        ! Write Snapshot to File
        call InitSnapshots(coord_temp, boundff)
        deallocate (istr)
    end do
    

    ! ****Create Folder Structure for PrePro & Solver Output**** !
    print *, 'Create Directories'
    call createDirectoriesInit()
    if (IV%SystemType == 'W')   then    ! AerOpt is executed from a Windows machine
                             
        call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')   ! Submits create directory file
            
    else                        ! AerOpt is executed from a Linux machine
                            
        call system('chmod a+x ./FileCreateDir.scr')    
        call system('./FileCreateDir.scr')    ! Submits create directory file
            
    end if
    

    ! ****Call 2D Preprocessor and pass on input parameters**** !
    print *, 'Start Preprocessing'
    if (IV%SystemType == 'Q') then
        allocate(character(len=61) :: pathLin_Prepro)
        pathLin_Prepro = '/eng/cvcluster/egnaumann/2DEngInlSim/PrePro/2DPreProcessorLin'
    else
        allocate(character(len=56) :: pathLin_Prepro)
        pathLin_Prepro = '/home/david.naumann/2DEngInlSim/PrePro/2DPreProcessorLin'
    end if
    pathWin = 'Flite2D\PreProcessing'   
    do ii = 1, IV%NoNests
    
        ! Determine correct String      
        call DetermineStrLen(istr, ii)
        ! write Inputfile
        call PreProInpFile()
        
        if (IV%SystemType == 'W') then
             
            allocate(character(len=29) :: strSystem)
            strSystem = pathWin
            
        else
            
            ! write command (for Linux)
            IntSystem = 10 + len(trim(IV%filename)) + len(istr) + len(pathLin_Prepro)
            allocate(character(len=IntSystem) :: strSystem)
            strSystem = pathLin_Prepro//' > '//trim(IV%filename)//istr//'.outpre'
            
        end if
        print *, 'Preprocessing Snapshot', ii
        print *, ' '
        call system(strSystem)   ! System operating command called to activate fortran       
        deallocate (istr)
        deallocate (strSystem)
    
    end do
    print *, 'Finished Preprocessing'
    
    
    ! ****Call 2D FLITE Solver and pass on input parameters**** !
    print *, 'Call FLITE 2D Solver'
    if (IV%SystemType /= 'B') then
        allocate(character(len=55) :: pathLin_Solver)
        pathLin_Solver = '/eng/cvcluster/egnaumann/2DEngInlSim/Solver/2DSolverLin'
    else
        allocate(character(len=50) :: pathLin_Solver)
        pathLin_Solver = '/home/david.naumann/2DEngInlSim/Solver/2DSolverLin'
    end if
    
    do ii = 1, IV%NoNests
            
        ! Determine correct String      
        call DetermineStrLen(istr, ii)                      
        ! Creates the input file including Solver Parameters and a second file including I/O filenames
        call WriteSolverInpFile()
        ! writes the batchfile to execute Solver on Cluster
        call writeBatchFile()
    
        ! Is AerOpt executed from Linux or Windows?                
        if (IV%SystemType == 'W')   then    ! AerOpt is executed from a Windows machine
            
            ! Transfer Files from Windows Machine onto Cluster
            call transferFilesWin()            
            call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')
            
            if (IV%runOnCluster == 'Y') then
                call Triggerfile()           ! Triggerfile for submission
            else
                call TriggerFile2()  ! Triggerfile for submission
            end if
                    
            ! Submits Batchfile via Putty
            call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Trigger.sh', 'putty')
                
        else    ! AerOpt is executed from a Linux machine
            
            ! Transfer Files in correct folder on Cluster
            call transferFilesLin()
            call system('chmod a+x ./FileCreateDir.scr')
            call system('./FileCreateDir.scr')
            
            if (IV%runOnCluster == 'Y') then
                call Triggerfile()     ! Triggerfile for submission
            else
                call TriggerFile2()               ! Triggerfile for submission
            end if
            
            ! Submits Batchfile
            call system('chmod a+x ./Trigger.sh')
            call system('./Trigger.sh')
                    
        end if
                
        deallocate(istr) 
                
    end do
    print *, 'Finished Submitting Jobs to FLITE 2D Solver'
    
    ! ****Wait & Check for FLITE Solver Output**** !
    print*, 'Start Sleep'
    jobcheck = 0
    waitTime = 0
    j = 1
    do while (jobcheck==0)
        
        ! Wait Function
        print*, 'Sleep'
        call SleepQQ(IV%delay*1000)
        print*, 'Wake Up - Check ', j
        j = j + 1
        
        ! Check Status of Simulation by checking the existence of all error files
        do ii = IV%NoNests, 1, -1
            
            ! Determine correct String      
            call DetermineStrLen(istr, ii)
            ! Creates File containing Linux commands to check for last file
            call CheckSimStatus()
            ! Submit File
            if (IV%SystemType == 'W')   then
                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'CheckStatus.scr', 'plink')
            else
                call system('chmod a+x ./CheckStatus.scr')
                call system('./CheckStatus.scr')
            end if
            ! Creates File to transfer response from Windows to Linux
            if (IV%SystemType == 'W')   then
                call CheckSimStatus2()
                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'CheckStatus.scr', 'psftp')
            end if
            
            open(1, file='check.txt')
            read(1,*) jobcheck
            close(1)
            deallocate(istr)
            if (jobcheck == 0) EXIT           
            
        end do
        
        waitTime = (IV%delay/3600.0) + waitTime
        if (waitTime > IV%waitMax) then
            STOP 'Cluster Simulation Time exceeded maximum waiting Time'
        end if
        
    end do
    print*, 'End Sleep - Jobs are finished'
    
!!!! Check amount of .rsd files & conversion!
    
    ! ****Optimize Mesh by the help of Cuckoo Search and POD**** !
    print *, 'Start Optmization'
    call SubOptimization(cond, MxDisp_Move, np)
    ! Output: Optimized mesh via Cuckoo Search and POD
    
    
    ! ****Generate Optimum Mesh and Safe in file**** !
    coord_temp = coord
    call SubGenerateInitialMeshes(coord_temp, connecf, boundf, coarse, connecc, Coord_CP, Rect, NestOpt)
    ! Output: Optimum Coordinates - 1 Mesh with moved boundaries based on optimum Control Point Coordinates
    
    ! Safe Optimum Geometry in Text File
    open(99, file= OutFolder//'/OptimumMesh.txt')         
    write(99,'(1I8)') np
    write(99,'(1I8)') ne
    write(99,'(2f12.7)') transpose(coord_temp)
11  format(3I8)
    close(99)
    
    end program AerOpt