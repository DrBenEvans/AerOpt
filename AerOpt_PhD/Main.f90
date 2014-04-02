program AerOpt
   
    ! ****Initializing Parameters and Implement Modules**** !
    use CreateInitialNests
    use GenerateInitialMeshes
    use Toolbox
    use Optimization
    use ReadData
    use InputData
    use CFD
    
    implicit none
    integer :: ii
    
    ! ****User Input****** !
    call SubInputData(IV)
    
    ! Check xmax, ymax, zmax & NoDim Input
    ! 1D: only x considered, 2D: x & y considered, 3D: x, y & z considered
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
    newdir = '2DEngInletSnapshots_'//IV%version//'_'//date(3:8)//'_'//time(1:4) !'2DEngInletSnapshots_1.6_140320_1150'
    
    
    ! ****Read Input Data(Fine Mesh, Coarse Mesh, CP Coordinates, Influence Box/Rectangle (IB)**** !
    print *, 'Start Read Data'
    call SubReadData()
    ! Output: Boundf, Coord, Connecf, Coord_CP
    
    
    ! ****Sub-Section: Create Initial Nests for the CFD Solver****** ! 
    ! ***********included in CreateInitialNests module************** !
    print *, 'Start LHS Sampling - Create Initial Nests'
    call SubCreateInitialNests()                !Sampling of initial points/nests via LHC    
    ! Output: InitialNests - Sampling Points for initial Nests
    
    !allocate(ArrayTemp(1000,(IV%NoDim*IV%NoCP)))
    open(29,file='Output_Data/InitialNests1000.txt')
    write(29, *) 'Initial Nests'
    write(29,'(1000f13.10)') InitialNests
    close(29) 
    !InitialNests = ArrayTemp(1:999,:)
    !deallocate(ArrayTemp)
    
!!!!!! IMPLEMENT double-check, wether Dimension of file and Input are compliant OR error check while Reading files
    
    ! ****Generate initial Meshes/Snapshots**** !
    call IdentifyBoundaryFlags()
    ! Output: Boundary Matrix incluing flags of adiabatic viscous wall, far field & engine inlet (boundf)
    allocate(RD%coord_temp(RD%np,IV%nodim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main "    
    do ii = 1, IV%NoNests
        print *, "Generating Mesh", ii, "/", IV%NoNests
        RD%coord_temp = RD%coord
        call SubGenerateInitialMeshes(InitialNests(ii,:))
        ! Output: new coordinates - Mesh with moved boundaries based on Initial Nest
        
!!!!! IMPLEMENT Mesh Quality Test

        ! Determine correct String      
        call DetermineStrLen(istr, ii) 
        ! Write Snapshot to File
        call InitSnapshots()
        deallocate (istr)
    end do
    deallocate(RD%coord_temp)

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
    
        call PreProcessing(ii)
    
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
            
       call Solver(ii) 
                
    end do
    print *, 'Finished Submitting Jobs to FLITE 2D Solver'
    
    
    ! ****Wait & Check for FLITE Solver Output**** !
    print*, 'Start Sleep'
    call Sleep()
    print*, 'End Sleep - Jobs are finished'
    
    if (IV%SystemType == 'W') then
        call TransferSolutionOutput()
    end if
    
    
    ! ****Check Simulation Results**** !
    print*, 'Start Check for Convergence'
    ii = 0
    call CheckforConvergence(ii)
     print*, 'All Solutions converged'
     
    
    ! ****Optimize Mesh by the help of Cuckoo Search and POD**** !
    print *, 'Start Optmization'
    call SubOptimization()
    ! Output: Optimized mesh via Cuckoo Search and POD
    
    
    ! ****Generate Optimum Mesh and Safe in file**** !
    allocate(RD%coord_temp(RD%np,IV%nodim),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main " 
    RD%coord_temp = RD%coord
    call SubGenerateInitialMeshes(NestOpt)
    ! Output: Optimum Coordinates - 1 Mesh with moved boundaries based on optimum Control Point Coordinates
    
    ! Safe Optimum Geometry in Text File
    open(99, file= OutFolder//'/OptimumMesh300.txt')         
    write(99,'(1I8)') RD%np
    write(99,'(1I8)') RD%ne
    write(99,'(2f12.7)') transpose(RD%coord_temp)
11  format(3I8)
    close(99)
    
    end program AerOpt