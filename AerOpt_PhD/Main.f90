program AerOpt
   
    ! ****Initializing Parameters and Implement Modules**** !
    use CreateSnapshots
    use GenerateMesh
    use Toolbox
    use Optimization
    use ReadData
    use InputData
    use CFD
    use FDGD
    
    implicit none

    print *, ''
    print *, '**************************************************************************'
    print *, '**                                                                      **'
    print *, '***********************  WELCOME TO THE AEROPT TOOL  *********************'     
    print *, '***********  AN AUTOMATED AERODYNAMIC OPTIMISATION SOFTWARE  *************'
    print *, '**                                                                      **'
    print *, '*********************   written by Dr. DAVID NAUMANN   *******************'
    print *, '*********************   supervised by Dr. BEN EVANS    *******************'
    print *, '**************************************************************************'
    print *, ''
    
    !! Technical Terms
    ! Geometry - A distinct Mesh/Shape
    ! Nest - Coordinates of all Control Points applied to one Geometry
    ! Snapshot - Can mean both the initial Geometry or initial Nest applied to construct the POD
    
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
    
    newdir = '2DEngInletSnapshots_'//IV%version//'_'//date(3:8)//'_'//time(1:4) 
    
    ! ****Read Input Data(Fine Mesh, Coarse Mesh, CP Coordinates, Influence Box/Rectangle (IB)**** !
    print *, 'Start Read Data'
    call SubReadData()
    ! Output: Boundf, Coord, Connecf, Coord_CP
    
    ! ****Create Folder Structure for PrePro & Solver Output**** !
    print *, 'Create Directories'
    call createDirectoriesInit()
    if (IV%SystemType == 'W' .and. IV%RunOnCluster == 'Y')   then    ! AerOpt is executed from a Windows machine connected to a Linux machine
        
        call createDirectoriesInit()
        call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')   ! Submits create directory file
     
    elseif (IV%SystemType == 'W') then ! AerOpt is executed from a Windows machine alone
        
        call createDirectoriesWindows()
        call system('FileCreateDir.bat')    ! Submits create directory file
        
    elseif (IV%SystemType == 'Q' .or. IV%SystemType == 'B') then     ! AerOpt is executed from a Linux machine
        
        call createDirectoriesInit()
        call system('chmod a+x ./FileCreateDir.scr')    
        call system('./FileCreateDir.scr')    ! Submits create directory file
            
    else       
        STOP 'INPUT ERROR: System Type selected does not exist! Program stopped.'       
    end if        
    
    ! *****Create Initial Nests for the Snapshots****** ! 
    print *, 'Start LHS Sampling - Create Initial Nests'
    call SubCreateSnapshots()    
    ! Output: Initial Nests - Sampling Points for Snapshots
  
! Same movement: Hard coded for testing Mesh Distortion for Aerofoil case - Control Nodes move equally
    if (IV%samemovement == .true.) then
        Snapshots(:,1) = IV%xmax
        Snapshots(:,2) = IV%xmax
        Snapshots(:,3) = IV%ymax
        Snapshots(:,4) = IV%ymax
    end if

    allocate(character(len=3) :: istr)
    write(istr, '(1f3.1)') IV%Ma
    open(29, file=newdir//'/Snapshots'//istr//'.txt', form='formatted',status='unknown')
    deallocate(istr)
    write(29, *) 'Snapshots'   
    write(29,'(<IV%NoSnap>f17.10)') Snapshots
    close(29)

!!!! IMPLEMENT double-check, wether Dimension of file and Input are compliant OR error check while Reading files

    ! ****Generate Full Fidelity Solutions of Snapshots**** !
    if (IV%MeshMovement == 1) then
        print *, 'Start PreMeshing'
        call PreMeshing()
        print *, 'End PreMeshing - All Area Coefficients calculated'
    end if
    print *, ''
    call SubCFD(1, IV%NoSnap, Snapshots, IV%NoSnap)
    call PostSolverCheck(IV%NoSnap, 0)
    
    ! ****Optimize Mesh by the help of Cuckoo Search and POD**** !
    call SubOptimization()
    ! Output: Optimized mesh via Modified Cuckoo Search and POD

    call timestamp()
    
    end program AerOpt