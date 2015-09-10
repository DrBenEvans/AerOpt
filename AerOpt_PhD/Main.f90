program AerOpt
   
    ! ****Initializing Parameters and Implement Modules**** !
    use CreateSnapshots
    use GenerateMesh
    use Toolbox
    use Optimization
    use ReadData
    use InputData
    use CFD
    
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
    ! 2D: x & y and  considered, 3D: x, y & z considered
    if (IV%NoDim == 2) then
        print *, 'zrange will be set to zero in 2 Dimensions'
        IV%zrange = 0.0
    end if
    
    ! Automatically generates a random initial number based on time and date
    call RANDOM_SEED
    
    ! Get Time and Date for File and Folder Name creation
    call DATE_AND_TIME(date, time)
    allocate(character(len=1) :: istr)
    write( istr, '(I1)' )  IV%NoDim
    
    newdir = 'AerOpt'//istr//'D_'//IV%version//'_'//date(3:8)//'_'//time(1:4) 
    deallocate(istr)
    
    ! ****Read Input Data(Fine Mesh, Coarse Mesh, CP Coordinates, Influence Box/Rectangle (IB)**** !
    print *, 'Start Read Data'
    if (IV%NoDim == 2) then
        call SubReadData()
    elseif (IV%NoDim == 3) then
        call SubReadData_3D()
    end if    
    ! Output: Boundf, Coord, Connec, Coord_CN
    
    ! ****Create Folder Structure for PrePro & Solver Output**** !
    print *, 'Create Directories'
    call CreateFolderStructure()      
    
    ! *****Create Initial Nests for the Snapshots****** ! 
    print *, 'Start LHS Sampling - Create Initial Nests'
    call SubCreateSnapshots()    
    ! Output: Initial Nests - Sampling Points for Snapshots

    allocate(character(len=3) :: istr)
    write(istr, '(1f3.1)') IV%Ma
    open(29, file=newdir//'/Snapshots'//istr//'.txt', form='formatted',status='unknown')
    deallocate(istr)
    write(29, *) 'Snapshots'   
    write(29,'(<IV%NoSnap>f17.10)') Snapshots
    close(29)

!!!! IMPLEMENT double-check, wether Dimension of file and Input are compliant OR error check while Reading files

    ! ****Generate Full Fidelity Solutions of Snapshots**** !
    print *, 'Start PreMeshing'
    if (IV%NoDim == 2) then
        call PreMeshing()
    elseif (IV%NoDim == 3) then
        call PreMeshing_3D()
    end if
    print *, 'End PreMeshing - All Area Coefficients calculated'
    print *, ''
    call SubCFD(1, IV%NoSnap, Snapshots, IV%NoSnap)
    call PostSolverCheck(IV%NoSnap, 0)
    
    ! ****Optimize Mesh by the help of Cuckoo Search and POD**** !
    call SubOptimization()
    ! Output: Optimized mesh via Modified Cuckoo Search and POD

    call timestamp()
    
    end program AerOpt