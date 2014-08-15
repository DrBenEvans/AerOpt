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
    
    newdir = '2DEngInletSnapshots_1.6_140319_1945' !'2DEngInletSnapshots_'//IV%version//'_'//date(3:8)//'_'//time(1:4)
    
    ! ****Read Input Data(Fine Mesh, Coarse Mesh, CP Coordinates, Influence Box/Rectangle (IB)**** !
    print *, 'Start Read Data'
    call SubReadData()
    ! Output: Boundf, Coord, Connecf, Coord_CP
    
    
    ! ****Create Folder Structure for PrePro & Solver Output**** !
    print *, 'Create Directories'
    call createDirectoriesInit()
    if (IV%SystemType == 'W')   then    ! AerOpt is executed from a Windows machine
                             
        call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')   ! Submits create directory file
            
    else                        ! AerOpt is executed from a Linux machine
                            
        call system('chmod a+x ./FileCreateDir.scr')    
        call system('./FileCreateDir.scr')    ! Submits create directory file
            
    end if


    ! **** Create Initial Nests for the Snapshots****** ! 
    print *, 'Start LHS Sampling - Create Initial Nests'
    call SubCreateSnapshots()    
    ! Output: Initial Nests - Sampling Points for Snapshots
    
    allocate(character(len=3) :: istr)
    write(istr, '(1f3.1)') IV%Ma
    open(29, file=newdir//'/Snapshots'//istr//'.txt', form='formatted',status='new')
    deallocate(istr)
    write(29, *) 'Snapshots'
    write(29,'(<IV%NoSnap>f13.10)') Snapshots
    close(29)
    
!!!! IMPLEMENT double-check, wether Dimension of file and Input are compliant OR error check while Reading files

    call IdentifyBoundaryFlags()
    ! Output: Boundary Matrix incluing flags of adiabatic viscous wall, far field & engine inlet (boundf)

    ! **** Generate Full Fidelity Solutions of Snapshots**** !
    call SubCFD(1, IV%NoSnap, Snapshots, IV%NoSnap)
    call PostSolverCheck(IV%NoSnap, 0)
    
    ! ****Optimize Mesh by the help of Cuckoo Search and POD**** !
    call SubOptimization()
    ! Output: Optimized mesh via Modified Cuckoo Search and POD
    
    
    ! ****Generate Optimum Mesh and Safe in file**** !
    allocate(RD%coord_temp(RD%np,IV%nodim),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main " 
    RD%coord_temp = RD%coord
    call SubGenerateMesh(NestOpt)
    ! Output: Optimum Coordinates - 1 Mesh with moved boundaries based on optimum Control Point Coordinates
    
    ! Safe Optimum Geometry in Text File
    call InitSnapshots(005)
    call timestamp()
    
    end program AerOpt