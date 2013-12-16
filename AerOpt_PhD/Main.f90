program AerOpt
   
    ! Initializing Parameters and Implement Modules  
    use CreateInitialNests
    use GenerateInitialMeshes
    use Toolbox
    use Optimization
    
    implicit none
    integer :: i, j, k                ! Simple Loop Variables

    ! User Input
    !pathexec = '/eng/cvcluster/egevansbj/codes/prepro/2Dprepro_duct'
    !pathexec2 = '/eng/cvcluster/egmanon/Duct2d'
    
    !real, parameter :: Mach = 0.5            ! Mach Number
    real, parameter :: xmax = 0.00         ! Maximum horizontal displacement of Control Nodes
    real, parameter :: ymax = 0.02         ! Maximum vertical displacement of Control Nodes
    real, parameter :: hMa = 0.5         ! Maximum vertical displacement of Control Nodes
    real, parameter :: zmax = 0.00         ! Maximum lateral displacement of Control Nodes
    real, parameter :: engFMF = 1.0        ! What is it?
    integer, parameter :: NoNests = 30     ! Number of Nests (Cuckoo Search)
    integer, parameter :: NoIter = -3      ! ?? Number of Iterations
    integer, parameter :: NoCP = 7         ! Number of Control Points
    integer, parameter :: NoDim = 2        ! Number of Dimensions
    integer, parameter :: NoG = 5            ! Number of Generations
    ! 1D: only x considered, 2D: x & y considered, 3D: x, y & z considered
     
    character, parameter :: runOnCluster = 'Y'
    character(len=2) :: istr
    character(len=1) :: ShowText
    
    ! Check xmax, ymax, zmax & NoDim Input
    select case (NoDim)
    case (1)
        if (ymax /= 0 .or. zmax /= 0) then
            print *, 'ymax and/or zmax remain unconsidered in 1 Dimension'
            print *, 'Input any and click enter to continue'
            read(*, *) Showtext
        endif
    case (2)
        if (zmax /= 0) then
            print *, 'zmax remains unconsidered in 2 Dimensions'
            print *, 'Input any and click enter to continue'
            read(*, *) Showtext
        endif
    end select
    
    ! Automatically generates a random initial number based on time and date
    call RANDOM_SEED    
    
    ! ****Sub-Section: Create Initial Nests for the CFD Solver****** ! 
    ! ***********included in CreateInitialNests module************** !  
    call SubCreateInitialNests(NoNests, NoDim, NoCP, xmax, ymax, zmax)                !Sampling of initial points/nests via LHC    
    ! Output: InitialNests - Sampling Points for initial Nests
    
    ! ****Read Input Data(Fine Mesh, Coarse Mesh, CP Coordinates, Influence Box/Rectangle (IB)**** !
    print *, 'Start Read Data'
    call ReadData(NoCP, NoDim)
    ! Output: Boundf, Coord, Connecf, Coord_CP
    
    !!!!!! IMPLEMENT double-check, wether Dimension of file and Input are compliant OR error check while Reading files
    
    ! ****Generate initial Meshes/Snapshots**** !
    allocate(coord_temp(np,NoDim))
    do i = 1, NoNests
        print *, "Generating Mesh", i, "/", NoNests
        coord_temp = coord
        call SubGenerateInitialMeshes(NoDim, NoCP, coord_temp, connecf, boundf, coarse, connecc, Coord_CP,Rect, InitialNests(i,:))
        ! Output: New Coordinates - 30 Snapshots with moved boundaries based on initial nests
         
        ! Safe Snapshot in Text File
        write( istr, '(I2)' )  i
        open(99, file='Output_Data/Snapshot'//istr//'.txt')         
        write(99,10) transpose(coord_temp)
10      format(2f12.7)        
        close(99)
    end do
    
    !!!!! IMPLEMENT Mesh Quality Test
    
    ! ****Optimize Mesh by the help of Cuckoo Search and POD**** !
    call SubOptimization(NoNests, NoCP, NoDim, cond, InitialNests, MxDisp_Move, np, xmax, hMa)
    ! Output: Optimized mesh via Cuckoo Search and POD
    
    
end program AerOpt