program AerOpt
   
    ! Initializing Parameters and Implement Modules  
    use CreateInitialNests
    use GenerateInitialMeshes
    use Toolbox
    
    implicit none
    integer :: i, j, k                ! Simple Loop Variables

    ! User Input
! Open and Read files
    !pathexec = '/eng/cvcluster/egevansbj/codes/prepro/2Dprepro_duct'
    !pathexec2 = '/eng/cvcluster/egmanon/Duct2d'
    
    real, parameter :: Ma = 0.5            ! Mach Number
    real, parameter :: xmax = 0.01         ! Maximum horizontal displacement of Control Nodes
    real, parameter :: ymax = 0.02         ! Maximum vertical displacement of Control Nodes
    real, parameter :: zmax = 0.00         ! Maximum lateral displacement of Control Nodes
    real, parameter :: engFMF = 1.0        ! What is it?
    integer, parameter :: NoNests = 30     ! Number of Nests (Cuckoo Search)
    integer, parameter :: NoIter = -3      ! ?? Number of Iterations
    integer, parameter :: NoCP = 7         ! Number of Control Points
    integer, parameter :: NoDim = 2        ! Number of Dimensions
    ! 1D: only x considered, 2D: x & y considered, 3D: x, y & z considered
    
    character (len=4), parameter :: filename = 'case'       ! filename of Meshfiles 
    character, parameter :: runOnCluster = 'Y'
    character(len=:), allocatable :: ShowText ! =       ! 
    character(len=:), allocatable :: TitleCaption ! = 'Dimension Check'   ! 
    
    ! Check xmax, ymax, zmax & NoDim Input
    select case (NoDim)
    case (1)
        if (ymax /= 0 .or. zmax /= 0) then
            allocate(character(60) :: Showtext)
            ShowText = 'ymax and/or zmax remain unconsidered in 1 Dimension'
            allocate(character(15) :: TitleCaption)
            TitleCaption = 'Dimension Check'
            call MsgBox(ShowText, TitleCaption)
        endif
    case (2)
        if (zmax /= 0) then
            allocate(character(60) :: Showtext)
            ShowText = 'zmax remains unconsidered in 2 Dimensions'
            allocate(character(15) :: TitleCaption)
            TitleCaption = 'Dimension Check'
            call MsgBox(ShowText, TitleCaption)
        endif
    end select
        
    !!!!!! PLOT initial Mesh: Actual and Coarse Mesh! --> MATLAB output file?
    
    ! Automatically generates a random initial number based on time and date
    call RANDOM_SEED    
    
    ! ****Sub-Section: Create Initial Nests for the CFD Solver****** ! 
    ! ***********included in CreateInitialNests module************** !  
    call SubCreateInitialNests(NoNests, NoDim, NoCP, xmax, ymax, zmax)                !Sampling of initial points/nests via LHC    
    ! Output: InitialNests - Sampling Points for initial Nests
    
    ! ****Generate initial Meshes**** !
    call ReadData(NoCP, NoDim)
    ! Output: Boundf, Coordf, Connecf, Coord_CP
    !!!!!! Implement double-check, wether Dimension of file and Input are compliant OR error check while Reading files
    
    do i = 1, NoNests
        print *, "Generating Mesh", i, "/", NoNests
        call SubGenerateInitialMeshes(NoDim, NoCP, coordf, connecf, boundf, Coord_CP) ! Output: 
    end do
    
    

        
end program AerOpt