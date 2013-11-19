program AerOpt
   
    ! Initializing Parameters and Implement Modules
    use CreateInitialNests
    use Toolbox
    integer :: i, j, k                ! Simple Loop Variables   
    
    ! User Input
!!! TO-DO: create text files including the input data of meshmat & datamat
! Open and Read files
    !pathexec = '/eng/cvcluster/egevansbj/codes/prepro/2Dprepro_duct'
    !pathexec2 = '/eng/cvcluster/egmanon/Duct2d'
    
    real, parameter :: Ma = 0.5            ! Mach Number
    real, parameter :: xmax = 0.02         ! Maximum horizontal displacement of Control Nodes
    real, parameter :: ymax = 0.02         ! Maximum vertical displacement of Control Nodes
    real, parameter :: zmax = 0.05         ! Maximum lateral displacement of Control Nodes
    real, parameter :: engFMF = 1.0        ! What is it?
    integer, parameter :: NoNests = 30     ! Number of Nests (Cuckoo Search)
    integer, parameter :: NoIter = -3      ! ?? Number of Iterations
    integer, parameter :: NoCP = 7         ! Number of Control Points
    integer, parameter :: NoDim = 2        ! Number of Dimensions
    character (len=4), parameter :: filename = 'case'   ! filename of Meshfiles 
    character, parameter :: runOnCluster = 'Y'
    integer, dimension(NoNests) :: rp
    integer, dimension(NoCP) :: ones
    
    !! PLOT initial Mesh: Actual and Coarse Mesh! --> MATLAB output file?
    
    !****Sub-Section: Create Initial Nests for the CFD Solver****** 
    !***********included in CreateInitialNests module**************
    
    call RANDOM_SEED    ! Automatically generates a random initial number based on time and date
    
    !call SubCreateInitialNests(NoNests, NoDim, NoCP, xmax, ymax, zmax)                !Sampling of initial points/nests via LHC    
    ! Output:
    call randperm(NoNests, rp)
    write (*,2) rp
    ones = (/ (0, i=1,NoCP) /)
    write (*,2) ones
2   format(1i3)
        
end program AerOpt