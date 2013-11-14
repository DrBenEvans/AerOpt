program AerOpt
   
    ! Initializing Parameters and Implement Modules
    use CreateInitialNests
    integer :: i, j, k, NoNests, NoIter, NoCtrPts, NoDim   ! Number of Nests, Iterations, Control Points, Dimensions
    real :: xmax, ymax, Ma, engFMF                  ! max Displacements, Mach number, ??
    character (len=4) filename
    character :: runOnCluster
    
    ! User Input
!!! TO-DO: create text files including the input data of meshmat & datamat
! Open and Read files
    !pathexec = '/eng/cvcluster/egevansbj/codes/prepro/2Dprepro_duct'
    !pathexec2 = '/eng/cvcluster/egmanon/Duct2d'
    
    Ma = 0.5            ! Mach Number
    xmax = 0.00         ! Maximum horizontal displacement of Control Nodes
    ymax = 0.02         ! Maximum vertical displacement of Control Nodes
    zmax = 0.05         ! Maximum lateral displacement of Control Nodes
    NoNests = 30        ! Number of Nests (Cuckoo Search)
    NoIter = -3         ! ?? Number of Iterations
    NoCtrPts = 7        ! Number of Control Points
    NoDim = 2           ! Number of Dimensions
    filename = 'case'
    runOnCluster = 'Y'
    engFMF = 1.0        ! What is it?
    
    !! PLOT initial Mesh: Actual and Coarse Mesh! --> MATLAB output file?
    
    !****Sub-Section: Create Initial Nests for the CFD Solver****** 
    !***********included in CreateInitialNests module**************
    call CreateRandomNumber
    print *, rn !Output: rn = Random Number
    
    call LHC(NoDim, NoCtrPts, xmax, ymax, zmax)                !Sampling of initial points/nests via LHC
    
    ! Output:
end program AerOpt