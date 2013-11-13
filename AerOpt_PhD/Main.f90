program AerOpt
   
    ! Initializing Parameters and Implement Modules
    use CreateInitialNests
    integer :: NoNests, NoIter
    real :: xmax, ymax, Ma, engFMF
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
    NoNests = 30        ! Number of Nests (Cuckoo Search)
    NoIter = -3         ! ?? Number of Iterations
    filename = 'case'
    runOnCluster = 'Y'
    engFMF = 1.0        ! Application?
    print *, Ma
    
    !! Plot initial Mesh: Actual and Coarse Mesh! --> MATLAB output file?
    call CreateRandomNumber
    
end program AerOpt