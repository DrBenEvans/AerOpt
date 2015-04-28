      
module InputData
    
    type InputVariablesData
    
        real :: Ma                      ! Mach number
        real, dimension(100) :: xrange      ! Range of horizontal displacement of Control Nodes
        real, dimension(100) :: yrange      ! Range of vertical displacement of Control Nodes 
		real, dimension(100) :: zrange      ! Range of lateral displacement of Control Nodes
        real, dimension(100) :: angle       ! Range of change in angle of Control Nodes
        real, dimension(50) :: CNconnecttrans    ! allows to connect Control Node movements, e.g. CNconnect(5) = 3 - CN 5 is connected to movement of CN 3
        real, dimension(50) :: CNconnectangle    ! allows to connect Control Node movements, e.g. CNconnect(5) = 3 - CN 5 is connected to movement of CN 3       
        real :: gamma                   ! Ratio of specific heats
        real :: R                       ! specific gas constant
        real :: Re                      ! Reynoldsnumber
        real :: AlphaInflowDirection    ! Angle/Direction of Flow within the Solver(0: Left to Right, 180: Right to Left)
        real :: Tamb					! ambient Temperature [K]
        real :: Pamb    				! ambient Pressure [Pa]
        real :: engFMF                  ! Solver variable - engines Front Mass Flow
        real :: Top2Low                 ! Fraction of Top to Low Nests
        integer :: NoNests              ! Number of Nests (Cuckoo Search)
        integer :: NoSnap               ! Number of initial Snapshots
        integer :: NoCN                 ! Number of Control Nodes
        integer :: NoDim                ! Number of Dimensions
        integer :: DoF                  ! Degrees of Freedom in the System        
        integer :: NoG                  ! Number of Generations
        integer :: NoPOMod              ! No of POD Modes considered
        integer :: NoLeviSteps          ! Number of Levy walks per movement
        integer :: NoIter               ! Batch File variable - Number of Iterations
        integer :: turbulencemodel      ! Turbulence Model applied in Solver
        logical :: constrain            ! Constrain: Include boundaries of design space for Levy Walk - 1:Yes 0:no
        logical :: AdaptSamp            ! Adaptive Sampling - T: Active
        integer :: delay                ! Delay per check in seconds
        integer :: waitMax              ! maximum waiting time in hours
        integer :: maxit                ! maximum Iterations for Solver File before stopping convergence try
        integer :: NoDelBP              ! Number of Points placed on boundars for Delaunay Triangulation
        real :: Aconst                  ! Levy Flight parameter (determined emperically)   
        character(len=20) :: filename   ! I/O file of initial Meshes for FLITE solver
        character :: runOnCluster       ! Run On Cluster or Run on Engine?
        character :: SystemType         ! Windows('W'), Cluster/QSUB ('Q') or HPCWales/BSUB ('B') System? (Cluster, HPCWales = Linux, Visual Studio = Windows)    
        character(len=20) :: UserName   ! Putty Username - Cluster: egnaumann
        character(len=20) :: Password   ! Putty Password
        character(len=3) :: version
        integer :: MeshMovement
        integer :: ObjectiveFunction
        logical :: Meshtest
        logical :: Pol                  ! POD using Polynomial
        logical :: multiquadric         ! RBF type for POD
        logical :: POD
		logical :: meanP
    
    end type InputVariablesData
    
    type(InputVariablesData) :: IV
    integer :: allocatestatus                                   ! Check Allocation Status for Large Arrays
    character(len=10) :: DataFolder = 'Input_Data'                ! Input Folder Name
    character(len=9) :: InFolder = 'Mesh_Data'                ! Input Folder Name
    character(len=11) :: OutFolder = 'Solver_Data'              ! Output Folder Name
    integer :: IntSystem                                        ! Length of System command string; used for character variable allocation
    integer :: maxDoF                                           ! maximum Degrees of Freedom available
    integer :: DoFtransmotion                                   ! Degrees of Freedom in the System considering translatoric movements
    integer :: av                                               ! Allocater Variable
    real :: waitTime                                            ! waiting time for Simulation Results
    integer :: jobcheck                                         ! Check Variable for Simulation 
    character(len=:), allocatable :: istr                       ! Number of I/O file
    character(len=29) :: pathWinPrePro                          ! Path to Windows preprocessor file
    character(len=24) :: pathWinSolver                          ! Path to Windows solver file
    character(len=:), allocatable :: pathLin_Prepro             ! Path to Linux preprocessor file
    character(len=:), allocatable :: pathLin_Solver             ! Path to Linux Solver file
    character(len=:), allocatable :: strSystem                  ! System Command string for communication with FLITE Solver
    character(len=8) :: date                                    ! Container for current date
    character(len=10) :: time                                   ! Container for current time
    character(len=35) :: newdir                                 ! Name of new folder for a new Solution generated by 2D solver
    double precision :: alpha, beta                             ! Matrix Multiplicaion LAPACK variables
    
contains
    
    subroutine SubInputData(IV)
    
        ! Variables
        implicit none
        type(InputVariablesData) :: IV
    
        ! Body of SubInputData
        namelist /InputVariables/ IV
    
        IV%Ma = 0.5  		            ! Mach number
        IV%Tamb = 30					! ambient Temperature [deg]
        IV%Pamb = 101325				! ambient Pressure [Pa]
        IV%R = 287                  	! specific gas constant
        IV%gamma = 1.4                  ! Ratio of specific heats
        IV%Re = 0.0						! Reynoldsnumber
        IV%xrange = 0.00                ! Maximum horizontal displacement of Control Nodes    
        IV%yrange = 0.00                ! Maximum vertical displacement of Control Nodes    
        IV%zrange = 0.00                ! Maximum lateral displacement of Control Nodes
        IV%angle = 0.00                 ! For Aerofoil Angle of Attack Test
        IV%CNconnecttrans = 0           ! allows to connect Control Node movements, e.g. CNconnect(5) = 3 - CN 5 is connected to movement of CN 3
        IV%CNconnectangle = 0           ! allows to connect Control Node movements, e.g. CNconnect(5) = 3 - CN 5 is connected to movement of CN 3
        IV%engFMF = 1.0			        ! engines Front Mass Flow(Solver variable)
        IV%AlphaInflowDirection = 0.0   ! Angle/Direction of Flow within the Solver(0: Left to Right, 180: Right to Left)
        IV%turbulencemodel = 0          ! Turbulence Model applied in Solver
        IV%Top2Low = 0.75		        ! Fraction of Top to Low Cuckoo Nests
        IV%NoSnap = 100                 ! Number of initial Snapshots
        IV%NoCN = 7			            ! Number of Control Nodes
        IV%NoDim = 2			        ! Number of Dimensions in Space 
        IV%DoF = 7                      ! Degrees of freedom in the system 
        IV%NoG = 100		            ! Number of Generations
        IV%NoPOMod = -1			        ! No of POD Modes considered 
        IV%NoLeviSteps = 100         	! Number of Levy walks per movement 
        IV%NoIter = -3               	! Batch File variable - Number of Iterations 
        IV%constrain = .TRUE.         	! Constrain: Include boundaries of design space for Levy Walk - 1:Yes 0:no
        IV%delay = 300               	! Sleep Time between check for Simulation Results in seconds
        IV%waitMax = 48			        ! maximum waiting time in hours
        IV%maxit = 30000                ! maximum Number of Iterations to reach convergence in Solver
        IV%Aconst = 0.01		        ! Levy Flight parameter (determined emperically)
        IV%filename = 'Snapshot'        ! I/O file of initial Meshes for FLITE solver
        IV%runOnCluster = 'Y'           ! Run On Cluster or Run on Engine?
        IV%SystemType = 'Q'             ! Windows('W'), Cluster/QSUB ('Q') or HPCWales/BSUB ('B') System? (Cluster, HPCWales = Linux, Visual Studio = Windows)
        IV%UserName = 'egnaumann'       ! Putty Username - Cluster: egnaumann
        IV%Password = 'Fleur666'        ! Putty Password
        IV%version = '1.8'
        IV%NoDelBP = 8                  ! Number of Delaunay Boundary Points       
        IV%ObjectiveFunction = 1        ! What is the optimisation target? 1 - Lift/Drag, 2 - Distortion, 3 - zero Lift Optimisation
        
        ! POD
        IV%POD = .false.                ! Activation of POD - TRUE is ACTIVE
        IV%Pol = .false.                ! Application of Polynomial?
        IV%multiquadric = .true.        ! using multiquadratic RBF function for POD 
        IV%AdaptSamp = .FALSE.          ! Adaptive Sampling - T: Active
        IV%meanP = .true.
        
        ! For Mesh Deformation
        IV%MeshMovement = 1
        IV%Meshtest = .true.
        
		! Read User Input Parameters
        open(1,file = DataFolder//'/AerOpt_InputParameters.txt',form='formatted',status='old')
        read(1,InputVariables)
        close(1)
        
        ! Derive maximum Degrees of Freedom
        maxDoF = 4*IV%NoCN			! 3 dimensions & angular motion per node
        
        ! Number of Nests per Generation 
        IV%NoNests = 10*IV%DoF     
        if (IV%NoNests > IV%NoSnap .or. IV%POD == .false.) then
            IV%NoNests = IV%NoSnap
        end if
        
        ! Path for PreProcessor
        if (IV%SystemType == 'Q') then
            allocate(character(len=61) :: pathLin_Prepro)
            pathLin_Prepro = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/PrePro/2DPreProcessorLin'
        else
            allocate(character(len=56) :: pathLin_Prepro)
            pathLin_Prepro = '/home/'//trim(IV%UserName)//'/AerOpt/PrePro/2DPreProcessorLin'
        end if
        pathWinPrePro = 'Executables\PreProcessing.exe' 
   
        ! Path for Solver
        if (IV%SystemType /= 'B') then
            allocate(character(len=55) :: pathLin_Solver)
            pathLin_Solver = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/Solver/2DSolverLin'
        else
            allocate(character(len=50) :: pathLin_Solver)
            pathLin_Solver = '/home/'//trim(IV%UserName)//'/AerOpt/Solver/2DSolverLin'
        end if
        pathWinSolver = 'Executables\2DSolver.exe'
        
    end subroutine SubInputData
    
end module InputData