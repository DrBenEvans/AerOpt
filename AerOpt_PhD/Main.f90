program AerOpt
   
    ! Initializing Parameters and Implement Modules  
    use CreateInitialNests
    use GenerateInitialMeshes
    use Toolbox
    use Optimization
    use ReadData
    
    implicit none
    integer :: i, j, k                ! Simple Loop Variables

    ! User Input
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
    integer, parameter :: NoPOMod = -1      ! No of POD Modes considered
    integer, parameter :: NoLeviSteps = 100 ! Number of Levy walks per movement
    logical, parameter :: constrain = 1
    real :: Aconst = 0.01                        ! Levy Flight parameter (determined emperical) 
    real, parameter :: p = 0.75             ! Fraction of Top to Low Nests
    ! 1D: only x considered, 2D: x & y considered, 3D: x, y & z considered
     
    character, parameter :: runOnCluster = 'Y'
    character, parameter :: IsLin = 'N'             ! Windows('N') or Linux('Y') System to run on? (Cluster = Linux, Visual Studio = Windows)
    character(len=:), allocatable :: istr, strPrepro
    character(len=1) :: ShowText
    character(len=51) :: pathLin
    character(len=8), parameter :: filename = 'Snapshot'
    character(len=7) :: pathWin
    integer :: PreInt
    
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
    print *, 'Start LHS Sampling - Create Initial Nests'
    call SubCreateInitialNests(NoNests, NoDim, NoCP, xmax, ymax, zmax)                !Sampling of initial points/nests via LHC    
    ! Output: InitialNests - Sampling Points for initial Nests
    
    ! ****Read Input Data(Fine Mesh, Coarse Mesh, CP Coordinates, Influence Box/Rectangle (IB)**** !
    print *, 'Start Read Data'
    call SubReadData(NoCP, NoDim)
    ! Output: Boundf, Coord, Connecf, Coord_CP
    
    !!!!!! IMPLEMENT double-check, wether Dimension of file and Input are compliant OR error check while Reading files
    
    ! ****Generate initial Meshes/Snapshots**** !
    allocate(coord_temp(np,NoDim))
    allocate(boundff(nbf,(NoDim+1)))
    boundff(:,1:2) = boundf
    do i = 1, NoNests
        print *, "Generating Mesh", i, "/", NoNests
        coord_temp = coord
        call SubGenerateInitialMeshes(NoDim, NoCP, coord_temp, connecf, boundf, coarse, connecc, Coord_CP,Rect, InitialNests(i,:))
        ! Output: New Coordinates - 30 Snapshots with moved boundaries based on initial nests
        
        call IdentifyBoundaryFlags()
        ! Output: Boundary Matrix incluing flags of adiabatic viscous wall, far field & engine inlet (boundff)
        
        ! Safe Snapshot in Data File
        if (i < 10) then
            allocate(character(len=1) :: istr)
            write( istr, '(I1)' )  i
        else
            allocate(character(len=2) :: istr)
            write( istr, '(I2)' )  i
        end if
        open(99, file='Output_Data/Snapshot'//istr//'.dat')
        write(99,*) 1
        write(99,*) 'David Naumann'
        write(99,*) 'NoTrgElem NoNodes NoBound'        
        write(99,*) ne, np, nbf
        write(99,*) 'Connectivities'
        do j = 1, ne
            write(99,*) j, connecf(j,:)
        end do
        write(99,*) 'Coordinates'
        do j = 1, np
            write(99,*) j, coord_temp(j,:)
        end do
        write(99,*) 'Boundary Faces'
        do j = 1, nbf
            write(99,*) boundff(j,:)
        end do
10      format(2f12.7)        
        close(99)
        deallocate (istr)
    end do
    
    ! ****Call 2D Preprocessor and pass on input parameters**** !
    print *, 'Start Preprocessing'
    pathLin = '/eng/cvcluster/egevansbj/codes/prepro/2Dprepro_duct'
    pathWin = ''   
    do i = 1, NoNests
        ! Determine correct String      
        if (i < 10) then
            allocate(character(len=1) :: istr)
            write( istr, '(I1)' )  i
        else
            allocate(character(len=2) :: istr)
            write( istr, '(I2)' )  i
        end if
        
        if (IsLin == 'N') then
            
            ! write Inputfile (for Windows)
            open(11, file='Input_Data\PreprocessingInput.txt')        
            write(11,*) 'Output_Data\', filename, istr, '.dat'
            write(11,*) 'f'
            write(11,*) 1
            write(11,*) 0
            write(11,*) 0
            write(11,*) 'Output_Data\', filename, istr, '.sol'
            close(11)
            allocate(character(len=29) :: strPrepro)
            strPrepro = 'Flite2D\PreProcessing'
            
        else
            
            ! write command (for Linux)
            PreInt = 59 + 3*len(filename) + 3*len(istr) + len(pathLin)
            allocate(character(len=PreInt) :: strPrepro)
            strPrepro = '/bin/echo -e "Output_Data/'//filename//istr//'.dat\nf\n1\n0\n0\n' &        ! Assemble system command string
            //filename//istr//'.sol\n" | '//pathLin//'/Aggl2d > '//filename//istr//'.outpre'
            
        end if
        print *, 'Preprocessing Snapshot', i
        print *, ' '
        call system(strPrepro)   ! System operating command called to activate fortran       
        deallocate (istr)
        deallocate (strPrepro)
    end do
    print *, 'Finished Preprocessing'
    
    !!!!! IMPLEMENT Mesh Quality Test
    
    ! ****Optimize Mesh by the help of Cuckoo Search and POD**** !
    print *, 'Start Optmization'
    call SubOptimization(NoNests, NoCP, NoDim, cond, InitialNests, MxDisp_Move, np, xmax, hMa, p, Aconst, NoPOMod, NoLeviSteps, NoG, constrain)
    ! Output: Optimized mesh via Cuckoo Search and POD
    
    coord_temp = coord
    call SubGenerateInitialMeshes(NoDim, NoCP, coord_temp, connecf, boundf, coarse, connecc, Coord_CP, Rect, NestOpt)
    ! Output: Optimum Coordinates - 1 Mesh with moved boundaries based on optimum Control Point Coordinates
     
    ! Safe Optimum Geometry in Text File
    open(99, file='Output_Data/OptimumMesh.txt')         
    write(99,'(1I8)') np
    write(99,'(1I8)') ne
    write(99,10) transpose(coord_temp)
11  format(3I8)
    close(99)
    
end program AerOpt