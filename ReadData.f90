module ReadData
    
    use InputData
    use Toolbox
    type ReadVariablesData
    
        integer :: ne, np, nbf, nbc                                             ! Number of elements, Nodes/Points & boundary faces
        integer, dimension(:), allocatable :: boundtype, boundpart              ! Boundary Information of Mesh: Type, Part of boundary    
        integer, dimension(:,:), allocatable :: connec, bound                   ! Connectivity & Boundary Matrix of Fine Mesh    
        double precision, dimension(:,:), allocatable :: coord                  ! Coordinates Matrix of Mesh
        double precision, dimension(:,:), allocatable :: Coord_CN               ! desired Coordinates of the Control Nodes
        double precision, dimension(:,:), allocatable :: coord_temp             ! Coordinates Matrix of Mesh for temporary Mesh Movement
        
    end type ReadVariablesData
    
    type(ReadVariablesData) :: RD
    double precision, dimension(:,:), allocatable :: ArrayTemp
    
contains
      
    subroutine SubReadData()
    ! Objective: Reads the Mesh data and Control Points Coordinates
    
        ! Variables
        implicit none
        integer :: i, NoHeadLines, clean, NoColumns1, NoColumns2
        character(len=5) :: container
        double precision, dimension(:), allocatable :: Input
        
        ! Body of ReadData
        NoColumns1 = 5
        NoColumns2 = 5
        clean = 4
        open(1, file= DataFolder//'/'//trim(IV%MeshFileName)//'.dat', form='formatted',status='old')
        write(*,*) 'File '//trim(IV%MeshFileName)//' opened...'
        write(*,*) "" 
        read(1,*) NoHeadLines
        do i = 1, NoHeadLines
            read(1,*) container
            if (container == 'clean') then
                NoColumns1 = IV%NoDim + 2
                NoColumns2 = IV%NoDim + 1
                clean = IV%NoDim + 1
            end if
        end do
        
        ! Read in Header parameters and allocate arrays
        read(1,*)        
        read(1,*) RD%ne, RD%np, RD%nbf
        allocate(RD%connec(RD%ne,IV%NoDim+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        allocate(RD%coord(RD%np,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        allocate(RD%bound(RD%nbf,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        allocate(RD%boundtype(RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        
        ! Connectivities
        allocate(Input(NoColumns1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        read(1,*) 
        do i = 1, RD%ne
            read(1, *) Input
            RD%connec(i,:) = Input(2:4)
        end do
        deallocate(Input)
        
        ! Coordinates
        allocate(Input(NoColumns2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        read(1,*) 
        do i = 1, RD%np
            read(1, *) Input
            RD%coord(i,:) = Input(2:3)          
        end do

        ! Trash
        if (clean == 4) then
            read(1,*) 
            do i = 1, RD%np
                read(1, *)
            end do 
        end if
        
        ! Boundary Faces
        read(1,*) 
        do i = 1, RD%nbf
            read(1, *) Input
            RD%bound(i,:) = Input(1:2)
            RD%boundtype(i) = Input(clean)
            ! RD%boundpart
        end do
        deallocate(Input)        
        close(1)
    
        allocate(RD%Coord_CN(IV%NoCN,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        open(2, file= DataFolder//'/Control_Nodes.txt', form='formatted',status='old')
        do i = 1, IV%NoCN
            read(2, *) RD%Coord_CN(i,:)
        end do
        close(2)
    
    end subroutine SubReadData
    
    subroutine SubReadData_3D()
    ! Objective: Reads the Mesh data and Control Points Coordinates
    
        ! Variables
        implicit none
        integer :: i, j, NoSurf
        double precision, dimension(:,:), allocatable :: Surf
        
        ! Body of ReadData
        open(1, file= DataFolder//'/'//trim(IV%MeshFileName)//'.plt', form='unformatted',status='old')
        write(*,*) 'File '//trim(IV%MeshFileName)//' opened...'
        write(*,*) ""      
        read(1) RD%ne, RD%np, RD%nbf
        
        allocate(RD%connec(RD%ne,IV%NoDim+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        allocate(RD%coord(RD%np,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        allocate(RD%bound(RD%nbf,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        allocate(RD%boundtype(RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        allocate(RD%boundpart(RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        read(1) RD%connec
        read(1) RD%coord
        read(1) RD%bound, RD%boundpart, RD%boundpart
        close(1)
        
        open(1, file= DataFolder//'/'//trim(IV%MeshFileName)//'.bco', form='formatted',status='old')
        write(*,*) 'File '//trim(IV%MeshFileName)//' opened...'
        write(*,*) ""      
        read(1,*)
        read(1,*) NoSurf
        read(1,*)
        allocate(Surf(NoSurf,2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        do i = 1, NoSurf
            read(1,*) Surf(i,:)
        end do
        close(1)
        do i = 1, RD%nbf
            do j = 1, NoSurf
                if (RD%boundpart(i) == Surf(j,1)) then
                    RD%boundtype(i) = Surf(j,2)
                end if
            end do
        end do
        
        allocate(RD%Coord_CN(IV%NoCN,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        open(2, file= DataFolder//'/Control_Nodes.txt', form='formatted',status='old')
        do i = 1, IV%NoCN
            read(2, *) RD%Coord_CN(i,:)
        end do
        close(2)
    
    end subroutine SubReadData_3D
    
    subroutine CreateFolderStructure()
    
        ! Variables
        implicit none
    
        ! Body of CreateFolderStructure
        call createDirectories()
        if (IV%SystemType == 'W')   then    ! AerOpt is executed from a Windows machine          
            if (IV%RunOnCluster == 'Y') then
                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'psftp')   ! Submits create directory file
            end if
            call system('Communication')    ! Submits create directory file
        else
            call system('chmod a+x ./Communication')    
            call system('./Communication')    ! Submits create directory file     
        end if  
    
    end subroutine CreateFolderStructure
    
    subroutine writeDatFile(ii)
    !Objective: Create Outputfile of each Snapshot as Input for the Pre Processor
    
        ! Variables
        implicit none
        integer :: ii, j
    
        ! Body of writeDatFile
        
        ! Determine correct String      
        call DetermineStrLen(istr, ii) 
        
        open(99, file= trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.dat', form='formatted',status='unknown')
        write(99,*) 2
        write(99,*) 'clean'
        write(99,*) 'David Naumann'
        write(99,*) 'NoElem NoNodes NoBound'        
        write(99,*) RD%ne, RD%np, RD%nbf
        write(99,*) 'Connectivities'
        do j = 1, RD%ne
            write(99,*) j, RD%connec(j,:)
        end do
        write(99,*) 'Coordinates'

        do j = 1, RD%np
            write(99,*) j, RD%coord_temp(j,:)
        end do
        write(99,*) 'Boundary Faces'
        do j = 1, RD%nbf
            write(99,*) RD%bound(j,:), RD%boundtype(j)
        end do
10      format(2f12.7)        
        close(99)
        deallocate (istr)
    
    end subroutine writeDatFile
    
    subroutine readDatFile(ii)
    
        ! Variables
        implicit none
        integer :: i, ii, NoHeadLines, clean, NoColumns1, NoColumns2
        character(len=5) :: container
        double precision, dimension(:), allocatable :: Input
    
        ! Body of readDatFile
        NoColumns1 = 5
        NoColumns2 = 5
        clean = 4
        
        ! Determine correct String      
        call DetermineStrLen(istr, ii) 
        open(1, file= trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.dat', form='formatted',status='unknown')
        read(1,*) NoHeadLines
        do i = 1, NoHeadLines
            read(1,*) container
            if (container == 'clean') then
                NoColumns1 = IV%NoDim + 2
                NoColumns2 = IV%NoDim + 1
                clean = IV%NoDim + 1
            end if
        end do
        deallocate(istr)

        ! Read in Header parameters and allocate arrays
        read(1,*)        
        read(1,*)
        allocate(RD%coord_temp(RD%np,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
         
        ! Connectivities
        allocate(Input(NoColumns1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        read(1,*) 
        do i = 1, RD%ne
            read(1, *) Input
        end do
        deallocate(Input)
         
        ! Coordinates
        allocate(Input(NoColumns2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "

        read(1,*) 
        do i = 1, RD%np
            read(1, *) Input
            RD%coord_temp(i,:) = Input(2:3)          
        end do

	 deallocate(Input)
	 close(1)
       
    end subroutine readDatFile
    
    subroutine writepltFile(ii)
    ! Objective: write 3D Mesh file (.plt) as Output
    
        ! Variables
        implicit none
        integer :: ii, j
    
        ! Body of writepltFile
        ! Determine correct String      
        call DetermineStrLen(istr, ii) 
        
        open(1, file= trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.plt', form='unformatted',status='unknown')    
        write(1) RD%ne, RD%np, RD%nbf
        write(1) RD%connec
        write(1) RD%coord_temp
        write(1) RD%bound, RD%boundpart, RD%boundpart
        close(1)
        deallocate (istr)
        
    end subroutine writepltFile
    
    subroutine writebcoFile(ii)
    
        ! Variables
        implicit none
        integer :: ii
        
        ! Body of writebcoFile
        ! Determine correct String      
        call DetermineStrLen(istr, ii) 
        call system('cp '//DataFolder//'/'//trim(IV%Meshfilename)//'.bco '//'"'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.bco"')
        deallocate(istr)
    
    end subroutine writebcoFile
    
    subroutine PreProInpFile()
    ! Objective: Create the Inputfile for the PreProcessor in Case of the Windows Executable
    
        ! Variables
        implicit none
    
        ! Body of PreProInpFile
        open(11, file=trim(IV%SimulationName)//'/'//InFolder//'/PreprocessingInput.txt', form='formatted',status='unknown')        
        write(11,'(A)') trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.dat'  ! Name of Input file
        write(11,*) 'f'                                     ! Hybrid Mesh?
        write(11,*) 1                                       ! Number of Grids
        write(11,*) 0                                       ! Directionality Parameters
        write(11,*) 0                                       ! Visualization Modes
        write(11,'(A)') trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.sol'  ! Output File name
        close(11)
    
    end subroutine PreProInpFile
    
    subroutine PreProInpFile_3D()
    ! Objective: Create the Inputfile for the PreProcessor in Case of the Windows Executable
    
        ! Variables
        implicit none
    
        ! Body of PreProInpFile
        open(11, file=trim(IV%SimulationName)//'/'//InFolder//'/PreprocessingInput.txt', form='formatted',status='unknown')        
        write(11,*) trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr  ! Name of Input files
        if (IV%Re == 0.0) then
            write(11,*) 1                       ! Inviscid
        else
            write(11,*) 2                       ! Viscid
        end if
        write(11,*) 1                           ! Number of Grids
        if (IV%Re == 0.0) then
            write(11,*) 'f'
        else
            write(11,*) 't'
        end if
        write(11,*) IV%NoProcessors             ! Number of Processors
        write(11,*) 'f'                         ! Rolling ground?
        write(11,*) 1                           ! Starting step
        write(11,*) 1                           ! Number of steps per cycle
        close(11)
    
    end subroutine PreProInpFile_3D
    
    subroutine WriteSolverInpFile()
    ! Objective: Create Solver Input files
    
        ! Variables
        implicit none
        character(len=7) :: strEFMF, strMa      ! String for Ma number & Engine Inlet Front Mass Flow Number
        character(len=2) :: strNI               ! String for Number of Iterations
        character(len=12) :: strRe              ! String for Reynoldsnumber
        character(len=7) :: strAlpha            ! String for Solver inflow direction
        character(len=5) :: strGamma            ! String for Gamma Ratio
        character(len=6) :: strmaxit            ! String for maximum Iterations
        character(len=1) :: strturbulencemodel  ! String for Turbulence model
        character(len=255) :: Output
    
        ! Body of WriteSolverInpFile
        write( strEFMF, '(F7.3)' )  IV%engFMF
        write( strMa, '(F7.3)' )  IV%Ma
        write( strNI, '(I2)' )  IV%NoIter
        write( strRe, '(F12.2)' )  IV%Re
        write( strAlpha, '(F7.2)' )  IV%AlphaInflowDirection
        write( strGamma, '(F5.2)' )  IV%gamma
        write( strmaxit, '(I6)' )  IV%maxit
        write( strturbulencemodel, '(I1)' )  IV%turbulencemodel
        
        
        ! Create Input File for Solver
        open(5, file=trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//'.inp', form='formatted',status='unknown')
        write(5,*) '&inputVariables'
        write(5,*) 'ivd%numberOfMGIterations = ' ,strNI , ','
        write(5,*) 'ivd%alpha = ' ,trim(strAlpha), ','
        write(5,*) 'ivd%MachNumber = ' ,trim(strMa), ','
        write(5,*) 'ivd%ReynoldsNumber = ', trim(strRe), ','
        write(5,*) 'ivd%gamma = ', trim(strGamma), ','
        write(5,*) 'ivd%turbulenceModel = ', trim(strturbulencemodel), ','
        write(5,*) 'ivd%writeToFileInterval = 100,'
        write(5,*) 'ivd%useTimeResidual = .false.,'
        write(5,*) 'ivd%enginesFrontMassFlow = ', trim(strEFMF), ','
        write(5,*) 'ivd%maxitt = ', trim(strmaxit), ','
        ! Others
        write(5,*) 'ivd%numberOfRelaxationSteps = 1, '
        write(5,*) 'ivd%useMatrixDissipation = .false.,'
        write(5,*) 'ivd%numberOfGridsToUse = 1,'
        write(5,*) 'ivd%viscosityScheme = 1,'
        write(5,*) 'ivd%boundaryTerm = 1,'
        write(5,*) 'ivd%CFLNumber = 1.0,'
        write(5,*) 'ivd%turbulenceCFLFactor = 1.0,'
        write(5,*) 'ivd%multigridScheme = 3,'
        write(5,*) 'ivd%tripFactor = 1.0,'
        write(5,*) 'ivd%dissipationScheme = 2,'
        write(5,*) 'ivd%coarseGriddissipationScheme = 2,'
        write(5,*) 'ivd%secondOrderDissipationFactor = 0.3,'
        write(5,*) 'ivd%fourthOrderDissipationFactor = 0.2,'
        write(5,*) 'ivd%coarseGridDissipationFactor = 0.5,'
        write(5,*) 'ivd%turbulenceSmoothingFactor = 0.0,'
        write(5,*) 'ivd%numberOfRSSteps = 0,'
        write(5,*) 'ivd%tripNodes(1) = 203,'
        write(5,*) 'ivd%tripNodes(2) = 252,'
        write(5,*) 'ivd%residualSmoothingFactor = 0.0,'
        write(5,*) 'ivd%numberOfPSSteps = 0,'
        write(5,*) 'ivd%prolongationSmoothingFactor = 0.0,'
        write(5,*) 'ivd%useDissipationWeighting = .false.,'
        write(5,*) 'ivd%normalComponentRelaxation = 1.0,'
        write(5,*) 'ivd%sizeOfSeparationField = 25,'
        write(5,*) 'ivd%numberOfTriperations = 0,'
        write(5,*) '/'
        close(5)
        
        ! Create Read File for Solver
        open(1, file=trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/SolverInput'//istr//'.sh', form='formatted',status='unknown')
        if (IV%SystemType == 'B') then
            write(1,*) InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//'.inp'            ! Control Filename
            write(1,*) InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.sol'      ! Computation Filename
            write(1,*) '' !trim(IV%filename), istr, '.unk'     ! Start-up File = previous Result File
            write(1,*) InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.unk'     ! Result filename
            write(1,*) InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.rsd'      ! Residual Filename
        elseif (IV%SystemType == 'Q' .or. IV%RunOnCluster == 'Y') then
            write(1,'(A)') trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//'.inp'    ! Control Filename
            write(1,'(A)') trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.sol'    ! Solution Filename
            !write(1,'(A)') '' trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk'    ! Start-up File = previous Result File
            write(1,*) ' '
            write(1,'(A)') trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.unk'    ! Result filename
            write(1,'(A)') trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.rsd'    ! Residual Filename
        elseif (IV%SystemType == 'W' .or. IV%SystemType == 'L') then
            write(1,'(A)') trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//'.inp'            ! Control Filename
            write(1,'(A)') trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.sol'      ! Computation Filename
            write(1,'(A)') '' !trim(IV%SimulationName), '/', InFolder, '/', trim(IV%filename), istr, '.unk'     ! Start-up File = previous Result File
            write(1,'(A)') trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.unk'     ! Result filename
            write(1,'(A)') trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.rsd'      ! Residual Filename
        end if
        close(1)
    
    end subroutine WriteSolverInpFile
 
    subroutine WriteSolverInpFile_3D()
    ! Objective: Create Solver Input files
    
        ! Variables
        implicit none
        character(len=7) :: strEFMF, strMa      ! String for Ma number & Engine Inlet Front Mass Flow Number
        character(len=2) :: strNI               ! String for Number of Iterations
        character(len=12) :: strRe              ! String for Reynoldsnumber
        character(len=7) :: strAlpha            ! String for Solver inflow direction
        character(len=5) :: strGamma            ! String for Gamma Ratio
        character(len=5) :: strmaxit            ! String for maximum Iterations
        character(len=1) :: strturbulencemodel  ! String for Turbulence model
        character(len=3) :: strProcessors       ! String for Turbulence model
    
        ! Body of WriteSolverInpFile
        write( strEFMF, '(F7.3)' )  IV%engFMF
        write( strMa, '(F7.3)' )  IV%Ma
        write( strNI, '(I2)' )  IV%NoIter
        write( strRe, '(F12.2)' )  IV%Re
        write( strAlpha, '(F7.2)' )  IV%AlphaInflowDirection
        write( strGamma, '(F5.2)' )  IV%gamma
        write( strmaxit, '(I5)' )  IV%maxit
        write( strturbulencemodel, '(I1)' )  IV%turbulencemodel
        write( strProcessors, '(I3)' )  IV%NoProcessors      
        
        ! Create Input File for Solver
        open(5, file=trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.inp', form='formatted',status='unknown')
        write(5,*) '&inputVariables'
        write(5,*) 'ivd%numberOfProcesses = ',trim(strProcessors) , ','
        write(5,'(A)') 'ivd%dataDirectory = "'//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/"'
        write(5,*) 'ivd%flowType = 0,'                              ! 0 = steady, 1 = unsteady, requires different types of objective functions
        write(5,*) 'ivd%numberOfMGIterations = ' ,strNI , ','
        write(5,*) 'ivd%alpha = ' ,trim(strAlpha), ','
        write(5,*) 'ivd%beta = ' ,trim(strAlpha), ','
        write(5,*) 'ivd%MachNumber = ' ,trim(strMa), ','
        write(5,*) 'ivd%ReynoldsNumber = ', trim(strRe), ','
        write(5,*) 'ivd%gamma = ', trim(strGamma), ','
        write(5,*) 'ivd%turbulenceModel = ', trim(strturbulencemodel), ','   !1 = Spalart-Allmaras, 2=K-epsilon, 3=SST
        write(5,*) 'ivd%restartNumber = 0,'                          !restart, if set to one, it asks for .rst_ files (previous .unk_)
        write(5,*) 'ivd%writeToFileInterval = 100,'
        write(5,*) 'ivd%useTimeResidual = .false.,'
        write(5,*) 'ivd%enginesFrontMassFlow = ', trim(strEFMF), ','
        write(5,*) 'ivd%datum(1) = 0.0,'                            !datum from which moments are computed (e.g. neutral point)
        write(5,*) 'ivd%datum(2) = 0.0,'
        write(5,*) 'ivd%datum(3) = 0.0,'
        write(5,*) 'ivd%groundplane = 0.0,'
        write(5,*) 'ivd%engineFlowType = 2,'                       !this section is just to do with engine inlets / exhausts
        write(5,*) 'ivd%engineBCRelaxation = 1.0,'
        write(5,*) 'ivd%numberOfEORelaxationSteps = 0,'
        write(5,*) 'ivd%movingwall = .true.,'                       !do you have wheels or rolling ground?
        write(5,*) '/'
        close(5)
        
        ! Create Read File for Solver
        open(1, file=trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/run.inp', form='formatted',status='unknown')
            write(1,'(A)') trim(IV%filename)//istr
        close(1)
    
    end subroutine WriteSolverInpFile_3D

    subroutine writeBatchFile(start, ending)
    
        ! Variables
        implicit none
        integer :: i, start, ending
        character(len=:), allocatable :: strProc, strwait, strStart, strEnd, istring, strCores
    
        ! Body of writeBatchFile
        call DetermineStrLen(strStart, start)
        call DetermineStrLen(strCores, (ending-start +1))
        call DetermineStrLen(strEnd, ending)
        call DetermineStrLen(strProc, IV%NoProcessors)
        call DetermineStrLen(strwait, IV%waitmax)  
        if (IV%SystemType /= 'B') then
            open(1, file= trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/batchfile', form='formatted',status='unknown')
            write(1,*) '#PBS -N ' ,trim(IV%filename), istr
            write(1,*) '#PBS -q oh'
            write(1,*) '#PBS -l nodes=', strProc
            write(1,*) '#PBS -l walltime=', strwait, ':00:00'
            write(1,*) '#PBS -l mem=1gb'
            write(1,*) '#PBS -m bea'
            write(1,*) '#PBS -M 717761@swansea.ac.uk'
            write(1,'(A)') pathSolver//' < '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/SolverInput'//istr//'.sh > out'
            write(1,*) 'cd ..'
            write(1,*) 'cd ..'
            write(1,'(A)') 'mv '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.rsd '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//OutFolder// '/'//trim(IV%filename)//istr//'.rsd'
            write(1,'(A)') 'mv '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.unk '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk'
        else
            open(1, file= trim(IV%SimulationName)//'/batchfile', form='formatted',status='unknown')
            write(1,'(A)') '#!/bin/bash --login'
            write(1,'(A)') '#SBATCH --job-name=GeomBatch'
            write(1,'(A)') '#SBATCH -o output.%J'
            write(1,'(A)') '#SBATCH --error=error.%J'
            write(1,'(A)') '#SBATCH -n '//strCores
            write(1,'(A)') '#maximum job time in D-HH:MM'
            write(1,'(A)') '#SBATCH -t '//strwait//':00:00'
            write(1,'(A)') '#SBATCH --mem-per-cpu=4000'
            write(1,'(A)') 'module purge'
            write(1,'(A)') 'module load parallel'
            write(1,'(A)') 'srun="srun -n1 -N1 --exclusive"'
            if (start /= 1) then
                write(1,'(A)') 'parallel="parallel -N 1 --delay .2 -j $SLURM_NTASKS --joblog parallel_joblog2 --resume"'
            else
                write(1,'(A)') 'parallel="parallel -N 1 --delay .2 -j $SLURM_NTASKS --joblog parallel_joblog --resume"'
            end if
            write(1,'(A)') '$parallel "$srun --input='//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//'{1}'//'/SolverInput'//'{1}'//'.sh '//pathSolver//'" ::: {'//strStart//'..'//strEnd//'}'
            write(1,'(A)') ' '
            do i = start, ending
                call DetermineStrLen(istring, i)
                write(1,'(A)') 'mv '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istring//'/'//trim(IV%filename)//istring//'.rsd '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istring//'.rsd'
                write(1,'(A)') 'mv '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istring//'/'//trim(IV%filename)//istring//'.unk '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istring//'.unk'
                deallocate(istring)
            end do
        end if
        close(1)
    
    end subroutine writeBatchFile
    
    subroutine writeBatchFile_3D()
    
        ! Variables
        implicit none
        character(len=:), allocatable :: strProc, strwait, strMemory
    
        ! Body of writeBatchFile
        call DetermineStrLen(strProc, IV%NoProcessors)
        call DetermineStrLen(strwait, IV%waitmax)
        call DetermineStrLen(strMemory, IV%Memory)
        open(1, file= trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/batchfile', form='formatted',status='unknown')  
        if (IV%SystemType /= 'B') then
            write(1,*) '#PBS -N ' ,trim(IV%filename), istr
            write(1,*) '#PBS -q oh'
            write(1,*) '#PBS -l nodes='//strProc
            write(1,*) '#PBS -l walltime='//strwait//':00:00'
            write(1,*) '#PBS -l mem='//strMemory//'gb'
            write(1,*) '#PBS -m bea'
            write(1,*) '#PBS -M 717761@swansea.ac.uk'
            write(1,'(A)') 'cd '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr
            write(1,'(A)') 'mpiexec '//pathSolver
        else
            write(1,*) '#BSUB -J ' ,trim(IV%filename), istr
            write(1,*) '#BSUB -o ' ,trim(IV%filename), istr, '.o'
            write(1,*) '#BSUB -e ' ,trim(IV%filename), istr, '.e'
            !write(1,*) '#BSUB -q <enter a queue>'
            write(1,*) '#BSUB -n '//strProc
            write(1,*) '#BSUB -W '//strwait//':00'
            write(1,'(A)') 'cd '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr
            write(1,'(A)') 'mpirun -machinefile $LSB_DJOB_HOSTFILE -n '//strProc//' '//pathSolver
            !write(1,*) '#BSUB -n 1'
            !#BSUB -R "span[ptile=12]"
            
        end if
        close(1)
    
    end subroutine writeBatchFile_3D
    
    subroutine generateUnkFile()
    
        ! Variables
        implicit none
        integer :: i
    
        ! Body of generateUnkFile
        do i = 1, IV%NoNests
            call DetermineStrLen(istr, i)
            open(1, file= 'Communication', form='formatted',status='unknown')
            write(1,'(A)') trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/plotreg.reg'
            write(1,'(A)') trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.res'
            write(1,'(A)') trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.unk'
            write(1,'(A)') 'f'
            write(1,'(A)') 'f'
            close(1)
            call system(trim(IV%filepath)//'/Executables/makeplot < Communication > /dev/null')
            call system('mv '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.rsd '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd')
            call system('mv '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.unk '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk')
            deallocate(istr)
        end do
        
    end subroutine generateUnkFile
    
    subroutine createDirectories()
    ! Objective: Create Directories from Windows in Linux
    
        ! Variables
        implicit none
        character(len=255) :: currentDir
        character(len=:), allocatable :: Output
        integer :: i
        
        ! Body of createDirectories
        call getcwd(currentDir)
        if (IV%SystemType == 'W') then
            open(1, file='Communication.bat', form='formatted',status='unknown')
        else
            open(1, file='Communication', form='formatted',status='unknown')
        end if
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,'(A)') 'cd '//trim(IV%filepath)
        write(1,*) 'mkdir ', trim(IV%SimulationName)
        write(1,*) 'cd ', trim(IV%SimulationName)
        write(1,*) 'mkdir ', InFolder
        write(1,*) 'mkdir ', OutFolder
        write(1,*) 'mkdir ', TopFolder
        write(1,*) 'cd ', InFolder
        do i = 1, IV%NoSnap
            call DetermineStrLen(istr, i)
            write(1,'(A)') 'mkdir '//trim(IV%filename)//istr
            deallocate(istr)
        end do
        close(1)
    
    end subroutine createDirectories
    
    subroutine transferFilesWin()
    ! Objective: transfer File to send files from Windows to Linux
    
        ! Variables
        implicit none
        character(len=255) :: currentDir
        integer :: i
        
        ! Body of createDirectories
        call getcwd(currentDir)
        open(1, file='Communication', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder
        ! Put Commmands to Transfer Data from Input_Data Folder on Windows to created Input Folder on Cluster
        write(1, '(A)') 'put "'//trim(currentDir)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/SolverInput'//istr//'.sh"'
        if (IV%NoDim == 2) then
            write(1, '(A)') 'put "'//trim(currentDir)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'.sol"'
        elseif (IV%NoDim == 2) then
            do i = 1, IV%NoProcessors
                call DetermineStrLen(istr, i)
                write(1, '(A)') 'put "'//trim(currentDir)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'.sol_'//istr//'"'
                deallocate(istr)
            end do
        end if
        write(1, '(A)') 'put "'//trim(currentDir)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//'.inp"'
        write(1, '(A)') 'put "'//trim(currentDir)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/batchfile.sh"'
        close(1)
    
    end subroutine transferFilesWin
    
    subroutine TriggerFile()
    ! Objectives: Triggerfile for Cluster Simulation & Parallelisation
    
        ! Variables
        implicit none
    
        ! Body of TriggerFile
        open(1, file='Communication', form='formatted',status='unknown')
        
        if (IV%SystemType == 'Q') then
            write(1,*) 'cd '
            write(1,*) 'cd ..'
            write(1,'(A)') 'cd '//trim(IV%filepath)//'/'//trim(IV%SimulationName)

            if (IV%SystemType /= 'B') then
                write(1,*) 'qsub batchfile'
            else
                write(1,*) 'sbatch batchfile'
            end if
        end if
        close(1)
          
    end subroutine TriggerFile
    
    subroutine TriggerFileQ()
    ! Objectives: Triggerfile for Cluster Simulation & Parallelisation
    
        ! Variables
        implicit none
    
        ! Body of TriggerFile
        open(1, file='Communication', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,'(A)') 'cd '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr 
        
        if (IV%SystemType /= 'B') then
            write(1,*) 'qsub batchfile'
        else
            write(1,*) 'sbatch batchfile'
        end if
        close(1)
          
    end subroutine TriggerFileQ
    
    subroutine TriggerFile2()
    ! Objectives: Triggerfile for noncluster Simulation
    
        ! Variables
        implicit none
    
        ! Body of TriggerFile
        open(1, file='Communication', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ..'
        write(1,*) 'cd ..'
        write(1,*) pathSolver
        close(1)
                        
    end subroutine TriggerFile2
    
    subroutine TransferSolutionOutput()
    
        ! Variables
        implicit none
        character(len=:), allocatable :: Output
        integer :: strOut      
    
        ! Body of TransferSolutionOutput
        open(1, file='Communication', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//OutFolder
        write(1, '(A)') 'get '//trim(IV%filename)//istr//'.unk' ! Put Commmands to Transfer Data from Output_Data Folder on Linux to Output_Folder on Windows
        write(1, '(A)') 'get '//trim(IV%filename)//istr//'.rsd'
        close(1)
        
    end subroutine TransferSolutionOutput
    
    subroutine DeleteLogFiles()
    
        ! Variables
        implicit none
    
        ! Body of DeleteLogFiles
        open(1, file='Communication', form='formatted',status='unknown')
        write(1,'(A)') 'chmod 777 '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/parallel_joblog*'
        write(1,'(A)') 'chmod 777 '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/error.*'
        write(1,'(A)') 'chmod 777 '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/output.*'
        write(1,'(A)') 'rm '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/error.*'
        write(1,'(A)') 'rm '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/output.*'
        write(1,'(A)') 'rm '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/parallel_joblog*'
        close(1)
    
    end subroutine DeletelogFiles

    subroutine DeleteErrorFiles()
    
        ! Variables
        implicit none
    
        ! Body of DeleteErrorFiles
        open(1, file='Communication', form='formatted',status='unknown')
        write(1,'(A)') 'find '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//' -name "*.e*" -type f -delete'
        write(1,'(A)') 'find '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//' -name "*.o*" -type f -delete'
        close(1)
    
    end subroutine DeleteErrorFiles
       
    subroutine moveTopNestFilesLin(i, idx)

        ! Variables
        implicit none
        integer :: i, idx
        character(len=255) :: Output
        character(len=:), allocatable :: NoGenstr
        character(len=:), allocatable :: indexStr
        
        ! Body of moveTopNestFiles
        call DetermineStrLen(istr, i)
        call DetermineStrLen(NoGenstr, OV%Gen) 
        call DetermineStrLen(indexStr, idx) 
        open(1, file='Communication', form='formatted',status='unknown')
        write(1, '(A)') 'mv '//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk "'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'_'//indexStr//'.unk"'
        write(1, '(A)') 'mv '//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd "'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'_'//indexStr//'.rsd"'
        if (IV%NoDim == 2) then
            write(1, '(A)') 'mv '//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.dat "'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'_'//indexStr//'.dat"'
        else
            write(1, '(A)') 'mv '//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.plt "'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'_'//indexStr//'.plt"'
        end if
        close(1)
        deallocate(istr)
        deallocate(NoGenstr)
        call system('chmod a+x Communication')
        call system('./Communication')    ! Submits Move file
        
    end subroutine moveTopNestFilesLin

    subroutine moveTopNestFilesLin_POD(i)

        ! Variables
        implicit none
        integer :: i
        character(len=255) :: Output
        character(len=:), allocatable :: NoGenstr
        
        ! Body of moveTopNestFiles
        call DetermineStrLen(istr, i)
        call DetermineStrLen(NoGenstr, OV%Gen) 
        open(1, file='Communication2', form='formatted',status='unknown')
	 write(1, '(A)') 'mv '//trim(IV%SimulationName)//'/POD_Pressure'//istr//'.txt "'//trim(IV%SimulationName)//'/'//TopFolder//'/POD_Pressure'//NoGenstr//'.txt"'
        close(1)
        deallocate(istr)
        deallocate(NoGenstr)
        call system('chmod a+x Communication2')
        call system('./Communication2')    ! Submits Move file
       
    end subroutine moveTopNestFilesLin_POD
    
    subroutine moveTopNestFilesWin(i, idx)

        ! Variables
        implicit none
        integer :: i, idx
        character(len=:), allocatable :: NoGenstr
        character(len=:), allocatable :: indexStr
        
        ! Body of moveTopNestFiles
        call DetermineStrLen(istr, i)
        call DetermineStrLen(NoGenstr, OV%Gen)
        call DetermineStrLen(indexStr, idx) 
        open(1, file='Communication.bat', form='formatted',status='unknown')
        write(1, '(A)') 'move '//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk "'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'_'//indexStr//'.unk"'
        write(1, '(A)') 'move '//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd "'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'_'//indexStr//'.rsd"'
        if (IV%NoDim == 2) then
            write(1, '(A)') 'move '//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.dat "'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'_'//indexStr//'.dat"'
        else
            write(1, '(A)') 'move '//trim(IV%SimulationName)//'\'//InFolder//'\'//trim(IV%filename)//istr//'\'//trim(IV%filename)//istr//'.plt "'//trim(IV%SimulationName)//'\'//TopFolder//'\'//trim(IV%filename)//NoGenstr//'.plt"'
        end if
        close(1)
        deallocate(istr)
        deallocate(NoGenstr)
        call system('Communication.bat >nul 2>&1')    ! Submits create directory file
        
    end subroutine moveTopNestFilesWin
    
    subroutine generateEnSightFileLin()
    
       ! Variables
        implicit none
        character(len=:), allocatable :: NoGenstr
        LOGICAL :: exe_exists

        ! Body of generateEnSightFile
        inquire(file=trim(IV%filepath)//'/Executables/engen_2D', exist=exe_exists)
        if(.not.exe_exists) then
            return
        end if

        call DetermineStrLen(NoGenstr, OV%Gen)
        call system('cp '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'_1.dat '//trim(IV%filepath)//'/'//trim(IV%filename)//NoGenstr//'.dat')
        call system('cp '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'_1.unk '//trim(IV%filepath)//'/'//trim(IV%filename)//NoGenstr//'.resp')
        open(1, file= 'Communication', form='formatted',status='unknown')
        write(1,'(A)') trim(IV%filename)//NoGenstr
        write(1,'(A)') trim(IV%filename)//NoGenstr
        write(1,'(A)') 'f'
        write(1,*) IV%Ma
        write(1,*) IV%Tamb
        write(1,*) IV%Pamb
        write(1,'(A)') '1017'
        close(1)
        call system(trim(IV%filepath)//'/Executables/engen_2D < Communication > /dev/null')
        call system('mv *ENSIGHT* '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//TopFolder)
        call system('rm '//trim(IV%filepath)//'/'//trim(IV%filename)//NoGenstr//'.dat')
        call system('rm '//trim(IV%filepath)//'/'//trim(IV%filename)//NoGenstr//'.resp')
        deallocate(NoGenstr)       
        
    end subroutine generateEnSightFileLin
    
    subroutine generateEnSightFileLin3D(i)
    
       ! Variables
        implicit none
        integer :: i
        character(len=:), allocatable :: NoGenstr
    
        ! Body of generateEnSightFile
        call DetermineStrLen(NoGenstr, OV%Gen)
        call DetermineStrLen(istr, i)
        call system('mv '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/base.plt '//trim(IV%filepath)//'/'//trim(IV%filename)//NoGenstr//'.plt')
        call system('mv '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'.unk '//trim(IV%filepath)//'/'//trim(IV%filename)//NoGenstr//'.unk')
        open(1, file= 'Communication', form='formatted',status='unknown')
        write(1,'(A)') trim(IV%filename)//NoGenstr
  ! Always hybrid?
        write(1,'(A)') 't'
        write(1,*) IV%Ma
        write(1,*) IV%Tamb
        write(1,*) IV%Pamb
        write(1,'(A)') '1017'
        close(1)
        call system(trim(IV%filepath)//'/Executables/engen_3D < Communication > /dev/null')
        call system('mv *ENSIGHT* '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//TopFolder)
        call system('rm '//trim(IV%filepath)//'/'//trim(IV%filename)//NoGenstr//'.plt')
        call system('mv '//trim(IV%filepath)//'/'//trim(IV%filename)//NoGenstr//'.unk '//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'.unk ')
        deallocate(NoGenstr)       
        
    end subroutine generateEnSightFileLin3D
        
    subroutine generateEnSightFileWin()
! ALSO 3D 
       ! Variables
        implicit none
        character(len=:), allocatable :: NoGenstr
        logical :: exe_exists
    
        ! Body of generateEnSightFile
        inquire(file=trim(IV%filepath)//'/Executables/engen_2D', exist=exe_exists)
        if(.not.exe_exists) then
            return
        end if

        call DetermineStrLen(NoGenstr, OV%Gen)
        open(1, file= 'Communication', form='formatted',status='unknown')
        write(1,'(A)') trim(IV%filename)//NoGenstr//'.dat'
        write(1,'(A)') trim(IV%filename)//NoGenstr//'.resp'
        write(1,'(A)') 'f'
        write(1,*) IV%Ma
        write(1,*) IV%Tamb
        write(1,*) IV%Pamb
        write(1,'(A)') '1017'
        close(1)
        open(1, file= 'Communication.bat', form='formatted',status='unknown')
        write(1, '(A)') 'copy '//trim(IV%SimulationName)//'\'//TopFolder//'\'//trim(IV%filename)//NoGenstr//'_1.dat "'//trim(IV%filename)//NoGenstr//'.dat"'
        write(1, '(A)') 'copy '//trim(IV%SimulationName)//'\'//TopFolder//'\'//trim(IV%filename)//NoGenstr//'_1.unk "'//trim(IV%filename)//NoGenstr//'.resp"'
        write(1,'(A)') 'Executables\engen_2D < Communication >nul 2>&1'
        write(1,'(A)') 'move *ENSIGHT* '//trim(IV%SimulationName)//'\'//TopFolder
        write(1,'(A)') 'del '//trim(IV%filename)//NoGenstr//'.dat'
        write(1,'(A)') 'del '//trim(IV%filename)//NoGenstr//'.resp'
        close(1)
        call system('Communication.bat >nul 2>&1')
        deallocate(NoGenstr)  
        
    end subroutine generateEnSightFileWin
    
    subroutine writeFitnessandNestoutput()
    
        ! Variables
        implicit none
        logical :: ex
    
        ! Body of writeFitnessandNestoutput
        allocate(character(len=3) :: istr)
        write(istr, '(1f3.1)') IV%Ma
        
        inquire(file=trim(IV%SimulationName)//'/Nests'//istr//'.txt', exist = ex)       
        if (ex == .true.) then
            open(19,file=trim(IV%SimulationName)//'/FitnessAll.txt',form='formatted',status='old',position='append') ! Comment for GUI
            open(29,file=trim(IV%SimulationName)//'/Nests'//istr//'.txt',form='formatted',status='old',position='append')
        else
            open(29,file=trim(IV%SimulationName)//'/Nests'//istr//'.txt', form='formatted',status='new')
            write(29, *) 'Snapshots'
            write(29,'(<IV%NoSnap>f17.10)') CS%Snapshots
            open(19,file=trim(IV%SimulationName)//'/FitnessAll.txt', form='formatted',status='new')                   ! Comment for GUI
        end if
        
        ! Store moved Nests in Output Analysis File
        write(29, *) 'Generation', OV%Gen
        write(29,'(<IV%NoNests>f17.10)') OV%Nests
        close(29)
      
        ! Store Fitness values
        write(19,'(1I4)',advance="no") OV%Gen
        write(19,'(<IV%NoNests>f17.10)') OV%Fi       
        if (IV%objectivefunction == 7) then
            write(19,'(<IV%NoNests>f17.10)') OV%Precoutput
        end if
        close(19)
        deallocate(istr)
        
    end subroutine writeFitnessandNestoutput
    
    end module ReadData
