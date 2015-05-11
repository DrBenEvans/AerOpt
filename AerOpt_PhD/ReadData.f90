module ReadData
    
    use InputData
    use Toolbox
    type ReadVariablesData
    
        integer :: ne, np, nbf, nbc, NoParts                                    ! Number of elements, Nodes/Points & boundary faces
        integer, dimension(:), allocatable :: boundtype, boundpart, MovingParts ! Boundary Information of Mesh: Type, Part of boundary    
        integer, dimension(:,:), allocatable :: connec, bound                   ! Connectivity & Boundary Matrix of Fine Mesh    
        double precision, dimension(:,:), allocatable :: coord                  ! Coordinates Matrix of Mesh
        double precision, dimension(:,:), allocatable :: Coord_CN               ! desired Coordinates of the Control Points
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
        character(len=200) :: MeshFileName
        character(len=5) :: container
        double precision, dimension(:), allocatable :: Input
        
        ! Body of ReadData
        NoColumns1 = 5
        NoColumns2 = 5
        clean = 4
        !write(*,'(A)',advance="no") "Enter Mesh Filename: "
        !read(*,'(A)') MeshFileName
        MeshFileName = 'Mesh.dat'
        open(1, file= DataFolder//'/'//trim(MeshFileName), form='formatted',status='old')
        write(*,*) 'File '//trim(MeshFileName)//' opened...'
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
        
        ! If only specified parts wish to be moved, include this back into the .dat file!
        !read(1,'(1I8)') RD%NoParts
        RD%NoParts = 0 ! Delete to specify parts
        if (RD%NoParts /= 0) then
            allocate(RD%MovingParts(RD%NoParts),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
            do i = 1, RD%NoParts
                read(1, *) RD%MovingParts(i)
            end do
        end if
        
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
        close(1)
    
        allocate(RD%Coord_CN(IV%NoCN,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        open(2, file= DataFolder//'/Control_Nodes.txt', form='formatted',status='old')
        do i = 1, IV%NoCN
            read(2, *) RD%Coord_CN(i,:)
        end do
        close(2)
    
    end subroutine SubReadData
    
    subroutine CreateFolderStructure()
    
        ! Variables
        implicit none
    
        ! Body of CreateFolderStructure
        call createDirectoriesInit()
        if (IV%SystemType == 'W' .and. IV%RunOnCluster == 'Y')   then    ! AerOpt is executed from a Windows machine connected to a Linux machine
        
            call createDirectoriesInit()
            call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')   ! Submits create directory file
            call createDirectoriesWindows()
            call system('FileCreateDir.bat')    ! Submits create directory file
     
        elseif (IV%SystemType == 'W') then ! AerOpt is executed from a Windows machine alone
        
            call createDirectoriesWindows()
            call system('FileCreateDir.bat')    ! Submits create directory file
        
        elseif (IV%SystemType == 'L') then ! AerOpt is executed from a Windows machine alone
        
            call createDirectoriesWindows()
            call system('chmod a+x ./FileCreateDir.bat')  
            call system('./FileCreateDir.bat')    ! Submits create directory file
        
        elseif (IV%SystemType == 'Q' .or. IV%SystemType == 'B') then     ! AerOpt is executed from a Linux machine
        
            call createDirectoriesInit()
            call system('chmod a+x ./FileCreateDir.scr')    
            call system('./FileCreateDir.scr')    ! Submits create directory file
            
        else       
            STOP 'INPUT ERROR: System Type selected does not exist! Program stopped.'       
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
        
        open(99, file= newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'.dat', form='formatted',status='unknown')
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
    
    subroutine PreProInpFile()
    ! Objective: Create the Inputfile for the PreProcessor in Case of the Windows Executable
    
        ! Variables
        implicit none
    
        ! Body of PreProInpFile
        open(11, file=newdir//'/'//InFolder//'/PreprocessingInput.txt', form='formatted',status='unknown')        
        write(11,*) newdir//'/'//InFolder, '/', trim(IV%filename), istr, '.dat'  ! Name of Input file
        write(11,*) 'f'                                     ! Hybrid Mesh?
        write(11,*) 1                                       ! Number of Grids
        write(11,*) 0                                       ! Directionality Parameters
        write(11,*) 0                                       ! Visualization Modes
        write(11,*) newdir//'/'//InFolder, '/', trim(IV%filename), istr, '.sol'  ! Output File name
        close(11)
    
    end subroutine PreProInpFile
    
    subroutine WriteSolverInpFile()
    ! Objective: Create Solver Input files
    
        ! Variables
        implicit none
        character(len=:), allocatable :: fname
        character(len=7) :: strEFMF, strMa      ! String for Ma number & Engine Inlet Front Mass Flow Number
        character(len=2) :: strNI               ! String for Number of Iterations
        character(len=12) :: strRe              ! String for Reynoldsnumber
        character(len=7) :: strAlpha            ! String for Solver inflow direction
        character(len=5) :: strGamma            ! String for Gamma Ratio
        character(len=5) :: strmaxit            ! String for maximum Iterations
        character(len=1) :: strturbulencemodel  ! String for Turbulence model
        character(len=255) :: Output
        integer :: fileLength
    
        ! Body of WriteSolverInpFile
        fileLength = len(IV%filename) + len(istr) + 15
        allocate(character(len=fileLength) :: fname)
        fname = newdir//'/'//InFolder//'/'//trim(IV%filename)//'.inp'
        write( strEFMF, '(F7.3)' )  IV%engFMF
        write( strMa, '(F7.3)' )  IV%Ma
        write( strNI, '(I2)' )  IV%NoIter
        write( strRe, '(F12.2)' )  IV%Re
        write( strAlpha, '(F7.2)' )  IV%AlphaInflowDirection
        write( strGamma, '(F5.2)' )  IV%gamma
        write( strmaxit, '(I5)' )  IV%maxit
        write( strturbulencemodel, '(I1)' )  IV%turbulencemodel
        
        
        ! Create Input File for Solver
        open(5, file=fname, form='formatted',status='unknown')
        write(5,*) '&inputVariables'
        write(5,*) 'ivd%numberOfMGIterations = ' ,strNI , ','
        write(5,*) 'ivd%numberOfGridsToUse = 1,'
        write(5,*) 'ivd%viscosityScheme = 1,'       ! how to turn viscosity off?
        write(5,*) 'ivd%boundaryTerm = 1,'          ! ??
        write(5,*) 'ivd%CFLNumber = 1.0,'
        write(5,*) 'ivd%turbulenceCFLFactor = 1.0,'
        write(5,*) 'ivd%alpha = ' ,trim(strAlpha), ','
        write(5,*) 'ivd%MachNumber = ' ,trim(strMa), ','
        write(5,*) 'ivd%numberOfRelaxationSteps = 1, '
        write(5,*) 'ivd%ReynoldsNumber = ', trim(strRe), ','
        write(5,*) 'ivd%gamma = ', trim(strGamma), ','
        write(5,*) 'ivd%turbulenceModel = ', trim(strturbulencemodel), ','
        write(5,*) 'ivd%multigridScheme = 3,'
        write(5,*) 'ivd%tripFactor = 1.0,'
        write(5,*) 'ivd%dissipationScheme = 2,'
        write(5,*) 'ivd%coarseGriddissipationScheme = 2,'
        write(5,*) 'ivd%secondOrderDissipationFactor = 0.2,'
        write(5,*) 'ivd%fourthOrderDissipationFactor = 0.2,'
        write(5,*) 'ivd%coarseGridDissipationFactor = 0.5,'
        write(5,*) 'ivd%turbulenceSmoothingFactor = 0.0,'
        write(5,*) 'ivd%numberOfRSSteps = 0,'
        write(5,*) 'ivd%tripNodes(1) = 1,'                ! ?? very specific?
        write(5,*) 'ivd%tripNodes(2) = 3,'                ! ?? very specific?
        write(5,*) 'ivd%residualSmoothingFactor = 0.0,'
        write(5,*) 'ivd%numberOfPSSteps = 0,'
        write(5,*) 'ivd%prolongationSmoothingFactor = 0.0,'
        write(5,*) 'ivd%useMatrixDissipation = .false.,'
        write(5,*) 'ivd%writeToFileInterval = 100,'
        write(5,*) 'ivd%useTimeResidual = .false.,'
        write(5,*) 'ivd%useDissipationWeighting = .false.,'
        write(5,*) 'ivd%normalComponentRelaxation = 1.0,'
        write(5,*) 'ivd%sizeOfSeparationField = 25,'
        write(5,*) 'ivd%numberOfTriperations = 0,'
        write(5,*) 'ivd%enginesFrontMassFlow = ', trim(strEFMF), ','
        write(5,*) 'ivd%maxitt = ', trim(strmaxit), ','
        write(5,*) 'ivd%numberOfDissipationLayers = 0,'
        write(5,*) 'ivd%HighOrder=.false.,'
        write(5,*) '/'
        close(5)
        
        ! Create Read File for Solver
        open(1, file=newdir//'/'//InFolder//'/SolverInput'//istr//'.sh', form='formatted',status='unknown')
        if (IV%SystemType == 'B') then
            write(1,*) trim(IV%filename), '.inp'            ! Control Filename
            write(1,*) trim(IV%filename), istr, '.sol'      ! Computation Filename
            write(1,*) '' !trim(IV%filename), istr, '.resp'     ! Start-up File = previous Result File
            write(1,*) trim(IV%filename), istr, '.resp'     ! Result filename
            write(1,*) trim(IV%filename), istr, '.rsd'      ! Residual Filename
        elseif (IV%SystemType == 'Q' .or. IV%RunOnCluster == 'Y') then
            Output = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//'.inp'    ! Control Filename
            write(1,'(A)') trim(Output)
            Output = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'.sol'    ! Control Filename
            write(1,'(A)') trim(Output)
            Output = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.resp'    ! Start-up File = previous Result File
            write(1,'(A)') '' !trim(Output) 
            Output = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.resp'    ! Result filename
            write(1,'(A)') trim(Output)
            Output = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd'    ! Residual Filename
            write(1,'(A)') trim(Output)
        elseif (IV%SystemType == 'W' .or. IV%SystemType == 'L') then
            write(1,*) newdir, '/', InFolder, '/', trim(IV%filename), '.inp'            ! Control Filename
            write(1,*) newdir, '/', InFolder, '/', trim(IV%filename), istr, '.sol'      ! Computation Filename
            write(1,*) '' !newdir, '/', InFolder, '/', trim(IV%filename), istr, '.resp'     ! Start-up File = previous Result File
            write(1,*) newdir, '/', InFolder, '/', trim(IV%filename), istr, '.resp'     ! Result filename
            write(1,*) newdir, '/', InFolder, '/', trim(IV%filename), istr, '.rsd'      ! Residual Filename
        end if
        close(1)
    
    end subroutine WriteSolverInpFile
    
    subroutine writeBatchFile()
    
        ! Variables
        implicit none
        character(len=255) :: Output
    
        ! Body of writeBatchFile
        open(1, file= newdir//'/'//InFolder//'/batchfile.sh', form='formatted',status='unknown')  
        if (IV%SystemType /= 'B') then
            write(1,*) '#PBS -N ' ,trim(IV%filename), istr
            write(1,*) '#PBS -q oh'
            write(1,*) '#PBS -l nodes=1:ppn=1'
            write(1,*) '#PBS -l walltime=24:00:00'
            write(1,*) '#PBS -l mem=1gb'
            write(1,*) '#PBS -m bea'
            write(1,*) '#PBS -M 717761@swansea.ac.uk'
            Output = pathSolver//' < /eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/SolverInput'//istr//'.sh'
            write(1,'(A)') trim(Output)
        else
            write(1,*) '#BSUB -J ' ,trim(IV%filename), istr
            write(1,*) '#BSUB -o ' ,trim(IV%filename), istr, '.o'
            write(1,*) '#BSUB -e ' ,trim(IV%filename), istr, '.e'
            !write(1,*) '#BSUB -q <enter a queue>'
            write(1,*) '#BSUB -n 1'
            write(1,*) '#BSUB -W 24:00'
            Output = pathSolver//' < /home/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/SolverInput'//istr//'.sh'
            write(1,'(A)') trim(Output)
            write(1,*) 'cd ..'
            write(1,*) 'mv ', InFolder, '/', trim(IV%filename), istr, '.rsd ', OutFolder, '/', trim(IV%filename), istr, '.rsd'
            write(1,*) 'mv ', InFolder, '/', trim(IV%filename), istr, '.resp ', OutFolder, '/', trim(IV%filename), istr, '.resp'
            write(1,*) '#BSUB -n 1'
            !#BSUB -R "span[ptile=12]"
            
        end if
        close(1)
    
    end subroutine writeBatchFile
    
    subroutine createDirectoriesInit()
    ! Objective: Create Directories from Windows in Linux
    
        ! Variables
        implicit none
        character(len=255) :: currentDir
        character(len=:), allocatable :: Output
        integer :: strOut
        
        ! Body of createDirectories
        call getcwd(currentDir)
        open(1, file='FileCreateDir.scr', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%UserName)
        write(1,*) 'chmod 777 AerOpt'
        write(1,*) 'cd AerOpt'
        write(1,*) 'mkdir ', newdir
        write(1,*) 'cd ', newdir
        write(1,*) 'mkdir ', InFolder
        write(1,*) 'mkdir ', OutFolder
        write(1,*) 'mkdir ', TopFolder
        write(1,*) 'cd ..'
        write(1,*) 'cd ..'
        write(1,*) 'chmod 711 AerOpt'
        close(1)
    
    end subroutine createDirectoriesInit
    
    subroutine createDirectoriesWindows()
    ! Objective: Create Directories in Windows
    
        ! Variables
        implicit none
        character(len=255) :: currentDir
        character(len=:), allocatable :: Output
        integer :: strOut
        
        ! Body of createDirectories
        call getcwd(currentDir)
        open(1, file='FileCreateDir.bat', form='formatted',status='unknown')
        write(1,*) 'mkdir ', newdir
        write(1,*) 'cd ', newdir
        write(1,*) 'mkdir ', InFolder
        write(1,*) 'mkdir ', OutFolder
        write(1,*) 'mkdir ', TopFolder
        write(1,*) 'cd ..'
        write(1,*) 'cd ..'
        close(1)
    
    end subroutine createDirectoriesWindows
    
    subroutine transferFilesWin()
    ! Objective: Create Directories from Windows in Linux
    
        ! Variables
        implicit none
        character(len=255) :: currentDir
        character(len=:), allocatable :: Output
        integer :: strOut
        
        ! Body of createDirectories
        call getcwd(currentDir)
        open(1, file='FileCreateDir.scr', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%UserName) , '/AerOpt/', newdir, '/', InFolder
        
        ! Put Commmands to Transfer Data from Input_Data Folder on Windows to created Input Folder on Cluster
        strOut = len(trim(currentDir)) + 33
        allocate(character(len=strOut) :: Output)
        Output = 'put "'//trim(currentDir)//'/'//newdir//'/'//InFolder//'/SolverInput'//istr//'.sh"'
        write(1, '(A)') Output       
        deallocate(Output)
                      
        strOut = len(trim(currentDir)) + 22
        allocate(character(len=strOut) :: Output)
        Output = 'put "'//trim(currentDir)//'/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'.sol"'
        write(1, '(A)') Output
        deallocate(Output)

        strOut = len(trim(currentDir)) + 23
        allocate(character(len=strOut) :: Output)
        Output = 'put "'//trim(currentDir)//'/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//'.inp"'
        write(1, '(A)') Output
        deallocate(Output)
        
        strOut = len(trim(currentDir)) + 25
        allocate(character(len=strOut) :: Output)
        Output = 'put "'//trim(currentDir)//'/'//newdir//'/'//InFolder//'/batchfile.sh"'
        write(1, '(A)') Output
        close(1)
    
    end subroutine transferFilesWin
    
    
    subroutine TriggerFile()
    ! Objectives: Triggerfile for Cluster Simulation & Parallelisation
    
        ! Variables
        implicit none
    
        ! Body of TriggerFile
        open(1, file='Trigger.sh', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%UserName), '/AerOpt/', newdir, '/', InFolder
        
        if (IV%SystemType /= 'B') then
            write(1,*) 'qsub batchfile.sh'
        else
            write(1,*) 'bsub < batchfile.sh'
        end if
        close(1)
                        
    end subroutine TriggerFile
    
    subroutine TriggerFile2()
    ! Objectives: Triggerfile for noncluster Simulation
    
        ! Variables
        implicit none
    
        ! Body of TriggerFile
        open(1, file='Trigger.sh', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ..'
        write(1,*) 'cd ..'
        write(1,*) pathSolver
        close(1)
                        
    end subroutine TriggerFile2
    
    subroutine CheckSimStatus()
    ! Objectives: writes file that checks for error files and stores response in check.txt file
    
        ! Variables
        implicit none
        character(len=255) :: Output
    
        ! Body of CheckSimStatus
        open(1, file='CheckStatus.scr', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%UserName), '/AerOpt/', newdir, '/', InFolder      
        write(1,*) '[ -e ', trim(IV%filename), istr, '.e* ] && echo 1 > check.txt || echo 0 > check.txt'
        if (IV%SystemType == 'Q')   then
            write(1,*) 'cd ..'
            write(1,*) 'cd ..'
            Output = 'mv /eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/check.txt check.txt'
            write(1, '(A)') trim(Output)
        elseif (IV%SystemType == 'B') then
            write(1,*) 'cd ..'
            write(1,*) 'cd ..'
            Output = 'mv /home/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/check.txt check.txt'
            write(1, '(A)') trim(Output)
        end if
        close(1)
    
    end subroutine CheckSimStatus
    
    subroutine CheckSimStatus2()
    ! Objectives: Transfers Check file to correct folder
    
        ! Variables
        implicit none
    
        ! Body of CheckSimStatus
        open(1, file='CheckStatus.scr', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'       
        write(1,*) 'cd ', trim(IV%UserName), '/AerOpt/', newdir, '/', InFolder
        write(1,*) 'get check.txt'       
        close(1)
    
    end subroutine CheckSimStatus2
    
    subroutine TransferSolutionOutput()
    
        ! Variables
        implicit none
        character(len=:), allocatable :: Output
        integer :: strOut      
    
        ! Body of TransferSolutionOutput
        open(1, file='FileCreateDir.scr', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%UserName), '/AerOpt/', newdir, '/', OutFolder
        
        ! Put Commmands to Transfer Data from Output_Data Folder on Linux to Output_Folder on Windows
        strOut = len(trim(IV%filename)) + len(istr) + 9
        allocate(character(len=strOut) :: Output)
        Output = 'get '//trim(IV%filename)//istr//'.resp'
        write(1, '(A)') Output       
        deallocate(Output)
        
        strOut = len(trim(IV%filename)) + len(istr) + 8
        allocate(character(len=strOut) :: Output)
        Output = 'get '//trim(IV%filename)//istr//'.rsd'
        write(1, '(A)') Output       
        deallocate(Output)
        
    end subroutine TransferSolutionOutput
    
    subroutine DeleteErrorFiles(i)
    
        ! Variables
        implicit none
        character(len=*) :: i
    
        ! Body of DeleteErrorFiles
        open(1, file='FileCreateDir.scr', form='formatted',status='unknown')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%UserName) , '/AerOpt/', newdir, '/', InFolder
        write(1,*) ' chmod 777 ', trim(IV%filename), i, '.e*'
        write(1,*) ' rm ', trim(IV%filename), i, '.e*'
        close(1)
    
    end subroutine DeleteErrorFiles
    
       
    subroutine moveTopNestFilesLin(i, NoGen)

        ! Variables
        implicit none
        integer :: i, NoGen
        character(len=255) :: Output
        character(len=:), allocatable :: NoGenstr
        
        ! Body of moveTopNestFiles
        call DetermineStrLen(istr, i)
        call DetermineStrLen(NoGenstr, NoGen)
        open(1, file='FileCreateDir.scr', form='formatted',status='unknown')
        
		Output = 'mv '//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.resp "'//newdir//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'.resp"'
        write(1, '(A)') trim(Output)
        Output = 'mv '//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd "'//newdir//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'.rsd"'
        write(1, '(A)') trim(Output)
        Output = 'mv '//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'.dat "'//newdir//'/'//TopFolder//'/'//trim(IV%filename)//NoGenstr//'.dat"'
        write(1, '(A)') trim(Output)
        
		!Output = 'mv '//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.resp "'//TopFolder//'/'//trim(IV%filename)//'_'//NoGenstr//'.resp"'
  !      write(1, '(A)') trim(Output)
  !      Output = 'mv '//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd "'//TopFolder//'/'//trim(IV%filename)//'_'//NoGenstr//'.rsd"'
  !      write(1, '(A)') trim(Output)
  !      Output = 'mv '//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'.dat "'//TopFolder//'/'//trim(IV%filename)//'_'//NoGenstr//'.dat"'
  !      write(1, '(A)') trim(Output)
        
        close(1)
        deallocate(istr)
        deallocate(NoGenstr)
        
    end subroutine moveTopNestFilesLin
    
    subroutine moveTopNestFilesWin(i, NoGen)

        ! Variables
        implicit none
        integer :: i, NoGen
        character(len=255) :: Output
        character(len=:), allocatable :: NoGenstr
        
        ! Body of moveTopNestFiles
        call DetermineStrLen(istr, i)
        call DetermineStrLen(NoGenstr, NoGen)
        open(1, file='FileCreateDir.bat', form='formatted',status='unknown')
        
        Output = 'move '//newdir//'\'//OutFolder//'\'//trim(IV%filename)//istr//'.resp "'//newdir//'\'//TopFolder//'\'//trim(IV%filename)//NoGenstr//'.resp"'
        write(1, '(A)') trim(Output)
        Output = 'move '//newdir//'\'//OutFolder//'\'//trim(IV%filename)//istr//'.rsd "'//newdir//'\'//TopFolder//'\'//trim(IV%filename)//NoGenstr//'.rsd"'
        write(1, '(A)') trim(Output)
        Output = 'move '//newdir//'\'//InFolder//'\'//trim(IV%filename)//istr//'.dat "'//newdir//'\'//TopFolder//'\'//trim(IV%filename)//NoGenstr//'.dat"'
        write(1, '(A)') trim(Output)
        
        !Output = 'move '//newdir//'\'//OutFolder//'\'//trim(IV%filename)//istr//'.resp "'//TopFolder//'\'//trim(IV%filename)//'_'//NoGenstr//'.resp"'
        !write(1, '(A)') trim(Output)
        !Output = 'move '//newdir//'\'//OutFolder//'\'//trim(IV%filename)//istr//'.rsd "'//TopFolder//'\'//trim(IV%filename)//'_'//NoGenstr//'.rsd"'
        !write(1, '(A)') trim(Output)
        !Output = 'move '//newdir//'\'//InFolder//'\'//trim(IV%filename)//istr//'.dat "'//TopFolder//'\'//trim(IV%filename)//'_'//NoGenstr//'.dat"'
        !write(1, '(A)') trim(Output)
        
        close(1)
        deallocate(istr)
        deallocate(NoGenstr)
        
    end subroutine moveTopNestFilesWin
    
end module ReadData