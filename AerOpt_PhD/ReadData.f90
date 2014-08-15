module ReadData
    
    use InputData
    use Toolbox
    type ReadVariablesData
    
        integer :: ne, np, nbf, nbc                             ! Number of elements, Nodes/Points & boundary faces
        integer, dimension(:,:), allocatable :: connecc         ! Connectivity Matrix of Coarse Mesh    
        integer, dimension(:,:), allocatable :: connecf, boundf ! Connectivity & Boundary Matrix of Fine Mesh    
        real, dimension(:,:), allocatable :: coord              ! Coordinates Matrix of Fine Mesh (includes coordinates of coarse mesh)
        real, dimension(:,:), allocatable :: coarse             ! includes element allocation of nodes to coarse triangles and Area Coefficients of each node
        real, dimension(:,:), allocatable :: Coord_CP           ! desired Coordinates of the Control Points
        real, dimension(:,:), allocatable :: Rect               ! Rectangle definition of 'Influence Box'
        real, dimension(:,:), allocatable :: coord_temp         ! Coordinates Matrix of Fine Mesh
    
    end type ReadVariablesData
    
    type(ReadVariablesData) :: RD
    real, dimension(:,:), allocatable :: ArrayTemp
    
contains
      
    subroutine SubReadData()
    ! Objective: Reads the Mesh data and Control Points Coordinates
    
        ! Variables
        implicit none
        integer :: i
    
        ! Body of ReadData
        open(1, file= InFolder//'/Mesh_fine.txt', form='formatted',status='old')
        read(1, 11) RD%ne
        allocate(RD%connecf(RD%ne,IV%NoDim+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        read(1, 11) RD%np
        allocate(RD%coord(RD%np,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        read(1, 11) RD%nbf
        allocate(RD%boundf(RD%nbf,IV%NoDim + 1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        do i = 1, RD%ne
            read(1, *) RD%connecf(i,:)
        end do
        do i = 1, RD%np
            read(1, *) RD%coord(i,:)
        end do
        do i = 1, RD%nbf
            read(1, *) RD%boundf(i,1:2)
        end do 
    11  format(1I8)        
        close(1)
    
        allocate(RD%Coord_CP(IV%NoCP,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        open(2, file= InFolder//'/Control_Nodes.txt', form='formatted',status='old')
        do i = 1, IV%NoCP
            read(2, *) RD%Coord_CP(i,:)
        end do
        close(2)
    
        allocate(RD%Rect(IV%NoCP,IV%NoDim*4),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        open(3, file= InFolder//'/Rectangles.txt', form='formatted',status='old')
        do i = 1, IV%NoCP
            read(3, *) RD%Rect(i,:)
        end do
        close(3)
    
        open(4, file= InFolder//'/Mesh_coarse.txt', form='formatted',status='old')
        read(4, *) RD%nbc
        allocate(RD%connecc(RD%nbc,IV%NoDim+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        allocate(RD%coarse(RD%np-RD%nbf, IV%NoDim+2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        do i = 1, (RD%np - RD%nbf)
            read(4, *) RD%coarse(i,:)
        end do
        do i = 1, RD%nbc
            read(4, *) RD%connecc(i,:)
        end do
        close(4)
    
    end subroutine SubReadData
    
    subroutine InitSnapshots(ii)
    !Objective: Create Outputfile of each Snapshot as Input for the Pre Processor
    
        ! Variables
        implicit none
        integer :: ii, j
    
        ! Body of InitSnapshots
        
        ! Determine correct String      
        call DetermineStrLen(istr, ii) 
        
        open(99, file= newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'.dat', form='formatted',status='unknown')
        write(99,*) 1
        write(99,*) 'David Naumann'
        write(99,*) 'NoTrgElem NoNodes NoBound'        
        write(99,*) RD%ne, RD%np, RD%nbf
        write(99,*) 'Connectivities'
        do j = 1, RD%ne
            write(99,*) j, RD%connecf(j,:)
        end do
        write(99,*) 'Coordinates'
        do j = 1, RD%np
            write(99,*) j, RD%coord_temp(j,:)*15.0
        end do
        write(99,*) 'Boundary Faces'
        do j = 1, RD%nbf
            write(99,*) RD%boundf(j,:)
        end do
10      format(2f12.7)        
        close(99)
        deallocate (istr)
    
    end subroutine InitSnapshots
    
    subroutine PreProInpFile()
    ! Objective: Create the Inputfile for the PreProcessor in Case of the Windows Executable
    
        ! Variables
        implicit none
    
        ! Body of PreProInpFile
        open(11, file=newdir//'/PreprocessingInput.txt', form='formatted',status='unknown')        
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
        character(len=255) :: Output
        integer :: fileLength
    
        ! Body of WriteSolverInpFile
        fileLength = len(IV%filename) + len(istr) + 15
        allocate(character(len=fileLength) :: fname)
        fname = newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'.inp'
        write( strEFMF, '(F7.3)' )  IV%engFMF
        write( strMa, '(F7.3)' )  IV%Ma
        write( strNI, '(I2)' )  IV%NoIter
        
        ! Create Input File for Solver
        open(5, file=fname, form='formatted',status='unknown')
        write(5,*) '&inputVariables'
        write(5,*) 'ivd%numberOfMGIterations = ' ,strNI , ','
        write(5,*) 'ivd%numberOfGridsToUse = 1,'
        write(5,*) 'ivd%viscosityScheme = 1,'
        write(5,*) 'ivd%boundaryTerm = 1,'
        write(5,*) 'ivd%CFLNumber = 1.0,'
        write(5,*) 'ivd%turbulenceCFLFactor = 1.0,'
        write(5,*) 'ivd%alpha = 180,'
        write(5,*) 'ivd%MachNumber = ' ,strMa, ','
        write(5,*) 'ivd%numberOfRelaxationSteps = 1, '
        write(5,*) 'ivd%ReynoldsNumber = 6500000.0,'
        write(5,*) 'ivd%gamma = 1.4,'
        write(5,*) 'ivd%turbulenceModel = 1,'
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
        write(5,*) 'ivd%useMatrixDissipation = .false.,'
        write(5,*) 'ivd%writeToFileInterval = 100,'
        write(5,*) 'ivd%useTimeResidual = .false.,'
        write(5,*) 'ivd%useDissipationWeighting = .false.,'
        write(5,*) 'ivd%normalComponentRelaxation = 1.0,'
        write(5,*) 'ivd%sizeOfSeparationField = 25,'
        write(5,*) 'ivd%numberOfTriperations = 0,'
        write(5,*) 'ivd%enginesFrontMassFlow = ', strEFMF, ','
        write(5,*) 'ivd%maxitt = 50000,'
        write(5,*) '/'
        close(5)
        
        ! Create Read File for Solver
        open(1, file=newdir//'/'//InFolder//'/SolverInput'//istr//'.sh', form='formatted',status='unknown')
        if (IV%SystemType == 'B') then
            write(1,*) trim(IV%filename), istr, '.inp'    ! Control Filename
            write(1,*) trim(IV%filename), istr, '.sol'    ! Computation Filename
            write(1,*) ''
            write(1,*) trim(IV%filename), istr, '.resp'  ! Result filename
            write(1,*) trim(IV%filename), istr, '.rsd'   ! Residual Filename
        else
            Output = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'.inp'    ! Control Filename
            write(1,'(A)') trim(Output)
            Output = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'.sol'    ! Control Filename
            write(1,'(A)') trim(Output)
            write(1,*) ''
            Output = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.resp'    ! Control Filename
            write(1,'(A)') trim(Output)
            Output = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd'    ! Control Filename
            write(1,'(A)') trim(Output)
        end if
        close(1)
    
    end subroutine WriteSolverInpFile
    
    subroutine writeBatchFile()
    
        ! Variables
        implicit none
        character(len=255) :: Output
    
        ! Body of writeBatchFile
        open(1, file= newdir//'/'//InFolder//'/batchfile'//istr//'.sh', form='formatted',status='unknown')  
        if (IV%SystemType /= 'B') then
            write(1,*) '#PBS -N ' ,trim(IV%filename), istr
            write(1,*) '#PBS -q oh'
            write(1,*) '#PBS -l nodes=1:ppn=1'
            write(1,*) '#PBS -l walltime=24:00:00'
            write(1,*) '#PBS -l mem=1gb'
            write(1,*) '#PBS -m bea'
            write(1,*) '#PBS -M 717761@swansea.ac.uk'
            Output = pathLin_Solver//' < /eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/SolverInput'//istr//'.sh'
            write(1,'(A)') trim(Output)
        else
            write(1,*) '#BSUB -J ' ,trim(IV%filename), istr
            write(1,*) '#BSUB -o ' ,trim(IV%filename), istr, '.o'
            write(1,*) '#BSUB -e ' ,trim(IV%filename), istr, '.e'
            !write(1,*) '#BSUB -q <enter a queue>'
            write(1,*) '#BSUB -n 1'
            write(1,*) '#BSUB -W 24:00'
            Output = pathLin_Solver//' < /home/'//trim(IV%UserName)//'/AerOpt/'//newdir//'/'//InFolder//'/SolverInput'//istr//'.sh'
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
        write(1,*) 'cd ..'
        write(1,*) 'cd ..'
        write(1,*) 'chmod 711 AerOpt'
        close(1)
    
    end subroutine createDirectoriesInit
    
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
        Output = 'put "'//trim(currentDir)//'/'//InFolder//'/SolverInput'//istr//'.txt"'
        write(1, '(A)') Output       
        deallocate(Output)
                      
        strOut = len(trim(currentDir)) + 22
        allocate(character(len=strOut) :: Output)
        Output = 'put "'//trim(currentDir)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.sol"'
        write(1, '(A)') Output
        deallocate(Output)
        
        strOut = len(trim(currentDir)) + 22
        allocate(character(len=strOut) :: Output)
        Output = 'put "'//trim(currentDir)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.dat"'
        write(1, '(A)') Output
        deallocate(Output)

        strOut = len(trim(currentDir)) + 23
        allocate(character(len=strOut) :: Output)
        Output = 'put "'//trim(currentDir)//'/'//InFolder//'/'//trim(IV%filename)//istr//'.inp"'
        write(1, '(A)') Output
        deallocate(Output)
        
        strOut = len(trim(currentDir)) + 25
        allocate(character(len=strOut) :: Output)
        Output = 'put "'//trim(currentDir)//'/'//InFolder//'/batchfile'//istr//'.sh"'
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
            write(1,*) 'qsub batchfile', istr, '.sh'
        else
            write(1,*) 'bsub < batchfile', istr, '.sh'
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
        write(1,*) pathLin_Solver
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
            Output = 'mv /home/'//trim(IV%UserName)//'/AerOpt'//newdir//'/'//InFolder//'/check.txt check.txt'
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
    
end module ReadData