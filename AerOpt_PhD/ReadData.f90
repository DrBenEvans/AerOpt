module ReadData
    
    use InputData
    integer :: ne, np, nbf, nbc                             ! Number of elements, Nodes/Points & boundary faces
    integer, dimension(:,:), allocatable :: connecc         ! Connectivity Matrix of Coarse Mesh    
    integer, dimension(:,:), allocatable :: connecf, boundf ! Connectivity & Boundary Matrix of Fine Mesh    
    real, dimension(:,:), allocatable :: coord              ! Coordinates Matrix of Fine Mesh (includes coordinates of coarse mesh)
    real, dimension(:,:), allocatable :: coarse             ! includes element allocation of nodes to coarse triangles and Area Coefficients of each node
    real, dimension(:,:), allocatable :: Coord_CP           ! desired Coordinates of the Control Points
    real, dimension(:,:), allocatable :: Rect               ! Rectangle definition of 'Influence Box'
    
contains
      
    subroutine SubReadData()
    ! Objective: Reads the Mesh data and Control Points Coordinates
    
        ! Variables
        implicit none
    
        ! Body of ReadData
        open(1, file= InFolder//'/Mesh_fine.txt')
        read(1, 11) ne
        allocate(connecf(ne,IV%NoDim+1))
        read(1, 11) np
        allocate(coord(np,IV%NoDim))
        read(1, 11) nbf
        allocate(boundf(nbf,IV%NoDim))
        do i = 1, ne
            read(1, *) connecf(i,:)
        end do
        do i = 1, np
            read(1, *) coord(i,:)
        end do
        do i = 1, nbf
            read(1, *) boundf(i,:)
        end do 
    11  format(1I8)        
        close(1)
    
        allocate(Coord_CP(IV%NoCP,IV%NoDim))
        open(2, file= InFolder//'/Control_Nodes.txt')
        do i = 1, IV%NoCP
            read(2, *) Coord_CP(i,:)
        end do
        close(2)
    
        allocate(Rect(IV%NoCP,IV%NoDim*4))
        open(3, file= InFolder//'/Rectangles.txt')
        do i = 1, IV%NoCP
            read(3, *) Rect(i,:)
        end do
        close(3)
    
        open(4, file= InFolder//'/Mesh_coarse.txt')
        read(4, *) nbc
        allocate(connecc(nbc,IV%NoDim+1))
        allocate(coarse(np-nbf, IV%NoDim+2))
        do i = 1, (np - nbf)
            read(4, *) coarse(i,:)
        end do
        do i = 1, nbc
            read(4, *) connecc(i,:)
        end do
        close(4)
    
    end subroutine SubReadData
    
    subroutine InitSnapshots(coord_temp, boundff)
    !Objective: Create Outputfile of each Snapshot as Input for the Pre Processor
    
        ! Variables
        implicit none
        real, dimension(np,IV%NoDim) :: coord_temp
        integer, dimension(nbf,(IV%NoDim+1)) :: boundff
    
        ! Body of InitSnapshots
        open(99, file= OutFolder//'/'//trim(IV%filename)//istr//'.dat')
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
            write(99,*) j, coord_temp(j,:)*15.0
        end do
        write(99,*) 'Boundary Faces'
        do j = 1, nbf
            write(99,*) boundff(j,:)
        end do
10      format(2f12.7)        
        close(99)
    
    end subroutine InitSnapshots
    
    subroutine PreProInpFile()
    ! Objective: Create the Inputfile for the PreProcessor in Case of the Windows Executable
    
        ! Variables
        implicit none
    
        ! Body of PreProInpFile
        open(11, file=InFolder//'/PreprocessingInput.txt')        
        write(11,*) OutFolder, '/', trim(IV%filename), istr, '.dat'  ! Name of Input file
        write(11,*) 'f'                                     ! Hybrid Mesh?
        write(11,*) 1                                       ! Number of Grids
        write(11,*) 0                                       ! Directionality Parameters
        write(11,*) 0                                       ! Visualization Modes
        write(11,*) OutFolder, '/', trim(IV%filename), istr, '.sol'  ! Output File name
        close(11)
    
    end subroutine PreProInpFile
    
    subroutine WriteSolverInpFile()
    ! Objective: Create Solver Input files
    
        ! Variables
        implicit none
        character(len=:), allocatable :: fname
        character(len=7) :: strEFMF, strMa      ! String for Ma number & Engine Inlet Front Mass Flow Number
        character(len=2) :: strNI               ! String for Number of Iterations
        integer :: fileLength
    
        ! Body of WriteSolverInpFile
        fileLength = len(IV%filename) + len(istr) + 15
        allocate(character(len=fileLength) :: fname)
        fname = InFolder//'/'//trim(IV%filename)//istr//'.inp'
        write( strEFMF, '(F7.3)' )  IV%engFMF
        write( strMa, '(F7.3)' )  IV%Ma
        write( strNI, '(I2)' )  IV%NoIter
        
        ! Create Input File for Solver
        open(5, file=fname)
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
        open(1, file=InFolder//'/SolverInput'//istr//'.txt')
        if (IV%SystemType == 'B') then
            write(1,*) trim(IV%filename), istr, '.inp'    ! Control Filename
            write(1,*) trim(IV%filename), istr, '.sol'    ! Computation Filename
            write(1,*) ''
            write(1,*) trim(IV%filename), istr, '.resp'  ! Result filename
            write(1,*) trim(IV%filename), istr, '.rsd'   ! Residual Filename
        else
            write(1,*) '2DEngInlSim/', newdir, '/', InFolder, '/', trim(IV%filename), istr, '.inp'    ! Control Filename
            write(1,*) '2DEngInlSim/', newdir, '/', InFolder, '/', trim(IV%filename), istr, '.sol'    ! Computation Filename
            write(1,*) ''
            write(1,*) '2DEngInlSim/', newdir, '/', OutFolder, '/', trim(IV%filename), istr, '.resp'  ! Result filename
            write(1,*) '2DEngInlSim/', newdir, '/', OutFolder, '/', trim(IV%filename), istr, '.rsd'   ! Residual Filename
        end if
        close(1)
    
    end subroutine WriteSolverInpFile
    
    subroutine writeBatchFile()
    
        ! Variables
        implicit none
        character(len=255) :: Output
    
        ! Body of writeBatchFile
        open(1, file= InFolder//'/batchfile'//istr//'.sh')  
        if (IV%SystemType /= 'B') then
            write(1,*) '#PBS -N ' ,trim(IV%filename), istr
            write(1,*) '#PBS -q oh'
            write(1,*) '#PBS -l nodes=1:ppn=1'
            write(1,*) '#PBS -l walltime=24:00:00'
            write(1,*) '#PBS -l mem=1gb'
            write(1,*) '#PBS -m bea'
            write(1,*) '#PBS -M 717761@swansea.ac.uk'
            Output = pathLin_Solver//' < /eng/cvcluster/egnaumann/2DEngInlSim/'//newdir//'/'//InFolder//'/SolverInput'//istr//'.txt'
            write(1,'(A)') trim(Output)
        else
            write(1,*) '#BSUB -J ' ,trim(IV%filename), istr
            write(1,*) '#BSUB -o ' ,trim(IV%filename), istr, '.o'
            write(1,*) '#BSUB -e ' ,trim(IV%filename), istr, '.e'
            !write(1,*) '#BSUB -q <enter a queue>'
            write(1,*) '#BSUB -n 1'
            write(1,*) '#BSUB -W 24:00'
            write(1,*) 'ssh cf-log-001'
            Output = pathLin_Solver//' < /home/david.naumann/2DEngInlSim/'//newdir//'/'//InFolder//'/SolverInput'//istr//'.txt'
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
        open(1, file='FileCreateDir.scr')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%defPath)
        write(1,*) 'cd ..'
        write(1,*) 'chmod 777 2DEngInlSim'
        write(1,*) 'cd 2DEngInlSim'
        write(1,*) 'mkdir ', newdir
        write(1,*) 'cd ', newdir
        write(1,*) 'mkdir ', InFolder
        write(1,*) 'mkdir ', OutFolder
        write(1,*) 'cd ..'
        write(1,*) 'cd ..'
        write(1,*) 'chmod 711 2DEngInlSim'
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
        open(1, file='FileCreateDir.scr')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%defPath) , '/', newdir, '/', InFolder
        
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
    
    subroutine transferFilesLin()
    ! Objective: Move files from Linux in Linux
    
        ! Variables
        implicit none
        character(len=255) :: currentDir
        character(len=255) :: Output
        
        ! Body of createDirectories
        call getcwd(currentDir)
        write(*,*) trim(currentDir)
        open(1, file='FileCreateDir.scr')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%defPath), '/', newdir, '/', InFolder

        Output = 'mv '//trim(currentDir)//'/'//InFolder//'/SolverInput'//istr//'.txt SolverInput'//istr//'.txt'
        write(1, '(A)') trim(Output)
        
        Output = 'mv '//trim(currentDir)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.sol '//trim(IV%filename)//istr//'.sol'
        write(1, '(A)') trim(Output)
        
        Output = 'mv '//trim(currentDir)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.dat '//trim(IV%filename)//istr//'.dat'
        write(1, '(A)') trim(Output)
        
        Output = 'mv '//trim(currentDir)//'/'//InFolder//'/'//trim(IV%filename)//istr//'.inp '//trim(IV%filename)//istr//'.inp'
        write(1, '(A)') trim(Output)
        
        Output = 'mv '//trim(currentDir)//'/'//InFolder//'/batchfile'//istr//'.sh batchfile'//istr//'.sh'
        write(1, '(A)') trim(Output)
        
        close(1)
    
    end subroutine transferFilesLin
    
    subroutine TriggerFile()
    ! Objectives: Triggerfile for Cluster Simulation & Parallelisation
    
        ! Variables
        implicit none
    
        ! Body of TriggerFile
        open(1, file='Trigger.sh')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%defPath), '/', newdir, '/', InFolder
        
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
        open(1, file='Trigger.sh')
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
        character(len=255) :: Output
    
        ! Body of CheckSimStatus
        open(1, file='CheckStatus.scr')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%defPath), '/', newdir, '/', InFolder      
        write(1,*) '[ -e ', trim(IV%filename), istr, '.e* ] && echo 1 > check.txt || echo 0 > check.txt'
        if (IV%SystemType == 'Q')   then
            write(1,*) 'cd ..'
            write(1,*) 'cd ..'
            Output = 'mv /eng/cvcluster/egnaumann/2DEngInlSim/'//newdir//'/'//InFolder//'/check.txt check.txt'
            write(1, '(A)') trim(Output)
        elseif (IV%SystemType == 'B') then
            write(1,*) 'cd ..'
            write(1,*) 'cd ..'
            Output = 'mv /home/david.naumann/2DEngInlSim/'//newdir//'/'//InFolder//'/check.txt check.txt'
            write(1, '(A)') trim(Output)
        end if
        close(1)
    
    end subroutine CheckSimStatus
    
    subroutine CheckSimStatus2()
    ! Objectives: Transfers Check file to correct folder
    
        ! Variables
        implicit none
    
        ! Body of CheckSimStatus
        open(1, file='CheckStatus.scr')
        write(1,*) 'cd '
        write(1,*) 'cd ..'       
        write(1,*) 'cd ', trim(IV%defPath), '/', newdir, '/', InFolder
        write(1,*) 'get check.txt'       
        close(1)
    
    end subroutine CheckSimStatus2
    
    subroutine TransferSolutionOutput()
    
        ! Variables
        character(len=:), allocatable :: Output
        integer :: strOut        
    
        ! Body of TransferSolutionOutput
        open(1, file='FileCreateDir.scr')
        write(1,*) 'cd '
        write(1,*) 'cd ..'
        write(1,*) 'cd ', trim(IV%defPath), '/', newdir, '/', OutFolder
        
        ! Put Commmands to Transfer Data from Output_Data Folder on Linux to Output_Folder on Windows
        strOut = len(trim(IV%filename)) + len(istr) + 9
        allocate(character(len=strOut) :: Output)
        Output = 'get '//trim(IV%filename)//istr//'.resp'
        write(1, '(A)') Output       
        deallocate(Output)
    
    end subroutine TransferSolutionOutput
    
end module ReadData