module CFD
    
    use InputData
    use ReadData
    use Toolbox
    use CreateSnapshots
    use GenerateMesh
    use FDGD
    
contains
    
    subroutine SubCFD(vecIndex, CN_CoordinatesArray, sizing)
    
        ! Variables
        implicit none
        integer :: i, sizing
        integer, dimension(sizing) :: vecIndex
        double precision, dimension(sizing, maxDoF) :: CN_CoordinatesArray

        ! Body of SubCFD
      
        ! ****Generate Meshes**** !
        allocate(RD%coord_temp(RD%np,IV%nodim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main "
        i = 1
        do while (i .le. size(vecIndex, dim = 1))
            print *, "Generating Mesh", vecIndex(i), "/", vecIndex(sizing)
            RD%coord_temp = RD%coord 
            call SubMovemesh(CN_CoordinatesArray(i,:))
            !Output: new coordinates - Mesh with moved boundaries based on Initial Nest
            
    !!!!! IMPLEMENT Mesh Quality Test
    
            ! Write Snapshot to File
            if (IV%NoDim == 2) then
                call writeDatFile(vecIndex(i))
            elseif (IV%NoDim == 3) then
                call writepltFile(vecIndex(i))
                call writebcoFile(vecIndex(i))
            end if
        i = i + 1
        end do
        deallocate(RD%coord_temp)
 
        if (IV%Meshtest == .true.) then
          pause
        end if
       
        ! ****call 2d preprocessor and pass on input parameters**** !
        print *, 'start preprocessing'
        i = 1
        do while (i .le. size(vecIndex, dim = 1)) 
            call preprocessing(vecIndex(i))
        i = i + 1
        end do
        print *, 'finished preprocessing'
        
        ! ****call 2d flite solver and pass on input parameters**** !
        print *, 'call flite 2d solver'
        if (IV%SystemType == 'Q') then
            i = 1
            do while (i .le. size(vecIndex, dim = 1))       
                call solver(vecIndex(i))
            i = i + 1
            end do
        else
            i = 1
            do while (i .le. size(vecIndex, dim = 1))
                ! Determine correct String      
                call DetermineStrLen(istr, vecIndex(i))
                call WriteSolverInpFile()
                deallocate(istr)
                i = i + 1
            end do
            call writeBatchFile(vecIndex(1), vecIndex(sizing))
            call Triggerfile()     ! Triggerfile for submission
            call system('chmod a+x ./Communication')
            call system('./Communication')
        end if       
        print *, 'finished submitting jobs to flite 2d solver' 
        
    end subroutine SubCFD
    
    subroutine PostSolverCheck(NoFiles, InitConv)
    
        ! Variables
        implicit none
        integer :: NoFiles, i, InitConv
  
        ! Body of PostSolverCheck
        ! ****Wait & Check for FLITE Solver Output**** !
        if (IV%runOnCluster == 'Y') then
            call Sleep(NoFiles)
        end if
        
        if (IV%NoDim == 3) then
            call generateUnkFile()
        end if
        
        if (IV%SystemType == 'W' .and. IV%runOnCluster == 'Y') then
            allocate(character(len=200) :: strSystem)
            do i = 1, NoFiles
                call DetermineStrLen(istr, i)
                call TransferSolutionOutput()
                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')   ! Submits transfersolution Output file                              
                strSystem = 'move '//trim(IV%filename)//istr//'.unk "'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk"'
                call system(trim(strSystem))
                strSystem = 'move '//trim(IV%filename)//istr//'.rsd "'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd"'
                call system(trim(strSystem))                
                deallocate(istr)
            end do
            deallocate (strSystem)
        end if
    
        ! ****Check Simulation Results**** !
        print *, ''
        print *, '*************************************'
        print *, '***  Start Check for Convergence  ***'
        print *, '*************************************'
        print *, ''
        call CheckforConvergence(NoFiles)
        print*, 'All Solutions checked for convergence'
        
        call SubDeleteFiles()
        
    end subroutine PostSolverCheck
    
    subroutine PostSolverCheckInit(NoFiles, InitConv)
    
        ! Variables
        implicit none
        integer :: NoFiles, i, InitConv
  
        ! Body of PostSolverCheck
        ! ****Wait & Check for FLITE Solver Output**** !
        if (IV%runOnCluster == 'Y') then
            call Sleep(NoFiles)
        end if
        
        if (IV%NoDim == 3) then
            call generateUnkFile()
        end if
        
        if (IV%SystemType == 'W' .and. IV%runOnCluster == 'Y') then
            allocate(character(len=200) :: strSystem)
            do i = 1, NoFiles
                call DetermineStrLen(istr, i)
                call TransferSolutionOutput()
                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')   ! Submits transfersolution Output file                              
                strSystem = 'move '//trim(IV%filename)//istr//'.unk "'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk"'
                call system(trim(strSystem))
                strSystem = 'move '//trim(IV%filename)//istr//'.rsd "'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd"'
                call system(trim(strSystem))                
                deallocate(istr)
            end do
            deallocate (strSystem)
        end if
    
        ! ****Check Simulation Results**** !
        print *, ''
        print *, '*************************************'
        print *, '***  Start Check for Convergence  ***'
        print *, '*************************************'
        print *, ''
        i = 0
        call CheckforConvergenceInit(i, InitConv, NoFiles)
        print*, 'All Solutions converged'
        
        call SubDeleteFiles()
        
    end subroutine PostSolverCheckInit
    
    subroutine PreProcessing(i)
    
        ! Variables
        implicit none
        integer :: i
    
        ! Body of PreProcessing
    
        ! Determine correct String      
        call DetermineStrLen(istr, i)
    
        ! write Inputfile
        if (IV%NoDim == 2) then
            call PreProInpFile()
        elseif (IV%NoDim == 3) then
            call PreProInpFile_3D()
        end if
        
        if (IV%SystemType == 'W') then
! Parallelise by sending as Job?             
            allocate(character(len=100) :: strSystem)
            strSystem = pathPrePro//' < '//trim(IV%SimulationName)//'\'//InFolder//'/PreprocessingInput.txt >nul 2>&1'
            
        else
            
            ! write command (for Linux)
            allocate(character(len=100) :: strSystem)
            strSystem = pathPrepro//' < '//trim(IV%SimulationName)//'/'//InFolder//'/PreprocessingInput.txt > /dev/null'
            
        end if
        print *, 'Preprocessing Geometry', i
        print *, ' '
        call system(trim(strSystem))   ! System operating command called to activate fortran 
        if (IV%NoDim == 3) then
            if (IV%SystemType == 'W') then
                call system('move '//trim(IV%filepath)//'/plotreg.reg "'//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/plotreg.reg"')
                call system('move '//trim(IV%filepath)//'/base.plt "'//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/base.plt"')
            else
                call system('mv '//trim(IV%filepath)//'/plotreg.reg "'//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/plotreg.reg"')
                call system('mv '//trim(IV%filepath)//'/base.plt "'//trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/base.plt"')
            end if
        end if
        deallocate (istr)
        deallocate (strSystem)
    
    end subroutine PreProcessing
    
    subroutine Solver(i)
    
        ! Variables
        implicit none
        integer :: i
    
        ! Body of Solver
        
        ! Determine correct String      
        call DetermineStrLen(istr, i)
        
        ! Creates the input file including Solver Parameters and a second file including I/O filenames
        if (IV%NoDim == 2) then
            call WriteSolverInpFile()
        elseif (IV%NoDim == 3) then
            call WriteSolverInpFile_3D()
        end if
        
        ! writes the batchfile to execute Solver on Cluster
        if (IV%NoDim == 2) then
            call writeBatchFile(i,i)
        elseif (IV%NoDim == 3) then
            call writeBatchFile_3D()
        end if
    
        ! Is AerOpt executed from Linux or Windows?                
        if (IV%SystemType == 'W')   then    ! AerOpt is executed from a Windows machine
            
            if (IV%runOnCluster == 'Y') then
                ! Transfer Files from Windows Machine onto Cluster
                call transferFilesWin()            
                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'psftp')
                call TriggerFileQ()           ! Triggerfile for submission
                ! Submits Batchfile via Putty
                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'putty')
            else
                print *, 'Solving Geometry', i
                allocate(character(len=200) :: strSystem)
                strSystem = pathSolver//' < '//trim(IV%SimulationName)//'/'//InFolder//'\'//trim(IV%filename)//istr//'/SolverInput'//istr//'.sh >nul 2>&1'
                call system(trim(strSystem))
                strSystem = 'move '//trim(IV%SimulationName)//'\'//InFolder//'\'//trim(IV%filename)//istr//'\'//trim(IV%filename)//istr//'.unk "'//trim(IV%SimulationName)//'\'//OutFolder//'\'//trim(IV%filename)//istr//'.unk"'
                call system(trim(strSystem))
                strSystem = 'move '//trim(IV%SimulationName)//'\'//InFolder//'\'//trim(IV%filename)//istr//'\'//trim(IV%filename)//istr//'.rsd "'//trim(IV%SimulationName)//'\'//OutFolder//'\'//trim(IV%filename)//istr//'.rsd"'
                call system(trim(strSystem))
                deallocate (strSystem)
            end if 
              
        elseif (IV%SystemType == 'L')   then    ! AerOpt is executed on a Linux machine
            
		    print *, 'Solving Geometry', i
            allocate(character(len=200) :: strSystem)
            strSystem = pathSolver//' < '//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/SolverInput'//istr//'.sh > /dev/null'
            call system(trim(strSystem))
            strSystem = 'mv '//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.unk "'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk"'
            call system(trim(strSystem))
            strSystem = 'mv '//trim(IV%SimulationName)//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.rsd "'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd"'
            call system(trim(strSystem))
            deallocate (strSystem)
            
        else    ! AerOpt is executed on a Linux cluster
            
            if (IV%runOnCluster == 'Y') then
                call TriggerFileQ()     ! Triggerfile for submission
            else
                call TriggerFile2()               ! Triggerfile for submission
            end if
            
            ! Submits Batchfile
            call system('chmod a+x ./Communication')
            call system('./Communication')
                    
        end if
        deallocate(istr)
    
    end subroutine Solver
    
    subroutine Sleep(NoFiles)
    
        ! Variables
        implicit none
        integer :: i, NoFiles, j
        integer, dimension(13) :: fileinfo
        integer, dimension(9) :: timeend
    
        ! Body of Sleep
        print*, 'Start Sleep'
        jobcheck = 0
        waitTime = 0
        j = 1
        do while (jobcheck==0)
        
            ! Wait Function
            print*, 'Sleep', IV%Ma
            call SleepQQ(IV%delay*1000)
            print*, 'Wake Up - Check ', j
            j = j + 1
          
            ! Check Status of Simulations by checking if new output files have been moved
            do i = 1, NoFiles
                
                ! Extract last modification time to check if file has been newly moved
                fileinfo = 0
                call DetermineStrLen(istr, i)
                call stat(trim(IV%filepath)//'/'//trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', fileinfo)
                call ltime(fileinfo(10), timeend)
                deallocate(istr)
                
                ! Check time of files generated
                jobcheck = 0
                    if (OV%timestart(3) == timeend(4)) then
                        if (OV%timestart(5) == timeend(3)) then  
                            if (OV%timestart(6) == timeend(2)) then
                                if (OV%timestart(7) == timeend(1)) then
                                    jobcheck = 1
                                elseif (OV%timestart(7) < timeend(1)) then
                                    jobcheck = 1
                                end if
                            elseif (OV%timestart(6) < timeend(2)) then
                                jobcheck = 1
                            end if
                        elseif (OV%timestart(5) < timeend(3)) then
                            jobcheck = 1
                        end if
                    elseif (OV%timestart(3) < timeend(4)) then
                        jobcheck = 1
                    end if
                if (jobcheck == 0) EXIT                           
            
            end do
        
            waitTime = (IV%delay/3600.0) + waitTime
            if (waitTime > IV%waitMax) then
                STOP 'Cluster Simulation Time exceeded maximum waiting Time'
            end if
        
        end do
        print*, 'End Sleep - Jobs are finished'
        
    end subroutine Sleep
    
    subroutine CheckforConvergence(NoFiles)
    
        ! Variables
        implicit none
        integer :: i, NoFiles
        logical :: Converge
    
        ! Body of CheckforConvergence
        OV%converged = 1
        Converge = .true.
        do i = 1, NoFiles

            call FileCheckConvergence(Converge, i)  
          
            ! All diverged Snapshots are pulled halfway to midpoint(no movement center)
            if (Converge == .false.) then              
                    print *, 'File', i, 'failed to converge and will be set to -150'
                    OV%converged(i) = 0
            end if
            Converge = .true.
            
        end do
        
    end subroutine CheckforConvergence
    
    recursive subroutine CheckforConvergenceInit(Iter, InitConv, NoFiles)
    
        ! Variables
        implicit none
        integer :: ii, InitConv, i, NoFiles
        integer, save :: NoConv
        integer,intent(inout) :: Iter
        logical :: Converge
        character(len=200) :: strCommand
        integer, dimension(:), allocatable :: DivNestPos, tempArray
        double precision, dimension(:), allocatable :: MidPoints
        
        allocate(DivNestPos(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CheckForConvergence "
        allocate(MidPoints(maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CheckforConvergence "
    
        ! Body of CheckforConvergence
        print *, 'Iteration', (Iter + 1)
        call timestamp()
		OV%converged = 1
        Converge = .true.
        NoConv = 0 
        do i = 1, NoFiles
    
            call FileCheckConvergence(Converge, i)  
          
            ! All diverged Snapshots are pulled halfway to midpoint(no movement center)
            if (Converge == .false.) then              
                    print *, 'File', i, 'failed to converge and will be resimulated'
                    NoConv = NoConv + 1
                    DivNestPos(NoConv) = i
                    MidPoints = CS%MxDisp(:,1) - (CS%MxDisp(:,1) - CS%MxDisp(:,2))/2.0  ! Midpoint calculation
                    if (InitConv == 0) then
                        CS%Snapshots(i,:) = CS%Snapshots(i,:) - ((CS%Snapshots(i,:) - MidPoints)/2.0)   ! Half way between current Nest and Midpoint
                    else
                        OV%Nests_Move(i,:) = OV%Nests_Move(i,:) - ((OV%Nests_Move(i,:) - MidPoints)/2.0) ! Half way between current Nest and Midpoint
                        OV%Nests(i,:) = OV%Nests(i,:) - ((OV%Nests(i,:) - MidPoints)/2.0)
                    end if
            end if
            Converge = .true.
            
        end do
        
        ! Resize DivNestPosArray
        allocate(tempArray(NoConv),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CheckForConvergence "
        tempArray = DivNestPos(1:NoConv)
        deallocate(DivNestPos)
        allocate(DivNestPos(NoConv),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CheckForConvergence "        
        DivNestPos = tempArray
        deallocate(tempArray)
        
        Iter = Iter + 1
        
        if (NoConv /= 0 .and. Iter < 3) then
            
            !!** Re-Do diverged solutions **!!                
            if (InitConv == 0) then
                call SubCFD(DivNestPos, CS%Snapshots(DivNestPos,:), NoConv)
            else
                call SubCFD(DivNestPos, OV%Nests(DivNestPos,:), NoConv)    
            end if
        
            call Sleep(NoFiles)
            
            call SubDeleteFiles()
            
            call CheckforConvergenceInit(Iter, InitConv, NoFiles)
            
        end if
        
        if (NoConv /= 0 .and. Iter == 3) then
            OV%converged(DivNestPos) = 0
            print *, 'Convergence of CFD Simulations could not be achieved. Diverged solutions will be excluded!'
        end if
        
    end subroutine CheckforConvergenceInit
    
    subroutine FileCheckConvergence(Converge, NoFile)
    
        ! Variables
        implicit none
        integer :: FileSize, LastLine, NoFile, j
        double precision, dimension(8) :: Input
        logical, intent(in out) :: Converge
    
        ! Body of FileCheckConvergence
        ! Determine correct String number
        call DetermineStrLen(istr, NoFile)
            
        ! Open .rsd file to check, if the last line contains 'Nan' solutions, which would mean convergence fail
        open(1, file=trim(IV%SimulationName)//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', form='formatted', STATUS="OLD")     
        inquire(1, size = FileSize) 
        if (IV%NoDim == 3) then
            LastLine = FileSize/175
        else
            if (IV%SystemType == 'W') then
                LastLine = FileSize/107
            else     
                LastLine = FileSize/106
            end if
        end if
        
        ! Read until last line
        do j = 1, (LastLine - 1)
            read(1, *) Input
        end do
        read(1, *) Input
        close(1)
            
        ! Convergence = false if last line contains 'NaN' or is below convergence criteria
        do j = 1, 8
            if (isnan(Input(j))) then
                Converge = .false.
                exit   
            end if       
        end do
        if (Input(2) > IV%NoIter) then
            Converge = .false.
        end if
                
        deallocate(istr)
            
    end subroutine FileCheckConvergence
    
    subroutine SubDeleteFiles()
    
        ! Variables
        implicit none
        
        ! Body of SubDeleteErrorFiles
        
        if (IV%SystemType == 'B') then
            call DeleteLogFiles()
        else
            call DeleteErrorFiles()
        end if
        if (IV%SystemType == 'W' .and. IV%runOnCluster == 'Y')   then    ! AerOpt is executed from a Windows machine           
            call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'psftp')
        elseif (IV%SystemType /= 'W')  then                
            call system('chmod a+x ./Communication')
            call system('./Communication')
        end if
    
    end subroutine SubDeleteFiles
    
    
    
    
    !subroutine Sleep(NoFiles)
    !
    !    ! Variables
    !    implicit none
    !    integer :: i, NoFiles, j
    !
    !    ! Body of Sleep
    !    print*, 'Start Sleep'
    !    jobcheck = 0
    !    waitTime = 0
    !    j = 1
    !    do while (jobcheck==0)
    !    
    !        ! Wait Function
    !        print*, 'Sleep', IV%Ma
    !        call SleepQQ(IV%delay*1000)
    !        print*, 'Wake Up - Check ', j
    !        j = j + 1
    !    
    !        ! Check Status of Simulation by checking the existence of all error files
    !
    !        do i = 1, NoFiles
    !
    !            ! Determine correct String      
    !            call DetermineStrLen(istr, i)
    !            ! Creates File containing Linux commands to check for last file
    !            call CheckSimStatus()
    !            ! Submit File
    !            if (IV%SystemType == 'W')   then
    !                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'plink')
    !            else
    !                call system('chmod a+x ./Communication')
    !                call system('./Communication')
    !            end if
    !            ! Creates File to transfer response from Windows to Linux
    !            if (IV%SystemType == 'W')   then
    !                call CheckSimStatus2()
    !                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'psftp')
    !            end if
    !        
    !            open(1, file='check.txt',form='formatted',status='old')
    !            read(1,*) jobcheck
    !            close(1)
    !            deallocate(istr)
    !            if (jobcheck == 0) EXIT           
    !        
    !        end do
    !    
    !        waitTime = (IV%delay/3600.0) + waitTime
    !        if (waitTime > IV%waitMax) then
    !            STOP 'Cluster Simulation Time exceeded maximum waiting Time'
    !        end if
    !    
    !    end do
    !    print*, 'End Sleep - Jobs are finished'
    !    
    !end subroutine Sleep
    !
end module CFD