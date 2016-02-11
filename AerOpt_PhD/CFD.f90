module CFD
    
    use InputData
    use ReadData
    use Toolbox
    use CreateSnapshots
    use GenerateMesh
    use FDGD
    
    contains
    
    subroutine SubCFD(Start, Ending, CN_CoordinatesArray, sizing)
    
        ! Variables
        implicit none
        integer :: Start, Ending, i, sizing
        double precision, dimension(sizing, maxDoF) :: CN_CoordinatesArray

        ! Body of SubCFD
        ! ****Generate Snapshots (initial Meshes)**** !
        allocate(RD%coord_temp(RD%np,IV%nodim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main "
        do i = Start, Ending          
            print *, "Generating Mesh", i, "/", Ending
            RD%coord_temp = RD%coord 
            call SubMovemesh(CN_CoordinatesArray(i - Start + 1,:))
            !Output: new coordinates - Mesh with moved boundaries based on Initial Nest
            
    !!!!! IMPLEMENT Mesh Quality Test
    
            ! Write Snapshot to File
            if (IV%NoDim == 2) then
                call writeDatFile(i)
            elseif (IV%NoDim == 3) then
                call writepltFile(i)
                call writebcoFile(i)
            end if
        end do
        deallocate(RD%coord_temp)
 
        if (IV%Meshtest == .true.) then
          pause
        end if
       
        ! ****Call 2D Preprocessor and pass on input parameters**** !
        print *, 'Start Preprocessing'
        do i = Start, Ending  
            call PreProcessing(i)  
        end do
        print *, 'Finished Preprocessing'
    
        ! ****Call 2D FLITE Solver and pass on input parameters**** !
        print *, 'Call FLITE 2D Solver'
        do i = Start, Ending        
           call Solver(i)               
        end do
        print *, 'Finished Submitting Jobs to FLITE 2D Solver' 
        
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
                strSystem = 'move '//trim(IV%filename)//istr//'.unk "'//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk"'
                call system(trim(strSystem))
                strSystem = 'move '//trim(IV%filename)//istr//'.rsd "'//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd"'
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
        call CheckforConvergence(i, InitConv, NoFiles)
        print*, 'All Solutions converged'
        
        ! Delete the Error files to allow sleep check in next generation
        do i = 1, NoFiles
            call SubDeleteErrorFiles(i)
        end do
        
    end subroutine PostSolverCheck
    
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
            strSystem = pathPrePro//' < '//newdir//'\'//InFolder//'/PreprocessingInput.txt >nul 2>&1'
            
        else
            
            ! write command (for Linux)
            allocate(character(len=100) :: strSystem)
            strSystem = pathPrepro//' < '//newdir//'/'//InFolder//'/PreprocessingInput.txt > /dev/null'
            
        end if
        print *, 'Preprocessing Geometry', i
        print *, ' '
        call system(trim(strSystem))   ! System operating command called to activate fortran 
        if (IV%NoDim == 3) then
            if (IV%SystemType == 'W') then
                call system('move '//trim(IV%filepath)//'/plotreg.reg "'//trim(IV%filepath)//'/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'/plotreg.reg"')
                call system('move '//trim(IV%filepath)//'/base.plt "'//trim(IV%filepath)//'/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'/base.plt"')
            else
                call system('mv '//trim(IV%filepath)//'/plotreg.reg "'//trim(IV%filepath)//'/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'/plotreg.reg"')
                call system('mv '//trim(IV%filepath)//'/base.plt "'//trim(IV%filepath)//'/'//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'/base.plt"')
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
            call writeBatchFile()
        elseif (IV%NoDim == 3) then
            call writeBatchFile_3D()
        end if
    
        ! Is AerOpt executed from Linux or Windows?                
        if (IV%SystemType == 'W')   then    ! AerOpt is executed from a Windows machine
            
            if (IV%runOnCluster == 'Y') then
                ! Transfer Files from Windows Machine onto Cluster
                call transferFilesWin()            
                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'psftp')
                call Triggerfile()           ! Triggerfile for submission
                ! Submits Batchfile via Putty
                call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'putty')
            else
                print *, 'Solving Geometry', i
                allocate(character(len=200) :: strSystem)
                strSystem = pathSolver//' < '//newdir//'/'//InFolder//'\'//trim(IV%filename)//istr//'/SolverInput'//istr//'.sh >nul 2>&1'
                call system(trim(strSystem))
                strSystem = 'move '//newdir//'\'//InFolder//'\'//trim(IV%filename)//istr//'\'//trim(IV%filename)//istr//'.unk "'//newdir//'\'//OutFolder//'\'//trim(IV%filename)//istr//'.unk"'
                call system(trim(strSystem))
                strSystem = 'move '//newdir//'\'//InFolder//'\'//trim(IV%filename)//istr//'\'//trim(IV%filename)//istr//'.rsd "'//newdir//'\'//OutFolder//'\'//trim(IV%filename)//istr//'.rsd"'
                call system(trim(strSystem))
                deallocate (strSystem)
            end if 
              
        elseif (IV%SystemType == 'L')   then    ! AerOpt is executed on a Linux machine
            
		    print *, 'Solving Geometry', i
            allocate(character(len=200) :: strSystem)
            strSystem = pathSolver//' < '//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'/SolverInput'//istr//'.sh > /dev/null'
            call system(trim(strSystem))
            strSystem = 'mv '//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.unk "'//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk"'
            call system(trim(strSystem))
            strSystem = 'mv '//newdir//'/'//InFolder//'/'//trim(IV%filename)//istr//'/'//trim(IV%filename)//istr//'.rsd "'//newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd"'
            call system(trim(strSystem))
            deallocate (strSystem)
            
        else    ! AerOpt is executed on a Linux cluster
            
            if (IV%runOnCluster == 'Y') then
                call Triggerfile()     ! Triggerfile for submission
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
        
            ! Check Status of Simulation by checking the existence of all error files

            do i = 1, NoFiles

                ! Determine correct String      
                call DetermineStrLen(istr, i)
                ! Creates File containing Linux commands to check for last file
                call CheckSimStatus()
                ! Submit File
                if (IV%SystemType == 'W')   then
                    call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'plink')
                else
                    call system('chmod a+x ./Communication')
                    call system('./Communication')
                end if
                ! Creates File to transfer response from Windows to Linux
                if (IV%SystemType == 'W')   then
                    call CheckSimStatus2()
                    call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'psftp')
                end if
            
                open(1, file='check.txt',form='formatted',status='old')
                read(1,*) jobcheck
                close(1)
                deallocate(istr)
                if (jobcheck == 0) EXIT           
            
            end do
        
            waitTime = (IV%delay/3600.0) + waitTime
            if (waitTime > IV%waitMax) then
                STOP 'Cluster Simulation Time exceeded maximum waiting Time'
            end if
        
        end do
        print*, 'End Sleep - Jobs are finished'
        
    end subroutine Sleep
    
    recursive subroutine CheckforConvergence(Iter, InitConv, NoFiles)
    
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
        
        if (NoConv /= 0) then
            
            !!** Re-Do diverged solutions **!!
            do ii = 1, NoConv
                
                if (InitConv == 0) then
                    call SubCFD(DivNestPos(ii), DivNestPos(ii), CS%Snapshots(DivNestPos(ii),:), 1)
                else
                    call SubCFD(DivNestPos(ii), DivNestPos(ii), OV%Nests(DivNestPos(ii),:), 1)    
                end if
                
                ! Delete the Error files to allow sleep check
                call SubDeleteErrorFiles(NoFiles - IV%NoSnap + DivNestPos(ii))
                
            end do
            
            call Sleep(NoFiles)
            
            Iter = Iter + 1
            
            if (Iter < 4) then
                call CheckforConvergence(Iter, InitConv, NoFiles)
            end if
            
            if (NoConv /= 0) then
                STOP 'Convergence of CFD Simulations could not be achieved. Check initial Mesh and/or Movement constraints!'
            end if
        end if
        
    end subroutine CheckforConvergence
    
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
        open(1, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', form='formatted', STATUS="OLD")     
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
            
        ! Convergence = false if last line contains 'NaN'
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
    
    subroutine SubDeleteErrorFiles(i)
    
        ! Variables
        integer :: i
        
        ! Body of SubDeleteErrorFiles
        ! Determine correct String      
        call DetermineStrLen(istr, i)
        
        call DeleteErrorFiles(istr)
        if (IV%SystemType == 'W' .and. IV%runOnCluster == 'Y')   then    ! AerOpt is executed from a Windows machine           
            call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Communication', 'psftp')
        elseif (IV%SystemType /= 'W')  then                
            call system('chmod a+x ./Communication')
            call system('./Communication')
        end if
                
        deallocate(istr)
    
    end subroutine SubDeleteErrorFiles
    
end module CFD