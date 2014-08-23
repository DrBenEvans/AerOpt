module CFD
    
    use InputData
    use ReadData
    use Toolbox
    use CreateSnapshots
    use GenerateMesh
    use FDGD
    
    contains
    
    subroutine SubCFD(Start, Ending, CN_Coordinates, size)
    
        ! Variables
        implicit none
        integer :: Start, Ending, i, size
        double precision, dimension(size, maxDoF) :: CN_Coordinates
        
    
        ! Body of SubCFD
        ! ****Generate Snapshots (initial Meshes)**** !
        allocate(RD%coord_temp(RD%np,IV%nodim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main "    
        do i = Start, Ending
            
            print *, "Generating Mesh", i, "/", Ending
            RD%coord_temp = RD%coord 
            if (IV%MeshGeneration == 'FDGD') then
                call SubFDGD(dble(CN_Coordinates(i - Start + 1,:)))
            else if (IV%MeshGeneration == 'RBF') then
                call SubGenerateMesh(CN_Coordinates(i - Start + 1,:))
            end if
            ! Output: new coordinates - Mesh with moved boundaries based on Initial Nest
  
    !!!!! IMPLEMENT Mesh Quality Test

            ! Write Snapshot to File
            call writeDatFile(i)

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
    
    subroutine PostSolverCheck(m, InitConv, SolutionNumber, Nests_Move, Nests)
    
        ! Variables
        implicit none
        integer :: m, ii, InitConv
        integer, optional :: SolutionNumber
        double precision, dimension(IV%NoNests,DoF), optional :: Nests_Move
        double precision, dimension(IV%NoNests,maxDoF), optional :: Nests
    
        ! Body of PostSolverCheck
        ! ****Wait & Check for FLITE Solver Output**** !
        call Sleep(m)
    
        if (IV%SystemType == 'W') then
            call TransferSolutionOutput()
        end if
    
        ! ****Check Simulation Results**** !
        print *, ''
        print *, '*************************************'
        print *, '***  Start Check for Convergence  ***'
        print *, '*************************************'
        print *, ''
        ii = 0
        call CheckforConvergence(ii, InitConv, SolutionNumber, Nests_Move, Nests)
        print*, 'All Solutions converged'
        
    end subroutine PostSolverCheck
    
    subroutine PreProcessing(i)
    
        ! Variables
        implicit none
        integer :: i
    
        ! Body of PreProcessing
    
        ! Determine correct String      
        call DetermineStrLen(istr, i)
    
        ! write Inputfile
        call PreProInpFile()
        
        if (IV%SystemType == 'W') then
             
            allocate(character(len=29) :: strSystem)
            strSystem = pathWin
            
        else
            
            ! write command (for Linux)
            allocate(character(len=100) :: strSystem)
            strSystem = pathLin_Prepro//' < '//newdir//'/PreprocessingInput.txt'//' > /dev/null'
            
        end if
        print *, 'Preprocessing Snapshot', i
        print *, ' '
        call system(trim(strSystem))   ! System operating command called to activate fortran       
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
        call WriteSolverInpFile()
        
        ! writes the batchfile to execute Solver on Cluster
        call writeBatchFile()
    
        ! Is AerOpt executed from Linux or Windows?                
        if (IV%SystemType == 'W')   then    ! AerOpt is executed from a Windows machine
            
            ! Transfer Files from Windows Machine onto Cluster
            call transferFilesWin()            
            call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')
            
            if (IV%runOnCluster == 'Y') then
                call Triggerfile()           ! Triggerfile for submission
            else
                call TriggerFile2()  ! Triggerfile for submission
            end if
                    
            ! Submits Batchfile via Putty
            call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'Trigger.sh', 'putty')
                
        else    ! AerOpt is executed from a Linux machine
            
            if (IV%runOnCluster == 'Y') then
                call Triggerfile()     ! Triggerfile for submission
            else
                call TriggerFile2()               ! Triggerfile for submission
            end if
            
            ! Submits Batchfile
            call system('chmod a+x ./Trigger.sh')
            call system('./Trigger.sh')
                    
        end if
                
        deallocate(istr)
    
    end subroutine Solver
    
    subroutine Sleep(m)
    
        ! Variables
        implicit none
        integer :: i, m, j
    
        ! Body of Sleep
        print*, 'Start Sleep'
        jobcheck = 0
        waitTime = 0
        j = 1
        do while (jobcheck==0)
        
            ! Wait Function
            print*, 'Sleep'
            call SleepQQ(IV%delay*1000)
            print*, 'Wake Up - Check ', j
            j = j + 1
        
            ! Check Status of Simulation by checking the existence of all error files

            do i = m, 1, -1

                ! Determine correct String      
                call DetermineStrLen(istr, i)
                ! Creates File containing Linux commands to check for last file
                call CheckSimStatus()
                ! Submit File
                if (IV%SystemType == 'W')   then
                    call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'CheckStatus.scr', 'plink')
                else
                    call system('chmod a+x ./CheckStatus.scr')
                    call system('./CheckStatus.scr')
                end if
                ! Creates File to transfer response from Windows to Linux
                if (IV%SystemType == 'W')   then
                    call CheckSimStatus2()
                    call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'CheckStatus.scr', 'psftp')
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
    
    recursive subroutine CheckforConvergence(Iter, InitConv, SolutionNumber, Nests_Move, Nests)
    
        ! Variables
        implicit none
        integer :: ii, InitConv, i
        integer, save :: NoConv
        integer,intent(inout) :: Iter
        logical :: Converge
        character(len=200) :: strCommand
        integer, dimension(:), allocatable :: DivNestPos
        double precision, dimension(:), allocatable :: MidPoints
        integer, optional :: SolutionNumber      
        double precision, dimension(IV%NoNests,DoF), optional :: Nests_Move
        double precision, dimension(IV%NoNests,maxDoF), optional :: Nests
        
        allocate(DivNestPos(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CheckForConvergence "
        allocate(MidPoints(IV%NoCP*IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CheckforConvergence "
    
        ! Body of CheckforConvergence
        print *, 'Iteration', (Iter + 1)
        Converge = .true.
        NoConv = 0
        if (InitConv == 0) then           
            do i = 1, IV%NoSnap
            
                call FileCheckConvergence(Converge, i)  
            
                ! All diverged Snapshots are pulled halfway to midpoint(no movement center)
                if (Converge == .false.) then
                
                        print *, 'File', i, 'failed to converge and will be resimulated'
                        NoConv = NoConv + 1
                        DivNestPos(NoConv) = i
                        MidPoints = MxDisp(:,1) - (MxDisp(:,1) - MxDisp(:,2))/2.0  ! Midpoint calculation                 
                        Snapshots(i,:) = Snapshots(i,:) - ((Snapshots(i,:) - MidPoints)/2.0)   ! Half way between current Nest and Midpoint
                    
                end if
                Converge = .true.
            
            end do          
        else           
            do i = SolutionNumber, (SolutionNumber - IV%NoNests), -1
            
                call FileCheckConvergence(Converge, i)  
            
                ! All diverged Nests are pulled halfway to midpoint(no movement center)
                if (Converge == .false.) then
                
                        print *, 'File', i, 'failed to converge and will be resimulated'
                        NoConv = NoConv + 1
                        DivNestPos(NoConv) = i
                        MidPoints = MxDisp(:,1) - (MxDisp(:,1) - MxDisp(:,2))/2.0  ! Midpoint calculation
                        Nests_Move(i,:) = Nests_Move(i,:) - ((Nests_Move(i,:) - MidPoints)/2.0) ! Half way between current Nest and Midpoint
                        Nests(i,:) = Nests(i,:) - ((Nests(i,:) - MidPoints)/2.0)
                                  
                    
                end if
                Converge = .true.
            
            end do            
        end if
        
        if (NoConv /= 0) then
            
            !!** Re-Do diverged solutions **!!
            do ii = 1, NoConv
                
                call SubCFD(DivNestPos(ii), DivNestPos(ii), Snapshots(DivNestPos(ii),:), 1)
                
                ! Determine correct String      
                call DetermineStrLen(istr, DivNestPos(ii))
        
                call DeleteErrorFiles(istr)
                if (IV%SystemType == 'W')   then    ! AerOpt is executed from a Windows machine           
                    call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')
                else
                    call system('chmod a+x ./FileCreateDir.scr')
                    call system('./FileCreateDir.scr')
                end if
                
                deallocate(istr)
                
            end do
            
            if (InitConv == 0) then
                call Sleep(IV%NoSnap)
            else
                call Sleep(SolutionNumber)
            end if
            
            Iter = Iter + 1
            
            if (Iter < 4) then
                call CheckforConvergence(Iter, InitConv)
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
        if (IV%SystemType == 'W') then
            open(1, file=OutFolder//'/'//trim(IV%filename)//istr//'.rsd', form='formatted', STATUS="OLD")
        else
            open(1, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', form='formatted', STATUS="OLD")
        end if       
        inquire(1, size = FileSize)           
        LastLine = FileSize/106
            
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
        
        deallocate(istr)
            
    end subroutine FileCheckConvergence
    
end module CFD