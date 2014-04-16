module CFD
    
    use InputData
    use ReadData
    use Toolbox
    use GenerateInitialMeshes
    use CreateInitialNests
    
contains
    
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
            IntSystem = 10 + len(trim(IV%filename)) + len(istr) + len(pathLin_Prepro)
            allocate(character(len=IntSystem) :: strSystem)
            strSystem = pathLin_Prepro//' > '//trim(IV%filename)//istr//'.outpre'
            
        end if
        print *, 'Preprocessing Snapshot', i
        print *, ' '
        call system(strSystem)   ! System operating command called to activate fortran       
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
            
            ! Transfer Files in correct folder on Cluster
            call transferFilesLin()
            call system('chmod a+x ./FileCreateDir.scr')
            call system('./FileCreateDir.scr')
            
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
    
    subroutine Sleep()
    
        ! Variables
        implicit none
        integer :: i
    
        ! Body of Sleep
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
            do i = IV%NoSnap, 1, -1
            
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
            
                open(1, file='check.txt')
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
    
    end subroutine Sleep
    
end module CFD