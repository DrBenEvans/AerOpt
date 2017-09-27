module Optimization
    
    use CreateSnapshots
    use Toolbox
    use InputData
    use ReadData
    use GenerateMesh
    use CFD
    use FDGD
    use Smoothing
           
    contains
    
    subroutine SubOptimization()
    
        ! Variables
        implicit none
    
        ! Body of SubOptimization
        if (IV%Optimiser == 1) then
            call MCSOptimization()
        elseif (IV%Optimiser == 2) then
            call DEOptimization()
        elseif (IV%Optimiser == 3) then
            call PSOOptimization()
        end if
    
    end subroutine SubOptimization
    
    subroutine DEOptimization()
    
        ! Variables
        implicit none
        double precision :: Ftemp, Fopt, temp, Fibefore, CR, F
        integer :: i, j, k, l, ii, randInt, G, c
        double precision, dimension(:), allocatable :: NormFact, Fcompare, Frange
        double precision, dimension(:,:), allocatable :: Snapshots_Move, newSnapshots, tempAgent_Move, tempAgent
        integer, dimension(:), allocatable :: ind_Fi, ReferenceAgentsIndex
        logical :: Converge, exist
        character(len=5) :: strNoSnap

        allocate(NormFact(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(tempAgent_Move(IV%NoNests,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Snapshots_Move(IV%NoSnap,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Nests_Move(IV%NoNests,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Nests(IV%NoNests,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(tempAgent(IV%NoNests,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(ReferenceAgentsIndex(3),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Precoutput(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation " 
        allocate(Frange(2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation " 
        allocate(OV%converged(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        ! Specific for Adaptive Sampling
        allocate(newSnapshots(2,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Fcompare(2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "    
        
        ! Body of SubOptimization
        print *, ''
        print *, '*************************************'
        print *, '******    Start DE Optmization   ****'
        print *, '**                                 **'
        print *, '**  written by Dr. DAVID NAUMANN   **'
        print *, '**  supervised by Dr. BEN EVANS    **'
        print *, '*************************************'
        print *, ''
        
        ! Initializing Number of Generations
        OV%Gen = 1
        
        ! SVD parameters
        alpha = 1.0
        beta = 0.0
        
        ! Evaluate Fitness of first Generation
	    call date_and_time ( values = OV%timestart )
        call SubCFD((/ (j, j=1, IV%NoSnap) /), CS%Snapshots, IV%NoSnap)
        call PostSolverCheckInit(IV%NoSnap, 0)
        
        ! Extract moving initial Nests
        j = 1
        do i = 1, size(CS%cond)            
            if (CS%cond(i) == -1) then
                Snapshots_Move(:,j) = CS%Snapshots(:,i)
                j = j + 1
            end if
        end do
                   
        ! Normalize Snapshots(move) between 0 and 1
        do i = 1, IV%DoF        
            NormFact(i) = CS%MxDisp_Move(i,1) - CS%MxDisp_Move(i,2)
            Snapshots_Move(:,i) = (Snapshots_Move(:,i)-CS%MxDisp_Move(i,2))/NormFact(i)
        end do

        ! Extract Pressure of Snapshots
        call InitiatePressure()
        
        ! Extract Initial Fitness and Transfer to general Fitness vector
        call TransferInitialFitness(Snapshots_Move)
        Fibefore = OV%Fi(1)
        allocate(ind_Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "        
        
        ! Write Output File for Analysis including Initial and all moved Nests of each Generation
        call writeFitnessandNestoutput()
        
        ! Initiate DE Optimisation parameters
        CR = 0.9
        F = 0.8
      
        !!*** Loop over all DE Generations - each Generation creates new vector Agents ***!!
        do G = 2, IV%NoG
        OV%Gen = G
            print *, ''
            print *, '************************'
            print *, 'Generation ', OV%Gen
            print *, '************************'
            print *, ''
            print *, IV%Ma, IV%NoCN
            call timestamp()
            call date_and_time ( values = OV%timestart )

            !!****** Adaptive Sampling - Start New Jobs (first and last Fitness) *******!!
            if (OV%Gen > 2 .and. OV%Gen < IV%NoG .and. IV%AdaptSamp == .true.) then
                
                print *, 'Adaptive Sampling - Start Part 1 / 2'
                ! Extract First and Last Nest
                newSnapshots(1,:) = OV%Nests(1,:)
                newSnapshots(2,:) = OV%Nests(IV%NoNests,:)
                Fcompare(1) = OV%Fi(1)
                Fcompare(2) = OV%Fi(IV%NoNests)
                ! Get High Fidelity Solution
                call SubCFD((/(IV%NoSnap + 1), (IV%NoSnap + 2)/), newSnapshots, 2)
                print *, 'Adaptive Sampling - Finish Part 1 / 2'
                
            end if
                 
            !!*** Loop over all Agents ***!!
            print *, ''
            print *, 'Modify Agents for Generation', OV%Gen
            print *, ''
            do ii = 1, IV%NoNests
                
                ! Extract 3 Agents at random distinct from each other and the current agent
                ReferenceAgentsIndex = 0.0
                do j = 1, 3
                    
                    exist = .true.
                    do while (exist == .true.)
                        
                        exist = .false.
                        call random_number(CS%rn)
                        do k = 1, j    
                            if (nint(1 + (IV%NoSnap-1)*CS%rn) == ReferenceAgentsIndex(k) .or. nint(1 + (IV%NoSnap-1)*CS%rn) == ii) then
                                exist = .true.
                                EXIT
                            end if
                        end do
                        
                    end do
                    ReferenceAgentsIndex(j) = nint(1 + (IV%NoSnap-1)*CS%rn)
                end do
                
                ! Get Random Integer value in the interval [1, DoF]
                call random_number(CS%rn)
                randInt = nint(1 + (IV%DoF-1)*CS%rn)
                
                ! Mutate every dimension of the agent
                do j = 1, IV%DoF
                    call random_number(CS%rn)
                    if (CS%rn < CR .or. j == randInt) then ! j = randInt ensures, that at least one dimension changes to not have the same Agent in the new generation
                       tempAgent_Move(ii,j) = OV%Nests_Move(ReferenceAgentsIndex(1),j) + F*(OV%Nests_Move(ReferenceAgentsIndex(2),j) - OV%Nests_Move(ReferenceAgentsIndex(3),j))
                    else
                        tempAgent_Move(ii,j) = OV%Nests_Move(ii,j)                      
                    end if
                end do
               
                ! Check if out of bounds
                if (IV%constrain == .true.) then                
                    do k = 1, IV%DoF
                        if (tempAgent_Move(ii,k) > 1) then
                            tempAgent_Move(ii,k) = 1
                        end if
                        if (tempAgent_Move(ii,k) < 0) then
                           tempAgent_Move(ii,k) = 0 
                        end if                   
                    end do
                end if
                 
                ! Refill tempNests and de-normalize
                l = 1
                do k = 1, size(CS%cond)            
                    if (CS%cond(k) == -1) then
                        tempAgent(ii,k) = tempAgent_Move(ii,l)*NormFact(l) + CS%MxDisp_Move(l,2)
                        l = l + 1
                    end if
                end do

                ! Evaluate Fitness
                if (IV%POD == .true.) then
                    ! calculate Objective Function via POD (Fitness, Nest Locations)
                    call getObjectiveFunction(.true., Ftemp, tempAgent(ii,:))               
                    ! Output: ONE Fitnessvalue(Fi)
                    
                    ! Check if new Fitness is better than current agent
                    if (Ftemp > OV%Fi(ii)) then
                        OV%Fi(ii) = Ftemp
                        OV%Nests(ii,:) = tempAgent(ii,:)
                        OV%Nests_Move(ii,:) = tempAgent_Move(ii,:)
                        OV%Precoutput(ii) = OV%Precovery
                     end if
                end if
                
            end do
            
            if (IV%POD == .false.) then
            
                !Generate Full Fidelity Solution      
                call SubCFD((/ (j, j=1, IV%NoNests) /),tempAgent, IV%NoNests)          
                ! Check ALL new Nests
                print *, 'Generation: ', OV%Gen
                call PostSolverCheck(IV%NoNests, 1)
 
                ! Evaluate Fitness of Full Fidelity Nest Solutions
                print *, 'Extract Pressure of Generation', OV%Gen
                if (IV%NoDim == 2) then
                    call ExtractPressure(1, IV%NoNests)
                elseif (IV%NoDim == 3) then
                    call ExtractPressure_3D(1, IV%NoNests)
                end if
                do ii = 1, IV%NoNests
                    call getObjectiveFunction(.false., Ftemp, NoSnapshot=ii)
              
                    ! Check if new Fitness is better than current fitnesss
                    if (Ftemp > OV%Fi(ii)) then
                        OV%Fi(ii) = Ftemp
                        OV%Nests(ii,:) = tempAgent(ii,:)
                        OV%Nests_Move(ii,:) = tempAgent_Move(ii,:)
                        OV%Precoutput(ii) = OV%Precovery
                     end if
                end do
                deallocate(OV%pressure)
                deallocate(OV%MaLocal)
                deallocate(OV%pTamb)
            end if
            
            !!*** Adaptive Sampling - Finish and Integrate New Jobs (first and last Fitness) ***!!
            if (OV%Gen > 2 .and. OV%Gen < IV%NoG .and. IV%AdaptSamp == .true.) then
                call AdaptiveSampling(newSnapshots, Fcompare)
            end if
            
            ! Re-order Fitness in ascending order
            ind_Fi = (/ (i, i=1,IV%NoNests) /)
            call QSort(OV%Fi,size(OV%Fi), 'y', ind_Fi)
            ! Change to Descending Order for Maximization Problem
            do j = 1, int(IV%NoNests/2.0)
                temp = OV%Fi(j)
                OV%Fi(j) = OV%Fi(IV%NoNests-j+1)
                OV%Fi(IV%NoNests-j+1) = temp
                temp = ind_Fi(j)
                ind_Fi(j) = ind_Fi(IV%NoNests-j+1)
                ind_Fi(IV%NoNests-j+1) = temp
            end do
            
            ! Re-order Nests for next Generation
            OV%Nests_Move = OV%Nests_Move(ind_Fi,:)
            OV%Nests = OV%Nests(ind_Fi,:)
            OV%Precoutput = OV%Precoutput(ind_Fi)
            tempAgent = tempAgent(ind_Fi,:)
            tempAgent_Move = tempAgent_Move(ind_Fi,:)
            
            ! Print out Fitness values
            print *, 'Current best solutions:' , OV%Fi(1:nint(IV%NoNests*IV%Low2Top))
            
            ! write Nest and Fitness Output into File
            call writeFitnessandNestoutput()
            
            ! Store Files of all nests
            do j = 1, IV%NoNests
                if (IV%SystemType == 'W') then
                    call moveTopNestFilesWin(ind_Fi(j), j)
                    call generateEnSightFileWin()
                else
                    call moveTopNestFilesLin(ind_Fi(j),j)
                    if (IV%NoDim == 2) then
                        call generateEnSightFileLin()
                    else
                        call generateEnSightFileLin3D(ind_Fi(j))
                    end if
                end if
            end do
            Fibefore = OV%Fi(1)
            
        end do
        
        print *, 'Finished Differential Evolution'                  
        print *, ''
        print *, '************************'
        print *, 'Last Generation'
        print *, '************************'
        print *, ''
        print *, 'Optimum Fitness found:', OV%Fi(1)
        print *, ''
        print *, 'Optimum Geometry:'
        print *, ''
        print *, OV%Nests(ind_Fi(1),:)
        
    end subroutine DEOptimization
        
    subroutine PSOOptimization()
    
        ! Variables
        implicit none
        double precision :: Fopt, temp, Fibefore, control, W, rn1, rn2, F
        integer :: i, j, k, l, ii, G, c, minNeighbourSize, NeighbourhoodSize
        double precision, dimension(:), allocatable :: NormFact, Fcompare, Neighbourbest, LocalFi, Wrange
        double precision, dimension(:,:), allocatable :: Snapshots_Move, newSnapshots, LocalNest, LocalNest_Move, PSOvelocity
        integer, dimension(:), allocatable :: ind_Fi, NeighboursIndex
        logical :: Converge, exist
        character(len=5) :: strNoSnap

        allocate(NormFact(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(LocalNest_Move(IV%NoNests,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Snapshots_Move(IV%NoSnap,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Nests_Move(IV%NoNests,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Nests(IV%NoNests,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(LocalNest(IV%NoNests,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(PSOvelocity(IV%NoNests,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(NeighboursIndex(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Neighbourbest(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Wrange(2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Precoutput(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%converged(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        
        ! Specific for Adaptive Sampling
        allocate(newSnapshots(2,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Fcompare(2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "    
        
        ! Body of SubOptimization
        print *, ''
        print *, '*************************************'
        print *, '******    Start Optmization PSO   *******'
        print *, '**                                 **'
        print *, '**  written by Dr. DAVID NAUMANN   **'
        print *, '**  supervised by Dr. BEN EVANS    **'
        print *, '*************************************'
        print *, ''
        
        ! Initializing Generation
        OV%Gen = 1
        
        ! SVD parameters
        alpha = 1.0
        beta = 0.0
        
        ! Evaluate Fitness of first Generation
        call SubCFD((/ (j, j=1, IV%NoSnap) /), CS%Snapshots, IV%NoSnap)      
        call PostSolverCheckInit(IV%NoSnap, 0)
        
        ! Extract moving initial Nests
        j = 1
        do i = 1, size(CS%cond)            
            if (CS%cond(i) == -1) then
                Snapshots_Move(:,j) = CS%Snapshots(:,i)
                j = j + 1
            end if
        end do
                   
        ! Normalize Snapshots(move) between 0 and 1
        do i = 1, IV%DoF        
            NormFact(i) = CS%MxDisp_Move(i,1) - CS%MxDisp_Move(i,2)
            Snapshots_Move(:,i) = (Snapshots_Move(:,i)-CS%MxDisp_Move(i,2))/NormFact(i)
        end do

        ! Extract Pressure of Snapshots
        call InitiatePressure()
        
        ! Extract Initial Fitness and Transfer to general Fitness vector
        call TransferInitialFitness(Snapshots_Move)
        Fibefore = OV%Fi(1)
        allocate(ind_Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(LocalFi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        
        
        ! Write Output File for Analysis including Initial and all moved Nests of each Generation
        call writeFitnessandNestoutput()
        
        ! Initialisation of Optimisation Parameters
        minNeighbourSize = nint(0.25*IV%NoSnap)
        Neighbourhoodsize = maxval((/1, minNeighbourSize/))
        W = 0.1
        c = 0
        control = 0.75
        Wrange = (/ -0.5, 0.9/)
        PSOvelocity = 0.0
        LocalNest_Move = OV%Nests_Move
        LocalNest = OV%Nests
        F = 0.8
        
        !!*** Loop over all PSO Generations - each Generation creates new Particles ***!!
        do G = 2, IV%NoG
        OV%Gen = G
            print *, ''
            print *, '************************'
            print *, 'Generation ', OV%Gen
            print *, '************************'
            print *, ''
            print *, IV%Ma, IV%NoCN
            call timestamp()
            call date_and_time ( values = OV%timestart ) 

            !!****** Adaptive Sampling - Start New Jobs (first and last Fitness) *******!!
            if (OV%Gen > 2 .and. OV%Gen < IV%NoG .and. IV%AdaptSamp == .true.) then
                
                print *, 'Adaptive Sampling - Start Part 1 / 2'
                ! Extract First and Last Nest
                newSnapshots(1,:) = OV%Nests(1,:)
                newSnapshots(2,:) = OV%Nests(IV%NoNests,:)
                Fcompare(1) = OV%Fi(1)
                Fcompare(2) = OV%Fi(IV%NoNests)
                ! Get High Fidelity Solution
                call SubCFD((/(IV%NoSnap + 1), (IV%NoSnap + 2)/), newSnapshots, 2)
                print *, 'Adaptive Sampling - Finish Part 1 / 2'
                
            end if
            
            !!*** Loop over all Agents ***!!
            do ii = 1, IV%NoNests
                
                ! Extract a random subset of all Agents to act as Neighbours
                NeighboursIndex = 0.0
                do j = 1, Neighbourhoodsize
                    
                    exist = .true.
                    do while (exist == .true.)
                        
                        exist = .false.
                        call random_number(CS%rn)
                        do k = 1, j    
                            if (nint(1 + (IV%NoSnap-1)*CS%rn) == NeighboursIndex(k) .or. nint(1 + (IV%NoSnap-1)*CS%rn) == ii) then
                                exist = .true.
                                EXIT
                            end if
                        end do
                        
                    end do
                    NeighboursIndex(j) = nint(1 + (IV%NoSnap-1)*CS%rn)
                end do
                
                ! Extract Neighbours best Position
                Neighbourbest = OV%Nests_Move(minval(NeighboursIndex(1:Neighbourhoodsize)),:)
                
                ! Compute new velocity
                do j = 1, IV%DoF
                    call random_number(CS%rn)
                    rn1 = CS%rn
                    call random_number(CS%rn)
                    rn2 = CS%rn
                    PSOvelocity(ii,j) = W*PSOvelocity(ii,j) + control*rn1*(OV%Nests_Move(ii,j) - LocalNest_Move(ii,j)) + control*rn2*(Neighbourbest(j) - LocalNest_Move(ii,j))
                end do
                
                ! Update Position
                LocalNest_Move(ii,:) = LocalNest_Move(ii,:) + PSOvelocity(ii,:)
               
                ! Check if out of bounds
                if (IV%constrain == .true.) then                
                    do k = 1, IV%DoF
                        if (LocalNest_Move(ii,k) > 1) then
                            LocalNest_Move(ii,k) = 1
                        end if
                        if (LocalNest_Move(ii,k) < 0) then
                           LocalNest_Move(ii,k) = 0 
                        end if                   
                    end do
                end if           
                
                ! Refill Nests
                l = 1
                do k = 1, size(CS%cond)            
                    if (CS%cond(k) == -1) then
                        LocalNest(ii,k) = LocalNest_Move(ii,l)*NormFact(l) + CS%MxDisp_Move(l,2)
                        l = l + 1
                    end if
                end do
                
                ! Evaluate Fitness
                if (IV%POD == .true.) then
                    ! calculate Objective Function via POD (Fitness, Nest Locations)
                    call getObjectiveFunction(.true., LocalFi(ii), LocalNest(ii,:))               
                    ! Output: ONE Fitnessvalue(Fi)
                    
                    ! Check if new Fitness is better than current local best
                    if (LocalFi(ii) > OV%Fi(ii)) then
                        OV%Fi(ii) = LocalFi(ii)
                        OV%Nests(ii,:) = LocalNest(ii,:)
                        OV%Nests_Move(ii,:) = LocalNest_Move(ii,:)
                        OV%Precoutput(ii) = OV%Precovery
                    end if
                end if
            
            end do
            
            if (IV%POD == .false.) then
                
                ! Store moved Nests in Output Analysis File for Restart
                open(29,file=newdir//'/Neststemp.txt',form='formatted',status='unknown')
                write(29,'(<IV%NoNests>f17.10)') OV%Nests
                close(29)
            
                !Generate Full Fidelity Solution      
                call SubCFD((/ (j, j=1, IV%NoNests) /),LocalNest, IV%NoNests)              
                ! Check ALL new Nests
                print *, 'Generation: ', OV%Gen
                call PostSolverCheck(IV%NoNests, 1)   
 
                ! Evaluate Fitness of Full Fidelity Nest Solutions
                print *, 'Extract Pressure of Generation', OV%Gen
                if (IV%NoDim == 2) then
                    call ExtractPressure(1, IV%NoNests)
                elseif (IV%NoDim == 3) then
                    call ExtractPressure_3D(1, IV%NoNests)
                end if
                do ii = 1, IV%NoNests
                    call getObjectiveFunction(.false., LocalFi(ii), NoSnapshot=ii)
              
                    ! Check if new Fitness is better than current local best
                    if (LocalFi(ii) > OV%Fi(ii)) then
                        OV%Fi(ii) = LocalFi(ii)
                        OV%Nests(ii,:) = LocalNest(ii,:)
                        OV%Nests_Move(ii,:) = LocalNest_Move(ii,:)
                        OV%Precoutput(ii) = OV%Precovery
                    end if
                end do
                deallocate(OV%pressure)
                deallocate(OV%MaLocal)
                deallocate(OV%pTamb)
            end if
                    
            !!*** Adaptive Sampling - Finish and Integrate New Jobs (first and last Fitness) ***!!
            if (OV%Gen > 2 .and. OV%Gen < IV%NoG .and. IV%AdaptSamp == .true.) then
                call AdaptiveSampling(newSnapshots, Fcompare)
            end if
            
            ! Re-order Fitness in ascending order
            ind_Fi = (/ (i, i=1,IV%NoNests) /)
            call QSort(OV%Fi,size(OV%Fi), 'y', ind_Fi)
            ! Change to Descending Order for Maximization Problem
            do j = 1, int(IV%NoNests/2.0)
                temp = OV%Fi(j)
                OV%Fi(j) = OV%Fi(IV%NoNests-j+1)
                OV%Fi(IV%NoNests-j+1) = temp
                temp = ind_Fi(j)
                ind_Fi(j) = ind_Fi(IV%NoNests-j+1)
                ind_Fi(IV%NoNests-j+1) = temp
            end do
            
            ! Re-order Nests for next Generation
            LocalFi = LocalFi(ind_Fi)
            OV%Precoutput = OV%Precoutput(ind_Fi)
            LocalNest = LocalNest(ind_Fi,:)
            OV%Nests = OV%Nests(ind_Fi,:)
            OV%Nests_Move = OV%Nests_Move(ind_Fi,:)
            LocalNest_Move = LocalNest_Move(ind_Fi,:)   
            
            ! Print out Fitness values
            print *, 'Current best solutions:' , OV%Fi(1:nint(IV%NoNests*IV%Low2Top))
            
            ! write Nest and Fitness Output into File
            call writeFitnessandNestoutput()
            
            ! Update Neighbourhood
            ! The greater W the less likely it is for the PSO to get stuck in a local optimum
            ! hence for c > 10 a change is needed and the velocity direction alternates +/- with the maximum W values
            ! However for it not to go just forth and back, the Wrange values are asymmetric
            ! Also, for it to not shoot past the boundary it is important to keep -1 < W < +1
            ! In close vicinity, if it is here and there finding a new value ( 5 < c < 10), W is reduced by W = W/2 so it can search in the region of the local optimum
            ! If it comes to close and plenty of changes occur (c < 2), you have to assume it approaches a local minimum and the W value needs to increase (W = 2*W)
            if (Fibefore /= OV%Fi(1)) then
                c = maxval((/0,c-1/))
                Neighbourhoodsize = maxval((/1, minNeighbourSize/))
                if (c < 2) then
                    W = 2*W
                elseif (c > 5) then
                    W = W/2
                end if
                if (W < Wrange(1)) then
                    W = Wrange(1)
                end if
                if (W > Wrange(2)) then
                    W = Wrange(2)
                end if
            else
                if (c > 10 .and. W > 0) then
                    W = Wrange(1)
                else
                    W = Wrange(2)
                end if
                c = c + 1
                Neighbourhoodsize = minval((/(Neighbourhoodsize + minNeighbourSize),(IV%NoNests-1)/))
            end if
            
            !Store all nests
            do j = 1, IV%NoNests
                if (IV%SystemType == 'W') then
                    call moveTopNestFilesWin(ind_Fi(j), j)
                    !call generateEnSightFileWin()
                else
                    call moveTopNestFilesLin(ind_Fi(j),j)
                    if (IV%NoDim == 2) then
                        call generateEnSightFileLin()
                    else
                        call generateEnSightFileLin3D(ind_Fi(j))
                    end if
                end if
            end do
            Fibefore = OV%Fi(1)
            
        end do
        
        print *, 'Finished Particle Swarm Optimization'                  
        print *, ''
        print *, '************************'
        print *, 'Last Generation'
        print *, '************************'
        print *, ''
        print *, 'Optimum Fitness found:', OV%Fi(1)
        print *, ''
        print *, 'Optimum Geometry:'
        print *, ''
        print *, OV%Nests(ind_Fi(1),:)
        
    end subroutine PSOOptimization
    
    subroutine MCSOptimization()
    
        ! Variables
        implicit none
        double precision :: Ac, Ftemp, Fopt, temp, Fibefore
        integer :: i, j, k, l, ii, NoSteps, NoTop, NoDiscard, randomNest, G
        double precision, dimension(:), allocatable :: NormFact, tempNests_Move, dist, tempNests, Fcompare, newMaLocal
        double precision, dimension(:,:), allocatable :: Snapshots_Move, newSnapshots, TopNest, TopNest_Move
        integer, dimension(:), allocatable :: ind_Fi, ind_Fitrack
        logical :: Converge
        character(len=5) :: strNoSnap

        allocate(NormFact(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(tempNests_Move(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(dist(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Snapshots_Move(IV%NoSnap,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Nests_Move(IV%NoNests,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Nests(IV%NoNests,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(tempNests(maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%Precoutput(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(OV%converged(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        
        ! Specific for Adaptive Sampling
        allocate(newSnapshots(2,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Fcompare(2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "    
        
        ! Body of SubOptimization
        print *, ''
        print *, '*************************************'
        print *, '****    Start MCS Optmization   *****'
        print *, '**                                 **'
        print *, '**  written by Dr. DAVID NAUMANN   **'
        print *, '**  supervised by Dr. BEN EVANS    **'
        print *, '**  supported by Dr. SEAN WALTON   **'
        print *, '*************************************'
        print *, ''
    
        ! Initializing optimization parameters
        NoDiscard = floor(IV%Low2Top*IV%NoNests) 
        NoTop = IV%NoNests - NoDiscard
        !IV%Aconst = (sqrt(real(IV%DoF))/IV%NoLeviSteps)*IV%Aconst
        tempNests = (/ (0, i=1,(maxDoF)) /)
        OV%Gen = 1
        
        ! SVD parameters
        alpha = 1.0
        beta = 0.0
       
        ! Specific Parameters required for Top Nest monitoring
        allocate(TopNest(NoTop,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(TopNest_Move(NoTop,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation " 
        
        ! Evaluate Fitness of first Generation
        call date_and_time ( values = OV%timestart )
        call SubCFD((/ (j, j=1,IV%NoSnap) /), CS%Snapshots, IV%NoSnap)
        call PostSolverCheckInit(IV%NoSnap, 0)
       
        ! Extract moving initial Nests
        j = 1
        do i = 1, size(CS%cond)            
            if (CS%cond(i) == -1) then
                Snapshots_Move(:,j) = CS%Snapshots(:,i)
                j = j + 1
            end if
        end do
                   
        ! Normalize Snapshots(move) between 0 and 1
        do i = 1, IV%DoF        
            NormFact(i) = CS%MxDisp_Move(i,1) - CS%MxDisp_Move(i,2)
            Snapshots_Move(:,i) = (Snapshots_Move(:,i)-CS%MxDisp_Move(i,2))/NormFact(i)
        end do

        ! Extract Pressure of Snapshots
        call InitiatePressure()
        
        ! Extract Initial Fitness and Transfer to general Fitness vector
        call TransferInitialFitness(Snapshots_Move)
        Fibefore = OV%Fi(1)
        allocate(ind_Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(ind_Fitrack(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "   
    
        ! Write Output File for Analysis including Initial and all moved Nests of each Generation
        call writeFitnessandNestoutput()
       
        !!*** Loop over all Cuckoo Generations - each Generation creates new Nests ***!!
        do G = 2, IV%NoG
        OV%Gen = G
        ind_Fitrack = (/ (i, i=1,IV%NoNests) /)
            print *, ''
            print *, '************************'
            print *, 'Generation ', OV%Gen
            print *, '************************'
            print *, ''
            print *, IV%Ma, IV%NoCN
            call timestamp()
            ! Store timestamp
            call date_and_time ( values = OV%timestart )
            
            !!****** Adaptive Sampling - Start New Jobs (first and last Fitness) *******!!
            if (OV%Gen > 2 .and. OV%Gen < IV%NoG .and. IV%AdaptSamp == .true.) then
                
                print *, 'Adaptive Sampling - Start Part 1 / 2'
                ! Extract First and Last Nest
                newSnapshots(1,:) = OV%Nests(1,:)
                newSnapshots(2,:) = OV%Nests(IV%NoNests,:)
                Fcompare(1) = OV%Fi(1)
                Fcompare(2) = OV%Fi(IV%NoNests)
                ! Get High Fidelity Solution
                call SubCFD((/(IV%NoSnap + 1), (IV%NoSnap + 2)/), newSnapshots, 2)
                print *, 'Adaptive Sampling - Finish Part 1 / 2'
                
            end if
            
            !!*** Loop over Discarded Nests ***!!
            print *, ''
            print *, 'Modify Discarded Cuckoos for Generation', OV%Gen
            print *, ''
            do ii = IV%NoNests, (NoTop + 1), -1
                
                ! Perform Random Walk using Levy Flight with a Cauchy Distribution
                Ac = IV%Aconst/(sqrt(OV%Gen*1.0)*100.0)
                call random_number(CS%rn)
                NoSteps = nint(log(CS%rn)*(-IV%NoLeviSteps))
                NoSteps = minval((/ NoSteps, IV%NoLeviSteps /))
                tempNests_Move = Ac*LevyWalk(NoSteps, IV%DoF) + OV%Nests_Move(ii,:)
               
                ! Check if out of bounds
                if (IV%constrain == .true.) then                
                    do k = 1, IV%DoF
                        if (tempNests_Move(k) > 1) then
                            tempNests_Move(k) = 1
                        end if
                        if (tempNests_Move(k) < 0) then
                           tempNests_Move(k) = 0
                        end if                   
                    end do
                end if
                 
                ! Refill tempNests
                l = 1
                do k = 1, size(CS%cond)            
                    if (CS%cond(k) == -1) then
                        tempNests(k) = (tempNests_Move(l) - 0.5)*NormFact(l)
                        l = l + 1
                    end if
                end do

                if (IV%POD == .true.) then
					call DetermineStrLen(istr, ii)
                    ! calculate Objective Function via POD (Fitness, Nest Locations)
                    call getObjectiveFunction(.true., OV%Fi(ii), tempNests, ii)               
                    ! Output: ONE Fitnessvalue(Fi)
					deallocate(istr)
                end if
               
                ! Embed temporary Nest into Nests
                OV%Nests_Move(ii,:) = tempNests_Move
                OV%Nests(ii,:) = tempNests
                
            end do        
            
            !!*** Loop over Top Nests ***!!
            print *, ''
            print *, 'Modify Top Cuckoos for Generation', OV%Gen
            print *, ''
            do ii = 1, NoTop
           
                ! Pick one of the Top Nests
                call random_number(CS%rn)
                randomNest = nint((1 + (NoTop - 1)*CS%rn))
                if (randomNest == ii) then  ! Same Nest
                
                    ! Perform Random Walk instead                   
                    Ac = IV%Aconst/(OV%Gen**2.0)
                    call random_number(CS%rn)
                    NoSteps = nint(log(CS%rn)*(-IV%NoLeviSteps))
                    NoSteps = minval((/ NoSteps, IV%NoLeviSteps /))
                    tempNests_Move = Ac*LevyWalk(NoSteps, IV%DoF) + OV%Nests_Move(ii,:)
                    
                else    ! Different Nest
                    
                    if (OV%Fi(ii) > OV%Fi(randomNest)) then
                        
                        ! Cross-bread Nests in Direction of Nest j by Golden Ratio
                        dist = OV%Nests_Move(randomNest,:) - OV%Nests_Move(ii,:)   ! Calculate Distance between Nests
                        dist = dist/(0.5*(1+sqrt(5.0)))                           ! Apply Golden Ratio
                        tempNests_Move = OV%Nests_Move(ii,:) + dist              ! Move Less Fit Nest
                        
                    elseif (OV%Fi(randomNest) > OV%Fi(ii)) then
                        
                        ! Cross-bread in Direction of randomNest by Golden Ratio
                        dist = OV%Nests_Move(ii,:) - OV%Nests_Move(randomNest,:)   ! Calculate Distance between Nests
                        dist = dist/(0.5*(1+sqrt(5.0)))                           ! Apply Golden Ratio
                        tempNests_Move = OV%Nests_Move(randomNest,:) + dist     ! Move Less Fit Nest
                        
                    else
                        
                        ! Fitness is the same: Cross-bread Half Way
                        dist = OV%Nests_Move(randomNest,:) - OV%Nests_Move(ii,:)   ! Calculate Distance between Nests
                        dist = dist*0.5                                         ! Apply Golden Ratio
                        tempNests_Move = OV%Nests_Move(ii,:) + dist              ! Move Less Fit Nest
                
                    end if
                   
                end if
                
                k = 1
                do while (k /= (IV%NoNests+1))
                    ! Check if out of bounds             
                    if (IV%constrain == .true.) then                
                        do k = 1, IV%DoF
                            if (tempNests_Move(k) > 1) then
                                tempNests_Move(k) = 1
                            end if
                            if (tempNests_Move(k) < 0) then
                                tempNests_Move(k) = 0
                            end if                   
                        end do
                    end if
                    
                    ! Check, if Nest already exists
                    do k = 1, IV%NoNests                  
                        if (all(tempNests_move == OV%Nests_move(k,:)) .or. all(tempNests_move == TopNest_move(ceiling(k*(1-IV%Low2Top)),:))) then
                            print*, 'Nest already exists'
                            ! Perform Random Walk instead                   
                            Ac = IV%Aconst/(OV%Gen**2.0)
                            call random_number(CS%rn)
                            NoSteps = nint(log(CS%rn)*(-IV%NoLeviSteps))
                            NoSteps = minval((/ NoSteps, IV%NoLeviSteps /))
                            tempNests_Move = Ac*LevyWalk(NoSteps, IV%DoF) + OV%Nests_Move(ii,:)
                            EXIT
                        end if
                    end do
                end do
                
                ! Refill tempNests and de-normalize
                l = 1
                do k = 1, size(CS%cond)            
                    if (CS%cond(k) == -1) then
                        tempNests(k) = tempNests_Move(l)*NormFact(l) + CS%MxDisp_Move(l,2)
                        l = l + 1
                    end if
                end do
                
                ! Store tempNests
                TopNest(ii,:) = tempNests
                TopNest_move(ii,:) = tempNests_move
                
                if (IV%POD == .true.) then
					call DetermineStrLen(istr, ii)
                    ! Update Objective Function via POD (Fitness, Nest Locations)
                    call getObjectiveFunction(.true., Ftemp, tempNests, ii)
                    ! Output: ONE Fitnessvalue(Fi)
                    deallocate(istr)
                    
                    ! Check if new Fitness is better than a Random Top Nest, If yes replace values
                    call random_number(CS%rn)
                    randomNest = nint((1 + (NoTop - 1)*CS%rn))
                    if (Ftemp > OV%Fi(randomNest)) then
                        OV%Nests_Move(randomNest,:) = tempNests_Move
                        OV%Nests(randomNest,:) = tempNests
                        OV%Fi(randomNest) = Ftemp
                        OV%Precoutput(randomNest) = OV%Precovery
                        ind_Fitrack(randomNest) = ii
                     end if
                end if
                
            end do
            
            if (IV%POD == .false.) then  
                
                ! Store moved Nests in Output Analysis File
                open(29,file=newdir//'/Neststemp.txt',form='formatted',status='unknown')
                write(29,'(<IV%NoNests>f17.10)') OV%Nests
                close(29)
            
                !Generate Full Fidelity Solution of new TopNest      
                call SubCFD((/ (j, j=1,NoTop) /), TopNest, NoTop)
                ! Full Fidelity Solutions of Discarded Nests       
                call SubCFD((/ (j, j=(NoTop + 1), IV%NoNests) /), OV%Nests((NoTop + 1):IV%NoNests,:), NoDiscard)            
                ! Check ALL new Nests
                print *, 'Generation: ', OV%Gen
                call PostSolverCheck(IV%NoNests, 1)   
 
                ! Evaluate Fitness of Full Fidelity Nest Solutions
                print *, 'Extract Pressure of Generation', OV%Gen
                if (IV%NoDim == 2) then
                    call ExtractPressure(1, IV%NoNests)
                elseif (IV%NoDim == 3) then
                    call ExtractPressure_3D(1, IV%NoNests)
                end if
                do ii = 1, NoTop
                    call getObjectiveFunction(.false., Ftemp, NoSnapshot=ii)
              
                    ! Check if new Fitness is better than a Random Top Nest, If yes replace values
                    call random_number(CS%rn)
                    randomNest = nint((1 + (NoTop - 1)*CS%rn))
                    if (Ftemp > OV%Fi(randomNest)) then
                        OV%Nests_Move(randomNest,:) = TopNest_Move(ii,:)
                        OV%Nests(randomNest,:) = TopNest(ii,:)
                        OV%Fi(randomNest) = Ftemp
                        OV%Precoutput(randomNest) = OV%Precovery
                        ind_Fitrack(randomNest) = ii
                    end if
                end do

                do ii = (NoTop + 1), IV%NoNests
                    call getObjectiveFunction(.false., OV%Fi(ii), NoSnapshot=ii)
                    OV%Precoutput(ii) = OV%Precovery
                end do
                deallocate(OV%pressure)
                deallocate(OV%MaLocal)
                deallocate(OV%pTamb)
            end if
            
            !!*** Adaptive Sampling - Finish and Integrate New Jobs (first and last Fitness) ***!!
            if (OV%Gen > 2 .and. OV%Gen < IV%NoG .and. IV%AdaptSamp == .true.) then
                call AdaptiveSampling(newSnapshots, Fcompare)
            end if
            
            ! Re-order Fitness in ascending order
            ind_Fi = (/ (i, i=1,IV%NoNests) /)
            call QSort(OV%Fi,size(OV%Fi), 'y', ind_Fi)
            ! Change to Descending Order for Maximization Problem
            do j = 1, int(IV%NoNests/2.0)
                temp = OV%Fi(j)
                OV%Fi(j) = OV%Fi(IV%NoNests-j+1)
                OV%Fi(IV%NoNests-j+1) = temp
                temp = ind_Fi(j)
                ind_Fi(j) = ind_Fi(IV%NoNests-j+1)
                ind_Fi(IV%NoNests-j+1) = temp
            end do
            
            ! Re-order Nests for next Generation
            OV%Nests_Move = OV%Nests_Move(ind_Fi,:)
            OV%Nests = OV%Nests(ind_Fi,:)
            OV%Precoutput = OV%Precoutput(ind_Fi)
            
            ! Print out Fitness values
            print *, 'Current best solutions:' , OV%Fi(1:nint(IV%NoNests*(1-IV%Low2Top)))
            
            ! write Nest and Fitness Output into File
            call writeFitnessandNestoutput()
            
            ! Store Files of Top 5 % fraction of Nests in TopFolder
            if (IV%POD == .false.) then          
                do j = 1, IV%NoNests
                    if (IV%SystemType == 'W') then
                        call moveTopNestFilesWin(ind_Fitrack(ind_Fi(j)), j)
                        call generateEnSightFileWin()
                    else
                        call moveTopNestFilesLin(ind_Fitrack(ind_Fi(j)), j)
                        if (IV%NoDim == 2) then
                            call generateEnSightFileLin()
                        else
                            call generateEnSightFileLin3D(ind_Fitrack(ind_Fi(j)))
                        end if
                    end if
                end do
	        else
                do j = 1, IV%NoNests
                    if (IV%SystemType == 'W') then
                        call moveTopNestFilesWin(ind_Fitrack(ind_Fi(j)), j)
                    else
                        call moveTopNestFilesLin_POD(ind_Fitrack(ind_Fi(j)))
                    end if
                end do
            end if
            Fibefore = OV%Fi(1)
            
        end do
        
        print *, 'Finished Cuckoo Search'                  
        print *, ''
        print *, '************************'
        print *, 'Last Generation'
        print *, '************************'
        print *, ''
        print *, 'Optimum Fitness found:', OV%Fi(1)
        print *, ''
        print *, 'Optimum Geometry:'
        print *, ''
        print *, OV%Nests(ind_Fi(1),:)
        
    end subroutine MCSOptimization
    
    subroutine InitiatePressure()
    
        ! Variables
        implicit none
    
        ! Body of InitiatePressure
        if (IV%ObjectiveFunction == 2 .or. IV%ObjectiveFunction == 7) then
            call getengineInlet() ! Get boundary nodes, that define the Engine Inlet Plane
            ! Output: Engine Inlet Nodes(engInNodes)     
        end if
        
        call timestamp()
        if (IV%POD == .true.) then
            call AllocateModesCoeff()
            ! Output: allocated modes and coeff based on the Number of POD Modes desired. If < 0, all Modes are considered.
            call POD()
            ! Output: Modes and Coefficients of POD  
            call ComputeRBFWeights()
            ! Output: Weights for RBF interpolation
        else
            if (IV%NoDim == 2) then
                call ExtractPressure(1, IV%NoSnap)
            elseif (IV%NoDim == 3) then
                call ExtractPressure_3D(1, IV%NoSnap)
            end if
        end if             
        call timestamp()  
        
    end subroutine InitiatePressure
    
    subroutine TransferInitialFitness(Snapshots_Move)
    
        ! Variables
        implicit none
        integer :: ii, j
        double precision :: Ftemp, temp
        integer, dimension(:), allocatable :: ind_Fi_initial
        double precision, dimension(:), allocatable :: Fi_initial, Precoutput_temp
        double precision, dimension(:,:), allocatable :: Snapshots_Move
        
        ! Body of TransferInitialFitness
        allocate(ind_Fi_initial(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Fi_initial(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Precoutput_temp(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        
        ! Determine Objective Function
        do ii = 1, IV%NoSnap
            call getObjectiveFunction(.false., Ftemp, CS%Snapshots(ii,:), NoSnapshot=ii)
            Fi_initial(ii) = Ftemp
            Precoutput_temp(ii) = OV%Precovery
        end do
        deallocate(OV%pressure)
        deallocate(OV%MaLocal)
        deallocate(OV%pTamb)
        ! Output: Fitness depend on user Input (objective Function)
       
        open(19,file=newdir//'/Fitness_0.txt', form='formatted',status='unknown')
        !open(19,file=TopFolder//'/Fitness_0.txt', form='formatted',status='unknown')
        write(19,'(1I1)',advance="no") 0
        write(19,'(1f25.10)',advance="no") Fi_initial(1)
        if (IV%objectivefunction == 7) then
            write(19,'(1f17.10)') OV%Precoutput(1)
        end if
        close(19)
      
        ! Pass on Top Snapshot parameters to Nests
        ind_Fi_initial = (/ (j, j=1,IV%NoSnap) /)
        call QSort(Fi_initial,size(Fi_initial), 'y', ind_Fi_initial) ! Result in Ascending order
        ! Change to Descending Order for Maximization Problem
        do j = 1, int(IV%NoSnap/2.0)
            temp = Fi_initial(j)
            Fi_initial(j) = Fi_initial(IV%NoSnap-j+1)
            Fi_initial(IV%NoSnap-j+1) = temp
            temp = ind_Fi_initial(j)
            ind_Fi_initial(j) = ind_Fi_initial(IV%NoSnap-j+1)
            ind_Fi_initial(IV%NoSnap-j+1) = temp
        end do
        OV%Fi = Fi_initial(1:IV%NoNests)
        OV%Nests_Move = Snapshots_Move(ind_Fi_initial(1:IV%NoNests),:)
        OV%Nests = CS%Snapshots(ind_Fi_initial(1:IV%NoNests),:)
        OV%Precoutput = Precoutput_temp(ind_Fi_initial(1:IV%NoNests))

        do j = 1, IV%NoNests
            if (IV%SystemType == 'W') then
               call moveTopNestFilesWin(ind_Fi_initial(j), j)
               call generateEnSightFileWin()
            else
                call moveTopNestFilesLin(ind_Fi_initial(j),j)
                if (IV%NoDim == 2) then
                   call generateEnSightFileLin()
                else
                    call generateEnSightFileLin3D(ind_Fi_initial(j))
                end if
            end if
        end do
        deallocate(ind_Fi_initial)
        deallocate(Fi_initial)
        deallocate(Precoutput_temp)
        
    end subroutine TransferInitialFitness
    
    subroutine ExtractPressure(Start, Ending)
        
        ! Variables
        implicit none
        integer :: Start, Ending, Length, i, k, j, xindex
        double precision :: Vamb, rho_amb, Rspec, gamma, qamb
        double precision, dimension(:), allocatable :: Output, Vx, Vy, rho, e
        character(len=:), allocatable :: istr

        ! Body of Extract Pressure
        Length = Ending - Start + 1
        allocate(rho(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(e(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(Vx(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(Vy(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(Output(6*RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(OV%pressure(RD%np, Length),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(OV%MaLocal(size(OV%engInNodes, dim = 1), Length),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(OV%pTamb(Length),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
       
        ! Precalculate ambient Parameters
        print *, 'Start Pressure Extraction'
        Vamb = IV%Ma*sqrt(IV%gamma*IV%R*IV%Tamb)          ! ambient velocity
        rho_amb = (IV%Pamb)/(IV%Tamb*IV%R)             ! ambient Density
            
        !Extract pressure of Snapshot Output file      
        do i = Start, Ending
           
            ! Determine correct String number
            call DetermineStrLen(istr, i)
           
            open(11, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk', form='formatted',status='old')     
            read(11, *) Output  ! index, rho, Vx, Vy, Vz, e
            
            k = 0
            do j = 1, (6*RD%np), 6
                k = k + 1
                rho(k) = Output(j+1)
                Vx(k) = Output(j+2)
                Vy(k) = Output(j+3)
                e(k) = Output(j+4)
                !turbulence(k) = Output(j+5)
            end do
                
            ! Calculate Pressure
            !pressure(:,(i - Start + 1)) = rho_amb*(IV%gamma - 1.0)*rho*(Vamb**2)*(e - 0.5*(Vx**2 + Vy**2)) ! Dimensional
            OV%pressure(:,(i - Start + 1)) = (IV%gamma - 1.0)*rho*((e - 0.5*(Vx**2 + Vy**2))) ! Non-dimensional  
            ! Old Bernoulli Equation to calculate non-dimensional pressure:  pressure(:,i) = e + (1.0/2.0)*(IV%Ma**2)*rho*(Vx*Vx + Vy*Vy) 
            
            ! Local Mach number non-dimensional
            OV%MaLocal(:,(i - Start + 1)) = sqrt(Vx(OV%engInNodes)**2 + Vy(OV%engInNodes)**2)/(sqrt(IV%gamma*OV%pressure(OV%engInNodes,(i - Start + 1))/rho(OV%engInNodes)))
            
            
! ambient pressure calculated based on Reynolds number and Ma and input Temperature            
            ! free-stream dynamic pressure
            if (IV%AlphaInflowDirection > 90) then
                xindex = maxloc(RD%Coord(:,1), dim = 1)
            else
                xindex = minloc(RD%Coord(:,1), dim = 1)
            end if
            OV%pTamb(i - Start + 1) = OV%pressure(xindex,(i - Start + 1))*(1.0+(IV%gamma-1.0)*0.5*IV%Ma**2)**(IV%gamma/(IV%gamma-1.0))
             
            close(11)
            deallocate(istr)
                
        end do
        print *,'Pressure of .unk files extracted'
        deallocate(Output)
        deallocate(Vx)
        deallocate(Vy)
        deallocate(e)
        deallocate(rho)

    end subroutine ExtractPressure
    
    subroutine ExtractPressure_3D(Start, Ending)
        
        ! Variables
        implicit none
        integer :: Start, Ending, Length, i, k, j, nop
        double precision :: Vamb, rho_amb, Rspec, gamma
        double precision, dimension(:), allocatable :: Vx, Vy, Vz, rho, e
        character(len=:), allocatable :: istr

        ! Body of Extract Pressure
        Length = Ending - Start + 1
        allocate(rho(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(e(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(Vx(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(Vy(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(Vz(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        allocate(OV%pressure(RD%np, Length),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ExtractPressure "
        
        ! Precalculate ambient Parameters
        print *, 'Start Pressure Extraction'
        Vamb = IV%Ma*sqrt(IV%gamma*IV%R*IV%Tamb)          ! ambient velocity
        rho_amb = (IV%Pamb)/(IV%Tamb*IV%R)             ! ambient Density
            
        !Extract pressure of Snapshot Output file      
        do i = Start, Ending
           
            ! Determine correct String number
            call DetermineStrLen(istr, i)
           
            open(11, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.unk', form='unformatted',status='old')     
            read(11) nop
            if (RD%np /= nop) then
                STOP 'ERROR: .unk file does not correspond with .plt file. Check file transfer'
            end if
            read(11) rho, Vx, Vy, Vz, e
                
            ! Calculate Pressure
            !OV%pressure(:,(i - Start + 1)) = (IV%gamma - 1.0)*rho*(e - 0.5*(Vx**2 + Vy**2 + Vz**2))*rho_amb*(Vamb**2) ! Dimensional
            OV%pressure(:,(i - Start + 1)) =  (IV%gamma - 1.0)*rho*(e - 0.5*(Vx**2 + Vy**2 + Vz**2)) ! Non-dimensional  
            ! Old Bernoulli Equation to calculate non-dimensional pressure:  pressure(:,i) = e + (1.0/2.0)*(IV%Ma**2)*rho*(Vx*Vx + Vy*Vy + Vz*Vz) 
            
            close(11)
                !open(11, file=newdir//'/'//OutFolder//'/pressure'//istr//'.txt', form='formatted',status='unknown')
                !write(11,'(1F25.15)') pressure(:,i)
                !close(11)
            deallocate(istr)
                
        end do
        print *,'Pressure of .unk files extracted'
        deallocate(Vx)
        deallocate(Vy)
        deallocate(Vz)
        deallocate(e)
        deallocate(rho)

    end subroutine ExtractPressure_3D
        
    subroutine POD()
        
        ! Variables
        implicit none
        integer :: i
        double precision, dimension(:,:), allocatable :: pressure2, modestemp, Malocal2

        ! Body of POD
        print *, 'Get POD modes and coefficients of Snapshots'
        allocate(OV%meanpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        
        ! Extract pressure of Snapshot Output file           
        if (IV%NoDim == 2) then
            call ExtractPressure(1, IV%NoSnap)
        elseif (IV%NoDim == 3) then
            call ExtractPressure_3D(1, IV%NoSnap)
        end if

        ! Exclude mean pressure from POD       
        if (IV%meanP == .true.) then
            do i = 1, IV%NoSnap
                OV%meanpressure = OV%meanpressure + OV%pressure(:,i)
            end do
            OV%meanpressure = (1/IV%NoSnap)*OV%meanpressure
            do i = 1, IV%NoSnap
                OV%pressure(:,i) = OV%pressure(:,i) - OV%meanpressure
            end do
        end if
        
        ! Perform Single Value Decomposition
        allocate(pressure2(RD%np, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD "
        pressure2 = OV%pressure
        call SVD(pressure2, size(OV%pressure, Dim = 1), size(OV%pressure, Dim = 2), modestemp)
        !deallocate(pressure2)
        print *, 'Finished SVD'
        
        OV%modes = modestemp
!!!!! OV%modes = modestemp(:,1:IV%NoPOMod)

        ! Matmul: coeff = modes'*pressure   coeff = matmul(transpose(modes),pressure)
        CALL DGEMM('T','N',size( OV%modes, dim = 2), size( OV%pressure, dim = 2), size( OV%modes, dim = 1),alpha,OV%modes,size( OV%modes, dim = 1),OV%pressure,size( OV%pressure, dim = 1),beta,OV%coeff,size( OV%coeff, dim = 1))
        
        deallocate(modestemp)
        
        !!** TESTING **!!    
        !Output: Modes and Coefficients of POD       
        !open(23,file=newdir//'/Coefficients.txt')
        !write(23,'(<IV%NoSnap>f20.7)') coeff           
        !close(23)
        !open(23,file=newdir//'/Modes.txt')
        !write(23,'(10f13.10)') modes(:,1:10)           
        !close(23)
        !pressure2(:,2) = matmul(modes,coeff(:,2))
        !open(23,file=newdir//'/Pressure_Reconstruct.txt')
        !write(23,'(1f25.10)') pressure2(:,2)           
        !close(23)
        
        if (IV%objectivefunction == 7 .or. IV%objectivefunction == 2) then
            ! Perform Single Value Decomposition
            allocate(Malocal2(size(OV%Malocal, dim = 1), size(OV%Malocal, dim = 2)),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD "
            Malocal2 = OV%Malocal
            call SVD(Malocal2, size(OV%Malocal, dim = 1), size(OV%Malocal, dim = 2), modestemp)
            print *, 'Finished SVD'
        
            OV%modes2 = modestemp
        
            ! Matmul: coeff = modes'*Malocal   coeff = matmul(transpose(modes),Malocal)
            CALL DGEMM('T','N',size( OV%modes2, dim = 2), size( OV%Malocal, dim = 2), size( OV%modes2, dim = 1),alpha,OV%modes2,size( OV%modes2, dim = 1),OV%Malocal,size( OV%Malocal, dim = 1),beta,OV%coeff2,size( OV%coeff2, dim = 1))
        end if
            
        print *, 'All Modes and Coefficients Calculated'

    end subroutine POD
        
    subroutine getDistortion(Distortion, NoSnapshot)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration
        
        ! Variables
        implicit none
        integer :: NoEngIN, NoSnapshot, j
        double precision, dimension(:), allocatable :: PlaneX, PlaneY, dP, h, Area_trap, Pmid_x, Pmid_y, Area, pT
        double precision :: Distortion, pTmean
 
        allocate(dP(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(pT(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(PlaneX(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(PlaneY(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(h(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area_trap(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_x(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_y(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area(size(OV%engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
     
        ! Body of getDistortionandPressureRecovery             
        NoEngIN = size(OV%engInNodes)
            
        ! Output: engInNodes
        PlaneX = RD%coord(OV%engInNodes,1)
        PlaneY = RD%coord(OV%engInNodes,2)
        
        !!*** Calculate coordinates of midpoints and afterwards the Area between them for average weighter pressure calculation ***!!
                                    
        ! Midpoint and Area Calculation
        do j = 1, (NoEngIN - 1)                       
            Pmid_x(j) = (PlaneX(j) + PlaneX(j+1))/2.0
            Pmid_y(j) = (PlaneY(j) + PlaneY(j+1))/2.0
        end do     
            
        do j = 1, (NoEngIN - 2)
            Area(j) = sqrt((Pmid_x(j) - Pmid_x(j+1))**2 + (Pmid_y(j) - Pmid_y(j+1))**2)
        end do  
        
        ! Total Pressure
        pT = OV%pressure(OV%engInNodes,NoSnapshot)*(1+(IV%gamma-1)*0.5*OV%MaLocal(:,NoSnapshot)**2)**(IV%gamma/(IV%gamma-1))
        
        ! Calculate Total Area Weighted Average Pressure
        pTmean = sum(pT(2:(NoEngIN-1))*Area, dim = 1)/sum(Area, dim = 1)
            
        ! Calculate Pressure Deviation
        dP = abs(pT - pTmean)           
            
        ! Determine Length and Height of Intercepting Plane
        do j = 1, (NoEngIN-1)
            h(j) = DistP2P(2, PlaneX(j), PlaneX(j+1), PlaneY(j), PlaneY(j+1))
        end do
            
        ! Apply Trapezoidal Rule to numerically integrate the Distortion
        Area_trap = h*(dP(1:(NoEngIN-1)) + dP(2:NoEngIN))/2.0
        Distortion = sum(Area_trap, dim = 1)/(pTmean*sum(h, dim = 1))
        Distortion = Distortion*(-1) ! To adapt to maximization Problem
        ! Output: Distortion constrained by a pressure recovery cap
  
    end subroutine getDistortion
        
    subroutine getengineInlet()
    ! Output: Identify all nodes positioned at the engine Inlet (engInNodes)
        
        ! Variables
        implicit none
        integer :: i,j
        integer, dimension(:,:), allocatable :: nodesall
        integer, dimension(:), allocatable :: nodesvec
        double precision, dimension(2) :: point
        
        allocate(nodesall(RD%nbf,2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getEngineInlet "
      
        ! Body of getengineInlet
        j = 0
        do i = 1, RD%nbf
            if (RD%boundtype(i) == 8) then
                j = j + 1
                nodesall(j,:) = RD%bound(i,:)           
            end if
        end do
        ! Outcome: A list of all nodes related to the engine Inlet, including possible doubling
    
        if (IV%NoDim == 2) then
            allocate(nodesvec(2*j))
            nodesvec = (/nodesall(1:j,1), nodesall(1:j,2)/)
        elseif (IV%NoDim == 3) then
            allocate(nodesvec(3*j))
            nodesvec = (/nodesall(1:j,1), nodesall(1:j,2), nodesall(1:j,3)/)
        end if
        call QSortInt(nodesvec, size(nodesvec), 'n')   
        call UniqueInt(nodesvec, size(nodesvec), OV%engInNodes)
        ! Output: unique vector engInNodes
   
    end subroutine getengineInlet
        
    subroutine SVD(A, M, N, U)

        ! Parameters
        implicit none
        integer :: M, N
        integer :: LDA, LDU, LDVT
        integer          LWMAX
        parameter        ( LWMAX = 1000000)

        ! Local Scalars
        integer          INFO, LWORK

        ! Local Arrays
        double precision, dimension(:, :), allocatable :: V2
        double precision, dimension(:,:), allocatable :: U, VT
        double precision, dimension (:), allocatable :: S
        double precision, intent(in) ::                             A( M, N )
        double precision ::                             WORK( LWMAX )
        integer, dimension(8*min(M,N)) :: IWORK    
 
        ! Executable Statements
        write(*,*)'DGESVD Program Results'

        ! Define Array Size
        allocate(V2(N,N),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD"
        LDA = M
        LDU = M
        LDVT = N
        if (M > N) then
            allocate(U(LDU,N),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD"
        else
            allocate(U(LDU,M),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD"
        end if
        allocate(VT(LDVT,N),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD"
        allocate(S(N))            
            
        ! Query the optimal workspace.
        LWORK = -1
        !call DGESDD( 'O', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
        !call DGESVD( 'N', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO ) ! Right Singular Vector
        call DGESVD( 'S', 'N', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO ) ! Left Singular Vector
        LWORK = min( LWMAX, int( WORK( 1 ) ) )

        ! Compute SVD.
        !call DGESDD( 'O', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
        !call DGESVD( 'N', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO ) ! Right Singular Vector
        call DGESVD( 'S', 'N', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO ) ! Left Singular Vector

        ! Check for convergence.
        if( INFO.GT.0 ) then
        write(*,*)'The algorithm computing SVD failed to converge.'
        stop
        end if

        V2 = transpose(VT)
            
    end subroutine SVD
    
    function LevyWalk(NoSteps, DoF)
    
        ! Variables
        implicit none
        integer :: median, scale, l, m, NoSteps, DoF
        double precision, parameter :: pi = 3.14159265359
        double precision, dimension(DoF) :: LevyWalk
        double precision, dimension(NoSteps) :: y
    
        ! Body of LevyWalk - Each Dimension walks
                
        do l = 1, DoF
      
            ! Cauchy distribution
            median = 0
            scale = 1
    
            call random_number(CS%rn)
            do m = 1, NoSteps
                y(m) = median + scale*tan(pi*CS%rn)
            end do
    
            LevyWalk(l) = sum(y, dim = 1);
        
        end do
    
    end function LevyWalk
    
    subroutine getDistortionPOD(tempNests, Distortion)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration      
        
        ! Variables
        implicit none
        integer :: NoEngIN, NoSnapshot, j
        double precision, dimension(:), allocatable :: PlaneX, PlaneY, dP, h, Area_trap, Pmid_x, Pmid_y, Area, pT
        double precision :: Distortion, pTmean
        double precision, dimension(maxDoF) :: tempNests
        double precision, dimension(:), allocatable :: newpressure, newMalocal 
 
        allocate(newpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "  
        allocate(dP(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(pT(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(PlaneX(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(PlaneY(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(h(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area_trap(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_x(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_y(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area(size(OV%engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        
        ! Body of getDistortionPOD
        NoEngIN = size(OV%engInNodes)
            
        ! Output: engInNodes
        PlaneX = RD%coord(OV%engInNodes,1)
        PlaneY = RD%coord(OV%engInNodes,2)
            
        !!*** Calculate coordinates of midpoints and afterwards the Area between them for average weighter pressure calculation ***!! 
                    
        ! Midpoint and Area Calculation
        do j = 1, (NoEngIN - 1)                       
            Pmid_x(j) = (PlaneX(j) + PlaneX(j+1))/2.0
            Pmid_y(j) = (PlaneY(j) + PlaneY(j+1))/2.0
        end do     
            
        do j = 1, (NoEngIN - 2)
            Area(j) = sqrt((Pmid_x(j) - Pmid_x(j+1))**2 + (Pmid_y(j) - Pmid_y(j+1))**2)
        end do
                
        ! Extract Pressure using POD
        call InterpolateCoefficients(tempNests, newpressure, newMalocal)
            
        ! Total Pressure
        pT = newpressure(OV%engInNodes)*(1+(IV%gamma-1)*0.5*newMalocal**2)**(IV%gamma/(IV%gamma-1))
        
        ! Calculate Total Area Weighted Average Pressure
        pTmean = sum(pT(2:(NoEngIN-1))*Area, dim = 1)/sum(Area, dim = 1)
            
        ! Calculate Pressure Deviation
        dP = abs(pT - pTmean)           
            
        ! Determine Length and Height of Intercepting Plane
        do j = 1, (NoEngIN-1)
            h(j) = DistP2P(2, PlaneX(j), PlaneX(j+1), PlaneY(j), PlaneY(j+1))
        end do
                
        ! Apply Trapezoidal Rule to numerically integrate the Distortion
        Area_trap = h*(dP(1:(NoEngIN-1)) + dP(2:NoEngIN))/2.0
        Distortion = sum(Area_trap, dim = 1)/(pTmean*sum(h, dim = 1))
        Distortion = Distortion*(-1) ! To adapt to maximization Problem
        ! Output: Distortion constrained by a pressure recovery cap
  
    end subroutine getDistortionPOD
    
    subroutine getDistortionandPressureRecoveryPOD(tempNests, Distortion,NoSnapshot)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration      
        
        ! Variables
        implicit none
        integer :: NoEngIN, NoSnapshot, j
        double precision, dimension(:), allocatable :: PlaneX, PlaneY, dP, h, Area_trap, Pmid_x, Pmid_y, Area, pT
        double precision :: Distortion, pTmean
        double precision, dimension(maxDoF) :: tempNests
        double precision, dimension(:), allocatable :: newpressure, newMalocal
 
        allocate(newpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD " 
        allocate(dP(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(pT(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(PlaneX(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(PlaneY(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(h(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area_trap(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_x(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_y(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area(size(OV%engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        
        ! Body of getDistortionPOD
        NoEngIN = size(OV%engInNodes)
            
        ! Output: engInNodes
        PlaneX = RD%coord(OV%engInNodes,1)
        PlaneY = RD%coord(OV%engInNodes,2)
            
        !!*** Calculate coordinates of midpoints and afterwards the Area between them for average weighter pressure calculation ***!! 
                    
        ! Midpoint and Area Calculation
        do j = 1, (NoEngIN - 1)                       
            Pmid_x(j) = (PlaneX(j) + PlaneX(j+1))/2.0
            Pmid_y(j) = (PlaneY(j) + PlaneY(j+1))/2.0
        end do     
            
        do j = 1, (NoEngIN - 2)
            Area(j) = sqrt((Pmid_x(j) - Pmid_x(j+1))**2 + (Pmid_y(j) - Pmid_y(j+1))**2)
        end do
                
        ! Extract Pressure using POD
        call InterpolateCoefficients(tempNests, newpressure, newMalocal)
            
        ! Total Pressure
        pT = newpressure(OV%engInNodes)*(1+(IV%gamma-1)*0.5*newMalocal**2)**(IV%gamma/(IV%gamma-1))
        
        ! Calculate Total Area Weighted Average Pressure
        pTmean = sum(pT(2:(NoEngIN-1))*Area, dim = 1)/sum(Area, dim = 1)
                   
        ! Check Pressure Recovery constraint
        OV%Precovery = (pTmean/OV%pTamb(NoSnapshot))
        
        ! Calculate Pressure Deviation
        dP = abs(pT - pTmean)           
            
        ! Determine Length and Height of Intercepting Plane
        do j = 1, (NoEngIN-1)
            h(j) = DistP2P(2, PlaneX(j), PlaneX(j+1), PlaneY(j), PlaneY(j+1))
        end do
                
        ! Apply Trapezoidal Rule to numerically integrate the Distortion
        Area_trap = h*(dP(1:(NoEngIN-1)) + dP(2:NoEngIN))/2.0
        Distortion = sum(Area_trap, dim = 1)/(pTmean*sum(h, dim = 1))
        Distortion = Distortion*(-1) ! To adapt to maximization Problem
        ! Output: Distortion constrained by a pressure recovery cap
  
        Distortion = Distortion + (OV%Precovery-1)
        
    end subroutine getDistortionandPressureRecoveryPOD
    
    subroutine ComputeRBFWeights()
    ! Objective: Compute Weights for the RBF interpolation scheme applied later in the POD reconstruction
    
        ! Variables
        implicit none
        double precision, dimension(:,:),allocatable :: RBFMatrix, A, b, coeff_temp, coeff_temp2
        double precision, dimension(:), allocatable :: PolBasis, x
        double precision :: z, ShapeParameter
        integer :: l, m, LWMAX, LWORK, Info
        parameter        ( LWMAX = 10000)
        double precision, dimension(:), allocatable ::  WORK, Ipiv
    
        ! LAPACK Parameters
        allocate(Ipiv(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        allocate(Work(LWMAX),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        LWORK = -1
        Ipiv = 0.0
        Info = 0
        ! Output Parameters
        allocate(OV%Weights(IV%NoSnap,IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "  
        ! Matrix Parameters
        allocate(RBFMatrix(IV%NoSnap,IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        allocate(PolBasis(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        if (IV%Pol == .true.) then
            allocate(OV%PolCoeff(IV%NoSnap),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
            allocate(b((IV%NoSnap + 1), 1),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights " 
            allocate(A((IV%NoSnap + 1),(IV%NoSnap + 1)),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
            allocate(x(IV%NoSnap + 1),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        else
            allocate(A(IV%NoSnap,IV%NoSnap),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
            allocate(b(IV%NoSnap, 1),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
            allocate(x(IV%NoSnap),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "      
        end if
        
        ! Body of CoefficientInterpolation
        ! A*x = b with A = [ RBFMatrix PolBasis, Polbasis^T 0] x = Weights of RBF and polynomial  b = [f 0] with f being the original coefficients of one Snapshot
        ! Details, see Computation Approach Book

        allocate(coeff_temp(size(OV%coeff, dim=1), size(OV%coeff,dim=2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        coeff_temp = OV%coeff
        coeff_temp = transpose(coeff_temp)
        
        if (IV%objectivefunction == 7 .or. IV%objectivefunction == 2) then
            allocate(OV%PolCoeff2(IV%NoSnap),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
            allocate(OV%Weights2(IV%NoSnap,IV%NoSnap),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights " 
            allocate(coeff_temp2(size(OV%coeff2, dim=1), size(OV%coeff2,dim=2)),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
            coeff_temp2 = OV%coeff2
            coeff_temp2 = transpose(coeff_temp2)
        end if 
    
        ! Generate Matrix of Radial Basis Functions (multiquadratic)
        ShapeParameter = 1.0
        do m = 1, IV%NoSnap
            do l = 1, IV%NoSnap
                z = sqrt(sum(((CS%Snapshots(m,:)-CS%Snapshots(l,:))**2), dim = 1))
                RBFMatrix(m,l) = Multiquadratic(z, ShapeParameter)
            end do 
        end do
        
        ! Design Polynomial Basis - for multiquadratic it is a '1' Matrix (see Computational Approach Book)
        PolBasis = (/ (1, m=1,IV%NoSnap) /)

        ! Assemble A matrix
        A(1:IV%NoSnap,1:IV%NoSnap) = RBFMatrix        
        if (IV%Pol == .true.) then
            A(1:IV%NoSnap,(IV%NoSnap+1)) = PolBasis
            A((IV%NoSnap+1),1:IV%NoSnap) = PolBasis
            A((IV%NoSnap+1),(IV%NoSnap+1)) = 0
        end if
      
        !Get Inverse of Matrix A via LAPACK Library
		STOP 'Ainverse does not work!'
        if (IV%Pol == .true.) then
            call dgetrf((IV%NoSnap+1), (IV%NoSnap+1), A, (IV%NoSnap+1), Ipiv, info)
            call dgetri((IV%NoSnap+1), A, (IV%NoSnap+1), Ipiv, WORK, LWORK, Info)
            LWORK = min( LWMAX, int( WORK( 1 ) ) )
            call dgetri((IV%NoSnap+1), A, (IV%NoSnap+1), Ipiv, WORK, LWORK, Info)
        else
            call dgetrf(IV%NoSnap, IV%NoSnap, A, IV%NoSnap, Ipiv, info)
            call dgetri(IV%NoSnap, A, IV%NoSnap, Ipiv, WORK, LWORK, Info)
            LWORK = min( LWMAX, int( WORK( 1 ) ) )
            call dgetri(IV%NoSnap, A, IV%NoSnap, Ipiv, WORK, LWORK, Info)
        end if

        ! Loop over all Snapshots to get the Weight Vector for each Snapshot
        do l = 1, IV%NoSnap
        
            ! b Solution/Residual Vector: Contains coefficents of ith Snapshot
            b(1:IV%NoSnap,1) = coeff_temp(:,l)
            if (IV%Pol == .true.) then
                b(1+IV%NoSnap,1) = 0
            end if
        
            ! x = A(-1)*b
            CALL DGEMM('N','N',size( A, dim = 1), size( b, dim = 2), size( A, dim = 2),alpha,A,size( A, dim = 1),b,size( b, dim = 1),beta,x,size( x, dim = 1))
        
            ! Output: RBF Weights and Polynomial Coefficient
            OV%Weights(l,:) = x(1:IV%NoSnap)          
            if (IV%Pol == .true.) then
                OV%PolCoeff(l) = x(1+IV%NoSnap)
            end if

            if (IV%objectivefunction == 7 .or. IV%objectivefunction == 2) then
                ! b Solution/Residual Vector: Contains coefficents of ith Snapshot
                b(1:IV%NoSnap,1) = coeff_temp2(:,l)
                if (IV%Pol == .true.) then
                    b(1+IV%NoSnap,1) = 0
                end if
                
                ! x = A(-1)*b
                CALL DGEMM('N','N',size( A, dim = 1), size( b, dim = 2), size( A, dim = 2),alpha,A,size( A, dim = 1),b,size( b, dim = 1),beta,x,size( x, dim = 1))
                
                ! Output: RBF Weights and Polynomial Coefficient
                OV%Weights2(l,:) = x(1:IV%NoSnap)          
                if (IV%Pol == .true.) then
                    OV%PolCoeff2(l) = x(1+IV%NoSnap)
                end if
            end if
            
        end do
        
        !open(23,file=newdir//'/Weights.txt')
        !write(23,'(10f55.10)') OV%Weights(1:10,:)           
        !close(23)
        !open(23,file=newdir//'/Polynomial.txt')
        !write(23,'(1f25.10)') OV%PolCoeff           
        !close(23)
        
        print *, 'All Weights for POD evaluated'
        
    end subroutine ComputeRBFWeights
    
    function Multiquadratic(x, ShapeParameter)
    ! Objective: Generate Function value of multiquadratic RBF function
    
        ! Variables
        implicit none
        double precision :: x, Multiquadratic, ShapeParameter
    
        ! Body of Multiquadratic
        
        if (IV%multiquadric == .true.) then
            Multiquadratic = sqrt(1+(x**2)/ShapeParameter)
        else
            Multiquadratic = exp(-(x**2)/ShapeParameter)
        end if
    
    end function Multiquadratic
    
    subroutine InterpolateCoefficients(newNest, newpressure, newMalocal)
    
        ! Variables
        implicit none
        double precision, dimension(RD%np) :: newpressure
        double precision, dimension(:), allocatable, optional :: newMalocal
        double precision, dimension(maxDoF) :: newNest
        double precision, dimension (:), allocatable :: RBFVector, NormFact
        double precision, dimension(:,:), allocatable :: newCoeff
        double precision :: x, RBFterm, ShapeParameter
        integer :: l, m
    
        allocate(RBFVector(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in InterpolateCoefficients "
        allocate(newCoeff(IV%NoSnap,1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in InterpolateCoefficients "
        allocate(NormFact(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "     
        
        ! Body of InterpolateCoefficients
        ShapeParameter = 1.0
        do l = 1, IV%NoSnap
            x = sqrt(sum(((newNest-CS%Snapshots(l,:))**2), dim = 1))
            RBFVector(l) = Multiquadratic(x, ShapeParameter)
        end do

		
		newCoeff = 0.0    
        do l = 1, IV%NoSnap
            newCoeff(:,1) = newCoeff(:,1) + OV%Weights(:,l)*RBFVector(l)
        end do
        
        if (IV%Pol == .true.) then
            newCoeff(:,1) = OV%PolCoeff + newCoeff(:,1)
        end if
            
        ! Matmul: var1 = modes*newCoeff   newpressure = matmul(modes,newCoeff(:,1))
        CALL DGEMM('N','N',size( OV%modes, dim = 1), size( newCoeff, dim = 2), size( OV%modes, dim = 2),alpha,OV%modes,size( OV%modes, dim = 1),newCoeff,size( newCoeff, dim = 1),beta,newpressure,size( newpressure, dim = 1))
       
        if (IV%meanP == .true.) then
            newpressure = newpressure + OV%meanpressure
        end if
        
        if (IV%objectivefunction == 7 .or. IV%objectivefunction == 2) then
             allocate(newMalocal(size(OV%engInNodes, dim = 1)),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
            
            do l = 1, IV%NoSnap
                RBFterm = 0
                do m = 1, IV%NoSnap
                    RBFterm = RBFterm + OV%Weights2(l,m)*RBFVector(m)
                end do

                if (IV%Pol == .true.) then
                    newCoeff(l,1) = OV%PolCoeff2(l) + RBFterm
                else
                    newCoeff(l,1) = RBFterm
                end if
            end do
        
            ! Matmul: newMalocal = modes*newCoeff   newMalocal = matmul(modes,newCoeff(:,1))
            CALL DGEMM('N','N',size( OV%modes2, dim = 1), size( newCoeff, dim = 2), size( OV%modes2, dim = 2),alpha,OV%modes2,size( OV%modes2, dim = 1),newCoeff,size( newCoeff, dim = 1),beta,newMalocal,size( newMalocal, dim = 1))
        end if
        
		open(23,file=newdir//'/POD_Pressure'//istr//'.txt')
		write(23,'(1f25.10)') newpressure           
		close(23)
        
    end subroutine InterpolateCoefficients
    
    subroutine PressInterp(DoF, newpressure, tempNests)
    ! Objective: Interpolation of Coefficients with Radial Basis Functions, based on normalized Gaussian RBF (see Hardy theory)
    
        ! Variables
        implicit none
        double precision, dimension(:,:), allocatable :: B_ar, CoeffVec, Lambda, newCoeff, InitNests, coeff_temp
        double precision, dimension(:), allocatable :: Init_Nests_temp, Ipiv
        integer, dimension(:), allocatable :: ind_IN
        double precision, dimension(maxDoF) :: tempNests
        double precision, dimension(RD%np) :: newpressure
        double precision :: a, b, d, a2
        integer :: DoF, i, l, m, n, LWMAX, LWORK, Info, tu
        parameter        ( LWMAX = 10000)
        double precision, dimension(:), allocatable ::  WORK
    
        allocate(Ipiv(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(Work(LWMAX),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(InitNests(IV%NoSnap,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(Init_Nests_temp(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(Lambda(size(OV%coeff, dim = 1),1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(CoeffVec(size(OV%coeff, dim = 1),1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(newCoeff(size(OV%coeff, dim=2),1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(ind_IN(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(coeff_temp(size(OV%coeff, dim=1), size(OV%coeff,dim=2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(B_ar(IV%NoSnap, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "

        
        ! Body of PressInterp
        ! Initial parameters for LAPACK
        LWORK = -1
        Ipiv = 0.0
        Info = 0
        
        ! Sort InitNests Matrix
        InitNests = CS%Snapshots
        coeff_temp = OV%coeff
        ind_IN = (/ (i, i=1,IV%NoSnap) /)      
        Init_Nests_temp = InitNests(:, (1+maxDoF-DoF)) ! only for QSort
        call QSort(Init_Nests_temp, size(InitNests, dim = 1), 'y', ind_IN) ! Index only required
        deallocate(Init_Nests_temp)
            
        InitNests = InitNests(ind_IN,:)
        coeff_temp = coeff_temp(:,ind_IN)
        coeff_temp = transpose(coeff_temp) 
        
        ! Compute Shape Parameter
        a = 0
        b = 0
        do l = 1, (IV%NoSnap - 1)
            do m = (l+1), IV%NoSnap
                a = a + sqrt(sum(((InitNests(l,:)-InitNests(m,:))**2), dim = 1))    ! Sum of all distances between points
                b = b + 1                                                           ! Amount of points considered  
            end do    
        end do
        d = a/b             ! Mean Distance between Points
        d = 1 !0.25*(d**2)
        
        
        ! Compute Coefficients for Interpolation
        do m = 1, IV%NoSnap
            do n = 1, IV%NoSnap       
                a = sqrt(sum(((InitNests(m,:)-InitNests(n,:))**2), dim = 1)) ! Distance
                B_ar(m,n) = 1.0/d*exp(-(a**2)/d)           
            end do     
        end do
            
        !Get Inverse of Matrix via LAPACK Library
        call dgetrf(IV%NoSnap, IV%NoSnap, B_ar, IV%NoSnap, Ipiv, info)
        call dgetri(IV%NoSnap, B_ar, IV%NoSnap, Ipiv, WORK, LWORK, Info)
        LWORK = min( LWMAX, int( WORK( 1 ) ) )
        call dgetri(IV%NoSnap, B_ar, IV%NoSnap, Ipiv, WORK, LWORK, Info)
              
        ! For each 'pressure field' f(:,k) ie vector of coefficients corresponding to mode k        
        newCoeff(:,1) = (/ (0, i=1,IV%NoSnap) /)
        do l = 1, size(coeff_temp, dim = 2)
                                  
            CoeffVec(:,1) = coeff_temp(:,l)
            ! Matmul: Lambda = B_ar*CoeffVec   Lambda = matmul(B_ar,CoeffVec)           
            CALL DGEMM('N','N',size( B_ar, dim = 1), size( CoeffVec, dim = 2), size( B_ar, dim = 2),alpha,B_ar,size( B_ar, dim = 1),CoeffVec,size( CoeffVec, dim = 1),beta,Lambda,size( Lambda, dim = 1))
            
            do m = 1, IV%NoSnap
                a = sqrt(sum(((tempNests-InitNests(m,:))**2), dim = 1)) ! Distance
                newCoeff(l,1) = newCoeff(l,1) + Lambda(m,1)*(1.0/d*exp(-(a**2)/d))
            end do

        end do        
        ! Matmul: var1 = modes*newCoeff   newpressure = matmul(modes,newCoeff(:,1))
        CALL DGEMM('N','N',size( OV%modes, dim = 1), size( newCoeff, dim = 2), size( OV%modes, dim = 2),alpha,OV%modes,size( OV%modes, dim = 1),newCoeff,size( newCoeff, dim = 1),beta,newpressure,size( newpressure, dim = 1))
        
    end subroutine PressInterp
    
    subroutine AllocateModesCoeff()
    
        ! Variables
        implicit none
    
        ! Body of AllocateModesCoeff
        if (IV%NoPOMod < 0 .OR. IV%NoPOMod > IV%NoSnap) then
            allocate(OV%modes(RD%np, IV%NoSnap),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            allocate(OV%coeff(IV%NoSnap, IV%NoSnap),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            if (IV%objectivefunction == 7 .or. IV%objectivefunction == 2) then
                allocate(OV%modes2(size(OV%engInNodes, dim = 1), IV%NoSnap),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
                allocate(OV%coeff2(IV%NoSnap, IV%NoSnap),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            end if
        else
            allocate(OV%modes(RD%np, IV%NoPOMod),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            allocate(OV%coeff(IV%NoSnap, IV%NoPOMod),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            if (IV%objectivefunction == 7 .or. IV%objectivefunction == 2) then
                allocate(OV%modes2(size(OV%engInNodes, dim = 1), IV%NoPOMod),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
                allocate(OV%coeff2(IV%NoSnap, IV%NoPOMod),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            end if
        end if
    
    end subroutine AllocateModesCoeff
    
    subroutine getObjectiveFunction(PODevaluation, Fi, tempNests, NoSnapshot)
    
        ! Variables
        implicit none
        double precision, dimension(maxDoF), optional :: tempNests
        integer, optional :: NoSnapshot
        double precision :: Fi
        logical :: PODevaluation
    
        ! Body of getObjectiveFunction
        if (PODevaluation == .true.) then
            if (IV%ObjectiveFunction == 2) then
                call getDistortionPOD(tempNests, Fi)
            elseif (IV%ObjectiveFunction == 1) then
                allocate(RD%coord_temp(RD%np,IV%NoDim),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main "
                RD%coord_temp = RD%coord 
                call SubMovemesh(tempNests)       
                call getLiftandDragPOD(tempNests, Fi)
                deallocate(RD%coord_temp)
            elseif (IV%ObjectiveFunction == 3) then
                allocate(RD%coord_temp(RD%np,IV%NoDim),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main "
                RD%coord_temp = RD%coord 
                call SubMovemesh(tempNests)       
                call getmaxLiftPOD(tempNests, Fi)
                deallocate(RD%coord_temp)
            elseif (IV%ObjectiveFunction == 6) then
                allocate(RD%coord_temp(RD%np,IV%NoDim),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main "
                RD%coord_temp = RD%coord 
                call SubMovemesh(tempNests)
                call getzeroLiftPOD(tempNests, Fi)
                deallocate(RD%coord_temp)
            elseif (IV%ObjectiveFunction == 7) then
                call getDistortionandPressureRecoveryPOD(tempNests, Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 8) then
                call getPressureDistributionPOD(Fi, tempNests)
            elseif (IV%ObjectiveFunction == 97) then
                call getAckley(Fi, tempNests)
            elseif (IV%ObjectiveFunction == 98) then
                call getDeJong(Fi, tempNests)
            elseif (IV%ObjectiveFunction == 99) then
                call getRosenbrock(Fi, tempNests)
            end if
        else
            if (OV%converged(NoSnapshot) == 0) then
                Fi = -150
                RETURN
            end if
            if (IV%ObjectiveFunction == 2) then
                call getDistortion(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 1) then
                call getLiftandDrag(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 3) then
                call getmaxLift(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 4) then
                call getminDrag(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 5) then
                call getDownForce(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 6) then
                call getzeroLift(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 7) then
                call getDistortionandPressureRecovery(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 8) then
                call getPressureDistribution(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 9) then
                call getL1AreaError(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 97) then
                call getAckley(Fi, tempNests)
            elseif (IV%ObjectiveFunction == 98) then
                call getDeJong(Fi, tempNests)
            elseif (IV%ObjectiveFunction == 99) then
                call getRosenbrock(Fi, tempNests)
            end if
        end if
    
    end subroutine getObjectiveFunction
    
    subroutine getmaxLiftPOD(tempNests, Fi)
    
        ! Variables
        implicit none
        double precision :: Lift, Fi, p, dx, norm, ralpha
        double precision, dimension(2) :: n, m, k, tangent, vec1, vec2
        integer :: intexit, i, j, NoIB
        double precision, dimension(:,:), allocatable :: tangentArray
        double precision, dimension(:), allocatable :: lengthArray
        double precision, dimension(:), allocatable :: newpressure
        double precision, dimension(maxDoF) :: tempNests
        double precision, parameter :: pi = 3.14159265359
  
        ! Body of getmaxLift
        NoIB = size(Innerbound)
        allocate(newpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(tangentArray(NoIB, IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(lengthArray(NoIB),stat=allocateStatus)
        
        ! Calculate Tangents and face lengths for each boundary point
        tangentArray = 0.0
        lengthArray = 0.0
        do i = 1, RD%nbf
            vec1 = RD%coord_temp(RD%bound(i,1),:)
            vec2 = RD%coord_temp(RD%bound(i,2),:)
            tangent = vec2 - vec1
            dx = DistP2P(IV%NoDim, vec1(1), vec2(1), vec1(2), vec2(2))
            intexit = 0
            do j = 1, NoIB
                if (RD%bound(i,1) == Innerbound(j)) then
                    tangentArray(j,:) = tangentArray(j,:) + tangent
                    lengthArray(j) = lengthArray(j) + dx
                    intexit = intexit + 1
                    if (intexit == 2) then
                        EXIT
                    end if
                end if
                if (RD%bound(i,2) == Innerbound(j)) then
                    tangentArray(j,:) = tangentArray(j,:) + tangent
                    lengthArray(j) = lengthArray(j) + dx
                    intexit = intexit + 1
                    if (intexit == 2) then
                        EXIT
                    end if
                end if
            end do
        end do
        lengthArray = lengthArray*0.5

        ! normalize tangents
        do i = 1, NoIB
            norm = sqrt(sum(tangentArray(i,:)**2))
            tangentArray(i,:) = tangentArray(i,:)/norm
        end do

        ! calculate Lift and Drag
        call InterpolateCoefficients(tempNests, newpressure)
        ralpha = IV%AlphaInflowDirection*Pi/180
        m = (/cos(ralpha), sin(ralpha)/)  ! tangent to flow direction
        k = (/- m(2), m(1)/)                                                ! normal of flow direction
        Lift = 0.0
        do i = 2, NoIB
            
            ! pressure portion
            n = (/ -tangentArray(i,2), tangentArray(i,1)/)                 ! normal to tangent, pointing into the geometry
            p = newpressure(Innerbound(i))
            dx = lengthArray(i)

            Lift = Lift - 2.0*p*sum(n*k)*dx

            ! skin friction(viscosity/boundary layer) portion
            ! later
            
        end do
        Fi = Lift
        
    end subroutine getmaxLiftPOD
    subroutine getLiftandDragPOD(tempNests, Fi)
    
        ! Variables
        implicit none
        double precision :: Lift, Drag, Fi, p, dx, norm, ralpha
        double precision, dimension(2) :: n, m, k, tangent, vec1, vec2
        integer :: intexit, i, j, NoIB
        double precision, dimension(:,:), allocatable :: tangentArray
        double precision, dimension(:), allocatable :: lengthArray
        double precision, dimension(:), allocatable :: newpressure
        double precision, dimension(maxDoF) :: tempNests
        double precision, parameter :: pi = 3.14159265359
  
        ! Body of getLiftandDrag
        NoIB = size(Innerbound)
        allocate(newpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(tangentArray(NoIB, IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(lengthArray(NoIB),stat=allocateStatus)
        
        ! Calculate Tangents and face lengths for each boundary point
        tangentArray = 0.0
        lengthArray = 0.0
        do i = 1, RD%nbf
            vec1 = RD%coord_temp(RD%bound(i,1),:)
            vec2 = RD%coord_temp(RD%bound(i,2),:)
            tangent = vec2 - vec1
            dx = DistP2P(IV%NoDim, vec1(1), vec2(1), vec1(2), vec2(2))
            intexit = 0
            do j = 1, NoIB
                if (RD%bound(i,1) == Innerbound(j)) then
                    tangentArray(j,:) = tangentArray(j,:) + tangent
                    lengthArray(j) = lengthArray(j) + dx
                    intexit = intexit + 1
                    if (intexit == 2) then
                        EXIT
                    end if
                end if
                if (RD%bound(i,2) == Innerbound(j)) then
                    tangentArray(j,:) = tangentArray(j,:) + tangent
                    lengthArray(j) = lengthArray(j) + dx
                    intexit = intexit + 1
                    if (intexit == 2) then
                        EXIT
                    end if
                end if
            end do
        end do
        lengthArray = lengthArray*0.5

        ! normalize tangents
        do i = 1, NoIB
            norm = sqrt(sum(tangentArray(i,:)**2))
            tangentArray(i,:) = tangentArray(i,:)/norm
        end do

        ! calculate Lift and Drag
        call InterpolateCoefficients(tempNests, newpressure)
        ralpha = IV%AlphaInflowDirection*Pi/180
        m = (/cos(ralpha), sin(ralpha)/)  ! tangent to flow direction
        k = (/- m(2), m(1)/)                                                ! normal of flow direction
        Lift = 0.0
        Drag = 0.0
        do i = 2, NoIB
            
            ! pressure portion
            n = (/ -tangentArray(i,2), tangentArray(i,1)/)                 ! normal to tangent, pointing into the geometry
            p = newpressure(Innerbound(i))
            dx = lengthArray(i)

            Lift = Lift - 2.0*p*sum(n*k)*dx
            Drag = Drag - 2.0*p*sum(n*m)*dx

            ! skin friction(viscosity/boundary layer) portion
            ! later
            
        end do
        Fi = Lift/Drag
        
    end subroutine getLiftandDragPOD
    
    subroutine getLiftandDrag(Fi, NoSnapshot)
    
        ! Variables
        implicit none
        integer :: FileSize, LastLine, NoSnapshot, j
        double precision, dimension(8) :: Input      
        double precision :: Lift, Drag, Fi
       
        ! Body of getLiftandDrag
        call DetermineStrLen(istr, NoSnapshot)
        open(11, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', form='formatted',status='old')
        deallocate(istr)
        inquire(11, size = FileSize)           
        if (IV%NoDim == 3) then
            LastLine = FileSize/175
        else
            if (IV%SystemType == 'W') then
                LastLine = FileSize/106 !107
            else     
                LastLine = FileSize/106
            end if
        end if
            
        ! Read until last line
        do j = 1, (LastLine - 1)
            read(11, *) Input
        end do
        read(11, *) Input
        close(11)
        Lift = Input(3)
        Drag = Input(4)
        Fi = Lift/Drag
        
    end subroutine getLiftandDrag
    
    subroutine getzeroLiftPOD(tempNests, Fi)
    
        ! Variables
        implicit none
        double precision :: Lift, Drag, Fi, p, dx, norm, ralpha
        double precision, dimension(2) :: n, m, k, tangent, vec1, vec2
        integer :: intexit, i, j, NoIB
        double precision, dimension(:,:), allocatable :: tangentArray
        double precision, dimension(:), allocatable :: lengthArray
        double precision, dimension(:), allocatable :: newpressure
        double precision, dimension(maxDoF) :: tempNests
        double precision, parameter :: pi = 3.14159265359
  
        ! Body of getzeroLiftPOD
        NoIB = size(Innerbound)
        allocate(newpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(tangentArray(NoIB, IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(lengthArray(NoIB),stat=allocateStatus)
        
        ! Calculate Tangents and face lengths for each boundary point
        tangentArray = 0.0
        lengthArray = 0.0
        do i = 1, RD%nbf
            vec1 = RD%coord_temp(RD%bound(i,1),:)
            vec2 = RD%coord_temp(RD%bound(i,2),:)
            tangent = vec2 - vec1
            dx = DistP2P(IV%NoDim, vec1(1), vec2(1), vec1(2), vec2(2))
            intexit = 0
            do j = 1, NoIB
                if (RD%bound(i,1) == Innerbound(j)) then
                    tangentArray(j,:) = tangentArray(j,:) + tangent
                    lengthArray(j) = lengthArray(j) + dx
                    intexit = intexit + 1
                    if (intexit == 2) then
                        EXIT
                    end if
                end if
                if (RD%bound(i,2) == Innerbound(j)) then
                    tangentArray(j,:) = tangentArray(j,:) + tangent
                    lengthArray(j) = lengthArray(j) + dx
                    intexit = intexit + 1
                    if (intexit == 2) then
                        EXIT
                    end if
                end if
            end do
        end do
        lengthArray = lengthArray*0.5

        ! normalize tangents
        do i = 1, NoIB
            norm = sqrt(sum(tangentArray(i,:)**2))
            tangentArray(i,:) = tangentArray(i,:)/norm
        end do

        ! calculate Lift and Drag
        call InterpolateCoefficients(tempNests, newpressure)
        ralpha = IV%AlphaInflowDirection*Pi/180
        k = (/- m(2), m(1)/)                                                ! normal of flow direction 
        do i = 2, NoIB
            
            ! pressure portion
            n = (/ -tangentArray(i,2), tangentArray(i,1)/)                 ! normal to tangent, pointing into the geometry
            p = newpressure(Innerbound(i))
            dx = lengthArray(i)

            Lift = Lift - 2.0*p*sum(n*k)*dx

            ! skin friction(viscosity/boundary layer) portion
            ! later
            
        end do
        Fi = abs(Lift)*(-1)
        
    end subroutine getzeroLiftPOD
    
    subroutine getzeroLift(Fi, NoSnapshot)
    
        ! Variables
        implicit none
        integer :: FileSize, LastLine, NoSnapshot, j
        double precision, dimension(8) :: Input      
        double precision :: Lift, Drag, Fi
       
        ! Body of getzeroLift
        call DetermineStrLen(istr, NoSnapshot)
        open(11, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', form='formatted',status='old')
        deallocate(istr)
        inquire(11, size = FileSize)           
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
            read(11, *) Input
        end do
        read(11, *) Input
        close(11)
        Lift = Input(3)
        Fi = abs(Lift)*(-1)
        
    end subroutine getzeroLift
    
    subroutine getDownForce(Fi, NoSnapshot)
    
        ! Variables
        implicit none
        integer :: FileSize, LastLine, NoSnapshot, j
        double precision, dimension(8) :: Input      
        double precision :: Lift, Fi
       
        ! Body of getLiftandDrag
        call DetermineStrLen(istr, NoSnapshot)
        open(11, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', form='formatted',status='old')
        deallocate(istr)
        inquire(11, size = FileSize)           
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
            read(11, *) Input
        end do
        read(11, *) Input
        close(11)
        Lift = Input(3)
        
        Fi = -Lift       
        
    end subroutine getDownForce
    
    subroutine getmaxLift(Fi, NoSnapshot)
    
        ! Variables
        implicit none
        integer :: FileSize, LastLine, NoSnapshot, j
        double precision, dimension(8) :: Input      
        double precision :: Lift, Fi
       
        ! Body of getzeroLift
        call DetermineStrLen(istr, NoSnapshot)
        open(11, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', form='formatted',status='old')
        deallocate(istr)
        inquire(11, size = FileSize)           
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
            read(11, *) Input
        end do
        read(11, *) Input
        close(11)
        Lift = Input(3)
        Fi = Lift
        
    end subroutine getmaxLift
    
    subroutine getminDrag(Fi, NoSnapshot)
    
        ! Variables
        implicit none
        integer :: FileSize, LastLine, NoSnapshot, j
        double precision, dimension(8) :: Input      
        double precision :: Drag, Fi
       
        ! Body of getzeroLift
        call DetermineStrLen(istr, NoSnapshot)
        open(11, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', form='formatted',status='old')
        deallocate(istr)
        inquire(11, size = FileSize)           
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
            read(11, *) Input
        end do
        read(11, *) Input
        close(11)
        Drag = Input(4)
        Fi = -Drag
        
    end subroutine getminDrag
    
    subroutine getDistortionandPressureRecovery(Distortion, NoSnapshot)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration
    
        ! Variables
        implicit none
        integer :: NoEngIN, NoSnapshot, j
        double precision, dimension(:), allocatable :: PlaneX, PlaneY, dP, h, Area_trap, Pmid_x, Pmid_y, Area, pT
        double precision :: Distortion, pTmean
 
        allocate(PlaneX(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(PlaneY(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(dP(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(h(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area_trap(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_x(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_y(size(OV%engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area(size(OV%engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(pT(size(OV%engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
     
        ! Body of getDistortionandPressureRecovery             
        NoEngIN = size(OV%engInNodes)
            
        ! Output: engInNodes
        PlaneX = RD%coord(OV%engInNodes,1)
        PlaneY = RD%coord(OV%engInNodes,2)
        
        !!*** Calculate coordinates of midpoints and afterwards the Area between them for average weighter pressure calculation ***!!
                                    
        ! Midpoint and Area Calculation
        do j = 1, (NoEngIN - 1)                       
            Pmid_x(j) = (PlaneX(j) + PlaneX(j+1))/2.0
            Pmid_y(j) = (PlaneY(j) + PlaneY(j+1))/2.0
        end do     
            
        do j = 1, (NoEngIN - 2)
            Area(j) = sqrt((Pmid_x(j) - Pmid_x(j+1))**2 + (Pmid_y(j) - Pmid_y(j+1))**2)
        end do  
        
        ! Total Pressure
        pT = OV%pressure(OV%engInNodes,NoSnapshot)*(1+(IV%gamma-1)*0.5*OV%MaLocal(:,NoSnapshot)**2)**(IV%gamma/(IV%gamma-1))
        
        ! Calculate Total Area Weighted Average Pressure
        pTmean = sum(pT(2:(NoEngIN-1))*Area, dim = 1)/sum(Area, dim = 1)
           
        ! Check Pressure Recovery constraint
        OV%Precovery = (pTmean/OV%pTamb(NoSnapshot))
            
        ! Calculate Pressure Deviation
        dP = abs(pT - pTmean)           
            
        ! Determine Length and Height of Intercepting Plane
        do j = 1, (NoEngIN-1)
            h(j) = DistP2P(2, PlaneX(j), PlaneX(j+1), PlaneY(j), PlaneY(j+1))
        end do
            
        ! Apply Trapezoidal Rule to numerically integrate the Distortion
        Area_trap = h*(dP(1:(NoEngIN-1)) + dP(2:NoEngIN))/2.0
        Distortion = sum(Area_trap, dim = 1)/(pTmean*sum(h, dim = 1))
        Distortion = Distortion*(-1) ! To adapt to maximization Problem
        ! Output: Distortion constrained by a pressure recovery cap
            
        Distortion = Distortion + (OV%Precovery-1)
        
    end subroutine getDistortionandPressureRecovery
    
    subroutine getL1AreaError(L1_Error, NoSnapshot)
    
        ! Variables
        implicit none
        integer :: i, j, k, nbp_target, maxbp, NoSnapshot, nibp
        logical :: ex
        double precision :: A, B, C, Atarget, Btarget, Ctarget, detAB, detCB, detAC, xtemp, ytemp, integral, integral_target, L1_Error
        double precision, dimension(:), allocatable :: x, y, xtarget, ytarget, xintersect, yintersect, xintersect2, yintersect2
        integer, dimension(:), allocatable :: ind, ind_target, indices, indices_target, order, marker
    
        ! Body of getL1AreaError
        
        ! Read in Target Reference Data
        nibp = size(InnerBound)
        inquire(file=DataFolder//'/Geom_target.txt', exist = ex)
        if (ex == .true.) then
            open(29,file=DataFolder//'/Geom_target.txt',form='formatted',status='old')
        else
            STOP "The file 'Geom_target' does not exist. Please generate!"
        end if
        read(29,*) nbp_target
        allocate(ind_target(nbp_target),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(xtarget(nbp_target),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(ytarget(nbp_target),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        do i = 1, nbp_target
            read(29,*) xtarget(i), ytarget(i)
        end do
        close(29)
        
        ! Allocate Arrays
        allocate(ind(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(x(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(y(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        
        ! Extract coordinates
        call readDatFile(NoSnapshot)
        x = RD%coord_temp(orderedBoundaryIndex,1)
        y = RD%coord_temp(orderedBoundaryIndex,2)
        deallocate(RD%coord_temp)
        
        if (nibp > nbp_target) then
            maxbp = nibp
        else
            maxbp = nbp_target
        end if
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(xintersect(maxbp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(yintersect(maxbp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        
        ! Identify intersections between target and current boundary
        k = 0
        do j = 1, nibp-1
            do i = 1, nbp_target-1
                
                ! Calculate intersection point
                A = y(j+1) - y(j)
                B = x(j) - x(j+1)
                C = A*x(j) + B*y(j)
                
                Atarget = ytarget(i+1) - ytarget(i)
                Btarget = xtarget(i) - xtarget(i+1)
                Ctarget = Atarget*xtarget(i) + Btarget*ytarget(i)
                
                detAB = A*Btarget - Atarget*B
                detCB = C*Btarget - Ctarget*B
                detAC = A*Ctarget - Atarget*C
                if (abs(detAB) < abs(10e-12)) then
                    ! Parallel
                    xtemp = -1000
                    ytemp = -1000
                else
                    ! Non-parallel
                    xtemp = detCB/detAB
                    ytemp = detAC/detAB
                end if
                
                ! Check, if intersection point lays on the line segment considered 
                if ((xtemp - min(x(j),x(j+1))) > -10e-12 .and. (xtemp - max(x(j),x(j+1))) < 10e-12 .and. (ytemp - min(y(j),y(j+1))) > -10e-12  .and. (ytemp - max(y(j),y(j+1))) < 10e-12)  then
                    if ((xtemp - min(xtarget(i),xtarget(i+1))) > -10e-12 .and. (xtemp - max(xtarget(i),xtarget(i+1))) < 10e-12  .and. (ytemp - min(ytarget(i),ytarget(i+1))) > -10e-12  .and. (ytemp - max(ytarget(i),ytarget(i+1))) < 10e-12 )  then
                        k = k + 1
                        ind_target(k) = i
                        ind(k) = j
                        xintersect(k) = xtemp
                        yintersect(k) = ytemp
                    end if
                end if
            end do
        end do

        if (k == 0) then
            
            ! Integrate over both boundaries using trapezoidal rule
            integral = 0
            integral_target = 0
            L1_error = 0
            
            do j = 1, nibp-1 
                integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
            end do
            integral = integral + 0.5*((x(1) - x(nibp))*(y(1) + y(nibp)))
            do j = 1, nbp_target-1 
                integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(ytarget(j+1) + ytarget(j)))
            end do
            integral_target = integral_target + 0.5*((xtarget(1) - xtarget(nbp_target))*(ytarget(1) + ytarget(nbp_target)))
            
        else
        
            ! Identify doubled points
            allocate(marker(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
            marker = 0
            j = k
            do i = 2, j
                if (xintersect(i) == xintersect(i-1) .and. yintersect(i) == yintersect(i-1)) then
                    marker(i) = 1
                    k = k - 1
                end if
            end do
        
            ! Array allocation with new array size j
            allocate(xintersect2(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "             
            allocate(yintersect2(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            allocate(indices(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            allocate(indices_target(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
        
            ! Exclude all points identified to be existing twice
            k = 0
            do i = 1, j           
                if (marker(i) == 0) then
                    k = k + 1
                    xintersect2(k) = xintersect(i)
                    yintersect2(k) = yintersect(i)
                    indices(k) = ind(i)
                    indices_target(k) = ind_target(i)
                end if
            end do
        
            ! Hand over of points to previous array name
            deallocate(xintersect)
            allocate(xintersect(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            xintersect = xintersect2
            deallocate(xintersect2)
            deallocate(yintersect)
            allocate(yintersect(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            yintersect = yintersect2
            deallocate(yintersect2)
            deallocate(ind)
            deallocate(ind_target)
        
            ! Integrate over both boundaries using trapezoidal rule
            integral = 0
            integral_target = 0
            L1_error = 0
            do i = 1, k-1
            
                integral = integral + 0.5*((x(indices(i)+1) - xintersect(i))*(yintersect(i) + y(indices(i)+1)))
                do j = indices(i)+1, indices(i+1)-1
                    integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
                end do
                integral = integral + 0.5*((xintersect(i+1) - x(indices(i+1)))*(yintersect(i+1) + y(indices(i+1))))
            
                integral_target = integral_target + 0.5*((xtarget(indices_target(i)+1) - xintersect(i))*(yintersect(i) + ytarget(indices_target(i)+1)))
                do j = indices_target(i)+1, indices_target(i+1)-1 
                    integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(ytarget(j+1) + ytarget(j)))
                end do
                integral_target = integral_target + 0.5*((xintersect(i+1) - xtarget(indices_target(i+1)))*(yintersect(i+1) + ytarget(indices_target(i+1))))
            
                ! L1 Absolute Error norm (root mean square)
                L1_error = L1_error + abs(integral - integral_target)
                integral = 0
                integral_target = 0
            
            end do
       
            integral = integral + 0.5*((x(indices(k)+1) - xintersect(k))*(yintersect(k) + y(indices(k)+1)))
            do j = indices(k)+1, nibp-1 
                integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
            end do
            integral = integral + 0.5*((x(1) - x(nibp))*(y(1) + y(nibp)))
            do j = 1, indices(1)-1
                integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
            end do
            integral = integral + 0.5*((xintersect(1) - x(indices(1)))*(yintersect(1) + y(indices(1))))
            
            integral_target = integral_target + 0.5*((xtarget(indices_target(k)+1) - xintersect(k))*(yintersect(k) + ytarget(indices_target(k)+1)))
            do j = indices_target(k)+1, nbp_target-1 
                integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(ytarget(j+1) + ytarget(j)))
            end do
            integral_target = integral_target + 0.5*((xtarget(1) - xtarget(nbp_target))*(ytarget(1) + ytarget(nbp_target)))
            do j = 1, indices_target(1)-1
                integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(ytarget(j+1) + ytarget(j)))
            end do
            integral_target = integral_target + 0.5*((xintersect(1) - xtarget(indices_target(1)))*(yintersect(1) + ytarget(indices_target(1))))
        
        end if
        
        ! L1 Absolute Error norm
        L1_error = L1_error + abs(integral - integral_target)
        L1_error = -L1_error
        
    end subroutine getL1AreaError
    
    subroutine getPressureDistribution(L1_error, NoSnapshot)
        
        ! Variables
        implicit none
        integer :: i, j, k, nbp_target, maxbp, NoSnapshot, nibp
        logical :: ex
        double precision :: A, B, C, Atarget, Btarget, Ctarget, detAB, detCB, detAC, xtemp, ytemp, integral, integral_target, L1_Error, p0
        double precision, dimension(:), allocatable :: x, y, xtarget, Cptarget, xintersect, yintersect, xintersect2, yintersect2
        integer, dimension(:), allocatable :: ind, ind_target, indices, indices_target, order, marker
    
        ! Body of getPressureDistribution
        nibp = size(InnerBound)

        ! Read in Target Reference Data         
        inquire(file=DataFolder//'/Cp_target_Onera206.txt', exist = ex)
        if (ex == .true.) then
            open(29,file=DataFolder//'/Cp_target_Onera206.txt',form='formatted',status='old')
        else
            STOP "The file 'Cp_target' does not exist. Please generate!"
        end if

        read(29,*) nbp_target
        allocate(ind_target(nbp_target),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(xtarget(nbp_target),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(Cptarget(nbp_target),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        do i = 1, nbp_target
            read(29,*) xtarget(i), Cptarget(i)
        end do
        close(29)
       
        ! Allocate Arrays
        allocate(ind(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(x(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(y(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "

        ! Extract x coordinate
        call readDatFile(NoSnapshot)
        x = RD%coord_temp(orderedBoundaryIndex,1)
        deallocate(RD%coord_temp)
        
        ! Calculate y(pressure coefficient)
        p0 = 0.178569528995839 !Ma2.0: 0.178569528995839  !Ma0.5: 2.857162211411720 !Ma0.729: 1.344037930148188 ! Soft code later, this is specific for Re = 6.5e7, T = 255.658 and l = 1 feet
        y = (OV%pressure(orderedBoundaryIndex, NoSnapshot)/p0 - 1)/(0.5*IV%gamma*IV%Ma**2)
          
        if (nibp > nbp_target) then
            maxbp = nibp
        else
            maxbp = nbp_target
        end if
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(xintersect(maxbp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(yintersect(maxbp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        
        ! Identify intersections between target and current boundary
        k = 0
        do j = 1, nibp-1
            do i = 1, nbp_target-1
                
                ! Calculate intersection point
                A = y(j+1) - y(j)
                B = x(j) - x(j+1)
                C = A*x(j) + B*y(j)
                
                Atarget = Cptarget(i+1) - Cptarget(i)
                Btarget = xtarget(i) - xtarget(i+1)
                Ctarget = Atarget*xtarget(i) + Btarget*Cptarget(i)
                
                detAB = A*Btarget - Atarget*B
                detCB = C*Btarget - Ctarget*B
                detAC = A*Ctarget - Atarget*C
                if (abs(detAB) < abs(10e-12)) then
                    ! Parallel
                    xtemp = -1000
                    ytemp = -1000
                else
                    ! Non-parallel
                    xtemp = detCB/detAB
                    ytemp = detAC/detAB
                end if
                
                ! Check, if intersection point lays on the line segment considered 
                if ((xtemp - min(x(j),x(j+1))) > -10e-12 .and. (xtemp - max(x(j),x(j+1))) < 10e-12 .and. (ytemp - min(y(j),y(j+1))) > -10e-12  .and. (ytemp - max(y(j),y(j+1))) < 10e-12)  then
                    if ((xtemp - min(xtarget(i),xtarget(i+1))) > -10e-12 .and. (xtemp - max(xtarget(i),xtarget(i+1))) < 10e-12  .and. (ytemp - min(Cptarget(i),Cptarget(i+1))) > -10e-12  .and. (ytemp - max(Cptarget(i),Cptarget(i+1))) < 10e-12 )  then
                        k = k + 1
                        ind_target(k) = i
                        ind(k) = j
                        xintersect(k) = xtemp
                        yintersect(k) = ytemp
                    end if
                end if
            end do
        end do
        
        if (k == 0) then
            
            ! Integrate over both boundaries using trapezoidal rule
            integral = 0
            integral_target = 0
            L1_error = 0
            
            do j = 1, nibp-1 
                integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
            end do
            integral = integral + 0.5*((x(1) - x(nibp))*(y(1) + y(nibp)))
            do j = 1, nbp_target-1 
                integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(Cptarget(j+1) + Cptarget(j)))
            end do
            integral_target = integral_target + 0.5*((xtarget(1) - xtarget(nbp_target))*(Cptarget(1) + Cptarget(nbp_target)))
            
        else
           
            ! Identify doubled points
            allocate(marker(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
            marker = 0
            j = k
            do i = 2, j
                if (xintersect(i) == xintersect(i-1) .and. yintersect(i) == yintersect(i-1)) then
                    marker(i) = 1
                    k = k - 1
                end if
            end do
        
            ! Array allocation with new array size j
            allocate(xintersect2(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "             
            allocate(yintersect2(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            allocate(indices(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            allocate(indices_target(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
        
            ! Exclude all points identified to be existing twice
            k = 0
            do i = 1, j           
                if (marker(i) == 0) then
                    k = k + 1
                    xintersect2(k) = xintersect(i)
                    yintersect2(k) = yintersect(i)
                    indices(k) = ind(i)
                    indices_target(k) = ind_target(i)
                end if
            end do
        
            ! Hand over of points to previous array name
            deallocate(xintersect)
            allocate(xintersect(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            xintersect = xintersect2
            deallocate(xintersect2)
            deallocate(yintersect)
            allocate(yintersect(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            yintersect = yintersect2
            deallocate(yintersect2)
            deallocate(ind)
            deallocate(ind_target)
        
            ! Integrate over both boundaries using trapezoidal rule
            integral = 0
            integral_target = 0
            L1_error = 0
            do i = 1, k-1
            
                integral = integral + 0.5*((x(indices(i)+1) - xintersect(i))*(yintersect(i) + y(indices(i)+1)))
                do j = indices(i)+1, indices(i+1)-1
                    integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
                end do
                integral = integral + 0.5*((xintersect(i+1) - x(indices(i+1)))*(yintersect(i+1) + y(indices(i+1))))
            
                integral_target = integral_target + 0.5*((xtarget(indices_target(i)+1) - xintersect(i))*(yintersect(i) + Cptarget(indices_target(i)+1)))
                do j = indices_target(i)+1, indices_target(i+1)-1 
                    integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(Cptarget(j+1) + Cptarget(j)))
                end do
                integral_target = integral_target + 0.5*((xintersect(i+1) - xtarget(indices_target(i+1)))*(yintersect(i+1) + Cptarget(indices_target(i+1))))
            
                ! L1 Absolute Error norm (root mean square)
                L1_error = L1_error + abs(integral - integral_target)
                integral = 0
                integral_target = 0
            
            end do
       
            integral = integral + 0.5*((x(indices(k)+1) - xintersect(k))*(yintersect(k) + y(indices(k)+1)))
            do j = indices(k)+1, nibp-1 
                integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
            end do
            integral = integral + 0.5*((x(1) - x(nibp))*(y(1) + y(nibp)))
            do j = 1, indices(1)-1
                integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
            end do
            integral = integral + 0.5*((xintersect(1) - x(indices(1)))*(yintersect(1) + y(indices(1))))
            
            integral_target = integral_target + 0.5*((xtarget(indices_target(k)+1) - xintersect(k))*(yintersect(k) + Cptarget(indices_target(k)+1)))
            do j = indices_target(k)+1, nbp_target-1 
                integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(Cptarget(j+1) + Cptarget(j)))
            end do
            integral_target = integral_target + 0.5*((xtarget(1) - xtarget(nbp_target))*(Cptarget(1) + Cptarget(nbp_target)))
            do j = 1, indices_target(1)-1
                integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(Cptarget(j+1) + Cptarget(j)))
            end do
            integral_target = integral_target + 0.5*((xintersect(1) - xtarget(indices_target(1)))*(yintersect(1) + Cptarget(indices_target(1))))
        
        end if
        
        ! L1 Absolute Error norm
        L1_error = L1_error + abs(integral - integral_target)
        L1_error = -L1_error
        
    end subroutine getPressureDistribution
    
    
    subroutine getPressureDistributionPOD(L1_Error, tempNests)
        
        ! Variables
        implicit none
        integer :: i, j, k, nbp_target, maxbp, nibp
        logical :: ex
        double precision :: A, B, C, Atarget, Btarget, Ctarget, detAB, detCB, detAC, xtemp, ytemp, integral, integral_target, L1_Error, p0
        double precision, dimension(:), allocatable :: x, y, xtarget, Cptarget, xintersect, yintersect, xintersect2, yintersect2
        integer, dimension(:), allocatable :: ind, ind_target, indices, indices_target, order, marker
        double precision, dimension(maxDoF) :: tempNests
        double precision, dimension(:), allocatable :: newpressure
 
        allocate(newpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD " 
    
        ! Body of getPressureDistribution
        nibp = size(InnerBound)

        ! Read in Target Reference Data         
        inquire(file=DataFolder//'/Cp_target.txt', exist = ex)
        if (ex == .true.) then
            open(29,file=DataFolder//'/Cp_target.txt',form='formatted',status='old')
        else
            STOP "The file 'Cp_target' does not exist. Please generate!"
        end if

        read(29,*) nbp_target
        allocate(ind_target(nbp_target),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(xtarget(nbp_target),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(Cptarget(nbp_target),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        do i = 1, nbp_target
            read(29,*) xtarget(i), Cptarget(i)
        end do
        close(29)
       
        ! Allocate Arrays
        allocate(ind(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(x(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(y(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "

        ! Extract x coordinate
        allocate(RD%coord_temp(RD%np,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main "
        RD%coord_temp = RD%coord 
        call SubMovemesh(tempNests)
        x = RD%coord_temp(orderedBoundaryIndex,1)
        deallocate(RD%coord_temp)
        
        ! Calculate y(pressure coefficient)
        p0 = 1.344037930148188 ! Soft code later, this is specific for Ma = 0.729, Re = 6.5e7, T = 255.658 and l = 1 feet
        
        ! Extract Pressure using POD
        call InterpolateCoefficients(tempNests, newpressure)
        
        y = (newpressure(orderedBoundaryIndex)/p0 - 1)/(0.5*IV%gamma*IV%Ma**2)
          
        if (nibp > nbp_target) then
            maxbp = nibp
        else
            maxbp = nbp_target
        end if
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
        allocate(xintersect(maxbp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(yintersect(maxbp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        
        ! Identify intersections between target and current boundary
        k = 0
        do j = 1, nibp-1
            do i = 1, nbp_target-1
                
                ! Calculate intersection point
                A = y(j+1) - y(j)
                B = x(j) - x(j+1)
                C = A*x(j) + B*y(j)
                
                Atarget = Cptarget(i+1) - Cptarget(i)
                Btarget = xtarget(i) - xtarget(i+1)
                Ctarget = Atarget*xtarget(i) + Btarget*Cptarget(i)
                
                detAB = A*Btarget - Atarget*B
                detCB = C*Btarget - Ctarget*B
                detAC = A*Ctarget - Atarget*C
                if (abs(detAB) < abs(10e-12)) then
                    ! Parallel
                    xtemp = -1000
                    ytemp = -1000
                else
                    ! Non-parallel
                    xtemp = detCB/detAB
                    ytemp = detAC/detAB
                end if
                
                ! Check, if intersection point lays on the line segment considered 
                if ((xtemp - min(x(j),x(j+1))) > -10e-12 .and. (xtemp - max(x(j),x(j+1))) < 10e-12 .and. (ytemp - min(y(j),y(j+1))) > -10e-12  .and. (ytemp - max(y(j),y(j+1))) < 10e-12)  then
                    if ((xtemp - min(xtarget(i),xtarget(i+1))) > -10e-12 .and. (xtemp - max(xtarget(i),xtarget(i+1))) < 10e-12  .and. (ytemp - min(Cptarget(i),Cptarget(i+1))) > -10e-12  .and. (ytemp - max(Cptarget(i),Cptarget(i+1))) < 10e-12 )  then
                        k = k + 1
                        ind_target(k) = i
                        ind(k) = j
                        xintersect(k) = xtemp
                        yintersect(k) = ytemp
                    end if
                end if
            end do
        end do
        
        if (k == 0) then
            
            ! Integrate over both boundaries using trapezoidal rule
            integral = 0
            integral_target = 0
            L1_error = 0
            
            do j = 1, nibp-1 
                integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
            end do
            integral = integral + 0.5*((x(1) - x(nibp))*(y(1) + y(nibp)))
            do j = 1, nbp_target-1 
                integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(Cptarget(j+1) + Cptarget(j)))
            end do
            integral_target = integral_target + 0.5*((xtarget(1) - xtarget(nbp_target))*(Cptarget(1) + Cptarget(nbp_target)))
            
        else
           
            ! Identify doubled points
            allocate(marker(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Pressure Distribution "
            marker = 0
            j = k
            do i = 2, j
                if (xintersect(i) == xintersect(i-1) .and. yintersect(i) == yintersect(i-1)) then
                    marker(i) = 1
                    k = k - 1
                end if
            end do
        
            ! Array allocation with new array size j
            allocate(xintersect2(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "             
            allocate(yintersect2(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            allocate(indices(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            allocate(indices_target(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
        
            ! Exclude all points identified to be existing twice
            k = 0
            do i = 1, j           
                if (marker(i) == 0) then
                    k = k + 1
                    xintersect2(k) = xintersect(i)
                    yintersect2(k) = yintersect(i)
                    indices(k) = ind(i)
                    indices_target(k) = ind_target(i)
                end if
            end do
        
            ! Hand over of points to previous array name
            deallocate(xintersect)
            allocate(xintersect(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            xintersect = xintersect2
            deallocate(xintersect2)
            deallocate(yintersect)
            allocate(yintersect(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressureDistribution "
            yintersect = yintersect2
            deallocate(yintersect2)
            deallocate(ind)
            deallocate(ind_target)
        
            ! Integrate over both boundaries using trapezoidal rule
            integral = 0
            integral_target = 0
            L1_error = 0
            do i = 1, k-1
            
                integral = integral + 0.5*((x(indices(i)+1) - xintersect(i))*(yintersect(i) + y(indices(i)+1)))
                do j = indices(i)+1, indices(i+1)-1
                    integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
                end do
                integral = integral + 0.5*((xintersect(i+1) - x(indices(i+1)))*(yintersect(i+1) + y(indices(i+1))))
            
                integral_target = integral_target + 0.5*((xtarget(indices_target(i)+1) - xintersect(i))*(yintersect(i) + Cptarget(indices_target(i)+1)))
                do j = indices_target(i)+1, indices_target(i+1)-1 
                    integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(Cptarget(j+1) + Cptarget(j)))
                end do
                integral_target = integral_target + 0.5*((xintersect(i+1) - xtarget(indices_target(i+1)))*(yintersect(i+1) + Cptarget(indices_target(i+1))))
            
                ! L1 Absolute Error norm (root mean square)
                L1_error = L1_error + abs(integral - integral_target)
                integral = 0
                integral_target = 0
            
            end do
       
            integral = integral + 0.5*((x(indices(k)+1) - xintersect(k))*(yintersect(k) + y(indices(k)+1)))
            do j = indices(k)+1, nibp-1 
                integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
            end do
            integral = integral + 0.5*((x(1) - x(nibp))*(y(1) + y(nibp)))
            do j = 1, indices(1)-1
                integral = integral + 0.5*((x(j+1) - x(j))*(y(j+1) + y(j)))
            end do
            integral = integral + 0.5*((xintersect(1) - x(indices(1)))*(yintersect(1) + y(indices(1))))
            
            integral_target = integral_target + 0.5*((xtarget(indices_target(k)+1) - xintersect(k))*(yintersect(k) + Cptarget(indices_target(k)+1)))
            do j = indices_target(k)+1, nbp_target-1 
                integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(Cptarget(j+1) + Cptarget(j)))
            end do
            integral_target = integral_target + 0.5*((xtarget(1) - xtarget(nbp_target))*(Cptarget(1) + Cptarget(nbp_target)))
            do j = 1, indices_target(1)-1
                integral_target = integral_target + 0.5*((xtarget(j+1) - xtarget(j))*(Cptarget(j+1) + Cptarget(j)))
            end do
            integral_target = integral_target + 0.5*((xintersect(1) - xtarget(indices_target(1)))*(yintersect(1) + Cptarget(indices_target(1))))
        
        end if
        
        ! L1 Absolute Error norm
        L1_error = L1_error + abs(integral - integral_target)
        L1_error = -L1_error
        
    end subroutine getPressureDistributionPod
  
    
    subroutine AdaptiveSampling(newSnapshots, Fcompare)
    
        ! Variables
        implicit none
        logical :: Converge
        integer :: NoConv, k
        double precision :: Ftemp
        double precision, dimension(:), allocatable :: Fcompare
        double precision, dimension(:,:), allocatable :: tempSnapshots, newSnapshots
        integer, dimension(:), allocatable :: ConvA
    
        allocate(ConvA(2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation " 
        allocate(tempSnapshots((IV%NoSnap + 2*(OV%Gen-2)),maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        
        ! Body of AdaptiveSampling
        print *, 'Adaptive Sampling - Start Part 2 / 2'
                
        ! Store Snapshots in temporary Container       
        tempSnapshots(1:IV%NoSnap,:) = CS%Snapshots
        deallocate(CS%Snapshots)
                
        ! Check if Jobs for new Snapshots are ready
        call Sleep(IV%NoSnap + 2)
                
        ! Convergence Check
        NoConv = 0
        ConvA = 0
        do k = 1, 2
            Converge = .true.
            call FileCheckConvergence(Converge, (IV%NoSnap + k))
            if (Converge == .true.) then ! If Converged
                        
                ! Include new Snapshot
                NoConv = NoConv + 1
                tempSnapshots((IV%NoSnap + NoConv),:) = newSnapshots(k,:)
                ConvA(k) = 1
                        
            else
                ! Diverging Snapshot
            end if
        end do
        print *, 'NoConv:', NoConv
                
        allocate(character(len=3) :: istr)
        write(istr, '(1f3.1)') IV%Ma
        open(19,file=newdir//'/Fitness'//istr//'.txt',form='formatted',status='old',position='append')
        if (ConvA(1) == 1) then ! If converged check fitness
            call getObjectiveFunction(.false., Ftemp, NoSnapshot=(IV%NoSnap + 1))
            write(19,'(1I3, 1f17.10)',advance="no") 0, Ftemp
            print *, 'Real Fitness best: ', Ftemp
            do k = 1, IV%NoNests
                if (Fcompare(1) == OV%Fi(k)) then  ! Replace POD fitness with real fitness if still the same value                       
                    print *, 'Comparison POD/Real Fitness best: ', OV%Fi(k), '/', Ftemp
                    OV%Fi(k) = Ftemp
                    EXIT
                end if
            end do
        else ! Not converged set fitness to worst fitness to exclude solution
            OV%Fi(1) = OV%Fi(IV%NoNests)
        end if
 
        if (ConvA(2) == 1) then ! If converged check fitness
            call getObjectiveFunction(.false., Ftemp, NoSnapshot=(IV%NoSnap + 2))
            write(19,'(1I3, 1f17.10)',advance="no") 1, Ftemp
            print *, 'Comparison POD/Real Fitness worst: ', Fcompare(2), '/', Ftemp 
        end if
        write(19,*) ' '
        close(19)
        deallocate(OV%pressure)
        deallocate(OV%MaLocal)
        deallocate(OV%pTamb)
        deallocate(istr)
                               
        ! Resize Snapshots Array to include new Snapshots
        IV%NoSnap = IV%NoSnap + NoConv
        allocate(CS%Snapshots(IV%NoSnap,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "              
        CS%Snapshots = tempSnapshots(1:IV%NoSnap,:)
        deallocate(tempSnapshots)
                
        ! Re-Do POD including new Snapshots
        deallocate(OV%modes)
        deallocate(OV%coeff)
        call AllocateModesCoeff()
        call POD()
        print *, 'Adaptive Sampling - Finish Part 2 / 2'
    
    end subroutine AdaptiveSampling
    
    subroutine getRosenbrock(Fi, x)
    
        ! Variables
        implicit none
        integer :: i, d
        double precision :: Fi
        double precision, dimension(maxDoF) :: x
    
        ! Body of getRosenbrock
        d = 10
        Fi = 0
        do i = 1, d-1
            Fi = Fi + ((1 - x(i))**2 + 100*(x(i+1) -  x(i)**2)**2)
        end do
        Fi = -Fi
        OV%Precovery = abs(sqrt(sum(x(1:d)**2)) - sqrt(real(d)))
    
    end subroutine getRosenbrock
    
    subroutine getDeJong(Fi, x)
    
        ! Variables
        implicit none
        integer :: i, d
        double precision :: Fi
        double precision, dimension(maxDoF) :: x
    
        ! Body of getDeJong
        d = 50
        Fi = 0
        do i = 1, d
            Fi = Fi + x(i)**2
        end do
        Fi = -Fi
        OV%Precovery = abs(sqrt(sum(x(1:d)**2)))
    
    end subroutine getDeJong
    
    subroutine getAckley(Fi, x)
    
        ! Variables
        implicit none
        integer :: i, d
        double precision :: Fi, a1, a2
        double precision, dimension(maxDoF) :: x
        double precision, parameter :: pi = 3.14159265359
        double precision, parameter :: e = 2.718281828459045
    
        ! Body of getAckley
        d = 50
        Fi = 0
        a1 = 0.0
        a2 = 0.0
        do i = 1, d
            a1 = a1 + x(i)**2
            a2 = a2 + cos(2*pi*x(i))
        end do
        Fi = -20*exp(-0.2*sqrt(a1/d)) - exp(a2/d) + (20 + e)
        Fi = -Fi
        OV%Precovery = abs(sqrt(sum(x(1:d)**2)))
    
    end subroutine getAckley
        
end module Optimization