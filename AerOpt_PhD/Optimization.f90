module Optimization
    
    use CreateSnapshots
    use Toolbox
    use InputData
    use ReadData
    use GenerateMesh
    use CFD
    use FDGD
    
    double precision, dimension(:,:), allocatable :: modes, coeff       ! Modes and Coefficient derived by the POD method
    double precision, dimension(:), allocatable :: meanpressure         ! Mean Pressure excluded of Snapshots in POD
    double precision, dimension(:), allocatable :: engInNodes           ! Engine inlet Nodes
    double precision, dimension(:,:), allocatable :: pressure           ! Pressure of All initial Snapshots
    double precision, dimension(:), allocatable :: Fi                   ! Vector with all current Fitness values
    double precision, dimension(:), allocatable :: NestOpt              ! Control Point Coordinates of Optimum Geometry
    double precision, dimension(:), allocatable :: PolCoeff             ! Polynomial Coefficients for POD RBF interpolation
    double precision, dimension(:,:), allocatable :: Weights            ! Weights for POD RBF interpolation
    double precision, dimension(:,:), allocatable :: Nests_Move, Nests              ! Nest regenerated with each generation
    logical :: InitConv                                                 ! Initial Convergence Check - TRUE: yes
        
contains
    
    subroutine SubOptimization()
    
        ! Variables
        implicit none
        double precision :: Ac, Ftemp, Fopt, temp
        integer :: i,j, k, l, ii, iii, NoSteps, NoTop, NoDiscard, randomNest, NoConv, NoTest
        double precision, dimension(:), allocatable :: NormFact, tempNests_Move, dist, tempNests, Fi_initial, Fcompare
        double precision, dimension(:,:), allocatable :: Snapshots_Move, tempSnapshots, newSnapshots, NestsTest, TopNest, TopNest_Move, tempNestA
        integer, dimension(:), allocatable :: ind_Fi, ind_Fi_initial, ConvA
        logical :: Converge
        character(len=5) :: strNoSnap
        !double precision, dimension(20, maxDOF) :: Nesting
        !double precision, dimension(RD%np) :: newpressure

        allocate(NormFact(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(tempNests_Move(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(dist(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Snapshots_Move(IV%NoSnap,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Nests_Move(IV%NoNests,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Nests(IV%NoNests,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(tempNests(maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(NestOpt(maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "     
        ! Specific for Adaptive Sampling
        allocate(newSnapshots(2,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Fcompare(2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(ConvA(2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "      
        
        ! Body of SubOptimization
        print *, ''
        print *, '*************************************'
        print *, '******    Start Optmization   *******'
        print *, '**                                 **'
        print *, '**  written by Dr. DAVID NAUMANN   **'
        print *, '**  supervised by Dr. BEN EVANS    **'
        print *, '**  supported by Dr. SEAN WALTON   **'
        print *, '*************************************'
        print *, ''
        
        ! Initializing general parameters
        NoDiscard = nint(IV%Top2Low*IV%NoNests)
        NoTop = IV%NoNests - NoDiscard
        IV%Aconst = (sqrt(real(IV%DoF))/IV%NoLeviSteps)*IV%Aconst
        tempNests = (/ (0, i=1,(maxDoF)) /)
        alpha = 1.0
        beta = 0.0
        
        ! Specific Parameters required for Top Nest monitoring
        allocate(TopNest(NoTop,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(TopNest_Move(NoTop,IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation " 
        
        ! Extract moving initial Nests
        j = 1
        do i = 1, size(cond)            
            if (cond(i) == -1) then
                Snapshots_Move(:,j) = Snapshots(:,i)
                j = j + 1
            end if
        end do
                   
        ! Normalize Snapshots(move) between 0 and 1
        do i = 1, IV%DoF        
            NormFact(i) = MxDisp_Move(i,1) - MxDisp_Move(i,2)
            Snapshots_Move(:,i) = Snapshots_Move(:,i)/NormFact(i) + 0.5
        end do
            
        ! Extract Pressure of Snapshots
        call timestamp()
        if (IV%POD == .true.) then
            call AllocateModesCoeff()
            ! Output: allocated modes and coeff based on the Number of POD Modes desired. If < 0, all Modes are considered.
            call POD()
            ! Output: Modes and Coefficients of POD  
            call ComputeRBFWeights()
            ! Output: Weights for RBF interpolation
        else
            call ExtractPressure(1, IV%NoSnap)
        end if             
        call timestamp()         
        if (IV%ObjectiveFunction == 2) then
            call getengineInlet() ! Get boundary nodes, that define the Engine Inlet Plane
            ! Output: Engine Inlet Nodes(engInNodes)
        end if
        
        ! Approximate Pressure Field of RandNest input file        
        !open(11,file=newdir//'/RandNest.txt')
        !read(11,'(20f17.10)') Nesting
        !close(11)
        !do i = 1, 20
        !    call DetermineStrLen(istr, i)
        !    call InterpolateCoefficients(Nesting(i,:), newPressure)
        !    open(23,file=newdir//'/NewPressure'//istr//'.txt')
        !    write(23,'(1f25.10)') newpressure           
        !    close(23)
        !    deallocate(istr)
        !end do
        !pause

        allocate(ind_Fi_initial(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(ind_Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Fi_initial(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
         
        ! Determine Objective Function
        do ii = 1, IV%NoSnap 
            call getObjectiveFunction(.false., Ftemp, Snapshots(ii,:), NoSnapshot=ii, NoGen=0)
            Fi_initial(ii) = Ftemp
        end do
        deallocate(pressure)
        ! Output: Distortion of all Snapshots as the Fitness (Fi)
      
        ! Pass on Top Snapshot parameters to Nests
        ind_Fi_initial = (/ (i, i=1,IV%NoSnap) /)
        call QSort(Fi_initial,size(Fi_initial), 'y', ind_Fi_initial) ! Result in Ascending order
        ! Change to Descending Order for Maximization Problem
        do j = 1, int(IV%NoSnap/2.0)
            temp = Fi_initial(j)
            Fi_initial(j) = Fi_initial(IV%NoNests-j+1)
            Fi_initial(IV%NoNests-j+1) = temp
            temp = ind_Fi_initial(j)
            ind_Fi_initial(j) = ind_Fi_initial(IV%NoNests-j+1)
            ind_Fi_initial(IV%NoNests-j+1) = temp
        end do
        Fi = Fi_initial(1:IV%NoNests)
        Nests_Move = Snapshots_Move(ind_Fi_initial(1:IV%NoNests),:)
        Nests = Snapshots(ind_Fi_initial(1:IV%NoNests),:)
        deallocate(Fi_initial)
        deallocate(ind_Fi_initial)
        
        ! Write Output File for Analysis including Initial and all moved Nests of each Generation
        allocate(character(len=3) :: istr)
        write(istr, '(1f3.1)') IV%Ma
        open(29,file=newdir//'/Nests'//istr//'.txt', form='formatted',status='unknown')   
        write(29, *) 'Snapshots'
        write(29,'(<IV%NoSnap>f17.10)') Snapshots
        close(29)
        open(19,file=newdir//'/Fitness'//istr//'.txt', form='formatted',status='unknown')
        write(19,*) 'Fitness'
        close(19)
        deallocate(istr)
       
        !!*** Loop over all Cuckoo Generations - each Generation creates new Nests ***!!
        do iii = 1, (IV%NoG - 1)
            print *, ''
            print *, '************************'
            print *, 'Generation ', iii
            print *, '************************'
            print *, ''
            print *, IV%Ma, IV%NoCN
            call timestamp()
            
            ! Re-order Fitness in ascending order
            ind_Fi = (/ (i, i=1,IV%NoNests) /)
            call QSort(Fi,size(Fi), 'y', ind_Fi)
            ! Change to Descending Order for Maximization Problem
            do j = 1, int(IV%NoNests/2.0)
                temp = Fi(j)
                Fi(j) = Fi(IV%NoNests-j+1)
                Fi(IV%NoNests-j+1) = temp
                temp = ind_Fi(j)
                ind_Fi(j) = ind_Fi(IV%NoNests-j+1)
                ind_Fi(IV%NoNests-j+1) = temp
            end do
            
            ! Re-order Nests for next Generation
            Nests_Move = Nests_Move(ind_Fi,:)
            Nests = Nests(ind_Fi,:)
            
            ! Store Fitness values
            allocate(character(len=3) :: istr)
            write(istr, '(1f3.1)') IV%Ma
            print *, 'Current best solutions:' , Fi(1:6)
            open(19,file=newdir//'/Fitness'//istr//'.txt',form='formatted',status='old',position='append')
            write(19,'(<IV%NoNests>f17.10)') Fi
            close(19) 
                    
            ! Store moved Nests in Output Analysis File
            open(29,file=newdir//'/Nests'//istr//'.txt',form='formatted',status='old',position='append')
            write(29, *) 'Generation', iii
            write(29,'(<IV%NoNests>f17.10)') Nests
            close(29)
            deallocate(istr)
          
            !!****** Adaptive Sampling - Start New Jobs (first and last Fitness) *******!!
            if (iii > 1 .and. iii < (IV%NoG - 1) .and. IV%AdaptSamp == .true.) then
                
                print *, 'Adaptive Sampling - Start Part 1 / 2'
                ! Extract First and Last Nest
                newSnapshots(1,:) = Nests(1,:)
                newSnapshots(2,:) = Nests(IV%NoNests,:)
                Fcompare(1) = Fi(1)
                Fcompare(2) = Fi(IV%NoNests)
                ! Get High Fidelity Solution
                call SubCFD((IV%NoSnap + 1), (IV%NoSnap + 2), newSnapshots, 2)
                print *, 'Adaptive Sampling - Finish Part 1 / 2'
                
            end if
            
            !!*** Loop over Discarded Nests ***!!
            print *, ''
            print *, 'Modify Discarded Cuckoos for Generation', (iii+1)
            print *, ''
            do ii = IV%NoNests, (NoTop + 1), -1
                
                ! Perform Random Walk using Levy Flight with a Cauchy Distribution
                Ac = IV%Aconst/((iii+1)**(1.0/2.0))
                call random_number(rn)
                NoSteps = nint(log(rn)*(-IV%NoLeviSteps))
                NoSteps = minval((/ NoSteps, IV%NoLeviSteps /))
                tempNests_Move = Ac*LevyWalk(NoSteps, IV%DoF) + Nests_Move(ii,:)
               
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
                do k = 1, size(cond)            
                    if (cond(k) == -1) then
                        tempNests(k) = (tempNests_Move(l) - 0.5)*NormFact(l)
                        l = l + 1
                    end if
                end do

                if (IV%POD == .true.) then
                    ! calculate Objective Function via POD (Fitness, Nest Locations)
                    call getObjectiveFunction(.true., Fi(ii), tempNests)               
                    ! Output: ONE Fitnessvalue(Fi)
                end if
               
                ! Embed temporary Nest into Nests
                Nests_Move(ii,:) = tempNests_Move
                Nests(ii,:) = tempNests
                
            end do

            ! Full Fidelity Solutions passed on to Solver
            if (IV%POD == .false.) then        
                call SubCFD((iii*IV%NoNests + NoTop + 1), ((iii + 1)*IV%NoNests), Nests((NoTop + 1):IV%NoNests,:), NoDiscard)             
            end if            
            
            !!*** Loop over Top Nests ***!!
            print *, ''
            print *, 'Modify Top Cuckoos for Generation', (iii+1)
            print *, ''
            do ii = 1, NoTop
           
                ! Pick one of the Top Nests
                call random_number(rn)
                randomNest = nint((1 + (NoTop - 1)*rn))
                if (randomNest == ii) then  ! Same Nest
                
                    ! Perform Random Walk instead                   
                    Ac = IV%Aconst/((iii+1)**2.0)
                    call random_number(rn)
                    NoSteps = nint(log(rn)*(-IV%NoLeviSteps))
                    NoSteps = minval((/ NoSteps, IV%NoLeviSteps /))
                    tempNests_Move = Ac*LevyWalk(NoSteps, IV%DoF) + Nests_Move(ii,:)
                    
                else    ! Different Nest
                    
                    if (Fi(ii) > Fi(randomNest)) then
                        
                        ! Cross-bread Nests in Direction of Nest j by Golden Ratio
                        dist = Nests_Move(randomNest,:) - Nests_Move(ii,:)   ! Calculate Distance between Nests
                        dist = dist/(0.5*(1+sqrt(5.0)))                           ! Apply Golden Ratio
                        tempNests_Move = Nests_Move(ii,:) + dist              ! Move Less Fit Nest
                        
                    elseif (Fi(randomNest) > Fi(ii)) then
                        
                        ! Cross-bread in Direction of randomNest by Golden Ratio
                        dist = Nests_Move(ii,:) - Nests_Move(randomNest,:)   ! Calculate Distance between Nests
                        dist = dist/(0.5*(1+sqrt(5.0)))                           ! Apply Golden Ratio
                        tempNests_Move = Nests_Move(randomNest,:) + dist     ! Move Less Fit Nest
                        
                    else
                        
                        ! Fitness is the same: Cross-bread Half Way
                        dist = Nests_Move(randomNest,:) - Nests_Move(ii,:)   ! Calculate Distance between Nests
                        dist = dist*0.5                                         ! Apply Golden Ratio
                        tempNests_Move = Nests_Move(ii,:) + dist              ! Move Less Fit Nest
                
                    end if
                   
                end if
                    
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
                
                ! Refill tempNests and de-normalize
                l = 1
                do k = 1, size(cond)            
                    if (cond(k) == -1) then
                        tempNests(k) = (tempNests_Move(l) - 0.5)*NormFact(l)
                        l = l + 1
                    end if
                end do
                
                ! Store tempNests
                TopNest(ii,:) = tempNests
                TopNest_move(ii,:) = tempNests_move
                
                if (IV%POD == .true.) then
                    ! Update Objective Function via POD (Fitness, Nest Locations)
                    call getObjectiveFunction(.true., Ftemp, tempNests)
                    ! Output: ONE Fitnessvalue(Fi)
                    
                    ! Check if new Fitness is better than a Random Top Nest, If yes replace values
                    call random_number(rn)
                    randomNest = nint((1 + (NoTop - 1)*rn))
                    if (Ftemp > Fi(randomNest)) then
                        Nests_Move(randomNest,:) = tempNests_Move
                        Nests(randomNest,:) = tempNests
                        Fi(randomNest) = Ftemp 
                     end if
                end if
                
            end do
             
            ! Store moved Nests in Output Analysis File
            open(39,file=newdir//'/TopNest.txt',form='formatted',status='unknown',position='append')
            write(39, *) 'Generation', (iii+1)
            write(39,'(<NoTop>f17.10)') TopNest
            close(39)
          
            ! Store moved Nests in Output Analysis File
            open(29,file=newdir//'/Neststemp.txt',form='formatted',status='unknown')
            write(29,'(<IV%NoNests>f17.10)') Nests
            close(29)
            
            if (IV%POD == .false.) then                
                !Generate Full Fidelity Solution of new TopNest      
                call SubCFD((iii*IV%NoNests + 1), (iii*IV%NoNests + NoTop), TopNest, NoTop)
                ! Check ALL new Nests
                print *, 'Generation: ', iii
                call PostSolverCheck(((iii + 1)*IV%NoNests), 1, Nests_Move, Nests)   
 
                ! Evaluate Fitness of Full Fidelity Nest Solutions
                print *, 'Extract Pressure of Generation', iii
                call ExtractPressure((iii*IV%NoNests + 1), ((iii + 1)*IV%NoNests))
                do ii = 1, NoTop
                    call getObjectiveFunction(.false., Ftemp, NoSnapshot=ii, NoGen=iii)
              
                    ! Check if new Fitness is better than a Random Top Nest, If yes replace values
                    call random_number(rn)
                    randomNest = nint((1 + (NoTop - 1)*rn))
                    if (Ftemp > Fi(randomNest)) then
                        Nests_Move(randomNest,:) = TopNest_Move(ii,:)
                        Nests(randomNest,:) = TopNest(ii,:)
                        Fi(randomNest) = Ftemp
                    end if
                end do

                do ii = (NoTop + 1), IV%NoNests
                    call getObjectiveFunction(.false., Fi(ii), NoSnapshot=ii, NoGen=iii)
                end do
                deallocate(pressure)
            end if

            !!*** Adaptive Sampling - Finish and Integrate New Jobs (first and last Fitness) ***!!
            if (iii > 1 .and. iii < (IV%NoG - 1) .and. IV%AdaptSamp == .true.) then
                print *, 'Adaptive Sampling - Start Part 2 / 2'
                
                ! Store Snapshots in temporary Container
                allocate(tempSnapshots((IV%NoSnap + 2*IV%NoG),maxDoF),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
                tempSnapshots(1:IV%NoSnap,:) = Snapshots
                deallocate(Snapshots)
                
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
                        
                        ! Delete Errorfiles if Snapshot diverged      
                        call DetermineStrLen(istr, (IV%NoSnap + k))        
                        call DeleteErrorFiles(istr)
                        if (IV%SystemType == 'W' .and. IV%runOnCluster == 'Y')   then    ! AerOpt is executed from a Windows machine           
                            call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')
                        else
                            call system('chmod a+x ./FileCreateDir.scr')
                            call system('./FileCreateDir.scr')
                        end if                
                        deallocate(istr)
                  
                    end if
                end do
                print *, 'NoConv:', NoConv
                
                allocate(character(len=3) :: istr)
                write(istr, '(1f3.1)') IV%Ma
                open(19,file=newdir//'/Fitness'//istr//'.txt',form='formatted',status='old',position='append')
                if (ConvA(1) == 1) then ! If converged check fitness
                    call getObjectiveFunction(.false., Ftemp, NoSnapshot=(IV%NoSnap + 1), NoGen=0)
                    write(19,'(1I3, 1f17.10)',advance="no") 0, Ftemp
                    print *, 'Real Fitness best: ', Ftemp
                    do k = 1, NoTop
                        if (Fcompare(1) == Fi(k)) then  ! Replace POD fitness with real fitness if still the same value                       
                            print *, 'Comparison POD/Real Fitness best: ', Fi(k), '/', Ftemp
                            Fi(k) = Ftemp
                        end if
                    end do
                else ! Not converged set fitness to worst fitness to exclude solution
                    Fi(1) = Fi(IV%NoNests)
                end if
 
                if (ConvA(2) == 1) then ! If converged check fitness
                    call getObjectiveFunction(.false., Ftemp, NoSnapshot=(IV%NoSnap + 2), NoGen=0)
                    write(19,'(1I3, 1f17.10)',advance="no") 1, Ftemp
                    print *, 'Comparison POD/Real Fitness worst: ', Fcompare(2), '/', Ftemp 
                end if
                close(19)
                deallocate(pressure)
                deallocate(istr)
                               
                ! Resize Snapshots Array to include new Snapshots
                IV%NoSnap = IV%NoSnap + NoConv
                allocate(Snapshots(IV%NoSnap,maxDoF),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "              
                Snapshots = tempSnapshots(1:IV%NoSnap,:)
                deallocate(tempSnapshots)
                
                ! Re-Do POD including new Snapshots
                deallocate(modes)
                deallocate(coeff)
                call AllocateModesCoeff()
                call POD()
                print *, 'Adaptive Sampling - Finish Part 2 / 2'
                
            end if
            write(19,*) ''
        
        end do
        print *, 'Finished Cuckoo Search'   
        
        ! Write Results in File
        allocate(character(len=3) :: istr)
        write(istr, '(1f3.1)') IV%Ma
        open(19,file=newdir//'/Fitness'//istr//'.txt',form='formatted',status='old',position='append')
        call QSort(Fi,size(Fi, dim = 1), 'y', ind_Fi)

        ! Change to Descending Order for Maximization Problem
        do j = 1, int(IV%NoNests/2.0)
            temp = Fi(j)
            Fi(j) = Fi(IV%NoNests-j+1)
            Fi(IV%NoNests-j+1) = temp
            temp = ind_Fi(j)
            ind_Fi(j) = ind_Fi(IV%NoNests-j+1)
            ind_Fi(IV%NoNests-j+1) = temp
        end do

        Fopt = Fi(1)
        write(19,'(<IV%NoNests>f17.10)') Fi
        close(19)
        
        ! write final Nests in File
        open(29,file=newdir//'/Nests'//istr//'.txt',form='formatted',status='old',position='append')
        write(29, *) 'Generation', iii
        write(29,'(<IV%NoNests>f17.10)') Nests
        close(29)
        deallocate(istr)
        
        NestOpt = Nests(ind_Fi(1),:)
        print *, ''
        print *, '************************'
        print *, 'Last Generation'
        print *, '************************'
        print *, ''
        print *, 'Optimum Fitness found:', Fopt
        print *, ''
        print *, 'Optimum Geometry:'
        print *, ''
        print *, NestOpt
        
    end subroutine SubOptimization
    
    subroutine ExtractPressure(Start, Ending)
        
        ! Variables
        implicit none
        integer :: Start, Ending, Length, i, k, j
        double precision :: Vamb, rho_amb, Rspec, gamma
        double precision, dimension(:), allocatable :: Output, Vx, Vy, rho, e
        character(len=:), allocatable :: istr

        ! Body of Extract Pressure
        Length = Ending - Start + 1
        allocate(rho(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        allocate(e(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        allocate(Vx(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        allocate(Vy(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        allocate(Output(6*RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        allocate(pressure(RD%np, Length),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        
        ! Precalculate ambient Parameters
        print *, 'Start Pressure Extraction'
        Vamb = IV%Ma*sqrt(IV%gamma*IV%R*IV%Tamb)          ! ambient velocity
        rho_amb = (IV%Pamb)/(IV%Tamb*IV%R)             ! ambient Density
            
        !Extract pressure of Snapshot Output file      
        do i = Start, Ending
           
            ! Determine correct String number
            call DetermineStrLen(istr, i)
           
            open(11, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.resp', form='formatted',status='old')     
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
            pressure(:,(i - Start + 1)) = (IV%gamma - 1.0)*rho*((e - 0.5*(Vx**2 + Vy**2))) ! Non-dimensional  
            ! Old Bernoulli Equation to calculate non-dimensional pressure:  pressure(:,i) = e + (1.0/2.0)*(IV%Ma**2)*rho*(Vx*Vx + Vy*Vy) 
            
            close(11)
            deallocate(istr)
                
        end do
        print *,'Pressure of .resp files extracted'
        deallocate(Output)
        deallocate(Vx)
        deallocate(Vy)
        deallocate(e)
        deallocate(rho)

    end subroutine ExtractPressure
    
    subroutine POD()
        
        ! Variables
        implicit none
        integer :: i
        double precision, dimension(:,:), allocatable :: pressure2, modestemp, V

        ! Body of POD
        print *, 'Get POD modes and coefficients of Snapshots'
        allocate(meanpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        
        ! Extract pressure of Snapshot Output file           
        call ExtractPressure(1, IV%NoSnap)

        ! Exclude mean pressure from POD       
        if (IV%meanP == .true.) then
            do i = 1, IV%NoSnap
                meanpressure = meanpressure + pressure(:,i)
            end do
            meanpressure = (1/IV%NoSnap)*meanpressure
            do i = 1, IV%NoSnap
                pressure(:,i) = pressure(:,i) - meanpressure
            end do
        end if
        
        ! Perform Single Value Decomposition
        allocate(pressure2(RD%np, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        pressure2 = pressure
        call SVD(pressure2, size(pressure, Dim = 1), size(pressure, Dim = 2), modestemp)
        !deallocate(pressure2)
        print *, 'Finished SVD'
                
        modes = modestemp
!!!!! modes = modestemp(:,1:IV%NoPOMod)

        ! Matmul: coeff = modes'*pressure   coeff = matmul(transpose(modes),pressure)
        CALL DGEMM('T','N',size( modes, dim = 2), size( pressure, dim = 2), size( modes, dim = 1),alpha,modes,size( modes, dim = 1),pressure,size( pressure, dim = 1),beta,coeff,size( coeff, dim = 1))
        
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
        !open(23,file=newdir//'/Pressure_Real.txt')
        !write(23,'(1f25.10)') pressure(:,2)           
        !close(23)
        
        print *, 'All Modes and Coefficients Calculated'

    end subroutine POD
        
    subroutine getDistortion(Distortion, NoSnapshot)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration
        
        ! Variables
        implicit none
        integer :: NoEngIN, NoSnapshot, j
        double precision, dimension(:), allocatable :: PlaneX, PlaneY, dPress, h, Area_trap, Pmid_x, Pmid_y, Area, Press_mid
        double precision :: Press_ave, L, Distortion
 
        allocate(PlaneX(size(engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(PlaneY(size(engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(dPress(size(engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(h(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area_trap(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_x(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Pmid_y(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Press_mid(size(engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
        allocate(Area(size(engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortion "
     
        ! Body of getDistortion
        NoEngIN = size(engInNodes)
            
        ! Output: engInNodes
        PlaneX = RD%coord(engInNodes,1)
        PlaneY = RD%coord(engInNodes,2)
        
        !!*** Calculate coordinates of midpoints and afterwards the Area between them ***!!
                 
        ! Assign Right boundary nodes as First Midpoint to act as boundaries
        Pmid_x(1) = PlaneX(1)
        Pmid_y(1) = PlaneY(1)
        ! Assign Right boundary nodes as Last Midpoint to act as boundaries
        Pmid_x(NoEngIN-1) = PlaneX(NoEngIN)
        Pmid_y(NoEngIN-1) = PlaneY(NoEngIN)                    
        ! Midpoint and Area Calculation
        do j = 2, (NoEngIN - 2)                       
            Pmid_x(j) = (PlaneX(j) + PlaneX(j+1))/2.0
            Pmid_y(j) = (PlaneY(j) + PlaneY(j+1))/2.0
        end do      
            
        do j = 1, (NoEngIN - 2)
            Area(j) = sqrt((Pmid_x(j) - Pmid_x(j+1))**2 + (Pmid_y(j) - Pmid_y(j+1))**2)
        end do
        
        ! Area Weighted Average Pressure (calculated based on the Areas)
        Press_mid = pressure(engInNodes(2:(NoEngIN-1)),NoSnapshot) ! Extract Pressure of middle engine Inlet Nodes
        Press_ave = sum(Press_mid*Area, dim = 1)/sum(Area, dim = 1)  
        ! Calculate Pressure Deviation
        ! print *, pressure(engInNodes,i)
        do j = 1, (NoEngIN)
            dPress(j) = abs(pressure(engInNodes(j),NoSnapshot) - Press_ave)            
        end do
                 
        ! Determine Length and Height of Intercepting Plane
        L = 0
        do j = 1, (NoEngIN-1)
            h(j) = DistP2P(2, PlaneX(j), PlaneX(j+1), PlaneY(j), PlaneY(j+1))
            L = L + h(j)
        end do
              
        ! Apply Trapezoidal Rule to numerically integrate the Distortion
        Area_trap = h*(dPress(1:(NoEngIN-1)) + dPress(2:NoEngIN))/2.0
        Distortion = sum(Area_trap, dim = 1)/(Press_ave*L)
        Distortion = Distortion*(-1) ! To adapt to maximization Problem
        ! Output: Distortion
  
    end subroutine getDistortion
        
    subroutine getengineInlet()
    ! Output: Identify all nodes positioned at the engine Inlet (engInNodes)
        
        ! Variables
        implicit none
        integer :: i,j
        integer, dimension(:,:), allocatable :: nodesall
        double precision, dimension(:), allocatable :: nodesvec
        double precision, dimension(2) :: point
        
        allocate(nodesall(RD%nbf,2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getEngineInlet "
      
        ! Body of getengineInlet
        j = 0
        do i = 1, RD%nbf
            if (RD%boundf(i,3) == 8) then
                j = j + 1
                nodesall(j,:) = RD%boundf(i,1:2)           
            end if
        end do
        ! Outcome: A list of all nodes related to the engine Inlet, including possible doubling
    
        allocate(nodesvec(2*j))
        nodesvec = (/nodesall(1:j,1), nodesall(1:j,2)/)
        call QSort(nodesvec, size(nodesvec), 'n')   
        call Unique(nodesvec, size(nodesvec), engInNodes)
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
    
            call random_number(rn)
            do m = 1, NoSteps
                y(m) = median + scale*tan(pi*rn)
            end do
    
            LevyWalk(l) = sum(y, dim = 1);
        
        end do
    
    end function LevyWalk
    
    subroutine getDistortionPOD(tempNests, Distortion)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration
        
        ! Variables
        implicit none
        integer :: NoEngIN, k
        double precision, dimension(:), allocatable :: PlaneX, PlaneY, dPress, h, Area_trap, Pmid_x, Pmid_y, Area, Press_mid
        double precision, dimension(maxDoF) :: tempNests
        double precision :: Distortion
        double precision :: Press_ave, L
        double precision, dimension(:), allocatable :: newpressure
 
        allocate(newpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(PlaneX(size(engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(PlaneY(size(engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(dPress(size(engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(h(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(Area_trap(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(Pmid_x(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(Pmid_y(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(Press_mid(size(engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        allocate(Area(size(engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDistortionPOD "
        
        ! Body of getDistortionPOD
        NoEngIN = size(engInNodes)
            
        ! Output: engInNodes
        PlaneX = RD%coord(engInNodes,1)
        PlaneY = RD%coord(engInNodes,2)
            
        ! Calculate coordinates of midpoints and afterwards the Area between them
                 
        ! Assign Right boundary nodes as First Midpoint to act as boundaries
        Pmid_x(1) = PlaneX(1)
        Pmid_y(1) = PlaneY(1)
                    
        ! Midpoint and Area Calculation
        do k = 1, (NoEngIN - 2)                       
            Pmid_x(k+1) = (PlaneX(k+1) + PlaneX(k+2))/2.0
            Pmid_y(k+1) = (PlaneY(k+1) + PlaneY(k+2))/2.0
            Area(k) = sqrt((Pmid_x(k) - Pmid_x(k+1))**2 + (Pmid_y(k) - Pmid_y(k+1))**2)
        end do
                    
        ! Assign Right boundary nodes as Last Midpoint to act as boundaries
        Pmid_x(NoEngIN-1) = PlaneX(NoEngIN)
        Pmid_y(NoEngIN-1) = PlaneY(NoEngIN)
                
        ! Area Weighted Average Pressure (calculated based on the Areas)
        call InterpolateCoefficients(tempNests, newpressure)
        Press_mid = newpressure(engInNodes(2:(NoEngIN-1))) ! Extract Pressure of middle engine Inlet Nodes
        Press_ave = sum(Press_mid*Area, dim = 1)/sum(Area, dim = 1)
            
        ! Calculate Pressure Deviation
        do k = 1, (NoEngIN)
            dPress(k) = abs(newpressure(engInNodes(k)) - Press_ave)            
        end do
                
        ! Determine Length and Height of Intercepting Plane
        L = 0
        do k = 1, (NoEngIN-1)
            h(k) = DistP2P(2, PlaneX(k), PlaneX(k+1), PlaneY(k), PlaneY(k+1))
            L = L + h(k)
        end do
                
        ! Apply Trapezoidal Rule to numerically integrate the Distortion
        Area_trap = h*(dPress(1:(NoEngIN-1)) + dPress(2:NoEngIN))/2.0
        Distortion = sum(Area_trap, dim = 1)/(Press_ave*L)
        Distortion = Distortion*(-1) ! To adapt to maximization Problem
        ! Output: Distortion

  
    end subroutine getDistortionPOD
    
    subroutine ComputeRBFWeights()
    ! Objective: Compute Weights for the RBF interpolation scheme applied later in the POD reconstruction
    
        ! Variables
        implicit none
        double precision, dimension(:,:),allocatable :: RBFMatrix, A, b, coeff_temp
        double precision, dimension(:), allocatable :: PolBasis, x, NormFact
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
        allocate(PolCoeff(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        allocate(Weights(IV%NoSnap,IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "  
        ! Matrix Parameters
        allocate(RBFMatrix(IV%NoSnap,IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        allocate(PolBasis(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        if (IV%Pol == .true.) then
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
        allocate(NormFact(IV%DoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        
        ! Body of CoefficientInterpolation
        ! A*x = b with A = [ RBFMatrix PolBasis, Polbasis^T 0] x = Weights of RBF and polynomial  b = [f 0] with f being the original coefficients of one Snapshot
        ! Details, see Computation Approach Book
        
        ! Normalize Snapshots between 0 and 10
        !do i = 1, DoF        
        !    NormFact(i) = MxDisp_Move(i,1) - MxDisp_Move(i,2)
        !    Snapshots(:,i) = (Snapshots(:,i)/NormFact(i) + 0.5)
        !end do

        allocate(coeff_temp(size(coeff, dim=1), size(coeff,dim=2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        coeff_temp = coeff
        coeff_temp = transpose(coeff_temp) 
    
        ! Generate Matrix of Radial Basis Functions (multiquadratic)
        ShapeParameter = 1.0
        do m = 1, IV%NoSnap
            do l = 1, IV%NoSnap
                z = sqrt(sum(((Snapshots(m,:)-Snapshots(l,:))**2), dim = 1))
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
            Weights(l,:) = x(1:IV%NoSnap)          
            if (IV%Pol == .true.) then
                PolCoeff(l) = x(1+IV%NoSnap)
            end if
        
        end do
        
        !open(23,file=newdir//'/Weights.txt')
        !write(23,'(10f55.10)') Weights(1:10,:)           
        !close(23)
        !open(23,file=newdir//'/Polynomial.txt')
        !write(23,'(1f25.10)') PolCoeff           
        !close(23)
        
        ! Denormalize Snapshots
        !do i = 1, DoF        
        !    NormFact(i) = MxDisp_Move(i,1) - MxDisp_Move(i,2)
        !    Snapshots(:,i) = (Snapshots(:,i) - 0.5)*NormFact(i)
        !end do
        
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
    
    subroutine InterpolateCoefficients(newNest, newpressure)
    
        ! Variables
        implicit none
        double precision, dimension(RD%np) :: newpressure
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
            x = sqrt(sum(((newNest-Snapshots(l,:))**2), dim = 1))
            RBFVector(l) = Multiquadratic(x, ShapeParameter)
        end do
    
        do l = 1, IV%NoSnap
            RBFterm = 0
            do m = 1, IV%NoSnap
                RBFterm = RBFterm + Weights(l,m)*RBFVector(m)
            end do

            if (IV%Pol == .true.) then
                newCoeff(l,1) = PolCoeff(l) + RBFterm
            else
                newCoeff(l,1) = RBFterm
            end if
        end do
        
        ! Matmul: var1 = modes*newCoeff   newpressure = matmul(modes,newCoeff(:,1))
        CALL DGEMM('N','N',size( modes, dim = 1), size( newCoeff, dim = 2), size( modes, dim = 2),alpha,modes,size( modes, dim = 1),newCoeff,size( newCoeff, dim = 1),beta,newpressure,size( newpressure, dim = 1))
       
        if (IV%meanP == .true.) then
            newpressure = newpressure + meanpressure
        end if
        
        open(23,file=newdir//'/NewPressure.txt')
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
        allocate(Lambda(size(coeff, dim = 1),1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(CoeffVec(size(coeff, dim = 1),1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(newCoeff(size(coeff, dim=2),1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(ind_IN(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(coeff_temp(size(coeff, dim=1), size(coeff,dim=2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(B_ar(IV%NoSnap, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "

        
        ! Body of PressInterp
        ! Initial parameters for LAPACK
        LWORK = -1
        Ipiv = 0.0
        Info = 0
        
        ! Sort InitNests Matrix
        InitNests = Snapshots
        coeff_temp = coeff
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
        CALL DGEMM('N','N',size( modes, dim = 1), size( newCoeff, dim = 2), size( modes, dim = 2),alpha,modes,size( modes, dim = 1),newCoeff,size( newCoeff, dim = 1),beta,newpressure,size( newpressure, dim = 1))
                
        open(23,file=newdir//'/NewPressureOld.txt')
        write(23,'(1f25.10)') newpressure           
        close(23)
        
    end subroutine PressInterp
    
    subroutine AllocateModesCoeff()
    
        ! Variables
        implicit none
    
        ! Body of AllocateModesCoeff
        if (IV%NoPOMod < 0 .OR. IV%NoPOMod > IV%NoSnap) then
            allocate(modes(RD%np, IV%NoSnap),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            allocate(coeff(IV%NoSnap, IV%NoSnap),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        else
            allocate(modes(RD%np, IV%NoPOMod),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            allocate(coeff(IV%NoSnap, IV%NoPOMod),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        end if
    
    end subroutine AllocateModesCoeff
    
    subroutine getObjectiveFunction(PODevaluation, Fi, tempNests, NoSnapshot, NoGen)
    
        ! Variables
        implicit none
        double precision, dimension(maxDoF), optional :: tempNests
        integer, optional :: NoSnapshot, NoGen
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
                call getzeroLiftPOD(tempNests, Fi)
            end if
        else
            if (IV%ObjectiveFunction == 2) then
                call getDistortion(Fi, NoSnapshot)
            elseif (IV%ObjectiveFunction == 1) then
                call getLiftandDrag(Fi, (NoGen*IV%NoNests + NoSnapshot))
            elseif (IV%ObjectiveFunction == 3) then
                call getzeroLift(Fi, (NoGen*IV%NoNests + NoSnapshot))
            end if
        end if
    
    end subroutine getObjectiveFunction
    
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
            vec1 = RD%coord_temp(RD%boundf(i,1),:)
            vec2 = RD%coord_temp(RD%boundf(i,2),:)
            tangent = vec2 - vec1
            dx = DistP2P(IV%NoDim, vec1(1), vec2(1), vec1(2), vec2(2))
            intexit = 0
            do j = 1, NoIB
                if (RD%boundf(i,1) == Innerbound(j)) then
                    tangentArray(j,:) = tangentArray(j,:) + tangent
                    lengthArray(j) = lengthArray(j) + dx
                    intexit = intexit + 1
                    if (intexit == 2) then
                        EXIT
                    end if
                end if
                if (RD%boundf(i,2) == Innerbound(j)) then
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
        if (IV%SystemType == 'W') then
            LastLine = FileSize/107
        else     
            LastLine = FileSize/106
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
            vec1 = RD%coord_temp(RD%boundf(i,1),:)
            vec2 = RD%coord_temp(RD%boundf(i,2),:)
            tangent = vec2 - vec1
            dx = DistP2P(IV%NoDim, vec1(1), vec2(1), vec1(2), vec2(2))
            intexit = 0
            do j = 1, NoIB
                if (RD%boundf(i,1) == Innerbound(j)) then
                    tangentArray(j,:) = tangentArray(j,:) + tangent
                    lengthArray(j) = lengthArray(j) + dx
                    intexit = intexit + 1
                    if (intexit == 2) then
                        EXIT
                    end if
                end if
                if (RD%boundf(i,2) == Innerbound(j)) then
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
        if (IV%SystemType == 'W') then
            LastLine = FileSize/107
        else     
            LastLine = FileSize/106
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
        
end module Optimization