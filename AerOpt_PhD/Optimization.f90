module Optimization
    
    use CreateInitialNests
    use Toolbox
    use InputData
    use ReadData
    use GenerateInitialMeshes
    use CFD
    double precision, dimension(:,:), allocatable :: modes, coeff       ! Modes and Coefficient derived by the POD method
    real, dimension(:), allocatable :: engInNodes                       ! Engine inlet Nodes
    double precision, dimension(:,:), allocatable :: pressure           ! Pressure of All initial Snapshots
    real, dimension(:), allocatable :: Fi                               ! Vector with all current Fitness values
    real, dimension(:), allocatable :: NestOpt                          ! Control Point Coordinates of Optimum Geometry
        
contains
    
    subroutine SubOptimization()
    
        ! Variables
        implicit none
        real :: Ac, Ftemp, Fopt
        integer :: NoSteps, i, j, k, NoCPdim, NoTop, NoDiscard, l, randomNest
        real, dimension(:), allocatable :: NormFact, tempNests_Move, dist, tempNests, Fi_initial
        real, dimension(:,:), allocatable :: InitialNests_Move, newNests_Move, newNests
        integer, dimension(:), allocatable :: ind_Fi, ind_Fi_initial
        logical :: oob
        character(len=5) :: strNoSnap

        allocate(NormFact(av*IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(tempNests_Move(av*IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(dist(av*IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(InitialNests_Move(IV%NoSnap,av*IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(newNests_Move(IV%NoNests,av*IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(newNests(IV%NoNests,IV%NoDim*IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(tempNests(IV%NoDim*IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(NestOpt(IV%NoDim*IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        
        ! Body of SubOptimization            
        NoCPdim = av*IV%NoCP
        NoDiscard = nint(IV%Top2Low*IV%NoNests)
        NoTop = IV%NoNests - NoDiscard
        IV%Aconst = (sqrt(real(NoCPDim))/IV%NoLeviSteps)*IV%Aconst
        tempNests = (/ (0, i=1,(IV%NoDim*IV%NoCP)) /)
        
        ! Extract moving initial Nests
        j = 1
        do i = 1, size(cond)            
            if (cond(i) == -1) then
                InitialNests_Move(:,j) = InitialNests(:,i)
                j = j + 1
            end if
        end do
            
        print *, 'Generation        1'
        print *, 'Get POD modes and coeff for initial nests'
        
        ! Normalize InitialNests between 0 and 1
        do i = 1, NoCPdim        
            NormFact(i) = MxDisp_Move(i,1) - MxDisp_Move(i,2)
            InitialNests_Move(:,i) = InitialNests_Move(:,i)/NormFact(i) + 0.5
        end do
            
        ! Allocate modes and coeff size based on the user input of Number of POD Modes desired. If < 0, all Modes are considered.
        if (IV%NoPOMod < 0 .OR. IV%NoPOMod > IV%NoSnap) then
            allocate(modes(RD%np, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            allocate(coeff(IV%NoSnap, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        IV%NoPoMod = IV%NoSnap
        else
            allocate(modes(RD%np, IV%NoPOMod),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
            allocate(coeff(IV%NoSnap, IV%NoPOMod),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        end if
        
        call timestamp()
        call POD()
        call timestamp()
        ! Output: Modes and Coefficients of POD       
                
        call getengineInlet() ! Get boundary nodes, that define the Engine Inlet Plane
        ! Output: Engine Inlet Nodes(engInNodes)
        
        allocate(ind_Fi_initial(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(ind_Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Fi(IV%NoNests),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        allocate(Fi_initial(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Optimisation "
        
        ! Determine Distortion
        call getDistortion(Fi_initial)
        deallocate(pressure)
        ! Output: Distortion of all Snapshots as the Fitness (Fi)
        
        ! Pass on Top Snapshot parameters to Nests 
        ind_Fi_initial = (/ (i, i=1,IV%NoSnap) /)
        call QSort(Fi_initial,size(Fi), 'y', ind_Fi_initial) ! Result in Ascending order
        Fi = Fi_initial(1:IV%NoNests)
        newNests_Move = InitialNests_Move(ind_Fi_initial(1:IV%NoNests),:)
        newNests = InitialNests(ind_Fi_initial(1:IV%NoNests),:)
        deallocate(Fi_initial)
        deallocate(ind_Fi_initial)
        
        ! Write Output File for Analysis including Initial and all moved Nests of each Generation
        write(strNoSnap, *) IV%NoSnap
        open(29,file='Output_Data/newNests'//trim(strNoSnap)//'.txt')
        write(29, *) 'Initial Snapshots'
        write(29,'(<IV%NoSnap>f13.10)') InitialNests
        
        ! Loop over all Cuckoo Generations - each Generation creates new Nests
        do i = 2, IV%NoG
            print *, 'Generation ', i
            call timestamp()
            
            ! Re-order Fitness in ascending order
            ind_Fi = (/ (i, i=1,IV%NoNests) /)
            call QSort(Fi,size(Fi), 'y', ind_Fi)
            
            !! Change to Descending Order in Case of Maximization Problem
            !do j = 1, nint(IV%NoNests/2.0)
            !    Fi(j) = Fi(IV%NoNests-j+1)
            !    ind_Fi(j) = ind_Fi(IV%NoNests-j+1)
            !end do
            
            ! Re-order newNests for next Generation
            newNests_Move = newNests_Move(ind_Fi,:)
            newNests = newNests(ind_Fi,:)
            
            ! Store Optimum Fitness value
            print *, 'Current best solution:' , Fi(1)
            open(19,file='Output_Data/Fitness1000.txt')
            write(19,'(1f13.10)') Fi(1)                    
                    
            ! Store moved Nests in Output Analysis File
            write(29, *) 'Generation', (i-1)
            write(29,'(1000f13.10)') newNests
            
            !!*** Loop over Discarded Nests ***!!
            print *, 'Modify Discarded Cuckoos'
            do j = 1, NoDiscard
                
                ! Perform Random Walk using Levy Flight with a Cauchy Distribution
                Ac = IV%Aconst/(i**(1.0/2.0))
                call random_number(rn)
                NoSteps = nint(log(rn)*(-IV%NoLeviSteps))
                NoSteps = minval((/ NoSteps, IV%NoLeviSteps /))
                tempNests_Move = Ac*LevyWalk(NoSteps, NoCPdim) + newNests_Move(IV%NoNests-j+1,:)
                
                ! Check if out of bounds
                oob = 0
                do k = 1, NoCPdim
                    if (tempNests_Move(k) <= 1 .and. tempNests_Move(k) >= 0) then
                    else
                         oob = 1
                    end if               
                end do
                
                ! Re-evaluate Fitness of moved Nests
                if (oob == 1 .and. IV%constrain == .true.) then                  
                    ! Do nothing if out of bounds
                else
                    
                    ! Refill tempNests
                    l = 1
                    do k = 1, size(cond)            
                        if (cond(k) == -1) then
                            tempNests(k) = (tempNests_Move(l) - 0.5)*NormFact(l)
                            l = l + 1
                        end if
                    end do
                    
                    ! Update Values (Fitness, Nest Locations)
                    call ReEvaluateDistortion(NoCPDim, tempNests, Fi(IV%NoNests-j+1))
                    ! Output: ONE Fitnessvalue(Fi)
                    newNests_Move(IV%NoNests-j+1,:) = tempNests_Move
                    newNests(IV%NoNests-j+1,:) = tempNests
                    
                end if
                         
            end do
            
            !!*** Loop over Top Nests ***!!
            print *, 'Modify Top Cuckoos'
            do j = 1, NoTop
                
                ! Pick one of the Top Nests
                call random_number(rn)
                randomNest = nint((1 + (NoTop - 1)*rn))
                if (randomNest == j) then  ! Same Nest
                
                    ! Perform Random Walk instead                   
                    Ac = IV%Aconst/(i**2.0)
                    call random_number(rn)
                    NoSteps = nint(log(rn)*(-IV%NoLeviSteps))
                    NoSteps = minval((/ NoSteps, IV%NoLeviSteps /))
                    tempNests_Move = Ac*LevyWalk(NoSteps, NoCPdim) + newNests_Move(j,:)
                    
                else    ! Different Nest
                    
                    if (Fi(j) < Fi(randomNest)) then
                        
                        ! Cross-bread Nests in Direction of Nest j by Golden Ratio
                        dist = newNests_Move(randomNest,:) - newNests_Move(j,:)   ! Calculate Distance between Nests
                        dist = dist/(0.5*(1+sqrt(5.0)))                           ! Apply Golden Ratio
                        tempNests_Move = newNests_Move(j,:) + dist              ! Move Less Fit Nest
                        
                    elseif (Fi(randomNest) < Fi(j)) then
                        
                        ! Cross-bread in Direction of randomNest by Golden Ratio
                        dist = newNests_Move(j,:) - newNests_Move(randomNest,:)   ! Calculate Distance between Nests
                        dist = dist/(0.5*(1+sqrt(5.0)))                           ! Apply Golden Ratio
                        tempNests_Move = newNests_Move(randomNest,:) + dist     ! Move Less Fit Nest
                        
                    else
                        
                        ! Fitness is the same: Cross-bread Half Way
                        dist = newNests_Move(randomNest,:) - newNests_Move(j,:)   ! Calculate Distance between Nests
                        dist = dist*0.5                                         ! Apply Golden Ratio
                        tempNests_Move = newNests_Move(j,:) + dist              ! Move Less Fit Nest
                
                    end if
                   
                end if
                    
                ! Check if out of bounds
                oob = 0
                do k = 1, NoCPdim
                    if (tempNests_Move(k) <= 1 .and. tempNests_Move(k) >= 0) then
                    else
                         oob = 1
                    end if               
                end do
                
                ! Re-evaluate Fitness of moved Nests
                if (oob == 1 .and. IV%constrain == .true.) then                    
                    ! Do nothing if out of bounds
                else
                    
                    ! Refill tempNests
                    l = 1
                    do k = 1, size(cond)            
                        if (cond(k) == -1) then
                            tempNests(k) = (tempNests_Move(l) - 0.5)*NormFact(l)
                            l = l + 1
                        end if
                    end do
                    
                    ! Update Values (Fitness, Nest Locations)
                    call ReEvaluateDistortion(NoCPDim, tempNests, Ftemp)
                    ! Output: ONE Fitnessvalue(Fi)
                    
                    ! Check if new Fitness is better than a Random Top Nest, If yes replace values
                    call random_number(rn)
                    randomNest = nint((1 + (NoTop - 1)*rn))
                    if (Ftemp < Fi(randomNest)) then
                        newNests_Move(IV%NoNests-j+1,:) = tempNests_Move
                        newNests(IV%NoNests-j+1,:) = tempNests
                        Fi(randomNest) = Ftemp
                    end if
            
                end if
                
            end do
                     
        end do
        print *, 'Finished Cuckoo Search'   
        
        ! Find Optimum Cuckoo and write in File
        call QSort(Fi,size(Fi, dim = 1), 'y', ind_Fi)
        Fopt = Fi(1)
        write(19,'(1f13.10)') Fi(1)
        close(19)
        
        !write final Nests in File
        write(29, *) 'Generation', (i-1)
        write(29,'(1000f13.10)') newNests
        close(29)
        
        NestOpt = newNests(ind_Fi(1),:)
        
        print *, 'Optimum Fitness found:', Fopt
        print *, 'Optimum Geometry:'
        print *, NestOpt
        
!!! Calculate Error??
        
    end subroutine SubOptimization
        
    subroutine POD()
        
        ! Variables
        implicit none
        real :: Vamb, rho_amb, Rspec, gamma
        real, dimension(:), allocatable :: Output
        real, dimension(:), allocatable :: Vx, Vy, rho, e
        double precision, dimension(:,:), allocatable :: pressure2, var1, var2, modestemp
        character(len=:), allocatable :: istr
        double precision, dimension(:,:), allocatable :: V, var3, ones      
        
        ! Body of POD
        allocate(var3(1, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        allocate(ones(RD%np,1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
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
        allocate(pressure(RD%np, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        
        ! Precalculat ambient Parameters
        Vamb = IV%Ma*sqrt(IV%gamma*IV%R*IV%Tamb)          ! ambient velocity
        rho_amb = (IV%Pamb)/(IV%Tamb*IV%R)             ! ambient Density
            
        !Extract pressure of Snapshot Output file
        do i = 1, IV%NoSnap
            
            ! Determine correct String number
            call DetermineStrLen(istr, i)
            
            !if (IV%SystemType == 'W') then
            !     call TransferSolutionOutput()
            !     call communicateWin2Lin(trim(IV%Username), trim(IV%Password), 'FileCreateDir.scr', 'psftp')
            !     call system('move '//trim(IV%filename)//istr//'.resp '//OutFolder//trim(IV%filename)//istr//'.resp')
            !end if
            !
            if (IV%SystemType == 'W') then
                open(11, file=OutFolder//'/'//trim(IV%filename)//istr//'.resp')
            else
                open(11, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.resp')
            end if
            
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
            pressure(:,i) = rho_amb*(IV%gamma - 1)*rho*(Vamb**2)*(e - 0.5*(Vx**2 + Vy**2))
            ! Old Bernoulli Equation to calculate non-dimensional pressure:  pressure(:,i) = e + (1.0/2.0)*(IV%Ma**2)*rho*(Vx*Vx + Vy*Vy) 
            
            close(11)
            deallocate(istr)
            print *,'Pressure Snapshot', i 
                
        end do
        deallocate(Output)
        deallocate(Vx)
        deallocate(Vy)
        deallocate(e)
        deallocate(rho)

        ! Perform Single Value Decomposition
        allocate(pressure2(RD%np, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        pressure2 = pressure
        call SVD(pressure2, size(pressure, Dim = 1), size(pressure, Dim = 2),V)
        deallocate(pressure2)
            
        ! Calculate PO modes and initial Coefficients 
        allocate(var1(RD%np, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        allocate(modestemp(RD%np, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        allocate(var2(RD%np, IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in POD "
        
        var1 = matmul(pressure,V)
        ones(:,1) = (/ (1, i=1,RD%np) /)
        var2 = var1*var1
        do j = 1, IV%NoSnap
            var3(1,j) = sqrt(sum(var2(:,j), dim = 1))
        end do
        var2 = matmul(ones,var3)
        modestemp = var1 / var2
        modes = modestemp
!!!!! modes = modestemp(:,1:IV%NoPOMod)
        coeff = matmul(transpose(modes),pressure)
        
        deallocate(var1)
        deallocate(var2)
        deallocate(modestemp)
            
        ! Output: Modes and Coefficients of POD       
        !open(23,file='Output_Data/Coefficients.txt')
        !write(23,'(30f12.7)') coeff            
        !close(23)
        print *, 'All Modes and Coefficients Calculated'
                      
    end subroutine POD
        
    subroutine getDistortion(Distortion)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration
        
        ! Variables
        implicit none
        integer :: NoEngIN
        real, dimension(:), allocatable :: PlaneX, PlaneY, dPress, h, Area_trap, Pmid_x, Pmid_y, Area, Press_mid
        real, dimension(IV%NoSnap) :: Distortion
        real :: Press_ave, L
 
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
        
        ! Calculate coordinates of midpoints and afterwards the Area between them
        do i = 1, IV%NoSnap
                 
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
            Press_mid = pressure(engInNodes(2:(NoEngIN-1)),i) ! Extract Pressure of middle engine Inlet Nodes
            Press_ave = sum(Press_mid*Area, dim = 1)/sum(Area, dim = 1)
            
            ! Calculate Pressure Deviation
            ! print *, pressure(engInNodes,i)
            do j = 1, (NoEngIN)
                dPress(j) = abs(pressure(engInNodes(j),i) - Press_ave)            
            end do
                
            ! Determine Length and Height of Intercepting Plane
            L = 0
            do j = 1, (NoEngIN-1)
                h(j) = DistP2P(2, PlaneX(j), PlaneX(j+1), PlaneY(j), PlaneY(j+1))
                L = L + h(j)
            end do
                
            ! Apply Trapezoidal Rule to numerically integrate the Distortion
            Area_trap = h*(dPress(1:(NoEngIN-1)) + dPress(2:NoEngIN))/2.0
            Distortion(i) = sum(Area_trap, dim = 1)/(Press_ave*L)
            ! Output: Distortion
                
        end do
  
    end subroutine getDistortion
        
    subroutine getengineInlet()
    ! Output: Identify all nodes positioned at the engine Inlet (engInNodes)
        
        ! Variables
        implicit none
        integer, dimension(:,:), allocatable :: nodesall
        real, dimension(:), allocatable :: nodesvec
        real, dimension(2) :: point
        
        allocate(nodesall(RD%nbf,2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getEngineInlet "
        
        ! Body of getengineInlet
        j = 0
        do i = 1, RD%nbf
            point(1) = (RD%coord(RD%boundf(i,1),1)*15 + RD%coord(RD%boundf(i,2),1)*15)/2.0
            point(2) = (RD%coord(RD%boundf(i,1),2)*15 + RD%coord(RD%boundf(i,2),2)*15)/2.0
            if (point(1) == 0 .and. point(2) < 1 .and. point(2) > 0) then
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
        
    subroutine SVD(A, M, N, V2)

        ! Parameters
        implicit none
        integer :: M, N
        integer :: LDA, LDU, LDVT
        integer          LWMAX
        parameter        ( LWMAX = 10000)

        ! Local Scalars
        integer          INFO, LWORK

        ! Local Arrays
        double precision, dimension(:, :), allocatable :: V2
        double precision, dimension(:,:), allocatable :: U, VT
        double precision, dimension (:), allocatable :: S
        double precision, intent(in) ::                             A( M, N )
        double precision ::                             WORK( LWMAX )
        integer, dimension(8*min(M,N)) :: IWORK

        !call PRINT_MATRIX( 'Initial Matrix A', M, N, A, LDA )      
 
        ! Executable Statements
        write(*,*)'DGESVD Program Results'

        ! Define Array Size
        allocate(V2(N,N),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD "
        LDA = M
        LDU = M
        LDVT = N
        if (M > 1000) then
            LDU = 1
            allocate(U(LDU,M),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD "
        else
            allocate(U(LDU,M),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD "
        end if
        allocate(VT(LDVT,N),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in SVD "
        allocate(S(N))            
            
        ! Query the optimal workspace.
        LWORK = -1
        !call DGESDD( 'O', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
        call DGESVD( 'N', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
        LWORK = min( LWMAX, int( WORK( 1 ) ) )

        ! Compute SVD.
        !call DGESDD( 'O', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
        call DGESVD( 'N', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )

        ! Check for convergence.
        if( INFO.GT.0 ) then
        write(*,*)'The algorithm computing SVD failed to converge.'
        stop
        end if
            
        !! Print singular values.
        !call PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
        !
        !! Print left singular vectors.
        !call PRINT_MATRIX( 'Left singular vectors (stored columnwise)', M, N, U, LDU )
        !
        !! Print right singular vectors.
        !call PRINT_MATRIX( 'Right singular vectors (stored rowwise)', N, N, VT, LDVT )
        !
        ! Print right singular vectors transposed.
        !call PRINT_MATRIX( 'Right singular vectors (stored rowwise), Transposed', N, N, transpose(VT), LDVT )
            
        !Output: V Matrix
        !open(23,file='Output_Data/VMatrix.txt')
        !write(23,'(30f12.7)') transpose(VT)            
        !close(23)

        V2 = transpose(VT)
            
    end subroutine SVD
    
    function LevyWalk(NoSteps, NoCPdim)
    
        ! Variables
        implicit none
        integer :: median, scale, l, m, NoSteps, NoCPdim
        real, parameter :: pi = 3.14159265359
        real, dimension(NoCPdim) :: LevyWalk
        real, dimension(NoSteps) :: y
    
        ! Body of LevyWalk - Each Dimension walks
                
        do l = 1, NoCPdim
      
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
    
    subroutine ReEvaluateDistortion(NoCPDim, tempNests, Distortion)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration
        
        ! Variables
        implicit none
        integer :: NoEngIN, NoCPDim
        real, dimension(:), allocatable :: PlaneX, PlaneY, dPress, h, Area_trap, Pmid_x, Pmid_y, Area, Press_mid
        real, dimension(IV%NoDim*IV%NoCP) :: tempNests
        real :: Distortion
        real :: Press_ave, L
        real, dimension(:), allocatable :: newpressure
 
        allocate(newpressure(RD%np),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        allocate(PlaneX(size(engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        allocate(PlaneY(size(engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        allocate(dPress(size(engInNodes)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        allocate(h(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        allocate(Area_trap(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        allocate(Pmid_x(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        allocate(Pmid_y(size(engInNodes)-1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        allocate(Press_mid(size(engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        allocate(Area(size(engInNodes)-2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReEvaluateDistortion "
        
        ! Body of getDistortion
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
        call PressInterp(NoCPDim, newpressure, tempNests)
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
        ! Output: Distortion

  
    end subroutine ReEvaluateDistortion
    
    subroutine PressInterp(NoCPDim, newpressure, tempNests)
    ! Objective: Interpolation of Coefficients with Radial Basis Functions, based on normalized Gaussian RBF (see Hardy theory)
    
        ! Variables
        implicit none
        double precision, dimension(:,:), allocatable :: B_ar
        real, dimension(:,:), allocatable :: InitNests, coeff_temp ! = InitialNests, but will be manipulated in this function --> New Name, so it will not effect the InitialNests Array
        real, dimension(:), allocatable :: Init_Nests_temp, Lambda, newCoeff, Ipiv
        integer, dimension(:), allocatable :: ind_IN
        real, dimension(IV%NoDim*IV%NoCP) :: tempNests
        real, dimension(RD%np) :: newpressure
        real :: a, b, d, a2
        integer :: NoCPDim, l, m, n, LWMAX, LWORK, Info
        parameter        ( LWMAX = 10000)
        double precision, dimension(:), allocatable ::  WORK
    
        allocate(Ipiv(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(Work(LWMAX),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(InitNests(IV%NoSnap,IV%NoDim*IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(Init_Nests_temp(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(Lambda(size(coeff, dim=2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(newCoeff(size(coeff, dim=2)),stat=allocateStatus)
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
        InitNests = InitialNests
        coeff_temp = coeff
        ind_IN = (/ (i, i=1,IV%NoSnap) /)
        Init_Nests_temp = InitNests(:, (1+IV%NoDim*IV%NoCP-NoCPDim))
        call QSort(Init_Nests_temp, size(InitNests, dim = 1), 'y', ind_IN)
        InitNests = InitNests(ind_IN,:)
        coeff_temp = coeff_temp(:,ind_IN)
        
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
        d = 0.25*(d**2)
        
        
        ! Compute Coefficients for Interpolation
        do m = 1, IV%NoSnap
            do n = 1, IV%NoSnap       
                a = sqrt(sum(((InitNests(m,:)-InitNests(n,:))**2), dim = 1)) ! Distance
                B_ar(m,n) = 1.0/d*exp(-(a**2)/d)           
            end do     
        end do
            
        !Get Inverse of Matrix via LAPACK Library         
        call dgetrf(IV%NoSnap, IV%NoSnap, B_ar, IV%NoNests, Ipiv, info)
        call dgetri(IV%NoSnap, B_ar, IV%NoSnap, Ipiv, WORK, LWORK, Info)
        LWORK = min( LWMAX, int( WORK( 1 ) ) )
        call dgetri(IV%NoSnap, B_ar, IV%NoSnap, Ipiv, WORK, LWORK, Info)
            
        ! For each 'pressure field' f(:,k) ie vector of coefficients corresponding to mode k    
        newCoeff = (/ (0, i=1,IV%NoSnap) /)
        do l = 1, size(coeff_temp, dim = 2)
            
            Lambda = matmul(B_ar,coeff_temp(l,:))
            
            do m = 1, IV%NoSnap
                a = sqrt(sum(((tempNests-InitNests(m,:))**2), dim = 1)) ! Distance
                newCoeff(l) = newCoeff(l) + Lambda(m)*(1.0/d*exp(-(a**2)/d))
            end do 
        end do
        newpressure = matmul(modes,newCoeff)
        
    end subroutine PressInterp
    
    recursive subroutine CheckforConvergence(Iter)
    
        ! Variables
        implicit none
        integer :: FileSize, LastLine, ii
        integer, save :: NoConv
        integer,intent(inout) :: Iter
        logical :: Converge
        real, dimension(8) :: Input
        character(len=200) :: strCommand
        integer, dimension(:), allocatable :: DivNestPos
        real, dimension(:), allocatable :: MidPoints
        
        allocate(DivNestPos(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
        allocate(MidPoints(IV%NoCP*IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in PressInterp "
    
        ! Body of CheckforConvergence
        print *, 'Iteration', (Iter + 1)
        Converge = .true.
        NoConv = 0
        do i = 1, IV%NoSnap
            
            ! Determine correct String number
            call DetermineStrLen(istr, i)
            
            ! Open .rsd file to check, if the last line contains 'Nan' solutions, which would mean convergence fail
            if (IV%SystemType == 'W') then
                open(1, file=OutFolder//'/'//trim(IV%filename)//istr//'.rsd', STATUS="OLD")
            else
                open(1, file=newdir//'/'//OutFolder//'/'//trim(IV%filename)//istr//'.rsd', STATUS="OLD")
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
            
            ! All diverged nests are pulled halfway to midpoint(no movement center)
            if (Converge == .false.) then
                
                    print *, 'File', i, 'failed to converge and will be resimulated'
                    NoConv = NoConv + 1
                    DivNestPos(NoConv) = i
                    MidPoints = MxDisp(:,1) - (MxDisp(:,1) - MxDisp(:,2))/2.0  ! Midpoint calculation
                    InitialNests(i,:) = InitialNests(i,:) - ((InitialNests(i,:) - MidPoints)/2.0)   ! Half way between current Nest and Midpoint
                    
            end if
            Converge = .true.
            deallocate(istr)
            
        end do
        
        if (NoConv /= 0) then
            
            !!! Re-Do Mesh of diverged Nest
            allocate(RD%coord_temp(RD%np,IV%nodim),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Main " 
            do ii = 1, NoConv
            
                print *, "Regenerating Mesh", DivNestPos(ii), "/", NoConv
                RD%coord_temp = RD%coord
                call SubGenerateInitialMeshes(InitialNests(DivNestPos(ii),:))
                ! Output: new coordinates - Mesh with moved boundaries based on Initial Nest
            
!!!!! implement mesh quality test
            
                ! Determine correct string      
                call DetermineStrLen(istr, DivNestPos(ii))
            
                ! write snapshot to file
                call InitSnapshots()
                deallocate (istr)
            end do
            deallocate(RD%coord_temp)
            
            !! PreProcessing
            do ii = 1, NoConv
                call PreProcessing(DivNestPos(ii))
            end do
            
            !! Solver
            do ii = 1, NoConv
                
                call Solver(DivNestPos(ii))
                
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
            
            call Sleep()
            
            Iter = Iter + 1
            
            if (Iter < 4) then
                call CheckforConvergence(Iter)
            end if
            
            if (NoConv /= 0) then
                STOP 'Convergence of CFD Simulations could not be achieved. Check initial Mesh and/or Movement constraints!'
            end if
        end if
        
    end subroutine CheckforConvergence
    
    
end module Optimization