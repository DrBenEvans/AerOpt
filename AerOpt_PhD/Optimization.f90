module Optimization
    
    use CreateInitialNests
    use Toolbox
    real, dimension(:,:), allocatable :: Mesh_opt                       ! The final optimum Mesh
    double precision, dimension(:,:), allocatable :: modes, coeff       ! Modes and Coefficient derived by the POD method
    real, dimension(:), allocatable :: engInNodes                       ! Engine inlet Nodes
    double precision, dimension(:,:), allocatable :: pressure           ! Pressure of All initial Snapshots
    real, dimension(:), allocatable :: Fi                               ! Vector with all current Fitness values
    real, dimension(:), allocatable :: NestOpt                                       ! Control Point Coordinates of Optimum Geometry
        
contains
    
    subroutine SubOptimization(NoNests, NoCP, NoDim, cond, InitialNests, MxDisp_Move, np, xmax, hMa, p, Aconst, NoPOMod, NoLeviSteps, NoG, constrain)
    
        ! Variables       
        real, dimension(av*NoCP,2) :: MxDisp_Move
        real, dimension(av*NoCP) :: NormFact
        real, dimension(NoNests,NoDim*NoCP) :: InitialNests, newNests
        integer, dimension(NoCP*NoDim) :: cond
        real, dimension(:,:), allocatable :: Mesh_new 
        real, dimension(NoNests,av*NoCP) :: InitialNests_Move, newNests_Move
        integer :: NoCPdim, NoTop, NoDiscard, l, randomNest
        integer, dimension(NoNests) :: ind_Fi
        real, dimension(av*NoCP) :: tempNests_Move, dist
        real, dimension(NoDim*NoCP) :: tempNests
        logical :: constrain, oob
        real :: Ftemp, Fopt
        allocate(NestOpt(NoDim*NoCP))
        
        ! Body of SubOptimization            
        NoCPdim = av*NoCP
        NoDiscard = nint(p*NoNests)
        NoTop = NoNests - NoDiscard
        Aconst = (sqrt(real(NoCPDim))/NoLeviSteps)*Aconst
        tempNests = (/ (0, i=1,(NoDim*NoCP)) /)
        
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
        if (NoPOMod < 0) then
            allocate(modes(np, NoNests))
            allocate(coeff(NoNests, NoNests))
        else
            allocate(modes(np, NoPOMod))
            allocate(coeff(NoNests, NoPOMod))
        end if
            
        call POD(NoNests, hMa, np)
        ! Output: Modes and Coefficients of POD
        
        allocate(Fi(NoNests))    
        call getengineInlet(NoDim) ! Get boundary nodes, that define the Engine Inlet Plane
        ! Output: Engine Inlet Nodes(engInNodes)
        
        call getDistortion(NoDim, NoNests, Fi) ! Determine Distortion
        ! Output: Distortion as the Fitness (Fi)
        
        ! Loop over all Cuckoo Generations - each Generation creates new Nests
        do i = 2, NoG
            print *, 'Generation ', i
                
            ind_Fi = (/ (i, i=1,NoNests) /)
            call QSort(Fi,size(Fi), 'y', ind_Fi) ! Result in Ascending order
            
            ! Change to Descending Order
            !do j = 1, nint(NoNests/2.0)
            !    Fi(j) = Fi(NoNests-j+1)
            !    ind_Fi(j) = ind_Fi(NoNests-j+1)
            !end do
            print *, 'Current best solution:' , Fi(1)
            
            newNests_Move = InitialNests_Move(ind_Fi,:)
            newNests = InitialNests(ind_Fi,:)
            !!!!! If necessary, record Initial Nests (newnests) --> Output Text file
            
            !!*** Loop over Discarded Nests ***!!
            print *, 'Modify Discarded Cuckoos'
            do j = 1, NoDiscard
                
                ! Perform Random Walk using Levy Flight with a Cauchy Distribution
                Ac = Aconst/(i**(1.0/2.0))
                call random_number(rn)
                NoSteps = nint(log(rn)*(-NoLeviSteps))
                NoSteps = minval((/ NoSteps, NoLeviSteps /))
                tempNests_Move = Ac*LevyWalk(NoSteps, NoCPdim) + newNests_Move(NoNests-j+1,:)
                
                ! Re-evaluate Fitness of moved Nests and check if out of bounds
                oob = 0
                do k = 1, NoCPdim
                    if (tempNests_Move(k) <= 1 .and. tempNests_Move(k) >= 0) then
                    else
                         oob = 1
                    end if               
                end do
                
                
                if (oob == 1 .and. constrain == 1) then                  
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
                    call ReEvaluateDistortion(NoDim, NoNests, NoCP, NoCPDim, np, tempNests, Fi(NoNests-j+1))
                    ! Output: ONE Fitnessvalue(Fi)
                    newNests_Move(NoNests-j+1,:) = tempNests_Move
                    newNests(NoNests-j+1,:) = tempNests
                    
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
                    Ac = Aconst/(i**2.0)
                    call random_number(rn)
                    NoSteps = nint(log(rn)*(-NoLeviSteps))
                    NoSteps = minval((/ NoSteps, NoLeviSteps /))
                    tempNests_Move = Ac*LevyWalk(NoSteps, NoCPdim) + newNests_Move(NoNests-j+1,:)
                    
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
                    
                ! Re-evaluate Fitness of moved Nests and check if out of bounds
                oob = 0
                do k = 1, NoCPdim
                    if (tempNests_Move(k) <= 1 .and. tempNests_Move(k) >= 0) then
                    else
                         oob = 1
                    end if               
                end do
                
                
                if (oob == 1 .and. constrain == 1) then                    
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
                    call ReEvaluateDistortion(NoDim, NoNests, NoCP, NoCPDim, np, tempNests, Ftemp)
                    ! Output: ONE Fitnessvalue(Fi)
                    
                    ! Check if new Fitness is better than a Random Top Nest, If yes replace values
                    call random_number(rn)
                    randomNest = nint((1 + (NoTop - 1)*rn))
                    if (Ftemp < Fi(randomNest)) then
                        newNests_Move(NoNests-j+1,:) = tempNests_Move
                        newNests(NoNests-j+1,:) = tempNests
                        Fi(randomNest) = Ftemp
                    end if
            
                end if
                
            end do
                     
        end do
        print *, 'Finished Cuckoo Search'
        
        ! Find Optimum Cuckoo
        call QSort(Fi,size(Fi, dim = 1), 'y', ind_Fi)
        Fopt = Fi(1)
        NestOpt = newNests(ind_Fi(1),:)
        
        print *, 'Optimum Fitness found:', Fopt
        print *, 'Optimum Geometry:'
        print *, NestOpt
        
        !!! Calculate Error??
        
    end subroutine SubOptimization
        
    subroutine POD(NoNests, hMa, np)
        
        ! Variables
        real, dimension(:,:), allocatable :: Output
        double precision, dimension(:,:), allocatable :: pressure2, var1, var2, modestemp
        character(len=1) :: istr1
        character(len=2) :: istr2
        double precision, dimension(NoNests,NoNests) :: V
        double precision, dimension(np,1) :: ones
        double precision, dimension(1, NoNests) :: var3
            
        
        ! Body of POD
        allocate(Output(np, 5))
        allocate(pressure(np, NoNests))
        allocate(pressure2(np, NoNests))
        allocate(var1(np, NoNests))
        allocate(modestemp(np, NoNests))
        allocate(var2(np, NoNests))
            
        !Extract pressure of Snapshot Output file
        do i = 1, NoNests
            
            if (i < 10) then
                write( istr1, '(I1)' )  i
                open(8, file='Cases/case'//istr1//'.txt')
            else
                write( istr2, '(I2)' )  i
                open(8, file='Cases/case'//istr2//'.txt')
            end if
                
            do k = 1, 5
                do j = 1, np
                    read(8, *) Output(j,k)  ! rho, Vx, Vy, Vz, e
                end do
            end do
                
            pressure(:,i) = Output(:,5) + (1.0/2.0)*(hMa**2)*Output(:,1)*(Output(:,2)*Output(:,2)+Output(:,3)*Output(:,3)) !! Bernoulli Equation to calculate non-dimensional pressure            
            !!!!! Implement new pressure Calculation
            close(8)
            print *,'Pressure Snapshot', i 
                
        end do         
            
        ! Perform Single Value Decomposition
        ! call SVDCMP(pressure,np,NoNests,np,NoNests, V)
        pressure2 = pressure
        call SVD(pressure2, size(pressure, Dim = 1), size(pressure, Dim = 2),V)
            
        ! Calculate PO modes and initial Coefficients            
        var1 = matmul(pressure,V)
        ones(:,1) = (/ (1, i=1,np) /)
        var2 = var1*var1
        do j = 1, NoNests
            var3(1,j) = sqrt(sum(var2(:,j), dim = 1))
        end do
        var2 = matmul(ones,var3)
        modestemp = var1 / var2
        modes = modestemp
        coeff = matmul(transpose(modes),pressure)
           
            
        ! Output: Modes and Coefficients of POD       
        open(23,file='G:/Coefficients.txt')
        write(23,'(30f12.7)') coeff            
        close(23)
        print *, 'All Modes and Coefficients Calculated'
                      
    end subroutine POD
        
    subroutine getDistortion(NoDim, NoNests, Distortion)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration
        
        ! Variables
        use ReadData
        integer :: NoEngIN
        real, dimension(size(engInNodes)) :: PlaneX, PlaneY, dPress
        real, dimension(size(engInNodes)-1) :: h, Area_trap, Pmid_x, Pmid_y
        real, dimension(size(engInNodes)-2) :: Area, Press_mid
        real, dimension(NoNests) :: Distortion
        real :: Press_ave, L
 
        
        ! Body of getDistortion
        NoEngIN = size(engInNodes)
            
        ! Output: engInNodes
        PlaneX = coord(engInNodes,1)
        PlaneY = coord(engInNodes,2)
        
        ! Calculate coordinates of midpoints and afterwards the Area between them
        do i = 1, NoNests
                 
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
            
            print *, 'Area:'
            print *, Area
            
            ! Area Weighted Average Pressure (calculated based on the Areas)
            Press_mid = pressure(engInNodes(2:(NoEngIN-1)),i) ! Extract Pressure of middle engine Inlet Nodes
            Press_ave = sum(Press_mid*Area, dim = 1)/sum(Area, dim = 1)
            
            ! Calculate Pressure Deviation
            print *, pressure(engInNodes,i)
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
        
    subroutine getengineInlet(NoDim)
    ! Output: Identify all nodes positioned at the engine Inlet (engInNodes)
        
        ! Variables
        use Toolbox
        use GenerateInitialMeshes

        integer, dimension(size(boundf, dim = 1),2) :: nodesall
        real, dimension(:), allocatable :: nodesvec
        real, dimension(2) :: point
    
        ! Body of getengineInlet
        j = 0
        do i = 1, size(boundf, dim = 1)
            point(1) = (coord(boundf(i,1),1)*15 + coord(boundf(i,2),1)*15)/2.0
            point(2) = (coord(boundf(i,1),2)*15 + coord(boundf(i,2),2)*15)/2.0
            if (point(1) == 0 .and. point(2) < 1 .and. point(2) > 0) then
                j = j + 1
                nodesall(j,:) = boundf(i,1:2)           
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
    
        ! f90 mkl_lapack95_lp64.lib mkl_intel_lp64.lib mkl_core.lib mkl_sequential.lib
   
        ! Parameters
        integer          LDA, LDU, LDVT
        integer          LWMAX
        parameter        ( LWMAX = 10000)

        ! Local Scalars
        integer          INFO, LWORK

        ! Local Arrays
        double precision, dimension(N, N) :: V2
        double precision, dimension(:,:), allocatable :: U, VT
        double precision, dimension (:), allocatable :: S
        double precision, intent(in) ::                             A( M, N )
        double precision ::                             WORK( LWMAX )

        !call PRINT_MATRIX( 'Initial Matrix A', M, N, A, LDA )      
 
        ! Executable Statements
        write(*,*)'DGESVD Program Results'

        ! Define Array Size
        LDA = M
        LDU = M
        LDVT = N
        if (M > 1000) then
        LDU = 1
        allocate(U(LDU,M))
        else
        allocate(U(LDU,M))
        end if
        allocate(VT(LDVT,N))
        allocate(S(N))            
            
        ! Query the optimal workspace.
        LWORK = -1
        call DGESVD( 'N', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
        LWORK = min( LWMAX, int( WORK( 1 ) ) )

        ! Compute SVD.
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
        open(23,file='Output_Data/VMatrixhigh.txt')
        write(23,'(30f12.7)') transpose(VT)            
        close(23)

        V2 = transpose(VT)
            
    end subroutine SVD
    
    function LevyWalk(NoSteps, NoCPdim)
    
        ! Variables
        integer :: median, scale
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
    
    subroutine ReEvaluateDistortion(NoDim, NoNests, NoCP, NoCPDim, np, tempNests, Distortion)
    ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration
        
        ! Variables
        use GenerateInitialMeshes
        integer :: NoEngIN
        real, dimension(size(engInNodes)) :: PlaneX, PlaneY, dPress
        real, dimension(size(engInNodes)-1) :: h, Area_trap, Pmid_x, Pmid_y
        real, dimension(size(engInNodes)-2) :: Area, Press_mid
        real, dimension(NoDim*NoCP) :: tempNests
        real :: Distortion
        real :: Press_ave, L
        real, dimension(np) :: newpressure
 
        
        ! Body of getDistortion
        NoEngIN = size(engInNodes)
            
        ! Output: engInNodes
        PlaneX = coord(engInNodes,1)
        PlaneY = coord(engInNodes,2)
            
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
        call PressInterp(NoNests, NoCP, NoDim, NoCPDim, np, newpressure, InitialNests, tempNests)
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
    
    subroutine PressInterp(NoNests, NoCP, NoDim, NoCPDim, np, newpressure, InitNests, tempNests)
    ! Objective: Interpolation of Coefficients with Radial Basis Functions, based on normalized Gaussian RBF (see Hardy theory)
    
        ! Variables
        real, dimension(NoNests,NoDim*NoCP) :: InitNests ! = InitialNests, but will be manipulated in this function --> New Name, so it will not effect the InitialNests Array
        real, dimension(NoNests) :: Init_Nests_temp, Lambda, newCoeff
        integer, dimension(NoNests) :: ind_IN
        real, dimension(size(coeff, dim=1), size(coeff,dim=2)) :: coeff_temp
        real, dimension(NoDim*NoCP) :: tempNests
        real :: a, b, d
        double precision, dimension(NoNests, NoNests) :: B_ar
        real, dimension(np) :: newpressure
    
        ! Body of PressInterp
        
        
        open(99, file='G:/InitialNests.txt')         
        write(99,10) transpose(InitNests)
10      format(14f12.7)        
        close(99)
        
        open(15, file='G:/tempNests.txt')         
        write(15,15) tempNests
15      format(1f12.7)        
        close(15)
        
        ! Sort InitNests Matrix
        coeff_temp = coeff
        ind_IN = (/ (i, i=1,NoNests) /)
        Init_Nests_temp = InitNests(:, (1+NoDim*NoCP-NoCPDim))
        call QSort(Init_Nests_temp, size(InitNests, dim = 1), 'y', ind_IN)
        InitNests = InitNests(ind_IN,:)
        coeff_temp = coeff_temp(:,ind_IN)
        
        ! Compute Shape Parameter
        a = 0
        b = 0
        do l = 1, (NoNests - 1)
            do m = (l+1), NoNests
                a = a + norm2(InitNests(l,:)-InitNests(m,:))        ! Sum of all distances between points
                b = b + 1                                           ! Amount of points considered  
            end do    
        end do
        d = a/b             ! Mean Distance between Points
        d = 0.25*(d**2)
        
        ! For each 'pressure field' f(:,k) ie vector of coefficients corresponding to mode k
        newCoeff = (/ (0, i=1,NoNests) /)
        do l = 1, size(coeff_temp, dim = 1)
            
            ! Compute Coefficients for Interpolation
            do m = 1, NoNests
                do n = 1, NoNests
                    a = norm2(InitNests(m,:)-InitNests(n,:))        ! Distance
                    B_ar(m,n) = 1.0/d*exp(-(a**2)/d)           
                end do     
            end do
            
            call Inverse(B_ar, B_ar, NoNests)
            Lambda = matmul(B_ar,coeff_temp(l,:))
            
            do m = 1, NoNests
                a = norm2(tempNests-InitNests(m,:))                 ! Distance
                newCoeff(l) = newCoeff(l) + Lambda(m)*(1.0/d*exp(-(a**2)/d))
            end do
        newpressure = matmul(modes,newCoeff)
        end do
    end subroutine PressInterp

end module Optimization