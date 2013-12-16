    module Optimization
    
        use CreateInitialNests
        use Toolbox
        real, dimension(:,:), allocatable :: Mesh_opt
        double precision, dimension(:,:), allocatable :: modes, coeff
        real, dimension(:), allocatable :: engInNodes
        double precision, dimension(:,:), allocatable :: pressure
        real, dimension(:), allocatable :: Fi
        
    contains
    
        subroutine SubOptimization(NoNests, NoCP, NoDim, cond, InitialNests, MxDisp_Move, np, xmax, hMa)
    
            ! Variables       
            real, dimension(av*NoCP,2) :: MxDisp_Move
            real, dimension(av*NoCP) :: NormFact
            real, dimension(NoNests,NoDim*NoCP) :: InitialNests
            integer, dimension(NoCP*NoDim) :: cond
            real, dimension(:,:), allocatable :: Mesh_new 
            real, dimension(NoNests,av*NoCP) :: InitialNests_Move 
            integer, parameter :: NoPOMod = -1
            integer, parameter :: NoLeviSteps = 100
            real :: A = 0.01
            real, parameter :: p = 0.75    
            integer :: NoCPdim
        
    
            ! Body of SubOptimization            
            NoCPdim = av*NoCP
            NoDiscard = aint(p*NoNests)
            NoTop = NoNests - NoDiscard
            A = (sqrt(real(NoDim))/NoLeviSteps)*A
        
            ! Extract moving initial Nests
            j = 1
            do i = 1, size(cond)            
                if (cond(i) == -1) then
                    InitialNests_Move(:,j) = InitialNests(:,i)
                    j = j + 1
                end if
            end do
            
            print *, 'Get PO modes and coeff for initial nests'
        
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
            
            call getengineInlet(NoDim) ! Get boundary nodes, that define the Engine Inlet Plane
            call getDistortion(NoDim, NoNests) ! Determine Distortion
            !Output: Distortion
            allocate(Fi(NoNests))
            Fi = Distortion
            
            print *, 'Generation 1'
    
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
            open(23,file='Output_Data\Coefficients.txt')
            write(23,'(30f12.7)') coeff            
            close(23)
            print *, 'All Modes and Coefficients Calculated'
                      
        end subroutine POD
        
        subroutine getDistortion(NoDim, NoNests)
        ! Objective: Determine the Distortion of each Snapshot by using the Area Weighted Average Pressure and Trapezoidal Numerical Integration
        
            ! Variables
            use GenerateInitialMeshes
            integer :: NoEngIN
            real, dimension(size(engInNodes)) :: PlaneX, PlaneY, dPress
            real, dimension(size(engInNodes)-1) :: h, Area_trap, Pmid_x, Pmid_y
            real, dimension(size(engInNodes)-2) :: Area, Press_mid
            real, dimension(NONests) :: Distortion
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
                    
                ! Midpoint and Area Calculation
                do j = 1, (NoEngIN - 2)                       
                    Pmid_x(j+1) = (PlaneX(j+1) + PlaneX(j+2))/2.0
                    Pmid_y(j+1) = (PlaneY(j+1) + PlaneY(j+2))/2.0
                    Area(j) = sqrt((Pmid_x(j) - Pmid_x(j+1))**2 + (Pmid_y(j) - Pmid_y(j+1))**2)
                end do
                    
                ! Assign Right boundary nodes as Last Midpoint to act as boundaries
                Pmid_x(NoEngIN-1) = PlaneX(NoEngIN)
                Pmid_y(NoEngIN-1) = PlaneY(NoEngIN)
                
                ! Area Weighted Average Pressure (calculated based on the Areas)
                print *, pressure(118,i), pressure(234,i), pressure(715,i)
                Press_mid = pressure(engInNodes(2:(NoEngIN-1)),i) ! Extract Pressure of middle engine Inlet Nodes
                Press_ave = sum(Press_mid*Area, dim = 1)/sum(Area, dim = 1)
            
                ! Calculate Pressure Deviation
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
            call QSort(nodesvec, size(nodesvec))   
            call Unique(nodesvec, size(nodesvec), engInNodes)
            ! Output: unique vector engInNodes
    
        end subroutine getengineInlet

    end module Optimization