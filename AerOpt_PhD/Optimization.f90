    module Optimization
    
        use CreateInitialNests
        use Toolbox
        real, dimension(:,:), allocatable :: Mesh_opt
        double precision, dimension(:,:), allocatable :: modes, coeff
        
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
        
            j = 1
            do i = 1, size(cond)            
                if (cond(i) == -1) then
                    InitialNests_Move(:,j) = InitialNests(:,i)
                    j = j + 1
                end if
            end do
            
            print *, 'Generation 1'
            print *, 'Get PO modes and coeff for initial nests'
        
            ! Normalize InitialNests between 0 and 1
            do i = 1, NoCPdim        
               NormFact(i) = MxDisp_Move(i,1) - MxDisp_Move(i,2)
               InitialNests_Move(:,i) = InitialNests_Move(:,i)/NormFact(i) + 0.5
            end do
            
            if (NoPOMod < 0) then
                allocate(modes(np, NoNests))
                allocate(coeff(NoNests, NoNests))
            else
                allocate(modes(np, NoPOMod))
                allocate(coeff(NoNests, NoPOMod))
            end if
            
            call POD(NoNests, hMa, np)
            ! Output: Modes and Coefficients of POD
            
            call getDistortion
            !Output: Distortion
    
        end subroutine SubOptimization
        
        subroutine POD(NoNests, hMa, np)
        
            ! Variables
            real, dimension(:,:), allocatable :: Output
            double precision, dimension(:,:), allocatable :: pressure, pressure2, var1, var2, modestemp
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
            !call SVDCMP(pressure,np,NoNests,np,NoNests, V)
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
                      
        end subroutine POD
        
        subroutine getDistortion()
        
            ! Variables
!1. bflag
!2. getDistortion 
        
            ! Body of getDistortion
            
        
        end subroutine getDistortion

    end module Optimization