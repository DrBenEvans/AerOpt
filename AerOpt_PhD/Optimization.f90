    module Optimization
    
        real, dimension(:,:), allocatable :: Mesh_opt       
        
    contains
    
        subroutine SubOptimization(NoNests, NoCP, NoDim, cond, InitialNests, av1, MxDisp_Move)
    
            ! Variables
        
            real, dimension(:,:), allocatable :: MxDisp_Move
            real, dimension(:,:), allocatable :: NormFact
            real, dimension(NoNests,NoDim*NoCP) :: InitialNests
            integer, dimension(NoCP*NoDim) :: cond
            real, dimension(:,:), allocatable :: Mesh_new 
            real, dimension(:,:), allocatable :: InitialNests_Move 
            integer, parameter :: NoPOMod = -1
            integer, parameter :: MaxSteps = 100
            real, parameter :: A = 0.01
            real, parameter :: p = 0.75    
            integer :: av1, NoCPdim
        
    
            ! Body of subOptimization
            allocate(MxDisp_Move(av1*NoCP,2))
            allocate(NormFact(av1*NoCP,1))
            allocate(InitialNests_Move(NoNests,av1*NoCP))
            
            NoCPdim = av1*NoCP
        
            j = 1
            do i = 1, size(cond)            
                if (cond(i) == -1) then
                    InitialNests_Move(:,j) = InitialNests(:,i)
                    j = j + 1
                end if
            end do
        
            ! Normalize InitialNests between 0 and 1
            do i = 1, NoCPdim        
               NormFact(i,1) = MxDisp_Move(i,1) - MxDisp_Move(i,2)
               InitialNests_Move(:,i) = InitialNests_Move(:,i)/NormFact(i,1) + 0.5
            end do
    
        end subroutine subOptimization

    end module Optimization