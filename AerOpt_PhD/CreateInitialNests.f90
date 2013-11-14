module CreateInitialNests
    
        real :: rn                                            ! Counts random numbers
        real, dimension(:,:), allocatable :: MxDisp           ! Matrix with min/max Displacements
        real, dimension(:,:), allocatable :: MxDisp_NonMov    ! Matrix with min/max Displacements
        integer, dimension(:), allocatable :: cond            ! Identifies zero/non-zero values

    contains
    
    subroutine CreateRandomNumber
 
        call RANDOM_SEED    ! Automatically generates a random initial number based on time and date
        call RANDOM_NUMBER (rn)
        
    end subroutine CreateRandomNumber
    
    subroutine LHC(NoDim, NoCtrPts, xmax, ymax, zmax)
    
        ! Variables
        integer :: i, j
        integer, dimension(NoCtrPts) :: ones
        allocate(MxDisp((NoCtrPts*NoDim),2))
        allocate(cond(NoCtrPts*NoDim))
        
        !!****Body of LHC****!!      
        do i = 1, NoCtrPts           
            ones(i) = 1
        end do
        
        ! Initialize min/max Displacement Matrix
        ! NoDim automatically defines the size of the Matrix
        MxDisp(:,1) = (/xmax*ones, ymax*ones, zmax*ones/)
        MxDisp(:,2) = (/xmax*(-1)*ones, ymax*(-1)*ones, zmax*(-1)*ones/)
        
        write(*,*), 'Full Matrix'
        print *, MxDisp
        
        !! Reduction of MxDisp to Nonzero Values - (only for LHS-routine to reduce processing time)
        cond = MxDisp(:,1) /= MxDisp(:,2)
        print *, cond
        
        ! Derive size for reduced Matrix
        j = 0
        do i = 1, size(cond)            
            if (cond(i) == -1) then               
                j = j + 1
            end if
        end do
        allocate(MxDisp_NonMov(j,2))
        
        ! Reduce Matrix
        ! MxDisp_NonMov = MxDisp(cond,:)        
        j = 1
        do i = 1, size(cond)            
            if (cond(i) == -1) then
                MxDisp_NonMov(j,:) = MxDisp(i,:)
                j = j + 1
            end if
        end do
                
        write(*,*), 'Reduced Matrix'
        print *, MxDisp_NonMov        

        ! Remember to empty MxDisp_NonMov to safe storage
    end subroutine LHC
    
end module CreateInitialNests