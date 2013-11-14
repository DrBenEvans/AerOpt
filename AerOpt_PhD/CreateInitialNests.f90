module CreateInitialNests
    
        real :: rn                                          ! Counts random numbers
        real, dimension(:,:), allocatable :: MxDisp           ! Matrix with min/max Displacements
        real, dimension(:,:), allocatable :: MxDisp_NonMov    ! Matrix with min/max Displacements

        !In order to get a different sequence each time, we initialize the
        !seed of the random number function with the sum of the current
        !hour, minute, and second.
        
    contains
    
    subroutine CreateRandomNumber
 
        call RANDOM_SEED
        call RANDOM_NUMBER (rn)
        
    end subroutine CreateRandomNumber
    
    subroutine LHC(NoDim, NoCtrPts, xmax, ymax, zmax)
    
        ! Variables
        integer, dimension(NoCtrPts) :: ones
        allocate(MxDisp((NoDim*2),NoCtrPts))
          
        ! Body of LHC       
        do i = 1, NoCtrPts           
            ones(i) = 1
        end do
! Look into Details of Row and Column Definition
        MxDisp(1,:) = xmax*ones
        MxDisp(2,:) = xmax*(-1)*ones
        MxDisp(3,:) = ymax*ones
        MxDisp(4,:) = ymax*(-1)*ones
        if (NoDim == 3) then  
            MxDisp(5,:) = zmax*ones
            MxDisp(6,:) = zmax*ones
        end if
        
    end subroutine LHC
    
end module CreateInitialNests