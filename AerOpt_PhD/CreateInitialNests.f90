module CreateInitialNests
    
        real rn                ! Counts random numbers
        integer*4 timeArray(3)    ! Holds the hour, minute, and second
        integer, dimension(:), allocatable :: put, get

        !In order to get a different sequence each time, we initialize the
        !seed of the random number function with the sum of the current
        !hour, minute, and second.
        
    contains
    
    subroutine CreateRandomNumber
 
        call RANDOM_SEED
        call RANDOM_NUMBER (rn)
        print *, rn
        
    end subroutine CreateRandomNumber 
    
end module CreateInitialNests