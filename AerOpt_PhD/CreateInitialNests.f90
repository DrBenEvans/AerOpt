module CreateInitialNests
    
        use Toolbox
        use InputData
        real :: rn                                            ! Counts random numbers
        real, dimension(:,:), allocatable :: MxDisp           ! Matrix with min/max Displacements
        real, dimension(:,:), allocatable :: MxDisp_Move      ! Matrix with only moving min/max Displacements
        integer, dimension(:), allocatable :: cond            ! Identifies zero/non-zero values
        real, dimension(:,:), allocatable :: InitialNests     ! Matrix containing the initial Nests
        integer :: av                                         ! Allocater Variable
        
    contains
    
    subroutine SubCreateInitialNests()
    
        ! Variables
        implicit none
        integer, dimension(IV%NoCP) :: ones                      ! Vector with Ones                                    
        allocate(MxDisp((IV%NoCP*IV%NoDim),2))
        allocate(cond(IV%NoCP*IV%NoDim))
        allocate(InitialNests(IV%NoNests,IV%NoDim*IV%NoCP))

        
        !!****Body of SubCreateInitialNests****!!      
        ones = (/ (1, i=1,IV%NoCP) /)
        
        ! Initialize min/max Displacement Matrix
        ! NoDim automatically defines the size of the Matrix
        MxDisp(:,1) = (/IV%xmax*ones, IV%ymax*ones, IV%zmax*ones/)
        MxDisp(:,2) = (/IV%xmax*(-1)*ones, IV%ymax*(-1)*ones, IV%zmax*(-1)*ones/)
        
        !**Reduction of MxDisp to Nonzero Values - (only for LHS-routine to reduce processing time)**!
        cond = MxDisp(:,1) /= MxDisp(:,2)
        
        ! Derive size for reduced Matrix
        av = 0
        if (IV%xmax /= 0.00) then
            av = av + 1
        end if
        if (IV%ymax /= 0.00 .and. IV%NoDim > 1) then
            av = av + 1
        end if
        if (IV%zmax /= 0.00 .and. IV%NoDim == 3) then
            av = av + 1
        end if
        allocate(MxDisp_Move(av*IV%NoCP,2))
        
        ! Reduce Matrix     
        j = 1
        do i = 1, size(cond)            
            if (cond(i) == -1) then
                MxDisp_Move(j,:) = MxDisp(i,:)
                j = j + 1
            end if
        end do
        !**END Reduction of...**!
               

        ! Execute Latin Hypercube Sampling with movable min/max Displacements
        call LHS(InitialNests, IV%NoNests, IV%NoCP)     
        !Output: InitialNests - an initial Sampling via LHS
        InitialNests(1,:) = (/ (0, i=1,(IV%NoCP*IV%NoDim)) /)

    end subroutine SubCreateInitialNests
    
    subroutine LHS(Sampling, NoSampPoints, NoPerm)
    
        ! Variables
        implicit none
        integer :: NoSampPoints, NoPerm
        integer :: ms
        real, dimension(:,:), allocatable :: Sampling
        real, dimension(:,:), allocatable :: Sampling_1D
        integer, dimension(:), allocatable :: zeros
        
        ! Body of LHS
        !!*** Based on the degrees of freedom(xmax, ymax & zmax definition ****!!
        !!*** & Dimension restriction) the LHS is performed 1,2 or 3 times. ***!!
        !!**** The size of MxDisp_Move indicates the degrees of freedom ****!!
        allocate(Sampling(NoSampPoints, IV%NoDim*NoPerm),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in LHS "
        allocate(Sampling_1D(NoSampPoints, NoPerm),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in LHS "
        allocate(zeros(NoSampPoints),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in LHS "
        zeros = (/ (0, i=1,NoSampPoints) /)
        ms = size(MxDisp_Move,1)
        if ( ms == NoPerm) then
            
            print *, 'Case 1'
            call LHS_1D(MxDisp_Move, Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
                       
            ! Based on the Number of Dimensions, NonMoving Columns need to be filled with zeros
            j = 1
            do i = 1, size(cond)            
                if (cond(i) == -1) then
                    Sampling(:,i) = Sampling_1D(:,j)
                    j = j + 1
                else
                    Sampling(:,i) =  zeros ! Non Moving
                end if
            end do            
            
        elseif (ms == NoPerm*2) then
            
            ! All possible combinations of Number of Dimensions and Moving/Non-Moving restrictions are considered.
            print *, 'Case 2, Part 1'
            call LHS_1D(MxDisp_Move(1:NoPerm,:), Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            if (IV%xmax == 0) then
                Sampling(:,(NoPerm+1):ms) = Sampling_1D ! Write Data in y-Matrix section
                do i = 1, NoPerm
                Sampling(:,i) = zeros ! Non Moving
                end do
            else
                Sampling(:,1:NoPerm) = Sampling_1D ! Write Data in x-Matrix section
            end if
            
            print *, 'Case 2,Part 2'
            call LHS_1D(MxDisp_Move((NoPerm+1):ms,:), Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            if (IV%xmax == 0) then     ! x: 0, y & z full
                Sampling(:,(2*NoPerm+1):3*NoPerm) = Sampling_1D ! Write Data in z-Matrix section                         
            elseif (IV%ymax == 0) then ! y: 0, x & z full
                do i = 1, NoPerm
                    Sampling(:,(NoPerm+i)) = zeros ! Non Moving
                end do
                Sampling(:,(2*NoPerm+1):3*NoPerm) = Sampling_1D ! Write Data in z-Matrix section
            elseif (IV%NoDim == 2) then ! just 2D - x & y full
                Sampling(:,(NoPerm+1):ms) = Sampling_1D ! Write Data in y-Matrix section
            else                     ! z: 0, x & y full  
                Sampling(:,(NoPerm+1):ms) = Sampling_1D ! Write Data in y-Matrix section
                do i = 1, NoPerm
                    Sampling(:,(2*NoPerm+i)) = zeros ! Non Moving
                end do
            end if    
                    
        elseif (ms == NoPerm*3) then
        
            ! Only one combination possible --> Fixed
            print *, 'Case 3, Part 1'
            call LHS_1D(MxDisp_Move(1:NoPerm,:), Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            Sampling(:,1:NoPerm) = Sampling_1D
            print *, 'Case 3, Part 2'
            call LHS_1D(MxDisp_Move((NoPerm+1):(2*NoPerm),:), Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            Sampling(:,(NoPerm+1):(NoPerm*2)) = Sampling_1D
            print *, 'Case 3, Part 3'
            call LHS_1D(MxDisp_Move((2*NoPerm+1):ms,:), Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            Sampling(:,(NoPerm*2+1):ms) = Sampling_1D
            
        else
            print *, 'Impossible Number of Dimensions.'
            stop
        end if
        !Output: Sampling - an initial Sampling vis LHS
        
    end subroutine LHS
    
    subroutine LHS_1D(MD_Move, Sampling_1D, NoSampPoints, NoPerm)
    
        ! Variables
        implicit none
        integer :: NoSampPoints, NoPerm
        real, dimension(NoSampPoints,NoPerm), intent(out) :: Sampling_1D
        integer, dimension(NoSampPoints) :: rp
        real, dimension(NoPerm,2) :: MD_Move 
        double precision, dimension(NoSampPoints) :: linSamp, linSamp2
        double precision :: max, min, first, last, dlS, ds, dlL1, dlL2, dL, dbound
    
        ! Body of LHS_1D
        !!************** A complicated way of applying LHS for a 1D sampling. *************!!
        !!********* The idea is to maximize the minimum distance between points. **********!!
        !!****** However, that only plays a role in case of non-uniform distribution. *****!!
        !! Here, a uniform distribution (linSpacing) by picking the Midpoints is selected. !!        
        max = MD_Move(1,1)
        min = MD_Move(1,2)
        linSamp = linSpacing(max, min, NoSampPoints) ! linear splitting of Design Space/Movement Domain
        
        !! Find maximum minimum distance between Points/Nests to redefine limits
        dlS = linSamp(1) - linSamp(2)
        ds = dlS/100000
        do i = 1, 100000
                
            first = max - ds*i
            last = min + ds*i
            linSamp = linSpacing(first, last, NoSampPoints)
                
            dlL1 = max - linSamp(1)
            dlL2 = (linSamp(2) - linSamp(1))/2
            if ((dlL1 + dlL2) > 10**(-6)) then
                dL = max - linSamp(1)
                exit
            end if
                
        end do
        
        ! Calculate Spacing with new limits
        first = max - dL
        last = min + dL
        linSamp = linSpacing(first, last, NoSampPoints)
        
        ! Simplified version: Redefine Limits(max/min) via Midpoint Calculation        
        linSamp2 = linSpacing(max, min, (NoSampPoints + 1))
        dbound = abs((linSamp2(1) - linSamp2(2))/2)
        max = max - dbound
        min = min + dbound
        linSamp2 = linSpacing(max, min, NoSampPoints) ! Calculate Spacing with new limits
        
        ! Test simplified version vs complicated version
        k = 0
        do i = 1, NoSampPoints
            if (abs(linSamp(i) - linSamp2(i))/IV%xmax < 1D-6) then
                    k = k + 1                
            end if
        end do
        if (k == NoSampPoints) then
            print *, 'Simplified Version WORKS'
        else
            print *, 'Simplified Version FAILED'
            print *, k
        end if
        
        ! Execute Random Permutation of distributed Nests for each Control Point
        do i = 1, NoPerm               
            call randperm(NoSampPoints, rp) ! see Subroutine, Output are integers
            do j = 1, NoSampPoints               
                Sampling_1D(j,i) = linSamp(rp(j)) ! Randomly permuted integers applied as indices (rp)                               
            end do                
        end do        
    
    end subroutine LHS_1D
    
end module CreateInitialNests