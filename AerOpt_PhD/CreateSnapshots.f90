module CreateSnapshots
    
        use Toolbox
        use InputData
        double precision :: rn                                            ! Counts random numbers
        double precision, dimension(:,:), allocatable :: MxDisp           ! Matrix with min/max Displacements
        double precision, dimension(:,:), allocatable :: MxDisp_Move      ! Matrix with only moving min/max Displacements
        integer, dimension(:), allocatable :: cond            ! Identifies zero/non-zero values
        double precision, dimension(:,:), allocatable :: Snapshots     ! Matrix containing the initial Nests
        
    contains
    
    subroutine SubCreateSnapshots()
    ! Objective: Create Nests (Locations of  for Snapshots
    
        ! Variables
        implicit none
        integer :: i, j
        integer, dimension(IV%NoCN) :: ones                      ! Vector with Ones                                    
        allocate(MxDisp(maxDoF,2))
        allocate(cond(maxDoF))
        allocate(Snapshots(IV%NoSnap,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CreateSnapshots "
        allocate(MxDisp_Move(IV%DoF,2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CreateSnapshots "
        
        
        !!****Body of SubCreateSnapshots****!!      
        ones = (/ (1, i=1,IV%NoCN) /)
        
        ! Initialize min/max Motion Matrix    
        MxDisp(:,1) = (/IV%xrange(1:IV%NoCN), IV%yrange(1:IV%NoCN), IV%zrange(1:IV%NoCN), IV%angle(1:IV%NoCN)/) ! min
        MxDisp(:,2) = (/IV%xrange((IV%NoCN+1):2*IV%NoCN), IV%yrange((IV%NoCN+1):2*IV%NoCN), IV%zrange((IV%NoCN+1):2*IV%NoCN), IV%angle((IV%NoCN+1):2*IV%NoCN)/) ! max
         
        !**Reduction of MxDisp to Nonzero Values - (only for LHS-routine to reduce processing time)**!
        cond = MxDisp(:,1) /= MxDisp(:,2)         
        j = 1
        do i = 1, size(cond)            
            if (cond(i) == -1) then
                MxDisp_Move(j,:) = MxDisp(i,:)
                j = j + 1
            end if
        end do               

        ! Execute Latin Hypercube Sampling with movable min/max Displacements
        call LHS(Snapshots, IV%NoSnap, IV%NoCN)     
        !Output: Snapshots - an initial Sampling via LHS      
        Snapshots(1,:) = 0.0

    end subroutine SubCreateSnapshots
    
    subroutine LHS(Sampling, NoSampPoints, NoPerm)
    
        ! Variables
        implicit none
        integer :: NoSampPoints, NoPerm
        double precision, dimension(NoSampPoints, maxDoF) :: Sampling
        double precision, dimension(:,:), allocatable :: Sampling_1D
        
        ! Body of LHS
        !!*** Based on the degrees of freedom(xmax, ymax & zmax definition ****!!
        !!*** & Dimension restriction) the LHS is performed 1,2 or 3 times. ***!!
        !!**** The size of MxDisp_Move indicates the degrees of freedom ****!!
        Sampling = 0.0
        allocate(Sampling_1D(NoSampPoints, NoPerm),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in LHS "
            
        print *, 'x Dimension'
        call LHS_1D(MxDisp(1:NoPerm,:), Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
        Sampling(:,1:NoPerm) = Sampling_1D ! Write Data in x-Matrix section
            
        print *, 'y Dimension'
        call LHS_1D(MxDisp((NoPerm+1):2*NoPerm,:), Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
        Sampling(:,(NoPerm+1):2*NoPerm) = Sampling_1D ! Write Data in y-Matrix section
            
        print *, 'z Dimension'
        call LHS_1D(MxDisp((2*NoPerm+1):3*NoPerm,:), Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
        Sampling(:,(2*NoPerm+1):3*NoPerm) = Sampling_1D ! Write Data in z-Matrix section  

        print *, 'angular Dimension'
        call LHS_1D(MxDisp((3*NoPerm+1):4*NoPerm,:), Sampling_1D, NoSampPoints, NoPerm) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
        Sampling(:,(3*NoPerm+1):4*NoPerm) = Sampling_1D ! Write Data in z-Matrix section 
        !Output: Sampling - an initial Sampling vis LHS
        
    end subroutine LHS
    
    subroutine LHS_1D(MD_Move, Sampling_1D, NoSampPoints, NoPerm)
    
        ! Variables
        implicit none
        integer :: NoSampPoints, NoPerm, i, j, k
        double precision, dimension(NoSampPoints,NoPerm), intent(out) :: Sampling_1D
        integer, dimension(NoSampPoints) :: rp
        double precision, dimension(NoPerm,2) :: MD_Move 
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
            if (abs(linSamp(i) - linSamp2(i)) < 1D-6) then
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
    
end module CreateSnapshots