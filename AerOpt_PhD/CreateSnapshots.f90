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
        double precision :: newscore, bestscore  
        double precision, dimension(IV%NoSnap,maxDoF) :: Snaptemp
        allocate(MxDisp(maxDoF,2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CreateSnapshots "
        allocate(cond(maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CreateSnapshots "
        allocate(Snapshots(IV%NoSnap,maxDoF),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CreateSnapshots "
        allocate(MxDisp_Move(IV%DoF,2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CreateSnapshots "
        
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
        call LHS(Snapshots, IV%NoSnap, IV%NoCN)     ! Output: Snapshots - an initial Sampling via LHS        
        bestscore = score(Snapshots)
        do i = 1, 1000
            call LHS(Snaptemp, IV%NoSnap, IV%NoCN)
            
            newscore = score(Snaptemp)
            if (newscore > bestscore) then
                Snapshots = Snaptemp
                bestscore = newscore
            end if
        end do
        
        ! Include initial Geometry as Snapshot
        Snapshots(1,:) = 0.0

    end subroutine SubCreateSnapshots
    
    subroutine LHS(Sampling, NoSampPoints, NoPerm)
    
        ! Variables
        implicit none
        integer :: NoSampPoints, NoPerm, i, j
        double precision, dimension(NoSampPoints, maxDoF) :: Sampling
        double precision, dimension(:,:), allocatable :: Sampling_1D
        double precision :: max, min
        integer, dimension(NoSampPoints) :: rp
        double precision, dimension(NoSampPoints) :: linSamp
        
        ! Body of LHS
        Sampling = 0.0
        do j = 1, maxDoF
            
            max = MxDisp(j,1) - (MxDisp(j,1) - MxDisp(j,2))/(2*NoSampPoints)
            min = MxDisp(j,2) + (MxDisp(j,1) - MxDisp(j,2))/(2*NoSampPoints)
            linSamp = linSpacing(max, min, NoSampPoints) ! equal spacing of Design Space/Movement Domain in 1D
            
            call randperm(NoSampPoints, rp) ! generates a vector of random integer values       
            do i = 1, IV%NoSnap             
                Sampling(i,j) = linSamp(rp(i)) ! Randomly permuted integers applied as indices (rp)                               
            end do
            
        end do
         
    end subroutine LHS
    
    function score(Sampling)
    
        ! Variables
        implicit none
        integer :: a, i, j, k
        double precision :: dist, score
        double precision, dimension(:), allocatable :: SnapDist
        double precision, dimension(IV%NoSnap, maxDoF) :: Sampling
        
        ! Body of score
        allocate(SnapDist(IV%NoSnap**2/2-(IV%NoSnap/2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in CreateSnapshots "
        a = 0
        dist = 0
        do i = 1, IV%NoSnap
            do j = i + 1, IV%NoSnap
                do k = 1, IV%NoCN                   
                    dist = dist + (Sampling(i,k) - Sampling(j,k))**2
                end do
                a = a + 1
                SnapDist(a) = sqrt(dist)
            end do
        end do
        score = minval(SnapDist)
        
    end function score
    
end module CreateSnapshots