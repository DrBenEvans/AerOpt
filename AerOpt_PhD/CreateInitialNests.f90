module CreateInitialNests
    
        use Toolbox   
        real :: rn                                            ! Counts random numbers
        real, dimension(:,:), allocatable :: MxDisp           ! Matrix with min/max Displacements
        real, dimension(:,:), allocatable :: MxDisp_Move      ! Matrix with only moving min/max Displacements
        integer, dimension(:), allocatable :: cond            ! Identifies zero/non-zero values
        real, dimension(:,:), allocatable :: InitialNests     ! Matrix containing the initial Nests
        integer :: av                                         ! Allocater Variable
        integer, dimension(:,:), allocatable :: boundff       ! Boundary faces including boundary flags
        
    contains
    
    subroutine SubCreateInitialNests(NoNests, NoDim, NoCP, xmax, ymax, zmax)
    
        ! Variables
        integer, dimension(NoCP) :: ones                      ! Vector with Ones                                    
        allocate(MxDisp((NoCP*NoDim),2))
        allocate(cond(NoCP*NoDim))
        allocate(InitialNests(NoNests,NoDim*NoCP))

        
        !!****Body of SubCreateInitialNests****!!      
        ones = (/ (1, i=1,NoCP) /)
        
        ! Initialize min/max Displacement Matrix
        ! NoDim automatically defines the size of the Matrix
        MxDisp(:,1) = (/xmax*ones, ymax*ones, zmax*ones/)
        MxDisp(:,2) = (/xmax*(-1)*ones, ymax*(-1)*ones, zmax*(-1)*ones/)
        
        !**Reduction of MxDisp to Nonzero Values - (only for LHS-routine to reduce processing time)**!
        cond = MxDisp(:,1) /= MxDisp(:,2)
        
        ! Derive size for reduced Matrix
        av = 0
        if (xmax /= 0.00) then
            av = av + 1
        end if
        if (ymax /= 0.00 .and. NoDim > 1) then
            av = av + 1
        end if
        if (zmax /= 0.00 .and. NoDim == 3) then
            av = av + 1
        end if
        allocate(MxDisp_Move(av*NoCP,2))
        
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
        call LHS(av, MxDisp_Move, NoNests, NoCP, NoDim, xmax, ymax, zmax)     
        !Output: InitialNests - an initial Sampling via LHS
        InitialNests(1,:) = (/ (0, i=1,(NoCP*NoDim)) /)

    end subroutine SubCreateInitialNests
    
    subroutine LHS(av, MxDisp_Move, NoNests, NoCP, NoDim, xmax, ymax, zmax)
    
        ! Variables      
        integer :: ms, av
        real, dimension(NoCP*av,2) :: MxDisp_Move
        real, dimension(NoNests,NoCP) :: InitialNests_1D
        integer, dimension(NoNests) :: zeros
        
        ! Body of LHS
        !!*** Based on the degrees of freedom(xmax, ymax & zmax definition ****!!
        !!*** & Dimension restriction) the LHS is performed 1,2 or 3 times. ***!!
        !!**** The size of MxDisp_Move indicates the degrees of freedom ****!!
        zeros = (/ (0, i=1,NoNests) /)
        ms = size(MxDisp_Move,1)
        if ( ms == NoCP) then
            
            print *, 'Case 1'
            call LHS_1D(MxDisp_Move, NoNests, NoCP,InitialNests_1D) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
                       
            ! Based on the Number of Dimensions, NonMoving Columns need to be filled with zeros
            j = 1
            do i = 1, size(cond)            
                if (cond(i) == -1) then
                    InitialNests(:,i) = InitialNests_1D(:,j)
                    j = j + 1
                else
                    InitialNests(:,i) =  zeros ! Non Moving
                end if
            end do            
            
        elseif (ms == NoCP*2) then
            
            ! All possible combinations of Number of Dimensions and Moving/Non-Moving restrictions are considered.
            print *, 'Case 2, Part 1'
            call LHS_1D(MxDisp_Move(1:NoCP,:), NoNests, NoCP,InitialNests_1D) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            if (xmax == 0) then
                InitialNests(:,(NoCP+1):ms) = InitialNests_1D ! Write Data in y-Matrix section
                do i = 1, NoCp
                InitialNests(:,i) = zeros ! Non Moving
                end do
            else
                InitialNests(:,1:NoCP) = InitialNests_1D ! Write Data in x-Matrix section
            end if
            
            print *, 'Case 2,Part 2'
            call LHS_1D(MxDisp_Move((NoCP+1):ms,:), NoNests, NoCP,InitialNests_1D) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            if (xmax == 0) then     ! x: 0, y & z full
                InitialNests(:,(2*NoCP+1):3*NoCP) = InitialNests_1D ! Write Data in z-Matrix section                         
            elseif (ymax == 0) then ! y: 0, x & z full
                do i = 1, NoCp
                    InitialNests(:,(NoCP+i)) = zeros ! Non Moving
                end do
                InitialNests(:,(2*NoCP+1):3*NoCP) = InitialNests_1D ! Write Data in z-Matrix section
            elseif (NoDim == 2) then ! just 2D - x & y full
                InitialNests(:,(NoCP+1):ms) = InitialNests_1D ! Write Data in y-Matrix section
            else                     ! z: 0, x & y full  
                InitialNests(:,(NoCP+1):ms) = InitialNests_1D ! Write Data in y-Matrix section
                do i = 1, NoCp
                    InitialNests(:,(2*NoCP+i)) = zeros ! Non Moving
                end do
            end if    
                    
        elseif (ms == NoCP*3) then
        
            ! Only one combination possible --> Fixed
            print *, 'Case 3, Part 1'
            call LHS_1D(MxDisp_Move(1:NoCP,:), NoNests, NoCP,InitialNests_1D) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            InitialNests(:,1:NoCP) = InitialNests_1D
            print *, 'Case 3, Part 2'
            call LHS_1D(MxDisp_Move((NoCP+1):(2*NoCP),:), NoNests, NoCP,InitialNests_1D) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            InitialNests(:,(NoCP+1):(NoCP*2)) = InitialNests_1D
            print *, 'Case 3, Part 3'
            call LHS_1D(MxDisp_Move((2*NoCP+1):ms,:), NoNests, NoCP,InitialNests_1D) ! Output: Initial_Nests_1D: Sampling Points for one Dimension
            InitialNests(:,(NoCP*2+1):ms) = InitialNests_1D
            
        else
            print *, 'Impossible Number of Dimensions.'
            stop
        end if
        !Output: InitialNests - an initial Sampling vis LHS
        
    end subroutine LHS
    
    subroutine LHS_1D(MD_Move, NoNests, NoCP,InitialNests_1D)
    
        ! Variables
        real, dimension(NoNests,NoCP), intent(out) :: InitialNests_1D
        integer, dimension(NoNests) :: rp
        real, dimension(NoCP,2) :: MD_Move 
        double precision, dimension(NoNests) :: linSamp, linSamp2
        double precision :: max, min, first, last, dlS, ds, dlL1, dlL2, dL, dbound
    
        ! Body of LHS_1D
        !!************** A complicated way of applying LHS for a 1D sampling. *************!!
        !!********* The idea is to maximize the minimum distance between points. **********!!
        !!****** However, that only plays a role in case of non-uniform distribution. *****!!
        !! Here, a uniform distribution (linSpacing) by picking the Midpoints is selected. !!        
        max = MD_Move(1,1)
        min = MD_Move(1,2)
        linSamp = linSpacing(max, min, NoNests) ! linear splitting of Design Space/Movement Domain
        
        !! Find maximum minimum distance between Points/Nests to redefine limits
        dlS = linSamp(1) - linSamp(2)
        ds = dlS/100000
        do i = 1, 100000
                
            first = max - ds*i
            last = min + ds*i
            linSamp = linSpacing(first, last, NoNests)
                
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
        linSamp = linSpacing(first, last, NoNests)
        
        ! Simplified version: Redefine Limits(max/min) via Midpoint Calculation        
        linSamp2 = linSpacing(max, min, (NoNests + 1))
        dbound = abs((linSamp2(1) - linSamp2(2))/2)
        max = max - dbound
        min = min + dbound
        linSamp2 = linSpacing(max, min, NoNests) ! Calculate Spacing with new limits
        
        ! Test simplified version vs complicated version
        k = 0
        do i = 1, NoNests
            if (abs(linSamp(i) - linSamp2(i))/xmax < 1D-6) then
                    k = k + 1                
            end if
        end do
        if (k == NoNests) then
            print *, 'Simplified Version WORKS'
        else
            print *, 'Simplified Version FAILED'
            print *, k
        end if
        
        ! Execute Random Permutation of distributed Nests for each Control Point
        do i = 1, NoCP               
            call randperm(NoNests, rp) ! see Subroutine, Output are integers
            do j = 1, NoNests               
                InitialNests_1D(j,i) = linSamp(rp(j)) ! Randomly permuted integers applied as indices (rp)                               
            end do                
        end do        
    
    end subroutine LHS_1D
    
    subroutine IdentifyBoundaryFlags()
    
        ! Variables
        use ReadData
        real, dimension(2) :: point
    
        ! Body of IdentifyBoundaryFlags
        do k = 1, nbf
            point(1) = (coord(boundf(k,1),1)*15 + coord(boundf(k,2),1)*15)/2.0
            point(2) = (coord(boundf(k,1),2)*15 + coord(boundf(k,2),2)*15)/2.0
            if (point(1) < 20 .and. point(1) > 0 .and. point(2) < 3 .and. point(2) > (-2)) then
                boundff(k,3) = 6    ! Adiabatic Viscous Wall
            elseif (point(1) == 0 .and. point(2) < 1 .and. point(2) > 0) then
                boundff(k,3) = 8    ! Engine Inlet
            else    
                boundff(k,3) = 3    ! Far Field
            end if       
        end do
    
    end subroutine IdentifyBoundaryFlags
    
end module CreateInitialNests