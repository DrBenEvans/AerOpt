module CreateInitialNests
    
        use Toolbox   
        real :: rn                                            ! Counts random numbers
        real, dimension(:,:), allocatable :: MxDisp           ! Matrix with min/max Displacements
        real, dimension(:,:), allocatable :: MxDisp_Move      ! Matrix with only moving min/max Displacements
        integer, dimension(:), allocatable :: cond            ! Identifies zero/non-zero values
        integer :: av                                         ! Allocater Variable       
        
    contains
    
    subroutine SubCreateInitialNests(NoNests, NoDim, NoCP, xmax, ymax, zmax)
    
        ! Variables
        integer, dimension(NoCP) :: ones
        allocate(MxDisp((NoCP*NoDim),2))
        allocate(cond(NoCP*NoDim))

        
        !!****Body of SubCreateInitialNests****!!      
        do i = 1, NoCP           
            ones(i) = 1
        end do
        
        ! Initialize min/max Displacement Matrix
        ! NoDim automatically defines the size of the Matrix
        MxDisp(:,1) = (/xmax*ones, ymax*ones, zmax*ones/)
        MxDisp(:,2) = (/xmax*(-1)*ones, ymax*(-1)*ones, zmax*(-1)*ones/)
        
        write(*,*), 'Full Matrix'
        print *, MxDisp
        
        !**Reduction of MxDisp to Nonzero Values - (only for LHS-routine to reduce processing time)**!
        cond = MxDisp(:,1) /= MxDisp(:,2)
        write(*,*), 'Zero/Nonzero'
        print *, cond
        
        ! Derive size for reduced Matrix
        av = 0
        if (xmax /= 0.00) then
            av = av + 1
        end if
        if (ymax /= 0.00) then
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
        
        write(*,*), 'Reduced Matrix'
        print *, MxDisp_Move        

        ! Execute Latin Hypercube Sampling with movable min/max Displacements
        call LHS(MxDisp_Move, NoNests, NoCP)     
        ! Output: FinalSamplingPoints - an initial Sampling
        
        ! Remember to empty MxDisp_Move to safe storage
    end subroutine SubCreateInitialNests
    
    subroutine LHS(MxDisp_Move, NoNests, NoCP)
    
        ! Variables      
        integer :: ms
        real, dimension(NoCP*av,2) :: MxDisp_Move
    
        ! Body of LHS
        !!*** Based on the Dimension, the LHS is performed 1,2 or 3 times. ***!!
        !!**** The size of MxDisp_Move indicates the number of Dimensions ****!!
        ms = size(MxDisp_Move,1)
        if ( ms == NoCP) then
            
            print *, 'Case 1'
            call LHS_1D(MxDisp_Move, NoNests, NoCP)
            
        elseif (ms == NoCP*2) then
            
            print *, 'Case 2'
            call LHS_1D(MxDisp_Move(1:NoCP,:), NoNests, NoCP)
            print *, 'Part 2'
            call LHS_1D(MxDisp_Move((NoCP+1):ms,:), NoNests, NoCP)
            
        elseif (ms == NoCP*3) then
            
            print *, 'Case 3'
            call LHS_1D(MxDisp_Move(1:NoCP,:), NoNests, NoCP)
            print *, 'Part 2'
            call LHS_1D(MxDisp_Move((NoCP+1):(2*NoCP),:), NoNests, NoCP)
            print *, 'Part 3'
            call LHS_1D(MxDisp_Move((2*NoCP+1):ms,:), NoNests, NoCP)
            
        else
            print *, 'Impossible Number of Dimensions.'
            stop
        end if
        !Output: FinalSamplingPoints - an initial Sampling
        
    end subroutine LHS
    
    subroutine LHS_1D(MD_Move, NoNests, NoCP)
    
        ! Variables
        integer, dimension(NoNests) :: rp
        real, dimension(NoCP,2) :: MD_Move
        double precision, dimension(NoNests) :: linSamp, linSamp2
        real, dimension(NoNests, NoCP) :: FinalSampPoints, FinalSampPoints2
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
                FinalSampPoints(j,i) = linSamp(rp(j)) ! Randomly permuted integers applied as indices (rp)                               
            end do                
        end do
        
        ! Output to Check
        write (*,1) transpose(FinalSampPoints)
1           format(7f10.5)
    
    end subroutine LHS_1D
    
    subroutine randperm(N, p)
    
!! Source: coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2006-03/msg00748.html
!! Based on Knuth's algorithm

        integer, intent(in) :: N
        integer, dimension(:), intent(out) :: p

        integer :: temp

        p = (/ (i, i=1,N) /)

        do j=N,2,-1

            call random_number(rn)
            k = floor(j*rn) + 1

            ! exchange p(k) and p(j)
            temp = p(k)
            p(k) = p(j)
            p(j) = temp

        end do

        end subroutine randperm
    
end module CreateInitialNests