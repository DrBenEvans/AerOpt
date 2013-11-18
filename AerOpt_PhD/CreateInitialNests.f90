module CreateInitialNests
    
        real :: rn                                            ! Counts random numbers
        real, dimension(:,:), allocatable :: MxDisp           ! Matrix with min/max Displacements
        real, dimension(:,:), allocatable :: MxDisp_NonMov    ! Matrix with min/max Displacements
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
        allocate(MxDisp_NonMov(av*NoCP,2))
        
        ! Reduce Matrix
        ! MxDisp_NonMov = MxDisp(cond,:)        
        j = 1
        do i = 1, size(cond)            
            if (cond(i) == -1) then
                MxDisp_NonMov(j,:) = MxDisp(i,:)
                j = j + 1
            end if
        end do
        !**END Reduction of...**!
        
        write(*,*), 'Reduced Matrix'
        print *, MxDisp_NonMov        

        ! Execute Latin Hypercube Sampling with min/max Displacements
        call LHS(MxDisp_NonMov, NoNests, NoCP)       
        ! Output:
        
        ! Remember to empty MxDisp_NonMov to safe storage
    end subroutine SubCreateInitialNests
    
    subroutine LHS(MxDisp_NonMov, NoNests, NoCP)
    
        ! Variables
        real, dimension(NoCP*av,2) :: MxDisp_NonMov
        double precision, dimension(NoNests) :: linSamp, linLoop
        real, dimension(NoNests, NoCP) :: FinalSampPoints
        integer :: ms
        integer, dimension(NoNests) :: rp
        double precision :: max, min, first, last, dx, dlS, ds, dlL_test, dlL1, dlL2, dL, a
    
        ! Body of LHS
        ms = size(MxDisp_NonMov,1)

        if ( ms == NoCP) then
            
            ! Implement Section as Sub and Implement it in other Cases 2 or 3 times respectively
            max = MxDisp_NonMov(1,1)
            min = MxDisp_NonMov(1,2)
            dx = max - min
            do i = 0, (NoNests - 1)
                a = i/ real(NoNests - 1)
                linSamp(i+1) = max - dx*a              
            end do
            dlS = linSamp(1) - linSamp(2)
            ds = dlS/100000
            do i = 1, 100000
                
                first = max - ds*i
                last = min + ds*i
                do j = 0, (NoNests - 1)
                    a = j/ real(NoNests - 1)  
                    linLoop(j+1) = first - (first - last)*a
                end do
                dlL1 = max - linLoop(1)
                dlL2 = (linLoop(3) - linLoop(2))/2
                dlL_test = (linLoop(2) - linLoop(1))/2
                if ((dlL1 + dlL2) > 10**(-6)) then
                    dL = max - linLoop(1)
                    exit
                end if
                
            end do
            
            first = max - dL
            last = min + dL
            do i = 0, (NoNests - 1)
                a = i/ real(NoNests - 1)  
                linLoop(i+1) = first - (first - last)*a
            end do
            
            print *, size(linLoop)
            print *, NoNests
            
            do i = 1, NoCP               
                call randperm(NoNests, rp)
                do j = 1, NoNests                   
                    FinalSampPoints(j,i) = linLoop(rp(j))
                end do
                
            end do
            
            print *, 'Case 1'
            print *, FinalSampPoints
            
        elseif (ms == NoCP*2) then
            print *, 'Case 2'
        elseif (ms == NoCP*3) then
            print *, 'Case 3'
        else
            print *, 'Impossible Number of Dimensions.'
            stop
        end if
    
    end subroutine LHS
    
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