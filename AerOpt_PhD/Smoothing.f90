    module Smoothing
    
    use ReadData
    use InputData
    use Toolbox
    use CreateSnapshots
    use FDGD
    
    contains
    
    subroutine SmoothingLinear(CNDisp)
    
        ! Variables
        implicit none
        integer :: ismooth, NoSweeps, NoPmove, NoIterFDGD, isweep, i, NoSweepsinit
        double precision, dimension(:), allocatable :: ddn_before, ddn_after, nx, ny, sx, sy, x, y, xnew, ynew, CNxfinal, CNyfinal, CNxsmooth, CNysmooth, beta1, beta2
        double precision, dimension(maxDoF) :: CNDisp
        double precision :: conv, initialResidual, res, convergence, smoothfactor
        logical ::  shapeenclosed
        
        ! Body of SubSmoothing
        convergence = -3
        initialResidual = 0.0
        conv = 0.0
        NoPmove = size(MovingGeomIndex, dim = 1)
        smoothfactor = 0.1
        NoSweepsinit = NoPmove*smoothfactor
        smoothfactor = 1
        shapeenclosed = 1
        allocate(x(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(y(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(xnew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ynew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNxfinal(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNyfinal(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNxsmooth(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNysmooth(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(nx(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ny(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sx(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sy(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(beta1(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(beta2(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        x = RD%coord(orderedBoundaryIndex,1)
        y = RD%coord(orderedBoundaryIndex,2)
        
        ! CN of smoothed surface
        CNxsmooth = x(CN_ind)
        CNysmooth = y(CN_ind)
        
        ! Final CN
        CNxfinal = x(CN_ind) + CNDisp(1:IV%NoCN)
        CNyfinal = y(CN_ind) + CNDisp((IV%NoCN+1):(2*IV%NoCN))
        
        do while (conv > convergence)            
            
            ! 1. Calculate Second Derivative before       
            allocate(ddn_before(NoPmove),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
            call CalcSecondDerivative(ddn_before, NoPmove, x, y, nx, ny, sx, sy, shapeenclosed, beta1, beta2)
            
            !2. Move Boundary Nodes (linear)
            call LinearMotion(x, y, CNxfinal, CNyfinal)
            
            do isweep = 1, NoSweeps
                !3. Calculate Second Derivative after
                allocate(ddn_after(NoPmove),stat=allocateStatus)
                if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
                call CalcSecondDerivative(ddn_after, NoPmove, x, y, nx, ny, sx, sy, shapeenclosed, beta1, beta2)
                
                !4. Apply smoothing
                do ismooth = 1, NoPmove
                    xnew(ismooth) = x(ismooth) + smoothfactor*(ddn_after(ismooth) - ddn_before(ismooth))*nx(ismooth)
                    ynew(ismooth) = y(ismooth) + smoothfactor*(ddn_after(ismooth) - ddn_before(ismooth))*ny(ismooth)  
                end do
                x = xnew
                y = ynew
            end do
            
            !5. Check for convergence 
            res = sqrt(sum((x(CN_ind) - CNxsmooth)**2 + (y(CN_ind) - CNysmooth)**2))
            if (initialResidual == 0) then
                initialResidual = res
            end if
            conv = log(res/initialResidual)/log(10.0)
            CNxsmooth = x(CN_ind) 
            CNysmooth = y(CN_ind)
        end do
        
        !6. if converged: Move CN locations and bound back onto desired position
        call LinearMotion(x, y, CNxfinal, CNyfinal)
        
        ! Hand over Coordinates to temporary coord
        RD%coord_temp(orderedBoundaryIndex,1) = x
        RD%coord_temp(orderedBoundaryIndex,2) = y
    
    end subroutine SmoothingLinear
    
    recursive subroutine SmoothingFDGD(CNDisp, NoIterFDGD)
    
        ! Variables
        implicit none
        integer :: ismooth, NoSweeps, NoPmove, NoIterFDGD, isweep, intersect, i, NoSweepsinit
        double precision, dimension(:), allocatable :: ddn_before, ddn_after, nx, ny, sx, sy, x, y, xnew, ynew, CNxsmooth, CNysmooth, beta1, beta2, smoothing
        double precision, dimension(maxDoF) :: CNDisp, zeros
        double precision :: conv, initialResidual, res, convergence, betamean, betamin, magnitude
        double precision, dimension(2) :: a, b
        logical :: shapeenclosed
        
        ! Body of SubSmoothing
        zeros = 0.0
        convergence = -3
        initialResidual = 0.0
        conv = 0.0
        NoPmove = size(MovingGeomIndex, dim = 1)
        NoSweepsinit = NoPmove*maxval(IV%smoothfactor(1:IV%NoCN))
        shapeenclosed = .false.
        allocate(x(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(y(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(xnew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ynew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNxsmooth(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNysmooth(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(nx(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ny(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sx(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sy(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ddn_before(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ddn_after(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(beta1(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(beta2(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(smoothing(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "

        ! Extract actual boundary order from boundary faces
        x = RD%coord(orderedBoundaryIndex,1)
        y = RD%coord(orderedBoundaryIndex,2)

        ! CN of smoothed surface
        CNxsmooth = x(CN_indordered)
        CNysmooth = y(CN_indordered)
 
        do while (conv > convergence)

            ! 1. Calculate Second Derivative before       
            call CalcSecondDerivative(ddn_before, NoPmove, x, y, nx, ny, sx, sy, shapeenclosed, beta1, beta2)
            
            !2. Pre-Processing of starting Geometry for FDGD
            call PreMeshingBoundary()
            
            !3. Set new DelaunayCoord
            call RelocateCN(CNDisp, NoIterFDGD)
            
            !4. Move Boundary Mesh via FDGD
            call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1)) 
            x = RD%coord_temp(orderedBoundaryIndex,1)
            y = RD%coord_temp(orderedBoundaryIndex,2)
            
            ! Calculate edge lengths of each element on left(beta1) and right(beta2) side
            call calcBeta(x, y, NoPmove, shapeenclosed, betamin, betamean, beta1, beta2)
    
            ! Check for badly conditioned elements
            if (betamin*20 < betamean) then
                
                ! Identify worst element and modify
                call findPoint(x, y, NoPmove, betamin, beta1, beta2)
                RD%coord_temp(orderedBoundaryIndex,1) = x
                RD%coord_temp(orderedBoundaryIndex,2) = y
            
                ! Move Boundary Mesh to initial position
                call PreMeshingBoundary()
                call RelocateCN(zeros, NoIterFDGD)
                call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1)) 
                x = RD%coord_temp(orderedBoundaryIndex,1)
                y = RD%coord_temp(orderedBoundaryIndex,2)
                
                ! Re-calculate second derivative after element modification
                call CalcSecondDerivative(ddn_before, NoPmove, x, y, nx, ny, sx, sy, shapeenclosed, beta1, beta2)
                
                ! Move Boundary Mesh back to desired position
                call RelocateCN(CNDisp, NoIterFDGD)
                call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1)) 
                x = RD%coord_temp(orderedBoundaryIndex,1)
                y = RD%coord_temp(orderedBoundaryIndex,2)
                
                ! Re-calculate edge lengths
                call calcBeta(x, y, NoPmove, shapeenclosed, betamin, betamean, beta1, beta2)
                
            end if
 
            ! Number of sweeps and the individual smoothing is pre-calculated and fixed during the smoothing
            ! This is not exact, but accurate enough
            nosweeps = floor(nosweepsinit*(betamean/betamin)**2)
            call calcSmoothing(x, y, NoPmove, IV%smoothfactor(1:IV%NoCN)/maxval(IV%smoothfactor(1:IV%NoCN)), smoothing, shapeenclosed)
                
            do isweep = 1, NoSweeps
                
                !7. Calculate Second Derivative after
                ddn_after = 0.0
                beta1 = 0.0
                beta2 = 0.0
                call CalcSecondDerivative(ddn_after, NoPmove, x, y, nx, ny, sx, sy, shapeenclosed, beta1, beta2)
                betamin = minval((/minval(beta1),minval(beta2)/))

                !9. Apply smoothing
                do ismooth = 1, NoPmove
                    xnew(ismooth) = x(ismooth) + betamin**2/2.0*smoothing(ismooth)*(ddn_after(ismooth) - ddn_before(ismooth))*nx(ismooth)
                    ynew(ismooth) = y(ismooth) + betamin**2/2.0*smoothing(ismooth)*(ddn_after(ismooth) - ddn_before(ismooth))*ny(ismooth)  
                end do
                x = xnew
                y = ynew
    
            end do
            RD%coord_temp(orderedBoundaryIndex,1) = x
            RD%coord_temp(orderedBoundaryIndex,2) = y
    
            !10. Check for convergence 
            res = sqrt(sum((x(CN_indordered) - CNxsmooth)**2 + (y(CN_indordered) - CNysmooth)**2))
            if (abs(initialResidual - res) < 10e-10) then
                conv = -4
            end if 
            if (initialResidual == 0) then
                initialResidual = res
            end if
            if (conv /= -4) then
                conv = log(res/initialResidual)/log(10.0)
            end if
                
            CNxsmooth = x(CN_indordered) 
            CNysmooth = y(CN_indordered)
        end do
        
        !8. if converged: Move CN locations and bound back onto desired position
        call PreMeshingBoundary()
        call RelocateCN(CNDisp, NoIterFDGD)
        call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))
                
        ! Check for valid background mesh
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        call CheckforIntersections(DelaunayCoordDomain, DelaunayElemDomain, intersect)
        
        ! Move Domain Nodes
        if (intersect == 1) then
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                CNDisp = CNDisp*(1.0/NoIterFDGD)
        else
            NoIterFDGD = NoIterFDGD*2
            call SmoothingFDGD(CNDisp, NoIterFDGD)
        end if
    
    end subroutine SmoothingFDGD
    
    
    
    
    
    
    subroutine CalcSecondDerivative(ddn, NoPmove, x, y, nx, ny, sx, sy, shapeenclosed, beta1, beta2)
    
        ! Variables
        implicit none
        integer :: ip, NoPmove, i
        double precision :: magnitude, dx1, dx2
        logical :: shapeenclosed
        double precision, dimension(:), allocatable :: x, y, ddn, sx, sy, nx, ny
        double precision, dimension(NoPmove) :: beta1, beta2
        
        if (shapeenclosed == .true.) then
            call CalcSecondDerivativeclosed(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)
        else
            call CalcSecondDerivativeopen(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)
        end if
    
    end subroutine CalcSecondDerivative
    
    subroutine CalcSecondDerivativeclosed(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)
    
        ! Variables
        implicit none
        integer :: ip, NoPmove, i
        double precision :: magnitude, dx1, dx2
        double precision, dimension(:), allocatable :: x, y, ddn, sx, sy, nx, ny
        double precision, dimension(NoPmove) :: beta1, beta2
    
        ! Body of CalcSecondDerivative                
        ! ComputeNormal
        nx(1) = -(y(2)-y(NoPmove))
        ny(1) = x(2)-x(NoPmove)
        do ip = 2, NoPmove-1
            nx(ip) = -(y(ip+1)-y(ip-1))
            ny(ip) = x(ip+1)-x(ip-1)
        end do
        nx(NoPmove) = -(y(1)-y(NoPmove-1))
        ny(NoPmove) = x(1)-x(NoPmove-1)

        ! Normalize
        do i = 1, NoPmove
            magnitude = sqrt(nx(i)**2 + ny(i)**2)
            nx(i) = nx(i)/magnitude
            ny(i) = ny(i)/magnitude
        end do
            
        ! Compute the s vectors
        sx = ny
        sy = -nx
             
        ! Compute Gradient
        dx1 = (x(1) - x(NoPmove))*sx(1) + (y(1) - y(NoPmove))*sy(1)
        dx2 = (x(2) - x(1))*sx(1) + (y(2) - y(1))*sy(1)
        ddn(1) = 2*(dx1*(x(2)*nx(1) + y(2)*ny(1)) - &
        (dx1+dx2)*(x(1)*nx(1) + y(1)*ny(1)) + &
        dx2*(x(NoPmove)*nx(1) + y(NoPmove)*ny(1)))/(dx1*dx2*(dx1 + dx2))
        beta1(1) = dx1
        beta2(1) = dx2
        do i = 2, NoPmove-1
            dx1 = (x(i) - x(i-1))*sx(i) + (y(i) - y(i-1))*sy(i)
            dx2 = (x(i+1) - x(i))*sx(i) + (y(i+1) - y(i))*sy(i)
            ddn(i) = 2*(dx1*(x(i+1)*nx(i) + y(i+1)*ny(i)) - &
            (dx1+dx2)*(x(i)*nx(i) + y(i)*ny(i)) + &
            dx2*(x(i-1)*nx(i) + y(i-1)*ny(i)))/(dx1*dx2*(dx1 + dx2))  
            beta1(i) = dx1
            beta2(i) = dx2
        end do
        dx1 = (x(NoPmove) - x(NoPmove-1))*sx(NoPmove) + (y(NoPmove) - y(NoPmove-1))*sy(NoPmove)
        dx2 = (x(1) - x(NoPmove))*sx(NoPmove) + (y(1) - y(NoPmove))*sy(NoPmove)
        ddn(NoPmove) = 2*(dx1*(x(1)*nx(NoPmove) + y(1)*ny(NoPmove)) - &
        (dx1+dx2)*(x(NoPmove)*nx(NoPmove) + y(NoPmove)*ny(NoPmove)) + &
        dx2*(x(NoPmove-1)*nx(NoPmove) + y(NoPmove-1)*ny(NoPmove)))/(dx1*dx2*(dx1 + dx2))
        beta1(NoPmove) = dx1
        beta2(NoPmove) = dx2
    
    end subroutine CalcSecondDerivativeclosed
    
    subroutine CalcSecondDerivativeopen(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)
    
        ! Variables
        implicit none
        integer :: ip, NoPmove, i
        double precision :: magnitude, dx1, dx2
        double precision, dimension(:), allocatable :: x, y, ddn, sx, sy, nx, ny
        double precision, dimension(NoPmove) :: beta1, beta2
    
        ! Body of CalcSecondDerivative                
        ! ComputeNormal
        do ip = 2, NoPmove-1
            nx(ip) = -(y(ip+1)-y(ip-1))
            ny(ip) = x(ip+1)-x(ip-1)
        end do
        nx(1) = nx(2)
        ny(1) = ny(2)
        nx(NoPmove) = nx(NoPmove-1)
        ny(NoPmove) = ny(NoPmove-1)

        ! Normalize
        do i = 1, NoPmove
            magnitude = sqrt(nx(i)**2 + ny(i)**2)
            nx(i) = nx(i)/magnitude
            ny(i) = ny(i)/magnitude
        end do
            
        ! Compute the s vectors
        sx = ny
        sy = -nx
             
        ! Compute Gradient
        ddn(1) = 0
        beta1(1) = (x(2) - x(1))*sx(1) + (y(2) - y(1))*sy(1)
        beta2(1) = beta1(1)
        do i = 2, NoPmove-1
            dx1 = (x(i) - x(i-1))*sx(i) + (y(i) - y(i-1))*sy(i)
            dx2 = (x(i+1) - x(i))*sx(i) + (y(i+1) - y(i))*sy(i)
            ddn(i) = 2*(dx1*(x(i+1)*nx(i) + y(i+1)*ny(i)) - &
            (dx1+dx2)*(x(i)*nx(i) + y(i)*ny(i)) + &
            dx2*(x(i-1)*nx(i) + y(i-1)*ny(i)))/(dx1*dx2*(dx1 + dx2))  
            beta1(i) = dx1
            beta2(i) = dx2
        end do
        beta1(NoPmove) = (x(NoPmove) - x(NoPmove-1))*sx(NoPmove) + (y(NoPmove) - y(NoPmove-1))*sy(NoPmove)
        beta2(NoPmove) = beta1(NoPmove)
        ddn(NoPmove) = 0
    
    end subroutine CalcSecondDerivativeopen
    
    subroutine calcBeta(x, y, NoPmove, shapeenclosed, betamin, betamean, beta1, beta2)
    
        ! Variables
        implicit none
        double precision :: betamin, betamean
        double precision, dimension(NoPmove) :: x, y, beta1, beta2
        integer :: NoPmove
        logical :: shapeenclosed
        
        ! Body of calcBeta
        if (shapeenclosed == .true.) then
            call calcBetaclosed(x, y, NoPmove, shapeenclosed, betamin, betamean, beta1, beta2)
        else
            call calcBetaopen(x, y, NoPmove, shapeenclosed, betamin, betamean, beta1, beta2)
        end if
    
    end subroutine calcBeta
    
    subroutine calcBetaclosed(x, y, NoPmove, shapeenclosed, betamin, betamean, beta1, beta2)
    
        ! Variables
        implicit none
        double precision :: betamin, betamean, magnitude
        logical :: shapeenclosed
        double precision, dimension(NoPmove) :: x, y, beta1, beta2, sx, sy
        integer :: NoPmove, i
        
        ! Body of calcBeta
        ! Compute tangential direction (s) and normalize
        sx(1) = x(2)-x(NoPmove)
        sy(1) = y(2)-y(NoPmove)
        do i = 2, NoPmove-1
            sx(i) = x(i+1)-x(i-1)
            sy(i) = y(i+1)-y(i-1)
        end do
        sx(NoPmove) = x(1)-x(NoPmove-1)
        sy(NoPmove) = y(1)-y(NoPmove-1)
               
        do i = 1, NoPmove ! Normalize
            magnitude = sqrt(sx(i)**2 + sy(i)**2)
            sx(i) = sx(i)/magnitude
            sy(i) = sy(i)/magnitude
        end do
            
        ! Calculate beta and derive the necessary number of sweeps
        beta1(1) = dot_product((/ x(1)-x(NoPmove), y(1)-y(NoPmove) /), (/sx(1), sy(1)/))
        beta2(1) = dot_product((/ x(2)-x(1), y(2)-y(1) /), (/sx(1), sy(1)/))
        do i = 2, NoPmove-1
            beta1(i) = dot_product((/x(i)-x(i-1), y(i)-y(i-1)/),(/sx(i), sy(i)/))
            beta2(i) = dot_product((/x(i+1)-x(i), y(i+1)-y(i)/),(/sx(i), sy(i)/))
        end do
        beta1(NoPmove) = dot_product((/x(NoPmove)-x(NoPmove-1), y(NoPmove)-y(NoPmove-1)/),(/sx(NoPmove), sy(NoPmove)/))
        beta2(NoPmove) = dot_product((/x(1)-x(NoPmove), y(1)-y(NoPmove)/),(/sx(NoPmove), sy(NoPmove)/))
        betamin = minval((/minval(beta1), minval(beta2)/))
        betamean = (sum(beta1)+sum(beta2))/(size(beta1, dim = 1)+size(beta2, dim = 1))
    
    end subroutine calcBetaclosed
    
    subroutine calcBetaopen(x, y, NoPmove, shapeenclosed, betamin, betamean, beta1, beta2)

        ! Variables
        implicit none
        double precision :: betamin, betamean, magnitude
        logical :: shapeenclosed
        double precision, dimension(NoPmove) :: x, y, beta1, beta2, sx, sy
        integer :: NoPmove, i

        ! Body of calcBeta
        ! Compute tangential direction (s) and normalize
        sx(1) = x(3)-x(1)
        sy(1) = y(3)-y(1)
        do i = 2, NoPmove-1
            sx(i) = x(i+1)-x(i-1)
            sy(i) = y(i+1)-y(i-1)
        end do
        sx(NoPmove) = x(NoPmove)-x(NoPmove-2)
        sy(NoPmove) = y(NoPmove)-y(NoPmove-2)

        do i = 1, NoPmove ! Normalize
            magnitude = sqrt(sx(i)**2 + sy(i)**2)
            sx(i) = sx(i)/magnitude
            sy(i) = sy(i)/magnitude
        end do

        ! Calculate beta and derive the necessary number of sweeps
        do i = 2, NoPmove-1
            beta1(i) = dot_product((/x(i)-x(i-1), y(i)-y(i-1)/),(/sx(i), sy(i)/))
            beta2(i) = dot_product((/x(i+1)-x(i), y(i+1)-y(i)/),(/sx(i), sy(i)/))
        end do
        beta1(1) = beta1(2)
        beta2(1) = beta2(2)
        beta1(NoPmove) = beta1(NoPmove-1)
        beta2(NoPmove) = beta2(NoPmove-1)
        betamin = minval((/minval(beta1), minval(beta2)/))
        betamean = (sum(beta1)+sum(beta2))/(size(beta1, dim = 1)+size(beta2, dim = 1))
    
    end subroutine calcBetaopen
    
    subroutine LinearMotion(x, y, CNxfinal, CNyfinal)
    
        ! Variables
        implicit none
        integer :: i, j
        double precision :: CNdist
        double precision, dimension(maxDoF) :: CNDisp
        double precision, dimension(:), allocatable :: x, y, CNxfinal, CNyfinal
    
        ! Body of LinearMotion
        do i = 1, IV%NoCN
            CNDisp(i) = CNxfinal(i) - x(CN_ind(i))
            CNDisp(i+IV%NoCN) = CNyfinal(i) - y(CN_ind(i))
            x(CN_ind(i)) = x(CN_ind(i)) + CNDisp(i)
            y(CN_ind(i)) = y(CN_ind(i)) + CNDisp(i+IV%NoCN)
        end do
        
        ! Move boundary points linearly
        do i = 1, IV%NoCN          
            if (i /= IV%NoCN) then
                CNdist = abs(x(CN_ind(i+1)) - x(CN_ind(i)))
                do j = CN_ind(i)+1, CN_ind(i+1)-1
                    y(j) = y(j) + CNDisp(i+IV%NoCN) - abs(x(CN_ind(i)) - x(j))*(1/CNdist)*CNDisp(i+IV%NoCN)
                    y(j) = y(j) + CNDisp(i+IV%NoCN+1) - abs(x(CN_ind(i+1)) - x(j))*(1/CNdist)*CNDisp(i+IV%NoCN+1)
                end do
                CNdist = abs(y(CN_ind(i+1)) - y(CN_ind(i)))
                do j = CN_ind(i)+1, CN_ind(i+1)-1
                    x(j) = x(j) + CNDisp(i) - abs(y(CN_ind(i)) - y(j))*(1/CNdist)*CNDisp(i)
                    x(j) = x(j) + CNDisp(i+1) - abs(y(CN_ind(i+1)) - y(j))*(1/CNdist)*CNDisp(i+1)
                end do
            else
                CNdist = x(CN_ind(1)) - x(CN_ind(i))
                do j = CN_ind(i)+1, CN_ind(i+1)-1
                    y(j) = y(j) + CNDisp(i+IV%NoCN) - abs(x(CN_ind(i)) - x(j))*(1/CNdist)*CNDisp(i+IV%NoCN)
                    y(j) = y(j) + CNDisp(1+IV%NoCN) - abs(x(CN_ind(1)) - x(j))*(1/CNdist)*CNDisp(1+IV%NoCN)
                end do                  
                CNdist = y(CN_ind(1)) - y(CN_ind(i))
                do j = CN_ind(i)+1, CN_ind(i+1)-1
                    x(j) = x(j) + CNDisp(i) - abs(y(CN_ind(i)) - y(j))*(1/CNdist)*CNDisp(i)
                    x(j) = x(j) + CNDisp(1) - abs(y(CN_ind(1)) - y(j))*(1/CNdist)*CNDisp(1)
                end do
            end if
        end do    


    end subroutine LinearMotion
    
    subroutine calcSmoothing(x, y, NoPmove, smoothfactor, smoothing, shapeenclosed)
    
        ! Variables
        implicit none
        integer :: NoPmove
        logical :: shapeenclosed
        double precision, dimension(IV%NoCN) :: smoothfactor
        double precision, dimension(NoPmove) :: x, y, smoothing
        
        ! Body of calcAlpha
        if (shapeenclosed == .true.) then
            smoothing = smoothingclosed(x, y, NoPmove, smoothfactor)
        else
            smoothing = smoothingopen(x, y, NoPmove, smoothfactor)
        end if
    
    end subroutine calcsmoothing
    
    function smoothingclosed(x, y, NoPmove, smoothfactor)
    
        ! Variables
        implicit none
        integer :: NoPmove, i, j
        double precision :: delta
        double precision, dimension(IV%NoCN) :: smoothfactor
        double precision, dimension(NoPmove) :: x, y, smoothingclosed
        double precision, dimension(:), allocatable :: dist
    
        ! Body of smoothingopen
        allocate(dist(NoPmove+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        dist(1) = 0.0
        do i = 2, NoPmove
            dist(i) = dist(i-1) + sqrt((x(i)-x(i-1))**2 + (y(i)-y(i-1))**2) 
        end do
        dist(NoPmove+1) = dist(NoPmove) + sqrt((x(1)-x(NoPmove))**2 + (y(1)-y(NoPmove))**2)
 
        smoothingclosed = 0.0
        do i = 1, IV%NoCN-1
            delta = (dist(CN_indordered(i+1))-dist(CN_indordered(i)))
            do j = CN_indordered(i)+1, CN_indordered(i+1)
                smoothingclosed(j) = smoothfactor(i)*(1 - (dist(j) - dist(CN_indordered(i)))/delta) + smoothfactor(i+1)*((dist(j) - dist(CN_indordered(i)))/delta)
            end do
        end do
        
        delta = (dist(NoPmove+1)-dist(CN_indordered(IV%NoCN)) + dist(CN_indordered(1)))
        do j = CN_indordered(IV%NoCN)+1, NoPmove
            smoothingclosed(j) = smoothfactor(IV%NoCN)*(1 - (dist(j) - dist(CN_indordered(IV%NoCN)))/delta) + smoothfactor(1)*((dist(j) - dist(CN_indordered(IV%NoCN)))/delta)
        end do
        do j = 1, CN_indordered(1)
            smoothingclosed(j) = smoothfactor(IV%NoCN)*(1 - (dist(j) + dist(NoPmove+1) - dist(CN_indordered(IV%NoCN)))/delta) + smoothfactor(1)*((dist(j) + dist(NoPmove+1) - dist(CN_indordered(IV%NoCN)))/delta)
        end do
    
    end function smoothingclosed
    
    function smoothingopen(x, y, NoPmove, smoothfactor)
    
        ! Variables
        implicit none
        integer :: NoPmove, i, j
        double precision :: delta
        double precision, dimension(IV%NoCN) :: smoothfactor
        double precision, dimension(NoPmove) :: x, y, smoothingopen
        double precision, dimension(:), allocatable :: dist
    
        ! Body of smoothingopen
        allocate(dist(NoPmove+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        do i = 2, NoPmove
            dist(i) = dist(i-1) + sqrt((x(i)-x(i-1))**2 + (y(i)-y(i-1))**2) 
        end do

        smoothingopen = 0.0
        do i = 1, IV%NoCN-1
            delta = (dist(CN_indordered(i+1))-dist(CN_indordered(i)))
            do j = CN_indordered(i)+1, CN_indordered(i+1)
                smoothingopen(j) = smoothfactor(i)*(1 - (dist(j) - dist(CN_indordered(i)))/delta) + smoothfactor(i+1)*((dist(j) - dist(CN_indordered(i)))/delta)
            end do
        end do
                
    end function smoothingopen
    
    subroutine findPoint(x, y, NoPmove, betamin, beta1, beta2)
    
        ! Variables
        implicit none
        double precision, dimension(2) :: PA, PB, PC, AC, AB, ABrot, PE
        double precision, dimension(2,2) :: T
        double precision :: CAB, ABC, gamma, betamin
        integer :: sign, minind, NoPmove
        double precision, dimension(NoPmove) :: x, y, beta1, beta2
        real, PARAMETER :: pi = 3.1415927
    
        ! Body of findPoint
        if (minval(beta1) == betamin) then
            minind = minloc(beta1, dim = 1)
            if (minind == 1) then
                PA = (/x(minind), y(minind)/)
                PB = (/x(NoPmove), y(NoPmove)/)
                PC = (/x(minind+1), y(minind+1)/)
            elseif (minind == NoPmove) then
                PA = (/x(minind), y(minind)/)              
                PB = (/x(minind-1), y(minind-1)/)
                PC = (/x(1), y(1)/)
            else
                PA = (/x(minind), y(minind)/)
                PB = (/x(minind-1), y(minind-1)/)
                PC = (/x(minind+1), y(minind+1)/)
            end if
            sign = -1
        else
            minind = minloc(beta2, dim = 1)
            if (minind == 1) then
                PA = (/x(minind), y(minind)/)
                PB = (/x(minind+1), y(minind+1)/)
                PC = (/x(NoPmove), y(NoPmove)/)
            elseif (minind == NoPmove) then
                PA = (/x(minind), y(minind)/)              
                PB = (/x(1), y(1)/)
                PC = (/x(minind-1), y(minind-1)/)
            else
                PA = (/x(minind), y(minind)/)
                PB = (/x(minind+1), y(minind+1)/)
                PC = (/x(minind-1), y(minind-1)/)              
            end if
            sign = 1
        end if
    
        AC = PC - PA
        AB = PB - PA    
        CAB = acos(dot_product(AC,AB)/(sqrt(AC(1)**2+AC(2)**2)*sqrt(AB(1)**2+AB(2)**2))) ! Angle CAB
        ABC = sign*(pi - CAB)/2 ! isoscele angle ABC
        ! Rotate Vector AB
        T(1,:) = (/cos(ABC), sin(ABC)/)
        T(2,:) = (/-sin(ABC), cos(ABC)/)
        ABrot = matmul(-AB,T)
        ! Find Intersection point E
        gamma = ((PA(2)-PB(2))*ABrot(1) - (PA(1)-PB(1))*ABrot(2))/ &
            (AC(1)*ABrot(2) - AC(2)*ABrot(1))
        PE = PA + gamma*AC
        
        if (sign == 1) then
            if (minind == 1) then
                x(NoPmove) = PE(1)
                y(NoPmove) = PE(2)
            else
                x(minind-1) = PE(1)
                y(minind-1) = PE(2)
            end if
        else
            if (minind == NoPmove) then
                x(1) = PE(1)
                y(1) = PE(2)
            else
                x(minind+1) = PE(1)
                y(minind+1) = PE(2)
            end if
        end if
        
    end subroutine findPoint
    
    end module Smoothing