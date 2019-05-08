    module Smoothing

    use ReadData
    use InputData
    use Toolbox
    use CreateSnapshots
    use FDGD

    contains

    recursive subroutine SmoothingLinear(CNDisp, NoMove, counter)

        ! Variables
        implicit none
        integer :: ismooth, NoSweeps, NoPmove, NoMove, isweep, intersect, i, NoSweepsinit, counter
        double precision, dimension(:), allocatable :: ddn_before, ddn_after, dn_before, dn_after, nx, ny, sx, sy, x, y, xnew, ynew, xbefore, ybefore, CNxfinal, CNyfinal, CNxsmooth, CNysmooth, beta1, beta2, smoothing, normx, normy
        double precision, dimension(maxDoF) :: CNDisp
        double precision :: conv, initialResidual, res, betamean, betamin, magnitude, N2CN
        double precision, dimension(2) :: a

        ! Body of SubSmoothing
        initialResidual = 0.0
        conv = 0.0
        NoPmove = size(InnerBound, dim = 1)
        N2CN = NoPmove/IV%NoCN
        NoSweepsinit = N2CN**2*maxval(IV%smoothfactor(1:IV%NoCN))
        allocate(x(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(y(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(xnew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ynew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(xbefore(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ybefore(NoPmove),stat=allocateStatus)
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
        allocate(dn_before(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(dn_after(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(beta1(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(beta2(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(smoothing(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(normx(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(normy(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNxfinal(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(CNyfinal(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "

        ! Extract actual boundary order from boundary faces
        x = RD%coord_temp(orderedBoundaryIndex,1)
        y = RD%coord_temp(orderedBoundaryIndex,2)

        ! Save for Intersection
        xbefore = x
        ybefore = y

        ! CN of smoothed surface
        CNxsmooth = x(CN_indordered)
        CNysmooth = y(CN_indordered)

        ! Final CN positions
        do i = 1, IV%NoCN
            CNxfinal(i) = x(CN_indordered(i)) + CNDisp(i)/NoMove
            CNyfinal(i) = y(CN_indordered(i)) + CNDisp(IV%NoCN+i)/NoMove
            if (IV%CNconnecttrans(i) /= 0) then
                CNxfinal(i) = x(CN_indordered(i)) + CNDisp(IV%CNconnecttrans(i))/NoMove
                CNyfinal(i) = y(CN_indordered(i)) + CNDisp(IV%NoCN+IV%CNconnecttrans(i))/NoMove
            end if
        end do

        !Pre-calculate smoothing and x/y-normalization value
        call calcSmoothing(x, y, NoPmove, IV%smoothfactor(1:IV%NoCN)/maxval(IV%smoothfactor(1:IV%NoCN)), smoothing)
        !call normXY(x, y, normx, normy, CNxfinal, CNyfinal, NoPmove)

        do while (conv > IV%smoothconvergence)

            ! Calculate Second Derivative before
            call CalcSecondDerivativeInitial(ddn_before, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)

            ! Calculate edge lengths of each element on left(beta1) and right(beta2) side
            call calcBeta(x, y, NoPmove, betamin, betamean, beta1, beta2)

            ! Move Boundary Nodes (linear)
            call quasi1DLinear(x, y, CNxfinal, CNyfinal, NoPmove)

            ! Check for badly conditioned elements
            if (betamin*20 < betamean) then

                !print*, 'small'
                !! Identify worst element and modify
                !call findPoint(x, y, NoPmove, betamin, beta1, beta2)
                !
                !! Move Boundary Nodes to initial position
                !CNxfinal = x(CN_indordered) - CNDisp(1:IV%NoCN)/NoMove
                !CNyfinal = y(CN_indordered) - CNDisp((IV%NoCN+1):(2*IV%NoCN))/NoMove
                !call quasi1DLinear(x, y, CNxfinal, CNyfinal, NoPmove)
                !
                !! Re-calculate second derivative after element modification
                !call CalcSecondDerivative(ddn_before, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)
                !
                !! Re-calculate edge lengths
                !call calcBeta(x, y, NoPmove, betamin, betamean, beta1, beta2)
                !
                !! Move Boundary Mesh back to desired position
                !CNxfinal = x(CN_indordered) + CNDisp(1:IV%NoCN)/NoMove
                !CNyfinal = y(CN_indordered) + CNDisp((IV%NoCN+1):(2*IV%NoCN))/NoMove
                !call quasi1DLinear(x, y, CNxfinal, CNyfinal, NoPmove)

            end if

            ! Number of sweeps is pre-calculated and fixed during the smoothing
            ! This is not exact, but accurate enough
            NoSweeps = floor(NoSweepsInit*(betamean/betamin)**2)
if (NoSweeps > NoSweepsInit*50) then
    print *, 'NoSweeps is too big. It will be reduced.'
    NoSweeps = NoSweepsInit*50
end if

            do isweep = 1, NoSweeps

                ! Calculate Second Derivative after
                ddn_after = 0.0
                call CalcSecondDerivative(ddn_after, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)

                ! Apply smoothing
                do ismooth = 1, NoPmove
            xnew(ismooth) = x(ismooth) + betamin**2/2.0*smoothing(ismooth)*(ddn_after(ismooth) - ddn_before(ismooth))*nx(ismooth) !*normx(ismooth) !-(dn_after(ismooth) - dn_before(ismooth)))*nx(ismooth)
                    ynew(ismooth) = y(ismooth) + betamin**2/2.0*smoothing(ismooth)*(ddn_after(ismooth) - ddn_before(ismooth))*ny(ismooth) !*normy(ismooth) !-(dn_after(ismooth) - dn_before(ismooth)))*ny(ismooth)
                end do
                x = xnew
                y = ynew

            end do
!open(1, file= 'XYcoord.dat', form='formatted', status = 'unknown')
!do isweep = 1, 100
!    write(1,'(2f22.15)') x(isweep), y(isweep)
!end do
!close(1)
            ! Check for convergence
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
print *, conv
            CNxsmooth = x(CN_indordered)
            CNysmooth = y(CN_indordered)
        end do

        ! If converged: Move CN locations and bound back onto desired position
        call quasi1DLinear(x, y, CNxfinal, CNyfinal, NoPmove)

        ! Hand over Coordinates to temporary coord
        RD%coord_temp(orderedBoundaryIndex,1) = x
        RD%coord_temp(orderedBoundaryIndex,2) = y

        ! Rotate
        do i = 1, IV%NoCN
            if (IV%angle(i) /= 0) then
                call AngleofAttack(CNDisp(4*i)/NoMove, i)
            end if
        end do

        ! Check for valid background mesh
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        call CheckforIntersections(DelaunayCoordDomain, DelaunayElemDomain, intersect, x, y)
        call CheckBoundIntersections(x, y, intersect)

        ! Move Domain Nodes
        if (intersect == 0) then
            if (counter < NoMove) then
                counter = counter + 1
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                call PreMeshingMid()
                call SmoothingLinear(CNDisp, NoMove, counter)
            else
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                if (NoMove > 1) then
                    call PreMeshingEnd()
                end if
            end if
        elseif (intersect == 1) then
            RD%coord_temp(orderedBoundaryIndex,1) = xbefore
            RD%coord_temp(orderedBoundaryIndex,2) = ybefore
            counter = 2*counter-1
            NoMove = NoMove*2
            call SmoothingLinear(CNDisp, NoMove, counter)
        else
            RD%coord_temp(orderedBoundaryIndex,1) = xbefore
            RD%coord_temp(orderedBoundaryIndex,2) = ybefore
            CNDisp = CNDisp/2.0
            call SmoothingLinear(CNDisp, NoMove, counter)
        end if

    end subroutine SmoothingLinear

    recursive subroutine SmoothingFDGD(CNDisp, NoMove, counter)

        ! Variables
        implicit none
        integer :: ismooth, NoSweeps, NoPmove, NoMove, isweep, intersect, i, NoSweepsinit, counter
        double precision, dimension(:), allocatable :: ddn_before, ddn_after, nx, ny, sx, sy, x, y, xnew, ynew, xbefore, ybefore, CNxsmooth, CNysmooth, beta1, beta2, smoothing
        double precision, dimension(maxDoF) :: CNDisp, zeros
        double precision :: conv, initialResidual, res, betamean, betamin, magnitude, N2CN
        double precision, dimension(2) :: a, b

        ! Body of SubSmoothing
        zeros = 0.0
        initialResidual = 0.0
        conv = 0.0
        NoPmove = size(InnerBound, dim = 1)
        N2CN = NoPmove/IV%NoCN
        NoSweepsinit = N2CN**2*maxval(IV%smoothfactor(1:IV%NoCN))
        allocate(x(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(y(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(xnew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ynew(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(xbefore(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ybefore(NoPmove),stat=allocateStatus)
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
        x = RD%coord_temp(orderedBoundaryIndex,1)
        y = RD%coord_temp(orderedBoundaryIndex,2)

        ! Save for Intersection
        xbefore = x
        ybefore = y

        ! CN of smoothed surface
        CNxsmooth = x(CN_indordered)
        CNysmooth = y(CN_indordered)

        do while (conv > IV%smoothconvergence)

            ! Calculate Second Derivative before
            call CalcSecondDerivativeInitial(ddn_before, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)

            ! Calculate edge lengths of each element on left(beta1) and right(beta2) side
            call calcBeta(x, y, NoPmove, betamin, betamean, beta1, beta2)

            ! Pre-Processing of starting Geometry for FDGD
            call PreMeshingBoundary()

            ! Set new DelaunayCoord
            call RelocateCN(CNDisp, NoMove, counter)

            ! Move Boundary Mesh via FDGD
            call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))
            x = RD%coord_temp(orderedBoundaryIndex,1)
            y = RD%coord_temp(orderedBoundaryIndex,2)

            ! Check for badly conditioned elements
            if (betamin*20 < betamean) then

                ! Identify worst element and modify
                call findPoint(x, y, NoPmove, betamin, beta1, beta2)
                RD%coord_temp(orderedBoundaryIndex,1) = x
                RD%coord_temp(orderedBoundaryIndex,2) = y

                ! Move Boundary Mesh to initial position
                call PreMeshingBoundary()
                call RelocateCN(zeros, NoMove, counter)
                call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))
                x = RD%coord_temp(orderedBoundaryIndex,1)
                y = RD%coord_temp(orderedBoundaryIndex,2)

                ! Re-calculate second derivative after element modification
                call CalcSecondDerivative(ddn_before, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)

                ! Re-calculate edge lengths
                call calcBeta(x, y, NoPmove, betamin, betamean, beta1, beta2)

                ! Move Boundary Mesh back to desired position
                call RelocateCN(CNDisp, NoMove, counter)
                call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))
                x = RD%coord_temp(orderedBoundaryIndex,1)
                y = RD%coord_temp(orderedBoundaryIndex,2)

            end if

            ! Number of sweeps and the individual smoothing is pre-calculated and fixed during the smoothing
            ! This is not exact, but accurate enough
            nosweeps = floor(nosweepsinit*(betamean/betamin)**2)
            call calcSmoothing(x, y, NoPmove, IV%smoothfactor(1:IV%NoCN)/maxval(IV%smoothfactor(1:IV%NoCN)), smoothing)

            do isweep = 1, NoSweeps

                ! Calculate Second Derivative after
                ddn_after = 0.0
                call CalcSecondDerivative(ddn_after, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)

                ! Apply smoothing
                do ismooth = 1, NoPmove
                    xnew(ismooth) = x(ismooth) + betamin**2/2.0*smoothing(ismooth)*(ddn_after(ismooth) - ddn_before(ismooth))*nx(ismooth)
                    ynew(ismooth) = y(ismooth) + betamin**2/2.0*smoothing(ismooth)*(ddn_after(ismooth) - ddn_before(ismooth))*ny(ismooth)
                end do
                x = xnew
                y = ynew

            end do
            RD%coord_temp(orderedBoundaryIndex,1) = x
            RD%coord_temp(orderedBoundaryIndex,2) = y

            ! Check for convergence
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

        !If converged: Move CN locations and bound back onto desired position
        call PreMeshingBoundary()
        call RelocateCN(CNDisp, NoMove, counter)
        call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))

        ! Rotate
        do i = 1, IV%NoCN
            if (IV%angle(i) /= 0) then
                call AngleofAttack(CNDisp(4*i)/NoMove, i)
            end if
        end do

        ! Check for valid background mesh
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        call CheckforIntersections(DelaunayCoordDomain, DelaunayElemDomain, intersect)
        call CheckBoundIntersections(x, y, intersect)

        ! Move Domain Nodes
        if (intersect == 0) then
            if (counter < NoMove) then
                counter = counter + 1
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                call PreMeshingMid()
                call SmoothingFDGD(CNDisp, NoMove, counter)
            else
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                if (NoMove > 1) then
                    call PreMeshingEnd()
                end if
                if (allocated(dRot) == .true.) then
                    deallocate(dRot)
                end if
            end if
        elseif (intersect == 1) then
            RD%coord_temp(orderedBoundaryIndex,1) = xbefore
            RD%coord_temp(orderedBoundaryIndex,2) = ybefore
            if (allocated(dRot) == .true.) then
                deallocate(dRot)
            end if
            NoMove = NoMove*2
            counter = 2*counter-1
            call SmoothingFDGD(CNDisp, NoMove, counter)
        else
            RD%coord_temp(orderedBoundaryIndex,1) = xbefore
            RD%coord_temp(orderedBoundaryIndex,2) = ybefore
            if (allocated(dRot) == .true.) then
                deallocate(dRot)
            end if
            CNDisp = CNDisp/2.0
            call SmoothingFDGD(CNDisp, NoMove, counter)
        end if

    end subroutine SmoothingFDGD




    subroutine CalcSecondDerivativeInitial(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)

        ! Variables
        implicit none
        integer :: ip, NoPmove, i
        double precision :: magnitude, dx1, dx2
        double precision, dimension(:), allocatable :: x, y, ddn, sx, sy, nx, ny
        double precision, dimension(NoPmove) :: beta1, beta2

        if (IV%shapeenclosed == .true.) then
            call CalcSecondDerivativeclosedInitial(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)
        else
            call CalcSecondDerivativeopenInitial(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)
        end if

    end subroutine CalcSecondDerivativeInitial

    subroutine CalcSecondDerivativeclosedInitial(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)

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

    end subroutine CalcSecondDerivativeclosedInitial

    subroutine CalcSecondDerivativeopenInitial(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)

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

    end subroutine CalcSecondDerivativeopenInitial

    subroutine CalcSecondDerivative(ddn, NoPmove, x, y, nx, ny, sx, sy, beta1, beta2)

        ! Variables
        implicit none
        integer :: ip, NoPmove, i
        double precision :: magnitude, dx1, dx2
        double precision, dimension(:), allocatable :: x, y, ddn, sx, sy, nx, ny
        double precision, dimension(NoPmove) :: beta1, beta2

        if (IV%shapeenclosed == .true.) then
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
        dx1 = beta1(1)
        dx2 = beta2(1)
        ddn(1) = 2*(dx1*(x(2)*nx(1) + y(2)*ny(1)) - &
        (dx1+dx2)*(x(1)*nx(1) + y(1)*ny(1)) + &
        dx2*(x(NoPmove)*nx(1) + y(NoPmove)*ny(1)))/(dx1*dx2*(dx1 + dx2))
        do i = 2, NoPmove-1
            dx1 = beta1(i)
            dx2 = beta2(i)
            ddn(i) = 2*(dx1*(x(i+1)*nx(i) + y(i+1)*ny(i)) - &
            (dx1+dx2)*(x(i)*nx(i) + y(i)*ny(i)) + &
            dx2*(x(i-1)*nx(i) + y(i-1)*ny(i)))/(dx1*dx2*(dx1 + dx2))
        end do
        dx1 = beta1(NoPmove)
        dx2 = beta2(NoPmove)
        ddn(NoPmove) = 2*(dx1*(x(1)*nx(NoPmove) + y(1)*ny(NoPmove)) - &
        (dx1+dx2)*(x(NoPmove)*nx(NoPmove) + y(NoPmove)*ny(NoPmove)) + &
        dx2*(x(NoPmove-1)*nx(NoPmove) + y(NoPmove-1)*ny(NoPmove)))/(dx1*dx2*(dx1 + dx2))

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
        do i = 2, NoPmove-1
            dx1 = beta1(i)
            dx2 = beta2(i)
            ddn(i) = 2*(dx1*(x(i+1)*nx(i) + y(i+1)*ny(i)) - &
            (dx1+dx2)*(x(i)*nx(i) + y(i)*ny(i)) + &
            dx2*(x(i-1)*nx(i) + y(i-1)*ny(i)))/(dx1*dx2*(dx1 + dx2))
        end do
        ddn(NoPmove) = 0

    end subroutine CalcSecondDerivativeopen

    subroutine calcBeta(x, y, NoPmove, betamin, betamean, beta1, beta2)

        ! Variables
        implicit none
        double precision :: betamin, betamean
        double precision, dimension(NoPmove) :: x, y, beta1, beta2
        integer :: NoPmove

        ! Body of calcBeta
        if (IV%shapeenclosed == .true.) then
            call calcBetaclosed(x, y, NoPmove, betamin, betamean, beta1, beta2)
        else
            call calcBetaopen(x, y, NoPmove, betamin, betamean, beta1, beta2)
        end if

    end subroutine calcBeta

    subroutine calcBetaclosed(x, y, NoPmove, betamin, betamean, beta1, beta2)

        ! Variables
        implicit none
        double precision :: betamin, betamean, magnitude
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
beta1(1) = 0.1
beta2(1) = 0.1
        betamin = minval((/minval(beta1), minval(beta2)/))
        betamean = (sum(beta1)+sum(beta2))/(size(beta1, dim = 1)+size(beta2, dim = 1))

    end subroutine calcBetaclosed

    subroutine calcBetaopen(x, y, NoPmove, betamin, betamean, beta1, beta2)

        ! Variables
        implicit none
        double precision :: betamin, betamean, magnitude
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

    subroutine calcSmoothing(x, y, NoPmove, smoothfactor, smoothing)

        ! Variables
        implicit none
        integer :: NoPmove
        double precision, dimension(IV%NoCN) :: smoothfactor
        double precision, dimension(NoPmove) :: x, y, smoothing

        ! Body of calcAlpha
        if (IV%shapeenclosed == .true.) then
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
        integer, dimension(:), allocatable :: CN_sort, sortind

        ! Order CN low to high
        allocate(CN_sort(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sortind(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        CN_sort = CN_indordered
        sortind = (/ (j, j=1,IV%NoCN) /)
        call QSortInt(CN_sort, size(CN_sort, dim = 1), 'y', sortind)

        ! Body of smoothingclosed
        allocate(dist(NoPmove+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        dist(1) = 0.0
        do i = 2, NoPmove
            dist(i) = dist(i-1) + sqrt((x(i)-x(i-1))**2 + (y(i)-y(i-1))**2)
        end do
        dist(NoPmove+1) = dist(NoPmove) + sqrt((x(1)-x(NoPmove))**2 + (y(1)-y(NoPmove))**2)

        smoothingclosed = 0.0
        do i = 1, IV%NoCN-1
            delta = (dist(CN_sort(i+1))-dist(CN_sort(i)))
            do j = CN_sort(i)+1, CN_sort(i+1)
                smoothingclosed(j) = smoothfactor(sortind(i))*(1 - (dist(j) - dist(CN_sort(i)))/delta) + smoothfactor(sortind(i+1))*((dist(j) - dist(CN_sort(i)))/delta)
            end do
        end do

        delta = (dist(NoPmove+1)-dist(CN_sort(IV%NoCN)) + dist(CN_sort(1)))
        do j = CN_sort(IV%NoCN)+1, NoPmove
            smoothingclosed(j) = smoothfactor(sortind(IV%NoCN))*(1 - (dist(j) - dist(CN_sort(IV%NoCN)))/delta) + smoothfactor(sortind(1))*((dist(j) - dist(CN_sort(IV%NoCN)))/delta)
        end do
        do j = 1, CN_sort(1)
            smoothingclosed(j) = smoothfactor(sortind(IV%NoCN))*(1 - (dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta) + smoothfactor(sortind(1))*((dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta)
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
        integer, dimension(:), allocatable :: CN_sort, sortind

        ! Order CN low to high
        allocate(CN_sort(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sortind(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        CN_sort = CN_indordered
        sortind = (/ (j, j=1,IV%NoCN) /)
        call QSortInt(CN_sort, size(CN_sort, dim = 1), 'y', sortind)

        ! Body of smoothingopen
        allocate(dist(NoPmove+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        do i = 2, NoPmove
            dist(i) = dist(i-1) + sqrt((x(i)-x(i-1))**2 + (y(i)-y(i-1))**2)
        end do

        smoothingopen = 0.0
        do i = 1, IV%NoCN-1
            delta = (dist(CN_sort(i+1))-dist(CN_sort(i)))
            do j = CN_sort(i)+1, CN_sort(i+1)
                smoothingopen(j) = smoothfactor(sortind(i))*(1 - (dist(j) - dist(CN_sort(i)))/delta) + smoothfactor(sortind(i+1))*((dist(j) - dist(CN_sort(i)))/delta)
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

    subroutine quasi1DLinear(x, y, CNxfinal, CNyfinal, NoPmove)

        ! Variables
        implicit none
        integer :: NoPmove
        double precision, dimension(:), allocatable :: x, y, CNxfinal, CNyfinal

        ! Body of quasi1DLinear
        if (IV%shapeenclosed == .true.) then
            call quasi1DLinearClosed(x, y, CNxfinal, CNyfinal, NoPmove)
        else
            call quasi1DLinearOpen(x, y, CNxfinal, CNyfinal, NoPmove)
        end if

    end subroutine quasi1DLinear

    subroutine quasi1DLinearOpen(x, y, CNxfinal, CNyfinal, NoPmove)

        ! Variables
        implicit none
        integer :: i, j, NoPmove
        double precision :: delta
        double precision, dimension(IV%NoCN,2) :: dCN
        double precision, dimension(:), allocatable :: x, y, CNxfinal, CNyfinal, dist
        integer, dimension(:), allocatable :: CN_sort, sortind

        ! Order CN low to high
        allocate(CN_sort(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sortind(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        CN_sort = CN_indordered
        sortind = (/ (j, j=1,IV%NoCN) /)
        call QSortInt(CN_sort, size(CN_sort, dim = 1), 'y', sortind)

        ! Body of quasi1DLinear
        ! Calculate required motion of CN
        do i = 1, IV%NoCN
            dCN(i,1) = CNxfinal(i) - x(CN_sort(i))
            dCN(i,2) = CNyfinal(i) - y(CN_sort(i))
        end do

        ! Calculate distance along boundary
        allocate(dist(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        dist = 0.0
        do i = 2, NoPmove
            dist(i) = dist(i-1) + sqrt((x(i)-x(i-1))**2 + (y(i)-y(i-1))**2)
        end do

        ! Move boundary points quasi-1D-linear
        do i = 1, IV%NoCN - 1
            delta = dist(CN_sort(i+1))-dist(CN_sort(i))
            do j = CN_sort(i)+1, CN_sort(i+1)
                y(j) = y(j) + dCN(i,2)*(1 - (dist(j) - dist(CN_sort(i)))/delta) + dCN(i+1,2)*((dist(j) - dist(CN_sort(i)))/delta)
                x(j) = x(j) + dCN(i,1)*(1 - (dist(j) - dist(CN_sort(i)))/delta) + dCN(i+1,1)*((dist(j) - dist(CN_sort(i)))/delta)
            end do
        end do

    end subroutine quasi1DLinearOpen

    subroutine quasi1DLinearClosed(x, y, CNxfinal, CNyfinal, NoPmove)

        ! Variables
        implicit none
        integer :: i, j, NoPmove
        double precision :: delta
        double precision, dimension(IV%NoCN,2) :: dCN
        double precision, dimension(:), allocatable :: x, y, CNxfinal, CNyfinal, dist
        integer, dimension(:), allocatable :: CN_sort, sortind

        ! Order CN low to high
        allocate(CN_sort(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sortind(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        CN_sort = CN_indordered
        sortind = (/ (j, j=1,IV%NoCN) /)
        call QSortInt(CN_sort, size(CN_sort, dim = 1), 'y', sortind)

        ! Body of quasi1DLinear
        ! Calculate required motion of CN
        do i = 1, IV%NoCN
            dCN(i,1) = CNxfinal(sortind(i)) - x(CN_sort(i))
            dCN(i,2) = CNyfinal(sortind(i)) - y(CN_sort(i))
        end do

        ! Calculate distance along boundary
        allocate(dist(NoPmove+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        dist = 0.0
        do i = 2, NoPmove
            dist(i) = dist(i-1) + sqrt((x(i)-x(i-1))**2 + (y(i)-y(i-1))**2)
        end do
        dist(NoPmove+1) = dist(NoPmove) + sqrt((x(1)-x(NoPmove))**2 + (y(1)-y(NoPmove))**2)

        ! Move boundary points quasi-1D-linear
        do i = 1, IV%NoCN - 1
            delta = dist(CN_sort(i+1))-dist(CN_sort(i))
            do j = CN_sort(i)+1, CN_sort(i+1)
                y(j) = y(j) + dCN(i,2)*(1 - (dist(j) - dist(CN_sort(i)))/delta) + dCN(i+1,2)*((dist(j) - dist(CN_sort(i)))/delta)
                x(j) = x(j) + dCN(i,1)*(1 - (dist(j) - dist(CN_sort(i)))/delta) + dCN(i+1,1)*((dist(j) - dist(CN_sort(i)))/delta)
            end do
        end do

        delta = (dist(NoPmove+1) + dist(CN_sort(1)) - dist(CN_sort(IV%NoCN)))
        do j = CN_sort(IV%NoCN)+1, NoPmove
            y(j) = y(j) + dCN(IV%NoCN,2)*(1 - (dist(j) - dist(CN_sort(IV%NoCN)))/delta) + dCN(1,2)*((dist(j) - dist(CN_sort(IV%NoCN)))/delta)
            x(j) = x(j) + dCN(IV%NoCN,1)*(1 - (dist(j) - dist(CN_sort(IV%NoCN)))/delta) + dCN(1,1)*((dist(j) - dist(CN_sort(IV%NoCN)))/delta)
        end do
        do j = 1, CN_sort(1)
            y(j) = y(j) + dCN(IV%NoCN,2)*(1 - (dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta) + dCN(1,2)*((dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta)
            x(j) = x(j) + dCN(IV%NoCN,1)*(1 - (dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta) + dCN(1,1)*((dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta)
        end do

    end subroutine quasi1DLinearClosed

    subroutine CalcFirstDerivative(dn, NoPmove, x, y, nx, ny, sx, sy)

        ! Variables
        implicit none
        integer :: ip, NoPmove, i
        double precision :: magnitude, dx1, dx2
        double precision, dimension(:), allocatable :: x, y, dn, sx, sy, nx, ny

        if (IV%shapeenclosed == .true.) then
            call CalcFirstDerivativeclosed(dn, NoPmove, x, y, nx, ny, sx, sy)
        else
            !call CalcFirstDerivativeopen(dn, NoPmove, x, y, nx, ny, sx, sy)
        end if

    end subroutine CalcFirstDerivative

    subroutine CalcFirstDerivativeclosed(dn, NoPmove, x, y, nx, ny, sx, sy)

        ! Variables
        implicit none
        integer :: ip, NoPmove, i
        double precision :: magnitude, dx1, dx2
        double precision, dimension(:), allocatable :: x, y, dn, sx, sy, nx, ny

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
        dn(1) = (dx1**2*(x(2)*nx(1) + y(2)*ny(1)) + &
        (dx2**2-dx1**2)*(x(1)*nx(1) + y(1)*ny(1)) - &
        dx2**2*(x(NoPmove)*nx(1) + y(NoPmove)*ny(1)))/(dx1*dx2*(dx1 + dx2))
        do i = 2, NoPmove-1
            dx1 = (x(i) - x(i-1))*sx(i) + (y(i) - y(i-1))*sy(i)
            dx2 = (x(i+1) - x(i))*sx(i) + (y(i+1) - y(i))*sy(i)
            dn(i) = (dx1**2*(x(i+1)*nx(i) + y(i+1)*ny(i)) + &
            (dx2**2-dx1**2)*(x(i)*nx(i) + y(i)*ny(i)) - &
            dx2**2*(x(i-1)*nx(i) + y(i-1)*ny(i)))/(dx1*dx2*(dx1 + dx2))
        end do
        dx1 = (x(NoPmove) - x(NoPmove-1))*sx(NoPmove) + (y(NoPmove) - y(NoPmove-1))*sy(NoPmove)
        dx2 = (x(1) - x(NoPmove))*sx(NoPmove) + (y(1) - y(NoPmove))*sy(NoPmove)
        dn(NoPmove) = (dx1**2*(x(1)*nx(NoPmove) + y(1)*ny(NoPmove)) + &
        (dx2**2-dx1**2)*(x(NoPmove)*nx(NoPmove) + y(NoPmove)*ny(NoPmove)) - &
        dx2**2*(x(NoPmove-1)*nx(NoPmove) + y(NoPmove-1)*ny(NoPmove)))/(dx1*dx2*(dx1 + dx2))

    end subroutine CalcFirstDerivativeclosed

    subroutine normXY(x, y, normx, normy, CNxfinal, CNyfinal, NoPmove)

        ! Variables
        implicit none
        integer :: NoPmove
        double precision, dimension(:), allocatable :: x, y, CNxfinal, CNyfinal, normx, normy

        ! Body of quasi1DLinear
        if (IV%shapeenclosed == .true.) then
            call normXYclosed(x, y, normx, normy, CNxfinal, CNyfinal, NoPmove)
        else
            call normXYopen(x, y, normx, normy, CNxfinal, CNyfinal, NoPmove)
        end if

    end subroutine normXY

    subroutine normXYopen(x, y, normx, normy, CNxfinal, CNyfinal, NoPmove)

        ! Variables
        implicit none
        integer :: i, j, NoPmove
        double precision :: delta, norm
        double precision, dimension(IV%NoCN,2) :: dCN
        double precision, dimension(:), allocatable :: x, y, CNxfinal, CNyfinal, dist, normx, normy
        integer, dimension(:), allocatable :: CN_sort, sortind

        ! Order CN low to high
        allocate(CN_sort(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sortind(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        CN_sort = CN_indordered
        sortind = (/ (j, j=1,IV%NoCN) /)
        call QSortInt(CN_sort, size(CN_sort, dim = 1), 'y', sortind)

        ! Body of quasi1DLinear
        ! Calculate required motion of CN
        do i = 1, IV%NoCN
            dCN(i,1) = abs(CNxfinal(sortind(i)) - x(CN_sort(i)))
            dCN(i,2) = abs(CNyfinal(sortind(i)) - y(CN_sort(i)))
        end do

        ! Calculate distance along boundary
        allocate(dist(NoPmove+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        dist = 0.0
        do i = 2, NoPmove
            dist(i) = dist(i-1) + sqrt((x(i)-x(i-1))**2 + (y(i)-y(i-1))**2)
        end do
        dist(NoPmove+1) = dist(NoPmove) + sqrt((x(1)-x(NoPmove))**2 + (y(1)-y(NoPmove))**2)

        ! Calculate y and x contribution based on distance and magnitude
        normy = 0.5
        normx = 0.5
        do i = 1, IV%NoCN - 1
            delta = dist(CN_sort(i+1))-dist(CN_sort(i))
            norm = dCN(i,2) + dCN(i,1) + dCN(i+1,2) + dCN(i+1,1)
            if (norm /= 0) then
                do j = CN_sort(i), CN_sort(i+1)
                    normy(j) = dCN(i,2)/norm*(1 - (dist(j) - dist(CN_sort(i)))/delta) + dCN(i+1,2)/norm*((dist(j) - dist(CN_sort(i)))/delta)
                    normx(j) = dCN(i,1)/norm*(1 - (dist(j) - dist(CN_sort(i)))/delta) + dCN(i+1,1)/norm*((dist(j) - dist(CN_sort(i)))/delta)
                end do
            end if
            if (normx(CN_sort(i)) == 0 .and. normy(CN_sort(i)) == 0) then
                normx(CN_sort(i)) = 0.5
                normy(CN_sort(i)) = 0.5
            end if
        end do
        if (normx(CN_sort(IV%NoCN)) == 0 .and. normy(CN_sort(IV%NoCN)) == 0) then
            normx(CN_sort(IV%NoCN)) = 0.5
            normy(CN_sort(IV%NoCN)) = 0.5
        end if

        ! Ensure Sum of normy + normx = 1
        normy = normy/(normy + normx)
        normx = normx/(normy + normx)

    end subroutine normXYopen

    subroutine normXYclosed(x, y, normx, normy, CNxfinal, CNyfinal, NoPmove)

        ! Variables
        implicit none
        integer :: i, j, NoPmove
        double precision :: delta, norm
        double precision, dimension(IV%NoCN,2) :: dCN
        double precision, dimension(:), allocatable :: x, y, CNxfinal, CNyfinal, dist, normx, normy
        integer, dimension(:), allocatable :: CN_sort, sortind

        ! Order CN low to high
        allocate(CN_sort(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(sortind(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        CN_sort = CN_indordered
        sortind = (/ (j, j=1,IV%NoCN) /)
        call QSortInt(CN_sort, size(CN_sort, dim = 1), 'y', sortind)

        ! Body of quasi1DLinear
        ! Calculate required motion of CN
        do i = 1, IV%NoCN
            dCN(i,1) = CNxfinal(sortind(i)) - x(CN_sort(i))
            dCN(i,2) = CNyfinal(sortind(i)) - y(CN_sort(i))
        end do

        ! Calculate distance along boundary
        allocate(dist(NoPmove+1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        dist = 0.0
        do i = 2, NoPmove
            dist(i) = dist(i-1) + sqrt((x(i)-x(i-1))**2 + (y(i)-y(i-1))**2)
        end do
        dist(NoPmove+1) = dist(NoPmove) + sqrt((x(1)-x(NoPmove))**2 + (y(1)-y(NoPmove))**2)

        ! Calculate y and x contribution based on distance and magnitude
        normy(CN_sort(1)) = 0.5
        normx(CN_sort(1)) = 0.5
        do i = 1, IV%NoCN - 1
            delta = dist(CN_sort(i+1))-dist(CN_sort(i))
            norm = dCN(i,2) + dCN(i,1) + dCN(i+1,2) + dCN(i+1,1)
            if (norm == 0) then
                normy((CN_sort(i)+1):CN_sort(i+1)) = 0.5
                normx((CN_sort(i)+1):CN_sort(i+1)) = 0.5
            else
                do j = CN_sort(i), CN_sort(i+1)
                    normy(j) = dCN(i,2)/norm*(1 - (dist(j) - dist(CN_sort(i)))/delta) + dCN(i+1,2)/norm*((dist(j) - dist(CN_sort(i)))/delta)
                    normx(j) = dCN(i,1)/norm*(1 - (dist(j) - dist(CN_sort(i)))/delta) + dCN(i+1,1)/norm*((dist(j) - dist(CN_sort(i)))/delta)
                end do
            end if
            if (normx(CN_sort(i)) == 0 .and. normy(CN_sort(i)) == 0) then
                normx(CN_sort(i)) = 0.5
                normy(CN_sort(i)) = 0.5
            end if
        end do

        delta = (dist(NoPmove+1) + dist(CN_sort(1)) - dist(CN_sort(IV%NoCN)))
        norm = dCN(IV%NoCN,2) + dCN(IV%NoCN,1) + dCN(1,2) + dCN(1,1)
        if (norm == 0) then
            normy((CN_sort(IV%NoCN)+1):NoPmove) = 0.5
            normx((CN_sort(IV%NoCN)+1):NoPmove) = 0.5
            normy(1:CN_sort(1)) = 0.5
            normx(1:CN_sort(1)) = 0.5
        else
            do j = CN_sort(IV%NoCN)+1, NoPmove
                normy(j) = dCN(IV%NoCN,2)/norm*(1 - (dist(j) - dist(CN_sort(IV%NoCN)))/delta) + dCN(1,2)/norm*((dist(j) - dist(CN_sort(IV%NoCN)))/delta)
                normx(j) = dCN(IV%NoCN,1)/norm*(1 - (dist(j) - dist(CN_sort(IV%NoCN)))/delta) + dCN(1,1)/norm*((dist(j) - dist(CN_sort(IV%NoCN)))/delta)
            end do
            do j = 1, CN_sort(1)
                normy(j) = dCN(IV%NoCN,2)/norm*(1 - (dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta) + dCN(1,2)/norm*((dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta)
                normx(j) = dCN(IV%NoCN,1)/norm*(1 - (dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta) + dCN(1,1)/norm*((dist(j) + dist(NoPmove+1) - dist(CN_sort(IV%NoCN)))/delta)
            end do
        end if
        if (normx(CN_sort(IV%NoCN)) == 0 .and. normy(CN_sort(IV%NoCN)) == 0) then
            normx(CN_sort(IV%NoCN)) = 0.5
            normy(CN_sort(IV%NoCN)) = 0.5
        end if

        ! Ensure Sum of normy + normx = 1
        normy = normy/(normy + normx)
        normx = normx/(normy + normx)

    end subroutine normXYclosed

    subroutine CheckBoundIntersections(x, y, intersect)

        ! Variables
        implicit none
        integer :: i, j, nibp, intersect
        double precision :: A, B, C, Atarget, Btarget, Ctarget, detAB, detCB, detAC, xtemp, ytemp
        double precision, dimension(:), allocatable :: x, y

        ! Body of getPressureDistribution
        nibp = size(InnerBound)

        ! Identify intersections between target and current boundary
        do j = 1, nibp-1
            do i = j+2, nibp-1

                ! Calculate intersection point
                A = y(j+1) - y(j)
                B = x(j) - x(j+1)
                C = A*x(j) + B*y(j)

                Atarget = y(i+1) - y(i)
                Btarget = x(i) - x(i+1)
                Ctarget = Atarget*x(i) + Btarget*y(i)

                detAB = A*Btarget - Atarget*B
                detCB = C*Btarget - Ctarget*B
                detAC = A*Ctarget - Atarget*C
                if (abs(detAB) < abs(10e-12)) then
                    ! Parallel
                    xtemp = -1000
                    ytemp = -1000
                else
                    ! Non-parallel
                    xtemp = detCB/detAB
                    ytemp = detAC/detAB
                end if

                ! Check, if intersection point lays on the line segment considered
                if ((xtemp - min(x(j),x(j+1))) > -10e-12 .and. (xtemp - max(x(j),x(j+1))) < 10e-12 .and. (ytemp - min(y(j),y(j+1))) > -10e-12  .and. (ytemp - max(y(j),y(j+1))) < 10e-12)  then
                    if ((xtemp - min(x(i),x(i+1))) > -10e-12 .and. (xtemp - max(x(i),x(i+1))) < 10e-12  .and. (ytemp - min(y(i),y(i+1))) > -10e-12  .and. (ytemp - max(y(i),y(i+1))) < 10e-12 )  then
                        print *, 'Intersection along Boundary identified. Movement will be halfed.'
                        intersect = 2
                        EXIT
                    end if
                end if
            end do
            if (intersect == 2) then
                EXIT
            end if
        end do

    end subroutine CheckBoundIntersections

    end module Smoothing
