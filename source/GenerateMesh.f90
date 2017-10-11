    module GenerateMesh
    use Toolbox
    use ReadData
    use InputData
    use CreateSnapshots
    use FDGD
    use Smoothing

    contains
    
    subroutine SubMovemesh(CN_CoordinatesMatrix)
    
        ! Variables
        implicit none
        integer :: size, ii, iii
        double precision, dimension(maxDoF) :: CN_CoordinatesMatrix

        ! Body of SubMovemesh        
        ii = 1
        iii = 1
        if (IV%NoDim == 2) then
! Add new features
            select case (IV%MeshMovement)
            case (1)
                call SmoothingLinear(CN_CoordinatesMatrix, ii, iii)
            case (2)
                call SmoothingFDGD(CN_CoordinatesMatrix, ii, iii)
            case (3)
                call SubRBF(CN_CoordinatesMatrix)
            case (4)
                call SubFDGD(CN_CoordinatesMatrix, ii, iii)
            case (5)
                call SubQuasi1DLinear(CN_CoordinatesMatrix, ii, iii)
            case default
                call SmoothingLinear(CN_CoordinatesMatrix, ii, iii)
            end select
        elseif (IV%NoDim == 3) then
            select case (IV%MeshMovement)
            case (1)
                call SmoothingLinear(CN_CoordinatesMatrix, ii, iii)
            case (2)
                call SmoothingFDGD(CN_CoordinatesMatrix, ii, iii)
            case (3)
                call SubRBF(CN_CoordinatesMatrix)
            case (4)
                call SubFDGD_3D(CN_CoordinatesMatrix, ii, iii)
            case (5)
                call SubQuasi1DLinear(CN_CoordinatesMatrix, ii, iii)
            case default
                call SmoothingLinear(CN_CoordinatesMatrix, ii, iii)
            end select
        end if
    
    end subroutine SubMovemesh
    
    recursive subroutine SubQuasi1DLinear(CNDisp, NoMove, counter)
    
        ! Variables
        implicit none
        integer :: NoMove, counter, intersect, NoPmove, i
        double precision, dimension(:), allocatable :: x, y, xbefore, ybefore, CNxfinal, CNyfinal
        double precision, dimension(maxDoF) :: CNDisp
        double precision, dimension(2) :: a
        
        ! Body of SubQuasi1DLinear
        NoPmove = size(InnerBound, dim = 1)
        allocate(x(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(y(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(xbefore(NoPmove),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in Smoothing "
        allocate(ybefore(NoPmove),stat=allocateStatus)
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
        
        ! Final CN positions
        do i = 1, IV%NoCN
            CNxfinal(i) = x(CN_indordered(i)) + CNDisp(i)/NoMove
            CNyfinal(i) = y(CN_indordered(i)) + CNDisp(IV%NoCN+i)/NoMove
            if (IV%CNconnecttrans(i) /= 0) then
                CNxfinal(i) = x(CN_indordered(i)) + CNDisp(IV%CNconnecttrans(i))/NoMove
                CNyfinal(i) = y(CN_indordered(i)) + CNDisp(IV%NoCN+IV%CNconnecttrans(i))/NoMove
            end if
        end do
        
        ! Perform linear motion
        call quasi1DLinear(x, y, CNxfinal, CNyfinal, NoPmove)
        
        ! Hand over Coordinates to temporary coord
        RD%coord_temp(orderedBoundaryIndex,1) = x
        RD%coord_temp(orderedBoundaryIndex,2) = y
        
        ! Rotate
        do i = 1, IV%NoCN
            if (IV%angle(i) /= 0) then
                call AngleofAttack(CNDisp(4*i)/NoMove,i)
            end if
        end do
        
        ! Check for valid background mesh
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        call CheckforIntersections(DelaunayCoordDomain, DelaunayElemDomain, intersect)

        ! Move Domain Nodes
        if (intersect == 1) then
            if (counter < NoMove) then
                counter = counter + 1
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))              
                call PreMeshingMid()
                call SubQuasi1DLinear(CNDisp, NoMove, counter)               
            else
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                if (NoMove > 1) then
                    call PreMeshingEnd()
                end if
            end if
        else
            RD%coord_temp(orderedBoundaryIndex,1) = xbefore
            RD%coord_temp(orderedBoundaryIndex,2) = ybefore
            counter = 2*counter-1
            NoMove = NoMove*2
            call SubQuasi1DLinear(CNDisp, NoMove, counter)
        end if
    
    end subroutine SubQuasi1DLinear

    subroutine SubRBF(NestDisp)

    ! Variables
    implicit none
    integer :: lb, i, k, j
    double precision, dimension(maxDoF) :: NestDisp
    double precision, dimension(:), allocatable :: dCP2N, CP_ind, size_ib
    integer, dimension(:,:), allocatable :: IB
    double precision :: c, w, dis, x1, y1, x2, y2, x3, y3, xp, yp

    ! Body of RBF
    allocate(dCP2N(RD%nbf),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in RBF "
    allocate(IB(RD%nbf, IV%NoCN),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in RBF "
    allocate(CP_ind(IV%NoCN),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in RBF "
    allocate(size_ib(IV%NoCN),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in RBF "       

    ! Relocating boundary nodes based on their distance to Control node and the Control Nodes Displacements (NestDisp)
    ! Method: Gaussian RBF function
    !c = 1.9
    do i = 1, IV%NoCN
        !c = DistP2P(IV%NoDim, IV%xrange(i), IV%xrange(IV%NoCN+i), IV%yrange(i), IV%yrange(IV%NoCN+i), IV%zrange(i), IV%zrange(IV%NoCN+i))/2.0
        do j = 1, size_ib(i)
            dis = RD%coord_temp(CP_ind(i),1) - RD%coord_temp(IB(j,i),1)   ! Distance of Control Node to a Node in the Influence Box         
            w = exp(-(dis**2)/(c**2))   ! Gaussian
            RD%coord_temp(IB(j,i),:) = (/(RD%coord_temp(IB(j,i),1) + w*NestDisp(i)),(RD%coord_temp(IB(j,i),2) + w*NestDisp(IV%NoCN+i))/)   ! New coordinates
        end do
    end do

    !! Move Domain Nodes
    !! Check for valid background mesh
    !call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))          
    !call CheckforIntersections(DelaunayCoordDomain, DelaunayElemDomain, intersect)
    !    
    !! Move Domain Nodes
    !if (intersect == 1) then
    !    !if (NoMove > 1) then              
    !    !    NoMove = NoMove - 1
    !    !    call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
    !    !    call PreMeshing3()
    !    !    call SubFDGD(NestDisp, NoMove)               
    !    !else
    !        call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
    !        NestDisp = NestDisp*(1.0/NoMove)
    !    !end if
    !else
    !    NoMove = NoMove*2
    !    call SubFDGD(NestDisp, NoMove)
    !end if

    end subroutine SubRBF

end module GenerateMesh