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
            if (IV%MeshMovement == 1) then
                call SmoothingLinear(CN_CoordinatesMatrix)
            else if (IV%MeshMovement == 2) then
                call SmoothingFDGD(CN_CoordinatesMatrix, ii)
            else if (IV%MeshMovement == 3) then
                call SubRBF(CN_CoordinatesMatrix)
            else if (IV%MeshMovement == 4) then            
                call SubFDGD(CN_CoordinatesMatrix, ii, iii)
            end if
        elseif (IV%NoDim == 3) then
            if (IV%MeshMovement == 1) then
                call SmoothingLinear(CN_CoordinatesMatrix)
            else if (IV%MeshMovement == 2) then
                call SmoothingFDGD(CN_CoordinatesMatrix, ii)
            else if (IV%MeshMovement == 3) then
                call SubRBF(CN_CoordinatesMatrix)
            else if (IV%MeshMovement == 4) then            
                call SubFDGD_3D(CN_CoordinatesMatrix, ii, iii)
            end if
        end if
    
    end subroutine SubMovemesh

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