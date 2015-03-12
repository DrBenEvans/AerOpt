    module GenerateMesh
    use Toolbox
    use ReadData
    use InputData
    use CreateSnapshots
    use FDGD

    contains
    
    subroutine SubMovemesh(CN_CoordinatesVector)
    
        ! Variables
        implicit none
        integer :: size
        double precision, dimension(maxDoF) :: CN_CoordinatesVector

        ! Body of SubMovemesh
        if (IV%MeshMovement == 1) then
            call SubFDGD(CN_CoordinatesVector)
        else if (IV%MeshMovement == 2) then
            call SubRBF(CN_CoordinatesVector)
        end if
    
    end subroutine SubMovemesh

    subroutine SubRBF(NestDisp)

    ! Variables
    implicit none
    integer :: lb, i, k, j
    double precision, dimension(IV%NoCP*IV%NoDim) :: NestDisp
    double precision, dimension(:), allocatable :: dCP2N, CP_ind, size_ib
    integer, dimension(:,:), allocatable :: IB
    double precision :: c, w, dis, x1, y1, x2, y2, x3, y3, xp, yp

    ! Body of RBF
    allocate(dCP2N(RD%nbf),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in RBF "
    allocate(IB(RD%nbf, IV%NoCP),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in RBF "
    allocate(CP_ind(IV%NoCP),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in RBF "
    allocate(size_ib(IV%NoCP),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in RBF "

    ! Find Real Control Nodes based on Coordinates - !! on long terms redundant only required to run once           
    ! Also Identify boundary nodes within the influence box of each Control Point
    !!!! NOTE: Hardwired for a 2D problem!
    do i = 1, IV%NoCP
        k = 0
        do j = 1, RD%nbf
            dCP2N(j) = DistP2P(IV%NoDim, RD%Coord_CP(i,1), RD%coord_temp(j,1), RD%Coord_CP(i, 2), RD%coord_temp(j,2))  ! Calculate Distances

            if (Rectcheck(RD%Rect(i,:), RD%coord_temp(j,:)) == 1) then ! Check if boundary node is in influence box - Note: first values in coordinates matrix are per default the boundary nodes
                k = k + 1
                if (i /= 6) then
                    IB(k,i) = j
                elseif (k > 15) then ! Bloodhound engine inlet specific. Dependant on the definition of the rectangles!
                    IB((k-15),i) = j
                end if                           
            end if 
        end do
        CP_ind(i) = minloc(dCP2N,dim=1) ! Returns the index(Position of Node in Coord Matrix) of the minimum Value          
        if (i /= 6) then
            size_ib(i) = k
        else
            size_ib(i) = k - 15
        end if

    end do       

    ! Relocating boundary nodes based on their distance to Control node and the Control Nodes Displacements (NestDisp)
    ! Method: Gaussian RBF function
    ! c = 1.9
    do i = 1, IV%NoCp
        c = abs((RD%Rect(i,3) - RD%Rect(i,1)))/2.0
        do j = 1, size_ib(i)
            dis = RD%coord_temp(CP_ind(i),1) - RD%coord_temp(IB(j,i),1)   ! Distance of Control Node to a Node in the Influence Box         
            w = exp(-(dis**2)/(c**2))   ! Gaussian
            RD%coord_temp(IB(j,i),:) = (/(RD%coord_temp(IB(j,i),1) + w*NestDisp(i)),(RD%coord_temp(IB(j,i),2) + w*NestDisp(IV%NoCP+i))/)   ! New coordinates
        end do
    end do

    if (IV%ObjectiveFunction == 2) then
        ! Relocation of the internal nodes(non-boundary)
        ! Method: Maintain the Area coefficient
        do i = 1, RD%np - RD%nbf
            x1 = RD%coord_temp(RD%connecc(RD%coarse(i,1),1),1)
            y1 = RD%coord_temp(RD%connecc(RD%coarse(i,1),1),2)
            x2 = RD%coord_temp(RD%connecc(RD%coarse(i,1),2),1)
            y2 = RD%coord_temp(RD%connecc(RD%coarse(i,1),2),2)
            x3 = RD%coord_temp(RD%connecc(RD%coarse(i,1),3),1)
            y3 = RD%coord_temp(RD%connecc(RD%coarse(i,1),3),2)

            xp = x1*RD%coarse(i,2) + x2*RD%coarse(i,3) + x3*RD%coarse(i,4)
            yp = y1*RD%coarse(i,2) + y2*RD%coarse(i,3) + y3*RD%coarse(i,4)
            RD%coord_temp((i+RD%nbf),:) = (/xp, yp/)
        end do
    else
        ! Move Domain Nodes
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        allocate(DelaunayElem(size(DelaunayElemDomain, dim = 1),size(DelaunayElemDomain, dim = 2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD "
        DelaunayElem = DelaunayElemDomain
        call RelocateMeshPoints(AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
        call CheckforIntersections()
        deallocate(DelaunayCoord)
        deallocate(DelaunayElem)
    end if

    end subroutine SubRBF

end module GenerateMesh