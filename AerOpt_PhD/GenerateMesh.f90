    module GenerateMesh
    use Toolbox
    use ReadData
    use InputData
    use CreateSnapshots

    contains

    subroutine SubGenerateMesh(NestDisp)

    ! Variables
    implicit none
    integer :: lb
    real, dimension(IV%NoCP*IV%NoDim) :: NestDisp
    real, dimension(:), allocatable :: dCP2N, CP_ind, size_ib
    integer, dimension(:,:), allocatable :: IB
    real :: c, w, dis, x1, y1, x2, y2, x3, y3, xp, yp

    ! Body of GenerateInitialMeshes
    allocate(dCP2N(RD%nbf),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
    allocate(IB(RD%nbf, IV%NoCP),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
    allocate(CP_ind(IV%NoCP),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
    allocate(size_ib(IV%NoCP),stat=allocateStatus)
    if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "

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
            dis = ((RD%coord_temp(CP_ind(i),1) - RD%coord_temp(IB(j,i),1)))*1.5   ! Distance of Control Node to a Node in the Influence Box         
            w = exp(-(dis**2)/(c**2))   ! Gaussian
            RD%coord_temp(IB(j,i),:) = (/(RD%coord_temp(IB(j,i),1) + w*NestDisp(i)),(RD%coord_temp(IB(j,i),2) + w*NestDisp(IV%NoCP+i))/)   ! New coordinates
        end do
    end do

    !!!!!! PLOT initial Nests --> MATLAB output file

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

    end subroutine SubGenerateMesh

    subroutine LocateNode()

    ! Variables
    double precision, dimension(3,3) :: A
    real :: D

    ! Body of LocateNode
    !A1 = /coord(



    end subroutine LocateNode

    subroutine IdentifyBoundaryFlags()

    ! Variables
    implicit none
    real, dimension(2) :: point

    ! Body of IdentifyBoundaryFlags
    do k = 1, RD%nbf
        point(1) = (RD%coord(RD%boundf(k,1),1)*15 + RD%coord(RD%boundf(k,2),1)*15)/2.0
        point(2) = (RD%coord(RD%boundf(k,1),2)*15 + RD%coord(RD%boundf(k,2),2)*15)/2.0
        if (point(1) < 20 .and. point(1) > 0 .and. point(2) < 3 .and. point(2) > (-2)) then
            RD%boundf(k,3) = 6    ! Adiabatic Viscous Wall
        elseif (point(1) == 0 .and. point(2) < 1 .and. point(2) > 0) then
            RD%boundf(k,3) = 8    ! Engine Inlet
        else    
            RD%boundf(k,3) = 3    ! Far Field
        end if       
    end do

    end subroutine IdentifyBoundaryFlags

end module GenerateMesh