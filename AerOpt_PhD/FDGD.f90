module FDGD
    
    use ReadData
    use InputData
    use Toolbox
    use CreateSnapshots
    
    double precision, dimension(:,:), allocatable :: AreaCoeffBound, AreaCoeffDomain
    double precision, dimension(:,:), allocatable :: DelaunayElem, DelaunayCoord
    double precision, dimension(:,:), allocatable :: DelaunayCoordBound
    integer :: NoElem
    
    contains
    
    subroutine SubFDGD(NestDisp)
    
        !**************************************************************************
        !***                    written by David Naumann                        ***
        !***                        Swansea University                          ***
        !**************************************************************************
        ! Objective:
        ! Part A
        ! 1. Create a Delaunay Triangulation with CN, interface points & domain boundary points.
        ! 2. Identify correct Delaunay Element and calculate Area Coefficients 
        !    of boundary points, that need to be moved.
        ! Part B
        ! 1. Move the mesh according to the CN movement
        ! 2. Finally, Validity check of Delaunay Mesh:
        !    Check for intersections of triangles after the movement
        !
        !**************************************************************************

        implicit none
        double precision, dimension(IV%NoCP*IV%NoDim) :: NestDisp
        integer :: i
        
        ! Body of SubFDGD
        
        ! Move Mesh of Boundary Nodes
        allocate(DelaunayCoord(size(DelaunayCoordBound, dim = 1), size(DelaunayCoordBound, dim = 2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        DelaunayCoord = DelaunayCoordBound
        ! Read in relocation of CP
        do i = 1, IV%NoCP
            DelaunayCoord(i,1) = DelaunayCoord(i,1) + NestDisp(i)
            DelaunayCoord(i,2) = DelaunayCoord(i,2) + NestDisp(i+IV%NoCP)
        end do
        call getDelaunayElem()
        call RelocateMeshPoints(AreaCoeffBound, size(AreaCoeffBound, dim = 1))
        call CheckforIntersections()
        deallocate(DelaunayCoord)
        deallocate(DelaunayElem)
        
        ! Move Mesh of Domain Nodes
        call getDelaunayCoordDomain()
        call getDelaunayElem()
        call RelocateMeshPoints(AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
        call CheckforIntersections()
        deallocate(DelaunayCoord)
        deallocate(DelaunayElem)
    
    end subroutine SubFDGD
    
    subroutine PreMeshing()
    
        ! Variables
        implicit none
        integer, dimension(:), allocatable ::  DomainIndex
        integer, dimension(:), allocatable :: MovingGeomIndex
    
        ! Body of PreMeshing
    
        ! Preparation for Boundary Movement
        call getDelaunayCoordBound(MovingGeomIndex)
        allocate(DelaunayCoord(size(DelaunayCoordBound, dim = 1), size(DelaunayCoordBound, dim = 2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        DelaunayCoord = DelaunayCoordBound
        call getDelaunayElem()
        call getAreaCoefficients(MovingGeomIndex, size(MovingGeomIndex), AreaCoeffBound)
        deallocate(DelaunayCoord)
        deallocate(DelaunayElem)
    
        ! Preparation for Domain Movement
        call getDelaunayCoordDomain()
        call getDelaunayElem()
        call getDomainIndex(DomainIndex)
        call getAreaCoefficients(DomainIndex, size(DomainIndex), AreaCoeffDomain)
        deallocate(DelaunayCoord)
        deallocate(DelaunayElem)
    
    end subroutine PreMeshing
    
    
    subroutine getDelaunayElem()
    ! Objective: Do the Triangulation and get the Delaunay Element Matrix
    
        ! Variables
        implicit none
        integer :: i, FileSize
        character(len=8) :: Filename = 'Delaunay'
    
        ! Body of getDelaunayElem
        open(1, file= Filename//'_nodes.txt', form='formatted', status = 'unknown')
        write(1,'(2f22.15)') transpose(DelaunayCoord)
        close(1)
    
        ! Performt Delaunay Triangulation via external executable
        allocate(character(len=100) :: strSystem,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        strSystem = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/DelaunayTriangulation < Delaunay/DelaunayInput.txt > /dev/null'
        call system(trim(strSystem))
        deallocate(strSystem)
        
        open(1, file= Filename//'_elements.txt', form='formatted',status='old')
        inquire(1, size = FileSize)
        NoElem = FileSize/32
        allocate(DelaunayElem(NoElem,IV%NoDim + 1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in BruteForce "
        do i = 1, NoElem
            read(1, *) DelaunayElem(i,:)
        end do        
        close(1)
    
    end subroutine getDelaunayElem
    
    
    subroutine getAreaCoefficients(RowIndex, NoP, AreaCoeff) 
    
        ! Variables
        implicit none
        double precision :: xa, ya, xb, yb, xc, yc, xp, yp, S1, S2, S3, S, e1, e2, e3
        integer :: NoNodesMove, i, j, NoP
        integer, dimension(NoP) :: RowIndex
        double precision, dimension(:,:), allocatable :: AreaCoeff
        
        ! Body of getAreaCoefficients        
           
        ! Identify correct Delaunay Triangle and store Area Coefficients
        allocate(AreaCoeff(NoP,5),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in BruteForce "
        do i = 1, NoP
            xp = RD%coord(RowIndex(i),1)
            yp = RD%coord(RowIndex(i),2)
            do j = 1, NoElem
                xa = DelaunayCoord(DelaunayElem(j,1),1)
                ya = DelaunayCoord(DelaunayElem(j,1),2)
                xb = DelaunayCoord(DelaunayElem(j,2),1)
                yb = DelaunayCoord(DelaunayElem(j,2),2)
                xc = DelaunayCoord(DelaunayElem(j,3),1)
                yc = DelaunayCoord(DelaunayElem(j,3),2)
                S1 = 0.5*((xp*yb + yp*xc + xb*yc) - (xp*yc + yp*xb + yb*xc))
                S2 = 0.5*((xa*yp + ya*xc + xp*yc) - (xa*yc + ya*xp + yp*xc))
                S3 = 0.5*((xa*yb + ya*xp + xb*yp) - (xa*yp + ya*xb + yb*xp))
                S = 0.5*((xa*yb + ya*xc + xb*yc) - (xa*yc + ya*xb + yb*xc))
                e1 = S1/S
                e2 = S2/s
                e3 = S3/S
                if (e1 >= 0 .and. e2 >= 0 .and. e3 >= 0) then
                    AreaCoeff(i,:) = (/ dble(RowIndex(i)), dble(j), e1, e2, e3/)
                    EXIT
                end if
            end do
        end do
        
! Redundant
        open(1, file= 'Area_Coefficients.txt',form='formatted',status='unknown')
        write(1,'(4f22.15)') transpose(AreaCoeff)
        close(1)
    
    end subroutine getAreaCoefficients
    
    
    subroutine CheckforIntersections()
    
        ! Variables
        implicit none
        double precision :: xa, ya, xb, yb, xc, yc, S
        integer :: j, NoElem
    
        ! Body of CheckforIntersections
        do j = 1, NoElem
            xa = DelaunayCoord(DelaunayElem(j,1),1)
            ya = DelaunayCoord(DelaunayElem(j,1),2)
            xb = DelaunayCoord(DelaunayElem(j,2),1)
            yb = DelaunayCoord(DelaunayElem(j,2),2)
            xc = DelaunayCoord(DelaunayElem(j,3),1)
            yc = DelaunayCoord(DelaunayElem(j,3),2)
            S = 0.5*((xa*yb + ya*xc + xb*yc) - (xa*yc + ya*xb + yb*xc))
            if (S < 0) then  ! If Area negative, the order of the points has changed (DelaunayElement order) and hence an intersection of elements took place
                print *, 'Intersection identified. Program paused.'
                pause
                strSystem = '2D_Delaunay.exe < DelaunayInput.txt'
                call system(trim(strSystem))
            end if
        end do
    
    end subroutine CheckforIntersections
    
    subroutine RelocateMeshPoints(AreaCoeff, NoP)
    
        ! Variables
        implicit none
        integer :: i, NoP
        double precision :: x1, y1, x2, y2, x3, y3, xp, yp, inp
        double precision, dimension(NoP,5) :: AreaCoeff
        
        ! Body of RelocateMeshPoints        
        do i = 1, NoP
            x1 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),1),1)
            y1 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),1),2)
            x2 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),2),1)
            y2 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),2),2)
            x3 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),3),1)
            y3 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),3),2)
            xp = x1*AreaCoeff(i,3) + x2*AreaCoeff(i,4) + x3*AreaCoeff(i,5)
            yp = y1*AreaCoeff(i,3) + y2*AreaCoeff(i,4) + y3*AreaCoeff(i,5)
            RD%coord_temp(AreaCoeff(i,1),:) = (/xp, yp/)
        end do
        
        !! Relocation of the internal nodes(non-boundary)
        !! Method: Maintain the Area coefficient
        !do i = 1, RD%np - RD%nbf
        !    x1 = RD%coord_temp(RD%connecc(RD%coarse(i,1),1),1)
        !    y1 = RD%coord_temp(RD%connecc(RD%coarse(i,1),1),2)
        !    x2 = RD%coord_temp(RD%connecc(RD%coarse(i,1),2),1)
        !    y2 = RD%coord_temp(RD%connecc(RD%coarse(i,1),2),2)
        !    x3 = RD%coord_temp(RD%connecc(RD%coarse(i,1),3),1)
        !    y3 = RD%coord_temp(RD%connecc(RD%coarse(i,1),3),2)
        !
        !    xp = x1*RD%coarse(i,2) + x2*RD%coarse(i,3) + x3*RD%coarse(i,4)
        !    yp = y1*RD%coarse(i,2) + y2*RD%coarse(i,3) + y3*RD%coarse(i,4)
        !    RD%coord_temp((i+RD%nbf),:) = (/xp, yp/)
        !end do
    
        
    end subroutine RelocateMeshPoints
    
    subroutine getDelaunayCoordDomain()
    
        ! Variables
        implicit none
    
        ! Body of getDelaunayCoordDomain
        integer, dimension(:), allocatable :: BoundIndex
    
        !Identify Boundary Indices
        call getBoundaryIndex(BoundIndex) 
        allocate(DelaunayCoord(size(BoundIndex),IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        DelaunayCoord = RD%Coord_temp(BoundIndex,:)
    
    end subroutine getDelaunayCoordDomain
    
    subroutine getDelaunayCoordBound(MovingGeomIndex)
    
        ! Variables
        implicit none
        integer :: i, j, k, l, nbp, ndbp, nibp, nobp, overlap, testx
        double precision :: maxx, maxy, minx, spacing
        double precision, dimension(:), allocatable :: dist, CP_ind
        integer, dimension(:), allocatable ::  BoundIndex, InnerBound, OuterBound, nodesvec, nodesvec2, NonMovingGeomIndex, distindex
        integer, dimension(:), allocatable :: MovingGeomIndex
        real, PARAMETER :: Pi = 3.1415927
        real, dimension(2,2) :: A
        real, dimension(2,1) :: b
        real, dimension(2) :: x, vector, vecNormal
        integer :: LWMAX, LWORK, Info
        parameter        ( LWMAX = 10000)
        double precision, dimension(:), allocatable ::  WORK, Ipiv
    
        ! Body of getDelaunayCoordBound
        ! LAPACK Parameters
        allocate(Ipiv(IV%NoSnap),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        allocate(Work(LWMAX),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ComputeRBFWeights "
        LWORK = -1
        Ipiv = 0.0
        Info = 0
        
        ! Number of Delaunay Boundary Points
        ndbp = 8                         
        allocate(nodesvec(2*RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(nodesvec2(2*RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        
        ! Separate internal and external boundary
        k = 0
        l = 0
        do i = 1, RD%nbf
            if (RD%boundf(i,3) == 3) then
                k = k + 1
                nodesvec(k) = RD%boundf(i,1)
                k = k + 1
                nodesvec(k) = RD%boundf(i,2)
            else
                l = l + 1
                nodesvec2(l) = RD%boundf(i,1)
                l = l + 1
                nodesvec2(l) = RD%boundf(i,2)
            end if
        end do
        call QSortInt(nodesvec(1:k), k, 'n') 
        call UniqueInt(nodesvec(1:k),k, OuterBound)
        call QSortInt(nodesvec2(1:l), l, 'n') 
        call UniqueInt(nodesvec2(1:l), l, InnerBound)
        nobp = size(OuterBound)
        nibp = size(InnerBound)
        
        ! Identify closest nodes to input CN coordinates in Mesh
        allocate(dist(nibp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(CP_ind(IV%NoCP),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "      
        do i = 1, IV%NoCP
            do j = 1, nibp
                dist(j) = DistP2P(IV%NoDim, RD%Coord_CP(i,1), RD%Coord(InnerBound(j),1), RD%Coord_CP(i, 2), RD%Coord(InnerBound(j),2))  ! Calculate Distances          
            end do
            CP_ind(i) = minloc(dist,dim=1) ! Returns the index(Position of Node in Coord Matrix) of the minimum Value 
        end do
        deallocate(dist)
        
        ! Extract moving boundary nodes
        nodesvec = 0
        nodesvec2 = 0
        k = 0
        l = 0
        overlap = 1
        if (RD%NoParts == 0) then ! If all parts move
            allocate(MovingGeomIndex(l),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
            MovingGeomIndex = InnerBound
            allocate(DelaunayCoordBound(IV%NoCP+ndbp,IV%NoDim),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        else    ! if only some parts of the boundary move
            do i = 1, RD%nbf
                do j = 1, RD%NoParts
                    if (RD%MovingParts(j) == RD%boundf(i,4)) then
                        k = k + 1
                        nodesvec(k) = RD%boundf(i,1)
                        k = k + 1
                        nodesvec(k) = RD%boundf(i,2)
                    else
                        l = l + 1
                        nodesvec2(l) = RD%boundf(i,1)
                        l = l + 1
                        nodesvec2(l) = RD%boundf(i,2)                   
                    end if
                end do
            end do
            call QSortInt(nodesvec(1:k), k, 'n') 
            call UniqueInt(nodesvec(1:k),k, MovingGeomIndex)
            call QSortInt(nodesvec2(1:l), l, 'n') 
            call UniqueInt(nodesvec2(1:l), l, NonMovingGeomIndex)
            ! Identify overlapping nodes between moving and non-moving parts to input into Delaunay Coordinates (to ensure smooth interfaces)
            do i = 1, size(MovingGeomIndex)
                do j = 1, size(NonMovingGeomIndex)
                    if (MovingGeomIndex(i) == NonMovingGeomIndex(i)) then
                        overlap = overlap + 1
                        nodesvec(overlap) = MovingGeomIndex(i)
                    end if
                end do        
            end do
            allocate(DelaunayCoordBound(IV%NoCP + ndbp + overlap,IV%NoDim),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        end if
        deallocate(nodesvec)
        deallocate(nodesvec2)
        deallocate(NonMovingGeomIndex)
                    
        ! Integrate CN coordinates into Delaunay Coordinates required for Triangulation
        DelaunayCoordBound(1:IV%NoCP,:) = RD%Coord(CP_ind,:)
        ! Input overlapping nodes into Delaunay Coordinates
        DelaunayCoordBound((IV%NoCP + 1):(IV%NoCP + overlap),:) = RD%Coord(nodesvec(1:overlap),:)
        overlap = overlap - 1
        
        ! Find starting point for Boundary Delaunay Coordinates
        allocate(dist(size(OuterBound)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(distindex(size(OuterBound)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        maxy = maxval(RD%Coord(OuterBound,2))
        minx = minval(RD%Coord(OuterBound,1))
        maxx = maxval(RD%Coord(OuterBound,1))
        do i = 1, nobp
            if (RD%Coord(OuterBound(i),1) == minx) then
                if (testx == 0) then
                    testx = 1
                else
                    testx = 2
                    EXIT
                end if
            end if
        end do
        if (testx == 2) then     ! it is rectangular or a halfcircle
            DelaunayCoordBound((IV%NoCP + overlap + 1),:) = (/minx, maxy/)
            spacing = (maxx-minx)/2
        else                    ! it is a circle and the starting point is arbitrary
            DelaunayCoordBound((IV%NoCP + overlap + 1),:) = (/RD%Coord(OuterBound(1),1), RD%Coord(OuterBound(1),2)/)
            spacing = Pi*(maxx-minx)/8
        end if
        
        ! Evenly distribute Delaunay Coordinates on the domain
        do i = 1, (ndbp - 1)        
            do j = 1, nobp
                dist(j) = DistP2P(IV%NoDim, DelaunayCoordBound((IV%NoCP + overlap + i),1), RD%Coord(Outerbound(j), 1), DelaunayCoordBound((IV%NoCP + overlap + i),2), RD%Coord(Outerbound(j), 2))  ! Calculate Distances          
            end do
            call QSort(dist, size(dist), 'y', distindex)
            do j = 1, nobp
                if (dist(j) > spacing) then
                    Exit
                end if
            end do
            if ( i == 1) then
                DelaunayCoordBound((IV%NoCP + overlap + i + 1),:) = RD%Coord(OuterBound(distindex(j)),:)
            else
                ! Check, if direction is correct
                vector = (/DelaunayCoordBound((IV%NoCP + overlap + i),1) - DelaunayCoordBound((IV%NoCP + overlap + i - 1),1), DelaunayCoordBound((IV%NoCP + overlap + i),2) - DelaunayCoordBound((IV%NoCP + overlap + i - 1),2)/)
                vecNormal = (/ - vector(1), vector(2)/)
                A(:,1) = vector
                A(:,2) = vecNormal
                !Get Inverse of Matrix via LAPACK Library
                call dgetrf(2, 2, A, 2, Ipiv, info)
                call dgetri(2, A, 2, Ipiv, WORK, LWORK, Info)
                LWORK = min( LWMAX, int( WORK( 1 ) ) )
                call dgetri(2,A, 2, Ipiv, WORK, LWORK, Info)
                do k = 1, nobp
                    b(:,1) = (/DelaunayCoordBound((IV%NoCP + overlap + i),1) - RD%Coord(OuterBound(distindex(j)),1), DelaunayCoordBound((IV%NoCP + overlap + i),2) - RD%Coord(OuterBound(distindex(j)),2)/)
                    ! x = A(-1)*b
                    CALL DGEMM('N','N',size( A, dim = 1), size( b, dim = 2), size( A, dim = 2),alpha,A,size( A, dim = 1),b,size( b, dim = 1),beta,x,size( x, dim = 1))
                    if (x(2) < 0) then
                        DelaunayCoordBound((IV%NoCP + overlap + i + 1),:) = RD%Coord(OuterBound(distindex(j)),:)
                        EXIT
                    else
                        j = j + 1
                    end if
                end do
            end if
        end do
        
    end subroutine getDelaunayCoordBound
    
    subroutine getBoundaryIndex(BoundIndex)
    
        ! Variables
        implicit none
        integer, dimension(:), allocatable ::  BoundIndex, nodesvec, nodesvec2
    
        ! Body of getBoundaryIndex
        allocate(nodesvec(2*RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(nodesvec2(2*RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        nodesvec = (/RD%boundf(:,1), RD%boundf(:,2)/)
        call QSortInt(nodesvec, size(nodesvec), 'n') 
        call UniqueInt(nodesvec, size(nodesvec), BoundIndex)
    
    end subroutine getBoundaryIndex
    
    subroutine getDomainIndex(DomainIndex)
    
        ! Variables
        implicit none
        integer, dimension(:), allocatable ::  BoundIndex, DomainIndex
        integer :: i, j
    
        ! Body of getDomainIndex        
        call getBoundaryIndex(BoundIndex)
        allocate(DomainIndex(RD%np - size(BoundIndex)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        j = 1
        do i = 1, RD%np
            if (BoundIndex(j) == i) then
                j = j + 1
            else
                DomainIndex(i - j + 1) = i
            end if
        end do
    
    end subroutine getDomainIndex
    
end module FDGD