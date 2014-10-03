module FDGD
    
    use ReadData
    use InputData
    use Toolbox
    use CreateSnapshots
    
    double precision, dimension(:,:), allocatable :: AreaCoeffBound, AreaCoeffDomain
    double precision, dimension(:,:), allocatable :: DelaunayElem, DelaunayCoord
    double precision, dimension(:,:), allocatable :: DelaunayCoordBound, DelaunayElemBound, DelaunayElemDomain
    integer, dimension(:), allocatable ::  InnerBound
    
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
        if (IV%alpha == 0) then
            do i = 1, IV%NoCP
                DelaunayCoord(i,1) = DelaunayCoord(i,1) + NestDisp(i)
                DelaunayCoord(i,2) = DelaunayCoord(i,2) + NestDisp(i+IV%NoCP)
            end do
        else           
            call AngleofAttack(NestDisp(4))
        end if

        ! Move Mesh of Boundary Nodes
        allocate(DelaunayElem(size(DelaunayElemBound, dim = 1),size(DelaunayElemBound, dim = 2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD "
        DelaunayElem = DelaunayElemBound
        call RelocateMeshPoints(AreaCoeffBound, size(AreaCoeffBound, dim = 1))       
        call CheckforIntersections()        
        deallocate(DelaunayCoord)
        deallocate(DelaunayElem)
      
        ! Move Mesh of Domain Nodes
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        allocate(DelaunayElem(size(DelaunayElemDomain, dim = 1),size(DelaunayElemDomain, dim = 2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD "
        DelaunayElem = DelaunayElemDomain
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
        allocate(DelaunayElemBound(size(DelaunayElem, dim = 1), size(DelaunayElem, dim = 2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        DelaunayElemBound = DelaunayElem
        call getAreaCoefficients(MovingGeomIndex, size(MovingGeomIndex), AreaCoeffBound)
        deallocate(DelaunayCoord)
        deallocate(DelaunayElem)

        ! Preparation for Domain Movement
        print *, 'Get all Boundary Nodes for Delaunay Triangulation in FDGD'
        call getDelaunayCoordDomain(RD%Coord, size(RD%Coord, dim = 1), size(RD%Coord, dim = 2))
        call getDelaunayElem()
        allocate(DelaunayElemDomain(size(DelaunayElem, dim = 1), size(DelaunayElem, dim = 2)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        DelaunayElemDomain = DelaunayElem
        call getDomainIndex(DomainIndex)
        call getAreaCoefficients(DomainIndex, size(DomainIndex), AreaCoeffDomain)
        deallocate(DelaunayCoord)
        deallocate(DelaunayElem)
  
    end subroutine PreMeshing
    
    
    subroutine getDelaunayElem()
    ! Objective: Do the Triangulation and get the Delaunay Element Matrix
    
        ! Variables
        implicit none
        integer :: i, FileSize, NoElem, temp
        double precision :: xa, ya, xb, yb, xc, yc, S
        character(len=8) :: Filename = 'Delaunay'
    
        ! Body of getDelaunayElem
        open(1, file= Filename//'_nodes.txt', form='formatted', status = 'unknown')
        write(1,'(2f22.15)') transpose(DelaunayCoord)
        close(1)
    
        ! Performt Delaunay Triangulation via external executable
        allocate(character(len=100) :: strSystem,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        if (IV%systemType == 'B') then
            strSystem = '/home/'//trim(IV%UserName)//'/AerOpt/DelaunayTriangulation < Delaunay/DelaunayInput.txt > /dev/null'
        else
            strSystem = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/DelaunayTriangulation < Delaunay/DelaunayInput.txt > /dev/null'
        end if
        call system(trim(strSystem))
        deallocate(strSystem)
        
        open(1, file= Filename//'_elements.txt', form='formatted',status='old')
        inquire(1, size = FileSize)
        NoElem = FileSize/31
        allocate(DelaunayElem(NoElem,IV%NoDim + 1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in BruteForce "
        do i = 1, NoElem
            read(1, *) DelaunayElem(i,:)
        end do        
        close(1)
        
        do i = 1, NoElem
            xa = DelaunayCoord(DelaunayElem(i,1),1)
            ya = DelaunayCoord(DelaunayElem(i,1),2)
            xb = DelaunayCoord(DelaunayElem(i,2),1)
            yb = DelaunayCoord(DelaunayElem(i,2),2)
            xc = DelaunayCoord(DelaunayElem(i,3),1)
            yc = DelaunayCoord(DelaunayElem(i,3),2)
            S = 0.5*((xa*yb + ya*xc + xb*yc) - (xa*yc + ya*xb + yb*xc))
            if (S <= 0) then  ! If Area negative, the order of the nodes is incorrect and positions will be swapped
                print *, 'Area negative'
                temp = DelaunayElem(i,1)
                DelaunayElem(i,1) = DelaunayElem(i,3)
                DelaunayElem(i,3) = temp
            end if
        end do
        
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
            do j = 1, size(DelaunayElem, dim = 1)
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
                e2 = S2/S
                e3 = S3/S
                if (e1 >= 0.0 .and. e2 >= 0.0 .and. e3 >= 0.0) then
                    AreaCoeff(i,:) = (/ dble(RowIndex(i)), dble(j), e1, e2, e3/)
                    EXIT
                end if
            end do
        end do
        
        ! Test Output
        !open(1, file= newdir//'/'//'Area_Coefficients.txt',form='formatted',status='unknown')
        !write(1,'(5f22.15)') transpose(AreaCoeff)
        !close(1)

    end subroutine getAreaCoefficients
    
    
    subroutine CheckforIntersections()
    
        ! Variables
        implicit none
        double precision :: xa, ya, xb, yb, xc, yc, S
        integer :: j, NoElem
    
        ! Body of CheckforIntersections
        NoElem = size(DelaunayElem, dim = 1)
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
   
    end subroutine RelocateMeshPoints
    
    subroutine getDelaunayCoordDomain(Coord, dim1, dim2)
    
        ! Variables
        implicit none  
        integer :: dim1, dim2
        double precision :: Coord(dim1, dim2)
        integer, dimension(:), allocatable :: BoundIndex
        
        ! Body of getDelaunayCoordDomain
        !Identify Boundary Indices
        call getBoundaryIndex(BoundIndex)
        allocate(DelaunayCoord(size(BoundIndex),IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        DelaunayCoord = Coord(BoundIndex,:)
    
    end subroutine getDelaunayCoordDomain
    
    subroutine getDelaunayCoordBound(MovingGeomIndex)
    
        ! Variables
        implicit none
        integer :: i, j, k, l, nbp, nibp, nobp, overlap, testx, testy, circ, lin
        double precision :: maxx, maxy, minx, miny, spacing
        double precision, dimension(:), allocatable :: dist, CP_ind
        integer, dimension(:), allocatable ::  OuterBound, nodesvec, nodesvec2,nodesvec3, nodesvec4, NonMovingGeomIndex, distindex, MovingGeomIndex
        real, PARAMETER :: Pi = 3.1415927
        logical :: mp
        integer :: LWMAX, LWORK, Info
        parameter        ( LWMAX = 10000)
        double precision, dimension(:), allocatable ::  WORK, Ipiv
        double precision, dimension(:,:), allocatable :: IdealCoord
    
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
        allocate(nodesvec3(k),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(nodesvec4(l),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        nodesvec3 = nodesvec
        nodesvec4 = nodesvec2
        call QSortInt(nodesvec3, k, 'n') 
        call UniqueInt(nodesvec3,k, OuterBound)
        call QSortInt(nodesvec4, l, 'n') 
        call UniqueInt(nodesvec4, l, InnerBound)
        deallocate(nodesvec3)
        deallocate(nodesvec4)
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
        overlap = 0
        mp = .false.
        if (RD%NoParts == 0) then ! If all parts move
            allocate(MovingGeomIndex(nibp),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
            MovingGeomIndex = InnerBound
            allocate(DelaunayCoordBound(IV%NoCP+IV%NoDelBP,IV%NoDim),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
 
        else    ! if only some parts of the boundary move
            do i = 1, RD%nbf
                do j = 1, RD%NoParts
                    if (RD%MovingParts(j) == RD%boundf(i,4)) then
                        k = k + 1
                        nodesvec(k) = RD%boundf(i,1)
                        k = k + 1
                        nodesvec(k) = RD%boundf(i,2)
                        mp = .true.
                        EXIT
                    end if                  
                end do
                if (mp == .false.) then
                    l = l + 1
                    nodesvec2(l) = RD%boundf(i,1)
                    l = l + 1
                    nodesvec2(l) = RD%boundf(i,2)
                else
                    mp = .false.
                end if
            end do
            allocate(nodesvec3(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
            allocate(nodesvec4(l),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
            nodesvec3 = nodesvec
            nodesvec4 = nodesvec2
            call QSortInt(nodesvec3, k, 'n') 
            call UniqueInt(nodesvec3,k, MovingGeomIndex)
            call QSortInt(nodesvec4, l, 'n') 
            call UniqueInt(nodesvec4, l, NonMovingGeomIndex)
            deallocate(nodesvec3)
            deallocate(nodesvec4)

            ! Identify overlapping nodes between moving and non-moving parts to input into Delaunay Coordinates (to ensure smooth interfaces)
            do i = 1, size(MovingGeomIndex)
                do j = 1, size(NonMovingGeomIndex)
                    if (MovingGeomIndex(i) == NonMovingGeomIndex(j)) then
                        overlap = overlap + 1
                        nodesvec(overlap) = MovingGeomIndex(i)
                    end if
                end do        
            end do
            deallocate(NonMovingGeomIndex)
            
            allocate(DelaunayCoordBound(IV%NoCP + IV%NoDelBP + overlap,IV%NoDim),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        end if
        
        overlap = overlap + 1                
        ! Integrate CN coordinates into Delaunay Coordinates required for Triangulation
        DelaunayCoordBound(1:IV%NoCP,:) = RD%Coord(InnerBound(CP_ind),:)
        ! Input overlapping nodes into Delaunay Coordinates
        DelaunayCoordBound((IV%NoCP + 1):(IV%NoCP + overlap),:) = RD%Coord(nodesvec(1:overlap),:)
        overlap = overlap - 1
        deallocate(nodesvec)
        deallocate(nodesvec2)
    
        ! Find Points for Boundary Delaunay Coordinates
        maxy = maxval(RD%Coord(OuterBound,2))
        miny = minval(RD%Coord(OuterBound,2))
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
        do i = 1, nobp
            if (RD%Coord(OuterBound(i),1) == miny) then
                if (testy == 0) then
                    testy = 1
                else
                    testy = 2
                    EXIT
                end if
            end if
        end do
        if (testx == 2 .and. testy == 2) then   ! it is rectangular
            DelaunayCoordBound((IV%NoCP + overlap + 1),:) = (/minx, maxy/)
            spacing = (maxx - minx)/(IV%NoDelBP/4)
            j = 0
            allocate(IdealCoord(IV%NoDelBP,IV%NoDim),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/(minx + i*spacing), maxy/)
                if (IdealCoord(j,1) == maxx) then
                    EXit
                end if
            end do
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/maxx, (maxy - i*spacing)/)
                if (IdealCoord(j,2) == miny) then
                    EXit
                end if
            end do
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/(maxx - i*spacing), miny/)
                if (IdealCoord(j,1) == minx) then
                    EXit
                end if
            end do
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/minx, (miny + i*spacing)/)
                if (IdealCoord(j,2) == maxy) then
                    EXit
                end if
            end do
            ! Identify real nodes to ideal Delaunay Boundary coordinates in Mesh
            allocate(dist(nobp),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh " 
            do i = 1, (IV%NoDelBP - 1)
                do j = 1, nobp
                    dist(j) = DistP2P(IV%NoDim, IdealCoord(i,1), RD%Coord(OuterBound(j),1), IdealCoord(i,2), RD%Coord(OuterBound(j),2))  ! Calculate Distances          
                end do
                DelaunayCoordBound((IV%NoCP + overlap + i + 1),:) = RD%Coord(OuterBound(minloc(dist,dim=1)),:)
            end do
            deallocate(dist)
        elseif (testx == 2) then  ! it is a halfcircle
            DelaunayCoordBound((IV%NoCP + overlap + 1),:) = (/minx, maxy/)
            circ = nint(0.61*(IV%NoDelBP + 2))
            spacing = sqrt(((maxy-miny)**2)/2 - ((maxy-miny)**2)/2*cos(Pi/(circ - 1.0))) ! HalfCircle
            call DistributeDomainDelaunayCoord(spacing, OuterBound, nobp, overlap, 1, (circ - 2), dble((/ 1.0, 0.5/)))
            DelaunayCoordBound((IV%NoCP + overlap + circ),:) = (/minx, miny/)
            lin = (IV%NoDelBP + 2) - circ
            spacing = (maxy - miny)/((lin - 1) + 0.01)  ! Line       
            call DistributeDomainDelaunayCoord(spacing, OuterBound, nobp, overlap, circ, (IV%NoDelBP - 1), dble((/ -1.0, 0.5/)))
        else ! it is a circle and the starting point is arbitrary
            DelaunayCoordBound((IV%NoCP + overlap + 1),:) = (/RD%Coord(OuterBound(maxloc(RD%Coord(OuterBound,2))),1), maxy /)
            spacing = sqrt(((maxy-miny)**2)/2 - (maxy-miny)**2)*cos(real(360/IV%NoDelBP))
            call DistributeDomainDelaunayCoord(spacing, OuterBound, nobp,  overlap, 1, (IV%NoDelBP - 1), dble((/ 0.0, -1.0/)))
        end if
    end subroutine getDelaunayCoordBound
    
    subroutine DistributeDomainDelaunayCoord(spacing, OuterBound, nobp, overlap, start, ending, direction)
    
        ! Variables
        implicit none
        integer :: nobp, nb, overlap, i, j, k, start, ending
        double precision :: spacing
        double precision, dimension(2,2) :: A
        double precision, dimension(2,1) :: b
        double precision, dimension(2) :: vector, vecNormal, direction
        double precision, dimension(2,1) :: x
        double precision, dimension(:), allocatable :: dist
        integer, dimension(nobp) ::  OuterBound, distindex
    
        ! Body of DistributeDomainDelaunayCoord
        allocate(dist(nobp),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        
        ! Evenly distribute Delaunay Coordinates on the domain
        do i = start, ending

            distindex = (/ (j, j=1,nobp) /)            
            do j = 1, nobp
                dist(j) = DistP2P(IV%NoDim, DelaunayCoordBound((IV%NoCP + overlap + i),1), RD%Coord(Outerbound(j), 1), DelaunayCoordBound((IV%NoCP + overlap + i),2), RD%Coord(Outerbound(j), 2))  ! Calculate Distances          
            end do
            call QSort(dist, size(dist), 'y', distindex)
            do j = 1, nobp
                if (dist(j) .ge. spacing) then
                    Exit
                end if
            end do
 
            if ( i == start) then ! right
                vector = direction
                !DelaunayCoordBound((IV%NoCP + overlap + i + 1),:) = RD%Coord(OuterBound(distindex(j)),:)
            else
                vector = (/DelaunayCoordBound((IV%NoCP + overlap + i),1) - DelaunayCoordBound((IV%NoCP + overlap + i - 1),1), DelaunayCoordBound((IV%NoCP + overlap + i),2) - DelaunayCoordBound((IV%NoCP + overlap + i - 1),2)/)
            end if
            
            ! Check, if direction is correct
            vecNormal = (/ - vector(2), vector(1)/)
            A(:,1) = vector
            A(:,2) = vecNormal
            A = inv(A)
            do k = 1, (nobp - j)                 
                b(:,1) = (/DelaunayCoordBound((IV%NoCP + overlap + i),1) - RD%Coord(OuterBound(distindex(j)),1), DelaunayCoordBound((IV%NoCP + overlap + i),2) - RD%Coord(OuterBound(distindex(j)),2)/)
               
                ! x = A(-1)*b
                x = matmul(A, b)
                if (x(1,1) < 0) then
                    DelaunayCoordBound((IV%NoCP + overlap + i + 1),:) = RD%Coord(OuterBound(distindex(j)),:)
                    EXIT
                else
                    j = j + 1
                end if
            end do
            
        end do
    
    end subroutine DistributeDomainDelaunayCoord
    
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
        integer :: i, j, k
    
        ! Body of getDomainIndex
        call getBoundaryIndex(BoundIndex)
        allocate(DomainIndex(RD%np - size(BoundIndex)),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        j = 1
        k = 1
        do i = 1, RD%np
            if (BoundIndex(j) == i) then ! if part of Bound, do not include in Domain
                if (j < size(BoundIndex)) then ! necessary to avoid Out Of Array Bounds request
                    j = j + 1
                else
                    k = 0
                end if
            else
                DomainIndex(i - j + k) = i
            end if
        end do

    end subroutine getDomainIndex
    
    subroutine AngleofAttack(alpha)
    
        ! Variables
        implicit none
        double precision, dimension(2,2) :: T
        double precision, dimension(2) :: P1, P2, P21
        double precision :: alpha, ralpha
        real, PARAMETER :: Pi = 3.1415927

        ! Body of AngleofAttack
        ralpha = alpha*Pi/180
        T(1,:) = (/cos(ralpha), sin(ralpha)/)
        T(2,:) = (/-sin(ralpha), cos(ralpha)/)
        P1 = (/DelaunayCoordBound(1,1), DelaunayCoordBound(1,2)/)
        P2 = (/DelaunayCoordBound(2,1), DelaunayCoordBound(2,2)/)
        P21 = P2 - P1
        DelaunayCoord(2,:) = matmul( T, P21) + P1
         
    end subroutine AngleofAttack
    
end module FDGD