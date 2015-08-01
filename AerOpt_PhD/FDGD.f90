module FDGD
    
    use ReadData
    use InputData
    use Toolbox
    use CreateSnapshots
    
    double precision, dimension(:,:), allocatable :: AreaCoeffBound, AreaCoeffDomain
    double precision, dimension(:), allocatable :: CN_ind, CN_indordered 
    double precision, dimension(:,:), allocatable :: DelaunayCoordBound, DelaunayCoordDomain
    integer, dimension(:,:), allocatable :: DelaunayElem, DelaunayElemBound, DelaunayElemDomain
    integer, dimension(:), allocatable ::  MovingGeomIndex, InnerBound ! Currently the same, but sort out in future
    integer :: overlap
    
    contains
    
    recursive subroutine SubFDGD(NestDisp, NoMove, counter)
    
        !**************************************************************************
        !***                    written by David Naumann                        ***
        !***                        Swansea University                          ***
        !**************************************************************************
        ! Objective:
        ! Part A - PreMeshing
        ! 1. Create a Delaunay Triangulation with CN, interface points & domain boundary points.
        ! 2. Identify correct Delaunay Element and calculate Area Coefficients of each boundary node
        ! 3. Create a Delaunay Triangulation with ALL boundary points (far field & geometry).
        ! 4. Identify correct Delaunay Element and calculate Area Coefficients of each domain node
        !
        ! Part B - Perform FDGD motion
        ! 1. Relocate Control Nodes according to new delta (change) given in Nests.
        ! 2. Move the boundary nodes according to the CN movement
        ! 3. Validity check of Delaunay Graph: Check for intersections of triangles
        !
        !**************************************************************************

        implicit none
        double precision, dimension(maxDoF) :: NestDisp, zeros
        integer :: i, intersect, NoMove, counter
        
        ! Move Control Nodes
        call RelocateCN(NestDisp, NoMove/counter)
        
        ! Move Boundary Nodes
        call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))       
        !call CheckforIntersections()        
         
        ! Check for valid background mesh
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        call CheckforIntersections(DelaunayCoordDomain, DelaunayElemDomain, intersect)

        ! Move Domain Nodes
        if (intersect == 1) then
            !if (counter < NoMove) then
            !    counter = counter + 1
            !    call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))              
            !    call PreMeshingMid()
            !    call SubFDGD(NestDisp, NoMove, counter)               
            !else
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                NestDisp = NestDisp*(1.0/NoMove)
            !end if
        else          
            zeros = 0.0
            call RelocateCN(zeros, counter)
            call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))  
            NoMove = NoMove*2
            call SubFDGD(NestDisp, NoMove, counter)
        end if

    end subroutine SubFDGD
    
    subroutine RelocateCN(NestDisp, counter)
    
        ! Variables
        implicit none
        double precision, dimension(maxDoF) :: NestDisp
        integer :: i
        integer :: counter

        ! Relocate CN (Control Nodes) applying rotative and translative motion
        do i = 1, IV%NoCN
            DelaunayCoordBound(i,1) = RD%Coord(InnerBound(CN_ind(i)),1) + (1.0/counter)*NestDisp(i)
            DelaunayCoordBound(i,2) = RD%Coord(InnerBound(CN_ind(i)),2) + (1.0/counter)*NestDisp(i+IV%NoCN)
            ! DelaunayCoord(i,2) = RD%Coord(InnerBound(CN_ind(i)),3) + NestDisp(i+2*IV%NoCN) for z
            if (IV%angle(i) /= 0) then
                DelaunayCoordBound(i,:) = AngleofAttack((1.0/counter)*NestDisp(i+3*IV%NoCN), i)
            end if
            if (IV%CNconnecttrans(i) /= 0) then
                DelaunayCoordBound(i,1) = RD%Coord(InnerBound(CN_ind(i)),1) + (1.0/counter)*NestDisp(IV%CNconnecttrans(i))
                DelaunayCoordBound(i,2) = RD%Coord(InnerBound(CN_ind(i)),2) + (1.0/counter)*NestDisp(IV%CNconnecttrans(i)+IV%NoCN)
                ! DelaunayCoord(i,2) = RD%Coord(InnerBound(CN_ind(i)),3) + NestDisp(IV%CNconnecttrans(i)+2*IV%NoCN) for z
            end if
        end do
        !DelaunayCoordcurrent = DelaunayCoord
    
    end subroutine RelocateCN
    
    subroutine PreMeshing()
    
        ! Variables
        implicit none
        integer, dimension(:), allocatable ::  DomainIndex
    
        ! Body of PreMeshing
    
        ! Preparation for Boundary Movement
        call getDelaunayCoordBound()
        call getDelaunayElem(DelaunayElemBound, DelaunayCoordBound)
        call getAreaCoefficients(DelaunayCoordBound, DelaunayElemBound, MovingGeomIndex, size(MovingGeomIndex), AreaCoeffBound)

        ! Preparation for Domain Movement
        print *, 'Get all Boundary Nodes for Delaunay Triangulation in FDGD'
        allocate(DelaunayCoordDomain(RD%nbf,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in BruteForce "
        call getDelaunayCoordDomain(RD%Coord, size(RD%Coord, dim = 1), size(RD%Coord, dim = 2))
        call getDelaunayElem(DelaunayElemDomain, DelaunayCoordDomain)
        call getDomainIndex(DomainIndex)
        call getAreaCoefficients(DelaunayCoordDomain, DelaunayElemDomain, DomainIndex, size(DomainIndex), AreaCoeffDomain)
        
        ! Further preparation work
        call OrderBoundary()

    end subroutine PreMeshing
    
    subroutine PreMeshingMid()
    
        ! Variables
        implicit none
        integer, dimension(:), allocatable ::  DomainIndex
    
        ! Body of PreMeshing
        call PreMeshingBoundary()

        ! Preparation for Domain Movement
        print *, 'Re-Do FDGD Delaunay Triangulation'
        call getDelaunayCoordDomain(RD%Coord, size(RD%Coord, dim = 1), size(RD%Coord, dim = 2))
        deallocate(DelaunayElemDomain)
        call getDelaunayElem(DelaunayElemDomain, DelaunayCoordDomain)
        call getDomainIndex(DomainIndex)
        deallocate(AreaCoeffDomain)
        call getAreaCoefficients(DelaunayCoordDomain, DelaunayElemDomain, DomainIndex, size(DomainIndex), AreaCoeffDomain)

    end subroutine PreMeshingMid
    
    subroutine PreMeshingBoundary()
    
        ! Variables
        implicit none
    
        ! Body of PreMeshing
        DelaunayCoordBound(1:IV%NoCN,:) = RD%Coord_temp(orderedBoundaryIndex(CN_indordered),:)
        deallocate(DelaunayElemBound,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        deallocate(AreaCoeffBound,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        call getDelaunayElem(DelaunayElemBound, DelaunayCoordBound)
        call getAreaCoefficients2(DelaunayCoordBound, DelaunayElemBound, orderedBoundaryIndex, size(orderedBoundaryIndex), AreaCoeffBound)

    end subroutine PreMeshingBoundary
    
    
    subroutine getDelaunayElem(DelaunayElem, DelaunayCoord)
    ! Objective: Do the Triangulation and get the Delaunay Element Matrix
    
        ! Variables
        implicit none
        integer :: i, FileSize, NoElem, temp
        double precision :: xa, ya, xb, yb, xc, yc, S
        character(len=8) :: Filename = 'Delaunay'
        integer, dimension(:,:), allocatable :: DelaunayElem
        double precision, dimension(:,:), allocatable :: DelaunayCoord
    
        ! Body of getDelaunayElem
        open(1, file= Filename//'_nodes.txt', form='formatted', status = 'unknown')
        write(1,'(2f22.15)') transpose(DelaunayCoord)
        close(1)
    
        ! Perform Delaunay Triangulation via external executable
        allocate(character(len=100) :: strSystem,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        if (IV%systemType == 'B') then
            strSystem = '/home/'//trim(IV%UserName)//'/AerOpt/Delaunay/DelaunayTriangulation < Delaunay/DelaunayInput.txt > /dev/null'
        elseif (IV%systemType == 'Q') then
            strSystem = '/eng/cvcluster/'//trim(IV%UserName)//'/AerOpt/Delaunay/DelaunayTriangulation < Delaunay/DelaunayInput.txt > /dev/null'
        elseif (IV%systemType == 'W') then
            strSystem = 'Executables\2D_Delaunay.exe < Executables/DelaunayInput.txt >nul 2>&1'
        elseif (IV%systemType == 'L') then
            strSystem = 'Executables/2D_Delaunay < Executables/DelaunayInput.txt > /dev/null'
        else
            STOP 'INPUT ERROR: System Type selected does not exist! Program stopped.'
        end if
        call system(trim(strSystem))
        deallocate(strSystem)
       
        open(1, file= Filename//'_elements.txt', form='formatted',status='old')
        inquire(1, size = FileSize)
        if (IV%SystemType == 'W') then          
            NoElem = FileSize/32
        else
            NoElem = FileSize/31
        end if
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
    
    
    subroutine getAreaCoefficients(DelaunayCoord, DelaunayElem, RowIndex, NoP, AreaCoeff) 
    
        ! Variables
        implicit none
        double precision :: xa, ya, xb, yb, xc, yc, xp, yp, S1, S2, S3, S, e1, e2, e3
        integer :: NoNodesMove, i, j, NoP,a ,b, c, d
        integer, dimension(NoP) :: RowIndex
        double precision, dimension(:,:), allocatable :: AreaCoeff
        double precision, intent(in), dimension(:, :) :: DelaunayCoord
        integer, intent(in), dimension(:, :) :: DelaunayElem
        
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
    
    subroutine getAreaCoefficients2(DelaunayCoord, DelaunayElem, RowIndex, NoP, AreaCoeff) 
    
        ! Variables
        implicit none
        double precision :: xa, ya, xb, yb, xc, yc, xp, yp, S1, S2, S3, S, e1, e2, e3
        integer :: NoNodesMove, i, j, NoP,a ,b, c, d
        integer, dimension(NoP) :: RowIndex
        double precision, dimension(:,:), allocatable :: AreaCoeff
        double precision, intent(in), dimension(:, :) :: DelaunayCoord
        integer, intent(in), dimension(:, :) :: DelaunayElem
        
        ! Body of getAreaCoefficients        
           
        ! Identify correct Delaunay Triangle and store Area Coefficients
        allocate(AreaCoeff(NoP,5),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in BruteForce "
        do i = 1, NoP
            xp = RD%coord_temp(RowIndex(i),1)
            yp = RD%coord_temp(RowIndex(i),2)
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
                if (e1 > -10e-12 .and. e2 > -10e-12 .and. e3 > -10e-12) then
                    AreaCoeff(i,:) = (/ dble(RowIndex(i)), dble(j), e1, e2, e3/)
                    EXIT
                end if
            end do
        end do
        
        ! Test Output
        !open(1, file= newdir//'/'//'Area_Coefficients.txt',form='formatted',status='unknown')
        !write(1,'(5f22.15)') transpose(AreaCoeff)
        !close(1)

    end subroutine getAreaCoefficients2
    
    subroutine CheckforIntersections(DelaunayCoord, DelaunayElem, intersect)
    
        ! Variables
        implicit none
        double precision :: xa, ya, xb, yb, xc, yc, S
        integer :: j, NoElem, intersect
        double precision, dimension(:,:) :: DelaunayCoord
        integer, dimension(:,:) :: DelaunayElem
    
        ! Body of CheckforIntersections
        intersect = 1
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
                print *, 'Intersection identified. Movement will be reduced.'
                intersect = 0
                EXIT
            end if
        end do
    
    end subroutine CheckforIntersections
    
    subroutine RelocateMeshPoints(DelaunayCoord, DelaunayElem, AreaCoeff, NoP)
    
        ! Variables
        implicit none
        integer :: i, NoP
        double precision :: x1, y1, x2, y2, x3, y3, xp, yp, inp
        double precision, dimension(NoP,5) :: AreaCoeff
        double precision, dimension(:,:) :: DelaunayCoord
        integer, dimension(:,:) :: DelaunayElem

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
        DelaunayCoordDomain = Coord(BoundIndex,:)
    
    end subroutine getDelaunayCoordDomain
    
    subroutine getDelaunayCoordBound()
    
        ! Variables
        implicit none
        integer :: i, j, k, l, nbp, nobp, testx, testy, circ, lin, nibp
        double precision :: maxx, maxy, minx, miny, spacing
        double precision, dimension(:), allocatable :: dist
        integer, dimension(:), allocatable ::  OuterBound, nodesvec, nodesvec2,nodesvec3, nodesvec4, NonMovingGeomIndex, distindex
        real, PARAMETER :: Pi = 3.1415927
        logical :: mp
        double precision, dimension(:,:), allocatable :: IdealCoord
    
        ! Body of getDelaunayCoordBound
        
        ! Number of Delaunay Boundary Points                       
        allocate(nodesvec(2*RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(nodesvec2(2*RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        
        ! Separate internal and external boundary
        k = 0
        l = 0
        do i = 1, RD%nbf
            if (RD%boundtype(i) == 3 .or. RD%boundtype(i) == 4) then
                k = k + 1
                nodesvec(k) = RD%bound(i,1)
                k = k + 1
                nodesvec(k) = RD%bound(i,2)
            else
                l = l + 1
                nodesvec2(l) = RD%bound(i,1)
                l = l + 1
                nodesvec2(l) = RD%bound(i,2)
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
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound "
        allocate(CN_ind(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound "      
        do i = 1, IV%NoCN
            do j = 1, nibp
                dist(j) = DistP2P(IV%NoDim, RD%Coord_CN(i,1), RD%Coord(InnerBound(j),1), RD%Coord_CN(i, 2), RD%Coord(InnerBound(j),2))  ! Calculate Distances          
            end do
            CN_ind(i) = minloc(dist,dim=1) ! Returns the index(Position of Node in Coord Matrix) of the minimum Value
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
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound "
            MovingGeomIndex = InnerBound
            allocate(NoNMovingGeomIndex(nobp),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound "
            NoNMovingGeomIndex = OuterBound 
        else    ! if only some parts of the boundary move
            do i = 1, RD%nbf
                do j = 1, RD%NoParts
                    if (RD%MovingParts(j) == RD%boundpart(i)) then
                        k = k + 1
                        nodesvec(k) = RD%bound(i,1)
                        k = k + 1
                        nodesvec(k) = RD%bound(i,2)
                        mp = .true.
                        EXIT
                    end if                  
                end do
                if (mp == .false.) then
                    l = l + 1
                    nodesvec2(l) = RD%bound(i,1)
                    l = l + 1
                    nodesvec2(l) = RD%bound(i,2)
                else
                    mp = .false.
                end if
            end do
            allocate(nodesvec3(k),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound "
            allocate(nodesvec4(l),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound "
            nodesvec3 = nodesvec
            nodesvec4 = nodesvec2
            call QSortInt(nodesvec3, k, 'n') 
            call UniqueInt(nodesvec3,k, MovingGeomIndex)
            call QSortInt(nodesvec4, l, 'n') 
            call UniqueInt(nodesvec4, l, NonMovingGeomIndex)
            deallocate(nodesvec3)
            deallocate(nodesvec4)
        end if
        
        ! Identify overlapping nodes between moving and non-moving parts to input into Delaunay Coordinates (to ensure smooth interfaces)
        do i = 1, size(MovingGeomIndex)
            do j = 1, size(NonMovingGeomIndex)
                if (MovingGeomIndex(i) == NonMovingGeomIndex(j)) then
                    do k = 1, IV%NoCN           ! Check, that identified Coordinate is not already picked as CN
                        if (MovingGeomIndex(i) == InnerBound(CN_ind(k))) then
                            testx = 2
                            EXIT
                        end if
                    end do
                    if (testx /= 2) then
                        overlap = overlap + 1
                        nodesvec(overlap) = MovingGeomIndex(i)
                    end if
                end if
            end do        
        end do
        deallocate(NonMovingGeomIndex)
        allocate(DelaunayCoordBound(IV%NoCN + IV%NoDelBP + overlap,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound"
        
                       
        ! Integrate CN coordinates into Delaunay Coordinates required for Triangulation
        DelaunayCoordBound(1:IV%NoCN,:) = RD%Coord(InnerBound(CN_ind),:)
        ! Input overlapping nodes into Delaunay Coordinates
        if (overlap /= 0) then
            DelaunayCoordBound((IV%NoCN + 1):(IV%NoCN + overlap),:) = RD%Coord(nodesvec(1:overlap),:)
        end if
        deallocate(nodesvec)
        deallocate(nodesvec2)
    
        ! Find Points for Boundary Delaunay Coordinates
        maxy = maxval(RD%Coord(OuterBound,2))
        miny = minval(RD%Coord(OuterBound,2))
        minx = minval(RD%Coord(OuterBound,1))
        maxx = maxval(RD%Coord(OuterBound,1))
        testx = 0
        do i = 1, nobp
            if (RD%Coord(OuterBound(i),1) == minx) then
                if (RD%Coord(OuterBound(i),2) == miny .or. RD%Coord(OuterBound(i),2) == maxy) then
                    testx = testx + 1
                end if
            elseif (RD%Coord(OuterBound(i),1) == maxx) then
                if (RD%Coord(OuterBound(i),2) == miny .or. RD%Coord(OuterBound(i),2) == maxy) then
                    testx = testx + 1
                end if
            end if
        end do
        if (testx == 4) then   ! it is rectangular
            DelaunayCoordBound((IV%NoCN + overlap + 1),:) = (/minx, maxy/)
            spacing = (maxx - minx)/(IV%NoDelBP/4)
            j = 0
            allocate(IdealCoord(IV%NoDelBP,IV%NoDim),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound"
            IdealCoord = 0.0
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/(minx + i*spacing), maxy/)
                if ((IdealCoord(j,1) - maxx) < 10e-8) then
                    EXit
                end if
            end do
            spacing = (maxy - miny)/(IV%NoDelBP/4)
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/maxx, (maxy - i*spacing)/)
                if ((IdealCoord(j,2) - miny) < 10e-8) then
                    EXit
                end if
            end do
            spacing = (maxx - minx)/(IV%NoDelBP/4)
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/(maxx - i*spacing), miny/)
                if ((IdealCoord(j,1) - minx) < 10e-8) then
                    EXit
                end if
            end do
            spacing = (maxy - miny)/(IV%NoDelBP/4)
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/minx, (miny + i*spacing)/)
                if ((IdealCoord(j,2) - maxy) < 10e-8) then
                    EXit
                end if
            end do
            ! Identify real nodes to ideal Delaunay Boundary coordinates in Mesh
            allocate(dist(nobp),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound " 
            do i = 1, (IV%NoDelBP - 1)
                do j = 1, nobp
                    dist(j) = DistP2P(IV%NoDim, IdealCoord(i,1), RD%Coord(OuterBound(j),1), IdealCoord(i,2), RD%Coord(OuterBound(j),2))  ! Calculate Distances          
                end do
                DelaunayCoordBound((IV%NoCN + overlap + i + 1),:) = RD%Coord(OuterBound(minloc(dist,dim=1)),:)
            end do
            deallocate(dist)
        elseif (testx == 2) then  ! it is a halfcircle
            DelaunayCoordBound((IV%NoCN + overlap + 1),:) = (/minx, maxy/)
            circ = nint(0.61*(IV%NoDelBP + 2))
            spacing = sqrt(((maxy-miny)**2)/2 - ((maxy-miny)**2)/2*cos(Pi/(circ - 1.0))) ! HalfCircle
            call DistributeDomainDelaunayCoord(spacing, OuterBound, nobp, overlap, 1, (circ - 2), dble((/ 1.0, 0.5/)))
            DelaunayCoordBound((IV%NoCN + overlap + circ),:) = (/minx, miny/)
            lin = (IV%NoDelBP + 2) - circ
            spacing = (maxy - miny)/((lin - 1) + 0.01)  ! Line       
            call DistributeDomainDelaunayCoord(spacing, OuterBound, nobp, overlap, circ, (IV%NoDelBP - 1), dble((/ -1.0, 0.5/)))
        else ! it is a circle and the starting point is arbitrary
            DelaunayCoordBound((IV%NoCN + overlap + 1),:) = (/RD%Coord(OuterBound(maxloc(RD%Coord(OuterBound,2))),1), maxy /)
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
                dist(j) = DistP2P(IV%NoDim, DelaunayCoordBound((IV%NoCN + overlap + i),1), RD%Coord(Outerbound(j), 1), DelaunayCoordBound((IV%NoCN + overlap + i),2), RD%Coord(Outerbound(j), 2))  ! Calculate Distances          
            end do
            call QSort(dist, size(dist), 'y', distindex)
            do j = 1, nobp
                if (dist(j) .ge. spacing) then
                    Exit
                end if
            end do
 
            if ( i == start) then ! right
                vector = direction
                !DelaunayCoordBound((IV%NoCN + overlap + i + 1),:) = RD%Coord(OuterBound(distindex(j)),:)
            else
                vector = (/DelaunayCoordBound((IV%NoCN + overlap + i),1) - DelaunayCoordBound((IV%NoCN + overlap + i - 1),1), DelaunayCoordBound((IV%NoCN + overlap + i),2) - DelaunayCoordBound((IV%NoCN + overlap + i - 1),2)/)
            end if
            
            ! Check, if direction is correct
            vecNormal = (/ - vector(2), vector(1)/)
            A(:,1) = vector
            A(:,2) = vecNormal
            A = inv(A)
            do k = 1, (nobp - j)                 
                b(:,1) = (/DelaunayCoordBound((IV%NoCN + overlap + i),1) - RD%Coord(OuterBound(distindex(j)),1), DelaunayCoordBound((IV%NoCN + overlap + i),2) - RD%Coord(OuterBound(distindex(j)),2)/)
               
                ! x = A(-1)*b
                x = matmul(A, b)
                if (x(1,1) < 0) then
                    DelaunayCoordBound((IV%NoCN + overlap + i + 1),:) = RD%Coord(OuterBound(distindex(j)),:)
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
        nodesvec = (/RD%bound(:,1), RD%bound(:,2)/)
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
    
    function AngleofAttack(alpha, CNindex)
    
        ! Variables
        implicit none
        integer :: CNindex
        double precision, dimension(2) :: AngleofAttack
        double precision, dimension(2,2) :: T
        double precision, dimension(2) :: P1, P2, P21
        double precision :: alpha, ralpha
        real, PARAMETER :: Pi = 3.1415927

        ! Body of AngleofAttack
        if (IV%CNconnectangle(CNindex) == 0) then
            STOP "No Control Node connection is given for a specified angle input!"
        end if
        ralpha = alpha*Pi/180
        T(1,:) = (/cos(ralpha), sin(ralpha)/)
        T(2,:) = (/-sin(ralpha), cos(ralpha)/)
        P1 = RD%Coord(InnerBound(CN_ind(IV%CNconnectangle(CNindex))),:)
        P2 = RD%Coord(InnerBound(CN_ind(CNindex)),:)
        P21 = P2 - P1
        AngleofAttack = matmul( T, P21) + P1
         
    end function AngleofAttack
    
    subroutine OrderBoundary()
    
        ! Variables
        implicit none
        integer :: indOrder, indjump, i, nmbf, j
        double precision, dimension(1000,2) :: jumped
        double precision, dimension(:), allocatable :: dist
    
        ! Body of OrderBoundary
        nmbf = size(MovingGeomIndex)
        allocate(orderedBoundaryIndex(nmbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        do j = 1, nmbf
            if (RD%boundtype(j) /= 3 .and. RD%boundtype(j) /= 4) then
                EXIT
            end if 
        end do
        orderedBoundaryIndex(1) = RD%bound(j,1)
        orderedBoundaryIndex(2) = RD%bound(j,2) 
        indOrder = 2
        indjump = 0
        do i = (j+1), nmbf+j
            call RecOrderBoundary(i, indOrder, jumped, indjump, orderedBoundaryIndex)
        end do
        
        ! Evaluate CN_indordered for ordered boundary
        ! Identify closest nodes to input CN coordinates in Mesh
        allocate(dist(nmbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound "
        allocate(CN_indordered(IV%NoCN),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound "      
        do i = 1, IV%NoCN
            do j = 1, nmbf
                dist(j) = DistP2P(IV%NoDim, RD%Coord_CN(i,1), RD%Coord(orderedBoundaryIndex(j),1), RD%Coord_CN(i, 2), RD%Coord(orderedBoundaryIndex(j),2))  ! Calculate Distances          
            end do
            CN_indordered(i) = minloc(dist,dim=1) ! Returns the index(Position of Node in Coord Matrix) of the minimum Value
        end do
        deallocate(dist)
        
    end subroutine OrderBoundary
    
    recursive subroutine RecOrderBoundary(i, indOrder, jumped, indjump, oBI)
    
        ! Variables
        implicit none
        integer :: indOrder, indjump, i, j
        logical :: match
        double precision, dimension(1000,2) :: jumped
        integer, dimension(*) :: oBI
    
        ! Body of RecOrderBoundary
        if (oBI(indOrder) == RD%bound(i,1)) then
            indOrder = indOrder + 1    
            oBI(indOrder) = RD%bound(i,2)
        else
            match = .false.
            do j = 1, indjump
                if (oBI(indOrder) == jumped(j,1)) then
                    indOrder = indOrder + 1  
                    oBI(indOrder) = jumped(j,2)
                    jumped(j,:) = jumped(indjump,:)
                    indjump = indjump - 1
                    match = .true.
                    exit
                end if
            end do
            if (match == .false.) then
                indjump = indjump + 1
                jumped(indjump,:) = RD%bound(i,:)
            else
                call RecOrderBoundary(i, indOrder, jumped, indjump, oBI)    
            end if
        end if
            
    end subroutine RecOrderBoundary
    
end module FDGD