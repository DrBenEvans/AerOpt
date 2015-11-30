module FDGD
    
    use ReadData
    use InputData
    use Toolbox
    use CreateSnapshots
    
    double precision, dimension(:,:), allocatable :: AreaCoeffBound, AreaCoeffDomain
    double precision, dimension(:), allocatable :: CN_ind, CN_indordered
    double precision, dimension(:,:), allocatable :: DelaunayCoordBound, DelaunayCoordDomain, dRot
    integer, dimension(1000,3) :: DelaunayElemBound
    integer, dimension(10000,3) :: DelaunayElemDomain
    integer, dimension(:), allocatable ::  InnerBound
    integer :: overlap, nbp
    integer, dimension(:), allocatable :: orderedBoundaryIndex   ! Vector containing the Index of boundary points in the correct order
    
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
        call RelocateCN(NestDisp, NoMove, counter)
        
        ! Move Boundary Nodes
        call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))       
        !call CheckforIntersections()        
         
        ! Rotate
        do i = 1, IV%NoCN
            if (IV%angle(i) /= 0) then
                call AngleofAttack(NestDisp(4*i)/NoMove, i)
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
                call SubFDGD(NestDisp, NoMove, counter)               
            else
                call RelocateMeshPoints(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                if (NoMove > 1) then
                    call PreMeshingEnd()
                end if
                if (allocated(dRot) == .true.) then
                    deallocate(dRot)
                end if
            end if
        else          
            zeros = 0.0
            call RelocateCN(zeros, NoMove, counter)
            call RelocateMeshPoints(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))  
            counter = 2*counter-1
            NoMove = NoMove*2
            if (allocated(dRot) == .true.) then
                deallocate(dRot)
            end if
            call SubFDGD(NestDisp, NoMove, counter)
        end if

    end subroutine SubFDGD
    
    recursive subroutine SubFDGD_3D(NestDisp, NoMove, counter)
    
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
        
        ! Rotate
        !do i = 1, IV%NoCN
        !    if (IV%angle(i) /= 0) then
        !        call AngleofAttack_3D(CNDisp(4*i)/NoMove*counter,i)
        !    end if
        !end do
        
        ! Move Control Nodes
        call RelocateCN_3D(NestDisp, NoMove, counter)
        
        ! Move Boundary Nodes
        call RelocateMeshPoints_3D(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))       
        !call CheckforIntersections()        
         
        ! Check for valid background mesh
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        call CheckforIntersections_3D(DelaunayCoordDomain, DelaunayElemDomain, intersect)

        ! Move Domain Nodes
        if (intersect == 1) then
            !if (counter < NoMove) then
            !    counter = counter + 1
            !    call RelocateMeshPoints_3D(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))              
            !    call PreMeshingMid()
            !    call SubFDGD_3D(NestDisp, NoMove, counter)               
            !else
                call RelocateMeshPoints_3D(DelaunayCoordDomain, DelaunayElemDomain, AreaCoeffDomain, size(AreaCoeffDomain, dim = 1))
                NestDisp = NestDisp*(1.0/NoMove)
            !end if
        else          
            zeros = 0.0
            call RelocateCN_3D(zeros, NoMove, counter)
            call RelocateMeshPoints_3D(DelaunayCoordBound, DelaunayElemBound, AreaCoeffBound, size(AreaCoeffBound, dim = 1))  
            NoMove = NoMove*2
            call SubFDGD_3D(NestDisp, NoMove, counter)
        end if

    end subroutine SubFDGD_3D

    subroutine RelocateCN(NestDisp, NoMove, counter)
    
        ! Variables
        implicit none
        double precision, dimension(maxDoF) :: NestDisp
        integer :: i, counter, NoMove

        ! Relocate CN (Control Nodes) applying rotative and translative motion
        if (allocated(dRot) == .true.) then
            do i = 1, IV%NoCN
                DelaunayCoordBound(i,1) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),1) + (real(counter)/NoMove)*NestDisp(i) + dRot(orderedBoundaryIndex(CN_indordered(i)),1)
                DelaunayCoordBound(i,2) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),2) + (real(counter)/NoMove)*NestDisp(i+IV%NoCN) + dRot(orderedBoundaryIndex(CN_indordered(i)),2)
                if (IV%CNconnecttrans(i) /= 0) then
                    DelaunayCoordBound(i,1) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),1) + (real(counter)/NoMove)*NestDisp(IV%CNconnecttrans(i))  + dRot(orderedBoundaryIndex(CN_indordered(i)),1)
                    DelaunayCoordBound(i,2) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),2) + (real(counter)/NoMove)*NestDisp(IV%CNconnecttrans(i)+IV%NoCN)  + dRot(orderedBoundaryIndex(CN_indordered(i)),2)
                end if
            end do
        else
            do i = 1, IV%NoCN
                DelaunayCoordBound(i,1) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),1) + (real(counter)/NoMove)*NestDisp(i)
                DelaunayCoordBound(i,2) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),2) + (real(counter)/NoMove)*NestDisp(i+IV%NoCN)
                if (IV%CNconnecttrans(i) /= 0) then
                    DelaunayCoordBound(i,1) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),1) + (real(counter)/NoMove)*NestDisp(IV%CNconnecttrans(i))
                    DelaunayCoordBound(i,2) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),2) + (real(counter)/NoMove)*NestDisp(IV%CNconnecttrans(i)+IV%NoCN)
                end if
            end do
        end if
    
    end subroutine RelocateCN
    
    subroutine RelocateCN_3D(NestDisp, NoMove, counter)
    
        ! Variables
        implicit none
        double precision, dimension(maxDoF) :: NestDisp
        integer :: i
        integer :: NoMove, counter

        ! Relocate CN (Control Nodes) applying rotative and translative motion
        do i = 1, IV%NoCN
            DelaunayCoordBound(i,1) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),1) + (real(counter)/NoMove)*NestDisp(i)
            DelaunayCoordBound(i,2) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),2) + (real(counter)/NoMove)*NestDisp(i+IV%NoCN)
            DelaunayCoordBound(i,3) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),3) + (real(counter)/NoMove)*NestDisp(i+2*IV%NoCN)
            if (IV%CNconnecttrans(i) /= 0) then
                DelaunayCoordBound(i,1) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),1) + (real(counter)/NoMove)*NestDisp(IV%CNconnecttrans(i))
                DelaunayCoordBound(i,2) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),2) + (real(counter)/NoMove)*NestDisp(IV%CNconnecttrans(i)+IV%NoCN)
                DelaunayCoordBound(i,3) = RD%Coord(orderedBoundaryIndex(CN_indordered(i)),3) + (real(counter)/NoMove)*NestDisp(IV%CNconnecttrans(i)+2*IV%NoCN)
            end if
        end do
    
    end subroutine RelocateCN_3D
    
    subroutine PreMeshing()
    
        ! Variables
        implicit none
        integer, dimension(:), allocatable ::  DomainIndex
        integer :: i
    
        ! Body of PreMeshing
        
        ! Preparation for Boundary Movement
        call getDelaunayCoordBound()

        ! Further preparation work
        call OrderBoundary()

        call getDelaunayElem2(DelaunayElemBound, size(DelaunayElemBound, dim = 1), size(DelaunayElemBound, dim = 2), DelaunayCoordBound)

        call getAreaCoefficients(DelaunayCoordBound, DelaunayElemBound, orderedBoundaryIndex, size(orderedBoundaryIndex), AreaCoeffBound)

        ! Preparation for Domain Movement
        print *, 'Get all Boundary Nodes for Delaunay Triangulation in FDGD'
        allocate(DelaunayCoordDomain(RD%nbf,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in BruteForce "
        call getDelaunayCoordDomain(RD%Coord, size(RD%Coord, dim = 1), size(RD%Coord, dim = 2))
        call getDelaunayElem2(DelaunayElemDomain, size(DelaunayElemDomain, dim = 1), size(DelaunayElemDomain, dim = 2), DelaunayCoordDomain)
        call getDomainIndex(DomainIndex)
        call getAreaCoefficients(DelaunayCoordDomain, DelaunayElemDomain, DomainIndex, size(DomainIndex), AreaCoeffDomain)
       
    end subroutine PreMeshing
    
    subroutine PreMeshing_3D()
    
        ! Variables
        implicit none
        integer, dimension(:), allocatable ::  DomainIndex
    
        ! Body of PreMeshing
    
        ! Preparation for Boundary Movement
        call getDelaunayCoordBound_3D()
        !call getDelaunayElem_3D(DelaunayElemBound, DelaunayCoordBound)
        call getAreaCoefficients_3D(DelaunayCoordBound, DelaunayElemBound, InnerBound, size(InnerBound), AreaCoeffBound)

        ! Preparation for Domain Movement
        print *, 'Get all Boundary Nodes for Delaunay Triangulation in FDGD'
        allocate(DelaunayCoordDomain(nbp,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in BruteForce "
        call getDelaunayCoordDomain(RD%Coord, size(RD%Coord, dim = 1), size(RD%Coord, dim = 2))
        !call getDelaunayElem_3D(DelaunayElemDomain, DelaunayCoordDomain)
        call getDomainIndex(DomainIndex)
        call getAreaCoefficients_3D(DelaunayCoordDomain, DelaunayElemDomain, DomainIndex, size(DomainIndex), AreaCoeffDomain)
        
        ! Further preparation work - How in 3D?
        !call OrderBoundary_3D()

    end subroutine PreMeshing_3D
    
    subroutine PreMeshingMid()
    
        ! Variables
        implicit none
        integer, dimension(:), allocatable ::  DomainIndex
    
        ! Body of PreMeshing
        call PreMeshingBoundary()

        ! Preparation for Domain Movement
        call getDelaunayCoordDomain(RD%Coord_temp, size(RD%Coord_temp, dim = 1), size(RD%Coord_temp, dim = 2))
        !deallocate(DelaunayElemDomain,stat=allocateStatus)
        !if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        deallocate(AreaCoeffDomain,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        call getDomainIndex(DomainIndex)
        call getDelaunayElem2(DelaunayElemDomain, size(DelaunayElemDomain, dim = 1), size(DelaunayElemDomain, dim = 2), DelaunayCoordDomain)
        call getAreaCoefficients2(DelaunayCoordDomain, DelaunayElemDomain, DomainIndex, size(DomainIndex), AreaCoeffDomain)

    end subroutine PreMeshingMid
    
    subroutine PreMeshingBoundary()
    
        ! Variables
        implicit none

        ! Body of PreMeshing
        DelaunayCoordBound(1:IV%NoCN,:) = RD%Coord_temp(orderedBoundaryIndex(CN_indordered),:)
        !deallocate(DelaunayElemBound,stat=allocateStatus)
        !if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        deallocate(AreaCoeffBound,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        call getDelaunayElem2(DelaunayElemBound, size(DelaunayElemBound, dim = 1), size(DelaunayElemBound, dim = 2), DelaunayCoordBound)
        call getAreaCoefficients2(DelaunayCoordBound, DelaunayElemBound, orderedBoundaryIndex, size(orderedBoundaryIndex), AreaCoeffBound)

    end subroutine PreMeshingBoundary
    
    subroutine PreMeshingEnd()
    
        ! Variables
        implicit none
        integer, dimension(:), allocatable ::  DomainIndex
    
        ! Body of PreMeshing
        DelaunayCoordBound(1:IV%NoCN,:) = RD%Coord(orderedBoundaryIndex(CN_indordered),:)
        !deallocate(DelaunayElemBound,stat=allocateStatus)
        !if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        deallocate(AreaCoeffBound,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        call getDelaunayElem2(DelaunayElemBound, size(DelaunayElemBound, dim = 1), size(DelaunayElemBound, dim = 2), DelaunayCoordBound)
        call getAreaCoefficients(DelaunayCoordBound, DelaunayElemBound, orderedBoundaryIndex, size(orderedBoundaryIndex), AreaCoeffBound)

        ! Preparation for Domain Movement
        call getDelaunayCoordDomain(RD%Coord, size(RD%Coord, dim = 1), size(RD%Coord, dim = 2))
        !deallocate(DelaunayElemDomain,stat=allocateStatus)
        !if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        deallocate(AreaCoeffDomain,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        call getDomainIndex(DomainIndex)
        call getDelaunayElem2(DelaunayElemDomain, size(DelaunayElemDomain, dim = 1), size(DelaunayElemDomain, dim = 2), DelaunayCoordDomain)
        call getAreaCoefficients(DelaunayCoordDomain, DelaunayElemDomain, DomainIndex, size(DomainIndex), AreaCoeffDomain)

        end subroutine PreMeshingEnd
        
    subroutine PreMeshingBoundary_3D()
    
        ! Variables
        implicit none
! CN_indordered, how to order in 3D   
        ! Body of PreMeshing
        if (IV%NoDim == 2) then
            DelaunayCoordBound(1:IV%NoCN,:) = RD%Coord_temp(orderedBoundaryIndex(CN_indordered),:)
        else
            DelaunayCoordBound(1:IV%NoCN,:) = RD%Coord_temp(InnerBound(CN_ind),:)
        end if
        !deallocate(DelaunayElemBound,stat=allocateStatus)
        !if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        deallocate(AreaCoeffBound,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        call getDelaunayElem2(DelaunayElemBound, size(DelaunayElemBound, dim = 1), size(DelaunayElemBound, dim = 2), DelaunayCoordBound)
        call getAreaCoefficients_3D(DelaunayCoordBound, DelaunayElemBound, orderedBoundaryIndex, size(orderedBoundaryIndex), AreaCoeffBound)

    end subroutine PreMeshingBoundary_3D
        
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
        
        open(2, file = 'DelaunayInput.txt', form='formatted', status = 'unknown')
        write(2,*) 'Delaunay'
        close(2)
        
        ! Perform Delaunay Triangulation via external executable
        allocate(character(len=100) :: strSystem,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        if (IV%systemType == 'B' .or. IV%systemtype == 'Q') then
            strSystem = trim(IV%filepath)//'/Executables/Delaunay_2D < DelaunayInput.txt > /dev/null'
        elseif (IV%systemType == 'W') then
            strSystem = 'Executables\Delaunay_2D.exe < DelaunayInput.txt >nul 2>&1'
        elseif (IV%systemType == 'L') then
            strSystem = 'Executables/Delaunay_2D < DelaunayInput.txt > /dev/null'
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
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in DelaunayElem "
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
                print *, 'Area negative when generating Delaunay Background mesh'
                temp = DelaunayElem(i,1)
                DelaunayElem(i,1) = DelaunayElem(i,3)
                DelaunayElem(i,3) = temp
            end if
        end do
        
    end subroutine getDelaunayElem
    
    subroutine getDelaunayElem_3D(DelaunayElem, DelaunayCoord)
    ! Objective: Do the Triangulation and get the Delaunay Element Matrix
    
        ! Variables
        implicit none
        integer :: i, FileSize, NoElem, temp
        double precision :: xa, ya, xb, yb, xc, yc, S
        integer, dimension(:,:), allocatable :: DelaunayElem
        double precision, dimension(:,:), allocatable :: DelaunayCoord
    
        ! Body of getDelaunayElem
        open(1, file= 'Delaunay.txt', form='formatted', status = 'unknown')
        write(1,'(3f22.15)') transpose(DelaunayCoord)
        close(1)
        open(1, file = 'DelaunayInput.txt', form='formatted', status = 'unknown')
        write(1,*) 'Delaunay.txt'
        close(1)
    
        ! Perform Delaunay Triangulation via external executable
        allocate(character(len=100) :: strSystem,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        if (IV%systemType == 'B' .or. IV%systemtype == 'Q') then
            strSystem = trim(IV%filepath)//'/Executables/Delaunay_3D < DelaunayInput.txt > /dev/null'
        elseif (IV%systemType == 'W') then
            strSystem = 'Executables\Delaunay_3D.exe < DelaunayInput.txt >nul 2>&1'
        elseif (IV%systemType == 'L') then
            strSystem = 'Executables/Delaunay_3D < DelaunayInput.txt > /dev/null'
        else
            STOP 'INPUT ERROR: System Type selected does not exist! Program stopped.'
        end if
        call system(trim(strSystem))
        deallocate(strSystem)
        open(1, file= 'Delaunay.tetra.txt', form='formatted',status='unknown')
        inquire(1, size = FileSize)
        if (IV%SystemType == 'W') then          
            NoElem = FileSize/42
        else
            NoElem = FileSize/41
        end if
        allocate(DelaunayElem(NoElem,IV%NoDim + 1),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in BruteForce "
        do i = 1, NoElem
            read(1, *) DelaunayElem(i,:)
        end do        
        close(1)
        
    end subroutine getDelaunayElem_3D
        
    subroutine getAreaCoefficients_3D(DelaunayCoord, DelaunayElem, RowIndex, NoP, AreaCoeff) 
    
        ! Variables
        implicit none
        double precision :: xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd, xp, yp, zp, S1, S2, S3, S4, S, e1, e2, e3, e4
        integer :: NoNodesMove, i, j, NoP
        integer, dimension(NoP) :: RowIndex
        double precision, dimension(:,:), allocatable :: AreaCoeff
        double precision, intent(in), dimension(:, :) :: DelaunayCoord
        integer, intent(in), dimension(:, :) :: DelaunayElem
        double precision, dimension(4,4) :: A
        
        ! Body of getAreaCoefficients        
           
        ! Identify correct Delaunay Triangle and store Area Coefficients
        allocate(AreaCoeff(NoP,6),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in BruteForce "
        do i = 1, NoP
            xp = RD%coord(RowIndex(i),1)
            yp = RD%coord(RowIndex(i),2)
            zp = RD%coord(RowIndex(i),3)
            do j = 1, size(DelaunayElem, dim = 1)
                xa = DelaunayCoord(DelaunayElem(j,1),1)
                ya = DelaunayCoord(DelaunayElem(j,1),2)
                za = DelaunayCoord(DelaunayElem(j,1),3)
                xb = DelaunayCoord(DelaunayElem(j,2),1)
                yb = DelaunayCoord(DelaunayElem(j,2),2)
                zb = DelaunayCoord(DelaunayElem(j,2),3)
                xc = DelaunayCoord(DelaunayElem(j,3),1)
                yc = DelaunayCoord(DelaunayElem(j,3),2)
                zc = DelaunayCoord(DelaunayElem(j,3),3)
                xd = DelaunayCoord(DelaunayElem(j,4),1)
                yd = DelaunayCoord(DelaunayElem(j,4),2)
                zd = DelaunayCoord(DelaunayElem(j,4),3)
                A(1,:) = (/ xa, ya, za, 1/)
                A(2,:) = (/ xb, yb, zb, 1/)
                A(3,:) = (/ xc, yc, zc, 1/)
                A(4,:) = (/ xd, yd, zd, 1/)
                S = (1.0/6.0)*Det(A,4)
                A(1,:) = (/ xp, yp, zp, 1/)
                A(2,:) = (/ xb, yb, zb, 1/)
                A(3,:) = (/ xc, yc, zc, 1/)
                A(4,:) = (/ xd, yd, zd, 1/)
                S1 = (1.0/6.0)*Det(A,4)
                A(1,:) = (/ xa, ya, za, 1/)
                A(2,:) = (/ xp, yp, zp, 1/)
                A(3,:) = (/ xc, yc, zc, 1/)
                A(4,:) = (/ xd, yd, zd, 1/)
                S2 = (1.0/6.0)*Det(A,4)
                A(1,:) = (/ xa, ya, za, 1/)
                A(2,:) = (/ xb, yb, zb, 1/)
                A(3,:) = (/ xp, yp, zp, 1/)
                A(4,:) = (/ xd, yd, zd, 1/)
                S3 = (1.0/6.0)*Det(A,4)
                A(1,:) = (/ xa, ya, za, 1/)
                A(2,:) = (/ xb, yb, zb, 1/)
                A(3,:) = (/ xc, yc, zc, 1/)
                A(4,:) = (/ xp, yp, zp, 1/)
                S4 = (1.0/6.0)*Det(A,4)
                e1 = S1/S
                e2 = S2/S
                e3 = S3/S
                e4 = S4/S
                if (e1 > -10e-10 .and. e2 > -10e-10 .and. e3 > -10e-10 .and. e4 > -10e-10) then
                    AreaCoeff(i,:) = (/ dble(RowIndex(i)), dble(j), e1, e2, e3, e4/)
                    EXIT
                end if
            end do
        end do
        
        ! Test Output
        !open(1, file= newdir//'/'//'Area_Coefficients.txt',form='formatted',status='unknown')
        !write(1,'(6f22.15)') transpose(AreaCoeff)
        !close(1)

    end subroutine getAreaCoefficients_3D
    
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
                if (e1 > -10e-10 .and. e2 > -10e-10 .and. e3 > -10e-10) then
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
                if (e1 > -10e-10 .and. e2 > -10e-10 .and. e3 > -10e-10) then
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
        !NoElem = size(DelaunayElem, dim = 1)
        j = 1
        do while (DelaunayElem(j,1) /= 0)        
            xa = DelaunayCoord(DelaunayElem(j,1),1)
            ya = DelaunayCoord(DelaunayElem(j,1),2)
            xb = DelaunayCoord(DelaunayElem(j,2),1)
            yb = DelaunayCoord(DelaunayElem(j,2),2)
            xc = DelaunayCoord(DelaunayElem(j,3),1)
            yc = DelaunayCoord(DelaunayElem(j,3),2)
            S = 0.5*((xa*yb + ya*xc + xb*yc) - (xa*yc + ya*xb + yb*xc))
            j = j + 1
            if (S < 0) then  ! If Area negative, the order of the points has changed (DelaunayElement order) and hence an intersection of elements took place
                print *, 'Intersection identified. Movement will be split up.'
                !call writeDatfile(10)
                !open(1, file= 'Delaunaynodes.dat', form='formatted', status = 'unknown')
                !write(1,'(2f22.15)') transpose(DelaunayCoord)
                !close(1)
                !open(1, file= 'DelaunayElem.dat', form='formatted', status = 'unknown')
                !write(1,'(3I8)') transpose(DelaunayElem)
                !close(1)
                intersect = 0
                EXIT
            end if
        end do
    
    end subroutine CheckforIntersections
    
    subroutine CheckforIntersections_3D(DelaunayCoord, DelaunayElem, intersect)
    
        ! Variables
        implicit none
        double precision :: xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd, S
        integer :: j, NoElem, intersect
        double precision, dimension(:,:) :: DelaunayCoord
        double precision, dimension(4,4) :: A
        integer, dimension(:,:) :: DelaunayElem
    
        ! Body of CheckforIntersections
        intersect = 1
        NoElem = size(DelaunayElem, dim = 1)
        do j = 1, NoElem
            xa = DelaunayCoord(DelaunayElem(j,1),1)
            ya = DelaunayCoord(DelaunayElem(j,1),2)
            za = DelaunayCoord(DelaunayElem(j,1),3)
            xb = DelaunayCoord(DelaunayElem(j,2),1)
            yb = DelaunayCoord(DelaunayElem(j,2),2)
            zb = DelaunayCoord(DelaunayElem(j,2),3)
            xc = DelaunayCoord(DelaunayElem(j,3),1)
            yc = DelaunayCoord(DelaunayElem(j,3),2)
            zc = DelaunayCoord(DelaunayElem(j,3),3)
            A(1,:) = (/ xa, ya, za, 1/)
            A(2,:) = (/ xb, yb, zb, 1/)
            A(3,:) = (/ xc, yc, zc, 1/)
            A(4,:) = (/ xd, yd, zd, 1/)
            S = (1/6)*Det(A,4)
            if (S < 0) then  ! If Area negative, the order of the points has changed (DelaunayElement order) and hence an intersection of elements took place
                print *, 'Intersection identified. Movement will be split up.'
                intersect = 0
                EXIT
            end if
        end do
    
    end subroutine CheckforIntersections_3D
    
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
    
    subroutine RelocateMeshPoints_3D(DelaunayCoord, DelaunayElem, AreaCoeff, NoP)
    
        ! Variables
        implicit none
        integer :: i, NoP
        double precision :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, xp, yp, zp, inp
        double precision, dimension(NoP,6) :: AreaCoeff
        double precision, dimension(:,:) :: DelaunayCoord
        integer, dimension(:,:) :: DelaunayElem

        ! Body of RelocateMeshPoints        
        do i = 1, NoP
            x1 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),1),1)
            y1 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),1),2)
            z1 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),1),3)
            x2 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),2),1)
            y2 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),2),2)
            z2 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),2),3)
            x3 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),3),1)
            y3 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),3),2)
            z3 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),3),3)
            x4 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),4),1)
            y4 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),4),2)
            z4 = DelaunayCoord(DelaunayElem(AreaCoeff(i,2),4),3)
            xp = x1*AreaCoeff(i,3) + x2*AreaCoeff(i,4) + x3*AreaCoeff(i,5) + x4*AreaCoeff(i,6)
            yp = y1*AreaCoeff(i,3) + y2*AreaCoeff(i,4) + y3*AreaCoeff(i,5) + y4*AreaCoeff(i,6)
            zp = z1*AreaCoeff(i,3) + z2*AreaCoeff(i,4) + z3*AreaCoeff(i,5) + z4*AreaCoeff(i,6)
            RD%coord_temp(AreaCoeff(i,1),:) = (/xp, yp, zp/)
        end do    
  
    end subroutine RelocateMeshPoints_3D
    
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
        integer :: i, j, k, l, nobp, testx, testy, circ, lin, nibp
        double precision :: maxx, maxy, minx, miny, spacing
        double precision, dimension(:), allocatable :: dist
        integer, dimension(:), allocatable ::  OuterBound, nodesvec, nodesvec2,nodesvec3, nodesvec4, distindex
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
        
        ! Identify overlapping nodes between moving and non-moving parts to input into Delaunay Coordinates (to ensure smooth interfaces)
        nodesvec = 0
        nodesvec2 = 0
        overlap = 0
        do i = 1, size(InnerBound)
            do j = 1, size(OuterBound)
                if (InnerBound(i) == OuterBound(j)) then
                    do k = 1, IV%NoCN           ! Check, that identified Coordinate is not already picked as CN
                        if (InnerBound(i) == InnerBound(CN_ind(k))) then
                            testx = 2
                            EXIT
                        end if
                    end do
                    if (testx /= 2) then
                        overlap = overlap + 1
                        nodesvec(overlap) = InnerBound(i)
                    end if
                end if
            end do        
        end do
        allocate(DelaunayCoordBound(IV%NoCN + IV%NoDelBP + overlap,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound"
        nbp = nibp + nobp - overlap ! Total number of boundary points
        
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
                if (abs(IdealCoord(j,1) - maxx) < 10e-8) then
                    EXit
                end if
            end do
            spacing = (maxy - miny)/(IV%NoDelBP/4)
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/maxx, (maxy - i*spacing)/)
                if (abs(IdealCoord(j,2) - miny) < 10e-8) then
                    EXit
                end if
            end do
            spacing = (maxx - minx)/(IV%NoDelBP/4)
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/(maxx - i*spacing), miny/)
                if (abs(IdealCoord(j,1) - minx) < 10e-8) then
                    EXit
                end if
            end do
            spacing = (maxy - miny)/(IV%NoDelBP/4)
            do i = 1, IV%NoDelBP
                j = j + 1
                IdealCoord(j,:) = (/minx, (miny + i*spacing)/)
                if (abs(IdealCoord(j,2) - maxy) < 10e-8) then
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
    
    subroutine getDelaunayCoordBound_3D()
    
        ! Variables
        implicit none
        integer :: i, j, k, l, nobp, testx, testy, circ, lin, nibp
        double precision :: maxx, maxy, minx, miny, spacing, maxz, minz
        double precision, dimension(:), allocatable :: dist
        integer, dimension(:), allocatable ::  OuterBound, nodesvec, nodesvec2,nodesvec3, nodesvec4, distindex
        real, PARAMETER :: Pi = 3.1415927
        double precision, dimension(:,:), allocatable :: IdealCoord
    
        ! Body of getDelaunayCoordBound
        
        ! Number of Delaunay Boundary Points                       
        allocate(nodesvec(3*RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        allocate(nodesvec2(3*RD%nbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        
        ! Separate internal and external boundary
        k = 0
        l = 0
        do i = 1, RD%nbf
            if (RD%boundtype(i) == 2 .or. RD%boundtype(i) == 3 .or. RD%boundtype(i) == 4) then
                k = k + 1
                nodesvec(k) = RD%bound(i,1)
                k = k + 1
                nodesvec(k) = RD%bound(i,2)
                k = k + 1
                nodesvec(k) = RD%bound(i,3)
            else
                l = l + 1
                nodesvec2(l) = RD%bound(i,1)
                l = l + 1
                nodesvec2(l) = RD%bound(i,2)
                l = l + 1
                nodesvec2(l) = RD%bound(i,3)
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
                dist(j) = DistP2P(IV%NoDim, RD%Coord_CN(i,1), RD%Coord(InnerBound(j),1), RD%Coord_CN(i, 2), RD%Coord(InnerBound(j),2), RD%Coord_CN(i, 3), RD%Coord(InnerBound(j),3))  ! Calculate Distances          
            end do
            CN_ind(i) = minloc(dist,dim=1) ! Returns the index(Position of Node in Coord Matrix) of the minimum Value
        end do
        deallocate(dist)
        
        ! Identify overlapping nodes between moving and non-moving parts to input into Delaunay Coordinates (to ensure smooth interfaces)
        ! There will be a range of nodes per overlap. How to reduce to sufficient amount? Otherwise system should work as it is.
        ! Find adjacent nodes -> Create array per "intercept plane" -> pick 2-4 nodes with max distant between them per plane
        nodesvec = 0
        nodesvec2 = 0
        overlap = 0
        do i = 1, size(InnerBound)
            do j = 1, size(OuterBound)
                if (InnerBound(i) == OuterBound(j)) then
                    do k = 1, IV%NoCN           ! Check, that identified Coordinate is not already picked as CN
                        if (InnerBound(i) == InnerBound(CN_ind(k))) then
                            testx = 2
                            EXIT
                        end if
                    end do
                    if (testx /= 2) then
                        overlap = overlap + 1
                        nodesvec(overlap) = InnerBound(i)
                    end if
                end if
            end do        
        end do
        allocate(DelaunayCoordBound(IV%NoCN + IV%NoDelBP + overlap,IV%NoDim),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in getDelaunayCoordBound"
        nbp = nibp + nobp - overlap ! Total number of boundary points
                       
        ! Integrate CN coordinates into Delaunay Coordinates required for Triangulation
        DelaunayCoordBound(1:IV%NoCN,:) = RD%Coord(InnerBound(CN_ind),:)
        ! Input overlapping nodes into Delaunay Coordinates
        if (overlap /= 0) then
            DelaunayCoordBound((IV%NoCN + 1):(IV%NoCN + overlap),:) = RD%Coord(nodesvec(1:overlap),:)
        end if
        deallocate(nodesvec)
        deallocate(nodesvec2)
    
        ! Find Points for Boundary Delaunay Coordinates
        maxz = maxval(RD%Coord(OuterBound,3))
        minz = minval(RD%Coord(OuterBound,3))
        maxy = maxval(RD%Coord(OuterBound,2))
        miny = minval(RD%Coord(OuterBound,2))
        minx = minval(RD%Coord(OuterBound,1))
        maxx = maxval(RD%Coord(OuterBound,1))
! Only corner coordinates allowed for the moment. Sensible to introduce possibility to choose more DN?
        DelaunayCoordBound((IV%NoCN + overlap + 1),:) = (/minx, miny, minz/)
        DelaunayCoordBound((IV%NoCN + overlap + 2),:) = (/minx, miny, maxz/)
        DelaunayCoordBound((IV%NoCN + overlap + 3),:) = (/minx, maxy, minz/)
        DelaunayCoordBound((IV%NoCN + overlap + 4),:) = (/minx, maxy, maxz/)
        DelaunayCoordBound((IV%NoCN + overlap + 5),:) = (/maxx, miny, minz/)
        DelaunayCoordBound((IV%NoCN + overlap + 6),:) = (/maxx, miny, maxz/)
        DelaunayCoordBound((IV%NoCN + overlap + 7),:) = (/maxx, maxy, minz/)
        DelaunayCoordBound((IV%NoCN + overlap + 8),:) = (/maxx, maxy, maxz/)
		
    end subroutine getDelaunayCoordBound_3D
    
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
        integer, dimension(:), allocatable ::  BoundIndex, nodesvec
    
        ! Body of getBoundaryIndex
        if (IV%NoDim == 3) then
            allocate(nodesvec(3*RD%nbf),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
            nodesvec = (/RD%bound(:,1), RD%bound(:,2), RD%bound(:,3)/)
        else
            allocate(nodesvec(2*RD%nbf),stat=allocateStatus)
            if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
            nodesvec = (/RD%bound(:,1), RD%bound(:,2)/)
        end if
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
    
    subroutine AngleofAttack(alpha, CNindex)
    
        ! Variables
        implicit none
        integer :: CNindex, nibp, i
        double precision, dimension(2,2) :: T
        double precision, dimension(2) :: P1, P2, P21, Pnew
        double precision :: alpha, ralpha
        real, PARAMETER :: Pi = 3.1415927

        ! Body of AngleofAttack
        nibp = size(InnerBound, dim = 1)
        if (allocated(dRot) == .true.) then
            deallocate(dRot)
        end if
        allocate(dRot(nibp, 2),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in MoveMesh "
        ralpha = alpha*Pi/180
        T(1,:) = (/cos(ralpha), sin(ralpha)/)
        T(2,:) = (/-sin(ralpha), cos(ralpha)/)
        P1 = RD%Coord(orderedBoundaryIndex(CN_indordered(CNindex)),:)
        do i = 1, nibp
            P2 = RD%Coord_temp(orderedBoundaryIndex(i),:)
            P21 = P2 - P1
            Pnew = matmul( T, P21) + P1
            dRot(orderedBoundaryIndex(i),:) = Pnew - P2
            RD%Coord_temp(orderedBoundaryIndex(i),:) = Pnew
        end do
         
    end subroutine AngleofAttack
    
    subroutine OrderBoundary()
    
        ! Variables
        implicit none
        integer :: indOrder, indjump, i, nmbf, j
        double precision, dimension(1000,3) :: jumped
        double precision, dimension(:), allocatable :: dist
    
        ! Body of OrderBoundary
        nmbf = size(InnerBound)
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
        !start = j+1
        !add = j
! DO WHILE
        do i = j+1, nmbf+j    ! TO-DO - introduce new variable to replace function of j
            call RecOrderBoundary(i, indOrder, jumped, indjump, orderedBoundaryIndex)
        end do       
        ! Check for all values being non-zero
        !test = .true.
        !do j = 1, nmbf
        !    if (orderedBoundaryIndex(j) == 0) then  ! If zero, means that
        !        test = .false.
        !        do k = 1, nmbf
        !            if (RD%boundtype(jumped(k,3)) /= 3 .and. RD%boundtype(jumped(k,3)) /= 4) then
        !                EXIT
        !            end if 
        !        end do
        !        start = start + j - 1
        !        add = add + k
        !        orderedBoundaryIndex(indOrder+1) = 0
        !        orderedBoundaryIndex(indOrder+2) = jumped(k,1)
        !        orderedBoundaryIndex(indOrder+3) = jumped(k,2)
        !        indOrder = indOrder + 3
        !        jumped(k,:) = jumped(indjump,:)
        !        indjump = indjump - 1
        !        EXIT
        !    end if
        !end do

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
    
    subroutine OrderBoundary_3D()
    
        ! Variables
        implicit none
        integer :: indOrder, indjump, i, nmbf, j
        double precision, dimension(1000,3) :: jumped
        double precision, dimension(:), allocatable :: dist
    
        ! Body of OrderBoundary
        nmbf = size(InnerBound)
        allocate(orderedBoundaryIndex(nmbf),stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in ReadData "
        do j = 1, nmbf
            if (RD%boundtype(j) /= 3 .and. RD%boundtype(j) /= 4) then
                EXIT
            end if 
        end do
        orderedBoundaryIndex(1) = RD%bound(j,1)
        orderedBoundaryIndex(2) = RD%bound(j,2)
        orderedBoundaryIndex(3) = RD%bound(j,3) 
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
                dist(j) = DistP2P(IV%NoDim, RD%Coord_CN(i,1), RD%Coord(orderedBoundaryIndex(j),1), RD%Coord_CN(i, 2), RD%Coord(orderedBoundaryIndex(j),2), RD%Coord_CN(i,3), RD%Coord(orderedBoundaryIndex(j),3))  ! Calculate Distances          
            end do
            CN_indordered(i) = minloc(dist,dim=1) ! Returns the index(Position of Node in Coord Matrix) of the minimum Value
        end do
        deallocate(dist)
        
    end subroutine OrderBoundary_3D
    
    recursive subroutine RecOrderBoundary(i, indOrder, jumped, indjump, oBI)
 ! How to identify correct order? What order is necessary for the smoothing?   
        ! Variables
        implicit none
        integer :: indOrder, indjump, i, j
        logical :: match
        double precision, dimension(1000,3) :: jumped
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
                jumped(indjump,:) = (/RD%bound(i,:), RD%boundtype(i)/)
            else
                call RecOrderBoundary(i, indOrder, jumped, indjump, oBI)    
            end if
        end if
            
    end subroutine RecOrderBoundary
    
    subroutine getDelaunayElem2(DelaunayElem, rows, columns, DelaunayCoord)
    ! Objective: Do the Triangulation and get the Delaunay Element Matrix
    
        ! Variables
        implicit none
        integer :: i, FileSize, NoElem, temp, rows, columns
        double precision :: xa, ya, xb, yb, xc, yc, S
        integer, dimension(rows,columns) :: DelaunayElem
        double precision, dimension(:,:), allocatable :: DelaunayCoord
   
        ! Body of getDelaunayElem
        open(1, file= 'Delaunay_nodes.txt', form='formatted', status = 'unknown')
        write(1,'(2f22.15)') transpose(DelaunayCoord)
        close(1, iostat = allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"

        open(2, file = 'DelaunayInput.txt', form='formatted', status = 'unknown')
        write(2,*) 'Delaunay'
        close(2)
 
        ! Perform Delaunay Triangulation via external executable
        allocate(character(len=100) :: strSystem,stat=allocateStatus)
        if(allocateStatus/=0) STOP "ERROR: Not enough memory in FDGD"
        if (IV%systemType == 'B' .or. IV%systemtype == 'Q') then
            strSystem = trim(IV%filepath)//'/Executables/Delaunay_2D < DelaunayInput.txt > /dev/null'
        elseif (IV%systemType == 'W') then
            strSystem = 'Executables\Delaunay_2D.exe < DelaunayInput.txt >nul 2>&1'
        elseif (IV%systemType == 'L') then
            strSystem = 'Executables/Delaunay_2D < DelaunayInput.txt > /dev/null'
        else
            STOP 'INPUT ERROR: System Type selected does not exist! Program stopped.'
        end if
        call system(trim(strSystem))
        deallocate(strSystem)
   
        open(3, file= 'Delaunay_elements.txt', form='formatted',status='old')
        inquire(3, size = FileSize)
        if (IV%SystemType == 'W') then          
            NoElem = FileSize/32
        else
            NoElem = FileSize/31
        end if
        DelaunayElem = 0
        do i = 1, NoElem
            read(3, *) DelaunayElem(i,:)
        end do        
        close(3)
          
        do i = 1, NoElem
            xa = DelaunayCoord(DelaunayElem(i,1),1)
            ya = DelaunayCoord(DelaunayElem(i,1),2)
            xb = DelaunayCoord(DelaunayElem(i,2),1)
            yb = DelaunayCoord(DelaunayElem(i,2),2)
            xc = DelaunayCoord(DelaunayElem(i,3),1)
            yc = DelaunayCoord(DelaunayElem(i,3),2)
            S = 0.5*((xa*yb + ya*xc + xb*yc) - (xa*yc + ya*xb + yb*xc))
            if (S <= 0) then  ! If Area negative, the order of the nodes is incorrect and positions will be swapped
                print *, 'Area negative when generating Delaunay Background mesh'
                temp = DelaunayElem(i,1)
                DelaunayElem(i,1) = DelaunayElem(i,3)
                DelaunayElem(i,3) = temp
            end if
        end do
        
    end subroutine getDelaunayElem2
    
end module FDGD