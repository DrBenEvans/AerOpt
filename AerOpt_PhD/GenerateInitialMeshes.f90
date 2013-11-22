module GenerateInitialMeshes
    
        use Toolbox
        integer :: ne, np, nbf            ! Number of elements, Nodes/Points & boundary faces
        real, dimension(:,:), allocatable :: coordf ! Coordinates Matrix of Fine Mesh
        integer, dimension(:,:), allocatable :: connecf, boundf ! Connectivity & Boundary Matrix of Fine Mesh
        real, dimension(:,:), allocatable :: Coord_CP
        real, dimension(:,:), allocatable :: Rect
            
contains
        
    subroutine SubGenerateInitialMeshes(NoDim, NoCP, coordf, connecf, boundf, Coord_CP, Rect)
        
        ! Variables
        real, dimension(np,NoDim) :: coordf
        integer, dimension(ne,NoDim+1) :: connecf   
        integer, dimension(nbf,NoDim) :: boundf
        real, dimension(NoCP, NoDim) :: Coord_CP
        real, dimension(nbf) :: dCP2N
        integer, dimension(NoCP) :: mind, av
        real, dimension(NoCP,NoDim*4) :: Rect
        integer, dimension(:, :), allocatable :: InfluenceBox
        integer, dimension(:), allocatable :: IB1, IB2, IB3, IB4, IB5, IB6, IB7
        integer :: lb

        ! Body of GenerateInitialMeshes
        
        ! Find Real Control Nodes based on Coordinates - on long terms redundant only required to run once           
        ! Also Identify boundary nodes within the influence box of each Control Point
        allocate(InfluenceBox(nbf, NoCP)) 
        do i = 1, NoCP
            k = 0
            do j = 1, nbf
                dCP2N(j) = DistP2P(NoDim, Coord_CP(i,1), coordf(j,1), Coord_CP(i, 2), coordf(j,2))  ! Calculate Distances
                
                if (Rectcheck(Rect(i,:), coordf(j,:)) == 1) then ! Check if boundary node is in influence box - Note: first values in coordinates matrix are per default the boundary nodes
                    k = k + 1
                    InfluenceBox(k,i) = j                    
                end if 
            end do
            mind(i) = minloc(dCP2N,dim=1) ! Returns the index(Position of Node in Coord Matrix) of the minimum Value          
            av(i) = k
        end do       
        !print *, mind
        
        ! Seperate Array to Influence Box of each Control Point
        allocate(IB1(av(1)));   IB1 = InfluenceBox(1:av(1),1)
        allocate(IB2(av(2)));   IB2 = InfluenceBox(1:av(2),2)
        allocate(IB3(av(3)));   IB3 = InfluenceBox(1:av(3),3)
        allocate(IB4(av(4)));   IB4 = InfluenceBox(1:av(4),4)
        allocate(IB5(av(5)));   IB5 = InfluenceBox(1:av(5),5)
        allocate(IB6(av(6)));   IB6 = InfluenceBox(1:av(6),6)
        allocate(IB7(av(7)));   IB7 = InfluenceBox(1:av(7),7)
        deallocate(InfluenceBox)
        !print *, IB1
        !print *, IB2
        
        
    end subroutine SubGenerateInitialMeshes
        
    subroutine ReadData(NoCP, NoDim)
    !! Objective: Reads the Mesh data and Control Points Coordinates
    
        ! Body of ReadData
        open(1, file='C:\Users\EG717761\Documents\Visual Studio 2010\Projects\AerOpt_PhD\Input_Data\Mesh_fine.txt')
        read(1, 1) ne
        allocate(connecf(ne,NoDim+1))
        read(1, 1) np
        allocate(coordf(np,NoDim))
        read(1, 1) nbf
        allocate(boundf(nbf,NoDim))
        do i = 1, ne
            read(1, *) connecf(i,:)
        end do
        do i = 1, np
            read(1, *) coordf(i,:)
        end do
        do i = 1, nbf
            read(1, *) boundf(i,:)
        end do
    1   format(1i9)   
        close(1)
    
        allocate(Coord_CP(NoCP,NoDim))
        open(2, file='C:\Users\EG717761\Documents\Visual Studio 2010\Projects\AerOpt_PhD\Input_Data\Control_Nodes.txt')
        do i = 1, NoCP
            read(2, *) Coord_CP(i,:)
        end do
        close(2)
        
        allocate(Rect(NoCP,NoDim*4))
        open(3, file='C:\Users\EG717761\Documents\Visual Studio 2010\Projects\AerOpt_PhD\Input_Data\Rectangles.txt')
        do i = 1, NoCP
            read(3, *) Rect(i,:)
        end do
        close(3)
    
    end subroutine ReadData
    
end module GenerateInitialMeshes