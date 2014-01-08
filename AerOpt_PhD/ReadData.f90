module ReadData
    
    integer :: ne, np, nbf, nbc                             ! Number of elements, Nodes/Points & boundary faces
    integer, dimension(:,:), allocatable :: connecc         ! Connectivity Matrix of Coarse Mesh    
    integer, dimension(:,:), allocatable :: connecf, boundf ! Connectivity & Boundary Matrix of Fine Mesh    
    real, dimension(:,:), allocatable :: coord              ! Coordinates Matrix of Fine Mesh (includes coordinates of coarse mesh)
    real, dimension(:,:), allocatable :: coarse             ! includes element allocation of nodes to coarse triangles and Area Coefficients of each node
    real, dimension(:,:), allocatable :: Coord_CP           ! desired Coordinates of the Control Points
    real, dimension(:,:), allocatable :: Rect               ! Rectangle definition of 'Influence Box'
    
contains
      
    subroutine SubReadData(NoCP, NoDim)
    ! Objective: Reads the Mesh data and Control Points Coordinates
    
    ! Body of ReadData
    open(1, file='Input_Data/Mesh_fine.txt')
    read(1, 11) ne
    allocate(connecf(ne,NoDim+1))
    read(1, 11) np
    allocate(coord(np,NoDim))
    read(1, 11) nbf
    allocate(boundf(nbf,NoDim))
    do i = 1, ne
        read(1, *) connecf(i,:)
    end do
    do i = 1, np
        read(1, *) coord(i,:)
    end do
    do i = 1, nbf
        read(1, *) boundf(i,:)
    end do 
11  format(1I8)        
    close(1)
    
    allocate(Coord_CP(NoCP,NoDim))
    open(2, file='Input_Data/Control_Nodes.txt')
    do i = 1, NoCP
        read(2, *) Coord_CP(i,:)
    end do
    close(2)
    
    allocate(Rect(NoCP,NoDim*4))
    open(3, file='Input_Data/Rectangles.txt')
    do i = 1, NoCP
        read(3, *) Rect(i,:)
    end do
    close(3)
    
    open(4, file='Input_Data/Mesh_coarse.txt')
    read(4, *) nbc
    allocate(connecc(nbc,NoDim+1))
    allocate(coarse(np-nbf, NoDim+2))
    do i = 1, (np - nbf)
        read(4, *) coarse(i,:)
    end do
    do i = 1, nbc
        read(4, *) connecc(i,:)
    end do
    close(4)
    
    end subroutine SubReadData
end module ReadData