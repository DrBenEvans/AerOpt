module GenerateInitialMeshes
    
        use Toolbox
        integer :: ne, np, nbf            ! Number of elements, Nodes/Points & boundary faces
        real, dimension(:,:), allocatable :: coordf ! Coordinates Matrix of Fine Mesh
        integer, dimension(:,:), allocatable :: connecf, boundf ! Connectivity & Boundary Matrix of Fine Mesh
        real, dimension(:,:), allocatable :: Coord_CP
            
contains
        
    subroutine SubGenerateInitialMeshes(NoDim, NoCP, coordf, connecf, boundf, Coord_CP)
        
        ! Variables
        real, dimension(np,NoDim) :: coordf
        integer, dimension(ne,NoDim+1) :: connecf   
        integer, dimension(nbf,NoDim) :: boundf
        real, dimension(NoCP, NoDim) :: Coord_CP
        real, dimension(nbf) :: dCP2N
        integer, dimension(NoCP) :: mind
        integer :: lb

        ! Body of GenerateInitialMeshes - long terms - redundant only required to run once           
        do i = 1, NoCP
            do j = 1, nbf
                dCP2N(j) = DistP2P(NoDim, Coord_CP(i,1), coordf(j,1), Coord_CP(i, 2), coordf(j,2))
            end do
            mind(i) = minloc(dCP2N,dim=1) ! Returns the index(Position of Node in Coord Matrix) of the minimum Value          
        end do
        
        print *, mind
        !!!!! DEALLOCATE big variables at the end!!!!
        
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
    
    end subroutine ReadData
    
end module GenerateInitialMeshes