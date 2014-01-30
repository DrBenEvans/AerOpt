module GenerateInitialMeshes
    
        use Toolbox
        use ReadData
        real, dimension(:,:), allocatable :: coord_temp  ! Coordinates Matrix of Fine Mesh
            
contains
        
    subroutine SubGenerateInitialMeshes(NoDim, NoCP, coord_temp, connecf, boundf, coarse, connecc, Coord_CP, Rect, NestDisp)
        
        ! Variables
        real, dimension(np,NoDim) :: coord_temp
        real, dimension(NoCP, NoDim) :: Coord_CP
        real, dimension(NoCP,NoDim*4) :: Rect
        real, dimension((np-nbf), NoDim+2) :: coarse
        real, dimension(nbf) :: dCP2N
        real, dimension(NoCP*NoDim) :: NestDisp
        integer, dimension(nbf, NoCP) :: IB
        integer, dimension(ne,NoDim+1) :: connecf
        integer, dimension(nbc,NoDim+1) :: connecc
        integer, dimension(nbf,NoDim) :: boundf
        integer, dimension(NoCP) :: CP_ind, size_ib
        real :: c, w, dis, x1, y1, x2, y2, x3, y3, xp, yp
        integer :: lb

        ! Body of GenerateInitialMeshes
        
        ! Find Real Control Nodes based on Coordinates - !! on long terms redundant only required to run once           
        ! Also Identify boundary nodes within the influence box of each Control Point
        ! NOTE: Hardwired for a 2D problem!
        do i = 1, NoCP
            k = 0
            do j = 1, nbf
                dCP2N(j) = DistP2P(NoDim, Coord_CP(i,1), coord_temp(j,1), Coord_CP(i, 2), coord_temp(j,2))  ! Calculate Distances
                
                if (Rectcheck(Rect(i,:), coord_temp(j,:)) == 1) then ! Check if boundary node is in influence box - Note: first values in coordinates matrix are per default the boundary nodes
                    k = k + 1
                    if (i /= 6) then
                        IB(k,i) = j
                    elseif (k > 15) then
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
        c = 1.9
        do i = 1, NoCp
            do j = 1, size_ib(i)
                dis = coord_temp(CP_ind(i),1) - coord_temp(IB(j,i),1)   ! Distance of Control Node to a Node in the Influence Box         
                w = exp(-(dis**2)/(c**2))   ! Gaussian
                coord_temp(IB(j,i),:) = (/(coord_temp(IB(j,i),1) + w*NestDisp(i)),(coord_temp(IB(j,i),2) + w*NestDisp(NoCP+i))/)   ! New coordinates
            end do
        end do
        
!!!!!! PLOT initial Nests --> MATLAB output file
        
        ! Relocation of the internal nodes(non-boundary)
        ! Method: Maintain the Area coefficient
        do i = 1, np - nbf
            x1 = coord_temp(connecc(coarse(i,1),1),1)
            y1 = coord_temp(connecc(coarse(i,1),1),2)
            x2 = coord_temp(connecc(coarse(i,1),2),1)
            y2 = coord_temp(connecc(coarse(i,1),2),2)
            x3 = coord_temp(connecc(coarse(i,1),3),1)
            y3 = coord_temp(connecc(coarse(i,1),3),2)
            
            xp = x1*coarse(i,2) + x2*coarse(i,3) + x3*coarse(i,4)
            yp = y1*coarse(i,2) + y2*coarse(i,3) + y3*coarse(i,4)
            coord_temp((i+nbf),:) = (/xp, yp/)
        end do
        
    end subroutine SubGenerateInitialMeshes
    
end module GenerateInitialMeshes