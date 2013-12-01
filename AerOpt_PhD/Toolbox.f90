    module Toolbox
    
        contains
    
        function linSpacing(max, min, NoInt)
        !! Objective: Generate a uniform, linear Distribution of NoInt values in the (max, min) Domain
        
            ! Variables
            integer :: a, NoInt
            double precision :: c, max, min
            double precision, dimension(NoInt) :: linSpacing
    
            ! Body of linSpacing
            do a = 0, (NoInt - 1)
                c = a/ real(NoInt - 1)  
                linSpacing(a+1) = max - (max - min)*c
            end do
    
        end function linSpacing
        
        function DistP2P (NoDim, Xa, Xb, Ya, Yb, Za, Zb)
        !! Objective: Calculate the Distance in any Dimension
        
            ! Variables
            optional :: Ya, Yb, Za, Zb
            real ::  DistP2P
               
            ! Body of DistP2P
            select case (NoDim)
            case (1)
                DistP2P = (Xa - Xb)
            case (2)
                DistP2P = sqrt((Xa - Xb)**2 + (Ya - Yb)**2)
            case (3)
                DistP2P = sqrt((Xa - Xb)**2 + (Ya - Yb)**2 + (Za - Zb)**2)
            end select
        
        end function DistP2P
        
        function RectCheck(r, p)
        !! Objective: Check of values to be positioned in a Rectangle
        
            ! Variables
            real, dimension(2) :: AB, BC, AP, BP, p
            real, dimension(8) :: r
            logical :: Rectcheck
            
            ! Body or RectCheck
            Rectcheck = 0
            AB = (/(r(3)-r(1)), (r(4) - r(2))/) ! B-A -> (x,y)
            BC = (/(r(5)-r(3)), (r(6) - r(4))/) ! C-B -> (x,y)
            AP = (/(p(1)-r(1)), (p(2) - r(2))/) ! P-A -> (x,y)
            BP = (/(p(1)-r(3)), (p(2) - r(4))/) ! P-B -> (x,y)
            
            if (0 <= dot_product(AB, AP) .and. &
                dot_product(AB, AP) <= dot_product(AB, AB) .and. &
                0 <= dot_product(BC, BP) .and. &
                dot_product(BC, BP) <= dot_product(BC, BC)) then
            Rectcheck = 1
            end if          
        
        end function RectCheck
        
    end module Toolbox