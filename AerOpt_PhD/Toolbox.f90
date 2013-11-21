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
        
        subroutine MsgBox(ShowText, TitleCaption)
        !! Objective: Output a MsgBox to better inform the User. Also breaks the program.
        
            ! Variables
            use ifwin 
            integer(SINT) :: ret
            character (len=60) :: ShowText
            character (len=15) :: TitleCaption  
             
            ! Body of MsgBox
            ret = MessageBox ( & GetForegroundWindow(), &   ! Handle to window 
            ShowText, &                                     ! Text (don't forget C-string) 
            TitleCaption, &                                 ! Caption for title bar 
            MB_ICONINFORMATION + MB_OK)                     ! Type flags
            
        end subroutine MsgBox
        
        function DistP2P (NoDim, Xa, Xb, Ya, Yb, Za, Zb)
        
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
        
        !subroutine RectCheck
        !    
        !
        !end subroutine RectCheck
    end module Toolbox