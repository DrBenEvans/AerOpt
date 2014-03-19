module Toolbox
            
contains
    
    function linSpacing(max, min, NoInt)
    !! Objective: Generate a uniform, linear Distribution of NoInt values in the (max, min) Domain
        
        ! Variables
        integer :: a, NoInt
        double precision :: c, max, min
        double precision, dimension(NoInt) :: linSpacing
    
        ! Body of linSpacing
        if (NoInt == 1) then
            linSpacing(1) = (max - min)/2
        else
            do a = 0, (NoInt - 1)
                c = a/ real(NoInt - 1)  
                linSpacing(a+1) = max - (max - min)*c
            end do
        end if
    
    end function linSpacing
        
    function DistP2P (NoDim, Xa, Xb, Ya, Yb, Za, Zb)
    ! Objective: Calculate the Distance in any Dimension
        
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
    ! Objective: Check of values to be positioned in a Rectangle
        
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

    Subroutine PRINT_MATRIX( DESC, M, N, A, LDA )
    ! Objective: Print a Matrix
        
        ! Variables
        character*(*) ::   DESC
        integer ::         M, N, LDA
        double precision :: A( LDA, * )
        integer ::        I, J
            
        ! Body of Print Matrix
        write(*,*)
        write(*,*) DESC
        do I = 1, M
            write(*,9998) ( A( I, J ), J = 1, N )
        end do

        9998 format( 11(:,1X,F6.2) )
        return
            
    End Subroutine PRINT_MATRIX
        
    recursive subroutine QSort(a,na, sel, ind)
    ! Objective: Sort/Order an Array in ascending order
    
        ! Dummy Arguments
        integer, intent(in) :: nA
        real, dimension(nA), intent(in out) :: A
        integer, dimension(nA), intent(in out), optional :: ind
        character(len=1), optional :: sel
 
        ! Local Variables
        integer :: left, right
        real :: random
        real :: pivot
        integer :: marker
        real :: temp
 
        if (nA > 1) then
 
            call random_number(random)
            pivot = A(int(random*real(nA-1))+1)   ! random pivor (not best performance, but avoids worst-case)
            left = 0
            right = nA + 1
 
            do while (left < right)
                right = right - 1
                do while (A(right) > pivot)
                    right = right - 1
                end do
                left = left + 1
                do while (A(left) < pivot)
                    left = left + 1
                end do
                
                ! Swap numbers
                if (left < right) then
                    temp = A(left)
                    A(left) = A(right)
                    A(right) = temp
                    if (sel == 'y') then
                        temp = ind(left)
                        ind(left) = ind(right)
                        ind(right) = temp
                    end if
                end if
            end do
 
            if (left == right) then
                marker = left + 1
            else
                marker = left
            end if
 
            call QSort(A(:marker-1),marker-1, sel, ind(:marker-1)) !recursive call -1
            call QSort(A(marker:),nA-marker+1, sel, ind(marker:)) ! recursive call +1
 
        end if
 
    end subroutine QSort

    subroutine Unique(vector_in, sz, vector_out)
    
        ! Variables
        integer :: sz
        real, dimension(sz), intent(in) :: vector_in
        real, dimension(:), allocatable, intent(out) :: vector_out
        integer, dimension(sz) :: marker
    
        ! Body of Unique
        j = 1
        marker(1) = 0
        do i = 2, sz
            if (vector_in(i) == vector_in(i-1)) then
                marker(i) = 1
            else
                marker(i) = 0
                j = j + 1
            end if   
        end do
        
        allocate(vector_out(j))
        j = 1
        do i = 1, sz            
            if (marker(i) == 0) then
                vector_out(j) = vector_in(i)
                j = j + 1
            end if
        end do
    
    end subroutine Unique
    
    subroutine randperm(N, p)
    
    !! Source: coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2006-03/msg00748.html
    !! Based on Knuth's algorithm

    integer, intent(in) :: N
    integer, dimension(:), intent(out) :: p

    integer :: temp

    p = (/ (i, i=1,N) /)

    do j=N,2,-1

        call random_number(rn)
        k = floor(j*rn) + 1

        ! exchange p(k) and p(j)
        temp = p(k)
        p(k) = p(j)
        p(j) = temp

    end do

    end subroutine randperm
    
    subroutine inverse(a,c,n)
    
        !============================================================
        ! Inverse matrix
        ! Method: Based on Doolittle LU factorization for Ax=b
        ! Alex G. December 2009
        ! http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
        !-----------------------------------------------------------
        ! input ...
        ! a(n,n) - array of coefficients for matrix A
        ! n      - dimension
        ! output ...
        ! c(n,n) - inverse matrix of A
        ! comments ...
        ! the original matrix a(n,n) will be destroyed 
        ! during the calculation
        !===========================================================
        implicit none 
        integer n
        double precision :: a(n,n)
        double precision, intent(out) :: c(n,n)
        double precision L(n,n), U(n,n), b(n), d(n), x(n)
        double precision coeff
        integer i, j, k

        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices
        L=0.0
        U=0.0
        b=0.0

        ! step 1: forward elimination
        do k=1, n-1
           do i=k+1,n
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                 a(i,j) = a(i,j)-coeff*a(k,j)
              end do
           end do
        end do

        ! Step 2: prepare L and U matrices 
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i=1,n
          L(i,i) = 1.0
        end do
        ! U matrix is the upper triangular part of A
        do j=1,n
          do i=1,j
            U(i,j) = a(i,j)
          end do
        end do

        ! Step 3: compute columns of the inverse matrix C
        do k=1,n
          b(k)=1.0
          d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
            d(i)=b(i)
            do j=1,i-1
              d(i) = d(i) - L(i,j)*d(j)
            end do
          end do
        ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
              x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
          end do
        ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
            c(i,k) = x(i)
          end do
          b(k)=0.0
        end do
        
    end subroutine inverse

    subroutine DetermineStrLen(istr, i)
    
        ! Variables
        character(len=:), allocatable :: istr
    
        ! Body of DetermineStrLen     
        if (i < 10) then
            allocate(character(len=1) :: istr)
            write( istr, '(I1)' )  i
        elseif (i < 100) then
            allocate(character(len=2) :: istr)
            write( istr, '(I2)' )  i
        elseif (i < 1000) then
            allocate(character(len=3) :: istr)
            write( istr, '(I3)' )  i
        else
            allocate(character(len=4) :: istr)
            write( istr, '(I4)' )  i
        end if
        
    end subroutine DetermineStrLen
    
    subroutine communicateWin2Lin(Username, Password, fname, channel)
    
        ! Variables
        integer :: IntSystem
        character(len=*) :: Username, Password, fname, channel
        character(len=:), allocatable :: strSyst
    
        ! Body of communicateWin2Lin
        if (channel == 'psftp') then
            IntSystem = 74 + len(Username) + len(Password) + len(fname)
            allocate(character(len=IntSystem) :: strSyst)
            strSyst = '"C:\Program Files (x86)\WinSCP\PuTTY\psftp" '//UserName//'@encluster.swan.ac.uk -pw '//Password//' -b '//fname
        else
            IntSystem = 79 + len(Username) + len(Password) + len(fname)
            allocate(character(len=IntSystem) :: strSyst)
            strSyst = '"C:\Program Files (x86)\WinSCP\PuTTY\'//channel//'" -ssh '//UserName//'@encluster.swan.ac.uk -pw '//Password//' -m '//fname
        end if
        call system(strSyst)
            
    end subroutine communicateWin2Lin      

end module Toolbox