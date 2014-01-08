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


        !******************************************************************************
        !  Given a matrix A, with logical dimensions M by N and physical dimensions MP by NP, this
        !  routine computes its singular value decomposition, A = U.W.V'. The matrix U replaces 
        !  A on output. The diagonal matrix of singular values W is output as a vector W. The matrix 
        !  V (not the Transpose V') is output as V. M must be greater or equal to N; if it is smaller,
        !  then A should be filled up to square with zero rows.
        !******************************************************************************
        !  Verified: working on 6 Feb 1998
        !
        ! NOTES: 
        ! 1. When used to solve Linear Equations with n equations and n unknowns,
        !     mp=np=m=n.
        ! 2. When finding inverses, n = soln and b = IdentityMatrix(n)
        !
        ! Modifications to Original Program:
        ! 1. Equalties/Inequalities are replaced by Differences and compared with EPS
        ! 2. The Matrix U in U.W.V' is actually stored in "a"
        !
        ! DEBUG:
        ! 1.  "IMPLICIT NONE" and "REAL(DBP)"
        ! 2.  Parameters might need to be larger
        !******************************************************************************
        !SUBROUTINE SVDCMP(a, m, n, mp, nnp, V1)
        !
        !    ! Variables
        !    implicit none      
        !    integer, parameter :: dbp = SELECTED_REAL_KIND (15,307)
        !    real(dbp)  :: EPS
        !    parameter (EPS=3.0d-15)       	!EPS is the relative precision
        !
        !    integer :: m,n,mp,nnp, i,j,l,k, its, nm, jj
        !    real(dbp) :: g, sscale, anorm, s, f,h,x,z,y,c
        !
        !    integer, parameter :: nmax = 10000        !Maximum anticipated value of N
        !    real(dbp) :: A(mp,nnp), w(nnp), v1(nnp,nnp), rv1(nmax)
        !
        !    
        !    ! Body of SVDCMP
        !    PRINT *, 'Precision chosen as ',eps
        !
        !    if (m.lt.n) then
        !    PRINT *, 'You must augment A with extra zero rows'
        !    call exit(10)
        !    ENDIF 
        !
        !            !Householder Reduction to bidiagonal form
        !    !(see Forsythe,Malcolm,Moler, "Computer Methods for Mathematical Computations"
        !   g=0.0d0
        !   sscale = 0.0d0
        !   anorm = 0.0d0
        !   do i = 1,n
        !      l = i + 1
        !      rv1(i) = sscale*g
        !      g = 0.0d0
        !      s = 0.0d0
        !      sscale = 0.0d0
        !      if (i.le.m) then
        !         do k = i,m
        !            sscale = sscale + dABS(a(k,i))
        !         end do       ! k loop
        !!         if (sscale.ne.0.0d0) then
			     !   if ( dabs(sscale-0.0d0).gt.EPS ) then
        !            do k = i,m
        !               a(k,i) = a(k,i) / sscale
        !               s = s + a(k,i)*a(k,i)
        !            end do    ! k loop
        !            f = a(i,i)
        !            g = - SIGN(SQRT(s),f)
        !            h = f*g - s
        !            a(i,i) = f - g
        !            if (i.ne.n) then
        !               do j = l,n
        !                  s = 0.0d0
        !                  do k = i,m
        !                     s = s + a(k,i)*a(k,j)
        !                  end do      ! k loop
        !                  f = s / h
        !                  do k = i, m 
        !                     a(k,j) = a(k,j) + f*a(k,i)
        !                  end do   ! k loop
        !               end do      ! j loop
        !            end if
        !            do k = i, m 
        !               a(k,i) = sscale * a(k,i)
        !            end do         ! k loop
        !         end if
        !      end if
        !
        !      w(i) = sscale * g
        !      g = 0.0d0
        !      s = 0.0d0
        !      sscale = 0.0d0
        !      if ((i.le.m).AND.(i.ne.n)) then
        !         do k = l, n
        !            sscale = sscale + dABS(a(i,k))
        !         end do         ! k loop
        !!         if (sscale.ne.0.0d0) then
			     !   if ( dabs(sscale-0.0d0).gt.EPS ) then
        !            do k = l, n
        !               a(i,k) = a(i,k) /sscale
        !               s = s + a(i,k) * a(i,k)
        !            end do      ! k loop 
        !            f = a(i,l) 
        !            g = - SIGN(SQRT(s),f)
        !            h = f * g - s
        !            a(i,l) = f - g
        !            do k = l, n
        !               rv1(k) = a(i,k) / h
        !            end do      ! k loop
        !            if (i.ne.m) then
        !               do j = l, m 
        !                  s = 0.0d0
        !                  do k = l, n 
        !                     s = s + a(j,k)*a(i,k)
        !                  end do   ! k loop
        !                  do k = l, n 
        !                     a(j,k) = a(j,k) + s*rv1(k)
        !                  end do   ! k loop
        !               end do      ! j loop
        !            end if
				    !    do k = l, n
        !               a(i,k) = sscale * a(i,k)
        !   	        end do
        !         end if
        !      end if
        !      anorm = MAX(anorm, (dABS(w(i)) + dABS(rv1(i))))
        !   end do
        !
        !! Accumulation of right-hand Transformations
        !   do i = n, 1, -1 
        !      if (i.lt.n) then
        !!         if (g.ne.0.0d0) then
			     !   if ( dabs(g-0.0d0).gt.EPS ) then
        !            do j = l, n       ! Double division to avoid possible overflow
        !               v1(j,i) = (a(i,j) / a(i,l)) / g
        !            end do      ! j loop
        !            do j = l, n
        !               s = 0.0d0
        !               do k = l, n
        !                  s = s + a(i,k)*v1(k,j)
        !               end do   ! k loop
        !               do k = l, n
        !                  v1(k,j) = v1(k,j) + s * v1(k,i)
        !               end do   ! k loop
        !   	        end do      ! j loop
        !         end if
        !         do j = l, n 
        !            v1(i,j) = 0.0d0
        !            v1(j,i) = 0.0d0
        !         end do
        !      end if
        !      v1(i,i) = 1.0d0
        !      g = rv1(i)
        !      l = i
        !   end do
        !
        !! Accumulation of left-hand Transformations
        !   do i = n, 1, -1
        !      l = 1 + i
        !      g = w(i)
        !      if (i.lt.n) then
        !         do j = l, n
        !            a(i,j) = 0.0d0
        !         end do
        !      end if
        !!      if (g.ne.0.0d0) then
        !      if ( dabs(g-0.0d0).gt.EPS ) then
        !         g = 1.0d0 / g
        !         if (i.ne.n) then
        !            do j = l,n 
        !               s = 0.0d0
        !               do k = l, m
        !                  s = s + a(k,i)*a(k,j)
        !               end do   ! k loop
        !               f = (s/a(i,i)) * g
        !               do k = i, m 
        !                  a(k,j) = a(k,j) + f * a(k,i)
        !               end do   ! k loop
        !            end do      ! j loop
        !         end if
        !         do j = i, m 
        !            a(j,i) = a(j,i) * g
        !         end do         ! j loop
        !      else
        !         do j = i, m
        !            a(j,i) = 0.0d0
        !         end do         ! j loop
        !      end if
        !      a(i,i) = a(i,i) + 1.0d0
        !   end do               ! i loop
        !
        !! Diagonalization of the bidigonal form
        !   do k = n, 1, -1                  !Loop over singular values
        !      do its = 1,30                 !Loop over allowed iterations
        !         do l = k, 1, -1            !Test for splitting
        !            nm = l - 1              ! Note that rv1(1) is always zero
        !!           if ( (dABS(rv1(l))+anorm) .eq. anorm ) GO TO 2
        !!          	if ( (dABS(w(nm))+anorm) .eq. anorm ) GO TO 1
        !            if ( dabs((dABS(rv1(l))+anorm) - anorm).lt.eps ) GO TO 2
        !  	        if ( dabs((dABS(w(nm))+anorm) - anorm).lt.eps ) GO TO 1
        !         end do      !  l loop
        !
        !1        c = 0.0d0                  ! Cancellation of rv1(l), if l>1 :
        !         s = 1.0d0
        !         do i = l, k
        !            f = s * rv1(i)
        !!            if ( (dABS(f)+anorm) .ne. anorm ) then
        !            if ( dabs( (dABS(f)+anorm) - anorm) .GT. eps ) then
        !
        !               g = w(i)
        !               h = SQRT(f*f + g*g)
        !               w(i) = h
        !               h = 1.0d0 / h
        !               c = g * h
        !               s = -f * h
        !               do j = 1, m
        !                  y = a(j,nm)
        !                  z = a(j,i)
        !                  a(j,nm) = (y*c) + (z*s)
        !                  a(j,i) = -(y*s) + (z*c)
        !               end do   ! j loop
        !            end if
        !         end do         ! i loop
        !2        z = w(k) 
        !         if (l .eq. k) then         ! convergence
				    !    if (z .lt. 0.0d0) then  ! Singular value is made non-negative
        !               w(k) = -z
        !               do j = 1,n
        !                  v1(j,k) = -v1(j,k)
        !               end do         ! j loop
	    		 !       end if
        !            GO TO 3
        !         end if
        !         if (its.eq.30) then
        ! 	        PRINT*, 'No Convergence in 30 iterations'
        !            call exit(10)
        !         ENDIF
        !         x = w(l)          ! Shift from bottom 2-by-2 minor
        !         nm = k - 1
        !         y = w(nm)
        !         g = rv1(nm)
        !         h = rv1(k)
        !         f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( 2.0d0*h*y)
        !         g = SQRT(f*f + 1.0d0)
        !         f = ( (x-z)*(x+z) + h*((y/(f+SIGN(g,f))) - h) ) / x
        !
        !! Next   QR Transformation
        !         c = 1.0d0
        !         s = 1.0d0
        !         do j = l, nm
        !  	        i = j + 1
        !            g = rv1(i)
        !            y = w(i)
        !            h = s*g
        !            g = c*g
        !            z = SQRT(f*f + h*h)
        !            rv1(j) = z
        !            c = f/z
        !            s = h/z
        !            f = (x*c) + (g*s)
        !            g = -(x*s) + (g*c)
        !            h = y*s
        !            y = y*c
        !            do jj = 1, n
        !               x = v1(jj,j)
        !               z = v1(jj,i)
        !               v1(jj,j) = (x*c) + (z*s)
        !               v1(jj,i) = -(x*s) + (z*c)
        !   	        end do
        !            z = SQRT(f*f + h*h)
        !            w(j) = z
        !!            if (z.ne.0.0d0) then
        !            if (  dabs(z-0.0d0).gt.eps  ) then
        !               z = 1.0d0 / z
        !               c = f*z
        !               s = h*z
        !            end if
        !            f = (g*c) + (y*s)
        !            x = -(g*s) + (y*c)
        !            do jj = 1, m
        !               y = a(jj,j)
        !               z = a(jj,i)
        !               a(jj,j) = (y*c) + (z*s)
        !               a(jj,i) = -(y*s) + (z*c)
        !            end do
        !         end do         ! j loop
        !         rv1(l) = 0.0d0
        !         rv1(k) = f
        !         w(k) = x
        !      end do            ! its loop
        !3  continue
        !   end do               ! k loop
        !   
        !    open(24,file='Output_Data\VMatrixhigh2.txt')
        !    write(24,'(30f12.7)') V1          
        !    close(24)
        !   return
        !END SUBROUTINE svdcmp
        
end module Toolbox