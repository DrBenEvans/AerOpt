module Toolbox

    type OptimizationVariables

        integer, dimension(:), allocatable :: engInNodes           ! Engine inlet Nodes
        integer(kind = 4), dimension(8) :: timestart
        integer :: Gen
        double precision :: Precovery
        double precision, dimension(:,:), allocatable :: modes, coeff, modes2, coeff2       ! Modes and Coefficient derived by the POD method
        double precision, dimension(:), allocatable :: meanpressure         ! Mean Pressure excluded of Snapshots in POD
        double precision, dimension(:), allocatable :: pTamb                ! total pressure of free-stream flow of each Nest
        double precision, dimension(:,:), allocatable :: pressure           ! Static Pressure at every mesh node for each Nest
        double precision, dimension(:,:), allocatable :: MaLocal            ! Local Mach number at every mesh node for each Nest
        double precision, dimension(:), allocatable :: Fi                   ! Vector with all current Fitness values
        double precision, dimension(:), allocatable :: converged            ! Vector flagging converged Snapshots
        double precision, dimension(:), allocatable :: Precoutput           ! Vector with all current total pressure values
        double precision, dimension(:), allocatable :: PolCoeff, PolCoeff2  ! Polynomial Coefficients for POD RBF interpolation
        double precision, dimension(:,:), allocatable :: Weights, Weights2  ! Weights for POD RBF interpolation
        double precision, dimension(:,:), allocatable :: Nests_Move, Nests              ! Nest regenerated with each generation
        logical :: InitConv                                                 ! Initial Convergence Check - TRUE: yes

    end type OptimizationVariables

    type(OptimizationVariables) :: OV


    type CreateSnapshotsVariables

        double precision :: rn                                            ! Counts random numbers
        double precision, dimension(:,:), allocatable :: MxDisp           ! Matrix with min/max Displacements
        double precision, dimension(:,:), allocatable :: MxDisp_Move      ! Matrix with only moving min/max Displacements
        integer, dimension(:), allocatable :: cond                          ! Identifies zero/non-zero values
        double precision, dimension(:,:), allocatable :: Snapshots      ! Matrix containing the initial Nests

    end type CreateSnapshotsVariables

    type(CreateSnapshotsVariables) :: CS

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
    implicit none
    integer :: NoDim
    double precision :: Xa, Xb
    double precision, optional :: Ya, Yb, Za, Zb
    double precision ::  DistP2P

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
    double precision, dimension(2) :: AB, BC, AP, BP, p
    double precision, dimension(8) :: r
    logical :: Rectcheck

    ! Body or RectCheck
    Rectcheck = 0
    AB = (/(r(3)-r(1)), (r(4) - r(2))/) ! B-A -> (x,y)
    BC = (/(r(5)-r(3)), (r(6) - r(4))/) ! C-B -> (x,y)
    AP = (/(p(1)-r(1)), (p(2) - r(2))/) ! P-A -> (x,y)
    BP = (/(p(1)-r(3)), (p(2) - r(4))/) ! P-B -> (x,y)

    if (0 <= dot_product(AB, AP) .and. &
    dot_product(AB, AP) <= dot_product(AB, AB) .and. &
0   <= dot_product(BC, BP) .and. &
    dot_product(BC, BP) <= dot_product(BC, BC)) then
        Rectcheck = 1
    end if

    end function RectCheck

    Subroutine PRINT_MATRIX( DESC, M, N, A, LDA )
    ! Objective: Print a Matrix

    ! Variables
    implicit none
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
    implicit none
    integer, intent(in) :: nA
    double precision, dimension(nA), intent(in out) :: A
    integer, dimension(nA), intent(in out), optional :: ind
    character(len=1), optional :: sel

    ! Local Variables
    integer :: left, right
    real :: random
    double precision :: pivot
    integer :: marker
    double precision :: temp

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

    recursive subroutine QSortInt(a,na, sel, ind)
    ! Objective: Sort/Order an Array in ascending order

    ! Dummy Arguments
    implicit none
    integer, intent(in) :: nA
    integer, dimension(nA), intent(in out) :: A
    integer, dimension(nA), intent(in out), optional :: ind
    character(len=1), optional :: sel

    ! Local Variables
    integer :: left, right
    real :: random
    double precision :: pivot
    integer :: marker
    double precision :: temp

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

    call QSortInt(A(:marker-1),marker-1, sel, ind(:marker-1)) ! recursive call -1
    call QSortInt(A(marker:),nA-marker+1, sel, ind(marker:))  ! recursive call +1

    end if

    end subroutine QSortInt

    subroutine Unique(vector_in, sz, vector_out)

    ! Variables
    implicit none
    integer :: sz, i, j
    double precision, dimension(sz), intent(in) :: vector_in
    double precision, dimension(:), allocatable, intent(out) :: vector_out
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

    subroutine UniqueInt(vector_in, sz, vector_out)

    ! Variables
    implicit none
    integer :: sz, i, j
    integer, dimension(sz), intent(in) :: vector_in
    integer, dimension(:), allocatable, intent(out) :: vector_out
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

    end subroutine UniqueInt

    subroutine randperm(N, p)

    !! Source: coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2006-03/msg00748.html
    !! Based on Knuth's algorithm
    implicit none
    integer, intent(in) :: N
    integer, dimension(:), intent(out) :: p
    double precision :: rn
    integer :: temp, j, k, i

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

    subroutine timestamp ( )

    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    May 31 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    31 May 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
      implicit none

      character ( len = 8 ) ampm
      integer ( kind = 4 ) d
      character ( len = 8 ) date
      integer ( kind = 4 ) h
      integer ( kind = 4 ) m
      integer ( kind = 4 ) mm
      character ( len = 9 ), parameter, dimension(12) :: month = (/ &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' /)
      integer ( kind = 4 ) n
      integer ( kind = 4 ) s
      character ( len = 10 ) time
      integer ( kind = 4 ) values(8)
      integer ( kind = 4 ) y
      character ( len = 5 ) zone

      call date_and_time ( date, time, zone, values )

      y = values(1)
      m = values(2)
      d = values(3)
      h = values(5)
      n = values(6)
      s = values(7)
      mm = values(8)

      if ( h < 12 ) then
        ampm = 'AM'
      else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h < 12 ) then
          ampm = 'PM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
        trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

      return

    end subroutine timestamp

    subroutine DetermineStrLen(istr, i)

    ! Variables
    implicit none
    integer :: i
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
    elseif (i < 10000) then
        allocate(character(len=4) :: istr)
        write( istr, '(I4)' )  i
    else
        allocate(character(len=5) :: istr)
        write( istr, '(I5)' )  i
    end if

    end subroutine DetermineStrLen

    subroutine communicateWin2Lin(Username, Password, fname, channel)

    ! Variables
    implicit none
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

    function Det(A, N)
    ! Objective: Calculates Determinant of Matrix

    ! Variables
    implicit none
    integer :: N, ID, j, icode
    double precision, dimension(N,N) :: A
    double precision :: Det
    integer, dimension(N) :: INDX(N)


    ! Body of Det
    call LUDCMP(A,N,INDX,ID,ICODE)  !see module lu.f90

    Det=DFLOAT(ID)

    do J=1, N
        Det=Det*A(J,J)
    end do

    if (abs(Det) < 1.5D-16) then
        Det = 0
    end if

    end function Det

    Subroutine LUDCMP(A,N,INDX,D,CODE)
    ! Objective: Performs the LU - Decomposition

    ! Variables
    implicit none
    integer, parameter :: NMAX = 100
    double precision, parameter :: TINY = 1.5D-17
    double precision ::  AMAX,DUM, SUM, A(N,N),VV(NMAX)
    integer :: CODE, D, INDX(N), I, J, K, IMAX, N

    ! Body of LU Decomposition
    D=1; CODE=0

    DO I=1,N
        AMAX=0.d0
        DO J=1,N
            IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
        END DO ! j loop
        IF(AMAX.LT.TINY) THEN
            CODE = 1
            RETURN
        END IF
        VV(I) = 1.d0 / AMAX
    END DO ! i loop

    DO J=1,N
        DO I=1,J-1
            SUM = A(I,J)
            DO K=1,I-1
                SUM = SUM - A(I,K)*A(K,J)
            END DO ! k loop
            A(I,J) = SUM
        END DO ! i loop
        AMAX = 0.d0
        DO I=J,N
            SUM = A(I,J)
            DO K=1,J-1
                SUM = SUM - A(I,K)*A(K,J)
            END DO ! k loop
            A(I,J) = SUM
            DUM = VV(I)*DABS(SUM)
            IF(DUM.GE.AMAX) THEN
                IMAX = I
                AMAX = DUM
            END IF
        END DO ! i loop

        IF(J.NE.IMAX) THEN
            DO K=1,N
                DUM = A(IMAX,K)
                A(IMAX,K) = A(J,K)
                A(J,K) = DUM
            END DO ! k loop
            D = -D
            VV(IMAX) = VV(J)
        END IF

        INDX(J) = IMAX
        IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

        IF(J.NE.N) THEN
            DUM = 1.d0 / A(J,J)
            DO I=J+1,N
                A(I,J) = A(I,J)*DUM
            END DO ! i loop
        END IF
    END DO ! j loop

    RETURN
    END subroutine LUDCMP

    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function inv(A) result(Ainv)

      implicit none
      double precision, dimension(:,:), intent(in) :: A
      double precision, dimension(size(A,1),size(A,2)) :: Ainv

      double precision, dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info

      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
    end function inv

    integer function fact(n)

        implicit none
        integer, intent(in) :: n
        integer :: i, a

        a = 1
        do i = 1, n
            a = a*i
        end do
        fact = a

    end function fact

    end module Toolbox
