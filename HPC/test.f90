program solve_poisson
    implicit none
    integer :: N, i
    real(8) :: a, b, h, max_error, residual
    real(8), allocatable :: x(:), u(:), f(:), c(:)
    real(8), allocatable :: lower(:), diag(:), upper(:), rhs(:) ! Tridiagonal

    !--- Paramètres du problème ---
    a = 0.0d0
    b = 1.0d0
    N = 100       ! nombre de points internes
    h = (b - a)/(N+1)

    allocate(x(0:N+1), u(0:N+1), f(1:N), c(1:N))
    allocate(lower(1:N-1))
    allocate(diag(1:N))
    allocate(upper(1:N-1))
    allocate(rhs(1:N))

    !--- Conditions aux limites ---
    u(0) = 0.0d0  ! ua
    u(N+1) = 0.0d0 ! ub

    !--- Maillage et fonctions ---
    do i = 0, N+1
        x(i) = a + i*h
    end do

    do i = 1, N
        f(i) = 1.0d0                  ! exemple: f(x) = 1
        c(i) = 0.0d0                  ! exemple: c(x) = 0
    end do

    !--- Construction de la matrice tridiagonale ---
    do i = 1, N
        diag(i) = 2.0d0 + h**2 * c(i)   ! diagonale
    end do
    do i = 1, N-1
        lower(i) = -1.0d0                 ! sous-diagonale
        upper(i) = -1.0d0                 ! sur-diagonale
    end do

    !--- Construction du côté droit ---
    do i = 1, N
        rhs(i) = h**2 * f(i)
    end do
    rhs(1) = rhs(1) + u(0)                ! prendre en compte la condition u(a)
    rhs(N) = rhs(N) + u(N+1)              ! prendre en compte la condition u(b)

    !--- Résolution tridiagonale par Thomas ---
    call thomas(lower, diag, upper, rhs, u(1:N), N)

    !--- Écriture dans un fichier pour plot ---
    open(unit=10, file='solution.dat', status='replace')
    do i = 0, N+1
        write(10,'(2F12.6)') x(i), u(i)
    end do
    close(10)

    max_error = maxval(abs(u(1:N) - 0.5d0*x(1:N)*(1.0d0-x(1:N))))
    print *, "Erreur maximale par rapport à la solution exacte :", max_error
    residual = 0.0d0
    do i = 1, N
        residual = max(residual, abs(-u(i-1) + (2 + h**2*c(i))*u(i) - u(i+1) - h**2*f(i)))
    end do
    print *, "Max residual =", residual

    !--- Affichage ---
    print *, "x_i      u_i"
    do i = 0, N+1
        print "(F8.4,2X,F8.4)", x(i), u(i)
    end do

    contains

    subroutine thomas(a,b,c,d,u,n)
        integer, intent(in) :: n
        real(8), intent(in) :: a(:), b(:), c(:)
        real(8), intent(inout) :: d(:)
        real(8), intent(out) :: u(:)
        real(8), allocatable :: c_star(:)
        integer :: i

        allocate(c_star(1:n-1))
        !--- Forward sweep ---
        c_star(1) = c(1)/b(1)
        d(1) = d(1)/b(1)
        do i = 2, n-1
            c_star(i) = c(i)/(b(i) - a(i-1)*c_star(i-1))
        end do
        do i = 2, n
            d(i) = (d(i) - a(i-1)*d(i-1)) / (b(i) - a(i-1)*c_star(i-1))
        end do

        !--- Back substitution ---
        u(n) = d(n)
        do i = n-1, 1, -1
            u(i) = d(i) - c_star(i)*u(i+1)
        end do
    end subroutine thomas

end program solve_poisson
