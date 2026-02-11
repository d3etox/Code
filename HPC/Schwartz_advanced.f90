program schwartz_methods
    use mpi
    implicit none

    integer :: ierr, rank, nprocs
    integer :: choice

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

    if (rank == 0) then
        print *, "========================"
        print *, " Methode de Schwarz 1D "
        print *, "========================"
        print *, "1. Resolution globale"
        print *, "2. Schwarz MPI"
        print *, "3. Schwartz N processus avec benchmark"
        print *, ""
        choice = 2
    end if

    call MPI_Bcast(choice, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    select case(choice)
    case(1)
        if (rank == 0) call solve_global()

    case(2)
        call solve_schwarz_mpi(rank, nprocs, MPI_COMM_WORLD)

    case(3)
        call solve_schwartz_nprocs(rank, nprocs, MPI_COMM_WORLD)
        
    case default
        if (rank == 0) print *, "Choix invalide"
    end select

    call MPI_Finalize(ierr)

contains
!=========================================================
subroutine solve_global()
    implicit none
    integer :: N, i
    real(8) :: a, b, h, max_error
    real(8), allocatable :: x(:), u(:), f(:), c(:)
    real(8), allocatable :: lower(:), diag(:), upper(:), rhs(:)

    a = 0.d0; b = 1.d0
    N = 200
    h = (b-a)/(N+1)

    allocate(x(0:N+1), u(0:N+1), f(1:N), c(1:N))
    allocate(lower(1:N-1), diag(1:N), upper(1:N-1), rhs(1:N))

    u(0) = 0.d0; u(N+1) = 0.d0

    do i = 0, N+1
        x(i) = a + i*h
    end do

    f = 1.d0
    c = 0.d0

    diag = 2.d0 + h*h*c
    lower = -1.d0
    upper = -1.d0
    rhs = h*h*f

    call solve_tridiag_thomas(lower, diag, upper, rhs, u(1:N), N)

    open(10,file="solution_global.dat",status="replace")
    do i = 0, N+1
        write(10,'(2F12.6)') x(i), u(i)
    end do
    close(10)

    max_error = maxval(abs(u(1:N) - 0.5d0*x(1:N)*(1.d0-x(1:N))))
    print *, "Erreur max globale :", max_error

    deallocate(x,u,f,c,lower,diag,upper,rhs)
end subroutine solve_global

!=========================================================
subroutine solve_schwarz_mpi(rank, nprocs, comm)
    use mpi
    implicit none   

    integer, intent(in) :: rank, nprocs, comm

    integer, parameter :: N = 400
    integer, parameter :: max_iter = 5000
    real(8), parameter :: tol = 1.d-8
    real(8), parameter :: L = 1.d0
    real(8), parameter :: overlap_ratio = 0.2d0

    integer :: i, iter, ierr
    integer :: Nloc, Nwork, overlap
    real(8) :: h, err, err_loc

    real(8), allocatable :: x(:), u(:), u_old(:), f(:)
    real(8), allocatable :: lower(:), diag(:), upper(:), rhs(:)

    integer :: left, right
    real(8) :: sendL, sendR, recvL, recvR
    character(len=10) :: rank_str

    ! -----------------------------
    ! Découpage du domaine
    ! -----------------------------
    h = L / (N+1)

    Nwork  = N / nprocs
    overlap = max(2, int(overlap_ratio * Nwork))
    Nloc = Nwork + 2*overlap

    allocate(x(0:N+1), u(0:Nloc+1), u_old(0:Nloc+1), f(1:Nloc))
    allocate(lower(1:Nloc-1), diag(1:Nloc), upper(1:Nloc-1), rhs(1:Nloc))

    u = 0.d0
    f = 1.d0

    lower = -1.d0
    upper = -1.d0
    diag  =  2.d0

    do i = 0, N+1
        x(i) = i*h
    end do

    ! Voisins MPI
    left  = rank - 1
    right = rank + 1
    if (left  < 0) left  = MPI_PROC_NULL
    if (right >= nprocs) right = MPI_PROC_NULL

    ! Conditions aux limites globales
    if (rank == 0)        u(0) = 0.d0
    if (rank == nprocs-1) u(Nloc+1) = 0.d0

    ! -----------------------------
    ! Itérations de Schwarz
    ! -----------------------------
    do iter = 1, max_iter

        u_old = u

        ! ---- échange des recouvrements ----
        sendL = u(overlap+1)
        sendR = u(Nloc-overlap)

        call MPI_SENDRECV(sendL,1,MPI_DOUBLE_PRECISION,left,0, &
                          recvR,1,MPI_DOUBLE_PRECISION,right,0,comm,MPI_STATUS_IGNORE,ierr)

        call MPI_SENDRECV(sendR,1,MPI_DOUBLE_PRECISION,right,1, &
                          recvL,1,MPI_DOUBLE_PRECISION,left,1,comm,MPI_STATUS_IGNORE,ierr)

        if (right /= MPI_PROC_NULL) u(Nloc+1) = recvR
        if (left  /= MPI_PROC_NULL) u(0)      = recvL

        ! ---- second membre ----
        rhs = h*h * f

        rhs(overlap+1)        = rhs(overlap+1)        + u(0)
        rhs(Nloc-overlap)     = rhs(Nloc-overlap)     + u(Nloc+1)

        ! ---- résolution sur partie interne ----
        call solve_tridiag_thomas( lower(overlap+1:Nloc-overlap-1), &
                     diag (overlap+1:Nloc-overlap),   &
                     upper(overlap+1:Nloc-overlap-1), &
                     rhs  (overlap+1:Nloc-overlap),   &
                     u    (overlap+1:Nloc-overlap),   &
                     Nloc-2*overlap )

        ! ---- critère de convergence ----
        err_loc = maxval(abs(u(overlap+1:Nloc-overlap) - u_old(overlap+1:Nloc-overlap)))
        call MPI_ALLREDUCE(err_loc, err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

    
        if (err < tol) exit

    end do

        ! Écrire seulement la partie locale "interne"
        write(rank_str,'(I0)') rank
        open(10,file='solution_rank'//trim(rank_str)//'.dat',status='replace')

        do i = overlap+1, Nloc-overlap
            write(10,'(2F12.6)') x(rank*Nwork + i - overlap), u(i)
        end do

        close(10)


    if (rank==0) print *, "Convergence Schwarz en", iter, "iterations"

    deallocate(u, u_old, f, lower, diag, upper, rhs)

end subroutine solve_schwarz_mpi


subroutine solve_schwartz_nprocs(rank, nprocs, comm)
    use mpi
    implicit none   

    integer, intent(in) :: rank, nprocs, comm

    integer, parameter :: N = 400
    integer, parameter :: max_iter = 5000
    real(8), parameter :: tol = 1.d-8
    real(8), parameter :: L = 1.d0
    real(8), parameter :: overlap_ratio = 0.2d0

    integer :: i, iter, ierr
    integer :: Nloc, Nwork, overlap, istart, iend, Nlast
    real(8) :: h, err, err_loc

    real(8), allocatable :: x(:), u(:), u_old(:), f(:)
    real(8), allocatable :: lower(:), diag(:), upper(:), rhs(:)

    integer :: left, right
    real(8) :: sendL, sendR, recvL, recvR
    character(len=10) :: rank_str

    ! -----------------------------
    ! Découpage du domaine
    ! -----------------------------
    h = L / (N+1)

    ! Points internes par rang
    Nwork = N / nprocs
    ! Dernier rang prend le reste si N n'est pas divisible
    if (rank == nprocs-1) Nlast = N - Nwork*(nprocs-1)
    overlap = max(2, int(overlap_ratio * Nwork))
    Nloc = Nwork + 2*overlap
    if (rank == nprocs-1) Nloc = Nlast + 2*overlap

    allocate(x(0:N+1), u(0:Nloc+1), u_old(0:Nloc+1), f(1:Nloc))
    allocate(lower(1:Nloc-1), diag(1:Nloc), upper(1:Nloc-1), rhs(1:Nloc))

    u = 0.d0
    f = 1.d0
    lower = -1.d0
    upper = -1.d0
    diag  =  2.d0

    do i = 0, N+1
        x(i) = i*h
    end do

    ! Voisins MPI
    left  = rank - 1
    right = rank + 1
    if (left  < 0) left  = MPI_PROC_NULL
    if (right >= nprocs) right = MPI_PROC_NULL

    ! Conditions aux limites globales
    if (rank == 0)        u(0) = 0.d0
    if (rank == nprocs-1) u(Nloc+1) = 0.d0

    ! -----------------------------
    ! Itérations de Schwarz
    ! -----------------------------
    do iter = 1, max_iter

        u_old = u

        ! ---- échange des recouvrements ----
        sendL = u(overlap+1)
        sendR = u(Nloc-overlap)

        call MPI_SENDRECV(sendL,1,MPI_DOUBLE_PRECISION,left,0, &
                          recvR,1,MPI_DOUBLE_PRECISION,right,0,comm,MPI_STATUS_IGNORE,ierr)

        call MPI_SENDRECV(sendR,1,MPI_DOUBLE_PRECISION,right,1, &
                          recvL,1,MPI_DOUBLE_PRECISION,left,1,comm,MPI_STATUS_IGNORE,ierr)

        if (right /= MPI_PROC_NULL) u(Nloc+1) = recvR
        if (left  /= MPI_PROC_NULL) u(0)      = recvL

        ! ---- second membre ----
        rhs = h*h * f
        rhs(overlap+1)        = rhs(overlap+1)        + u(0)
        rhs(Nloc-overlap)     = rhs(Nloc-overlap)     + u(Nloc+1)

        ! ---- résolution sur partie interne ----
        call solve_tridiag_thomas( lower(overlap+1:Nloc-overlap-1), &
                     diag (overlap+1:Nloc-overlap),   &
                     upper(overlap+1:Nloc-overlap-1), &
                     rhs  (overlap+1:Nloc-overlap),   &
                     u    (overlap+1:Nloc-overlap),   &
                     Nloc-2*overlap )

        ! ---- critère de convergence ----
        err_loc = maxval(abs(u(overlap+1:Nloc-overlap) - u_old(overlap+1:Nloc-overlap)))
        call MPI_ALLREDUCE(err_loc, err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

        if (rank==0 .and. mod(iter,100)==0) print '(A,I5,A,E12.5)', "Iter ", iter, " err = ", err
        if (err < tol) exit

    end do

    ! -----------------------------
    ! Écriture des fichiers locaux
    ! -----------------------------
    write(rank_str,'(I0)') rank
    open(10,file='solution_rank'//trim(rank_str)//'.dat',status='replace')

    ! indices globaux
    istart = rank*Nwork + 1
    iend   = istart + (Nloc - 2*overlap) - 1
    if (rank == nprocs-1) iend = N  ! dernier rang ajuste à N

    do i = overlap+1, Nloc-overlap
        write(10,'(2F12.6)') x(istart + i - overlap - 1), u(i)
    end do

    close(10)

    if (rank==0) print *, "Convergence Schwarz en", iter, "iterations"

    deallocate(u, u_old, f, lower, diag, upper, rhs)

end subroutine solve_schwartz_nprocs

!=========================================================
subroutine solve_tridiag_thomas(a,b,c,d,u,n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: a(:), b(:), c(:)
    real(8), intent(inout) :: d(:)
    real(8), intent(out) :: u(:)
    real(8), allocatable :: c_star(:)
    integer :: i

    allocate(c_star(1:n-1))

    c_star(1) = c(1)/b(1)
    d(1) = d(1)/b(1)

    do i = 2, n-1
        c_star(i) = c(i)/(b(i) - a(i-1)*c_star(i-1))
    end do

    do i = 2, n
        d(i) = (d(i) - a(i-1)*d(i-1)) / (b(i) - a(i-1)*c_star(i-1))
    end do

    u(n) = d(n)
    do i = n-1, 1, -1
        u(i) = d(i) - c_star(i)*u(i+1)
    end do

    deallocate(c_star)
end subroutine solve_tridiag_thomas

end program schwartz_methods