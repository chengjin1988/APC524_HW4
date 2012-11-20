!*****************************************
!           APC524  Homework 4
!          MPI implementation
!*****************************************

program heat_mpi
    implicit none
    include 'mpif.h'
    real(kind=8), parameter :: pi = 2*acos(0.0)
    real(kind=8), parameter :: k = 1.0, time=10.0*pi**2.0/k
    real(kind=8), parameter :: factr=0.9
    character(len=100) :: arg
    integer :: nx, ny, it, nt, i, j
    real(kind=8) :: dt, dx, x, current_x, dfactr, averageT
    real(kind=8), allocatable :: T(:,:), deltaT(:,:), buff(:,:)

    real(kind=8) :: tstart, tend, telapsed

! MPI parameters
    integer :: rank, nproc, ierror, tag1, tag2, status
    integer :: rank_left, rank_right, irank

!    open(unit=1, file='final.dat')
    call mpi_init(ierror)
    call mpi_comm_size(mpi_comm_world, nproc, ierror)
    call mpi_comm_rank(mpi_comm_world, rank, ierror)
! Read input from the command line
    if (iargc() /= 1) then
       write(*,'(a)') "usage:  ./heat_mpi <nx>"
       stop
    else
       call getarg(1,arg)
       read(arg,*) ny
    endif

! Set up dt and dx
    dx = pi/ny
    dt = factr*(dx**2.0)/(4.0*k)
    nt =  floor(time/dt)
    dfactr = dt*k/(dx*dx)

! Allocate memory
    nx = ny/nproc                      ! domain decomposition
    allocate(T(nx+2, ny+1), deltaT(nx+2, ny+1), buff(nx+2, ny+1))
! Setup the initial condition for each domain

    if ( rank == 0 ) then
       call cpu_time(tstart)
    end if

    current_x = rank*(pi/nproc)
    do i = 1, nx+2
       x = current_x + dx*(i-1)
       T(i, 1) = (cos(x))**2                 ! boundary condition at y = 0
       T(i, ny+1) = (sin(x))**2              ! boundary condition at y = pi
    enddo
    T(:, 2:ny) = 0

! Start integrating
    rank_left = rank - 1
    rank_right = rank + 1
    if (rank == 0)    rank_left = nproc - 1
    if (rank == nproc - 1)    rank_right = 0

    do it = 1, nt+1

       ! DO LAPLACIAN
       ! boundary condition in y
       deltaT(:,1) = 0
       deltaT(:,ny+1) = 0

       ! interior points in each domain
       do j = 2, ny
          do i = 2, nx+1
             deltaT(i, j) = dfactr*( T(i-1,j) + T(i+1, j) + T(i, j-1) + T(i, j+1) - 4*T(i,j) )
          enddo
       enddo
       ! END LAPLACIAN

       ! DO EULER IN TIME
       ! for interior points in each domain
       do j = 1, ny+1
           do i = 2, nx+1
              T(i,j) = T(i,j) + deltaT(i,j)
           enddo
       enddo

       ! for the boundary points in each domain by MPI_ISEND and MPI_IRECV

       ! send out data
       call MPI_SEND(T(2,1:ny+1), ny+1, MPI_DOUBLE_PRECISION, rank_left, tag1, MPI_COMM_WORLD, ierror)
       call MPI_SEND(T(nx+1,1:ny+1), ny+1, MPI_DOUBLE_PRECISION, rank_right, tag2, MPI_COMM_WORLD, ierror)
       ! receive data
       call MPI_RECV(T(1, 1:ny+1), ny+1, MPI_DOUBLE_PRECISION, rank_left, tag2, MPI_COMM_WORLD, status, ierror)
       call MPI_RECV(T(nx+2, 1:ny+1), ny+1, MPI_DOUBLE_PRECISION, rank_right, tag1, MPI_COMM_WORLD, status, ierror)
    enddo

    dx = pi/ny
    if ( rank/= 0 ) then         ! all the rest processors send to processor 0
       call MPI_SEND(T(:,:), (nx+2)*(ny+1), MPI_DOUBLE_PRECISION, 0, tag1, MPI_COMM_WORLD, ierror)
    else                         ! irank = 0 collects the data to evaluate average Temp
       averageT = averageT+sum(T(1:nx,:))
       do irank = 1, nproc-1
          call MPI_RECV(buff(:,:), (nx+2)*(ny+1), MPI_DOUBLE_PRECISION, irank, tag1, MPI_COMM_WORLD, status, ierror)
          averageT = averageT+sum(buff(1:nx, :))
       enddo
       averageT = averageT/(4*nx*(ny+1))
    endif

    if ( rank == 0 ) then
       call cpu_time(tend)
       telapsed = tend-tstart
       write(*,'(a,i5, a, i5)') "nx = ", nx, "   ny = ", ny
       write(*,'(a, e12.6, a, e12.6)') "dx = ", dx, "    factor = ", factr
       write(*,'(a, e12.6)') "time = ", time
       write(*,'(a, e12.6)') "volume averaged temperature =  ", averageT
       write(*,'(a, e14.6)') "Elapsed CPU time = ", telapsed
    endif
    do i = 2, nx+1
       do j = 1, ny+1
          x = dx*(i-1) + current_x
          write(*, '(e14.6, e14.6, e14.6)') x, dx*(j-1), T(i, j)
       enddo
    enddo


    call MPI_FINALIZE(ierror)
end program

