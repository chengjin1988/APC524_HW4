!*****************************************
!           APC524  Homework 4
!          serial implementation
!*****************************************

program heat_serial
    implicit none
    real(kind=8), parameter :: pi = 2*acos(0.0)
    real(kind=8), parameter :: k = 1.0, time=10.0*pi**2.0/k
    real(kind=8), parameter :: factr=0.9
    character(len=100) :: arg
    integer :: n, nx, it, nt, i, j
    real(kind=8) :: dt, dx, x, dfactr
    real(kind=8), allocatable :: T(:,:), deltaT(:,:)

    real(kind=8) :: tstart, tend, telapsed
! Read input from the command line
    n = iargc()
    if (n /= 1) then
       write(*,'(a)') "usage:  ./heat_serial <nx>"
       stop
    else
       call getarg(1,arg)
       read(arg,*) nx
    endif

! Set up dt and dx
    dx = pi/nx
    dt = factr*(dx**2.0)/(4.0*k)
    nt =  floor(time/dt)
    dfactr = dt*k/(dx*dx)

! Allocate memory
    allocate(T(nx, nx+1), deltaT(nx, nx+1))

! Setup the initial condition
    call cpu_time(tstart)
    do i = 1, nx
       x = dx*(i-1)
       T(i, 1) = (cos(x))**2                 ! boundary condition at y = 0
       T(i, nx+1) = (sin(x))**2              ! boundary condition at y = pi
    enddo
    do j = 2, nx
       do i = 1, nx
          T(i,j) = 0
       enddo
    enddo
! Start integrating
    do it = 1, nt+1

       ! DO LAPLACIAN

       ! boundary condition in x
       do j = 2, nx
          deltaT(1,j) = dfactr*(T(nx,j) + T(2,j) + T(1, j-1) + T(1, j+1) - 4*T(1, j))
          deltaT(nx, j) = dfactr*(T(nx-1,j) + T(1,j) + T(nx, j-1) + T(nx, j+1) - 4*T(nx, j))
       enddo

       ! interior points
       do j = 2, nx
          do i = 2, nx-1
             deltaT(i, j) = dfactr*( T(i-1,j) + T(i+1, j) + T(i, j-1) + T(i, j+1) - 4*T(i,j) )
          enddo
       enddo
       ! END LAPLACIAN

       ! DO EULER
       do j = 2, nx
           do i = 1, nx
              T(i,j) = T(i,j) + deltaT(i,j)
           enddo
       enddo
    enddo

    call cpu_time(tend)
    telapsed = tend - tstart

    write(*,'(a,i5, a, e12.6)') "nx = ", nx, "    time = ", time
    write(*,'(a, e12.6, a, e12.6)') "dx = ", dx, "    factor = ", factr
    write(*,'(a, e12.6)') "volume averaged temperature =  ", sum(T)/(nx*(nx+1))
    write(*,'(a, e14.6)') "Elapsed CPU time = ", telapsed

    do i = 1, nx
       do j = 1, nx+1
          write(*, '(e14.6, e14.6, e14.6)') dx*(i-1), dx*(j-1), T(i, j)
       enddo
    enddo
    close(1)
end program

