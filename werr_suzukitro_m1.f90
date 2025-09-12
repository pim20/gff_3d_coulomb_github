program suzuki_trotter_m1_grid
  implicit none
  double precision :: beta, r(3), rprime(3), val
  integer :: i, j, k, nx, ny, nz
  double precision :: xmax, ymax, zmax, dx, dy, dz
  double precision :: x, y, z

  ! parameters
  beta = 1.0d0
  r = (/0.0d0, 0.0d0, 0.0d0/)    ! fix r at origin

  nx = 10   ! number of x points
  ny = 10   ! number of y points
  nz = 10   ! number of z points

  xmax = 3.0d0
  ymax = 3.0d0
  zmax = 3.0d0

  dx = xmax/(nx-1)
  dy = ymax/(ny-1)
  dz = zmax/(nz-1)

  print *, "# x   y   z   K^(1)(0,(x,y,z);beta)"

  do i=1,nx
     x = (i-1)*dx
     do j=1,ny
        y = (j-1)*dy
        do k=1,nz
           z = (k-1)*dz
           rprime = (/x, y, z/)
           call kernel_m1(beta, r, rprime, val)
           print *, x, y, z, val
        end do
     end do
  end do

contains

  ! Potential W(r) = erf(r)/r
  double precision function potential(r)
    implicit none
    double precision, intent(in) :: r(3)
    double precision :: rr
    rr = sqrt(r(1)**2 + r(2)**2 + r(3)**2)
    if (rr .eq. 0d0) then
       potential = 2d0/sqrt(atan(1.0d0))   ! limit r→0
    else
       potential = erf(rr)/rr
    endif
  end function potential

  ! Free particle propagator
  double precision function k0(r, rprime,beta)
    implicit none
    double precision, intent(in) :: r(3), rprime(3), beta
    double precision :: dr2, pi,a0,b0
    pi= 4.0d0*atan(1.0d0)
    dr2 = (r(1)-rprime(1))**2 + (r(2)-rprime(2))**2 + (r(3)-rprime(3))**2
    a0=sqrt((1.0d0/(2.0d0*beta*pi))**3)
    b0= 1.0d0/(2.0d0*beta)
    k0 = a0* exp(-dr2*b0)
  end function k0

  ! Kernel for m=1
  subroutine kernel_m1(beta, r, rprime, val)
    implicit none
    double precision, intent(in) :: beta, r(3), rprime(3)
    double precision, intent(out):: val
    val = k0(r, rprime, beta) * exp(-beta * potential(rprime))
  end subroutine kernel_m1

end program suzuki_trotter_m1_grid
