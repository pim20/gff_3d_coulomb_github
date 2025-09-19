module suzukit_1d
implicit none
contains
!this module is 1 dimension suzuki trotter

subroutine setupgrid_call()
  implicit none
  real(8):: a,b,dx
  integer:: n
  b=2.0d0
  a=-2.0d0
  n= 12
  dx= (b-a)/(n-1)
  !print*, "dx", dx
end subroutine

subroutine setup_grid(x, n,a,b)
  implicit none
  real(8), intent(in) :: a,b
  integer, intent(in) :: n
  real(8), intent(out) :: x(n)
  real(8):: dx
  integer :: i
  dx= (b-a)/(n-1)
  do i = 1, n
    x(i) = a+ (i-1)* dx
  end do
end subroutine
subroutine gfree(beta, x, xprime,  gf0)
  implicit none
  real(8), intent(in):: x, xprime, beta
  real(8), intent(out):: gf0
  real(8):: a0,b0,pi
  !m=1, hbar=1
  pi= 4.0d0*atan(1.0d0)
  a0=sqrt(1.0d0/(2.0d0*beta*pi))
  b0= 1.0d0/(2.0d0*beta)
  gf0 = a0*exp(-b0*(x-xprime)**2)
end subroutine

subroutine wpot(beta, nwpot,x, wx)
  implicit none
  real(8), intent(in):: x, beta
  real(8), intent(out):: wx
  integer, intent(in):: nwpot
  select case(nwpot)
    case(1)
      wx=0.0d0
    case(2)
      wx= 2.0d0
  end select

end subroutine

subroutine suzuki_trotterm1(beta,nwpot, x, xprime, val)
  implicit none
  real(8), intent(in):: beta, x, xprime
  real(8), intent(out):: val
  real(8):: gf0, wx
  integer, intent(in)::nwpot
  call gfree(beta, x, xprime, gf0)
  call wpot(beta,nwpot,x, wx)
  val = gf0* exp(-wx)

end subroutine

end module
