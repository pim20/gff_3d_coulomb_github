module suzukit_1d
implicit none
contains
!this module is 1 dimension suzuki trotter
!_______________________________________________________________________________
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
!_______________________________________________________________________________
subroutine multiplication(a , b, multval)
  implicit none
  real(8), intent(in):: a,b
  real(8), intent(out):: multval
  multval = a*b
end subroutine



!_______________________________________________________________________________
!this is for setting the grid
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
!_______________________________________________________________________________
subroutine trapezoid(a,b,fx, summ)
  implicit none
  real(8), intent(in)::a,b
  real(8), intent(in):: fx(:)
  real(8), intent(out):: summ
  real(8):: h
  integer:: i,n
  n = size(fx)
  h= (b-a)/(n-1)
  summ= fx(1)+fx(n)
  do i =2, n-1
  	summ =summ+ 2.0d0*fx(i)
  end do

  summ= summ*h/2.0d0
end subroutine

!_______________________________________________________________________________
!this is 1D code for green's function of H0
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
!_______________________________________________________________________________
!This subroutine writes the different form of interaction potential(W) used
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
!_______________________________________________________________________________
!this subroutine is for potential part, writing w as exp(-w)
subroutine expow(beta, nwpot, x, expwval)
  implicit none
  real(8), intent(in):: beta, x
  real(8), intent(out):: expwval
  integer, intent(in):: nwpot
  real(8):: wx
  call wpot(beta, nwpot, x,wx)
  expwval = exp(-beta*(wx))

end subroutine




!_______________________________________________________________________________
!this subroutine is for suzuki trotter but m=1 term
subroutine suzuki_trotterm1(beta,nwpot, x, xprime, val)
  implicit none
  real(8), intent(in):: beta, x, xprime
  real(8), intent(out):: val
  real(8):: gf0, wx, expwval
  integer, intent(in)::nwpot
  call gfree(beta, x, xprime, gf0)
  call expow(beta, nwpot, x, expwval)
  call multiplication(gf0 , expwval, val)
  !val = gf0* exp(-(beta*wx))

end subroutine
!_______________________________________________________________________________
!this subroutine for for suzuki trotter m=2
! for the kernel type \bra{x} exp(-beta(H0+W))\ket{xprime}
!applying suzuki trotter, here m=2
!\bra{x}exp{-beta/2(H0)} exp(-beta/2(W)) exp{-beta/2(H0)} exp(-beta/2(W))\ket{xprime}
! introducing identity in between |int dx'' \ket{x''} \bra{x''}
!the terms become int dx'' \bra{x} exp{-beta/2(H0)} \ket{x''} exp(-beta/2(W(x'')))*
! \bra{x''}  exp{-beta/2(H0)} \ket{x'} exp(-beta/2(W(x')))
subroutine suzuki_trotterm2(beta, nwpot, ngridpt, dx, xgrid, valm2)
  implicit none
  integer(8), intent(in) :: ngridpt
  real(8), intent(in)    :: xgrid(ngridpt)
  real(8), intent(in)    :: beta, dx
  integer, intent(in)    :: nwpot
  real(8), intent(out)   :: valm2(ngridpt, ngridpt)

  integer :: i, j, k
  real(8) :: x, xprime, xdprime
  real(8) :: gf1, gf2, expwxp, expwxd, integrand, sum

  do i = 1, ngridpt
    x = xgrid(i)
    do j = 1, ngridpt
      xprime = xgrid(j)

      ! potential factor at x'
      call expow(beta/2.d0, nwpot, xprime, expwxp)

      sum = 0.d0
      do k = 1, ngridpt
        xdprime = xgrid(k)

        call gfree(beta/2.d0, x, xdprime, gf1)
        call gfree(beta/2.d0, xdprime, xprime, gf2)
        call expow(beta/2.d0, nwpot, xdprime, expwxd)

        integrand = gf1 * expwxd * gf2
        if (k == 1 .or. k == ngridpt) then
          sum = sum + 0.5d0 * integrand
        else
          sum = sum + integrand
        end if
      end do

      valm2(i,j) = expwxp * sum * dx
    end do
  end do
end subroutine suzuki_trotterm2

  !call gfree(beta, x, xdprime, gfxxd0)
  !call gfree(beta, xdprime, xprime, gfxdxp0)
  !call expow(beta, nwpot, xdprime, expwvalxdp)
  !call expow(beta, nwpot, xprime, expwvalxp)
end module
