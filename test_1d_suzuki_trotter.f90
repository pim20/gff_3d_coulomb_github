program test_suzukitrit_1d
use suzukit_1d
implicit none

call test_suzuki_trotter_m2()
!call test_suzuki_trotter_m3()

contains
!_______________________________________________________________________________
!This subroutine test e^{-beta(H0+W)} suzuki trotter m=1 (first order)
!W is taken zero and constant represented by (nwpot)
subroutine test_suzuki_trotter_m1()
  implicit none
  real(8)::beta, val, x, xprime, a,b,dx
  real(8), allocatable::xgrid(:)
  integer:: nwpot, n, i,j
  nwpot=2
  beta = 1.0d0
  n=100
  a=-2.0d0
  b=2.0d0
  open(unit=1, file= "testsuzukit1d_wconst.txt", status= "replace", action= "write")
  allocate(xgrid(n))
  call setup_grid(xgrid, n,a,b,dx)
  do i =1, n
    x =xgrid(i)
    do j =1,n
      xprime = xgrid(j)
      call suzuki_trotterm1(beta,nwpot, x, xprime, val)
      write(1,*) x, xprime, val
    end do
  end do
  !call setupgrid_call()
  deallocate(xgrid)
  close(1)
  end subroutine test_suzuki_trotter_m1
!_______________________________________________________________________________
! This subroutine test e^{-beta(H0+W)} suzuki trotter m=2 (Second order)

subroutine test_suzuki_trotter_m2()
  implicit none
  real(8) :: beta, x, xprime
  integer :: nwpot, n,i,j
  real(8) :: dx,a,b
  real(8), allocatable :: xgrid(:)
  real(8), allocatable :: valm2(:,:)
!
  beta  = 1.0d0
  nwpot =  1
  n     =  100
  a     =  -2.0d0
  b     =  2.0d0
  dx    =  (b-a)/(n-1)

  allocate(xgrid(n))
  allocate(valm2(n,n))
  !do i = 1,n
  !  xgrid(i) = a+(i-1)*dx
  !end do
  call suzuki_trotterm2(beta, nwpot, n, a,b,xgrid, valm2)
  do i =1, n
    x = xgrid(i)
    do j =1, n
      xprime= xgrid(j)
      write(*,*) x, xprime,valm2(i,j)
    end do
  end do

  deallocate(xgrid)
  deallocate(valm2)

end subroutine
!_______________________________________________________________________________
subroutine test_suzuki_trotter_m3()
  implicit none
  real(8) :: beta, x, xprime
  integer :: nwpot, n,i,j
  real(8) :: dx,a,b
  real(8), allocatable :: xgrid(:)
  real(8), allocatable :: valm3(:,:)
!
  beta  = 1.0d0
  nwpot =  1
  n     =  100
  a     =  -2.0d0
  b     =  2.0d0
  dx    =  (b-a)/(n-1)

  allocate(xgrid(n))
  allocate(valm3(n,n))
  !do i = 1,n
  !  xgrid(i) = a+(i-1)*dx
  !end do
  call suzuki_trotterm3(beta, nwpot, n, a,b, xgrid,valm3)
  do i =1, n
    x = xgrid(i)
    do j =1, n
      xprime= xgrid(j)
      write(*,*) x, xprime,valm3(i,j)
    end do
  end do

  deallocate(xgrid)
  deallocate(valm3)

end subroutine




end program
