program test_suzukitrit_1d
use suzukit_1d
implicit none
!-------------------------------------------------------------------------------
!This test file is for different step of suzuki trotter that is m
!nwpot is the case for interaction potential(W)
!nwpot =1 is W=0
!nwpot =1 is W= const
!-------------------------------------------------------------------------------

!call test_suzuki_trotter_m1()
!call test_suzuki_trotter_m2()
!call test_suzuki_trotter_m3()
call test_gfree()

contains
!_______________________________________________________________________________
!This subroutine test e^{-beta(H0+W)} suzuki trotter m=1 (first order)
!W is taken zero and constant represented by (nwpot)
subroutine test_suzuki_trotter_m1()
  implicit none
  real(8)::beta, val, x, xprime, a,b,dx
  real(8), allocatable::xgrid(:)
  integer:: nwpot, n, i,j
  nwpot=1
  beta = 1.0d0
  n=100
  a=-2.0d0
  b=2.0d0
  open(unit=1, file= "test_st_m1_w0.txt", status= "replace", action= "write")
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
  integer :: ierr
!
  beta  = 1.0d0
  nwpot =  1
  n     =  100
  a     =  -2.0d0
  b     =  2.0d0
  dx    =  (b-a)/(n-1)

  allocate(xgrid(n))
  allocate(valm2(n,n))
  open(unit=4, file="test_st_m2_w0.txt", status="replace", action="write")
!  if (ierr /= 0) then
!    print *, "Error opening file"
!    stop
!  end if
  call suzuki_trotterm2(beta, nwpot, n, a,b,xgrid, valm2)
  do i =1, n
    x = xgrid(i)
    do j =1, n
      xprime= xgrid(j)
      write(4,*) x, xprime,valm2(i,j)
    end do
  end do

  deallocate(xgrid)
  deallocate(valm2)
  close(4)
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
  open(unit=3, file= "test_st_m3_w0.txt", status= "replace", action= "write")
  call suzuki_trotterm3(beta, nwpot, n, a,b, xgrid,valm3)
  do i =1, n
    x = xgrid(i)
    do j =1, n
      xprime= xgrid(j)
      write(3,*) x, xprime,valm3(i,j)
    end do
  end do

  deallocate(xgrid)
  deallocate(valm3)
 close(3)
end subroutine

!_______________________________________________________________________________
subroutine test_gfree()
  implicit none
  real(8)::beta, val, x, xprime, a,b,dx
  real(8)::gf0
  real(8), allocatable::xgrid(:)
  integer:: nwpot, n, i,j
  beta  = 1.0d0
  n     =  100
  a     =  -2.0d0
  b     =  2.0d0
  !dx    =  (b-a)/(n-1)
  allocate(xgrid(n))
  call setup_grid(xgrid, n,a,b,dx)
  open(unit=6, file= "test_gfree_func.txt", status ="replace", action ="write")
  do i =1, n
    x= xgrid(i)
    do j =1, n
      xprime = xgrid(j)
      call gfree(beta, x, xprime,  gf0)
      write(6,*) x, xprime, gf0
    end do
  end do
  close(6)


end subroutine test_gfree


subroutine compare_st_mval_w0()
  implicit none
  real(8)  :: beta, x, xprime
  real(8)  :: gf0
  real(8)  :: a,b
  real(8)  ::

end subroutine


end program
