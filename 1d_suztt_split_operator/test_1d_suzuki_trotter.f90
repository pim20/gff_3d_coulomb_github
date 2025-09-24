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
!call test_gfree()
call compare_st_mval_w0()

contains
!_______________________________________________________________________________
!This subroutine test e^{-beta(H0+W)} suzuki trotter m=1 (first order)
!W is taken zero and constant represented by (nwpot)
subroutine test_suzuki_trotter_m1()
  implicit none
  real(8)::beta, val, x, xprime, a,b,dx
  real(8), allocatable::xgrid(:), valm1(:,:)
  integer:: nwpot, n, i,j
  nwpot=1
  beta = 1.0d0
  n=100
  a=-2.0d0
  b=2.0d0
  open(unit=1, file= "test_st_m1_w0.txt", status= "replace", action= "write")
  allocate(xgrid(n))
  allocate(valm1(n,n))
  !call setup_grid(xgrid, n,a,b,dx)
  call suzuki_trotterm1(beta,nwpot, n, a,b, xgrid, valm1)
  do i =1, n
    x =xgrid(i)
    do j =1,n
      xprime = xgrid(j)
      write(1,*) x, xprime, valm1(i,j)
    end do
  end do
  !call setupgrid_call()
  deallocate(xgrid)
  deallocate(valm1)
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
  real(8)  :: beta
  integer  :: n
  real(8), allocatable:: xm(:), xp(:)
  real(8), allocatable:: val_m1(:,:), val_m2(:,:), val_m3(:,:), val_gf(:,:)
  real(8):: diff1, diff2, diff3,sum1, sum2, sum3
  real(8)  :: gf0
  real(8)  :: a,b
  integer:: i,j
  ! Same grid size as in test runs
  n = 100

  allocate(xm(n), xp(n))
  allocate(val_m1(n,n), val_m2(n,n), val_m3(n,n), val_gf(n,n))
  !real(8)  ::
  ! read the file m=1
  open(unit=10, file="test_st_m1_w0.txt", status="old", action="read")
  do i = 1, n
     do j = 1, n
        read(10,*) xm(i), xp(j), val_m1(i,j)
     end do
  end do
  close(10)

  !read m=2
  open(unit=11, file= "test_st_m2_w0.txt", status ="old", action= "read")
  do i= 1, n
    do j=1, n
      read(11,*) xm(i), xp(j), val_m2(i,j)
    end do
  end do
  close(11)

  ! Read m=3
    open(unit=12, file="test_st_m3_w0.txt", status="old", action="read")
    do i = 1, n
       do j = 1, n
          read(12,*) xm(i), xp(j), val_m3(i,j)
       end do
    end do
    close(12)

    ! Read free kernel
    open(unit=13, file="test_gfree_func.txt", status="old", action="read")
    do i = 1, n
       do j = 1, n
          read(13,*) xm(i), xp(j), val_gf(i,j)
       end do
    end do
    close(13)
    !now for the diff

    sum1 =0.0d0
    sum2 =0.0d0
    sum3 =0.0d0
    open(unit =15, file= "compare_diff_mval_with_gfree.txt", status="replace", action="write")
    do i = 1, n
     do j = 1, n
        diff1 = abs(val_m1(i,j) - val_gf(i,j))
        diff2 = abs(val_m2(i,j) - val_gf(i,j))
        diff3 = abs(val_m3(i,j) - val_gf(i,j))

        sum1 = sum1+diff1
        sum2 = sum2+diff2
        sum3 = sum3+diff3
        write(15,*) xm(i), xp(j), diff1, diff2, diff3
     end do
   end do

    close(15)

    print *, "Sum(|m1 - gfree|) = ", sum1
    print *, "Sum(|m2 - gfree|) = ", sum2
    print *, "Sum(|m3 - gfree|) = ", sum3
    deallocate(xm, xp, val_m1, val_m2, val_m3, val_gf)


end subroutine


end program
