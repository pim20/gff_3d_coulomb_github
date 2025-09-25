program test_suzukitrit_1d
  use suzukit_1d
  use st_params
  implicit none
!-------------------------------------------------------------------------------
! Test file for different Suzuki–Trotter steps (m = 1,2,3)
! nwpot controls the interaction potential (W):
!   nwpot = 1  W = 0
!   nwpot = 2  W = const (example case)
!-------------------------------------------------------------------------------

  call test_suzuki_trotter_m1()
  call test_suzuki_trotter_m2()
   !call test_suzuki_trotter_m3()
  !call test_gfree()
  !call compare_st_mval_w0()

contains
!===============================================================================
! m = 1 Suzuki–Trotter test
!===============================================================================
subroutine test_suzuki_trotter_m1()
  use st_params
  implicit none
  real(8) :: x, xprime
  real(8), allocatable :: xgrid(:), valm1(:,:)
  integer :: i,j

  open(unit=1, file= "test_st_m1_w0_11ab.txt", status="replace", action="write")
  open(unit=31, file= "test_st_m1_w0_integrated_11ab.txt", status="replace", action="write")
  allocate(xgrid(n), valm1(n,n))

  call suzuki_trotterm1(beta, nwpot, n, a, b, xgrid, valm1)

  do i = 1, n
    x = xgrid(i)
    write(31,*) xgrid(i), xgrid(i), valm1(i,i)
    do j = 1, n
      xprime = xgrid(j)
      write(1,*) x, xprime, valm1(i,j)
    end do
  end do

  deallocate(xgrid, valm1)
  close(1)
  close(31)
end subroutine test_suzuki_trotter_m1

!===============================================================================
! m = 2 Suzuki–Trotter test
!===============================================================================
subroutine test_suzuki_trotter_m2()
  use st_params
  implicit none
  real(8) :: x, xprime
  real(8), allocatable :: xgrid(:), valm2(:,:)
  integer :: i,j

  allocate(xgrid(n), valm2(n,n))
  open(unit=4, file="test_st_m2_w0_11ab.txt", status="replace", action="write")
  open(unit=41, file="test_st_m2_w0_integrated_11ab.txt", status="replace", action="write")
  call suzuki_trotterm2(beta, nwpot, n, a, b, xgrid, valm2)

  do i = 1, n
    x = xgrid(i)
    write(41,*) xgrid(i), xgrid(i), valm2(i,i)
    do j = 1, n
      xprime = xgrid(j)
      write(4,*) x, xprime, valm2(i,j)
    end do
  end do

  deallocate(xgrid, valm2)
  close(4)
  close(41)
end subroutine test_suzuki_trotter_m2

!===============================================================================
! m = 3 Suzuki–Trotter test
!===============================================================================
subroutine test_suzuki_trotter_m3()
  use st_params
  implicit none
  real(8) :: x, xprime
  real(8), allocatable :: xgrid(:), valm3(:,:)
  integer :: i,j

  allocate(xgrid(n), valm3(n,n))
  open(unit=3, file="test_st_m3_w0_11ab.txt", status="replace", action="write")

  call suzuki_trotterm3(beta, nwpot, n, a, b, xgrid, valm3)

  do i = 1, n
    x = xgrid(i)
    do j = 1, n
      xprime = xgrid(j)
      write(3,*) x, xprime, valm3(i,j)
    end do
  end do

  deallocate(xgrid, valm3)
  close(3)
end subroutine test_suzuki_trotter_m3

!===============================================================================
! Free propagator (analytic kernel) test
!===============================================================================
subroutine test_gfree()
  use st_params
  implicit none
  real(8) :: x, xprime, gf0, dx
  real(8), allocatable :: xgrid(:)
  integer :: i,j

  allocate(xgrid(n))
  call setup_grid(xgrid, n, a, b, dx)

  open(unit=6, file= "test_gfree_func_11ab.txt", status="replace", action="write")
  do i = 1, n
    x = xgrid(i)
    do j = 1, n
      xprime = xgrid(j)
      call gfree(beta, x, xprime, gf0)
      write(6,*) x, xprime, gf0
    end do
  end do
  close(6)

  deallocate(xgrid)
end subroutine test_gfree

!===============================================================================
! Compare ST approximations (m=1,2,3) with analytic free propagator
!===============================================================================
subroutine compare_st_mval_w0()
  use st_params
  implicit none
  real(8), allocatable :: xm(:), xp(:)
  real(8), allocatable :: val_m1(:,:), val_m2(:,:), val_m3(:,:), val_gf(:,:)
  real(8) :: diff1, diff2, diff3, sum1, sum2, sum3
  integer :: i,j

  allocate(xm(n), xp(n))
  allocate(val_m1(n,n), val_m2(n,n), val_m3(n,n), val_gf(n,n))

  ! Read m=1
  open(unit=10, file="test_st_m1_w0_11ab.txt", status="old", action="read")
  do i = 1, n
     do j = 1, n
        read(10,*) xm(i), xp(j), val_m1(i,j)
     end do
  end do
  close(10)

  ! Read m=2
  open(unit=11, file="test_st_m2_w0_11ab.txt", status="old", action="read")
  do i = 1, n
     do j = 1, n
        read(11,*) xm(i), xp(j), val_m2(i,j)
     end do
  end do
  close(11)

  ! Read m=3
  open(unit=12, file="test_st_m3_w0_11ab.txt", status="old", action="read")
  do i = 1, n
     do j = 1, n
        read(12,*) xm(i), xp(j), val_m3(i,j)
     end do
  end do
  close(12)

  ! Read free kernel
  open(unit=13, file="test_gfree_func_11ab.txt", status="old", action="read")
  do i = 1, n
     do j = 1, n
        read(13,*) xm(i), xp(j), val_gf(i,j)
     end do
  end do
  close(13)

  ! Differences along diagonal
  sum1 = 0.0d0
  sum2 = 0.0d0
  sum3 = 0.0d0
  open(unit=15, file="compare_diff_mval_with_gfree_11ab.txt", status="replace", action="write")
  do i = 1, n
     diff1 = abs(val_m1(i,i) - val_gf(i,i))
     diff2 = abs(val_m2(i,i) - val_gf(i,i))
     diff3 = abs(val_m3(i,i) - val_gf(i,i))
     sum1 = sum1 + diff1
     sum2 = sum2 + diff2
     sum3 = sum3 + diff3
     write(15,*) xm(i), xp(i), diff1, diff2, diff3
  end do
  close(15)

  print *, "Sum(|m1 - gfree|) = ", sum1
  print *, "Sum(|m2 - gfree|) = ", sum2
  print *, "Sum(|m3 - gfree|) = ", sum3

  deallocate(xm, xp, val_m1, val_m2, val_m3, val_gf)
end subroutine compare_st_mval_w0

end program test_suzukitrit_1d
