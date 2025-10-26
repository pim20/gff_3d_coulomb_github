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

  !call test_suzuki_trotter_m1()
  call test_suzuki_trotter_m2()
   !call test_suzuki_trotter_m3()
  !call test_gfree()
  !call compare_st_mval_w0()
  !call gf_inetgrate()
  !call test_suzuki_trotter_m2_simpson()

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

  open(unit=1, file= "test_st_m1_w0_10ab.txt", status="replace", action="write")
  open(unit=31, file= "test_st_m1_w0_integrated_10ab.txt", status="replace", action="write")
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
  real(8):: start, finish
  !call cpu_time(start)
  allocate(xgrid(n), valm2(n,n))
  open(unit=4, file="test_st_m2_w0_10ab.txt", status="replace", action="write")
  open(unit=41, file="test_st_m2_integrated.txt", status="replace", action="write")
  call suzuki_trotterm2(beta, nwpot, n, a, b, xgrid, valm2)

  do i = 1, n
    x = xgrid(i)
    write(41,*) xgrid(i), xgrid(i), valm2(i,i)
    do j = 1, n
      xprime = xgrid(j)
      write(4,*) x, xprime, valm2(i,j)
    end do
  end do
  !call cpu_time(finish)

  !print *, 'CPU time used (seconds): ', finish - start
  deallocate(xgrid, valm2)
  close(4)
  close(41)
end subroutine test_suzuki_trotter_m2
!===============================================================================
! m = 2 Suzuki–Trotter test, simpson rule
!===============================================================================
subroutine test_suzuki_trotter_m2_simpson()
  use st_params
  implicit none
  real(8) :: x, xprime
  real(8), allocatable :: xgrid(:), valm2(:,:)
  integer :: i,j
  real(8):: start, finish
  call cpu_time(start)
  allocate(xgrid(n), valm2(n,n))
  open(unit=4, file="simp_st_m2_w0_10ab_step_500.txt", status="replace", action="write")
  open(unit=41, file="simp_st_m2_w0_integrated_10ab_step_500.txt", status="replace", action="write")
  call suzuki_trotterm2_simpson(beta, nwpot, n, a, b, xgrid, valm2)

  do i = 1, n
    x = xgrid(i)
    write(41,*) xgrid(i), xgrid(i), valm2(i,i)
    do j = 1, n
      xprime = xgrid(j)
      write(4,*) x, xprime, valm2(i,j)
    end do
  end do
  call cpu_time(finish)

  print *, 'CPU time used (seconds): ', finish - start
  deallocate(xgrid, valm2)
  close(4)
  close(41)
end subroutine test_suzuki_trotter_m2_simpson

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
  open(unit=3, file="test_st_m3_w0_10ab.txt", status="replace", action="write")
  open(unit=51, file= "test_st_m3_w0_integrated_10ab.txt", status="replace", action="write")

  call suzuki_trotterm3(beta, nwpot, n, a, b, xgrid, valm3)

  do i = 1, n
    x = xgrid(i)
    write(51,*) xgrid(i), xgrid(i), valm3(i,i)
    do j = 1, n
      xprime = xgrid(j)
      write(3,*) x, xprime, valm3(i,j)
    end do
  end do

  deallocate(xgrid, valm3)
  close(3)
  close(51)
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

  open(unit=6, file= "test_gfree_func_10ab.txt", status="replace", action="write")
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
  real(8) :: diff1_int, diff2_int, diff3_int
  integer :: i,j

  allocate(xm(n), xp(n))
  allocate(val_m1(n,n), val_m2(n,n), val_m3(n,n), val_gf(n,n))

  ! Read m=1
  open(unit=10, file="test_st_m1_w0_10ab.txt", status="old", action="read")
  do i = 1, n
     do j = 1, n
        read(10,*) xm(i), xp(j), val_m1(i,j)
     end do
  end do
  close(10)

  ! Read m=2
  open(unit=11, file="test_st_m2_w0_10ab.txt", status="old", action="read")
  do i = 1, n
     do j = 1, n
        read(11,*) xm(i), xp(j), val_m2(i,j)
     end do
  end do
  close(11)

  ! Read m=3
  open(unit=12, file="test_st_m3_w0_10ab.txt", status="old", action="read")
  do i = 1, n
     do j = 1, n
        read(12,*) xm(i), xp(j), val_m3(i,j)
     end do
  end do
  close(12)

  ! Read free kernel
  open(unit=13, file="test_gfree_func_10ab.txt", status="old", action="read")
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
  open(unit=15, file="compare_diff_mval_with_gfree_10ab_integrated.txt", status="replace", action="write")
  do i = 1, n
     diff1_int = abs(val_m1(i,i) - val_gf(i,i))
     diff2_int = abs(val_m2(i,i) - val_gf(i,i))
     diff3_int = abs(val_m3(i,i) - val_gf(i,i))
     sum1 = sum1 + diff1_int
     sum2 = sum2 + diff2_int
     sum3 = sum3 + diff3_int
     write(15,*) xm(i), xp(i), xm(i)-xp(i), diff1_int, diff2_int, diff3_int
  end do
  close(15)

  print *, "Sum(|m1 - gfree|) = ", sum1
  print *, "Sum(|m2 - gfree|) = ", sum2
  print *, "Sum(|m3 - gfree|) = ", sum3

  open(unit=16, file="compare_diff_mval_with_gfree_10ab.txt", status="replace", action="write")
  do i = 1, n
    do j =1,n
     diff1 = abs(val_m1(i,j) - val_gf(i,j))
     diff2 = abs(val_m2(i,j) - val_gf(i,j))
     diff3 = abs(val_m3(i,j) - val_gf(i,j))
     !sum1 = sum1 + diff1
     !sum2 = sum2 + diff2
     !sum3 = sum3 + diff3
     write(16,*) xm(i), xp(j), xm(i)-xp(j), diff1, diff2, diff3
   end do
  end do
  close(16)

  deallocate(xm, xp, val_m1, val_m2, val_m3, val_gf)
end subroutine compare_st_mval_w0
!-------------------------------------------------------------------------------
!subroutine for gfree integrating
!-------------------------------------------------------------------------------
subroutine gf_inetgrate()
  use st_params
  implicit none
  real(8):: x
  integer:: i
  real(8), allocatable:: xgrid(:), gfgrid(:)
  real(8):: dx,summ,gf0
  !n=200
  !beta=1.0d0
  !a=-10.0d0
  !b=10.0d0

  allocate(xgrid(n), gfgrid(n))
  call setup_grid(xgrid, n, a, b, dx)
  summ=0.0d0
  do i =1, n
    x = xgrid(i)
    call gfree(beta, x, 2.0d0,  gf0)
    gfgrid(i)=gf0
  end do
  call trapezoid(a,b,gfgrid, summ)
  print *, "summ", summ
  deallocate(xgrid, gfgrid)
  !call gfree(beta, x, 2.0d0,  gf0)
end subroutine
!!------------------------------------------------------------------------------
!for exponential function
!-------------------------------------------------------------------------------







!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
end program test_suzukitrit_1d
