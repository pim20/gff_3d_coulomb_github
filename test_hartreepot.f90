program test_hartree
implicit none
call test_neon_hartree()

contains
  !-----------------------------------------------------------------------------
  subroutine test_neon_hartree()
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: n, i, ios
  character(len=*), parameter :: infile='electron_ne.txt'
  character(len=*), parameter :: outfile='neon_hartree.txt'
  real(dp), allocatable :: r(:), rho(:), vh(:)
  real(dp) :: dr, pi, EH, Qtot

  pi = 4.0_dp*atan(1.0_dp)

  ! --- Count number of data rows ---
  n = 0
  open(10,file=infile,status='old')
  do
    read(10,*,iostat=ios)
    if (ios /= 0) exit
    n = n + 1
  end do
  close(10)

  if (n <= 0) then
    print *, 'No data in ', infile
    return
  end if

  allocate(r(n), rho(n), vh(n))

  ! --- Read data into arrays ---
  open(11,file=infile,status='old')
  do i = 1, n
    read(11,*) r(i), rho(i)
  end do
  close(11)

  ! --- Call your Hartree routine ---
  call hartree_radial(n, r, rho, vh)
  ! --- Total charge (electron count) ---
 Qtot = 0.0_dp
 do i = 1, n-1
   dr = r(i+1)-r(i)
   Qtot = Qtot + 4.0_dp*pi * 0.5_dp * &
          ( r(i)**2 * rho(i) + r(i+1)**2 * rho(i+1) ) * dr
 end do

  ! --- Write results for plotting ---
  open(12,file=outfile,status='replace')
  write(12,'(A)') '# i r rho V_H'
  do i = 1, n
    write(12,'(I6,3ES20.10)') i, r(i), rho(i), vh(i)
  end do
  close(12)

  ! --- Compute Coulomb (Hartree) energy ---
  EH = 0.0_dp
  do i = 1, n-1
    dr = r(i+1)-r(i)
    EH = EH + 2.0_dp*pi * 0.5_dp * &
         ( r(i)**2 * rho(i) * vh(i) + r(i+1)**2 * rho(i+1) * vh(i+1) ) * dr
  end do
  print '(A,F16.8)', 'Total electronic charge = ', Qtot
  print '(A,F16.8)', 'Coulomb (Hartree) energy = ', EH
  print *, 'Results written to ', outfile

  deallocate(r, rho, vh)
end subroutine test_neon_hartree

  !-----------------------------------------------------------------------------
subroutine test_exp_hartree()
  implicit none
  integer, parameter :: n = 1000
  double precision :: r(n), rho(n), v(n)
  integer :: i
  double precision :: rmax, dr

  rmax = 20.0d0
  dr = rmax/(n-1)
  do i = 1, n
    r(i) = (i-1)*dr
    ! example: exponential density normalized arbitrarily
    rho(i) = exp(-r(i))
  end do

  call hartree_radial(n, r, rho, v)

  ! print a few values
  do i = 1, 10, 1
     print '(I4, 2F12.6)', i, r(i), v(i)
  end do
end subroutine
  !_____________________________________________________________________________


  subroutine test_gaussian_hartree()
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: n = 200
    real(dp), parameter :: rmax = 10.0_dp
    real(dp), parameter :: alpha = 0.5_dp   ! Gaussian width parameter
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp), parameter :: Qtarget = 1.0_dp ! normalize to unit charge
    real(dp) :: dr, rho0
    real(dp), allocatable :: r(:), rho(:), vnum(:), vanal(:)
    integer :: i

    allocate(r(n), rho(n), vnum(n), vanal(n))
    dr = rmax/(n-1)

    ! Normalization: rho(r) = rho0 * exp(-alpha r^2), with total charge Qtarget
    rho0 = Qtarget * alpha**1.5_dp / (pi**1.5_dp)

    ! Build grid and density
    do i = 1, n
      r(i) = (i-1)*dr
      rho(i) = rho0 * exp(-alpha * r(i)**2)
    end do

    ! Call your subroutine to compute numerical Hartree
    call hartree_radial(n, r, rho, vnum)

    ! Analytic Hartree potential: V(r) = Q * erf(sqrt(alpha) * r) / r
    do i = 1, n
      if (r(i) > 1.0d-12) then
        vanal(i) = Qtarget * erf(sqrt(alpha)*r(i)) / r(i)
      else
        vanal(i) = Qtarget * 2.0_dp*sqrt(alpha)/sqrt(pi)  ! r=0 limit
      end if
    end do

    open(unit =1, file ="gaussian_hartree.txt", status = "replace", action ="write")
    ! Print a few sample values
    write(1,*) '   i        r        V_num        V_analytic     diff'
    do i = 1, n
      write(1,*) i, r(i), vnum(i), vanal(i), vnum(i)-vanal(i)
    end do

    ! Maximum error across grid
    print *, 'Max abs error = ', maxval(abs(vnum-vanal))

    deallocate(r, rho, vnum, vanal)
  end subroutine test_gaussian_hartree





subroutine hartree_radial(n, r, rho, vhart)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: r(*)       ! r(1..n) increasing (may have r(1)=0.0)
  double precision, intent(in) :: rho(*)     ! rho(1..n)
  double precision, intent(out) :: vhart(*)  ! output: Hartree potential at r(i)
  integer :: i
  double precision :: pi
  double precision, allocatable :: Q(:)   ! enclosed charge Q(r) = 4*pi * integral_0^r r'^2 rho dr'
  double precision, allocatable :: S(:)   ! outer integral S(r) = 4*pi * integral_r^inf r' rho dr'
  double precision :: dr, term
  pi = 4.0d0*atan(1.0d0)

  if (n .le. 0) return

  allocate(Q(n), S(n))

  ! --- Compute Q cumulatively with trapezoid on r^2 * rho
  Q(1) = 0.0d0
  do i = 1, n-1
    dr = r(i+1) - r(i)
    if (dr .lt. 0.0d0) then
      print *, 'hartree_radial: error - r must be increasing'
      stop
    end if
    ! trapezoid for integral of r'^2 * rho from r(i) to r(i+1)
    term = 0.5d0*( r(i)**2 * rho(i) + r(i+1)**2 * rho(i+1) ) * dr
    Q(i+1) = Q(i) + 4.0d0*pi * term
  end do

  ! --- Compute S from the tail: S(n) = 0 if we assume rho beyond r(n) ~ 0
  S(n) = 0.0d0
  do i = n-1, 1, -1
    dr = r(i+1) - r(i)
    ! trapezoid for integral of r' * rho from r(i) to r(i+1)
    term = 0.5d0*( r(i) * rho(i) + r(i+1) * rho(i+1) ) * dr
    S(i) = S(i+1) + 4.0d0*pi * term
  end do

  ! --- Combine to get V_H
  do i = 1, n
    if (r(i) > 0.0d0) then
      vhart(i) = Q(i)/r(i) + S(i)
    else
      ! r = 0 limit: Q(r)/r -> 0, so just S(1)
      vhart(i) = S(i)
    end if
  end do

  deallocate(Q, S)
end subroutine hartree_radial

end program test_hartree
