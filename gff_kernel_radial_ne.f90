!_______________________________________________________________________________
! ============================================================================
! kernel_radial.f90  (standalone, no external libs)
! Compute K_beta(r, r') using symmetric Suzuki–Trotter.
! - Reads radial density file: rho_r.dat  (two columns: r  rho(r))
! - Builds Hartree potential V_H(r) from rho(r)_r
! - Propagates u=r*psi with Crank–Nicolson for kinetic step
! - Outputs a kernel column: K_beta(r, r') for a chosen r' grid index
!
! ============================================================================

program kernel_radial
  implicit none

contains
!_____________________________Density inputs____________________________________
!-------------------------------------------------------------------------------
  subroutine read_rho_radial(fname, rgrid, rho)
    implicit none
    character(len=*), intent(in) :: fname
    real(8), intent(in)  :: rgrid(:)
    real(8), intent(out) :: rho(:)
    integer :: ios, i
    real(8) :: rr, rv
    rho = 0.0d0
    open(10, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'WARNING: ', trim(fname),' not found. Using mock Gaussian density.'
      do i = 1, size(rgrid)
        rho(i) = 10.0d0 * exp(- (rgrid(i)/1.0d0)**2)
      end do
      return
    end if

     do
       read(10, *, iostat = ios) rr, rv
       if (ios /= 0) exit
       call grid_nearest(rr,rv,rgrid,rho)
     end do
    close(10)
  end subroutine read_rho_radial
!-------------------------------------------------------------------------------
  subroutine grid_nearest(rr, rv, rgrid, rho)
  implicit none
  real(8), intent(in):: rr, rv
  real(8), intent(in):: rgrid(:)
  real(8), intent(out)::rho(:)
  integer:: i
  i = int(rr / (rgrid(2)-rgrid(1))) + 1
  i = max(1, min(i, size(rgrid)))
  rho(i) = rv

  end subroutine
!____________________________Potential__________________________________________
!-------------------------------------------------------------------------------
subroutine compute_potential_hartree(r, rho, VH)
implicit none
real(8), intent(in):: r(:), rho(:)
real(8), intent(out):: VH(:)




end subroutine
!_______________________________________________________________________________
subroutine errorfunctionbyr(r, erfr)
implicit none
real(8), intent(in):: r
real(8), intent(out):: erfr
!give errorfunction
erfr = erf(r)/r

end subroutine

subroutine free_kernel(r, rprime, gffree)
implicit none
real(8), intent(in):: r, rprime
real(8), intent(out):: gffree

end subroutine

subroutine gfree(x,xprime,beta,g0)
  implicit none
  real(8), intent(in)::x,xprime,beta
  real(8), intent(out)::g0
  real(8):: a0, b0,pi
  !m=1, hbar=1
  pi= 4.0d0*atan(1.0d0)
  a0=sqrt(1.0d0/(2.0d0*beta*pi))
  b0= 1.0d0/(2.0d0*beta)
  g0 = a0*exp(-b0*(x-xprime)**2)
end subroutine gfree

end program kernel_radial
