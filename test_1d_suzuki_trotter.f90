program test_suzukitrit_1d
use suzukit_1d
implicit none
real(8)::beta, val, x, xprime, a,b
real(8), allocatable::xgrid(:)
integer:: nwpot, n, i,j
nwpot=1
beta = 1.0d0
n=100
a=-2.0d0
b=2.0d0
open(unit=1, file= "testsuzukit1d.txt", status= "replace", action= "write")
allocate(xgrid(n))
call setup_grid(xgrid, n,a,b)
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
end program
