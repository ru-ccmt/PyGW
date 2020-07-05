
subroutine cum_simps(np,n,x,f,g)
  ! !INPUT/OUTPUT PARAMETERS:
  !   np : order of fitting polynomial (in,integer)
  !   n  : number of points (in,integer)
  !   x  : abscissa array (in,real(n))
  !   f  : function array (in,real(n))
  !   g  : integrated function (out,real(n))
  ! !DESCRIPTION:
  !   Calculates the integrals $g(x_i)$ of a function $f$ defined on a set of
  !   points $x_i$ as follows
  !   $$ g(x_i)=\int_0^{x_i} f(x)\,dx. $$
  !   This is performed by piecewise fitting the function to polynomials of order
  !   $n_p-1$ and performing the integrations analytically.
  !
  ! !REVISION HISTORY:
  !   Created May 2002 (JKD)
  !EOP
  !BOC
  implicit none
  ! arguments
  integer, intent(in) :: np
  integer, intent(in) :: n
  real(8), intent(in) :: x(n)
  real(8), intent(in) :: f(n)
  real(8), intent(out) :: g(n)
  ! local variables
  integer i,i0,npo2
  ! automatic arrays
  real(8) c(np)
  ! external functions
  real(8) polynom
  external polynom
  if (n.lt.2) then
     write(*,*)
     write(*,'("Error(fint): n < 2 : ",I8)') n
     write(*,*)
     stop
  end if
  if (np.lt.2) then
     write(*,*)
     write(*,'("Error(fint): np < 2 : ",I8)') np
     write(*,*)
     stop
  end if
  if (n.lt.np) then
     write(*,*)
     write(*,'("Error(fint): n < np : ",2I8)') n,np
     write(*,*)
     stop
  end if
  npo2=np/2
  g(1)=0.d0
  do i=2,n
     if (i.le.npo2) then
        i0=1
     else if (i.gt.n-npo2) then
        i0=n-np+1
     else
        i0=i-npo2
     end if
     g(i)=polynom(-1,np,x(i0),f(i0),c,x(i))+g(i0)
  end do
  return
end subroutine cum_simps
!EOC

real(8) function polynom(m,np,xa,ya,c,x)
  implicit none
  integer, intent(in) :: m
  integer, intent(in) :: np
  real(8), intent(in) :: xa(np)
  real(8), intent(in) :: ya(np)
  real(8), intent(out) :: c(np)
  real(8), intent(in) :: x
  ! !LOCAL VARIABLES:
  integer i,j
  real(8) x1,t1,sum
  if (np.le.0) then
     write(6,*)
     write(6,'("Error(polynom): np <= 0 : ",I8)') np
     write(6,*)
     print *, "ERROR: polynom", "np<=0"
  end if
  if (m.gt.np) then
     polynom=0.d0
     return
  end if
  x1=xa(1)
  ! find the polynomial coefficients in divided differences form
  c(:)=ya(:)
  do i=2,np
     do j=np,i,-1
        c(j)=(c(j)-c(j-1))/(xa(j)-xa(j+1-i))
     end do
  end do
  ! special case m=0
  if (m.eq.0) then
     sum=c(1)
     t1=1.d0
     do i=2,np
        t1=t1*(x-xa(i-1))
        sum=sum+c(i)*t1
     end do
     polynom=sum
     return
  end if
  ! convert to standard form
  do j=1,np-1
     do i=1,np-j
        c(np-i)=c(np-i)-(xa(np-i-j+1)-x1)*c(np-i+1)
     end do
  end do
  if (m.gt.0) then
     ! take the m'th derivative
     do j=1,m
        do i=m+1,np
           c(i)=c(i)*dble(i-j)
        end do
     end do
     t1=c(np)
     do i=np-1,m+1,-1
        t1=t1*(x-x1)+c(i)
     end do
     polynom=t1
  else
     ! find the integral
     t1=c(np)/dble(np)
     do i=np-1,1,-1
        t1=t1*(x-x1)+c(i)/dble(i)
     end do
     polynom=t1*(x-x1)
  end if
  return
end function polynom
!EOC

