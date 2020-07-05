real*8 function polynom(m,np,xa,ya,c,x)
  implicit none
  integer, intent(in) :: m     ! which derivative 1,2,.... Here 0 means integral
  integer, intent(in) :: np    ! how many points to use for the polynomial fit
  real(8), intent(in) :: xa(np)! 
  real(8), intent(in) :: ya(np)! 
  real(8), intent(out):: c(np) ! temporary array
  real(8), intent(in) :: x     ! the value of x at the point we need derivative
  ! local variables
  integer i,j
  real(8) x1,t1,sum
  if (np <= 0) then
     write(6,*)
     write(6,'("Error(polynom): np <= 0 : ",I8)') np
     write(6,*)
     call exit(1)
  end if
  if (m > np) then
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

subroutine derv(ucpl,ucl,rr,npt,np)
  REAL*8, intent(out):: ucpl(npt)
  REAL*8, intent(in) :: ucl(npt), rr(npt)
  INTEGER, intent(in):: npt
  INTEGER, intent(in):: np  ! default np=4
  !f2py integer optional, intent(in) :: np=4
  real(8), external :: polynom
  ! locals
  INTEGER :: ir, ir0
  INTEGER :: np_2
  REAL*8  :: c(np)
  np_2 = np/2
  ucpl(1)=0.0d0
  do ir=2,npt
     if(ir <= np_2) then
        ir0 = 1
     elseif(ir > npt-np_2)then
        ir0 = npt-np+1
     else
        ir0 = ir-np_2
     endif
     ucpl(ir) = polynom(1,np,rr(ir0),ucl(ir0),c,rr(ir))
  enddo
end subroutine derv

