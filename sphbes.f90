
subroutine sphbes(n,x,sj,sjp)
  !   Calculates the spherical bessel function $j_l(x)$ and its first derivative $j'_l(x)$ 
  !   at the point $x\le 0$, for l=0,1,...n       
  implicit none
  integer, intent(in) :: n    ! Maximum l to be calculated
  real(8), intent(in) :: x    ! The argument at which the bessel functions are calculated
  real(8), intent(out) :: sj(*) ! values of spherical bessel function
  real(8), intent(out) :: sjp(*)! values of first derivative of the spherical bessel function
  !
  integer(4) :: l
  real(8) :: factor  ! Multiplying factor between the spherical bessel func and the regular bessel func of half integer order
  real(8) :: rj(n+1) ! Array containing the regular bessel function of order i-1/2
  real(8) :: rjp(n+1)! Array containing the first derivative of the regular bessel function of order i-1/2
  real(8), parameter :: rtpio2 = 1.25331413731550025120788d+0
  external bessj
  intrinsic sqrt      

  if(n.lt.0.or.x.lt.0.)then
     write(6,*) "ERROR in sphbes:bad arguments"
     write(6,*)'  input: x = ',x
     write(6,*)'  input: n = ',n
     stop 'bad arguments in sphbes'
  endif
      
  if(x.gt.0.0d0)then
     call bessj(x,n,rj,rjp)
     
     factor=rtpio2/sqrt(x)
     do l=1,n+1
        sj(l)=factor*rj(l)
        sjp(l)=factor*rjp(l)-sj(l)/(2.0d0*x)
     enddo
  else
     do l=1,n+1
        sj(l)=0.0d0
        sjp(l)=0.0d0
     enddo
     sj(1)=1.0d0
     sjp(2)=1.0d0/3.0d0 
  endif
end subroutine sphbes

subroutine bessj(x,nl,rj,rjp)
  !Returns the regular Bessel functions of first kind
  !$\texttt{rj(l)}=J_{l-\frac{1}{2}}$ and their derivatives
  !$\texttt{rjp(l)}=J'_{l-\frac{1}{2}}$, for positive \texttt{x} and for
  !$1\le l \le nl+1$, using Steed's method.
  implicit none
  real(8), intent(in) :: x
  integer, intent(in) :: nl
  real(8), intent(out) :: rj(*)            
  real(8), intent(out) :: rjp(*)            
  !    
  integer :: i, isign, l 
  real(8) :: xnu, b, c, d, del, f, fact, gam, h, p, q, rjl, rjl1, rjmu, rjp1, rjpl, rjtemp, w, xi, xi2, xmu, xmu2
  integer, parameter :: maxit=10000
  real(8), parameter :: eps=1.e-16
  real(8), parameter :: fpmin=1.e-30
  real(8), parameter :: pi=3.141592653589793d+0
  ! Original subroutine: bessjy.for (c) copr. 1986-92 numerical recipes software &124i..
  if(x.le.0..or.nl.lt.0) stop 'bad arguments in bessj'
  xnu=dble(nl)+0.5d0
  xmu=0.5d0
  xmu2=xmu*xmu
  xi=1.d0/x
  xi2=2.d0*xi
  !     The Wronskian
  w=xi2/pi
  !     Evaluate the continued fraction expansion for J'_(nl-1/2)/J_(nl-1/2)
  !     by the modified Lentz's method. isign keeps track of sign changes in
  !     the denominator
  isign=1
  h=xnu*xi
  if(h.lt.fpmin)h=fpmin
  b=xi2*xnu
  d=0.d0
  c=h
  do i=1,maxit
     b=b+xi2
     d=b-d
     if(abs(d).lt.fpmin) d=fpmin
     c=b-1.d0/c
     if(abs(c).lt.fpmin) c=fpmin
     d=1.d0/d
     del=c*d
     h=del*h
     if(d.lt.0.d0) isign=-isign
     if(abs(del-1.d0).lt.eps) goto 1
  enddo
!11   continue
  stop 'x too large in bessjy; try asymptotic expansion'
1 continue
  !     Initialize J and J' for downward recurrence
  rjl=isign*fpmin
  rjpl=h*rjl
  !     Store values for later rescaling      
  rjl1=rjl
  rjp1=rjpl
  !     Downward recurrence (unnormalized)
  fact=xnu*xi
  do l=nl,0,-1
     rjtemp=fact*rjl+rjpl
     fact=fact-xi
     rjpl=fact*rjtemp-rjl
     rjl=rjtemp
  enddo
  if(rjl.eq.0.d0)rjl=eps
  f=rjpl/rjl
  !     Equation 6.7.3 Numerical Recipies in Fortran
  p=-.5d0*xi
  q=1.d0
  !     Equation 6.7.6 Numerical Recipies in Fortran
  gam=(p-f)/q
  !     Equation 6.7.7 Numerical Recipies in Fortran
  rjmu=sqrt(w/((p-f)*gam+q))
  rjmu=sign(rjmu,rjl)
  !     Scale original J and J'
  fact=rjmu/rjl
  rj(nl+1)=rjl1*fact
  rjp(nl+1)=rjp1*fact
  fact=xnu*xi
  !     Downward recurence
  do l=nl,1,-1
     rj(l)=fact*rj(l+1)+rjp(l+1)
     fact=fact-xi
     rjp(l)=fact*rj(l)-rj(l+1)
  enddo ! l
end subroutine bessj


subroutine ylm(y,v,lmax)
  !  This subroutine calculates spherical harmonics up to l=lmax for a
  ! given vector in cartesian coordinates
  !
  !  The spherical harmonics (Condon and Shortley convention)
  !  $Y_{0,0},Y_{1,-1},Y_{1,0},Y_{1,1},Y_{2,-2} ...Y_{LMAX,LMAX}$
  !           for vector V (given in Cartesian coordinates)
  !           are calculated. In the Condon Shortley convention the
  !           spherical harmonics are defined as
  !
  !\begin{equation}                            
  ! Y_{l,m}=(-1)^m \sqrt{\frac{1}{2\pi}} P_l^m(cos(\theta))e^{im\phi}
  !\end{equation}
  !
  !\noindent
  !where $P_l^m(cos(\theta))$ is the normalized Associated Legendre
  !           function. Thus,
  !\begin{equation}                            
  ! Y_{l,-m}=(-1)^m Y^*_{l,m}
  !\end{equation}
  !
  !The output is writen to the vector Y(:) such that
  ! Y(k)$=Y_{l,m}$ with $k=l^2+l+m+1$, thus,
  ! Y(1)$=Y_{0,0}$, Y(2)$=Y_{1,-1}$, Y(3)$=Y_{1,0}$, Y(4)$=Y_{1,1}$...
  !Y(lmax*lmax+1)$=Y_{lmax,-lmax}... $Y((lmax+1)*(lmax+1))$=Y_{lmax,lmax}$
  !
  !   The basic algorithm used to calculate the spherical
  !   harmonics for vector V is as follows:
  !
  !\begin{subequations}
  !\begin{align}
  ! Y_{0,0}=&\sqrt{\frac{1}{4\pi}}\\
  ! Y_{1,0}=&\sqrt{\frac{3}{4\pi}}cos(\theta)\\
  ! Y_{1,1}=&-\sqrt{\frac{3}{8\pi}}sin(\theta)e^{i\phi}\\
  ! Y_{1,-1}=&-Y_{1,1}\\
  ! Y_{l,l}=&-\sqrt{\frac{2l+1}{2l}}sin(\theta)e^{i\phi}Y_{l-1,l-1}\\
  ! Y_{l,m}=&\sqrt{\frac{(2l-1)(2l+1)}{(l-m)(l+m)}}%
  !cos(\theta)Y_{l-1,m}-\sqrt{\frac{(l-1+m)(l-1-m)(2l+1)}{(2l-3)(l-m)(l+m)}}%
  !Y_{l-2,m}
  !\end{align}
  !\end{subequations}
  !
  implicit none
  integer,intent(in):: lmax  ! maximum l for which the spherical  harmonics are calculated
  real(8),intent(in) :: v(3) ! vector (cartesian coordinates) for which  the spherical harmonics are calculated
  complex*16,intent(out) :: y((lmax+1)*(lmax+1)) ! the values of the spherical harmonics dimension (lmax+1)^2
  !
  integer ::  i2l, i4l2, index, index2, l, m, msign
  real(8) ::  a, b, c, ab, abc, abmax, abcmax
  real(8) ::  d4ll1c, d2l13, pi
  real(8) ::  costh, sinth, cosph, sinph
  real(8) ::  temp1, temp2, temp3
  real(8) ::  yllr, yll1r, yl1l1r, ylmr
  real(8) ::  ylli, yll1i, yl1l1i, ylmi
  pi = (4.0d+0)*atan(1.0d+0)
  !        y(0,0)
  yllr = 1.0d+0/sqrt(4.0d+0*pi)
  ylli = 0.0d+0
  y(1) = cmplx(yllr,ylli,8)
  !        continue only if spherical harmonics for (l .gt. 0) are desired
  if (lmax .le. 0) return 
  !        calculate sin(phi), cos(phi), sin(theta), cos(theta)
  !        theta, phi ... polar angles of vector v
  abmax  = max(abs(v(1)),abs(v(2)))
  if (abmax .gt. 0.0d+0) then
     a = v(1)/abmax
     b = v(2)/abmax
     ab = sqrt(a*a+b*b)
     cosph = a/ab
     sinph = b/ab
  else
     cosph = 1.0d+0
     sinph = 0.0d+0
  endif
  abcmax = max(abmax,abs(v(3)))
  if (abcmax .gt. 0.0d+0) then
     a = v(1)/abcmax
     b = v(2)/abcmax
     c = v(3)/abcmax
     ab = a*a + b*b
     abc = sqrt(ab + c*c)
     costh = c/abc
     sinth = sqrt(ab)/abc
  else
     costh = 1.0d+0
     sinth = 0.0d+0
  endif
  !        y(1,0)
  y(3) = cmplx(sqrt(3.0d+0)*yllr*costh,0.0d+0,8)
  !        y(1,1) ( = -conjg(y(1,-1)))
  temp1 = -sqrt(1.5d+0)*yllr*sinth
  y(4) = cmplx(temp1*cosph,temp1*sinph,8)
  y(2) = -conjg(y(4))
  !
  do l = 2, lmax
     index  = l*l+1
     index2 = index + 2*l
     msign  = 1 - 2*mod(l,2)
     !        yll = y(l,l) = f(y(l-1,l-1)) ... formula 1
     yl1l1r = dble(y(index-1))
     yl1l1i = aimag(y(index-1))
     temp1 = -sqrt(dble(2*l+1)/dble(2*l))*sinth
     yllr = temp1*(cosph*yl1l1r - sinph*yl1l1i)
     ylli = temp1*(cosph*yl1l1i + sinph*yl1l1r)
     y(index2) = cmplx(yllr,ylli,8)
     y(index)  = msign*conjg(y(index2))
     index2 = index2 - 1
     index  = index  + 1
     !        yll1 = y(l,l-1) = f(y(l-1,l-1)) ... formula 2
     !               (the coefficient for y(l-2,l-1) in formula 2 is zero)
     temp2 = sqrt(dble(2*l+1))*costh
     yll1r = temp2*yl1l1r
     yll1i = temp2*yl1l1i
     y(index2) = cmplx(yll1r,yll1i,8)
     y(index)  = -msign*conjg(y(index2))
     index2 = index2 - 1
     index  = index  + 1
     !
     i4l2 = index2 - 4*l + 2
     i2l  = index2 - 2*l
     d4ll1c = costh*sqrt(dble(4*l*l-1))
     d2l13  = -sqrt(dble(2*l+1)/dble(2*l-3))
     do m = l - 2, 0, -1
        !        ylm = y(l,m) = f(y(l-2,m),y(l-1,m)) ... formula 2
        temp1 = 1.0d+0/sqrt(dble((l+m)*(l-m)))
        temp2 = d4ll1c*temp1
        temp3 = d2l13*sqrt(dble((l+m-1)*(l-m-1)))*temp1
        ylmr = temp2*dble(y(i2l))  + temp3*dble(y(i4l2))
        ylmi = temp2*aimag(y(i2l)) + temp3*aimag(y(i4l2))
        y(index2) = cmplx(ylmr,ylmi,8)
        y(index)  = msign*conjg(y(index2))
        !
        msign  = -msign
        index2 = index2 - 1
        index  = index  + 1
        i4l2   = i4l2   - 1
        i2l    = i2l    - 1
     enddo  ! m
  enddo ! l
end subroutine ylm
