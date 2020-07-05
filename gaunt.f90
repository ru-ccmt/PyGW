SUBROUTINE cmp_all_Gaunt(cgcoef, maxj)
  IMPLICIT NONE
  INTEGER, intent(in) :: maxj
  REAL*8, intent(out), dimension( (maxj+1)*(maxj+2)*(maxj+3)*(16*maxj*maxj+29*maxj+10)/60 ) :: cgcoef
  ! External function
  REAL*8 :: Gaunt
  ! Temporaries
  REAL*8, parameter :: pi = 3.14159265358979d0
  INTEGER :: m1, m2, m3, l1, l2, l3, i
  REAL*8  :: c
  !ntot=(maxj+1)*(maxj+2)*(maxj+3)*(16*maxj*maxj+29*maxj+10)/60
  i=0      
  do l1=0,maxj
     do l2=0,l1
        do l3=l1-l2,l1+l2
           do m1=-l1,l1
              do m2=0,l2
                 m3=m1+m2
                 i=i+1
                 if(mod(l1+l2+l3,2).eq.0)then
                    if(iabs(m3).le.l3)then
                       cgcoef(i) = Gaunt(l3,m3,l2,m2,l1,m1)
                    else
                       cgcoef(i)=0.0d0
                    endif
                 else
                    cgcoef(i)=0.0d0
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo
END SUBROUTINE cmp_all_Gaunt

REAL*8 function getcgcoef(l1,l2,l3,m1,m2,cgcoef)
!This subroutine gets the gaunt coefficient $G^{l3,m1+m2}_{l1,l2,m1,m2}$
!from the vector \verb"cgcoef" by:
!\begin{equation}
!G^{l3,m1+m2}_{l1,l2,m1,m2}=\alpha \verb"cgcoef"(i)
!\end{equation}
!calculating i by:
!\begin{equation}
!\begin{aligned}
!i=&\tfrac{1}{60}(16l^2-3l-3)(l+2)(l+1)l+\tfrac{1}{3}ll'(l'+1)(4l'-1)+%
!\tfrac{1}{6}l'(l'-1)(4l'+7)+\\
!&+(2l+1)(l'+1)(L-l+l')+(l'+1)(m+l)+m'+l'+1 
!\end{aligned}
!\end{equation}
!where
!
!\begin{subequations}
!\begin{align}
!l&=l1 & l'&=l2 & m&=m1 & m'&=m2 & \alpha&=1 &\text{ if $l1 \ge l2$} \\
!l&=l2 & l'&=l1 & m&=m2 & m'&=m1 & \alpha&=1 &\text{ if $l1 < l2$} \\
! &    &   &    & m&=m1 & m'&=m2 & \alpha&=1 &\text{ if $m2 \ge 0$} \\
! &    &   &    & m&=-m1 & m'&=-m2 & \alpha&=(-1)^{l+l'-L} &\text{ if $m2 < 0$}
!\end{align}
!\end{subequations}
!
  implicit none
  REAL*8, intent(in) :: cgcoef(*)
  integer(4), intent(in) :: l1, l2, l3, m1, m2
  ! !LOCAL VARIABLES:
  integer(4) :: j1,j2,mj1,mj2
  integer(4) :: par,ing
  integer(4) :: ind1,ind2,ind3,ind4
  real(8) :: fact
  logical :: trcond
  intrinsic mod
  intrinsic abs
  par=mod(abs(l1+l2-l3),2)    
  fact=1.0d0
  trcond=(abs(m1+m2).le.l3).and.(abs(l1-l2).le.l3).and.(l1+l2.ge.l3)
  if(trcond)then
     if(l1.lt.l2)then
        j1=l2
        mj1=m2
        j2=l1
        mj2=m1
     else
        j1=l1
        mj1=m1
        j2=l2
        mj2=m2
     endif
     if(mj2.lt.0)then
        mj2=-mj2
        mj1=-mj1
        fact=(-2.0d0*par+1.0d0)
     endif
     ind1=(16*j1*j1-3*j1-3)*(j1+2)*(j1+1)*j1/60
     ind2=j1*j2*(j2+1)*(4*j2-1)/3
     ind3=j2*(j2-1)*(4*j2+7)/6
     ind4=(2*j1+1)*(j2+1)*(l3-j1+j2)
     ing=ind1+ind2+ind3+ind4+(j2+1)*(mj1+j1)+mj2+j2+1 
     getcgcoef=fact*cgcoef(ing)   
  else
     getcgcoef=0.0d0
  endif
end function getcgcoef



REAL*8 function Gaunt(l1, m1, l2, m2, l3, m3)
  IMPLICIT NONE
  INTEGER, intent(in) :: l1, m1, l2, m2, l3, m3
  REAL*8, PARAMETER   :: pi = 3.14159265358979d0
  REAL*8 :: l1_, l2_, l3_, mm1_, m2_, m3_, zero
  ! function calls
  REAL*8             :: f3j
  REAL*8             :: mone
  l1_ = l1;   l2_ = l2;   l3_ = l3
  mm1_ = -m1; m2_ = m2; m3_ = m3
  zero = 0
  ! Calculates <Y_{l1m1}|Y_{l2m2}|Y_{l3m3}>
  if (l1.LT.0 .OR. l2.LT.0 .OR. l3.LT.0) print *, "Quantum number l must be non-negative!"
  Gaunt = mone(m1)*sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*pi))*f3j(l1_,zero,l2_,zero,l3_,zero)*f3j(l1_,mm1_,l2_,m2_,l3_,m3_)
  return
END function Gaunt

REAL*8 function iFactorial(j)
  IMPLICIT NONE
  INTEGER, intent(in) :: j
  INTEGER :: i
  REAL*8 :: x
  if (j<0) print *, "iFactorial defined only for non-negative numbers!"
  x=1
  iFactorial = x
  if (j.eq.1) return
  DO i=2,j
     x = x*i
  END DO
  iFactorial = x
  return
end function iFactorial

REAL*8 function dFactorial(x)
  IMPLICIT NONE
  REAL*8, intent(in) :: x
  REAL*8, PARAMETER :: spi2 = 0.8862269254527579
  REAL*8 :: y, r
  r=1
  y=x
  DO WHILE(y.gt.1.0)
     r = r * y
     y = y -1.
  ENDDO
  IF (abs(y-0.5).LT.1e-10) r = r*spi2
  dFactorial = r
  return
END function dFactorial

REAL*8 function mone(i)
  INTEGER, intent(in) :: i
  mone = 1 - 2*MOD(abs(i),2)
  return
end function mone

REAL*8 function Delta(j1, j2, j)
  IMPLICIT NONE
  REAL*8, intent(in) :: j1, j2, j
  ! function calls
  REAL*8 :: dFactorial
  Delta = sqrt(dFactorial(j1+j2-j)*dFactorial(j1-j2+j)*dFactorial(-j1+j2+j)/dFactorial(j1+j2+j+1))
  return
END function Delta

REAL*8 function f3j(j1, m1, j2, m2, j3, m3)
  IMPLICIT NONE
  REAL*8, intent(in) :: j1, j2, j3, m1, m2, m3
  INTEGER            :: tmin, tmax, t
  REAL*8             :: sum, v1, v2, dn
  ! function calls
  REAL*8             :: dFactorial
  REAL*8             :: iFactorial
  REAL*8             :: Delta
  REAL*8             :: mone
  f3j=0
  IF (abs(m1+m2+m3) .GT. 1e-10) return
  IF (abs(j1-j2)-1e-14 .GT. j3 .OR. j3 .GT. j1+j2+1e-14) return
  if (abs(m1) .GT. j1 .OR. abs(m2) .GT. j2 .OR. abs(m3) .GT. j3) return
  tmin = INT(max(max(0.0,j2-j3-m1),j1-j3+m2)+1e-14)
  tmax = INT(min(min(j1+j2-j3,j1-m1),j2+m2)+1e-14)
  sum=0
  DO t=tmin, tmax
     v1 = dFactorial(j3-j2+m1+t)*dFactorial(j3-j1-m2+t)
     v2 = dFactorial(j1+j2-j3-t)*dFactorial(j1-m1-t)*dFactorial(j2+m2-t)
     sum = sum + mone(t)/(iFactorial(t)*v1*v2)
  END DO
  dn = dFactorial(j1+m1)*dFactorial(j1-m1)*dFactorial(j2+m2)*dFactorial(j2-m2)*dFactorial(j3+m3)*dFactorial(j3-m3)
  f3j = mone(INT(j1-j2-m3))*Delta(j1,j2,j3)*sqrt(dn)*sum
  return
END function f3j


REAL*8 function ClebschG(j,m,j1,m1,j2,m2)
  IMPLICIT NONE
  REAL*8, intent(in) :: j,m,j1,m1,j2,m2
  INTEGER            :: tmin, tmax, t
  REAL*8             :: sum, v1, v2
  ! function calls
  REAL*8             :: iFactorial
  REAL*8             :: dFactorial
  REAL*8             :: mone
  REAL*8             :: Delta

  ClebschG = 0
  IF (m1+m2 .NE. m) return
  tmin = INT(max(max(0.0,j2-j-m1),j1-j+m2)+1e-14)
  tmax = INT(min(min(j1+j2-j,j1-m1),j2+m2)+1e-14)
  sum=0;
  DO t=tmin, tmax
     v1 = sqrt((2*j+1)*dFactorial(j1+m1)*dFactorial(j1-m1)*dFactorial(j2+m2)*dFactorial(j2-m2)*dFactorial(j+m)*dFactorial(j-m))
     v2 = iFactorial(t)*dFactorial(j1+j2-j-t)*dFactorial(j1-m1-t)*dFactorial(j2+m2-t)*dFactorial(j-j2+m1+t)*dFactorial(j-j1-m2+t)
     sum = sum + mone(t)*v1/v2
  END DO
  ClebschG = sum*Delta(j1,j2,j)
  return
END function ClebschG
