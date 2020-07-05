!  Recursion for Pade coefficient (J.Serene)
!****************************
SUBROUTINE Padecof(Pt,gn,zn,nn)
!****************************
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Pt(nn)
  COMPLEX*16, intent(in)  :: gn(nn),zn(nn)
  INTEGER, intent(in)     :: nn
  ! local
  COMPLEX*16  :: p(nn,nn)
  INTEGER :: i, j
  p=0.0
  do j=1,nn
     p(1,j) = gn(j)
  enddo
  do j=2,nn
     do i=2,j
        p(i,j)=(p(i-1,i-1)-p(i-1,j))/(zn(j)-zn(i-1))/p(i-1,j)
     enddo
  enddo
  do j=1,nn
     Pt(j)=p(j,j)
  enddo
  return
end SUBROUTINE Padecof

!  Calculation of a Green's function for a given pade-coeff-p(i,j)
!  on the real axis e=e+i0
!*****************************
SUBROUTINE PadeG(Gx,x,zn,Pt,nn,m)
!*****************************
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Gx(m)
  COMPLEX*16, intent(in)  :: x(m), zn(nn),Pt(nn)
  INTEGER, intent(in)     :: nn, m
  ! locals
  COMPLEX*16 :: aw(0:nn),bw(0:nn)
  INTEGER :: i, j
  aw(0)=(0.d0,0.d0)
  aw(1)=Pt(1)
  bw(0)=(1.d0,0.d0)
  bw(1)=(1.d0,0.d0)
  do j=1,m
     do i=1,nn-1
        aw(i+1)=aw(i)+(x(j)-zn(i))*Pt(i+1)*aw(i-1)
        bw(i+1)=bw(i)+(x(j)-zn(i))*Pt(i+1)*bw(i-1)
     enddo
     Gx(j)=aw(nn)/bw(nn)
  enddo
  return
end SUBROUTINE PadeG

!subroutine setpatrd(npar,z,f,a)
!  ! a == Pt, f==gn, z==zn
!  !
!  ! this subroutine calculate coeffients in Pade's approximation 
!  ! using Thiele's reciprocal difference method as described in 
!  !   H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977) 
!  ! or
!  !   K.-H. Lee and K. J. Chang, Phys. Rev. B 54 R8285 (1996) 
!  !
!  ! Input:
!  !  npar       --- order of the Pade's approximation, the number of poles is equal to npar/2
!  !  z(1..npar) --- complex points, on which, the values of the function to be fitted are given
!  !  f(1..npar)
!  !  Output: 
!  !  a(1..npar) --- coefficients of Pade's approximation
!  ! 
!  implicit none 
!  integer(4),intent(in)::npar
!  complex(8),intent(in)::z(npar),f(npar)
!  complex(8),intent(out)::a(npar) 
!
!  integer(4)::n,p, i
!  complex(8)::g(npar,npar)
!  g=0.d0
!  g(1:npar,1)=f(1:npar) 
!
!  do i=1,npar
!     write(6,'(I2,1x,2F21.16,1x,2F20.16)') i, z(i), f(i)
!  enddo
!  
!  do p=2,npar
!     do n=p,npar
!        g(n,p)=(g(p-1,p-1)-g(n,p-1))/((z(n)-z(p-1))*g(n,p-1))
!     enddo
!
!     do i=1,npar
!        write(6,'(I2,1x,2F20.16)') i, g(i,p)
!     enddo
!     
!  enddo
!  do n=1,npar
!     a(n)=g(n,n)
!  enddo
!
!  do i=1,npar
!     write(6,'(I2,1x,2F21.16,1x,2F20.16)') i, z(i), a(i)
!  enddo
!  
!end subroutine setpatrd

subroutine init_c(x,y,c,m)
  ! The algorithm is described in PRB 61, 5147 (2000). It is equivalent to pade
  ! but using matrix inversion method
  ! This subroutine calculates the values of the \texttt{2n+2} parameters $c_k$ 
  ! of a function of the form:
  !
  !\begin{equation}\label{init_c-1}     
  !f(x,\{c\})=\frac{P(x,\{c\})}{Q(x,\{c\})}
  !\end{equation}
  !
  !where 
  !
  !\begin{equation}\label{init_c-2}     
  !P(x,\{c\})=\sum\limits_{k=1}^{n+1}{c_{k}x^{k-1}},
  !\end{equation}
  !
  !\begin{equation}\label{init_c-3}     
  !Q(x,\{c\})=1+\sum\limits_{k=n+2}^{2n+2}{c_k x^{k-n-1}}
  !\end{equation}
  !
  !by adjusting them to a set of \texttt{2n+2} $(x,y)$-pairs.
  !
  implicit none
  integer,    intent(in) :: m
  complex*16, intent(in) :: x(m)
  complex*16, intent(in) :: y(m)
  complex*16,  intent(out):: c(m)
  !
  integer(4) :: i, j, n, info
  integer(4) :: ip(1:m)
  complex*16 :: xtn
  complex*16 :: amat(m,m)
  complex*16, parameter :: czero = (0.0d0,0.0d0)      
  complex*16, parameter :: cone  = (1.0d0,0.0d0)      
  external zgetrf
  external zgetrs
  !
  n=m/2-1
  c(1:m) = y(1:m)
  !
  do i=1,m
     xtn = cone
     do j=1,n+1
        amat(i,j) = xtn              ! amat[i,j] = x[i]^(j-1)          j=1,..m/2-1
        xtn = xtn * x(i)   
        amat(i,j+n+1) = -y(i)*xtn    ! amat[i,j+m/2] = -y[i] * x[i]^j  j=1,..m/2-1
     enddo
  enddo
  call zgetrf(m,m,amat,m,ip,info)   ! creates LU factorization of matrix amat
  if (info.ne.0) then
     write(6,*) 'acont: zgetrf, info =',info
     stop
  endif
  call zgetrs('n',m,1,amat,m,ip,c,m,info) ! c <= amat^{-1} * c
  if (info.ne.0) then
     write(6,*) 'acont: zgetrs, info =',info
     stop
  endif
  return
end subroutine init_c


subroutine nllsq(x,y,ndata,a,ma,chisq)
  ! This subroutine fits the set for data points 
  !\texttt{x(1:ndata)}, \texttt{y(1:ndata)}, and a nonlinear function
  !dependent on \texttt{ma} coefficients \texttt{a(1:ma)} using the
  !Levenberg-Marquardt method.
  !
  implicit none
  integer,    intent(in)    :: ndata     ! Number of datapoints to be fitted
  real(8),    intent(in)    :: x(ndata)  ! abscisas of the datapoints
  complex*16, intent(in)    :: y(ndata)  ! data value of the datapoint
  integer,    intent(in)    :: ma        ! Number of coefficients
  real(8),    intent(inout) :: a(ma)
  real(8),    intent(out) :: chisq
  !
  integer :: it
  real(8) :: alambda, alambdaold, chisqold, deltach
  real(8) :: covar(ma,ma), alpha(ma,ma)
  logical :: converg
  real(8), parameter :: tol = 1.0d-10
  external mrqmin
  external stdesc
  !
  converg = .false.
  call stdesc(x,y,ndata,a,ma,chisqold)
  alambda    = -0.1d0
  alambdaold =  1.0d0
  it=0
  do while (.not.converg)
     it=it+1
     call mrqmin(x,y,ndata,a,ma,covar,alpha,ma,chisq,alambda)
     deltach=abs(chisq-chisqold)
     if(alambda.lt.alambdaold)then
        converg = (deltach.lt.tol).or.(it.gt.1000)
        chisqold=chisq
     else  
        converg = (it.gt.1000)
     endif
     alambdaold=alambda
  enddo
  if (it.gt.1000) write(*,*)'WARNING: chisq did not converge'
  alambda = 0.0d0
  return
end subroutine nllsq

subroutine ratfun(x,a,y,dyda,hess,n)
  ! Given the abscissa \texttt{x}, and the parameters \texttt{a(1:4n+4)}
  !this subroutine returns in \texttt{y} the value of the function:
  !
  !\begin{equation}\label{ratfun-1}     
  !f(x,\{c\})=\frac{P(x,\{c\})}{Q(x,\{c\})}
  !\end{equation}
  !
  !where 
  !\begin{equation}\label{ratfun-2}     
  !P(x,\{c\})=\sum\limits_{k=1}^{n+1}{c_{k}x^{k-1}},
  !\end{equation}
  !
  !\begin{equation}\label{ratfun-3}     
  !Q(x,\{c\})=1+\sum\limits_{k=n+2}^{2n+2}{c_k x^{k-n-1}}
  !\end{equation}
  !
  !and $c_k=\texttt{a(k)}+i\texttt{a(k+2n+2)}$. In \texttt{dyda} the partial
  !derivatives of the function with respect to the parameters:
  !
  !\begin{equation}\label{ratfun-4}     
  !\frac{\partial f(x,\{c\})}{\partial c_j}=\left\{
  !\begin{array}{lc}
  !{\displaystyle \frac{x^{j-1}}{Q(x,\{c\})}}& 1 \leq j \leq n+1 \vspace{3mm}\\
  !{\displaystyle -\frac{P(x,\{c\})}{[Q(x,\{c\})]^2}x^{j-n-1}}& n+2 %
  !\leq j \leq 2n+2 \\
  !\end{array}
  !\right.
  !\end{equation}
  ! 
  ! The algorithm for calculating the derivatives is:
  !
  !\begin{subequations}\label{ratfun-5}
  !\begin{align}
  !\frac{\partial f(x,\{c\})}{\partial c_1}=&
  !\frac{1}{Q(x,\{c\})}\\
  !\frac{\partial f(x,\{c\})}{\partial c_{n+2}}=&
  !-\frac{f(x,\{c\})}{Q(x,\{c\})}\\
  !\frac{\partial f(x,\{c\})}{\partial c_j}=&\left\{
  !\begin{array}{lc}
  !{\displaystyle \frac{\partial f(x,\{c\})}{\partial c_{j-1}}x}& 2 \leq j \leq n+1 \\
  !{\displaystyle -f(x,\{c\})\frac{\partial f(x,\{c\})}{\partial%
  !c_{j-n-1}}}& n+3 \leq j \leq 2n+2 \\
  !\end{array}
  !\right.
  !\end{align}
  !\end{subequations}
  !
  !
  ! The hessian is returned in the array \texttt{hess}, and is calculated
  !as:
  !\begin{equation}\label{ratfun-6}     
  !\frac{\partial^2 f(x,\{c\})}{\partial c_k\partial c_j}=\left\{
  !\begin{array}{lc}
  ! 0 & 1 \leq j,k \leq n+1 \\
  !{\displaystyle -\frac{x^{j-1+k-n-1}}{[Q(x,\{c\})]^2}}& %
  !1 \leq j \leq n+1 < k \leq 2n+2 \vspace{3mm}\\
  !{\displaystyle 2\frac{P(x,\{c\})}{[Q(x,\{c\})]^3}x^{j+k-2n-2}}& n+2 %
  !\leq j,k \leq 2n+2 \\
  !\end{array}
  !\right.
  !\end{equation}
  !
  ! which are calculated as:
  !
  !\begin{equation}\label{ratfun-7}     
  !\frac{\partial^2 f(x,\{c\})}{\partial c_k\partial c_j}=\left\{
  !\begin{array}{lc}
  ! 0 & 1 \leq j,k \leq n+1 \vspace{3mm}\\
  !{\displaystyle -\frac{\partial f(x,\{c\})}{\partial c_{j}}\frac{\partial%
  !f(x,\{c\})}{\partial c_{k-n-1}}}&1 \leq j \leq n+1 < k \leq 2n+2 
  !\vspace{3mm}\\
  !{\displaystyle 2f(x,\{c\})\frac{\partial f(x,\{c\})}{\partial%
  !c_{j-n-1}}\frac{\partial f(x,\{c\})}{\partial%
  !c_{k-n-1}}}& n+2 %
  !\leq j,k \leq 2n+2 \\
  !\end{array}
  !\right.
  !\end{equation}
  !
  ! The derivatives with respect to the real parameters are easily
  !calculated by using
  !
  !\begin{equation}\label{ratfun-8}     
  !\frac{\partial c_k}{\partial a_j}=\delta_{k,j}+i\delta_{k,j-2n-2}
  !\end{equation}
  ! and the chain rule.
  !\begin{equation}\label{ratfun-9}     
  !\frac{\partial f(x,\{c\})}{\partial a_j}=\sum_{k=1}^{2*n+2}%
  !{\frac{\partial f(x,\{c\})}{\partial c_k}\frac{\partial c_k}{\partial
  !a_j}}
  !\end{equation}
  !
  implicit none
  integer,    intent(in)  :: n
  complex*16, intent(in)  :: x
  real(8),    intent(in)  :: a(4*n+4)
  complex*16, intent(out) :: y
  complex*16, intent(out) :: dyda(4*n+4), hess(4*n+4,4*n+4)
  !
  integer :: i, j, l, m
  complex*16 :: q, p, ca(2*n+2)
  complex*16, parameter :: czero = (0.0d0,0.0d0)      
  complex*16, parameter :: cone  = (1.0d0,0.0d0)      
  complex*16, parameter :: cim   = (0.0d0,1.0d0)      
  !     Initializations
  dyda(:)=czero
  !     Convert the real coefficient to complex ones.      
  do i=1,2*n+2
     ca(i)=cmplx(a(i),a(2*n+2+i),8)
  enddo
  !     calculate the qominator and perator of y
  q=czero
  p=czero
  do i=1,n
     q = (q + ca(2*n+3-i)) * x
     p = (p + ca(n+2-i)) * x
  enddo
  q = (q + ca(n+2))* x + cone
  p = p + ca(1)
  y = p / q
  !     Calculate the partial derivatives with respect to the complex
  !     coefficients ca(:) (equal to the derivative with respect to the real
  !     part of the coefficient
  dyda(1) = 1 / q
  do i=1, n
     dyda(i+1) = dyda(i) * x
     dyda(n+1+i) = - y * dyda(i+1)
  enddo
  dyda(2*n+2) = - y * dyda(n+1) * x
  !     Convert to the derivatives with respect to the real coeficients a(:)
  !     For the real part of the coefficient it is done, just extend to the
  !     imaginary part.  
  do i=1,2*n+2
     dyda(2*n+2+i) = dyda(i) * cim
  enddo
  !     calculate the Hessian
  do i=1, n+1
     l= i + n + 1
     do j=1, n+1
        m = j + n + 1
        hess(i,j) = czero
        hess(l,j) = -dyda(i) * dyda(j)
        hess(i,m) = hess(l,j)
        hess(l,m) = -y * dyda(i) * dyda(j)
     enddo
  enddo
  !     extend to the imaginary part of the coefficients
  do i=1,2*n+2
     do j=1,2*n+2
        hess(2*n+2+i,j) = hess (i,j) * cim      
        hess(i,2*n+2+j) = hess (i,j) * cim      
        hess(2*n+2+i,2*n+2+j)= - hess(i,j)
     enddo
  enddo
  return
end  subroutine ratfun


subroutine pd_funval(y,x,a,n,m)
  ! Given the abscissa \texttt{x}, and the parameters \texttt{a(1:4n+4)}
  !this subroutine returns in \texttt{y} the value of the function:
  !
  !\begin{equation}\label{ratfun-1}     
  !f(x,\{c\})=\frac{P(x,\{c\})}{Q(x,\{c\})}
  !\end{equation}
  !
  !where 
  !\begin{equation}\label{ratfun-2}     
  !P(x,\{c\})=\sum\limits_{k=1}^{n+1}{c_{k}x^{k-1}},
  !\end{equation}
  !
  !\begin{equation}\label{ratfun-3}     
  !Q(x,\{c\})=1+\sum\limits_{k=n+2}^{2n+2}{c_k x^{k-n-1}}
  !\end{equation}
  !
  !and $c_k=\texttt{a(k)}+i\texttt{a(k+2n+2)}$.
  implicit none
  integer, intent(in)  :: n, m
  real(8), intent(in)  :: x(2*m)
  real(8), intent(in)  :: a(4*n+4)
  real(8), intent(out) :: y(2*m)
  !
  integer :: i, j
  complex*16 :: q, p, yc, xc, ca(2*n+2)
  complex*16, parameter :: czero = (0.0d0,0.0d0)      
  complex*16, parameter :: cone  = (1.0d0,0.0d0)      
  complex*16, parameter :: cim   = (0.0d0,1.0d0)      
  !     Convert the real coefficient to complex ones.
  do i=1,2*n+2
     ca(i)=cmplx(a(i),a(2*n+2+i),8)
  enddo
  do j=1,m
     xc = x(j) + x(j+m)*cim
     q=czero
     p=czero
     do i=1,n
        q = (q + ca(2*n+3-i)) * xc
        p = (p + ca(n+2-i)) * xc
     enddo
     q = (q + ca(n+2))* xc + cone
     p = p + ca(1)
     yc = p / q
     y(j)   = dble(yc)
     y(j+m) = aimag(yc)
  enddo

  return
end subroutine pd_funval


subroutine pd_funvalc(yc,dyc, xc,ca, n, m)
  ! Given the abscissa \texttt{x}, and the parameters \texttt{a(1:2n+2)}
  !this subroutine returns in \texttt{y} the value of the function:
  !
  !\begin{equation}\label{acrgn-1}     
  !f(x,\{c\})=\frac{P(x,\{c\})}{Q(x,\{c\})}
  !\end{equation}
  !
  !where 

  !\begin{equation}\label{acrgn-2}     
  !P(x,\{c\})=\sum\limits_{k=1}^{n+1}{c_{k}x^{k-1}},
  !\end{equation}
  !
  !\begin{equation}\label{acrgn-3}     
  !Q(x,\{c\})=1+\sum\limits_{k=n+2}^{2n+2}{c_k x^{k-n-1}}
  !\end{equation}
  !
  !and $c_k=\texttt{a(k)}+i\texttt{a(k+2n+2)}$. In \texttt{dy} the 
  !derivative of the function with respect to \texttt{x}:
  !
  !\begin{equation}\label{acrgn-4}     
  !\frac{\partial f(x,\{c\})}{\partial x}=\frac{1}{Q(x,\{c\})}\frac{\partial %
  !P(x,\{c\})}{\partial x}-\frac{P(x,\{c\})}{[Q(x,\{c\})]^2}\frac{\partial %
  !Q(x,\{c\})}{\partial x} %
  !\end{equation}
  ! 
  ! with:
  !
  !\begin{subequations}\label{acrgn-5}
  !\begin{align}
  !\frac{\partial P(x,\{c\})}{\partial x}=&\sum\limits_{k=2}^{n+1}{%
  !(k-1)c_{k}x^{k-2}}\\
  !\frac{\partial Q(x,\{c\})}{\partial x}=&
  !\sum\limits_{k=n+2}^{2n+2}{(k-n-1)c_k x^{k-n-2}}\\
  !\end{align}
  !\end{subequations}
  !
  implicit none
  integer,    intent(in)  :: n, m      ! number of fitting parameters and number of frequencies to evaluate on
  complex*16, intent(in)  :: xc(m)     ! x==frequencies
  complex*16, intent(in)  :: ca(2*n+2) ! ca==fitting coeffients
  complex*16, intent(out) :: yc(m)     ! fitted selfec at real frequency 
  complex*16, intent(out) :: dyc(m)    ! first order derivative of selfec with respect to frequency
  !
  integer(4) :: i, j
  complex*16 :: q, p, dq, dp
  complex(8), parameter :: czero = (0.0d0,0.0d0)
  complex(8), parameter :: cone  = (1.0d0,0.0d0)
  complex(8), parameter :: cim   = (0.0d0,1.0d0)
  do j=1,m
     q=czero
     p=czero
     dq=czero
     dp=czero
     do i=n,1,-1
        q = (q + ca(n+2+i)) * xc(j)
        p = (p + ca(i+1))   * xc(j)
     enddo
     q = (q + ca(n+2)) * xc(j) + cone
     p = p + ca(1)
     do i=n,2,-1
        dq = (dq + ca(n+2+i)*dble(i+1))*xc(j)
        dp = (dp + ca(i+1)*dble(i)) * xc(j)
     enddo
     dq = dq + ca(n+2) + 2.0d0*ca(n+3)*xc(j)
     dp = dp + ca(2)
     yc(j) = p / q
     dyc(j) = dp / q - yc(j) * dq / q
  enddo
end subroutine pd_funvalc



subroutine pd_jacoby(dyda, x,a,n,m)
  ! Given the abscissa \texttt{x}, and the parameters \texttt{a(1:4n+4)}
  !this subroutine returns in \texttt{y} the value of the function:
  !
  !\begin{equation}\label{ratfun-1}     
  !f(x,\{c\})=\frac{P(x,\{c\})}{Q(x,\{c\})}
  !\end{equation}
  !
  !where 
  !\begin{equation}\label{ratfun-2}     
  !P(x,\{c\})=\sum\limits_{k=1}^{n+1}{c_{k}x^{k-1}},
  !\end{equation}
  !
  !\begin{equation}\label{ratfun-3}     
  !Q(x,\{c\})=1+\sum\limits_{k=n+2}^{2n+2}{c_k x^{k-n-1}}
  !\end{equation}
  !
  !and $c_k=\texttt{a(k)}+i\texttt{a(k+2n+2)}$. In \texttt{dyda} the partial
  !derivatives of the function with respect to the parameters:
  !
  !\begin{equation}\label{ratfun-4}     
  !\frac{\partial f(x,\{c\})}{\partial c_j}=\left\{
  !\begin{array}{lc}
  !{\displaystyle \frac{x^{j-1}}{Q(x,\{c\})}}& 1 \leq j \leq n+1 \vspace{3mm}\\
  !{\displaystyle -\frac{P(x,\{c\})}{[Q(x,\{c\})]^2}x^{j-n-1}}& n+2 %
  !\leq j \leq 2n+2 \\
  !\end{array}
  !\right.
  !\end{equation}
  ! 
  ! The algorithm for calculating the derivatives is:
  !
  !\begin{subequations}\label{ratfun-5}
  !\begin{align}
  !\frac{\partial f(x,\{c\})}{\partial c_1}=&
  !\frac{1}{Q(x,\{c\})}\\
  !\frac{\partial f(x,\{c\})}{\partial c_{n+2}}=&
  !-\frac{f(x,\{c\})}{Q(x,\{c\})}\\
  !\frac{\partial f(x,\{c\})}{\partial c_j}=&\left\{
  !\begin{array}{lc}
  !{\displaystyle \frac{\partial f(x,\{c\})}{\partial c_{j-1}}x}& 2 \leq j \leq n+1 \\
  !{\displaystyle -f(x,\{c\})\frac{\partial f(x,\{c\})}{\partial%
  !c_{j-n-1}}}& n+3 \leq j \leq 2n+2 \\
  !\end{array}
  !\right.
  !\end{align}
  !\end{subequations}
  !
  implicit none
  integer, intent(in)  :: n,m
  real(8), intent(in)  :: x(2*m)
  real(8), intent(in)  :: a(4*n+4)
  real(8), intent(out) :: dyda(2*m,4*n+4)
  !
  integer :: i, j
  complex*16 :: y(m)
  complex*16 :: q, p, xc, ca(2*n+2)
  complex*16 :: cdyda(4*n+4)
  complex*16, parameter :: czero = (0.0d0,0.0d0)      
  complex*16, parameter :: cone  = (1.0d0,0.0d0)      
  complex*16, parameter :: cim   = (0.0d0,1.0d0)      
  !     Convert the real coefficient to complex ones.      
  do i=1,2*n+2
     ca(i)=cmplx(a(i),a(2*n+2+i),8)
  enddo
  !     Initializations
  !     calculate the qominator and perator of y
  do j=1,m
     xc = x(j) + x(j+m)*cim
     cdyda(:)=czero
     q=czero
     p=czero
     do i=1,n
        q = (q + ca(2*n+3-i)) * xc
        p = (p + ca(n+2-i)) * xc
     enddo
     q = (q + ca(n+2))* xc + cone
     p = p + ca(1)
     y(j) = p / q
     !     Calculate the partial derivatives with respect to the complex
     !     coefficients ca(:) (equal to the derivative with respect to the real
     !     part of the coefficient
     cdyda(1) = 1 / q
     do i=1, n
        cdyda(i+1) = cdyda(i) * xc
        cdyda(n+1+i) = - y(j) * cdyda(i+1)
     enddo
     cdyda(2*n+2) = - y(j) * cdyda(n+1) * xc
     !     Convert to the derivatives with respect to the real coeficients a(:)
     !     For the real part of the coefficient it is done, just extend to the
     !     imaginary part.  
     do i=1,2*n+2
        cdyda(2*n+2+i) = cdyda(i) * cim
     enddo
     dyda(j,:)   = dble(cdyda(:))
     dyda(j+m,:) = aimag(cdyda(:))
  enddo
  return
end subroutine pd_jacoby
