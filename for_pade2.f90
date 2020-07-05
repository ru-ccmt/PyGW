subroutine mrqmin(x,y,ndata,a,ma,covar,alpha,nca,chisq,alamda)
  !Levenberg-Marquardt method, attempting to reduce the value
  !$\chi^2$ of a fit between a set of data points
  !\texttt{x(1:ndata)}, \texttt{y(1:ndata)}, and a nonlinear function
  !dependent on \texttt{ma} coefficients \texttt{a(1:ma)}. The program 
  !returns current best-fit values for the
  !parameters \texttt{a(1:ma)}, and $\chi^2 = \texttt{chisq}$. The
  !arrays \texttt{covar(1:nca,1:nca)}, \texttt{alpha(1:nca,1:nca)}
  !with physical dimension \texttt{nca} ($\geq$ the number of fitted
  !parameters) are used as working space during most iterations.
  ! On the first call provide an initial guess for the
  !parameters \texttt{a}, and set $\texttt{alamda}<0$ for
  !initialization (which then sets $\texttt{alamda}=.001$). If
  !a step succeeds \texttt{chisq} becomes smaller and \texttt{alamda}
  !decreases by a factor of 10. If a step fails \texttt{alamda} grows
  !by a factor of 10. You must call this routine repeatedly until
  !convergence is achieved. Then, make one final call with
  !$\texttt{alamda}=0$, so that \texttt{covar(1:ma,1:ma)} returns the
  !covariance matrix, and \texttt{alpha} the curvature matrix.
  !(Parameters held fixed will return zero covariances.)
  implicit none
  integer(4), intent(in) :: ndata ! Number of datapoints to be fitted
  integer(4), intent(in) :: ma    ! Number of coefficients
  integer(4), intent(in) :: nca   ! Size of working-space arrays
  real(8),    intent(in) :: x(*)  ! abscisas of the datapoints
  complex*16, intent(in) :: y(*)  ! data value of the datapoint
  real(8),    intent(inout) :: alamda, chisq
  real(8),    intent(inout) :: a(ma)
  real(8),    intent(inout) :: covar(nca,nca), alpha(nca,nca)
  !
  integer :: j, k, l
  real(8) :: ochisq
  integer, parameter :: max = 50
  real(8) :: atry(max), beta(max), da(max)
  !
  save ochisq,atry,beta,da
  external covsrt
  external gaussj
  external mrqcof
  ! Original subroutine: mrqmin.for (c) copr. 1986-92 numerical recipes
  ! software &680i..
  ! Last modified: 7th. Jul. 2005 by RGA
  if(alamda.lt.0.0d0)then ! Initialization.
     alamda=1.0d-3
     call mrqcof(x,y,ndata,a,ma,alpha,beta,nca,chisq)
     ochisq=chisq
     do j=1,ma
        atry(j)=a(j)
     enddo
  endif
  !     Alter linearized fitting matrix, by augmenting diagonal elements.
  do j=1,ma
     do k=1,ma
        covar(j,k)=alpha(j,k)
     enddo ! k
     covar(j,j)=alpha(j,j)*(1.0d0+alamda)
     da(j)=beta(j)
  enddo ! j
  call gaussj(covar,ma,nca,da,1,1)
  !      Once converged, evaluate covariance matrix. 
  !      
  !      if(alamda.eq.0.0d0)then
  !        call covsrt(covar,nca,ma,ia,mfit)
  !        
  !       Spread out alpha to its full size too.
  !        
  !        call covsrt(alpha,nca,ma,ia,mfit) 
  !        return
  !      endif
  j=0
  !      
  !     Did the trial succeed?
  ! 
  do l=1,ma
     atry(l)=a(l)+da(l)
  enddo ! l

  call mrqcof(x,y,ndata,atry,ma,covar,da,nca,chisq)

  if(chisq.lt.ochisq) then     ! Success, accept the new solution.
     alamda=0.1d0*alamda
     ochisq=chisq
     do j=1,ma
        do k=1,ma
           alpha(j,k)=covar(j,k)
        enddo ! k
        beta(j)=da(j)
     enddo ! j
     do l=1,ma
        a(l)=atry(l)
     enddo ! l
  else                        ! Failure, increase alamda and return.
     alamda=1.0d+1*alamda
     chisq=ochisq
  endif
  return
end subroutine mrqmin
      
subroutine stdesc(x,y,ndata,a,ma,varsq)
  ! This program perform one iteration of the steepest descent method.
  implicit none
  integer,    intent(in) :: ndata ! Number of datapoints to be fitted
  integer,    intent(in) :: ma    ! Number of coefficients
  real(8),    intent(in) :: x(*)  ! abscisas of the datapoints
  complex*16, intent(in) :: y(*)  ! data value of the datapoint
  real(8), intent(inout) :: varsq
  real(8), intent(inout) :: a(ma)
  !
  integer(4) :: j, i, it
  real(8) :: varsqold, num, denomi, denom, lambda, deltach
  real(8) :: atemp(ma)
  real(8), allocatable :: beta(:), alpha(:,:)
  logical :: conv
  real(8), parameter :: tolstd = 1.0d-3
  external mrqcof
  !
  allocate(beta(1:ma))
  allocate(alpha(1:ma,1:ma))
  varsqold=0.0d0
  conv=.false.
  it=0
  do while (.not.conv)
     it=it+1
     call mrqcof(x,y,ndata,a,ma,alpha,beta,ma,varsq)
     num=0.0d0
     denom=0.0d0
     do i=1,ma
        num=num+beta(i)*beta(i)
        denomi=0.0d0
        do j=1,ma
           denomi=denomi+alpha(i,j)*beta(j)
        enddo
        denom=denom+denomi*beta(i)
     enddo
     lambda=num/denom
     do j=1,ma 
        a(j)=a(j)+lambda*beta(j)
     enddo
     deltach=abs(varsq-varsqold)
     varsqold=varsq
     atemp(1:ma)=a(1:ma)
     conv=((deltach.lt.tolstd).or.(it.gt.200))
  enddo
  if(it.gt.200)write(*,*)'WARNING: stepest descent did not converge after',it, 'iterations'
  deallocate(beta)
  deallocate(alpha)
  return
end subroutine stdesc

subroutine mrqcof(x,y,ndata,a,ma,alpha,beta,nca,chisq)
  !  Used by \texttt{mrqmin} to evaluate the linearized 
  !  fitting matrix \texttt{alpha}, and vector \texttt{beta} as:
  !
  ! \begin{equation}\label{mrqcof-01}
  ! \begin{aligned}
  !\beta_k &\equiv - \frac{1}{2}\frac{\partial\chi^2}{\partial a_k} &
  !\alpha_{kl}&\equiv \frac{1}{2}\frac{\partial^2\chi^2}{\partial
  !a_k\partial a_l} , 
  !\end{aligned}
  !\end{equation}
  !and calculate $\chi^2$. 
  implicit none
  integer(4), intent(in) :: ndata ! Number of datapoints to be fitted
  integer(4), intent(in) :: ma    ! Number of coefficients
  integer(4), intent(in) :: nca   ! Size of working-space arrays
  real(8),    intent(in) :: x(*)  ! abscisas of the datapoints
  complex*16, intent(in) :: y(*)  ! data value of the datapoint
  real(8),    intent(inout) :: chisq
  real(8),    intent(inout) :: a(ma), beta(ma), alpha(nca,nca)
  integer(4) :: i, j, k, l, m
  complex*16 :: wt, dy, ymod, cx
  complex*16 :: dyda(ma), d2yda(ma,ma) 
  external ratfun
  !
  ! Original subroutine: mrqcof.for (c) copr. 1986-92 numerical recipes
  ! software &681i..
  ! Last modified: 7th. Jul. 2005 by RGA
  !
  do j=1,ma
     !      
     !       Initialize (symmetric) alpha, beta. 
     !
     do k=1,j
        alpha(j,k)=0. 
     enddo ! k
     beta(j)=0. 
  enddo ! j
  chisq=0.0d0 
  !      
  !     Summation loop over all data. 
  !
  do i=1,ndata 
     cx=cmplx(0.0d0,x(i),8)
     call ratfun(cx,a,ymod,dyda,d2yda,ma/4-1)
     dy=y(i)-ymod 
     do l=1,ma 
        wt=dyda(l) 
        do m=1,l 
           alpha(l,m)=alpha(l,m)+real(wt*conjg(dyda(m)))
        enddo
        beta(l)=beta(l)+real(dy*conjg(wt))
     enddo ! l
     chisq = chisq + dble(dy*conjg(dy)) ! And find \chi^2. 
  enddo ! i
  chisq=chisq/ndata
  !     Fill in the symmetric side. 
  do j=2,ma
     do k=1,j-1 
        alpha(k,j)=alpha(j,k) 
     enddo ! k
  enddo ! j 
  return 
end subroutine mrqcof

      
subroutine gaussj(a,n,np,b,m,mp)
  ! Linear equation solution by Gauss-Jordan elimination, equation
  !\ref{gaussj-1} below.
  !
  !\begin{equation}\label{gaussj-1}
  !\begin{aligned}
  !  \begin{bmatrix}
  !    a_{11}&a_{12}&a_{13}&a_{14}\\
  !    a_{21}&a_{22}&a_{23}&a_{24}\\
  !    a_{31}&a_{32}&a_{33}&a_{34}\\
  !    a_{41}&a_{42}&a_{43}&a_{44}
  !  \end{bmatrix}\cdot &
  !  \begin{bmatrix}
  !    \begin{pmatrix}
  !      x_{11}\\
  !      x_{21}\\
  !      x_{31}\\
  !      x_{41}
  !    \end{pmatrix}\sqcup
  !    \begin{pmatrix}
  !      x_{12}\\
  !      x_{22}\\
  !      x_{32}\\
  !      x_{42}
  !    \end{pmatrix}\sqcup
  !    \begin{pmatrix}
  !      x_{13}\\
  !      x_{23}\\
  !      x_{33}\\
  !      x_{43}
  !    \end{pmatrix}\sqcup
  !    \begin{pmatrix}
  !      y_{11}&y_{12}&y_{13}&y_{14}\\
  !      y_{21}&y_{22}&y_{23}&y_{24}\\
  !      y_{31}&y_{32}&y_{33}&y_{34}\\
  !      y_{41}&y_{42}&y_{43}&y_{44}
  !    \end{pmatrix}
  !  \end{bmatrix}\\
  != & \begin{bmatrix}
  !    \begin{pmatrix}
  !      b_{11}\\
  !      b_{21}\\
  !      b_{31}\\
  !      b_{41}
  !    \end{pmatrix}\sqcup
  !    \begin{pmatrix}
  !      b_{12}\\
  !      b_{22}\\
  !      b_{32}\\
  !      b_{42}
  !    \end{pmatrix}\sqcup
  !    \begin{pmatrix}
  !      b_{13}\\
  !      b_{23}\\
  !      b_{33}\\
  !      b_{43}
  !    \end{pmatrix}\sqcup
  !    \begin{pmatrix}
  !      1 & 0 & 0 & 0 \\
  !      0 & 1 & 0 & 0 \\
  !      0 & 0 & 1 & 0 \\
  !      0 & 0 & 0 & 1 \\
  !    \end{pmatrix}
  !  \end{bmatrix}
  !\end{aligned}
  !\end{equation}
  ! \texttt{a(1:n,1:n)}
  ! is an imput matrix stored in an array of physical dimensions \texttt{np}
  ! by \texttt{np}.\texttt{b(1:n,1:m)} is an input matrix containing the
  ! \texttt{m} right-hand side vectors, stored in an array of physical
  ! dimensions \texttt{np} by \texttt{mp}. On output, \texttt{a(1:n,1:n)} is
  ! replaced by its matrix inverse, and \texttt{b(1:n,1:m)} is replaced by
  ! the corresponding set of solution vectors.
  implicit none
  integer, intent(in)    :: m, mp, n, np
  real(8), intent(inout) :: a(np,np), b(np,mp)
  integer(4) :: i, j, k, l, ll, icol, irow
  integer(4) :: indxc(1:np) ! Used for bookkeeping on the pivoting
  integer(4) :: indxr(1:np) ! Used for bookkeeping on the pivoting
  integer(4) :: ipiv(1:np)  ! Used for bookkeeping on the pivoting
  real(8) :: big, dum, pivinv
  ! Original subroutine: gaussj.for (c) copr. 1986-92 numerical recipes
  ! software &30i..
  ! Last modified: 7th. Jul. 2005 by RGA
  do j=1,n
     ipiv(j)=0
  enddo ! j
  !     This is the main loop over the columns to be reduced
  do i=1,n
     big=0.0d0
     !        This is the outer loop of the search for a pivot element
     do j=1,n
        if(ipiv(j).ne.1)then
           do  k=1,n
              if (ipiv(k).eq.0) then
                 if (abs(a(j,k)).ge.big)then
                    big=abs(a(j,k))
                    irow=j
                    icol=k
                 endif
              else if (ipiv(k).gt.1) then
                 stop 'singular matrix in gaussj'
              endif
           enddo ! k
        endif
     enddo ! j
     ipiv(icol)=ipiv(icol)+1
     !        
     !       We now have the pivot element, so we interchange rows, if needed,
     !       to put the pivot element on the diagonal. The columns are not
     !       physically interchanged, only relabeled: indxc(i), the column of
     !       the ith pivot element, is the ith column that is reduced, while
     !       indxr(i) is the row in which that pivot   element was originally
     !       located. If indxr(i) \neq indxc(i) there is an implied column
     !       interchange. With this form of bookkeeping, the solution b's will
     !       end up in the correct order, and the inverse matrix will be
     !       scrambled by columns.
     !                 
     if (irow.ne.icol) then
        do l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        enddo ! l
        do l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        enddo ! l
     endif
     indxr(i)=irow
     indxc(i)=icol
     !       We are now ready to divide the pivot row by the pivot element,
     !       located at irow and icol.
     if (a(icol,icol).eq.0.0d0) stop 'singular matrix in gaussj'
     pivinv=1./a(icol,icol)
     a(icol,icol)=1.
     do l=1,n
        a(icol,l)=a(icol,l)*pivinv
     enddo ! l
     do l=1,m
        b(icol,l)=b(icol,l)*pivinv
     enddo ! l
     !       Next, we reduce the rows...
     do ll=1,n
        if(ll.ne.icol)then   ! ... except for the pivor one, of course.
           dum=a(ll,icol)
           a(ll,icol)=0.0d0
           do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
           enddo ! l
           do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
           enddo ! l
        endif
     enddo ! ll
  enddo ! i
  !      
  !     This is the end of the main loop over columns of the reduction. It
  !     only remains to unscramble the solution in view of the column
  !     interchanges. We do this by interchanging pairs of columns in the
  !     reverse order that the permutation was bilt up.
  !         
  do l=n,1,-1
     if(indxr(l).ne.indxc(l))then
        do k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
        enddo ! k
     endif
  enddo ! l
  return
end subroutine gaussj
