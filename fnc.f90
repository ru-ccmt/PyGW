subroutine fermi(y, x, n)
  real*8, intent(out) :: y(n)
  real*8, intent(in)  :: x(n)
  integer,intent(in)  :: n
  do i=1,n
     y(i) = 1.d0/(exp(x(i))+1.d0)
  end do
end subroutine fermi

subroutine gauss(y, x, s, n)
  real*8, intent(out) :: y(n)
  real*8, intent(in)  :: x(n)
  real*8, intent(in)  :: s
  integer,intent(in)  :: n
  !
  real*8, parameter :: sr2pi = 2.506628274631d0
  do i=1,n
     y(i) = 1.d0/(s*sr2pi)*exp(-x(i)*x(i)/(2.d0*s*s))
  enddo
end subroutine gauss

subroutine cfermi(y, x, n)
  real*8, intent(out) :: y(n)
  real*8, intent(in)  :: x(n)
  integer,intent(in)  :: n
  do i=1,n
     if (x(i)>10.) then
        y(i) = 0.d0
     else if (x(i) < -10.d0) then
        y(i) = 1.d0
     else
        y(i) = 1.d0/(exp(x(i))+1.d0)
     endif
  end do
end subroutine cfermi
subroutine cgauss(y, x, s, n)
  real*8, intent(out) :: y(n)
  real*8, intent(in)  :: x(n)
  real*8, intent(in)  :: s
  integer,intent(in)  :: n
  !
  real*8, parameter :: sr2pi = 2.506628274631
  do i=1,n
     if ( abs(x(i))>10.d0*s) then
        y(i) = 0.d0
     else
        y(i) = 1.d0/(s*sr2pi)*exp(-x(i)*x(i)/(2*s*s))
     endif
  end do
end subroutine cgauss


subroutine fr_convolution(res, iom, enk, mwm, omega, womeg, nom)
  !  We are computing  
  !       sigma(iom) = -1/beta \sum_{iOm} mwm[iOm]/(iom-eps+iOm)
  !   Because mwm[-iOm]=mwm[iOm] and at T=0, we can rewrite
  !       sigma(iom) = 1/pi Integrate[ mwm[Om]*(eps-iom)/((eps-iom)^2 + Om^2) , {Om,0,Infinity} ]
  !   Finally, we notice that
  !                         Integrate[ 1/((eps-iom)^2 + Om^2) , {Om,0,Infinity} ] = pi*sign(eps)/(2*(eps-iom))
  !   therefore
  !       sigma(iom) = (eps-iom)/pi Integrate[ (mwm[Om]-mwm[iom])/((eps-iom)^2 + Om^2) , {Om,0,Infinity} ] + mwm[iom] * sign(eps)/2
  !
  complex*16, intent(out) :: res
  integer,    intent(in)  :: iom, nom
  real*8,     intent(in)  :: enk
  real*8,     intent(in)  :: omega(nom), womeg(nom)
  complex*16, intent(in)  :: mwm(nom)
  !
  real(8),    parameter  :: pi = 3.14159265358979d0
  complex*16 :: eps_om, eps_om2, sc0
  integer    :: i
  
  eps_om = cmplx(enk, -omega(iom), 8)
  eps_om2 = eps_om * eps_om
  sc0 = 0.d0
  do i=1,nom
     if (i==iom) cycle
     sc0 = sc0 + (mwm(i)-mwm(iom)) * womeg(i) / (omega(i)**2 + eps_om2)
  enddo
  res = eps_om*sc0/pi
  if (abs(mwm(iom))>1e-16) then
     res = res + mwm(iom)/pi * atan(omega(nom)/eps_om) !0.5d0*pi*sign(1.d0, enk)
  endif
end subroutine fr_convolution

subroutine fr_convolution2(res, enk, mwm, mwm_iom, omega, omega_iom, womeg, nom)
  !  We are computing  
  !       sigma(iom) = -1/beta \sum_{iOm} mwm[iOm]/(iom-eps+iOm)
  !   Because mwm[-iOm]=mwm[iOm] and at T=0, we can rewrite
  !       sigma(iom) = 1/pi Integrate[ mwm[Om]*(eps-iom)/((eps-iom)^2 + Om^2) , {Om,0,Infinity} ]
  !   Finally, we notice that
  !                         Integrate[ 1/((eps-iom)^2 + Om^2) , {Om,0,Infinity} ] = pi*sign(eps)/(2*(eps-iom))
  !   therefore
  !       sigma(iom) = (eps-iom)/pi Integrate[ (mwm[Om]-mwm[iom])/((eps-iom)^2 + Om^2) , {Om,0,Infinity} ] + mwm[iom] * sign(eps)/2
  !
  complex*16, intent(out) :: res
  complex*16, intent(in)  :: mwm_iom
  integer,    intent(in)  :: nom
  real*8,     intent(in)  :: enk
  real*8,     intent(in)  :: omega(nom), womeg(nom)
  real*8,     intent(in)  :: omega_iom
  complex*16, intent(in)  :: mwm(nom)
  !
  real(8),    parameter  :: pi = 3.14159265358979d0
  complex*16 :: eps_om, eps_om2, sc0
  integer    :: i
  
  eps_om = cmplx(enk, -omega_iom, 8)
  eps_om2 = eps_om * eps_om
  sc0 = 0.d0
  do i=1,nom
     !if (i==iom) cycle
     sc0 = sc0 + (mwm(i)-mwm_iom) * womeg(i) / (omega(i)**2 + eps_om2)
  enddo
  res = eps_om*sc0/pi
  if (abs(mwm_iom)>1e-16) then
     res = res + mwm_iom/pi * atan(omega(nom)/eps_om) !0.5d0*pi*sign(1.d0, enk)
  endif
end subroutine fr_convolution2

subroutine many_fr_convolutions(sC, enk, Ul, omega, womeg, nom, nl, nb2)
  implicit none
  integer,    intent(in) :: nom, nl, nb2
  complex*16, intent(out):: sC(nl,nb2,nom)
  complex*16, intent(in) :: Ul(nl,nom) ! nl=22,nom=32
  real*8,     intent(in) :: enk(nb2) ! nb2=44
  real*8,     intent(in) :: omega(nom), womeg(nom)
  !
  integer    :: il, ie2, iom
  do il=1,nl
     do iom=1,nom
        do ie2=1,nb2
           call fr_convolution2(sC(il,ie2,iom), enk(ie2), Ul(il,:), Ul(il,iom), omega, omega(iom), womeg, nom)
        enddo
     enddo
  enddo
end subroutine many_fr_convolutions

subroutine all_fr_convolutions(sc_p, enk, mwm, omega, womeg, nb1,nb2,nom)
  implicit none
  integer,    intent(in) :: nom, nb1, nb2
  complex*16, intent(out):: sc_p(nb1,nom)
  complex*16, intent(in) :: mwm(nom,nb1,nb2)
  real*8,     intent(in) :: enk(nb2)
  real*8,     intent(in) :: omega(nom), womeg(nom)
  !
  complex*16 :: sc, sc_sum
  integer    :: ie1, ie2, iom
  do ie1=1,nb1
     do iom=1,nom
        sc_sum = 0.d0
        do ie2=1,nb2
           call fr_convolution2(sc, enk(ie2), mwm(:,ie1,ie2), mwm(iom,ie1,ie2), omega, omega(iom), womeg, nom)
           sc_sum = sc_sum + sc
        enddo
        sc_p(ie1,iom) = sc_sum
     enddo
  enddo
end subroutine all_fr_convolutions

subroutine few_fr_convolutions(res, enk, mwm, mwm_iom, omega, omega_iom, womeg, nom, nb2)
  implicit none
  integer,    intent(in) :: nom, nb2
  complex*16, intent(out):: res
  complex*16, intent(in) :: mwm(nb2,nom), mwm_iom(nb2)
  real*8,     intent(in) :: enk(nb2)
  real*8,     intent(in) :: omega(nom), womeg(nom), omega_iom
  !
  integer    :: ie2
  complex*16 :: sc
  res = 0.d0
  do ie2=1,nb2
     call fr_convolution2(sc, enk(ie2), mwm(ie2,:), mwm_iom(ie2), omega, omega_iom, womeg, nom)
     res = res + sc
  enddo
end subroutine few_fr_convolutions
    
