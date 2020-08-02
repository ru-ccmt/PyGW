
integer Function get_lm(l,m)
  implicit none
  integer, intent(in) :: l, m
  get_lm = l*(l+1) + m + 1
  return
end Function get_lm


subroutine calcmmatvv_MT(mt12, n1,n2, alfa, beta, gama, nLO_at, mult, iul_ul, iul_udl, iudl_ul, iudl_udl, iulol_ul, iulol_udl, iul_ulol, iudl_ulol, iulol_ulol, nat, lmax, nbmax, ntnt, ndf, nbmax_g, nLOmax, ntnt_g, ndf_g, nlo, lomx, lomax)
  ! This subroutine calculates < psi_{n1}| -i*nabla | psi_{n2}> -- the contribution of the MT-Spheres of the bands n1 and n2
  implicit none
  complex*16, intent(out):: mt12(3)
  integer,    intent(in) :: nat, lmax, nbmax, ntnt, ndf, nbmax_g, nLOmax, ntnt_g, ndf_g, lomx, nlo, lomax
  integer,    intent(in) :: n1, n2, mult(nat)
  complex*16, intent(in) :: alfa(nbmax,ntnt,ndf)
  complex*16, intent(in) :: beta(nbmax,ntnt,ndf)
  complex*16, intent(in) :: gama(nbmax_g,nLOmax,ntnt_g,ndf_g)
  integer,    intent(in) :: nLO_at(0:lomax,nat)
  real(8),    intent(in) :: iul_ul(2,0:lmax-1,nat), iul_udl(2,0:lmax-1,nat), iudl_ul(2,0:lmax-1,nat), iudl_udl(2,0:lmax-1,nat)
  real(8),    intent(in) :: iulol_ul(2,nlo,0:lomx,nat), iulol_udl(2,nlo,0:lomx,nat), iul_ulol(2,nlo,0:lomx,nat), iudl_ulol(2,nlo,0:lomx,nat), iulol_ulol(2,nlo,nlo,0:lomx,nat)
  !
  integer    :: iat, ieq, idf
  integer    :: l, m, lm, l1m1, l1mm1, l1m
  integer    :: ilo, jlo
  real(8)    :: denom, flm1, flm3, flm5 
  complex*16 :: xpy1, xpy2, xmy1, xmy2, z1, z2, pxpy, pxmy, pz
  complex*16 :: aa, ab, ba, bb, ag, bg, ga, gb
  complex*16 :: czero, imag
  intrinsic conjg
  intrinsic sqrt
  interface
     integer Function get_lm(l,m)
       integer, intent(in) :: l, m
     end Function get_lm
  end interface
  !--------------------------------------------------------------------!
  !  contribution from MT spheres : < psi_{n1}| i*nabla | psi_{n2}>    !
  !--------------------------------------------------------------------!
  !  DICTIONARY
  ! iul_ul[1,2]     < u_{l+1}  | d/dr-l/r | u_l  >, < u_l  | d/dr + (l+2)/r | u_{l+1}  >
  ! iul_udl[1,2]    < u_{l+1}  | d/dr-l/r |udot_l>, < u_l  | d/dr + (l+2)/r |udot_{l+1}>
  ! iudl_ul[1,2]    <udot_{l+1}| d/dr-l/r | u_l  > ,<udot_l| d/dr + (l+2)/r | u_{l+1}  >
  ! iudl_udl[1,2]   <udot_{l+1}| d/dr-l/r |udot_l>, <udot_l| d/dr + (l+2)/r |udot_{l+1}>
  ! iulol_ul[1,2]   < ulo_{l+1}| d/dr-l/r | u_l  >, < ulo_l| d/dr + (l+2)/r | u_{l+1}  >
  ! iulol_udl[1,2]  < ulo_{l+1}| d/dr-l/r |udot_l>, < ulo_l| d/dr + (l+2)/r |udot_{l+1}>
  ! iul_ulol[1,2]   < u_{l+1}  | d/dr-l/r | ulo_l>, < u_l  | d/dr + (l+2)/r | ulo_{l+1}>
  ! iudl_ulol[1,2]  <udot_{l+1}| d/dr-l/r | ulo_l>, <udot_l| d/dr + (l+2)/r | ulo_{l+1}>
  ! iulol_ulol[1,2] < ulo_{l+1,ilo}| d/dr-l/r | ulo_{l,jlo}>, < ulo_{l,ilo}| d/dr + (l+2)/r | ulo_{l+1,jlo}>
  !
  czero = cmplx(0.d0,0.d0,8)
  imag  = cmplx(0.d0,1.d0,8)
  pxpy = czero
  pxmy = czero
  pz   = czero
  idf = 0
  
  do iat = 1,nat
     do ieq = 1,mult(iat)
        idf = idf+1
        !$OMP PARALLEL DO PRIVATE(lm,l,m,ilo,jlo,denom,l1mm1,l1m,l1m1,flm1,flm3,flm5,aa,ab,ba,bb,ag,bg,ga,gb,xpy1,xpy2,xmy1,xmy2,z1,z2)&
        !$OMP& REDUCTION(+:pxpy,pxmy,pz)&
        !$OMP& SCHEDULE(STATIC)
        do lm=1,lmax*lmax
           l = int(aint(sqrt(lm-1.d0)+1e-6))
           m = lm-1-l*(l+1)
           denom = dble((2*l+1)*(2*l+3))
           l1mm1 = get_lm(l+1,m-1) ! (l+1,m-1) == (l+1)*(l+2) + m  = lm + 2*l+1
           l1m   = get_lm(l+1,m)   ! (l+1, m ) == (l+1)*(l+2) + m + 1 = lm + 2*l+2
           l1m1  = get_lm(l+1,m+1) ! (l+1,m+1) == (l+1)*(l+2) + m + 2
           flm1 = -sqrt(dble((l+m+1)*(l+m+2))/denom)  !  == -a(l,m)
           flm3 =  sqrt(dble((l-m+1)*(l-m+2))/denom)  !  ==  a(l,-m)
           flm5 =  sqrt(dble((l-m+1)*(l+m+1))/denom)  !  ==  f(l,m)
           !  xpy1 = < psi_{n1}(l+1,m+1) | (d/dr-l/r) | psi_{n2}(l,m) >
           aa = conjg(alfa(n1,l1m1,idf))*alfa(n2,lm,idf)  !  A[n1](l+1,m+1)^*  A[n2](l,m)
           ab = conjg(alfa(n1,l1m1,idf))*beta(n2,lm,idf)  !  A[n1](l+1,m+1)^*  B[n2](l,m)
           ba = conjg(beta(n1,l1m1,idf))*alfa(n2,lm,idf)  !  B[n1](l+1,m+1)^*  A[n2](l,m)
           bb = conjg(beta(n1,l1m1,idf))*beta(n2,lm,idf)  !  B[n1](l+1,m+1)^*  B[n2](l,m)
           !  < A[n1](l+1,m+1) u_{l+1} + B[n1](l+1,m+1) udot_{l+1}| (d/dr-l/r) | A[n2](l,m) u_l + B[n2](l,m) udot_l >
           xpy1 = aa * iul_ul(1,l,iat) + ab * iul_udl(1,l,iat) + ba * iudl_ul(1,l,iat) + bb * iudl_udl(1,l,iat)
           !  xpy2 = < psi_{n1}(l,m)  | d/dr + (l+2)/r | psi_{n2}(l+1,m-1) > 
           aa = conjg(alfa(n1,lm,idf))*alfa(n2,l1mm1,idf) !  A[n1](l,m)^* A[n2](l+1,m-1)
           ab = conjg(alfa(n1,lm,idf))*beta(n2,l1mm1,idf) !  A[n1](l,m)^* B[n2](l+1,m-1)
           ba = conjg(beta(n1,lm,idf))*alfa(n2,l1mm1,idf) !  B[n1](l,m)^* A[n2](l+1,m-1)
           bb = conjg(beta(n1,lm,idf))*beta(n2,l1mm1,idf) !  B[n1](l,m)^* B[n2](l+1,m-1)
           !  < A[n1](l,m) u_l + B[n1](l,m) udot_l  | d/dr + (l+2)/r | A[n2](l+1,m-1) u_{l+1} + B[n2](l+1,m-1) udot_{l+1} >
           xpy2 = aa * iul_ul(2,l,iat) + ab * iul_udl(2,l,iat) + ba * iudl_ul(2,l,iat) + bb * iudl_udl(2,l,iat)
           !  xmy1 = < psi_{n1}(l+1,m-1) | d/dr-l/r | psi_{n2}(l,m)  >
           aa = conjg(alfa(n1,l1mm1,idf))*alfa(n2,lm,idf) !  A[n1](l+1,m-1)^* A[n2](l,m)
           ab = conjg(alfa(n1,l1mm1,idf))*beta(n2,lm,idf) !  A[n1](l+1,m-1)^* B[n2](l,m)
           ba = conjg(beta(n1,l1mm1,idf))*alfa(n2,lm,idf) !  B[n1](l+1,m-1)^* A[n2](l,m)
           bb = conjg(beta(n1,l1mm1,idf))*beta(n2,lm,idf) !  B[n1](l+1,m-1)^* B[n2](l,m)
           !  < A[n1](l+1,m-1) u_{l+1} + B[n1](l+1,m-1) udot _{l+1} | d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l >
           xmy1 = aa * iul_ul(1,l,iat) + ab * iul_udl(1,l,iat) + ba * iudl_ul(1,l,iat) + bb * iudl_udl(1,l,iat)
           !  xmy2 = < psi_{n1}(l,m) | d/dr + (l+2)/r | psi_{n2}(l+1,m+1) >
           aa = conjg(alfa(n1,lm,idf))*alfa(n2,l1m1,idf) ! A[n1](l,m)^* A[n2](l+1,m+1)
           ab = conjg(alfa(n1,lm,idf))*beta(n2,l1m1,idf) ! A[n1](l,m)^* B[n2](l+1,m+1)
           ba = conjg(beta(n1,lm,idf))*alfa(n2,l1m1,idf) ! B[n1](l,m)^* A[n2](l+1,m+1)
           bb = conjg(beta(n1,lm,idf))*beta(n2,l1m1,idf) ! B[n1](l,m)^* B[n2](l+1,m+1)
           !  < A[n1](l,m) u_l +B[n1](l,m) udot_l | d/dr + (l+2)/r | A[n2](l+1,m+1) u_{l+1}  + B[n2](l+1,m+1) udot_{l+1}>
           xmy2 = aa * iul_ul(2,l,iat) + ab * iul_udl(2,l,iat) + ba * iudl_ul(2,l,iat) + bb * iudl_udl(2,l,iat)
           !  z1 = < psi_{n1}(l+1,m) | d/dr-l/r | psi_{n2}(l,m) >
           aa = conjg(alfa(n1,l1m,idf))*alfa(n2,lm,idf) !  A[n1](l+1,m)^* A[n2](l,m)
           ab = conjg(alfa(n1,l1m,idf))*beta(n2,lm,idf) !  A[n1](l+1,m)^* B[n2](l,m)
           ba = conjg(beta(n1,l1m,idf))*alfa(n2,lm,idf) !  B[n1](l+1,m)^* A[n2](l,m)
           bb = conjg(beta(n1,l1m,idf))*beta(n2,lm,idf) !  B[n1](l+1,m)^* B[n2](l,m)
           !  < A[n1](l+1,m) u_{l+1} + B[n1](l+1,m) udot_{l+1}| d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l > 
           z1 = aa * iul_ul(1,l,iat) + ab * iul_udl(1,l,iat) + ba * iudl_ul(1,l,iat) + bb * iudl_udl(1,l,iat)
           !  z2 = < psi_{n1}(l,m) | d/dr + (l+2)/r | psi_{n2}(l+1,m) >
           aa = conjg(alfa(n1,lm,idf))*alfa(n2,l1m,idf) !  A[n1](l,m)^* A[n2](l+1,m)
           ab = conjg(alfa(n1,lm,idf))*beta(n2,l1m,idf) !  A[n1](l,m)^* B[n2](l+1,m)
           ba = conjg(beta(n1,lm,idf))*alfa(n2,l1m,idf) !  B[n1](l,m)^* A[n2](l+1,m)
           bb = conjg(beta(n1,lm,idf))*beta(n2,l1m,idf) !  B[n1](l,m)^* B[n2](l+1,m)
           !  < A[n1](l,m) u_l +B[n1](l,m) udot_l | d/dr + (l+2)/r | A[n2](l+1,m) u_{l+1} + B[n2](l+1,m) udot_{l+1} > 
           z2 = aa * iul_ul(2,l,iat) + ab * iul_udl(2,l,iat) + ba * iudl_ul(2,l,iat) + bb * iudl_udl(2,l,iat)
           !  the loop for the Local orbitals on the right side
           do jlo = 1,nLO_at(l,iat) 
              !  xpy1 += < psi_{n1}(l+1,m+1) | d/dr-l/r | C[n2](l,m) ulo_l>
              ag = conjg(alfa(n1,l1m1,idf))*gama(n2,jlo,lm,idf) ! A[n1](l+1,m+1)^* C[n2](l,m)
              bg = conjg(beta(n1,l1m1,idf))*gama(n2,jlo,lm,idf) ! B[n1](l+1,m+1)^* C[n2](l,m)
              !  < A[n1](l+1,m+1) u_{l+1} + B[n1](l+1,m+1) udot_{l+1} | d/dr-l/r | C[n2](l,m) ulo_l>
              xpy1 = xpy1 + ag * iul_ulol(1,jlo,l,iat) + bg * iudl_ulol(1,jlo,l,iat)
              !  xmy1 += < psi_{n1}(l+1,m-1) | d/dr-l/r | C[n2](l,m) ulo_l>
              ag = conjg(alfa(n1,l1mm1,idf))*gama(n2,jlo,lm,idf) !  A[n1](l+1,m-1)^* C[n2](l,m)
              bg = conjg(beta(n1,l1mm1,idf))*gama(n2,jlo,lm,idf) !  B[n1](l+1,m-1)^* C[n2](l,m)
              !  < A[n1](l+1,m-1) u_{l+1} + B[n1](l+1,m-1) udot_{l+1}| d/dr-l/r | C[n2](l,m) ulo_l>
              xmy1 = xmy1 + ag * iul_ulol(1,jlo,l,iat) + bg * iudl_ulol(1,jlo,l,iat)
              !  z1 += < psi_{n1}(l+1,m) | d/dr-l/r | C[n2](l,m) ulo_l> 
              ag = conjg(alfa(n1,l1m,idf))*gama(n2,jlo,lm,idf)  !  A[n1](l+1,m)^* C[n2](l,m)
              bg = conjg(beta(n1,l1m,idf))*gama(n2,jlo,lm,idf)  !  B[n1](l+1,m)^* C[n2](l,m)
              !  < A[n1](l+1,m) u_{l+1} + B[n1](l+1,m) udot_{l+1} | d/dr-l/r | C[n2](l,m) ulo_l> 
              z1 = z1 + ag * iul_ulol(1,jlo,l,iat) + bg * iudl_ulol(1,jlo,l,iat)
           enddo

           do jlo = 1,nLO_at(l+1,iat) 
              !  xpy2 += < psi_{n1}(l,m) | d/dr + (l+2)/r | C[n2](l+1,m-1) ulo_{l+1}>
              ag = conjg(alfa(n1,lm,idf))*gama(n2,jlo,l1mm1,idf)  !  A[n1](l,m)^*  C[n2](l+1,m-1)
              bg = conjg(beta(n1,lm,idf))*gama(n2,jlo,l1mm1,idf)  !  B[n1](l,m)^*  C[n2](l+1,m-1)
              ! < A[n1](l,m) u_l  + B[n1](l,m) udot_l | d/dr + (l+2)/r | C[n2](l+1,m-1) ulo_{l+1}>
              xpy2 = xpy2 + ag * iul_ulol(2,jlo,l,iat) + bg * iudl_ulol(2,jlo,l,iat)
              !  xmy2 += < psi_{n1}(l,m) | d/dr + (l+2)/r | C[n2](l+1,m+1) ulo_{l+1}>
              ag = conjg(alfa(n1,lm,idf))*gama(n2,jlo,l1m1, idf)  !  A[n1](l,m)^*  C[n2](l+1,m+1)
              bg = conjg(beta(n1,lm,idf))*gama(n2,jlo,l1m1, idf)  !  B[n1](l,m)^*  C[n2](l+1,m+1)
              !  < A[n1](l,m) u_l + B[n1](l,m) udot_l | d/dr + (l+2)/r | C[n2](l+1,m+1) ulo_{l+1}>
              xmy2 = xmy2 + ag * iul_ulol(2,jlo,l,iat) + bg * iudl_ulol(2,jlo,l,iat)
              !KH : Introducing a bug to be compatible with pygab
              !xmy2 = xmy1 + ag * iul_ulol(2,jlo,l,iat) + bg * iudl_ulol(2,jlo,l,iat)

              !  z2 += < psi_{n1}(l,m) | d/dr + (l+2)/r | C[n2](l+1,m) ulo_{l+1}>
              ag = conjg(alfa(n1,lm,idf))*gama(n2,jlo,l1m,idf)  ! A[n1](l,m)^*  C[n2](l+1,m)
              bg = conjg(beta(n1,lm,idf))*gama(n2,jlo,l1m,idf)  ! B[n1](l,m)^*  C[n2](l+1,m)
              !  < A[n1](l,m) u_l + B[n1](l,m) udot_l | d/dr + (l+2)/r | C[n2](l+1,m) ulo_{l+1}>
              z2 = z2 + ag * iul_ulol(2,jlo,l,iat) + bg * iudl_ulol(2,jlo,l,iat)
           enddo
           !  the loop for the Local orbitals on the left side 
           do ilo = 1,nLO_at(l+1,iat)
              !  xpy1 += < C[n1](l+1,m+1) ulo_{l+1}| d/dr-l/r | psi_{n2}(l,m) >   
              ga = conjg(gama(n1,ilo,l1m1,idf))*alfa(n2, lm,idf)  ! C[n1](l+1,m+1)^*  A[n2](l,m)
              gb = conjg(gama(n1,ilo,l1m1,idf))*beta(n2, lm,idf)  ! C[n1](l+1,m+1)^*  B[n2](l,m)
              !  < C[n1](l+1,m+1) ulo_{l+1}| d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l >
              xpy1 = xpy1 + ga * iulol_ul(1,ilo,l,iat) + gb * iulol_udl(1,ilo,l,iat)
              !  xmy1 +=  C[n1](l+1,m-1) ulo_{l+1}| d/dr-l/r | psi_{n2}(l,m) >
              ga = conjg(gama(n1,ilo,l1mm1,idf))*alfa(n2,lm,idf)  !  C[n1](l+1,m-1)^*  A[n2](l,m)
              gb = conjg(gama(n1,ilo,l1mm1,idf))*beta(n2,lm,idf)  !  C[n1](l+1,m-1)^*  B[n2](l,m)
              !  < C[n1](l+1,m-1) ulo_{l+1}| d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l >
              xmy1 = xmy1 + ga * iulol_ul(1,ilo,l,iat) + gb * iulol_udl(1,ilo,l,iat)
              !  z1 += < C[n1](l+1,m) ulo_{l+1}| d/dr-l/r | psi_{n2}(l,m) >
              ga = conjg(gama(n1,ilo,l1m,  idf))*alfa(n2,lm,idf)  !  C[n1](l+1,m)^*  A[n2](l,m)
              gb = conjg(gama(n1,ilo,l1m,  idf))*beta(n2,lm,idf)  !  C[n1](l+1,m)^*  B[n2](l,m)
              !  < C[n1](l+1,m) ulo_{l+1}| d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l >
              z1 = z1 + ga * iulol_ul(1,ilo,l,iat) + gb * iulol_udl(1,ilo,l,iat)
           enddo
           do ilo = 1,nLO_at(l,iat)
              !  xpy2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | psi_{n2}(l+1,m-1) >
              ga = conjg(gama(n1,ilo,lm,idf))*alfa(n2,l1mm1,idf) !  C[n1](l,m)^*  A[n2](l+1,m-1)
              gb = conjg(gama(n1,ilo,lm,idf))*beta(n2,l1mm1,idf) !  C[n1](l,m)^*  B[n2](l+1,m-1)
              !  < C[n1](l,m) ulo_l| d/dr + (l+2)/r | A[n2](l+1,m-1) u_{l+1} + B[n2](l+1,m-1) udot_{l+1} >
              xpy2 = xpy2 + ga * iulol_ul(2,ilo,l,iat) + gb * iulol_udl(2,ilo,l,iat)
              !  xmy2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | psi_{n2}(l+1,m+1) >
              ga = conjg(gama(n1,ilo,lm,idf))*alfa(n2,l1m1, idf)  !  C[n1](l,m)^*  A[n2](l+1,m+1)
              gb = conjg(gama(n1,ilo,lm,idf))*beta(n2,l1m1, idf)  !  C[n1](l,m)^*  B[n2](l+1,m+1)
              !  < C[n1](l,m) ulo_l| d/dr + (l+2)/r | A[n2](l+1,m+1) u_{l+1} +  B[n2](l+1,m+1) udot_{l+1} >
              xmy2 = xmy2 + ga * iulol_ul(2,ilo,l,iat) + gb * iulol_udl(2,ilo,l,iat)
              ! z2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | psi_{n2}(l+1,m) >
              ga = conjg(gama(n1,ilo,lm,idf))*alfa(n2,l1m,idf)  !  C[n1](l,m)^*  A[n2](l+1,m)
              gb = conjg(gama(n1,ilo,lm,idf))*beta(n2,l1m,idf)  !  C[n1](l,m)^*  B[n2](l+1,m)
              !  < C[n1](l,m) ulo_l| d/dr + (l+2)/r | A[n2](l+1,m) u_{l+1} + B[n2](l+1,m) udot_{l+1} >
              z2 = z2 + ga * iulol_ul(2,ilo,l,iat) + gb * iulol_udl(2,ilo,l,iat)
           enddo
           do jlo = 1,nLO_at(l,iat)
              do ilo = 1, nLO_at(l+1,iat)
                 ! xpy1 += < C[n1](l+1,m+1) ulo_{l+1}| d/dr-l/r | C[n2](l,m) ulo_l>
                 xpy1 = xpy1 + conjg(gama(n1,ilo,l1m1, idf))*gama(n2,jlo,lm,   idf) * iulol_ulol(1,ilo,jlo,l,iat)
                 ! xmy1 += < C[n1](l+1,m-1) ulo_{l+1}| d/dr-l/r | C[n2](l,m) ulo_l>
                 xmy1 = xmy1 + conjg(gama(n1,ilo,l1mm1,idf))*gama(n2,jlo,lm,   idf) * iulol_ulol(1,ilo,jlo,l,iat)
                 ! z1 += < C[n1](l+1,m) ulo_{l+1}| d/dr-l/r | C[n2](l,m) ulo_l>
                 z1   = z1   + conjg(gama(n1,ilo,l1m,  idf))*gama(n2,jlo,lm,   idf) * iulol_ulol(1,ilo,jlo,l,iat)
                 ! xpy2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | C[n2](l+1,m-1) ulo_{l+1}>
                 xpy2 = xpy2 + conjg(gama(n1,jlo,lm,   idf))*gama(n2,ilo,l1mm1,idf) * iulol_ulol(2,jlo,ilo,l,iat)
                 ! xmy2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | C[n2](l+1,m+1) ulo_{l+1}>
                 xmy2 = xmy2 + conjg(gama(n1,jlo,lm,   idf))*gama(n2,ilo,l1m1, idf) * iulol_ulol(2,jlo,ilo,l,iat)
                 ! z2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | C[n2](l+1,m) ulo_{l+1}>
                 z2   = z2   + conjg(gama(n1,jlo,lm,   idf))*gama(n2,ilo,l1m,  idf) * iulol_ulol(2,jlo,ilo,l,iat)
                 ! iulol_ulol[1,2] < ulo_{l+1,ilo}| d/dr-l/r | ulo_{l,jlo}>, < ulo_{l,ilo}| d/dr + (l+2)/r | ulo_{l+1,jlo}>

                 !write(6,'(A,4I4,1x,12F12.7)') 'sixth:', l, m, jlo, ilo, xpy1, xmy1, z1, xpy2, xmy2, z2
              enddo
           enddo
           ! xpy1 = < psi_{n1}(l+1,m+1)| d/dr-l/r       | psi_{n2}(l,m)   >
           ! xmy1 = < psi_{n1}(l+1,m-1)| d/dr-l/r       | psi_{n2}(l,m)   >
           ! z1   = < psi_{n1}(l+1,m)  | d/dr-l/r       | psi_{n2}(l,m)   >
           ! xpy2 = < psi_{n1}(l,m)    | d/dr + (l+2)/r | psi_{n2}(l+1,m-1) >
           ! xmy2 = < psi_{n1}(l,m)    | d/dr + (l+2)/r | psi_{n2}(l+1,m+1) >
           ! z2   = < psi_{n1}(l,m)    | d/dr + (l+2)/r | psi_{n2}(l+1,m) >
           ! flm1   == -a(l,m)
           ! flm3   ==  a(l,-m)
           ! flm5   ==  f(l,m)
           ! pxpy <= -a(l,m)  < psi_{n1}(l+1,m+1)| d/dr-l/r | psi_{n2}(l,m)> + a(l,-m) < psi_{n1}(l,m) | d/dr+(l+2)/r | psi_{n2}(l+1,m-1) >
           ! pxmy <=  a(l,-m) < psi_{n1}(l+1,m-1)| d/dr-l/r | psi_{n2}(l,m)> - a(l,m)  < psi_{n1}(l,m) | d/dr+(l+2)/r | psi_{n2}(l+1,m+1) >
           ! pz   <=  f(l,m)  < psi_{n1}(l+1,m)  | d/dr-l/r | psi_{n2}(l,m)> + f(l,m)  < psi_{n1}(l,m) | d/dr+(l+2)/r | psi_{n2}(l+1,m) >)
           pxpy = pxpy + flm1 * xpy1 + flm3 * xpy2
           pxmy = pxmy + flm3 * xmy1 + flm1 * xmy2
           pz   = pz   + flm5 * ( z1 + z2 )
           !enddo
           !enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  enddo
  mt12(1) = -(pxpy+pxmy)*imag/2.d0
  mt12(2) = (-pxpy+pxmy)/2.d0
  mt12(3) = -pz*imag
  !write(6,'(I3,1x,I3,1x,A,F16.12,1x,F16.12,2x,F16.12,1x,F16.12,2x,F16.12,1x,F16.12)') n2, n1, 'mt12=', mt12(:)
end subroutine calcmmatvv_MT

subroutine calcmmatvv_MT2(mt12, calfa_n1, cbeta_n1, cgama_n1, alfa_n2, beta_n2, gama_n2, nLO_at, mult, iul_ul, iul_udl, iudl_ul, iudl_udl, iulol_ul, iulol_udl, iul_ulol, iudl_ulol, iulol_ulol, nat, lmax, ntnt, ndf, nLOmax, ntnt_g, ndf_g, nlo, lomx, lomax)
  ! This subroutine calculates < psi_{n1}| -i*nabla | psi_{n2}> -- the contribution of the MT-Spheres of the bands n1 and n2
  implicit none
  complex*16, intent(out):: mt12(3)
  integer,    intent(in) :: nat, lmax, ntnt, ndf, nLOmax, ntnt_g, ndf_g, lomx, nlo, lomax
  integer,    intent(in) :: mult(nat)
  complex*16, intent(in) :: calfa_n1(ntnt,ndf), cbeta_n1(ntnt,ndf), cgama_n1(nLOmax,ntnt_g,ndf_g)
  complex*16, intent(in) :: alfa_n2(ntnt,ndf), beta_n2(ntnt,ndf), gama_n2(nLOmax,ntnt_g,ndf_g)
  integer,    intent(in) :: nLO_at(0:lomax,nat)
  real(8),    intent(in) :: iul_ul(2,0:lmax-1,nat), iul_udl(2,0:lmax-1,nat), iudl_ul(2,0:lmax-1,nat), iudl_udl(2,0:lmax-1,nat)
  real(8),    intent(in) :: iulol_ul(2,nlo,0:lomx,nat), iulol_udl(2,nlo,0:lomx,nat), iul_ulol(2,nlo,0:lomx,nat), iudl_ulol(2,nlo,0:lomx,nat), iulol_ulol(2,nlo,nlo,0:lomx,nat)
  !
  integer    :: iat, ieq, idf
  integer    :: l, m, lm, l1m1, l1mm1, l1m
  integer    :: ilo, jlo
  real(8)    :: denom, flm1, flm3, flm5 
  complex*16 :: xpy1, xpy2, xmy1, xmy2, z1, z2, pxpy, pxmy, pz
  complex*16 :: aa, ab, ba, bb, ag, bg, ga, gb
  complex*16 :: czero, imag
  intrinsic conjg
  intrinsic sqrt
  interface
     integer Function get_lm(l,m)
       integer, intent(in) :: l, m
     end Function get_lm
  end interface
  !--------------------------------------------------------------------!
  !  contribution from MT spheres : < psi_{n1}| i*nabla | psi_{n2}>    !
  !--------------------------------------------------------------------!
  !  DICTIONARY
  ! iul_ul[1,2]     < u_{l+1}  | d/dr-l/r | u_l  >, < u_l  | d/dr + (l+2)/r | u_{l+1}  >
  ! iul_udl[1,2]    < u_{l+1}  | d/dr-l/r |udot_l>, < u_l  | d/dr + (l+2)/r |udot_{l+1}>
  ! iudl_ul[1,2]    <udot_{l+1}| d/dr-l/r | u_l  > ,<udot_l| d/dr + (l+2)/r | u_{l+1}  >
  ! iudl_udl[1,2]   <udot_{l+1}| d/dr-l/r |udot_l>, <udot_l| d/dr + (l+2)/r |udot_{l+1}>
  ! iulol_ul[1,2]   < ulo_{l+1}| d/dr-l/r | u_l  >, < ulo_l| d/dr + (l+2)/r | u_{l+1}  >
  ! iulol_udl[1,2]  < ulo_{l+1}| d/dr-l/r |udot_l>, < ulo_l| d/dr + (l+2)/r |udot_{l+1}>
  ! iul_ulol[1,2]   < u_{l+1}  | d/dr-l/r | ulo_l>, < u_l  | d/dr + (l+2)/r | ulo_{l+1}>
  ! iudl_ulol[1,2]  <udot_{l+1}| d/dr-l/r | ulo_l>, <udot_l| d/dr + (l+2)/r | ulo_{l+1}>
  ! iulol_ulol[1,2] < ulo_{l+1,ilo}| d/dr-l/r | ulo_{l,jlo}>, < ulo_{l,ilo}| d/dr + (l+2)/r | ulo_{l+1,jlo}>
  !
  czero = cmplx(0.d0,0.d0,8)
  imag  = cmplx(0.d0,1.d0,8)
  pxpy = czero
  pxmy = czero
  pz   = czero
  idf = 0
  do iat = 1,nat
     do ieq = 1,mult(iat)
        idf = idf+1
        !$OMP PARALLEL DO PRIVATE(lm,l,m,ilo,jlo,denom,l1mm1,l1m,l1m1,flm1,flm3,flm5,aa,ab,ba,bb,ag,bg,ga,gb,xpy1,xpy2,xmy1,xmy2,z1,z2)&
        !$OMP& REDUCTION(+:pxpy,pxmy,pz)&
        !$OMP& SCHEDULE(STATIC)
        do lm=1,lmax*lmax
           l = int(aint(sqrt(lm-1.d0)+1e-6))
           m = lm-1-l*(l+1)
           denom = dble((2*l+1)*(2*l+3))
           l1mm1 = get_lm(l+1,m-1) ! (l+1,m-1) == (l+1)*(l+2) + m  = lm + 2*l+1
           l1m   = get_lm(l+1,m)   ! (l+1, m ) == (l+1)*(l+2) + m + 1 = lm + 2*l+2
           l1m1  = get_lm(l+1,m+1) ! (l+1,m+1) == (l+1)*(l+2) + m + 2
           flm1 = -sqrt(dble((l+m+1)*(l+m+2))/denom)  !  == -a(l,m)
           flm3 =  sqrt(dble((l-m+1)*(l-m+2))/denom)  !  ==  a(l,-m)
           flm5 =  sqrt(dble((l-m+1)*(l+m+1))/denom)  !  ==  f(l,m)
           !  xpy1 = < psi_{n1}(l+1,m+1) | (d/dr-l/r) | psi_{n2}(l,m) >
           aa = calfa_n1(l1m1,idf)*alfa_n2(lm,idf)  !  A[n1](l+1,m+1)^*  A[n2](l,m)
           ab = calfa_n1(l1m1,idf)*beta_n2(lm,idf)  !  A[n1](l+1,m+1)^*  B[n2](l,m)
           ba = cbeta_n1(l1m1,idf)*alfa_n2(lm,idf)  !  B[n1](l+1,m+1)^*  A[n2](l,m)
           bb = cbeta_n1(l1m1,idf)*beta_n2(lm,idf)  !  B[n1](l+1,m+1)^*  B[n2](l,m)
           !  < A[n1](l+1,m+1) u_{l+1} + B[n1](l+1,m+1) udot_{l+1}| (d/dr-l/r) | A[n2](l,m) u_l + B[n2](l,m) udot_l >
           xpy1 = aa * iul_ul(1,l,iat) + ab * iul_udl(1,l,iat) + ba * iudl_ul(1,l,iat) + bb * iudl_udl(1,l,iat)
           !  xpy2 = < psi_{n1}(l,m)  | d/dr + (l+2)/r | psi_{n2}(l+1,m-1) > 
           aa = calfa_n1(lm,idf)*alfa_n2(l1mm1,idf) !  A[n1](l,m)^* A[n2](l+1,m-1)
           ab = calfa_n1(lm,idf)*beta_n2(l1mm1,idf) !  A[n1](l,m)^* B[n2](l+1,m-1)
           ba = cbeta_n1(lm,idf)*alfa_n2(l1mm1,idf) !  B[n1](l,m)^* A[n2](l+1,m-1)
           bb = cbeta_n1(lm,idf)*beta_n2(l1mm1,idf) !  B[n1](l,m)^* B[n2](l+1,m-1)
           !  < A[n1](l,m) u_l + B[n1](l,m) udot_l  | d/dr + (l+2)/r | A[n2](l+1,m-1) u_{l+1} + B[n2](l+1,m-1) udot_{l+1} >
           xpy2 = aa * iul_ul(2,l,iat) + ab * iul_udl(2,l,iat) + ba * iudl_ul(2,l,iat) + bb * iudl_udl(2,l,iat)
           !  xmy1 = < psi_{n1}(l+1,m-1) | d/dr-l/r | psi_{n2}(l,m)  >
           aa = calfa_n1(l1mm1,idf)*alfa_n2(lm,idf) !  A[n1](l+1,m-1)^* A[n2](l,m)
           ab = calfa_n1(l1mm1,idf)*beta_n2(lm,idf) !  A[n1](l+1,m-1)^* B[n2](l,m)
           ba = cbeta_n1(l1mm1,idf)*alfa_n2(lm,idf) !  B[n1](l+1,m-1)^* A[n2](l,m)
           bb = cbeta_n1(l1mm1,idf)*beta_n2(lm,idf) !  B[n1](l+1,m-1)^* B[n2](l,m)
           !  < A[n1](l+1,m-1) u_{l+1} + B[n1](l+1,m-1) udot _{l+1} | d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l >
           xmy1 = aa * iul_ul(1,l,iat) + ab * iul_udl(1,l,iat) + ba * iudl_ul(1,l,iat) + bb * iudl_udl(1,l,iat)
           !  xmy2 = < psi_{n1}(l,m) | d/dr + (l+2)/r | psi_{n2}(l+1,m+1) >
           aa = calfa_n1(lm,idf)*alfa_n2(l1m1,idf) ! A[n1](l,m)^* A[n2](l+1,m+1)
           ab = calfa_n1(lm,idf)*beta_n2(l1m1,idf) ! A[n1](l,m)^* B[n2](l+1,m+1)
           ba = cbeta_n1(lm,idf)*alfa_n2(l1m1,idf) ! B[n1](l,m)^* A[n2](l+1,m+1)
           bb = cbeta_n1(lm,idf)*beta_n2(l1m1,idf) ! B[n1](l,m)^* B[n2](l+1,m+1)
           !  < A[n1](l,m) u_l +B[n1](l,m) udot_l | d/dr + (l+2)/r | A[n2](l+1,m+1) u_{l+1}  + B[n2](l+1,m+1) udot_{l+1}>
           xmy2 = aa * iul_ul(2,l,iat) + ab * iul_udl(2,l,iat) + ba * iudl_ul(2,l,iat) + bb * iudl_udl(2,l,iat)
           !  z1 = < psi_{n1}(l+1,m) | d/dr-l/r | psi_{n2}(l,m) >
           aa = calfa_n1(l1m,idf)*alfa_n2(lm,idf) !  A[n1](l+1,m)^* A[n2](l,m)
           ab = calfa_n1(l1m,idf)*beta_n2(lm,idf) !  A[n1](l+1,m)^* B[n2](l,m)
           ba = cbeta_n1(l1m,idf)*alfa_n2(lm,idf) !  B[n1](l+1,m)^* A[n2](l,m)
           bb = cbeta_n1(l1m,idf)*beta_n2(lm,idf) !  B[n1](l+1,m)^* B[n2](l,m)
           !  < A[n1](l+1,m) u_{l+1} + B[n1](l+1,m) udot_{l+1}| d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l > 
           z1 = aa * iul_ul(1,l,iat) + ab * iul_udl(1,l,iat) + ba * iudl_ul(1,l,iat) + bb * iudl_udl(1,l,iat)
           !  z2 = < psi_{n1}(l,m) | d/dr + (l+2)/r | psi_{n2}(l+1,m) >
           aa = calfa_n1(lm,idf)*alfa_n2(l1m,idf) !  A[n1](l,m)^* A[n2](l+1,m)
           ab = calfa_n1(lm,idf)*beta_n2(l1m,idf) !  A[n1](l,m)^* B[n2](l+1,m)
           ba = cbeta_n1(lm,idf)*alfa_n2(l1m,idf) !  B[n1](l,m)^* A[n2](l+1,m)
           bb = cbeta_n1(lm,idf)*beta_n2(l1m,idf) !  B[n1](l,m)^* B[n2](l+1,m)
           !  < A[n1](l,m) u_l +B[n1](l,m) udot_l | d/dr + (l+2)/r | A[n2](l+1,m) u_{l+1} + B[n2](l+1,m) udot_{l+1} > 
           z2 = aa * iul_ul(2,l,iat) + ab * iul_udl(2,l,iat) + ba * iudl_ul(2,l,iat) + bb * iudl_udl(2,l,iat)
           !  the loop for the Local orbitals on the right side
           do jlo = 1,nLO_at(l,iat) 
              !  xpy1 += < psi_{n1}(l+1,m+1) | d/dr-l/r | C[n2](l,m) ulo_l>
              ag = calfa_n1(l1m1,idf)*gama_n2(jlo,lm,idf) ! A[n1](l+1,m+1)^* C[n2](l,m)
              bg = cbeta_n1(l1m1,idf)*gama_n2(jlo,lm,idf) ! B[n1](l+1,m+1)^* C[n2](l,m)
              !  < A[n1](l+1,m+1) u_{l+1} + B[n1](l+1,m+1) udot_{l+1} | d/dr-l/r | C[n2](l,m) ulo_l>
              xpy1 = xpy1 + ag * iul_ulol(1,jlo,l,iat) + bg * iudl_ulol(1,jlo,l,iat)
              !  xmy1 += < psi_{n1}(l+1,m-1) | d/dr-l/r | C[n2](l,m) ulo_l>
              ag = calfa_n1(l1mm1,idf)*gama_n2(jlo,lm,idf) !  A[n1](l+1,m-1)^* C[n2](l,m)
              bg = cbeta_n1(l1mm1,idf)*gama_n2(jlo,lm,idf) !  B[n1](l+1,m-1)^* C[n2](l,m)
              !  < A[n1](l+1,m-1) u_{l+1} + B[n1](l+1,m-1) udot_{l+1}| d/dr-l/r | C[n2](l,m) ulo_l>
              xmy1 = xmy1 + ag * iul_ulol(1,jlo,l,iat) + bg * iudl_ulol(1,jlo,l,iat)
              !  z1 += < psi_{n1}(l+1,m) | d/dr-l/r | C[n2](l,m) ulo_l> 
              ag = calfa_n1(l1m,idf)*gama_n2(jlo,lm,idf)  !  A[n1](l+1,m)^* C[n2](l,m)
              bg = cbeta_n1(l1m,idf)*gama_n2(jlo,lm,idf)  !  B[n1](l+1,m)^* C[n2](l,m)
              !  < A[n1](l+1,m) u_{l+1} + B[n1](l+1,m) udot_{l+1} | d/dr-l/r | C[n2](l,m) ulo_l> 
              z1 = z1 + ag * iul_ulol(1,jlo,l,iat) + bg * iudl_ulol(1,jlo,l,iat)
           enddo

           do jlo = 1,nLO_at(l+1,iat) 
              !  xpy2 += < psi_{n1}(l,m) | d/dr + (l+2)/r | C[n2](l+1,m-1) ulo_{l+1}>
              ag = calfa_n1(lm,idf)*gama_n2(jlo,l1mm1,idf)  !  A[n1](l,m)^*  C[n2](l+1,m-1)
              bg = cbeta_n1(lm,idf)*gama_n2(jlo,l1mm1,idf)  !  B[n1](l,m)^*  C[n2](l+1,m-1)
              ! < A[n1](l,m) u_l  + B[n1](l,m) udot_l | d/dr + (l+2)/r | C[n2](l+1,m-1) ulo_{l+1}>
              xpy2 = xpy2 + ag * iul_ulol(2,jlo,l,iat) + bg * iudl_ulol(2,jlo,l,iat)
              !  xmy2 += < psi_{n1}(l,m) | d/dr + (l+2)/r | C[n2](l+1,m+1) ulo_{l+1}>
              ag = calfa_n1(lm,idf)*gama_n2(jlo,l1m1, idf)  !  A[n1](l,m)^*  C[n2](l+1,m+1)
              bg = cbeta_n1(lm,idf)*gama_n2(jlo,l1m1, idf)  !  B[n1](l,m)^*  C[n2](l+1,m+1)
              !  < A[n1](l,m) u_l + B[n1](l,m) udot_l | d/dr + (l+2)/r | C[n2](l+1,m+1) ulo_{l+1}>
              xmy2 = xmy2 + ag * iul_ulol(2,jlo,l,iat) + bg * iudl_ulol(2,jlo,l,iat)
              !KH : Introducing a bug to be compatible with pygab
              !xmy2 = xmy1 + ag * iul_ulol(2,jlo,l,iat) + bg * iudl_ulol(2,jlo,l,iat)

              !  z2 += < psi_{n1}(l,m) | d/dr + (l+2)/r | C[n2](l+1,m) ulo_{l+1}>
              ag = calfa_n1(lm,idf)*gama_n2(jlo,l1m,idf)  ! A[n1](l,m)^*  C[n2](l+1,m)
              bg = cbeta_n1(lm,idf)*gama_n2(jlo,l1m,idf)  ! B[n1](l,m)^*  C[n2](l+1,m)
              !  < A[n1](l,m) u_l + B[n1](l,m) udot_l | d/dr + (l+2)/r | C[n2](l+1,m) ulo_{l+1}>
              z2 = z2 + ag * iul_ulol(2,jlo,l,iat) + bg * iudl_ulol(2,jlo,l,iat)
           enddo
           !  the loop for the Local orbitals on the left side 
           do ilo = 1,nLO_at(l+1,iat)
              !  xpy1 += < C[n1](l+1,m+1) ulo_{l+1}| d/dr-l/r | psi_{n2}(l,m) >   
              ga = cgama_n1(ilo,l1m1,idf)*alfa_n2(lm,idf)  ! C[n1](l+1,m+1)^*  A[n2](l,m)
              gb = cgama_n1(ilo,l1m1,idf)*beta_n2(lm,idf)  ! C[n1](l+1,m+1)^*  B[n2](l,m)
              !  < C[n1](l+1,m+1) ulo_{l+1}| d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l >
              xpy1 = xpy1 + ga * iulol_ul(1,ilo,l,iat) + gb * iulol_udl(1,ilo,l,iat)
              !  xmy1 +=  C[n1](l+1,m-1) ulo_{l+1}| d/dr-l/r | psi_{n2}(l,m) >
              ga = cgama_n1(ilo,l1mm1,idf)*alfa_n2(lm,idf)  !  C[n1](l+1,m-1)^*  A[n2](l,m)
              gb = cgama_n1(ilo,l1mm1,idf)*beta_n2(lm,idf)  !  C[n1](l+1,m-1)^*  B[n2](l,m)
              !  < C[n1](l+1,m-1) ulo_{l+1}| d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l >
              xmy1 = xmy1 + ga * iulol_ul(1,ilo,l,iat) + gb * iulol_udl(1,ilo,l,iat)
              !  z1 += < C[n1](l+1,m) ulo_{l+1}| d/dr-l/r | psi_{n2}(l,m) >
              ga = cgama_n1(ilo,l1m,idf)*alfa_n2(lm,idf)  !  C[n1](l+1,m)^*  A[n2](l,m)
              gb = cgama_n1(ilo,l1m,idf)*beta_n2(lm,idf)  !  C[n1](l+1,m)^*  B[n2](l,m)
              !  < C[n1](l+1,m) ulo_{l+1}| d/dr-l/r | A[n2](l,m) u_l + B[n2](l,m) udot_l >
              z1 = z1 + ga * iulol_ul(1,ilo,l,iat) + gb * iulol_udl(1,ilo,l,iat)
           enddo
           do ilo = 1,nLO_at(l,iat)
              !  xpy2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | psi_{n2}(l+1,m-1) >
              ga = cgama_n1(ilo,lm,idf)*alfa_n2(l1mm1,idf) !  C[n1](l,m)^*  A[n2](l+1,m-1)
              gb = cgama_n1(ilo,lm,idf)*beta_n2(l1mm1,idf) !  C[n1](l,m)^*  B[n2](l+1,m-1)
              !  < C[n1](l,m) ulo_l| d/dr + (l+2)/r | A[n2](l+1,m-1) u_{l+1} + B[n2](l+1,m-1) udot_{l+1} >
              xpy2 = xpy2 + ga * iulol_ul(2,ilo,l,iat) + gb * iulol_udl(2,ilo,l,iat)
              !  xmy2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | psi_{n2}(l+1,m+1) >
              ga = cgama_n1(ilo,lm,idf)*alfa_n2(l1m1, idf)  !  C[n1](l,m)^*  A[n2](l+1,m+1)
              gb = cgama_n1(ilo,lm,idf)*beta_n2(l1m1, idf)  !  C[n1](l,m)^*  B[n2](l+1,m+1)
              !  < C[n1](l,m) ulo_l| d/dr + (l+2)/r | A[n2](l+1,m+1) u_{l+1} +  B[n2](l+1,m+1) udot_{l+1} >
              xmy2 = xmy2 + ga * iulol_ul(2,ilo,l,iat) + gb * iulol_udl(2,ilo,l,iat)
              ! z2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | psi_{n2}(l+1,m) >
              ga = cgama_n1(ilo,lm,idf)*alfa_n2(l1m,idf)  !  C[n1](l,m)^*  A[n2](l+1,m)
              gb = cgama_n1(ilo,lm,idf)*beta_n2(l1m,idf)  !  C[n1](l,m)^*  B[n2](l+1,m)
              !  < C[n1](l,m) ulo_l| d/dr + (l+2)/r | A[n2](l+1,m) u_{l+1} + B[n2](l+1,m) udot_{l+1} >
              z2 = z2 + ga * iulol_ul(2,ilo,l,iat) + gb * iulol_udl(2,ilo,l,iat)
           enddo
           do jlo = 1,nLO_at(l,iat)
              do ilo = 1, nLO_at(l+1,iat)
                 ! xpy1 += < C[n1](l+1,m+1) ulo_{l+1}| d/dr-l/r | C[n2](l,m) ulo_l>
                 xpy1 = xpy1 + cgama_n1(ilo,l1m1, idf)*gama_n2(jlo,lm,   idf) * iulol_ulol(1,ilo,jlo,l,iat)
                 ! xmy1 += < C[n1](l+1,m-1) ulo_{l+1}| d/dr-l/r | C[n2](l,m) ulo_l>
                 xmy1 = xmy1 + cgama_n1(ilo,l1mm1,idf)*gama_n2(jlo,lm,   idf) * iulol_ulol(1,ilo,jlo,l,iat)
                 ! z1 += < C[n1](l+1,m) ulo_{l+1}| d/dr-l/r | C[n2](l,m) ulo_l>
                 z1   = z1   + cgama_n1(ilo,l1m,  idf)*gama_n2(jlo,lm,   idf) * iulol_ulol(1,ilo,jlo,l,iat)
                 ! xpy2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | C[n2](l+1,m-1) ulo_{l+1}>
                 xpy2 = xpy2 + cgama_n1(jlo,lm,   idf)*gama_n2(ilo,l1mm1,idf) * iulol_ulol(2,jlo,ilo,l,iat)
                 ! xmy2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | C[n2](l+1,m+1) ulo_{l+1}>
                 xmy2 = xmy2 + cgama_n1(jlo,lm,   idf)*gama_n2(ilo,l1m1, idf) * iulol_ulol(2,jlo,ilo,l,iat)
                 ! z2 += < C[n1](l,m) ulo_l| d/dr + (l+2)/r | C[n2](l+1,m) ulo_{l+1}>
                 z2   = z2   + cgama_n1(jlo,lm,   idf)*gama_n2(ilo,l1m,  idf) * iulol_ulol(2,jlo,ilo,l,iat)
                 ! iulol_ulol[1,2] < ulo_{l+1,ilo}| d/dr-l/r | ulo_{l,jlo}>, < ulo_{l,ilo}| d/dr + (l+2)/r | ulo_{l+1,jlo}>

                 !write(6,'(A,4I4,1x,12F12.7)') 'sixth:', l, m, jlo, ilo, xpy1, xmy1, z1, xpy2, xmy2, z2
              enddo
           enddo
           ! xpy1 = < psi_{n1}(l+1,m+1)| d/dr-l/r       | psi_{n2}(l,m)   >
           ! xmy1 = < psi_{n1}(l+1,m-1)| d/dr-l/r       | psi_{n2}(l,m)   >
           ! z1   = < psi_{n1}(l+1,m)  | d/dr-l/r       | psi_{n2}(l,m)   >
           ! xpy2 = < psi_{n1}(l,m)    | d/dr + (l+2)/r | psi_{n2}(l+1,m-1) >
           ! xmy2 = < psi_{n1}(l,m)    | d/dr + (l+2)/r | psi_{n2}(l+1,m+1) >
           ! z2   = < psi_{n1}(l,m)    | d/dr + (l+2)/r | psi_{n2}(l+1,m) >
           ! flm1   == -a(l,m)
           ! flm3   ==  a(l,-m)
           ! flm5   ==  f(l,m)
           ! pxpy <= -a(l,m)  < psi_{n1}(l+1,m+1)| d/dr-l/r | psi_{n2}(l,m)> + a(l,-m) < psi_{n1}(l,m) | d/dr+(l+2)/r | psi_{n2}(l+1,m-1) >
           ! pxmy <=  a(l,-m) < psi_{n1}(l+1,m-1)| d/dr-l/r | psi_{n2}(l,m)> - a(l,m)  < psi_{n1}(l,m) | d/dr+(l+2)/r | psi_{n2}(l+1,m+1) >
           ! pz   <=  f(l,m)  < psi_{n1}(l+1,m)  | d/dr-l/r | psi_{n2}(l,m)> + f(l,m)  < psi_{n1}(l,m) | d/dr+(l+2)/r | psi_{n2}(l+1,m) >)
           pxpy = pxpy + flm1 * xpy1 + flm3 * xpy2
           pxmy = pxmy + flm3 * xmy1 + flm1 * xmy2
           pz   = pz   + flm5 * ( z1 + z2 )
           !enddo
           !enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
  enddo
  mt12(1) = -(pxpy+pxmy)*imag/2.d0
  mt12(2) = (-pxpy+pxmy)/2.d0
  mt12(3) = -pz*imag
  !write(6,'(I3,1x,I3,1x,A,F16.12,1x,F16.12,2x,F16.12,1x,F16.12,2x,F16.12,1x,F16.12)') n2, n1, 'mt12=', mt12(:)
end subroutine calcmmatvv_MT2

subroutine calcmmatvv_is(mpw, i,j, kl, Ak, indgk, nvk, ipwint, gindex, i_g0, k2cartes, n1, n2, n3, ngk, ng, nbmax, ngi)
  ! Computing interstitial contribution of
  !          < psi_{i}| (-i*\nabla) | psi_{j}>_I
  !                       = \sum_{G1,G2} A*_{G1,i} A_{G2,j} <G1+k| -i*\nabla |G2+k>_I =
  !                       = \sum_{G1,G2} A*_{G1,i} A_{G2,j} (G2+k) <G1|G2>_I
  implicit none
  complex*16, intent(out) :: mpw(3)
  integer,    intent(in)  :: i, j, nvk, n1, n2, n3, nbmax, ngi, ng, ngk
  real(8),    intent(in)  :: k2cartes(3,3)
  real(8),    intent(in)  :: kl(3)                 ! k-point in semi-cartesian form
  integer,    intent(in)  :: gindex(ng,3)          ! all G-points 
  integer,    intent(in)  :: indgk(ngk)
  integer,    intent(in)  :: i_g0(2*n1+1,2*n2+1,2*n3+1)
  complex*16, intent(in)  :: ipwint(ng)
  complex*16, intent(in)  :: Ak(nbmax,ngi)
  !
  real(8) :: kcar(3)
  integer :: ig2, ig1, iGl2(3), iGl1(3), iG(3), idg
  real(8) :: Gc2(3), kpG2(3)
  complex*16 :: intmom, modmom
  !
  kcar = matmul(k2cartes, kl) ! change k-point to cartesian coordinates
  !print *, 'nvk=', nvk, 'ngk=', ngk, 'i=', i, 'j=', j
  !do ig2=1,nvk
  !   write(6,'(I4,1x,I4,1x,2F10.6,2x,2F10.6)') i, ig2, Ak(i,ig2), Ak(j,ig2)
  !enddo
  !write(6,*)
  
  mpw = 0.d0
  do ig2=1,nvk
     iGl2 = gindex(indgk(ig2)+1,:)
     Gc2 = matmul(k2cartes,iGl2)
     kpG2(:) = kcar(:) + Gc2(:)

     !write(6,'(I4,1x,3I3,1x,3F10.4,1x,3F10.4)') ig2, iGl2, Gc2, kpG2
     intmom = 0.d0
     do ig1=1,nvk
        iGl1 = gindex(indgk(ig1)+1,:)
        iG = iGl2-iGl1
        idg = i_g0(iG(1)+n1+1,iG(2)+n2+1,iG(3)+n3+1)
        intmom = intmom + ipwint(idg+1) * conjg(Ak(i,ig1))
        !write(6,'(A,I4,1x,3I3,1x,I4,1x,2F10.6,1x,2F10.6)') '     ', ig1, iG, idg, ipwint(idg+1), Ak(i,ig1)
     enddo
     modmom = Ak(j,ig2)*intmom
     mpw(:) = mpw(:) + kpG2(:)*modmom
  enddo
  !write(6,'(I3,1x,I3,1x,A,F16.12,1x,F16.12,2x,F16.12,1x,F16.12,2x,F16.12,1x,F16.12)') j, i, 'mpw=', mpw(:)
end subroutine calcmmatvv_is


subroutine calcmmatcv(mcv, iat, iv, corind, alfa,beta,gama, iul_ucl,iudl_ucl,iucl_ul,iucl_udl,iulol_ucl,iucl_ulol,nLO_at, ncg,nat,nloat,nbmax,ndf,ntnt,nbmax_g,nLOmax,ntnt_g,ndf_g,ncg_a,lomax)
  ! This subroutine calculates the momentum matrix element
  ! <psi| -i\nabla | psi>
  ! between the core and the valence states.
  IMPLICIT NONE
  integer,    intent(in) :: iv, iat, ncg, ntnt, nbmax, ndf, nbmax_g, nLOmax, ntnt_g, ndf_g, nloat, lomax, nat, ncg_a
  complex*16, intent(out):: mcv(ncg,3)
  integer,    intent(in) :: corind(ncg,5)
  real*8,     intent(in) :: iul_ucl(2,ncg_a), iudl_ucl(2,ncg_a), iucl_ul(2,ncg_a), iucl_udl(2,ncg_a)
  real*8,     intent(in) :: iulol_ucl(2,nloat,ncg_a), iucl_ulol(2,nloat,ncg_a)
  complex*16, intent(in) :: alfa(nbmax,ntnt,ndf)
  complex*16, intent(in) :: beta(nbmax,ntnt,ndf)
  complex*16, intent(in) :: gama(nbmax_g,nLOmax,ntnt_g,ndf_g)
  integer,    intent(in) :: nLO_at(0:lomax,nat)
  !
  interface
     integer Function get_lm(l,m)
       integer, intent(in) :: l, m
     end Function get_lm
  end interface
  !
  complex*16 :: mxcv, mxvc, mycv, myvc, mzcv, mzvc
  integer :: lm, l1m1, l1mm1, l1m,lm1m,lm1mm1,lm1m1
  integer :: iatt, l, m, ilo, ic, idf, icg
  real(8) :: denom1, denom2
  real(8) :: flm1, flm2, flm3, flm4, flm5, flm6
  complex*16 :: xpy1,xpy2,xmy1,xmy2,z1,z2
  complex*16 :: pxpy,pxmy,pz, imag
  !
  !  All matrix elements we need here:
  ! iul_ucl  [1-2,ic]    = [ iul1ucl  ,iulucl1  ] = [<   u^v_{l+1}| d/dr -l/r |u^c_l>, <   u^v_{l-1}| d/dr +(l+1)/r |u^c_l>]
  ! iudl_ucl [1-2,ic]    = [ iudl1ucl ,iudlucl1 ] = [<udot^v_{l+1}| d/dr -l/r |u^c_l>, <udot^v_{l-1}| d/dr +(l+1)/r |u^c_l>]
  ! iucl_ul  [1-2,ic]    = [ iucl1ul  ,iuclul1  ] = [< u^c_l | d/dr -(l-1)/r|   u^v_{l-1}>, < u^c_l | d/dr +(l+2)|   u^v_{l+1}>]
  ! iucl_udl [1-2,ic]    = [ iucl1udl ,iucludl1 ] = [< u^c_l | d/dr -(l-1)/r|udot^v_{l-1}>, < u^c_l | d/dr +(l+2)|udot^v_{l+1}>]
  !
  ! iulol_ucl[1-2,ic,ilo]= [ iulol1ucl,iulolucl1] = [<ulo_{l+1}| d/dr - l/r |u^c_l>, <ulo_{l-1}| d/dr+(l+1)/r |u^c_l>]
  ! iucl_ulol[1-2,ic,ilo]= [ iucl1ulol,iuclulol1] = [< u^c_l | d/dr -(l-1)/r|   ulo_{l-1}>, < u^c_l | d/dr +(l+2)|   ulo_{l+1}>]
  ! iucl_ucl [1-2,ic,jc] = [ iucl1ucl ,iuclucl1 ] = [< u^c_l | d/dr -(l-1)/r|   u^c_{l-1}>, < u^c_l | d/dr +(l+2)/r |u^c_{l+1}>]
  imag  = cmplx(0.d0,1.d0,8)
  do icg = 1, ncg
     xpy1 = 0.d0
     xmy1 = 0.d0
     z1   = 0.d0
     iatt = corind(icg,1)+1
     if (iatt .ne. iat) cycle  !  only over one atom at a time
     idf = corind(icg,2)+1
     ic  = corind(icg,3)+1
     l  = corind(icg,4)
     m  = corind(icg,5)
     denom1 = dble((2*l+1)*(2*l+3))
     denom2 = dble((2*l-1)*(2*l+1)) 
     lm    = get_lm(l,m)     ! (l,  m  ) == l*(l+1) + m + 1
     l1mm1 = get_lm(l+1,m-1) ! (l+1,m-1) == (l+1)*(l+2) + m  = lm + 2*l+1
     l1m   = get_lm(l+1,m)   ! (l+1, m ) == (l+1)*(l+2) + m + 1 = lm + 2*l+2
     l1m1  = get_lm(l+1,m+1) ! (l+1,m+1) == (l+1)*(l+2) + m + 2
     lm1mm1= get_lm(l-1,m-1) ! 
     lm1m  = get_lm(l-1,m)
     lm1m1 = get_lm(l-1,m+1)
     !
     flm1 = -sqrt(dble((l+m+1)*(l+m+2))/denom1) ! -a(l,m)
     flm2 = sqrt(dble((l-m-1)*(l-m))/denom2)    !  a(l-1,-m-1)
     flm3 = sqrt(dble((l-m+1)*(l-m+2))/denom1)  !  a(l,-m)
     flm4 = -sqrt(dble((l+m-1)*(l+m))/denom2)   ! -a(l-1,m-1) 
     flm5 = sqrt(dble((l-m+1)*(l+m+1))/denom1)  !  f(l,m)
     flm6 = sqrt(dble((l-m)*(l+m))/denom2)      !  f(l-1,m)
     !       calculate the matrix elements <core/p/valence>
     !  -(l-1)<= m-1 <= l-1 or 0<=l-1+m-1
     if ((l > 0) .and. (l+m-2 >= 0)) then
        ! xpy1 = < u^c_l | d/dr -(l-1)/r |   alfa[iv](l-1,m-1) u^v_{l-1} + beta[iv](l-1,m-1) udot^v_{l-1} + gama[iv](l-1,m-1) ulo_{l-1} >
        xpy1 = alfa(iv,lm1mm1,idf) * iucl_ul(1,ic) + beta(iv,lm1mm1,idf) * iucl_udl(1,ic)
        do ilo=1,nLO_at(l-1,iat) 
           xpy1 = xpy1 + gama(iv,ilo,lm1mm1,idf) * iucl_ulol(1,ilo,ic)
        enddo
     endif
     ! -(l-1)<= m+1 <= l-1  or l-m-2>=0
     if ((l > 0) .and. (l-m-2 >= 0)) then
        ! xmy1 = < u^c_l | d/dr -(l-1)/r| alfa[iv](l-1,m+1) u^v_{l-1} + beta[iv](l-1,m+1) udot_{l-1} + gama[iv](l-1,m+1) ulo_{l-1} >
        xmy1 = alfa(iv,lm1m1,idf) * iucl_ul(1,ic) + beta(iv,lm1m1,idf) * iucl_udl(1,ic)
        do ilo=1,nLO_at(l-1,iat)
           xmy1 = xmy1 + gama(iv,ilo,lm1m1,idf) * iucl_ulol(1,ilo,ic)
        enddo
     endif
     if ((l > 0) .and. (l+m-1 >= 0) .and. (l-m-1 >= 0)) then
        ! z1 = < u^c_l | d/dr -(l-1)/r| alfa[iv](l-1,m) u^v_{l-1} + beta[iv](l-1,m) udot^v_{l-1} + gama[iv](l-1,m) ulo_{l-1} >
        z1 = alfa(iv,lm1m,idf) * iucl_ul(1,ic) + beta(iv,lm1m,idf) * iucl_udl(1,ic)
        do ilo=1,nLO_at(l-1,iat) 
           z1 = z1 + gama(iv,ilo,lm1m,idf) * iucl_ulol(1,ilo,ic) 
        enddo
     endif
     ! xpy2 = < u^c_l | d/dr +(l+2)/r | alfa[iv](l+1,m-1) u^v_{l+1} + beta[iv](l+1,m-1) udot_{l+1} + gama[iv](l+1,m-1) ulo_{l+1}>
     xpy2 = alfa(iv,l1mm1,idf) * iucl_ul(2,ic) + beta(iv,l1mm1,idf) * iucl_udl(2,ic)
     ! xmy2 = < u^c_l | d/dr +(l+2)/r | alfa[iv](l+1,m+1) u^v_{l+1}> + beta[iv](l+1,m+1) udot^v_{l+1} + gama[iv](l+1,m+1) ulo_{l+1} >
     xmy2 = alfa(iv,l1m1,idf) * iucl_ul(2,ic) + beta(iv,l1m1,idf) * iucl_udl(2,ic)
     !  z2  = < u^c_l | d/dr +(l+2)/r | alfa[iv](l+1,m) u^v_{l+1}> + beta[iv](l+1,m) udot^v_{l+1} + gama[iv](l+1,m) ulo_{l+1} >
     z2 = alfa(iv,l1m,idf) * iucl_ul(2,ic) + beta(iv,l1m,idf) * iucl_udl(2,ic)
     do ilo=1,nLO_at(l+1,iat) 
        xpy2 = xpy2 + gama(iv,ilo,l1mm1,idf)*iucl_ulol(2,ilo,ic)
        xmy2 = xmy2 + gama(iv,ilo,l1m1, idf)*iucl_ulol(2,ilo,ic)
        z2   = z2   + gama(iv,ilo,l1m,  idf)*iucl_ulol(2,ilo,ic)
     enddo
     ! pxpy =-a(l-1, m-1) < u^c_l | d/dr -(l-1)/r | psi[iv](l-1,m-1)> + a(l,-m) < u^c_l | d/dr +(l+2)/r | psi[iv](l+1,m-1) >
     ! pxmy = a(l-1,-m-1) < u^c_l | d/dr -(l-1)/r | psi[iv](l-1,m+1)> - a(l, m) < u^c_l | d/dr +(l+2)/r | psi[iv](l+1,m+1) >
     !  pz  =    f(l-1,m) < u^c_l | d/dr -(l-1)/r | psi[iv](l-1, m) > + f(l, m) < u^c_l | d/dr +(l+2)/r | psi[iv](l+1, m ) >
     pxpy = flm4 * xpy1 + flm3 * xpy2
     pxmy = flm2 * xmy1 + flm1 * xmy2
     pz   = flm6 * z1   + flm5 * z2 
     !
     !write(6,'(A,2F14.10,1x,A,2F14.10,1x,A,2F14.10)') 'a=', alfa(iv,l1mm1,idf), 'b=', beta(iv,l1mm1,idf), 'iucl2=', iucl_ul(2,ic), iucl_udl(2,ic)
     !do ilo=1,nLO_at(l+1,iat) 
     !   write(6,'(A,2F14.10,1x,A,F14.10)') 'c=', gama(iv,ilo,l1mm1,idf), 'iucl2=', iucl_ulol(2,ilo,ic)
     !enddo
     !write(6,'(A,4I3,6F16.12,A,6F16.12)') 'xpy2,xmy2,z2', icg, ic, l, m, xpy2, xmy2, z2, ' : ', pxpy, pxmy, pz
     !
     mxcv = -(pxpy+pxmy)*imag/2.d0
     mycv = -(pxpy-pxmy)/2.d0
     mzcv = -imag*pz
     !
     xpy2 = 0.d0
     xmy2 = 0.d0
     z2   = 0.d0
     !       calculate the matrix elements <valence/p/core>
     ! xpy1 = < alfa[iv](l+1,m+1) u^v_{l+1} + beta[iv](l+1,m+1) udot_{l+1} + gama[iv](l+1,m+1) ulo_{l+1} | d/dr -l/r |u^c_l>
     xpy1 = conjg(alfa(iv,l1m1,idf))*iul_ucl(1,ic) + conjg(beta(iv,l1m1,idf))*iudl_ucl(1,ic)
     ! xmy1 = < alfa[iv](l+1,m-1) u^v_{l+1} + beta[iv](l+1,m-1) udot_{l+1} + gama[iv](l+1,m-1) ulo_{l+1} | d/dr -l/r |u^c_l>
     xmy1 = conjg(alfa(iv,l1mm1,idf))*iul_ucl(1,ic) + conjg(beta(iv,l1mm1,idf))*iudl_ucl(1,ic)
     ! z1   = < alfa[iv](l+1,m  ) u^v_{l+1} + beta[iv](l+1,m  ) udot_{l+1} + gama[iv](l+1,m  ) ulo_{l+1} | d/dr -l/r |u^c_l>
     z1 = conjg(alfa(iv,l1m,idf))*iul_ucl(1,ic) + conjg(beta(iv,l1m,idf))*iudl_ucl(1,ic)
     do ilo=1,nLO_at(l+1,iat)  
        xpy1 = xpy1 + conjg(gama(iv,ilo,l1m1, idf))*iulol_ucl(1,ilo,ic)
        xmy1 = xmy1 + conjg(gama(iv,ilo,l1mm1,idf))*iulol_ucl(1,ilo,ic)
        z1   = z1   + conjg(gama(iv,ilo,l1m,  idf))*iulol_ucl(1,ilo,ic)
     enddo
     if ((l > 0) .and. (l-m-2 >= 0)) then
        ! xpy2 = < alfa[iv](l-1,m+1) u^v_{l-1} + beta[iv](l-1,m+1) udot_{l-1} + gama[iv](l-1,m+1) ulo_{l-1} | d/dr +(l+1)/r |u^c_l>
        xpy2 = conjg(alfa(iv,lm1m1,idf))*iul_ucl(2,ic) + conjg(beta(iv,lm1m1,idf))*iudl_ucl(2,ic)
        do ilo=1,nLO_at(l-1,iat) 
           xpy2 = xpy2 + conjg(gama(iv,ilo,lm1m1,idf))*iulol_ucl(2,ilo,ic)
        enddo
     endif
     ! We need -(l-1) <= m-1 <= l-1 or 0<=l-1+m-1
     if ((l > 0) .and. (l+m-2 >= 0)) then
        ! xmy2 = < alfa[iv](l-1,m-1) u^v_{l-1} + beta[iv](l-1,m-1) udot_{l-1} + gama[iv](l-1,m-1) ulo_{l-1} | d/dr +(l+1)/r |u^c_l>
        xmy2 = conjg(alfa(iv,lm1mm1,idf))*iul_ucl(2,ic) + conjg(beta(iv,lm1mm1,idf))*iudl_ucl(2,ic)
        do ilo=1,nLO_at(l-1,iat) 
           xmy2 = xmy2 + conjg(gama(iv,ilo,lm1mm1,idf))*iudl_ucl(2,ic)
        enddo
     endif
     ! We need -(l-1)<=m<=l-1 => 0<=l+m-1 and l-m-1>=0
     if ((l > 0) .and. (l+m-1 >= 0) .and. (l-m-1 >= 0)) then
        ! z2 = < alfa[iv](l-1,m) u^v_{l-1} + beta[iv](l-1,m) udot_{l-1} + gama[iv](l-1,m) ulo_{l-1} | d/dr +(l+1)/r |u^c_l>
        z2 = conjg(alfa(iv,lm1m,idf))*iul_ucl(2,ic) + conjg(beta(iv,lm1m,idf))*iudl_ucl(2,ic)
        do ilo=1,nLO_at(l-1,iat) 
           z2 = z2 + conjg(gama(iv,ilo,lm1m,idf))*iulol_ucl(2,ilo,ic)
        enddo
     endif
     ! pxpy =-a(l,m)  < psi[iv](l+1,m+1) | d/dr -l/r |u^c_l> + a(l-1,-m-1) < psi[iv](l-1,m+1) | d/dr +(l+1)/r |u^c_l>
     ! pxmy = a(l,-m) < psi[iv](l+1,m-1) | d/dr -l/r |u^c_l> - a(l-1,m-1)  < psi[iv](l-1,m-1) | d/dr +(l+1)/r |u^c_l>
     !  pz  = f(l,m)  < psi[iv](l+1,  m) | d/dr -l/r |u^c_l> + f(l-1,  m)  < psi[iv](l-1,  m) | d/dr +(l+1)/r |u^c_l>
     pxpy = flm1 * xpy1 + flm2 * xpy2
     pxmy = flm3 * xmy1 + flm4 * xmy2
     pz   = flm5 * z1   + flm6 * z2
     !
     !write(6,'(A,4I3,6F16.12,A,6F16.12)') 'xpy3,xmy3,z3', icg, ic, l, m, xpy1, xmy1, z1, ' : ', pxpy, pxmy, pz
     !
     mxvc = -(pxpy+pxmy)*imag/2.d0
     myvc = -(pxpy-pxmy)/2.d0
     mzvc = -imag*pz
     mcv(icg,1) = (mxcv+conjg(mxvc))/2.d0
     mcv(icg,2) = (mycv+conjg(myvc))/2.d0
     mcv(icg,3) = (mzcv+conjg(mzvc))/2.d0
     !write(6,'(A,4I3,6F16.12,A,6F16.12)') 'xpy4,xmy4,z4', icg, ic, l, m, mxvc, myvc, mzvc, ' : ', mcv(icg,:)
     !write(6,'(I3,1x,I3,1x,F16.12,1x,F16.12,2x,F16.12,1x,F16.12,2x,F16.12,1x,F16.12,2x,F16.12)') iv, icg, mcv(icg,:), real(sum(mcv(icg,:)*conjg(mcv(icg,:)))/3.d0)
  enddo ! icg
end subroutine calcmmatcv
