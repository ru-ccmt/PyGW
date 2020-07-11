subroutine mt_vxcnn(Vxcmt, irk, uxcu, alfa, beta, gama, cgcoef, lmxc, lxcm, ibgw, nbgw, mult, iatnr, nLO_at, lmax, nLOmax, nat, nt, ndf, nbnd, lomax, maxlxc, nrmax2)
  implicit none
  complex*16,intent(out):: Vxcmt(nbgw,nbgw)
  real(8),   intent(in) :: uxcu(nat,maxlxc,nt*nt,nrmax2)                ! radial matrix elements of Vxc
  complex*16,intent(in) :: alfa(nbnd,nt*nt,ndf), beta(nbnd,nt*nt,ndf), gama(nbnd,nLOmax,nt*nt,ndf)  ! <psi_{KS}|u_alpha>
  real*8,    intent(in) :: cgcoef(*)                                     ! gaunt coefficients
  integer,   intent(in) :: lmxc(2,maxlxc,nat), lxcm(nat), iatnr(nat)       ! index arrays for storing uxcu
  integer,   intent(in) :: ibgw, nbgw                                    ! cutoff in band space
  integer,   intent(in) :: mult(nat)  
  integer,   intent(in) :: nLO_at(lomax+1,nat)
  integer,   intent(in) :: nat, nt, ndf, nLOmax, lomax, nbnd
  integer,   intent(in) :: lmax, maxlxc, nrmax2, irk
  !
  interface
     REAL*8 function getcgcoef(l1,l2,l3,m1,m2,cgcoef)
       REAL*8, intent(in) :: cgcoef(*)
       integer(4), intent(in) :: l1, l2, l3, m1, m2
     end function getcgcoef
  end interface
  !
  integer    :: idf, iat, ieq, ie1, ie2, isgn, ir1, ir2, ir12, nr1, nr2
  integer    :: Lb, mb, lxc, Lx, l1, m1, l2, m2, l12, l1m1, l2m2
  real(8)    :: uxcuw
  complex*16 :: io1, imag, angint, af, caf
  !
  Vxcmt = 0.d0
  imag = cmplx(0.d0,1.d0,8)
  !
  idf=0
  do iat = 1, nat                        !*  Loop over inequivalent atoms:
     do ieq = 1, mult(iat)               !*  Loop over equivalent atoms:
        idf = idf + 1
        do l2=0,lmax                     !*  loop over l2
           do l1=0,lmax                  !*  Loop over l1
              if(l1 > lomax) then 
                 nr1 = 2
              else
                 nr1 = nLO_at(l1+1,iat) + 2
              endif
              if(l2 > lomax) then 
                 nr2 = 2
              else 
                 nr2 = nLO_at(l2+1,iat) + 2 
              endif
              l12 = l1+1 + l2*(lmax+1)             !* just combined index of (l1, l2)
              do m1=-l1,l1                         !* loop over m1
                 l1m1 = l1*l1+l1+m1+1              !* combined (l1,m1) index
                 do lxc = 1, lxcm(iat)             !*  Loop over vxc components: its L value first
                    Lx = lmxc(1,lxc,iat)           !*
                    Lb = abs(Lx)                   !*  L,m = (Lb,mb) value of the Vxc
                    mb = lmxc(2,lxc,iat)
                    io1 = 1.d0
                    if (Lx < 0) io1 = imag
                    if ( .not.( (abs(Lb-l2) <= l1) .and. (l1 <= Lb+l2)) ) cycle !* Triangular rule for sum of Lb+l2 = l1 in quantum vector form
                    do isgn = 1, -1, -2 
                       mb = mb * isgn
                       m2 = m1 - mb                 !* <m1 | V_mb | m2>
                       if( abs(m2) > l2 ) cycle 
                       l2m2 = l2*l2 + l2 + m2 + 1   !* combined index for (l2,m2)
                       !* Calculate the angular integral  <Y_{l1m1}|V_{lb,mb}|Y_{l2,m2}>
                       angint = getcgcoef(Lb,l2,l1,mb,m2,cgcoef) ! gaunt coefficients == <Y_{l1m1}|Y_{lb,mb}|Y_{l2,m2}>
                       if (abs(angint).lt.1.0d-12) cycle 
                       if ( iatnr(iat) > 0) then    !* this means this atom has cubic environment, and we transformed potential to cubic
                          if ( MOD(Lb,2) == 1 ) then!* means Lb in [3,7,9]
                             if (mb < 0) then       !* In cubic symmetry, we can only have l=(0,3,4,6,7,8,9,10...)
                                angint = -angint * imag
                             else
                                angint = angint * imag
                             endif
                          endif
                       else                         !* non-cubic environment on this atom
                          if (mb < 0) then
                             angint = angint * io1/sqrt(2.0d0)
                          elseif (mb > 0) then
                             angint = angint * io1 * (-1.0)**mb*sign(1,Lx)/sqrt(2.0d0)
                          endif
                       endif
                       !write(6,'(A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,2F14.9)') 'lb=',Lb,'mb=',mb,'l2=',l2,'l1=',l1,'m2=',m2,'ang=',angint
                       !* finished angular part
                       do ir2 = 1, nr2       !* over all radial functions at this (l2,m2)
                          do ie2 = ibgw+1, nbgw     !* Over the gw bands, i.e., not necessary all bands participate in GW calculation
                             select case (ir2)  !* bands to radial functions transformation af = A(ie,r2l2m2)
                             case(1)
                                af = alfa(ie2,l2m2,idf)
                             case(2)
                                af = beta(ie2,l2m2,idf)
                             case default
                                af = gama(ie2,ir2-2,l2m2,idf)
                             end select
                             !write(6,'(A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,2F14.9,1x,A,2F14.9,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x)') 'irk=',irk,'l2m2=',l2m2,'ie=',ie,'ir2=',ir2,'af=',af, 'ang=', angint, 'lb=',Lb,'l2=',l2,'l1=',l1,'mb=',mb,'m2=',m2
                             do ir1 = 1, nr1    !* over all radial functions at this (l1,m1)
                                do ie1 = ibgw+1, nbgw     !* Over the gw bands, i.e., not necessary all bands participate in GW calculation
                                   select case(ir1)!* bands to radial functions transformation caf = A(ie,r1l1m1)
                                   case(1)
                                      caf = conjg(alfa(ie1,l1m1,idf))
                                   case(2)
                                      caf = conjg(beta(ie1,l1m1,idf))
                                   case default
                                      caf = conjg(gama(ie1,ir1-2,l1m1,idf))
                                   end select
                                   ir12 = (ir2-1)*nr1 + ir1
                                   ! <u_{l1,m1}| V_{lxc} | u_{l2,m2}>
                                   uxcuw = uxcu(iat,lxc,l12,ir12)
                                   Vxcmt(ie1,ie2) = Vxcmt(ie1,ie2) + caf * uxcuw * af * angint
                                   !write(6,'(A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,2F14.9,1x,A,2F14.9,1x,A,F14.9)') 'irk=',irk,'iat=',iat,'ieq=',ieq,'l2=',l2,'l1=',l1,'m1=',m1,'lxc=',lxc,'isgn=',isgn,'ie=',ie,'ir2=',ir2,'ir1=',ir1,'vxcn=',vxcmt(ie), 'caf=', caf, 'uxcu=', uxcuw
                                enddo
                             enddo
                          enddo
                       enddo ! ie
                       !
                       if(mb.eq.0) exit   ! we do now want to have two copies of 0, because we have and 0*sign == [0,-0]
                    enddo ! isign
                 enddo ! lxc
              enddo ! m1
              !do ie=ibgw+1,nbgw
              !   write(6,'(A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,2F14.9)') 'irk=',irk,'l2=',l2,'l1=',l1,'ie=',ie,'vxcn=',Vxcmt(ie)
              !enddo
           enddo ! l1
        enddo ! i2
     enddo ! ieq
  enddo ! iat
end subroutine mt_vxcnn

!subroutine interst_vxcnn(Vxci, Vxsc, Aeigc, jipw, istpw, nbmax, nvk, nksxc, npw2apw)
!  implicit none
!  complex*16,intent(out) :: Vxci(nbmax,nbmax)   ! exchange-correlation potential
!  complex*16, intent(in) :: Aeigc(nbmax,nvk)    ! this is the eigenvector withouth local orbitals.
!  integer,    intent(in) :: jipw(nvk,nvk)       ! jipw is index for difference between two G points from vector file
!  complex*16, intent(in) :: Vxsc(nksxc)
!  complex*16, intent(in) :: istpw(npw2apw,nksxc)
!  integer,    intent(in) :: nbmax, nvk, nksxc, npw2apw
!  !
!  complex*16 :: inti(nvk,nvk), tmat(nbmax,nvk)
!  complex*16 :: alpha, beta
!  integer    :: ikxc
!  !
!  Vxci(:,:) = 0.d0
!  do ikxc=1,nksxc
!     call get_inti(inti, jipw, istpw(:,ikxc), nvk,  npw2apw)
!     ! Vxc_{ij} = \sum_{G1,G2,K} (Aeig_{i,G1})* V_{xc,K}<G1| e^{iK\vr}|G2> Aeig_{j,G2}
!     alpha = Vxsc(ikxc)
!     beta = 0.d0
!     print *, 'nvk=', nvk, 'nbmax=', nbmax
!     call zgemm('c','n', nbmax,   nvk, nvk, alpha, Aeigc, nbmax, inti,  nvk,   beta, tmat, nbmax )
!     alpha = 1.d0
!     beta  = 1.d0
!     print *, 'nvk=', nvk, 'nbmax=', nbmax
!     call zgemm('n','t', nbmax, nbmax, nvk, alpha, tmat,  nbmax, Aeigc, nbmax, beta,  Vxci, nbmax )
!  enddo
!end subroutine interst_vxcnn
                    
subroutine intstipw(istpw, ksxc, gindex, timat, tau, pos, vmt, Rmt, mult, k2cartes, npw2apw,  nksxc, nsym, ng, nat, ndf)
  ! This routine computes integral over interstitials for difference of two vectors iK and iG, where
  !   iK is one of the vectors listed in Vxc, and iG is one reciprocal vector up to some cutoff (used in vector file).
  !  We compute
  !                 Integrate[ (iK - iG)\vr, {\vr in interstitial}]
  ! However, we average over all members of the star corresponding to iK, i.e.,
  !          
  !               \sum_{iK members of the star} Integrate[ (iK - iG)\vr, {\vr in interstitial}] /{#group_operations}
  implicit none
  interface
     COMPLEX*16 Function int1ipw(iG, k2cartes, vmt, Rmt, pos, mult, nat, ndf)
       integer, intent(in) :: iG(3)
       real(8), intent(in) :: k2cartes(3,3)
       REAL*8,  intent(in) :: vmt(nat)
       REAL*8,  intent(in) :: Rmt(nat)
       REAL*8,  intent(in) :: pos(3,ndf)
       INTEGER, intent(in) :: mult(nat)
       INTEGER, intent(in) :: nat, ndf
     end Function int1ipw
  end interface
  complex*16,intent(out) :: istpw(npw2apw,nksxc)
  integer,    intent(in) :: ksxc(nksxc,3)
  integer,    intent(in) :: timat(nsym,3,3)
  real*8,     intent(in) :: tau(nsym,3), k2cartes(3,3)
  integer,    intent(in) :: gindex(ng,3)          ! all G-points 
  integer,    intent(in) :: npw2apw, nksxc, nsym, ng
  real*8,     intent(in) :: vmt(nat), Rmt(nat), pos(3,ndf)
  integer,    intent(in) :: mult(nat), nat, ndf
  !
  real*8, parameter :: two_pi = 6.28318530717959d0
  integer    :: ik, i, isym, iKs(3), iG(3), iG_i(3)!, iK_sym(3)
  !real*8     :: ph
  complex*16 :: int1!, phim
  integer    :: iKsym(3,nsym)
  real*8     :: phs(nsym)
  complex*16 :: imag, csum, ph_im(nsym)
  imag = cmplx( 0.d0, 1.d0 , kind=kind(1.0d0))
  
  do ik=1,nksxc                                           ! over all G vectors in Vxc
     iKs = ksxc(ik,:)                                     ! G vector in Vxc
     !
     do isym=1,nsym
        iKsym(:,isym) = matmul(iKs, timat(isym,:,:) )     ! all members of the star
     end do
     phs(:) = matmul(iKs, transpose(tau) )                ! phase due to translation for all members of the star
     ph_im(:) = exp( phs(:)*imag*two_pi)                  ! the corresponding complex phase
     do i=1,npw2apw                                       ! over all G vectors
        iG_i = gindex(i,:)
        csum = 0.d0
        do isym=1,nsym                                    ! all members of the star
           iG = iKsym(:,isym) - iG_i(:)                   ! iG_new = g_symmetry K_{xc} - G
           int1 = int1ipw(iG, k2cartes, vmt, Rmt, pos, mult, nat, ndf) ! integral over interstitials for this iG_new
           csum = csum + int1 * ph_im(isym)               ! now also ading the phase
        enddo
        istpw(i,ik) = csum/nsym
     end do
  end do
  !
!     do isym=1,nsym                                      ! over all members of the star
!        iK_sym = matmul(iKs, timat(isym,:,:) )           ! member of the star
!        ph     = dot_product(iKs,   tau(isym,:) )*two_pi ! phase for this member of the star : K_{xc}*g_symmetry_tau
!        phim   = cmplx( cos(ph), -sin(ph) , kind=kind(1.0d0)) ! exp(-K_{xc}*g_symmetry_tau)
!        do i=1,npw2apw                                   ! over all G vectors
!           iG = iK_sym - gindex(i,:)                     ! iG_new = g_symmetry K_{xc}-G
!           int1 = int1ipw(iG, k2cartes, vmt, Rmt, pos, mult, nat, ndf) ! integral over interstitials for this iG_new
!           istpw(i,ik) = istpw(i,ik) + int1 * phim       ! now also ading the phase
!        end do
!     end do
!  end do
!  istpw(:,:) = istpw(:,:)/nsym
end subroutine intstipw

COMPLEX*16 Function int1ipw(iG, k2cartes, vmt, Rmt, pos, mult, nat, ndf)
  ! This function calculated integral of the plane wave in the interstitials
  !
  !  Integrate[ exp(i G\vr), {\vr in interstitials}]
  !
  IMPLICIT NONE
  integer, intent(in) :: iG(3)
  real(8), intent(in) :: k2cartes(3,3)
  REAL*8,  intent(in) :: vmt(nat)
  REAL*8,  intent(in) :: Rmt(nat)
  REAL*8,  intent(in) :: pos(3,ndf)
  INTEGER, intent(in) :: mult(nat)
  INTEGER, intent(in) :: nat, ndf
  !
  real*8, parameter :: two_pi = 6.28318530717959d0
  INTEGER    :: idf, iat, ieq!, ki(3)
  REAL*8     :: ak, x, j1, intmod, pk, Gc(3)
  COMPLEX*16 :: integc
  COMPLEX*16 :: phase
  !
  Gc = matmul(k2cartes, iG)
  ak = sqrt(dot_product(Gc,Gc))
  !print *, 'ak=', ak, 'iG=', iG, 'Gc=', Gc
  if (ak < 1e-10) then
     integc = 0.d0
     do iat=1,nat
        integc = integc + vmt(iat)*mult(iat)
     end do
     int1ipw = 1.d0 - integc
  else
     integc = 0.d0
     !
     !ph = exp( matmul(iG,pos) * imag )
     !idf = 1
     !do iat=1,nat
     !   phase = sum(ph[idf:idf+mult(iat)-1])
     !   x = ak*Rmt(iat)
     !   j1 = (sin(x)/x-cos(x))/x             ! spherical_bessel == (sin(x)/x-cos(x))/x
     !   integc = integc + phase * vmt(iat) * 3*j1/x
     !   idf = idf + mult(iat)
     !end do
     !
     idf = 0
     do iat=1,nat
        phase = 0.d0
        do ieq=1,mult(iat)
           idf = idf + 1
           pk = dot_product(iG, pos(:,idf))*two_pi
           !pk = (iG(1)*pos(1,idf) + iG(2)*pos(2,idf) + iG(3)*pos(3,idf))*two_pi
           phase = phase + cmplx( cos(pk), sin(pk) , kind=kind(1.0d0))
           !print *, 'iat=', iat, 'pk=', pk, 'phase=', phase
        end do
        x = ak * Rmt(iat)                    ! kr[iat] = |kc|*Rmt[iat]
        j1 = (sin(x)/x-cos(x))/x             ! spherical_bessel == (sin(x)/x-cos(x))/x
        intmod = vmt(iat) * 3*j1/x           ! intmod[iat] = j1*3*vmt/kr
        integc = integc + intmod * phase
        !print *, 'x=', x, 'j1=', j1, 'intmod=', intmod, 'phase=', phase, 'integc=', integc
     end do
     !print *, 'Finally=', integc
     int1ipw = -integc
  end if
  return
end Function int1ipw

subroutine two_gs_to_one(jipw, gindex, indgk, i_g0, nvk, ng, n1, n2, n3)
  ! This subroutine calculates index array jimp[iv1,iv2] which finds index of two reciprocal vectors G1-G2
  ! in fixed basis (pw.gindex), where  G1 and G2 are reciprocal vectors from Hamiltonian
  ! corresponding to a particular irreducible k-point.
  implicit none
  integer, intent(out) :: jipw(nvk,nvk)  ! jipw is index for difference between two G points from vector file
  integer, intent(in)  :: gindex(ng,3)          ! all G-points 
  integer, intent(in)  :: indgk(nvk)     ! be careful. In python indgk has larger dimension. Just take out the slice of it.
  integer, intent(in)  :: i_g0(2*n1+1,2*n2+1,2*n3+1)
  integer, intent(in)  :: nvk, ng, n1, n2, n3
  !
  integer :: iGc1(3), iGc2(3), iG(3)
  integer :: iv1, iv2!, ii
  !
  !print *, 'n1=', n1, 'n2=', n2, 'n3=', n3
  jipw(:,:) = -1
  do iv1 = 1,nvk
     do iv2 = 1,nvk
        iGc1 = gindex(indgk(iv1)+1,:)
        iGc2 = gindex(indgk(iv2)+1,:)
        iG = iGc1 - iGc2
        jipw(iv1,iv2) = i_g0(iG(1)+n1+1,iG(2)+n2+1,iG(3)+n3+1)
        !write(6,'(I4,1x,I4,3x,3I3,3x,3I3,3x,I3)') iv1, iv2, iGc1, iGc2, jipw(iv1,iv2)
     enddo
  enddo
end subroutine two_gs_to_one

subroutine get_inti(inti, jipw, istpw, nvk, npw2apw)
  ! Now we calculate interstitail integral for the difference G1-G2 (with G1 and G2 from Hamiltonian basis)
  !       inti[1,2] = Integrate[ (iK - i(G1-G2))\vr, {\vr in interstitial}]  = <G1|e^{iK\vr}|G2>
  implicit none
  complex*16, intent(out) :: inti(nvk,nvk)        ! output
  integer,    intent(in)  :: jipw(nvk,nvk)        ! jipw is index for difference between two G points from vector file
  complex*16, intent(in)  :: istpw(npw2apw)       ! interstitial integral Integrate[ (iK - iG)\vr, {\vr in interstitial}] where K is averaged over all members of the star
  integer,    intent(in)  :: nvk, npw2apw
  !
  integer :: iv1, iv2
  !
  inti = 0.d0
  do iv1=1,nvk
     do iv2=1,nvk
        if ( jipw(iv1,iv2) >= 0 ) then
           inti(iv1,iv2) = istpw(jipw(iv1,iv2)+1)
        end if
     end do
  end do
end subroutine get_inti

                    
!subroutine inter_vxcnn(gindex, indgk, nvk)
!  integer, intent(in) :: gindex(ng,3)          ! all G-points 
!  integer, intent(in) :: indgk(ngk)
!  integer, intent(in) :: nvk, ng
!  !     --------------------
!  !     Interstitial region:
!  !     --------------------
!  nvk=nv(irk) 
!  allocate(jipw(nvk,nvk),inti(nvk,nvk),tmat1(nvk,ibgw:nbgw),tvec1(nvk),tvec2(nvk))
!  do iv1=1, nvk
!     if(indgkir(iv1,irk).gt.npw) then
!        write(msg,*) "indgkir(iv1,irk) > npw","iv1,irk,indgkir=",iv1,irk,indgkir(iv1,irk)
!        call outerr(sname,msg)
!     endif
!     do iv2=1, nvk
!        if(indgkir(iv2,irk).gt.npw) then
!           write(msg,*)  "indgkir(iv2,irk) > npw","  iv2,irk,indgkir=",iv2,irk,indgkir(iv2,irk)
!           call outerr(sname,msg)
!        endif
!        do i=1,3
!           ik(i) = gindex(i,indgkir(iv1,irk))-gindex(i,indgkir(iv2,irk))
!        enddo
!        if((minval(ik).ge.ngmin).and.(maxval(ik).le.ngmax)) then
!           jipw(iv1,iv2)=ig0(ik(1),ik(2),ik(3))
!        else
!           jipw=0
!           write(6,'(a,i6,2x,2i4,2x,3i4)') 'WARNING: jipw=0 in w2k_calcvxcnn',jipw(iv1,iv2),ngmin,ngmax,ik
!        endif
!     enddo ! iv2
!  enddo ! iv1
!  
!  !
!  !     Loop over the mixed basis functions:
!  !
!  do ikxc = 1, nksxc
!     do iv1=1, nvk   !* loop over G 
!        do iv2=1, nvk   !* loop over G'
!           if(jipw(iv1,iv2).gt.0)then
!              inti(iv1,iv2)=istpw(ikxc,jipw(iv1,iv2))
!           else  
!              inti(iv1,iv2)=czero
!           endif
!        enddo ! iv2
!     enddo ! iv1
!     
!     call cpu_time(time1) 
!     call zgemm('n','n',nvk,nbgw-ibgw+1,nvk,vxcs(ikxc,isp),inti,nvk,zzk(:,ibgw:nbgw),maxngk,czero,tmat1,nvk)
!     call cpu_time(time2) 
!     time_lapack=time_lapack+time2-time1
!
!     !print *, 'shape(tvec1)=', shape(tvec1)
!     !print *, 'shape(tvec2)=', shape(tvec2)
!     !print *, 'nvk=', nvk, 'ie=', ibgw, nbgw, 'irk=', irk
!     !print *, 'shape(sumiv2)=', shape(sumiv2)
!     !print *, 'shape(vxci)=', shape(vxci)
!     !print *, 'shape(zzk)=', shape(zzk)
!     !print *, 'shape(tmat)=', shape(tmat1)
!     
!     do ie=ibgw,nbgw
!        tvec1(1:nvk)=zzk(1:nvk,ie)
!        tvec2(1:nvk)=tmat1(1:nvk,ie)
!        !KH           !sumiv2(ie)=zdotc(nvk,tvec1,1,tvec2,1)
!        !sumiv2(ie) = dot_product(tvec1, tvec2)
!        sumiv2(ie) = mzdotc(nvk, tvec1, tvec2)
!        vxci(ie,irk)=vxci(ie,irk)+sumiv2(ie)
!     enddo ! ie
!  enddo ! ikxc
!  deallocate(jipw,inti,tmat1,tvec1,tvec2) 
!end subroutine inter_vxcnn

REAL*8 function getcgcoef(l1,l2,l3,m1,m2,cgcoef)
!!! <Y_{l3,m1+m2}|Y_{l1,m1}Y_{l2,m2}>
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
  if (trcond) then
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


subroutine calc_minm(mmat, ql, nstart,nend,mstart,mend, alfa,beta,gama,alfp,betp,gamp, s3r, pos, mult, nmix, big_l, nLO_at, ncore, cgcoef, lmax, mbsiz, nbmax, mbmax, nat, ntnt, nLOmax, ndf, n_mix, nhow_many_fnc, max_nmix, lomax, nbmax_g,mbmax_g,ntnt_g,ndf_g)
  !This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$ of equation \ref{Mdef} for n, m both belongings to LAPW states 
  !ql(i)=dble(qlist(i,iq))/dble(idvq)
  implicit none
  complex(8), intent(out):: mmat(nstart:nend,mbsiz,mstart:mend)
  real(8),    intent(in) :: ql(3)
  integer,    intent(in) :: nstart,nend,mstart,mend  !! ranges of bands n and m, respectively
  complex*16, intent(in) :: alfa(nbmax,ntnt,ndf)
  complex*16, intent(in) :: beta(nbmax,ntnt,ndf)
  complex*16, intent(in) :: gama(nbmax_g,nLOmax,ntnt_g,ndf_g)
  complex*16, intent(in) :: alfp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: betp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: gamp(mbmax_g,nLOmax,ntnt_g,ndf_g)
  real(8),    intent(in) :: s3r(n_mix,nhow_many_fnc,nhow_many_fnc,nat)   ! s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} > 
  real(8),    intent(in) :: pos(3,ndf)
  integer,    intent(in) :: mult(nat)
  integer,    intent(in) :: nmix(nat)          ! how many product functions per atom
  integer,    intent(in) :: big_l(max_nmix,nat)! L for product basis
  integer,    intent(in) :: nLO_at(lomax+1,nat), ncore(nat)
  real*8,     intent(in) :: cgcoef(*)          ! gaunt coefficients
  integer,    intent(in) :: lmax, mbsiz, nbmax, mbmax, ntnt, nLOmax, ndf, n_mix, nhow_many_fnc, nat, max_nmix, lomax
  integer,    intent(in) :: nbmax_g,mbmax_g,ntnt_g,ndf_g
  !
  interface
     REAL*8 function getcgcoef(l1,l2,l3,m1,m2,cgcoef)
       REAL*8, intent(in) :: cgcoef(*)
       integer(4), intent(in) :: l1, l2, l3, m1, m2
     end function getcgcoef
  end interface
  !
  real(8),    parameter  :: pi = 3.14159265358979d0
  logical :: PRINT
  integer :: iat, lb, mb, idf, ie1, ie2, ieq
  complex*16  :: na,nb,nc(nLOmax), angint, phs, sm, two_pi_i, phs_angint
  integer :: lc10, lc20, ilo1, ilo2, l2min, l2max, la1, lb1, lc1, la2, lb2, lc2, im, irm
  integer :: l1, m1, lm1, l2, m2, lm2
  !
  PRINT = .False.
  two_pi_i = 2.0d0 * pi * cmplx(0.d0,1.d0, 8)
  mmat(:,:,:) = cmplx(0.d0, 0.d0, 8)
  !
  if (PRINT) write(6,*) 'calcminm'
  nc = 0
  idf = 0
  im  = 0
  ! put together : alfas(ie1,lmall) == [alfa(ie1,lm1,idf), beta(ie1,lm1,idf),gama(ie1,ilo1,lm1,idf)]
  ! put together : alfps(ie2,lmall) == [alfp(ie1,lm1,idf), betp(ie1,lm1,idf),gamp(ie1,ilo1,lm1,idf)]
  ! output       : alfps * MM * alfas.H
  ! where        : MM(lmall,lmall)
  do iat = 1, nat               !! Loop over atoms
     do ieq = 1, mult(iat)      !! Loop over equivalent atoms
        idf = idf + 1
        phs = exp( -dot_product( ql, pos(:,idf) )*two_pi_i )
        do irm = 1, nmix(iat)   !!  Loop over mixed functions:
           lb = big_l(irm,iat)
           do mb = -lb,lb
              im = im + 1       !! index of the product basis functions
              do l1 = 0, lmax              ! index of the single-particle function u_{lm}
                 l2min = abs(lb-l1)        ! l2 is bounded by triangular rule of <Y_{l1m1}|Y_{l2,m2}|Y_{lb,mb}>
                 l2max = min(lb+l1,lmax)
                 la1 = ncore(iat) + l1 + 1 ! this is index for the matrix element s3r previously computed s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} >
                 lb1 = la1 + lmax + 1      ! why extra 1 I do not understand!!!
                 lc10 = ncore(iat) + 2*(lmax+1) + sum(nLO_at(:l1,iat))  ! where local orbitals start
                 do l2 = l2min,l2max           !! loop over l2:
                    la2 = l2 + ncore(iat) + 1  
                    lb2 = la2 + lmax + 1
                    lc20 = ncore(iat) + 2*(lmax+1) + sum(nLO_at(:l2,iat))  ! where local orbitals start
                    do m1 = -l1,l1             
                       lm1 = l1*l1 + l1 + m1 + 1   ! index for Y_{l1,m1}
                       m2 = m1 - mb                ! This comes from matrix element <Y_{l1m1}|Y_{l2,m2}|Y_{lb,mb}>, which requires m1==mb+m2
                       if ( abs(m2) .gt. l2) cycle 
                       lm2 = l2*l2 + l2 + m2 + 1   ! index for Y_{l2,m2}
                       angint = getcgcoef(l2,lb,l1,m2,mb, cgcoef)  ! gaunt coefficients = <Y_{l1m1}|Y_{l2,m2}|Y_{lb,mb}>* = <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>
                       if ( abs(angint) .lt. 1.0d-8) cycle
                       phs_angint = phs*angint
                       !write(6,'(I4,1x,4I3,2F16.10,2x,2I3)') im, l1, l2, m1, m2, phs_angint, lc10, lc20
                       do ie1 = nstart+1,nend                        !! Loop over eigenfunctions at k
                          ! s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} > 
                          na = s3r(irm,la1,la2,iat) * alfa(ie1,lm1,idf) + s3r(irm,lb1,la2,iat) * beta(ie1,lm1,idf) ! na=alfa(ie1,l1m1)<u_{a1,l1}|u_{a2,l2}|u_{irm,lb}>   +beta(ie1,l1m1)<udot_{a1,l1}|u_{a2,l2}|u_{irm,lb}>
                          nb = s3r(irm,la1,lb2,iat) * alfa(ie1,lm1,idf) + s3r(irm,lb1,lb2,iat) * beta(ie1,lm1,idf) ! nb=alfa(ie1,l1m1)<u_{a1,l1}|udot_{a2,l2}|u_{irm,lb}>+beta(ie1,l1m1)<udot_{a1,l1}|udot_{a2,l2}|u_{irm,lb}>
                          do ilo1 = 1,nLO_at(l1+1,iat) 
                             lc1 = lc10 + ilo1
                             na = na + s3r(irm,lc1,la2,iat) * gama(ie1,ilo1,lm1,idf) ! na += <u_{irm,lb}| u_{lo,a1,l1}    u_{a2,l2} > * gama(ie1,ilo,l1m1)
                             nb = nb + s3r(irm,lc1,lb2,iat) * gama(ie1,ilo1,lm1,idf) ! nb += <u_{irm,lb}| u_{lo,a1,l1} udot_{a2,l2}|> * gama(ie1,ilo,l1m1)
                          enddo
                          !write(6,'(I3,1x,A,2F16.10,1x,A,2F16.10)') ie1, 'na=',  na, 'nb=', nb
                          ! na = < u_{irm,lb} |    u_{a2,l2} psi_{ie1} > = < u_{irm,lb} | (alfa u_{a1,l1} + beta udot_{a1,l1} + gama ulo_{a1,l1}) *    u_{a2,l2} >
                          ! nb = < u_{irm,lb} | udot_{a2,l2} psi_{ie1} > = < u_{irm,lb} | (alfa u_{a1,l1} + beta udot_{a1,l1} + gama ulo_{a1,l1}) * udot_{a2,l2} >
                          ! nc = < u_{irm,lb} |  ulo_{a2,l2} psi_{ie1} > = < u_{irm,lb} | (alfa u_{a1,l1} + beta udot_{a1,l1} + gama ulo_{a1,l1}) *  ulo_{a2,l2} >
                          do ilo2 = 1,nLO_at(l2+1,iat)
                             lc2 = lc20 + ilo2
                             nc(ilo2) = s3r(irm,la1,lc2,iat) * alfa(ie1,lm1,idf) + s3r(irm,lb1,lc2,iat) * beta(ie1,lm1,idf)
                             do ilo1 = 1,nLO_at(l1+1,iat) 
                                lc1 = lc10 + ilo1
                                nc(ilo2) = nc(ilo2) + s3r(irm,lc1,lc2,iat) * gama(ie1,ilo1,lm1,idf)
                             enddo
                          enddo
                          ! sm = < u_{irm,lb} | psi^*_{ie2} psi_{ie1} > = < u_{irm,lb} | (alfp * u_{a2,l2} + betp * udot_{a2,l2} + gamp * ulo_{a2,l2}) psi_{ie1} >
                          do ie2 = mstart+1,mend
                             sm = alfp(ie2,lm2,idf) * na + betp(ie2,lm2,idf) * nb
                             do ilo2 = 1,nLO_at(l2+1,iat)
                                sm = sm + gamp(ie2,ilo2,lm2,idf) * nc(ilo2)
                             enddo
                             mmat(ie1,im,ie2) = mmat(ie1,im,ie2) + phs_angint * sm
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              if (PRINT) then
                 do ie1=nstart+1,nend
                    do ie2 = mstart+1,mend
                       write(6,'(A,3I4,1x,2F16.10)') 'mm=', im,ie2,ie1, mmat(ie1,im,ie2)
                    end do
                 end do
              endif
           enddo
        enddo
     enddo
  enddo
end subroutine calc_minm


subroutine calc_minm2(mmat, ql, nstart,nend,mstart,mend, alfa,beta,gama,alfp,betp,gamp, s3r, pos, mult, nmix, big_l, nLO_at, ncore, cgcoef, lmax, mbsiz, nbmax, mbmax, nat, ntnt, nLOmax, ndf, n_mix, nhow_many_fnc, max_nmix, lomax, nbmax_g,mbmax_g,ntnt_g,ndf_g)
  !This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$ of equation \ref{Mdef} for n, m both belongings to LAPW states 
  !ql(i)=dble(qlist(i,iq))/dble(idvq)
  implicit none
  complex(8), intent(out):: mmat(nend-nstart,mbsiz,mend-mstart) 
  real(8),    intent(in) :: ql(3)
  integer,    intent(in) :: nstart,nend,mstart,mend  !! ranges of bands n and m, respectively
  complex*16, intent(in) :: alfa(nbmax,ntnt,ndf)
  complex*16, intent(in) :: beta(nbmax,ntnt,ndf)
  complex*16, intent(in) :: gama(nbmax_g,nLOmax,ntnt_g,ndf_g)
  complex*16, intent(in) :: alfp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: betp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: gamp(mbmax_g,nLOmax,ntnt_g,ndf_g)
  real(8),    intent(in) :: s3r(n_mix,nhow_many_fnc,nhow_many_fnc,nat)   ! s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} > 
  real(8),    intent(in) :: pos(3,ndf)
  integer,    intent(in) :: mult(nat)
  integer,    intent(in) :: nmix(nat)          ! how many product functions per atom
  integer,    intent(in) :: big_l(max_nmix,nat)! L for product basis
  integer,    intent(in) :: nLO_at(lomax+1,nat), ncore(nat)
  real*8,     intent(in) :: cgcoef(*)          ! gaunt coefficients
  integer,    intent(in) :: lmax, mbsiz, nbmax, mbmax, ntnt, nLOmax, ndf, n_mix, nhow_many_fnc, nat, max_nmix, lomax
  integer,    intent(in) :: nbmax_g,mbmax_g,ntnt_g,ndf_g
  !
  interface
     REAL*8 function getcgcoef(l1,l2,l3,m1,m2,cgcoef)
       REAL*8, intent(in) :: cgcoef(*)
       integer(4), intent(in) :: l1, l2, l3, m1, m2
     end function getcgcoef
  end interface
  !
  real(8), parameter  :: pi = 3.14159265358979d0
  logical   :: PRINT
  integer   :: iat, lb, mb, idf, ie1, ie2, ieq
  complex*16:: angint, phs, two_pi_i, phs_angint, ccsm
  integer   :: lc10, lc20, ilo1, ilo2, l2min, l2max, la1, lb1, lc1, la2, lb2, lc2, im, irm, im_st, im1, im1_size
  integer   :: l1, m1, lm1, l2, m2, lm2
  complex*16:: aa, ab, ba, bb, ac(nLOmax), ca(nLOmax), bc(nLOmax), cb(nLOmax), cc(nLOmax,nLOmax), cone, czero
  integer   :: ncl2, ic2, ie1p, cind(nLOmax,(lmax+1)**2)
  complex*16, allocatable :: na(:,:,:), nb(:,:,:), nc(:,:,:), mmt2(:,:,:), talfp(:,:), tbetp(:,:), tgamp(:,:)
  !
  allocate( talfp((lmax+1)**2,mend-mstart), tbetp((lmax+1)**2,mend-mstart) )
  !
  PRINT = .False.
  two_pi_i = 2.0d0 * pi * cmplx(0.d0,1.d0, 8)
  cone = cmplx(1.d0, 0.d0,8)
  czero = cmplx(0.d0,0.d0,8)
  mmat(:,:,:) = cmplx(0.d0, 0.d0, 8)
  !
  if (PRINT) write(6,*) 'calcminm', ' nLOmax=', nLOmax
  
  idf = 0
  im  = 0
  im_st=0
  do iat = 1, nat               !! Loop over atoms
     ncl2 = 0
     do l2=0,lmax
        do m2=-l2,l2
           lm2 = l2*l2 + l2 + m2 + 1
           do ilo2 = 1,nLO_at(l2+1,iat)
              ncl2 = ncl2 + 1
              cind(ilo2,lm2) = ncl2
              !print *, 'l=', l2, 'm=', m2, 'lm=', lm2, 'ilo=', ilo2, 'nl=', ncl2
           enddo
        enddo
     enddo
     if (ncl2 > 0) then
        allocate( tgamp(ncl2,mend-mstart) )
        tgamp = 0.d0
     endif
     
     do ieq = 1, mult(iat)      !! Loop over equivalent atoms
        idf = idf + 1
        !
        do ie2 = 1,mend-mstart
           do lm2=1,(lmax+1)**2
              talfp(lm2,ie2) = alfp(ie2+mstart,lm2,idf)
              tbetp(lm2,ie2) = betp(ie2+mstart,lm2,idf)
           enddo
        enddo
        if (ncl2 > 0 ) then
           do ie2 = 1,mend-mstart
              do l2=0,lmax
                 do m2=-l2,l2
                    lm2 = l2*l2 + l2 + m2 + 1
                    do ilo2 = 1,nLO_at(l2+1,iat)
                       tgamp(cind(ilo2,lm2),ie2) = gamp(ie2+mstart,ilo2,lm2,idf)
                    enddo
                 enddo
              enddo
           enddo
        endif
        !
        phs = exp( -dot_product( ql, pos(:,idf) )*two_pi_i )
        ! determine the size of this block of product-basis, which we want to calculate
        im1_size = 0
        do irm = 1, nmix(iat)   !!  Loop over mixed functions:
           lb = big_l(irm,iat)
           im1_size = im1_size + 2*lb+1
        enddo
        
        allocate( na(nend-nstart,im1_size,(lmax+1)**2), nb(nend-nstart,im1_size,(lmax+1)**2) )
        na = 0.d0
        nb = 0.d0
        if (ncl2 > 0) then
           allocate( nc(ncl2,nend-nstart,im1_size) )
           nc = 0.d0
        endif
        
        im1 = 0
        do irm = 1, nmix(iat)   !!  Loop over mixed functions:
           lb = big_l(irm,iat)
           do mb = -lb,lb
              im = im + 1       !! index of the product basis functions
              im1 = im1 + 1
              !print *, im_st + im1, im
              do l1 = 0, lmax              ! index of the single-particle function u_{lm}
                 l2min = abs(lb-l1)        ! l2 is bounded by triangular rule of <Y_{l1m1}|Y_{l2,m2}|Y_{lb,mb}>
                 l2max = min(lb+l1,lmax)
                 la1 = ncore(iat) + l1 + 1 ! this is index for the matrix element s3r previously computed s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} >
                 lb1 = la1 + lmax + 1      ! why extra 1 I do not understand!!!
                 lc10 = ncore(iat) + 2*(lmax+1) + sum(nLO_at(:l1,iat))  ! where local orbitals start
                 do l2 = l2min,l2max           !! loop over l2:
                    la2 = l2 + ncore(iat) + 1  
                    lb2 = la2 + lmax + 1
                    lc20 = ncore(iat) + 2*(lmax+1) + sum(nLO_at(:l2,iat))  ! where local orbitals start
                    do m1 = -l1,l1             
                       lm1 = l1*l1 + l1 + m1 + 1   ! index for Y_{l1,m1}
                       m2 = m1 - mb                ! This comes from matrix element <Y_{l1m1}|Y_{l2,m2}|Y_{lb,mb}>, which requires m1==mb+m2
                       if ( abs(m2) .gt. l2) cycle 
                       lm2 = l2*l2 + l2 + m2 + 1   ! index for Y_{l2,m2}
                       angint = getcgcoef(l2,lb,l1,m2,mb, cgcoef)  ! gaunt coefficients = <Y_{l1m1}|Y_{l2,m2}|Y_{lb,mb}>* = <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>
                       if ( abs(angint) .lt. 1.0d-8) cycle
                       phs_angint = phs*angint
                       
                       aa = phs_angint * s3r(irm,la1,la2,iat) ! aa = <u_{irm,lb}|    u_{a1,l1}    u_{a2,l2}> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       ab = phs_angint * s3r(irm,lb1,la2,iat) ! ab = <u_{irm,lb}| udot_{a1,l1}    u_{a2,l2}> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       ba = phs_angint * s3r(irm,la1,lb2,iat) ! ba = <u_{irm,lb}|    u_{a1,l1} udot_{a2,l2}> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       bb = phs_angint * s3r(irm,lb1,lb2,iat) ! bb = <u_{irm,lb}| udot_{a1,l1} udot_{a2,l2}> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       do ilo1 = 1,nLO_at(l1+1,iat) 
                          lc1 = lc10 + ilo1
                          ac(ilo1) = phs_angint * s3r(irm,lc1,la2,iat) ! ac = <u_{irm,lb}| u_{lo,a1,l1}    u_{a2,l2} > <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                          bc(ilo1) = phs_angint * s3r(irm,lc1,lb2,iat) ! bc = <u_{irm,lb}| u_{lo,a1,l1} udot_{a2,l2}|> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       enddo
                       do ilo2 = 1,nLO_at(l2+1,iat)
                          lc2 = lc20 + ilo2
                          ca(ilo2) = phs_angint * s3r(irm,la1,lc2,iat) ! ca = < u_{irm,lb} |    u_{a1,l1} ulo_{a2,l2} > <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                          cb(ilo2) = phs_angint * s3r(irm,lb1,lc2,iat) ! cb = < u_{irm,lb} | udot_{a1,l1} ulo_{a2,l2} > <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                          do ilo1 = 1,nLO_at(l1+1,iat)
                              lc1 = lc10 + ilo1
                             cc(ilo1,ilo2) = phs_angint * s3r(irm,lc1,lc2,iat) ! cc = < u_{irm,lb} | ulo_{a1,l1}) ulo_{a2,l2} > <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                          enddo
                       enddo
                       !write(6,'(I4,1x,4I3,2F16.10,2x,2I3)') im, l1, l2, m1, m2, phs_angint, lc10, lc20
                       ! na = < u_{irm,lb} |    u_{a2,l2} psi_{ie1} >e^{-iqr} = < u_{irm,lb} | (alfa u_{a1,l1} + beta udot_{a1,l1} + gama ulo_{a1,l1}) *    u_{a2,l2} >e^{-iqr}
                       ! nb = < u_{irm,lb} | udot_{a2,l2} psi_{ie1} >e^{-iqr} = < u_{irm,lb} | (alfa u_{a1,l1} + beta udot_{a1,l1} + gama ulo_{a1,l1}) * udot_{a2,l2} >e^{-iqr}
                       ! nc = < u_{irm,lb} |  ulo_{a2,l2} psi_{ie1} >e^{-iqr} = < u_{irm,lb} | (alfa u_{a1,l1} + beta udot_{a1,l1} + gama ulo_{a1,l1}) *  ulo_{a2,l2} >e^{-iqr}
                       do ie1=1,nend-nstart                        !! Loop over eigenfunctions at k
                          ie1p = ie1+nstart
                          na(ie1,im1,lm2) = na(ie1,im1,lm2) + aa * alfa(ie1p,lm1,idf) + ab * beta(ie1p,lm1,idf) 
                          nb(ie1,im1,lm2) = nb(ie1,im1,lm2) + ba * alfa(ie1p,lm1,idf) + bb * beta(ie1p,lm1,idf) 
                          do ilo1 = 1,nLO_at(l1+1,iat)
                             na(ie1,im1,lm2) = na(ie1,im1,lm2) + ac(ilo1) * gama(ie1p,ilo1,lm1,idf) ! na += <u_{irm,lb}| u_{lo,a1,l1}    u_{a2,l2} > * gama(ie1,ilo,l1m1)<Y|Y|Y>e^{-iqr}
                             nb(ie1,im1,lm2) = nb(ie1,im1,lm2) + bc(ilo1) * gama(ie1p,ilo1,lm1,idf) ! nb += <u_{irm,lb}| u_{lo,a1,l1} udot_{a2,l2}|> * gama(ie1,ilo,l1m1)<Y|Y|Y>e^{-iqr}
                          enddo
                          do ilo2 = 1,nLO_at(l2+1,iat)
                             ccsm = 0.d0
                             do ilo1 = 1,nLO_at(l1+1,iat) 
                                ccsm = ccsm + cc(ilo1,ilo2) * gama(ie1p,ilo1,lm1,idf)
                             enddo
                             ic2 = cind(ilo2,lm2)
                             nc(ic2,ie1,im1) = nc(ic2,ie1,im1) + ccsm + ca(ilo2) * alfa(ie1p,lm1,idf) + cb(ilo2) * beta(ie1p,lm1,idf)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              
              !! sm = < u_{irm,lb} | psi^*_{ie2} psi_{ie1} ><Y|Y|Y>e^{-iqr}
              !do ie1= 1,nend-nstart
              !   do ie2 = 1,mend-mstart
              !      sm = 0.d0
              !      do l2=0,lmax
              !         do m2=-l2,l2
              !            lm2 = l2*l2 + l2 + m2 + 1
              !            !do lm2=1,(lmax+1)**2
              !            sm = sm + alfp(ie2+mstart,lm2,idf) * na(ie1,lm2)
              !            sm = sm + betp(ie2+mstart,lm2,idf) * nb(ie1,lm2)
              !            do ilo2 = 1,nLO_at(l2+1,iat)
              !               sm = sm + gamp(ie2+mstart,ilo2,lm2,idf) * nc(ilo2,lm2,ie1)
              !            enddo
              !         enddo
              !      enddo
              !      mmat(ie1,ie2,im) = sm
              !   enddo
              !enddo

              !mmt2 = matmul( na, talfp)             ! na(ie1,lm2)*talfp(lm2,ie2)
              !mmat(:,:,im) = mmt2(:,:)
              !mmt2 = matmul( nb, tbetp)             ! nb(ie1,lm2)*talfp(lm2,ie2)
              !mmat(:,:,im) = mmat(:,:,im) + mmt2(:,:)
              !if (ncl2 > 0) then
              !   mmt2 = matmul( transpose(nc), tgamp)               ! nc(ic2,ie1)*tgamp(ic2,ie2)
              !   mmat(:,:,im) = mmat(:,:,im) + mmt2(:,:)
              !endif
              
           enddo
        enddo

        allocate( mmt2(nend-nstart,im1_size,mend-mstart) )
        ! mmat(ie1,im,ie2) = na(ie1,im,lm)*talfp(lm,ie2) + nb(ie1,im,lm)*tbetp(lm,ie2) + nc(ic,ie1,im)*tgamp(ic,ie2)
        call zgemm('n','n', (nend-nstart)*im1_size, mend-mstart, (lmax+1)**2, cone, na,(nend-nstart)*im1_size, talfp,(lmax+1)**2, czero, mmt2,(nend-nstart)*im1_size)
        call zgemm('n','n', (nend-nstart)*im1_size, mend-mstart, (lmax+1)**2, cone, nb,(nend-nstart)*im1_size, tbetp,(lmax+1)**2, cone,  mmt2,(nend-nstart)*im1_size)
        if (ncl2 > 0) then
           call zgemm('t','n', (nend-nstart)*im1_size, mend-mstart, ncl2, cone, nc,ncl2, tgamp,ncl2, cone,  mmt2,(nend-nstart)*im1_size)
        endif
        mmat(:,(im_st+1):(im_st+im1_size),:) = mmt2(:,:,:)
        deallocate( mmt2 )
        
        if (PRINT) then
           do im1=1,im1_size
              im = im_st + im1
              do ie1=1,nend-nstart
                 do ie2 = 1,mend-mstart
                    write(6,'(A,3I4,1x,2F16.10)') 'mm=', im,ie2,ie1, mmat(ie1,im,ie2)
                 end do
              end do
           end do
        endif
        
        im_st = im_st + im1_size
        
        deallocate( na, nb )
        if (ncl2 > 0) deallocate( nc )
     enddo
     if (ncl2 > 0) deallocate( tgamp )
  enddo
  deallocate( talfp, tbetp )
end subroutine calc_minm2


subroutine calc_minm_MT(mmat, ql, nstart,nend,mstart,mend, alfa,beta,gama,alfp,betp,gamp, s3r, pos, mult, nmix, big_l, nLO_at, ncore, cgcoef, lmax, loctmatsize, nbmax, mbmax, nat, ntnt, nLOmax, ndf, n_mix, nhow_many_fnc, max_nmix, lomax, nbmax_g,mbmax_g,ntnt_g,ndf_g)
  !This subroutine calculates the matrix elements $M^i_{nm}(\vec{k},\vec{q})$
  !                               mmat(ie1,ie2,im) = < u^{product}_{im,lb} | psi^*_{ie2,k-q} psi_{ie1,k} > e^{-iqr}
  !                                                                psi^*_{ie2,k-q} e^{i(k-q)*r} psi_{ie1,k} e^{-ik*r} 
  implicit none
  complex*16, intent(out):: mmat(nend-nstart,mend-mstart,loctmatsize)
  real(8),    intent(in) :: ql(3)
  integer,    intent(in) :: nstart,nend,mstart,mend  !! ranges of bands n and m, respectively
  complex*16, intent(in) :: alfa(nbmax,ntnt,ndf)
  complex*16, intent(in) :: beta(nbmax,ntnt,ndf)
  complex*16, intent(in) :: gama(nbmax_g,nLOmax,ntnt_g,ndf_g)
  complex*16, intent(in) :: alfp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: betp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: gamp(mbmax_g,nLOmax,ntnt_g,ndf_g)
  real(8),    intent(in) :: s3r(n_mix,nhow_many_fnc,nhow_many_fnc,nat)   ! s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} > 
  real(8),    intent(in) :: pos(3,ndf)
  integer,    intent(in) :: mult(nat)
  integer,    intent(in) :: nmix(nat)          ! how many product functions per atom
  integer,    intent(in) :: big_l(max_nmix,nat)! L for product basis
  integer,    intent(in) :: nLO_at(lomax+1,nat), ncore(nat)
  real*8,     intent(in) :: cgcoef(*)          ! gaunt coefficients
  integer,    intent(in) :: lmax, loctmatsize, nbmax, mbmax, ntnt, nLOmax, ndf, n_mix, nhow_many_fnc, nat, max_nmix, lomax
  integer,    intent(in) :: nbmax_g,mbmax_g,ntnt_g,ndf_g
  !
  interface
     REAL*8 function getcgcoef(l1,l2,l3,m1,m2,cgcoef)
       REAL*8, intent(in) :: cgcoef(*)
       integer(4), intent(in) :: l1, l2, l3, m1, m2
     end function getcgcoef
  end interface
  !
  real(8), parameter  :: pi = 3.14159265358979d0
  logical   :: PRINT
  integer   :: iat, lb, mb, idf, ie1, ie2, ieq
  real(8)   :: angint
  complex*16:: phs, two_pi_i, phs_angint, ccsm
  integer   :: lc10, lc20, ilo1, ilo2, l2min, l2max, la1, lb1, lc1, la2, lb2, lc2, im, irm
  integer   :: l1, m1, lm1, l2, m2, lm2
  complex*16:: aa, ab, ba, bb, ac(nLOmax), ca(nLOmax), bc(nLOmax), cb(nLOmax), cc(nLOmax,nLOmax), cone, czero
  integer   :: ncl2, ci2, ie1p, cind(nLOmax,(lmax+1)**2)
  complex*16, allocatable :: na(:,:), nb(:,:), nc(:,:), talfp(:,:), tbetp(:,:), tgamp(:,:), mmt2(:,:)
  !
  allocate( na(nend-nstart,(lmax+1)**2), nb(nend-nstart,(lmax+1)**2) )
  allocate( talfp((lmax+1)**2,mend-mstart), tbetp((lmax+1)**2,mend-mstart) )
  allocate( mmt2(nend-nstart,mend-mstart) )
  !
  PRINT = .False.
  two_pi_i = 2.0d0 * pi * cmplx(0.d0,1.d0, 8)
  cone = cmplx(1.d0, 0.d0, 8)
  czero = cmplx(0.d0,0.d0, 8)
  mmat(:,:,:) = cmplx(0.d0, 0.d0, 8)
  !
  if (PRINT) write(6,*) 'calcminm', ' nLOmax=', nLOmax
  
  idf = 0
  im  = 0
  do iat = 1, nat               !! Loop over atoms
     ncl2 = 0
     do l2=0,lmax
        do m2=-l2,l2
           lm2 = l2*l2 + l2 + m2 + 1
           do ilo2 = 1,nLO_at(l2+1,iat)
              ncl2 = ncl2 + 1
              cind(ilo2,lm2) = ncl2
           enddo
        enddo
     enddo
     if (ncl2 > 0) then
        allocate( tgamp(ncl2,mend-mstart), nc(ncl2, nend-nstart) )
        tgamp = czero
        nc    = czero
     endif
     
     do ieq = 1, mult(iat)      !! Loop over equivalent atoms
        idf = idf + 1
        !
        do ie2 = 1,mend-mstart
           do lm2=1,(lmax+1)**2
              talfp(lm2,ie2) = alfp(ie2+mstart,lm2,idf)
              tbetp(lm2,ie2) = betp(ie2+mstart,lm2,idf)
           enddo
        enddo
        if (ncl2 > 0 ) then
           do ie2 = 1,mend-mstart
              do l2=0,lmax
                 do m2=-l2,l2
                    lm2 = l2*l2 + l2 + m2 + 1
                    do ilo2 = 1,nLO_at(l2+1,iat)
                       tgamp(cind(ilo2,lm2),ie2) = gamp(ie2+mstart,ilo2,lm2,idf)
                    enddo
                 enddo
              enddo
           enddo
        endif
        !
        phs = exp( -dot_product( ql, pos(:,idf) )*two_pi_i )
        !
        do irm = 1, nmix(iat)   !!  Loop over mixed functions:
           lb = big_l(irm,iat)
           do mb = -lb,lb
              im = im + 1       !! index of the product basis functions
              na = 0.d0
              nb = 0.d0
              if (ncl2 > 0) nc = 0.d0
              do l1 = 0, lmax              ! index of the single-particle function u_{lm}
                 l2min = iabs(lb-l1)       ! l2 is bounded by triangular rule of <Y_{l1m1}|Y_{l2,m2}|Y_{lb,mb}>
                 l2max = min(lb+l1,lmax)
                 la1 = ncore(iat) + l1 + 1 ! this is index for the matrix element s3r previously computed s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} >
                 lb1 = la1 + lmax + 1      ! why extra 1 I do not understand!!!
                 lc10 = ncore(iat) + 2*(lmax+1) + sum(nLO_at(:l1,iat))  ! where local orbitals start
                 do l2 = l2min,l2max           !! loop over l2:
                    la2 = l2 + ncore(iat) + 1  
                    lb2 = la2 + lmax + 1
                    lc20 = ncore(iat) + 2*(lmax+1) + sum(nLO_at(:l2,iat))  ! where local orbitals start
                    do m1 = -l1,l1             
                       lm1 = l1*l1 + l1 + m1 + 1   ! index for Y_{l1,m1}
                       m2 = m1 - mb                ! This comes from matrix element <Y_{l1m1}|Y_{l2,m2}|Y_{lb,mb}>, which requires m1==mb+m2
                       if ( abs(m2) .gt. l2) cycle 
                       lm2 = l2*l2 + l2 + m2 + 1   ! index for Y_{l2,m2}
                       angint = getcgcoef(l2,lb,l1,m2,mb, cgcoef)  ! gaunt coefficients = <Y_{l1m1}|Y_{l2,m2}|Y_{lb,mb}>* = <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>
                       if ( abs(angint) < 1.0d-8) cycle
                       phs_angint = phs*angint
                       !
                       aa = phs_angint * s3r(irm,la1,la2,iat) ! aa = <u_{irm,lb}|    u_{a1,l1}    u_{a2,l2}> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       ab = phs_angint * s3r(irm,lb1,la2,iat) ! ab = <u_{irm,lb}| udot_{a1,l1}    u_{a2,l2}> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       ba = phs_angint * s3r(irm,la1,lb2,iat) ! ba = <u_{irm,lb}|    u_{a1,l1} udot_{a2,l2}> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       bb = phs_angint * s3r(irm,lb1,lb2,iat) ! bb = <u_{irm,lb}| udot_{a1,l1} udot_{a2,l2}> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       do ilo1 = 1,nLO_at(l1+1,iat) 
                          lc1 = lc10 + ilo1
                          ac(ilo1) = phs_angint * s3r(irm,lc1,la2,iat) ! ac = <u_{irm,lb}| u_{lo,a1,l1}    u_{a2,l2} > <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                          bc(ilo1) = phs_angint * s3r(irm,lc1,lb2,iat) ! bc = <u_{irm,lb}| u_{lo,a1,l1} udot_{a2,l2}|> <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                       enddo
                       do ilo2 = 1,nLO_at(l2+1,iat)
                          lc2 = lc20 + ilo2
                          ca(ilo2) = phs_angint * s3r(irm,la1,lc2,iat) ! ca = < u_{irm,lb} |    u_{a1,l1} ulo_{a2,l2} > <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                          cb(ilo2) = phs_angint * s3r(irm,lb1,lc2,iat) ! cb = < u_{irm,lb} | udot_{a1,l1} ulo_{a2,l2} > <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                          do ilo1 = 1,nLO_at(l1+1,iat)
                              lc1 = lc10 + ilo1
                             cc(ilo1,ilo2) = phs_angint * s3r(irm,lc1,lc2,iat) ! cc = < u_{irm,lb} | ulo_{a1,l1}) ulo_{a2,l2} > <Y_{lb,mb}|Y_{l1,m1}Y*_{l2,m2}>*e^{-iqr}
                          enddo
                       enddo
                       ! na = < u_{irm,lb} |    u_{a2,l2} psi_{ie1} >e^{-iqr} = < u_{irm,lb} | (alfa u_{a1,l1} + beta udot_{a1,l1} + gama ulo_{a1,l1}) *    u_{a2,l2} >e^{-iqr}
                       ! nb = < u_{irm,lb} | udot_{a2,l2} psi_{ie1} >e^{-iqr} = < u_{irm,lb} | (alfa u_{a1,l1} + beta udot_{a1,l1} + gama ulo_{a1,l1}) * udot_{a2,l2} >e^{-iqr}
                       ! nc = < u_{irm,lb} |  ulo_{a2,l2} psi_{ie1} >e^{-iqr} = < u_{irm,lb} | (alfa u_{a1,l1} + beta udot_{a1,l1} + gama ulo_{a1,l1}) *  ulo_{a2,l2} >e^{-iqr}
                       do ie1=1,nend-nstart                        !! Loop over eigenfunctions at k
                          ie1p = ie1+nstart
                          na(ie1,lm2) = na(ie1,lm2) + aa * alfa(ie1p,lm1,idf) + ab * beta(ie1p,lm1,idf) 
                          nb(ie1,lm2) = nb(ie1,lm2) + ba * alfa(ie1p,lm1,idf) + bb * beta(ie1p,lm1,idf) 
                          do ilo1 = 1,nLO_at(l1+1,iat)
                             na(ie1,lm2) = na(ie1,lm2) + ac(ilo1) * gama(ie1p,ilo1,lm1,idf) ! na += <u_{irm,lb}| u_{lo,a1,l1}    u_{a2,l2} > * gama(ie1,ilo,l1m1)<Y|Y|Y>e^{-iqr}
                             nb(ie1,lm2) = nb(ie1,lm2) + bc(ilo1) * gama(ie1p,ilo1,lm1,idf) ! nb += <u_{irm,lb}| u_{lo,a1,l1} udot_{a2,l2}|> * gama(ie1,ilo,l1m1)<Y|Y|Y>e^{-iqr}
                          enddo
                          do ilo2 = 1,nLO_at(l2+1,iat)
                             ccsm = 0.d0
                             do ilo1 = 1,nLO_at(l1+1,iat) 
                                ccsm = ccsm + cc(ilo1,ilo2) * gama(ie1p,ilo1,lm1,idf)
                             enddo
                             ci2 = cind(ilo2,lm2)
                             nc(ci2,ie1) = nc(ci2,ie1) + ccsm + ca(ilo2) * alfa(ie1p,lm1,idf) + cb(ilo2) * beta(ie1p,lm1,idf)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              
              !! sm = < u_{irm,lb} | psi^*_{ie2} psi_{ie1} >e^{-iqr}
              call zgemm('n','n', nend-nstart,mend-mstart,(lmax+1)**2, cone, na,nend-nstart, talfp,(lmax+1)**2, czero, mmt2,nend-nstart)
              call zgemm('n','n', nend-nstart,mend-mstart,(lmax+1)**2, cone, nb,nend-nstart, tbetp,(lmax+1)**2, cone,  mmt2,nend-nstart)
              if (ncl2 > 0) then
                 call zgemm('t','n', nend-nstart,mend-mstart,ncl2, cone, nc,ncl2, tgamp,ncl2, cone,  mmt2,nend-nstart)
              endif
              mmat(:,:,im) = mmt2(:,:)
              if (.False.) then
                 do ie1=1,nend-nstart
                    do ie2 = 1,mend-mstart
                       write(6,'(A,3I4,1x,2F16.10)') 'mm=', im,ie2,ie1, mmat(ie1,ie2,im)
                    end do
                 end do
              endif
           enddo
        enddo
     enddo
     if (ncl2 > 0) then
        deallocate( tgamp, nc )
     endif
  enddo
  deallocate( mmt2 )
  deallocate( talfp, tbetp )
  deallocate( na, nb )
end subroutine calc_minm_MT

subroutine calc_minm_IS(mmat, nstart,nend,mstart,mend, Aeigk,Aeigq, iumklap, mpwipw,  nvi,nvj, indgk_i, indgk_j, indgq_inverse, gindex, i_g0, Vol, ngq, ngq_barc, ng, n1, n2, n3, ngi, ngj, inv_size, nbmax_i, nbmax_j)
  ! Calculates overlap between product basis and two single-particle eigenvectors from vector file in the interstitials.
  !    Given the product basis in the interstitials is 1/sqrt(O) * e^{i*Gp*r}  (Gp goes up to ngq),
  !    it computes the integral
  !                              1/sqrt(O)*< e^{i*iGp*r}|psi_{ie2}^* |psi_{ie1}>_{Interstitials}
  !    which is evaluated like   1/sqrt(O)*< e^{i*iGp*r}|e^{i(-Gj+Gi)*r > Aeig^*(Gj,ie2)*Aeig(Gi,ie1) or
  !          mmat(ie1,ie2,iGp) = 1/sqrt(O) * mpwipw(iGp,Gi-Gj) * Aeig^*(Gj,ie2)*Aeig(Gi,ie1)
  implicit none
  integer,    intent(in) :: nvi, nvj, ngq, ngq_barc, ng, n1, n2, n3, ngi, ngj, inv_size, nbmax_i, nbmax_j
  integer,    intent(in) :: nstart,nend, mstart,mend     ! ranges of bands for k and k-q point
  complex*16, intent(out):: mmat(nend-nstart,mend-mstart,ngq)
  integer,    intent(in) :: iumklap(3)
  complex*16, intent(in) :: mpwipw(ngq,ngq_barc)        ! 1/sqrt(O)<G|K>_{inter}
  integer,    intent(in) :: indgk_i(ngi), indgk_j(ngj)! reciprocal vectors indices from vector file for ki(==k) and kj(==k+q)
  integer,    intent(in) :: indgq_inverse(inv_size)     ! inverse index of indgq for product basis
  integer,    intent(in) :: gindex(ng,3)                ! all G-points
  integer,    intent(in) :: i_g0(2*n1+1,2*n2+1,2*n3+1)  ! inverse index that translates from G point to index in gindex
  complex*16, intent(in) :: Aeigk(nbmax_i,ngi)          ! eigenvector from vector file for k-point
  complex*16, intent(in) :: Aeigq(nbmax_j,ngj)          ! eigenvector from vector file for q-point
  real(8),    intent(in) :: Vol
  !
  integer   :: iGv(3)
  integer   :: idg, iGp, ie1, ie2, ivi, ivj, ndim, mdim, ierr, iK
  real(8)   :: sqvi
  complex*16:: cone, czero
  integer,    allocatable :: iKpw(:,:)    
  complex*16, allocatable :: tmat(:,:), tmat2(:,:), mnn(:,:)
  logical   :: PRINT
  PRINT = .False.
  !
  cone  = cmplx(1.d0,0.d0,8)
  czero = cmplx(0.d0,0.d0,8)
  sqvi  = 1.d0/sqrt(Vol)

  if (PRINT) then
     print *, 'ngq=', ngq, 'ngq_barc=', ngq_barc
     print *, 'mpwipw='
     do iGp = 1,ngq
        do iK = 1,ngq
           write(6,'(I4,1x,I4,1x,2F12.8)') iGp, iK, mpwipw(iGp, iK)
        enddo
     enddo
  endif

  ndim = nend - nstart
  mdim = mend - mstart
  allocate( iKpw(nvi,nvj), tmat(nvj,nvi), tmat2(nvj,ndim), mnn(mdim,ndim), stat=ierr)
  if (ierr.ne.0) then
     print *, "Fail to allocate work arrays in calc_minm_IS"
  end if
  
  iKpw(:,:)=0
  do ivj = 1,nvj 
     do ivi = 1,nvi
        iGv = gindex(indgk_i(ivi)+1,:) - gindex(indgk_j(ivj)+1,:) + iumklap(:)  ! G_ki - G_kj + umklap : reciprocal vectors from vector file
        idg = 0
        if ( abs(iGv(1))<=n1 .and. abs(iGv(2))<=n2 .and. abs(iGv(3))<=n3 ) idg = i_g0(iGv(1)+n1+1,iGv(2)+n2+1,iGv(3)+n3+1)+1  ! is this difference not too large, and is in our fixed mesh?
        if (idg > 0 .and. idg<=inv_size) iKpw(ivi,ivj) = indgq_inverse(idg)+1     ! now we get the index in ngq_barc, for which mpwipw = 1/sqrt(olap)<G|I|K> was computed. This is index for K
        !write(6,'(I3,1x,I3,3x,3I3,2x,I4,2x,I7,2x,3I3,2x,3I3,2x,3I3)') ivj, ivi, iGv, idg, iKpw(ivi,ivj), gindex(indgk_i(ivi)+1,:), gindex(indgk_j(ivj)+1,:), iumklap(:)
     enddo 
  enddo
  tmat = czero
  do iGp = 1,ngq                          ! over G in the product-interstitial basis
     do ivi = 1,nvi                       ! G1 from vector file at ki==k
        do ivj = 1,nvj                    ! G2 from vector file at kj==k+q
           if ( iKpw(ivi,ivj).gt.0 ) then ! G1-G2 is iKpw in pw.indgq[:,iKpw]
              tmat(ivj,ivi) = mpwipw(iGp, iKpw(ivi,ivj)) ! 1/sqrt(O)*<G|I|K> = 1/sqrt(O)*Int[e^{i(G_{ki}-G_{kj}-iGp)*r},inter] = 1/sqrt(O)*<iGp|I|G_{ki}-G_{kj}> 
              !write(6,'(I4,1x,I4,1x,I4,1x,2F12.8)') iGp, ivi, ivj, tmat(ivj,ivi)
           endif
        enddo
     enddo

     if (PRINT) then
        print *, 'tmat='
        do ivi=1,nvi
           do ivj=1,nvj
              write(6,'(I3,1x,I3,1x,I3,3x,2F12.8)') iGp, ivi, ivj, tmat(ivj,ivi)
           enddo
        enddo
     endif
     
     ! Computes the integral : 1/sqrt(O)*< e^{i*iGp*r}|psi_{ie2,k+q}^* |psi_{ie1,k}>_{Interstitials}
     !
     !  tmat2(nvj,ie1) = tmat(nvj,nvi) * Aeigk.T(nvi,nstart+ie1)
     call zgemm('n','t', nvj,ndim,nvi,  cone, tmat,nvj, Aeigk(nstart+1:nend,:nvi),ndim, czero, tmat2,nvj )
     !  mnn(ie2,ie1) = Aeigq(mstart+ie2+1,nvj)^* * tmat(nvj,nvi) * Aeig.T(nvi,nstart+ie1)
     call zgemm('n','n', mdim,ndim,nvj, cone, Aeigq(mstart+1:mend,:nvj),mdim, tmat2,nvj, czero, mnn,mdim)
     do ie1 = 1,nend-nstart
        do ie2 = 1,mend-mstart
           mmat(ie1,ie2,iGp) = mnn(ie2,ie1)*sqvi
        enddo
     enddo
     
     if (PRINT) then
        print *, 'mnn=', 'sqvi=', sqvi
        do ie1=1,ndim
           do ie2=1,mdim
              write(6,'(I3,1x,I3,1x,I3,3x,2F12.8)') iGp, ie1, ie2, mmat(ie1,ie2,iGp) !mnn(ie2,ie1)
           enddo
        enddo
     endif
     
  enddo
  deallocate(iKpw,tmat,tmat2,mnn)
end subroutine calc_minm_IS

subroutine calc_minc(mmat,kl,nstart,nend,cstart,cend, alfa,beta,gama,s3r,corind, pos, mult, nmix, big_l, nLO_at, ncore, cgcoef, lmax, loctmatsize, nbmax, nat, ntnt, nLOmax, ndf, n_mix, nhow_many_fnc, max_nmix, lomax, nbmax_g,ntnt_g,ndf_g,ncg)
  ! calculates matrix elements <Product_Basis| psi_{i1,k} psi_{icore}^*> e^{-ikr}, where i1\in[nstart,nend] and icore\in[cstart,cent]
  !
  implicit none
  integer,    intent(in) :: nstart,nend, cstart,cend, loctmatsize
  integer,    intent(in) :: lmax,nbmax,ntnt,nLOmax,ndf,n_mix,nhow_many_fnc,nat,max_nmix,lomax,ncg!,lmixmax
  integer,    intent(in) :: nbmax_g,ntnt_g,ndf_g
  complex*16, intent(out):: mmat(nend-nstart,cend-cstart,loctmatsize)
  real(8),    intent(in) :: kl(3)  ! k in semi-cartesian coordinates
  complex*16, intent(in) :: alfa(nbmax,ntnt,ndf)
  complex*16, intent(in) :: beta(nbmax,ntnt,ndf)
  complex*16, intent(in) :: gama(nbmax_g,nLOmax,ntnt_g,ndf_g)
  real(8),    intent(in) :: s3r(n_mix,nhow_many_fnc,nhow_many_fnc,nat)   ! s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} > 
  real(8),    intent(in) :: pos(3,ndf)
  integer,    intent(in) :: mult(nat)
  integer,    intent(in) :: nmix(nat)          ! how many product functions per atom
  integer,    intent(in) :: big_l(max_nmix,nat)! L for product basis
  integer,    intent(in) :: nLO_at(lomax+1,nat), ncore(nat)
  real*8,     intent(in) :: cgcoef(*)          ! gaunt coefficients
  integer,    intent(in) :: corind(ncg,5)
  !
  interface
     REAL*8 function getcgcoef(l1,l2,l3,m1,m2,cgcoef)
       REAL*8, intent(in) :: cgcoef(*)
       integer(4), intent(in) :: l1, l2, l3, m1, m2
     end function getcgcoef
  end interface
  !
  real(8), parameter  :: pi = 3.14159265358979d0
  !
  integer    :: lb, mb, iat, idf, ic, l1, m1, icg, ie1, ie1p, ilo1, im, irm, l2, l1min, l1max, la1, lb1, lc1, lc10, lm1, m2, ieq, im0
  complex*16 :: aa, bb, cc(nLOmax)
  real(8)    :: angint
  complex*16 :: two_pi_i, csm, phs, phs_angint
  logical, parameter   :: EXPECTED_WRONG_BUT_EQUAL_TO_PYGAP = .False.
  integer, allocatable :: locmixind3(:,:)
  !
  allocate( locmixind3(max_nmix,ndf) )
  locmixind3(:,:) = 0
  
  idf=0
  im=0
  do iat=1,nat
     do ieq=1,mult(iat)
        idf=idf+1
        do irm=1,nmix(iat)
           lb = big_l(irm,iat)
           locmixind3(irm,idf) = im
           im = im + 2*lb+1
        enddo
     enddo
  enddo
  

  two_pi_i = 2.0d0 * pi * cmplx(0.d0,1.d0, 8)
  mmat = 0.d0

  !$OMP PARALLEL DO PRIVATE(icg,iat,idf,ic,l2,m2,phs,irm,lb,im0,mb,im,l1min,l1max,l1,la1,lb1,lc10,m1,angint,lm1,phs_angint,aa,bb,cc,lc1,ie1,ie1p,csm,ilo1)&
  !$OMP& SHARED(mmat)&
  !$OMP& SCHEDULE(STATIC)
  do icg = 1, ncg                                 !! loop over core states, including l,m in core
     if( icg <= cstart .or. icg > cend) cycle
     iat = corind(icg,1)+1
     idf = corind(icg,2)+1
     ic  = corind(icg,3)+1
     l2  = corind(icg,4)
     m2  = corind(icg,5)
     !print *, 'iat=', iat, 'idf=', idf, 'ic=', ic, 'l2=', l2, 'm2=', m2
     phs = exp( -dot_product( kl, pos(:,idf) )*two_pi_i )  ! because the core state does not have e^{(ik-iq)r} phase factor. Maybe we should assume it has?
     do irm = 1, nmix(iat)                        !! Loop over product functions
        lb = big_l(irm,iat)                       ! L of the product function
        im0 = locmixind3(irm,idf)
        do mb = -lb,lb                            ! M of the product function
           im = im0 + lb + mb + 1
           l1min = iabs(lb-l2)        ! l1 is bounded by triangular rule of <Y_{l2m2}|Y_{l1,m1}|Y_{lb,mb}>
           l1max = min(lb+l2,lmax)
           do l1 = l1min,l1max           ! loop over l1
              la1 = l1 + ncore(iat) + 1  ! index of the valence state u_{l1}
              lb1 = la1 + lmax + 1       ! index of the valence state udot_{l1} 
              lc10 = ncore(iat) + 2*(lmax+1) + sum(nLO_at(:l1,iat))  ! where local orbitals start for u_lo_{l1}
              if (EXPECTED_WRONG_BUT_EQUAL_TO_PYGAP) then
                 m1 = m2 - mb            ! This comes from matrix element <Y_{l2m2}|Y_{l1,m1}|Y_{lb,mb}>, which requires m2==mb+m1
                 angint = getcgcoef(l1,lb,l2,m1,mb, cgcoef)  ! gaunt coefficients = <Y_{l2m2}|Y_{l1,m1}|Y_{lb,mb}>* = <Y_{lb,mb}|Y_{l2,m2}Y*_{l1,m1}>
              else
                 ! we need <Y_{lb,mb}| Y*_{l2,m2} Y_{l1,m1}>
                 m1 = m2 + mb
                 angint = getcgcoef(l2,lb,l1,m2,mb, cgcoef)  ! gaunt coefficients = <Y_{lb,mb}| Y*_{l2,m2} Y_{l1,m1}>== < Y_{l1,m1}| Y_{lb,mb} Y_{l2,m2}>^*
              endif
              if ( abs(m1).gt.l1) cycle 
              lm1 = l1*l1 + l1 + m1 + 1                   ! index for Y_{l1,m1}
              
              if ( abs(angint) < 1.0d-8) cycle
              phs_angint = phs*angint
              !<Product_Basis| psi_{icore}^* psi_{i1,k} >
              aa =  phs_angint * s3r(irm,ic,la1,iat) ! aa = <u_{irm,lb}| uc_{a1,l2}    u_{a1,l1}> <Y_{lb,mb}|Y*_{l2,m2}Y_{l1,m1}>*e^{-ikR}
              bb =  phs_angint * s3r(irm,ic,lb1,iat) ! bb = <u_{irm,lb}| uc_{a1,l2} udot_{a1,l1}> <Y_{lb,mb}|Y*_{l2,m2}Y_{l1,m1}>*e^{-ikR}
              do ilo1 = 1,nLO_at(l1+1,iat)
                 lc1 = lc10 + ilo1
                 cc(ilo1) = phs_angint * s3r(irm,ic,lc1,iat) ! cc = <u_{irm,lb}| uc_{a1,l2} ulo_{a1,l1}> <Y_{lb,mb}|Y*_{l2,m2}Y_{l1,m1}>*e^{-ikR}
              enddo
              do ie1 = 1,nend-nstart
                 ie1p = ie1+nstart
                 csm = aa * alfa(ie1p,lm1,idf) + bb * beta(ie1p,lm1,idf) ! <u_{irm,lb}| uc_{a1,l2}      u_{a1,l1}> <Y_{lb,mb}| Y*_{l2,m2} Y_{l1,m1}> * e^{-ikR}*alfa(ie1,l1,m1)
                 do ilo1 = 1,nLO_at(l1+1,iat)                            !+<u_{irm,lb}| uc_{a1,l2}   udot_{a1,l1}> <Y_{lb,mb}| Y*_{l2,m2} Y_{l1,m1}> * e^{-ikR}*beta(ie1,l1,m1)
                    csm = csm + cc(ilo1) * gama(ie1p,ilo1,lm1,idf)       !+<u_{irm,lb}| uc_{a1,l2}ulo_{ilo,a1,l1}> <Y_{lb,mb}| Y*_{l2,m2} Y_{l1,m1}> * e^{-ikR}*gama(ie1,l1,m1)
                 enddo
                 mmat(ie1,icg-cstart,im) = mmat(ie1,icg-cstart,im) + csm
              enddo
           enddo
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  deallocate( locmixind3 )
  if (.False.) then
     print *, 'calc_minc'
     do im=1,loctmatsize
        do ie1=1,nend-nstart
           do icg=1,cend-cstart
              write(6,'(I3,1x,I3,1x,I3,1x,2F14.10)') im, ie1+nstart, icg+cstart, mmat(ie1,icg,im)  ! icg==ie2
           enddo
        enddo
     enddo
  endif
end subroutine calc_minc

subroutine calc_minc2(mmat,kmql,cstart,cend, mstart,mend,alfp,betp,gamp,s3r,corind, pos, mult, nmix, big_l, nLO_at, ncore, cgcoef, lmax, loctmatsize, mbmax, nat, ntnt, nLOmax, ndf, n_mix, nhow_many_fnc, max_nmix, lomax, mbmax_g,ntnt_g,ndf_g,ncg)
  ! calculates matrix elements <Product_Basis| psi_{icore} psi_{i2,k-q}^* >e^{i(k-q)r}, where i2\in[mstart,mend] and icore\in[cstart,cent]
  implicit none
  integer,    intent(in) :: mstart,mend, cstart,cend, loctmatsize
  integer,    intent(in) :: lmax,mbmax,ntnt,nLOmax,ndf,n_mix,nhow_many_fnc,nat,max_nmix,lomax,ncg!,lmixmax
  integer,    intent(in) :: mbmax_g,ntnt_g,ndf_g
  complex*16, intent(out):: mmat(cend-cstart,mend-mstart,loctmatsize)
  real(8),    intent(in) :: kmql(3)  ! k in semi-cartesian coordinates
  complex*16, intent(in) :: alfp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: betp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: gamp(mbmax_g,nLOmax,ntnt_g,ndf_g)
  real(8),    intent(in) :: s3r(n_mix,nhow_many_fnc,nhow_many_fnc,nat)   ! s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} > 
  real(8),    intent(in) :: pos(3,ndf)
  integer,    intent(in) :: mult(nat)
  integer,    intent(in) :: nmix(nat)          ! how many product functions per atom
  integer,    intent(in) :: big_l(max_nmix,nat)! L for product basis
  integer,    intent(in) :: nLO_at(lomax+1,nat), ncore(nat)
  real*8,     intent(in) :: cgcoef(*)          ! gaunt coefficients
  integer,    intent(in) :: corind(ncg,5)
  !
  interface
     REAL*8 function getcgcoef(l1,l2,l3,m1,m2,cgcoef)
       REAL*8, intent(in) :: cgcoef(*)
       integer(4), intent(in) :: l1, l2, l3, m1, m2
     end function getcgcoef
  end interface
  !
  real(8), parameter  :: pi = 3.14159265358979d0
  !
  integer    :: lb, mb, iat, idf, ic, l1, m1, icg, ie2, ie2p, ilo1, im, irm, l2, l1min, l1max, la1, lb1, lc1, lc10, lm1, m2, ieq, im0
  complex*16 :: aa, bb, cc(nLOmax)
  real(8)    :: angint
  complex*16 :: two_pi_i, csm, phs, phs_angint
  logical, parameter   :: EXPECTED_WRONG_BUT_EQUAL_TO_PYGAP = .False.
  integer, allocatable :: locmixind3(:,:)
  !
  allocate( locmixind3(max_nmix,ndf) )
  locmixind3(:,:) = 0
  
  idf=0
  im=0
  do iat=1,nat
     do ieq=1,mult(iat)
        idf=idf+1
        do irm=1,nmix(iat)
           lb = big_l(irm,iat)
           locmixind3(irm,idf) = im
           im = im + 2*lb+1
        enddo
     enddo
  enddo
  

  two_pi_i = 2.0d0 * pi * cmplx(0.d0,1.d0, 8)
  mmat = 0.d0
  !$OMP PARALLEL DO PRIVATE(icg,iat,idf,ic,l2,m2,phs,irm,lb,im0,mb,im,l1min,l1max,l1,la1,lb1,lc10,m1,angint,lm1,phs_angint,aa,bb,cc,lc1,ie2,ie2p,csm,ilo1)&
  !$OMP& SHARED(mmat)&
  !$OMP& SCHEDULE(STATIC)
  do icg = 1, ncg                                 !! loop over core states, including l,m in core
     if( icg <= cstart .or. icg > cend) cycle
     iat = corind(icg,1)+1
     idf = corind(icg,2)+1
     ic  = corind(icg,3)+1
     l2  = corind(icg,4) ! == lc
     m2  = corind(icg,5) ! == mc
     !print *, 'iat=', iat, 'idf=', idf, 'ic=', ic, 'l2=', l2, 'm2=', m2
     phs = exp( dot_product( kmql, pos(:,idf) )*two_pi_i )  ! because the core state does not have e^{(ik-iq)r} phase factor. Maybe we should assume it has?
     do irm = 1, nmix(iat)                        !! Loop over product functions
        lb = big_l(irm,iat)                       ! L of the product function
        im0 = locmixind3(irm,idf)
        do mb = -lb,lb                            ! M of the product function
           im = im0 + lb + mb + 1
           l1min = iabs(lb-l2)        ! l1 is bounded by triangular rule of <Y_{l2m2}|Y_{l1,m1}|Y_{lb,mb}>
           l1max = min(lb+l2,lmax)
           do l1 = l1min,l1max           ! loop over l1
              la1 = l1 + ncore(iat) + 1  ! index of the valence state u_{l1}
              lb1 = la1 + lmax + 1       ! index of the valence state udot_{l1} 
              lc10 = ncore(iat) + 2*(lmax+1) + sum(nLO_at(:l1,iat))  ! where local orbitals start for u_lo_{l1}
              if (EXPECTED_WRONG_BUT_EQUAL_TO_PYGAP) then
                 ! we need <Y_{lb,mb}| Y*_{lc,mc} Y_{l1,m1}>, here m2=mc and l2=lc
                 m1 = m2 + mb
                 angint = getcgcoef(l2,lb,l1,m2,mb, cgcoef)  ! gaunt coefficients = <Y_{lb,mb}| Y*_{lc,mc} Y_{l1,mv+mb}>== < Y_{l1,mc+mb}| Y_{lb,mb} Y_{lc,mc}>^*
              else
                 m1 = m2 - mb            ! This comes from matrix element <Y_{l2m2}|Y_{l1,m1}|Y_{lb,mb}>, which requires m2==mb+m1
                 angint = getcgcoef(l1,lb,l2,m1,mb, cgcoef)  ! gaunt coefficients = <Y_{l2,m1+mb}|Y_{l1,m1} Y_{lb,mb}>* = <Y_{lb,mb} | Y_{l2,m2} Y*_{l1,m1}>
              endif
              if ( abs(m1).gt.l1) cycle 
              lm1 = l1*l1 + l1 + m1 + 1                   ! index for Y_{l1,m1}
              
              if ( abs(angint) < 1.0d-8) cycle
              phs_angint = phs*angint
              !<Product_Basis| psi_{icore} psi_{i2,k-q}^* >
              aa =  phs_angint * s3r(irm,ic,la1,iat) ! aa = <u_{irm,lb}| uc_{a1,l2}    u_{a1,l1}> <Y_{lb,mb}|Y_{l2,m2} Y*_{l1,m1}>*e^{i(k-q)R}
              bb =  phs_angint * s3r(irm,ic,lb1,iat) ! bb = <u_{irm,lb}| uc_{a1,l2} udot_{a1,l1}> <Y_{lb,mb}|Y_{l2,m2} Y*_{l1,m1}>*e^{i(k-q)R}
              do ilo1 = 1,nLO_at(l1+1,iat)
                 lc1 = lc10 + ilo1
                 cc(ilo1) = phs_angint * s3r(irm,ic,lc1,iat) ! cc = <u_{irm,lb}| uc_{a1,l2} ulo_{a1,l1}> <Y_{lb,mb}|Y_{l2,m2} Y*_{l1,m1}>*e^{i(k-q)R}
              enddo
              do ie2 = 1,mend-mstart
                 ie2p = ie2+mstart
                 csm = aa * alfp(ie2p,lm1,idf) + bb * betp(ie2p,lm1,idf) ! <u_{irm,lb}| uc_{a1,l2}      u_{a1,l1}> <Y_{lb,mb}| Y_{lc,mc} Y*_{l1,m1}> * e^{i(k-q)R}*alfp(ie2,l1,m1)
                 do ilo1 = 1,nLO_at(l1+1,iat)                            !+<u_{irm,lb}| uc_{a1,l2}   udot_{a1,l1}> <Y_{lb,mb}| Y_{lc,mc} Y*_{l1,m1}> * e^{i(k-q)R}*betp(ie2,l1,m1)
                    csm = csm + cc(ilo1) * gamp(ie2p,ilo1,lm1,idf)       !+<u_{irm,lb}| uc_{a1,l2}ulo_{ilo,a1,l1}> <Y_{lb,mb}| Y_{lc,mc} Y*_{l1,m1}> * e^{i(k-q)R}*gamp(ie2,l1,m1)
                 enddo
                 mmat(icg-cstart,ie2,im) = mmat(icg-cstart,ie2,im) + csm
              enddo
           enddo
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  deallocate( locmixind3 )
  if (.False.) then
     print *, 'calc_minc2'
     do im=1,loctmatsize
        do ie2=1,mend-mstart
           do icg=1,cend-cstart
              write(6,'(I3,1x,I3,1x,I3,1x,2F14.10)') im, ie2+mstart, icg+cstart, mmat(icg,ie2,im)  ! icg==ie2
           enddo
        enddo
     enddo
     print *, 'calc_minc2_end'
  endif
end subroutine calc_minc2

subroutine calc_minc3(mmat,kmql,cstart,cend, mstart,mend,alfp,betp,gamp,s3r,corind, pos, mult, nmix, big_l, nLO_at, ncore, cgcoef, lmax, loctmatsize, mbmax, nat, ntnt, nLOmax, ndf, n_mix, nhow_many_fnc, max_nmix, lomax, mbmax_g,ntnt_g,ndf_g,ncg)
  ! calculates matrix elements <Product_Basis| psi_{icore} psi_{i2,k-q}^* >e^{i(k-q)r}, where i2\in[mstart,mend] and icore\in[cstart,cent]
  implicit none
  integer,    intent(in) :: mstart,mend, cstart,cend, loctmatsize
  integer,    intent(in) :: lmax,mbmax,ntnt,nLOmax,ndf,n_mix,nhow_many_fnc,nat,max_nmix,lomax,ncg!,lmixmax
  integer,    intent(in) :: mbmax_g,ntnt_g,ndf_g
  complex*16, intent(out):: mmat(cend-cstart,mend-mstart,loctmatsize)
  real(8),    intent(in) :: kmql(3)  ! k in semi-cartesian coordinates
  complex*16, intent(in) :: alfp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: betp(mbmax,ntnt,ndf)
  complex*16, intent(in) :: gamp(mbmax_g,nLOmax,ntnt_g,ndf_g)
  real(8),    intent(in) :: s3r(n_mix,nhow_many_fnc,nhow_many_fnc,nat)   ! s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} > 
  real(8),    intent(in) :: pos(3,ndf)
  integer,    intent(in) :: mult(nat)
  integer,    intent(in) :: nmix(nat)          ! how many product functions per atom
  integer,    intent(in) :: big_l(max_nmix,nat)! L for product basis
  integer,    intent(in) :: nLO_at(lomax+1,nat), ncore(nat)
  real*8,     intent(in) :: cgcoef(*)          ! gaunt coefficients
  integer,    intent(in) :: corind(ncg,5)
  !
  interface
     REAL*8 function getcgcoef(l1,l2,l3,m1,m2,cgcoef)
       REAL*8, intent(in) :: cgcoef(*)
       integer(4), intent(in) :: l1, l2, l3, m1, m2
     end function getcgcoef
  end interface
  !
  real(8), parameter  :: pi = 3.14159265358979d0
  !
  integer    :: lb, mb, iat, idf, ic, l1, m1, icg, ie2, ie2p, ilo1, im, irm, l2, l1min, l1max, la1, lb1, lc1, lc10, lm1, m2, ieq, im0
  complex*16 :: aa, bb, cc(nLOmax)
  real(8)    :: angint
  complex*16 :: two_pi_i, csm, phs, phs_angint
  logical, parameter   :: EXPECTED_WRONG_BUT_EQUAL_TO_PYGAP = .False.
  logical, parameter   :: orig = .True.
  integer, allocatable :: locmixind3(:,:)
  complex*16, allocatable :: aax(:,:), bbx(:,:), ccx(:,:,:)
  !
  complex*16, allocatable :: res(:,:)
  complex*16 :: czero, cone
  integer :: nb, lmx
  !
  allocate( locmixind3(max_nmix,ndf) )
  locmixind3(:,:) = 0
  
  idf=0
  im=0
  do iat=1,nat
     do ieq=1,mult(iat)
        idf=idf+1
        do irm=1,nmix(iat)
           lb = big_l(irm,iat)
           locmixind3(irm,idf) = im
           im = im + 2*lb+1
        enddo
     enddo
  enddo

  if (.not. orig) then
     allocate( aax((lmax+1)**2,loctmatsize) )
     allocate( bbx((lmax+1)**2,loctmatsize) )
     allocate( ccx(nLOmax,(lmax+1)**2,loctmatsize) )
     nb = mend-mstart
     allocate( res(nb,loctmatsize) )
  endif
  
  two_pi_i = 2.0d0 * pi * cmplx(0.d0,1.d0, 8)
  mmat = 0.d0
  do icg = 1, ncg                                 !! loop over core states, including l,m in core
     if( icg <= cstart .or. icg > cend) cycle
     iat = corind(icg,1)+1
     idf = corind(icg,2)+1
     ic  = corind(icg,3)+1
     l2  = corind(icg,4) ! == lc
     m2  = corind(icg,5) ! == mc
     !print *, 'iat=', iat, 'idf=', idf, 'ic=', ic, 'l2=', l2, 'm2=', m2

     if (.not. orig) then
        aax = 0.d0
        bbx = 0.d0
        ccx = 0.d0
     endif
     
     phs = exp( dot_product( kmql, pos(:,idf) )*two_pi_i )  ! because the core state does not have e^{(ik-iq)r} phase factor. Maybe we should assume it has?
     do irm = 1, nmix(iat)                        !! Loop over product functions
        lb = big_l(irm,iat)                       ! L of the product function
        im0 = locmixind3(irm,idf)
        do mb = -lb,lb                            ! M of the product function
           im = im0 + lb + mb + 1

           l1min = iabs(lb-l2)        ! l1 is bounded by triangular rule of <Y_{l2m2}|Y_{l1,m1}|Y_{lb,mb}>
           l1max = min(lb+l2,lmax)
           do l1 = l1min,l1max           ! loop over l1
              la1 = l1 + ncore(iat) + 1  ! index of the valence state u_{l1}
              lb1 = la1 + lmax + 1       ! index of the valence state udot_{l1} 
              lc10 = ncore(iat) + 2*(lmax+1) + sum(nLO_at(:l1,iat))  ! where local orbitals start for u_lo_{l1}
              if (EXPECTED_WRONG_BUT_EQUAL_TO_PYGAP) then
                 ! we need <Y_{lb,mb}| Y*_{lc,mc} Y_{l1,m1}>, here m2=mc and l2=lc
                 m1 = m2 + mb
                 angint = getcgcoef(l2,lb,l1,m2,mb, cgcoef)  ! gaunt coefficients = <Y_{lb,mb}| Y*_{lc,mc} Y_{l1,mv+mb}>== < Y_{l1,mc+mb}| Y_{lb,mb} Y_{lc,mc}>^*
              else
                 m1 = m2 - mb            ! This comes from matrix element <Y_{l2m2}|Y_{l1,m1}|Y_{lb,mb}>, which requires m2==mb+m1
                 angint = getcgcoef(l1,lb,l2,m1,mb, cgcoef)  ! gaunt coefficients = <Y_{l2,m1+mb}|Y_{l1,m1} Y_{lb,mb}>* = <Y_{lb,mb} | Y_{l2,m2} Y*_{l1,m1}>
              endif
              if ( abs(m1).gt.l1) cycle 
              lm1 = l1*l1 + l1 + m1 + 1                   ! index for Y_{l1,m1}
              
              if ( abs(angint) < 1.0d-8) cycle
              phs_angint = phs*angint
              !<Product_Basis| psi_{icore} psi_{i2,k-q}^* >
              if (orig) then
                 aa =  phs_angint * s3r(irm,ic,la1,iat) ! aa = <u_{irm,lb}| uc_{a1,l2}    u_{a1,l1}> <Y_{lb,mb}|Y_{l2,m2} Y*_{l1,m1}>*e^{i(k-q)R}
                 bb =  phs_angint * s3r(irm,ic,lb1,iat) ! bb = <u_{irm,lb}| uc_{a1,l2} udot_{a1,l1}> <Y_{lb,mb}|Y_{l2,m2} Y*_{l1,m1}>*e^{i(k-q)R}
                 do ilo1 = 1,nLO_at(l1+1,iat)
                    lc1 = lc10 + ilo1
                    cc(ilo1) = phs_angint * s3r(irm,ic,lc1,iat) ! cc = <u_{irm,lb}| uc_{a1,l2} ulo_{a1,l1}> <Y_{lb,mb}|Y_{l2,m2} Y*_{l1,m1}>*e^{i(k-q)R}
                 enddo
                 do ie2 = 1,mend-mstart
                    ie2p = ie2+mstart
                    csm = aa * alfp(ie2p,lm1,idf) + bb * betp(ie2p,lm1,idf) ! <u_{irm,lb}| uc_{a1,l2}      u_{a1,l1}> <Y_{lb,mb}| Y_{lc,mc} Y*_{l1,m1}> * e^{i(k-q)R}*alfp(ie2,l1,m1)
                    do ilo1 = 1,nLO_at(l1+1,iat)                            !+<u_{irm,lb}| uc_{a1,l2}   udot_{a1,l1}> <Y_{lb,mb}| Y_{lc,mc} Y*_{l1,m1}> * e^{i(k-q)R}*betp(ie2,l1,m1)
                       csm = csm + cc(ilo1) * gamp(ie2p,ilo1,lm1,idf)       !+<u_{irm,lb}| uc_{a1,l2}ulo_{ilo,a1,l1}> <Y_{lb,mb}| Y_{lc,mc} Y*_{l1,m1}> * e^{i(k-q)R}*gamp(ie2,l1,m1)
                    enddo
                    mmat(icg-cstart,ie2,im) = mmat(icg-cstart,ie2,im) + csm
                 enddo
              else
                 aax(lm1,im) =  aax(lm1,im) + phs_angint * s3r(irm,ic,la1,iat) ! aa = <u_{irm,lb}| uc_{a1,l2}    u_{a1,l1}> <Y_{lb,mb}|Y_{l2,m2} Y*_{l1,m1}>*e^{i(k-q)R}
                 bbx(lm1,im) =  bbx(lm1,im) + phs_angint * s3r(irm,ic,lb1,iat) ! bb = <u_{irm,lb}| uc_{a1,l2} udot_{a1,l1}> <Y_{lb,mb}|Y_{l2,m2} Y*_{l1,m1}>*e^{i(k-q)R}
                 do ilo1 = 1,nLO_at(l1+1,iat)
                    lc1 = lc10 + ilo1
                    ccx(ilo1,lm1,im) = ccx(ilo1,lm1,im) + phs_angint * s3r(irm,ic,lc1,iat) ! cc = <u_{irm,lb}| uc_{a1,l2} ulo_{a1,l1}> <Y_{lb,mb}|Y_{l2,m2} Y*_{l1,m1}>*e^{i(k-q)R}
                 enddo
              endif
           enddo
        enddo
     enddo

     
     if (.not. orig) then
        czero = cmplx(0.d0,0.d0,8)
        cone  = cmplx(1.d0,0.d0,8)
        lmx = (lmax+1)**2
        
        call zgemm('n','n', nb,loctmatsize,lmx, cone, alfp((mstart+1),1,idf),mbmax, aax,lmx, czero, res, nb)
        call zgemm('n','n', nb,loctmatsize,lmx, cone, betp((mstart+1),1,idf),mbmax, bbx,lmx, cone,  res, nb)
        if (mbmax_g == mbmax) then ! there are also local orbitals
           call zgemm('n','n', nb,loctmatsize,lmx*nLOmax, cone, gamp((mstart+1),1,1,idf),mbmax, ccx,lmx*nLOmax, cone, res, nb)
        endif
        mmat(icg-cstart,:,:) = res(:,:)
     endif
  enddo
  deallocate( locmixind3 )
  if (.not. orig) then
     deallocate( aax, bbx, ccx )
     deallocate( res )
  endif

  if (.False.) then
     print *, 'calc_minc'
     do im=1,loctmatsize
        do ie2=1,mend-mstart
           do icg=1,cend-cstart
              write(6,'(I3,1x,I3,1x,I3,1x,2F14.10)') im, ie2+mstart, icg+cstart, mmat(icg,ie2,im)  ! icg==ie2
           enddo
        enddo
     enddo
  endif
end subroutine calc_minc3
