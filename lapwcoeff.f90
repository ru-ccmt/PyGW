!subroutine gap2_set_lapwcoef(alm,blm,clm, kl, gindex, indgk, abcelo, rmt, pos, mult, umt, &
!     & rotloc, trotij, Vol, iop, ifrot, k2cartes, nLO_at, nlo, lapw, nLOmax, nat, nvk, nt, ngk, ng, ndf, lomax, nloat)
subroutine gap2_set_lapwcoef(alm,blm,clm, kl, indgk, iop, ifrot, nvk, &
  & gindex, abcelo, rmt, pos, mult, umt, rotloc, trotij, Vol, k2cartes, nLO_at, nlo, lapw, nLOmax, nat, nt, ngk, ng, ndf, lomax, nloat)
  
  ! This subroutine calculates the coefficients Alm, Blm and Clm of the
  !augmentation for LAPW wave functions, according to the formulas:
  !
  !\begin{subequations}
  !\begin{align}
  !A^{\alpha}_{lm}(\vec{k}+\vec{G})=&\frac{4\pi}{\sqrt{\Omega}}R_{MT_{\alpha}}^2 %
  !   i^l Y^*_{lm}(\hat{k+G})e^{i(\vec{k}+\vec{G})\cdot\vec{r}_{\alpha}}a_l \\    
  !B^{\alpha}_{lm}(\vec{k}+\vec{G})=&\frac{4\pi}{\sqrt{\Omega}}R_{MT_{\alpha}}^2 %
  !   i^l Y^*_{lm}(\hat{k+G})e^{i(\vec{k}+\vec{G})\cdot\vec{r}_{\alpha}}b_l \\    
  !C^{\alpha}_{lm}(\vec{k}+\vec{G})=&\frac{4\pi}{\sqrt{\Omega}}R_{MT_{\alpha}}^2 %
  !   i^l Y^*_{lm}(\hat{k+G})e^{i(\vec{k}+\vec{G})\cdot\vec{r}_{\alpha}}c_l     
  !\end{align}
  !\end{subequations}
  !use bands,     only: nv
  !use eigenvec,  only: alm,blm,clm
  !use kpoints,   only: idvk, klist, kpirind
  !use lapwlo,    only: nt,lapw,nlo_tot,abcelo,umt,lmax,lomax,nlo_at,nlomax
  !use struk,     only: mult, pos, rmt, vi, nat,rotij,rotloc
  !use recipvec,  only: gindex, indgk, rk
  !use task,      only: casename
  ! !INPUT PARAMETERS:
  implicit none
  complex*16, intent(out) :: alm(ngk,nt*nt,ndf)
  complex*16, intent(out) :: blm(ngk,nt*nt,ndf)
  complex*16, intent(out) :: clm(ngk,nLOmax,nt*nt,ndf)
  !
  real(8), intent(in) :: kl(3)                 ! k-point in semi-cartesian form
  integer, intent(in) :: gindex(ng,3)          ! all G-points 
  integer, intent(in) :: indgk(ngk)
  real(8), intent(in) :: rmt(nat)
  real(8), intent(in) :: umt(nat,nt,5)  ! umt[iat,lmax,:]=[p,pe,dp,dpe,pei]
  integer, intent(in) :: iop   !  =1 --- alfa, 2 --- alfap
  logical, intent(in) :: ifrot ! true if local rotations are taken  into account   
  real(8), intent(in) :: k2cartes(3,3)
  integer, intent(in) :: mult(nat) 
  real(8), intent(in) :: abcelo(nat,lomax+1,nloat,3)
  real(8), intent(in) :: Vol, pos(3,ndf)
  real(8), intent(in) :: rotloc(nat,3,3)
  real(8), intent(in) :: trotij(ndf,3,3)
  integer, intent(in) :: nLO_at(lomax+1,nat), nlo(lomax+1,nat), lapw(nt,nat)
  integer, intent(in) :: nat, nvk, nt, ngk, ng, ndf, nLOmax, lomax, nloat
  !
  real(8), allocatable    :: al(:,:,:)   ! coefficient al(kn) for atom i
  real(8), allocatable    :: bl(:,:,:)   ! coefficient bl(kn) for atom i
  real(8), allocatable    :: kpg(:,:)
  complex*16, allocatable :: yl(:,:)     ! temporary storage of   spherical harmonics for local orbitals
  complex*16, allocatable :: phs(:)
  !
  real(8),    parameter :: pi = 3.14159265358979d0
  complex*16, parameter :: imag = (0.0d+0,1.0d+0) ! Imaginary unit.
  complex*16, parameter :: czero = 0.0d+0         !  complex zero.
  !
  integer :: n, ig, iat, l,m, idf, ieq, ilm, jlo,ilo0, i1, i2, lmax
  integer :: igvec(3), iGc(3)
  real(8) :: rk(ngk)
  real(8) :: rkn, rmt2, rhoatm, arg, p, pe, dp, dpe
  real(8) :: dfj(0:nt-1), fj(0:nt-1), Gc(3), kcar(3)
  real(8) :: vec(3), rotv1(3), rotv2(3)
  real(8) :: alo,blo,clo,two_pi
  complex*16:: cil, comc, cil_rhoatm
  logical, parameter :: ldbg=.false.
  !
  external lohns
  external sphbes
  external ylm
  lmax = nt-1
  two_pi = 2.0d+0*pi

  if(ldbg) write(6,*) "------------ start set_lapwcoef --------"
  
  ! change k-point to cartesian coordinates
  kcar = matmul(k2cartes, kl) ! This is just 2pi/a*Identity for orhogonal systems
  !
  allocate( al(nvk,0:nt,nat), bl(nvk,0:nt,nat) )
  allocate( yl((nt+1)*(nt+1),ngk) )
  allocate( phs(ngk), kpg(1:3,1:ngk) )
  al=0.d0
  bl=0.d0
  yl=0.d0
  phs=0.d0
  kpg(1:3,1:ngk) = 0.0d0
  
  !write(6,'(A,3F12.7,1x,A,3F12.7)') 'kvec=', kcar(:), 'kl=', kl(:)
  rk(1:ngk) = 0.0d0
  do ig=1,ngk
     iGc = gindex(indgk(ig)+1,:)
     Gc = matmul(k2cartes,iGc)   ! This is just 2pi/a*Idenity for orthogonal systems
     kpg(:,ig) = kcar + Gc
     rk(ig) = sqrt( sum(kpg(1:3,ig)**2) )
     !write(6,'(A,I3,1x,A,3I3,1x,A,3F12.7,1x,A,F10.7)') 'ig=', ig, 'iG=', iGc(:), 'G=', Gc(:), '|k+G|=', rk(ig)
  enddo

  do n = 1,nvk
     rkn = rk(n)
     do iat = 1, nat
        rmt2 = rmt(iat)**2
        call sphbes(nt-1,rmt(iat)*rkn,fj,dfj)
        do l = 0,lmax 
           p   = umt(iat,l+1,1) ! p
           pe  = umt(iat,l+1,2) ! pe
           dp  = umt(iat,l+1,3) ! dp
           dpe = umt(iat,l+1,4) ! dpe
           if ( lapw(l+1,iat) == 1) then
              al(n,l,iat) = rkn*dfj(l)*pe - fj(l)*dpe
              bl(n,l,iat) = fj(l)*dp - rkn*dfj(l)*p
           else
              al(n,l,iat) = fj(l)/p/rmt2
              bl(n,l,iat) = 0.d0
           endif
           !if (ldbg) write(6,'(A,I4,1x,A,I2,1x,A,I3,1x,A,F16.8,1x,A,F16.8,1x,A,2F16.8)') 'n=', n, 'iat=', iat, 'l=', l, 'al=', al(n,l,iat), 'bl=', bl(n,l,iat), 'fj=', fj, dfj
           !if (ldbg) write(6,'(A,I4,1x,A,I2,1x)') 'n=', n, 'iat=', iat
           !if (ldbg) write(6,'(A,I4,1x,A,I2,1x,A,I3,1x)') 'n=', n, 'iat=', iat, 'l=', l
           !if (ldbg) write(6,'(A,I4,1x,A,I2,1x,A,I3,1x,A,F16.8)') 'n=', n, 'iat=', iat, 'l=', l, 'al=', al(n,l,iat)
           !if (ldbg) write(6,'(A,I4,1x,A,I2,1x,A,I3,1x,A,F16.8,1x,A,F16.8,1x)') 'n=', n, 'iat=', iat, 'l=', l, 'al=', al(n,l,iat), 'bl=', bl(n,l,iat)
           if (ldbg) write(6,'(A,I4,1x,A,I2,1x,A,I3,1x,A,F16.8,1x,A,F16.8,1x,A)') 'n=', n, 'iat=', iat, 'l=', l, 'al=', al(n,l,iat), 'bl=', bl(n,l,iat), 'fj='
           !if (ldbg) write(6,'(A,I4,1x,A,I2,1x,A,I3,1x,A,F16.8,1x,A,F16.8,1x,A,F16.8)') 'n=', n, 'iat=', iat, 'l=', l, 'al=', al(n,l,iat), 'bl=', bl(n,l,iat), 'fj=', fj
           !if (ldbg) write(6,'(A,I4,1x,A,I2,1x,A,I3,1x,A,F16.8,1x,A,F16.8,1x,A,8F16.8)') 'n=', n, 'iat=', iat, 'l=', l, 'al=', al(n,l,iat), 'bl=', bl(n,l,iat), 'fj=', fj(:4), dfj(:4)
        enddo
     enddo
  enddo

  if (ldbg) then
     print *, 'kvec=', kl
     print *, 'pos=', pos
     idf = 0
     do iat = 1, nat
        print *, 'rotloc=', rotloc(iat,:,:)
        do ieq = 1, mult(iat)
           idf = idf + 1
           print *, 'rotij =', transpose(trotij(idf,:,:))
        end do
     end do
  end if
  
  alm=0.d0
  blm=0.d0
  clm=0.d0
  idf = 0
  do iat = 1, nat
     rhoatm = 4.0d+0*pi*sqrt(1.d0/Vol)*rmt(iat)*rmt(iat)
     do ieq = 1, mult(iat)
        idf = idf + 1
        do ig = 1, ngk
           igvec(1:3) = gindex(indgk(ig)+1,:)
           arg = two_pi * dot_product(igvec(:)+kl(:), pos(:,idf))
           phs(ig) = cmplx(dcos(arg),dsin(arg),8)
           vec(1:3) = kpg(1:3,ig)
           if (ifrot) then
              rotv1 = matmul(vec,trotij(idf,:,:))
              rotv2 = matmul(rotloc(iat,:,:), rotv1)
              call ylm(yl(1,ig),rotv2,nt)
           else
              call ylm(yl(1,ig),vec,nt)
           endif
           if (ldbg) write(6,'(3I3,1x,F12.8,3x,2F12.8,3x,3F12.8,3x,6F10.6)') igvec(:), arg/two_pi, phs(ig), kpg(:,ig), yl(:3,ig)
        enddo
        cil = (1.0d+0,0.0d+0)
        ilm = 0
        do l = 0, lmax
           cil_rhoatm = cil * rhoatm
           do m = -l, l
              ilm = ilm + 1
              do ig = 1, nvk
                 comc = cil_rhoatm * phs(ig) * conjg(yl(ilm,ig))
                 if (iop.ne.1) comc = conjg(comc)
                 alm(ig,ilm,idf) = al(ig,l,iat) * comc
                 blm(ig,ilm,idf) = bl(ig,l,iat) * comc
              enddo ! ig
              do ig = nvk+1,ngk ! over nlo_tot
                 alm(ig,ilm,idf) = 0.d0
                 blm(ig,ilm,idf) = 0.d0
                 if (nLOmax > 0 .and. l<=lomax) clm(ig,:,ilm,idf) = 0.d0
              enddo ! ig 
              if( l > lomax) cycle 
              ilo0 = lapw(l+1,iat)
              do jlo = lapw(l+1,iat), nLO_at(l+1,iat)
                 call lohns (i1,i2, iat,l,jlo-ilo0+1,  mult,nlo,lomax,nat)
                 alo = abcelo(iat,l+1,jlo+1,1)
                 blo = abcelo(iat,l+1,jlo+1,2)
                 clo = abcelo(iat,l+1,jlo+1,3)
                 if(ldbg) then
                    write(6,'(A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,F10.6,1x,A,F10.6,1x,A,F10.6,1x)') 'iat=', iat, 'l=', l,&
                         & 'ilm=', ilm,'jlo=', jlo, 'alo=', alo, 'blo=', blo, 'clo=', clo
                 endif
                 do ig = nvk+i1,nvk+i2
                    comc = cil_rhoatm * phs(ig) * conjg(yl(ilm,ig))
                    if(iop.ne.1) comc = conjg(comc)
                    alm(ig,ilm,idf) = alo * comc
                    blm(ig,ilm,idf) = blo * comc
                    if(jlo >= 1) clm(ig,jlo,ilm,idf) = clo * comc
                    if (ldbg) then 
                       if(jlo.eq.0) then 
                          write(6,'(A,i6,i6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6)') 'a,b', ig,ilm,alm(ig,ilm,idf),&
                               & blm(ig,ilm,idf)
                       else
                          write(6,'(A,i6,i6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6)') 'abc', ig,ilm,alm(ig,ilm,idf),&
                               & blm(ig,ilm,idf),clm(ig,jlo,ilm,idf) 
                       endif
                    endif
                 enddo ! ig
              enddo ! jlo
           enddo ! m
           cil = cil*imag
        enddo ! l
     enddo ! ieq 
  enddo !iat
  deallocate(al,bl,yl,phs,kpg)
  !write(6,*) 'almblm', 'kl=', kl
  !do idf=1,ndf
  !   do ig=1,ngk
  !      do ilm=1,(lmax+1)**2
  !         write(6,'(A,I2,1x,A,I2,1x,A,I3,1x,A,I3,1x,A,2F14.9,1x,A,2F14.9)') 'irk=', 0, 'idf=', idf, 'i=', ig, 'j=', ilm, 'alm=', alm(ig,ilm,idf), 'blm=', blm(ig,ilm,idf)
  !      end do
  !   end do
  !enddo
  return
end subroutine gap2_set_lapwcoef

subroutine dmft1_set_lapwcoef(alm,blm,clm, debug, iop, ifrot, kl, kirr,timat_ik,tau_ik, indgk,  nvk, gindex, abcelo, rmt, pos, mult, umt, &
     & rotloc, trotij, tauij, Vol, k2cartes, nLO_at, nlo, lapw, nLOmax, &
     & nat, nt, ngk, ng, ndf, lomax, nloat)
  ! timat_ik should be transposed of what we have in Python
  !
  ! This subroutine calculates the coefficients Alm, Blm and Clm of the
  !augmentation for LAPW wave functions, according to the formulas:
  !
  !\begin{subequations}
  !\begin{align}
  !A^{\alpha}_{lm}(\vec{k}+\vec{G})=&\frac{4\pi}{\sqrt{\Omega}}R_{MT_{\alpha}}^2 %
  !   i^l Y^*_{lm}(\hat{k+G})e^{i(\vec{k}+\vec{G})\cdot\vec{r}_{\alpha}}a_l \\    
  !B^{\alpha}_{lm}(\vec{k}+\vec{G})=&\frac{4\pi}{\sqrt{\Omega}}R_{MT_{\alpha}}^2 %
  !   i^l Y^*_{lm}(\hat{k+G})e^{i(\vec{k}+\vec{G})\cdot\vec{r}_{\alpha}}b_l \\    
  !C^{\alpha}_{lm}(\vec{k}+\vec{G})=&\frac{4\pi}{\sqrt{\Omega}}R_{MT_{\alpha}}^2 %
  !   i^l Y^*_{lm}(\hat{k+G})e^{i(\vec{k}+\vec{G})\cdot\vec{r}_{\alpha}}c_l     
  !\end{align}
  !\end{subequations}
  !use bands,     only: nv
  !use eigenvec,  only: alm,blm,clm
  !use kpoints,   only: idvk, klist, kpirind
  !use lapwlo,    only: nt,lapw,nlo_tot,abcelo,umt,lmax,lomax,nlo_at,nlomax
  !use struk,     only: mult, pos, rmt, vi, nat,rotij,rotloc
  !use recipvec,  only: gindex, indgk, rk
  !use task,      only: casename
  ! !INPUT PARAMETERS:
  implicit none
  complex*16, intent(out) :: alm(ngk,nt*nt,ndf)
  complex*16, intent(out) :: blm(ngk,nt*nt,ndf)
  complex*16, intent(out) :: clm(ngk,nLOmax,nt*nt,ndf)
  !
  real(8), intent(in) :: kl(3), kirr(3) ! irreducible k-point in semi-cartesian form and the corresponding irreducible point
  integer, intent(in) :: gindex(ng,3)   ! all G-points 
  integer, intent(in) :: indgk(ngk)
  real(8), intent(in) :: rmt(nat)
  real(8), intent(in) :: umt(nat,nt,5)  ! umt[iat,lmax,:]=[p,pe,dp,dpe,pei]
  integer, intent(in) :: iop   !  =1 --- alfa, 2 --- alfap
  logical, intent(in) :: ifrot, debug ! true if local rotations are taken  into account   
  real(8), intent(in) :: k2cartes(3,3)
  integer, intent(in) :: mult(nat) 
  real(8), intent(in) :: abcelo(nat,lomax+1,nloat,3)
  real(8), intent(in) :: Vol, pos(3,ndf)
  real(8), intent(in) :: rotloc(nat,3,3)
  real(8), intent(in) :: trotij(ndf,3,3)
  integer, intent(in) :: nLO_at(lomax+1,nat), nlo(lomax+1,nat), lapw(nt,nat)
  integer, intent(in) :: timat_ik(3,3)
  real(8), intent(in) :: tau_ik(3), tauij(ndf,3)
  integer, intent(in) :: nat, nvk, nt, ngk, ng, ndf, nLOmax, lomax, nloat
  !
  real(8), allocatable    :: al(:,:,:)   ! coefficient al(kn) for atom i
  real(8), allocatable    :: bl(:,:,:)   ! coefficient bl(kn) for atom i
  complex*16, allocatable :: yl(:,:)     ! temporary storage of   spherical harmonics for local orbitals
  complex*16, allocatable :: phs(:)
  !
  real(8),    parameter :: pi = 3.14159265358979d0
  complex*16, parameter :: imag = (0.0d+0,1.0d+0) ! Imaginary unit.
  complex*16, parameter :: czero = 0.0d+0         !  complex zero.
  !
  integer :: n, ig, iat, l,m, idf, ieq, ilm, jlo,ilo0, i1, i2, lmax
  integer :: igvec(3)!, iGc(3)
  real(8) :: bk3(3), bkrot(3), pos_idf(3), rotij(3,3), rotloc_x_k2cartes(3,3), rotloc_x_k2cartes_x_rotij(3,3)!, rk(ngk)
  real(8) :: rkn, rmt2, rhoatm, arg, p, pe, dp, dpe
  real(8) :: dfj(0:nt-1), fj(0:nt-1)!, kcar(3), Gc(3)
  real(8) :: vec(3)!, rotv1(3), rotv2(3)
  real(8) :: alo,blo,clo,two_pi
  complex*16:: cil, comc, cil_rhoatm
  real(8) :: pos_first(3)
  logical, parameter :: ldbg=.false.
  !
  external lohns
  external sphbes
  external ylm
  lmax = nt-1
  two_pi = 2.0d+0*pi

  if(ldbg) write(6,*) "------------ start set_lapwcoef --------"
  
  ! change k-point to cartesian coordinates
  !kcar = matmul(k2cartes, kl) ! This is just 2pi/a*Identity for orhogonal systems
  !
  allocate( al(nvk,0:nt,nat), bl(nvk,0:nt,nat) )
  allocate( yl((nt+1)*(nt+1),ngk) )
  allocate( phs(ngk) )
  al=0.d0
  bl=0.d0
  yl=0.d0
  phs=0.d0

  do n = 1,nvk
     vec(:) = matmul(k2cartes, gindex(indgk(n)+1,:) + kirr(:))
     rkn = sqrt( sum( vec(1:3)**2 ) )
     do iat = 1, nat
        rmt2 = rmt(iat)**2
        call sphbes(nt-1,rmt(iat)*rkn,fj,dfj)
        do l = 0,lmax 
           p   = umt(iat,l+1,1) ! p
           pe  = umt(iat,l+1,2) ! pe
           dp  = umt(iat,l+1,3) ! dp
           dpe = umt(iat,l+1,4) ! dpe
           if ( lapw(l+1,iat) == 1) then
              al(n,l,iat) = rkn*dfj(l)*pe - fj(l)*dpe
              bl(n,l,iat) = fj(l)*dp - rkn*dfj(l)*p
           else
              al(n,l,iat) = fj(l)/p/rmt2
              bl(n,l,iat) = 0.d0
           endif
           !if (debug) write(6,'(A,I4,1x,A,I2,1x,A,I3,1x,A,F16.10,1x,A,F16.10,1x,A,2F16.10)') 'n=', n, 'iat=', iat, 'l=', l , 'al=', al(n,l,iat), 'bl=', bl(n,l,iat), 'fj=', fj(l), dfj(l)
        enddo
     enddo
  enddo
  alm=0.d0
  blm=0.d0
  clm=0.d0
  idf = 0
  do iat = 1, nat
     rhoatm = 4.0d+0*pi*sqrt(1.d0/Vol)*rmt(iat)*rmt(iat)
     rotloc_x_k2cartes = matmul(rotloc(iat,:,:), k2cartes)
     pos_first(:) = pos(:,idf+1)
     do ieq = 1, mult(iat)
        idf = idf + 1
        rotij = transpose(trotij(idf,:,:))
        rotloc_x_k2cartes_x_rotij = matmul( rotloc_x_k2cartes, rotij )
        pos_idf = matmul(pos_first, rotij) + tauij(idf,:)  ! instead of pos(:,idf) we use an equivalent expression, up to a lattice vector. This will respect symmetry exactly.
        do ig = 1, ngk
           igvec(1:3) = gindex(indgk(ig)+1,:)
           ! note that : kl ~ dot_product(timat_ik, kirr) up to an reciprocal vector
           bk3(:) = igvec(:)+kirr(:)
           bkrot(:) = matmul(timat_ik, bk3(:) )
           arg  =  two_pi * dot_product( bkrot(:), pos_idf ) + two_pi * dot_product(tau_ik(:), bk3(:) )
           ! old : arg = two_pi * dot_product(igvec(:)+kl(:), pos(:,idf))
           phs(ig) = cmplx(dcos(arg),dsin(arg),8)
           if (ifrot) then
              vec = matmul(rotloc_x_k2cartes_x_rotij, bkrot )
              call ylm(yl(1,ig),vec,nt)
           else
              vec(1:3) = matmul(k2cartes, igvec + kl)   ! This is just 2pi/a*Idenity for orthogonal systems
              call ylm(yl(1,ig),vec,nt)
           endif
           !if (debug) write(6,'(I4,1x,I4,1x,3I3,1x,F12.8,3x,2F12.8,3x,6F10.6)') ig, indgk(ig)+1, igvec(:), arg, phs(ig), yl(:3,ig)
        enddo
        cil = (1.0d+0,0.0d+0)
        ilm = 0
        do l = 0, lmax
           cil_rhoatm = cil * rhoatm
           do m = -l, l
              ilm = ilm + 1
              do ig = 1, nvk
                 comc = cil_rhoatm * phs(ig) * conjg(yl(ilm,ig))
                 if (iop.ne.1) comc = conjg(comc)
                 alm(ig,ilm,idf) = al(ig,l,iat) * comc
                 blm(ig,ilm,idf) = bl(ig,l,iat) * comc
                 !if (debug) then
                 !   write(6,'(A,i6,i6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6)') 'a0b', ig,ilm,alm(ig,ilm,idf),blm(ig,ilm,idf)
                 !endif
              enddo ! ig
              do ig = nvk+1,ngk ! over nlo_tot
                 alm(ig,ilm,idf) = 0.d0
                 blm(ig,ilm,idf) = 0.d0
                 if (nLOmax > 0 .and. l<=lomax) clm(ig,:,ilm,idf) = 0.d0
              enddo ! ig 
              if( l > lomax) cycle 
              ilo0 = lapw(l+1,iat)
              do jlo = lapw(l+1,iat), nLO_at(l+1,iat)
                 call lohns (i1,i2, iat,l,jlo-ilo0+1,  mult,nlo,lomax,nat)
                 alo = abcelo(iat,l+1,jlo+1,1)
                 blo = abcelo(iat,l+1,jlo+1,2)
                 clo = abcelo(iat,l+1,jlo+1,3)
                 !if(debug) then
                 !   write(6,'(A,I2,1x,A,I2,1x,A,I2,1x,A,I2,1x,A,F10.6,1x,A,F10.6,1x,A,F10.6,1x)') 'iat=', iat, 'l=', l,&
                 !        & 'ilm=', ilm,'jlo=', jlo, 'alo=', alo, 'blo=', blo, 'clo=', clo
                 !endif
                 do ig = nvk+i1,nvk+i2
                    comc = cil_rhoatm * phs(ig) * conjg(yl(ilm,ig))
                    if(iop.ne.1) comc = conjg(comc)
                    alm(ig,ilm,idf) = alo * comc
                    blm(ig,ilm,idf) = blo * comc
                    if(jlo >= 1) clm(ig,jlo,ilm,idf) = clo * comc
                    !if (debug) then 
                    !   if(jlo.eq.0) then 
                    !      write(6,'(A,i6,i6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6)') 'a,b', ig,ilm,alm(ig,ilm,idf),&
                    !           & blm(ig,ilm,idf)
                    !   else
                    !      write(6,'(A,i6,i6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6,f12.6)') 'abc', ig,ilm,alm(ig,ilm,idf),&
                    !           & blm(ig,ilm,idf),clm(ig,jlo,ilm,idf) 
                    !   endif
                    !endif
                 enddo ! ig
              enddo ! jlo
           enddo ! m
           cil = cil*imag
        enddo ! l
     enddo ! ieq 
  enddo !iat
  deallocate(al,bl,yl,phs)!,kpg)
  return
end subroutine dmft1_set_lapwcoef

subroutine lohns(i1,i2, iat,l0,jlo, mult,nlo,lomax,nat)
  ! lohns calculates indices as they appear in the Hamiltonian within lapw1 step.
  ! Can copy this routine from cdftKS or src_lapw1
  implicit none
  integer, intent(out):: i1, i2
  integer, intent(in) :: iat, l0
  integer, intent(in) :: mult(nat)
  integer, intent(in) :: nlo(lomax+1,nat), lomax, nat
  integer     :: i,l, jlo, jlo1
  i2 = 1
  do i = 1, iat - 1
     do l = 0, lomax
        do jlo1 = 1,nlo(l+1,i)
           i2 = i2 + (2*l+1)*mult(i)
        enddo
     enddo
  enddo
  do l = 0, l0-1
     do jlo1 = 1,nlo(l+1,iat)
        i2 = i2 + (2*l+1)*mult(iat)
     enddo
  enddo
  do jlo1 =1,jlo
     i1 = i2
     i2 = i2 + (2*l+1)*mult(iat)
  enddo
  i2 = i2 - 1
  return
end subroutine lohns


