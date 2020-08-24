
subroutine Ewald_Summation(sgm, pq, lam_max, pos, br2, rbas, alat, Vol, ortho, r_cf, g_cf, eta,  ndf)
  ! coul_strcnst
  ! Ewald sumation for the Coulomb repulsion
  implicit none
  complex*16, intent(out):: sgm((lam_max+1)*(lam_max+1),(ndf+1)*ndf/2)
  real(8),    intent(in) :: pq(3)         ! q in cartesian coordinates
  real(8),    intent(in) :: r_cf, g_cf, eta ! cutoffs in real and momentum space, eta==Ewalds parameter.
  real(8),    intent(in) :: pos(3,ndf)    ! atom's positions
  real(8),    intent(in) :: alat(3)       ! (a,b,c) size of the unit cell
  real(8),    intent(in) :: br2(3,3), rbas(3,3) ! reciprocal and real space unit cell vectors
  real(8),    intent(in) :: Vol           ! volume of the unit cell
  logical,    intent(in) :: ortho         ! Is the unit cell orthogoal or not
  integer,    intent(in) :: ndf, lam_max  ! num. atoms, maximum l-used
  !
  interface
     integer Function genr_str_size(Rmax,rshift,rbs)
       ! Finds number of lattice vectors to be included in the calculation of the structure constants. We use cutoff Rmax.
       real(8), intent(in)     :: Rmax      ! Maximum radius to search 
       real(8), intent(in)     :: rshift(3) ! Shift of the origin
       real(8), intent(in)     :: rbs(3,3)  ! Bravais lattice basis
     end Function genr_str_size
  end interface
  !
  real(8),    parameter  :: pi = 3.14159265358979d0
  real(8),    allocatable:: rstr(:,:)
  integer,    allocatable:: indx(:)
  complex*16, allocatable:: ylam(:), stmp1(:), stmp2(:), cm(:)
  complex*16 :: exp_q_aa, ttmp, exp_g_aa, pref_gausg_exp_g_aa, gtolam, imag
  real(8)    :: rleng, gleng, gausr, gausg, gammaor, pref, q_x_aa, g_x_aa, gamlam
  real(8)    :: rbs(3,3), r_aa(3), rp_aa(3), mq(3), gpq_aa(3), g(3)
  integer    :: np, ng, idf, jdf, i1, ijdf, ilm, l, m
  !
  imag = cmplx(0.d0, 1.d0, 8)
  allocate( ylam((lam_max+1)**2), stmp1((lam_max+1)**2), stmp2((lam_max+1)**2), cm((lam_max+1)**2) )
  rbs=transpose(rbas)
  mq(:) = -pq(:)
  do idf=1,ndf
     do jdf=idf,ndf
        if(ortho) then
           r_aa(:) = (pos(:,idf)-pos(:,jdf)) * alat(:)   ! vector between two different atoms in cartesian coordinates
        else
           r_aa = matmul(pos(:,idf)-pos(:,jdf), rbas)    ! vector between two different atoms in cartesian coordinates
        endif
        !! Calculate all real space lattice vectors R, such that R+r_aa < r_cf
        np = genr_str_size(r_cf, r_aa, rbs)            ! np is how many were found, looking for |r_aa+R| < r_cf
        allocate( rstr(4,np), indx(np) )                 ! space for these vectors, and its sortings
        call genr_str(rstr, indx, np, r_cf, r_aa, rbs) ! now find all of these vectors R+r_aa

        !write(6,'(A,I4,2x,A)') 'np=', np, 'rbs='
        !do i=1,3
        !   write(6,'(10x,3F10.6)') rbs(i,:)
        !enddo
        
        stmp1 = 0.d0 ! Initialize the temporal storage of the lattice sums
        !! Calculate the sum over the real space lattice
        do i1=1,np
           rp_aa(1:3) = rstr(1:3, indx(i1))
           rleng      = rstr(4,   indx(i1))
           
           !! calculate the values of the spherical harmonics at rp_aa
           call ylm(ylam, rp_aa, lam_max)
           q_x_aa   = dot_product(rp_aa, mq)
           exp_q_aa = cmplx (cos(q_x_aa), sin(q_x_aa), 8)
           gausr   = exp(-rleng**2/eta**2)
           gammaor = sqrt(pi)*erfc(rleng/eta)/rleng
           
           !write(6,'(A,I4,F12.5,5x,3F12.5,5x,2F12.5)') 'cc', i1, rleng, rp_aa(:), gausr, gammaor
           
           cm(1) = gammaor * exp_q_aa * ylam(1)
           do l = 1,lam_max
              gammaor = gammaor * (l-0.5d0)/rleng + rleng**(l-2) * gausr / ( eta**(2*l-1) )
              ttmp = gammaor * exp_q_aa
              do m = -l,l
                 ilm = l*l + l + m + 1
                 cm(ilm) =  ttmp * ylam(ilm)
              enddo ! m
           enddo ! l
           stmp1(:) = stmp1(:) + cm(:)
           !if (i1 .eq. np) write(6,'(A,I4,1x,F10.4,2x,6F12.7)') 'yy', i1, rleng, stmp1(17), stmp1(21), stmp1(25)
        enddo ! i1
        !! calculate the vectors for the sum in reciprocal space
        ng = genr_str_size(g_cf, pq, br2)             ! looking for |q+G| < g_cf
        deallocate( rstr, indx )
        allocate( rstr(4,ng), indx(ng) )
        call genr_str(rstr, indx, ng, g_cf, pq, br2)  ! generates reciprocal G for |q+G| < g_cf
        
        stmp2 = 0.d0
        !! Calculate the reciprocal lattice sum
        pref = 4.0d0*pi * sqrt(pi) / Vol
        do i1=1,ng
           gpq_aa(1:3) = rstr(1:3, indx(i1))  ! this is q+G
           g(1:3)      = gpq_aa(1:3) - pq(1:3)! this is G
           gleng       = rstr(4, indx(i1))    ! |q+G|
           !!calculate the values of the spherical harmonics at gpq_aa
           call ylm(ylam, gpq_aa,lam_max)     ! Y_lm(G+q)
           g_x_aa = dot_product( g, r_aa )     ! G * r_aa
           exp_g_aa = cmplx( cos(g_x_aa), sin(g_x_aa), 8)
           !
           gausg = exp( -(eta * gleng)**2/4.d0 )
           pref_gausg_exp_g_aa = pref * gausg * exp_g_aa
           !
           gtolam = 1.0d0 / gleng**2

           !write(6,'(A,I4,F12.5,5x,3F12.5,5x,2F12.5)') 'zzc', i1, gleng, gausg, g_x_aa, dble(gtolam), exp_g_aa
           
           cm(1) = pref_gausg_exp_g_aa * gtolam * ylam(1)
           do l = 1,lam_max
              gtolam = gtolam * (-gleng * imag/ 2.0d0)
              ttmp = pref_gausg_exp_g_aa * gtolam
              do m = -l,l
                 ilm = l*l + l + m + 1
                 cm(ilm) = ttmp * ylam(ilm)
              enddo ! mu
           enddo ! l
           stmp2(:) = stmp2(:) + cm(:)

           !if (i1 .eq. ng) write(6,'(A,I4,1x,F10.4,2x,6F12.7)') 'zz', i1, gleng, stmp2(17), stmp2(21), stmp2(25)
        enddo !i1
        deallocate( rstr, indx )
        !
        ijdf = idf + jdf*(jdf-1)/2 ! combined index for (idf,jdf). Notice idf <= jdf
        gamlam = sqrt(pi)
        sgm(1,ijdf) = ( stmp1(1) + stmp2(1) ) / gamlam

        do l = 1,lam_max
           gamlam = gamlam * (l-0.5d0)
           do m = -l,l
              ilm = l*l + l + m + 1
              sgm(ilm,ijdf) = ( stmp1(ilm) + stmp2(ilm) ) / gamlam
           enddo
        enddo
        if(idf .eq. jdf) sgm(1,ijdf) = sgm(1,ijdf) - 1.0d0/(eta*pi)

        gamlam = sqrt(pi)
        !write(6,'(A,2F12.7,2x,2F12.7,2x,2F12.7)') 'xx', stmp1(1)/gamlam, stmp2(1)/gamlam, sgm(1,1)
        
     enddo ! jdf
  enddo ! idf
  deallocate( ylam, stmp1, stmp2, cm )
  if(.False.) then
     write(6,*) 'sgm='
     do idf=1,ndf
        do jdf=idf,ndf
           ijdf=idf+jdf*(jdf-1)/2
           do l=0,lam_max
              do m=-l,l
                 ilm = l*l + l + m + 1
                 if (abs(sgm(ilm,ijdf)) > 1e-5) then
                    if (abs(aimag(sgm(ilm,ijdf))) > 1e-10) then
                       write(6,'(4i5,2e12.4)') m,l,idf,jdf,sgm(ilm,ijdf)
                    else
                       write(6,'(4i5,e12.4)') m,l,idf,jdf,dble(sgm(ilm,ijdf))
                    endif
                 endif
              enddo
           enddo
        enddo
     enddo
  endif
  return
end subroutine Ewald_Summation


subroutine MT_Coulomb(vmat, vq, iat, loctmatsize, big_l, nmix, im_start, tilg, rrint, djmm, rtlij, sgm, mult, nat, ndf, max_nmix, dimdj, nijrm, tsize, lam_max)
  ! coul_setvm0
  implicit none
  integer,   intent(in) :: nat, ndf, max_nmix, tsize, lam_max, nijrm, dimdj, loctmatsize
  complex*16,intent(out):: Vmat(loctmatsize,loctmatsize)
  real(8),   intent(in) :: vq(3)              ! q in non-cartesian
  integer,   intent(in) :: iat
  integer,   intent(in) :: big_l(max_nmix,nat)! L for product basis
  integer,   intent(in) :: nmix(nat)          ! how many product functions per atom
  integer,   intent(in) :: im_start(nat+1)    ! where does im (index for product functions) start at this atom == idf
  real(8),   intent(in) :: tilg(tsize)        ! \tilde{g}_{lm,l'm'}=\sqrt{4\pi}(-1)^{l}\sqrt{\frac{(l+l'+m+m')!(l+l'-m-m')!}{(2l+1)(2l'+1)[2(l+l')+1]!(l+m)!(l-m)!(l'+m')!(l'-m')!}}
  real(8),   intent(in) :: rrint(nat,nijrm)   ! rrint[iat, (im,jm)] = <u_{jm,L}| (r_<)^L / (r_>)^{L+1} | u_{im,L} >
  complex*16,intent(in) :: djmm(dimdj,ndf)    ! Rotation matrices for spherical harmonics D^j_{mm′}
  real(8),   intent(in) :: rtlij(max_nmix,max_nmix,nat,nat) ! < r^L_{im} | u_{im,iat} > * <r^L_{jm}| u_{jm,jat}>
  complex*16,intent(in) :: sgm((lam_max+1)*(lam_max+1),(ndf+1)*ndf/2)
  integer,   intent(in) :: mult(nat)          ! how many atoms of this sort
  !
  interface
     real(8) function gettildeg(l1,l2,m1,m2,tilg,tsize)
       integer, intent(in) :: l1, l2, m1, m2
       integer, intent(in) :: tsize
       real(8), intent(in) :: tilg(tsize)
     end function gettildeg
     complex*16 function getdjmm(idf,l,m1,m2, djmm, ndf, dimdj)
       integer,    intent(in) :: ndf, dimdj
       complex*16, intent(in) :: djmm(dimdj,ndf)
       integer,    intent(in) :: idf ! index of atoms for which the rotation matrix is retrieved 
       integer,    intent(in) :: l, m1, m2        
     end function getdjmm
  end interface
  !
  real(8),    parameter  :: pi = 3.14159265358979d0
  logical :: q_is_zero
  integer :: ieq, jeq, maxbigl, im, jm, iirm, jjrm, ijrm, minu, ijdf, irm, jrm, it, jat, idf, jdf!, minu1
  integer :: l1, l2, m1, m2, m3, m4, lm34, iidf, jjdf
  real(8) :: fourpi, tg
  complex*16 :: dj3, dj4, stc, sum34
  fourpi = 4.d0*pi

  q_is_zero = .False.
  if ( sum(abs(vq)) < 1e-7 ) then ! q==0
     q_is_zero  = .True. 
  endif

  idf = 0
  do it=1,iat-1
     idf = idf + mult(it)
  enddo
  maxbigl = maxval(big_l)
  Vmat(:,:) = 0.d0
  
  im = im_start(iat)
  do ieq = 1, mult(iat)        !* Loop over equivalent atoms:
     idf = idf+1
     ijdf=(idf*(idf+1))/2      !* ijdf for jdf==idf
     do irm = 1, nmix(iat)     !* Loop over mixed functions:
        l1 = big_l(irm,iat)
        do m1 = -l1,l1
           im = im + 1
           jdf=0
           do jat=1,nat
              jm = im_start(jat)
              do jeq=1,mult(jat)
                 jdf=jdf+1
                 iidf = min(idf,jdf)
                 jjdf = max(idf,jdf)
                 ijdf = iidf + jjdf*(jjdf-1)/2  ! This is how we stored Ewalds sum
                 do jrm = 1, nmix(jat)
                    l2 = big_l(jrm,jat)
                    do m2 = -l2,l2
                       jm = jm + 1
                       if (.not. (q_is_zero .and. l1.eq.0 .and. l2.eq.0) ) then
                          sum34 = cmplx(0.d0,0.d0,8)
                          do m3 = -l1,l1    ! -m3 is m1 if atom1's coordinate axis coincides with the global axis. In general, there is extra rotation.
                             do m4 = -l2,l2 !  m4 is m2 if atom2's coordinate axis coincides with the global axis.
                                if ( idf <= jdf) then
                                   lm34 = (l1+l2+1)*(l1+l2) + m4 + m3 + 1  ! Index for the Ewald sumation, l=l1+l2, m=m3+m4
                                   stc = sgm(lm34, ijdf)       ! Ewald summation \sum_R Y_{l1+l2,m3+m4} e^{-iq*R}/R
                                   minu = 1-2*mod(abs(m3),2)   ! == (-1)**m3 : This minus is from the two center expansion, which has both Y_{l3,m3}^* and Y_{l4,m4}^* conjugated, and for delta function, we need Y_{l3,m3}=(-1)^m3 Y_{l3,-m3}
                                   !if (minu .ne. minu1) print *, 'ERROR', m3, minu,minu1
                                   !write(6,'(I4,1x,I4,1x,2F16.8)') ijdf, lm34, stc
                                else
                                   lm34 = (l1+l2+1)*(l1+l2) -m4 - m3 + 1  ! Here we need conjugated Ewald sum of e^{-iq*(-R)}/R Y_{l1+l2,m3+m4} = [e^{-iq*R}/R Y_{l1+l2,-m3-m4} ]^*  (-1)^{m3+m4}
                                   stc = conjg( sgm(lm34, ijdf) )         ! need conjugated Ewald sum
                                   minu = 1-2*mod(abs(m4+l1+l2),2)        ! == (-1)**(m4+l1+l2) : now that we have extra (-1)^(m3+m4) and extra (-1)^m3, we should have (-1)^m4 only
                                   !if (minu .ne. minu1) print *, 'ERROR', m4+l1+l2, minu, minu1
                                   !write(6,'(I4,1x,I4,1x,2F16.8)') ijdf, lm34, stc
                                endif
                                tg = gettildeg(l1,l2,m3,m4,tilg,tsize)
                                dj3 = conjg( getdjmm(idf, l1, -m3, m1, djmm, ndf, dimdj) )  ! rotating Y_{l1,m3}^* into global coordinate system Y_{l1,m1}. Se explanation for -m3 below. Comes from two center expansion.
                                dj4 =        getdjmm(jdf, l2,  m4, m2, djmm, ndf, dimdj)    ! rotating Y_{l2,m4} into global coordinate system Y_{l2,m2}.
                                sum34 = sum34 + minu * stc * tg * dj3 * dj4
                             enddo ! m4
                          enddo ! m3
                          ! rtlij = < r^L_{im} | u_{im,iat} > * <r^L_{jm}| u_{jm,jat}>
                          ! Vmat = rtlij * Ewald_sum * (-1)^m1 * tildeg(l1,l2,-m1,m2). We also take into account the rotation when global and local axis are not alligned.
                          Vmat(im,jm) = rtlij(jrm,irm,jat,iat) * sum34
                          !write(6,'(I2,1x,I3,1x,2I2,1x,I2,1x,I3,1x,2I2,1x,2F15.10,1x,F16.8)') idf, irm, l1, m1, jdf, jrm, l2, m2, sum34, rtlij(jrm,irm,jat,iat) 
                          if ( (m1.eq.m2) .and. (l1.eq.l2) .and. (idf.eq.jdf) ) then
                             ! Contribution for the same atom, we use standard one-center Laplace expansion
                             iirm = min(irm,jrm)
                             jjrm = max(irm,jrm)
                             ijrm = jjrm*(jjrm-1)/2 + iirm  ! this is how integrals are stored.
                             !  Vmt += 4*pi * <u_{jm,L,idf}| (r_<)^L / (r_>)^{L+1} | u_{im,L,idf} >/(2*l1+1)
                             Vmat(im,jm) = Vmat(im,jm) + fourpi * rrint(iat,ijrm)/dble(2*l1+1)
                             ! Here we do not need to care to rotate to the global coordinate systm, because the expansion is valid also in the local coordinate system.
                             ! When we peformed the Ewald sum, we have done it in the global system, while all atom-centered integrations are done in local coordinate system, hence the rotation was necessary.
                          endif
                       endif ! if((iq.ne.1).or.(l1.ne.0).or.(l2.ne.0))
                    enddo ! m2
                 enddo ! jrm
              enddo ! jeq
           enddo ! jat
        enddo ! m1  
     enddo ! irm
  enddo ! ieq
end subroutine MT_Coulomb


subroutine cmp_pw_olap(olap, indgq, ipwint, gindex, i_g0, ngq, ng, maxngq_barc, n1, n2, n3)
  ! Calculates overlap between plane waves from indgq
  implicit none
  complex*16, intent(out):: olap(ngq,ngq)
  integer,   intent(in)  :: indgq(maxngq_barc)  !
  integer,   intent(in)  :: gindex(ng,3)       ! all G-points
  complex*16, intent(in) :: ipwint(ng)
  integer, intent(in)    :: i_g0(2*n1+1,2*n2+1,2*n3+1)
  integer,    intent(in) :: ngq, ng, maxngq_barc, n1, n2, n3
  !
  integer :: ipw, jpw, iG(3), idg
  do ipw=1,ngq
     olap(ipw,ipw) = ipwint(1)
     do jpw=ipw+1,ngq
        iG  = gindex(indgq(ipw)+1,:) - gindex(indgq(jpw)+1,:)
        idg = i_g0(iG(1)+n1+1,iG(2)+n2+1,iG(3)+n3+1)
        olap(ipw,jpw) = ipwint(idg+1)
        olap(jpw,ipw) = conjg(ipwint(idg+1))
     enddo
  enddo
end subroutine cmp_pw_olap

subroutine cmp_pw_olap2(olap, indgq, ipwint, gindex, i_g0, ngq, ngq_barc, ng, maxngq_barc, n1, n2, n3)
  ! Calculates overlap between plane waves from indgq
  implicit none
  complex*16,intent(out) :: olap(ngq,ngq_barc)
  integer,   intent(in)  :: indgq(maxngq_barc)  !
  integer,   intent(in)  :: gindex(ng,3)       ! all G-points
  complex*16, intent(in) :: ipwint(ng)
  integer, intent(in)    :: i_g0(2*n1+1,2*n2+1,2*n3+1)
  integer,    intent(in) :: ngq, ngq_barc, ng, maxngq_barc, n1, n2, n3
  !
  integer :: ipw, jpw, iGv(3), idg
  olap = 0.d0
  do ipw=1,ngq_barc
     do jpw=1,ngq
        iGv  = gindex(indgq(ipw)+1,:) - gindex(indgq(jpw)+1,:)   ! G_{ipw}-G_{jpw}
        if ( abs(iGv(1))<=n1 .and. abs(iGv(2))<=n2 .and. abs(iGv(3))<=n3 ) idg = i_g0(iGv(1)+n1+1,iGv(2)+n2+1,iGv(3)+n3+1)+1
        if (idg>0) olap(jpw,ipw) = ipwint(idg)  ! olap[jpw,ipw]=Integrate[ (e^{i(G[ipw]-G[jpw])r},{r in interstitials}]/V_cell = <G_{jpw}|G_{ipw}>_{interstitials}
     enddo
  enddo
end subroutine cmp_pw_olap2
  
subroutine Mixed_Coulomb(vtemp, vq, iat, Q_Vq, jlam, gindex, indgq, gqlen, gqlen0, G_unique, ngq_barc, sloctmatsize, big_l, k2cartes, rotloc, trotij, mult, pos, Vol, nmix, nat, ng, ngs, ngs2, ndf, maxngq_barc)
  !   Calculates the matrix element between an atomic product function and a plane wave for the interstitial
  !                    <e^{(q+G)\vr}|V_{coul}|u_{irm}>  = 4*pi/|q+G|^2 e^{-i(q+G)_R }<q+G|u_{irm}>
  !   coul_setvm0
  implicit none
  integer,   intent(in) :: nat, ndf, ng, ngs, ngs2, sloctmatsize, maxngq_barc
  complex*16,intent(out):: vtemp(sloctmatsize,ngq_barc)
  real(8),   intent(in) :: vq(3)              ! q in non-cartesian
  integer,   intent(in) :: iat, ngq_barc
  logical,   intent(in) :: Q_Vq               ! Q_Vq==True : computing V_Mixed, Q_Vq==False : computing mpwmix
  integer,   intent(in) :: gindex(ng,3)       ! all G-points
  integer,   intent(in) :: indgq(maxngq_barc) !
  real(8),   intent(in) :: gqlen(maxngq_barc) ! notice that this is the length with repetitions
  real(8),   intent(in) :: gqlen0(ngs)        ! this contains unique lengths only
  integer,   intent(in) :: big_l(nmix), nmix  !
  real(8),   intent(in) :: jlam(nmix,ngs)     ! <j_L((q+G)r)|u_{mix,L}>=Int[j_L((q+G)r) u_{mix,L}/r r^2,r]
  integer,   intent(in) :: mult(nat)          !
  real(8),   intent(in) :: pos(3,ndf)         ! 
  real(8),   intent(in) :: k2cartes(3,3)      ! transform to cartesian
  real(8),   intent(in) :: rotloc(nat,3,3), trotij(ndf,3,3)
  real(8),   intent(in) :: Vol                ! Vol unit cell
  integer,   intent(in) :: G_unique(ngs2)     ! index of unique lengths
  !
  real(8),    parameter  :: pi = 3.14159265358979d0
  complex*16, allocatable :: sph(:,:)
  integer    :: ieq, idf, it, maxbigl, im, ipw0, irm, l1, m1, ilm_sph, ipw, ig, ig0, ig1
  real(8)    :: rotation1(3,3), rotation(3,3), Gq(3)
  real(8)    :: fourpi, prefac
  complex*16 :: imag, two_pi_i, cval, expqr, prefac2
  !
  fourpi = 4.d0*pi
  imag = cmplx(0.d0, 1.d0, 8)
  two_pi_i = 2.0d0*pi*imag
  maxbigl = maxval(big_l)
  allocate( sph( (maxbigl+1)*(maxbigl+1), ngq_barc) )
  vtemp(:,:) = 0.d0
  !
  ipw0 = 1
  if ( sum(abs(vq)) < 1e-7 ) then ! q==0
     ipw0=2 ! should not use G==0 when q==0, i.e., the diverging term
  endif
  prefac = fourpi/sqrt(Vol)
  !
  idf = 0
  do it=1,iat-1
     idf = idf + mult(it)
  enddo
  !
  im=0
  do ieq = 1, mult(iat)        ! Loop over equivalent atoms:
     idf = idf+1               ! Index for this atom in the long list of all atoms
     ! First computes atomic centered spherical harmonics Y_lm(q+G) for all G's
     rotation1 = matmul(transpose(trotij(idf,:,:)), k2cartes) ! total rotation on this atom is  rotloc*rotij*k2cartes (G+q)_{semi-cartesian}
     rotation  = matmul( rotloc(iat,:,:), rotation1 )
     do ig=1,ngq_barc
        Gq = matmul(rotation, vq(:) + gindex(indgq(ig)+1, :) )
        !write(6,'(A,I3,1x,A,I4,1x,A,3F10.6)') 'idf=', idf, 'ig=', ig, 'G+q=', Gq(:)
        call ylm(sph(:,ig), Gq, maxbigl)  ! Spherical harmonics on this atom, properly rotated to the global coordinate system
     enddo
     do irm = 1, nmix                           ! Loop over all atomic product functions on this atom
        l1 = big_l(irm)                         ! L for the product function
        do m1 = -l1,l1                          ! M for the product function
           im = im + 1                          ! index of the product function
           ilm_sph = l1*(l1+1) + m1 + 1         ! index of the spherical harmonics
           prefac2 = prefac * (-imag)**l1
           do ipw=1,ngs
              ! from ig0 to ig1 the length of the vector |q+G| is the same
              ig0 = max(G_unique(ipw)+1,ipw0)   ! when q==0, we remove G=0 divergent term
              if (ipw < ngs) then
                 ig1 = G_unique(ipw+1)          ! the entry at which this next length starts
              else
                 ig1 = ngq_barc                 ! or until the last entry
              endif
              if (Q_Vq) then
                 cval = (fourpi/gqlen0(ipw)**2)  * prefac2 * jlam(irm,ipw)  ! 4*pi/|q+G|^2 <q+G|u_{irm}>
              else
                 cval =                            prefac2 * jlam(irm,ipw)  !              <q+G|u_{irm}>
              endif
              do ig = ig0,ig1
                 expqr = exp( -two_pi_i * dot_product(pos(:,idf), gindex(indgq(ig)+1,:) ) ) ! e^{i (q+G) R_iat }
                 vtemp(im,ig) = sph(ilm_sph, ig) * expqr * cval   ! 4*pi/|q+G|^2 e^{-i(q+G)_R }<q+G|u_{irm}>
                 !write(6,'(I4,1x,I4,1x,I2,1x,I2,1x,I4,1x,2F10.5,1x,2F10.5,1x,2F16.12)') ig, im,  l1,m1,indgq(ig)+1,cval,sph(ilm_sph, ig),   vtemp(im,ig)
              enddo
           enddo ! ipw
        enddo ! m1  
     enddo ! irm
  enddo ! ieq
  deallocate( sph )
end subroutine Mixed_Coulomb


integer Function genr_str_size(Rmax,rshift,rbs)
  ! Finds number of lattice vectors to be included in the
  ! calculation of the structure constants. We use cutoff Rmax.
  implicit none
  real(8), intent(in)     :: Rmax      ! Maximum radius to search 
  real(8), intent(in)     :: rshift(3) ! Shift of the origin
  real(8), intent(in)     :: rbs(3,3)  ! Bravais lattice basis
  !
  integer(4) :: i,imax,ir,i1,i2,i3
  real(8) :: rleng              ! minimum length of the basis vectors.
  real(8) :: ivec(3), lrbs(3), rps(3)     ! vector belonging to the real space lattice and length of the basis vectors
  !
  do i=1,3
     lrbs(i) = sqrt(dot_product(rbs(:,i),rbs(:,i))) ! length of the three basis vectors
  enddo
  imax  = idint( Rmax/minval(lrbs) ) + 1          ! maximum needed index that we need to iterate
  ir=0
  do i1=-imax,imax
     do i2=-imax,imax
        do i3=-imax,imax
           ivec = (/i1,i2,i3/)
           rps(:) = matmul(rbs,ivec) + rshift(:)
           rleng = sqrt(dot_product(rps,rps))
           if ((rleng <= Rmax).and.(rleng > 1.0d-6)) ir=ir+1
        enddo
     enddo
  enddo
  genr_str_size = ir
  return
end Function genr_str_size

subroutine genr_str(rstr, indx, nr, Rmax, rshift, rbs)
  ! Generates the indexes of the lattice vectors to be included in the
  ! calculation of the structure constants. We use cutoff Rmax.
  implicit none
  real(8), intent(out)  :: rstr(4,nr)
  integer, intent(out)  :: indx(nr)
  integer, intent(in)   :: nr        !number of vectors
  real(8), intent(in)   :: Rmax      ! Maximum radius to search 
  real(8), intent(in)   :: rshift(3) ! Shift of the origin
  real(8), intent(in)   :: rbs(3,3)  ! Bravais lattice basis
  !
  integer(4) :: i,imax,ir,i1,i2,i3
  real(8) :: rleng              ! minimum length of the basis vectors.
  real(8) :: ivec(3), lrbs(3), rps(3)     ! vector belonging to the real space lattice and length of the basis vectors
  !
  do i=1,3
     lrbs(i) = sqrt(dot_product(rbs(:,i),rbs(:,i))) ! length of the three basis vectors
  enddo
  imax  = idint( Rmax/minval(lrbs) ) + 1          ! maximum needed index that we need to iterate
  ir=0
  do i1=-imax,imax
     do i2=-imax,imax
        do i3=-imax,imax
           ivec = (/i1,i2,i3/)
           rps(:) = matmul(rbs,ivec) + rshift(:)
           
           rleng = sqrt(dot_product(rps,rps))
           if ((rleng <= Rmax).and.(rleng > 1.0d-6)) then
              ir=ir+1
              if (ir > nr) then
                 write(6,*) 'ERROR genr_str in for_Coulomb.f90 : nr < ir. Call genr_str_size first and make sure there is enough space allocated.'
              endif
              rstr(1:3,ir) = rps(1:3)
              rstr(4,  ir) = rleng
           endif
        enddo
     enddo
  enddo
  if (.not.(nr==ir)) write(6,'(A,I4,A,I4,A)') 'ERROR genr_str in for_Coulomb.f90 : nr=',nr, 'ir=',ir,'. Call genr_str_size first and make sure there is enough space allocated.'

  indx = 0
  call qsort4(nr,rstr,indx)
end subroutine genr_str

subroutine pw_get_short_length(ngs, gind4, glen, nk)
  ! This function find how many unique lengths are present
  implicit none
  integer, intent(out)    :: ngs, gind4(nk)
  REAL(8), intent(in)     :: glen(nk)
  INTEGER, intent(in)     :: nk
  !
  integer    :: ipw, igl
  real(8)    :: glprev
  !
  glprev=0.d0 ! Notice that we remove G=0 by setting glprev to zero. Make it negative if we want to inlcude G=0
  igl=0
  do ipw=1,nk
     if ( abs(glen(ipw)-glprev) > 1e-10 ) then
        igl=igl+1          ! how many different lengths
        glprev = glen(ipw)
     endif
     gind4(ipw) = igl      ! gind4(ipw_long) = index_to_unique
  enddo
  ngs=igl
end subroutine pw_get_short_length

subroutine pw_get_short_and_phase(glen0, phase, ngs, gind4, gind, glen, pos, mult, nk, nat, ndf)
  ! This function calculated integral of the plane wave in the interstitials
  ! for all vectors in gind array with length glen:
  !  Integrate[ exp(i G\vr), {\vr in interstitials}]
  ! 
  ! We could use  Function int1ipw from for_vxcnn.f90, which does the same thing, but for
  ! a single vector
  implicit none
  real(8), intent(out)    :: glen0(ngs)
  complex*16, intent(out) :: phase(ndf,ndf,ngs)
  INTEGER, intent(in)     :: gind(nk,3), gind4(nk)
  REAL(8), intent(in)     :: glen(nk)
  INTEGER, intent(in)     :: mult(nat)
  REAL(8), intent(in)     :: pos(3,ndf)
  integer, intent(in)     :: ngs        ! length of the short arrays
  INTEGER, intent(in)     :: nk, nat, ndf
  !
  real(8),    parameter  :: pi = 3.14159265358979d0
  !
  integer    :: ipw, igl, jgl, iat, ieq, jat, jeq, idf, jdf
  real(8)    :: dpos(3)
  complex*16 :: expg, two_pi_i
  !
  glen0=0.d0    ! glen0 stores only different lengths  (much smaller than number of all G-points)
  jgl=0
  do ipw=1,nk
     !print *, 'ipw=', ipw, 'jgl=', jgl
     if (gind4(ipw) .ne. jgl) then
        jgl = gind4(ipw)
        glen0(jgl) = glen(ipw) ! fancy way of setting smaller mesh of glen0 with inequivalent lengths only
        !write(6,'(A,I3,A,F10.6)') 'glen0[',jgl,']=',glen0(jgl)
     endif
  enddo
  if (jgl.ne.ngs) then
     write(6,'(A,I4,A,I4)') 'ERROR in pw_get_short_and_phase for_Coulomb.f90 ngs=', ngs, 'while jgl=', jgl
  endif
  !
  two_pi_i = 2.0d0*pi*cmplx(0.d0,1.d0,8)
  phase(:,:,:) = cmplx(0.d0,0.d0,8)
  idf=0
  do iat = 1, nat
     do ieq = 1, mult(iat)
        idf=idf+1
        jdf=0    
        do jat=1,nat
           do jeq=1,mult(jat)
              jdf=jdf+1
              if (idf .eq. jdf) then
                 do ipw = 2, nk                 ! over all G's, excepy G=0
                    igl = gind4(ipw)            ! igl is index in the shorter array
                    phase(idf,jdf,igl) = phase(idf,jdf,igl) + 1.d0 ! phase for the entire set of G vectors which have equal length!
                 enddo ! ipw
              else
                 dpos(1:3) = pos(1:3,idf)-pos(1:3,jdf)
                 do ipw = 2, nk                 ! over all G's, excepy G=0
                    expg = exp( -dot_product( gind(ipw,:), dpos) * two_pi_i )  ! exp(-(G*(R_i-R_j)*2*pi*1j)
                    igl = gind4(ipw)            ! igl is index in the shorter array
                    phase(idf,jdf,igl) = phase(idf,jdf,igl) + expg      ! phase for the entire set of G vectors which have equal length!
                 enddo ! ipw
              endif
              !do igl=1,ngs
              !   write(6,'(3I5,2x,2F12.4)') idf, jdf, igl, phase(idf,jdf,igl)
              !end do
           enddo ! jeq
        enddo ! jat
     enddo ! ieq
  enddo ! iat
end subroutine pw_get_short_and_phase


!---------------------------------------------------------------
! Sort a single-precision real array by index, with a fuzzy equality test
subroutine qsort4(array_size,value4,index)
  integer, intent(in)    :: array_size
  integer, intent(inout) :: index(array_size)
  real(8), intent(in)    :: value4(4,array_size)
  integer, parameter     :: QSORT_THRESHOLD=8
  include "qsort_inline.inc"
contains
  ! set up initial index:
  subroutine init()
    integer :: i
    do i=1,array_size
       index(i)=i
    end do
  end subroutine init
  ! swap indices a,b
  subroutine swap(a,b)
    integer, intent(in) :: a,b
    integer :: hold
    hold=index(a)
    index(a)=index(b)
    index(b)=hold
  end subroutine swap
  ! circular shift-right by one:
  subroutine rshift(left,right)
    implicit none
    integer, intent(in) :: left, right
    integer :: hold, i
    hold=index(right)
    ! This sytnax is valid, but has poor optimization in GFortran:
    ! index(left+1:right)=index(left:right-1)
    do i=right,left+1,-1
      index(i)=index(i-1)
    end do
    index(left)=hold
  end subroutine rshift
  logical function less_than(a,b)
    integer, intent(in) :: a,b
    real(8), parameter :: small=1.0e-6
    if ( abs(value4(4,index(a))-value4(4,index(b))) < small ) then
       less_than = index(a) < index(b)
    else
       less_than = value4(4,index(a)) < value4(4,index(b))
    end if
  end function less_than
end subroutine qsort4

subroutine spher_bessel(sj,l,x,nx)
  !   Calculates the spherical bessel function $j_l(x)$ for many x(nx) and l<=lmax
  implicit none
  integer, intent(in)  :: l, nx      ! Maximum l to be calculated
  real(8), intent(in)  :: x(nx)      ! The argument at which the bessel functions are calculated
  real(8), intent(out) :: sj(l+1,nx) ! values of spherical bessel function
  !
  real(8) :: sjp(l+1)                ! values of first derivative of the spherical bessel function
  integer :: ix
  !
  do ix=1,nx
     call sphbes(l, x(ix), sj(:,ix), sjp)
  end do
end subroutine spher_bessel

real(8) function gettildeg(l1,l2,m1,m2,tilg,tsiz)
  integer, intent(in) :: l1, l2, m1, m2
  integer, intent(in) :: tsiz
  real(8), intent(in) :: tilg(tsiz)
  !
  integer :: j1, j2, mj1, mj2
  integer :: index1, index2, index3, index4
  integer :: par, fact, tind
  par = mod(abs(l1-l2),2)
  !! set the indexes j1, j2, mj1, mj2, and the multiplication factor
  ! first assign (j1,j2)=(l1,l2) so that j1 >= j2. Notice that function
  ! is simmetric with respect to (l1,m1) and (l2,m2)
  if(l1 < l2) then
     j1 = l2
     j2 = l1
     mj1 = m2
     mj2 = m1
     fact = 1-2*par  !!set the multiplication factor to (-1)^(l1-l2)
  else
     j1 = l1
     j2 = l2
     mj1 = m1
     mj2 = m2
     fact = 1
  endif
  if(mj2 < 0) then
     mj1 = -mj1
     mj2 = -mj2
  endif
  !! get the position of the value searched in the vector tilg
  index1 = j2 + mj2 + 1
  index2 = (j2+1)*(mj1+j1)
  index3 = (j2+1)*j1*j2 + (j2-1)*j2/2
  index4 = (j1+1)*(j1+2)*(3*j1-1)*j1/12
  tind = index1 + index2 + index3 + index4
  gettildeg = fact * tilg(tind)
end function gettildeg

complex*16 function getdjmm(idf,l,m1,m2, djmm, ndf, dimdj)
  ! Extracts the value of $D^j_{m_1m_2}$ from the djmm vector by using:
  !\begin{equation}
  !D^j_{m_1m_2}=\left\{
  !\begin{array}{ll}
  !\texttt{djmm(i)} & m_2 \ge 0 \\
  !(-1)^{m_1-m_2}\texttt{djmm(i)}^* & m_2 < 0
  !\end{array}
  !\right.
  !\end{equation}
  !where the index is given by
  !\begin{equation}
  !\texttt{i}=\tfrac{1}{6}l(l+1)(4j-1)+(l+1)(\textrm{sgn}(m_2)m_1+l)+|m_2|+1
  !\end{equation}
  !use rotylm, only: djmm
  implicit none
  integer,    intent(in) :: ndf, dimdj
  complex*16, intent(in) :: djmm(dimdj,ndf)
  integer,    intent(in) :: idf ! index of atoms for which the rotation matrix is retrieved 
  integer,    intent(in) :: l, m1, m2        
  !
  integer :: i,par
  !
  intrinsic iabs
  intrinsic sign
  !     first check that the input parameters are physical
  if ( (iabs(m2) <= l) .and. (iabs(m1) <= l) ) then
     !       calculate the index of the djmm vector
     i = l*(l+1)*(4*l-1)/6+(l+1)*(sign(1,m2)*m1+l)+iabs(m2)+1
     if(m2 >= 0)then
        getdjmm = djmm(i,idf)  !  For m2>=0 D^j_m1m2=djmm(i)
     else
        !         For m2<0 D^j_m1m2=(-1)^(m1-m2)djmm^*(i)
        par = mod(abs(m1-m2),2)
        getdjmm = conjg(djmm(i,idf)) * (1-2*par) ! (-1)^(m1-m2)* conj()
     endif
     return
  else
     !       Error handling
     write(6,*) "unphysical values as input",'l =',l,'m1 =',m1,'m2 =',m2  
     call exit(1)
  endif
end function getdjmm

