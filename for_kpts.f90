
subroutine k2indx(indx, k, div)
  IMPLICIT NONE
  integer, intent(out):: indx
  integer, intent(in) :: k(3), div(3)
  indx = ( k(1)*div(2) + k(2) )*div(3) + k(3)
end subroutine k2indx

subroutine indx2k(k, indx, div)
  IMPLICIT NONE
  integer, intent(out):: k(3)
  integer, intent(in) :: indx, div(3)
  !
  integer :: i1 
  k(3) = mod(indx,  div(3))
  i1 = indx / div(3)
  k(2) = mod(i1, div(2))
  k(1) = i1 / div(2)
end subroutine indx2k

subroutine k_m_q_id(kqid, div, nkp)
  IMPLICIT NONE
  integer, intent(in) :: div(3), nkp
  integer, intent(out):: kqid(nkp,nkp)
  !
  integer :: i, iq, ik, kq, indx
  integer :: qr(3), kr(3), kmq(3)
  !
  if (nkp .ne. div(1)*div(2)*div(3)) then
     print *, 'ERROR nkp=', nkp, 'while div*div*div=', div(1)*div(2)*div(3)
  end if
  do iq=1,nkp
     call indx2k(qr, iq-1, div)
     do ik=1,nkp
        call indx2k(kr, ik-1, div)
        do i=1,3
           kq = kr(i)-qr(i)
           if (kq < 0) kq = kq + div(i)
           kmq(i) = mod(kq, div(i))
        enddo
        call k2indx(indx, kmq, div)
        kqid(ik,iq) = indx
     end do
  end do
end subroutine k_m_q_id

subroutine ngq_size(ngq, ngq_barc, ngqlen, q_c, k_c, maxlen_mb, maxlen_coul, nq, nk)
  IMPLICIT NONE
  INTEGER, intent(out) :: ngq(nq), ngq_barc(nq), ngqlen(nq)
  REAL*8, intent(in)   :: q_c(nq,3)
  REAL*8, intent(in)   :: k_c(nk,3)
  REAL*8, intent(in)   :: maxlen_mb, maxlen_coul
  INTEGER, intent(in)  :: nq, nk
  !
  REAL*8  :: akq, aG, aG_prev
  REAL*8  :: Gpq_len(nk)
  INTEGER :: indxc(nk)
  INTEGER :: iq, ik, n, i0, i
  REAL*8, parameter :: eps = 1e-6
  ngq = 0
  ngq_barc = 0
  ngqlen = 0
  do iq=1,nq
     n = 0
     do ik=1,nk
        akq = norm2(k_c(ik,:)+q_c(iq,:))
        if ( akq < maxlen_mb  ) ngq(iq) = ngq(iq) + 1
        if ( akq < maxlen_coul) then
           n = n + 1             ! n == ngq_barc(iq)
           Gpq_len(n) = akq      ! remember the length of this vector, so that we can sort them below
        endif
     end do
     ngq_barc(iq) = n
     indxc=0
     !call hpsort_eps_epw(n, Gpq_len, indxc, eps) ! sorting all possible lengths of |q+G|
     call qsort(n, Gpq_len, indxc) ! sorting all possible lengths of |q+G|
     ! Now finding out how many of |q+G| are different, i.e., unique
     aG_prev = -1.d0
     do i0=1,n
        i = indxc(i0)
        !aG = Gpq_len(i0)  ! this array is already sorted in hpsort
        aG = Gpq_len(i)  ! with qsort it is not yet sorted
        if ( abs(aG-aG_prev) > 1e-6 ) then
           ngqlen(iq) = ngqlen(iq) + 1   ! how many unique lengths of |q+G|
           aG_prev = aG
           !print *, 'iq=', iq, 'i0=', i0, 'i=', i, 'aG=', aG, 'ngqlen=', ngqlen(iq)
        endif
     end do
  end do
end subroutine ngq_size

subroutine ngq_sort(indgq, indgqlen, G_unique, q_c, k_c, maxlen_coul, ngqlen, maxngq_barc, maxngqlen, nq, nk)
  IMPLICIT NONE
  INTEGER, intent(out) :: indgq(nq, maxngq_barc)
  REAL*8,  intent(out) :: indgqlen(nq, maxngq_barc)
  INTEGER, intent(out) :: G_unique(nq, maxngqlen)
  REAL*8,  intent(in)  :: q_c(nq,3)
  REAL*8,  intent(in)  :: k_c(nk,3)
  REAL*8,  intent(in)  :: maxlen_coul
  INTEGER, intent(in)  :: ngqlen(nq)
  INTEGER, intent(in)  :: maxngqlen, maxngq_barc, nq, nk
  !
  REAL*8  :: Gpq_len(maxngq_barc)
  INTEGER :: which_ik(maxngq_barc), indxc(maxngq_barc)
  INTEGER :: iq, ik, n, i, i0, ii
  REAL*8  :: akq, aG_prev, aG
  REAL*8, parameter :: eps = 1e-6
  !
  G_unique(:,:)=0
  do iq=1,nq
     n = 0
     do ik=1,nk
        akq = norm2(k_c(ik,:)+q_c(iq,:))
        if ( akq < maxlen_coul) then
           n = n + 1            ! n == ngq_barc(iq)
           Gpq_len(n) = akq     ! remember the length of this vector, so that we can sort them below
           which_ik(n) = ik     ! remember which G is this
           !write(6,'(A,3F12.6,2x,3F12.6,2x,F14.8,3x,2I4)') 'G,q=', k_c(ik,:), q_c(iq,:), akq, n, ik
        endif
     end do
     indxc=0
     !call hpsort_eps_epw(n, Gpq_len, indxc, eps) ! sorting all possible lengths of |q+G|
     call qsort(n, Gpq_len, indxc) ! sorting all possible lengths of |q+G|
     aG_prev = -1.d0
     ii = 0
     do i0=1,n
        i = indxc(i0)
        !aG = Gpq_len(i0)   ! this array is already sorted in hpsort
        aG = Gpq_len(i)    ! with qsort it is not yet sorted
        indgq   (iq,i0) = which_ik(i)-1 ! index for python is smaller for 1 compared to fortran
        indgqlen(iq,i0) = aG
        if ( abs(aG-aG_prev) > 1e-6 ) then
           ii = ii + 1
           G_unique(iq,ii) = i0-1  ! index for python is smaller for 1 compared to fortran
           aG_prev = aG
        endif
     end do
     if (ngqlen(iq) .ne. ii) then
        WRITE(6,*) 'ERROR in ngq_sort: for_kpts.f90: ngqlen=', ngqlen(iq),'while ii=', ii, 'should be equal!'
     endif
  end do
end subroutine ngq_sort

subroutine pw_integ(ipwint, gind, glen, vmt, Rmt, pos, mult, nk, nat, ndf)
  ! This function calculated integral of the plane wave in the interstitials
  ! for all vectors in gind array with length glen:
  !  Integrate[ exp(i G\vr), {\vr in interstitials}]
  ! 
  ! We could use  Function int1ipw from for_vxcnn.f90, which does the same thing, but for
  ! a single vector
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: ipwint(nk)
  INTEGER, intent(in)     :: gind(nk,3)
  REAL*8, intent(in)      :: glen(nk)
  REAL*8, intent(in)      :: vmt(nat)
  REAL*8, intent(in)      :: Rmt(nat)
  INTEGER, intent(in)     :: mult(nat)
  REAL*8, intent(in)      :: pos(3,ndf)
  INTEGER, intent(in)     :: nk, nat, ndf
  !
  real(8),    parameter  :: pi = 3.14159265358979d0
  INTEGER    :: ippw, i, idf, iat, ki(3), ieq
  REAL*8     :: ak, x, j1, intmod!, pk
  COMPLEX*16 :: integc
  COMPLEX*16 :: phase, two_pi_i
  two_pi_i = 2.0d0*pi*cmplx(0.d0,1.d0,8)
  ipwint = 0
  do ippw=1,nk
     !i = indx(ippw)+1
     i = ippw
     ki = gind(i,:)  ! ki=(i0,i1,i2)
     ak = glen(i)    ! ak=|kc|
     !print *, ippw, i, ki, ak
     if (ak < 1e-10) then
        integc = 0
        do iat=1,nat
           integc = integc + vmt(iat)*mult(iat)
        end do
        ipwint(ippw) = 1.d0 - integc
     else
        integc = 0
        idf = 0
        do iat=1,nat
           phase = 0
           do ieq=1,mult(iat)
              idf = idf + 1
              !pk = dot_product(pos(:,idf),ki)*2.d0*pi
              !phase = phase + dcmplx( cos(pk), sin(pk) )
              phase = phase + exp(two_pi_i * dot_product(pos(:,idf),ki))
           end do
           x = ak * Rmt(iat)                    ! kr[iat] = |kc|*Rmt[iat]
           j1 = (sin(x)/x-cos(x))/x             ! spherical_bessel == (sin(x)/x-cos(x))/x
           intmod = vmt(iat) * 3*j1/x           ! intmod[iat] = j1*3*vmt/kr
           integc = integc + intmod * phase
        end do
        ipwint(ippw) = -integc
     end if
  end do
end subroutine pw_integ

subroutine Count_number_PW(npw2apw, ngindx, ng, kmax, gmax, br2)
  implicit none
  integer, intent(out) :: npw2apw, ngindx
  integer, intent(in)  :: ng(3)
  real*8,  intent(in)  :: br2(3,3)
  real*8,  intent(in)  :: kmax, gmax
  !
  real*8  :: ki(3), kc(3), kk
  integer :: i1, i2, i3
  npw2apw = 0
  ngindx = 0
  do i1=-ng(1),ng(1)
     do i2=-ng(2),ng(2)
        do i3=-ng(3),ng(3)
           ki(1) = i1
           ki(2) = i2
           ki(3) = i3
           kc = matmul(br2,ki)
           kk = sqrt(dot_product(kc,kc))
           if ( kk < 2.*kmax+4) npw2apw = npw2apw + 1
           if ( kk < gmax ) ngindx = ngindx + 1
        enddo
     enddo
  enddo
end subroutine Count_number_PW

subroutine Generate_PW(glen, G_c, gindex, ng, ngindx, gmax, pia, br2, ortho, is_CXZ, Q_sort)
  implicit none
  real*8,  intent(out) :: glen(ngindx)
  real*8,  intent(out) :: G_c(ngindx,3)
  integer, intent(out) :: gindex(ngindx,3)
  !
  integer, intent(in)  :: ngindx
  integer, intent(in)  :: ng(3)
  real*8,  intent(in)  :: br2(3,3)
  real*8,  intent(in)  :: gmax, pia
  logical, intent(in)  :: ortho, is_CXZ, Q_sort
  !
  real*8  :: ki(3), kc(3), kk
  integer :: i1, i2, i3, ig, iig
  integer, allocatable :: indx(:)
  real*8,  allocatable :: tGc(:,:), tglen(:)
  integer, allocatable :: tgindex(:,:)
  ig = 0
  do i1=-ng(1),ng(1)
     do i2=-ng(2),ng(2)
        do i3=-ng(3),ng(3)
           ki(1) = i1
           ki(2) = i2
           ki(3) = i3
           kc = matmul(br2,ki)
           kk = sqrt(dot_product(kc,kc))
           if ( kk < gmax ) then
              ig = ig + 1
              glen(ig) = kk
              G_c(ig,:) = kc(:)
              if (ortho) then
                 gindex(ig,1) = nint( kc(1)/pia +0.1 )
                 gindex(ig,2) = nint( kc(2)/pia +0.1 )
                 gindex(ig,3) = nint( kc(3)/pia +0.1 )
              else
                 if (is_CXZ) then
                    gindex(ig,:) = (/i1+i3,i2,i3-i1/)
                 else
                    gindex(ig,:) = (/i1,i2,i3/)
                 endif
              endif
           endif
        enddo
     enddo
  enddo
  if (ig .ne. ngindx) then
     write(6,*) 'ERROR in Generate_PW : for_kpts.f90 number of plane waves=', ig, 'should be ', ngindx
     call exit(1)
  endif
  if (Q_sort) then
     allocate( indx(ngindx), tglen(ngindx), tgindex(ngindx,3), tGc(ngindx,3) )
     indx=0
     call qsort(ngindx,glen,indx)
     do ig=1,ngindx
        iig = indx(ig)
        tglen(ig) = glen(iig)
        tgindex(ig,:) = gindex(iig,:)
        tGc(ig,:) = G_c(iig,:)
     enddo
     glen(:) = tglen(:)
     gindex(:,:) = tgindex(:,:)
     G_c(:,:) = tGc(:,:)
     deallocate(indx, tglen, tgindex, tGc)
  endif
end subroutine Generate_PW

!---------------------------------------------------------------
! Sort a single-precision real array by index, with a fuzzy equality test
subroutine qsort(array_size,value,index)
  integer, intent(in)    :: array_size
  integer, intent(inout) :: index(array_size)
  real(8), intent(in)    :: value(array_size)
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
    if ( abs(value(index(a))-value(index(b))) < small ) then
       less_than = index(a) < index(b)
    else
       less_than = value(index(a)) < value(index(b))
    end if
  end function less_than
end subroutine qsort

subroutine hpsort_eps_epw (n, ra, ind, eps)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  !use kinds, ONLY : DP
  implicit none  
  integer, parameter:: dp=kind(0.d0)
  !-input/output variables
  integer, intent(in)   :: n  
  real(DP), intent(in)  :: eps
  integer :: ind (n)  
  real(DP) :: ra (n)
  !-local variables
  integer :: i, ir, j, l, iind  
  real(DP) :: rra  
!
  ! initialize index array
  IF (ind (1) .eq.0) then  
     DO i = 1, n  
        ind (i) = i  
     ENDDO
  ENDIF
  ! nothing to order
  IF (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
    ! still in hiring phase
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    ELSE  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       IF ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  

  END DO sorting    
contains 
  !  internal function 
  !  compare two real number and return the result
  logical function hslt( a, b )
    REAL(DP) :: a, b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt
  !
end subroutine hpsort_eps_epw

!program main
!  INTEGER,parameter :: n=10
!  REAL*8  :: ra(n)
!  INTEGER :: ind(n)
!  INTEGER :: i
!  REAL*8  :: eps
!  
!  ra(1) = 3.0
!  ra(2) = 1.0
!  ra(3) = 0.1
!  ra(4) = 0.2
!  ra(5) = 1.0
!  ra(6) = 10.
!  ra(7) = 11.
!  ra(8) = 0.01
!  ra(9) = 0.02
!  ra(10)= 0.0
!  ind(:) = 0
!
!  eps = 1e-10
!  call hpsort_eps_epw (n, ra, ind, eps)
!
!  do i=1,n
!     print *, i, ind(i), ra(i)
!  enddo
!  
!end program main

!program main
!  integer :: i, nkp
!  integer :: div(3)
!  data div / 10, 10, 10/
!  integer :: kp(3)
!  integer, allocatable :: kqid(:,:)
!  
!  nkp  = div(1)*div(2)*div(3)
!
!  do ik=1,nkp
!     call indx2k(kp, ik-1, div)
!     print *, ik-1, kp
!  enddo
!  
!  allocate( kqid(nkp,nkp) )
!
!  call k_m_q_id(kqid, div, nkp)
!
!  do ik=1,nkp
!     do iq=1,nkp
!        print *, ik, iq, kqid(ik,iq)
!     end do
!  end do
!  deallocate( kqid )
!end program main
