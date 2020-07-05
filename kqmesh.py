# @Copyright 2020 Kristjan Haule

from scipy import *
import itertools
from timeit import default_timer as timer

import for_kpts as fkp
from inout import *

class KQmesh:
    """This class generates the k- and q-point grids and initilizes tetrahedral integration. 
    The result for the k-mesh has to be identical to that of the WIEN2k code. 
    The q-mesh is the same but without k0shift. 
    The tetrahedra are different since the q-vector breaks their symmetry all of them have to be taken into account separatedly.
    The geometrical weights of the set of q points depends on the corresponding k-point...
    """
    def k2indx(self, k):
        return (k[0]*self.ndiv[1] + k[1])*self.ndiv[2] + k[2]
    def indx2k(self, ind):
        k2 = ind %  self.ndiv[2]
        i1 = ind / self.ndiv[2]
        k1 = i1 % self.ndiv[1]
        k0 = i1 / self.ndiv[1]
        return array([k0,k1,k2])
    
    def __init__(self, nkdivs, k0shift, strc, latgen, fout):
        self.ndiv = nkdivs
        self.shift = k0shift
        self.aaa = array([strc.a, strc.b, strc.c])
        self.gbas = latgen.gbas
        
        if latgen.ortho or strc.lattice[1:3]=='CXZ':
            self.k2icartes = array(round_(latgen.br2*self.aaa/(2*pi)),dtype=int)
            self.k2cartes = 2*pi/self.aaa*identity(3) # == latgen.pia
        else:
            self.k2icartes = identity(3,dtype=int)
            self.k2cartes  = latgen.br2
            
        print >> fout, 'Transformation to semi-cartesian coordinates (k2icartes)='
        for i in range(3):
            print >> fout, ('%8.5f '*3) % tuple(self.k2icartes[i,:])
        print >> fout, 'Transformation from semi-cartesian to cartesian (k2cartes)='
        for i in range(3):
            print >> fout, ('%8.5f '*3) % tuple(self.k2cartes[i,:])
       
        
        nkp = nkdivs[0]*nkdivs[1]*nkdivs[2]
        self.sym = []
        for isym in range(strc.Nsym):
            t1   = dot(latgen.tizmat[isym,:,:].T, k0shift[:3])      # symmetry operation on k0shift
            toff = array([t1[i] % (2*nkdivs[i]) for i in range(3)]) # negative i numbers become n-i
            diff = abs(toff - k0shift[:3])                          # how much does symmetry operation change k0shift
            if sum(diff % 2) < 1e-6:
                self.sym.append( isym )
            
        if len(self.sym) < strc.Nsym:
            print >> fout, 'WARNING in KQmesh.init : The k-point offset selected reduces the symmetry group of the mesh, resulting in a larger number of irred. k-points'
        
        tm1 = timer()
        
        # This part calculates a set of reduced k-points for the integration
        # with the tetrahedron method. The k-points are reduced by the symmetry and 
        # a weight is put to each irreducible k-point to tell how many points it 
        # represents before the symmetry operation.
        divsh = 2
        self.kii_ind = -ones(nkp, dtype=int)
        self.weight=[]
        kk = zeros(3,dtype=int)
        self.iksym = -ones(nkp, dtype=int)
        for k_id,iii in enumerate(itertools.product(range(nkdivs[0]),range(nkdivs[1]),range(nkdivs[2]))):
            #k_id = self.k2indx(iii) # gives unique integer index, corresponding to three integeres for k-point (iii[0]/N,iii[1]/N,iii[2]/N)
            # can be calculated either with k2indx or just iterated
            if self.kii_ind[k_id] < 0: # we did not came across this k_id yet
                nirkp = len(self.weight)
                self.kii_ind[k_id] = nirkp
                kk = divsh*array(iii) + k0shift[:3]   # Notice that iii = (kk-k0shift)/divsh
                #print 'kk=', kk, 'iii=', iii, 'k_id=', k_id
                starr = set([])
                for i,isym in enumerate(self.sym):
                    # here kk is simple lattice 1BZ point, like kk=(i/N,j/N,k/N)
                    # t1 is now another k-point of the star, but is not necessary in the  1BZ
                    t1 = dot(latgen.tizmat[isym,:,:].T, kk)
                    # t1 is now brought to the 1BZ, and becomes t2. t2 is equivalent to t1.
                    t2 = array([t1[j] % (divsh*nkdivs[j]) for j in range(3)])
                    # this is now translated into lattice notation of kk, i.e., (i/N,j/N,k/N). As above iii=(kk-k0shift)/divsh
                    iii_s = (t2 - k0shift[:3])/divsh
                    ks_id = self.k2indx(iii_s)
                    #print 't1=', t1, 't2=', t2, 'ii_s=', ii_s
                    starr.add(ks_id)
                    if self.kii_ind[ks_id] >= 0:
                        #if kpav[ks_id] != kpav[k_id]:
                        if self.kii_ind[ks_id] != self.kii_ind[k_id]:
                            #print 'WARNING : previous kpav['+str(ks_id)+']='+str(kpav[ks_id])+' new kpav['+str(k_id)+']='+str(kpav[k_id])
                            print 'WARNING : previous kii_ind['+str(ks_id)+']='+str(self.kii_ind[ks_id])+' new kii_ind['+str(k_id)+']='+str(self.kii_ind[k_id])
                            #kpav[ks_id] = kpav[k_id]
                            self.kii_ind[ks_id] = self.kii_ind[k_id]
                    else:
                        #kpav[ks_id] = kpav[k_id]
                        self.kii_ind[ks_id] = self.kii_ind[k_id]
                    if self.iksym[ks_id]<0: # The group operation is not unique. We would rather save the first one and not the last one
                        self.iksym[ks_id] = isym
                    #print ' %2d %2d' % (i, ks_id), ' ks=(%3d%3d%3d)'%tuple(iii_s), 't1=(%3d%3d%3d)'%tuple(t1), 't2=(%3d%3d%3d)'%tuple(t2), 'kpav=', kpav, 'star=', starr
                # weight of each  k-point depends on how many k-points are in each star.
                # len(starr) gives all unique k-points we found that are in this particular star.
                self.weight.append( len(starr) )
                #print 'weight['+str(nirkp)+']=', weight[-1]
                #print '-'*64
            #else:
                #print 'iii=', iii, 'k_id=', k_id
        
        tm2 = timer()
        print >> fout, 'Generating k-point mesh'
        print >> fout, '## KQm:init: t(find irred k-points)=%14.9f' % (tm2-tm1,)
        
        # Now we construct the inverse index, which will give first long index for one member of the star
        self.k_ind = zeros(len(self.weight),dtype=int)
        for i in range(len(self.kii_ind)-1,-1,-1):
            self.k_ind[self.kii_ind[i]] = i

        # common multiple of all divisions, including in the presence of the shift
        if sum(abs(array(self.shift[:3]))) > 1e-3:
            # because shift is 1/2 of div, we need to increase common multiple
            self.LCM = reduce(lcm, [self.ndiv[i] if shift[i]==0 else self.ndiv[i]*2 for i in range(3)] )
        else:
            self.LCM = reduce(lcm, self.ndiv)
            
        tm3 = timer()
        print >> fout, '## KQm:init: t(inverse kii_ind)    =%14.9f' % (tm3-tm2,)

        print >> fout, 'Legend for k-mesh below:'
        print >> fout, '   kii_ind -- index to the corresponding irreducible k-point, i.e., the unique index for irreducible k-point'
        print >> fout, '   readkp  -- index to the same corresponding irreduible k-point, but now in terms of its position in this list'
        print >> fout, '   iksym   -- symmetry operation which converts irreducible k-point to the reducible k-point'
        print >> fout, '   following three columns -- (b1,b2,b2) fractional coordinates along the reciprocal lattice vectors. They only go from 0 to 1.'
        print >> fout, '   i   kii_ind[i] readkp[i]   iksym[i]'
        kshft = array(k0shift[:3])
        for i in range(len(self.kii_ind)):
            kk = (self.indx2k(i) + kshft/float(divsh)) * 1./nkdivs
            print >> fout, '%4d  %4d       %4d        %4d' % (i, self.kii_ind[i], self.k_ind[self.kii_ind[i]], self.iksym[i]), '  ', ('%7.4f '*3) % tuple(kk)

            
        #for i in range(len(self.k_ind)):
        #    print >> fout, 'k_ind['+str(i)+']='+str(self.k_ind[i])
        
        # ikpredin == kii_ind
        print >> fout, 'Total  weight of k-points=', sum(self.weight), ' expected to  be', nkdivs[0]*nkdivs[1]*nkdivs[2]
        print >> fout, 'Legend for the irreducible k-points list below:'
        print >> fout, '  kpt w/o shft -- k-point in lattice representation, which needs to be divided by '+str(self.LCM)+' to get k-point coordinate'
        print >> fout, '  kpt w shft   -- k-point in the same lattice representation, but with possible shift for half latice unit'
        print >> fout, '  wgh          -- weight of each k-point'
        print >> fout, '  kpt cartesian1-- semi cartesian k-point representation without common factor LCM='+str(self.LCM)
        print >> fout, '  kpt cartesian2-- semi cartesian representation including LCM'
        print >> fout, '  kpt cartesian3-- cartesian representation including LCM and proper length of the three reciprocal vectors'
        print >> fout, ' ik  kpt w/o shft     kpt w shft          wgh    kpt-cartesian1  kpt-cartesian2  kpts-cartesian3'
        for i in range(len(self.weight)):
            iii = self.indx2k(self.k_ind[i])
            kl = array(iii) + array(k0shift[:3])/float(divsh)
            _klist0_ = array(round_(kl*self.LCM/self.ndiv), dtype=int)
            _klist_ = dot(self.k2icartes, _klist0_ )
            #kc = dot(self.k2icartes, kl)/self.ndiv
            kc = _klist_/float(self.LCM)
            kcc = dot(self.k2cartes, kc)
            print >> fout, '%3d' %(i,), '[%3d,%3d,%3d]' % tuple(iii), '   ', ('[%4.1f,%4.1f,%4.1f]') % tuple(kl), '%5d' % (self.weight[i],), '   ', _klist_, kc, kcc

        Print_ALL = False
        if Print_ALL:
            print >> fout, 'All k-points'
            print >> fout,  ' ik irrk kpt_w/o_shft      kpt_w_shft            kpt_cartesian'
            for ik,iii in enumerate(itertools.product(range(self.ndiv[0]),range(self.ndiv[1]),range(self.ndiv[2]))):
                kl = array(iii) + array(k0shift[:3])/float(divsh)
                kc = dot(self.k2icartes, kl)/self.ndiv
                kcc = dot(self.k2cartes, kc)
                print >> fout, '%3d %3d' %(ik,self.kii_ind[ik]), '[%3d,%3d,%3d]' % tuple(iii), '   ', ('[%4.1f,%4.1f,%4.1f]') % tuple(kl), '   ', kc, kcc
            
        tm4 = timer()
        print >> fout, '## KQm:init: t(printing)           =%14.9f' % (tm4-tm3,)
            
        # ikpredin() == kii_ind[]
        # readkp()   == k_ind[kii_ind[i]]
        # ikpid()    == inverse(k_ind)
        # ikp()      == self.indx2k(k_ind[i])
        # klist()    == array(self.indx2k(k_ind[i]))*divsh + array(k0shift[:3])
        ### generic k-point can be obtained by
        #      isym = iksym[i]
        #      k_i = IZMAT[isym] * k_irr
        #
        # ikpid(redkp[cornid]) = kii_ind[i]

    def tetra(self, latgen, fout):
        """The tetrahedra are tried to be related by the symmetry. The tetrahedron is defined
           by the irreducible k-point on the vertices. If the vertices of two tetrahedra goes
           into the same irreducible point set, the weight of this tetrahedron will be added 
           once to reflect this relation. After this, we can just do the integration over the 
           irreducible tetrahedron and multiply it by a weight before summation. 
        """
        div = array(self.ndiv, dtype=int)
        divsh = 2
        
        tm1 = timer()
         
        # This part selects the tetrahedra so that they share the shortest 
        # diagonal of the containing cube. The tetrahedra will be divided according
        # to this diagonal to lower the possible error of linearization method.
        p0 = [[0,0,0,0,1,1,1,1],
              [0,0,1,1,0,0,1,1],
              [0,1,0,1,0,1,0,1]]
        p0 = array(p0)
        diag = zeros(4)
        # calculate main diagonals
        for i in range(4):
            sm = dot( [(p0[k,i]-p0[k,7-i])/float(div[k]) for k in range(3)], self.gbas )
            diag[i] = sum(sm**2)
        # find smallest diagonal
        self.mnd = argmin(diag)
        
        print >> fout, 'diag=', diag, 'armin(diag)=', self.mnd
        
        tet0=[[[0,0,0], [0,0,1], [0,1,1], [1,1,1]], 
              [[0,0,0], [0,1,1], [0,1,0], [1,1,1]],
              [[0,0,0], [0,1,0], [1,1,0], [1,1,1]],
              [[0,0,0], [1,1,0], [1,0,0], [1,1,1]],
              [[0,0,0], [1,0,0], [1,0,1], [1,1,1]],
              [[0,0,0], [1,0,1], [0,0,1], [1,1,1]]]
        tet0 = array(tet0)
        tet = zeros(shape(tet0),dtype=int)
        # rotate tetraedra
        if self.mnd==0:
            tet[:,:,:] =   tet0[:,:,:]
        elif self.mnd==1:
            tet[:,:,0] =   tet0[:,:,0]
            tet[:,:,1] =   tet0[:,:,2]
            tet[:,:,2] = 1-tet0[:,:,1]
        elif self.mnd==2:
            tet[:,:,0] =   tet0[:,:,0]
            tet[:,:,1] = 1-tet0[:,:,2]
            tet[:,:,2] =   tet0[:,:,1]
        elif self.mnd==3:
            tet[:,:,0] =   tet0[:,:,0]
            tet[:,:,1] = 1-tet0[:,:,1]
            tet[:,:,2] = 1-tet0[:,:,2]

        # atet and tetc are similar, except tetc contains all k-point and all possible tetrahedra
        self.tetc = zeros((4,div[0]*div[1]*div[2]*6),dtype=int, order='F')
        ctet = [0,0,0,0]
        wtet=[]
        all_tet={}
        for ik,iii in enumerate(itertools.product(range(div[0]),range(div[1]),range(div[2]))):
            orig = array(iii,dtype=int)
            #print 'iii=', iii
            for t in range(6):
                for i in range(4):
                    #print 'orig=', orig, 'tet=', tet[t,i,:], 'div=', div
                    corn = (orig+tet[t,i,:]) % div
                    #print 'corn=', corn
                    #cornid = self.k2indx(corn)
                    cornid = fkp.k2indx(corn,div)
                    ctet[i] = self.kii_ind[cornid] # ireducible k-points for each corner.
                    self.tetc[i,6*ik+t]= cornid    # reducible k-points in all four corners
                #print '     ', t, '%3d%3d%3d%3d' % tuple(array(ctet)+1)
                ktet = tuple(sorted(ctet))
                #print 'ktet=', ktet, 'ctet=', ctet
                if all_tet.has_key(ktet):
                    it = all_tet[ktet]
                    wtet[it] += 1
                else:
                    it = len(wtet)     # how many unique tetra we already have
                    wtet.append(1)     # now we have one more  with  weight 1
                    all_tet[ktet] = it # at which index is this tetrahedron saved
                
                # here we know that ik==self.k2indx(orig)
                
                #print 'it=', it, 'wtet=', wtet, 'all_tet=', all_tet
        #self.atet = [[] for i in range(len(all_tet))]
        #for ktet in all_tet.keys():
        #    it = all_tet[ktet]
        #    self.atet[it] = ktet
        self.atet = zeros((4,len(all_tet)),order='F',dtype=int)
        for ktet in all_tet.keys():
            it = all_tet[ktet]
            self.atet[:,it] = ktet
        
        tm2 = timer()
        print >> fout, '## KQm:tetra:t(find irr tetrahedra)=%14.9f' % (tm2-tm1,)
        
        #print 'tetrahedron   tindx  wgh'
        #for ktet in all_tet.keys():
        #    it = all_tet[ktet]
        #    print ktet, "%4d %4d" % (it, wtet[it])

        self.ntet = 6*div[0]*div[1]*div[2]
        vtet = 1/float(self.ntet)
        self.wtet = array(wtet,dtype=float)/float(self.ntet)
        # this is not done in gap2. just here.
        #print 'wtet=', self.wtet
        
        print >> fout, 'Total weight of tetrahedra =', sum(self.wtet), 'expected weight=', 1 # self.ntet
        print >> fout, 'indx  tetrahedron  wgh'
        #for it,ktet in enumerate(self.atet):
        #    print >> fout, '%4d' % (it,), '(%2d,%2d,%2d,%2d)' % ktet, self.wtet[it]
        for it in range(len(self.wtet)):
            ktet = tuple(self.atet[:,it])
            print >> fout, '%4d' % (it,), '(%2d,%2d,%2d,%2d)' % ktet, (self.wtet[it]*self.ntet)
        
        # (outet(j,it),j=1,4) => self.atet[4,ntet]
        # wtet()              => self.wtet[ntet]
        
        print >> fout, 'Number of symmetry operations: nsym='+str(len(self.sym))
        print >> fout, 'Total number of k-points: nkp='+str(div[0]*div[1]*div[2])
        print >> fout, 'Number of irred. k-points:  len(weight)='+str(len(self.weight))
        print >> fout, 'Total number of tetrahedra: ntet='+str(self.ntet)
        print >> fout, 'Number of inequivalent tetrahedra: len(atet)='+str(shape(self.atet)[1])

        tm3 = timer()
        print >> fout, '## KQm:tetra t(print terahedra)    =%14.9f' % (tm3-tm2,)
        
        #### Here working with difference of two momenta: k-q
        nkp = div[0]*div[1]*div[2]
        shift = array(self.shift[:3])


        self.kirlist0 = zeros((len(self.weight),3),dtype=int)
        for i in range(len(self.weight)):
            iii = self.indx2k(self.k_ind[i])
            self.kirlist0[i,:] = ( iii*self.LCM + shift*self.LCM/divsh )/self.ndiv   # lattice vectors
            
        self.klist = zeros((nkp,3),dtype=int)
        for ik in range(nkp):
            iii = self.indx2k(ik)
            self.klist[ik,:]  = ( iii*self.LCM + shift*self.LCM/divsh )/self.ndiv  # maybe we can eliminate now klist. WARNING: This should be different for othogonal and .not.orthogonal
            
        tm4 = timer()
        print >> fout, '## KQm:tetra t(generate k-list)    =%14.9f' % (tm4-tm3,)

        if sum(abs(shift)) <= 1e-3:
            self.qlist = copy(self.klist)
            self.LCMq = self.LCM
        else:
            self.LCMq = reduce(lcm, div)
            self.qlist = zeros((nkp,3),dtype=int)
            for iq in range(nkp):
                iii = self.indx2k(iq)
                self.qlist[iq,:] = iii*(self.LCMq/div[:])
        
        tm5 = timer()
        print >> fout, '## KQM:tetra t(generate q-list)    =%14.9f' % (tm5-tm4,)
        
        #vt = 1./self.ntet
        
        #print >> fout, 'k-points q-points k_cartesian'
        #for ik in range(nkp):
        #    print >> fout, ik, self.klist[ik,:], self.qlist[ik,:]#, kc/self.LCM

        tm6 = timer()
        print >> fout, '## KQm:tetra t(print k and q)      =%14.9f' % (tm6-tm5,)
            
        # Given a q vector in the submesh coordinates, this code will give you the
        # data about how one tetrahedron is related with another one by this q vector.
        # outet has nothing to do with the q vector, it just gives the number k-points on 
        # the vertices of each tetrahedron. While sib(:) tells us how one tetrahedron is 
        # related to another by the q-vector. If sib(1)=10, it means the first tetrahedron
        # is related to the tenth tetrahedron by this q-vector. The tetrahedron is not 
        # reduced by symmetry in this subroutine.
        #print 'kqtetra='
        if FORT:
            self.kqid = fkp.k_m_q_id(div,nkp)
        else:
            self.kqid = zeros((nkp,nkp),dtype=int)
            for iq,jjj in enumerate(itertools.product(range(div[0]),range(div[1]),range(div[2]))):
                qr = array(jjj,dtype=int)
                for ik,iii in enumerate(itertools.product(range(div[0]),range(div[1]),range(div[2]))):
                    kr = array(iii,dtype=int)
                    kmq = (kr-qr) % div
                    ikq = self.k2indx(kmq)
                    self.kqid[ik,iq] = ikq

        # remember:
        #            ik = 0...nkp and t=0...6
        #            linkt[iq, 6*ik+t] = 6*self.kqid[ik,iq]+t
        #
        #print >> fout, 'tgenq:'
        #for iq in range(nkp):
        #    print >> fout, '%4d ' % (iq+1,), '%3d%3d%3d' % tuple(self.qlist[iq])
        #    for ik in range(nkp):
        #        for t in range(6):
        #            print >> fout, '          %4d  %4d' % (6*ik+t+1, 6*self.kqid[ik,iq]+t+1)
        
                    
            
        tm7 = timer()
        print >> fout, '## KQm:tetra t(get kqid)           =%14.9f' % (tm7-tm6,)

        print >> fout, ' idvk=', self.LCM, ' idvq=', self.LCMq
        if latgen.ortho or strc.lattice[1:3]=='CXZ':
            # just for debugging
            self.kirlist = copy(self.kirlist0)
            print >> fout, 'Changing klist and kirlist to cartesian coordinates'
            for ik in range(len(self.kirlist)):
                self.kirlist[ik,:] = dot(self.k2icartes, self.kirlist[ik,:] )
            for ik in range(len(self.klist)):
                self.klist[ik,:] = dot(self.k2icartes, self.klist[ik,:] )
            self.qlistc = zeros(shape(self.qlist),dtype=int)
            for iq in range(nkp):
                self.qlistc[iq,:] = dot(self.k2icartes, self.qlist[iq,:])
        else:
            self.kirlist = self.kirlist0
            
        if True:
            print >> fout, '#k-mesh was transformed to cartesian=', (latgen.ortho or strc.lattice[1:3]=='CXZ')
            print >> fout, '#Irreducible k-points with LCM=', self.LCM
            print >> fout, ' ik kpt  cartesian    wgh   kpt cartesian/LCM'
            for ik in range(len(self.kirlist)):
                kc = self.kirlist[ik,:]
                kcr = kc/float(self.LCM)
                kcc = dot(self.k2cartes, kcr)
                print >> fout, '%3d' %(ik,), '[%3d,%3d,%3d]' % tuple(kc), ' %5d' % (self.weight[ik],), '  ', kcr, kcc
            
            print >> fout, '# All k-points with LCM=', self.LCM
            print >> fout, '  ik kpt cartesian     kpt cartesian'
            for ik in range(nkp):
                kc = self.klist[ik,:]
                kcr = kc/float(self.LCM)
                kcc = dot(self.k2cartes, kcr)
                print >> fout, '%3d' %(ik,), '[%3d,%3d,%3d]' % tuple(kc), '   ', kcr, kcc
                
            print >> fout, '# All q-points LCM=', self.LCMq
            print >> fout, ' iq     q_cartesian                                   iq      idvq'
            for iq in range(nkp):
                qc = self.qlistc[iq,:]
                qcr = qc/float(self.LCMq)
                qcc = dot(self.k2cartes, qcr)
                print >> fout, '%3d  %15.10f%15.10f%15.10f' % (iq,qcc[0],qcc[1],qcc[2]), ' %3d%3d%3d' % tuple(qc), '  ', self.LCMq
            
