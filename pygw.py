#!/usr/bin/env python
# @Copyright 2020 Kristjan Haule

Parallel = True
if Parallel :
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    msize = comm.Get_size()
    mrank = comm.Get_rank()
    master=0
else:
    msize = 1
    mrank = 0
    master = 0

ForceReal = False       # This forces Coulomb interaction to be truncated to real part only, and is not exact
Real_sqrt_olap = True   # This is a good way of ensuring that basis in the interstitials is as close to plane-waves as possible, but orthogonal
Real_sqrt_Vcoul = False # This similarly ensures that sqrt(V_c) is as close as possible to V_c in terms of phase. But it can make the matrix slightly bigger
DMFT1_like = False      # This forces the reducible k-points to be handled like in dmft1 code. Somehow, this has difficulty for more than one atom per unit cell.
Save_eps = False
NEW_KCW = True
Debug_Print = False

from scipy import *
import re
import os
import os,psutil,sys
import itertools
import glob
from scipy import linalg
from scipy import optimize
from pylab import linalg as la
from timeit import default_timer as timer
from fractions import gcd
from scipy import special
import numpy.polynomial.legendre as legendre
import numpy.polynomial.laguerre as laguerre

import gwienfile as w2k
import gaunt as gn          # gaunt coefficients
import for_kpts as fkp      # 
import for_tetrahedra as ft #
import radials as rd        # radial wave functions
import radd                 # radial derivative of a function
import fnc                  # fermi and gauss functions
import lapwcoeff as lapwc   # computes lapw coefficients alm,blm,clm
import for_vxcnn as fvxcn   # matrix elements for Vxcn
import for_Coulomb as fCl
import for_q_0 as f_q_0

from inout import *
from kqmesh import *
import mcommon as mcmn

    
class ProductBasis:
    nspin_mb=1 # option for constructing mixed basis in the spin-polarized case (nspin=2) 
               # 2 or 1 -- 2==use both spin up and down radial wave functions, 1==use only spin up radial wave functions
    mb_ludot = False  # control whether considering \dot{u}_l when setting up the mixed basis 
    mb_lcore = True   # control whether considering core states when building the mixed basis 
    #lmbmax   = 2     # maximum l of radial functions considered when building the mixed basis set
    #lblmax   = 0     # the maximum L for the mixed basis (big l)
    #mb_emin = -1e10   # the energy cut-off for mixed basis funcitons, any radial functions  
    #mb_emax = 20.0    # with energy lower than mb_emin are dropped in constructing mixbasis 
    #
    def __init__(self, io_data, strc, in1, radf, core, nspin, fout):
        (self.kmr, self.pwm, self.lmbmax, self.wftol, self.lblmax, self.mb_emin, self.mb_emax) = io_data # decoding input files from input

        # maximum l at each atom is encoded in lmbmax, and each digit correspond to an atom, i.e., 34 means lmax=3 and lmax=4 at first
        # and second atom, respectively. If the number is negative, we use some default value, depending on Z.
        # If the number is positive, we use the same lmax=lmbmax at each atom.
        self.lmbmax_at = zeros(strc.nat,dtype=intc)
        if self.lmbmax > 10:
            ii, p = self.lmbmax, []
            while (ii>10):
                p.append( ii % 10 )
                ii = ii/10
            p.append(ii)
            p = p[::-1]
            for iat in range(min(strc.nat,len(p))):
                self.lmbmax_at[iat] = p[iat]
        else:  
            if self.lmbmax < 0:
                lvmax_at = zeros(strc.nat, dtype=int)
                for iat in range(strc.nat):
                    znuc= strc.Znuc[iat]
                    if znuc <= 2.0 :
                        lvmax_at[iat] = 0
                    elif znuc > 2.0 and znuc <= 18.0:
                        lvmax_at[iat] = 1
                    elif znuc > 18.0 and znuc <= 54.0:
                        lvmax_at[iat] = 2
                    else:
                        lvmax_at[iat] = 3
                self.lmbmax_at[:] = lvmax_at[:] + abs(self.lmbmax) 
            elif self.lmbmax > 0:
                self.lmbmax_at[:] = self.lmbmax
        # similarly decoding lblmax_at to self.lblmax_at
        self.lblmax_at = zeros(strc.nat,dtype=intc)
        if self.lblmax > 10:
            ii, p = self.lblmax, []
            while (ii>10):
                p.append( ii % 10 )
                ii = ii/10
            p.append(ii)
            p = p[::-1]
            for iat in range(min(strc.nat,len(p))):
                self.lblmax_at[iat] = p[iat]
        else:
            if self.lblmax <= 10:
                self.lblmax_at = self.lmbmax_at*2
            else:
                self.lblmax_at[:] = self.lblmax
      
        self.lmbmax = max(self.lmbmax_at)
        self.lblmax = max(self.lblmax_at)
        print >> fout, '  Atom-specific lmbmax==(products made from u_l with l<lmbmax) and lblmax==(u_l1*u_l2->w_L with  L<lblmax):'
        print >> fout, '     iat      lmbmax   lblmax'
        for iat in range(strc.nat):
            print >> fout, '    %4d %8d %8d' % (iat,self.lmbmax_at[iat],self.lblmax_at[iat])
        
        #print >> fout, 'lmbmax= ', self.lmbmax
        #print >> fout, 'lblmax= ', self.lblmax 
        #
        #     Initial estimation of the maximum number of mixed basis functions       
        lomaxp1 = shape(in1.nlo)[0]
        nmixmax = (in1.nt + lomaxp1 * in1.nlomax) * (self.lmbmax+1)*(self.lmbmax+1)*self.nspin_mb
        
        self.ncore = zeros(strc.nat, dtype=int)
        self.big_l=[]
        self.ul_product=[]
        self.us_product=[]
        # calculate the radial part of the mixed basis functions  (mb_setumix)
        for iat in range(strc.nat):
            rp, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
            
            _big_l_=[]
            _ul_product_=[]
            _us_product_=[]
            
            print >> fout, ' '*10+'Product basis functions for atom '+strc.aname[iat]+'\n'
            nrf = self.count_nrf(iat,in1,core)
            #print >> fout, "# rad. func. to be considered on the atom %2d is %3d " % (iat+1, nrf)
            nup = nrf*(nrf+1)/2*self.nspin_mb
            #print >> fout, "mb_setuprod: # of uprod for atom %2d is %3d " % (iat,nup)
            ulprod = []
            usprod = []
            eles = []
            for isp in range(min(self.nspin_mb,nspin)):
                ul_all = zeros((nrf,npt))
                us_all = zeros((nrf,npt))
                l_all = zeros(nrf,dtype=int)
                #### collect_u
                irf = 0 
                # core states
                for ic,ecore in enumerate(core.eig_core[isp][iat]):
                    if ecore > self.mb_emin:
                        l_all[irf] = core.l_core[iat][ic]
                        ul_all[irf,:] = core.ul_core[isp][iat][ic][:]
                        us_all[irf,:] = core.us_core[isp][iat][ic][:]
                        irf = irf + 1
                self.ncore[iat] = irf
                for l in range(self.lmbmax_at[iat]+1):
                    l_all[irf] = l 
                    ul_all[irf,:] = radf.ul[isp,iat,l,:npt]
                    us_all[irf,:] = radf.us[isp,iat,l,:npt]
                    irf = irf + 1 
                    if self.mb_ludot:
                        l_all[irf] = l
                        ul_all[irf,:] = radf.udot [isp,iat,l,:npt]
                        us_all[irf,:] = radf.usdot[isp,iat,l,:npt]
                        irf = irf + 1
                    
                    if (l < lomaxp1):
                        for ilo in in1.nLO_at_ind[iat][l]:
                            #print 'in local orbital ilo=', ilo, 'l=', l, 'E=', in1.Elo[0,iat,l,ilo], 'cutoff=', self.mb_emax
                            if in1.Elo[0,iat,l,ilo] >= self.mb_emax : continue
                            l_all[irf] = l
                            ul_all[irf,:] = radf.ulo [isp,iat,l,ilo,:npt]
                            us_all[irf,:] = radf.uslo[isp,iat,l,ilo,:npt]
                            irf = irf + 1
                # construct all possible products of radial functions
                for jrf in range(nrf):
                    l2, u2, us2 = l_all[jrf], ul_all[jrf,:], us_all[jrf,:]
                    for irf in range(jrf+1):
                        l1, u1, us1 = l_all[irf], ul_all[irf,:],  us_all[irf,:]
                        eles.append( (l1,l2) )
                        ulprod.append( u1[:]*u2[:]/rp[:] )
                        usprod.append( us1[:]*us2[:]/rp[:] )
            # calculate the overlap matrix of product functions
            olap = zeros((len(eles),len(eles)))
            for i in range(len(eles)):
                for j in range(i+1):
                    olap[i,j] = rd.rint13g(strc.rel, ulprod[i], usprod[i], ulprod[j], usprod[j], dh, npt, strc.r0[iat])
                    olap[j,i] = olap[i,j]

            #for i in range(len(olap)):
            #    for j in range(i,len(olap)):
            #        print "%2d %2d  %12.7f " % (i,j,olap[i,j])


            nl = zeros(self.lblmax_at[iat]+1, dtype=int)
            for i in range(len(eles)):
                l_max = eles[i][0]+eles[i][1]
                l_min = abs(eles[i][0]-eles[i][1])
                for l in range(self.lblmax_at[iat]+1):
                    if l_min <= l and l <= l_max:
                        nl[l] += 1

            #print 'nl=', nl
            # diagonalize the overlap matrix of the product functions for each  L-block
            for L in range(self.lblmax_at[iat]+1):
                if nl[L] == 0: continue
                acceptable_func = filter(lambda i: ( L >= abs(eles[i][0]-eles[i][1]) and L <= eles[i][0]+eles[i][1] ), range(len(eles)))
                #print >> fout, ' '*14+'for L=',L,'acceptable_func=', acceptable_func
                # prepare those functions that enter this L
                ul_prod = zeros( (npt,len(acceptable_func)) )
                us_prod = zeros( (npt,len(acceptable_func)) )
                for i,i1 in enumerate(acceptable_func):
                    ul_prod[:,i] = ulprod[i1][:]
                    us_prod[:,i] = usprod[i1][:]
                
                # the L-block overlap matrix (uml)
                uml = zeros( (nl[L],nl[L]) )
                # generate the L-block overlap matrix
                #print >> fout, (' '*14)+' - generate the L-block overlap matrix'
                for i,i1 in enumerate(acceptable_func):
                    uml[i,i] = olap[i1,i1]
                    #print 'uml['+str(i)+','+str(i)+']='+str(uml[i,i])
                    for j,i2 in enumerate(acceptable_func[i+1:]):
                        #print 'i,k=', i, i+j+1
                        uml[i ,i+j+1] = olap[i1,i2]
                        uml[i+j+1,i ] = olap[i1,i2] 
                #print 'uml=', uml
                #print 'ind=', ind
                #print 'acceptable_func=', acceptable_func
                w, Ud = linalg.eigh(uml)

                finite_contribution = filter(lambda i: w[i]>self.wftol, range(len(acceptable_func)) )
                #print >> fout, ' '*14+'finite_contribution=', finite_contribution
                Udn = zeros( (len(acceptable_func),len(finite_contribution)) )
                for i,j in enumerate(finite_contribution):
                    #print 'Udn['+str(i)+']= Ud['+str(j)+']'
                    Udn[:,i] = Ud[:,j]
                
                #print 'w=', w, 'v=', Ud
                ul_prodn = dot(ul_prod, Udn)
                us_prodn = dot(us_prod, Udn)

                #print >> fout, ' '*14+'shape(ul_prodn)=', shape(ul_prodn)
                #  proper normalization of the resulting radial product wave functions
                #print >> fout, " - normalize the rad mixed func"
                norms=[]
                for i,i1 in enumerate(finite_contribution):
                    norm = rd.rint13g(strc.rel, ul_prodn[:,i], us_prodn[:,i], ul_prodn[:,i], us_prodn[:,i], dh, npt, strc.r0[iat])
                    norms.append(norm)
                    ul_prodn[:,i] *= 1./sqrt(norm)
                    us_prodn[:,i] *= 1./sqrt(norm)
                    _big_l_.append(L)
                    _ul_product_.append(ul_prodn[:,i])
                    _us_product_.append(us_prodn[:,i])

                print >> fout, (' '*5)+'L=%2d Nr. of products:%4d Nr. of basis functions%4d' % (L,  len(acceptable_func), len(finite_contribution) )  
                print >> fout, ' '*14+'functions have norm', norms

            self.big_l.append(_big_l_)
            self.ul_product.append(_ul_product_)
            self.us_product.append(_us_product_)
            
            nwf = sum([(2*L+1) for L in self.big_l[iat]])
            print >> fout, (' '*54)+'----'
            print >> fout, (' '*5)+'Total number of radial functions'+' '*17+'%4d    Maximum L %4d'  % ( len(self.big_l[iat]), max(self.big_l[iat]) )
            print >> fout, (' '*5)+'Total number of basis functions '+' '*17+'%4d' % (nwf,)

        print >> fout, '  maxbigl=', max([max(self.big_l[iat]) for iat in range(len(self.big_l))])
        #     Calculate the total number of mixed wave functions (including M)    
        #     = size of the local part of the matrices
        ww = [sum([(2*L+1)*strc.mult[iat] for L in self.big_l[iat]]) for iat in range(strc.nat)]
        loctmatsize = sum(ww)

        lmixmax    = max(ww)
        print >> fout, ' Max. nr. of MT-sphere wavefunctions per atom %6d' % (lmixmax,)
        print >> fout, ' Total  nr. of MT-sphere wavefunctions        %6d' % (loctmatsize,)
        #
        # set an array that stores the general index of the mixed
        # function for a given mixed function of a given atom
        #
        ndf = sum(strc.mult)
        self.locmixind = zeros((ndf,lmixmax),dtype=int)
        nmixlm = zeros(ndf, dtype=intc)
        self.atm = []
        self.Prod_basis=[]
        self.iProd_basis={}
        idf=0
        imix=0
        for iat in range(strc.nat):
            for ieq in range(strc.mult[iat]):
                self.atm.append( iat )
                im = 0
                for irm in range(len(self.big_l[iat])):
                    L = self.big_l[iat][irm]
                    for M in range(-L,L+1):
                        self.locmixind[idf,im] = imix
                        self.iProd_basis[(idf,irm,L,M)] = imix
                        self.Prod_basis.append( (idf,irm,L,M) )
                        im += 1
                        imix += 1
                nmixlm[idf]=im
                idf += 1

        print >> fout, ' List of all product basis functions:  total#:  ', imix
        for i in range(len(self.Prod_basis)):
            (idf,irm,L,M) = self.Prod_basis[i]
            print >> fout, '  %3d %s%2d %s%3d %s%2d %s%3d' % (i, 'idf=', idf, 'irm=', irm, 'L=', L, 'M=', M)


        #print 'locmixind=', self.locmixind
        #print 'nmixlm=', nmixlm

    def count_nrf(self,iat,in1,core):
      nrf = 0
      # core states
      for ic,ecore in enumerate(core.eig_core[0][iat]):
          if ecore > self.mb_emin:
              nrf = nrf + 1 
      # normal LAPW radial functions 
      if self.mb_ludot: # whether to consider u_dot
        nrf += (self.lmbmax_at[iat]+1)*2
      else:
        nrf += self.lmbmax_at[iat]+1

      # LO basis
      lomaxp1 = shape(in1.nlo)[0]
      for l in range(min(lomaxp1,self.lmbmax_at[iat]+1)): 
          for ilo in in1.nLO_at_ind[iat][l]:
              if in1.Elo[0,iat,l,ilo] < self.mb_emax :
                  nrf += 1
      return nrf
    
    def cmp_matrices(self, strc, in1, radf, core, nspin, fout):
        """Calculated all matrices related to radial mixbasis functions
         
           Some general matrices (independent of eigen-vectors ) related to
           bare Coulomb interaction:
             - rtl
             - rrint
             - tilg
        """
        print >> fout, 'set_mixbasis: calc mixmat'
        import cum_simps as cs

        cfein = 1/137.0359895**2 if strc.rel else 1e-22

        N = max([len(self.big_l[iat]) for iat in range(strc.nat)])
        self.rtl   = zeros((strc.nat,N))          # < r^L | u_{im,at} >
        self.rrint = zeros((strc.nat,N*(N+1)/2))  # <u_{jm,L,at}| (r_<)^L / (r_>)^{L+1} | u_{im,L,at} >/
        for iat in range(strc.nat):
            #npt = strc.nrpt[iat]
            #dh  = log(strc.rmt[iat]/strc.r0[iat])/(npt - 1)      # logarithmic step for the radial mesh
            #dd = exp(dh)
            #rp = strc.r0[iat]*dd**range(npt)            
            rp, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
            #rm   = rp[::-1]
            for irm in range(len(self.big_l[iat])):
                L = self.big_l[iat][irm]
                r_to_L = rp**(L+1) 
                # Int[ u_i(r) r^{L+1}, {r,0,R}]
                rxov = rd.rint13g(strc.rel, self.ul_product[iat][irm], self.us_product[iat][irm], r_to_L, r_to_L, dh, npt, strc.r0[iat])
                if (rxov < 0):
                  self.rtl[iat,irm] = -rxov
                  self.ul_product[iat][irm] *= -1
                  self.us_product[iat][irm] *= -1
                else:  
                  self.rtl[iat,irm] = rxov
                
        for iat in range(strc.nat):
            #npt = strc.nrpt[iat]
            #dh  = log(strc.rmt[iat]/strc.r0[iat])/(npt - 1)      # logarithmic step for the radial mesh
            #dd = exp(dh)
            #rp = strc.r0[iat]*dd**range(npt)
            rp, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
            rm   = rp[::-1]
            for irm in range(len(self.big_l[iat])):
                L = self.big_l[iat][irm]
                r_to_L = rp**(L+1) 
                # double radial integrals <u_j| (r_<)^L/(r_>)^{L+1} | u_i>
                # the oustide loop for double  radial integrals
                u_i = self.ul_product[iat][irm] + cfein*self.us_product[iat][irm]
                u_i_r_to_L = cs.cum_simps(6, rp, u_i[:]*r_to_L[:])  # Int[ u_i(r)*r * r^L, {r,0,rx}]
                w = u_i[:]*rp[:]/r_to_L[:]
                u_i_m_to_L = cs.cum_simps(6, rm, w[::-1])[::-1]     # -Int[ u_i(r)*r / r^{L+1}, {r,R,rx}]
                
                for jrm in range(irm,len(self.big_l[iat])):
                    # the inside loop for double radial integrals
                    if self.big_l[iat][jrm] != L: continue
                    u_j = self.ul_product[iat][jrm] + cfein*self.us_product[iat][jrm]
                    #   Int[ u_j(rx)*( rx /rx^{L+1} * Int[ u_i(r)*r * r^L, {r,0,rx}] - rx^{L+1} * Int[ u_i(r)*r / r^{L+1}, {r,R,rx}] ), {rx,0,R}]
                    # = Int[ u_j(rx)*rx * Int[ u_i(r)*r * (r_<)^L/(r_>)^{L+1}, {r,0,R}], {rx,0,R}]
                    uij = u_j[:]*(rp[:]/r_to_L[:]*u_i_r_to_L[:] - r_to_L[:]*u_i_m_to_L[:])
                    ijrm = irm + jrm*(jrm+1)/2 # store into this index
                    self.rrint[iat,ijrm] = cs.cum_simps(6,rp,uij)[-1]
                
        if Debug_Print:
            for iat in range(strc.nat):
                print >> fout, (' '*5)+'Radial integrals'
                print >> fout, (' '*13)+'N  L     <r^(L+1)|v_L>'+' '*9+'<r^(L+2)|v_L> '
                for irm in range(len(self.big_l[iat])):
                    L = self.big_l[iat][irm]
                    print >> fout, ((' '*10)+'%4d%3d%18.10f ') % (irm,L,self.rtl[iat,irm])
                    
            for iat in range(strc.nat):
                print >> fout, (' '*5)+'Double radial integrals'
                print >> fout, (' '*13)+'N1  N2  L     <v_i| (r_<)^L/(r_<)^{L+1} |v_j>'
                for irm in range(len(self.big_l[iat])):
                    L = self.big_l[iat][irm]
                    for jrm in range(irm,len(self.big_l[iat])):
                        if self.big_l[iat][jrm] != L: continue
                        ijrm = irm + jrm*(jrm+1)/2
                        print >> fout, ((' '*10)+'%4d%4d %3d%18.10f%5d') % (irm,jrm,L,self.rrint[iat,ijrm],ijrm)

        lomaxp1 = shape(in1.nlo)[0]
        how_many_fnc = zeros((nspin,strc.nat),dtype=int)
        for isp in range(nspin):
            for iat in range(strc.nat):
                # first collect all core radial wavefunctions
                how_many_functions  = len(core.l_core[iat])
                # ul functions
                for l in range(in1.nt):
                    how_many_functions += 1
                # energy derivative ul_dot functions
                for l in range(in1.nt):
                    how_many_functions += 1
                # LO local orbitals
                for l in range(lomaxp1):
                    for ilo in in1.nLO_at_ind[iat][l]:
                        how_many_functions += 1
                how_many_fnc[isp,iat] = how_many_functions
        
        n_mix = max([len(self.big_l[iat]) for iat in range(strc.nat)])
        #self.s3r = zeros( (nspin,strc.nat,amax(how_many_fnc),amax(how_many_fnc),n_mix) )
        self.s3r = zeros( (n_mix,amax(how_many_fnc),amax(how_many_fnc),strc.nat,nspin), order='F' )
        orb_info = [' core','     ',' dot ',' lo  ']
        for isp in range(nspin):
            for iat in range(strc.nat):
                #npt = strc.nrpt[iat]
                #dh  = log(strc.rmt[iat]/strc.r0[iat])/(npt - 1)      # logarithmic step for the radial mesh
                #dd = exp(dh)
                #rp = strc.r0[iat]*dd**range(npt)
                rp, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
                
                rwf_all = []
                lrw_all = []
                # first collect all core radial wavefunctions
                for ic,lc in enumerate(core.l_core[iat]):
                    lrw_all.append((lc,0))
                    rwf_all.append( (core.ul_core[isp][iat][ic][:], core.us_core[isp][iat][ic][:]) )   # note that these wave functions are normalized to sqrt(occupancy/degeneracy)
                # ul functions
                for l in range(in1.nt):
                    rwf_all.append( (radf.ul[isp,iat,l,:npt],radf.us[isp,iat,l,:npt]) )
                    lrw_all.append((l,1))
                # energy derivative ul_dot functions
                for l in range(in1.nt):
                    rwf_all.append( (radf.udot[isp,iat,l,:npt],radf.usdot[isp,iat,l,:npt]) )
                    lrw_all.append((l,2))
                # LO local orbitals
                for l in range(lomaxp1):
                    for ilo in in1.nLO_at_ind[iat][l]:
                        rwf_all.append( (radf.ulo[isp,iat,l,ilo,:npt],radf.uslo[isp,iat,l,ilo,:npt]) )
                        lrw_all.append((l,3))
                # now performing the integration
                for ir2 in range(len(lrw_all)):
                    l2, it2 = lrw_all[ir2]
                    a2, b2 = rwf_all[ir2]
                    for ir1 in range(ir2+1):
                        l1, it1 = lrw_all[ir1]
                        a1, b1 = rwf_all[ir1]
                        a3 = a1[:] * a2[:] / rp[:]
                        b3 = b1[:] * b2[:] / rp[:]
                        for irm in range(len(self.big_l[iat])):  # over all product functions
                            L = self.big_l[iat][irm]             # L of the product function
                            if L < abs(l1-l2) or L > (l1+l2): continue # triangular inequality violated
                            ## s3r == < u-product-basis_{irm}| u_{a1,l1} u_{a2,l2} > 
                            rint = rd.rint13g(strc.rel, self.ul_product[iat][irm], self.us_product[iat][irm], a3, b3, dh, npt, strc.r0[iat])
                            #self.s3r[isp,iat,ir2,ir1,irm] = rint
                            #self.s3r[isp,iat,ir1,ir2,irm] = rint
                            self.s3r[irm,ir1,ir2,iat,isp] = rint
                            self.s3r[irm,ir2,ir1,iat,isp] = rint
                if Debug_Print:
                    print >> fout, (' '*5)+'Integrals <v_(NL) | u_(l1)*u_(l2)> for atom  %10s' % (strc.aname[iat],)
                    print >> fout, (' '*13)+'N   L   l1  u_   l2 u_        <v | u*u>'
                    for irm in range(len(self.big_l[iat])):  # over all product functions
                        L = self.big_l[iat][irm]             # L of the product function
                        for ir2 in range(len(lrw_all)):
                            l2, it2 = lrw_all[ir2]
                            for ir1 in range(ir2+1):
                                l1, it1 = lrw_all[ir1]
                                #if self.s3r[isp,iat,ir2,ir1,irm]!=0:
                                if self.s3r[irm,ir1,ir2,iat,isp]!=0:
                                    #print >> fout, ((' '*10)+('%4d'*3)+'%s%4d%s%19.11e') % (irm,L,l1,orb_info[it1],l2,orb_info[it2],self.s3r[isp,iat,ir2,ir1,irm])
                                    print >> fout, ('%4d%4d'+(' '*2)+('%4d'*3)+'%s%4d%s%19.11e') % (ir1,ir2,irm,L,l1,orb_info[it1],l2,orb_info[it2],self.s3r[irm,ir1,ir2,iat,isp])
                            
    def get_djmm(self,idf,l,m1,m2):
      if (abs(m2) > l) or (abs(m1) > l):
        return 0
      #  calculate the index of the djmm vector
      sgn_m2 = 1 if m2>=0 else -1  # note that this is like sign(m2), except at m2=0
      idx = l*(l+1)*(4*l-1)/6 + (l+1)*(sgn_m2*m1+l)+abs(m2)
      if m2 >= 0:
          # For m2>=0 D^j_m1m2=djmm(i)
          #print 'extracting l=', l, 'm1=', m1, 'm2=', m2, 'idx=', idx, 'val=', self.djmm[idf,idx]
          return self.djmm[idx,idf]
      else:
          # For m2<0 D^j_m1m2=(-1)^(m1-m2)djmm^*(i)
          m_1_m1m2 = 1-2*(abs(m1-m2) % 2)  # (-1)^(m1-m2)
          val = m_1_m1m2 * conj(self.djmm[idx,idf])
          #print 'extracting l=', l, 'm1=', m1, 'm2=', m2, 'idx=', idx, 'val=', val
          return val
    
    def generate_djmm(self,strc,latgen_trotij,in1_nt,fout):
        """
         Generates the rotation matrices for spherical harmonics $D^j_{mm'}$ up
        to \texttt{maxbigl} once
         $D^1_{mm'}$ is known using the inverse Clebsch-Gordan series:
        
        \begin{equation}
        D^j_{\mu m}=\sum\limits_{\mu_1=-1}^1{\sum\limits_{m_1=1}^1{C(1,j-1,j;%
        \mu_1,\mu-\mu_1)C(1,j-1,j;m_1,m-m_1)D^1_{\mu_1m_1}D^{j-1}_{\mu-\mu_1,m-m_1}}}
        \end{equation}.
        
         Since the program already calculates the Gaunt coefficients
        $G^{LM}_{l_1m_1,l_2m_2}$ we replace the product of Clebsch-Gordan
        coefficients using:
        
         \begin{equation}
        C(1,j-1,j;\mu_1,\mu-\mu_1)C(1,j-1,j;m_1,m-m_1)=\sqrt{\frac{4\pi(2j+1)}%
        {3(2j-1)}}\frac{G^{j\mu}_{1\mu_1,j-1\mu-\mu_1}G^{jm}_{1m_1,j-1m-m_1}}%
        {G^{j0}_{10,j-10}}
        \end{equation}
        """
        #
        #     calculate Gaunt coefficients
        #
        maxnt = max(in1_nt, 2*(self.lmbmax+1))
        import gaunt as gn
        self.cgcoef = gn.cmp_all_gaunt(maxnt)
        # print l1,l2,l3,m1,m2,m3, gn.getcgcoef(l1,l2,l3,m1,m2,cgcoef)
        
        iateq_ind=[]
        for iat in range(len(strc.mult)):
            iateq_ind += [iat]*strc.mult[iat]
        
        # Calculate the dimension of djmm
        #maxbigl = amax(self.big_l)
        maxbigl = max([max(self.big_l[iat]) for iat in range(len(self.big_l))])   # bug jul.7 2020
        dimdj=(maxbigl+1)*(maxbigl+2)*(4*maxbigl+3)/6
        ndf = len(iateq_ind)
        #print 'ndf=', ndf, 'dimdj=', dimdj, 'iateq_ind=', iateq_ind
        # Allocate djmm and initialize it to zero
        self.djmm = zeros((dimdj,ndf), dtype=complex, order='F')

        # Transform a rotation matrix from cartesian $(x,y,z)$ to spherical basis $(Y_{1,-11},Y_{1,0},Y_{1,1})$.
        s2 = 1/sqrt(2.)
        C2Sph = array([[ s2, -s2*1j,   0],
                       [ 0,    0,      1],
                       [-s2, -s2*1j,   0]], dtype=complex)
        #  Loop over all atoms
        for idf in range(len(iateq_ind)): # atom index counting all atoms
            iat = iateq_ind[idf]          # atom index counting only equivalent
            self.djmm[0,idf] = 1.0
            # Calculate the rotation matrix in the cartesian basis
            # rotcart=transpose(rotloc x rotij)
            rotcart = dot(strc.rotloc[iat], latgen_trotij[idf,:,:].T)
            #print 'rotcart=', rotcart, 'rotloc=', strc.rotloc[iat], 'rotij=', latgen_trotij[idf,:,:].T
            # Transform the rotation matrix to the spherical basis
            rotsph = dot( dot( C2Sph, rotcart), C2Sph.conj().T )
            # Obtain the rotation matrix for j=1 D^1_{m,m')=transpose(rotsph)
            #print 'rotsph=', rotsph
            
            #print 'rotsph=', rotsph
            for mu in [-1,0,1]:
                for m in [0,1]:
                    idx = 2*(mu+1)+m+1  # index for l=1
                    self.djmm[idx,idf] = rotsph[m+1,mu+1]
                    if Debug_Print:
                        print >> fout, ' djmm[idf=%2d,l=%2d,mu=%3d,m=%2d,idx=%4d]=%16.10f %16.10f' % (idf+1,1,mu,m,idx+1,self.djmm[idx,idf].real,self.djmm[idx,idf].imag)
            #  Obtain the rotation matrix for j> 1 by recursive relation
            for l in range(2,maxbigl+1):
                sql = sqrt(4*pi*(2*l+1)/(2*l-1.0)/3.)
                prefac = sql/gn.getcgcoef(1,l-1,l,0,0,self.cgcoef)
                #print '%2d %12.7f' % (l, prefac)
                _djmm_ = zeros( (2*l+1,l+1), dtype=complex )
                for mu in range(-l,l+1):
                    for mu1 in [-1,0,1]:
                        mu2 = mu-mu1
                        #print 'mu2=', mu2
                        if abs(mu2) > l-1: continue
                        cg1 = gn.getcgcoef(1,l-1,l,mu1,mu2,self.cgcoef)
                        #print '%3d %3d %3d %3d %12.7f' % (l, mu, mu1, mu2, cg1)
                        for m in range(l+1):
                            for m1 in [-1,0,1]:
                                m2 = m-m1
                                if abs(m2) > l-1: continue
                                dj1 = self.get_djmm(idf,1,mu1,m1)
                                dj2 = self.get_djmm(idf,l-1,mu2,m2)
                                cg2 = gn.getcgcoef(1,l-1,l,m1,m2,self.cgcoef)
                                #print '%3d %3d %3d %3d %3d  %12.7f %12.7f  %12.7f %12.7f  %12.7f' % (l, mu, mu1, m, m1, dj1.real, dj1.imag, dj2.real, dj2.imag, cg2)
                                _djmm_[l+mu,m] += cg1*cg2*dj1*dj2
                # pack the current results into array self.djmm
                for mu in range(-l,l+1):
                    for m in range(l+1):
                        idx = l*(l+1)*(4*l-1)/6+(l+1)*(mu+l)+m
                        self.djmm[idx,idf] = _djmm_[l+mu,m]*prefac
                        if Debug_Print:
                            print >> fout, ' djmm[idf=%2d,l=%2d,mu=%3d,m=%2d,idx=%4d]=%16.10f %16.10f' % (idf+1,l,mu,m,idx+1,self.djmm[idx,idf].real,self.djmm[idx,idf].imag)


def Check_Equal_k_lists(klist, kqm, fout):
    klst = array(klist)
    Equal_k_lists = True
    if len(klist) == len(kqm.kirlist):
        for ik in range(len(klist)):
            if sum(abs(klst[ik,:] - kqm.kirlist[ik,:]/float(kqm.LCM))) > 1e-5:
                Equal_k_lists = False
    else:
        Equal_k_lists = False
        
    if not Equal_k_lists:
        print >> fout, 'Num. irr. k-points=', len(kqm.kirlist), 'while in energy file irr k-points=', len(klist)
        print >> fout, 'ERROR: the irreducible k-mesh generated here is inconsistent with that used in WIEN2k!!!'
        print >> fout, ' - there may be some bugs in libbzint '
        print >> fout, ' - use KS-vectors in the full BZ (check the flag \'-sv\' in gap2_init) to avoid such problems'
        print >> fout, 'k-points from vector file:'
        for ik in range(len(klist)):
            print >> fout, klist[ik]
        print >> fout, 'k-points generated here:'
        for ik in range(len(kqm.kirlist)):
            print >> fout, kqm.kirlist[ik,:]/float(kqm.LCM)
        sys.exit(1)
    

def Convert2CubicPotential(Vlm, lmxc, strc):
    c_kub = zeros((11,11))
    c_kub[ 0, 0] = 1.0
    c_kub[ 3, 2] = 1.0
    c_kub[ 4, 0] = 0.5*sqrt(7./3.)
    c_kub[ 4, 4] = 0.5*sqrt(5./3.)
    c_kub[ 6, 0] = 0.5*sqrt(0.50)
    c_kub[ 6, 2] = 0.25*sqrt(11.0)
    c_kub[ 6, 4] =-0.5*sqrt(7.0/2.0)
    c_kub[ 6, 6] =-0.25*sqrt(5.0)
    c_kub[ 7, 2] = 0.5*sqrt(13./6.)
    c_kub[ 7, 6] = 0.5*sqrt(11./16.)
    c_kub[ 8, 0] = 0.125*sqrt(33.)
    c_kub[ 8, 4] = 0.25*sqrt(7./3.)
    c_kub[ 8, 8] = 0.125*sqrt(65./3.)
    c_kub[ 9, 2] = 0.25*sqrt(3.)
    c_kub[ 9, 4] = 0.5*sqrt(17./6.)
    c_kub[ 9, 6] =-0.25*sqrt(13.)
    c_kub[ 9, 8] =-0.5*sqrt(7./6.)
    c_kub[10, 0] = 0.125*sqrt(65./6.)
    c_kub[10, 2] = 0.125*sqrt(247./6.)
    c_kub[10, 4] =-0.25*sqrt(11./2.)
    c_kub[10, 6] = 0.0625*sqrt(19./3.)
    c_kub[10, 8] =-0.125*sqrt(187./6.)
    c_kub[10,10] =-0.0625*sqrt(85.)
    sq2 = sqrt(2.0)
    for iat in range(strc.nat):
        #print 'lmxc=', lmxc[iat]
        if strc.iatnr[iat]>0:  # this means this atom has cubic environment, and we will transform potential to cubic
            lxc = 0
            while lxc < len(lmxc[iat]):
                lx,mx = lmxc[iat][lxc]
                #print 'iat=', iat, 'lxc=', lxc, 'l=', lx, 'm=', mx
                if lx==0:
                    if mx == 0:
                        lxc += 1
                    else:
                        break
                elif lx==-3:
                    if mx == 2:
                        Vlm[iat][lxc,:] *= -1.0/sq2
                        lxc += 1
                    else:
                        break
                elif lx in (4,6,-7,-9):
                    la = abs(lx)
                    c1, c2 = c_kub[la,mx], c_kub[la,mx+4]
                    if mx == 0:
                        sq1 = 1.0
                    else:
                        sq1 = sq2
                    tmp  =  Vlm[iat][lxc,:] * c1 + Vlm[iat][lxc+1,:] * c2
                    Vlm[iat][lxc  ,:] =  tmp * (c1/sq1)
                    Vlm[iat][lxc+1,:] =  tmp * (c2/sq2)

                    #print 'c1=', c1, 'c2=', c2, 'l=', lx, 'm=', mx, 'V[0]=', Vlm[iat][lxc,-1], Vlm[iat][lxc+1,-1]
                    lxc=lxc+2
                elif lx in (8,10):
                    if mx == 0:
                        sq1 = 1.0
                    else:
                        sq1 = sq2
                    la = abs(lx)
                    c1, c2, c3 = c_kub[la,mx], c_kub[la,mx+4], c_kub[la,mx+8]
                    tmp = Vlm[iat][lxc,:] * c1 + Vlm[iat][lxc+1,:] * c2 + Vlm[iat][lxc+2,:] * c3
                    Vlm[iat][lxc  ,:] = tmp * c1/sq1
                    Vlm[iat][lxc+1,:] = tmp * c2/sq2
                    Vlm[iat][lxc+2,:] = tmp * c3/sq2
                    lxc += 3
                else:
                    break
            
                    
            
class KohnShamSystem:
    """ This class set up Kohn-Sham reference system, including core
        stat`<es, KS energies and orbitals, KS exchange-correlation potential 
    """
    PRINT = False
    def __init__(self, io, kqm, case, in1, strc, core, radf, nspin, fout):
        
        # Reading w2k energy files and its KS-eigenvalues
        spflag = ['up','dn'] if nspin==2 else ['']
        (self.klist, wegh, Ebnd, self.hsrws) = w2k.Read_energy_file(case+'.energy'+spflag[0], strc, fout)
        band_max = min(map(len,Ebnd))
        self.Ebnd = zeros( (nspin,len(Ebnd),band_max) )
        for ik in range(len(Ebnd)):
            self.Ebnd[0,ik,:] = Ebnd[ik][:band_max]
        if nspin==2:
            (self.klist, wegh, Ebnd, self.hsrws) = w2k.Read_energy_file(case+'.energy'+spflag[1], strc, fout)
            for ik in range(len(Ebnd)):
                self.Ebnd[1,ik,:] = Ebnd[ik][:band_max]
        # converting to Hartrees
        self.Ebnd *= Ry2H  # convert bands to Hartree
        # number of basis functions excluding local orbitals
        self.nv = array(self.hsrws)-in1.nlo_tot

        # Setting up several cutoff-parameters, which depend on KS-energies
        self.set_band_par(io['emax_pol'], io['emax_sc'], io['iop_core'], core, nspin, fout)

        # Recompute Fermi energy if neeeded
        if io['efermi'] >= 1e-2: # Recompute the Fermi energy
            #if False:
            #    (EF, Eg, evbm, ecbm, eDos) = self.Fermi(kqm.atet, kqm.wtet, core.nval, nspin)
            #else:
            (EF, Eg, evbm, ecbm, eDos) = mcmn.calc_Fermi(self.Ebnd[0], kqm.atet, kqm.wtet, core.nval, nspin)
            
        else:
            print >> fout, ' Use the Fermi energy from case.ingw'
            evbm = max( filter(lambda x: x < EF, self.Ebnd.ravel()) )
            ecbm = min( filter(lambda x: x > EF, self.Ebnd.ravel()) )
            Eg = ecbm - evbm
            eDos = sum([ft.dostet(EF, self.Ebnd, kqm.atet, kqm.wtet) for isp in range(nspin)])*2.0/nspin
        # Printing results of Fermi energy and gaps
        if Eg >= 0:
            print >> fout, '\n'+'-'*32+'\nFermi: Insulating, KS E_Fermi[eV]=%-12.6f Gap[eV]=%-12.6f  EVBM[eV]=%-12.6f  ECBM[eV]=%-12.6f' % (EF*H2eV, Eg*H2eV, evbm*H2eV, ecbm*H2eV)
        else:
            print >> fout, '\n'+'-'*32+'\nFermi: Metallic, KS E_Fermi[eV]=%-12.6f  DOS[E_f]=%-12.6f' % (EF*H2eV, eDos)
        print >> fout, '-'*32
        self.EFermi = EF
        # We want the Fermi energy to be at zero
        # Now we change the band energies so that EF=0
        self.Ebnd -= EF
        if len(core.corind)>0:
            nsp = len(core.eig_core)       # bug jul.7 2020
            nat = len(core.eig_core[0])    # bug jul.7 2020
            #nsp,nat,nc = shape(core.eig_core)
            for isp in range(nsp):
                for iat in range(nat):
                    core.eig_core[isp][iat][:] = array(core.eig_core[isp][iat][:])*Ry2H - EF
        EF = 0.0
        self.EF = EF
        self.Eg = Eg
        print >> fout, 'Set EF to ', self.EF

        if False:
            print >> fout, 'KS: bande'
            for ik in range(shape(self.Ebnd)[1]):
                for ie in range(shape(self.Ebnd)[2]):
                    print >> fout, '%3d %3d %18.12f' % (ik,ie,self.Ebnd[0,ik,ie])
        
        
        # Setting some extra cutoff-parameters
        self.set_band_par2(io['ibgw'], io['nbgw'], io['emingw'], io['emaxgw'], core, nspin, fout)
        
        # Calculating many radial integrals needed for evaluating the Coulomb repulsion
        self.core_valence_integrals(case, in1, strc, core, radf, nspin, fout)
        self.valence_valence_integrals(case, in1, strc, core, radf, nspin, fout)
        print >> fout, '\n'+'-'*32

    def VectorFileRead(self, case, strc, latgen, kqm, pw, fout):
        def Shift_to_1BZ(k):
            " This gives k_{1BZ}-k"
            if k==0:
                return 0
            elif k>0:
                return -int(floor(k))
            elif k<0:
                return int(ceil(-k))
            
        nsp, nkp, band_max = shape(self.Ebnd)
        (heads, all_Gs, self.all_As, all_Ek) = w2k.Read_vector_file(case, strc, fout)
        
        self.indgkir = []
        if self.PRINT:
            print >> fout, 'indgkir: only some G-vectors generated above are used in the vector file. indgkir is index to this subset of G-vectors'
        for irk in range(len(all_Ek)):
            diffE = sum(abs(all_Ek[irk][:band_max]*Ry2H-self.EFermi - self.Ebnd[0,irk,:]))
            k, kname, wgh, ios, n0, nb = heads[irk]
            diffk = sum(abs(k-self.klist[irk]))
            if diffk > 1e-3:
                print >> fout, 'WARNING it seems k-points from vector file and energy file are not compatible. diffk=', diffk
            if diffE > 1e-3:
                print >> fout, 'WARNING it seems band energies from vector file and energy file are not compatible. diffE=', diffE
            
            nG,n3 = shape(all_Gs[irk])
            ind = zeros(nG,dtype=int)
            for i in range(nG):
                iG = all_Gs[irk][i]          # reciprocal vector read from vector file
                ind[i] = pw.ig0[tuple(iG)]   # what is the index of this G in my gvec (previously generated fixed reciprocal mesh)?
                if self.PRINT:
                    print >> fout, 'irk=%3d i=%3d indgkir=%3d' % (irk+1, i+1, ind[i]+1)
            self.indgkir.append(ind)
            # zzk == all_As[irk][ib,hsrows]
            # zzkall(:,:,irk,isp)=zzk
        
        # In the next several lines, we will find permutation of the reciprocal vectors needed to transform
        # eigenvectors from an irreducible k-point mesh to generic all-k point mesh
        #Gbas = array(round_(kqm.gbas * kqm.aaa), dtype=int) # This transforms from lattice to cartesian coordinate system
        #Gbas = kqm.k2icartes   # This transforms from lattice to semi-cartesian coordinate system (or does not transform when ortho=False)
        
        ankp = kqm.ndiv[0]*kqm.ndiv[1]*kqm.ndiv[2]          # all k-points, reducible as well
        #print >> fout, 'Gbas=', Gbas
        if True: #self.PRINT:
            print >> fout, 'Umklap shifts for all k-points. (If any component is not 0 or 1, gap2 code would give different result)'
            for ik in range(ankp):
                irk = kqm.kii_ind[ik]
                #write(6,'(I4,1x,I4,1x,A,3I3,A,2x,A,3I3,A)') ikp, kpirind(ikp), '(', ikvec(:), ')', '(', irkvec(:), ')'
                #print >> fout, '%4d %4d' % (ik+1, irk+1), tuple(kqm.klist[ik,:]), tuple(kqm.kirlist0[irk,:])
                if kqm.k_ind[irk] != ik : # this k-point is rreducible
                    isym = kqm.iksym[ik]       # which group operation was used to get this momentum vector
                    jkvec = dot(kqm.kirlist0[irk], latgen.tizmat[isym])
                    #g0 = map(lambda x: int( (1-sign(x))/2. ), jkvec)        # is equal to 1 if jkvec is negative
                    g0 = map(lambda x: Shift_to_1BZ(x/kqm.LCM), jkvec)
                    g0 = dot(kqm.k2icartes, g0)
                    #if latgen.ortho or strc.lattic[:3]=='CXZ':
                    #    g0 = dot(Gbas, g0)
                    print >> fout, '%3d' % (ik+1,), '%3d'*3 % tuple(g0)

        if self.PRINT:
            print >> fout, 'vector_G'
        self.phase_arg = []
        self.indgk = []   # permutation index for all k-points
        for ik in range(ankp): # over all k-points, i.e, k_generic
            irk = kqm.kii_ind[ik]
            if kqm.k_ind[irk] == ik : # this k-point is irreducible, hence the order of reciprocal vectors stays the same
                self.indgk.append( self.indgkir[irk] )
                self.phase_arg.append( 0.0 )
            else:
                Gs = all_Gs[irk]
                isym = kqm.iksym[ik]       # which group operation was used to get this momentum vector
                # Notice that k_generic*timat = k_irr + delta_G*timat   or   k_generic = k_irr*timat + delta_G
                # where k_generic and k_irr are in the 1-BZ, while k_irr*timat is generically not,
                #    and needs delta_G, which is called g0 here.
                # Note that g0 is here computed in lattice coordinates when ortho=True,
                #    and needs to be converted to cartesian below.
                
                ## timat*(k+G) == k_irr + G_irr
                ## tizmat*(k0+G0) == k_irr0 + G_irr0
                jkvec = dot(kqm.kirlist0[irk], latgen.tizmat[isym])  # generick k-point, but possibly outside 1BZ. Note that this is done in integer (lattice) coordinates, not semi-cartesians.
                #g0 = map(lambda x: int( (1-sign(x))/2. ), jkvec)    # BUG WAS HERE: This shift delta_G to bring k-point into 1BZ
                g0 = map(lambda x: Shift_to_1BZ(x/kqm.LCM), jkvec)   # This shift delta_G to bring k-point into 1BZ, i.e., k_{1BZ}-k
                g0 = dot(kqm.k2icartes, g0)                          # now transforming to semi-cartesian coordinates (from lattice integer coordinates)
                tsymat = strc.timat[isym]  # symmetry operation that transforms k_generic[ik] = k_irr[irk]*timat + delta_G
                # Notice that this transforms the vector directly to cartesian coordinate system, not lattice system, like latgen.tizmat[isym]
                
                G_p_k = all_Gs[irk][:,:] + kqm.klist[ik,:]/float(kqm.LCM) # k+G in semi-cartesian
                self.phase_arg.append( dot(G_p_k, strc.tau[isym,:]) )     # positions are in the same system, i.e., semi-cartesian : (k+G)*tau[isym]
                #phase = exp(-2*pi*self.phase_arg[-1]*1j)
                
                ind = zeros(shape(Gs)[0],dtype=int)
                # We have <k+G| H |k+G> given in terms of G's in indgkir order (as read from vector file).
                # We know that  G + k_irr = G + (k_generic-delta_G)*timat. If we change the order of G's, such that
                # G_new*timat = G, than we have  (G_new*timat + (k_generic-delta_G)*timat)*r = (G_new+k_generic-delta_G)*timat*r
                # and because timat is symmetry operation timat*r should be equivalent to r.
                # Now we see that G_new-delta_G will work for generik k-vector in the same way as G works for 
                # irreducible k_irr. Notice that timat^2=1 hence G_new = G*timat, and we need G*timat-delta_G
                iG_transformed = dot(all_Gs[irk], tsymat)   # all G-vectors in cartesian corrdinate systems are transformed at once here. This is equivalent to: fortran_timat[isym].G
                iG_transformed -= g0                        # We need G_new - delta_G, possible umklap shift
                for i in range(shape(Gs)[0]):               # over all reciprocal vectors in the plane wave basis
                    iG = iG_transformed[i]                  # for generic k-point, this is the correct set of G's, i.e, fortran_timat[isym].(k+G) == G_{transformed}+k_1BZ
                    ind[i] = pw.ig0[tuple(iG)]              # this gives index of each G-vector, and hence gives the necessary permutation of the eigenvector
                    if self.PRINT:
                        print >> fout, '%s%3d '*6 % ('ik=', ik+1, 'irk=', irk+1, 'isym=', isym+1, 'indgkir=', self.indgkir[irk][i]+1, 'ig=', i+1, 'indgk=', ind[i]+1), '%s%12.7f' %('arg=', args[i])
                self.indgk.append(ind)                      # permutation of the eigenvectors, which should work for generic k-point
                #print 'rred ik=', ik, 'with irred ik=', irk, 'and shift g0=', g0, 'jkvec=', jkvec, 'isym=', isym, 'and timat=', (tsymat.T).tolist()

        #print 'indgk=', len(self.indgk), map(shape,self.indgk)
        print >> fout, '-'*32, '\n'

    def Vxc(self, case, in1, strc, radf, fout):
        (lmxc, Vxclm, self.ksxc, self.Vxcs) = w2k.Read_xc_file(case, strc, fout)
        for iat in range(strc.nat): Vxclm[iat] *= Ry2H
        self.Vxcs  *= Ry2H
        Convert2CubicPotential(Vxclm, lmxc, strc)
        
        maxlxc = max([len(lmxc[iat]) for iat in range(strc.nat)])
        # create fortran equivalent of lmxc
        self.lxcm_f = zeros(strc.nat, dtype=int)
        self.lmxc_f = zeros((2,maxlxc,strc.nat), dtype=int, order='F')  
        for iat in range(strc.nat):
            self.lxcm_f[iat] = len(lmxc[iat])    # bug jul.7 2020
            for lxc in range(len(lmxc[iat])):
                self.lmxc_f[:,lxc,iat] = lmxc[iat][lxc][:]
        
        # number of functions at each (iat,l)
        nrf = [[ 2+len(in1.nLO_at_ind[iat][l]) if l < len(in1.nLO_at_ind[iat]) else 2 for l in range(in1.nt)] for iat in range(strc.nat)]   
        nrmax = max(map(max,nrf))
        self.uxcu = zeros( (strc.nat, maxlxc, in1.nt**2, nrmax**2) )
        #print 'shape(self.uxcu)=', shape(self.uxcu)
        
        isp=0
        # computing matrix elements <l2,m2| V_{bl,m} |l1,m1>
        for iat in range(strc.nat):
            #npt = strc.nrpt[iat]
            #dh  = log(strc.rmt[iat]/strc.r0[iat])/(npt - 1)      # logarithmic step for the radial mesh
            #dd = exp(dh)
            #rx = strc.r0[iat]*dd**range(npt)
            rx, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
            
            ars = [[] for l in range(in1.nt)]
            for l in range(len(in1.nLO_at_ind[iat])):
                ars[l] = [radf.ul[isp,iat,l,:]/rx, radf.udot [isp,iat,l,:]/rx]
                ars[l] += [radf.ulo[isp,iat,l,ilo,:]/rx for ilo in in1.nLO_at_ind[iat][l]]
            for l in range(len(in1.nLO_at_ind[iat]),in1.nt):
                ars[l] = [radf.ul[isp,iat,l,:]/rx, radf.udot [isp,iat,l,:]/rx]
            
            brs = [[] for l in range(in1.nt)]
            for l in range(len(in1.nLO_at_ind[iat])):
                brs[l] = [radf.us[isp,iat,l,:]/rx, radf.usdot [isp,iat,l,:]/rx]
                brs[l] += [radf.uslo[isp,iat,l,ilo,:]/rx for ilo in in1.nLO_at_ind[iat][l]]
            for l in range(len(in1.nLO_at_ind[iat]),in1.nt):
                brs[l] = [radf.us[isp,iat,l,:]/rx, radf.usdot [isp,iat,l,:]/rx]

            #for l in range(in1.nt):
            #    print 'nrf=', len(ars[l]), nrf[iat][l]
            for lxc in range(len(lmxc[iat])):
                lx,mx = lmxc[iat][lxc]
                lb = abs(lx)
                for l2 in range(in1.nt):                        # right radial function
                    for l1 in range(in1.nt):                    # left radial function
                        if not( abs(lb-l2) <= l1  and l1 <= (lb+l2) ) and not( abs(lb-l1) <= l2 and l2 <= (lb+l1) ): continue
                        nr1 = len(ars[l1])
                        for ir2 in range(len(ars[l2])):  # how many functions we have at l2?
                            a3 = ars[l2][ir2] * Vxclm[iat][lxc]
                            b3 = brs[l2][ir2] * Vxclm[iat][lxc]
                            for ir1 in range(len(ars[l1])): # how many functions we have at l1?
                                #print 'iat=', iat, 'lxc=', lxc, 'l2=', l2, 'l1=', l1, 'ir2=', ir2, 'ir1=', ir1, 'res=', rr
                                self.uxcu[iat, lxc, l2*in1.nt+l1, ir2*nr1+ir1 ] = rd.rint13g(strc.rel, ars[l1][ir1], brs[l1][ir1], a3, b3, dh, npt, strc.r0[iat])
        if self.PRINT:
            print >> fout, '-'*32
            print >> fout, 'matrix elements of Vxc'
            print >> fout, '-'*32
            for iat in range(strc.nat):
                for lxc in range(len(lmxc[iat])):
                    lb = abs(lmxc[iat][lxc][0])
                    for l2 in range(in1.nt):
                        nr2 = nrf[iat][l2]
                        for l1 in range(in1.nt):
                            if not( abs(lb-l2) <= l1  and l1 <= (lb+l2) ) and not( abs(lb-l1) <= l2 and l2 <= (lb+l1) ): continue
                            nr1 = nrf[iat][l1]
                            for ir2 in range(nr2):
                                for ir1 in range(nr1):
                                    print >> fout, 'iat=%2d lxc=%2d l2=%2d l1=%2d ir2=%2d ir1=%2d uxcu=%14.10f' % (iat+1, lxc+1, l2, l1, ir2+1, ir1+1, self.uxcu[iat, lxc, l2*in1.nt+l1, ir2*nr1+ir1 ] )
                                
                                
    def set_band_par2(self, io_ibgw, io_nbgw, io_emingw, io_emaxgw, core, nspin, fout):
        nbnd = shape(self.Ebnd)[2]
        nkp = shape(self.Ebnd)[1]
        nomax_numin = [0,10000]                        # [index of the last valence band, index of the first conduction band]
        for isp in range(nspin):
            nocc_at_k = [len(filter(lambda x: x<0, self.Ebnd[isp,ik,:])) for ik in range(nkp)] # how many occuiped bands at each k-point
            nomax = max(nocc_at_k)-1                 # index of the last valence band
            numin = min(nocc_at_k)                   # index of the first conduction band
            nomax_numin[0] = max(nomax_numin[0],nomax)
            nomax_numin[1] = min(nomax_numin[1],numin)
        self.nomax_numin = nomax_numin
        print >> fout, ' Highest occupied band: ', self.nomax_numin[0]
        print >> fout, ' Lowest unoccupied band:', self.nomax_numin[1]
        # Set the total number of bands considered in the summation over states
        # for the calculations the exchange (x) and correlation self-energies
        self.nbands_x = self.ncg_x + self.nomax_numin[0]+1 # only occupied bands are needed
        self.nbands_c = self.ncg_c + self.nbmaxsc     # some unoccupied bands are also needed
        print >> fout, ' Number of bands considered in Sx(nbands_x):', self.nbands_x
        print >> fout, ' Number of bands considered in Sc(nbands_c):', self.nbands_c
        if io_ibgw < 0:
            nocc_at_k = [[len(filter(lambda x: x<io_emingw, self.Ebnd[isp,ik,:])) for ik in range(nkp)] for isp in range(nspin)]# how many bands below io['emingw'] at each k-point
            self.ibgw = min(map(min,nocc_at_k))
            #print 'ibgw=', self.ibgw,  'nocc_at_k=', nocc_at_k
        else:
            self.ibgw = io_ibgw
        if io_nbgw <= 0:
            nocc_at_k = [[len(filter(lambda x: x<io_emaxgw, self.Ebnd[isp,ik,:])) for ik in range(nkp)] for isp in range(nspin)]# how many bands below io['emingw'] at each k-point
            self.nbgw = max(map(max,nocc_at_k))
            #print [[filter(lambda x: x<io_emaxgw, self.Ebnd[isp,ik,:]) for ik in range(nkp)] for isp in range(nspin)]
            #print 'nbgw=', self.nbgw, 'nocc_at_k=', nocc_at_k
        else:
            self.nbgw = min(io_nbgw,nbnd)

        if self.ibgw > self.nomax_numin[1]:
            print >> fout, '*KohnShamSystem: WARNING - range of gw bands!! ibgw=', self.ibgw, 'numin=', self.nomax_numin[1]
            print >> fout, "*Now we will set ibgw to 1 "
            self.ibgw = 0
        
        if self.nbgw <= self.nomax_numin[0]:
            print >> fout, "*KohnShamSystem: WARNING - range of gw bands!! nbgw,nomax=", self.nbgw, self.nomax_numin[0]
            print >> fout, "*Now we will set nbgw to nbmax"
            self.nbgw = nbnd

        #print >> fout, ' Index of bands considered in GW calculations: (%d,%d)' %(self.ibgw,self.nbgw), 'in energy range',io_emingw, 'to', io_emaxgw, 'Hartree'

        nbandsgw = self.nbgw - self.ibgw
        nvelgw = core.nval - 2*self.ibgw
        print >> fout, ' Nr. of bands (nbmax)               :', nbnd
        print >> fout, ' Nr. of bands used in P (nbmaxpol)   :', self.nbmaxpol
        print >> fout, ' Nr. of bands for Sc(nbmaxsc)       :', self.nbmaxsc
        print >> fout, ' Nr. of gw bands (nbandsgw)         :', nbandsgw
        print >> fout, ' Range of GW bands (ibgw,nbgw)      :', self.ibgw,self.nbgw, 'with energy range', io_emingw, 'to', io_emaxgw, 'Hartree'
        print >> fout, ' Nr. val. electrons (nvel)          :', int(core.nval)
        print >> fout, ' Nr. val. electrons in GW (nvelgw)  :',int(nvelgw)
        print >> fout, ' Nr. core states(ncg)                     :', self.ncg
        print >> fout, ' Nr. core states for exchange selfx(ncg_x):', self.ncg_x
        print >> fout, ' Nr. core states for correlat selfc(ncg_c):', self.ncg_c
        print >> fout, ' Nr. core states for polariz  eps  (ncg_p):', self.ncg_p
        print >> fout, '-'*32
            
        
    def set_band_par(self, io_emax_pol, io_emax_sc, io_iop_core, core, nspin, fout):
        band_max = shape(self.Ebnd)[2] # this is max number of bands which are present at all k-points. Some k-points might have more bands, but not all points have them.
        # here we set nbmaxpol, which determines the number of bands used for polmat calculations
        if io_emax_pol < 0:
            self.nbmaxpol = band_max
        else:
            self.nbmaxpol = max([len(filter(lambda x: x < io_emax_pol, self.Ebnd[0,ik,:])) for ik in range(len(self.Ebnd[0]))])
            if nspin==2:
                nbmaxpol = max([len(filter(lambda x: x < io_emax_pol, self.Ebnd[1,ik,:])) for ik in range(len(self.Ebnd[1]))])
                self.nbmaxpol = max(self.nbmaxpol,nbmaxpol)
            
        if io_emax_sc < 0:
            self.nbmaxsc = band_max
        else:
            self.nbmaxsc = max([len(filter(lambda x: x < io_emax_sc, self.Ebnd[0,ik,:])) for ik in range(len(self.Ebnd[0]))])
            if nspin==2:
                nbmaxsc = max([len(filter(lambda x: x < io_emax_sc, self.Ebnd[1,ik,:])) for ik in range(len(self.Ebnd[1]))])
                self.nbmaxsc = max(self.nbmaxsc,nbmaxsc)
        
        #  Set the number of core states considered in the summation over states
        ncg = len(core.corind) # number of all core states
        self.ncg_x = ncg       # for exchange
        if io_iop_core == 0:
            self.ncg_c, self.ncg_p  = ncg, ncg
        elif io_iop_core == 1:
            self.ncg_p, self.ncg_c  = 0, ncg
        else:
            self.ncg_c, self.ncg_p = 0, 0
        self.ncg = ncg
        
    def core_valence_integrals(self, case, in1, strc, core, radf, nspin, fout):
        # momradintc
        ilocals={}
        lomaxp1 = shape(in1.nlo)[0]
        # iul_ucl   = [<   u^v_{l+1}| d/dr -l/r |u^c_l>, <   u^v_{l-1}| d/dr +(l+1)/r |u^c_l>]
        # iudl_ucl  = [<udot^v_{l+1}| d/dr -l/r |u^c_l>, <udot^v_{l-1}| d/dr +(l+1)/r |u^c_l>]
        # iulol_ucl = [<ulo_{l+1}| d/dr - l/r |u^c_l>, <ulo_{l-1}| d/dr+(l+1)/r |u^c_l>]
        # iucl_ul   = [< u^c_l | d/dr -(l-1)/r|   u^v_{l-1}>, < u^c_l | d/dr +(l+2)|   u^v_{l+1}>]
        # iucl_udl  = [< u^c_l | d/dr -(l-1)/r|udot^v_{l-1}>, < u^c_l | d/dr +(l+2)|udot^v_{l+1}>]
        # iucl_ulol = [< u^c_l | d/dr -(l-1)/r|   ulo_{l-1}>, < u^c_l | d/dr +(l+2)|   ulo_{l+1}>]
        # iucl_ucl  = [< u^c_l | d/dr -(l-1)/r|   u^c_{l-1}>, < u^c_l | d/dr +(l+2)/r |u^c_{l+1}>] 
        #
        # <u^v_{l+1}|d/dr  -  l/r|u^c_l> = Integrate[ {r*d/dr(u_{c,l}/r) - l   u_{c,l}/r }* u_{v,l+1} , {r,0,R_MT}] = Integrate[ (d/dr(u_{c,l}) - u_{c,l}/r - l   u_{c,l}/r) * u_{v,l+1}, {r,0,R_MT}]
        # <u^v_{l-1}|d/dr+(l+1)/r|u^c_l> = Integrate[ {r*d/dr(u_{c,l}/r)+(l+1) u_{c,l}/r }* u_{v,l-1} , {r,0,R_MT}] = Integrate[ (d/dr(u_{c,l}) - u_{c,l}/r+(l+1) u_{c,l}/r) * u_{v,l-1}, {r,0,R_MT}]
        self.iul_ucl  = [[[] for i in range(strc.nat)] for j in range(nspin)]
        self.iudl_ucl = [[[] for i in range(strc.nat)] for j in range(nspin)]
        self.iulol_ucl  = [[[] for i in range(strc.nat)] for j in range(nspin)]
        # <u^c_l|d/dr-(l-1)/r|u^v_{l-1}> = Integrate[ u_{c,l} *{ r*d/dr(u_{v,l-1}/r) -(l-1)*u_{v,l-1}/r }, {r,0,R_MT}] = Integrate[ u_{c,l} *{ d/dr(u_{v,l-1}) - u_{v,l-1}/r -(l-1)u_{v,l-1}/r } , {r,0,R_MT}]
        # <u^c_l|d/dr+(l+2)/r|u^v_{l+1}> = Integrate[ u_{c,l} *{ r*d/dr(u_{v,l+1}/r) +(l+2)*u_{v,l+1}/r }, {r,0,R_MT}] = Integrate[ u_{c,l} *{ d/dr(u_{v,l+1}) - u_{v,l+1}/r +(l+2)u_{v,l+1}/r } , {r,0,R_MT}]
        # Note that major (u_v) and minor component (u_{minor}) for valence electrons are related by derivative:
        #                           r*d/dr(u_v/r) = (d u_v/dr - u_v/r) = 2*u_{minor} 
        self.iucl_ul  = [[[] for i in range(strc.nat)] for j in range(nspin)]
        self.iucl_udl = [[[] for i in range(strc.nat)] for j in range(nspin)]
        self.iucl_ulol = [[[] for i in range(strc.nat)] for j in range(nspin)]
        # <u^c_l| d/dr - (l-1)/r |u^c_{l-1}> = Integrate[ u_{core,l} * { r d/dr (u_{core,l-1}/r) - (l-1) * u_{core,l-1}/r }, {r,0,R_MT}]
        # <u^c_l| d/dr + (l+2)/r |u^c_{l+1}> = Integrate[ u_{core,l} * { r d/dr (u_{core,l+1}/r) + (l+2) * u_{core,l+1}/r }, {r,0,R_MT}]
        for isp in range(nspin):
            for iat in range(strc.nat):
                rx, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
                _iul_ucl_   = zeros((2,len(core.l_core[iat])),order='F')
                _iudl_ucl_  = zeros((2,len(core.l_core[iat])),order='F')
                _iucl_ul_   = zeros((2,len(core.l_core[iat])),order='F')
                _iucl_udl_  = zeros((2,len(core.l_core[iat])),order='F')
                al_lo = reduce(lambda x,y: x+y, in1.nLO_at_ind[iat])
                nloat = 1
                if al_lo: nloat = max(al_lo) + 1
                _iulol_ucl_ = zeros((2,nloat,len(core.l_core[iat])),order='F')
                _iucl_ulol_ = zeros((2,nloat,len(core.l_core[iat])),order='F')
                #_iucl_ucl_ = zeros((2,len(core.l_core[iat]),len(core.l_core[iat])))
                for ic,lc in enumerate(core.l_core[iat]):
                    ucl  = core.ul_core[isp][iat][ic][:npt]
                    ucor = ucl/rx
                    ucpl = radd.derv(ucl,rx) - ucor  # ucpl == r d(u/r)/dr = du/dr - u/r
                    cfxs = [lc, -(lc+1)]
                    for dl in [1,-1]:
                        # Integrate[ {r*d/dr(u_{c,l}/r) -cfx * u_{c,l}/r }* u_{v,l+-1} , {r,0,R_MT}]
                        lv = lc + dl
                        il = (1-dl)/2
                        cfx = cfxs[il]
                        if lv >= 0:
                            #  iul_ucl  = [<   u^v_{l+1}| d/dr -l/r |u^c_l>, <   u^v_{l-1}| d/dr +(l+1)/r |u^c_l>]
                            #  iudl_ucl = [<udot^v_{l+1}| d/dr -l/r |u^c_l>, <udot^v_{l-1}| d/dr +(l+1)/r |u^c_l>]
                            #           = int(u_{l+1} d u_{core_l}/dr r^2 dr) - l * int(u_{l+1} u_{core_l} r dr)
                            ul   = radf.ul  [isp,iat,lv]
                            udl  = radf.udot[isp,iat,lv]
                            _iul_ucl_[il,ic]  = rd.rint13_nr(ul, ucpl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ul, ucor, dh, npt, strc.r0[iat])
                            #  iudl1ucl = <udot_{l+1}| d/dr -l/r | u^c_l> = int(udot_{l+1} d u_{core_l}/dr r^2 dr) - l * int(udot_{l+1} u_{core_l} r dr)
                            _iudl_ucl_[il,ic] = rd.rint13_nr(udl, ucpl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(udl, ucor, dh, npt, strc.r0[iat])
                            #  Now for local orbitals, l< lomax
                            if lc+1 < lomaxp1:
                                for ilo in in1.nLO_at_ind[iat][lv]:
                                    # iulol_ucl = [<ulo_{l+1}| d/dr-l/r |u^c_l>, <ulo_{l-1}| d/dr+(l+1)/r |u^c_l>]
                                    # iulol1ucl = int(ulo_{l+1} ducore_l/dr r^2 dr) - l*int(ulo_{l+1} ucore_l r dr)
                                    ulol = radf.ulo [isp,iat,lv,ilo]
                                    #ilocals[('iulol_ucl',il,isp,iat,ic,ilo)] = rd.rint13_nr(ulol, ucpl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ulol, ucor, dh, npt, strc.r0[iat])
                                    _iulol_ucl_[il,ilo,ic] = rd.rint13_nr(ulol, ucpl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ulol, ucor, dh, npt, strc.r0[iat])

                    cfxs = [lc-1, -(lc+2)]
                    for dl in [-1,1]:
                        # iucl_ul = [<u^c_l|d/dr -(l-1)| u^v_{l-1}>, <u^c_l|d/dr +(l+2)| u^v_{l+1}>] = Integrate[ u_{c,l} *{ r*d/dr(u_{v,l-+1}/r) -cfx * u_{v,l-+1}/r }, {r,0,R_MT}]
                        # iucl_udl= [<u^c_l|d/dr -(l-1)| udot^v_{l-1}>, <u^c_l|d/dr +(l+2)| udot^v_{l+1}>]
                        # iucl_ulol=[<u^c_l|d/dr -(l-1)|    ulo_{l-1}>, <u^c_l|d/dr +(l+2)|    ulo_{l+1}>]
                        lv = lc + dl
                        il = (dl+1)/2
                        cfx = cfxs[il]
                        if lv >= 0:
                            upl  = 2.0*radf.us   [isp,iat,lv,:npt]    # this is equivalent to r*d/dr u_l
                            udpl = 2.0*radf.usdot[isp,iat,lv,:npt]    # this is equivalent to r*d/dr udot_l
                            uor  =     radf.ul   [isp,iat,lv,:npt]/rx # this is u_l/r
                            udor =     radf.udot [isp,iat,lv,:npt]/rx # this is udot_l/r
                            #  iucl1ul = int( u_{core,l} d/dr u_{lv} r^2 dr) - cfx * int( u_{core_l} u_{lv} r dr)
                            _iucl_ul_[il,ic]  = rd.rint13_nr(ucl, upl,  dh, npt, strc.r0[iat]) - cfx*rd.rint13_nr(ucl, uor,  dh, npt, strc.r0[iat])
                            #  iul1udl = int( u_{core,l} d/dr udot_{lv} r^2 dr) - cfx * int( u_{core_l} udot_{lv} r dr)
                            _iucl_udl_[il,ic] = rd.rint13_nr(ucl, udpl, dh, npt, strc.r0[iat]) - cfx*rd.rint13_nr(ucl, udor, dh, npt, strc.r0[iat])
                            #  And now for local orbitals l< lomax
                            if lv < lomaxp1:
                                for ilo in in1.nLO_at_ind[iat][lv]:
                                    ulopl = 2.0*radf.uslo[isp,iat,lv,ilo,:npt]   # this is equivalent to d*d/dr u_{lo}
                                    uloor =     radf.ulo[isp,iat,lv,ilo,:npt]/rx # this is u_{lo}/r
                                    # iucl1ulol = int( u_{core,l} d/dr ulo_{lv} r^2 dr) - cfx * int( u_{core,l} ulo_{lv} r dr)
                                    #ilocals[('iucl_ulol',il,isp,iat,ic,ilo)] = rd.rint13_nr(ucl, ulopl, dh, npt, strc.r0[iat]) - cfx*rd.rint13_nr(ucl, uloor, dh, npt, strc.r0[iat])
                                    _iucl_ulol_[il,ilo,ic] = rd.rint13_nr(ucl, ulopl, dh, npt, strc.r0[iat]) - cfx*rd.rint13_nr(ucl, uloor, dh, npt, strc.r0[iat])
                    
                    # Now the core core elements
                    # iucl_ucl = [< u^c_l |d/dr-(l-1)/r|u^c_{l-1}>, < u^c_l |d/dr+(l+2)/r|u^c_{l+1}>] = Integrate[ u_{core,l} * { r d/dr (u_{core,l-+1}/r) - cfx * u_{core,l-+1}/r }, {r,0,R_MT}]
                    cfxs = [lc-1, -(lc+2)]
                    for jc,cl in enumerate(core.l_core[iat]):
                        for dl in [-1,1]:
                            if cl == lc + dl:
                                il = (dl+1)/2
                                cfx = cfxs[il]
                                ul  = core.ul_core[isp][iat][jc][:npt]
                                uor = ul/rx
                                upl = radd.derv(ul,rx) - uor  # upl == r d(u/r)/dr = du/dr - u/r
                                #_iucl_ucl_[il,ic,jc] = rd.rint13_nr(ucl, upl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ucl, uor, dh, npt, strc.r0[iat])
                                ilocals[('iucl_ucl',il,isp,iat,ic,jc)] = rd.rint13_nr(ucl, upl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ucl, uor, dh, npt, strc.r0[iat])
                self.iul_ucl  [isp][iat] = _iul_ucl_
                self.iudl_ucl [isp][iat] = _iudl_ucl_
                self.iulol_ucl[isp][iat] = _iulol_ucl_
                self.iucl_ul  [isp][iat] = _iucl_ul_
                self.iucl_udl [isp][iat] = _iucl_udl_
                self.iucl_ulol[isp][iat] = _iucl_ulol_
                #
        if Debug_Print:
            for isp in range(nspin):
                for iat in range(strc.nat):
                    for ic,lc in enumerate(core.l_core[iat]):
                        print >> fout, 'iul1ucl (%2d,%2d)=%10.7f' %(iat,ic,self.iul_ucl [isp][iat][0,ic])
                        print >> fout, 'iudl1ucl(%2d,%2d)=%10.7f' %(iat,ic,self.iudl_ucl[isp][iat][0,ic])
                        print >> fout, 'iucl1ul (%2d,%2d)=%10.7f' %(iat,ic,self.iucl_ul [isp][iat][0,ic])
                        print >> fout, 'iucl1udl(%2d,%2d)=%10.7f' %(iat,ic,self.iucl_udl[isp][iat][0,ic])
                        #     
                        print >> fout, 'iulucl1 (%2d,%2d)=%10.7f' %(iat,ic,self.iul_ucl [isp][iat][1,ic])
                        print >> fout, 'iudlucl1(%2d,%2d)=%10.7f' %(iat,ic,self.iudl_ucl[isp][iat][1,ic])
                        print >> fout, 'iuclul1 (%2d,%2d)=%10.7f' %(iat,ic,self.iucl_ul [isp][iat][1,ic])
                        print >> fout, 'iucludl1(%2d,%2d)=%10.7f' %(iat,ic,self.iucl_udl[isp][iat][1,ic])
                        #
                        for ilo in in1.nLO_at_ind[iat][lc+1]:
                            print >> fout, 'iulol1ucl(%2d,%2d,%2d)=%10.7f' % (iat,ic,ilo,self.iulol_ucl[isp][iat][0,ilo,ic])
                        for ilo in in1.nLO_at_ind[iat][lc-1]:
                            print >> fout,  'iulolucl1(%2d,%2d,%2d)=%10.7f' % (iat,ic,ilo,self.iulol_ucl[isp][iat][1,ilo,ic])
                        for ilo in in1.nLO_at_ind[iat][lc-1]:
                            print >> fout, 'iucl1ulol(%2d,%2d,%2d)=%10.7f' % (iat,ic,ilo,self.iucl_ulol[isp][iat][0,ilo,ic])
                        for ilo in in1.nLO_at_ind[iat][lc+1]:
                            print >> fout, 'iuclulol1(%2d,%2d,%2d)=%10.7f' % (iat,ic,ilo,self.iucl_ulol[isp][iat][1,ilo,ic])
                        for jc,cl in enumerate(core.l_core[iat]):
                            if ilocals.has_key(('iucl_ucl',0,isp,iat,ic,jc)):
                                print >> fout, 'iucl1ucl(%2d,%2d,%2d)=%10.7f' % (iat,ic,jc, ilocals[('iucl_ucl',0,isp,iat,ic,jc)])
                            if ilocals.has_key(('iucl_ucl',1,isp,iat,ic,jc)):
                                print >> fout, 'iuclucl1(%2d,%2d,%2d)=%10.7f' % (iat,ic,jc, ilocals[('iucl_ucl',1,isp,iat,ic,jc)])
                
                #iucl1ucl == iucl_ucl[0]
                #iuclucl1 == iucl_ucl[1]
                #iul1ucl  [isp][iat] = _iul_ucl_[0,:]
                #iudl1ucl [isp][iat] = _iudl_ucl_[0,:]
                #iulol1ucl[isp][iat] = _iulol_ucl_[0,:]
                #iulucl1  [isp][iat] = _iul_ucl_[1,:]
                #iudlucl1 [isp][iat] = _iudl_ucl_[1,:]
                #iulolucl1[isp][iat] = _iulol_ucl_[1,:]
                #iucl1ul  [isp][iat] = _iucl_ul_[0,:]
                #iucl1udl [isp][iat] = _iucl_udl_[0,:]
                #iucl1ulol[isp][iat] = _iucl_ulol_[0,:]
                #iuclul1  [isp][iat] = _iucl_ul_[1,:]
                #iucludl1 [isp][iat] = _iucl_udl_[1,:]
                #iuclulol1[isp][iat] = _iucl_ulol_[1,:]
                #call momradintv(iat,isp)
    def valence_valence_integrals(self, case, in1, strc, core, radf, nspin, fout):
        # momradintv
        self.iul_ul   = zeros((2,in1.nt-1,strc.nat,nspin),order='F')
        self.iul_udl  = zeros((2,in1.nt-1,strc.nat,nspin),order='F')
        self.iudl_ul  = zeros((2,in1.nt-1,strc.nat,nspin),order='F')
        self.iudl_udl = zeros((2,in1.nt-1,strc.nat,nspin),order='F')
             
        self.ilocals={}
        for isp in range(nspin):
            for iat in range(strc.nat):
                #npt = strc.nrpt[iat]
                #dh  = log(strc.rmt[iat]/strc.r0[iat])/(npt - 1)      # logarithmic step for the radial mesh
                #dd = exp(dh)
                #rx = strc.r0[iat]*dd**range(npt)
                rx, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
                for l in range(in1.nt-1):
                    l1s = [l+1,l]
                    l2s = [l,l+1]
                    cfxs = [l,-(l+2)]
                    for il in [0,1]:
                        l1, l2, cfx = l1s[il], l2s[il], cfxs[il]
                        ul   =   radf.ul   [isp,iat,l1]
                        udl  =   radf.udot [isp,iat,l1]
                        upl  = 2*radf.us   [isp,iat,l2] # note 2*us == d/dr ul(r) - ul(r)/r = r d/dr( ul(r)/r )
                        udpl = 2*radf.usdot[isp,iat,l2]
                        uor  =   radf.ul   [isp,iat,l2]/rx
                        udor =   radf.udot [isp,iat,l2]/rx
                        #  iul_ul  = < u_{l1}  | d/dr - cfx/r |  u_{l2} > = int(u_{l1} d/dr(u_l2) r^2 dr) - cfx*int(u_{l1} u_l2 r dr)
                        self.iul_ul[il,l,iat,isp]  = rd.rint13_nr(ul, upl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ul, uor, dh, npt, strc.r0[iat])
                        #  iul_udl = < u_{l1}  | d/dr - cfx/r |udot_{l2}> = int(u_{l1} d udot_l2/dr r^2 dr) - cfx*int(u_{l1} udot_l2 r dr)
                        self.iul_udl[il,l,iat,isp] = rd.rint13_nr(ul, udpl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ul, udor, dh, npt, strc.r0[iat])
                        #  iudl_ul = <udot_{l1}| d/dr - cfx/r |  u_{l2} > = int(udot_{l1} du_l2/dr r^2 dr) - cfx*int(udot_{l1} u_l2 r dr)
                        self.iudl_ul[il,l,iat,isp] = rd.rint13_nr(udl, upl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(udl, uor, dh, npt, strc.r0[iat])
                        #  iudl_udl= <udot_{l1}| d/dr - cfx/r |udot_{l2}> = int(udot_{l1} d udot_l2/dr r^2 dr) - cfx*int(udot_{l1} udot_l2 r dr)
                        self.iudl_udl[il,l,iat,isp]= rd.rint13_nr(udl, udpl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(udl, udor, dh, npt, strc.r0[iat])
                        for ilo in in1.nLO_at_ind[iat][l1]:
                            ulol = radf.ulo[isp,iat,l1,ilo]
                            #  iulol_ul = < ulo_{l1}| d/dr - cfx/r | u_{l2} > = int(ulo_{l1} du_l2/dr r^2 dr) - cfx*int(ulo_{l1} u_l2 r dr)
                            self.ilocals[('iulol_ul',il,isp,iat,l,ilo)] = rd.rint13_nr(ulol, upl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ulol, uor, dh, npt, strc.r0[iat])
                            #  iulol_udl= < ulo_{l1}| d/dr - cfx/r |udot_{l2}> = int(ulo_{l1} d udot_l2/dr r^2 dr) - cfx*int(ulo_{l1} udot_l2 r dr)
                            self.ilocals[('iulol_udl',il,isp,iat,l,ilo)]= rd.rint13_nr(ulol, udpl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ulol, udor, dh, npt, strc.r0[iat])
                        #  And now for local orbitals, only for l<= lomax
                        for jlo in in1.nLO_at_ind[iat][l2]:
                            ulopl = 2*radf.uslo[isp,iat,l2,jlo]
                            uloor = radf.ulo[isp,iat,l2,jlo]/rx
                            #  iul_ulol = < u_{l1}  | d/dr - cfx/r | ulo_{l2}> = int(u_{l1} d ulo_l2/dr r^2 dr) - cfx*int(u_{l1} ulo_l2 r dr)
                            self.ilocals[('iul_ulol',il,isp,iat,l,jlo)] = rd.rint13_nr(ul, ulopl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ul, uloor, dh, npt, strc.r0[iat])
                            # iudl_ulol = <udot_{l1}| d/dr - cfx/r | ulo_{l2}> = int(udot_{l1} d ulo_l2/dr r^2 dr) - cfx*int(udot_{l1} ulo_l2 r dr)
                            self.ilocals[('iudl_ulol',il,isp,iat,l,jlo)] = rd.rint13_nr(udl, ulopl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(udl, uloor, dh, npt, strc.r0[iat])
                            for ilo in in1.nLO_at_ind[iat][l1]:
                                ulol = radf.ulo[isp,iat,l1,ilo]
                                # iulol_ulol = < ulo_{l1,ilo}| d/dr - cfx/r | ulo_{l2,jlo}> = int(ulo_{l1} d ulo_l2/dr r^2 dr) - cfx*int(ulo_{l1} ulo_l2 r dr)
                                self.ilocals[('iulol_ulol',il,isp,iat,l,jlo,ilo)] = rd.rint13_nr(ulol, ulopl, dh, npt, strc.r0[iat]) - cfx * rd.rint13_nr(ulol, uloor, dh, npt, strc.r0[iat])


        # printing
        if Debug_Print:
            for isp in range(nspin):
                for iat in range(strc.nat):
                    rx, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
                    for l in range(in1.nt-1):
                        print >> fout, 'iul1ul  (%2d,%2d)=%10.7f' % (iat,l,self.iul_ul  [0,l,iat,isp]), ' iulul1  (%2d,%2d)=%10.7f' % (iat,l,self.iul_ul[1,l,iat,isp])
                        print >> fout, 'iul1udl (%2d,%2d)=%10.7f' % (iat,l,self.iul_udl [0,l,iat,isp]), ' iuludl1 (%2d,%2d)=%10.7f' % (iat,l,self.iul_udl[1,l,iat,isp])
                        print >> fout, 'iudl1ul (%2d,%2d)=%10.7f' % (iat,l,self.iudl_ul [0,l,iat,isp]), ' iudlul1 (%2d,%2d)=%10.7f' % (iat,l,self.iudl_ul[1,l,iat,isp])
                        print >> fout, 'iudl1udl(%2d,%2d)=%10.7f' % (iat,l,self.iudl_udl[0,l,iat,isp]), ' iudludl1(%2d,%2d)=%10.7f' % (iat,l,self.iudl_udl[1,l,iat,isp])
            # printing local orbitals
            for iky in self.ilocals.keys():
                if iky[1]==0:
                    name = re.sub('_','1',iky[0])
                else:
                    name = re.sub('_','',iky[0])+'1'
                print >> fout, '%-9s%s=%10.7f' % (name,str(tuple(iky[2:])),self.ilocals[iky])
        # DICTIONARY
        #iul1ul    ,iulul1    = iul_ul[0,1]
        #iul1udl   ,iuludl1   = iul_udl[0,1]
        #iudl1ul   ,iudlul1   = iudl_ul[0,1]
        #iudl1udl  ,iudludl1  = iudl_udl[0,1]
        #iulol1ul  ,iulolul1  = iulol_ul[0,1] 
        #iulol1udl ,iuloludl1 = iulol_udl[0,1]
        #iul1ulol  ,iululol1  = iul_ulol[0,1]
        #iudl1ulol ,iudlulol1 = iudl_ulol[0,1]
        #iulol1ulol,iulolulol1= iulol_ulol[0,1]
    
    def Give_fortran_ilocals(self, isp,strc,in1):
        #for l in range(in1.nt-1):
        #    l1s = [l+1,l]
        #    l2s = [l,l+1]
        #    for il in [0,1]:
        #        l1, l2 = l1s[il], l2s[il]
        #        for jlo in in1.nLO_at_ind[iat][l2]:
        #            for ilo in in1.nLO_at_ind[iat][l1]:
        #                print 'l=',l, 'il=',il, 'jlo=',jlo, 'ilo=',ilo, '     l1=', l1, 'l2=', l2
        #
        #iat=0
        #print 'nLO_at=', in1.nLO_at[:,iat]
        #lomax = shape(in1.nLO_at)[0]-1
        #print 'lomax=', lomax
        #for l in range(in1.nt-1):
        #    for jlo in range(in1.nLO_at[l,iat]):
        #        for ilo in range(in1.nLO_at[l+1,iat]):
        #            print 'need l=', l, 'il=', 0, 'ilo=', ilo+1, 'jlo=', jlo+1
        #            print 'need l=', l, 'il=', 1, 'jlo=', jlo+1, 'ilo=', ilo+1
        kys = self.ilocals.keys()
        lomax, nlomax = 1, 1
        if kys:
            lomax = max([kys[i][4] for i in range(len(kys))])+1
            nlomax = max([kys[i][5] for i in range(len(kys))])
            #print 'nlomax=', nlomax, 'lomax=', lomax

        iulol_ul   = zeros((2,nlomax,lomax,strc.nat), order='F')
        iulol_udl  = zeros((2,nlomax,lomax,strc.nat), order='F')
        iul_ulol   = zeros((2,nlomax,lomax,strc.nat), order='F')
        iudl_ulol  = zeros((2,nlomax,lomax,strc.nat), order='F')
        iulol_ulol = zeros((2,nlomax,nlomax,lomax,strc.nat), order='F')
        for k in self.ilocals.keys():
            if k[0]=='iul_ulol' and k[2]==isp:
                (il,_isp_,iat,l,jlo) = k[1:]
                iul_ulol[il,jlo-1,l,iat] = self.ilocals[k]
            elif k[0]=='iudl_ulol' and k[2]==isp:
                (il,_isp_,iat,l,jlo) = k[1:]
                iudl_ulol[il,jlo-1,l,iat] = self.ilocals[k]
            elif k[0]=='iulol_ul' and k[2]==isp:
                (il,_isp_,iat,l,jlo) = k[1:]
                iulol_ul[il,jlo-1,l,iat] = self.ilocals[k]
            elif k[0]=='iulol_udl' and k[2]==isp:
                (il,_isp_,iat,l,jlo) = k[1:]
                iulol_udl[il,jlo-1,l,iat] = self.ilocals[k]
            elif k[0]=='iulol_ulol' and k[2]==isp:
                (il,_isp_,iat,l,jlo,ilo) = k[1:]
                iulol_ulol[il,ilo-1,jlo-1,l,iat] = self.ilocals[k]
                #print 'setting il=', il, 'l=', l, 'ilo=', ilo, 'jlo=', jlo, ' = ', self.ilocals[k]
        return (iulol_ul, iulol_udl, iul_ulol, iudl_ulol, iulol_ulol)

            

class PlaneWaves:
    """Generates fixed basis of plane waves G, which will be used in subsequent calculations.
    """
    def __init__(self, hsrws, kmr, pwm, case, strc, in1, latgen, kqm, debug, fout):
        tm1 = timer()
        maxngk = max(hsrws)
        kxcmax = w2k.Read_xc_file_dimension(case, strc, fout)
        
        kmax = in1.rkmax/min(strc.rmt)
        ng = zeros(3,dtype=int)
        for i in range(3):
            bleng = sqrt(sum(latgen.br2[i,:]**2))
            ng[i] = int(kmr * kmax * pwm/bleng) + 1
        
        print >> fout, 'kxcmax=', kxcmax, 'ng1,ng2,ng3=', ng
        igm = array([max(4*ng[i],2*kxcmax+ng[i]) for i in range(3)],dtype=int)
        
        ng = map(int, igm/2)
        npw = (2*ng[0]+1)*(2*ng[1]+1)*(2*ng[2]+1)
        
        maxgxc = 2.*kmax+4
        maxgcoul = kmr*kmax*(pwm+2.)
        gmax = max(maxgxc, maxgcoul)
        
        tm2 = timer()
        print >> fout, '## PlaneWave t(read_xc)            =%14.9f' % (tm2-tm1,)
        print >> fout, 'gmax=', gmax, 'br2=', latgen.br2
        
        ortho = (latgen.ortho or strc.lattice[1:3]=='CXZ')
        if FORT:
            #glen=[]    # |G|
            #gindex=[]  #  G in terms of (i0,i1,i2)
            #G_c=[]     #  G in cartesian
            self.npw2apw, ngindx = fkp.count_number_pw(ng,kmax,gmax,latgen.br2)
            is_CXZ = (strc.lattice[1:3]=='CXZ')
            self.glen, self.G_c, self.gindex = fkp.generate_pw(ng,ngindx,gmax,latgen.pia,latgen.br2,ortho,is_CXZ,True)
        else:
            glen=[]    # |G|
            gindex=[]  #  G in terms of (i0,i1,i2)
            G_c=[]     #  G in cartesian
            self.npw2apw = 0
            for ii,ki in enumerate(itertools.product(range(-ng[0],ng[0]+1),range(-ng[1],ng[1]+1),range(-ng[2],ng[2]+1))):
                kc = dot(latgen.br2,ki)
                kk = sqrt(sum(kc**2))
                if kk < 2.*kmax+4: self.npw2apw += 1
                if kk < gmax :
                    glen.append( kk )
                    if ortho:
                        gindex.append( map(int, round_(kc/latgen.pia)) )
                    else:
                        if strc.lattice[1:3]=='CXZ':
                            gindex.append( [ki[0]+ki[2],ki[1],ki[2]-ki[0]] )
                        else:
                            gindex.append( ki )
                    G_c.append( kc )
            G_c = array(G_c)
            gindex = array(gindex,dtype=int)
            
            # Now we will sort plane waves according to their length
            indx = argsort(glen, kind='stable')  # just obtaining index to the sorted sequence
            
            # here we rearange all arrays so that they are in sorted order
            self.glen = zeros(shape(glen))           # |G|
            self.G_c  = zeros(shape(G_c))            # \vG in cartesian
            self.gindex = zeros(shape(gindex),dtype=int) # \vG in integer
            #self.ig0={}
            for i0,i in enumerate(indx):
                self.gindex[i0,:] = gindex[i]
                self.glen[i0]     = glen[i]
                self.G_c[i0]      = G_c[i]
                #iG = gindex[i]
                #self.ig0[tuple(iG)] = i0
        
        
        tm3 = timer()
        print >> fout, '## PlaneWave t(gen_ipw_fixed_basis)=%14.9f' % (tm3-tm2,)
        #####
        if False:
            ft = open('python_sorted_index.dat', 'w')
            for i in range(len(self.glen)):
                print >> ft, '%4d  %4d%4d%4d  %15.10f  %15.10f%15.10f%15.10f' % (i+1, self.gindex[i,0], self.gindex[i,1], self.gindex[i,2], self.glen[i], self.G_c[i,0], self.G_c[i,1], self.G_c[i,2])
            ft.close()
            
        # This is just temporary to debug the code, we want to be completely compatible with the fortran cdoe.
        if False:  ##### ATTENTION ????
            fix = open('sorted_index.dat')
            self.glen=[]
            self.gindex=[]
            self.G_c=[]
            for ii,line in enumerate(fix):
                iG = map(int,line.split()[1:4])
                Gcc = iG * latgen.pia
                kk = sqrt(sum(Gcc**2))
                self.glen.append( kk )
                self.gindex.append( iG )
                self.G_c.append( Gcc )
            self.glen = array(self.glen)
            self.G_c = array(self.G_c)
            self.gindex = array(self.gindex, dtype=int)
            
        ##### 
        self.ig0={}
        for i in range(ngindx):
            self.ig0[ tuple(self.gindex[i]) ] = i
        
        # so that we can keep using short glen, G_c and gindex
        glen = self.glen      # |G|
        G_c  = self.G_c       # \vG in cartesian
        gindex = self.gindex  # \vG in integer
        
        if Debug_Print:
            print >> fout, 'Sorted Gs'
            for i0 in range(len(gindex)):
                print >> fout, '%4d ' % (i0+1,), '(%3d,%3d,%3d)' % tuple(self.gindex[i0]), '%9.5f' % (self.glen[i0],), '[%9.5f,%9.5f,%9.5f]' % tuple(self.G_c[i0])
        
        tm4 = timer()
        print >> fout, '## PlaneWave t(sort_ipw)           =%14.9f' % (tm4-tm3,)
        
        # This part calculates the integral of a plane wave with wave vector k
        #     1/V * Integral[ e^{i*k*r} , {interstitial}]
        # belonging to the reciprocal Bravais lattice in the
        # interstitial region by the difference between the integral over the whole
        # unit cell and the Muffin Tin spheres
        # First, get vmt == the volume ration for each atom type
        self.vmt = array([4*pi*strc.rmt[iat]**3/3. * 1/latgen.Vol for iat in range(strc.nat)])
        mult = array(strc.mult)
        Rmt = array(strc.rmt)
        if FORT:
            self.ipwint = fkp.pw_integ(gindex,glen,self.vmt,Rmt,strc.vpos,mult)
        else:
            vmt_3 = self.vmt*3.0
            #
            self.ipwint = zeros(len(gindex),dtype=complex) # result= Integrate[e^{i*k*r},{r in interstitials}]_{cell,interstitial}/V_{cell}
            for i in range(len(gindex)):
                ki, ak = gindex[i], glen[i]  # ki=(i0,i1,i2); ak=|kc|
                if ak < 1e-10:
                    self.ipwint[i] = 1 - sum(self.vmt*mult)
                else:
                    kr  = ak * Rmt                       # kr[iat] = |kc|*Rmt[iat]
                    j1 = special.spherical_jn(1,kr)      # spherical_bessel == (sin(x)/x-cos(x))/x
                    intmod = vmt_3 * j1 / kr             # intmod[iat] = j1*3*vmt/kr
                    integc = 0j
                    for iat in range(strc.nat):
                        ekr = exp(2*pi*1j*dot(array(strc.pos[iat]),ki))  # ekr[ieq]
                        integc += intmod[iat] * sum(ekr)                 # sum(ekr)==sum(phase over all equivalent atoms)
                    self.ipwint[i] = -integc
        if debug:
            print >> fout, 'ipwint: '
            for i in range(len(gindex)):
                ki, ak = gindex[i], glen[i]  # ki=(i0,i1,i2); ak=|kc|
                print >> fout, '%4d ' % (i+1,), '%3d'*3 % tuple(ki), ' %16.11f%16.11f' %(self.ipwint[i].real, self.ipwint[i].imag)
        
        tm5 = timer()
        print >> fout, '## PlaneWave t(integral_ipw)       =%14.9f' % (tm5-tm4,)
        # First determine the maximum number of plane waves in variuos parts of the calculation
        maxlen_mb   = kmr * kmax      # |G+q| should be smaller than maxlen_mb for Mixbasis (maxngq)
        maxlen_coul = maxlen_mb * pwm # |G+q| should be smaller than maxlen_coul for bare Coulomb
        q_c = dot(kqm.qlist, latgen.br2.T)/float(kqm.LCMq)  # q in cartesian coordinates
        if FORT: # faster code written in Fortran
            self.ngq,self.ngq_barc,self.ngqlen = fkp.ngq_size(q_c,G_c,maxlen_mb,maxlen_coul)
        else:    # slower but equivalent Python code
            self.ngq = zeros(len(kqm.qlist),dtype=int) # plane waves for Mixbasis
            self.ngq_barc = zeros(len(kqm.qlist),dtype=int) # plane waves for bare Coulomb
            for iq in range(len(kqm.qlist)):
                #print '%3d ' % (iq+1,), ('%15.10f'*3) % tuple(q_c[iq]), ' %3d%3d%3d' % tuple(kqm.qlist[iq,:])
                kpq = [linalg.norm(G_c[i,:] + q_c[iq,:]) for i in range(len(G_c))]
                self.ngq[iq] = sum(1 for akq in kpq if akq<maxlen_mb)
                self.ngq_barc[iq] = sum(1 for akq in kpq if akq<maxlen_coul)
        maxngq = max(self.ngq)               # max number of plane waves for Mixbasis across all q-points
        maxngq_barc = max(self.ngq_barc)     # max number of plane waves for bare Coulomb across all q-points
        if maxngq == len(gindex):
            print >> fout, 'WARNING !! maxngq = npw !!!'
        if maxngq_barc == len(gindex):
            print >> fout, 'WARNING!! maxngq_barc = npw !!!'
        
        if FORT:  # faster code written in Fortran
            maxngqlen   = max(self.ngqlen)
            self.indgq, self.gqlen, self.G_unique = fkp.ngq_sort(q_c,G_c,maxlen_coul,self.ngqlen,maxngq_barc,maxngqlen)
        else:     # slower but equivalent Python code
            self.indgq    = zeros((len(kqm.qlist),maxngq_barc),dtype=int)  # index to the plane wave G for which |G+q| is beyond certain cutoff for bare Coulomb
            self.gqlen = zeros((len(kqm.qlist),maxngq_barc))               # the value of |G+q| for plane waves beyond certain cutoff for bare Coulomb
            G_unique_ = [[] for iq in range(len(kqm.qlist))]               # many of |G+q| are equal, and we want to store only unique terms. This is just index to the first occurence with unique length
            
            Gpq_len = zeros(maxngq_barc)
            which_ik= zeros(maxngq_barc, dtype=int)
            for iq in range(len(kqm.qlist)):
                kpq = [linalg.norm(G_c[i,:] + q_c[iq,:]) for i in range(len(G_c))]
                n = 0
                for i in range(len(G_c)):
                    if (kpq[i] < maxlen_coul):
                        Gpq_len[n] = kpq[i]
                        which_ik[n] = i
                        n += 1
                indxc = argsort(Gpq_len[:n], kind='stable')
                aG_prev = -1
                for i0,i in enumerate(indxc):
                    aG = Gpq_len[i]
                    self.indgq[iq,i0] = which_ik[i]
                    self.gqlen[iq,i0] = aG
                    if abs(aG-aG_prev)>1e-6:
                        G_unique_[iq].append(i0)
                        aG_prev = aG
            self.ngqlen = [len(G_unique_[iq]) for iq in range(len(kqm.qlist))]
            maxngqlen   = max(self.ngqlen)
            # Now creating array from list of lists for G_unique, to be compatible with above Fortran code
            self.G_unique = zeros((len(G_unique_),maxngqlen),dtype=int)
            for iq in range(len(G_unique_)):
                self.G_unique[iq,:self.ngqlen[iq]] = G_unique_[iq][:]



        tm7 = timer()
        print >> fout, '## PlaneWave t(max |k+q|)          =%14.9f' % (tm7-tm5,)

        print >> fout, 'Nr. of IPW(npw):',npw
        print >> fout, 'Nr. of pws for lapw overlaps (npw2apw):',self.npw2apw
        print >> fout, 'Max nr. of IPW for APW (maxngk):',maxngk
        print >> fout, 'Max nr. of IPW for Mixbasis (maxngq) :',maxngq
        print >> fout, 'Max nr. of IPW for bare Coulomb (maxngq_barc) :',maxngq_barc

        
        print >> fout, 'maxlen_mb=', maxlen_mb, ' maxlen_coul=', maxlen_coul
        print >> fout, 'gqlen and G_unique'
        for iq in range(len(kqm.qlist)):
            print >> fout, '%5d' %(iq+1,), 'q=', q_c[iq,:], 'ngq=%5d ngq_barc=%5d' % (self.ngq[iq], self.ngq_barc[iq])
            if debug:
                for i in range(self.ngq_barc[iq]):
                    ii = self.indgq[iq,i]
                    print >> fout, '   %5d %15.10f' % (i+1,self.gqlen[iq,i]), '%4d' % (ii+1,), '%15.10f'*3 % tuple(G_c[ii,:]), gindex[ii,:]
            if Debug_Print:
                for i in range(self.ngqlen[iq]):
                    i0 = self.G_unique[iq,i]
                    ii = self.indgq[iq,i0]
                    print >> fout, '    --- %3d %15.10f' % (i0, self.gqlen[iq,i0]), ii, gindex[ii,:], G_c[ii,:]
                
    def convert_ig0_2_array(self):
        ivects = self.ig0.keys()
        avects = array( ivects, dtype=int)
        mx, my, mz = max(abs(avects[:,0])), max(abs(avects[:,1])), max(abs(avects[:,2]))
        iag0 = -ones((2*mx+1,2*my+1,2*mz+1), dtype=int)
        for ik in ivects:
            iag0[ik[0]+mx, ik[1]+my, ik[2]+mz] = self.ig0[ik]
        return iag0

    def inverse_indgq(self, iq):
        # this is indggq
        max_len = max(self.indgq[iq,:self.ngq_barc[iq]])+1
        indgq_inverse = -ones(max_len, dtype=int)
        for ibasis in range(self.ngq_barc[iq]):
            ig = self.indgq[iq,ibasis]   # ig == index in gvec
            indgq_inverse[ig] = ibasis # if we have index in gvec, we can index in basis in mesh up to ngq_barc
        return indgq_inverse
    
class FrequencyMesh:
    def __init__(self, iopfreq, nomeg, omeg_min, omeg_max, iopMultiple, fout):
        #iopMultiple=20
        fginfo = ['Equally spaced mesh', 'Grid for Gauss-Laguerre quadrature', 'Grid for double Gauss-Legendre quadrature,', 'Grid of Tan-mesh for convolution', 'Using SVD basis and Tan-mesh for convolution']
        self.iopfreq = iopfreq
        print >> fout, 'Frequnecy grid for convolution: iopfreq='+str(iopfreq)+' om_max='+str(omeg_max)+' om_min='+str(omeg_min)+' nom='+str(nomeg)+': ', fginfo[iopfreq-1]
        if iopfreq == 1:     # Equaly spaced mesh (for tests or for emac calcultions, not for integration)
            self.omega = linspace(omeg_min, omeg_max, nomeg)
            self.womeg = zeros(nomeg)
        elif iopfreq == 2:   # Grid for Gauss-Laguerre quadrature
            self.omega, wu = laguerre.laggauss(nomeg)
            self.womeg = wu * exp(self.omega)
        elif iopfreq == 3:   # Double Gauss-Legendre quadrature from 0 to omegmax and from omegmax to infinity 
            n = nomeg/2
            u, wu = legendre.leggauss(n)
            self.omega, self.womeg = zeros(2*n), zeros(2*n)
            self.omega[:n] = 0.5*omeg_max*(u+1)
            self.womeg[:n] = 0.5*omeg_max*wu
            u  = u[::-1]   # turn them around
            wu = wu[::-1]  # turn them around
            self.omega[n:] = 2*omeg_max/(u+1)
            self.womeg[n:] = 2*omeg_max*wu/(u+1)**2
        elif iopfreq == 4:   # tan mesh
            # good combination : omeg_max = 20, omeg_min = 0.02, nomeg=32
            om0, dom0 = mcmn.Give2TanMesh(omeg_min, omeg_max, nomeg)
            n = len(om0)/2
            self.omega, self.womeg = om0[n:], dom0[n:]
            # another iopMultiple-times more precise mesh for self-energy integration
            minx = omeg_min/(1 + 0.2*(iopMultiple-1.))
            om1, dom1 = mcmn.Give2TanMesh(minx, omeg_max, nomeg*iopMultiple)
            n = len(om1)/2
            self.omega_precise, self.womeg_precise = om1[n:], dom1[n:]
            #self.omega_precise, self.womeg_precise = self.omega, self.womeg
            #print >> fout, 'Frequnecy with tan mesh for convolution nom=', len(self.omega)
        elif iopfreq == 5:
            # good combination : omeg_max = 20, omeg_min = 0.02, nomeg=32
            om0, dom0 = mcmn.Give2TanMesh(omeg_min, omeg_max, nomeg)
            n = len(om0)/2
            self.omega, self.womeg = om0[n:], dom0[n:]
        else:
          print >> fout, 'ERROR: wrong iopfreq=', iopfreq, ' should be between 1...3'
          
        #print >> fout, 'frequencies and weights'
        for i in range(len(self.omega)):
            print >> fout, '%3d  x_i=%16.10f w_i=%16.10f' % (i+1, self.omega[i], self.womeg[i])
            
class Kweights:
    def __init__(self, io, ks, kqm, fout):
        # setkiw(io, ks, kqm, fout):
        """Computes tetrahedron weights for Green's function like objects
        """
        wkir = array(kqm.weight)
        (nspin,nkir,nbnd) = shape(ks.Ebnd)
        ankp = kqm.ndiv[0]*kqm.ndiv[1]*kqm.ndiv[2]
        _ankp_ = float(ankp)
        
        if io.iop_bzint == 0:   # perform BZ integration over tetrahedra in the IBZ
            isp=0
            self.kiw = ft.intw(ks.Ebnd[isp], ks.EF, kqm.atet, kqm.wtet, io.iop_bcor==1)
            self.kwfer = zeros(shape(self.kiw))
            if (ks.Eg <= 0): # metallic
                self.kwfer = ft.intwsurf(ks.Ebnd[isp], ks.EF, kqm.atet, kqm.wtet)
            
        elif io.iop_bzint == 1:  # Set the weights by the smearing
            (ns,nkir,nbnd) = shape(ks.Ebnd)
            isp=0
            self.kwfer = zeros((nkir,nbnd))
            self.kiw   = zeros((nkir,nbnd))
            wkir_nkp = wkir/_ankp_
    
            self.kiw = cfermi(ks.Ebnd[isp].flatten()/io.esmear).reshape((nkir,nbnd))
            if (ks.Eg <= 0): # metallic
                self.kwfer = cgauss(ks.Ebnd[isp].flatten(),io.esmear).reshape((nkir,nbnd))
    
            # Finding bands which are either not fully empty (nomx) or not fully occupied (numn)
            nomx, numn = ks.nomax_numin
            nomx = max([len(filter(lambda x: x>1e-4,   self.kiw[ik,:])) for ik in range(len(self.kiw))]) # band at which the weight is very small
            numn = min([nbnd-len(filter(lambda x: x<1-1e-4, self.kiw[ik,:])) for ik in range(len(self.kiw))]) # band which is not fully occupied, and weight deviates from 1 for tiny amount
                
            for ik in range(nkir):
                self.kiw  [ik,:] *= wkir_nkp[ik]
                self.kwfer[ik,:] *= wkir_nkp[ik]
                
            print >> fout, 'nomx=', nomx, 'numn=', numn
            if nomx > ks.nomax_numin[0]:
                print >> fout, " nomax in terms of the occupation:", nomx
                ks.nomax_numin[0] = nomx
            
            if numn < ks.nomax_numin[1]:
                print >> fout, " numin in terms of the occupation:", numn 
                ks.nomax_numin[1] = numin
        else:
            print >> fout,  "ERROR: unsupported option iop_bzint=", io.iop_bzint
            sys.exit(1)
            
        
        print >> fout, 'tetrahedron weight: kiw kwfer:'
        for ik in range(nkir):
            for ib in range(nbnd):
                if self.kiw[ik,ib]!=0 or self.kwfer[ik,ib]!=0:
                    print >> fout, '%3d %3d %16.10f   %16.10f' % (ik+1, ib+1, self.kiw[ik,ib], self.kwfer[ik,ib])
    
        # set up the band-independent k-weight in the IBZ                
        self.kwt_ibz = self.kiw[:,0] # should be equal to wkir/ankp
        
        # kiw and kwfer correspond to the weight for an irreducible point.
        # If we work with all (even reducible) k-points, the weight has to be divided by wkir
        for ik in range(nkir):            
            self.kiw  [ik,:] *= 1.0/wkir[ik]
            self.kwfer[ik,:] *= 1.0/wkir[ik]
            
        # set up the band-independent k-weight in the RBZ
        kwt_bz = zeros(ankp)
        for ik in range(ankp):
            irk = kqm.kii_ind[ik]
            kwt_bz[ik] = self.kwt_ibz[irk]  # should be equal to 1/ankp
        
        print >> fout, 'tetrahedron weight: kiw kwfer:'
        for ik in range(nkir):
            for ib in range(nbnd):
                if self.kiw[ik,ib]!=0 or self.kwfer[ik,ib]!=0:
                    print >> fout, '%3d %3d %16.10f   %16.10f' % (ik+1, ib+1, self.kiw[ik,ib], self.kwfer[ik,ib])
        
        #print >> fout, 'k-points weights for all points'
        #for ik in range(ankp):
        #    print >> fout, '%3d %10.6f' % (ik, kwt_bz[ik])
    
#def ImproveDegEigenvectors(eigvals, sgi):
#    smalle = 1e-10
#    N = len(eigvals)
#    its=0
#    while its < N:
#        ite =  its+1
#        while (ite<N and abs(eigvals[ite]-eigvals[its])<smalle): ite += 1 # looping through almost equal eigenvalues
#        # In case of exactly two equal eigenvalues, we can make many entries exactly vanish. We will find the best linear combination of the two
#        # eigenvectors, such that most entries in one eigenvector are zero
#        if (ite-its == 2 ):
#            #print 'deg=', its,ite, eigvals[its:ite]
#            rt = [sgi[i,(its+1)]/sgi[i,its] if abs(sgi[i,its])>1e-50 else rand() for i in range(N)]
#            #rt = sgi[:,(its+1)]/sgi[:,its] # ratio between all components of the two eigenvectors
#            indx = numpy.argsort(rt)       # sort all ratio components
#            #print 'sorted=', [rt[indx[i]] for i in range(len(indx))]
#            # find which ratio occurs most often, i.e., has largest degeneracy in ratio vector
#            max_deg=1
#            max_pair=()
#            jts = 0
#            while jts < N:
#                jte = jts+1
#                while (jte < N and abs(rt[indx[jte]]-rt[indx[jts]])<1e-6 ): jte += 1 # loop over all components of the eigenvector, which have degenerate ratio
#                if (jte-jts > max_deg): 
#                    max_deg = jte-jts   # this is the largest degeneracy up to now, with max_deg components degenrate
#                    #print 'num=', jte-jts, jts,jte, [rt[indx[j]] for j in range(jts,jte)]
#                    max_pair = (jts,jte)
#                jts=jte
#            #print 'max_pair=', max_pair, sorted([indx[i] for i in range(max_pair[0],max_pair[1])])
#            if max_pair: # we found degeneracy, and is saved in max_pair. Degenrate components are: sorted([indx[i] for i in range(max_pair[0],max_pair[1])])
#                r = rt[indx[max_pair[0]]] # this is the ratio of all those degenrate components
#                a, b = r/sqrt(1+r**2), 1.0/sqrt(1+r**2) # coefficients of the linear combination
#                #print 'r=', r, 'a=', a, 'b=', b
#                vits=  b*sgi[:,its] + a*sgi[:,its+1]
#                vit = -a*sgi[:,its] + b*sgi[:,its+1]
#                #print 'old[it] =', ('%10.6f '*N) % tuple(sgi[:,its+1].real)
#                #print 'old[its]=', ('%10.6f '*N) % tuple(sgi[:,its].real)
#                #print 'new[it] =', ('%10.6f '*N) % tuple(vit.real)
#                #print 'new[its]=', ('%10.6f '*N) % tuple(vits.real)
#                #print 'orthogonality=', [dot(vits,vits), dot(vits,vit), dot(vit,vits), dot(vit,vit)]
#                sgi[:,its]   = vits
#                sgi[:,its+1] = vit
#        its = ite
#    return sgi

    
class MatrixElements2Band:
    DEBUG = False
    def __init__(self, io, pw, strc, latgen):
        # Instead of dictionary pw.ig0, which converts integer-G_vector into index in fixed basis, we will use fortran compatible integer array
        self.i_g0 = pw.convert_ig0_2_array() # Instead of dictionary pw.ig0, which converts integer-G_vector into index in fixed basis, we will use fortran compatible integer array
        
        alat = array([strc.a, strc.b, strc.c])
        self.calc_tildeg(2*(io.lmbmax+1))
        self.eta = self.optimal_eta(latgen.rbas, latgen.br2)
        self.r_cf = 2 * self.r_cutoff(io.stctol, self.eta, 10, alat)
        self.g_cf = 2 * self.g_cutoff(io.stctol, self.eta, 10, latgen.pia, latgen.Vol)

    
    def Vxc(self, strc, in1, latgen, kqm, ks, radf, pw, pb, fout):
        t_lapwcoef, t_lapack, t_muffint, t_interst = 0, 0, 0, 0
        print >> fout, ' calc vxcnn'
        if FORT:
            #  interstitial integral Integrate[ (iK - iG)\vr, {\vr in interstitial}] where K is averaged over all members of the star
            #  Here K is read from V_{xc} and G is from the fixed basis with length smaller than pw.npw2apw
            istpw = fvxcn.intstipw(ks.ksxc, pw.gindex, strc.timat, strc.tau, strc.vpos, pw.vmt, strc.rmt, strc.mult, kqm.k2cartes, pw.npw2apw)
        else:
            #  We will compute :
            #  iK_sym[isym,j,ik] = dot(ks.ksxc[ik,:], timat[isym][:,j]) = dot(timat[isym].T[j,:], ks.ksxc.T[:,ik]) 
            imat  = zeros((strc.Nsym,3,3),dtype=int)
            for isym in range(strc.Nsym):
                imat[isym,:,:] = strc.timat[isym,:,:].T
            imat = reshape(imat,(strc.Nsym*3,3))
            iK_symat = reshape( la.matmul(imat, ks.ksxc.T), (strc.Nsym,3,len(ks.ksxc)) ) # iK_symat[isym,:,ik]
            iK_tsymat = zeros( (strc.Nsym,len(ks.ksxc),3), dtype=int )                       # iK_tsymat[isym,ik,3]
            for isym in range(strc.Nsym):
                iK_tsymat[isym,:,:] = iK_symat[isym,:,:].T
            #
            # phasemat[ik,isym] = dot(ks.ksxc[ik,:],strc.tau[isym,:])*2*pi  = dot(ks.ksxc[ik,:],strc.tau.T[:,isym])*2*pi 
            phasemat = la.matmul(ks.ksxc, strc.tau.T)*(2*pi*1j)
            phimat = exp(-phasemat)       # phasemat[ik,isym], phimat[ik,isym]
            #
            istpw = zeros( (pw.npw2apw, len(ks.ksxc) ), dtype=complex, order='F')
            for ik in range(len(ks.ksxc)): # over all G vectors in Vxc
                for ig in range(pw.npw2apw): # over all G vectors we are interested in
                    istpw[ig,ik] = sum([fvxcn.int1ipw(iK_tsymat[isym,ik,:]-pw.gindex[ig,:], kqm.k2cartes, pw.vmt, strc.rmt, strc.vpos, strc.mult) * phimat[ik,isym] for isym in range(strc.Nsym)])/float(strc.Nsym)
                #to_print = istpw[ik,:10].real
                #print '%3d' % (ik,), '%12.7f'*10 % tuple(to_print)
            #for ik in range(len(ks.ksxc)):
            #    print '%3d' % (ik,), '%12.7f'*10 % tuple(istpw[:10,ik].real)

        #i_g0 = pw.convert_ig0_2_array()
        #nbmax = min([shape(ks.all_As[irk])[0] for irk in range(len(ks.all_As))])
        nbmax = ks.nbgw  # the band cutoff for this calculation
        Vxct = zeros((len(kqm.kirlist),nbmax,nbmax),dtype=complex)
        isp=0
        for irk in range(len(kqm.kirlist)):
            kl = array(kqm.kirlist[irk,:])/float(kqm.LCM)  # k in semi-cartesian form
            ik = kqm.k_ind[irk]   # index in all-kpoints, not irreducible
            tm10 = timer()
            # Get alm,blm,clm coefficients
            if DMFT1_like:
                isym = kqm.iksym[ik]
                timat_ik, tau_ik = strc.timat[isym].T, strc.tau[isym,:]
                alm,blm,clm = lapwc.dmft1_set_lapwcoef(False,1, True, kl, kl, timat_ik,tau_ik, ks.indgkir[irk], ks.nv[irk], pw.gindex, radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.tauij, latgen.Vol,  kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax)
            else:    
                alm,blm,clm = lapwc.gap2_set_lapwcoef(kl, ks.indgk[ik], 1, True, ks.nv[irk], pw.gindex, radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.Vol, kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax)
            tm11 = timer()
            
            (ngk,ntnt,ndf) = shape(alm)
            (ngk,nLOmax,ntnt,ndf) = shape(clm)
            # eigen-vector
            Aeig = array(ks.all_As[irk][:nbmax,:],dtype=complex)
            if False:
                print 'alm,blm'
                for idf in range(ndf):
                    for ig in range(ngk):
                        for ilm in range(ntnt):
                            print 'irk=%2d idf=%2d i=%3d j=%3d alm=%14.9f%14.9f blm=%14.9f%14.9f' % (irk+1, idf+1, ig+1, ilm+1, alm[ig,ilm,idf].real, alm[ig,ilm,idf].imag, blm[ig,ilm,idf].real, blm[ig,ilm,idf].imag)
                print 'eigenvector', shape(ks.all_As[irk])
                for ib in range(nbmax):
                    for j in range(ngk):
                        print 'irk=%2d ib=%3d j=%3d A=%14.9f%14.9f' % (irk+1,ib+1,j+1,Aeig[ib,j].real,Aeig[ib,j].imag)
            
            # And now change alm,blm,clm to band basis, which we call alfa,beta,gama
            alfa = reshape( la.matmul(Aeig, reshape(alm, (ngk,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
            beta = reshape( la.matmul(Aeig, reshape(blm, (ngk,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
            if in1.nlomax > 0:
                gama = reshape( la.matmul(Aeig, reshape(clm, (ngk,ntnt*ndf*nLOmax)) ), (nbmax,nLOmax,ntnt,ndf) )
            else:
                gama = zeros((nbmax,nLOmax,ntnt,ndf),dtype=complex,order='F')
            
            tm12 = timer()
            t_lapwcoef += tm11-tm10
            t_lapack   += tm12-tm11
            
            # The muffin-thin part of <psi_{k}|V_{xc}|psi_{k}>
            Vxcmt = fvxcn.mt_vxcnn(irk,ks.uxcu,alfa,beta,gama,pb.cgcoef,ks.lmxc_f,ks.lxcm_f,ks.ibgw,ks.nbgw,strc.mult,strc.iatnr,in1.nLO_at,in1.lmax,in1.nt)
            tm13 = timer()
            
            if self.DEBUG:
                print 'alfabeta'
                for idf in range(ndf):
                    for l in range(in1.lmax):
                        for m in range(-l,l+1):
                            lm = l*l + l + m 
                            for ie in range(nbmax):
                                print ('%s%2d '+'%s%3d '*2+ '%s%14.9f%14.9f '*3) % ('irk=',irk+1,'lm=',lm+1,'ie=',ie+1, 'alf=', alfa[ie,lm,idf].real, alfa[ie,lm,idf].imag, 'bet=', beta[ie,lm,idf].real, beta[ie,lm,idf].imag, 'gam=', gama[ie,0,lm,idf].real, gama[ie,0,lm,idf].imag)
            
            # This part calculates index array jimp[iv1,iv2] which finds index of two reciprocal vectors G1-G2
            # in fixed basis (pw.gindex), where  G1 and G2 are reciprocal vectors from Hamiltonian
            # corresponding to a particular irreducible k-point.
            nvk = ks.nv[irk]
            if FORT:
                jipw = fvxcn.two_gs_to_one(pw.gindex,ks.indgk[ik][:nvk],self.i_g0)
            else:
                jipw = -ones((nvk,nvk),dtype=int)  # jipw is index for difference between two G points from vector file
                for i,(iv1,iv2) in enumerate(itertools.product(range(nvk),range(nvk))):
                    iGc1 = pw.gindex[ks.indgk[ik][iv1]]
                    iGc2 = pw.gindex[ks.indgk[ik][iv2]]
                    diG = tuple(iGc1-iGc2)
                    if pw.ig0.has_key(diG):
                        jipw[iv1,iv2] = pw.ig0[diG]
                
            # Now we calculate interstitail integral for the difference G1-G2 (with G1 and G2 from Hamiltonian basis)
            # and K from V_{xc}:
            #       inti[1,2] = Integrate[ (iK - i(G1-G2))\vr, {\vr in interstitial}]  = <G1|e^{iK\vr}|G2>
            Aeigc = Aeig[:,:nvk] # this is the eigenvector withouth local orbitals.
            # interstitial exchange-correlation
            Vxci = zeros((nbmax,nbmax),dtype=complex)
            for ikxc in range(len(ks.ksxc)):
                if FORT:
                    inti = fvxcn.get_inti(jipw, istpw[:,ikxc])
                else:
                    inti = zeros((nvk,nvk),dtype=complex, order='F')
                    for i,(iv1,iv2) in enumerate(itertools.product(range(nvk),range(nvk))):
                        if jipw[iv1,iv2] >= 0:
                            inti[iv1,iv2] = istpw[jipw[iv1,iv2],ikxc]
                
                # Vxc_{ij} = \sum_{G1,G2,K} (Aeig_{i,G1})* V_{xc,K}<G1| e^{iK\vr}|G2> Aeig_{j,G2}
                tmat1 = la.matmul( conj(Aeigc), inti ) * ks.Vxcs[ikxc]
                Vxci += la.matmul( tmat1, Aeigc.T)
                
            tm14 = timer()
            t_muffint  += tm13-tm12
            t_interst += tm14-tm13
            
            ns,ne = ks.ibgw,ks.nbgw
            Vxct[irk,ns:ne,ns:ne] = Vxcmt[ns:ne,ns:ne] + Vxci[ns:ne,ns:ne]
            if self.DEBUG:
                for i in range(len(Vxcmt)):
                    print '%s%3d %s%14.9f%14.9f' % ('ie=', i+1, 'vxcmt=', Vxcmt[i,i].real, Vxcmt[i,i].imag)
                for ie in range(nbmax):
                    print '%s%3d %s%3d %s%14.9f %14.9f' % ('irk=', irk, 'ie=', ie, 'Vxci=', Vxci[ie,ie].real, Vxci[ie,ie].imag)

        print >> fout, '## Vxc:      t(lapwcoef)           =%14.9f' % t_lapwcoef
        print >> fout, '## Vxc:      t(lapack)             =%14.9f' % t_lapack
        print >> fout, '## Vxc:      t(muffin)             =%14.9f' % t_muffint
        print >> fout, '## Vxc:      t(interst)            =%14.9f' % t_interst
        
        
        for irk in range(len(kqm.kirlist)):
            for i in range(nbmax):
                print >> fout, '%s%3d %s%3d %s%14.9f %14.9f' % ('irk=', irk+1, 'ie=', i+1, 'Vxc=', Vxct[irk,i,i].real, Vxct[irk,i,i].imag)
        return Vxct

    def optimal_eta(self, rbas, br2):
        """ This function calculates the optimal value of eta for the lattice
            summations needed to obtain the structure constants.
        """
        lrbs = [ linalg.norm(rbas[i,:]) for i in range(3)]
        lgbs = [ linalg.norm(br2 [:,i]) for i in range(3)]
        return sqrt(2*min(lrbs)/min(lgbs))

    def calc_tildeg(self, lmax):
        """ Calculates $\tilde{g}_{lm,l'm'}$ according to equation \ref{tildea},
             \begin{equation}\label{calctilg}
               \tilde{g}_{lm,l'm'}=\sqrt{4\pi}(-1)^{l}\sqrt{\frac{(l+l'+m+m')!(l+l'-m-m')!}%
                {(2l+1)(2l'+1)[2(l+l')+1]!(l+m)!(l-m)!(l'+m')!(l'-m')!}}
              \end{equation}
             needed for the calculation of the structure constants, for l and l' = 0
        ...  lmax and stores them in memory.
        """

        self.tilg = zeros((lmax+1)*(lmax+2)*(lmax+3)*(3*lmax+2)/12)
        i=0
        for l1 in range(lmax+1):
            p0 = (-1)**l1 * 8*pi*sqrt(pi)
            for l2 in range(l1+1):
                denom = (2*l1+1)*(2*l2+1)*(2*(l1+l2)+1)
                p1 = p0 / sqrt(denom)                
                for m1 in range(-l1,l1+1):
                    for m2 in range(l2+1):
                        jm1 = l1+m1
                        jm2 = l2+m2
                        jm3 = l1-m1
                        jm4 = l2-m2
                        combj = special.binom( jm1+jm2, jm1 ) * special.binom( jm3+jm4, jm3 )
                        self.tilg[i] = p1*sqrt(combj)
                        i+=1

        #for i in range(len(self.tilg)):
        #    print '%s%5d %s%20.10f' % ('i=', i, 'tildg=', self.tilg[i])
    
    def r_cutoff(self, tol, eta, lambdamax, alat):
        """ Estimates the cutoff radius of the sums in real space for the calculation of the structure constants by the solving the equation:
       \begin{equation}
       \mathfrak{E}_{R,\lambda}^{\textrm{tol}}=\left\{%
       \begin{array}{ll}
       \frac{4\pi}{(\lambda -2)\Gamma(\lambda+\tfrac{1}{2})}%
       \left(\frac{\Gamma[\tfrac{\lambda}{2}+\tfrac{3}{2},\left(\tfrac{R_c}{\eta}\right)^2]}%
       {\eta^{\lambda-2}}-\frac{\Gamma[\lambda+\tfrac{1}{2},\left(\tfrac{R_c}{\eta}\right)^2]}%
       {R_c^{\lambda-2}}\right)&\lambda \neq 2\\
       \frac{4\pi}{\Gamma(\tfrac{5}{2})}\left[\tfrac{\eta}{R_c}%
       \Gamma[3,\left(\tfrac{R_c}{\eta}\right)^2]-\Gamma[\tfrac{5}{2},\left(\tfrac{R_c}{\eta}\right)^2]\right]&
       \lambda=2\\
       \end{array}
       \right.
       \end{equation}
        and taking the maximum value of $R_c$ obtained for $\lambda = 1...$ \verb lambdamax.
        """
        rnot = max(alat)
        rct = 50 * ones(lambdamax+1)
        eps = zeros(lambdamax+1)
        which_l = range(lambdamax+1)
        which_l.remove(2)
        four_pi = 4*pi
        ls = arange(lambdamax+1)           # all possible l's
        ls[2] = 0 # jus so that we do not divide by 0
        etal = 4*pi/(eta**(ls-2) * (ls-2)) # 4*pi/(eta**(l-2)*(l-2))
        ls[2] = 2
        gmms = special.gamma(ls+0.5)       # 
        gmns = special.gamma((ls+1)/2.)    #
        for i in range(1,101):
            x = i/2.0
            x2 = x**2
            gaml32 = special.gammaincc( 3,   x2 ) * 2
            gmm = gmms[2] # special.gamma( 2.5 )
            gaml12 = special.gammaincc( 2.5, x2 ) * gmm
            prefac = four_pi/gmm
            eps[2]= abs( prefac * ( gaml32 / x - gaml12 ) )
            if (eps[2] < tol) and (x < rct[2]) :
                rct[2] = x
            #print '%s%3d %s%10.4f %s%3d %s%12.4f %s%12.4f %s%12.4f %s%12.4f %s%12.9f %s%12.9f' % ('i=',i,'x=',x,'l1=',2,'g32=',gaml32,'g12=',gaml12,'gm=',gmms[2],'pref=',prefac,'eps=',eps[2],'rct=',rct[2])
            for l1 in which_l:
                gaml32 = special.gammaincc( (l1+1)/2., x2 ) * gmns[l1]
                gaml12 = special.gammaincc(  l1 + 0.5, x2 ) * gmms[l1]
                eps[l1] = abs( etal[l1] * ( gaml32 - gaml12/x**(l1-2) ) / gmms[l1] )
                if (eps[l1] < tol) and (x < rct[l1]) :
                    rct[l1] = x
                
                #print '%s%3d %s%10.4f %s%3d %s%12.4f %s%12.4f %s%12.4f %s%12.4f %s%12.9f %s%12.9f' % ('i=',i,'x=',x,'l1=',l1,'g32=',gaml32,'g12=',gaml12,'gm=',gmms[l1],'pref=',prefac,'eps=',eps[l1],'rct=',rct[l1])
        return max(rct) * eta
    
    def g_cutoff(self, tol, eta, lambdamax, pia, Vol):
        """ Estimates the cutoff radius of the sums in reciprocal space for the
         calculation of the structure constants by the solving the equation:
         \begin{equation}
         \mathfrak{E}_{G,\lambda}^{\textrm{tol}}=\frac{8(\pi)^{\frac{5}{2}}}{\Omega%
         \Gamma(\lambda+\tfrac{1}{2})\eta^{\lambda+1}}%
         \Gamma\left[\tfrac{\lambda+1}{2},\left(\tfrac{\eta G_c}{2}\right)^2\right]
         \end{equation}
          and taking the maximum value of $G_c$ obtained for $\lambda = 1...$ \verb lambdamax.
        """
        rnot = max(pia)
        rct = 50 * ones(lambdamax+1)
        eps = zeros(lambdamax+1)
        ls = arange(lambdamax+1)
        gmms = special.gamma(ls+0.5)       # 
        prefac = 8 * pi**2.5 /( Vol * eta**(ls+1) * gmms )
        gmns = special.gamma((ls+1)/2.)    #
        for i in range(1,101):
            x = i/2.0
            for l1 in range(lambdamax+1):
                gaml32 = special.gammaincc( (l1+1)/2., x**2 ) * gmns[l1]
                eps[l1] = abs( prefac[l1]*gaml32 )
                if (eps[l1] < tol) and (x < rct[l1]) :
                    rct[l1] = x
        return max(rct)*2/eta

    def frcut_coul(self, rcut_coul, alat, ndivs, Vol):
       if rcut_coul < 0.0:    # truncated/screened Coulomb interaction 
           if int(rcut_coul) == -2:
               return max(alat[:]*ndivs[:])/2
           else:
               vtot = Vol*ndivs[0]*ndivs[1]*ndivs[2]
               return (vtot/(4*pi/3))**(1/3.)
           print >> fout, 'set default rcut_coul to ', rcut_coul
       else:
           return rcut_coul

    def Zero_Momentum_Coulomb(self, kmax, mbsize, ortho, strc, latgen, pb, fout):
        gmax = 10 * pb.kmr * kmax
        ng = [ int(gmax*latgen.pia[i])+1 for i in range(3)]
        np_tmp, ngindx = fkp.count_number_pw(ng,kmax,gmax,latgen.br2)
        glen, G_c, gindex = fkp.generate_pw(ng,ngindx,gmax,latgen.pia,latgen.br2,ortho,strc.lattice[1:3]=='CXZ',True)
        
        ngs, gind4 = fCl.pw_get_short_length(glen)
        print >> fout, 'gmax=', gmax, 'ng=', ng, 'ngindx=', ngindx, 'len(gindex)=', len(gindex), 'ngs=', ngs
        
        glen0,phase = fCl.pw_get_short_and_phase(ngs,gind4,gindex,glen,strc.vpos,strc.mult)
        
        sing = {}
        irms = [[] for iat in range(strc.nat)]
        for iat in range(strc.nat):
            #npt = strc.nrpt[iat]
            #dh  = log(strc.rmt[iat]/strc.r0[iat])/(npt - 1)      # logarithmic step for the radial mesh
            #dd = exp(dh)
            #rx = strc.r0[iat]*dd**range(npt)
            rx, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
            
            sinf = zeros((ngs,npt))
            for ig in range(ngs):   # now iterating only over G's of different length
                # each |G| length appears only once
                sinf[ig,:] = sin( glen0[ig] * rx )  # sinf[ir] = sin(G*r_i)
                
            for irm in range(len(pb.big_l[iat])):  # over all product functions
                L = pb.big_l[iat][irm]             # L of the product function
                if L==0:                           # only L==0 is important for q==0
                    au, bu = pb.ul_product[iat][irm], pb.us_product[iat][irm]
                    # now iterating only over G's of different length                            
                    sing[(iat,irm)] = [ rd.rint13g(strc.rel, au, bu, sinf[ig], sinf[ig], dh, npt, strc.r0[iat]) for ig in range(ngs) ]
            
                    irms[iat].append(irm) # which product basis have L==0

        const = 16*pi**2/latgen.Vol
        G_2_4 = const/glen0**4  # const/G^4
        sinsing={}
        for (iat,irm) in sing.keys():
            for (jat,jrm) in sing.keys():
                sinsing[(iat,jat,irm,jrm)] = G_2_4[:] * sing[(iat,irm)][:] * sing[(jat,jrm)][:]  # <u_i|sin(G*r)><u_j|sin(G*r)>*const/G^4

        
        Vmat = zeros((mbsize,mbsize), dtype=complex)
        for idf in range(len(pb.atm)):
            iat = pb.atm[idf]
            for irm in irms[iat]:
                im = pb.iProd_basis[(idf,irm,0,0)]
                for jdf in range(len(pb.atm)):
                    jat = pb.atm[jdf]
                    for jrm in irms[jat]:
                        jm = pb.iProd_basis[(jdf,jrm,0,0)]
                        vmat = sum(phase[idf,jdf,:] * sinsing[(iat,jat,irm,jrm)][:])
                        Vmat[im,jm] = vmat
                        
        return Vmat

    def Calc_Olap_PW_Product(self, iat, gqlen, Lmax, rx, dh, strc, pb):
        # calcjlam
        # First we compute spherical bessel functions for all plane waves up to large cutoff ( indgq[iq,:] ). Different lengths are storeed in gqlen, which we are using here.
        npt = len(rx)
        Kqr = outer(gqlen, rx)                           # Kq[ig]*rx[ir]
        jln = reshape( fCl.spher_bessel(Lmax, Kqr.flatten()), (Lmax+1,len(gqlen),npt) )    # all spherical bessels: jln[L,iq,ir] = j_L( |G+q|*r )
        for ir in range(npt):
            jln[:,:,ir] *= rx[ir]     # jln[L,iq,ir] = j_L( |K+g|*r ) * r
        jlam = zeros( ( len(pb.big_l[iat]), len(gqlen) ) )   # <j_L((q+K)r)|u_{mix,L}>=Int[j_L((q+G)r) u_{mix,L}/r r^2,r]
        for irm in range(len(pb.big_l[iat])):  # over all product functions
            L = pb.big_l[iat][irm]             # L of the product function
            au, bu = pb.ul_product[iat][irm], pb.us_product[iat][irm]
            for ig in range(len(gqlen)):
                jlam[irm,ig] = rd.rint13g(strc.rel, au, bu, jln[L,ig], jln[L,ig], dh, npt, strc.r0[iat])  # <j_L((q+G)r)|u_{mix,L}>=Int[j_L((q+G)r) u_{mix,L}/r r^2,r]
        return jlam

    def OrthogonalizedBasisInterstitials(self, iq, pw):
        ### mpwipw
        ### diagipw
        if False:
            print 'pw.indgq'
            for ipw in range(pw.ngq[iq]):
                idG = pw.gindex[pw.indgq[iq,ipw],:]
                print '%4d%4d' % (ipw+1,pw.indgq[iq,ipw]+1), ' ', ('%3d'*3) % tuple(idG), '  %14.10f%14.10f' % (pw.ipwint[ipw].real,pw.ipwint[ipw].imag)
            print 'pw.indgq_done'
        
        if FORT:
            # Computes overlap between plane waves : olap[jpw,ipw] = <G_{jpw}|G_{ipw}>_{interstitials}/V_{cell}
            olap = fCl.cmp_pw_olap2(pw.indgq[iq], pw.ipwint, pw.gindex, self.i_g0, pw.ngq[iq], pw.ngq[iq])
        else:
            # Computes overlap between plane waves : olap[jpw,ipw] = <G_{jpw}|G_{ipw}>_{interstitials}/V_{cell}
            olap = zeros( (pw.ngq[iq],pw.ngq[iq]), dtype=complex)
            for ipw in range(pw.ngq[iq]):
               olap[ipw,ipw] = pw.ipwint[0]
               for jpw in range(ipw+1,pw.ngq[iq]):
                   idG = pw.gindex[pw.indgq[iq,ipw],:] - pw.gindex[pw.indgq[iq,jpw],:] # G_{ipw}-G_{jpw}
                   idg = pw.ig0[tuple(idG)]                                            # index of G_{ipw}-G_{jpw}
                   olap[jpw,ipw] =      pw.ipwint[idg]                                 # olap(jpw,ipw) = <G_{jpw}|G_{ipw}>=Integrate[e^{(G_{ipw}-G_{jpw})r},{r in interstitials}]
                   olap[ipw,jpw] = conj(pw.ipwint[idg])                                # olap(ipw,jpw) = <G_{ipw}|G_{jpw}>
                   #print '%4d%4d' % (pw.indgq[iq,ipw]+1,pw.indgq[iq,jpw]+1), ('%3d'*3) % tuple(idG), '%4d' % (idg+1,)

        if sum(abs(olap.imag)) < 1e-7:  # overlap is real
            olap = array(olap.real,dtype=float)
        
        if False:
            print 'olap_PW'
            for ipw in range(pw.ngq[iq]):
                for jpw in range(pw.ngq[iq]):
                    print '%3d %3d %17.12f%17.12f%17.12f%17.12f' % (ipw+1,jpw+1,olap[ipw,jpw].real,-olap[ipw,jpw].imag,olap[jpw,ipw].real,-olap[jpw,ipw].imag)
            
        # diagonalizing overlap basis for interstitial
        epsipw, sgi = linalg.eigh(olap)   # sgi[pw.ngq[iq],pw.ngq[iq]]
        
        if False:
            print 'epsipw:'
            for i in range(len(epsipw)):
                print '%2d %s%16.12f ' % (i,'e=',epsipw[i])
                #print ('%10.6f '*len(epsipw)) % tuple(sgi[:,i].real)
        

        if False: ##### ATTENTION : ONLY FOR DEBUGGING
            dat = loadtxt('si_eigvec_fort.dat')
            sgi_fort = dat[:,0::2] + dat[:,1::2]*1j
            sgi = sgi_fort
            
        # Finally, taking 1/sqrt(olap) 
        if Real_sqrt_olap:
            # Now we will try to compute real 1/sqrt(O) = U 1/sqrt(eig) U^H
            sgiH = la.matmul(sgi, dot(diag(1/sqrt(abs(epsipw))), conj(sgi.T)) )
        else:
            #sgi = ImproveDegEigenvectors(epsipw, sgi)
            sgiH = conj(sgi.T)
            for i in range(len(epsipw)):
                sgiH[i,:] *= 1/sqrt(abs(epsipw[i]))

        if False:
            print 'sgi='
            for ipw in range(pw.ngq[iq]):
                for jpw in range(pw.ngq[iq]):
                    print '%3d %3d %17.12f%17.12f%17.12f%17.12f' % (ipw+1,jpw+1,sgiH[ipw,jpw].real,-sgiH[ipw,jpw].imag,sgiH[jpw,ipw].real,-sgiH[jpw,ipw].imag)
        
        if FORT:
            # coul_mpwipw
            # note that this matrix of overlap tmat[i,j]=<G_i|G_j> is non-diagonal, and is of shape tmat[ngq,ngq_barc], where ngq_barc>ngq.
            # Here ngq is the size of the basis in the interstitials, while ngq_barc is some high-G cutoff.
            tmat = fCl.cmp_pw_olap2(pw.indgq[iq], pw.ipwint, pw.gindex, self.i_g0, pw.ngq[iq], pw.ngq_barc[iq]) 
        else:
            tmat = zeros((pw.ngq[iq],pw.ngq_barc[iq]),dtype=complex)
            for ipw in range(pw.ngq_barc[iq]):
              for jpw in range(pw.ngq[iq]):
                idG = pw.gindex[pw.indgq[iq,ipw],:] - pw.gindex[pw.indgq[iq,jpw],:]  # G_{ipw}-G_{jpw}
                idg = pw.ig0[tuple(idG)]                                             # index of G_{ipw}-G_{jpw}
                tmat[jpw,ipw] = pw.ipwint[idg]                                       # tmat[jpw,ipw] = <G_{jpw}|G_{ipw}>_{int}=Integrate[e^{i*(G_{ipw}-G_{jpw}r)},{r in interstitial}]
                #print ipw, jpw, idG, idg, tmat[ipw,jpw]
        # mpwipw[i,j] = <G_i,G_j>
        # where G_j has much larger dimension that G_i, because it needs to allow the possibility of G_j to be any combination of reciprocal vectors from KS-vector file.
        mpwipw = dot(sgiH,tmat)   # mpwipw[pw.ngq[iq],pw.ngq_barc[iq]]

        if False:
            print 'mpwipw'
            for jpw in range(pw.ngq[iq]):
                for ipw in range(pw.ngq_barc[iq]):
                    print '%3d %4d %16.10f %16.10f' % (jpw+1,ipw+1,mpwipw[jpw,ipw].real,mpwipw[jpw,ipw].imag)
            print 'end_mpwipw'
        return mpwipw
    
    def fwi0(self, latgen, pb):
        loctmatsize = len(pb.Prod_basis)
        wi0 = zeros(loctmatsize)
        fct = sqrt(4*pi/latgen.Vol) 
        for idf in range(len(pb.atm)):
            iat = pb.atm[idf]
            for irm in range(len(pb.big_l[iat])):  # over all product functions
                L = pb.big_l[iat][irm]             # L of the product function
                if L==0:                           # only L==0 is important for q==0
                    im = pb.iProd_basis[(idf,irm,0,0)]
                    wi0[im] = fct * pb.rtl[iat,irm]  # 4*pi/Vol*< r^0 | u_{im,at} > because rtl=< r^L | u_{im,at} >
                    #print 'idf=', idf, 'irm=', irm, 'L=', L, 'im=', im, 'wi0=', wi0[im]
        return wi0
                
    def Coulomb_from_PW(self, iq, io, strc, in1, latgen, kqm, ks, radf, pw, pb, fout):
        """ Calculating Coulomb by  
                <u_{product} | K+q > 4*pi/(K+q)^2 < K+q| u_{product}> 
            where  u_{product} is the product basis defined in enire space, i.e.,
            both in MT and interstitials. 
            In the MT the plane wave is expaned 
                  <r|K+q> = 4*pi*(i)^l j_l(|q+K|*r) Y_{lm}(\vr) Y_{lm}(q+K)
            and the bessel functions are integrated in MT. 

            Also calculates mpwipw[G,K] == 1/sqrt(O) * <e^{iG*r}| e^{i*K*r}>
            Maybe we should compute mpwiw outside this, because it is used to define 
            othogonalized basis in the interstitials
        """
        ortho = (latgen.ortho or strc.lattice[1:3]=='CXZ')
        alat = array([strc.a, strc.b, strc.c])

        vq = kqm.qlistc[iq,:]/float(kqm.LCMq)
        
        if False:  ##### ATTENTION ONLY FOR DEBUGGING?????
            dat = array(loadtxt('indgq_fort.dat').T, dtype=int)
            pw.indgq[iq,:pw.ngq_barc[iq]] = dat[1,:]
            #for ig in range(pw.ngq[iq]):
            #    print ig, pw.indgq[iq,ig], dat[1,ig]
        
        gqlen = array([ pw.gqlen[iq, pw.G_unique[iq,i]] for i in range(pw.ngqlen[iq]) ])  # lengths of unique |q+G|
        Lmax =  max( [ max(pb.big_l[iat]) for iat in range(strc.nat)])                    # maximum possible value of L in product basis
        
        loctmatsize = len(pb.Prod_basis)
        mbsize =  loctmatsize + pw.ngq[iq]
        im_start = zeros(strc.nat+1, dtype=int)                                           # first index for product basis on each atom
        idf = 0
        for iat in range(strc.nat):
            im_start[iat] = pb.iProd_basis[(idf,0,0,0)]                                   # first index of product basis on each atom
            idf += strc.mult[iat]
        im_start[strc.nat] = loctmatsize                                                  # the last index for product basis on the last atom
        # mpwmix : <u_{product_basis_everywhere}| e^{i*K*r}>
        mpwmix = zeros((mbsize,pw.ngq_barc[iq]), dtype=complex)
        for iat in range(strc.nat):
            rx, dh, npt = strc.radial_mesh(iat)
            jlam = self.Calc_Olap_PW_Product(iat, gqlen, Lmax, rx, dh, strc, pb)
            istr, iend = im_start[iat], im_start[iat+1]
            # Now computing matrix elements <e^{(q+G)\vr}|V_{coul}|u_{irm}>  = 4*pi/|q+G|^2 e^{-i(q+G)_R }<q+G|u_{irm}>
            mpwmix[istr:iend,:] = fCl.mixed_coulomb(vq,iat+1,False,jlam,pw.gindex,pw.indgq[iq],pw.gqlen[iq],gqlen,pw.G_unique[iq],pw.ngq_barc[iq],iend-istr,pb.big_l[iat],kqm.k2cartes,strc.rotloc,latgen.trotij,strc.mult,strc.vpos,latgen.Vol)
        mpwmix = conj(mpwmix)  # <u_{product_basis}| K>
        self.mpwipw = self.OrthogonalizedBasisInterstitials(iq, pw)  # 1/sqrt(Olap)*<G|K>
        mpwmix[loctmatsize:mbsize,:] = self.mpwipw
        
        #if (vq==0): mpwmix[:,0] += wi0[:]. But this is not used anyway. So should be irrelevant.
        
        gqlen = pw.gqlen[iq,:pw.ngq_barc[iq]]
        if sum(abs(vq)) < 1e-7: gqlen[0] = 1e100
        Q2 = gqlen**2
        Vq = 4*pi/Q2

        Vmm = conj(mpwmix)
        Vmm *= Vq          # Vmm[i,j] = conj(mpwmix[i,j]) * Vq[j]
        Vmat = dot(mpwmix, Vmm.T)


        self.ev, self.Vmat = linalg.eigh(Vmat)
        
        #self.wi0 = zeros(mbsize)
        #self.wi0[:loctmatsize] = me.wi0(latgen, pb)
        #self.wi0[loctmatsize:mbsize] = mpwipw[:pw.ngq[iq]]
        #print 'wi0'
        #for i in range(len(wi0)):
        #    print >> fout, '%3d %16.9f%16.9f' % (i, wi0[i].real, wi0[i].imag)
        print >> fout, 'mpwmix loctmatsize=', loctmatsize
        for j in range(pw.ngq_barc[iq]):
            for i in range(mbsize):
                print >> fout, '%3d %3d %16.9f%16.9f' % (j+1, i+1, mpwmix[i,j].real, mpwmix[i,j].imag)
        #print >> fout, '** barc eigenvalues **'
        #for i in range(len(ev)):
        #    print >> fout, '%5d%18.9f' % (i+1, ev[-1-i])
        
    def Coulomb(self, iq, io, strc, in1, latgen, kqm, ks, radf, pw, pb, fout):
        """ Calculating Coulomb in the interstitail and mixed term, we use
                <u_{product} | K+q > 4*pi/(K+q)^2 < K+q| u_{product}> 
            where  u_{product} is the product basis defined in enire space.
            In the MT part we use the real space expression for the Coulomb,
            which is cutoff-free, but it requires two center Laplace expansion.
            Also needs Ewalds summation and real space integration.

            Also calculates mpwipw[G,K] == 1/sqrt(O) * <e^{iG*r}| e^{i*K*r}>
            Maybe we should compute mpwiw outside this, because it is used to define 
            othogonalized basis in the interstitials
        """
        t_sph_bess, t_mixed, t_MT = 0,0,0
        
        ortho = (latgen.ortho or strc.lattice[1:3]=='CXZ')
        alat = array([strc.a, strc.b, strc.c])
        
        tm9 = timer()
        #rcut_coul = self.frcut_coul(io.rcut_coul, alat, kqm.ndiv, latgen.Vol)
        
        vq = kqm.qlistc[iq,:]/float(kqm.LCMq)
        
        ndf = sum(strc.mult)
        lam_max = 4*(io.lmbmax+1)
        loctmatsize = len(pb.Prod_basis)
        mbsize =  loctmatsize + pw.ngq[iq]
        kmax = in1.rkmax/min(strc.rmt)
        print >> fout, 'mbsize=', mbsize, 'loctmatsize=', loctmatsize, 'kmax=', kmax
        
        if iq==0:
            Vmat = self.Zero_Momentum_Coulomb(kmax, mbsize, ortho, strc, latgen, pb, fout)
        else:
            Vmat = zeros((mbsize,mbsize), dtype=complex)
        
        if False:  ##### ATTENTION ONLY FOR DEBUGGING?????
            dat = array(loadtxt('indgq_fort.dat').T, dtype=int)
            pw.indgq[iq,:pw.ngq_barc[iq]] = dat[1,:]
            #for ig in range(pw.ngq[iq]):
            #    print ig, pw.indgq[iq,ig], dat[1,ig]
        
        tm10 = timer()
        self.mpwipw = self.OrthogonalizedBasisInterstitials(iq, pw)
        
        tm16 = timer()
        cq = dot(kqm.k2cartes, vq)
        sgm = fCl.ewald_summation(cq, lam_max, strc.vpos, latgen.br2, latgen.rbas, alat, latgen.Vol, ortho, self.r_cf, self.g_cf, self.eta) # coul_strcnst
        tm17 = timer()
        
        gqlen = array([ pw.gqlen[iq, pw.G_unique[iq,i]] for i in range(pw.ngqlen[iq]) ])  # lengths of unique |q+G|
        Lmax =  max( [ max(pb.big_l[iat]) for iat in range(strc.nat)])                    # maximum possible value of L in product basis
        nmix = array([ len(pb.big_l[iat]) for iat in range(strc.nat) ], dtype=int)        # nmix is the length of product basis on each atom
        max_nmix = max( nmix )
        im_start = zeros(strc.nat+1, dtype=int)                                           # first index for product basis on each atom
        big_l = zeros( (max_nmix,strc.nat), dtype=int, order='F' )                        # big_l in fortran array form
        idf = 0
        for iat in range(strc.nat):
            im_start[iat] = pb.iProd_basis[(idf,0,0,0)]                                   # first index of product basis on each atom
            nm = len(pb.big_l[iat])                                                       # length of the atomic product basis on this atom
            big_l[:nm,iat] = pb.big_l[iat][:]                                             # saving big_l to firtran-like array
            idf += strc.mult[iat]                                                         # what is the index of atom in array of all atoms, including equivalent
        im_start[strc.nat] = loctmatsize                                                  # the last index for product basis on the last atom
        
        #tm18 = timer()
        rtlij = zeros((max_nmix,max_nmix,strc.nat,strc.nat),order='F')
        for iat in range(strc.nat):
            nim = len(pb.big_l[iat])
            for jat in range(strc.nat):
                njm = len(pb.big_l[jat])
                rtlij[:njm,:nim,jat,iat] = outer(pb.rtl[jat,:njm] , pb.rtl[iat,:nim])   # bug jul.7 2020
        
        #tm19 = timer()
        Vmatit = zeros((loctmatsize,pw.ngq_barc[iq]), dtype=complex)
        
        for iat in range(strc.nat):
            rx, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
            
            tm20 = timer()
            ## First we compute spherical bessel functions for all plane waves up to large cutoff ( indgq[iq,:] ).
            # Different lengths are storeed in gqlen, which we are using here.
            # Next we compute matrix elements with the atomic centered product basis functions
            jlam = self.Calc_Olap_PW_Product(iat, gqlen, Lmax, rx, dh, strc, pb)
            tm21 = timer()
            
            istr, iend = im_start[iat], im_start[iat+1]
            
            # Now computing matrix elements <e^{(q+G)\vr}|V_{coul}|u_{irm}>  = 4*pi/|q+G|^2 e^{-i(q+G)_R }<q+G|u_{irm}>
            Vmatit[istr:iend,:] = fCl.mixed_coulomb(vq,iat+1,True,jlam,pw.gindex,pw.indgq[iq],pw.gqlen[iq],gqlen,pw.G_unique[iq],pw.ngq_barc[iq],iend-istr,pb.big_l[iat],kqm.k2cartes,strc.rotloc,latgen.trotij,strc.mult,strc.vpos,latgen.Vol)
            tm22 = timer()
            # Here is the muffin-thin part
            Vmat[:loctmatsize,:loctmatsize] += fCl.mt_coulomb(vq,iat+1,loctmatsize,big_l,nmix,im_start,self.tilg,pb.rrint,pb.djmm,rtlij,sgm,strc.mult,lam_max)
            tm23 = timer()
            
            t_sph_bess += tm21-tm20
            t_mixed += tm22-tm21
            t_MT += tm23-tm22

        #tm24 = timer()
        # This transforms into the plane-wave-product basis
        mat2 = dot(self.mpwipw, Vmatit.T)      # mat2[ pw.ngq[iq], loctmatsize]
        # Setting mixed- contribution to final storage
        Vmat[loctmatsize:mbsize,:loctmatsize] = mat2[:,:]
        Vmat[:loctmatsize,loctmatsize:mbsize] = conj(mat2.T)
        #tm25 = timer()
        
        # Now diagonal plane-wave part, i.e, plane-wave-product basis: mpwipw[i,j] * Vq[j] * mpwipw.H[j,k]
        gqlen = pw.gqlen[iq,:pw.ngq_barc[iq]]
        if sum(abs(vq)) < 1e-7: gqlen[0] = 1e100
        
        Q2 = gqlen**2
        Vq = 4*pi/Q2
        # Vmm[i,k] = mpwipw[i,j] * Vq[j] * mpwipw.H[j,k]
        Vmm = zeros((pw.ngq[iq],pw.ngq_barc[iq]), dtype=complex)
        for k in range(pw.ngq[iq]):
            Vmm[k,:] = Vq[:]*conj(self.mpwipw[k,:])
        m2 = dot(self.mpwipw, Vmm.T)
        Vmat[loctmatsize:mbsize,loctmatsize:mbsize] = m2[:,:]
        self.loctmatsize = loctmatsize
            
        tm11 = timer()
        print >> fout, '## Coulomb: t(optimal_etas_cutoffs)=%14.9f' % (tm10-tm9,)
        print >> fout, '## Coulomb: t(Ewald)               =%14.9f' % (tm17-tm16,)
        print >> fout, '## Coulomb: t(spherical_bessel)    =%14.9f' % t_sph_bess
        print >> fout, '## Coulomb: t(Muffin-thin)         =%14.9f' % t_MT
        print >> fout, '## Coulomb: t(mixed_MT_PW)         =%14.9f' % t_mixed

        self.ev, self.Vmat = linalg.eigh(Vmat)

        if False: ### ATTENTION  Only for debugging
            dat = loadtxt('CoulEigvec.dat')
            Vm = dat[:,0::2] + dat[:,1::2]*1j
            self.Vmat[:,:] = Vm[:,:]  # This should make Coulomb diagonalized matrix exactly the same as in fortran
        
        if iq==0:
            wi0 = zeros(mbsize, dtype=complex)
            wi0[:loctmatsize] = self.fwi0(latgen, pb)             # <Chi_product| 1 >
            wi0[loctmatsize:mbsize] = self.mpwipw[:pw.ngq[iq],0]  # <G_orthohonal_basis|G=0>
            wi0new = dot(wi0, conj(self.Vmat))                    # V^\dagger < Chi_product| 1> = overlap of singular vector with unity in space
            
            self.immax = argmax(abs(wi0new))                      # which singular eigenvector V_{l,:} has largest overlap with unity?
            print >> fout, '- Maximum singular eigenvector **', self.immax, abs(wi0new[self.immax]), self.ev[self.immax]
            
            alpha = (latgen.Vol/(6*pi**2))**(1./3.)
            ankp = kqm.ndiv[0]*kqm.ndiv[1]*kqm.ndiv[2]
            ifst = 1
            expQ = exp(-alpha*Q2[ifst:])
            sf1 = sum(expQ/gqlen[ifst:])/ankp
            sf2 = sum(expQ/Q2[ifst:])/ankp
            for iq in range(1,len(kqm.qlist)):
                Q1 = pw.gqlen[iq,:pw.ngq_barc[iq]]
                Q2 = Q1**2
                expQ = exp(-alpha*Q2)
                f1 = sum(expQ/Q1)/ankp
                f2 = sum(expQ/Q2)/ankp
                sf1 += f1
                sf2 += f2
            intf  = latgen.Vol/(4*pi**2)
            intf1 = intf/alpha
            intf2 = intf*sqrt(pi/alpha)
            self.singc1 = intf1-sf1   # see Eq.B6 (page 364) in gap2 paper (Eq.B3..Eq.B6)
            self.singc2 = intf2-sf2   # also Ref.61 in gap2 paper
            
        if False:
            print >> fout, 'Vmat='
            for im in range(mbsize):
                for jm in range(mbsize):
                    if (Vmat[im,jm] != 0 ):
                        if (im < loctmatsize and jm < loctmatsize):
                            (idf,irm,iL,iM) = pb.Prod_basis[im]
                            (jdf,jrm,jL,jM) = pb.Prod_basis[jm]
                            print >> fout, ('%3d '*2+'%2d'*3+'%3d'+','+'%2d'*3+'%3d'+' %14.9f%14.9f') % (im+1,jm+1,idf+1,iL,iM,irm+1, jdf+1,jL,jM,jrm+1, Vmat[im,jm].real, Vmat[im,jm].imag)
                        elif (im < loctmatsize):
                            (idf,irm,iL,iM) = pb.Prod_basis[im]
                            ii = pw.indgq[iq,jm-loctmatsize]
                            jG = pw.gindex[ii,:]
                            print >> fout, ('%3d '*2+'%2d'*3+'%3d'+','+'   '+'%2d'*3+' %14.9f%14.9f') % (im+1,jm+1,idf+1,iL,iM,irm+1, jG[0],jG[1],jG[2], Vmat[im,jm].real, Vmat[im,jm].imag)
                        elif (jm < loctmatsize):
                            (jdf,jrm,jL,jM) = pb.Prod_basis[jm]
                            ii = pw.indgq[iq,im-loctmatsize]
                            iG = pw.gindex[ii,:]
                            print >> fout, ('%3d '*2+'%2d'*3+'   '+','+'%2d'*3+'%3d'+' %14.9f%14.9f') % (im+1,jm+1,iG[0],iG[1],iG[2], jdf+1,jL,jM,jrm+1, Vmat[im,jm].real, Vmat[im,jm].imag)
                        else:
                            ii = pw.indgq[iq,im-loctmatsize]
                            iG = pw.gindex[ii,:]
                            ii = pw.indgq[iq,jm-loctmatsize]
                            jG = pw.gindex[ii,:]
                            print >> fout, ('%3d '*2+'%2d'*3+'   '+','+'   '+'%2d'*3+' %14.9f%14.9f') % (im+1,jm+1,iG[0],iG[1],iG[2], jG[0],jG[1],jG[2], Vmat[im,jm].real, Vmat[im,jm].imag)
        #print >> fout, '** barc eigenvalues **'
        #for i in range(len(self.ev)):
        #    print >> fout, '%5d%18.9f' % (i+1, self.ev[-1-i])

    def Coul_setev(self, iq, fout, evtol=1e-8):
        # evtol -- eigenvalues bigger than this cutoff will be kept
        mbsize = len(self.ev)
        if Real_sqrt_Vcoul:
            #??### LAST CHANGE:
            matsize = mbsize
            self.barcev = self.ev
            # Here we literaly compute sqrt(V_coulomb) using eigenvectors and eigenvalues
            # The hope is that such procedure leads to more stable and real values of V_coulomb, which do not have arbitrary complex phase due to diagonalization
            self.barcvm = la.matmul(self.Vmat, dot(diag(sqrt(abs(self.ev))), conj(self.Vmat.T)) )
            #print 'Imaginary part of Vcoul=', sum(abs(self.barcvm.imag))/(mbsize*mbsize)
            
            print >> fout, "  - Old/New basis set size =", mbsize, matsize
            print >> fout, '** barc eigenvalues **'
            for im in range(matsize-1,-1,-1):
                print >> fout, '%3d %14.8f' % (matsize-im, self.barcev[im]),'v=', ('%13.8f'*10) % tuple(self.barcvm[:10,im].real)
        else:
            im_kept = self.ev > evtol
            if iq==0:
                if io.iop_coul_x == 0:      # even if the eigenvalue is larger than cutooff, we will still remove it
                    im_kept[self.immax] = False
                else:                       # now if io.iop_coul_x !=0 then we will always keep this eigenvector
                    im_kept[self.immax] = True
            
            matsize = sum(1 for i in im_kept if i)
            
            self.barcev = zeros(matsize)
            self.barcvm = zeros((mbsize,matsize), dtype=complex, order='F')
            
            im=matsize-1
            for jm in range(mbsize):
                if im_kept[jm]:
                    self.barcev[im] = self.ev[jm]
                    self.barcvm[:,im] = self.Vmat[:,jm]*sqrt(abs(self.ev[jm]))
                    im -= 1
                    
            print >> fout, "  - Old/New basis set size =", mbsize, matsize
            print >> fout, '** barc eigenvalues **'
            for im in range(matsize):
                print >> fout, '%3d %14.8f' % (im+1, self.barcev[im]),'v=', ('%13.8f'*10) % tuple(self.barcvm[:10,im].real)

            
    def calc_minm(self, ik, iq, band_limits, mode, strc, in1, latgen, kqm, ks, radf, pw, pb, core, rest, t_times, fout, PRINT=False):
        #print 'mode=', mode, 'ik=', ik, 'band_limits=', band_limits
        (isp, indggq, nmix, max_nmix, big_l, ql) = rest
        (mbsize,matsize) = shape(self.barcvm)
        (nst, nend, mst, mend, cst, cend) = band_limits
        #(nst, nend, mst, mend, cst, cend) = (ks.ibgw, ks.nbgw, 0, nomx+1, 0, ks.ncg_x)
        irk = kqm.kii_ind[ik]
        jk = kqm.kqid[ik,iq]  # index of k-q
        jrk = kqm.kii_ind[jk]
        kpts = (ik,jk,iq)
        
        t0 = timer()
        # First create alpha, beta, gamma for k
        Aeigk = array(ks.all_As[irk], dtype=complex)   # eigenvector from vector file
        kil = array(kqm.klist[ik,:])/float(kqm.LCM)  # k in semi-cartesian form

        if kqm.k_ind[irk] != ik:                       # not irreducible
            Aeigk *= exp(-2*pi*1j * ks.phase_arg[ik][:])  # adding phase : zzk[1:ngk,ib] = phase[1:ngk] * zzk[1:ngk,ib]
        if DMFT1_like:
            kirr = array(kqm.kirlist[irk,:])/float(kqm.LCM)  # k in semi-cartesian form
            isym = kqm.iksym[ik]
            timat_ik, tau_ik = strc.timat[isym].T, strc.tau[isym,:]
            alm,blm,clm = lapwc.dmft1_set_lapwcoef(False, 1, True, kil, kirr,timat_ik,tau_ik, ks.indgkir[irk], ks.nv[irk], pw.gindex, radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.tauij, latgen.Vol, kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax)
        else:
            alm,blm,clm = lapwc.gap2_set_lapwcoef(kil, ks.indgk[ik], 1, True, ks.nv[irk], pw.gindex, radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.Vol, kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax)
            
        (ngi,nLOmax,ntnt,ndf) = shape(clm)
        (ngi,ntnt,ndf) = shape(alm)
        (nbmax,ngi2) = shape(Aeigk)
        # And now change alm,blm,clm to band basis, which we call alfa,beta,gama
        alfa = reshape( la.matmul(Aeigk, reshape(alm, (ngi,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
        beta = reshape( la.matmul(Aeigk, reshape(blm, (ngi,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
        if in1.nlomax > 0:
            gama = reshape( la.matmul(Aeigk, reshape(clm, (ngi,ntnt*ndf*nLOmax)) ), (nbmax,nLOmax,ntnt,ndf) )
        else:
            gama = zeros((1,1,1,1),dtype=complex,order='F') # can not have zero size. 

        # And next create alpha, beta, gamma for k+q
        Aeigq = array( conj( ks.all_As[jrk] ), dtype=complex)  # eigenvector from vector file
        kjl = array(kqm.klist[jk,:])/float(kqm.LCM)          # k+q in semi-cartesian form
        if kqm.k_ind[jrk] != jk:                             # the k-q-point is reducible, eigenvector needs additional phase
            Aeigq *= exp( 2*pi*1j * ks.phase_arg[jk][:] )
        if DMFT1_like:
            kjrr = array(kqm.kirlist[jrk,:])/float(kqm.LCM)
            isym = kqm.iksym[jk]
            timat_ik, tau_ik = strc.timat[isym].T, strc.tau[isym,:]
            alm,blm,clm = lapwc.dmft1_set_lapwcoef(False, 2, True, kjl, kjrr,timat_ik,tau_ik, ks.indgkir[jrk], ks.nv[jrk], pw.gindex, radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.tauij, latgen.Vol, kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax)
        else:
            alm,blm,clm = lapwc.gap2_set_lapwcoef(kjl, ks.indgk[jk], 2, True, ks.nv[jrk], pw.gindex, radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.Vol, kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax)
            
        (ngj,nLOmax,ntnt,ndf) = shape(clm)
        (ngj,ntnt,ndf) = shape(alm)
        (nbmax,ngj2) = shape(Aeigq)
        # And now change alm,blm,clm to band basis, which we call alfa,beta,gama
        alfp = reshape( la.matmul(Aeigq, reshape(alm, (ngj,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
        betp = reshape( la.matmul(Aeigq, reshape(blm, (ngj,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
        if in1.nlomax > 0:
            gamp = reshape( la.matmul(Aeigq, reshape(clm, (ngj,ntnt*ndf*nLOmax)) ), (nbmax,nLOmax,ntnt,ndf) )
        else:
            gamp = zeros((1,1,1,1),dtype=complex,order='F') # can not have zero size
        
        abc_lapw = (alfa,beta,gama,alfp,betp,gamp)
        Aeigs = (Aeigk,Aeigq)
        
        if False:
            print 'alfa,beta,gama='
            for ie in range(shape(alfa)[0]):
                for lm in range(shape(alfa)[1]):
                    print 'ie=%3d lm=%3d alfa=%14.10f%14.10f beta=%14.10f%14.10f gama=%14.10f%14.10f' % (ie+1,lm+1,alfa[ie,lm,0].real, alfa[ie,lm,0].imag, beta[ie,lm,0].real, beta[ie,lm,0].imag, gama[ie,0,lm,0].real, gama[ie,0,lm,0].imag)
            print 'alfp,betp,gamp='
            for ie in range(shape(alfp)[0]):
                for lm in range(shape(alfp)[1]):
                    print 'ie=%3d lm=%3d alfa=%14.10f%14.10f beta=%14.10f%14.10f gama=%14.10f%14.10f' % (ie+1,lm+1,alfp[ie,lm,0].real, alfp[ie,lm,0].imag, betp[ie,lm,0].real, betp[ie,lm,0].imag, gamp[ie,0,lm,0].real, gamp[ie,0,lm,0].imag)

        t1 = timer()
        t_times[0] += t1-t0

        # Next computing overlap of two Kohn-Sham orbitals and product basis
        s3r = pb.s3r[:,:,:,:,isp]
        ## The muffin-tin part : mmat[ie1,ie2,im] = < u^{product}_{im,lb} | psi^*_{ie2} psi_{ie1} > e^{-iqr}
        #    where ie1=[nst,nend] and ie2=[mst,mend]
        mmat_mt = fvxcn.calc_minm_mt(ql,nst,nend,mst,mend, alfa,beta,gama,alfp,betp,gamp,s3r,strc.vpos,strc.mult,nmix,big_l,in1.nLO_at,pb.ncore,pb.cgcoef,in1.lmax,self.loctmatsize)
        ## The interstitial part : mmat(ie1,ie2,im) = 1/sqrt(O)*< e^{i*iGp*r}|psi_{ie2}^* |psi_{ie1}>_{Interstitials}
        iumklap = array(round_( kqm.klist[ik]/float(kqm.LCM) - kqm.klist[jk]/float(kqm.LCM) - ql ),dtype=int)
        nvi, nvj = ks.nv[irk], ks.nv[jrk]
        # The interstitial part : mmat[ie1,ie2,im] = < u^{product}_{im,lb} | psi^*_{ie2} psi_{ie1} > = 1/sqrt(O)*< e^{i*iG_{im}*r}|psi_{ie2,k+q}^* |psi_{ie1,k}>_{Interstitials}
        #    where ie1=[nst,nend] and ie2=[mst,mend]
        mmat_it = fvxcn.calc_minm_is(nst,nend,mst,mend,Aeigk,Aeigq,iumklap,self.mpwipw,nvi,nvj,ks.indgk[ik],ks.indgk[jk],indggq,pw.gindex,self.i_g0,latgen.Vol)
        # Combine MT and interstitial part together
        mmat = concatenate( (mmat_mt, mmat_it), axis=-1 )
        t2 = timer()
        # Now compute also product of sqrt(Vcoulomb)*mmat
        minm = la.matmul( mmat, conj(self.barcvm) )
        #print 'mmat_mt=', isfortran(mmat_mt), 'mmat_it=', isfortran(mmat_it), 'mmat=', isfortran(mmat), 'minm=', isfortran(minm), 'barcvm=', isfortran(self.barcvm)
        t3 = timer()
        t_times[1] += t2-t1
        t_times[2] += t3-t2

        minc = None
        if cend > 0:  # at least one core state, hence compute overlap of Kohn-Sham+core orbitals and product basis
            t1 = timer()
            if mode=='selfe':
                # computing <Product_Basis| psi_{i1,k} psi_{icore}^*> with i1=[nst,nend] and icore=[cst,cend], where i1 is occupied
                mmat_c = fvxcn.calc_minc(kil,nst,nend,cst,cend,alfa,beta,gama,s3r,core.corind,strc.vpos,strc.mult,nmix,big_l,in1.nLO_at,pb.ncore,pb.cgcoef,in1.lmax,self.loctmatsize)
            else:
                # computing <Product_Basis| psi_{icore} psi_{i2,k-q}^* > with i2=[mst,mend] and icore=[cst,cent], where i2 is empty
                mmat_c = fvxcn.calc_minc2(kjl,cst,cend,mst,mend,alfp,betp,gamp,s3r,core.corind,strc.vpos,strc.mult,nmix,big_l,in1.nLO_at,pb.ncore,pb.cgcoef,in1.lmax,self.loctmatsize)
            t2 = timer()
            # Now compute also
            minc = la.matmul( mmat_c, conj(self.barcvm[:self.loctmatsize,:]) )
            #minc = dot( mmat_c, conj(self.barcvm[:self.loctmatsize,:]) )
            t3 = timer()
            t_times[3] += t2-t1
            t_times[4] += t3-t2
            #if mode == 'polarization':
            #    print 'time=', t2-t1
            #    for im in range(self.loctmatsize):
            #        for ic in range(cend-cst):
            #            for ie2 in range(mend-mst):
            #                print '%4d %4d %4d %20.14f %20.14f' % (im, ic, ie2, mmat_c[ic,ie2,im].real, mmat_c[ic,ie2,im].imag)
            #    
            #    sys.exit(0)
            
        #if PRINT: # Conclusion: you will need to check mmat, which is compatible, and can not check minm, because it is not unique
        #    print >> fout, 'mmat=', 'ik=', ik
        #    for imix in range(shape(mmat)[2]):
        #        for ie1 in range(shape(mmat)[0]):
        #            for ie2 in range(shape(mmat)[1]):
        #                print >> fout, '%3d %3d %3d   %12.8f%12.8f' % (imix+1, ie1+1, ie2+1, mmat[ie1,ie2,imix].real, mmat[ie1,ie2,imix].imag)

            
        return (minm, minc)
    
    def calc_selfx(self, iq, strc, in1, latgen, kqm, ks, radf, pw, pb, core, kw, fout):
        isp = 0
        
        indggq = pw.inverse_indgq(iq)
        # nomx+1 is the number of occupied valence bands, nomx is the last valence band
        nomx, numin = ks.nomax_numin
        # Here we need both all occuped valence bands (nomx+1) and also all the core states. Hence nbands = nomx+1 + ncore
        
        nirkp = len(kqm.weight)
        (mbsize,matsize) = shape(self.barcvm)
        nmix = array([ len(pb.big_l[iat]) for iat in range(strc.nat) ], dtype=int)        # nmix is the length of product basis on each atom
        max_nmix = max( nmix )
        big_l = zeros( (max_nmix,strc.nat), dtype=int, order='F' )                        # big_l in fortran array form
        for iat in range(strc.nat):
            big_l[:len(pb.big_l[iat]),iat] = pb.big_l[iat][:]                             # saving big_l to fortran-like array

        ql = kqm.qlistc[iq,:]/float(kqm.LCMq)
        
        print >> fout
        
        selfx = zeros((nirkp,ks.nbgw-ks.ibgw))
        t_selfx = 0
        t_times = zeros(5)
        for irk in range(nirkp):
            kl = array(kqm.kirlist[irk,:])/float(kqm.LCM)  # k in semi-cartesian form
            ik = kqm.k_ind[irk]   # index in all-kpoints, not irreducible
            jk = kqm.kqid[ik,iq]  # index of k-q
            jrk = kqm.kii_ind[jk]
            kpts = (ik,jk,iq)
            
            band_limits = (ks.ibgw, ks.nbgw, 0, nomx+1, 0, ks.ncg_x)
            rest = (isp, indggq, nmix, max_nmix, big_l, ql)
            minm, minc = self.calc_minm(ik, iq, band_limits, 'selfe', strc, in1, latgen, kqm, ks, radf, pw, pb, core, rest, t_times, fout)
            (nst, nend, mst, mend, cst, cend)  = band_limits
            
            t2 = timer()
            #   ms[ie1,ie2] = sum_{im} minm[ie1,ie2,im] * conj(minm[ie1,ie2,im])
            ms = sum(minm * conj(minm), axis=2).real
            #   msw[ie1] = sum_{ie2} ms[ie1,ie2]*kwgh[ie2]
            kwgh = kw.kiw[jrk, mst:mend] # tetrahedral weights for occupied valence bands
            sx = -dot(ms, kwgh)
            t3 = timer()
            t_selfx += t3-t2

            if cend > 0:  # at least one core state, hence compute overlap of Kohn-Sham+core orbitals and product basis
                t2 = timer()
                #   ms[ie1,ie2] = sum_{im} minm[ie1,ie2,im] * conj(minm[ie1,ie2,im])
                ms = sum(minc * conj(minc), axis=2).real
                # msw[ie1] = sum_{ie2} ms[ie1,ie2]*kwgh[ie2]
                kwgh = kw.kwt_ibz[jrk]*ones(ks.ncg_x)
                sx += -dot(ms, kwgh)
                t3 = timer()
                t_selfx += t3-t2
                
            if iq==0:  # correction for q==0
                ### THIS NEEDS TO BE STUDIES. IT DOES NOT MAKE SENSE FOR METALS!
                cst = 4*pi/latgen.Vol*self.singc2*len(kqm.qlist)
                sx[:nomx+1-nst] += -cst * kw.kiw[jrk, nst:nomx+1] # tetrahedral weights for occupied valence bands
            
            selfx[irk,:] = sx[:]
            
            if False:
                print >> fout, 'selfx:minm'
                for ie1 in range(nend-nst):
                    for ie2 in range(mend-mst):
                        print >> fout, '%3d %3d' % (ie1, ie2), '%12.8f'*matsize % tuple(minm[ie1,ie2,:].real)
                        print >> fout, '%3d %3d' % (ie1, ie2), '%12.8f'*matsize % tuple(minm[ie1,ie2,:].imag)
                print >> fout, 'minc'
                if cend > 0 :
                    for ie1 in range(nend-nst):
                        for ie2 in range(cend-cst):
                            print >> fout, '%3d %3d' % (ie1, ie2), '%12.8f'*matsize % tuple(minm[ie1,ie2,:].real)
                            print >> fout, '%3d %3d' % (ie1, ie2), '%12.8f'*matsize % tuple(minm[ie1,ie2,:].imag)

            for i in range(len(sx)):
                print >> fout, 'dSigx[irk=%3d,ie=%3d]=%16.12f' % (irk, i+ks.ibgw, sx[i])
            
        print >> fout, '## selfx  : t(prep_minm) [iq=%-3d]  =%14.9f' % (iq,t_times[0])
        print >> fout, '## selfx  : t(minm)      [iq=%-3d]  =%14.9f' % (iq,t_times[1])
        print >> fout, '## selfx  : t(minm*sV)   [iq=%-3d]  =%14.9f' % (iq,t_times[2])
        print >> fout, '## selfx  : t(minc)      [iq=%-3d]  =%14.9f' % (iq,t_times[3])
        print >> fout, '## selfx  : t(minc*sV)   [iq=%-3d]  =%14.9f' % (iq,t_times[4])
        print >> fout, '## selfx  : t(cmp_selfx) [iq=%-3d]  =%14.9f' % (iq,t_selfx)
        return selfx
    
    def calc_head(self, strc, in1, latgen, kqm, ks, radf, pw, pb, core, kw, fr, kcw, io_iop_drude, io_eta_head, dUl, Ul, fout):
        ### calcomat
        isp = 0
        iq = 0
        indggq = pw.inverse_indgq(iq)
        # nomx+1 is the number of occupied valence bands, nomx is the last valence band
        nomx, numin = ks.nomax_numin
        # Here we need both all occuped valence bands (nomx+1) and also all the core states. Hence nbands = nomx+1 + ncore
        
        nirkp = len(kqm.weight)
        (mbsize,matsize) = shape(self.barcvm)
        nmix = array([ len(pb.big_l[iat]) for iat in range(strc.nat) ], dtype=int)        # nmix is the length of product basis on each atom
        max_nmix = max( nmix )
        big_l = zeros( (max_nmix,strc.nat), dtype=int, order='F' )                        # big_l in fortran array form
        for iat in range(strc.nat):
            big_l[:len(pb.big_l[iat]),iat] = pb.big_l[iat][:]                             # saving big_l to fortran-like array

        ql = kqm.qlistc[iq,:]/float(kqm.LCMq)

        isp=0
        (iulol_ul, iulol_udl, iul_ulol, iudl_ulol, iulol_ulol) = ks.Give_fortran_ilocals(isp,strc,in1)
        iul_ul   = ks.iul_ul[:,:,:,isp]   #= zeros((2,in1.nt-1,strc.nat,nspin),order='F')
        iul_udl  = ks.iul_udl[:,:,:,isp]  #= zeros((2,in1.nt-1,strc.nat,nspin),order='F')
        iudl_ul  = ks.iudl_ul[:,:,:,isp]  #= zeros((2,in1.nt-1,strc.nat,nspin),order='F')
        iudl_udl = ks.iudl_udl[:,:,:,isp] #= zeros((2,in1.nt-1,strc.nat,nspin),order='F')
        
        nst,nend,mst,mend = 0,nomx+1,numin,ks.nbmaxpol
        ncg = len(core.corind)
        mmatvv = zeros((nirkp,mend-mst,nend-nst,3),dtype=complex)
        mmatcv = zeros((nirkp,mend-mst,ncg,3),dtype=complex)

        t_coeff = 0
        t_matvv, t_matcv = 0, 0
        for irk in range(nirkp):
            kl = array(kqm.kirlist[irk,:])/float(kqm.LCM)  # k in semi-cartesian form
            ik = kqm.k_ind[irk]   # index in all-kpoints, not irreducible
            jk = kqm.kqid[ik,iq]  # index of k-q
            jrk = kqm.kii_ind[jk]
            kpts = (ik,jk,iq)

            t0 = timer()
            # First create alpha, beta, gamma for k
            Aeigk = array(ks.all_As[irk], dtype=complex)   # eigenvector from vector file
            kil = array(kqm.klist[ik,:])/float(kqm.LCM)  # k in semi-cartesian form
            if DMFT1_like:
                isym = kqm.iksym[ik]
                timat_ik, tau_ik = strc.timat[isym].T, strc.tau[isym,:]
                alm,blm,clm = lapwc.dmft1_set_lapwcoef(False, 1, False, kil, kil, timat_ik,tau_ik, ks.indgkir[irk], ks.nv[irk], pw.gindex, radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.tauij, latgen.Vol, kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax)
            else:
                alm,blm,clm = lapwc.gap2_set_lapwcoef(kil, ks.indgk[ik], 1, False, ks.nv[irk], pw.gindex, radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.Vol, kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax)
            
            (ngi,nLOmax,ntnt,ndf) = shape(clm)
            (ngi,ntnt,ndf) = shape(alm)
            (nbmax,ngi2) = shape(Aeigk)
            # And now change alm,blm,clm to band basis, which we call alfa,beta,gama
            alfa = reshape( la.matmul(Aeigk, reshape(alm, (ngi,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
            beta = reshape( la.matmul(Aeigk, reshape(blm, (ngi,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
            if in1.nlomax > 0:
                gama = reshape( la.matmul(Aeigk, reshape(clm, (ngi,ntnt*ndf*nLOmax)) ), (nbmax,nLOmax,ntnt,ndf) )
            else:
                gama = zeros((1,1,1,1),dtype=complex,order='F') # can not have zero size. 
            
            t1 = timer()
            tw1,tw2 = 0, 0
            t_coeff += t1-t0
            for ie2 in range(mst,mend):
                t3 = timer()
                alfa_n2, beta_n2 = alfa[ie2,:,:],  beta[ie2,:,:]
                if in1.nlomax > 0:
                    gama_n2 = gama[ie2,:,:,:]
                else:
                    gama_n2 = gama[0,:,:,:]
                for ie1 in range(nst,nend):
                    #  < psi_{ie1}| (-i*\nabla) | psi_{ie2}>_MT
                    if False:
                        p_mt = f_q_0.calcmmatvv_mt(ie1+1,ie2+1,alfa,beta,gama,in1.nLO_at,strc.mult,iul_ul,iul_udl,iudl_ul,iudl_udl,iulol_ul,iulol_udl,iul_ulol,iudl_ulol,iulol_ulol)
                    else:
                        calfa_n1, cbeta_n1 = conj(alfa[ie1,:,:]), conj(beta[ie1,:,:])
                        if in1.nlomax > 0:
                            cgama_n1 = conj(gama[ie1,:,:,:])
                        else:
                            cgama_n1 = gama[0,:,:,:]
                        p_mt = f_q_0.calcmmatvv_mt2(calfa_n1,cbeta_n1,cgama_n1,alfa_n2,beta_n2,gama_n2,in1.nLO_at,strc.mult,iul_ul,iul_udl,iudl_ul,iudl_udl,iulol_ul,iulol_udl,iul_ulol,iudl_ulol,iulol_ulol)
                    
                    #  < psi_{ie1}| (-i*\nabla) | psi_{ie2}>_I = \sum_{G1,G2} A*_{G1,ie1} A_{G2,ie2} <G1+k| -i*\nabla |G2+k>_I =
                    p_in = f_q_0.calcmmatvv_is(ie1+1,ie2+1,kil, Aeigk, ks.indgk[ik], ks.nv[irk], pw.ipwint,pw.gindex,self.i_g0,kqm.k2cartes)
                    if (ie1 == ie2):
                        p12 = (p_mt+p_in).real
                    else:
                        #  < psi_{ie2}| (-i*\nabla) | psi_{ie1}>_MT
                        p_mth = f_q_0.calcmmatvv_mt(ie2+1,ie1+1,alfa,beta,gama,in1.nLO_at,strc.mult,iul_ul,iul_udl,iudl_ul,iudl_udl,iulol_ul,iulol_udl,iul_ulol,iudl_ulol,iulol_ulol)
                        #  < psi_{ie2}| (-i*\nabla) | psi_{ie1}>_I = \sum_{G1,G2} A*_{G1,ie2} A_{G2,ie1} <G1+k| -i*\nabla |G2+k>_I =
                        p_inh = f_q_0.calcmmatvv_is(ie2+1,ie1+1,kil, Aeigk, ks.indgk[ik], ks.nv[irk], pw.ipwint,pw.gindex,self.i_g0,kqm.k2cartes)
                        p12 = (p_mt+p_in+conj(p_mth+p_inh))/2.
                    mmatvv[irk,ie2-mst,ie1-nst,:] = p12[:]
                    p2 = dot(p12, conj(p12)).real/3.
                    #print '%3d %3d' % (ie2+1, ie1+1), ('%16.12f %16.12f  '*3+'  %16.12f') % (p12[0].real, p12[0].imag, p12[1].real, p12[1].imag, p12[2].real, p12[2].imag, p2)
                t4 = timer()
                t_matvv += t4-t3
                if ncg > 0:
                    #  ( < psi_{ie2}| (-i*\nabla) | u_core>_MT + < u_core | (-i*\nabla) | psi_{ie2}>^*_MT )/2
                    ic_start=0
                    for iat in range(strc.nat):
                        if not core.l_core[iat]: # this atom has no core
                            continue
                        ic_end = ic_start + len(core.l_core[iat])
                        iul_ucl  = ks.iul_ucl[isp][iat]
                        iudl_ucl = ks.iudl_ucl[isp][iat]
                        iucl_ul  = ks.iucl_ul[isp][iat]
                        iucl_udl = ks.iucl_udl[isp][iat]
                        iulol_ucl = ks.iulol_ucl[isp][iat]
                        iucl_ulol = ks.iucl_ulol[isp][iat]
                        mmatcv[irk,ie2-mst,ic_start:ic_end,:] += f_q_0.calcmmatcv(iat+1,ie2+1,core.corind,alfa,beta,gama,iul_ucl,iudl_ucl,iucl_ul,iucl_udl,iulol_ucl,iucl_ulol,in1.nLO_at)
                        ic_start = ic_end
                    if ncg != ic_end:
                        print 'ERROR : missmatch with core states ic_end='+str(ic_end)+' while ncg='+str(ncg)
                        sys.exit(1)
                    for icg in range(ncg):
                        p12 = mmatcv[irk,ie2-mst,icg,:]
                        p2 = dot(p12, conj(p12)).real/3.
                        #print >> fout, '%3d %3d' % (ie2+1, icg+1), ('%16.12f %16.12f  '*3+'  %16.12f') % (p12[0].real, p12[0].imag, p12[1].real, p12[1].imag, p12[2].real, p12[2].imag, p2)
                t5 = timer()
                t_matcv += t5-t4
        print >> fout, '## cal_head t(ABC coeff)           =%14.9f' % t_coeff
        print >> fout, '## cal_head t(matvv)               =%14.9f' % t_matvv
        print >> fout, '## cal_head t(matcv)               =%14.9f' % t_matcv

        ### calchead
        t6 = timer()
        # Dielectric constant at Gamma-point
        c0_head = 0.   # 4*pi/Vol * Nspin* \sum_{i,k} 1/3*<psi_i| p |psi_i>^2 * delta(e_{k,i}-EF)
        ankp = kqm.ndiv[0]*kqm.ndiv[1]*kqm.ndiv[2]
        fspin = 2. # because we do not  sum  over spins
        (nb1,nb2,nkp,nom_nl) = shape(kcw) # the last index could be iomega or svd l
        head = zeros(nom_nl, dtype=complex)
        for irk in range(nirkp):
            ik = kqm.k_ind[irk]   # index in all-kpoints, not irreducible
            nvbm, ncbm = ks.nomax_numin[0]+1, ks.nomax_numin[1]
            kwt = kqm.weight[irk]
            coef = 4.*pi*kwt*fspin/latgen.Vol
            if (ks.Eg <= 0): # metallic
                # correction for metallic state
                for ie in range(ncbm,nvbm):
                    # 4*pi/Vol * Nspin*kw * 1/3*sum_{i = cross_EF} <psi_i| p |psi_i>^2
                    pnmkq2 = vdot(mmatvv[irk,ie-mst,ie-nst,:],mmatvv[irk,ie-mst,ie-nst,:]).real/3.0
                    # print 'pnmkq2=', pnmkq2, 'kwfer=', kw.kwfer[irk,ie]
                    c0_head += coef * pnmkq2 * kw.kwfer[irk,ie]
            
            #print >> fout, 'irk=', irk, 'c0_head=', c0_head
            
            for icg in range(ks.ncg_p):
                # core contibution
                iat,idf,ic,l,m = core.corind[icg,:]
                edif = ks.Ebnd[isp,irk,ncbm:ks.nbmaxpol] - core.eig_core[isp][iat][ic]
                mm = mmatcv[irk,(ncbm-mst):(ks.nbmaxpol-mst),icg,:]  # mm[irk,icg][ie,3] == < u_core| -i*\nabla | psi_{irk,ie}>
                mm2 = sum(mm*conj(mm), axis=1).real/3.               # mm2[ie] = 1/3. * sum_{i} |mm[ie,i=1:3]|^2=1/3(|<u_core|p_x|psi_{irk,ie}>|^2+|<u_core| p_y|psi_{irk,ie}>|^2+|<u_core| p_z|psi_{irk,ie}>|^2)
                termcv = mm2/edif**2
                _kcw_ = kcw[icg,:(ks.nbmaxpol-ncbm),ik,:]  # Polarization_weight[icg,ik][ie,iom]
                #_kcw_ = kcw[icg,ncbm:ks.nbmaxpol,ik,:]  # Polarization_weight[icg,ik][ie,iom]
                head -= coef*conj(dot(termcv,_kcw_))    # head[iom] -= 4*pi/Vol*Nspin \sum_{k,i,j} 1/3|<u_i|p|u_j>|^2/(E_{k,i}-E_{k,j})^2 * Polarization_weight[i,j,k,omega]
                        
            
            # valence contribution
            mm2 = sum(abs(mmatvv[irk,(ncbm-mst):(ks.nbmaxpol-mst),0:(nvbm),:])**2,axis=2)/3.
            termvv = zeros((ks.nbmaxpol-ncbm,nvbm))
            for ie2 in range(ks.nbmaxpol-ncbm):
                for ie1 in range(nvbm):
                    edif = ks.Ebnd[isp,irk,ncbm+ie2]-ks.Ebnd[isp,irk,ie1]
                    if abs(edif) < 1e-5:
                        #print >> fout, 'WARNING in calc_head : degenerate CB and VB'
                        termvv[ie2,ie1] = 0
                    else:
                        termvv[ie2,ie1] = mm2[ie2,ie1]/edif**2
            
            for iom_il in range(nom_nl):
                _kcw_ = kcw[ks.ncg_p:(ks.ncg_p+nvbm),:(ks.nbmaxpol-ncbm),ik,iom_il]  # Polarization_weight[ik,iom][ie1,ie2]
                #_kcw_ = kcw[ks.ncg_p:(ks.ncg_p+nvbm),ncbm:ks.nbmaxpol,ik,iom].real  # Polarization_weight[ik,iom][ie1,ie2]
                sm = sum(termvv * _kcw_.T) # sm = sum_{i,j} 1/3|<u_i|p|u_j>|^2/(E_{k,i}-E_{k,j})^2 * Polarization_weight[i,j,k,omega]
                head[iom_il] -= coef*sm       # head -= 4*pi/Vol*Nspin \sum_{k,i,j} 1/3|<u_i|p|u_j>|^2/(E_{k,i}-E_{k,j})^2 * Polarization_weight[i,j,k,omega]
                #print '%3d %s %16.10f %16.10f' % (iom+1, 'dhead=', sm, -coef*sm)
            
            #print >> fout, 'irk=', irk, 'head=', head
            if False:
                iom=2
                _kcw_ = kcw[ks.ncg_p:(ks.ncg_p+nvbm),:(ks.nbmaxpol-ncbm),ik,iom]  # [ie1,ie2]
                #_kcw_ = kcw[ks.ncg_p:(ks.ncg_p+nvbm),ncbm:ks.nbmaxpol,ik,iom].real  # [ie1,ie2]
                sm0 = sum(termvv * _kcw_.T)                                          # [ie1,ie2]
                sm1 = sum(termvv)
                sm2 = sum(_kcw_)
                print >> fout, 'sm0=', sm0
                print >> fout, 'sm1=', sm1
                print >> fout, 'sm2=', sm2
        
        if ks.Eg <= 0 and io_iop_drude==1 :
            # Add plasmon contribution
            print >> fout, " Intraband contribution : Calc. plasmon freq. (eV):", sqrt(c0_head)*H2eV
            wpl2 = c0_head

            # on imaginary axis
            if dUl is not None:
                head += dot( wpl2/(fr.omega + io_eta_head)**2, dUl)
            else:
                head += wpl2/(fr.omega + io_eta_head)**2
            # on real axis
            #head -= wpl2/(fr.omega * (fr.omega + io_eta_head*1j))
        t7 = timer()
        print >> fout, '## calcomat t(dielectric head)     =%14.9f' % (t7-t6)

        if dUl is not None:
            head_om = dot(head,Ul)
            for iom in range(len(head_om)):
                print >> fout, '%3d head= %16.10f%16.10f' % (iom+1, head_om[iom].real+1.0, head_om[iom].imag)
        else:
            for iom_il in range(nom_nl):
                print >> fout, '%3d head= %16.10f%16.10f' % (iom_il+1, head[iom_il].real+1.0, head[iom_il].imag)
        # Note that in this implementation we need to add 1.0 to head in frequency domain
        return (head, mmatcv, mmatvv, mst)
    
    def calc_eps(self, iq, head_quantities, strc, in1, latgen, kqm, ks, radf, pw, pb, core, kw, fr, kcw, ddir, dUl, Ul, fout):
        isp = 0
        fspin = 2.
        coefw = 2.*sqrt(pi/latgen.Vol)
        if iq==0:
            (head, mmatcv, mmatvv, mst) = head_quantities
            # Converting head from svd to frequency
            if Ul is not None:
                head = dot(head,Ul)
        
        ql = kqm.qlistc[iq,:]/float(kqm.LCMq)
        indggq = pw.inverse_indgq(iq)
        nmix = array([ len(pb.big_l[iat]) for iat in range(strc.nat) ], dtype=int)        # nmix is the length of product basis on each atom
        max_nmix = max( nmix )
        big_l = zeros( (max_nmix,strc.nat), dtype=int, order='F' )                        # big_l in fortran array form
        for iat in range(strc.nat):
            big_l[:len(pb.big_l[iat]),iat] = pb.big_l[iat][:]                             # saving big_l to fortran-like array
        
        (mbsize,matsize) = shape(self.barcvm)
        (nb1_kcw,nb2_kcw,nkp_kcw,nom_nil) = shape(kcw)
        
        #len(fr.omega)
        eps   = zeros((matsize,matsize,nom_nil), dtype=complex, order='F')
        epsw1 = zeros((matsize,nom_nil), dtype=complex)
        #epsw2 = zeros((matsize,len(fr.omega)), dtype=complex) == conj(epsw1)
        
        t_times= zeros(10)
        for ik,irk in enumerate(kqm.kii_ind):
            coef = -fspin 
            jk = kqm.kqid[ik,iq]
            jrk = kqm.kii_ind[jk]
            nvbm, ncbm = ks.nomax_numin[0]+1, ks.nomax_numin[1]
            # set the local array for band energies  
            enk = zeros( ks.ncg_p + ks.nbmaxpol)
            for icg in range(ks.ncg_p):
                iat,idf,ic = core.corind[icg,:3]
                enk[icg] = core.eig_core[isp][iat][ic]
            enk[ks.ncg_p:(ks.ncg_p+ks.nbmaxpol)] = ks.Ebnd[isp,irk,:ks.nbmaxpol] 
            # find the index for the band whose energy is higher than eminpol  
            
            band_limits = (0, nvbm, ncbm, ks.nbmaxpol, 0, ks.ncg_p)
            rest = (isp, indggq, nmix, max_nmix, big_l, ql)
            minm0, minc0 = self.calc_minm(ik, iq, band_limits, 'polarization', strc, in1, latgen, kqm, ks, radf, pw, pb, core, rest, t_times, fout)
            if minc0 is not None:
                minm = concatenate( (minc0, minm0), axis=0 ) # minc0[core-states,unoccupied-bands,product-basis], minm0[occupied-bands,unoccupied-bands,product-basis]
            else:
                minm = minm0
            
            nb1,nb2,matsiz = shape(minm)
            Ndouble_bands = (ks.ncg_p+nvbm)*(ks.nbmaxpol-ncbm)
            
            if (matsiz != matsize):
                print 'Error : matsize=', matsize, 'and matsiz=', matsiz
            if False:
                print >> fout, 'shape(minm)=', shape(minm), 'shape(minc0)=', shape(minc0), 'shape(minm0)=', shape(minm0)
                print >> fout, 'ik=', ik, 'calceps:minm'
                for imix in range(matsiz):
                    for ie1 in range(nb1):
                        for ie2 in range(nb2):
                            print >> fout, '%3d %3d %3d %16.12f%16.12f' % (imix+1, ie1+1, ie2+ncbm+1, minm[ie1,ie2,imix].real, minm[ie1,ie2,imix].imag)
                print >> fout, 'calceps:minm_end'
                
            t10 = timer()
            if iq==0:
                # this could be done for irreducible k-points only, hence this can be optimized
                # Calculate pm(ie12) = pm(ie1,ie2) = p_{ie1,ie2}/(e_2-e_1)
                # Needed for the wings
                mmc = sum(mmatcv[irk,(ncbm-mst):(ks.nbmaxpol-mst),:,:],axis=2)*(coefw/sqrt(3.)) # ~<psi_{i1}|p^2|psi_{i2}>
                mmv = sum(mmatvv[irk,(ncbm-mst):(ks.nbmaxpol-mst),:,:],axis=2)*(coefw/sqrt(3.))
                pm = zeros((ks.ncg_p+nvbm)*(ks.nbmaxpol-ncbm), dtype=complex) # ~<psi_{i1}|p^2|psi_{i2}>/dE
                ie12=0
                for ie1 in range(ks.ncg_p):        # core-occupied
                    for ie2 in range(ncbm,ks.nbmaxpol): # empty
                        edif = enk[ie1] - enk[ie2+ks.ncg_p]
                        if abs(edif) > 1e-10:
                            pm[ie12] = mmc[ie2-ncbm,ie1]/edif
                        ie12 += 1
                for ie1 in range(nvbm):        # valence-occupied
                    for ie2 in range(ncbm,ks.nbmaxpol): # empty
                        edif = enk[ie1+ks.ncg_p] - enk[ie2+ks.ncg_p]
                        if abs(edif) > 1e-10:
                            pm[ie12] = mmv[ie2-ncbm,ie1]/edif
                        ie12 += 1
                
                if False:
                    print 'ik=', ik+1, 'pm='
                    for ie in range(ie12):
                        print >> fout, '%4d  %18.14f%18.14f' % (ie+1, pm[ie].real, pm[ie].imag)
                    sys.exit(0)
            t11 = timer()
            t_times[5] += t11-t10
            
            if False:
                minm2 = reshape(minm, (nb1*nb2,matsiz))
                print >> fout, 'minm2='
                for ie12,(ie1,ie2) in enumerate(itertools.product(range(ks.ncg_p+nvbm),range(ks.nbmaxpol-ncbm))):
                    for im in range(matsiz):
                        print >> fout, ie12, im, minm2[ie12,im]
                sys.exit(0)
                
            if ForceReal:
                "This is dangerous: We are truncating minm matrix elements to real components only"
                print 'Imaginary part of minm=', sum(abs(minm.imag))/(matsiz*nb1*nb2), ' is set to zero!'
                minm = array(minm.real, dtype=float)
                
            t12 = timer()
            tmat  = zeros( (matsiz,Ndouble_bands), dtype=minm.dtype)
            
            # minm[ie12,im] = < u^{product}_{im} | psi^*_{ie2} psi_{ie1} >
            minm2 = reshape(minm, (nb1*nb2,matsiz))
            #print 'shape(minm2)=', shape(minm2), 'matsiz=', matsiz, 'Ndouble_bands=', Ndouble_bands
            cminm2 = conj(minm2)
            minm2Tc = minm2.T*coef
            t_times[6] += timer()-t12
            
            #nvbm, ncbm = ks.nomax_numin[0]+1, ks.nomax_numin[1]
            #kcw = zeros( (ks.ncg+nvbm,ks.nbmaxpol-ncbm,nkp,nom) )
            #kcw[:(ks.ncg_p+nvbm),:(ks.nbmaxpol-ncbm),ik,iom]
            
            for iom_il in range(nom_nil):
                t13 = timer()
                # tmat = < u^{product}_{im} | psi^*_{ie2} psi_{ie1} > F[ie12,om]
                tmat = minm2Tc * kcw[:(ks.ncg_p+nvbm),:(ks.nbmaxpol-ncbm),ik,iom_il].ravel() #  tmat[matsiz,ie12] = transpose(minm2[ie12,matsiz])*_kcw_[ie12]
                #tmat = minm2Tc * kcw[:(ks.ncg_p+nvbm),ncbm:ks.nbmaxpol,ik,iom].ravel() #  tmat[matsiz,ie12] = transpose(minm2[ie12,matsiz])*_kcw_[ie12]
                t14 = timer()
                eps[:,:,iom_il] += la.matmul( tmat, cminm2 )
                #eps[:,:,iom_il] += matmul( tmat, cminm2 )
                #eps[:,:,iom_il] += dot( tmat, cminm2 )
                #linalg.blas.zgemm(1.0, tmat, cminm2, beta=1.0, c=eps[:,:,iom_il], overwrite_c=True)
                #print 'minm2=', isfortran(minm2), 'tmat=', isfortran(tmat), 'eps=', isfortran(eps[:,:,iom]), 'minm=', isfortran(minm),'minm0=',isfortran(minm0),'minc0=',isfortran(minc0)
                
                t15 = timer()
                if iq==0: 
                    wtmp = dot(tmat,pm)
                    epsw1[:,iom_il] += wtmp
                    #epsw2[:,iom_il] += conj(wtmp)
                    #print 'ik=%3d  iom=%3d  wtmp^2=%18.14f' % (ik+1,iom+1, sum(abs(wtmp)**2) )
                t16 = timer()
                t_times[7] += t14-t13
                t_times[8] += t15-t14
                t_times[9] += t16-t15
                
        if Save_eps:
            save(ddir+'/eps.'+str(iq), eps)    
            
        if iq==0:
            if Ul is not None:
                # converting epsw1 to frequency, if in svd basis
                epsw1 = dot(epsw1,Ul)
            epsw2 = copy(conj(epsw1))
        
        epsd = diagonal(eps, axis1=0, axis2=1).T.copy()  # epsd[i,iom_il] = eps[i,i,iom_il]
        if Ul is not None:
            epsd = dot(epsd,Ul) # convert to frequency when in svd basis
        epw1=0
        for iom in range(len(fr.omega)):
            #epst = sum([eps[i,i,iom].real for i in range(matsiz)])
            epst = sum(epsd[:,iom].real)
            if iq==0:
                epw1 = sum(abs(epsw1[:,iom])**2)
            print >> fout, 'Tr(epsilon[%11.6f]) = %18.14f  epsw1^2=%18.14f' % (fr.omega[iom], epst, epw1)
        
        t17 = timer()
        PRINT = False
        if Ul is not None:
            # converting the intire eps to frequency, when in svd basis
            eps = dot(eps,Ul)
        t18 = timer()
        
        emac = zeros((len(fr.omega),2),dtype=complex)
        Id = identity(shape(eps)[0])
        for iom in range(len(fr.omega)):
            # force hermicity, which should be obeyed
            eps[:,:,iom] = (eps[:,:,iom] + conj(eps[:,:,iom].T))/2.
            
            # above we just computed V*polarization, while epsilon = 1+VP
            eps[:,:,iom] += Id
            
            if PRINT:
                eigs = linalg.eigvalsh(eps[:,:,iom])
                print >> fout, 'iom=%3d  eps Eigenvalues=' % (iom+1,)
                for i in range(len(eigs)-1,-1,-1):
                    print >> fout, '%3d %14.10f' % (i+1, eigs[i])
            # this is W/V = (1+VP)^{-1}
            eps[:,:,iom] = linalg.inv(eps[:,:,iom])
            
            if iq==0:
                bw1 =      dot(eps[:,:,iom], epsw1[:,iom])
                w2b = conj(dot(eps[:,:,iom], epsw2[:,iom]))
                ww = dot(epsw2[:,iom],bw1)
                
                emac[iom,0] = 1.0+head[iom]-ww
                emac[iom,1] = 1.0+head[iom]

                head[iom] = 1/(1+head[iom]-ww)
                
                print >> fout, 'iom=%2d ww=%13.10f new head=%16.10f emac=[%16.10f,%16.10f]' % (iom+1,ww.real,head[iom].real,emac[iom,0].real, emac[iom,1].real)
                
                epsw1[:,iom] = -head[iom]*bw1[:]
                epsw2[:,iom] = -head[iom]*w2b[:]
                
                eps[:,:,iom] += head[iom]*tensordot(bw1, w2b, axes = 0)

                head[iom] -= 1.0
            else:
                head, epsw1, epsw2 = None, None, None
            
            # Now we compute 1/(1+VP)-1
            eps[:,:,iom] -= Id
            
            wst = sum([eps[i,i,iom].real for i in range(matsiz)])
            print >> fout, 'Tr(1/eps[%3d]-1) = %18.14f' % (iom+1, wst)
        t19 = timer()

        if dUl is not None:
            # converting back to svd basis, if svd basis is available
            eps   = dot(eps,dUl)
            if iq==0:
                epsw1 = dot(epsw1,dUl)
                epsw2 = dot(epsw2,dUl)
                head  = dot(head,dUl)
            
        t20 = timer()
        print >> fout, 'shape(minm)=', shape(minm), '=(nb1,nb2,matsiz)'
        print >> fout, '## calc_eps: t(prep_minm) [iq=%-3d] =%14.9f' % (iq,t_times[0])
        print >> fout, '## calc_eps: t(minm)      [iq=%-3d] =%14.9f' % (iq,t_times[1])
        print >> fout, '## calc_eps: t(minm*sV)   [iq=%-3d] =%14.9f' % (iq,t_times[2])
        print >> fout, '## calc_eps: t(minc)      [iq=%-3d] =%14.9f' % (iq,t_times[3])
        print >> fout, '## calc_eps: t(minc*sV)   [iq=%-3d] =%14.9f' % (iq,t_times[4])
        print >> fout, '## calc_eps: t(wings)     [iq=%-3d] =%14.9f' % (iq,t_times[5])
        print >> fout, '## calc_eps: t(minm2Tc)   [iq=%-3d] =%14.9f' % (iq,t_times[6])
        print >> fout, '## calc_eps: t(minm*kcw)  [iq=%-3d] =%14.9f' % (iq,t_times[7])
        print >> fout, '## calc_eps: t(tmat*minm) [iq=%-3d] =%14.9f' % (iq,t_times[8])
        print >> fout, '## calc_eps: t(wings2)    [iq=%-3d] =%14.9f' % (iq,t_times[9])
        print >> fout, '## calc_eps: t(svd2iom)   [iq=%-3d] =%14.9f' % (iq,t18-t17)
        print >> fout, '## calc_eps: t(inverse)   [iq=%-3d] =%14.9f' % (iq,t19-t18)
        print >> fout, '## calc_eps: t(iom2svd)   [iq=%-3d] =%14.9f' % (iq,t20-t19)
        return (eps, epsw1, epsw2, head)

    def calc_selfc(self, iq, eps, epsw1, epsw2, head, strc, in1, latgen, kqm, ks, radf, pw, pb, core, kw, fr, ddir, dUl, Ul, fout, PRINT=True):
        isp = 0
        ql = kqm.qlistc[iq,:]/float(kqm.LCMq)
        indggq = pw.inverse_indgq(iq)
        nmix = array([ len(pb.big_l[iat]) for iat in range(strc.nat) ], dtype=int)        # nmix is the length of product basis on each atom
        max_nmix = max( nmix )
        big_l = zeros( (max_nmix,strc.nat), dtype=int, order='F' )                        # big_l in fortran array form
        for iat in range(strc.nat):
            big_l[:len(pb.big_l[iat]),iat] = pb.big_l[iat][:]                             # saving big_l to fortran-like array
        
        if iq==0:
            coefs2 = self.singc2*4*pi/latgen.Vol
            coefs1 = self.singc1*sqrt(4*pi/latgen.Vol)
        wkq = 1.0/len(kqm.qlist)
        
        nirkp = len(kqm.weight)
        
        if False and (dUl is None):
            for iom in range(len(fr.omega)):
                ee = linalg.eigvalsh(eps[:,:,iom])
                print >> fout, 'iom=%3d  eps Eigenvalues=' % (iom+1,)
                for i in range(len(ee)):
                    print >> fout, '%3d %14.10f' % (i+1,ee[i])
            (matsiz,matsize,nom) = shape(eps)
            print >> fout, 'eps='
            for iom in range(len(fr.omega)):
                for i in range(matsiz):
                    for j in range(matsize):
                        print >> fout, '%4d %4d %4d %18.14f%18.14f' % (iom+1, i+ks.ibgw+1, j+1, eps[i,j,iom].real, eps[i,j,iom].imag)

        (matsiz1,matsiz2,nom_nil) = shape(eps)
        sc_p = zeros( (nirkp, ks.nbgw-ks.ibgw, len(fr.omega) ), dtype=complex )
        t_times = zeros(8)
        for irk in range(nirkp):
            kl = array(kqm.kirlist[irk,:])/float(kqm.LCM)  # k in semi-cartesian form
            ik = kqm.k_ind[irk]   # index in all-kpoints, not irreducible
            jk = kqm.kqid[ik,iq]  # index of k-q
            jrk = kqm.kii_ind[jk]
            kpts = (ik,jk,iq)
            
            band_limits = (ks.ibgw, ks.nbgw, 0, ks.nbands_c-ks.ncg_c, 0, ks.ncg_c)
            rest = (isp, indggq, nmix, max_nmix, big_l, ql)
            minm0, minc0 = self.calc_minm(ik, iq, band_limits, 'selfe', strc, in1, latgen, kqm, ks, radf, pw, pb, core, rest, t_times, fout)
            (nst, nend, mst, mend, cst, cend)  = band_limits
            nmdim = (band_limits[1]-band_limits[0])*(band_limits[3]-band_limits[2]+band_limits[5]-band_limits[4])
            if minc0 is not None:
                minm = concatenate( (minc0, minm0), axis=1 )  # putting core and valence part together shape(minim)=(nb1,nb2,matsiz), where nb1=nbgw-ibgw, nb2=core+valence
            else:
                minm = minm0
            nb1,nb2,matsiz = shape(minm)
            minm2 = reshape(minm, (nb1*nb2,matsiz)) # now minm2[nb1*nb2,matsiz]
            cminm2 = conj(minm2)
            
            t_i1 = timer()
            mwm = zeros((nom_nil,nb1,nb2), dtype=complex)
            for iom_il in range(nom_nil):
                ## mwm[iom,ib1,ib2] = wkq * \sum_{im1,im2} minm2[nb1*nb2,im1]^* eps[im1,im2,iom] minm2[nb1*nb2,im2]
                mw  = la.matmul(cminm2,eps[:,:,iom_il]) # mw[nb1*nb2,matsiz] = minm2*[nb1*nb2,matsiz] * eps[matsiz,matsiz]
                mwm0 = reshape(wkq*sum(mw * minm2, axis=-1),(nb1,nb2))   # sum(mw[nb1*nb2,matsiz] * minm[nb1*nb2,matsiz],axis=1)
                mwm[iom_il,:,:] = mwm0
                if (iq==0):
                    # mwm[iom,ib1,ib2] += coefs2 head[iom] + coefs1 \sum_{im} minm[ib1,ib2,im] eps2[im,iom] + coefs1 \sum_{im} minm[ib1,ib2,im]^* epsw1[im,iom]
                    for ie1 in range(ks.ibgw, ks.nbgw):
                        ie2 = ie1 + ks.ncg_c
                        d2 = coefs2*head[iom_il] + coefs1*dot(minm[ie1-ks.ibgw,ie2,:],epsw2[:,iom_il]) + coefs1*vdot(minm[ie1-ks.ibgw,ie2,:],epsw1[:,iom_il])
                        mwm[iom_il,ie1-ks.ibgw,ie2] += d2
            t_i2 = timer()
            t_times[5] += t_i2-t_i1
            
            save(ddir+'/mwm.'+str(iq)+'.'+str(irk), mwm)
            
            t_i3 = timer()
            t_times[6] += t_i3-t_i2
            
            #if dUl is not None:
            #    enk = zeros( ks.nbands_c)
            #    for icg in range(ks.ncg_c):
            #        iat,idf,ic = core.corind[icg,:3]
            #        enk[icg] = core.eig_core[isp][iat][ic]
            #    enk[ks.ncg_c:ks.nbands_c] = ks.Ebnd[isp,jrk,:(ks.nbands_c-ks.ncg_c)]   # careful, we need G(k+q,:) here, not G(k,:)
            #    
            #    (nl,nom) = shape(Ul)
            #    # Cnl[l,ie2,iom] = 1/pi*Integrate[ Ul[l,iOm]*(enk[ie2]-iom)/((enk[ie2]-iom)^2 + Om^2),{Om,0,Infinity}]
            #    Cnl = fnc.many_fr_convolutions(enk,Ul,fr.omega,fr.womeg)
            #    (nl,nb2,nom) = shape(Cnl)  # note: nb2==ks.nbands_c
            #    Cnl = reshape(Cnl,(nl*nb2,nom))
            #    mwm2 = zeros( (ks.nbgw-ks.ibgw, nl*nb2), dtype=complex )
            #    for ie1 in range(ks.nbgw-ks.ibgw):
            #        mwm2[ie1,:] = reshape( mwm[:,ie1,:], nl*nb2 )
            #    sc_p[irk,:,:] += dot(mwm2,Cnl)
            #    
            #else:
            #    sc_p[irk,:,:] += mcmn.Compute_selfc_inside(iq, irk, ks.Ebnd[isp], mwm, fr.omega, fr.womeg, kqm, ks.ncg_c, core, Ul, fout, PRINT=True)
            sc_p[irk,:,:] += mcmn.Compute_selfc_inside(iq, irk, ks.Ebnd[isp], mwm, fr, kqm, ks.ncg_c, core, Ul, fout)
                
            t_i4 = timer()
            t_times[7] += t_i4-t_i3

            PRINT=False
            if PRINT:
                for ie1 in range(ks.nbgw-ks.ibgw):
                    for iom in range(len(fr.omega)):
                        print >> fout, 'dSigc[irk=%3d,ie=%3d,iom=%3d]=%16.12f%16.12f' % (irk,ie1+ks.ibgw,iom,sc_p[irk,ie1,iom].real,sc_p[irk,ie1,iom].imag)

            
        print >> fout, '## calc_selfc: t(prep_minm)[iq=%-3d]=%14.9f' % (iq,t_times[0])
        print >> fout, '## calc_selfc: t(minm)    [iq=%-3d] =%14.9f' % (iq,t_times[1])
        print >> fout, '## calc_selfc: t(minm*sV) [iq=%-3d] =%14.9f' % (iq,t_times[2])
        print >> fout, '## calc_selfc: t(minc)    [iq=%-3d] =%14.9f' % (iq,t_times[3])
        print >> fout, '## calc_selfc: t(minc*sV) [iq=%-3d] =%14.9f' % (iq,t_times[4])
        print >> fout, '## calc_selfc: t(mwm)     [iq=%-3d] =%14.9f' % (iq,t_times[5])
        print >> fout, '## calc_selfc: t(enk      [iq=%-3d] =%14.9f' % (iq,t_times[6])
        print >> fout, '## calc_selfc: t(convol)  [iq=%-3d] =%14.9f' % (iq,t_times[7])
        return (sc_p)
    
    
def Polarization_weights(iq, ks, kqm, core, fr, iop_bzintq, sgnfrq, dUl, fout, isp=0):
    nomx, numin = ks.nomax_numin
    ankp = len(kqm.kii_ind)
    nbnd = ks.ncg_p + ks.nbmaxpol
    enk = zeros((nbnd, ankp), order='F')
    for icg in range(ks.ncg_p):
      iat, idf, ic = core.corind[icg][0:3]
      enk[icg,:] = core.eig_core[isp][iat][ic]
    for ik in range(ankp):
        irk = kqm.kii_ind[ik]
        enk[ks.ncg_p:(ks.ncg_p+ks.nbmaxpol),ik] = ks.Ebnd[isp,irk,:ks.nbmaxpol]
    
    # iop_bzintq = io.iop_bzintq # could be [-1:precise_numeric_integration,0:analytic_defaul,1:discrete_sum_of_poles_for_testing]
    # sgnfrq = io.fflg
    #omgmax_ht = 4.0   # (only relevant for iop_bzintq==-2) the upper bound for the HT integration
    #nomg_ht = 1000    # (only relevant for iop_bzintq==-2) the number of freq points used for the Hilbert transform (HT)
    #beta = 100.       # (only relevant for iop_bzintq== 1)
    #kcw = ft.bz_calcqdepw(iop_bzintq, sgnfrq, enk, ks.EF, kqm.kqid[:,iq], ks.ncg_p, nomx, numin, fr.omega, kqm.tetc, omgmax_ht, nomg_ht, beta)
    
    if dUl is not None:
        kcw = ft.bz_calcqdepw_par3(enk, ks.EF, kqm.kqid[:,iq], ks.ncg_p, nomx, numin, dUl, fr.omega, kqm.tetc)
    else:
        kcw = ft.bz_calcqdepw_par2(enk, ks.EF, kqm.kqid[:,iq], ks.ncg_p, nomx, numin, fr.omega, kqm.tetc)
        #kcw = ft.bz_calcqdepw_par(enk, ks.EF, kqm.kqid[:,iq], ks.ncg_p, nomx, numin, fr.omega, kqm.tetc)
        
    if False:
        save('enk',enk)
        save('kqid', kqm.kqid)
        save('tetc', kqm.tetc)
        rest = {'EF': ks.EF, 'ncg': ks.ncg_p, 'nomx': nomx, 'numin':numin}
        print rest
        import pickle
        f = open("rest.pkl","wb")
        pickle.dump(rest,f)
        f.close()
        save('kcw.'+str(iq), kcw)
        
    if False:
        fl = sorted(glob.glob('kcw.*'))
        if len(fl)>0:
            n = int(fl[-1].split('.')[-1])
        else:
            n=0
        fnm = 'kcw.'+str(n+1)
        fo = open(fnm, 'w')
        print >> fo, 'nk=', shape(kcw)[2], 'shape(kcw)=', shape(kcw), 'kcw='
        (nb1,nb2,nkp,nom) = shape(kcw)
        for ik in range(nkp):
            for ib in range(ks.ncg_p+nomx+1):
                for jb in range(nb2):
                    print >> fo, '%4d %4d %4d' % (ik+1, ib+1, jb+1), ('%17.13f'*nom) % tuple(kcw[ib,jb,ik,:].real)
        fo.close()
        sys.exit(0)
    
    return kcw

def svd_functions(iomega, wiomeg, om_min, om_max, om_nom, svd_cutoff, fout):
    print >> fout, '******* SVD of Matsubara Kernel for frequency dependence *******'
    print >> fout, 'svd_cutoff=', svd_cutoff, 'real axis mesh has min='+str(om_min)+'H max='+str(om_max)+'H and nom='+str(om_nom)
    # The real frequency mesh
    om, dom = mcmn.Give2TanMesh(om_min,om_max,om_nom)  
    # The imaginary frequency mesh for positive and negative frequency
    iom  = hstack( (-iomega[::-1], iomega[:]) )
    diom = hstack( ( wiomeg[::-1], wiomeg[:]) )
    # bosonic kernel for real part
    Ker = zeros((len(iom),len(om)))
    num = om*dom/pi
    om2 = om**2
    for i,iw in enumerate(iom):
        Ker[i,:] = num*sqrt(diom[i])/(om2+iom[i]**2)
    u, s, vh = linalg.svd(Ker, full_matrices=True)
    n = where(s < svd_cutoff)[0][0]
    print >> fout, 'singular eigenvalues kept='
    for i in range(n):
        print >> fout, "%2d %16.10f" % (i+1, s[i])
    Ul = transpose(u[:,:n])
    Ul *= 1/sqrt(diom)
    N2 = len(iom)/2
    dUl = (Ul[:,N2:]*diom[N2:]*2.0).T  # This is used to compute coefficients cl. Needs factor of two, because we integrate only over positive Matsubara frequencies
    return (Ul[:,N2:], dUl)
    

class G0W0:
    ddir = 'data'
    def __init__(self):
        if mrank==master and not os.path.isdir(self.ddir):
            os.mkdir(self.ddir)
            
    def PrintRadials(self, radf, strc, in1):
        isp=0
        lomaxp1 = shape(in1.nlo)[0]
        for iat in range(strc.nat):
            rx, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
            for l in range(in1.nt):
                frf = open(self.ddir+'/RadialFunctions.'+str(iat)+'.'+str(l), 'w')
                print >> frf, '#   u   udot',
                if l < lomaxp1:
                    for ilo in in1.nLO_at_ind[iat][l]:
                        print >>frf, '   ulo[ilo='+str(ilo)+']',
                print >> frf
                for ir,r in enumerate(rx):
                    print >> frf, r, radf.ul[isp,iat,l,ir], radf.udot[isp,iat,l,ir],
                    if l < lomaxp1:
                        for ilo in in1.nLO_at_ind[iat][l]:
                            print >> frf, radf.ulo[isp,iat,l,ilo,ir],
                    print >> frf
        
    def Compute(self, io, ps):
        (case, nspin, fout) = (io.case, io.nspin, io.out)
        
        #print >> fout, 'mem-usage[io]=', ps.memory_info().rss*b2MB,'MB'
        
        print >> fout, '-'*32, "Set Struct", '-'*32
        fout.flush()
        strc = w2k.Struct(case, fout)
        strc.debugprint(fout)
        
        print >> fout, '-'*32, "Lattice generation", '-'*32
        fout.flush()
        latgen = w2k.Latgen(strc, fout)

        
        print >> fout, '-'*32
        fr = FrequencyMesh(io.iopfreq, io.nomeg, io.omegmin, io.omegmax, io.iopMultiple, fout)
        print >> fout, '-'*32

        if io.iopfreq==5: # This means SVD is on
            om_min = io.omegmin if io.omegmin > 1e-3 else 1e-3
            om_max = io.omegmax if io.omegmax > 10 else 10
            om_nom = io.nomeg*3 if io.nomeg > 30 else 30*3
            (Ul, dUl) = svd_functions(fr.omega, fr.womeg, om_min, om_max, om_nom, io.svd_cutoff, fout)
            save(self.ddir+'/Ul', Ul)
            save(self.ddir+'/dUl', dUl)
        else:
            Ul, dUl = None, None
        
        #print >> fout, 'mem-usage[strc,latgen]=', ps.memory_info().rss*b2MB,'MB'
        
        print >> fout, '-'*32, "w2k_readin1", '-'*32
        fout.flush()
        in1 = w2k.In1File(case, strc, fout, io.lomax)
        # call set_lapwlo         !* set up LAPW+lo basis
        (Elapw, Elo) = w2k.get_linearization_energies(case, in1, strc, nspin, fout)
        in1.Add_Linearization_Energy(Elapw, Elo)
         
        #print >> fout, 'mem-usage[in1]=', ps.memory_info().rss*b2MB,'MB'
        
        Vr = w2k.Read_Radial_Potential(case, strc.nat, nspin, strc.nrpt, fout)
        radf = w2k.RadialFunctions(in1,strc,Elapw,Elo,Vr,nspin,fout)

        if mrank==master:
            self.PrintRadials(radf, strc, in1)
        
        del Vr
        radf.get_ABC(in1, strc, fout)

        #print >> fout, 'mem-usage[Vr,radf]=', ps.memory_info().rss*b2MB,'MB'
        
        print >> fout,'-'*32, 'w2k_readcore'
        fout.flush()
        wcore = w2k.CoreStates(case, strc, nspin, fout)
        
        #print >> fout, 'mem-usage[wcore]=', ps.memory_info().rss*b2MB,'MB'
        
        print >> fout, '-'*32, "set_mixbasis", '-'*32
        fout.flush()
        pb = ProductBasis( (io.kmr, io.pwm, io.lmbmax, io.wftol, io.lblmax, io.mb_emin, io.mb_emax), strc, in1, radf, wcore, nspin,  fout)
        pb.cmp_matrices(strc, in1, radf, wcore, nspin, fout)
        pb.generate_djmm(strc,latgen.trotij,in1.nt,fout)

        #print >> fout, 'mem-usage[ProductBasis]=', ps.memory_info().rss*b2MB,'MB'
        
        divs = io.nkdivs
        #divs = [15,15,15]
        kqm = KQmesh(divs, io.k0shift, strc, latgen, fout)
        kqm.tetra(latgen, strc, fout)

        #print >> fout, 'mem-usage[KQmesh]=', ps.memory_info().rss*b2MB,'MB'
        
        io_data={ 'emax_pol': io.emax_pol, 'emax_sc':  io.emax_sc, 'iop_core': io.iop_core, 'efermi': io.efermi, 'ibgw': io.ibgw, 'nbgw': io.nbgw, 'emingw': io.emingw, 'emaxgw': io.emaxgw}
        ks = KohnShamSystem(io_data, kqm, case, in1, strc, wcore, radf, nspin, fout)

        #print >> fout, 'mem-usage[KohnSham]=', ps.memory_info().rss*b2MB,'MB'
        
        Check_Equal_k_lists(ks.klist, kqm, fout)
        
        pw = PlaneWaves(ks.hsrws, io.kmr, io.pwm, case, strc, in1, latgen, kqm, False, fout)
        
        #print >> fout, 'mem-usage[PlaneWaves]=', ps.memory_info().rss*b2MB,'MB'
        
        #  determine minimum rmt
        self.rmtmin = min(strc.rmt)
        print >> fout, 'minimum Rmt=', self.rmtmin
        print >> fout, 'volume fraction of each MT-sphere is', pw.vmt
        fout.flush()
        
        ks.VectorFileRead(case, strc, latgen, kqm, pw, fout)

        del ks.klist
        
        ks.Vxc(case, in1, strc, radf, fout)

        #print >> fout, 'mem-usage[ks.VectorFileRead,ks.Vxc]=', ps.memory_info().rss*b2MB,'MB'
        
        
        #print >> fout, 'mem-usage[FrequencyMesh]=', ps.memory_info().rss*b2MB,'MB'
        
        kw = Kweights(io, ks, kqm, fout)
        
        #print >> fout, 'mem-usage[Kweights]=', ps.memory_info().rss*b2MB,'MB'
        
        me = MatrixElements2Band(io, pw, strc, latgen)
        Vxct = me.Vxc(strc, in1, latgen, kqm, ks, radf, pw, pb, fout)
        if mrank==master:
            save(self.ddir+'/Vxct', Vxct)
        
        fout.flush()
        
        del ks.uxcu
        del ks.Vxcs
        del ks.ksxc
        del ks.lxcm_f
        del ks.lmxc_f

        #prepare_Wannier(case, strc, in1, latgen, kqm, ks, radf, pw, io.rmax, fout)
        
        #print >> fout, 'mem-usage[MatrixElements2Bands,Vxct]=', ps.memory_info().rss*b2MB,'MB'
        
        sigx = zeros( (len(kqm.weight),ks.nbgw-ks.ibgw) )
        sigc = zeros( (len(kqm.weight),ks.nbgw-ks.ibgw,len(fr.omega)), dtype=complex )
        t_Coul, t_wgh, t_selfx, t_setev, t_calc_eps, t_calc_sfc = 0, 0, 0, 0, 0, 0

        #print >> fout, 'mem-usage[sigx,sigc]=', ps.memory_info().rss*b2MB,'MB'

        iqs,iqe,sendcounts,displacements = mcmn.mpiSplitArray(mrank, msize, len(kqm.qlist) )
        print >> fout, 'processor rank=', mrank, 'will do', range(iqs,iqe)
        fout.flush()
        #iqs,iqe = 0, len(kqm.qlist)
        #iqs,iqe = len(kqm.qlist)-1,len(kqm.qlist)
        for iq in range(iqs,iqe):
            t1 = timer()
            # Calculates matrix elements of the Coulomb interaction
            me.Coulomb(iq, io, strc, in1, latgen, kqm, ks, radf, pw, pb, fout)
            t2 = timer()
            #print >> fout, 'mem-usage[Coulomb(iq='+str(iq)+']=', ps.memory_info().rss*b2MB,'MB'
            #me.Coulomb_from_PW(iq, io, strc, in1, latgen, kqm, ks, radf, pw, pb, fout)

            # removes eigenvalues of the Coulomb interaction which are very small
            me.Coul_setev(iq, fout, 1e-7)
            t3 = timer()

            # calculates the exchange self-energy
            sx = me.calc_selfx(iq, strc, in1, latgen, kqm, ks, radf, pw, pb, wcore, kw, fout)
            sigx += sx
            t4 = timer()

            # removes more eigenvalues of the Coulomb interaction, which are smaller than io.barcevtol
            me.Coul_setev(iq, fout, io.barcevtol)
            t5 = timer()

            # Using tetrahedron method, computes polarization in the band basis (Lindhardt formula), but is integrated over momentum
            print >> fout, '***** iq=', iq
            kcw = Polarization_weights(iq, ks, kqm, wcore, fr, io.iop_bzintq, io.fflg, dUl, fout)

            t6 = timer()
            #print >> fout, 'mem-usage[sigx,Polarization_weights(iq='+str(iq)+']=', ps.memory_info().rss*b2MB,'MB'
            
            if iq==0:  # expansion around q=0 of V*P
                # Note : head_quantities = (head, mmatcv, mmatvv, mst)
                head_quantities = me.calc_head(strc, in1, latgen, kqm, ks, radf, pw, pb, wcore, kw, fr, kcw, io.iop_drude, io.eta_head, dUl, Ul, fout)
            else:
                head_quantities = None

            # Calculates V*(1/epsilon-1) = W-V
            (eps, epsw1, epsw2, head) = me.calc_eps(iq, head_quantities, strc, in1, latgen, kqm, ks, radf, pw, pb, wcore, kw, fr, kcw, self.ddir, dUl, Ul, fout)
            del kcw

            t7 = timer()
            #print >> fout, 'mem-usage[eps(iq='+str(iq)+')]=', ps.memory_info().rss*b2MB,'MB'
            # Calculates correlation self-energy
            sc = me.calc_selfc(iq, eps, epsw1, epsw2, head, strc, in1, latgen, kqm, ks, radf, pw, pb, wcore, kw, fr, self.ddir, dUl, Ul, fout)
            sigc += sc
            t8 = timer()
            #print >> fout, 'mem-usage[sigc(iq='+str(iq)+')]=', ps.memory_info().rss*b2MB,'MB'
            fout.flush()
            
            t_Coul += t2-t1
            t_setev += t3-t2+t5-t4
            t_wgh += t6-t5
            t_selfx += t4-t3
            t_calc_eps += t7-t6
            t_calc_sfc += t8-t7

        print >> fout, 'q-loop finished'
        fout.flush()
        
        if Parallel:
            sigx = comm.reduce(sigx, op=MPI.SUM, root=master)
            sigc = comm.reduce(sigc, op=MPI.SUM, root=master)
            
        if mrank==master:
            save(self.ddir+'/Sigmax', sigx)
            save(self.ddir+'/Sigmac', sigc)
            save(self.ddir+'/omega', fr.omega)
            save(self.ddir+'/womeg', fr.womeg)
            if fr.iopfreq == 4:
                save(self.ddir+'/omega_precise', fr.omega_precise)
                save(self.ddir+'/womeg_precise', fr.womeg_precise)
            
            #print >> fout, 'mem-usage[after q-look]=', ps.memory_info().rss*b2MB,'MB'
            
            for irk in range(len(kqm.weight)):
                for ie1 in range(ks.ibgw, ks.nbgw):
                    print >> fout, ' Sigx[irk=%3d,ie=%3d]=%16.12f' % (irk, ie1, sigx[irk,ie1-ks.ibgw].real)
            for irk in range(len(kqm.weight)):
                for ie1 in range(ks.ibgw, ks.nbgw):
                    for iom in range(len(fr.omega)):
                        print >> fout, ' Sigc[irk=%3d,ie=%3d,iom=%3d]=%16.12f%16.12f' % (irk, ie1, iom, sigc[irk,ie1-ks.ibgw,iom].real, sigc[irk,ie1-ks.ibgw,iom].imag)
            
            
            print >> fout, '## Coulomb:     t(Coulomb)         =%14.9f' % (t_Coul,)
            print >> fout, '## Coulomb:     t(setev)           =%14.9f' % (t_setev,)
            print >> fout, '## calc_selfx:  t(selfx)           =%14.9f' % (t_selfx,)
            print >> fout, '## eps weights: t(kcw)             =%14.9f' % (t_wgh,)
            print >> fout, '## eps calc_eps:t(calc_eps)        =%14.9f' % (t_calc_eps,)
            print >> fout, '## calc_selfc:  t(selfc)           =%14.9f' % (t_calc_sfc,)
            
            nbnd, nom = ks.nbgw-ks.ibgw, len(fr.omega)
            print >> fout, "AnalizeSpectra : calceqp"
            if io.iop_es >=0 :
                print >> fout, "# Parameters used:"
                print >> fout, "#  Analytic continuation (iop_ac) =", io.iop_ac
                print >> fout, "#  Fermi level shift    (iop_es)  =", io.iop_es
                print >> fout, "#  Nr.freq points  (nomeg)        =", nom
                print >> fout, "#  Number of AC poles (npar_ac/2) =", io.npar_ac/2

            EF_qp = 0
            isp=0
            bands = copy(ks.Ebnd[isp])
            save(self.ddir+'/KS_qp', bands[:,ks.ibgw:ks.nbgw]) 
            # quasiparticle energies for G0W0 scheme
            (eqp, eqp_im) = mcmn.Compute_quasiparticles(bands[:,ks.ibgw:ks.nbgw], bands[:,ks.ibgw:ks.nbgw], sigc, sigx, Vxct[:,ks.ibgw:ks.nbgw,ks.ibgw:ks.nbgw], fr.omega, (io.iop_ac,io.iop_es,io.iop_gw0,io.npar_ac,io.iop_rcf), isp, fout, PRINT=True)
            # the fermi energy for G0W0 scheme
            (EF, Eg, evbm, ecbm, eDos) = mcmn.calc_Fermi(eqp, kqm.atet, kqm.wtet, wcore.nval-ks.ibgw*2, io.nspin)
            print >> fout, ':E_FERMI_QP(eV)=  %12.4f' % (EF*H2eV,)
            if Eg > 0:
                print >> fout, ':BandGap_QP(eV)=  %12.4f' % (Eg*H2eV,)
            else:
                print >> fout, ':DOS_at_Fermi_QP= %12.4f' % (eDos,)
            print >> fout, 'Fermi: evbm=%12.4f  ecbm=%12.4f ' % (evbm*H2eV, ecbm*H2eV)
            save(self.ddir+'/GW_qp', eqp-EF)
            
            # First analyzing Kohn-Sham bands
            (nomax,numin) = mcmn.Band_Analys(bands[:,ks.ibgw:ks.nbgw], ks.EF, nbnd, 'KS', kqm, fout)
            # Next analyzing G0W0 bands
            (nomax,numin) = mcmn.Band_Analys(eqp, EF, nbnd, 'GW', kqm, fout)

if __name__ == '__main__':
    #pid = os.getpid()
    #ps = psutil.Process(pid)
    io = InOut("gw.inp", "pygw.out", mrank==master)
    #io = InOut("gw.inp", 'pygw.'+str(mrank), True)
    wgw = G0W0()
    #wgw.Compute(io, ps)
    wgw.Compute(io, None)
    if mrank==master:
        print 'PYGW DONE'

