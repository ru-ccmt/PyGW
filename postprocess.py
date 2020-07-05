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

from scipy import *
from scipy import linalg
from timeit import default_timer as timer
from scipy import optimize
import sys
import os

#from pylab import *
import gwienfile as w2k

from inout import *
from kqmesh import *
from mcommon import *

import os,psutil
#from pympler import asizeof

ddir = 'data'

class LightFrequencyMesh:
    def __init__(self, io, omega, womeg, fout):
        fginfo = ['Equally spaced mesh', 'Grid for Gauss-Laguerre quadrature', 'Grid for double Gauss-Legendre quadrature,', 'Grid of Tan-mesh for convolution', 'Using SVD basis and Tan-mesh for convolution']
        self.iopfreq = io.iopfreq
        self.omega, self.womeg = omega, womeg
        if io.iopfreq == 4:
            # another iopMultiple-times more precise mesh for self-energy integration
            minx = io.omegmin/(1 + 0.2*(io.iopMultiple-1.))
            om1, dom1 = Give2TanMesh(minx, io.omegmax, io.nomeg*io.iopMultiple)
            n = len(om1)/2
            self.omega_precise, self.womeg_precise = om1[n:], dom1[n:]
        print >> fout, 'Frequnecy grid for convolution: iopfreq='+str(io.iopfreq)+' om_max='+str(io.omegmax)+' om_min='+str(io.omegmin)+' nom='+str(io.nomeg)+': ', fginfo[io.iopfreq-1]
        for i in range(len(self.omega)):
            print >> fout, '%3d  x_i=%16.10f w_i=%16.10f' % (i+1, self.omega[i], self.womeg[i])

class SCGW0:
    def __init__(self, io):
        self.sigx  = load(ddir+'/Sigmax.npy')
        self.sigc  = load(ddir+'/Sigmac.npy')
        self.Vxct  = load(ddir+'/Vxct.npy')
        self.omega = load(ddir+'/omega.npy')
        self.womeg = load(ddir+'/womeg.npy')
        self.Ul, self.dUl = None, None
        self.iopfreq = io.iopfreq
        if io.iopfreq == 5:
            self.Ul = load(ddir+'/Ul' )
            self.dUl= load(ddir+'/dUl')
        self.fr = LightFrequencyMesh(io, self.omega, self.womeg, io.out)
        
    def ReadKSEnergy(self, case, nspin, core, io, fout):
        # Reading w2k energy files and its KS-eigenvalues
        spflag = ['up','dn'] if nspin==2 else ['']
        (self.klist, self.wegh, Ebnd, self.hsrws) = w2k.Read_energy_file(case+'.energy'+spflag[0], strc, fout)
        band_max = min(map(len,Ebnd))
        self.Ebnd = zeros( (nspin,len(Ebnd),band_max) )
        for ik in range(len(Ebnd)):
            self.Ebnd[0,ik,:] = Ebnd[ik][:band_max]
        if nspin==2:
            (self.klist, self.wegh, Ebnd, self.hsrws) = w2k.Read_energy_file(case+'.energy'+spflag[1], strc, fout)
            for ik in range(len(Ebnd)):
                self.Ebnd[1,ik,:] = Ebnd[ik][:band_max]
        # converting to Hartrees
        self.Ebnd *= Ry2H  # convert bands to Hartree
        #######
        # Recompute the Fermi energy, if neeeded
        if io['efermi'] >= 1e-2: # Recompute the Fermi energy
            (EF, Eg, evbm, ecbm, eDos) = calc_Fermi(self.Ebnd[0], kqm.atet, kqm.wtet, core.nval, nspin)
            print >> fout, 'Fermi energy was recomputed and set to ', EF*H2eV
        else:
            print >> fout, ' Use the Fermi energy from case.ingw'
            evbm = max( filter(lambda x: x < EF, self.Ebnd.flatten()) )
            ecbm = min( filter(lambda x: x > EF, self.Ebnd.flatten()) )
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
            nsp,nat,nc = shape(core.eig_core)
            for isp in range(nsp):
                for iat in range(nat):
                    core.eig_core[isp][iat][:] = array(core.eig_core[isp][iat][:])*Ry2H - EF
        self.EF = 0.0
        self.Eg = Eg
        print >> fout, 'Set EF to ', self.EF
        
        ncg = len(core.corind) # number of all core states
        if io['iop_core'] in [0,1]:
            self.ncg_c = ncg
        else:
            self.ncg_c = 0
        print >> fout, 'ncg_c=', self.ncg_c
        
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
        if io['ibgw'] < 0:
            nocc_at_k = [[len(filter(lambda x: x<io['emingw'], self.Ebnd[isp,ik,:])) for ik in range(nkp)] for isp in range(nspin)]# how many bands below io['emingw'] at each k-point
            self.ibgw = min(map(min,nocc_at_k))
            #print 'ibgw=', self.ibgw,  'nocc_at_k=', nocc_at_k
        else:
            self.ibgw = io['ibgw']
        if self.ibgw > self.nomax_numin[1]:
            print >> fout, 'KohnShamSystem: WARNING - range of gw bands!! ibgw=',self.ibgw,'numin=',self.nomax_numin[1]
            print >> fout, '*Now we will set ibgw to 0'
            self.ibgw = 0
        print >> fout, 'ibgw=', self.ibgw
        
        
    def Compute_selfc(self, bands, core, kqm, fout, PRINT=True):
        nirkp = len(kqm.weight)
        Nallkp = len(kqm.qlist)*nirkp
        
        ikas,ikae,sendcounts,displacements = mpiSplitArray(mrank, msize, Nallkp )
        print >> fout, 'processor rank=', mrank, 'will do', range(ikas,ikae)
        
        ik0 = ikas % nirkp
        iq0 = ikas / nirkp
        mwm = load(ddir+'/mwm.'+str(iq0)+'.'+str(ik0)+'.npy')
        (nom, nb1, nb2) = shape(mwm)
        
        t_read, t_cmp = 0.0, 0.0
        
        sigc = zeros( (nirkp, nb1, len(self.omega) ), dtype=complex )
        
        for i in range(ikas,ikae):
            irk = i % nirkp
            iq  = i / nirkp
            t1 = timer()
            mwm = load(ddir+'/mwm.'+str(iq)+'.'+str(irk)+'.npy')
            t2 = timer()
            sigc[irk,:,:] += Compute_selfc_inside(iq, irk, bands, mwm, self.fr, kqm, self.ncg_c, core, self.Ul, fout, PRINT)
            t3 = timer()
            t_read += t2-t1
            t_cmp  += t3-t2

        print >> fout, '## Compute_selfc : t_read    =%10.5f' % (t_read,)
        print >> fout, '## Compute_selfc : t_compute =%10.5f' % (t_cmp,)

        if Parallel:
            sigc = comm.allreduce(sigc, op=MPI.SUM)
        
        if PRINT:
            for irk in range(nirkp):
                for ie1 in range(nb1):
                    for iom in range(nom):
                        print >> fout, 'Sigc[irk=%3d,ie=%3d,iom=%3d]=%16.12f%16.12f' % (irk+1, ie1+1, iom+1, sigc[irk,ie1,iom].real, sigc[irk,ie1,iom].imag)                
        return sigc
    
    def calceqp(self, io, strc, kqm, nval, fout):
        (nirkp, nbnd, nom) = shape(self.sigc)
        print >> fout, "SCGW0 : calceqp"
        #anc_type=['old-fashioned Pade with n='+str(io.npar_ac-1),'modified Pade with '+str(io.npar_ac)+' coefficients','Simple quasiparticle approximation']
        print >> fout, "# Parameters used:"
        print >> fout, "#  Analytic continuation (iop_ac) =", io.iop_ac, 'i.e., '+io.anc_type[io.iop_ac]
        print >> fout, "#  Fermi level shift    (iop_es)  =", io.iop_es
        print >> fout, "#  Nr.freq points  (nomeg)        =", len(self.omega)
        print >> fout, "#  Number of AC poles (npar_ac/2) =", io.npar_ac/2
        isp=0
        
        EF_qp = 0
        bands = copy(self.Ebnd[isp])
        
        if not os.path.isfile(ddir+'/KSenk') and mrank==master:
            save(ddir+'/KSenk', bands[:,self.ibgw:nbnd+self.ibgw])
        
        # quasiparticle energies for G0W0 scheme
        nst,nend = self.ibgw,self.ibgw+nbnd
        (self.eqp0, eqp_im) = Compute_quasiparticles(bands[:,nst:nend], self.Ebnd[isp][:,nst:nend], self.sigc, self.sigx, self.Vxct[:,nst:nend,nst:nend], self.omega, (io.iop_ac,io.iop_es,io.iop_gw0,io.npar_ac,io.iop_rcf), isp, fout, PRINT=True)

        # the fermi energy for G0W0 scheme
        print >> fout, 'total nval=', core.nval, 'but with ibgw=', self.ibgw, 'the resulting nval=', core.nval-self.ibgw*2
        (EF, Eg, evbm, ecbm, eDos) = calc_Fermi(self.eqp0, kqm.atet, kqm.wtet, core.nval-self.ibgw*2, io.nspin)
        print >> fout, ':E_FERMI_QP(eV)=  %12.4f' % (EF*H2eV,)
        if io.iop_esgw0 == 1:
            self.eqp0 -= EF # shift bands so that EF=0
            EF_qp += EF # remember what shift we applied
            EF = 0      # now set EF to zero
        if Eg > 0:
            print >> fout, ':BandGap_QP(eV)=  %12.4f' % (Eg*H2eV,)
        else:
            print >> fout, ':DOS_at_Fermi_QP= %12.4f' % (eDos,)
        print >> fout, 'Fermi: evbm=%12.4f  ecbm=%12.4f ' % (evbm*H2eV, ecbm*H2eV)
        # First analyzing Kohn-Sham bands
        (nomax,numin) = Band_Analys(bands, self.EF, nbnd, 'KS', kqm, fout)
        # Next analyzing G0W0 bands
        (nomax,numin) = Band_Analys(self.eqp0, EF, nbnd, 'GW', kqm, fout)
        if False and mrank==master:
            save(ddir+'/GW_qp', self.eqp0)
            save(ddir+'/KS_qp', bands[:,nst:nend])

        eqp = copy(self.eqp0)
        #if (False):
        if (nomax >= numin): # metallic
            print >> fout, 'metallic bands, we will not consider GW0 scheme'
            return None
        else:                # insulating k-dependent gap
            if (nomax < numin): # insulating
                Egk = copy(bands[:,numin]-bands[:,nomax])
            else:
                Egk = copy(bands[:,numin])
            
            mix = io.mix_sc
            for isc in range(io.nmax_sc):
                bands[:,nst:nend] = copy(eqp[:,:])
                #bands[:,:nbnd] = bands[:,:nbnd]*(1-mix) + eqp[:,:]*mix
                
                if (nomax < numin): # insulating
                    Egk_new = copy(bands[:,numin]-bands[:,nomax])
                else:
                    Egk_new = copy(bands[:,numin])
                
                ediff = max(abs(Egk-Egk_new))
                Egk = Egk_new
                print >> fout, '#scgw: isc=', isc, 'ediff=', ediff, 'Egk=', (Egk*H2eV).tolist()
                #print  '#scgw: isc=', isc, 'ediff=', ediff, 'Egk=', Egk.tolist()
                io.out.flush()
                if ediff < io.eps_sc: break

                # Recompute correlation self-energy using quasiparticle's green's function
                sigc = self.Compute_selfc(bands, core, kqm, fout, False)
                
                # Compute the new quasiparticle energies
                (eqp, eqp_im) = Compute_quasiparticles(bands[:,nst:nend], self.Ebnd[isp][:,nst:nend], sigc, self.sigx, self.Vxct[:,nst:nend,nst:nend], self.omega, (io.iop_ac,-1,io.iop_gw0,io.npar_ac,io.iop_rcf), isp, fout, PRINT=False)
                
                # and recompute the Fermi energy on this quasiparticle bands
                (EF, Eg, evbm, ecbm, eDos) = calc_Fermi(eqp, kqm.atet, kqm.wtet, core.nval, io.nspin)
                print >> fout, ':E_FERMI_QP(eV)=  %12.4f' % (EF*H2eV,)
                if io.iop_esgw0 == 1:
                    eqp -= EF # shift bands so that EF=0
                    EF_qp += EF # remember what shift we applied
                    EF = 0      # now set EF to zero
                if Eg > 0:
                    print >> fout, ':BandGap_QP(eV)=  %12.4f' % (Eg*H2eV,)
                else:
                    print >> fout, ':DOS_at_Fermi_QP= %12.4f' % (eDos,)
                print >> fout, 'Fermi: evbm=%12.4f  ecbm=%12.4f ' % (evbm*H2eV, ecbm*H2eV)
                print >> fout, 'eferqp0=', EF_qp
                
            if ediff >= 5e-3:
                print >> fout, 'WARNING : GW0 did not converge. Will not analyze'
                return None
            else:
                Band_Analys(eqp, EF, nbnd, 'GW0', kqm, fout)
                return eqp
    
class QPs:
    def __init__(self, rbas, tizmat, rmax, nkp, fout):
        " Creates the star of the space group. I think it might not work in non-symorphic groups, becuase of the screw axis is not used"
        for itt in range(5):    
            rvec = zeros(3)
            for i in range(3):
                rvec[i] = linalg.norm(rbas[i,:])

            nr = array(map(int,1./rvec * rmax))*2
            self.rbas = rbas
            #print 'nr=', nr, 'rmax=', rmax
            rlen=[]
            rind=[]
            for ii,ir in enumerate(itertools.product(range(-nr[0],nr[0]+1),range(-nr[1],nr[1]+1),range(-nr[2],nr[2]+1))):
                rvec = dot(ir,rbas)
                rr = linalg.norm(rvec)
                if (rr <= rmax):
                    rlen.append(rr)
                    rind.append(ir)
            indx = argsort(rlen) # kind='stable')  # just obtaining index to the sorted sequence
            # rearange arrays so that they are sorted
            self.rlen = zeros(shape(rlen))           # |R|
            self.rind = zeros(shape(rind), dtype=int)           # \vR in cartesian
            for i0,i in enumerate(indx):
                self.rlen[i0]   = rlen[i]
                self.rind[i0,:] = rind[i]
            invrind = -ones((2*nr[0]+1,2*nr[1]+1,2*nr[2]+1),dtype=int16)
            for i,r in enumerate(self.rind):
                invrind[r[0]+nr[0],r[1]+nr[1],r[2]+nr[2]] = i
            
            self.slen=[0.0]  # length of vector Rm in for each star
            self.rst = zeros((len(self.rlen),2),dtype=int)  # contains stars of distance R
            self.rst[0,1] = len(tizmat)                # R=0 is nsym degenerate
            ist=0                                 # R>0 will start with ist=1
            for ippw,r in enumerate(self.rind):
                if ippw==0: continue
                #print 'ippw=', ippw, 'r=', r, 'rst['+str(ippw)+',0]=',self.rst[ippw,0]
                if self.rst[ippw,0]==0:              # this r did not occur yet, hence it belongs to a new star
                    ist += 1                    # index of the new star
                    self.slen.append(self.rlen[ippw])
                    self.rst[ippw,0] = ist           # remember which star this r corresponds to
                    #print 'ist=', ist, 'rst['+str(ippw)+',0]='+str(ist), 'r=', r
                    for isym in range(len(tizmat)):  # now we go over all group operations, and generate all members of this star
                        r_star = dot(tizmat[isym],r) # new member of the star
                        jppw = invrind[r_star[0]+nr[0],r_star[1]+nr[1],r_star[2]+nr[2]] # should exist, where?
                        if jppw >= 0:
                            #print 'ist=', ist, 'rst['+str(jppw)+',0]='+str(ist), ' r_star=', r_star
                            self.rst[jppw,0] = ist  # this is in the same star
                            self.rst[jppw,1] += 1   # and how many times the same vector appears in this star, i.e., degeneracy
            self.nst = ist+1
            print >> fout, 'Number of k-points=', nkp, 'Number of stars=', self.nst
            if self.nst > nkp*1.3 :
                break
            rmax = rmax*1.2
            print >> fout, 'Since the number of stars should be substantially bigger than number of k-points, I need to increase rmax. New rmax=', rmax
            
    def ReadGap_qp(self,fname):
        fi = open(fname, 'r')
        dat = fi.next().split()
        (ib0_kip,ib1_kip,nkp1,nsp_qp) = map(int,dat[:4])
        eferqp1 = float(dat[4])
        klist1 = zeros((nkp1,3),dtype=int)
        kvecs1 = zeros((nkp1,3))
        eks1 = zeros((nsp_qp,nkp1,ib1_kip-ib0_kip+1))
        eqp1 = zeros((nsp_qp,nkp1,ib1_kip-ib0_kip+1))
        for is1 in range(nsp_qp):
            nvbm = 0
            ncbm = ib1_kip
            for ik in range(nkp1):
                dat = map(int,fi.next().split())
                ikvec, idv = dat[2:5], dat[5]
                klist1[ik] = ikvec
                kvecs1[ik] = array(ikvec)/float(idv)
                io, iu = ib0_kip-1, ib1_kip-1
                for ib in range(ib0_kip-1,ib1_kip):
                    line = fi.next()
                    #print '"'+line[4:24]+'"', '"'+line[24:44]+'"'
                    ii, eks1[is1,ik,ib], eqp1[is1,ik,ib] = int(line[:4]), float(line[4:24]), float(line[24:44])
                    if eks1[is1,ik,ib] < 0:
                        if io<ib : io=ib
                    else:
                        if iu>ib : iu=ib
                    #print 'ik=', ik, 'ib=', ib, 'eks1=', eks1[is1,ik,ib], 'eqp1=', eqp1[is1,ik,ib]
                fi.next()
                if io > nvbm : nvbm = io
                if iu < ncbm : ncbm = iu
            #print 'nvbm,ncbm=', nvbm+1,ncbm+1
            if nvbm >= ncbm:
                #print 'nvbm >= ncbm', nvbm+1, ncbm+1
                #print '  nvbm = ncbm-1 is forced'
                nvbm = ncbm - 1
        fi.close()
        return (kvecs1,eks1,eqp1)
    
    def Pickett_Interpolate(self, kvecs1, eks1, eqp1, kvecs2, eks2, fout):
        """ This interpolation algorithm is described in PRB 38, 2721 (1988).
        """
        nkp1, nb1 = shape(eks1)
        nkp2, nb2 = shape(eks2)
        dek1 = eqp1-eks1 # this is what we interpolate, i.e., only the difference between QP-bands and KS bands
        
        print >> fout, 'number of stars=', self.nst, 'number of computed k-points=', nkp1
        
        den = float(self.rst[0,1]) # number of all group operations, i.e., den = nsymp = rst[0,1]
        smat1 = zeros((nkp1,self.nst), dtype=complex) # "star function" of the input bands, small k-point mesh
        for ir,r in enumerate(self.rind):
            ist, pref = self.rst[ir,:]
            smat1[:,ist] += exp(2*pi*dot(kvecs1,r)*1j) * pref/den # This is the "star function" Sm[ik,ist]
        
        smat2 = zeros((nkp2,self.nst), dtype=complex)   # "star function" on the dense mesh of k-points
        for ir,r in enumerate(self.rind):
            ist, pref = self.rst[ir,:]
            smat2[:,ist] += exp(2*pi*dot(kvecs2,r)*1j) * pref/den # The star function Sm[k,ist]

        c1,c2 = 0.25,0.25      # The two coefficients mentioned in the paper as C1 and C2
        rho = zeros(self.nst)  # We use here rho = 1 - 2*c1*R^2 + c1^2*R^4 + c2*R^6
        rho[0]=1.
        rmin = self.slen[1]
        for ist in range(self.nst):
            x2 = (self.slen[ist]/rmin)**2
            x6 = x2**3
            rho[ist] = (1-c1*x2)**2 + c2*x6
        
        # Now we start solving the equations in the paper
        sm2 = zeros((nkp1-1,self.nst),dtype=complex) # sm2 <- Sm[k_i]-Sm[k_n] in the paper
        nb = min(nb1,nb2)                            # number of input energies.
        dele = zeros((nkp1-1,nb))
        
        for ik in range(nkp1-1):
            sm2[ik,:]  = smat1[ik+1,:]  - smat1[0,:]       #  sm2[ik,istar] = Sm[k_i]-Sm[k_0]
            dele[ik,:] = dek1[ik+1,:nb] - dek1[0,:nb]      #  dele <- e[k_j]-e[k_0] in the paper
        
        h = zeros((nkp1-1,nkp1-1),dtype=complex)          # H_ij in the paper
        for ik in range(nkp1-1):
            for jk in range(nkp1-1):
                h[ik,jk] += sum(sm2[ik,:]*conj(sm2[jk,:])/rho[:])
        
        Hinv = linalg.inv(h)
        Hf = dot(conj(sm2.T), Hinv)
        for ist in range(self.nst):
            Hf[ist,:] *= 1/rho[ist]
        coef = dot( Hf, dele )
        coef[0,:] = dek1[0,:nb] - dot(smat1[0,:],coef) # epsilon_m in the paper
        
        # dek2[ik,nb] = smat2.T[ik,ist] * coef[ist,ib]
        dek2 = dot(smat2, coef)     # this is the resulting energy on the dense grid : e_m * S_m(k) in the paper
        eqp2 = eks2 + dek2.real     # finally, adding back the Kohn-Sham energy
        return eqp2

def toSmallerArray(ebnd2, nbs, nbe):
    eks2 = zeros((len(ebnd2),nbe-nbs))
    for ik in range(len(ebnd2)):
        eks2[ik,:] = array(ebnd2[ik][nbs:nbe])
    return eks2



def SaveBandPlot(filename, bands, klist2, knames):
    def Path_distance(klist):
        return cumsum( [0]+[linalg.norm(klist[ik+1,:]-klist[ik,:]) for ik in range(len(klist)-1)] )
    
    def PrintLegend(knames):
        leg = '{'
        for ik in range(len(knames)):
            name = knames[ik].strip()
            if name:
                leg += str(ik)+':"'+name+'",'
        leg += '}'
        return leg
    
    nk, nb = shape(bands)
    fgw = open(filename, 'w')
    # k-path distance
    xc = Path_distance( array(klist2) )
    print >> fgw, '# leg=' + PrintLegend(knames)
    for ik in range(len(klist2)):
        print >> fgw, '%10.6f  ' % (xc[ik],), ('%14.8f '*nb) % tuple(bands[ik,:]*H2eV)
    fgw.close()
    
if __name__ == '__main__':
    
    band_energy_file = ''
    if len(sys.argv)>1 and os.path.isfile(sys.argv[1]):
        band_energy_file = sys.argv[1]
        if len(sys.argv)<3:
            print 'When you give energy file for interpolation, you must also give as additional argument the Fermi energy for this bands'
            sys.exit(0)
        else:
            band_EFermi = float(sys.argv[2])*Ry2H

    io = InOut("gw.inp", "pypp.out", mrank==master)

    strc = w2k.Struct(io.case, io.out)
    latgen = w2k.Latgen(strc, io.out)
    latgen.Symoper(strc, io.out)
    
    kqm = KQmesh(io.nkdivs, io.k0shift, strc, latgen, io.out)
    kqm.tetra(latgen, io.out)
    core = w2k.CoreStates(io.case, strc, io.nspin, io.out)
    io_data={ 'emax_pol': io.emax_pol, 'emax_sc':  io.emax_sc, 'iop_core': io.iop_core, 'efermi': io.efermi, 'ibgw': io.ibgw, 'nbgw': io.nbgw, 'emingw': io.emingw, 'emaxgw': io.emaxgw, 'ibgw':io.ibgw, 'emingw':io.emingw, 'emaxgw':io.emaxgw}
    ansp = SCGW0(io)
    ansp.ReadKSEnergy(io.case, io.nspin, core, io_data, io.out)

    if band_energy_file and mrank==master: # Can produce the band plot
        eks1 = load(ddir+'/KS_qp.npy')
        nkp1 = shape(eks1)[0]

        (klist2, wegh, ebnd2, hsrws, knames) = w2k.Read_energy_file(band_energy_file, strc, io.out, give_kname=True)

        klist1 = kqm.kirlist/float(kqm.LCM)
        qps = QPs(latgen.rbas, latgen.tizmat, io.rmax, nkp1, io.out )
        
        nbs = ansp.ibgw # ths first band to use in interpolation
        nbe = min([len(ebnd2[ik]) for ik in range(len(ebnd2))]) # the last band to use
        nbe = min(shape(eks1)[1]+ansp.ibgw, nbe)
        
        #print 'nbs=', nbs, 'nbe=', nbe
        eks2 = toSmallerArray(ebnd2, nbs, nbe )*Ry2H - band_EFermi #ansp.EFermi
        
        kvecs2 = cart2int(klist2,strc,latgen)
        kvecs1 = cart2int(klist1,strc,latgen)
        
        qpfile = ddir+'/GW_qp.npy'
        if os.path.isfile(qpfile):
            eqp1 = load(ddir+'/GW_qp.npy')
            eqp2 = qps.Pickett_Interpolate(kvecs1, eks1, eqp1, kvecs2, eks2, io.out)
            SaveBandPlot(ddir+'/GW_bands.dat', eqp2, klist2, knames)
            SaveBandPlot(ddir+'/KS_bands.dat', eks2, klist2, knames)

    eqpn = ansp.calceqp(io, strc, kqm, core.nval, io.out)

    if band_energy_file and mrank==master: # Can produce band plot
        eqp02 = qps.Pickett_Interpolate(kvecs1, eks1, ansp.eqp0, kvecs2, eks2, io.out)
        SaveBandPlot(ddir+'/G0W0_bands.dat', eqp02, klist2, knames)
        
    if eqpn is not None and mrank==master:
        save(ddir+'/GW0_qp', eqpn)
        if band_energy_file: # Can produce band plot
            eqpn2 = qps.Pickett_Interpolate(kvecs1, eks1, eqpn, kvecs2, eks2, io.out)
            SaveBandPlot(ddir+'/GW0_bands.dat', eqpn2, klist2, knames)
