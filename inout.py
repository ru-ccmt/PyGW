from scipy import *
import re
import os

Ry2H = 0.5         # Rydberg2Hartree
H2eV = 27.2113961  # Hartree2eV
Ry2eV= 13.60569805 # 
FORT = True        # use fortran modules
sr2pi = sqrt(2*pi)
b2MB = 1./1024**2

import numpy
if "lcm" not in dir(numpy):
    def lcm(a, b):
        if a > b:
            greater = a
        else:
            greater = b
        while True:
            if greater % a == 0 and greater % b == 0:
                lcm = greater
                break
            greater += 1
        return lcm

def fermi(x):
    return 1/(exp(x)+1.)
def gauss(x, s):
    return 1/(s*sr2pi)*exp(-x**2/(2*s**2))
def _cfermi_(x):
    if x > 10.:
        return 0.0
    elif x< -10:
        return 1.0
    else:
        return fermi(x)
cfermi = vectorize(_cfermi_)
def _cgauss_(x,s):
    if abs(x) > 10*s:
        return 0.0
    else:
        return gauss(x,s)
cgauss = vectorize(_cgauss_, excluded=['s'])

class InOut:
    def  __init__(self, gwinp, gwout, PRINT=True):
        if PRINT:
            self.out = open(gwout,'w')
        else:
            self.out = open(os.devnull, 'w')
        self.float = "[-+]?\d*\.?\d*[e|E]?\d*"
        self.inpdata = open(gwinp,'r').readlines()
        self.ParsAll(gwinp)
        #del self.inpdata
        
    def Pars(self, pattern, default=''):
        for line in self.inpdata:
            m = re.search(pattern,line)
            if m is not None:
                return m.group(1)
        return default
    
    def ParsBlock(self, pattern):
        for ist,line in enumerate(self.inpdata):
            m = re.search(pattern,line)
            if m is not  None:
                break
        block_data=[]
        if ist < len(self.inpdata)-1:
            iend = ist+1
            while (iend < len(self.inpdata) ):
                m = re.search('\%', self.inpdata[iend])
                if m is not None:
                    break
                no_comments = ((self.inpdata[iend]).split('#'))[0]
                block_data.append( no_comments.split('|') )
                iend += 1
        return block_data
        
    def ParsAll(self, gwinp):
        
        print >>  self.out, '*'*80
        print >> self.out, '*',' '*28, 'GW Program : PyGW',' '*28, '*'
        print >> self.out, '*'*80
        print >> self.out
        print >> self.out, '-'*32, ' readingw ', '-'*32
        
        self.case = self.Pars("CaseName\s*= \"(\w*)\"")
        
        print >> self.out, 'CaseName=', self.case
        
        #self.iop_scratch = self.Pars("UseScratch\s*=\s*(\d)")
        #if(self.iop_scratch == '0'):
        #    print >> self.out, " -- Do not use scratch space"
        #elif (self.iop_scratch == '1'):
        #    print >> self.out, " -- Use different scratch files for different proc."
        #else:
        #    print >> self.out, " -- Use same scratch files for different proc."
        
        if (not os.path.exists(self.case+".struct")):
          print >> self.out, 'ERROR: structure file '+self.case+'.struct does not exist'
          sys.exit(0)
        
        #self.savdir = self.Pars("SavDir\s*=\s*\"(.*)\"")
        #print >> self.out, "SavDir: ", self.savdir
        
        #self.taskname = self.Pars("Task\s*=\s*([\w|']*)")
        #print >> self.out, "Task: ", self.taskname
        
        self.nspin = int( self.Pars("nspin\s*=\s*(\d)") )
        print >> self.out, 'nspin=', self.nspin
        fspin = 2./self.nspin
        if(self.nspin == 2):
            self.spflag = ['up','dn']
        else:
            self.spflag = ['']
        
        self.iop_core =  0
        print >> self.out, 'iop_core=', self.iop_core, ' -- using all core states in this calculation'
        self.lomax =  int( self.Pars("LOmax\s*=\s*(\d+)") )      # Number of local orbitals
        print >> self.out, 'LOmax=', self.lomax, '-- maximum for local orbitals, needed to be specified for wien2k vector/energy file'
        
        #self.lsymvector = self.Pars("SymVector\s*=\s*([T|F])")   #     whether using vector files taking symmetry into account
        #if self.lsymvector=='T':
        #    print >> self.out, "Use symmetrized eigenvector file"
        #else:
        #    print >> self.out, "Use non-symmetrized eigenvector file"
        
        # Set the window (in Ry) for the number of unoccupied bands for which GW band correction is to be calculated
        self.emaxgw = float( self.Pars("emaxgw\s*=\s*("+self.float+")",1e4) )
        self.emingw = float( self.Pars("emingw\s*=\s*("+self.float+")",-1e4) ) 
        print >> self.out, '(emingw,emaxgw)=('+str(self.emingw)+','+str(self.emaxgw)+') energy range in Ry (determines active bands) for self-energy calculation (only for external legs)'
        # Convert to Hartree unit
        self.emaxgw *= Ry2H
        self.emingw *= Ry2H

        self.mb_emin = float( self.Pars("MB_emin\s*=\s*("+self.float+")",-1e10) )
        self.mb_emax = float( self.Pars("MB_emax\s*=\s*("+self.float+")",20) )
        #self.mb_emin *= Ry2H we actually use Ry for linearization energies
        #self.mb_emax *= Ry2H
        
        # Block size used for Minm matrix operations
        #print >> self.out, "Options related to Minm:"
        #self.mblksiz = int( self.Pars("Minm_mblksiz\s*=\s*(\d+)",10) )
        #print >> self.out, "block size for m-index(mblksiz):", self.mblksiz
        
        # Read the energy cut-off for polarization matrix and correlation selfenergies 
        #  when emaxpol/emaxsc is negative, all unoccupied bands are used  
        #self.eminpol = float( self.Pars("eminpol\s*=\s*("+self.float+")",-1e10) )
        self.emax_pol = float( self.Pars("emaxpol\s*=\s*("+self.float+")",-1e10) )
        self.emax_sc  = float( self.Pars("emaxsc\s*=\s*("+self.float+")",-1e10) )
        #self.emin_sc  = float( self.Pars("eminsc\s*=\s*("+self.float+")",-1e10) )
        print >> self.out, 'emaxpol=', self.emax_pol, ' -- upper energy cutoff (Ry) for computing polarization',
        if self.emax_pol < 0 :
            print >> self.out, '-- We take all available bands in computing polarization'
        else:
            print >> self.out
        print >> self.out, 'emaxsc =', self.emax_sc, '-- upper cutoff energy (Ry) in computing dynamic self-energy in internal loop',
        if self.emax_sc < 0:
            print >> self.out, '-- We take all available bands in computing dynamic self-energy'
        else:
            print >> self.out
        
        #self.eminpol *= Ry2H # convert to unit of Ha. 
        self.emax_pol*= Ry2H # convert to unit of Ha. 
        self.emax_sc *= Ry2H # convert to unit of Ha. 
        #self.eminsc  *= Ry2H # convert to unit of Ha.
        
        #self.nvel = float( self.Pars("nvel\s*=\s*("+self.float+")") )
        #print >> self.out, 'nvel=',  self.nvel
        
        #self.core_ortho = self.Pars("Core_ortho",'F')
        #print >> self.out, 'core_ortho=', self.core_ortho
        
        self.lcmplx = self.Pars("ComplexVector\s*=\s*([T|F])",'F')
        print >> self.out, 'Complex or real KS vectors?',
        if self.lcmplx=='T':
            print >> self.out, " -- Complex Vector"
        else:
            print >> self.out, " -- Real Vector"
            
        # Parameters related to self-consistency
        self.eps_sc = float( self.Pars("eps_sc\s*=\s*("+self.float+")",1e-4) )
        self.nmax_sc = int(  self.Pars("nmax_sc\s*=\s*(\d+)",20) )
        self.mix_sc = float( self.Pars("mix_sc\s*=\s*("+self.float+")",1.0) )
        # iop_vxc : control how to treat vxc 
        #  0 -- calculate from vxc data read from the wien2k case.r2v file 
        #  1 -- directly read from an external file
        #self.iop_vxc = int( self.Pars("iop_vxc\s*=\s*(\d+)",0) )
        
        self.iop_bzint  = 0 
        self.iop_bzintq = 0
        
        self.eta_head = float( self.Pars("eta_head\s*=\s*("+self.float+")",0.01) )
        print >> self.out, 'eta_head=', self.eta_head, '-- broadening of plasmon in the head'
        #self.ztol_sorteq = 0.01
        #self.tol_taylor = 10.
        self.esmear = 0.01
        #self.eta_freq = 0.01
        #self.n_gauq = 8
        #self.ztol_vol =  1e-10
        self.iop_bcor = 0
        
        self.efermi = float( self.Pars("EFermi\s*=\s*("+self.float+")",1e4) )
        if self.efermi < 1e2:
            self.efermi *= Ry2H
            print >> self.out, "Read Fermi energy from the gw input:", self.efermi
        
        #self.iop_metallic = 0 # it seems we always expect insulator
        #self.spinmom = 0
        
        mbd = self.ParsBlock('^\%BZConv')
        #self.bzcon = mbd[0][0].strip().strip('"')
        self.fdep  = mbd[0][1].strip().strip('"')
        
        #print >> self.out, 'bzcon=', self.bzcon
        print >> self.out, 'Note: Using tetrahedron method on imaginary axis'
        if self.fdep == 'nofreq':
            self.fflg=1
        elif self.fdep == 'refreq':
            self.fflg=2
        elif self.fdep=='imfreq' or fdep=='IMFREQ':
            self.fflg=3
        else:
            print >> self.out, "WARNING: unsupported option for fdep!"
            print >> self.out, "--Taking default value: imfreq"
            self.fdep = 'imfreq'
            self.fflg=3
        #print >> self.out, "fdep=", self.fdep
        
        self.rmax=40.0

        mbd = self.ParsBlock('^\%FourIntp')
        if mbd:
            rmax = float(mbd[0][0])
        else:
            rmax = 40.0
            
        #mbd = self.ParsBlock('^\%kMeshIntp')
        #if mbd:
        #    self.iop_kip    = int(mbd[0][0])
        #    self.eqptag_kip = mbd[0][1].strip().strip('"')
        #else:
        #    self.iop_kip=0
        #    self.eqptag_kip=''
        #print >> self.out, 'iop_kip=', self.iop_kip, 'eqptag_kip=', self.eqptag_kip
        
        self.iop_drude = 1
        self.omega_plasma = -1.0
        self.omega_plasma /= H2eV
        
        self.iop_epsw = 0
        self.iop_mask_eps = 0
        self.q0_eps = 1/sqrt(3.0)
        
        mbd = self.ParsBlock('^\%FreqGrid')
        self.iopfreq, self.nomeg, self.omegmax, self.omegmin, self.iopMultiple = 4, 32, 20., 0.02, 1
        if mbd:
            self.iopfreq = int(mbd[0][0])
            self.nomeg   = int(mbd[0][1])
            self.omegmax = float(mbd[0][2])
            self.omegmin = 0
            if self.iopfreq == 1 or self.iopfreq>3:
                self.omegmin = float(mbd[0][3])
            self.svd_cutoff = 1e-10
            if self.iopfreq == 5:
                if len(mbd[0])>4:
                    self.svd_cutoff = float(mbd[0][4])
                if self.svd_cutoff > 1e-5:
                    print >> self.out, 'WARNING: svd_cutoff is large ==', self.svd_cutoff, 'This would be very imprecise calculation.'
                    print >> self.out, 'WARNING: Setting cutoff to 1e-10'
                    self.svd_cutoff = 1e-10
            if self.iopfreq == 4:
                if len(mbd[0])>4:
                    self.iopMultiple  = int(mbd[0][4])
                if self.iopMultiple < 1 or self.iopMultiple > 1000:
                    print >> self.out, 'WARNING: mesh for convolution should have iopMultiple more points, where specified iopMultiple='+str(self.iopMultiple)
                    self.iopMultiple = 1
                    print >> self.out, 'WARNING: this value of iopMultiple does not make sense, hence setting it to '+str(self.iopMultiple)
                
        print >> self.out, 'FreqGrid:'
        fginfo = ['Equally spaced mesh', 'Grid for Gauss-Laguerre quadrature', 'Grid for double Gauss-Legendre quadrature,', 'Grid of Tan-mesh for convolution', 'Using SVD basis and Tan-mesh for convolution']
        print >> self.out, '  iopfreq= %4d' % (self.iopfreq,), '-- '+ fginfo[self.iopfreq-1]
        print >> self.out, '  nomeg  = %4d' % (self.nomeg,)  , '-- number of Matsubara frequency points'
        print >> self.out, '  omegmax= '+str(self.omegmax), '-- upper frequency cutoff in Hartree'
        print >> self.out, '  omegmin= '+str(self.omegmin), '-- the low energy cutoff in Hartree'
        if self.iopfreq == 4:
            print >> self.out, '  iopMultiple='+str(self.iopMultiple), '-- how many more points should be used in frequency convolution'
        
        self.nproc_col, self.nproc_row = 0, 1
        
        mbd = self.ParsBlock('^\%kmesh')
        self.nkdivs=[1,1,1]
        self.k0shift=[0,0,0]
        if mbd:
            self.nkdivs = [int(mbd[0][0]), int(mbd[0][1]), int(mbd[0][2])]
            self.k0shift = [int(mbd[1][0]), int(mbd[1][1]), int(mbd[1][2])]
        print >> self.out, 'kmesh:'
        print >> self.out, '  nkdivs =', self.nkdivs, ' -- how many k-points along each reciprocal vector'
        print >> self.out, '  deltak =', self.k0shift, ' -- shift of k-mesh'

        mbd = self.ParsBlock('^\%SelfEnergy')
        self.fnpol_ac, self.iop_es, self.iop_ac = 2, 0, 1
        if mbd:
            self.npol_ac = int(mbd[0][0])
            self.iop_es  = int(mbd[0][1])
            self.iop_ac  = int(mbd[0][2])
        
        self.iop_esgw0 = 1 # whether shift the Fermi energy during self-consistent GW0
        #self.iop_gw0   = 1 # how the do GW0 self-consistent iteration
        self.iop_gw0 = int( self.Pars("iop_gw0\s*=\s*(\d)", 1) )
        
        self.iop_rcf   = 0.8 # real-frequency cutoff for using conventional pade
        self.npar_ac=2*self.npol_ac
        if self.nomeg < self.npar_ac:
            print >> self.out, 'WARNING: not enough freq for analytic continuation'
            print >> self.out, '  - npar_ac .gt. nomeg'
            print >> self.out, '  - npar_ac,nomeg=', self.npar_ac, self.nomeg
            print >> self.out, '  - reset npar_ac =nomeg'
            self.npar_ac = self.nomeg
            self.npol_ac = self.npar_ac/2
        
        self.anc_type=['old-fashioned Pade with n='+str(self.npar_ac-1),'modified Pade (Rojas, Godby and Needs) with '+str(self.npar_ac)+' coefficients','Simple quasiparticle approximation']
        print >> self.out, 'SelfEnergy: (analytic continuation of self-energy)'
        print >> self.out, '  npol_ac='+str(self.npol_ac), ', iop_es='+str(self.iop_es), ', iop_ac='+str(self.iop_ac)
        print >> self.out, '  iop_ac =', self.iop_ac, 'i.e., analytic continuation type is '+self.anc_type[self.iop_ac]
        print >> self.out, '  number of AC poles (npol_ac) =', self.npar_ac/2
        if self.iop_ac==0:
            print >> self.out, '  iop_rcf ='+str(self.iop_rcf), ' -- above this energy (in Ry) we use modified Pade (Rojas, Godby and Needs) and below this energy the original pade'
        #print >> self.out, '- Nr. of poles used in analytic continuation:', self.npol_ac
        #print >> self.out, '- Option for calculating selfenergy(iop_es): ', self.iop_es
        #if self.iop_es == 0:
        #    print >> self.out, "  -- perturbative calculation"
        #else:
        #    print >> self.out, "  -- iterative calculation"
        
        #print >> self.out, '- Option of analytic continuation (iop_ac):', self.iop_ac
        #if self.iop_ac == 1:
        #    print >> self.out, "  -- RGN method(Rojas, Godby and Needs)"
        #else:
        #    print >> self.out, "  -- Pade's approximation "
        
        mbd = self.ParsBlock('^\%MixBasis')
        if mbd:
            self.kmr = float(mbd[0][0])
            self.lmbmax, self.wftol, self.lblmax = int(mbd[1][0]), float(mbd[1][1]), int(mbd[1][2])
        else:
            self.kmr = 1.0
            self.lmbmax = 3
            self.lblmax = self.lmbmax*2
            self.wftol = 1e-4
        print >> self.out, 'Product (mixed) basis parameters'
        print >> self.out, '  kmr=', self.kmr, '-- Interstitial: Maximum |G| for plane waves. Note RKmax_{PB}=RKmax_{in1}*kmr'
        print >> self.out, '  lmbmax=', self.lmbmax, '-- MT-Spheres: Maximum l for products'
        print >> self.out, '  wftol=', self.wftol, '-- Linear dependence tolerance, i.e., disregard basis functions with smaller singular values'
        print >> self.out, '  MB_emin=', self.mb_emin, '-- low-energy cutoff (in linearization energy) for radial functions included in the product basis'
        print >> self.out, '  MB_emax=', self.mb_emax, '-- high-energy cutoff (in linearization energy) for radial functions included in the product basis'
        self.nspin_mb = 1
        self.ibgw = -1   # default stating index for GW calculation. -1 means it will be set later.
        self.nbgw = -1   # default last band index for GW calculation, -1 means it will be set later
        self.barcevtol = float( self.Pars("barcevtol\s*=\s*("+self.float+")",-1e-10) )
        #print >> self.out, 'barcevtol=', self.barcevtol
        
        self.lvorb = False
        
        #     Read the parameters for the Bare coulomb potential
        self.pwm, self.stctol =  2.0, 1e-15
        mbd = self.ParsBlock('^\%BareCoul')
        self.pwm    = float(mbd[0][0])
        self.stctol = float(mbd[0][1])
        
        # set the trancation radius for the bare Coulomb interaction, needed for finite systems
        self.iop_coulvm, self.iop_coul_x, self.iop_coul_c = 0, 0, 0
        self.rcut_coul = -1.0
        print >> self.out, 'Parameters for Coulomb matrix:'
        print >> self.out, "  pwm="+str(self.pwm), "-- Maximum |G| in kmr units "
        print >> self.out, "  stctol="+str(self.stctol)+ "-- Error tolerance for struc. const in Ewald summation"
        #print >> self.out, "  Coulomb interaction for exchange",    self.iop_coul_x 
        #print >> self.out, "  Coulomb interaction for correlation", self.iop_coul_c 
        #print >> self.out, '-'*55
        
        #mbd = self.ParsBlock('^\%gw')
        #self.iop_sxc, self.iop_vxc = 0, 0
        #if mbd:
        #    self.iop_sxc = int(mbd[0][0])
        #    self.iop_vxc = int(mbd[0][1])
        #print >> self.out, 'gw: sxc=%d vxc=%d' % (self.iop_sxc, self.iop_vxc)
        self.out.flush()
