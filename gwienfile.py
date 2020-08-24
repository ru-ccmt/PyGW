#!/usr/bin/env python
# @Copyright 2020 Kristjan Haule

'''
Classes to handle reading of WIEN2K files.
'''
import sys, re, os, operator
from numpy import zeros, arange
from scipy import *
from scipy import linalg
import rdVec

def vectorso_exists(case):
    filename = env.SCRATCH+"/"+case+".vectorso"
    return os.path.isfile(filename) and os.path.getsize(filename) > 0


def in1c_exists(case):
    filename = case+".in1c"
    return os.path.isfile(filename) and os.path.getsize(filename) > 0


class Struct:
    def __init__(self, case_file, fh_info=sys.stdout):
        self.case_file = case_file
        self.parse()
        self.w2kCorrect()

    def parse(self):
        '''The file is structured for reading by fortran code,
        so data is positioned by line and by column.'''
        f = open(self.case_file + '.struct', 'r')
        self.title = f.next().strip()

        line = f.next()
        self.lattice = line[0:4].strip()
        self.nat = int(line[27:30])
        self.latticename = line[30:].strip()

        self.mode = f.next()[13:17]

        line = f.next()
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma = [float(line[i:i+10]) for i in range(0,60,10)]

        self.iatom  = []
        self.mult   = []
        self.isplit = []
        self.aname  = []
        self.nrpt   = []
        self.r0     = []
        self.rmt    = []
        self.Znuc   = []
        self.pos    = []
        self.rotloc = []
        self.iatnr  = [] # iatnr[iat]>0 => cubic, iatnr[iat]<0 => non-cubic
        
        for iat in range(self.nat):
            line = f.next()
            self.iatnr.append( int(line[4:8]) )
            pos = [[float(line[col:col+10]) for col in 12+13*arange(3)]]
            
            line = f.next()
            mult = int(line[15:17])
            self.mult.append(mult)

            self.isplit.append( int(line[34:37]) )

            for mu in range(mult-1):
                line = f.next()
                pos.append( [float(line[col:col+10]) for col in 12+13*arange(3)] )

            self.pos.append(pos)
                
            line = f.next()
            self.aname.append( line[0:10].strip() )
            self.nrpt.append( int(line[15:20]) )
            self.r0.append( float(line[25:35]) )
            self.rmt.append( float(line[40:50]) )
            self.Znuc.append( int(float(line[55:65])) )
            
            #rmt[iat] = self.r0[iat]*exp(dh[iat]*(self.nrpt[iat]-1))               # calc. the muffin tin radius                                                                     

            rt=[]
            for i in range(3):
                line = f.next()
                rt.append( [float(line[col:col+10]) for col in 20+10*arange(3)] )
            self.rotloc.append(array(rt).T)

        self.aZnuc=[]
        for iat in range(self.nat):
            for mu in range(self.mult[iat]):
                self.aZnuc.append(self.Znuc[iat])

        line = f.next()
        dt = line.split()
        if len(line)<2 or dt[1][:6]!='NUMBER':
            f.close()
            return
        
        self.Nsym = int(dt[0])
        self.timat  = zeros((self.Nsym,3,3),dtype=int) # it is transpose compared to w2k fortran convention
        self.tau    = zeros((self.Nsym,3))
        
        for isym in range(self.Nsym):
            for j in range(3):
                line = f.next()
                #print int(line[0:2]),int(line[2:4]),int(line[4:6]),float(line[6:16])
                self.timat[isym,j,:] = [int(line[0:2]),int(line[2:4]),int(line[4:6])]
                self.tau[isym,j]     = float(line[6:16])
            #print self.timat[isym,:,:]
            #print self.tau[isym]
            ii = int(f.next())
            if (ii != isym+1):
                print 'WARNING : issues with reading symmetry operations in struct file at isym=', isym+1, 'ii=', ii
                f.close()
                return
        f.close()
        if self.mode == 'RELA':
            self.rel = True
        elif self.mode == 'NREL':
            self.rel = False
        else:
            print 'ERROR : struct file '+strc.case_file+'.struct seems to have wrong mode '+strc.mode
        #print >> fout, 'RELA=', strc.mode

        # create array of pos, which is list of lists, and not necessary equal length
        #n_eqm = max([len(self.pos[iat]) for iat in range(self.nat)])
        #self.vpos = zeros((self.nat,n_eqm,3))
        #for iat in range(len(self.pos)):
        #    ne = len(self.pos[iat])
        #    self.vpos[iat,:ne,:] = strc.pos[iat][:]
        self.mult = array(self.mult)
        self.rotloc = array(self.rotloc)
        ndf = sum(self.mult)
        self.vpos = zeros((3,ndf),order='F')
        idf = 0
        for iat in range(self.nat):
          for ieq in range(self.mult[iat]):
            self.vpos[:,idf] = self.pos[iat][ieq][:]
            idf += 1

        # We here rearange symmetry operations such that idenity appears first. This is just for the convenience
        if sum(abs(self.timat[0,:,:]-identity(3)))>0 or sum(abs(self.tau[0,:]))>0:
            #print  'Idenity did not appear first in the structure file. We will rearange symmetry operations so that identity is first'
            iisym = None
            Id = identity(3)
            for isym in range(len(self.timat)):
                if sum(abs(self.timat[isym,:,:]-Id))==0 and sum(abs(self.tau[0,:]))==0:
                  iisym = isym
                  break
            #print 'Idenity found at place', iisym
            if iisym is None:
                print 'ERROR : we could not find Idenity in the list of symmetries'
            else:
                self.timat[0,:,:], self.timat[iisym,:,:] = copy(self.timat[iisym,:,:]), copy(self.timat[0,:,:])
                self.tau[0,:], self.tau[iisym,:] = copy(self.tau[iisym,:]), copy(self.tau[0,:])
            #print self.timat[0,:,:], self.timat[iisym,:,:]

            
    def w2kCorrect(self):
      if self.alpha == 0:
        self.alpha = 90.0
      if self.beta  == 0:
        self.beta = 90.0
      if self.gamma == 0 :
        self.gamma = 90.0

      if self.lattice.strip() == 'H':
          self.gamma = 120.
          

          
    def debugprint(self, f=sys.stdout):
        print >> f, '***** Structure File Contents *****'
        print >> f, 'title =', self.title
        print >> f, 'Number of sorts, nat =', self.nat
        print >> f, 'unit cell parameters (a,b,c) =', self.a, self.b, self.c
        print >> f, 'angles (alpha, beta, gamma) =', self.alpha, self.beta, self.gamma


        printlist = [
            (self.aname,  'Atom name, aname ='),
            (self.iatnr,  'Atom index, cubic/noncubic environment ='),
            (self.mult,   'Number of atoms of this type, mult ='),
            (self.isplit, 'Symmetry of the atom, isplit ='),
            (self.nrpt,   'Number of radial points, npt ='),
            (self.r0,     'First point in radial mesh, r0 ='),
            (self.rmt,    'Muffin tin radius, rmt ='),
            (self.Znuc,   'Atomic number, Znuc ='),
            ]

        for i in range(self.nat):
            print >> f, '---------- atom type', i, '------------'

            for var,docstr in printlist:
                print >> f, docstr, var[i]

            print >> f, 'Position(s) inside the unit cell:'
            for m in range(self.mult[i]):
                print >> f, '    pos =', self.pos[i][m]

            print >> f, 'Local rotation matrix:'
            for row in self.rotloc[i,:,:]:
                print >> f, row

        #print >> f, '***********************************'

    def flat(self, notflat):
        '''Return a flat view of given data as a list.
        Example: if w.mult = [2,4] and w.aname = ['V', 'O']
        w.flatten(w.aname) -> ['V', 'V', 'O', 'O', 'O', 'O']'''
        if notflat is self.pos:
            listoflists = self.pos
        else:
            listoflists = [[elem]*mult for elem,mult in zip(notflat, self.mult)]
        return reduce(operator.add, listoflists)

    def radial_mesh(self, iat):
        npt = self.nrpt[iat]
        dh  = log(self.rmt[iat]/self.r0[iat])/(npt - 1)      # logarithmic step for the radial mesh
        dd = exp(dh)
        return( self.r0[iat]*dd**range(npt), dh, npt )

    
class Latgen:
    """
      Given wien2k structure file creates lattice BR2 and gbas, rbas, Volume
    """
    def __init__(self, strc, fout):
        self.pia   = array([2.0*pi/strc.a, 2.0*pi/strc.b, 2.0*pi/strc.c])
        self.alpha = [strc.alpha*pi/180., strc.beta*pi/180., strc.gamma*pi/180.]
        self.ortho = False
        
        self.br2 = zeros((3,3))
        if strc.lattice[:1]=='H': # hexagonal
          print >> fout, 'hexagonal lattice'
          self.br2[0,0] = 2.0/sqrt(3.0)
          self.br2[0,1] = 1.0/sqrt(3.0)
          self.br2[1,1] = 1.0
          self.br2[2,2] = 1.0
          self.rvfac = 2.0/sqrt(3.0)
          self.ortho = False
          for j in range(3):
            self.br2[j,:] *= self.pia[j]
        elif strc.lattice[:1] in ['S','P']: # primitive or simple
          print >> fout, 'primitive or simple lattice'
          self.ortho = True
          for i in range(3):
            if( abs(self.alpha[i]-pi/2.) > 0.0001):
              print >> fout, 'alpha['+str(i)+'] not equal 90'
              self.ortho = False
              # includes triclinic, monoclinic, simple orthorhombic, tetragonal and cubic
          sinbc = sin(self.alpha[0])
          cosab = cos(self.alpha[2])
          cosac = cos(self.alpha[1])
          cosbc = cos(self.alpha[0])
          wurzel=sqrt(sinbc**2-cosac**2-cosab**2+2*cosbc*cosac*cosab)
          self.br2[0,0]= sinbc/wurzel*self.pia[0]
          self.br2[0,1]= (-cosab+cosbc*cosac)/(sinbc*wurzel)*self.pia[1]
          self.br2[0,2]= ( cosbc*cosab-cosac)/(sinbc*wurzel)*self.pia[2]
          self.br2[1,1]=  self.pia[1]/sinbc
          self.br2[1,2]= -self.pia[2]*cosbc/sinbc
          self.br2[2,2]=  self.pia[2]
          self.rvfac = 1.0/wurzel
        elif strc.lattice[:1] == 'F': # face centered
          print >> fout, 'face centered lattice'
          self.br2[0,0] = -1.0
          self.br2[1,0] =  1.0
          self.br2[2,0] =  1.0
          self.br2[0,1] =  1.0
          self.br2[1,1] = -1.0
          self.br2[2,1] =  1.0
          self.br2[0,2] =  1.0
          self.br2[1,2] =  1.0
          self.br2[2,2] = -1.0
          self.rvfac = 4.0
          self.ortho = True
          for j in range(3):
            self.br2[j,:] *= self.pia[j]
        elif strc.lattice[:1] == 'B': # body centered
          print >> fout, 'body centered lattice'
          self.br2[0,0] = 0.0
          self.br2[1,0] = 1.0
          self.br2[2,0] = 1.0
          self.br2[0,1] = 1.0
          self.br2[1,1] = 0.0
          self.br2[2,1] = 1.0
          self.br2[0,2] = 1.0
          self.br2[1,2] = 1.0
          self.br2[2,2] = 0.0
          self.rvfac = 2.0
          self.ortho = True
          for j in range(3):
            self.br2[j,:] *= self.pia[j]
        elif strc.lattice[:1] == 'R': # rhombohedral
          print >> fout, 'rhombohedral lattice'
          self.br2[0,0] =  1.0/sqrt(3.0)
          self.br2[0,1] =  1.0/sqrt(3.0)
          self.br2[0,2] = -2.0/sqrt(3.0)
          self.br2[1,0] = -1.0
          self.br2[1,1] =  1.0
          self.br2[2,0] =  1.0
          self.br2[2,1] =  1.0
          self.br2[2,2] =  1.0
          self.rvfac = 6.0/sqrt(3.0)
          self.ortho = False
          for j in range(3):
            self.br2[j,:] *= self.pia[j]
        elif strc.lattice[:1] == 'C': # base centered
          print >> fout, 'base centered lattice'
          if strc.lattice[1:3] == 'XZ':
            ix=0
            iy=1
            iz=2
          elif strc.lattice[1:3] == 'YZ':
            ix=1
            iy=0
            iz=2
          elif strc.lattice[1:3] == 'XY':
            ix=0
            iy=2
            iz=1
          if( abs(self.alpha[iz]-pi/2.0) < 0.0001 ):
            #   orthorombic case
            self.br2[ix,ix] =  1.0
            self.br2[ix,iz] =  1.0
            self.br2[iy,iy] =  1.0
            self.br2[iz,ix] = -1.0
            self.br2[iz,iz] =  1.0
            self.rvfac = 2.0
            self.ortho = True
            for j in range(3):
              self.br2[j,:] *= self.pia[j]
          else:
            #  monoclinic case
            print 'alpha['+str(iz)+'] not equal 90 degrees'
            sinab = sin(self.alpha[iz])
            cosab = cos(self.alpha[iz])
            self.br2[ix,ix] =  self.pia[ix]/sinab
            self.br2[ix,iy] = -self.pia[iy]*cosab/sinab
            self.br2[ix,iz] =  self.pia[ix]/sinab
            self.br2[iy,iy] =  self.pia[iy]
            self.br2[iz,ix] = -self.pia[iz]
            self.br2[iz,iz] =  self.pia[iz]
            self.rvfac = 2.0/sinab
            self.ortho= False
        else:
          print >> fout, 'ERROR wrong lattice=', strc.lattice
          sys.exit(1)
        
        #  define inverse of cellvolume
        vi = self.rvfac/ (strc.a * strc.b * strc.c)
        self.Vol = 1./vi
        # Calculate the basis vectors of the real lattice
        self.gbas = self.br2[:,:] / (2.0*pi)
        self.rbas = linalg.inv(self.gbas)
        for i in range(3):
            for j in range(3):
                if abs(self.rbas[i,j])<1e-14:
                    self.rbas[i,j]=0
        print >> fout, 'Ortho=', self.ortho
        print >> fout, 'BR2=', self.br2
        print >> fout, 'Unit cell volume=', self.Vol
        print >> fout, 'gbas=', self.gbas
        print >> fout, 'rbas=', self.rbas
        self.Rotdef( strc, fout)
        self.Symoper( strc, fout)
        
    def Rotdef(self, strc, fout):
        """   define rotation matrices if required
        Finds the symmetry operations, which transform equivalent
        atoms into each other
        the operation rotij[idf,3,3], tauij[idf,3] transforms
        a position of a "not equivalent" atom to the position of
        a corresponding equivalent atom (idf).
        """
        Nat_all = sum([len(strc.pos[i]) for i in range(len(strc.pos))])
        #print 'Nat_all=', Nat_all
        
        self.trotij = zeros((Nat_all,3,3))
        self.tauij  = zeros((Nat_all,3))
        
        nnat = 0
        for iat in range(strc.nat):
            poss = strc.pos[iat]
            pos0 = poss[0]
            for ieq in range(len(poss)):
                cpos = poss[ieq]
                #print >> fout, iat, ieq, cpos
                Found_sym = False
                for isym in range(strc.Nsym):
                    sym_pos = dot(strc.timat[isym,:,:],pos0)
                    sym_pos += strc.tau[isym,:]
                    for i in range(3):
                        sym_pos[i] = (sym_pos[i] % 1.0)
                    delta_pos = abs(sym_pos[:]-cpos[:])
                    #print >> fout, iat, ieq, isym, 'delta_pos=', delta_pos
                    if sum(delta_pos) < 3e-4:
                        Found_sym = True
                        break
                    # Not done yet, some lattices require half translation
                    # Check now for centered lattices
                    delta_cpos = array([(delta_pos[i]+0.5) % 1.0 for i in range(3)])
                    if strc.lattice[:1]=='B': #  bc
                        if sum(abs(delta_cpos)) < 3e-4:
                            Found_sym = True
                            break
                    elif strc.lattice[:1]=='F': #  fc
                        delta = delta_cpos[:]
                        for i in range(3):
                            delta[i] = delta_pos[i]
                            if sum(abs(delta)) < 3e-4:
                                Found_sym = True
                                break
                        if Found_sym:
                            break
                    elif strc.lattice[:1]=='C':
                        if strc.lattice[1:3]=='XY':
                            delta = [delta_cpos[0], delta_cpos[1], delta_pos[2] ]
                        elif strc.lattice[1:3]=='XZ':
                            delta = [delta_cpos[0], delta_pos[1], delta_cpos[2] ]
                        elif strc.lattice[1:3]=='YZ':
                            delta = [delta_pos[0], delta_cpos[1], delta_cpos[2] ]
                        if sum(abs(delta)) < 3e-4:
                            Found_sym = True
                            break
                if Found_sym:        
                    self.trotij[nnat,:,:] = strc.timat[isym,:,:]
                    self.tauij[nnat,:] = strc.tau[isym,:]
                else:
                    print >> fout, "ERROR in rotdef: no sym operation found for iat=",iat,"ieq=",ieq,"idf",nnat
                    print "ERROR in rotdef: no sym operation found for iat=",iat,"ieq=",ieq,"idf",nnat
                    sys.exit(1)
                # correct_rotij
                #  Redefines local rotation matix from struct-file with
                #  unitary transformation  u(-1) * S * u  for non-orthogonal lattice
                if not self.ortho and strc.lattice[1:3]!='CXZ':
                    # Note  in dmft1 and lapw2 we used : rotij_cartesian = BR1 * rotij * BR1inv^{-1}
                    #  but here we use the convention as in lapw1 : rotij_cartesian = BR2 * rotij * ((BR2^T)^{-1})^T
                    #  Note also that in this way we rotate real space vectors, not momentum space vectors. The latter are rotated with inverse of BR2 matrices (see Symoper)
                    rtij = transpose( dot( dot(self.gbas,self.trotij[nnat,:,:].T), self.rbas ) )
                    self.trotij[nnat,:,:] = rtij
                    
                print >> fout, 'For atom '+strc.aname[iat]+' ieq='+str(ieq)+' at pos=',cpos, 'found tauij=', self.tauij[nnat,:], 'rotij.T=', (self.trotij[nnat,:,:]).tolist()
                nnat += 1
        
    def Symoper(self,  strc, fout):
        """ Redefine symmetry matrices such that they will work in  lattice system
            rather than cartesian system.
            Example: bcc structure. We will use 1BZ k-point in lattice vectors written as k_l=(i/N, j/N, k/N).
                     For bcc structure the same k-point in cartesian coordinates should be k_c=( (i+j)/N, (i+k)/N, (j+k)/N )
                     For bcc structure BR2 = [[0,1,1],[1,0,1],[1,1,0]], hence
                     BR2*k_l = k_c
                     We could apply symmetry operations on k_c as R*k_c. Alternatively,  we can use
                     (BR2)^{-1}*R*BR2*k_l = rbas*R*gbas
        """
        if self.ortho:  # I suspect that here we should have "or strc.lattic[:3]=='CXZ'", hence I suspect that this does not work for CXZ!
            # If the system is orthogonal, than we can work with integer k-vectors even in lattice coordinate system
            # and we do not need to work in cartesian. 
            self.tizmat = zeros(shape(strc.timat),dtype=int)  # this is like timat matrix, but tizmat is used in integer representation, while timat has to be used only in semi-cartesian coordinates.
            self.iztau = zeros(shape(strc.tau),dtype=float)
            for isym in range(strc.Nsym):
                ttmp = round_(transpose( dot( dot(self.rbas, strc.timat[isym,:,:].T), self.gbas ) ))
                self.tizmat[isym,:,:] = ttmp[:,:]
                self.iztau[isym,:] = dot(self.rbas,strc.tau[isym,:])
                #  Note that this is very different than in dmft1 or dmft2
        else:
            # we are forced to stary in cartesian system
            #print 'WARNING : Check Symoper if it really works. I think it does not'
            self.tizmat = copy(strc.timat)
        
class In1File:
    def __init__(self, case_file, strc, fout, lomax=-1):
        if in1c_exists(case_file):
            fin1 = open(case_file+'.in1c','r')
        else:
            fin1 = open(case_file+'.in1', 'r')
        lines = fin1.readlines()

        if (lomax<0):
            # In origonal wien2k lomax is seto to 3, but the GW part changes that to lomax=10
            # Here we determine from in1 file what might lomax be, but be careful because in initializing GW we usualy set lomax->10
            # This has consequence for reading case.energy file, to find energy for radial functions. They are stored in backes of lomax.
            # first pass through the file to determine lomax
            lomax = 3
            iread = 2
            for iat in range(strc.nat):
                line_dat = lines[iread].split()
                iread += 1
                Ei0, nlr, iapw = float(line_dat[0]), int(line_dat[1]), int(line_dat[2])  # default parameters ei,iapw and number of exceptions nlr
                dat=[]
                for j in range(nlr):
                    line_dat = lines[iread].split()
                    iread += 1
                    dat.append( [int(line_dat[0]), float(line_dat[1]), float(line_dat[2]), line_dat[3], int(line_dat[4])] )
                for j in range(nlr):
                    l, iapw = dat[j][0], dat[j][4]
                    l_prev = dat[j-1][0] if j>0 else -1
                    if l==l_prev or iapw==1: # LO or lo
                        if l >= lomax:
                            lomax = l+1 
        print >> fout, 'lomax set to', lomax, ' make sure that this is  compatible with wien2k case.energy file'
        
        Ef = 0.5
        str_float = "[-+]?\d*\.?\d*[e|E]?\d*"
        m = re.search('EF=('+str_float+')',lines[0])
        if m is not None:
            Ef = float(m.group(1))
        line_dat=lines[1].split()
        self.rkmax, self.lmax, self.lnsmax = float(line_dat[0]), int(line_dat[1]), int(line_dat[2])
        self.nt = self.lmax + 1
        if self.lnsmax==0: self.lnsmax = 2
        #kmax = rkmax/rmtmin
        print >> fout, 'rkmax=', self.rkmax, 'lmax=', self.lmax, 'lnsmax=', self.lnsmax, 'nat=', strc.nat
        self.lapw = ones((self.lmax+1,strc.nat), dtype=int, order='F')
        self.nlo  = zeros((lomax+1,strc.nat), dtype=int, order='F')
        self.loor  = ones((lomax+1,strc.nat), dtype=int, order='F')
        self.nlo_tot = 0
        iread = 2
        for iat in range(strc.nat):
            line_dat = lines[iread].split()
            iread += 1
            Ei0, nlr, iapw = float(line_dat[0]), int(line_dat[1]), int(line_dat[2])  # default parameters ei,iapw and number of exceptions nlr
            if abs(Ei0-0.3)<1e-6: Ei0=Ef-0.2
            if (iapw == 1):
                self.lapw[:,iat] = 0                                     # default value for lapw(l,iat) = False
                print >> fout, 'atom=', strc.aname[iat], 'E=', Ei0, ' apw with n.cases=', nlr
            else:
                self.lapw[:,iat] = 1
                print >> fout, 'atom=', strc.aname[iat], 'E=', Ei0, 'lapw with n.cases=', nlr

            self.loor[:,iat] = 0
            self.nlo[:,iat] = 0
            # Read  exceptional cases. 
            #   1) If the same l appears twice, the second one is a LO
            #         set loor to true, nlo must be incremented by 1
            #         and g_tot by (2l+1)*mult[iat]
            #   2) If it is an APW+lo (iapw == 1): lapw(l,iat) is set to false, and
            #         g_tot has to be incremented, nlo =1 , elo = ei...
            dat=[]
            for j in range(nlr):
                line_dat = lines[iread].split()
                iread += 1
                dat.append( [int(line_dat[0]), float(line_dat[1]), float(line_dat[2]), line_dat[3], int(line_dat[4])] )
            for j in range(nlr):
                l, Ei, dE, scanflag, iapw = dat[j][0], dat[j][1], dat[j][2], dat[j][3], dat[j][4]
                l_prev = dat[j-1][0] if j>0 else -1
                this_is_LO  = (l==l_prev) # current must be an LO
                if this_is_LO: 
                    self.loor[l,iat] = 1
                    self.nlo[l,iat] += 1  # this is ilo in fortran
                    self.nlo_tot += (2*l+1)*strc.mult[iat]
                    print >> fout, (' '*10+'l=%2d    %15s') % (l, 'Local Orbital')
                else:
                    if (iapw==1):
                        self.lapw[l,iat]=0
                        self.nlo[l,iat] =1   # this is ilo in fortran
                        self.nlo_tot += (2*l+1)*strc.mult[iat]
                        print >> fout, (' '*10+'l=%2d    %15s') % (l, 'APW+lo')
                    else:
                        print >> fout, (' '*10+'l=%2d    %15s') % (l, 'LAPW')
                    
        
        nLO_at = zeros((strc.nat,lomax+1), dtype=int)
        self.nLO_at_ind = []
        for iat in range(strc.nat):
            _nlo_at_ind_=[]
            for l in range(lomax+1):
                if self.lapw[l,iat]:  # LAPW+LO, hence all orbitals are actually local orbitals
                    nLO_at[iat,l] = self.nlo[l,iat]
                    _nlo_at_ind_.append( range(1,self.nlo[l,iat]+1) )
                else:                 # this is APW+lo, hence one is not LO, but corresponds to lo in APW+lo
                    nLO_at[iat,l] = self.nlo[l,iat]-1
                    _nlo_at_ind_.append( range(1,self.nlo[l,iat]) )# the first orbital is not LO
            self.nLO_at_ind.append(_nlo_at_ind_)

        self.nLO_at = zeros((lomax+1,strc.nat), dtype=int, order='F')
        self.nLO_at[:,:] = copy(nLO_at.T)
        print >> fout, 'nLO_at_ind=', self.nLO_at_ind
        #print 'lapw=', self.lapw
        #print 'nlo=', self.nlo
        print >> fout, 'nLO_at=', nLO_at, '==', [[ len(self.nLO_at_ind[iat][l]) for l in range(lomax+1)] for iat in range(strc.nat)]
        self.nlomax = amax(nLO_at)

        
        self.l_newlo=0
        if self.nlomax>1:
            self.l_newlo = 1
            nloat = amax(self.nlo)
            print >> fout, "More than one LO's are detected: set l_newlo to", nloat

        print >> fout, "Max. number of LO's per l per atom  ( nlomax):", self.nlomax

            
        print >> fout, "Information on Local Orbitals (LO or lo)  "
        print >> fout, "  lomax=",lomax
        for l in range(lomax):
            print >> fout, "l=%2d LO? #lo #LO        " % (l,),
        print >> fout
        for iat in range(strc.nat):
            for l in range(lomax):
              print >> fout, "    %4d%4d%4d        " % (self.loor[l,iat],self.nlo[l,iat],self.nLO_at[l,iat]),
            print >> fout
        print >> fout, ' Total number of local orbitals (nlo_tot)= ', self.nlo_tot

    def Add_Linearization_Energy(self, Elapw, Elo):
        self.Elapw = Elapw
        self.Elo = Elo
        #El = zeros((lmax+1,strc.nat))
        #Elo = zeros((lomax+1,nloat,strc.nat))
        #umt    = zeros((6,lmax+1,strc.nat)) # collective storage of e,p,pe,dp,dpe,pei 
                                            # umt(0,:,:) -->   e(:,:)  -- linearization energy E_l for atom iat
                                            # umt(1,:,:) -->   p(:,:)  -- u_l(r,E)  at r = RMT(iat)
                                            # umt(2,:,:) -->   pe(:,:) -- ud(r,E) at r = RMT(iat)
                                            # umt(3,:,:) -->   dp(:,:) -- derivative of u_l(r,E) to r at r=RMT(iat)
                                            # umt(4,:,:) -->   dpe(:,:) --derivative of ud_l(r,E) to r at r=RMT(iat) 
                                            # umt(5,:,:) -->   pei(:,:) -- norm of ud_l(r,E) integrated over MT sphere 


def get_linearization_energies(case, in1, strc, nspin, fout):
    def is_float(str):
        try:
            float(str)
            return True
        except ValueError:
            return False    
    Elapw = zeros((nspin,strc.nat,in1.nt))
    #lomaxp1 = shape(in1.nlo)[1]
    #nloat = amax(in1.nlo)
    lomaxp1, nloat = shape(in1.nlo)[0], amax(in1.nlo)
    #print 'lomax=', lomaxp1-1, 'nloat=', nloat
    
    Elo   = zeros((nspin,strc.nat,lomaxp1,nloat))
    if nspin==2:
        spflag = ['up','dn']
    else:
        spflag = ['']
    for isp in range(nspin):
        fi = open(case+'.energy'+spflag[isp], 'r')
        for iat in range(strc.nat):
            line1 = fi.next()
            line2 = fi.next()
            possible_spacing = [9,12]
            spacing_type = sign(in1.l_newlo)  # spacing = if l_newlo==0 9 else 12
            
            # checking if this is really how linearization energies are written?
            correct_spacing=True
            spacing = possible_spacing[spacing_type] 
            for d in [line1[spacing*i:spacing*(i+1)] for i in range(len(line1)/spacing)]: 
                if not is_float(d):
                    correct_spacing=False
                    break
            if not correct_spacing:
                spacing_type = 1-spacing_type
                spacing = possible_spacing[spacing_type] 
            
            _elapw_ = map(float,[line1[spacing*i:spacing*(i+1)] for i in range(len(line1)/spacing)])
            _elo_   = map(float,[line2[spacing*i:spacing*(i+1)] for i in range(len(line2)/spacing)])
            Elapw[isp,iat,:] = _elapw_[:in1.nt]

            _nloat_ = int( len(_elo_)/lomaxp1 )
            if _nloat_ * lomaxp1 != len(_elo_):
                print 'WARNING lomax+1=', lomaxp1, 'and number of energies in case.energy is ', len(_elo_), 'which is not compatible'
            _elo_ = _elo_[:lomaxp1*_nloat_]
            _elo_ = reshape(_elo_, (_nloat_,lomaxp1) ).T

            # previoulsy wrong version
            #Elo[isp,iat,l,ilo] = _elo_[:,:nloat]
            for l in range(lomaxp1):
                if in1.lapw[l,iat]: # This is LAPW+LO, hence local orbitals start at ilo=1 and not at ilo=0
                    Elo[isp,iat,l,1:nloat] = _elo_[l,0:(nloat-1)]
                else:               # This is APW+lo, hence first local orbital is actually not LO
                    Elo[isp,iat,l,:] = _elo_[l,:nloat]

            
            for l in range(in1.nt):
                if not in1.lapw[l,iat]:
                    Elapw[isp,iat,l] -= 200
                #if in1.lapw[l,iat]:
                #    umt[0,l,iat,isp] = Elapw[l,iat,isp]
                #else:
                #    umt[0,l,iat,isp] = Elapw[l,iat,isp]-200
                #if l < lomaxp1:
                #    abcelo[4,:,l,iat,isp] = Elo[l,:,iat,isp]
            
            print >> fout, '-'*80
            print >> fout, 'Linearization energy (read from energy file) for spin='+str(isp+1)+' atom '+strc.aname[iat]+':'
            for l in range(in1.nt):
                print >> fout, 'Eapw[l='+str(l)+']=',Elapw[isp,iat,l]
            print >> fout
            for l in range(lomaxp1):
                for ilo in range(nloat):
                    if Elo[isp,iat,l,ilo]  < 90000.0:
                        print >> fout, 'Elo[ilo='+str(ilo)+',l='+str(l)+']=', Elo[isp,iat,l,ilo]
            print >> fout
    return (Elapw, Elo)
                        
def Read_Radial_Potential(case, nat, nspin, nrpt, fout):
    nrad = max(nrpt)
    fid = open(case+'.vsp','r')
    Vr = zeros((nspin,nat,nrad))
    for isp in range(nspin):
        line = fid.next()
        line = fid.next()
        line = fid.next()
        for iat in range(nat):
            npt = nrpt[iat]
            line = fid.next()
            dat = line.split()
            if int(dat[2]) != iat+1:
                print 'WARN: ATOMNUMBER seems wrong in reading '+case+'.vsp file iat=', iat+1, 'and dat=', dat
            line = fid.next()
            dat = line.split()
            num_lm = int(dat[3])
            line = fid.next()
            line = fid.next()
            line = fid.next()
            dat = line.split()
            r_l,r_m = int(dat[3]), int(dat[5])
            line = fid.next()
            
            ipt = 0
            while (fid):
                line = fid.next()
                for i in range(4):
                    Vr[isp,iat,ipt] = float(line[(3+i*19):(3+(i+1)*19)])
                    ipt += 1
                    if ipt >= npt: break
                if ipt >= npt : break

            for i in range(6):
                line = fid.next()
    fid.close()
    print >> fout, 'potential read from file '+case+'.vsp'
    return Vr



#def Read_energy_file(spflag, case, strc, fout, give_kname=False,Extended=False):
#    fi = open(case+'.energy'+spflag, 'r')
#    filename = case+'.energy'+spflag
def Read_energy_file(filename, strc, fout, give_kname=False,Extended=False):
    fi = open(filename, 'r')
    lines = fi.readlines()
    fi.close()
    iln = 2*strc.nat  # bug jul.7.2020
    if give_kname: knames=[]
    Ebnd=[]
    klist=[]
    wegh=[]
    hsrws=[]
    for ik in range(1000000):
        if iln >= len(lines): break
        line = lines[iln]; iln += 1
        if Extended:
            kks = [line[27*i:27*(i+1)] for i in range(3)]
            kname = line[81:91]
            hsr = [line[91+6*i:91+6*(i+1)] for i in range(2)]
            ws = line[103:108]
        else:
            kks = [line[19*i:19*(i+1)] for i in range(3)]
            kname = line[57:67]
            hsr = [line[67+6*i:67+6*(i+1)] for i in range(2)]
            ws = line[79:84]
        kp = map(float,kks)
        hsrows, Ne = int(hsr[0]), int(hsr[1])
        wgh = float(ws)
        klist.append(kp)     # kvecs2
        if give_kname: knames.append(kname)
        wegh.append(wgh)     # wk2
        hsrws.append(hsrows) # ngkir
        #print kp, kname, hsrows, Ne, wgh
        Ebands=[]
        for ib in range(Ne):
            line = lines[iln]; iln += 1
            dd = line.split()
            ibp1, ee = int(dd[0]), float(dd[1])
            Ebands.append(ee)
            #print ib+1, ibp1, ee
        Ebnd.append( Ebands )  # bande
        # He multiplies Ebnd by 0.5 because of the units!
    if give_kname:
        return (klist, wegh, Ebnd, hsrws, knames)
    else:
        return (klist, wegh, Ebnd, hsrws)

def Read_xc_file_dimension(case, strc, fout, ReadAll=False):
    fi = open(case+'.vxc', 'r')
    lines = fi.readlines()
    fi.close()
    iln = 4
    for iat in range(strc.nat):
        npt = strc.nrpt[iat]
        lcmax = int(lines[iln][15:18]); iln += 3
        Vxc = zeros(npt)
        for lxc in range(lcmax):
            lm = int(lines[iln][15:18]), int(lines[iln][23:25])
            #print lm
            iln += 2
            ipt = 0
            for il in range(10000):
                line = lines[iln]; iln+=1
                if ReadAll:
                    for i in range(4):
                        w = float(line[(3+19*i):(3+19*(i+1))])
                        Vxc[ipt]=w; ipt += 1
                        if ipt >= npt: break
                else:
                    ipt += 4
                if ipt>=npt: break
            iln += 2
        iln += 5   # bug jul. 7, 2020
    iln += 1       # bug jul. 7, 2020
    nksxc = int(lines[iln][13:19])
    iln += 1
    kxcmax = 0
    for ikxc in range(nksxc):
        line = lines[iln]; iln += 1
        ksxc = map(int, [line[3+5*i:3+5*(i+1)] for i in range(3)])
        if ReadAll:
            w = float(line[18:(18+19)]) + float(line[(18+19):(18+2*19)])*1j
        kxcmax = max( kxcmax, max(map(abs,ksxc)) )
    return kxcmax

def Read_xc_file(case, strc, fout):
    fi = open(case+'.vxc', 'r')
    lines = fi.readlines()
    fi.close()
    iln = 4
    #lxcm=[]
    lmxc=[]
    Vxclm=[]
    for iat in range(strc.nat):
        npt = strc.nrpt[iat]
        ncm = int(lines[iln][15:18]); iln += 3
        #lxcm.append( ncm )
        Vxc = zeros((ncm,npt))
        _lmxc_=[]
        for lxc in range(ncm):
            lm = [int(lines[iln][15:18]), int(lines[iln][23:25])]
            _lmxc_.append(lm)
            iln += 2
            ipt = 0
            for il in range(10000):
                line = lines[iln]; iln+=1
                for i in range(4):
                    w = float(line[(3+19*i):(3+19*(i+1))])
                    Vxc[lxc,ipt]=w; ipt += 1
                    if ipt >= npt: break
                if ipt>=npt: break
            iln += 2
        iln += 5     # bug jul.7 2020
        lmxc.append( _lmxc_ ) # lmxc[iat][lxc][0-1]
        Vxclm.append( Vxc )  # Vxclm[iat][lxc,ir]
    iln += 1         # bug jul.7 2020
    nksxc = int(lines[iln][13:19])
    iln += 1
    kxcmax = 0
    Vxcs = zeros((nksxc),dtype=complex)
    ksxc = zeros((nksxc,3),dtype=int)
    for ikxc in range(nksxc):
        line = lines[iln]; iln += 1
        ksxc[ikxc,:] = map(int, [line[3+5*i:3+5*(i+1)] for i in range(3)])
        Vxcs[ikxc] = float(line[18:(18+19)]) + float(line[(18+19):(18+2*19)])*1j
        #kxcmax = max( kxcmax, max(map(abs,ksxc)) )

    print >> fout, "*** Number of lm component for each atom "
    print >> fout, "   iat   lxcm(iat)"  
    for iat in range(strc.nat):
        print >> fout, '%5d %10d' % (iat, len(lmxc[iat]))
    print >> fout, "# Number of Institial plane waves from vxc=", nksxc
    return (lmxc, Vxclm, ksxc, Vxcs)

def Read_vector_file(case, strc, fout):
    vectortype = float
    vecread    = rdVec.fvread3
    vecwrite   = rdVec.fvwrite3
    so = ''
    if os.path.isfile(case+".inso") and os.path.getsize(case+".inso")>0:
        print >> fout, 'Found '+case+'.inso file, hence assuming so-coupling exists. Switching -so switch!'
        so = 'so'
        vectortype=complex
        
    if os.path.isfile(case+".in1c") and os.path.getsize(case+".in1c")>0:
        print >> fout, 'Found '+case+'.in1c file, hence assuming complex eigenvector. Switching to complex!'
        so = ''
        vectortype=complex
    
    if vectortype==complex:
        vecread = rdVec.fvread3c
        vecwrite = rdVec.fvwrite3c
    
    maxkpoints = 10000
    # opens vector file
    heads=[]
    all_Gs=[]
    all_As=[]
    all_Ek=[]
    fname, tape = case+'.vector'+so, 9
    Elinear = rdVec.fvopen(tape, fname, strc.nat)
    print >> fout, 'linearization energy=', Elinear
    for ik in range(maxkpoints):
        # Reads vector file
        head = rdVec.fvread1(tape)
        (k, kname, wgh, ios, n0, nb) = head
        if ios!=0: break # vector file is finished, no more k-points
        print >> fout, 'k=', k, 'nrows=', n0, 'nbands=', nb
        heads.append(head)
        # Reciprocal vectors
        Gs = rdVec.fvread2(tape, n0)
        all_Gs.append(Gs.T)
        # Reading KS eigensystem
        As=zeros((nb,n0), dtype=vectortype)
        Ek=zeros(nb, dtype=float)
        for i in range(nb):
            (num, ek, A) = vecread(tape, n0)
            As[i,:] = A  # KS eigenvector
            Ek[i] = ek   # KS eigenvalue
        all_As.append(As)
        all_Ek.append(Ek)
    rdVec.fvclose(tape)
    return (heads, all_Gs, all_As, all_Ek)
    
    
class RadialFunctions:
    def __init__(self, in1, strc, Elapw, Elo, Vr, nspin, fout):
        
        self.umt = zeros((nspin,strc.nat,in1.nt,5))
        
        import radials as rd
        Vr *= 0.5 # Converting from Rydberg to Hartree
        
        nrad = max(strc.nrpt)
        self.ul    = zeros((nspin,strc.nat,in1.nt,nrad))
        self.us   = zeros((nspin,strc.nat,in1.nt,nrad))
        self.udot = zeros((nspin,strc.nat,in1.nt,nrad))
        self.usdot= zeros((nspin,strc.nat,in1.nt,nrad))
        
        store_nodel = zeros((strc.nat,in1.nt),dtype=int)
        store_nodeu = zeros((strc.nat,in1.nt),dtype=int)
        
        for isp in range(nspin):
            for iat in range(strc.nat):
                print >> fout, "-----  Calc radfunc for atom ",iat
                print >> fout, "-- (L)APW states --"
                print >> fout, (' '*10+'POTENTIAL PARAMETERS FOR JATOM=%-3d name=%-10s') % (iat,strc.aname[iat])
                print >> fout, ' '*11+'L'+(' '*7)+'U(R)'+(' '*11)+"U'(R)",(' '*9)+'DU/DE'+(' '*10)+"DU'/DE"+(' '*9)+"NORM-U'"
                # val/cond states:  calculate  the radial functions u, udot and ulo
                npt = strc.nrpt[iat]
                dh  = log(strc.rmt[iat]/strc.r0[iat])/(npt - 1)      # logarithmic step for the radial mesh
                for l in range(in1.nt):
                    eh = 0.5*Elapw[isp,iat,l]  # converting to Hartrees
                    #       Calculate the function at el (u_l(r,E_l)):
                    a,b,nodes,uv,duv = rd.outwin(strc.rel, Vr[isp,iat,:], strc.r0[iat], dh, npt, eh, float(l), strc.Znuc[iat])
                    #     Calculate |u_l(r,E_l)|^2:
                    ovlp = rd.rint13g(strc.rel, a, b, a, b, dh, npt, strc.r0[iat])
                    #     Calculate the normalization factor:
                    trx = 1/sqrt(ovlp)
                    #     Store the normalized values at the boundaries in p and dp:
                    uv *= trx
                    duv *= trx
                    #     Normalize the function u_l:
                    a *= trx
                    b *= trx
                    
                    #print >> fout, 'uv[is='+str(isp)+',iat='+str(iat)+',l='+str(l)+']=', uv, 'duv=', duv
                    
                    dele = 2e-3 # the up and downward energy-shift in Hartrees
                    #delei = 0.25/dele
                    #  Calculate u_l(r,E_1=E_l-\Delta E)
                    ac,bc,nodel,uvc,duvc = rd.outwin(strc.rel, Vr[isp,iat,:], strc.r0[iat], dh, npt, eh-dele, float(l), strc.Znuc[iat])
                    ovlp = rd.rint13g(strc.rel, ac, bc, ac, bc, dh, npt, strc.r0[iat])
                    trx = 1/sqrt(ovlp)
                    uvc *= trx
                    duvc *= trx
                    ac *= trx
                    bc *= trx
                    
                    ae,be,nodeu,uve,duve = rd.outwin(strc.rel, Vr[isp,iat,:], strc.r0[iat], dh, npt, eh+dele, float(l), strc.Znuc[iat])
                    ovlp = rd.rint13g(strc.rel, ae, be, ae, be, dh, npt, strc.r0[iat])
                    trx = 1/sqrt(ovlp)
                    uve  *= trx
                    duve *= trx
                    ae *= trx
                    be *= trx
                    
                    uve = (uve - uvc)*(0.25/dele)
                    duve= (duve-duvc)*(0.25/dele)
                    ae = (ae-ac)*(0.25/dele)
                    be = (be-bc)*(0.25/dele)
                    
                    #Insure ortogonalization
                    # Calculate <u_l|udot_l>
                    cross = rd.rint13g(strc.rel, a, b, ae, be, dh, npt, strc.r0[iat])
                    if( cross > 0.05): print >> fout, 'For l='+str(l)+' correction='+str(-cross)+' overlap='+str(ovlp)
                    # Set orthogonalized udot_l
                    ae -= cross * a
                    be -= cross * b
                    uve  -= cross * uv  # pe
                    duve -= cross * duv # dpe
                    
                    # Store the orthogonalized values at the boundaries in pe and dpe:
                    self.umt[isp,iat,l,0] = uv   # p
                    self.umt[isp,iat,l,2] = duv  # dp
                    self.umt[isp,iat,l,1] = uve  # pe
                    self.umt[isp,iat,l,3] = duve # dpe
                    # Calculate |udot_l(r,E_l)|^2 and store it in  pei(j,iat)
                    self.umt[isp,iat,l,4] = rd.rint13g(strc.rel, ae, be, ae, be, dh, npt, strc.r0[iat])
                    
                    self.ul   [isp,iat,l,:npt] = a[:]
                    self.us   [isp,iat,l,:npt] = b[:]
                    self.udot [isp,iat,l,:npt] = ae[:]
                    self.usdot[isp,iat,l,:npt] = be[:]
                    
                    print >> fout, (' '*10+'%2d '+'%14.6f '*5+' '*5+'%2d '*3) % (l,self.umt[isp,iat,l,0],self.umt[isp,iat,l,2],self.umt[isp,iat,l,1],self.umt[isp,iat,l,3],self.umt[isp,iat,l,4],nodel,nodes,nodeu)
                    
                    store_nodel[iat,l] = nodel
                    store_nodeu[iat,l] = nodeu
        
        #  Calculate the radial function for local orbitals.
        nsp,nat,lomaxp1,nloat = shape(Elo)
        self.umtlo = zeros((nspin,strc.nat,lomaxp1,nloat,4))
        self.ulo  = zeros((nspin,strc.nat,lomaxp1,nloat,nrad))
        self.uslo = zeros((nspin,strc.nat,lomaxp1,nloat,nrad))
        #self.nLO_at = zeros((strc.nat,self.lomax+1), dtype=intc)
        for isp in range(nspin):
            for iat in range(strc.nat):
                print >> fout, (' '*10+'LOCAL ORBITAL POTENTIAL PARAMETERS FOR JATOM=%-3d name=%-10s') % (iat,strc.aname[iat])
                print >> fout, (' '*11)+'L'+(' '*7)+'U(R)'+(' '*10)+"U'(R)"+(' '*9)+'<u|u_LO>'+(' '*6)+'<u_dot|u_LO>'+' num. nodes'+ ' Elo'
                npt = strc.nrpt[iat]
                dh  = log(strc.rmt[iat]/strc.r0[iat])/(npt - 1)      # logarithmic step for the radial mesh
                for l in range(lomaxp1):
                    for ilo in in1.nLO_at_ind[iat][l]:
                        el = 0.5*Elo[isp,iat,l,ilo]   # in Hartrees
                        #print >> fout, 'Calculating function at l=', l, 'ilo=', ilo, 'Elo=', el*2.0
                        #       Calculate the function at el (u_l(r,Elo_l)):
                        a,b,nodes,uv,duv = rd.outwin(strc.rel, Vr[isp,iat,:], strc.r0[iat], dh, npt, el, float(l), strc.Znuc[iat])
                        #     Calculate |u_l(r,E_l)|^2:
                        ovlp = rd.rint13g(strc.rel, a, b, a, b, dh, npt, strc.r0[iat])
                        #     Calculate the normalization factor:
                        trx = 1/sqrt(ovlp)
                        #     Normalize the function u_l:
                        uv  *= trx   # plo
                        duv *= trx   # dplo
                        a *= trx
                        b *= trx
                        #   Calculate pi12lo(l,iat)=<u_l(r,E_l)|u_l(r,E_{lo})>:
                        pi12lo = rd.rint13g(strc.rel, self.ul[isp,iat,l,:npt], self.us[isp,iat,l,:npt], a, b, dh, npt, strc.r0[iat])
                        #   Calculate pe12lo(l,iat)=<\dot{u}_l(r,E_l)|u_l(r,E_{lo})>:
                        pe12lo = rd.rint13g(strc.rel, self.udot[isp,iat,l,:npt], self.usdot[isp,iat,l,:npt], a, b, dh, npt, strc.r0[iat])
                        
                        self.umtlo[isp,iat,l,ilo,0] = uv      # plo
                        self.umtlo[isp,iat,l,ilo,1] = duv     # dplo
                        self.umtlo[isp,iat,l,ilo,2] = pi12lo  # pi12lo
                        self.umtlo[isp,iat,l,ilo,3] = pe12lo  # pe12lo
                        #   Store the radial wave function :
                        self.ulo [isp,iat,l,ilo,:npt] = a[:]
                        self.uslo[isp,iat,l,ilo,:npt] = b[:]
                        
                        print >> fout, ((' '*10)+'%2d '+('%14.6f'*4)+(' '*5)+'  %2d   '+' %10.6f') % (l,uv,duv,pi12lo,pe12lo,nodes,el*2)
                for l in range(lomaxp1):
                    if len(in1.nLO_at_ind[iat][l])>0:
                        print >> fout, ' '*10+'LOCAL ORBITALS OVERLAPS:'
                        print >> fout, ' '*10+'jlo   klo      <u|u>           E[jlo]          E[klo]'
                    for jlo in in1.nLO_at_ind[iat][l]:
                        for klo in in1.nLO_at_ind[iat][l]:
                            pilolo = rd.rint13g(strc.rel, self.ulo[isp,iat,l,jlo,:], self.uslo[isp,iat,l,jlo,:], self.ulo[isp,iat,l,klo,:], self.uslo[isp,iat,l,klo,:], dh, npt, strc.r0[iat])
                            print >> fout, ((' '*7)+('%5d '*2)+('%15.6f '*3)) % (jlo,klo,pilolo,Elo[isp,iat,l,jlo],Elo[isp,iat,l,klo])
                    
    def get_ABC(self, in1, strc, fout):
        """Calculate the cofficients A,B,C for (L)APW and the local orbitals"""
        cutoff = 200
        #nspin, nat, nt, x = shape(self.umt)
        nspin, nat, lomaxp1, nloat, x = shape(self.umtlo)
        self.abcelo = zeros((nspin,nat,lomaxp1,nloat,3))
        for isp in range(nspin):
            for iat  in range(nat):
                Rmt = strc.rmt[iat]
                for l in range(lomaxp1):
                    #   P(l)   = ul_Rmt(1,l,jatom)                                                                                                                                                 
                    #   PE(l)  = ul_Rmt(2,l,jatom)                                                                                                                                                 
                    #   DP(l)  = dul_Rmt(1,l,jatom)                                                                                                                                                
                    #   DPE(l) = dul_Rmt(2,l,jatom)                                                                                                                                                
                    #   PLO(l) = ul_Rmt(2+jlo,l,jatom)                                                                                                                                             
                    #   DPLO(l)= dul_Rmt(2+jlo,l,jatom)                                                                                                                                            
                    #                                                                                                                                                                              
                    p   = self.umt[isp,iat,l,0] # p
                    pe  = self.umt[isp,iat,l,1] # pe
                    dp  = self.umt[isp,iat,l,2] # dp
                    dpe = self.umt[isp,iat,l,3] # dpe
                    pei = self.umt[isp,iat,l,4] # pei
                    if in1.lapw[l,iat]==0:  # APW+lo, hence ony u and dotu are used
                        alonorm=sqrt(1 + (p/pe)**2 * pei)
                        alo = 1/alonorm
                        blo = -p/(pe*alonorm)
                        clo = 0
                        #
                        ilo=0
                        self.abcelo[isp,iat,l,ilo,0] = alo
                        self.abcelo[isp,iat,l,ilo,1] = blo
                        self.abcelo[isp,iat,l,ilo,2] = clo
                        print >> fout, 'lo coefficient: iat=%-2d l=%-2d ilo=%-1d lapw=%1d a=%-12.7f b=%-12.7f c=%-12.7f' % (iat,l,ilo,in1.lapw[l,iat],self.abcelo[isp,iat,l,ilo,0],self.abcelo[isp,iat,l,ilo,1],self.abcelo[isp,iat,l,ilo,2])
                    for ilo in in1.nLO_at_ind[iat][l]:
                        plo    = self.umtlo[isp,iat,l,ilo,0]  # plo
                        dplo   = self.umtlo[isp,iat,l,ilo,1]  # dplo
                        pi12lo = self.umtlo[isp,iat,l,ilo,2]  # pi12lo
                        pe12lo = self.umtlo[isp,iat,l,ilo,3]  # pe12lo
                        if in1.lapw[l,iat]:   # This is LAPW+LO
                            #  We construct the LO orbtial as   u_new = ALO*u + BLO*dotu + CLO*u_LO                                                                                                        
                            #    and require       u_new(R) = 0                                                                                                                                            
                            #                  du_new/dr(R) = 0                                                                                                                                            
                            #         and     <u_new|u_new> = 1                                                                                                                                            
                            #  which leads to:                                                                                                                                                             
                            # xac =  (u_LO*ddotu-du_LO*dotu)*R^2                                                                                                                                           
                            # xbc = -(u_LO*du - du_LO*u)*R^2                                                                                                                                               
                            # clo =   1/sqrt( 1 + xac*(xac + 2*<u_LO|u>)+xbc*(xbc*<dotu|dotu>+2*<dotu|u_LO> )                                                                                              
                            # alo = xac/sqrt( 1 + xac*(xac + 2*<u_LO|u>)+xbc*(xbc*<dotu|dotu>+2*<dotu|u_LO> )                                                                                              
                            # blo = xbc/sqrt( 1 + xac*(xac + 2*<u_LO|u>)+xbc*(xbc*<dotu|dotu>+2*<dotu|u_LO> )                                                                                              
                            xac = (plo*dpe - dplo*pe)*Rmt**2
                            xbc = -(plo*dp - dplo*p)*Rmt**2
                            clo = 1/sqrt(1 + xac*(xac + 2*pi12lo) + xbc*(xbc*pei + 2*pe12lo))
                            clo = min(clo,cutoff)
                            alo = clo*xac
                            blo = clo*xbc
                            #print >> fout, '%s%12.7f '*9 % ('debug_lo p=', p, 'dp=', dp, 'plo=', plo, 'dplo=',  dplo, 'pe=', pe, 'dpe=', dpe, 'pi12lo=', pi12lo, 'pe12lo=', pe12lo, 'rm=', Rmt)
                        else: #  must be APW+lo+LO
                            # We construct the LO orbital as u_new = ALO*u + CLO*u_LO                                                                                                                   
                            #   and require       u_new(R) = 0                                                                                                                                          
                            #       and      <u_new|u_new> = 1                                                                                                                                          
                            # which leads to:                                                                                                                                                           
                            # xac = sqrt(1 + (u/u_LO)**2 - 2*(u/u_LO)*<u|u_LO>)                                                                                                                         
                            # ALO = 1/sqrt(1 + (u/u_LO)**2 - 2*(u/u_LO)*<u|u_LO>)                                                                                                                       
                            # CLO = -u/u_LO /sqrt(1 + (u/u_LO)**2 - 2*(u/u_LO)*<u|u_LO>)                                                                                                                
                            xbc=-p/plo
                            xac=sqrt(1+xbc**2+2*xbc*pi12lo)
                            alo = 1/xac
                            blo = 0
                            clo = xbc/xac
                        self.abcelo[isp,iat,l,ilo,0] = alo
                        self.abcelo[isp,iat,l,ilo,1] = blo
                        self.abcelo[isp,iat,l,ilo,2] = clo
                        print >> fout, 'lo coefficient: iat=%-2d l=%-2d ilo=%-1d lapw=%1d a=%-12.7f b=%-12.7f c=%-12.7f' % (iat,l,ilo,in1.lapw[l,iat],self.abcelo[isp,iat,l,ilo,0],self.abcelo[isp,iat,l,ilo,1],self.abcelo[isp,iat,l,ilo,2])

class CoreStates:
    def __init__(self, case, strc, nspin, fout):
        self.l2kappa={'S ': -1, 'P ': -2, 'PP': 1, 'D ': -3, 'DD': 2, 'F ': -4, 'FF': 3}
        fin = open(case+'.core','r')
        print >> fout, 'read core states occupuation information!'
        self.occ_inc = [[] for iat in range(strc.nat)]
        ncore = zeros(strc.nat,dtype=int)
        for iat in range(strc.nat):
            dat = fin.next().split()
            norb, shift, iprint = int(dat[0]), float(dat[1]), int(dat[2])
            ncore[iat] = norb
            for iorb in  range(norb):
                dat = fin.next().split(',')
                qn_n, qn_kappa, occ = int(dat[0]), int(dat[1]), int(dat[2].split()[0])
                self.occ_inc[iat].append( occ )
                if (occ < 1e-2):
                    ncore[iat] -= 1
                #print 'n='+str(qn_n)+' kappa=%2d' % (qn_kappa,)+' occ='+str(occ)
            print >> fout, 'Number of core states at '+str(iat)+'-th atom:', ncore[iat]
        line = fin.next()
        ncoremax = max(ncore)
        self.eig_core = [[[] for i in range(strc.nat)] for j in range(nspin)]
        self.l_core =  [[] for j in range(strc.nat)]
        self.corind = []
        n_sym_kap_ocm=[[] for iat in range(strc.nat)]

        ########
        #import radials as rd        # radial wave functions
        if ncoremax != 0:
            # core states present
            isp=0   # spin up
            for iat in range(strc.nat):
                line = fin.next()
                line = fin.next()
                t_atomname, t_norb = line[20:30], int(line[41:43])
                if (t_norb != len(self.occ_inc[iat])):
                    print >> fout, 'WARNING : reading '+case+'.core : norb read on top=', len(self.occ_inc[iat]),  'while t_norb=', t_norb
                print >> fout, 'N t  nm oc kp l Ene'
                for iorb in range(t_norb):
                    line = fin.next()
                    t_n, t_symb, t_nm, t_ec = int(line[1:2]), line[2:4], line[9:12], float(line[20:40])
                    if t_norb==1 and abs(t_ec) < 1e-5: continue
                    if self.occ_inc[iat][iorb] < 1e-2:
                        print >> fout, ' the core state ', t_nm, ' is frozen'
                        continue
                    kappa = self.l2kappa[t_symb]
                    lc = abs(kappa)-(1-sign(kappa))/2
                    self.l_core[iat].append( lc )
                    self.eig_core[isp][iat].append( t_ec )  # *0.5 Ry2Hartree
                    n_sym_kap_ocm[iat].append( (t_n,t_nm,kappa,abs(kappa)*2) )
                    Ry2H = 0.5  
                    print >> fout, t_n,t_symb,t_nm, self.occ_inc[iat][iorb], '%2d %1d' % (kappa,lc), t_ec*Ry2H
            if nspin==2:
                isp=1
                for iat in range(strc.nat):
                    line = fin.next()
                    t_atomname, t_norb = line[20:30], int(line[41:43])
                    for iorb in range(t_norb):
                        line = fin.next()
                        t_n, t_symb, t_nm, t_ec = int(line[1:2]), line[2:4], line[9:12], float(line[20:40])
                        if self.occ_inc[iat][iorb] < 1e-2:
                            continue
                        self.eigcore[isp][iat].append( t_ec) # *0.5 Ry2Hartree
            
            nclm_max  = max([ sum([2*lc+1 for lc in self.l_core[iat]]) for iat in range(strc.nat)])
            #l_core_max = max([ max(self.l_core[iat]) for iat in range(strc.nat)])
            l_core_max = max([ max(self.l_core[iat]) if self.l_core[iat] else 0 for iat in range(strc.nat)])
            print >> fout, 'l_core_max=', l_core_max, 'nclm_max=', nclm_max
            
            # self.ul_core[isp][iat][ic][ipt]
            #self.ul_core =  [[[[] for i in range(ncore[iat])] for iat in range(strc.nat)]]
            #self.us_core =  [[[[] for i in range(ncore[iat])] for iat in range(strc.nat)]]
            self.ul_core =  [[[] for iat in range(strc.nat)] for isp in range(nspin)]        # bug jul.7 2020
            self.us_core =  [[[] for iat in range(strc.nat)] for isp in range(nspin)]        # bug jul.7 2020
            for  isp in range(nspin):
                for iat in range(strc.nat):
                    #######
                    #rp, dh, npt = strc.radial_mesh(iat)  # logarithmic radial mesh
                    npt = strc.nrpt[iat]
                    norb = int(fin.next().split()[0])
                    print >> fout, ' total number of core wave functions from '+case+'.core = ', norb
                    for iorb in range(norb):
                        line = fin.next()
                        nm = line[17:20]
                        print >> fout, ' read orbital '+nm+' with occ=',self.occ_inc[iat][iorb]
                        ulc = zeros((2,npt))
                        for idot in [0,1]:
                            ipt = 0
                            while(fin):
                                line = fin.next()
                                for i in range(4):
                                    ulc[idot,ipt] = float(line[3+19*i:3+19*(i+1)])
                                    ipt += 1
                                    if ipt >= npt: break
                                if ipt >= npt: break
                        if self.occ_inc[iat][iorb] < 1e-2:
                            print >> fout, "  exclude core state ", nm 
                        else:
                            ic = len(self.ul_core[isp][iat])
                            lc = self.l_core[iat][ic]
                            t_n, t_nm, kappa, max_occ = n_sym_kap_ocm[iat][ic]
                            if (nm != t_nm):
                                print >> fout, 'WARNING reading '+case+'.core at iat='+str(iat)+' iorb='+str(iorb)+' nm='+nm+' t_nm='+t_nm
                            norm = sqrt(0.5*max_occ/(2*lc+1.0))
                            self.ul_core[isp][iat].append( ulc[0,:]*norm )
                            self.us_core[isp][iat].append( ulc[1,:]*norm )
                            #####
                            #norm0 = rd.rint13g(strc.rel, ulc[0,:], ulc[1,:], ulc[0,:], ulc[1,:], dh, npt, strc.r0[iat])
                            print >> fout, '    and the wave function for state ', nm,'. It is renormalized with coefficient %10.7f to account for proper degeneracy.' % (norm**2,)
            fin.close()
        else:
            nclm_max = 0
            
        ntot =  0
        ne_core = 0
        for iat in range(strc.nat):
            ntot = ntot + strc.Znuc[iat]*strc.mult[iat]
            occoreat = 0
            for ic in range(len(n_sym_kap_ocm[iat])):
                t_n, t_nm, kappa, max_occ = n_sym_kap_ocm[iat][ic]
                occoreat += max_occ
            ne_core += strc.mult[iat]*occoreat
        nval = ntot - ne_core
        nelfc = 0 # for now just neglect that, but maybe you should revisit this!!
        print >> fout
        print >> fout, ' '*10+'Total number of electrons:      %8.4f' % (ntot,)
        print >> fout, ' '*10+'Number of valence electrons:    %8.4f' % (nval,)
        print >> fout, ' '*10+'Number of core electrons:       %8.4f' % (ne_core,)
        print >> fout, ' '*10+'Number of frozen core electrons:%8.4f' % (nelfc,)
        self.ntot, self.nval, self.ncore = ntot, nval, ne_core
        
        # setup corind
        iateq_ind=[]
        for iat in range(len(strc.mult)):
            iateq_ind += [iat]*strc.mult[iat]
            
        ndf = len(iateq_ind)
        self.nclm = zeros(ndf, dtype=intc)
        self.clmind=zeros((ndf,nclm_max), dtype=intc)
        remember_weight=[]
        for idf,iat in enumerate(iateq_ind): # idf -- atom index counting all atoms, iat -- atom index counting only equivalent
            iclm = 0
            for ic,lc in enumerate(self.l_core[iat]):
                norm2 = 0.5*n_sym_kap_ocm[iat][ic][3]/(2*lc+1.0)
                for mc in range(-lc,lc+1):
                    self.clmind[idf,iclm] = len(self.corind)
                    self.corind.append( [iat,idf,ic,lc,mc] )
                    remember_weight.append( norm2 )
                    iclm=iclm+1
            self.nclm[idf] = iclm

        print >> fout, ' '*10+'Total Number of Core states:     ', sum(self.nclm), len(self.corind)
        self.corind = array(self.corind, dtype=int)
        if True:
            print >> fout, 'corrind index='
            for ic in range(len(self.corind)):
                print >> fout, 'iat=%3d idf=%3d ic=%3d lc=%3d mc=%3d' % tuple(self.corind[ic]), ' with weight=sqrt(%10.7f)' % (remember_weight[ic],)
        
        #print self.nclm
        #print self.clmind
        #print 'corind=', self.corind
        
class SCF2:
    def __init__(self, case, spin=''):
        self.case = case
        self.extn = 'scf2'+spin  # spin can be '', 'up' or 'dn'
        self.parse()

    def filename(self):
        return self.case + '.' + self.extn

    def parse(self):
        f = open(self.filename(), 'r')
        lines = f.readlines()
        f.close()

        # use [-1] to take values from last iteration
        self.NOE = [float(line[38:]) for line in lines if re.match(r':NOE', line)][-1]
        self.EF  = [float(line[38:]) for line in lines if re.match(r':FER', line)][-1]

    def debugprint(self, f=sys.stdout):
        print >> f, 'From file `%s`', self.filename()
        print >> f, 'NOE =', self.NOE
        print >> f, 'EF  =', self.EF, 'Ry'
