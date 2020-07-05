from scipy import *
from scipy import optimize
from numpy import linalg
from scipy import interpolate

from inout import *
import for_tetrahedra as ft
import fnc
import for_pade as fpade

from pylab import *

def Give2TanMesh(x0,L,Nw):
    def fun(x,x0,L,Nw):
        "x[0]=d, x[1]=w"
        d=x[0]
        w=x[1]
        #print 'd=', d, 'w=', w
        return array([L-w/tan(d), x0-w*tan(pi/(2*Nw)-d/Nw) ])
    tnw = tan(pi/(2*Nw))
    if x0 > L*0.25*Nw*tnw**2:
        x0 = L*0.25*Nw*tnw**2-1e-15
    xi=x0/L
    d0 = Nw/2.*(tnw-sqrt(tnw**2 - 4*xi/Nw))
    w0 = L*d0
    sol=optimize.root(fun, [d0,w0], args=(x0,L,Nw) )
    (d,w) = sol.x
    xt = linspace(0.0,1.0,2*Nw)*(pi-2*d) - pi/2 + d
    om = w*tan(xt)
    dom = (w*(pi-2*d)/(2.*Nw))/cos(xt)**2
    return om,dom


def calc_Fermi(Ebnd, kqm_atet, kqm_wtet, nval, nspin):
    eDos = 1.0
    if nspin==1 and (nval % 2 == 0):
        # just simply counting bands, and setting EF in the middle
        nvm   = int( nval/2.0 + 0.51 )
        evbm  = max( Ebnd[:,nvm-1] )
        ecbm  = min( Ebnd[:,nvm])
        EF  = 0.5*(evbm+ecbm)
        ocint = ft.idos(EF, Ebnd, kqm_atet, kqm_wtet)*2.0/nspin
        if ecbm >= evbm and abs(ocint-nval) < 1e-6:
            Eg = ecbm - evbm
            eDos = 0
    if eDos != 0: # the simple method did not succeed
        evbm = min(Ebnd.ravel()) # minimum of energy
        ecbm = max(Ebnd.ravel()) # maximum of energy
        ocmax = sum([ft.idos(ecbm, Ebnd, kqm_atet, kqm_wtet) for isp in range(nspin)])*2.0/nspin
        if ocmax <= nval:
            print 'ERROR in fermi: not enough bands : %s%-10.4f %s%-10.2f %s%-10.2f' % ('emax=', ecbm, 'ocmax=', ocmax, 'nel= ', nval )
            sys.exit(1)
        
        EF = optimize.brentq(lambda x: sum([ft.idos(x, Ebnd, kqm_atet, kqm_wtet) for isp in range(nspin)])*2.0/nspin-nval, evbm, ecbm)
        eDos = sum([ft.dostet(EF, Ebnd, kqm_atet, kqm_wtet) for isp in range(nspin)])*2.0/nspin
        
        # For insulator (including semiconductor, set fermi energy as the middle of gap
        if eDos < 1e-4:
            evbm = max( filter(lambda x: x < EF, Ebnd.ravel()) )
            ecbm = min( filter(lambda x: x > EF, Ebnd.ravel()) )
            EF = (evbm + ecbm)/2.
            Eg = ecbm - evbm
        else: 
            Eg = -eDos
            evbm, ecbm = EF, EF
    return (EF, Eg, evbm, ecbm, eDos)
 
def cart2int(klist,strc,latgen):
    """ Cartesian to integer representation of k-points"""
    if latgen.ortho:
        rbas = zeros((3,3))
        alat = array([strc.a, strc.b, strc.c])
        for j in range(3):
            rbas[j,:] = latgen.rbas[j,:]/alat[j] 
        iklist = dot(klist,rbas)
        return iklist
    else:
        return klist

def Compute_selfc_inside(iq, irk, bands, mwm, fr, kqm, ncg_c, core, Ul, fout, PRINT=False):
    def freqconvl(iom, enk, mwm, omega, womeg):
        """ We are computing frequency convolution
                sigma(iom) = -1/beta \sum_{iOm} mwm[iOm]/(iom-eps+iOm)
            Because mwm[-iOm]=mwm[iOm] and at T=0, we can rewrite
                sigma(iom) = 1/pi Integrate[ mwm[Om]*(eps-iom)/((eps-iom)^2 + Om^2) , {Om,0,Infinity} ]
            Finally, we notice that
                                  Integrate[ 1/((eps-iom)^2 + Om^2) , {Om,0,Infinity} ] = pi*sign(eps)/(2*(eps-iom))
            therefore
                sigma(iom) = (eps-iom)/pi Integrate[ (mwm[Om]-mwm[iom])/((eps-iom)^2 + Om^2) , {Om,0,Infinity} ] + mwm[iom] * sign(eps)/2
        """
        if FORT:
            return fnc.fr_convolution(iom+1, enk, mwm, omega, womeg)
        else:
            eps_om = enk - omega[iom]*1j
            sc0 = sum( (mwm[:]-mwm[iom]) * womeg / (omega**2 + eps_om**2) )
            return eps_om*sc0/pi + 0.5 * mwm[iom] * sign(enk)
    
    (nom, nb1, nb2) = shape(mwm)
    omega, womeg = fr.omega, fr.womeg
    
    ik = kqm.k_ind[irk]   # index in all-kpoints, not irreducible
    jk = kqm.kqid[ik,iq]  # index of k+q
    jrk = kqm.kii_ind[jk] # index of k+q in irreducible list

    #print 'Compute_selfc_inside:shape(bands)=', shape(bands), 'jrk=', jrk, 'nb2=', nb2
    isp=0
    enk = zeros( nb2 )
    for icg in range(ncg_c):
        iat,idf,ic = core.corind[icg,:3]
        enk[icg] = core.eig_core[isp][iat][ic]
    enk[ncg_c:nb2] = bands[jrk,:(nb2-ncg_c)]   # careful, we need G(k+q,:) here, not G(k,:)

    if Ul is not None:
        (nl,nom) = shape(Ul)
        # Cnl[l,ie2,iom] = 1/pi*Integrate[ Ul[l,iOm]*(enk[ie2]-iom)/((enk[ie2]-iom)^2 + Om^2),{Om,0,Infinity}]
        Cnl = fnc.many_fr_convolutions(enk, Ul, omega, womeg) # frequency convolution of all svd functions
        Cnl = reshape(Cnl,(nl*nb2,nom))
        mwm2 = zeros( (nb1, nl*nb2), dtype=complex )
        for ie1 in range(nb1):
            mwm2[ie1,:] = reshape( mwm[:,ie1,:], nl*nb2 )
        sc_p = dot(mwm2,Cnl)
    else:
        if not (fr.iopfreq==4 and len(fr.omega_precise)!=len(fr.omega)):
            sc_p = fnc.all_fr_convolutions(enk, mwm, omega, womeg)
        else:
            sc_p = zeros( (nb1,nom), dtype=complex )
            mwmp = zeros( (nb2,len(fr.omega_precise)), dtype=complex)
            for_asymptotic = (omega**2+1.0)
            for_asymptotic_precise = (fr.omega_precise**2+1.0)
            zeros_nb2 = zeros(nb2, dtype=complex)
            for ie1 in range(nb1):
                for ie2 in range(nb2):
                    # Interpolating W on more dense mesh for frequency convolution
                    # Note that we interpolate W*(om^2+1), which goes to constant at large om, and is easier to interpolate properly
                    mwmr = interpolate.CubicSpline(omega, for_asymptotic*mwm[:,ie1,ie2].real, bc_type=((2,0),(1,0)), extrapolate=True)
                    if sum(abs(mwm[:,ie1,ie2].imag))/nom > 1e-10:
                        mwmi = interpolate.CubicSpline(omega, for_asymptotic*mwm[:,ie1,ie2].imag, bc_type=((2,0),(1,0)), extrapolate=True)
                        mwmp[ie2,:] = (mwmr(fr.omega_precise) + mwmi(fr.omega_precise)*1j)/for_asymptotic_precise
                    else:
                        mwmp[ie2,:] = mwmr(fr.omega_precise)/for_asymptotic_precise
                for iom in range(nom):
                    # We now always subtract this peace, because treating some frequencies and bands differently seems to cause problems.
                    mwm_iom = mwm[iom,ie1,:]
                    sc_p[ie1,iom] = fnc.few_fr_convolutions(enk, mwmp, mwm_iom, fr.omega_precise, omega[iom], fr.womeg_precise)
                    #if irk==3:
                    #    rr = 0
                    #    for ie2 in range(nb2):
                    #        eps_om = enk[ie2] - omega[iom]*1j
                    #        to_sum = (mwmp[ie2,:]-mwm_iom[ie2]) / (fr.omega_precise**2 + eps_om**2)
                    #        sc0 = sum( to_sum * fr.womeg_precise)
                    #        cc = eps_om*sc0/pi + mwm_iom[ie2]/pi * arctan(omega[-1]/(eps_om)) #pi*0.5*sign(enk[ie2])
                    #        rr += cc 
                    #        print 'ie2=', ie2, 'en=', enk[ie2], 'iom=', iom, 'om=', omega[iom], 'cont=', cc
                    #        if ( abs(mwm_iom[ie2])>1e-10 ):
                    #            fo = open('case_to_study.dat', 'w')
                    #            print >> fo, '# ', eps_om
                    #            for i in range(len(fr.omega_precise)):
                    #                print >> fo, fr.omega_precise[i], mwmp[ie2,i].real
                    #            fo.close()
                    #        plot(fr.omega_precise, to_sum, 'o-')
                    #        plot(fr.omega_precise, to_sum*fr.womeg_precise*100, 's-')
                    #        show()
                    #    print 'res=', abs(rr-sc_p[ie1,iom]), rr
            #print 'iq=', iq, 'irk=', irk, 'convolution finished'
    return sc_p

def Compute_quasiparticles(bands, Ebnd, sigc, sigx, Vxct, omega, (iop_ac, iop_es, iop_gw0, npar_ac, iop_rcf), isp, fout, PRINT):
    (nirkp, nbnd, nom) = shape(sigc)
    eqp    = zeros(shape(bands))
    eqp_im = zeros(shape(bands))
    if (PRINT): print >> fout, 'Quasiparticle energies in eV'
    lwarn = True
    for irk in range(nirkp):
        for ie in range(nbnd):
            enk0 = Ebnd[irk,ie]
            vxc_nk = Vxct[irk,ie,ie].real
            enk = bands[irk,ie]
            # enk is the energy at which the self-energy is calculated
            debug = '.0.0' if (irk==0 and ie==0) else ''
            sig,dsig,apar = AnalyticContinuation(iop_ac, omega, sigc[irk,ie,:], npar_ac, enk, fout, iop_rcf, debug)
            # quasiparticle residue, but not at zero freqeucny, but at the energy of the band!
            z_nk = 1/(1-dsig.real)
            # self-energy at the energy of the band. This is what they believe is the best quasiparticle approximation
            s_nk = sig.real + sigx[irk,ie]
            if (z_nk > 1.0 or z_nk < 0):
                if (lwarn):
                    print >> fout, 'WARNING : Z_nk at energy e_k=', enk*H2eV, 'eV is unphysical', z_nk, 'irk=', irk, 'ie=', ie, 's_nk=', s_nk*H2eV, 'ds/dw=', dsig.real
                    lwarn = False  # we will warn only once
                z_nk = 1.0
                dsig = 0.0
            if iop_es == -1: # self-consistent GW0
                if iop_gw0==1: # this is default
                    delta = s_nk - vxc_nk + enk0 - enk  # this is used in self-consistent GW0
                elif iop_gw0==2:
                    delta = z_nk*(s_nk-vxc_nk) + enk0 - enk
                else:
                    delta = z_nk*(s_nk-vxc_nk + enk0 - enk)
            elif iop_es == 0: # this is used in G0W0
                delta = z_nk*(s_nk-vxc_nk)
            else:
                print >> fout, 'Not implemented here'
                sys.exit(1)
            eqp[irk,ie] = enk + delta
            eqp_im[irk,ie] = sig.imag*z_nk

            #print >> fout, 'iop_eps=%2d iop_gw0=%2d delta=%16.10f' % (iop_es, iop_gw0, delta)
            #print 'eqp['+str(irk)+','+str(ie)+']='+str(eqp[irk,ie]), 'and enk=', enk0
            if PRINT:
                print >> fout, 'eqp[irk=%3d,ie=%3d]=%16.10f and enk0=%16.10f enk=%16.10f znk=%10.8f snk=%16.10f vxc=%16.10f' % (irk,ie,eqp[irk,ie]*H2eV,enk0*H2eV,enk*H2eV,z_nk,s_nk.real*H2eV, vxc_nk.real*H2eV)
    return (eqp, eqp_im)
    
def Band_Analys(bande, EF, nbmax, titl, kqm, fout):
    bands = copy(bande[:,:nbmax])
    nirkp = shape(bands)[0]
    print >> fout, '-'*60
    print >> fout, '  '+titl+' Band Analysis'
    print >> fout, '-'*60
    print >> fout, '  Range of bands considered: %5d %5d' % (0,nbmax)
    print >> fout, '  EFermi[eV]= %10.4f' % (EF*H2eV,)
    
    if max(bands.ravel()) < EF or min(bands.ravel())>EF:
        print >> fout, 'WARNING from bandanaly:  - Fermi energy outside the energy range of bande!'
        print >> fout, 'minimal energy=', min(bands)*H2eV, 'max energy=', max(bands)*H2eV, 'EF=', EF*H2eV

    
    nocc_at_k = [len(filter(lambda x: x<EF, bands[ik,:])) for ik in range(nirkp)] # how many occuiped bands at each k-point
    nomax = max(nocc_at_k)-1                 # index of the last valence band
    numin = min(nocc_at_k)                   # index of the first conduction band

    ikvm = argmax(bands[:,nomax]) # maximum of the valence band
    ikcm = argmin(bands[:,numin]) # minimum of the conduction band

    Qmetal =  (nomax >= numin)
    if Qmetal:
        print >> fout, ' Valence and Conductance bands overlap: metallic!'
    if Qmetal:
        evbm = EF
    else:
        evbm = bands[ikvm,nomax]
    print >> fout, '  Band index for VBM and CBM=%4d %4d' % (nomax+1, numin+1)
    
    bands = (bands - evbm)*H2eV
        
    egap1 =     bands[ikcm,numin] - bands[ikvm,nomax]  # the smallest indirect gap between KS bands
    egap2 = min(bands[ikvm,numin:])-bands[ikvm,nomax]  # the direct gap starting from valence band
    egap3 = bands[ikcm,numin] - max(bands[ikcm,:(nomax+1)])# the direct gap starting from conduction band
    #print 'egap1=', egap1, 'egap=', egap2, 'egap3=', egap3
    
    if ikvm==ikcm: # direct gap
        print >> fout, ':BandGap_'+titl+' = %12.3feV' % egap1
        kp = kqm.kirlist[ikvm,:]/float(kqm.LCM)
        print >> fout, ('  Direct gap at k=  %8.3f%8.3f%8.3f') % tuple(kp), 'ik='+str(ikvm+1)
    else:
        print >> fout, (':BandGap_'+titl+' = %12.3f%12.3f%12.3f eV') % (egap1,egap2,egap3)
        kv = kqm.kirlist[ikvm,:]/float(kqm.LCM)
        kc = kqm.kirlist[ikcm,:]/float(kqm.LCM)
        print >> fout, '  Indirect gap, k(VBM)=%8.3f%8.3f%8.3f' % tuple(kv), 'ik='+str(ikvm+1)
        print >> fout, '                k(CBM)=%8.3f%8.3f%8.3f' % tuple(kc), 'ik='+str(ikcm+1)
    print >> fout, 'Range of each band with respect to VBM (eV):'
    print >> fout, ('%5s'+'%12s'*3) % ('n ','Bottom','Top','Width')
    for i in range(shape(bands)[1]):
        ebmin = min(bands[:,i])
        ebmax = max(bands[:,i])
        print >> fout, '%5d%12.3f%12.3f%12.3f' % (i+1, ebmin, ebmax, ebmax-ebmin)
    return (nomax,numin)



def padeMatrix(z,f,N,verbose=False):
    """
    Input variables:
    z       - complex. points in the complex plane.
    f       - complex. Corresponding (Green's) function in the complex plane. 
    N       - int. Number of Pade coefficients to use
    verbose - boolean. Determine if to print solution information
    Returns the obtained Pade coefficients.
    """
    # number of input points
    M = len(z)
    r = N/2
    y = f*z**r
    A = ones((M,N),dtype=complex)
    for i in range(M):
        A[i,:r] = z[i]**(arange(r))
        A[i,r:] = -f[i]*z[i]**(arange(r)) 
    # Calculated Pade coefficients
    # rcond=-1 means all singular values will be used in the solution.
    sol = linalg.lstsq(A, y, rcond=-1)
    # Pade coefficents
    x = sol[0]
    if verbose:
        print 'error_2= ',linalg.norm(dot(A,x)-y)
        print 'residuals = ', sol[1]
        print 'rank = ',sol[2]
        print 'singular values / highest_singlular_value= ',sol[3]/sol[3][0]
    return x

def epade(z,x):
    """
    Input variables:
    z - complex. Points where continuation is evaluated.
    x - complex. Pade approximant coefficient.
    
    Returns the value of the Pade approximant at the points z.
    """
    r = len(x)/2
    numerator = zeros(len(z),dtype=complex256)
    denomerator = zeros(len(z),dtype=complex256)
    for i in range(r):
        numerator += x[i]*z**i
        denomerator += x[r+i]*z**i
    denomerator += z**r
    return numerator/denomerator

def AnalyticContinuation(iop_ac,omg,sc,npar,enk, fout, iop_rcf=0.5, debug=''):
    # Analytic continuation of self-energy, and its evaluation at energy enk. Both the value and derivat
    # iop_ac==0  -- Thiele's reciprocal difference method as described in
    #       H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977). This is usual Pade, known in DMFT
    # iop_ac==1  -- Rojas, Godby and Needs (PRL 74, 1827 (1996). This is just a fit to imaginary axis data
    #
    #
    # We first take only a few values of self-energy and frequency, which will be used for Pade
    def f_unc(x, *argv):
        return fpade.pd_funval(x,argv)
    
    def f_jacob(x, *argv):
        return fpade.pd_jacoby(x,argv)
    #
    #if iop_ac==0 or npar==len(omg):
    if iop_ac==0 and abs(enk)/Ry2H < iop_rcf:
        # n = npar-1
        n = min(32, len(omg))
        z_p = omg[:n]*1j
        y_p = sc[:n]
        if True:
            apar = fpade.padecof(y_p, z_p)
            yc = fpade.padeg([enk], z_p, apar)
            dyc = (fpade.padeg([enk+1e-3], z_p, apar)-fpade.padeg([enk-1e-3], z_p, apar))/2e-3
        else:
            apar = padeMatrix(z_p,y_p, 2*len(z_p))
            yc  =  epade(array([enk]), apar)
            dyc = (epade(array([enk+1e-3]), apar) - epade(array([enk-1e-3]), apar))/2e-3
            
        #print >> fout, 'classic Pade: s['+str(enk*H2eV)+']='+str(yc[0]*H2eV)
    elif iop_ac==1 or (iop_ac==0 and abs(enk)/Ry2H > iop_rcf):
        ## Here we fit all values of self-energy at Matsubara pointz iw==z to the following rational function
        #          (c[1] + c[2] z)
        #  f(z) = ----------------------
        #          1 + c[3] z + c[4] z^2
        #  which is a rational function with two poles.
        #  It is also equivalent to pade of the level 4, i.e.,
        #  f(z) = a[1]/(1+a[2]z/(1+a[3]z/(1+a[4]z))), where non-linear relation between a[i] and c[i] exists
        #  
        iwpa = arange(npar)*int((len(omg)-1.)/(npar-1.))  # linear mesh of only a few points, equidistant.
        iwpa[-1] = len(omg)-1                             # always take the last point
        cx = zeros(len(iwpa), dtype=complex)  # imaginary frequency at selected points
        cy = zeros(len(iwpa), dtype=complex)  # self-energy at selected points
        for ip,iw in enumerate(iwpa):
            cx[ip] = omg[iw]*1j
            cy[ip] = sc[iw]
        apar = fpade.init_c(cx, cy)           # just a way to find a good starting point for minimization
        n = len(apar)                         # how many complex coefficints we want to fit
        anl = hstack( (apar.real,apar.imag) ) # we stack these complex coefficients (starting point coefficients) into real array
        if True:
            x = hstack( (-omg.imag, omg.real) )  # take all frequency as x values : x = i*omg
            y = hstack( ( sc.real, sc.imag) )    # take all self-energy points as y values: y = [real,imag]
            # now fit all self-energy points to rational function
            #   sigma(iom) = P(iom)/Q(iom), where
            #   P(iom) = sum_{k=0,n  } c_{k}   (iom)^k
            #   Q(iom) = sum_{k=1,n+1} c_{k+n} (iom)^k
            # and determine c_k by minimization of chi2 using Levenberg-Marquardt algorithm.
            popt, pcov = optimize.curve_fit( f_unc, x, y, p0=anl, jac=f_jacob, method='lm')
            apar = popt[:n] + popt[n:]*1j # these are ck coefficients
            # Evaluate self-energy at energy of the band ek
            yc,dyc = fpade.pd_funvalc([enk],apar)
        else:
            chisq=0.0
            chisq = fpade.nllsq(omg, sc, anl) 
            apar = anl[:n] + anl[n:]*1j
            yc,dyc = fpade.pd_funvalc([enk],apar)
        #print >> fout, 'modified Pade: s['+str(enk*H2eV)+']='+str(yc[0]*H2eV)
    else:
        s0 = sc[0].real
        dsdw = polyfit([0,omg[0],omg[1]], [0.0,sc[0].imag,sc[1].imag], 1)[0]
        yc = [s0 + dsdw*enk]
        dyc = [dsdw]
        apar = [s0,dsdw]
        #print >> fout, 'Simple quasiparticle approximation: s['+str(enk*H2eV)+']='+str(yc[0]*H2eV)
        
    #print >> fout, 'sig[e=%10.6f]=%10.6f %10.6f' % (enk,yc.real, yc.imag)
    #print ' cx, cy, apar'
    #for ip in range(len(iwpa)):
    #    print '%21.16f%21.16f  %20.16f%20.16f  %20.16f%20.16f' % (cx[ip].real, cx[ip].imag, cy[ip].real, cy[ip].imag, apar[ip].real, apar[ip].imag)
    if debug:
        romega = hstack( (-omg[::-1],omg) )
        if iop_ac==0 and abs(enk)/Ry2H < iop_rcf:
            #print 'iop_ac=', iop_ac, 'abs(enk)/Ry2H=', abs(enk)/Ry2H, 'iop_rcf=', iop_rcf
            if True:
                yre = fpade.padeg(romega, z_p, apar)
                yim = fpade.padeg(romega*1j, z_p, apar)
            else:
                yre = epade(romega, apar)
                yim = epade(romega*1j, apar)
        elif iop_ac==1 or (iop_ac==0 and abs(enk)/Ry2H > iop_rcf):
            yre,dyre = fpade.pd_funvalc(romega, apar)
            x = hstack( (-romega.imag, romega.real) )
            _yim_ = f_unc(x, *popt)
            yim = _yim_[:len(romega)] + _yim_[len(romega):]*1j
        else:
            yre = s0 + dsdw*romega
            yim = s0 + dsdw*romega*1j
        fo = open('sigma_ancont'+debug, 'w')
        print >> fo, '# enk=', enk*H2eV, 'and result is', yc[0]*H2eV
        for i in range(len(romega)):
            print >> fo, romega[i]*H2eV, yre[i].real*H2eV, yre[i].imag*H2eV, yim[i].real*H2eV, yim[i].imag*H2eV
        fo.close()
        fo = open('sigma_data'+debug, 'w')
        for i in range(len(omg)):
            print >> fo, omg[i]*H2eV, sc[i].real*H2eV, sc[i].imag*H2eV
        fo.close()
    return yc[0], dyc[0], apar
    
def mpiSplitArray(mrank,msize,leng):
    def SplitArray(irank,msize,leng):
        if leng % msize==0:
            pr_proc = int(leng/msize)
        else:
            pr_proc = int(leng/msize+1)
        if (msize<=leng):
            iqs,iqe = min(irank*pr_proc,leng) , min((irank+1)*pr_proc,leng)
        else:
            rstep=(msize+1)/leng
            if irank%rstep==0 and irank/rstep<leng:
                iqs = irank/rstep
                iqe = iqs+1
            else:
                if irank/rstep<leng:
                    iqs = irank/rstep
                    iqe = irank/rstep
                else:
                    iqs = leng-1
                    iqe = leng-1
        return iqs,iqe
    #print 'mrank=', mrank, 'msize=', msize, 'leng=', leng
    sendcounts=[]
    displacements=[]
    for irank in range(msize):
        iqs,iqe = SplitArray(irank,msize,leng)
        sendcounts.append((iqe-iqs))
        displacements.append(iqs)
    iqs,iqe = SplitArray(mrank,msize,leng)
    return iqs,iqe, array(sendcounts,dtype=int), array(displacements,dtype=int)


