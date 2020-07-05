def prepare_Wannier(case, strc, in1, latgen, kqm, ks, radf, pw, rmax, fout):
    isp=0
    ndk = array(kqm.ndiv,dtype=int)
    rlen0=[]
    rind0=[]
    for ii,ir in enumerate(itertools.product(arange(0,ndk[0]),arange(0,ndk[1]),arange(0,ndk[2]))):
        rvec = dot(ir,latgen.rbas)
        rlen0.append( linalg.norm(rvec) )
        rind0.append(ir)
    indx = argsort(rlen0) # kind='stable')  # just obtaining index to the sorted sequence
    # rearange arrays so that they are sorted
    rlen = zeros(shape(rlen0))  # |R|
    rind = zeros(shape(rind0))  # \vR in integer lattice coordinates
    for i0,i in enumerate(indx):
        rlen[i0]   = rlen0[i]
        rind[i0,:] = rind0[i]
        #print i0, rlen[i0], rind[i0,:]/float(LCM)

    #rmax *= 0.5
    rvec = zeros(3)
    for i in range(3):
        rvec[i] = linalg.norm(latgen.rbas[i,:])
    nr = array(map(int,1./rvec * rmax))*2
    
    #print 'nr=', nr, 'rmax=', rmax
    rlen0=[]
    rind0=[]
    for ii,ir in enumerate(itertools.product(range(-nr[0],nr[0]+1),range(-nr[1],nr[1]+1),range(-nr[2],nr[2]+1))):
        rvec = dot(ir,latgen.rbas)
        rr = linalg.norm(rvec)
        if (rr <= rmax):
            rlen0.append(rr)
            rind0.append(ir)
    indx = argsort(rlen0) # kind='stable')  # just obtaining index to the sorted sequence
    # rearange arrays so that they are sorted
    rlen2 = zeros(shape(rlen0))           # |R|
    rind2 = zeros(shape(rind0), dtype=int)           # \vR in cartesian
    for i0,i in enumerate(indx):
        rlen2[i0]   = rlen0[i]
        rind2[i0,:] = rind0[i]
    
    


        
    #r2icartes = linalg.inv(kqm.k2icartes)
    #if latgen.ortho or latgen.lattice[1:3]=='CXZ': # R to cartesian coordinates
    #    for ir in range(len(rind)):
    #        rind[ir,:] = dot(r2icartes, rind[ir,:] )

    klist = mcmn.cart2int(array(kqm.klist)/float(kqm.LCM),strc,latgen)

    
    for ik in range(1,len(klist)):
        kil = array(klist[ik,:]) #/float(kqm.LCM)
        kr = 2*pi*dot(rind,kil)
        expkr = sum(exp(1j*kr))
        if abs(expkr) > 1e-13:
            print 'WARNING \sum_R exp(i*k*R) is not 0 for k!=Gamma'
            print kil, expkr
            
    bq = zeros(6,dtype=int)
    for iq,q in enumerate(kqm.qlist):
        if sum(abs(q))==1:
            if q[0]==1 and q[1]==0 and q[2]==0:
                bq[0] = iq
            elif q[0]==0 and q[1]==1 and q[2]==0:
                bq[2] = iq
            elif q[0]==0 and q[1]==0 and q[2]==1:
                bq[4] = iq
        elif q[0]==(kqm.LCMq-1) and q[1]==0 and q[2]==0:
            bq[1] = iq
        elif q[0]==0 and q[1]==(kqm.LCMq-1) and q[2]==0:
            bq[3] = iq
        elif q[0]==0 and q[1]==0 and q[2]==(kqm.LCMq-1):
            bq[5] = iq
    print bq
    t_times = zeros(3)
    
    iateq_ind=[]
    for iat in range(len(strc.mult)):
        iateq_ind += [iat]*strc.mult[iat]
    
    nst, nend = ks.ibgw, ks.nbgw
    nend=5
    print 'nst=', nst, 'nend=', nend
        
    A = zeros((len(klist),nend-nst,nend-nst),dtype=complex)
    t0 = timer()
    for ik in range(len(kqm.klist)):
        #print 'nst=', nst, 'nend=', nend
        #(nst, nend, mst, mend, cst, cend) = band_limits
        #(nst, nend, mst, mend, cst, cend) = (ks.ibgw, ks.nbgw, 0, nomx+1, 0, ks.ncg_x)
        
        irk = kqm.kii_ind[ik]
        
        # First create alpha, beta, gamma for k
        Aeigk = array(ks.all_As[irk][nst:nend,:], dtype=complex)   # eigenvector from vector file
        if kqm.k_ind[irk] != ik:                       # not irreducible
            Aeigk *= exp(-2*pi*1j * ks.phase_arg[ik][:])  # adding phase : zzk[1:ngk,ib] = phase[1:ngk] * zzk[1:ngk,ib]
        kil = array(kqm.klist[ik,:])/float(kqm.LCM)  # k in semi-cartesian form
        alm,blm,clm = lapwc.set_lapwcoef(kil, pw.gindex, ks.indgk[ik], radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.Vol, 1, True, kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax, ks.nv[irk])
        
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
        if ik==0:
            how_big = [{} for n in range(nend-nst)]
            for idf in range(ndf):
                for ilm in range(ntnt):
                    for n in range(nend-nst):
                        how_big[n][(0,idf,ilm)] = abs(vdot(alfa[n,ilm,idf],alfa[n,ilm,idf]).real)
                        how_big[n][(1,idf,ilm)] = abs(vdot(beta[n,ilm,idf],alfa[n,ilm,idf]).real)
                        lc = int(sqrt(ilm))
                        iat = iateq_ind[idf]
                        for ilo in range(in1.nLO_at[lc,iat]):
                            how_big[n][(2,idf,ilm,ilo)] = abs(vdot(gama[n,ilo,ilm,idf],gama[n,ilo,ilm,idf]).real)
            orbs=[]
            for n in range(nend-nst):
                # for each band we want to find a localized orbital that has the largest overlap.
                # if such orbital was already selected because it has large overlap with some other band, we go down the list and accept orbital with the second largest overlap
                key_how_big = sorted(how_big[n], key=how_big[n].get, reverse=True)
                which_orbital = key_how_big[:nend-nst]
                
                for i in range(len(which_orbital)):
                    if which_orbital[i] not in orbs:
                        orbs.append( which_orbital[i] )
                        break
                for w in which_orbital:
                    print 'ibnd=', n, 'orbital=', w, 'v=', how_big[n][w]
            print 'selected orbitals : ', orbs
    
        for i,w in enumerate(orbs):
            idf, ilm = w[1:3]
            if w[0]==0:
                A[ik,:,i] = conj(alfa[:,ilm,idf])
            elif w[0]==1:
                A[ik,:,i] = conj(beta[:,ilm,idf])
            elif w[0]==2:
                ilo = w[3]
                A[ik,:,i] = conj(gama[:,ilo,ilm,idf])
    t1 = timer()
    t_times[0] += t1-t0

    print 'shape(Ebnd)=', shape(ks.Ebnd)
    Hk = zeros(shape(A), dtype=complex)
    for ik in range(len(klist)):
        u, s, vh = linalg.svd(A[ik,:,:])
        print ik, 'singular-values:', s.tolist()
        Uk = dot(u,vh)
        irk = kqm.kii_ind[ik]
        enk0 = ks.Ebnd[isp,irk,nst:nend]
        Hk[ik,:,:] = dot( conj(Uk.T), dot(diag(enk0), Uk) )

        #print ik, 'min=', min(enk0.ravel()), 'max=', max(enk0.ravel())
        
        #chk = linalg.eigvalsh(Hk[ik,:,:])
        #print 'chk=', chk.tolist(), enk0



    c1,c2 = 0.25,0.25      # The two coefficients mentioned in the paper as C1 and C2
    #rho = zeros(len(rind2))  # We use here rho = 1 - 2*c1*R^2 + c1^2*R^4 + c2*R^6
    x2 = (rlen2/rlen2[1])**2
    rho = (1-c1*x2)**2 + c2*x2**3
    
    smat1 = zeros((len(klist),len(rind2)),dtype=complex)
    for ik in range(len(klist)):
        smat1[ik,:] = exp( 2*pi*1j*dot(rind2,klist[ik]) )
    
    sm2 = zeros((len(klist)-1,len(rind2)),dtype=complex) # sm2 <- Sm[k_i]-Sm[k_n] in the paper
    nb = nend-nst                            # number of input energies.
    dele = zeros((len(klist)-1,nb,nb),dtype=complex)
    
    for ik in range(len(klist)-1):
        sm2[ik,:]  = smat1[ik+1,:]  - smat1[0,:]   #  sm2[ik,istar] = Sm[k_i]-Sm[k_0]
        dele[ik,:,:] = Hk[ik+1,:,:] - Hk[0,:,:]    #  dele <- e[k_j]-e[k_0] in the paper
    
    h = zeros((len(klist)-1,len(klist)-1),dtype=complex)          # H_ij in the paper
    for ik in range(len(klist)-1):
        for jk in range(len(klist)-1):
            h[ik,jk] += sum(sm2[ik,:]*conj(sm2[jk,:])/rho[:])
    
    Hinv = linalg.inv(h)
    Hf = dot(conj(sm2.T), Hinv)
    for ist in range(len(rho)):
        Hf[ist,:] *= 1/rho[ist]
    coef = tensordot( Hf, dele, axes=1 )
    
    coef[0,:,:] = Hk[0,:,:] - tensordot(smat1[0,:],coef,axes=1) # epsilon_m in the paper

    
    band_energy_file = case+'.energy_band'
    (klist2, wegh, ebnd2, hsrws, knames) = w2k.Read_energy_file(band_energy_file, strc, fout, give_kname=True)
    kvecs2 = mcmn.cart2int(klist2,strc,latgen)

    smat2 = zeros((len(kvecs2),len(rind2)),dtype=complex)
    for ik in range(len(kvecs2)):
        smat2[ik,:] = exp( 2*pi*1j*dot(rind2,kvecs2[ik]) )
    
    Hk2 = tensordot(smat2, coef, axes=1)     # this is the resulting energy on the dense grid : e_m * S_m(k) in the paper

    enks = zeros((len(kvecs2),nend-nst))
    for ik,k in enumerate(kvecs2):
        enk = linalg.eigvalsh(Hk2[ik])
        enks[ik,:] = enk
        print ik, 'k=', k, enks[ik,:].tolist()
    
    e2bnd = zeros((len(ebnd2),4))
    for ik in range(len(ebnd2)):
        e2bnd[ik,:] = ebnd2[ik][nst:nst+4]
    e2bnd -= 0.0079145774

    from pylab import *
    plot(enks[:,0]*2)
    plot(enks[:,1]*2)
    plot(enks[:,2]*2)
    plot(enks[:,3]*2)
    plot(e2bnd[:,0], ':')
    plot(e2bnd[:,1], ':')
    plot(e2bnd[:,2], ':')
    plot(e2bnd[:,3], ':')
    
    show()
    
    sys.exit(0)

    

        

    
    Hr = zeros((len(rind),nend-nst,nend-nst),dtype=complex)
    nkp = len(klist)
    for ir in range(len(rind)):
        kr = dot(klist,rind[ir]) # need to add pos*k
        expkr = exp(-1j*2*pi*kr)
        Hr[ir,:,:] = dot(Hk.T,expkr).T/nkp

    
    
    #sys.exit(0)
    for ik,k in enumerate(klist):
        kr = dot(rind,k)/float(kqm.LCM) # need to add pos*k
        expkr = exp(1j*2*pi*kr)
        Hk = dot(Hr.T,expkr).T
        chk = linalg.eigvalsh(Hk)
        irk = kqm.kii_ind[ik]
        diff = chk - ks.Ebnd[isp,irk,nst:nend]
        print 'chk[',k/float(kqm.LCM),']=', max(abs(diff)), chk.tolist(), ks.Ebnd[isp,irk,nst:nend].tolist()
        
    
    #klist1 = kqm.kirlist/float(kqm.LCM)
    #kvecs2 = klist2
    enks = zeros((len(kvecs2),nend-nst))
    for ik,k in enumerate(kvecs2):
        kr = dot(rind,k)
        expkr = exp(1j*2*pi*kr)
        Hk = dot(Hr.T,expkr).T
        enk = linalg.eigvalsh(Hk)
        enks[ik,:] = enk
        print ik, 'k=', k, enks[ik,:].tolist()
        
        
        #jk = kqm.kqid[ik,iq]  # index of k+q
        #jrk = kqm.kii_ind[jk]
        #
        ## And next create alpha, beta, gamma for k+q
        #Aeigq = array( conj( ks.all_As[jrk] ), dtype=complex)  # eigenvector from vector file
        #if kqm.k_ind[jrk] != jk:                             # the k-q-point is reducible, eigenvector needs additional phase
        #    Aeigq *= exp( 2*pi*1j * ks.phase_arg[jk][:] )
        #kjl = array(kqm.klist[jk,:])/float(kqm.LCM)          # k+q in semi-cartesian form
        #alm,blm,clm = lapwc.set_lapwcoef(kjl, pw.gindex, ks.indgk[jk], radf.abcelo[isp], strc.rmt, strc.vpos, strc.mult, radf.umt[isp], strc.rotloc, latgen.trotij, latgen.Vol, 2, True, kqm.k2cartes, in1.nLO_at, in1.nlo, in1.lapw, in1.nlomax, ks.nv[jrk])
        #(ngj,nLOmax,ntnt,ndf) = shape(clm)
        #(ngj,ntnt,ndf) = shape(alm)
        #(nbmax,ngj2) = shape(Aeigq)
        ## And now change alm,blm,clm to band basis, which we call alfa,beta,gama
        #alfp = reshape( la.matmul(Aeigq, reshape(alm, (ngj,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
        #betp = reshape( la.matmul(Aeigq, reshape(blm, (ngj,ntnt*ndf)) ), (nbmax,ntnt,ndf) )
        #if in1.nlomax > 0:
        #    gamp = reshape( la.matmul(Aeigq, reshape(clm, (ngj,ntnt*ndf*nLOmax)) ), (nbmax,nLOmax,ntnt,ndf) )
        #else:
        #    gamp = zeros((1,1,1,1),dtype=complex,order='F') # can not have zero size
        #
        #abc_lapw = (alfa,beta,gama,alfp,betp,gamp)

    e2bnd = zeros((len(ebnd2),4))
    for ik in range(len(ebnd2)):
        e2bnd[ik,:] = ebnd2[ik][nst:nst+4]
    e2bnd -= 0.0079145774
    plot(enks[:,0]*2)
    plot(enks[:,1]*2)
    plot(enks[:,2]*2)
    plot(enks[:,3]*2)
    plot(e2bnd[:,0], ':')
    plot(e2bnd[:,1], ':')
    plot(e2bnd[:,2], ':')
    plot(e2bnd[:,3], ':')
    
    show()
    
    print 'time=', t_times
