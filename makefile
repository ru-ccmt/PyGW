CMP = /usr/bin/f2py --opt='-fast' --link-lapack  --fcompiler=intelem
CMPM = /usr/bin/f2py --opt='-fast' --link-lapack  --fcompiler=intelem
FORT = ifort
FFLAGS = -fast

#CMP = f2py --opt='-O3 -ffree-line-length-none'  --fcompiler=gnu95
#CMPM = f2py --opt='-O3 -ffree-line-length-none -fopenmp'  --fcompiler=gnu95
#FORT = gfortran
#FFLAGS = -O2 -ffree-line-length-none -fopenmp

modules = fnc.so cum_simps.so gaunt.so rdVec.so radials.so radd.so for_kpts.so for_tetrahedra.so for_vxcnn.so lapwcoeff.so for_Coulomb.so for_q_0.so for_pade.so
objects = sphbes.o for_tetrahedra2.o for_tetrahedra3.o for_pade2.o

all : $(objects) $(modules) 

lapwcoeff.so : lapwcoeff.f90 sphbes.o
	$(CMP) -c $< sphbes.o -m lapwcoeff

for_Coulomb.so : for_Coulomb.f90 sphbes.o
	$(CMP) -c $< sphbes.o -m for_Coulomb

for_tetrahedra.so : for_tetrahedra.f90 for_tetrahedra2.o for_tetrahedra3.o
	$(CMPM) -c $< for_tetrahedra2.o for_tetrahedra3.o -m for_tetrahedra

for_q_0.so : for_q_0.f90
	$(CMPM) -c for_q_0.f90 -m for_q_0 

for_vxcnn.so : for_vxcnn.f90
	$(CMPM) -c for_vxcnn.f90 -m for_vxcnn

for_pade.so : for_pade.f90 for_pade2.o
	$(CMP) -c $< for_pade2.o -m for_pade

clean :
	rm -f $(modules) $(objects) *.pyc *.mod *.o

%.so : %.f90
	$(CMP) -c $< -m $(*F)

%.o : %.f90
	$(FORT) -fPIC -c $(FFLAGS) $<
