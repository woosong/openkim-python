from kimservice import *
from ase.lattice.cubic import FaceCenteredCubic 
from kimcalculator import *

N = 2

ar = FaceCenteredCubic('Ar', directions=[[1,0,0],[0,1,0],[1,1,1]], size=[N,N,N])
print "Made crystal"
ar.set_calculator(KIM_Calculator("ex_model_Ar_P_LJ"))

print "Calculating energy: "
print ar.get_potential_energy()
#print ar.get_forces()
