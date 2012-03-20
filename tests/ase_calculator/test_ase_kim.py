from kimservice import *
import ase
from ase.lattice.cubic import FaceCenteredCubic, SimpleCubic 
from ase.visualize import view
from kimcalculator import *
from numpy import *
from ase.calculators.emt import EMT

N = 12 
ar = FaceCenteredCubic('Ar', pbc=[(1,1,1)], directions=[[1,0,0],[2,1,0],[1,1,1]], size=[N,N,N])
print ar.get_cell()
#view(ar) 

calc1 = KIM_Calculator("ex_model_Ar_P_LJ", True)
ar.set_calculator(calc1)
print "energy = ", ar.get_potential_energy()

calc2 = KIM_Calculator("ex_model_Ar_P_LJ", False)
ar.set_calculator(calc2)
print "energy = ", ar.get_potential_energy()


