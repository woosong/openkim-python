from kimservice import *
import ase
from ase.lattice.cubic import FaceCenteredCubic, SimpleCubic 
from ase.visualize import view
from kimcalculator import *
from numpy import *
from ase.calculators.emt import EMT

N = 20 
ar = FaceCenteredCubic('Ar', pbc=[(1,0,1)], directions=[[1,0,0],[0,1,0],[1,1,1]], size=[N,N,N/2])
calc1 = KIM_Calculator("ex_model_Ar_P_LJ")
ar.set_calculator(calc1)
print "energy = ", ar.get_potential_energy()

