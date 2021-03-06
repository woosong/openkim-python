from kimservice import *
import ase
from ase.lattice.cubic import FaceCenteredCubic, SimpleCubic 
from ase.visualize import view
from kimcalculator import *
from numpy import *
from ase.calculators.emt import EMT

N = 10 
ar = FaceCenteredCubic('Ar', pbc=[(1,1,1)], directions=[[1,0,0],[0,1,0],[1,1,1]], size=[N,N,N])

calc1 = KIMCalculator("ex_model_Ar_P_LJ")
ar.set_calculator(calc1)
print "energy = ", ar.get_potential_energy()

