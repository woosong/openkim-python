#
#  KIM Calculator for ASE
#  
#  We make use of the SWIG python KIM interface to make use of KIM models  
#  through ASE
#  
#  Authors: Matt Bierbaum
#           YJ Chen
#
#############################################################################
#
#  Notes for work in progress:
#      1. name pointers to KIM with "km_" 
#      2. need info on neighbor list
#      3. for certain models need particleSize and particleCharge 
#  

from kimservice import *
import neighborlist as kimnl


class KIM_Calculator:

    def __init__(self,modelname):
        # correctly set up the .kim file so that is matches that
        # of modelname, i.e. will run any types
        self.tname = "test_ase_kim"
        self.testf = open("test_ase_kim.kim", "rt")
        self.teststring = "".join(self.testf.readlines())

        self.modelname = modelname
        status, self.pkim = KIM_API_init_str(self.teststring, self.modelname)
        if KIM_STATUS_OK > status:
            KIM_API_report_error('KIM_API_init',status)
        
        # We need to create a string that matches the capabilities
        # of the test that we wish to use.  We don't need to grab 
        # anything, really, but accomodate anything.  
        #status, self.pmdl = KIM_API_model_info(self.modelname)        

        self.km_numberOfAtoms  = None
        self.km_energy         = None
        self.km_forces         = None
        self.km_particleEnergy = None
        self.km_virial         = None
        self.km_particleVirial = None
        self.km_hessian        = None

    def update(self,atoms):
        """
        here we connect the KIM pointers to values in the ase atoms class
        """
        natoms = atoms.get_number_of_atoms()
        ntypes = len(set(atoms.get_atomic_numbers()))
        
        if self.km_numberOfAtoms != natoms:
            KIM_API_free(self.pkim)
            status, self.pkim = KIM_API_init_str(self.teststring, self.modelname)
            KIM_API_allocate(self.pkim, natoms, ntypes)
            kimnl.initialize(self.pkim)
            KIM_API_model_init(self.pkim)
 
        # we proceed to set the values of the standard things each model and atom class has
        self.km_numberOfAtoms      = KIM_API_get_data_ulonglong(self.pkim, "numberOfParticles")
        self.km_numberAtomTypes    = KIM_API_get_data_int(self.pkim, "numberParticleTypes")
        self.km_atomTypes          = KIM_API_get_data_int(self.pkim, "particleTypes")
        self.km_coordinates        = KIM_API_get_data_double(self.pkim, "coordinates")
        self.km_numberOfAtoms[0]   = natoms 
        self.km_numberAtomTypes[0] = ntypes
        self.km_coordinates[:]     = atoms.get_positions().flatten()
        
        # fill the proper chemical identifiers 
        symbols = atoms.get_chemical_symbols()
        for i in range(natoms):
            self.km_atomTypes[i] = KIM_API_get_partcl_type_code(self.pkim, symbols[i])
        
        # build the neighborlist (not a cell-based, type depends on model)
        #kimnl.build_neighborlist_allall(self.pkim)

        # check if model calculates energies, forces, virials
        if KIM_API_get_compute(self.pkim,"energy"):
            self.km_energy = KIM_API_get_data_double(self.pkim, "energy")
        if KIM_API_get_compute(self.pkim,"forces"):
            self.km_forces = KIM_API_get_data_double(self.pkim, "forces")
        if KIM_API_get_compute(self.pkim, "particleEnergy"):
            self.km_particleEnergy = KIM_API_get_data_double(self.pkim, "particleEnergy")
        if KIM_API_get_compute(self.pkim, "virial"):
            self.km_virial = KIM_API_get_data_double(self.pkim, "virial")
        """
        if KIM_API_get_compute(self.pkim, "particleVirial"):
            self.km_particleVirial = KIM_API_get_data_double(self.pkim, "particleVirial")
        if KIM_API_get_compute(self.pkim, "hessian"):
            self.km_hessian = KIM_API_get_data_double(self.pkim, "hessian")
        """

        KIM_API_model_compute(self.pkim)

        
    def get_potential_energy(self,atoms):
        self.update(atoms)
        if self.km_energy is not None:
            return self.km_energy
        else:
            raise SupportError("energy") 
    
    def get_potential_energies(self,atoms):
        self.update(atoms)
        if self.km_particleEnergy is not None:
            particleEnergies = self.km_particleEnergy.reshape((self.km_numberOfAtoms,3))
            return particleEnergies
        else:
            raise SupportError("potential energies") 
     
    def get_forces(self,atoms):
        self.update(atoms)
        if self.km_forces is not None:
            forces = self.km_forces.reshape((self.km_numberOfAtoms,3))
            return forces
        else:
            raise SupportError("forces") 
    
    def get_stress(self,atoms):
        self.update(atoms)
        if self.km_virial is not None:
            return self.km_virial
        else:
            raise SupportError("stress")

    def get_stresses(self,atoms):
        self.update(atoms)
        if self.km_particleVirial is not None:
            return self.km_particleVirial
        else:
            raise SupportError("stress per particle")

    def get_hessian(self,atoms):
        self.update(atoms) 
        if self.km_hessian is not None:
            return self.km_hessian
        else:
            raise SupportError("hessian") 


class SupportError(Exception):

    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)+" computation not supported by model"
