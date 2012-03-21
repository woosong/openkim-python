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
#      3. should write a calculation_required function
#       

from kimservice import *
import neighborlist as kimnl
from kimscriptmaker import makekimscript, checkIndex
import numpy
import time

class KIM_Calculator:

    def __init__(self,modelname,ucell=True):
        # correctly set up the .kim file so that is matches that
        # of modelname, i.e. will run any types
        self.testname = "test_"+modelname
        self.teststring = None

        self.modelname = modelname 
        status, km_pmdl = KIM_API_model_info(modelname)
        if KIM_STATUS_OK > status:
            KIM_API_report_error('KIM_API_model_info',status)
            raise InitializationError(self.modelname)
        km_pmdl = 0

        self.ucell = ucell

        # initialize pointers for kim
        self.km_numberOfAtoms  = None
        self.km_particleCharge = None
        self.km_energy         = None
        self.km_forces         = None
        self.km_particleEnergy = None
        self.km_virial         = None
        self.km_particleVirial = None
        self.km_hessian        = None
        
        # initialize ase atoms specifications 
        self.pbc = None
        self.cell = None
        self.cell_orthogonal = None

    def makeTestString(self,atoms):
        """
        makes string if it doesn't exist, if exists just keeps it as is
        """
        if self.teststring is None or self.cell_BC_changed(atoms):
            self.teststring = makekimscript(self.modelname,self.testname,atoms)


    def cell_BC_changed(self,atoms):
        """
        this function is to check whether BC has changed and cell orthogonality has changed
        because we might want to change neighbor list generator method
        """
        cell_orthogonal = ((abs(numpy.dot(atoms.get_cell()[0],atoms.get_cell()[1])) + abs(numpy.dot(atoms.get_cell()[0],atoms.get_cell()[2])) + abs(numpy.dot(atoms.get_cell()[1],atoms.get_cell()[2])))<10**(-8))
        if (self.pbc != atoms.get_pbc()).any() or self.cell_orthogonal != cell_orthogonal:
            return True
        else:
            return False

    def calculation_required2(self,atoms):
        """
        this checks whether or not the atoms configuration has changed and we need to recalculate..
        """
        return (self.km_energy is None or \
                (self.km_numberOfAtoms != atoms.get_number_of_atoms()) or \
                (self.km_atomTypes[:] != atoms.get_atomic_numbers()).any() or \
                (self.km_coordinates[:] != atoms.get_positions().flatten()).any() or \
                (self.pbc != atoms.get_pbc()).any() or \
                (self.cell != atoms.get_cell()).any())


    def update(self,atoms):
        """
        here we connect the KIM pointers to values in the ase atoms class
        """

        # here we only reinitialize the model if the number of Atoms / types of atoms have changed, or if the model is uninitialized
        natoms = atoms.get_number_of_atoms()
        ntypes = len(set(atoms.get_atomic_numbers()))

        if self.km_numberOfAtoms != natoms or self.km_numberAtomTypes != ntypes or cell_BC_changed:
            self.makeTestString(atoms)
            status, self.pkim = KIM_API_init_str(self.teststring, self.modelname)
            if KIM_STATUS_OK > status:
                KIM_API_report_error('KIM_API_init',status)
                raise InitializationError(self.modelname)
        
            KIM_API_allocate(self.pkim, natoms, ntypes)
            kimnl.initialize(self.pkim)
            KIM_API_model_init(self.pkim)
    
            # get pointers to model inputs
            self.km_numberOfAtoms      = KIM_API_get_data_ulonglong(self.pkim, "numberOfParticles")
            self.km_numberAtomTypes    = KIM_API_get_data_int(self.pkim, "numberParticleTypes")
            self.km_atomTypes          = KIM_API_get_data_int(self.pkim, "particleTypes")
            self.km_coordinates        = KIM_API_get_data_double(self.pkim, "coordinates")
            if checkIndex(self.pkim,"particleCharge") >= 0:
                self.km_particleCharge = KIM_API_get_data_double(self.pkim,"particleCharge")
            if checkIndex(self.pkim, "particleSize") >= 0:
                self.km_particleSize = KIM_API_get_data_double(self.pkim,"particleSize")
                print "WARNING: ASE atoms do not have particle sizes"
    
            # check what the model calculates and get model outputs
            if checkIndex(self.pkim,"energy") >= 0:
                self.km_energy = KIM_API_get_data_double(self.pkim, "energy")
            if checkIndex(self.pkim,"forces") >= 0:
                self.km_forces = KIM_API_get_data_double(self.pkim, "forces")
            if checkIndex(self.pkim, "particleEnergy") >= 0:
                self.km_particleEnergy = KIM_API_get_data_double(self.pkim, "particleEnergy")
            if checkIndex(self.pkim, "virial") >=0:
                self.km_virial = KIM_API_get_data_double(self.pkim, "virial")
            if checkIndex(self.pkim, "particleVirial") >=0:
                self.km_particleVirial = KIM_API_get_data_double(self.pkim, "particleVirial")
            if checkIndex(self.pkim, "hessian")>=0:
                self.km_hessian = KIM_API_get_data_double(self.pkim, "hessian")
                        
        if self.calculation_required2(atoms):
            # if the calculation is required we proceed to set the values of the standard things each model and atom class has
            self.km_numberOfAtoms[0]   = natoms 
            self.km_numberAtomTypes[0] = ntypes
            self.km_coordinates[:]     = atoms.get_positions().flatten()
            if self.km_particleCharge is not None:
                km_particleCharge[:] = atoms.get_charges()       
 
            # fill the proper chemical identifiers 
            symbols = atoms.get_chemical_symbols()
            for i in range(natoms):
                self.km_atomTypes[i] = KIM_API_get_partcl_type_code(self.pkim, symbols[i])
            
            # build the neighborlist (not a cell-based, type depends on model)
            #print KIM_API_get_NBC_method(self.pkim)
            kimnl.set_cell(atoms.get_cell().flatten(), atoms.get_pbc().flatten().astype('int8'))
            neigh_start = time.time()
            kimnl.build_neighborlist(self.pkim)
            #print "neigh time = ", time.time() - neigh_start

            calc_start = time.time()
            KIM_API_model_compute(self.pkim)
            #print "calc time = ", time.time() - calc_start 

            # set the ase atoms stuff to current configuration
            self.pbc = atoms.get_pbc()
            self.cell = atoms.get_cell()
            self.cell_orthogonal = ((abs(numpy.dot(atoms.get_cell()[0],atoms.get_cell()[1])) + abs(numpy.dot(atoms.get_cell()[0],atoms.get_cell()[2])) + abs(numpy.dot(atoms.get_cell()[1],atoms.get_cell()[2])))<10**(-8))

        
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

    def __del__(self):
        """ 
        Garbage collects the KIM API objects automatically
        """
        pass
        #KIM_API_free(self.pkim)

class SupportError(Exception):

    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)+" computation not supported by model"

class InitializationError(Exception):
            
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)+" initialization failed"
