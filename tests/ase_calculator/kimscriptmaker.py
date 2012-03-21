from kimservice import *
import numpy

#
# KIM script maker to make a .kim-like string on the fly and 
# stores the .kim file in the appropriate directory
#
# Please refer to standard.kim for OPENKIM related standards
#

def makekimscript(modelname,testname,atoms):
    """
    input: this file takes the modelname and testname and makes an appropriate test string for the model
    and also returns it as a string
    """
    kimstr = "TEST_NAME :=" +testname+"\n"

    # to get model units, inputs, outputs, options we call KIM_API_model_info
    status, km_pmdl = KIM_API_model_info(modelname) 
    # this needs to be a pointer    
 
    # BASE UNIT LINES
    unit_length = KIM_API_get_unit_length(km_pmdl)
    unit_energy = KIM_API_get_unit_energy(km_pmdl)    
    unit_charge = KIM_API_get_unit_charge(km_pmdl)
    unit_temperature = KIM_API_get_unit_temperature(km_pmdl)
    unit_time = KIM_API_get_unit_time(km_pmdl)    

    kimstr += "Unit_length :=" + unit_length +"\n"
    kimstr += "Unit_energy :=" + unit_energy +"\n"
    kimstr += "Unit_charge :=" + unit_charge +"\n"
    kimstr += "Unit_temperature :=" + unit_temperature +"\n"
    kimstr += "Unit_time :=" + unit_time +"\n"


    # SUPPORTED_ATOM/PARTICLE_TYPES
    
    kimstr += "SUPPORTED_AOM/PARTICLES_TYPES: \n"    

    # check ASE atoms class for which atoms it has
    acodes = set(atoms.get_atomic_numbers())
    asymbols = set(atoms.get_chemical_symbols())


    for code, symbol in zip(list(acodes),list(asymbols)):
        kimstr += symbol +" spec "+str(code) + "\n"   

    #print "put atoms and symbols in string"

    # CONVENTIONS
    kimstr += "CONVENTIONS:\n"
    
    # note: by default the convention for python is Zero-based lists
    kimstr += "ZeroBasedLists  flag\n"
    
    # Neighbor Access methods

    # note: KIM_API_get_index crashes python if the "name" doesn't exist
    # threfore we catch it with a try, except clause in checkIndex
    index1 = checkIndex(km_pmdl,"Neigh_IterAccess")
    index2 = checkIndex(km_pmdl,"Neigh_LocaAccess")
    index3 = checkIndex(km_pmdl, "Neigh_BothAccess")
    # index should be non-negative and which comes first would be preferable for the model 

    maxindex =  max([index1, index2, index3])
    #print maxindex

    if maxindex==index1 and index1>=0:
        kimstr += "Neigh_IterAccess  flag\n"
    elif maxindex==index2 and index2>=0:
        kimstr += "Neigh_LocaAccess  flag\n"
    elif maxindex==index3 and index3>=0:
        kimstr += "Neigh_BothAccess  flag\n"
    else:
        kimstr += "Neigh_LocaAccess  flag\n"
        print "WARNING: No Neighbor Access Method specified in Model"    

    #print "put neighbor access method in string"

    # Neighbor list and Boundary Condition (NBC) methods
    # here we can list all of the NBC methods because it will check against the model to find one that matches
  
    pbc = atoms.get_pbc()
    cell = atoms.get_cell()    
 
    cell_orthogonal = ((abs(numpy.dot(cell[0],cell[1])) + abs(numpy.dot(cell[0],cell[2])) + abs(numpy.dot(cell[1],cell[2])))<12**(-8))

    """
    index1 = checkIndex(km_pmdl,"NEIGH_PURE_H")
    index2 = checkIndex(km_pmdl, "NEIGH_PURE_F")
    index3 = checkIndex(km_pmdl,"NEIGH_RVEC_F")
    index4 = checkIndex(km_pmdl, "MI_OPBC_H")
    index5 = checkIndex(km_pmdl,"MI_OPBC_F")
    index6 = checkIndex(km_pmdl, "CLUSTER")

    maxindex = max([index1,index2,index3,index4,index5,index6])    
    """    

    if pbc.any():
        kimstr += "NEIGH_RVEC_F flag \n"
        if cell_orthogonal:
            kimstr += "MI_OPBC_H flag \n"
            kimstr += "MI_OPBC_F flag \n"
    else:
        kimstr += "NEIGH_PURE_H flag \n"
        kimstr += "NEIGH_PURE_F flag \n"
        kimstr += "CLUSTER flag \n"        
    # to check if something is there we should use KIM_API_get_index(km_pmdl, "string")
    
    # MODEL_INPUT section
    kimstr += "MODEL_INPUT:\n"
    if checkIndex(km_pmdl,"numberOfParticles")>=0:
        kimstr +="numberOfParticles  integer  none  []\n"
    if checkIndex(km_pmdl, "numberContributingParticles")>=0:
        kimstr +="numberContributingParticles  integer  none  []\n"
    if checkIndex(km_pmdl, "numberParticleTypes")>=0:
        kimstr +="numberParticleTypes  integer  none  []\n"
    if checkIndex(km_pmdl, "particleTypes")>=0:
        kimstr +="particleTypes  integer  none  [numberOfParticles]\n"
    if checkIndex(km_pmdl, "coordinates")>=0:
       kimstr +="coordinates  real*8  length  [numberOfParticles,3]\n" 
    if checkIndex(km_pmdl, "particleCharge")>=0:
        kimstr +="particleCharge  real*8  charge  [numberOfParticles]\n"
    if checkIndex(km_pmdl, "particleSize")>=0:
        kimstr +="particleSize  real*8  length  [numberOfParticles]\n"
    # confused about the get_neigh and neighObject pointer.  this is test dependent.  do we include this?   
    # how do we decide whether or not to include this? Include for now. decide later
    #################################################  
    if checkIndex(km_pmdl, "get_neigh")>=0:
        kimstr += "get_neigh  method  none []\n"
    if checkIndex(km_pmdl, "neighObject")>=0:
        kimstr += "neighObject  pointer  none  []\n"
    #############################################
    if checkIndex(km_pmdl, "process_dEdr")>=0:
        kimstr +="process_dEdr  method  none  []\n"
    if checkIndex(km_pmdl, "process_d2Edr2")>=0:
        kimstr +="process_d2Edr2  method  none  []\n"
    if checkIndex(km_pmdl, "boxSideLengths")>=0:
        kimstr +="boxSideLengths  real*8  length  [3]\n"
      
    # MODEL_OUTPUT section
    kimstr += "MODEL_OUTPUT: \n"  
    if checkIndex(km_pmdl, "compute")>=0:
        kimstr += "compute  method  none  []\n"
    if checkIndex(km_pmdl, "reinit") >= 0:
        kimstr += "reinit  method  none  []\n"
    if checkIndex(km_pmdl, "destroy") >= 0:
        kimstr += "destroy  method  none  []\n"
    if checkIndex(km_pmdl, "cutoff") >= 0:
        kimstr += "cutoff  real*8  length  []\n"
    if checkIndex(km_pmdl, "energy") >= 0:
        kimstr += "energy  real*8  energy  []\n"
    if checkIndex(km_pmdl, "forces") >= 0:
        kimstr += "forces  real*8  force  [numberOfParticles,3]\n"
    if checkIndex(km_pmdl, "particleEnergy") >=0 :
        kimstr += "particleEnergy  real*8  energy  [numberOfParticles]\n"
    if checkIndex(km_pmdl, "virial") >= 0:
        kimstr += "virial  real*8  energy  [6]\n"
    if checkIndex(km_pmdl, "particleVirial") >=0 :
        kimstr += "particleVirial  real*8  energy  [numberOfParticles,6]\n"
    if checkIndex(km_pmdl, "hessian") >= 0:
        kimstr += "hessian  real*8  pressure  [numberOfParticles,numberOfParticles,3,3]\n"
    
    return kimstr


def checkIndex(pkim,variablename):
    try:
        index = KIM_API_get_index(pkim,variablename)
    except:
        index = -1 
    return index

