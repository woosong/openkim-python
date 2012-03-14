import numpy
import pylab
from kimservice import *

FCCSPACING = 5.260
NCELLSPERSIDE = 2
DIM = 3
ATYPES = 1
NCLUSTERATOMS = (4*(NCELLSPERSIDE*NCELLSPERSIDE*NCELLSPERSIDE) + 6*(NCELLSPERSIDE*NCELLSPERSIDE) + 3*(NCELLSPERSIDE) + 1)

#
# create_FCC_configuration function
#
#  creates a cubic configuration of FCC atoms with lattice spacing
#  `FCCspacing' and `nCellsPerSide' cells along each direction.
#
#  With periodic==0. this will result in a total number of atoms equal to
#  4*(nCellsPerSide)**3 + 6*(nCellsPerSide)**2 + 3*(nCellsPerSide) + 1
#
#  With periodic==1 this will result in a total number of atoms equal to
#  4*(nCellsPerSide)**3
#
#  Returns the Id of the atom situated in the middle of the configuration
#  (this atom will have the most neighbors.)
#
def create_FCC_configuration(FCCspacing, nCellsPerSide, periodic, coords):
    # local variables 
    FCCshifts = numpy.zeros((4,3), dtype=numpy.float64)
    latVec = numpy.zeros(3, dtype=numpy.float64)
    MiddleAtomId = 0

    # create a cubic FCC cluster 
    FCCshifts[0,0] = 0.0
    FCCshifts[0,1] = 0.0
    FCCshifts[0,2] = 0.0
    
    FCCshifts[1,0] = 0.5*FCCspacing
    FCCshifts[1,1] = 0.5*FCCspacing
    FCCshifts[1,2] = 0.0
    
    FCCshifts[2,0] = 0.5*FCCspacing
    FCCshifts[2,1] = 0.0
    FCCshifts[2,2] = 0.5*FCCspacing
    
    FCCshifts[3,0] = 0.0
    FCCshifts[3,1] = 0.5*FCCspacing
    FCCshifts[3,2] = 0.5*FCCspacing
   
    a = 0
    for i in range(nCellsPerSide):
        latVec[0] = i*FCCspacing
        for j in range(nCellsPerSide):
            latVec[1] = j*FCCspacing
            for k in range(nCellsPerSide):
                latVec[2] = k*FCCspacing
                for m in range(4):
                    coords[a*DIM:(a+1)*DIM] = latVec + FCCshifts[m]
                    if i == nCellsPerSide/2 and j == nCellsPerSide/2 and k == nCellsPerSide/2 and m == 1:
                        MiddleAtomId = a
                    a += 1
            if not periodic:
                # add in the remaining three faces
                # pos-x face
                latVec[0] = NCELLSPERSIDE*FCCspacing
                latVec[1] = i*FCCspacing
                latVec[2] = j*FCCspacing
                coords[a*DIM:(a+1)*DIM] = latVec
                a += 1
                coords[a*DIM:(a+1)*DIM] = latVec + FCCshifts[3]
                a += 1
                # pos-y face
                latVec[0] = i*FCCspacing
                latVec[1] = NCELLSPERSIDE*FCCspacing
                latVec[2] = j*FCCspacing
                coords[a*DIM:(a+1)*DIM] = latVec
                a += 1
                coords[a*DIM:(a+1)*DIM] = latVec + FCCshifts[2]
                a += 1
                # pos-z face
                latVec[0] = i*FCCspacing
                latVec[1] = j*FCCspacing
                latVec[2] = NCELLSPERSIDE*FCCspacing
                coords[a*DIM:(a+1)*DIM] = latVec
                a += 1
                coords[a*DIM:(a+1)*DIM] = latVec + FCCshifts[1]
                a += 1
        if not periodic:
            # add in the remaining three faces
            # pos-x face
            latVec[0] = i*FCCspacing
            latVec[1] = NCELLSPERSIDE*FCCspacing
            latVec[2] = NCELLSPERSIDE*FCCspacing
            coords[a*DIM:(a+1)*DIM] = latVec
            a += 1
            # pos-y face
            latVec[0] = NCELLSPERSIDE*FCCspacing
            latVec[1] = i*FCCspacing
            latVec[2] = NCELLSPERSIDE*FCCspacing
            coords[a*DIM:(a+1)*DIM] = latVec
            a += 1
            # pos-z face
            latVec[0] = NCELLSPERSIDE*FCCspacing
            latVec[1] = NCELLSPERSIDE*FCCspacing
            latVec[2] = i*FCCspacing
            coords[a*DIM:(a+1)*DIM] = latVec
            a += 1
    if not periodic:
        for n in range(DIM):
            coords[a*DIM + n] = NCELLSPERSIDE*FCCspacing
        a += 1

    return MiddleAtomId

testname = "test_python_openkim"
testf = open("test_python_openkim.kim", "rt")
teststring = "".join(testf.readlines())
modelname = raw_input("Please enter a valid KIM model name:")

#status, pkim = KIM_API_init(testname, modelname)
status, pkim = KIM_API_init_str(teststring, modelname)
if KIM_STATUS_OK > status:
    KIM_API_report_error('KIM_API_init',status)
try:
    KIM_API_allocate(pkim, NCLUSTERATOMS, ATYPES)
    status = KIM_API_model_init(pkim)
    if KIM_STATUS_OK > status:
        raise kimservice.error("KIM_API_model_init")
    numberOfAtoms = KIM_API_get_data_ulonglong(pkim, "numberOfParticles")
    numberAtomTypes = KIM_API_get_data_int(pkim, "numberParticleTypes")
    atomTypes = KIM_API_get_data_int(pkim, "particleTypes")
    coordinates = KIM_API_get_data_double(pkim, "coordinates")
    cutoff = KIM_API_get_data_double(pkim, "cutoff")
    energy = KIM_API_get_data_double(pkim, "energy")
    forces = KIM_API_get_data_double(pkim, "forces")

    # Set values
    numberOfAtoms[0] = NCLUSTERATOMS
    numberAtomTypes[0] = ATYPES
   
    atypecode = KIM_API_get_partcl_type_code(pkim, "Ar")
    
    for i in range(numberOfAtoms[0]):
        atomTypes[i] = atypecode
    
    MiddleAtomId = create_FCC_configuration(FCCSPACING, NCELLSPERSIDE, 0, coordinates)
    KIM_API_model_compute(pkim)
     
    print "computed"
except error:
    KIM_API_report_error(error.message,errno)

# print results to screen
print "--------------------------------------------------------------------------------"
print "This is Test          : %s" % testname
print "Results for KIM Model : %s" % modelname
print "Forces:"
print "Atom     X                        Y                        Z"
for i in range(numberOfAtoms[0]):
    print "%2i   %25.15e%25.15e%25.15e" % (i, forces[i*DIM + 0], forces[i*DIM + 1], forces[i*DIM + 2])
print ""
print "Energy = %20.15e" %  energy[0]

try:
    KIM_API_model_destroy(pkim)
    KIM_API_free(pkim)
except error:
    KIM_API_report_error(error.message,errno)

