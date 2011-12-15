import numpy
import numpy.linalg as la
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

def set_NeighborList(pkim, coords, numberOfAtoms, cutoff):
    cutoff2 = cutoff*2
    neighbors = []
    fneighbors = []
    for i in range(numberOfAtoms):
        neighbors.append([])
        fneighbors.append([])
    for i in range(numberOfAtoms):
        for j in range(i+1, numberOfAtoms):
            dx = coords[i*DIM:i*DIM+DIM] - coords[j*DIM:j*DIM+DIM]
            dist = (dx**2).sum()
            if (dist < cutoff2):
                # zero-based
                fneighbors[j].append((i, dx))
                # full list
                neighbors[i].append((j,dx))

    # merge list so that for each i, [neighbors(j>i), neighbors(j<i)]
    for i in range(numberOfAtoms):
        fneighbors[i] = neighbors[i] + fneighbors[i]

    #drop one for "zero"
    HalfNNeighbors = numpy.array([len(l) for l in neighbors], dtype='int32')
    NNeighbors = numpy.array([len(l) for l in fneighbors], dtype='int32')
    neighborList = []
    RijList = []
    for l in fneighbors:
        neighs = [n for n,r in l]
        Rijs = [r for n,r in l]
        neighborList += neighs
        RijList += Rijs
    neighborList = numpy.array(neighborList, dtype='int32')
    RijList = numpy.array(RijList, dtype='double').flatten()
    return NNeighbors, HalfNNeighbors, neighborList, RijList

testname = "test_python_openkim_neigh"
testf = open("test_python_openkim_neigh.kim", "rt")
teststring = "".join(testf.readlines())
modelname = raw_input("Please enter a valid KIM model name:")

#status, pkim = KIM_API_init(testname, modelname)
status, pkim = KIM_API_init_str(teststring, modelname)
if KIM_STATUS_OK > status:
    KIM_API_report_error('KIM_API_init',status)

cnt = 1
try:
    KIM_API_allocate(pkim, NCLUSTERATOMS, ATYPES)

    numberOfAtoms = KIM_API_get_data_ulonglong(pkim, "numberOfAtoms")
    numberAtomTypes = KIM_API_get_data_int(pkim, "numberAtomTypes")
    numberContributingAtoms=KIM_API_get_data_int(pkim,"numberContributingAtoms")
    atomTypes = KIM_API_get_data_int(pkim, "atomTypes")
    coordinates = KIM_API_get_data_double(pkim, "coordinates")
    cutoff = KIM_API_get_data_double(pkim, "cutoff")
    energy = KIM_API_get_data_double(pkim, "energy")
    forces = KIM_API_get_data_double(pkim, "forces")

    # Set values 
    numberOfAtoms[0] = NCLUSTERATOMS
    numberContributingAtoms[0]=NCLUSTERATOMS
    numberAtomTypes[0] = ATYPES

    status = KIM_API_model_init(pkim)
    if KIM_STATUS_OK > status:
        raise kimservice.error("KIM_API_model_init")
    import neighborlist
    status = neighborlist.set_kim_periodic_full_neigh(pkim)
    status = neighborlist.set_kim_periodic_half_neigh(pkim)

    KIM_API_print(pkim)

    atypecode = KIM_API_get_aTypeCode(pkim, "Ar")
    for i in range(numberOfAtoms[0]):
        atomTypes[i] = atypecode
    
    MiddleAtomId = create_FCC_configuration(FCCSPACING, NCELLSPERSIDE, 0, coordinates)
   
    NNeighbors, HalfNNeighbors, neighborList, RijList = set_NeighborList(pkim, coordinates, numberOfAtoms[0], cutoff[0]*2)
    neighborlist.set_neigh_object(pkim, NNeighbors, HalfNNeighbors, neighborList, RijList)

    KIM_API_model_compute(pkim)

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
    neighborlist.free_neigh_object(pkim)
    KIM_API_model_destroy(pkim)
    KIM_API_free(pkim)
except error:
    KIM_API_report_error(error.message,errno)


