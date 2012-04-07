import numpy
import numpy.linalg as la
from kimservice import *
import kimneighborlist
import virial

# Copper
FCCSPACING = 3.715
NCELLSPERSIDE = 2
DIM = 3
ATYPES = 1

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

testname = "test_python_openkim_virial"
testf = open("test_python_openkim_virial.kim", "rt")
teststring = "".join(testf.readlines())
modelname = raw_input("Please enter a valid KIM model name:")

coords_dump = numpy.fromfile('cu_coordinates.dump', sep=' ')
coords_dump = coords_dump.reshape((len(coords_dump)/3, 3))
NCLUSTERATOMS = len(coords_dump)

#status, pkim = KIM_API_init(testname, modelname)
status, pkim = KIM_API_init_str(teststring, modelname)
if KIM_STATUS_OK > status:
    KIM_API_report_error('KIM_API_init',status)

cnt = 1
try:
    KIM_API_allocate(pkim, NCLUSTERATOMS, ATYPES)
    #get everything pointers KIM
    numberOfAtoms = KIM_API_get_data_ulonglong(pkim, "numberOfParticles")
    numberAtomTypes = KIM_API_get_data_int(pkim, "numberParticleTypes")
    numberContributingAtoms=KIM_API_get_data_int(pkim,"numberContributingParticles")
    atomTypes = KIM_API_get_data_int(pkim, "particleTypes")
    coordinates = KIM_API_get_data_double(pkim, "coordinates")
    cutoff = KIM_API_get_data_double(pkim, "cutoff")
    energy = KIM_API_get_data_double(pkim, "energy")
    forces = KIM_API_get_data_double(pkim, "forces")
    virialGlobal = KIM_API_get_data_double(pkim, "virial")
    #virialPerAtom = KIM_API_get_data_double(pkim, "particleVirial")

    # Set values
    numberOfAtoms[0] = NCLUSTERATOMS
    numberContributingAtoms[0]=NCLUSTERATOMS
    numberAtomTypes[0] = ATYPES

    if KIM_STATUS_OK > status:
        raise kimservice.error("KIM_API_model_init")
#    status = kimneighborlist.set_kim_periodic_full_neigh(pkim)
#    status = kimneighborlist.set_kim_periodic_half_neigh(pkim)
    kimneighborlist.nbl_initialize(pkim)

    status = virial.virial_init(pkim)
    status = virial.set_virial(pkim)

    status = KIM_API_model_init(pkim)

    atypecode = KIM_API_get_partcl_type_code(pkim, "Ar")
    for i in range(numberOfAtoms[0]):
        atomTypes[i] = atypecode
    
    coordinates[:] = coords_dump.flatten()[:]

    KIM_API_print(pkim)

#    NNeighbors, HalfNNeighbors, neighborList, RijList = set_NeighborList(pkim, coordinates, numberOfAtoms[0], cutoff[0]*4)
#    kimneighborlist.set_neigh_object(pkim, NNeighbors, HalfNNeighbors, neighborList, RijList)
    kimneighborlist.nbl_build_neighborlist(pkim)

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

print virialGlobal#, virialPerAtom

try:
    kimneighborlist.nbl_cleanup(pkim)
    KIM_API_model_destroy(pkim)
    KIM_API_free(pkim)
except error:
    KIM_API_report_error(error.message,errno)


