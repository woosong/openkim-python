%module kimneighborlist
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}

%{
#include "neighborlist.h"
%}
%apply (int DIM1, int *IN_ARRAY1)    {(int sz1, int* NNeighbors)}
%apply (int DIM1, int *IN_ARRAY1)    {(int sz2, int* HalfNNeighbors)}
%apply (int DIM1, int *IN_ARRAY1)    {(int sz3, int* neighborList)}
%apply (int DIM1, double *IN_ARRAY1) {(int sz4, double* RijList)}
%apply (int DIM1, double *IN_ARRAY1) {(int S1, double* Cell)}
%apply (int DIM1, int *IN_ARRAY1)    {(int S2, int* PBC)}
%apply (int DIM1, int *IN_ARRAY1)    {(int N1, int* Ghosts)}
%include "neighborlist.h"

