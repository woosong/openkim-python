%module neighborlist
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
%apply (int DIM1, int *IN_ARRAY1) {(int sz1, int* NNeighbors)}
%apply (int DIM1, int *IN_ARRAY1) {(int sz2, int* HalfNNeighbors)}
%apply (int DIM1, int *IN_ARRAY1) {(int sz3, int* neighborList)}
%apply (int DIM1, double *IN_ARRAY1) {(int sz4, double* RijList)}
%include "neighborlist.h"

