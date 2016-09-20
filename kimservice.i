/* ==============================================================
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common Development
// and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
//
// CDDL HEADER END
//
// ==============================================================
//
// KIM API :: python interface
//  
// Provides a swig-ed python interface to the Knowledge-base of
// Interatomic Models (KIM) API 
//  
// Authors:
//     Woosong Choi    
//     Matt Bierbaum   
//     YJ Chen
//
// ==================================================================*/


/* 
    Defines the modules name: kimservice
    next lines allow for numpy to work properly through
    double *data types.
*/
%module kimservice
%{
#define SWIG_FILE_WITH_INIT
%}

%{
static PyObject *kimerror = NULL;
static PyIntObject *kimerrno = NULL;
%}
%pythoncode %{
error=_kimservice.error
errno=_kimservice.errno
%}

/* all functions that return an error are treated as exceptions */
%typemap(argout) int *error {
    // TODO find out why this breaks
    //kimerrno->ob_ival = *$1;
    if (KIM_STATUS_OK > *$1) {
        PyErr_SetString(kimerror, "$symname");
        return NULL;
    }
} 
%typemap(in,numinputs=0) int *error (int temp) {
    $1 = &temp;
}

%include "numpy.i"
%init %{
import_array();

    // also init the error object
    kimerror = PyErr_NewException("kimservice.error", NULL, NULL);
    Py_INCREF(kimerror);
    PyModule_AddObject(m, "error", kimerror);

    kimerrno = (PyIntObject *)PyInt_FromLong(KIM_STATUS_OK);
    Py_INCREF(kimerrno);
    PyModule_AddObject(m, "errno", (PyObject *)kimerrno);
%}
%apply (double* IN_ARRAY1, int DIM1) {(double* seq, int n)};

/* 
    Individually wrap the initialization functions (and others that handle
    modification of pointers) so that they work well
*/
%{
#include "KIM_API_C.h"
#include "KIM_API_status.h"
typedef void* voidp;
int KIM_API_file_init_python(void **pkimmdl, char *testname, char *modelname) {
    return KIM_API_file_init((void*)pkimmdl, testname, modelname);
};
int KIM_API_init_str_python(void **pkimmdl, char *teststring, char *modelname) {
    return KIM_API_string_init((void*)pkimmdl, teststring, modelname);
};
int KIM_API_model_info_python(void **pmdl, char *modelname){
    return KIM_API_model_info((void*)pmdl, modelname);
};
char *const KIM_API_get_model_kim_str_python(const char *modelname){
    char *kimString;
    KIM_API_get_model_kim_str(modelname, &kimString);
    return kimString;
};
void KIM_API_free_python(void **pkimmdl, int *error) {
    KIM_API_free((void*)pkimmdl, error);
};
const char *KIM_API_get_NBC_method_python(void *pkimmdl){
    const char *method;
    KIM_API_get_NBC_method(pkimmdl, &method);
    return method;
}
int KIM_API_get_num_model_species_python(void *pkimmdl){
    int num, len;
    KIM_API_get_num_model_species(pkimmdl, &num, &len);
    return num;
}
const char *KIM_API_get_model_species_python(void *pkimmdl, int ind){
    const char *species;
    KIM_API_get_model_species(pkimmdl, ind, &species);
    return species;
}
int KIM_API_get_num_free_params_python(void *pkimmdl){
    int nfp, msl;
    KIM_API_get_num_free_params(pkimmdl, &nfp, &msl);
    return nfp;
}
int KIM_API_get_rank_python(void *pkimmdl, const char *nm, int *error){
    return (int)KIM_API_get_rank(pkimmdl, nm, error);
}
int KIM_API_get_shape_python(void *pkimmdl, const char *nm, int *error){
    int shape;
    KIM_API_get_shape(pkimmdl, nm, &shape, error);
    return shape;
}
const char *KIM_API_get_free_parameter_python(void *pkimmdl, const int index){
    const char *param;
    KIM_API_get_free_parameter(pkimmdl, index, &param);
    return param;
}
/* 
    templated functions to set the necessary data pointers, use these at
    the bottom of this file to define links for all data types
*/
template<typename T>
void KIM_API_get_data_dtype(T **data, int *n, void *kimmdl,char *nm,int *error) {
    *data = (T *)KIM_API_get_data(kimmdl, nm, error);
    *n = (int)KIM_API_get_size(kimmdl,nm,error);
};

template<typename T>
int KIM_API_set_data_dtype(void *kimmdl, char *nm,int size, T *dt) {
    return KIM_API_set_data(kimmdl, nm, (intptr_t)size, (void *)dt);
};

int KIM_API_set_data_single_intptr_tp(void *kimmdl,char *nm, intptr_t *p) {
    return KIM_API_set_data(kimmdl, nm, (intptr_t)1, (void *)p);
};
%}

template<typename T> 
void KIM_API_get_data_dtype(T **data, int *n, void *kimmdl,char *nm,int *error); 
template<typename T> 
int KIM_API_set_data_dtype(void *kimmdl, char *nm,int size, T *dt); 


/* 
    SWIG the entirety of the KIM API, leaving the last few functions
    that require special attention (numpy arrays, void **) to be 
    address in the subsequent lines
*/
%include "cpointer.i"
%typedef void* voidp;
%pointer_functions(voidp, voidpp);
%include "KIM_API_C.h"
%include "KIM_API_status.h"

/*  These functions until EOF require additional attention to work */
int KIM_API_file_init_python(void **pkimmdl, char *testkimfile, char *modelname);
%pythoncode %{
def KIM_API_file_init(testkimfile, modelname):
    ppkim = new_voidpp()
    status = KIM_API_file_init_python(ppkim, testkimfile, modelname)
    pkim = voidpp_value(ppkim)
    return (status, pkim) 
%}
int KIM_API_init_str_python(void **pkimmdl, char *teststring, char *modelname);
%pythoncode %{
def KIM_API_init_str(teststring, modelname):
    ppkim = new_voidpp()
    status = KIM_API_init_str_python(ppkim, teststring, modelname)
    pkim = voidpp_value(ppkim)
    return (status, pkim) 
%}
char *const KIM_API_get_model_kim_str_python(const char *modelname);
%pythoncode %{
def KIM_API_get_model_kim_str(modelname):
    return KIM_API_get_model_kim_str_python(modelname)
%}
int KIM_API_model_info_python(void **pmdl, char *modelname);
%pythoncode %{
def KIM_API_model_info(modelname):
    ppmdl = new_voidpp()
    status = KIM_API_model_info_python(ppmdl, modelname)
    pmdl = voidpp_value(ppmdl)
    return (status, pmdl)
%}
void KIM_API_free_python(void **pkimmdl, int *error);
%pythoncode %{
def KIM_API_free(pkim):
    ppkim = new_voidpp()
    voidpp_assign(ppkim, pkim)
    return KIM_API_free_python(ppkim)
%}
const char *KIM_API_get_NBC_method_python(void *pkimmdl);
%pythoncode %{
def KIM_API_get_NBC_method(pkim):
    return KIM_API_get_NBC_method_python(pkim)
%}
int KIM_API_get_num_model_species_python(void *pkimmdl);
%pythoncode %{
def KIM_API_get_num_model_species(pkim):
    return KIM_API_get_num_model_species_python(pkim)
%}
const char *KIM_API_get_model_species_python(void *pkimmdl, int ind);
%pythoncode %{
def KIM_API_get_model_species(pkim, ind):
    return KIM_API_get_model_species_python(pkim, ind)
%}
int KIM_API_get_num_free_params_python(void *pkimmdl);
%pythoncode %{
def KIM_API_get_num_free_params(pkim):
    return KIM_API_get_num_free_params_python(pkim)
%}
int KIM_API_get_rank_python(void *pkimmdl, const char *nm, int *error);
%pythoncode %{
def KIM_API_get_rank(pkim, nm):
    return int(KIM_API_get_rank_python(pkim, nm))
%}
int KIM_API_get_shape_python(void *pkimmdl, const char *nm, int *error);
%pythoncode %{
def KIM_API_get_shape(pkim, nm):
    return KIM_API_get_shape_python(pkim, nm)
%}
const char *KIM_API_get_free_parameter_python(void *pkimmdl, const int index);
%pythoncode %{
def KIM_API_get_free_parameter(pkim, index):
    return KIM_API_get_free_parameter_python(pkim, index)
%}

/*
    the library of templated function calls to access and modify
    data that lives in the KIM API
*/
%apply (double** ARGOUTVIEW_ARRAY1, int* DIM1) {(double **data, int *n)}
%template(KIM_API_get_data_double) KIM_API_get_data_dtype<double>;
%apply (int** ARGOUTVIEW_ARRAY1, int* DIM1) {(int **data, int *n)}
%template(KIM_API_get_data_int) KIM_API_get_data_dtype<int>;
%apply (unsigned int** ARGOUTVIEW_ARRAY1, int* DIM1) {(unsigned int **data, int *n)}
%template(KIM_API_get_data_uint) KIM_API_get_data_dtype<unsigned int>;
%apply (long ** ARGOUTVIEW_ARRAY1, int* DIM1) {(long **data, int *n)}
%template(KIM_API_get_data_long) KIM_API_get_data_dtype<long>;
%apply (unsigned long ** ARGOUTVIEW_ARRAY1, int* DIM1) {(unsigned long **data, int *n)}
%template(KIM_API_get_data_ulong) KIM_API_get_data_dtype<unsigned long>;
%apply (unsigned long long ** ARGOUTVIEW_ARRAY1, int* DIM1) {(unsigned long long **data, int *n)}
%template(KIM_API_get_data_ulonglong) KIM_API_get_data_dtype<unsigned long long>;

%apply (int DIM1, int *IN_ARRAY1) {(int size, int *dt)}
%template(KIM_API_set_data_int) KIM_API_set_data_dtype<int>;
%apply (int DIM1, unsigned int *IN_ARRAY1) {(int size, unsigned int *dt)}
%template(KIM_API_set_data_uint) KIM_API_set_data_dtype<unsigned int>;
%apply (int DIM1, double *IN_ARRAY1) {(int size, double *dt)}
%template(KIM_API_set_data_double) KIM_API_set_data_dtype<double>;
%apply (int DIM1, long *IN_ARRAY1) {(int size, long *dt)}
%template(KIM_API_set_data_long) KIM_API_set_data_dtype<long>;
%apply (int DIM1, unsigned long *IN_ARRAY1) {(int size, unsigned long *dt)}
%template(KIM_API_set_data_ulong) KIM_API_set_data_dtype<unsigned long>;
int KIM_API_set_data_single_intptr_tp(void *kimmdl,char *nm,intptr_t *p); 

%pythoncode %{
import inspect

def KIM_API_report_error(usermsg, ier):
    istack = inspect.stack()[1]
    errline = istack[2]
    errfile = istack[1]
    return _kimservice.KIM_API_report_error(errline, errfile, usermsg, ier)
%}
