%module kimservice
%{
#define SWIG_FILE_WITH_INIT
%}

%{
static PyObject *kimerror;
static PyIntObject *kimerrno;
%}
%pythoncode %{
error=_kimservice.error
errno=_kimservice.errno
%}

%typemap(argout) int *error {
    kimerrno->ob_ival = *$1;
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

%{
#include "KIM_API_C.h"
#include "KIM_API_status.h"
typedef void* voidp;
int KIM_API_init_python(void **pkimmdl, char *testname, char *modelname) {
    return KIM_API_init((void*)pkimmdl, testname, modelname);
};
int KIM_API_init_python_str(void **pkimmdl, char *teststring, char *modelname) {
    return KIM_API_string_init((void*)pkimmdl, teststring, modelname);
};
void KIM_API_free_python(void **pkimmdl, int *error) {
    KIM_API_free((void*)pkimmdl, error);
};

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


%include "cpointer.i"
%typedef void* voidp;
%pointer_functions(voidp, voidpp);
%include "KIM_API_C.h"
%include "KIM_API_status.h"

int KIM_API_init_python(void **pkimmdl, char *testname, char *modelname);
%pythoncode %{
def KIM_API_init(testname, modelname):
    ppkim = new_voidpp()
    status = KIM_API_init_python(ppkim, testname, modelname)
    pkim = voidpp_value(ppkim)
    return (status, pkim) 
%}
int KIM_API_init_python_str(void **pkimmdl, char *teststring, char *modelname);
%pythoncode %{
def KIM_API_init_str(teststring, modelname):
    ppkim = new_voidpp()
    status = KIM_API_init_python_str(ppkim, teststring, modelname)
    pkim = voidpp_value(ppkim)
    return (status, pkim) 
%}
void KIM_API_free_python(void **pkimmdl, int *error);
%pythoncode %{
def KIM_API_free(pkim):
    ppkim = new_voidpp()
    voidpp_assign(ppkim, pkim)
    return KIM_API_free_python(ppkim)
%}

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
