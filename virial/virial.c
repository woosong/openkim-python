#include "KIM_API_C.h"
#include "KIM_API_status.h"
#include "virial.h"
#include <stdlib.h>
#include <string.h>

void process_dEdr(void** ppkim,double *de,double *r,double ** pdx,int *i,int *j,int *ier);
void process_d2Edr2(void** ppkim,double *d2e,double **r,double ** pdx,int **i,int **j,int *ier);

static double *virialGlobal = NULL;
static double *virialPerAtom = NULL;
static double *hessian = NULL;
static int *numberOfAtoms = NULL;
static int virialGlobal_flag, virialPerAtom_flag, hessian_flag;
static int virialGlobal_need2add = 1;
static int virialPerAtom_need2add = 1;
static int hessian_need2add = 1;

int virial_init(void *pkim){
    virialGlobal_flag = 0;
    virialPerAtom_flag =0;
    hessian_flag = 0;

    int ierGlobal,ierPerAtom,ierHessian,kimerr;
    int i;
    virialGlobal  = (double *) KIM_API_get_data(pkim, "virial",&ierGlobal);
    virialPerAtom = (double *) KIM_API_get_data(pkim, "particleVirial",&ierPerAtom);
    hessian       = (double *) KIM_API_get_data(pkim, "hessian", &ierHessian);
    numberOfAtoms = (int *) KIM_API_get_data(pkim, "numberOfParticles",&kimerr);
    //halfNeighbors = !pkim->requiresFullNeighbors();

    if (kimerr !=KIM_STATUS_OK) return kimerr;
    if (ierGlobal == KIM_STATUS_OK && virialGlobal != NULL) {
        virialGlobal_flag = KIM_API_get_compute(pkim, "virial", &kimerr);
        if (virialGlobal_flag==1 && virialGlobal_need2add) *virialGlobal =0.0;
    }

    if (ierPerAtom == KIM_STATUS_OK && virialPerAtom != NULL) {
        virialPerAtom_flag = KIM_API_get_compute(pkim, "particleVirial", &kimerr);
        if (virialPerAtom_flag==1 && virialPerAtom_need2add) {
            for (i =0;i<(*numberOfAtoms)*6 ;i++) virialPerAtom[i]=0.0;
        }
    }

    if (ierHessian == KIM_STATUS_OK && hessian != NULL) {
        hessian_flag = KIM_API_get_compute(pkim, "hessian", &kimerr);
        if (hessian_flag==1 && hessian_need2add) {
            for (i=0; i<(*numberOfAtoms)*(*numberOfAtoms)*9; i++) hessian[i] = 0.0;
        }
    }
    return kimerr;
}

int set_virial(void* pkim) {
    return KIM_API_set_data(pkim, "process_dEdr", 1, (void*) &process_dEdr);
}

int set_hessian(void* pkim) {
    return KIM_API_set_data(pkim, "process_d2Edr2", 1, (void*) &process_d2Edr2);
}

void process_dEdr(void **ppkim,double *de,double *r,double ** pdx,int *i,int *j,int *ier){
   *ier=KIM_STATUS_FAIL;
   double vir[6],v;
   double *dx = *pdx;
   //v=(*de)/((*r)*(*r))/3.0;
   //v=(*de)/((*r)*(*r));

   v=-(*de)/(*r);
 

   if (virialGlobal_flag ==1 && virialGlobal_need2add) {
       static int numruns1 =0;
       numruns1++;
       vir[0] = v * dx[0] * dx[0];
       vir[1] = v * dx[1] * dx[1];
       vir[2] = v * dx[2] * dx[2];
       vir[3] = v * dx[1] * dx[2];
       vir[4] = v * dx[0] * dx[2];
       vir[5] = v * dx[0] * dx[1];
       virialGlobal[0] += vir[0];
       virialGlobal[1] += vir[1];
       virialGlobal[2] += vir[2];
       virialGlobal[3] += vir[3];
       virialGlobal[4] += vir[4];
       virialGlobal[5] += vir[5];
    }
    if (virialPerAtom_flag==1 && virialPerAtom_need2add ){
       static int numruns =0;
       numruns++;
       vir[0] =0.5 * v * dx[0] * dx[0];
       vir[1] =0.5 * v * dx[1] * dx[1];
       vir[2] =0.5 * v * dx[2] * dx[2];
       vir[3] =0.5 * v * dx[1] * dx[2];
       vir[4] =0.5 * v * dx[0] * dx[2];
       vir[5] =0.5 * v * dx[0] * dx[1];
       virialPerAtom[(*i)*6 + 0] += vir[0];
       virialPerAtom[(*i)*6 + 1] += vir[1];
       virialPerAtom[(*i)*6 + 2] += vir[2];
       virialPerAtom[(*i)*6 + 3] += vir[3];
       virialPerAtom[(*i)*6 + 4] += vir[4];
       virialPerAtom[(*i)*6 + 5] += vir[5];

       virialPerAtom[(*j)*6 + 0] += vir[0];
       virialPerAtom[(*j)*6 + 1] += vir[1];
       virialPerAtom[(*j)*6 + 2] += vir[2];
       virialPerAtom[(*j)*6 + 3] += vir[3];
       virialPerAtom[(*j)*6 + 4] += vir[4];
       virialPerAtom[(*j)*6 + 5] += vir[5];
     }
   
  
    *ier = KIM_STATUS_OK;
}


void process_d2Edr2(void **ppkim,double *de,double **r,double ** pdx,int **i,int **j,int *ier){
   *ier=KIM_STATUS_FAIL;
   double vir[6],v;
   double *dx = *pdx;
   //v=(*de)/((*r)*(*r))/3.0;
   //v=(*de)/((*r)*(*r));

   v=-(*de)/(**r);
 
/*
   if (hessian_flag ==1 && hessian_need2add) {
 static int numruns1 =0;
 numruns1++;
       vir[0] = v * dx[0] * dx[0];
       vir[1] = v * dx[1] * dx[1];
       vir[2] = v * dx[2] * dx[2];
       vir[3] = v * dx[1] * dx[2];
       vir[4] = v * dx[0] * dx[2];
       vir[5] = v * dx[0] * dx[1];
       virialGlobal[0] += vir[0];
       virialGlobal[1] += vir[1];
       virialGlobal[2] += vir[2];
       virialGlobal[3] += vir[3];
       virialGlobal[4] += vir[4];
       virialGlobal[5] += vir[5];
    }
    if (virialPerAtom_flag==1 && virialPerAtom_need2add ){
static int numruns =0;
numruns++;
       vir[0] =0.5 * v * dx[0] * dx[0];
       vir[1] =0.5 * v * dx[1] * dx[1];
       vir[2] =0.5 * v * dx[2] * dx[2];
       vir[3] =0.5 * v * dx[1] * dx[2];
       vir[4] =0.5 * v * dx[0] * dx[2];
       vir[5] =0.5 * v * dx[0] * dx[1];
       virialPerAtom[(*i)*6 + 0] += vir[0];
       virialPerAtom[(*i)*6 + 1] += vir[1];
       virialPerAtom[(*i)*6 + 2] += vir[2];
       virialPerAtom[(*i)*6 + 3] += vir[3];
       virialPerAtom[(*i)*6 + 4] += vir[4];
       virialPerAtom[(*i)*6 + 5] += vir[5];

       virialPerAtom[(*j)*6 + 0] += vir[0];
       virialPerAtom[(*j)*6 + 1] += vir[1];
       virialPerAtom[(*j)*6 + 2] += vir[2];
       virialPerAtom[(*j)*6 + 3] += vir[3];
       virialPerAtom[(*j)*6 + 4] += vir[4];
       virialPerAtom[(*j)*6 + 5] += vir[5];
     }
   */
  
    *ier = KIM_STATUS_OK;
}
