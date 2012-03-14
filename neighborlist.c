#include "KIM_API_status.h"
#include "KIM_API_C.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Define neighborlist structure */
typedef struct
{
   int iteratorId;
   int* NNeighbors;
   int* HalfNNeighbors;
   int* neighborList;
   double* RijList;
} NeighList;

/* 
   The c method to return the neighbor pointer, regardless of if it 
   came from python or the c libraries
*/
int get_neigh(void* kimmdl, int *mode, int *request, int* atom,
              int* numnei, int** nei1atom, double** Rij);

int initialize(void* kimmdl) {
    int status = KIM_API_set_data(kimmdl, "get_neigh", 1, (void*)&get_neigh);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__, "get numberOfParticles", status);
    
    int *natoms = (int*)KIM_API_get_data(kimmdl, "numberOfParticles", &status);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__, "get numberOfParticles", status);

    NeighList *nl = (NeighList*)malloc(sizeof(NeighList));
    nl->NNeighbors = (int*)malloc(sizeof(int)*(*natoms));
    nl->HalfNNeighbors = (int*)malloc(sizeof(int)*(*natoms));
    nl->neighborList = (int*)malloc(sizeof(int)*(*natoms));
    nl->RijList = (double*)malloc(sizeof(int)*(*natoms));
    status = KIM_API_set_data(kimmdl, "neighObject", 1, nl);
 
    return status;
}

void safefree(void *ptr) {
    if (ptr != NULL)
        free(ptr);
}

int free_neigh_object(void* pkim) {
    NeighList *nl;
    int status;

    nl = (NeighList*) KIM_API_get_data(pkim, "neighObject", &status);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

    if (nl==NULL)
        return KIM_STATUS_OK;

    safefree(nl->NNeighbors);
    safefree(nl->HalfNNeighbors);
    safefree(nl->neighborList);
    safefree(nl->RijList);
    safefree(nl);

    status = KIM_API_set_data(pkim, "neighObject", 1, NULL);
    return status;
}


/*
   Use this method to produce a neighbor list in python 
   and set it through the KIM API calls
*/
int set_neigh_object(void* kimmdl, int sz1, int* NNeighbors, 
                                   int sz2, int* HalfNNeighbors, 
                                   int sz3, int* neighborList, 
                                   int sz4, double* RijList) 
{
    //int set_neigh_object(void* kimmdl, int* NNeighbors, int* neighborList, double* RijList) {
    NeighList *nl;
    int status;
    nl = (NeighList*) KIM_API_get_data(kimmdl, "neighObject", &status);
    if (KIM_STATUS_OK == status) {
        free_neigh_object(kimmdl);
    };
    nl = (NeighList *)malloc(sizeof(NeighList));

    nl->iteratorId = -1;
    /*
    nl->NNeighbors = NNeighbors;
    nl->neighborList = neighborList;
    nl->RijList = RijList;
    */
    nl->NNeighbors = (int*)malloc(sizeof(int)*sz1);
    nl->HalfNNeighbors = (int*)malloc(sizeof(int)*sz2);
    nl->neighborList = (int*)malloc(sizeof(int)*sz3);
    nl->RijList = (double*)malloc(sizeof(double)*sz4);
    memcpy(nl->NNeighbors, NNeighbors, sizeof(int)*sz1);
    memcpy(nl->HalfNNeighbors, HalfNNeighbors, sizeof(int)*sz2);
    memcpy(nl->neighborList, neighborList, sizeof(int)*sz3);
    memcpy(nl->RijList, RijList, sizeof(double)*sz4);
    status = KIM_API_set_data(kimmdl, "neighObject", 1, nl);
    return status;
}


// Cluster neighbor locator might need redesign due to these static variables
static int nclusteratoms = 0;
static int ndim = 3;

int get_neigh(void* kimmdl, int *mode, int *request, int* atom,
              int* numnei, int** nei1atom, double** Rij)
{
   /* local variables */
   intptr_t* pkim = *((intptr_t**) kimmdl);
   int atomToReturn;
   int status;
   int* numberOfAtoms;
   char* method;
   NeighList* nl;
   int i;

   int cluster  = 0;
   int ishalf   = 0;

   /* initialize numnei */
   *numnei = 0;

   /* unpack neighbor list object */
   numberOfAtoms = (int*) KIM_API_get_data(pkim, "numberOfParticles", &status);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

   nl = (NeighList*) KIM_API_get_data(pkim, "neighObject", &status);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

   /* figure out the neighbor locator type */
   method = KIM_API_get_NBC_method(pkim, &status);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_NBC_method", status);

   if (!strcmp(method, "CLUSTER"))
        cluster = 1;
   if (!strcmp(method, "NEIGH_PURE_H"))
        ishalf = 1;
   if (!strcmp(method, "MI_OPBC_H"))
        ishalf = 1;

   safefree(method);

   /* check mode and request */
   if (0 == *mode) /* iterator mode */
   {
      if (0 == *request) /* reset iterator */
      {
         (*nl).iteratorId = -1;
         return KIM_STATUS_NEIGH_ITER_INIT_OK;
      }
      else if (1 == *request) /* increment iterator */
      {
         (*nl).iteratorId++;
         if ((*nl).iteratorId >= *numberOfAtoms)
         {
            return KIM_STATUS_NEIGH_ITER_PAST_END;
         }
         else
         {
            atomToReturn = (*nl).iteratorId;
         }
      }
      else /* invalid request value */
      {
         KIM_API_report_error(__LINE__, __FILE__,"Invalid request in get_cluster_neigh", KIM_STATUS_NEIGH_INVALID_REQUEST);
         return KIM_STATUS_NEIGH_INVALID_REQUEST;
      }
   }
   else if (1 == *mode) /* locator mode */
   {
      if ((*request >= *numberOfAtoms) || (*request < 0)) /* invalid id */
      {
         KIM_API_report_error(__LINE__, __FILE__,"Invalid atom ID in get_cluster_neigh", KIM_STATUS_PARTICLE_INVALID_ID);
         return KIM_STATUS_PARTICLE_INVALID_ID;
      }
      else
      {
         atomToReturn = *request;
      }
   }
   else /* invalid mode */
   {
      KIM_API_report_error(__LINE__, __FILE__,"Invalid mode in get_cluster_neigh", KIM_STATUS_NEIGH_INVALID_MODE);
      return KIM_STATUS_NEIGH_INVALID_MODE;
   }

   /* set the returned atom */
   *atom = atomToReturn;

   /* set the returned number of neighbors for the returned atom */
   *numnei = ishalf ? (*nl).HalfNNeighbors[*atom] : (*nl).NNeighbors[*atom];

   if (!cluster)
   {
        /* set the location for the returned neighbor list */
        *nei1atom = (*nl).neighborList;

        /* set the pointer to Rij to appropriate value */
        *Rij = (*nl).RijList;

        for(i = 0; i < *atom; i++) {
            *nei1atom += (*nl).NNeighbors[i];
            *Rij += (*nl).NNeighbors[i]*3;
        }
   }
   else
   {
        /* set the location for the returned neighbor list */
        *nei1atom = &((*nl).neighborList[(*atom)*nclusteratoms]);

        /* set the pointer to Rij to appropriate value */
        *Rij = &((*nl).RijList[(*atom)*ndim*nclusteratoms]);
   }

   return KIM_STATUS_OK;
}



/*===========================================================================
    A standard library of neighbor locators for use with 
    the python KIM API.  All of these were implemented in the
    ex_test_Ar_* files in the v1.0.0 API source
  ==========================================================================*/
int build_neighborlist_allall(void *kimmdl)
{
    int status;
    int* numberOfParticles;
    double* coords;
    double* cutoff;
    char *method;
    NeighList *nl;
 
    int periodic = 0;
    int ishalf   = 0;
    int rijlist  = 0;

    /* get the preferred neighborlist style from the API */
    method = KIM_API_get_NBC_method(kimmdl, &status);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_NBC_method", status);

    if (strcmp(method, "NEIGH_PURE_H") == 0)
        ishalf = 1;
    if (!strcmp(method, "MI_OPBC_F"))
        periodic = 1;
    if (!strcmp(method, "NEIGH_RVEC_F"))
        rijlist = 1;
    if (!strcmp(method, "MI_OPBC_H")){
        ishalf = 1;
        periodic = 1;
    }    
    safefree(method);

    /* get the data necessary for the neighborlist */
    KIM_API_getm_data(kimmdl, &status, 3*3,
            "numberOfParticles",    &numberOfParticles,     1,
            "coordinates",          &coords,                1,
            "cutoff",               &cutoff,                1);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);

    /* get the box side lengths if applicable */
    double *boxSideLengths = 0;
    if (periodic){
        boxSideLengths = (double*)KIM_API_get_data(kimmdl, "boxSideLengths", &status);
        if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__, "boxside", status);
    }

    /* unpack neighbor list object */
    nl = (NeighList*) KIM_API_get_data(kimmdl, "neighObject", &status);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

    int i, j, k, a;
    double dx[3];
    double r2;
    double rcut;
    double rcut2;

    rcut  = *cutoff;
    rcut2 = rcut*rcut;

    double **temp_neigh = (double**)malloc(sizeof(double*)*(*numberOfParticles));
    double **temp_rij   = (double**)malloc(sizeof(double*)*(*numberOfParticles));
    int total = 0;

    for (i=0; i<*numberOfParticles; ++i){
        temp_neigh[i] = (double*)malloc(sizeof(double)*(*numberOfParticles));
        temp_rij[i]   = (double*)malloc(sizeof(double)*(*numberOfParticles)*3);

        a = 0;
        for (j=0; j<*numberOfParticles; ++j){
            r2 = 0.0;

            if (i != j){
                
                for (k=0; k<3; k++){
                    dx[k] = coords[3*j+k] - coords[3*i+k];
                    
                    if (periodic){
                        while (abs(dx[k]) > 0.5*boxSideLengths[k])
                            dx[k] = dx[k] - boxSideLengths[k] * (dx[k] > 0 ? 1: -1);
                    }
    
                    r2 += dx[k]*dx[k];
                }
                
                if (r2 < 4*rcut2){ //4 for a buffer
                    // atom j is a neighbor of atom i
                    if (ishalf && j<i)
                        temp_neigh[i][a] = j; 
                    else 
                        temp_neigh[i][a] = j; 
                    if (rijlist){
                        for (k=0; k<3; k++)
                            temp_rij[i][3*a+k] = dx[k]; 
                    }
                    a++;
                }
            }
        }
        // atom i has a-1 neighbors
        if (ishalf)
            nl->HalfNNeighbors[i] = a;
        else
            nl->NNeighbors[i] = a;
        
        total += a;
    }

    safefree(nl->neighborList);
    safefree(nl->RijList);
    
    nl->neighborList = (int*)malloc(sizeof(int)*total);
    nl->RijList = (double*)malloc(sizeof(double)*total*3);

    int iter = 0;
    for (i=0; i<*numberOfParticles; i++){
        for (j=0; j<nl->NNeighbors[i]; j++){
            nl->neighborList[iter] = temp_neigh[i][j];
            
            for (k=0; k<3; k++)
                nl->RijList[3*iter+k] = temp_rij[i][3*j+k];
            iter++;
        }
        safefree(temp_neigh[i]);
        safefree(temp_rij[i]);
    }

    safefree(temp_neigh);
    safefree(temp_rij);

    return status;
}

int build_neighborlist_cell(void *kimmdl)
{
    return 0;
}
