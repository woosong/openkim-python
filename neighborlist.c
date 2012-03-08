#include "KIM_API_status.h"
#include "KIM_API_C.h"
#include <stdlib.h>
#include <string.h>

/* Define neighborlist structure */
typedef struct
{
   int iteratorId;
   int* NNeighbors;
   int* HalfNNeighbors;
   int* neighborList;
   double* RijList;
} NeighList;

int get_periodic_neigh(void* kimmdl, int *mode, int *request, int* atom,
                       int* numnei, int** nei1atom, double** Rij, int ishalf);

int get_periodic_neigh_half(void* kimmdl, int *mode, int *request, int* atom,
                       int* numnei, int** nei1atom, double** Rij);

int get_periodic_neigh_full(void* kimmdl, int *mode, int *request, int* atom,
                       int* numnei, int** nei1atom, double** Rij);

int get_cluster_neigh(void* kimmdl, int *mode, int *request, int* atom,
                      int* numnei, int** nei1atom, double** Rij);

int set_kim_periodic_full_neigh(void* kimmdl) {
    int status = KIM_API_set_data(kimmdl, "get_full_neigh", 1, (void*)&get_periodic_neigh_full);
    return status;
}
int set_kim_periodic_half_neigh(void* kimmdl) {
    int status = KIM_API_set_data(kimmdl, "get_half_neigh", 1, (void*)&get_periodic_neigh_half);
    return status;
}
int set_kim_cluster_full_neigh(void* kimmdl) {
    int status = KIM_API_set_data(kimmdl, "get_full_neigh", 1, (void*)&get_cluster_neigh);
    return status;
}
int set_kim_cluster_half_neigh(void* kimmdl) {
    int status = KIM_API_set_data(kimmdl, "get_half_neigh", 1, (void*)&get_cluster_neigh);
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

int set_neigh_object(void* kimmdl, int sz1, int* NNeighbors, int sz2, int* HalfNNeighbors, int sz3, int* neighborList, int sz4, double* RijList) {
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

int get_periodic_neigh_half(void* kimmdl, int *mode, int *request, int* atom,
        int* numnei, int** nei1atom, double** Rij) {
   return get_periodic_neigh(kimmdl, mode, request, atom, numnei, nei1atom, Rij,1);
}

int get_periodic_neigh_full(void* kimmdl, int *mode, int *request, int* atom,
        int* numnei, int** nei1atom, double** Rij) {
   return get_periodic_neigh(kimmdl, mode, request, atom, numnei, nei1atom, Rij,0);
}

/* Define prototypes */
int get_periodic_neigh(void* kimmdl, int *mode, int *request, int* atom,
                       int* numnei, int** nei1atom, double** Rij, int ishalf)
{
   /* local variables */
   intptr_t* pkim = *((intptr_t**) kimmdl);
   int atomToReturn;
   int status;
   int* numberOfAtoms;
   NeighList* nl;
    int i;

   /* initialize numnei */
   *numnei = 0;

   /* unpack neighbor list object */
   numberOfAtoms = (int*) KIM_API_get_data(pkim, "numberOfAtoms", &status);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

   nl = (NeighList*) KIM_API_get_data(pkim, "neighObject", &status);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

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
         KIM_API_report_error(__LINE__, __FILE__,"Invalid request in get_periodic_neigh", KIM_STATUS_NEIGH_INVALID_REQUEST);
         return KIM_STATUS_NEIGH_INVALID_REQUEST;
      }
   }
   else if (1 == *mode) /* locator mode */
   {
       // one based fix?
       // use zero based
      if ((*request >= *numberOfAtoms) || (*request < 0)) /* invalid id */
      {
         KIM_API_report_error(__LINE__, __FILE__,"Invalid atom ID in get_periodic_neigh", KIM_STATUS_PARTICLE_INVALID_ID);
         return KIM_STATUS_PARTICLE_INVALID_ID;
      }
      else
      {
         atomToReturn = *request;
      }
   }
   else /* invalid mode */
   {
      KIM_API_report_error(__LINE__, __FILE__,"Invalid mode in get_periodic_neigh", KIM_STATUS_NEIGH_INVALID_MODE);
      return KIM_STATUS_NEIGH_INVALID_MODE;
   }

   /* set the returned atom */
   *atom = atomToReturn;

   /* set the returned number of neighbors for the returned atom */
   *numnei = ishalf ? (*nl).HalfNNeighbors[*atom] : (*nl).NNeighbors[*atom];

   /* set the location for the returned neighbor list */
   *nei1atom = (*nl).neighborList;

   /* set the pointer to Rij to appropriate value */
   *Rij = (*nl).RijList;

   for(i = 0; i < *atom; i++) {
       *nei1atom += (*nl).NNeighbors[i];
       *Rij += (*nl).NNeighbors[i]*3;
    }

   return KIM_STATUS_OK;
}

// Cluster neighbor locator might need redesign due to these static variables
static int nclusteratoms = 0;
static int ndim = 3;

int get_cluster_neigh(void* kimmdl, int *mode, int *request, int* atom,
                      int* numnei, int** nei1atom, double** Rij)
{
   /* local variables */
   intptr_t* pkim = *((intptr_t**) kimmdl);
   int atomToReturn;
   int status;
   int* numberOfAtoms;
   NeighList* nl;

   /* initialize numnei */
   *numnei = 0;

   /* unpack neighbor list object */
   numberOfAtoms = (int*) KIM_API_get_data(pkim, "numberOfAtoms", &status);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

   nl = (NeighList*) KIM_API_get_data(pkim, "neighObject", &status);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

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
   *numnei = (*nl).NNeighbors[*atom];

   /* set the location for the returned neighbor list */
   *nei1atom = &((*nl).neighborList[(*atom)*nclusteratoms]);

   /* set the pointer to Rij to appropriate value */
   *Rij = &((*nl).RijList[(*atom)*ndim*nclusteratoms]);

   return KIM_STATUS_OK;
}
