#include "KIM_API_status.h"
#include "KIM_API_C.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "cvec.h"
#include "neighborlist.h"


/* Define neighborlist structure */
typedef struct
{
   int ready;
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
int nbl_build_neighborlist_cell(void *kimmdl);
int nbl_build_neighborlist_cell_rvec(void *kimmdl);
int nbl_build_neighborlist_allall(void *kimmdl);
int nbl_free_neigh_object(void *kimmdl);

int nbl_get_neigh(void* kimmdl, int *mode, int *request, int* atom,
              int* numnei, int** nei1atom, double** Rij);


int nbl_initialize(void *kimmdl){
    int status;

    nbl_free_neigh_object(kimmdl);               

    // setup a blank neighborlist
    NeighList *nl = (NeighList*)malloc(sizeof(NeighList));
    nl->NNeighbors     = (int*)malloc(sizeof(int)*1);
    nl->HalfNNeighbors = (int*)malloc(sizeof(int)*1);
    nl->neighborList   = (int*)malloc(sizeof(int)*1);
    nl->RijList        = (double*)malloc(sizeof(double)*1);
    nl->iteratorId     = -1;
    nl->ready          = 0;

    status = KIM_API_set_data(kimmdl, "neighObject", 1, nl);
    status = KIM_API_set_data(kimmdl, "get_neigh", 1, (void*)&nbl_get_neigh);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__, "get numberOfParticles", status);

    return status;
}

void safefree(void *ptr) {
    if (ptr != NULL)
        free(ptr);
}

int nbl_free_neigh_object(void* pkim) {
    NeighList *nl;
    int status;

    nl = (NeighList*) KIM_API_get_data(pkim, "neighObject", &status);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

    if (nl==NULL)
        return KIM_STATUS_OK;

    nl->ready = 0;
    safefree(nl->NNeighbors);
    safefree(nl->HalfNNeighbors);
    safefree(nl->neighborList);
    safefree(nl->RijList);

    return status;
}

int nbl_cleanup(void *pkim){
    nbl_free_neigh_object(pkim);

    /* still debating if this is necessary.  valgrind gives errors
       when it is present, so leaving out for now.  also give a 40byte leak
       when not turned on.  the former seems much worse */
/*    int status;
    NeighList *nl = (NeighList*) KIM_API_get_data(pkim, "neighObject", &status);
    safefree(nl);*/
    return KIM_STATUS_OK;
}

/*
   Use this method to produce a neighbor list in python 
   and set it through the KIM API calls
*/

int nbl_set_neigh_object(void* kimmdl, int sz1, int* NNeighbors, 
                                       int sz2, int* HalfNNeighbors, 
                                       int sz3, int* neighborList, 
                                       int sz4, double* RijList) 
{
    //int set_neigh_object(void* kimmdl, int* NNeighbors, int* neighborList, double* RijList) {
    NeighList *nl;
    int status;

    nl = (NeighList*) KIM_API_get_data(kimmdl, "neighObject", &status);
    if (KIM_STATUS_OK == status) {
        nbl_free_neigh_object(kimmdl);
    };

    if (nl == NULL)
        nbl_initialize(kimmdl);

    nl->ready          = 1; 
    nl->iteratorId     = -1;
    nl->NNeighbors     = (int*)malloc(sizeof(int)*sz1);
    nl->HalfNNeighbors = (int*)malloc(sizeof(int)*sz2);
    nl->neighborList   = (int*)malloc(sizeof(int)*sz3);
    nl->RijList        = (double*)malloc(sizeof(double)*sz4);
    
    memcpy(nl->NNeighbors, NNeighbors, sizeof(int)*sz1);
    memcpy(nl->HalfNNeighbors, HalfNNeighbors, sizeof(int)*sz2);
    memcpy(nl->neighborList, neighborList, sizeof(int)*sz3);
    memcpy(nl->RijList, RijList, sizeof(double)*sz4);
   
    return status;
}


// Cluster neighbor locator might need redesign due to these static variables
static int nclusteratoms = 0;
static int ndim = 3;

int nbl_get_neigh(void* kimmdl, int *mode, int *request, int* atom,
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

   /* if the user got here by mistake, correct it */
   if (nl == NULL)
        nbl_initialize(kimmdl);
   if (!nl->ready)
        nbl_build_neighborlist(kimmdl);

   /* figure out the neighbor locator type */
   method = KIM_API_get_NBC_method(pkim, &status);
   if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_NBC_method", status);
   if (!strcmp(method, "CLUSTER"))      cluster = 1;
   if (!strcmp(method, "NEIGH_PURE_H")) ishalf = 1;
   if (!strcmp(method, "MI_OPBC_H"))    ishalf = 1;
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
    the python KIM API.  
  ==========================================================================*/
double cellf[9] = {1.0, 0.0, 0.0, /*|*/ 0.0, 1.0, 0.0, /*|*/ 0.0, 0.0, 1.0};
double cellr[9] = {1.0, 0.0, 0.0, /*|*/ 0.0, 1.0, 0.0, /*|*/ 0.0, 0.0, 1.0};
int    pbc[3]   = {0, 0, 0};
int    ortho[3] = {1, 1, 1};
int    init     = 0;

inline double det(double a11, double a12, double a21, double a22){
    return (a11*a22) - (a12*a21);
}

inline double dot3(double *a, double *b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
} 

inline double det3(double *mat){
    double dd = (mat[0]*mat[4]*mat[8] - mat[0]*mat[5]*mat[7] -
                 mat[1]*mat[3]*mat[8] + mat[1]*mat[5]*mat[6] +
                 mat[2]*mat[3]*mat[7] - mat[2]*mat[4]*mat[6]);
    return dd;
}

void inverse(double *mat, double *inv){
    inv[0] = det(mat[4], mat[5], mat[7], mat[8]);
    inv[1] = det(mat[2], mat[1], mat[8], mat[7]);
    inv[2] = det(mat[1], mat[2], mat[4], mat[5]);
    inv[3] = det(mat[5], mat[3], mat[8], mat[6]);
    inv[4] = det(mat[0], mat[3], mat[6], mat[8]);
    inv[5] = det(mat[2], mat[0], mat[5], mat[3]);
    inv[6] = det(mat[3], mat[4], mat[6], mat[7]);
    inv[7] = det(mat[1], mat[0], mat[7], mat[6]);
    inv[8] = det(mat[0], mat[1], mat[3], mat[4]);

    int i = 0;
    double dd = det3(mat);
    for (i=0; i<9; i++)
        inv[i] /= dd;
}

void transpose(double *m){
    int i,j;
    double trans[9];
    memcpy(trans, m, sizeof(double)*9);
    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            m[3*i+j] = trans[3*j+i];
}


void is_orthogonal(double *mat, int *ortho){
    int i, j;
    
    for (i=0; i<3; i++){
        int or = 1;
        for (j=0; j<3; j++){
            if (abs(mat[3*i+j]) > 1e-8 && i != j)
                or = 0;
        }
        ortho[i] = or;
    }
}

void printmat(double *m){
    int i,j;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            printf(" %f ", m[3*i+j]);
        }
        printf("\n");
    }
}

int nbl_set_cell(int S1, double* Cell, int S2, int* PBC){
    if (S1 != 9) printf("WARNING: cell size greater than 9.  Taking first 9 elements\n");
    if (S2 != 3) printf("WARNING: pbc size greater than 3.  Taking first 3 elements\n");
 
    memcpy(pbc, PBC, sizeof(int)*3);
    memcpy(cellf, Cell, sizeof(double)*9);

    is_orthogonal(cellf, ortho);
    transpose(cellf);    
    inverse(cellf, cellr);        
    
    init = 1;
    return 0;
}

inline double fmod1(double a){
    double m = a - ((int)a);
    if (m >= -1e-8)
        return m;
    return m + 1.0;
}

void transform(double coords[3], double cell[9], double out[3]){
   int i;
    for (i=0; i<3; i++)
        out[i] = dot3(&cell[3*i], coords);
}

/* =====================================================
   call this, and it will decide which to use
   ====================================================*/
int nbl_build_neighborlist(void *kimmdl){
    int rvec;
    int status;
    char *method;

    method = KIM_API_get_NBC_method(kimmdl, &status);
    rvec = strcmp(method,"NEIGH_RVEC_F");
    if (KIM_STATUS_OK > status)
        KIM_API_report_error(__LINE__, __FILE__,"get_NBC_method", status);
    safefree(method);

    /* if the user didn't specify the box, don't die, just
       use the all neighbor list, can't assume 
       anything about periodic conditions or cell */
    if (init == 0)
        return nbl_build_neighborlist_allall(kimmdl);

    /* otherwise, the cell works for all types */
    if (!rvec)
        return nbl_build_neighborlist_cell_rvec(kimmdl);
    return nbl_build_neighborlist_cell(kimmdl);
}


/*======================================================
  the simple all-all neighbor list
  works for all neighborlists including rvec
  ======================================================*/
int nbl_build_neighborlist_allall(void *kimmdl)
{
    int status;
    int* numberOfParticles;
    double* coords;
    double* ncoords;
    double* cutoff;
    char *method;
    double *boxSideLengths;
    NeighList *nl;

    double tt[3];
 
    int ishalf   = 0;
    int rijlist  = 0;
    int mi_opbc  = 0;
    
    /* get the preferred neighborlist style from the API */
    method = KIM_API_get_NBC_method(kimmdl, &status);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_NBC_method", status);
    if (!strcmp(method, "NEIGH_PURE_H"))  ishalf  = 1;
    if (!strcmp(method, "MI_OPBC_F"))     mi_opbc = 1;
    if (!strcmp(method, "NEIGH_RVEC_F"))  rijlist = 1;
    if (!strcmp(method, "MI_OPBC_H")){    ishalf  = 1;  mi_opbc = 1; }    
    safefree(method);

    /* get the data necessary for the neighborlist */
    KIM_API_getm_data(kimmdl, &status, 4*3,
            "numberOfParticles",    &numberOfParticles,     1,
            "coordinates",          &coords,                1,
            "cutoff",               &cutoff,                1,
            "boxSideLengths",       &boxSideLengths,        mi_opbc);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);

    if (mi_opbc){
        boxSideLengths[0] = cellf[0];
        boxSideLengths[1] = cellf[4];  
        boxSideLengths[2] = cellf[8];
    }

    nl = (NeighList*) KIM_API_get_data(kimmdl, "neighObject", &status);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

    /* if the user got here by mistake, correct it */
    if (nl == NULL)
        nbl_initialize(kimmdl);
 
    /* reset the NeighList object to refill */
    nbl_free_neigh_object(kimmdl);
    
    nl->iteratorId     = -1;
    nl->NNeighbors     = (int*)malloc(sizeof(int)*(*numberOfParticles));
    nl->HalfNNeighbors = (int*)malloc(sizeof(int)*(*numberOfParticles));

    /* begin the actual neighborlist calculation */
    int i, j, k;
    int neighbors;
    double dx[3];
    double ds[3];
    double r2;
    double rcut;
    double rcut2;

    rcut  = *cutoff;
    rcut2 = rcut*rcut;

    ncoords = (double*)malloc(sizeof(double)*(*numberOfParticles)*3);
    for (i=0; i<*numberOfParticles; i++){
        transform(&coords[i*3], cellr, &ncoords[i*3]);
        for (j=0; j<3; j++){
            if (pbc[j] == 1)  
                ncoords[3*i+j] = fmod1(ncoords[3*i+j]);
        }
    }

    cvec *temp_neigh = (cvec*)malloc(sizeof(cvec));
    dvec *temp_rij   = (dvec*)malloc(sizeof(dvec));
    cvec_init(temp_neigh, 4);
    dvec_init(temp_rij,   4);

    int total = 0;
 
    for (i=0; i<*numberOfParticles; i++){
        neighbors = 0;
        for (j=(i+1)*ishalf; j<*numberOfParticles; j++){

            if (rijlist){
                for (tt[0]=-1*pbc[0]; tt[0]<=1*pbc[0]; tt[0]++){
                for (tt[1]=-1*pbc[1]; tt[1]<=1*pbc[1]; tt[1]++){
                for (tt[2]=-1*pbc[2]; tt[2]<=1*pbc[2]; tt[2]++){

                for (k=0; k<3; k++){
                    dx[k] = ncoords[3*j+k] - ncoords[3*i+k];
                    dx[k] += tt[k];
                }
               
                transform(dx, cellf, ds); 

                r2 = 0.0;
                for (k=0; k<3; k++)    
                    r2 += ds[k]*ds[k];
                
                if (r2 > 1e-10 && r2 < rcut2){ //4 for a buffer
                    cvec_insert_back(temp_neigh, j); 
                    if (rijlist){
                        for (k=0; k<3; k++)
                            dvec_insert_back(temp_rij, ds[k]); 
                    }
                    neighbors++;
                }
                } } }
            }
            else {
                for (k=0; k<3; k++){
                    dx[k] = ncoords[3*j+k] - ncoords[3*i+k];
                    if (pbc[k] != 0)
                        dx[k] += fmod1(dx[k]); 
                }
               
                transform(dx, cellf, ds); 

                r2 = 0.0;
                for (k=0; k<3; k++)    
                    r2 += ds[k]*ds[k];
                
                if (r2 > 1e-10 && r2 < rcut2){ //4 for a buffer
                    cvec_insert_back(temp_neigh, j); 
                    if (rijlist){
                        for (k=0; k<3; k++)
                            dvec_insert_back(temp_rij, ds[k]); 
                    }
                    neighbors++;
                }
            }
        }
        nl->HalfNNeighbors[i] = neighbors;
        nl->NNeighbors[i] = neighbors;
        
        total += neighbors;
    }

    nl->neighborList = (int*)malloc(sizeof(int)*total);
    nl->RijList      = (double*)malloc(sizeof(double)*total*3);

    memcpy(nl->neighborList, temp_neigh->array, sizeof(int)*total);
    if (rijlist != 0) memcpy(nl->RijList, temp_rij->array, sizeof(double)*total*3);

    /* cleanup all the memory usage */
    cvec_destroy(temp_neigh);
    dvec_destroy(temp_rij);
    safefree(temp_neigh);
    safefree(temp_rij);
    safefree(ncoords);
    
    nl->ready = 1;
    return status;
}


//========================================================================
// the more advanced cell neighbor locator
// this version works for rvec only - it would be inefficient 
// to combine this with the other
//========================================================================
inline void coords_to_index(double *x, int *size, int *index, double *max, double *min){   
    int i;
    for (i=0; i<3; i++)
        index[i] = (int)(((x[i]-min[i])/(max[i]-min[i]) - 1e-14) * size[i]);
}

inline int mod_rvec(int a, int b, int p, int *image){
    *image = 1;
    if (b==0) {if (a==0) *image=0; return 0;}
    if (p != 0){
        if (a>b)  return a-b-1;
        if (a<0)  return a+b+1;
    } else {
        if (a>b)  return b;
        if (a<0)  return 0;
    }
    *image = 0;
    return a;
}


int nbl_build_neighborlist_cell_rvec(void *kimmdl)
{
    int i, j, k;
    int status;
    int* numberOfParticles;
    double* coords;
    double* ncoords; 
    double* cutoff;
    NeighList *nl;
 
    /* get the data necessary for the neighborlist */
    KIM_API_getm_data(kimmdl, &status, 3*3,
            "numberOfParticles",    &numberOfParticles,     1,
            "coordinates",          &coords,                1,
            "cutoff",               &cutoff,                1);
    if (KIM_STATUS_OK > status) 
        KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);

    nl = (NeighList*) KIM_API_get_data(kimmdl, "neighObject", &status);
    if (KIM_STATUS_OK > status) 
        KIM_API_report_error(__LINE__, __FILE__,"get_data", status);

    /* if the user got here by mistake, correct it */
    if (nl == NULL)  nbl_initialize(kimmdl);

    /* reset the NeighList object to refill */
    nbl_free_neigh_object(kimmdl);
    nl->iteratorId     = -1;
    nl->NNeighbors     = (int*)malloc(sizeof(int)*(*numberOfParticles));
    nl->HalfNNeighbors = (int*)malloc(sizeof(int)*(*numberOfParticles));

    /* begin the actual neighborlist calculation */
    int copies[3];
    double R, R2, dist;

    int box = ortho[0]*ortho[1]*ortho[2];
    R  = *cutoff;
    R2 = (*cutoff)*(*cutoff);

    /* convert the coordinates to [0.0,1.0] for peiodic sides
       and find the side length of those that are not */
    double min[3] = {0.0,0.0,0.0};
    double max[3] = {1.0,1.0,1.0};
    ncoords = (double*)malloc(sizeof(double)*(*numberOfParticles)*3);
    for (i=0; i<*numberOfParticles; i++){
        transform(&coords[i*3], cellr, &ncoords[i*3]);
        for (j=0; j<3; j++){
            if (pbc[j] == 0){
                if (max[j] < ncoords[3*i+j]) max[j] = ncoords[3*i+j];
                if (min[j] > ncoords[3*i+j]) min[j] = ncoords[3*i+j];
            }
            else 
                ncoords[3*i+j] = fmod1(ncoords[3*i+j]);
        }
    }

    /* all of this to make the perfect box size for the cells */
    /* create unit vectors along the cell, column vectors */
    double span[9];
    for (i=0; i<3; i++){
        double len = 0.0;
        for (j=0; j<3; j++)
            len += cellf[3*j+i]*cellf[3*j+i];
        for (j=0; j<3; j++) 
            span[3*j+i] = cellf[3*j+i]/sqrt(len);
    }

    /* the smallest box needs to include the slanty-ness */
    double factor[3];
    for (i=0; i<3; i++){
        double tcos, tsin;
        double sinmin = 1.0;
        for (j=0; j<3; j++){
            if (i != j){
                tcos = span[3*0+i]*span[3*0+j] 
                     + span[3*1+i]*span[3*1+j] 
                     + span[3*2+i]*span[3*2+j];
                tsin = sqrt(1 - tcos*tcos);
                if (tsin < sinmin) sinmin = tsin;
            }
        }
        factor[i] = sinmin;
    }

    /* make the cell box */
    double length[3];
    int size[3];
    int size_total = 1;
    
    /* figure out how big the cells need to be */
    for (i=0; i<3; i++){
        double rtemp = sqrt(cellf[3*0+i]*cellf[3*0+i] 
                          + cellf[3*1+i]*cellf[3*1+i] 
                          + cellf[3*2+i]*cellf[3*2+i]);
        if (box == 0) rtemp *= factor[i]/sqrt(3.0);
        size[i] = (int)((rtemp / R) * (max[i]-min[i]));
        copies[i] = 1;

        /* if the cell is smaller than cutoff, make sure
           there are more copies to be periodic */
        if (size[i] <= 0){
            size[i] = 1;
            copies[i] = (int)(R/rtemp) + 1;
        }

        size_total *= size[i];
    }

    if (size_total > 1000000000){
        printf("%i x %i x %i ?! You thinks me a bit ambitious.\n", size[0], size[1], size[2]);
        printf("Check if you have particles at 1e10 or try implementing a hash cell list.\n");
        printf("Attempting to proceed...\n");
    }

    /* initialize the cells for the neighbors */ 
    cvec **cells = (cvec**)malloc(sizeof(cvec*)*size_total);
    for (i=0; i<size_total; i++){
        cells[i] = (cvec*)malloc(sizeof(cvec));
        cvec_init(cells[i], 128);
    }

    /* convert all particles to cells */
    int index[3];
    for (i=0; i<*numberOfParticles; i++){
        coords_to_index(&ncoords[3*i], size, index, max, min);
        int t = index[0] + index[1]*size[0] + index[2]*size[0]*size[1];
        cvec_insert_back(cells[t], i);
    }

    cvec *temp_neigh = (cvec*)malloc(sizeof(cvec));
    dvec *temp_rij   = (dvec*)malloc(sizeof(dvec));
    cvec_init(temp_neigh, (*numberOfParticles));
    dvec_init(temp_rij,   (*numberOfParticles));
   
    int total = 0;
    int neighbors;

    int tt[3];
    int tix[3];
    int image[3];

    double dx[3];
    double ds[3];

    for (i=0; i<*numberOfParticles; i++){
        neighbors = 0;
        /* translate back to a cell id */
        coords_to_index(&ncoords[3*i], size, index, max, min);

        /* loop over the neighboring cells, more than one if the box is small */
        for (tt[0]=-copies[0]; tt[0]<=copies[0]; tt[0]++){
        for (tt[1]=-copies[1]; tt[1]<=copies[1]; tt[1]++){
        for (tt[2]=-copies[2]; tt[2]<=copies[2]; tt[2]++){
            int goodcell = 1;    
            /* wrap the cell around as necessary 
               and if it isn't periodic, don't count the wrap around */
            for (j=0; j<3; j++){
                tix[j] = mod_rvec(index[j]+tt[j],size[j]-1,pbc[j],&image[j]);
                if (pbc[j] < image[j])
                    goodcell=0;
            }

            /* GTG? */
            if (goodcell){
                int ind = tix[0] + tix[1]*size[0] + tix[2]*size[0]*size[1];
                cvec *cell = cells[ind];

                /* for every particle in that cell */
                for (j=0; j<cell->elems; j++){
                    int n = cvec_at(cell, j); 

                    /* find the distance vecs in our 1x1x1 box */
                    dist = 0.0;
                    for (k=0; k<3; k++){
                        dx[k] = ncoords[3*n+k] - ncoords[3*i+k];
                
                        /* subtract off the image number in [0,1] box */
                        if (image[k])
                            dx[k] += tt[k];
                    }
       
                    /* take back to the curvy cell */
                    transform(dx, cellf, ds);

                    /* find the distance in real space */
                    for (k=0; k<3; k++) 
                        dist += ds[k]*ds[k];

                    /* do we have a neighbor and are not ourselves */
                    if (dist > 1e-10 && dist < R2){ 
                        cvec_insert_back(temp_neigh, n);
                        for (k=0; k<3; k++)
                            dvec_insert_back(temp_rij, ds[k]);
                        neighbors++;
                    }
                }
            }
        } } } 
        /* update that particles list */
        nl->HalfNNeighbors[i] = neighbors;
        nl->NNeighbors[i] = neighbors;
        
        total += neighbors;
    }

    nl->neighborList = (int*)malloc(sizeof(int)*total);
    nl->RijList      = (double*)malloc(sizeof(double)*total*3);

    memcpy(nl->neighborList, temp_neigh->array, sizeof(int)*total);
    memcpy(nl->RijList, temp_rij->array, sizeof(double)*total*3);

    /* cleanup all the memory usage */
    cvec_destroy(temp_neigh);
    dvec_destroy(temp_rij);
    for (i=0; i<size_total; i++){
        cvec_destroy(cells[i]); 
        safefree(cells[i]);
    }
    safefree(cells);
    safefree(temp_neigh);
    safefree(temp_rij);
    safefree(ncoords);

    nl->ready = 1; 
    return status;
}


//========================================================================
// the less advanced cell neighbor locator
// works for everything that is not rvec
//========================================================================
inline int mod(int a, int b, int p){
    if (p != 0){
        if (b==0) return 0;
        if (a>b)  return a-b-1;
        if (a<0)  return a+b+1;
    } else {
        if (b==0) return 0;
        if (a>b)  return b;
        if (a<0)  return 0;
    }
    return a;
}

cvec *get_relevant_cells(int *index, int *size, int *hash){
    cvec *list = (cvec*)malloc(sizeof(cvec));
    cvec_init(list, 27);
    int i,j,k;
 
    for (i=index[0]-1; i<=index[0]+1; i++){
    for (j=index[1]-1; j<=index[1]+1; j++){
    for (k=index[2]-1; k<=index[2]+1; k++){
        int ind = mod(i,size[0]-1,pbc[0]) 
                + mod(j,size[1]-1,pbc[1])*size[0]
                + mod(k,size[2]-1,pbc[2])*size[0]*size[1];
        if (hash[ind] == 0){ hash[ind] = 1; cvec_insert_back(list, ind); }
    } } }
 
    return list;       
}


int nbl_build_neighborlist_cell(void *kimmdl)
{
    int status;
    int* numberOfParticles;
    double* coords;
    double* ncoords; 
    double* cutoff;
    char *method;
    double *boxSideLengths;
    NeighList *nl;
 
    int ishalf   = 0;
    int rijlist  = 0;
    int mi_opbc  = 0;
    
    /* get the preferred neighborlist style from the API */
    method = KIM_API_get_NBC_method(kimmdl, &status);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_NBC_method", status);
    if (!strcmp(method, "NEIGH_PURE_H"))  ishalf  = 1;
    if (!strcmp(method, "MI_OPBC_F"))     mi_opbc = 1;
    if (!strcmp(method, "NEIGH_RVEC_F"))  rijlist = 1;
    if (!strcmp(method, "MI_OPBC_H")){    ishalf  = 1;   mi_opbc = 1; }    
    safefree(method);

    /* get the data necessary for the neighborlist */
    KIM_API_getm_data(kimmdl, &status, 4*3,
            "numberOfParticles",    &numberOfParticles,     1,
            "coordinates",          &coords,                1,
            "cutoff",               &cutoff,                1,
            "boxSideLengths",       &boxSideLengths,        mi_opbc);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);

    if (mi_opbc){
        boxSideLengths[0] = cellf[0];
        boxSideLengths[1] = cellf[4];  
        boxSideLengths[2] = cellf[8];
    }

    nl = (NeighList*) KIM_API_get_data(kimmdl, "neighObject", &status);
    if (KIM_STATUS_OK > status) KIM_API_report_error(__LINE__, __FILE__,"get_data", status);
    /* if the user got here by mistake, correct it */
    if (nl == NULL)
        nbl_initialize(kimmdl);
 
    /* reset the NeighList object to refill */
    nbl_free_neigh_object(kimmdl);
    nl->iteratorId     = -1;
    nl->NNeighbors     = (int*)malloc(sizeof(int)*(*numberOfParticles));
    nl->HalfNNeighbors = (int*)malloc(sizeof(int)*(*numberOfParticles));

    /* begin the actual neighborlist calculation */
    int i, j, k, l;
    int neighbors;
    double dx[3];
    double ds[3];
    double R, R2, dist;

    int box = ortho[0]*ortho[1]*ortho[2];
    R  = *cutoff;
    R2 = (*cutoff)*(*cutoff);

    double min[3] = {0.0,0.0,0.0};
    double max[3] = {1.0,1.0,1.0};
    ncoords = (double*)malloc(sizeof(double)*(*numberOfParticles)*3);
    for (i=0; i<*numberOfParticles; i++){
        transform(&coords[i*3], cellr, &ncoords[i*3]);
        for (j=0; j<3; j++){
            if (pbc[j] == 0){
                if (max[j] < ncoords[3*i+j]) max[j] = ncoords[3*i+j];
                if (min[j] > ncoords[3*i+j]) min[j] = ncoords[3*i+j];
            }
            else 
                ncoords[3*i+j] = fmod1(ncoords[3*i+j]);
        }
    }
 
    /* all of this to make the perfect box size for the cells */
    // create unit vectors along the cell, column vectors
    double span[9];
    for (i=0; i<3; i++){
        double len = 0.0;
        for (j=0; j<3; j++)
            len += cellf[3*j+i]*cellf[3*j+i];
        for (j=0; j<3; j++) 
            span[3*j+i] = cellf[3*j+i]/sqrt(len);
    }

    double factor[3];
    for (i=0; i<3; i++){
        double tcos, tsin;
        double sinmin = 1.0;
        for (j=0; j<3; j++){
            if (i != j){
                tcos = span[3*0+i]*span[3*0+j] + span[3*1+i]*span[3*1+j] + span[3*2+i]*span[3*2+j];
                tsin = sqrt(1 - tcos*tcos);
                if (tsin < sinmin) sinmin = tsin;
            }
        }
        factor[i] = sinmin;
    }
    
    /* make the cell box */
    int size[3];
    int size_total = 1;
    
    for (i=0; i<3; i++){
        double rtemp = sqrt(cellf[3*0+i]*cellf[3*0+i] 
                          + cellf[3*1+i]*cellf[3*1+i] 
                          + cellf[3*2+i]*cellf[3*2+i]);

        if (box == 0) rtemp *= factor[i]/sqrt(3.0);
        size[i] = (int)((rtemp / R)*(max[i]-min[i])) + 1;
        size_total *= size[i];
    }
   
    if (size_total > 1000000000){
        printf("%i x %i x %i ?! You thinks me a bit ambitious.\n", size[0], size[1], size[2]);
        printf("Check if you have particles at 1e10 or try implementing a hash cell list.\n");
        printf("Attempting to proceed...\n");
    }
 
    cvec **cells = (cvec**)malloc(sizeof(cvec*)*size_total);
    int *hash    =   (int*)malloc(sizeof(int)*size_total);
    for (i=0; i<size_total; i++){
        cells[i] = (cvec*)malloc(sizeof(cvec));
        hash[i]  = 0;
        cvec_init(cells[i], 4);
    }

    int index[3];
    for (i=0; i<*numberOfParticles; i++){
        coords_to_index(&ncoords[3*i], size, index, max, min);
        cvec_insert_back(cells[index[0] + index[1]*size[0] + index[2]*size[0]*size[1]], i);
    }

    cvec *temp_neigh = (cvec*)malloc(sizeof(cvec));
    dvec *temp_rij   = (dvec*)malloc(sizeof(dvec));
    cvec_init(temp_neigh, 4);
    dvec_init(temp_rij,   4);
   
    int total = 0;

    for (i=0; i<*numberOfParticles; i++){
        neighbors = 0;
        coords_to_index(&ncoords[3*i], size, index, max, min);
        cvec *ncells = get_relevant_cells(index, size, hash);

        for (l=0; l<ncells->elems; l++){
            int cellindex = cvec_at(ncells, l);
            cvec *cell = cells[cellindex];
            hash[cellindex] = 0;

            for (j=0; j<cell->elems; j++){
                int n = cvec_at(cell, j); 
                if ((i != n) && ((ishalf==0) || ((ishalf!=0) && (n > i)))){
                    dist = 0.0;
                    for (k=0; k<3; k++){
                        dx[k] = ncoords[3*n+k] - ncoords[3*i+k];
                    
                        if (pbc[k] != 0)
                            dx[k] = fmod1(dx[k]);
                    }
               
                    transform(dx, cellf, ds);
                    for (k=0; k<3; k++) 
                        dist += ds[k]*ds[k];

                    if (dist < R2){ //4 for a buffer
                        // atom j is a neighbor of atom i
                        cvec_insert_back(temp_neigh, n);
                        if (rijlist){
                            for (k=0; k<3; k++)
                                dvec_insert_back(temp_rij, ds[k]);
                        }
                        neighbors++;
                    }
                }
            }
        }
        nl->HalfNNeighbors[i] = neighbors;
        nl->NNeighbors[i] = neighbors;
        
        total += neighbors;

        cvec_destroy(ncells);
        safefree(ncells);
    }

    nl->neighborList = (int*)malloc(sizeof(int)*total);
    nl->RijList      = (double*)malloc(sizeof(double)*total*3);

    memcpy(nl->neighborList, temp_neigh->array, sizeof(int)*total);
    if (rijlist)
        memcpy(nl->RijList,      temp_rij->array,   sizeof(double)*total*3);

    /* cleanup all the memory usage */
    cvec_destroy(temp_neigh);
    dvec_destroy(temp_rij);
    for (i=0; i<size_total; i++){
        cvec_destroy(cells[i]); 
        safefree(cells[i]);
    }
    safefree(cells);
    safefree(temp_neigh);
    safefree(temp_rij);
    safefree(ncoords);
    safefree(hash);

    nl->ready = 1; 
    return status;
}

