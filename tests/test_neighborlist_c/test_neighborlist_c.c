/*
    This is a copy of the KIM API example ex_test_Ar_free_cluster_CLUSTER_C
    the changed sections as noted, and incorporate the neighborlist from
    up a directory
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "KIM_API_C.h"
#include "KIM_API_status.h"
#include "../../neighborlist.h"

#define FCCSPACING    5.260
#define NCELLSPERSIDE 2
#define DIM           3
#define ATYPES        1
#define NCLUSTERATOMS (4*(NCELLSPERSIDE*NCELLSPERSIDE*NCELLSPERSIDE) + 6*(NCELLSPERSIDE*NCELLSPERSIDE) + 3*(NCELLSPERSIDE) + 1)

/* Define prototypes */
static void create_FCC_configuration(double FCCspacing, int nCellsPerSide, int periodic,
                                     double *coords, int *MiddleAtomId);

double mcell[9] = {2e1,0,0,0,2e1,0,0,0,2e1};
int mpbc[3] = {1,1,1};

/* Main program */
int main(int argc, char* argv[])
{

   /* Local variable declarations */
   int i;


   /* KIM variable declarations */
   char testname[] = "test_neighborlist_c";
   char modelname[] = "ex_model_Ar_P_LJ";
   void* pkim;
   int status;
   int partcl_type_code;

   /* model inputs */
   int* numberOfParticles;
   int* numberParticleTypes;
   int* particleTypes;
   double* coords;
   /* model outputs */
   double* cutoff;
   double* energy;
   int middleDum;

   int lenteststring = 80*100;
   char kimfile[1024];
   char teststring[lenteststring];

   sprintf(kimfile, "%s.kim", testname);
   
   FILE *fkim = fopen(kimfile, "r");
   if (fkim)
       fread(teststring, sizeof(char), lenteststring, fkim); 
   else {
       printf("could not open %s\n", kimfile);
       exit(1);
   }
   fclose(fkim);

   /* Initialize the KIM Model */
   status = KIM_API_string_init(&pkim, teststring, modelname);
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_init", status);
      exit(1);
   }

   /* Allocate memory via the KIM system */
   KIM_API_allocate(pkim, NCLUSTERATOMS, ATYPES, &status);
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_allocate", status);
      exit(1);
   }

   /* call Model's init routine */
   status = KIM_API_model_init(pkim);
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_init", status);
      exit(1);
   }

   /* initialize the neighborlist */
   nbl_initialize(pkim);
   nbl_set_cell(9, mcell, 3, mpbc);

   /* Unpack data from KIM object */
   KIM_API_getm_data(pkim, &status, 6*3,
                     "numberOfParticles",   &numberOfParticles,   1,
                     "numberParticleTypes", &numberParticleTypes, 1,
                     "particleTypes",       &particleTypes,       1,
                     "coordinates",         &coords,              1,
                     "cutoff",              &cutoff,              1,
                     "energy",              &energy,              1);
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", status);
      exit(1);
   }

   /* Set values */
   *numberOfParticles   = NCLUSTERATOMS;
   *numberParticleTypes = ATYPES;
   partcl_type_code = KIM_API_get_partcl_type_code(pkim, "Ar", &status);
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_partcl_type_code", status);
      exit(1);
   }
   for (i = 0; i < *numberOfParticles; ++i)
   {
      particleTypes[i] = partcl_type_code;
   }

   /* set up the cluster atom positions */
   create_FCC_configuration(FCCSPACING, NCELLSPERSIDE, 0, coords, &middleDum);

   /* actually build the neighborlist */
   nbl_build_neighborlist(pkim);

   /* Call model compute */
   KIM_API_model_compute(pkim, &status);
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_compute", status);
      exit(1);
   }

   /* print results to screen */
   printf("--------------------------------------------------------------------------------\n");
   printf("This is Test          : %s\n",testname);
   printf("Results for KIM Model : %s\n",modelname);
   printf("Energy        = %25.15e\n", *energy);

   /* don't forget to destroy and deallocate */
   nbl_cleanup(pkim);

   KIM_API_model_destroy(pkim, &status);
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_destory", status);
      exit(1);
   }
   KIM_API_free(&pkim, &status);
   if (KIM_STATUS_OK > status)
   {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_free", status);
      exit(1);
   }

   /* everything is great */
   return 0;
}

/******************************************************************************
 *
 * create_FCC_configuration function
 *
 *  creates a cubic configuration of FCC atoms with lattice spacing
 *  `FCCspacing' and `nCellsPerSide' cells along each direction.
 *
 *  With periodic==0. this will result in a total number of atoms equal to
 *  4*(nCellsPerSide)**3 + 6*(nCellsPerSide)**2 + 3*(nCellsPerSide) + 1
 *
 *  With periodic==1 this will result in a total number of atoms equal to
 *  4*(nCellsPerSide)**3
 *
 *  Returns the Id of the atom situated in the middle of the configuration
 *  (this atom will have the most neighbors.)
 *
 ******************************************************************************/
static void create_FCC_configuration(double FCCspacing, int nCellsPerSide, int periodic,
                                     double *coords, int* MiddleAtomId)
{
   /* local variables */
   double FCCshifts[4][3];
   double latVec[3];
   int a;
   int i;
   int j;
   int k;
   int m;
   int n;

   /* create a cubic FCC cluster */
   FCCshifts[0][0] = 0.0;            FCCshifts[0][1] = 0.0;            FCCshifts[0][2] = 0.0;
   FCCshifts[1][0] = 0.5*FCCspacing; FCCshifts[1][1] = 0.5*FCCspacing; FCCshifts[1][2] = 0.0;
   FCCshifts[2][0] = 0.5*FCCspacing; FCCshifts[2][1] = 0.0;            FCCshifts[2][2] = 0.5*FCCspacing;
   FCCshifts[3][0] = 0.0;            FCCshifts[3][1] = 0.5*FCCspacing; FCCshifts[3][2] = 0.5*FCCspacing;

   *MiddleAtomId = 0; /* Always put middle atom as #0 */
   a = 1;            /* leave space for middle atom as atom #0 */
   for (i = 0; i < nCellsPerSide; ++i)
   {
      latVec[0] = ((double) i)*FCCspacing;
      for (j = 0; j < nCellsPerSide; ++j)
      {
         latVec[1] = ((double) j)*FCCspacing;
         for (k = 0; k < nCellsPerSide; ++k)
         {
            latVec[2] = ((double) k)*FCCspacing;
            for (m = 0; m < 4; ++m)
            {
               for (n = 0; n < DIM; ++n)
               {
                  coords[a*DIM + n] = latVec[n] + FCCshifts[m][n];
               }
               if ((nCellsPerSide/2 == i) && (nCellsPerSide/2 == j) &&
                   (nCellsPerSide/2 == k) && (1 == m))
               {
                  /* put middle atom as atom #0 */
                  for (n = 0; n < DIM; ++n)
                  {
                     coords[n] = latVec[n] + FCCshifts[m][n];
                  }
                  a--;
               }
               a++;
            }
         }
         if (!periodic)
         {
            /* add in the remaining three faces */
            /* pos-x face */
            latVec[0] = NCELLSPERSIDE*FCCspacing;
            latVec[1] = ((double) i)*FCCspacing;
            latVec[2] = ((double) j)*FCCspacing;
            for (n = 0; n < DIM; ++n)
            {
               coords[a*DIM + n] = latVec[n];
            }
            a++;
            for (n = 0; n < DIM; ++n)
            {
               coords[a*DIM + n] = latVec[n] + FCCshifts[3][n];
            }
            a++;
            /* pos-y face */
            latVec[0] = ((double) i)*FCCspacing;
            latVec[1] = NCELLSPERSIDE*FCCspacing;
            latVec[2] = ((double) j)*FCCspacing;
            for (n = 0; n < DIM; ++n)
            {
               coords[a*DIM + n] = latVec[n];
            }
            a++;
            for (n = 0; n < DIM; ++n)
            {
               coords[a*DIM + n] = latVec[n] + FCCshifts[2][n];
            }
            a++;
            /* pos-z face */
            latVec[0] = ((double) i)*FCCspacing;
            latVec[1] = ((double) j)*FCCspacing;
            latVec[2] = NCELLSPERSIDE*FCCspacing;
            for (n = 0; n < DIM; ++n)
            {
               coords[a*DIM + n] = latVec[n];
            }
            a++;
            for (n = 0; n < DIM; ++n)
            {
               coords[a*DIM + n] = latVec[n] + FCCshifts[1][n];
            }
            a++;
         }
      }
      if (!periodic)
      {
         /* add in the remaining three edges */
         latVec[0] = ((double) i)*FCCspacing;
         latVec[1] = NCELLSPERSIDE*FCCspacing;
         latVec[2] = NCELLSPERSIDE*FCCspacing;
         for (n = 0; n < DIM; ++n)
         {
            coords[a*DIM + n] = latVec[n];
         }
         a++;
         latVec[0] = NCELLSPERSIDE*FCCspacing;
         latVec[1] = ((double) i)*FCCspacing;
         latVec[2] = NCELLSPERSIDE*FCCspacing;
         for (n = 0; n < DIM; ++n)
         {
            coords[a*DIM + n] = latVec[n];
         }
         a++;
         latVec[0] = NCELLSPERSIDE*FCCspacing;
         latVec[1] = NCELLSPERSIDE*FCCspacing;
         latVec[2] = ((double) i)*FCCspacing;
         for (n = 0; n < DIM; ++n)
         {
            coords[a*DIM + n] = latVec[n];
         }
         a++;
      }
   }
   if (!periodic)
   {
      /* add in the remaining corner */
      for (n = 0; n < DIM; ++n)
      {
         coords[a*DIM + n] = NCELLSPERSIDE*FCCspacing;
      }
      a++;
   }

   return;
}
