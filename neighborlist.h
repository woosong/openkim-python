#ifndef _NEIGHBORLIST_H_
#define _NEIGHBORLIST_H_

/* 
    these are listed in the order they should be called with
    a choice between set_neigh_object and calculate_neigh
*/
int initialize(void* kimmdl); 
int set_neigh_object(void* kimmdl, int sz1, int* NNeighbors, 
                                   int sz2, int* HalfNNeighbors, 
                                   int sz3, int* neighborList, 
                                   int sz4, double* RijList); 
int build_neighborlist_allall(void *kimmdl);
int build_neighborlist_cell(void *kimmdl);
int free_neigh_object(void* kimmdl); 
#endif
