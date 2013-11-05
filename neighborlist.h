#ifndef _NEIGHBORLIST_H_
#define _NEIGHBORLIST_H_

/* 
    these are listed in the order they should be called with
    a choice between set_neigh_object and calculate_neigh
*/
int nbl_initialize(void* kimmdl); 
int nbl_set_cell(int S1, double* Cell, int S2, int* PBC);
int nbl_set_ghosts(int N1, int* Ghosts, int ishalf);
int nbl_build_neighborlist(void *kimmdl);
int nbl_set_neigh_object(void* kimmdl, int sz1, int* NNeighbors, 
                                   int sz2, int* HalfNNeighbors, 
                                   int sz3, int* neighborList, 
                                   int sz4, double* RijList); 
int nbl_cleanup(void* kimmdl); 

#endif
