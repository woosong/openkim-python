#ifndef _NEIGHBORLIST_H_
#define _NEIGHBORLIST_H_
int set_neigh_object(void* kimmdl, int sz1, int* NNeighbors, int sz2, int* neighborList, int sz3, double* RijList);  
int set_kim_periodic_full_neigh(void* kimmdl); 
int set_kim_periodic_half_neigh(void* kimmdl); 
int set_kim_cluster_full_neigh(void* kimmdl); 
int set_kim_cluster_half_neigh(void* kimmdl); 

#endif
