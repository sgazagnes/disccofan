#ifndef ATTRIBUTES_H
#define ATTRIBUTES_H

//#include "types.h"

void init_aux_data_store(AuxDataStore *store, size_t size_item, ulong size_array);
void clear_aux_data_store(AuxDataStore *store);
void *get_new_aux_data(AuxDataStore *store);
void realloc_store(AuxDataStore *store);

void *new_area_data(AuxDataStore *store);
void init_area_data(void *areaattr, bool init, ulong x, ulong y, ulong z, value gval);
void delete_area_data(void *areaattr);
void add_to_area_data(void *areaattr, void *data);
void merge_area_data(void *areaattr, void *childattr);
void merge_to_area_data(AuxDataStore *store, void **thisattr, void *areaattr, void *childattr); 
void clone_area_data(AuxDataStore *store, void **thisattr, void *areaattr);
double area_attribute(void *areaattr);

void *new_extent_data( AuxDataStore *store);
void init_extent_data(void *rectattr, bool init, ulong x, ulong y, ulong z, value gval);
void delete_extent_data(void *rectattr);
void add_to_extent_data(void *rectattr, void *data);
void merge_extent_data(void *rectattr, void *childattr);
void merge_to_extent_data( AuxDataStore *store, void **thisattr, void *rectattr, void *childattr);
void clone_extent_data( AuxDataStore *store, void **thisattr, void *rectattr);
double extent_rectarea_attribute(void *rectattr);
double extent_rectdiag_attribute(void *rectattr);
double extent_x_attribute(void *rectattr);
double extent_y_attribute(void *rectattr);
double extent_z_attribute(void *rectattr);

void *new_mean_data( AuxDataStore *store);
void init_mean_data(void *meanattr, bool init, ulong x, ulong y, ulong z, value gval);
void delete_mean_data(void *meanattr);
void add_to_mean_data(void *meanattr, void *data);
void merge_mean_data(void *meanattr, void *childattr);
void merge_to_mean_data( AuxDataStore *store, void **thisattr, void *meanattr, void *childattr);
void clone_mean_data( AuxDataStore *store, void **thisattr, void *meanattr);
double mean_area_attribute(void *meanattr);
double mean_x_attribute(void *meanattr);
double mean_y_attribute(void *meanattr);
double mean_z_attribute(void *meanattr);

void *new_wmean_data( AuxDataStore *store);
void init_wmean_data(void *wmeanattr, bool init, ulong x, ulong y, ulong z, value gval);
void delete_wmean_data(void *wmeanattr);
void add_to_wmean_data(void *wmeanattr, void *data);
void merge_wmean_data(void *wmeanattr, void *childattr);
void merge_to_wmean_data( AuxDataStore *store, void **thisattr, void *wmeanattr, void *childattr);
void clone_wmean_data( AuxDataStore *store, void **thisattr, void *wmeanattr);
double wmean_totalflux_attribute(void *wmeanattr);
double wmean_x_attribute(void *wmeanattr);
double wmean_y_attribute(void *wmeanattr);
double wmean_z_attribute(void *wmeanattr);


void *new_inertia_data(AuxDataStore *store);
void init_inertia_data(void *inertiaattr, bool init, ulong x, ulong y, ulong z, value gval);
void delete_inertia_data(void *inertiaattr);
void add_to_inertia_data(void *inertiaattr, void *data);
void merge_inertia_data(void *inertiaattr, void *childattr);
void merge_to_inertia_data(AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr);
void clone_inertia_data(AuxDataStore *store, void **thisattr, void *inertiaattr);
void *inertia_arr(void *inertiaattr);
void *inertia_attribute_arr(void *inertiaattr);
double inertia_area_attribute(void *inertiaattr);
double  inertia_trace_attribute(void *inertiaattr);
double  inertia_elon_attribute(void *inertiaattr);
double  inertia_flat_attribute(void *inertiaattr);
double  inertia_spar_attribute(void *inertiaattr);
double  inertia_ncom_attribute(void *inertiaattr);


void init_attrib_array(void *attribute, value *gvals, bool border[6], ulong *dims, ulong attr_off[3], ulong lwb, ulong upb, ulong size_att);


/*void *(*new_aux_data)(AuxDataStore *, void *);
void (*init_aux_data)(void *, bool , ulong , ulong,  ulong, value);
void (*delete_aux_data)(void *);
void (*add_to_aux_data)(void *, void *);
void (*merge_aux_data)(void *, void *);
void (*merge_to_aux_data)(AuxDataStore *, void **, void *, void *);
void (*clone_aux_data)(AuxDataStore *, void **, void *);
void (*create_mpi_aux_data)(void);*/
#endif 
