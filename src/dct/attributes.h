#ifndef ATTRIBUTES_H_
#define ATTRIBUTES_H_

//#include "types.h"

void init_aux_data_store(AuxDataStore *store, size_t size_item, ulong size_array);
void clear_aux_data_store(AuxDataStore *store);
void *get_new_aux_data(AuxDataStore *store);
void realloc_store(AuxDataStore *store);

void *new_area_data(AuxDataStore *store, void *data);
void init_area_data(void *areaattr, bool init, ulong x, ulong y, ulong z);

void *load_area_data( AuxDataStore *store, void *initval);
void delete_area_data(void *areaattr);
void add_to_area_data(void *areaattr, void *data);
void merge_area_data(void *areaattr, void *childattr);
void merge_to_area_data(AuxDataStore *store, void **thisattr, void *areaattr, void *childattr); 
void clone_area_data(AuxDataStore *store, void **thisattr, void *areaattr);
double area_attribute(void *areaattr);

void *new_encl_rect_data(AuxDataStore *store, void *data);
void init_encl_rect_data(void *rectattr, bool init, ulong x, ulong y, ulong z);
void *load_encl_rect_data( AuxDataStore *store, void *initval);
void delete_encl_rect_data(void *rectattr);
void add_to_encl_rect_data(void *rectattr, void *data);
void merge_encl_rect_data(void *rectattr, void *childattr);
void merge_to_encl_rect_data(AuxDataStore *store, void **thisattr, void *rectattr, void *childattr);
void clone_encl_rect_data( AuxDataStore *store, void **thisattr, void *rectattr);
double encl_rect_area_attribute(void *rectattr);
double encl_rect_diag_attribute(void *rectattr);

void *new_inertia_data(AuxDataStore *store, void *data);
void init_inertia_data(void *inertiaattr, bool init, ulong x, ulong y, ulong z);
void *load_inertia_data( AuxDataStore *store, void *initval);
void delete_inertia_data(void *inertiaattr);
void add_to_inertia_data(void *inertiaattr, void *data);
void merge_inertia_data(void *inertiaattr, void *childattr);
void merge_to_inertia_data(AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr);
void clone_inertia_data(AuxDataStore *store, void **thisattr, void *inertiaattr);
double inertia_attribute(void *inertiaattr);
double inertia_div_a2_attribute(void *inertiaattr);
double mean_x_attribute(void *inertiaattr);
double mean_y_attribute(void *inertiaattr);
double mean_z_attribute(void *inertiaattr);

void init_attrib_array(void *attribute, bool border[6], ulong *dims, ulong attr_off[3], ulong lwb, ulong upb, ulong size_att);

void *(*new_aux_data)(AuxDataStore *, void *);
void (*init_aux_data)(void *,  bool , ulong , ulong,  ulong);
void (*delete_aux_data)(void *);
void (*add_to_aux_data)(void *, void *);
void (*merge_aux_data)(void *, void *);
void (*merge_to_aux_data)(AuxDataStore *, void **, void *, void *);
void (*clone_aux_data)(AuxDataStore *, void **, void *);
void (*create_mpi_aux_data)(void);


#endif 
