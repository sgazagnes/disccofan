#ifndef ATTRIBUTES_H_
#define ATTRIBUTES_H_

//#include "types.h"

void init_aux_data_store(AuxDataStore *store, size_t size_item, ulong size_array);
void clear_aux_data_store(AuxDataStore *store);
void *get_new_aux_data(AuxDataStore *store);
void realloc_store(AuxDataStore *store);

void *new_area_data(AuxDataStore *store, double *init);
void *load_area_data( AuxDataStore *store, double *init);
void delete_area_data(void *areaattr);
void add_to_area_data(void *areaattr, double *init);
void merge_area_data(void *areaattr, void *childattr);
void merge_to_area_data(AuxDataStore *store, void **thisattr, void *areaattr, void *childattr); 
void clone_area_data(AuxDataStore *store, void **thisattr, void *areaattr);
double area_attribute(void *areaattr);

void *new_encl_rect_data(AuxDataStore *store, double *init);
void *load_encl_rect_data( AuxDataStore *store, double *init);
void delete_encl_rect_data(void *rectattr);
void add_to_encl_rect_data(void *rectattr, double *init);
void merge_encl_rect_data(void *rectattr, void *childattr);
void merge_to_encl_rect_data(AuxDataStore *store, void **thisattr, void *rectattr, void *childattr);
void clone_encl_rect_data( AuxDataStore *store, void **thisattr, void *rectattr);
double encl_rect_area_attribute(void *rectattr);
double encl_rect_diag_attribute(void *rectattr);

void *new_inertia_data(AuxDataStore *store, double *init);
void *load_inertia_data( AuxDataStore *store, double *init);
void delete_inertia_data(void *inertiaattr);
void add_to_inertia_data(void *inertiaattr, double *init);
void merge_inertia_data(void *inertiaattr, void *childattr);
void merge_to_inertia_data(AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr);
void clone_inertia_data(AuxDataStore *store, void **thisattr, void *inertiaattr);
double inertia_attribute(void *inertiaattr);
double inertia_div_a2_attribute(void *inertiaattr);
double mean_x_attribute(void *inertiaattr);
double mean_y_attribute(void *inertiaattr);
double mean_z_attribute(void *inertiaattr);

void *new_inertiaeor_data(AuxDataStore *store, double *init);
void *load_inertiaeor_data( AuxDataStore *store, double *init);
void delete_inertiaeor_data(void *inertiaattr);
void add_to_inertiaeor_data(void *inertiaattr, double *init);
void merge_inertiaeor_data(void *inertiaattr, void *childattr);
void merge_to_inertiaeor_data(AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr);
void clone_inertiaeor_data(AuxDataStore *store, void **thisattr, void *inertiaattr);
void *inertiaeor_attribute_arr(void *inertiaattr);
void correct_inertiaeor_data(void *inertiaattr, void *childattr);
double inertiaeor_attribute(void *inertiaattr);

/*void *new_eor_data(AuxDataStore *store, double *init);
void *load_eor_data( AuxDataStore *store, double *init);
void delete_eor_data(void *inertiaattr);
void add_to_eor_data(void *inertiaattr, double *init);
void merge_eor_data(void *inertiaattr, void *childattr);
void merge_to_eor_data(AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr);
void clone_eor_data(AuxDataStore *store, void **thisattr, void *inertiaattr);
void *eor_attribute(void *inertiaattr);*/


void *(*new_aux_data)(AuxDataStore *, double*);
void (*delete_aux_data)(void *);
void (*add_to_aux_data)(void *, double*);
void (*merge_aux_data)(void *, void *);
void (*merge_to_aux_data)(AuxDataStore *, void **, void *, void *);
void (*clone_aux_data)(AuxDataStore *, void **, void *);
void (*create_mpi_aux_data)(void);
void *(*read_aux_file_binary)(AuxDataStore *, FILE *);
void (*write_aux_file_binary)(FILE *, void *);

#endif /* DECISION3D_H_ */
