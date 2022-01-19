#include "types.h"
#include "attributes.h"
#include "checks.h"
#include "filter.h"
#include "eispack.h"

void init_aux_data_store(AuxDataStore *store, size_t size_item, ulong size_array){
  store->size_item  = size_item;
  store->size_alloc = size_array;
  store->item_curr  = 0;
  store->data = calloc(size_array, size_item);
  check_alloc(store->data, 100);
}

void clear_aux_data_store(AuxDataStore *store){ 
  free(store->data);
  free(store);
}

void *get_new_aux_data(AuxDataStore *store){
  
  char *p;
  if (store->item_curr >= store->size_alloc){
    debug(" Reallocation of attribute storage (may cause errors),size item %ld",  store->size_item);
    store->size_alloc = store->size_alloc > 1 ? store->size_alloc * 1.5 : store->size_alloc*100 ;
    store->data       = realloc(store->data, store->size_alloc * store->size_item);
    check_alloc(store->data, 101);
  }
  p = (char *)(store->data + (ulong) store->size_item * store->item_curr);
  store->item_curr++;
  return((void *) p);    
}

void realloc_store(AuxDataStore *store){  
  if (store->item_curr < store->size_alloc){
    store->data       = realloc(store->data, store->item_curr * store->size_item);
    check_alloc(store->data, 101);
  }
}


/****** Typedefs and functions for area attributes ******************************/

void *new_area_data(AuxDataStore *store, void *data){
   AreaData *areadata;
   areadata = store ? get_new_aux_data(store) :  calloc(1, sizeof(AreaData));
   check_alloc(areadata, 102);
   return(areadata);
} /* new_area_data */

void init_area_data(void *areaattr, bool init, ulong x, ulong y, ulong z, value gval){
  AreaData *areadata = areaattr;
  areadata->area = (ulong) init;
}

void *load_area_data( AuxDataStore *store, void *initval){
   AreaData *areadata;
   double *values = (double *) initval;
   ulong area = (ulong) (*values);
   areadata = store ? get_new_aux_data(store) :  calloc(1, sizeof(AreaData));
   check_alloc(areadata, 103);
   areadata->area = area;
   return(areadata);
} /* new_area_data */

void delete_area_data(void *areaattr){
  free(areaattr);
} /* delete_area_data */


void add_to_area_data(void *areaattr, void *data){
   AreaData *areadata = areaattr;
   areadata->area++;
} /* add_to_area_data */

void merge_area_data(void *areaattr, void *childattr){
   AreaData *areadata = areaattr;
   AreaData *childdata = childattr;

   areadata->area += childdata->area;
} /* merge_area_data */

void merge_to_area_data( AuxDataStore *store, void **thisattr, void *areaattr, void *childattr){
  AreaData *thisdata = *thisattr;
  AreaData *areadata = areaattr;
  AreaData *childdata = childattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(AreaData));
    check_alloc(thisdata, 104);
    *thisattr = thisdata;
  }
  thisdata->area = areadata->area + childdata->area;
} /* merge_to_area_data */

void clone_area_data( AuxDataStore *store, void **thisattr, void *areaattr){

  AreaData *thisdata = *thisattr;
  AreaData *areadata = areaattr;
  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(AreaData));
    check_alloc(thisdata, 105);
    *thisattr = thisdata;
  }

  thisdata->area = areadata->area;
} /* clone_area_data */

double area_attribute(void *areaattr){
   AreaData *areadata = areaattr;
   double area;

   area = areadata->area;
   return(area);
} /* area_attribute */

/****** Typedefs and functions for minimum enclosing rectangle attributes *******/

void *new_encl_rect_data( AuxDataStore *store, void *data){
   EnclRectData *rectdata;
   rectdata = store ? get_new_aux_data(store) :  calloc(1, sizeof(EnclRectData));
   check_alloc(rectdata, 106);
   return(rectdata);
} /* new_encl_rect_data */

void init_encl_rect_data(void *rectattr, bool init, ulong x, ulong y, ulong z, value gval){
  EnclRectData *rectdata = rectattr;
  // double *initcoord = (double *) data;
  if(init){
    rectdata->minX = rectdata->maxX = x;
    rectdata->minY = rectdata->maxY = y;
    rectdata->minZ = rectdata->maxZ = z;
  } else {
    rectdata->minX = rectdata->minY = rectdata->minZ = ULONG_MAX;
    rectdata->maxX = rectdata->maxY = rectdata->maxZ = 0;
  }

}

void *load_encl_rect_data( AuxDataStore *store, void *data){
   EnclRectData *rectdata;
    double *initcoord = (double *) data; 
   rectdata = store ? get_new_aux_data(store) :  calloc(1, sizeof(EnclRectData));
   check_alloc(rectdata, 107);   
   rectdata->minX = (ulong) initcoord[0];
   rectdata->maxX = (ulong) initcoord[1];
   rectdata->minY = (ulong) initcoord[2];
   rectdata->maxY = (ulong) initcoord[3];
   rectdata->minZ = (ulong) initcoord[4];
   rectdata->maxZ = (ulong) initcoord[5];
   return(rectdata);
} /* new_encl_rect_data */


void delete_encl_rect_data(void *rectattr)
{
  free(rectattr);
} /* delete_encl_rect_data */

void add_to_encl_rect_data(void *rectattr, void *data){
   EnclRectData *rectdata = rectattr;
   double *values = (double *) data; 

   rectdata->minX = MIN(rectdata->minX, (ulong) values[0]);
   rectdata->minY = MIN(rectdata->minY, (ulong) values[1]);
   rectdata->minZ = MIN(rectdata->minZ, (ulong) values[2]);
   rectdata->maxX = MAX(rectdata->maxX, (ulong) values[0]);
   rectdata->maxY = MAX(rectdata->maxY, (ulong) values[1]);
   rectdata->maxZ = MAX(rectdata->maxZ, (ulong) values[2]);
} /* add_to_encl_rect_data */

void merge_encl_rect_data(void *rectattr, void *childattr){
   EnclRectData *rectdata = rectattr;
   EnclRectData *childdata = childattr;

   rectdata->minX = MIN(rectdata->minX, childdata->minX);
   rectdata->minY = MIN(rectdata->minY, childdata->minY);
   rectdata->minZ = MIN(rectdata->minZ, childdata->minZ);
   rectdata->maxX = MAX(rectdata->maxX, childdata->maxX);
   rectdata->maxY = MAX(rectdata->maxY, childdata->maxY);
   rectdata->maxZ = MAX(rectdata->maxZ, childdata->maxZ);

} /* merge_encl_rect_data */

void merge_to_encl_rect_data( AuxDataStore *store, void **thisattr, void *rectattr, void *childattr){
  EnclRectData *thisdata = *thisattr;
  EnclRectData *rectdata = rectattr;
  EnclRectData *childdata = childattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(EnclRectData));
    check_alloc(thisdata, 108);
    *thisattr = thisdata;
  }
  thisdata->minX = MIN(rectdata->minX, childdata->minX);
  thisdata->minY = MIN(rectdata->minY, childdata->minY);
  thisdata->minZ = MIN(rectdata->minZ, childdata->minZ);
  thisdata->maxX = MAX(rectdata->maxX, childdata->maxX);
  thisdata->maxY = MAX(rectdata->maxY, childdata->maxY);
  thisdata->maxZ = MAX(rectdata->maxZ, childdata->maxZ);
} /* merge_to_encl_rect_data */

void clone_encl_rect_data( AuxDataStore *store, void **thisattr, void *rectattr){
  EnclRectData *thisdata = *thisattr;
  EnclRectData *rectdata = rectattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(EnclRectData));
    check_alloc(thisdata, 109);
    *thisattr = thisdata;
  }
  thisdata->minX = rectdata->minX;
  thisdata->minY = rectdata->minY;
  thisdata->minZ = rectdata->minZ;
  thisdata->maxX = rectdata->maxX;
  thisdata->maxY = rectdata->maxY;
  thisdata->maxZ = rectdata->maxZ;
} /* clone_encl_rect_data */

double encl_rect_area_attribute(void *rectattr){
  EnclRectData *rectdata = rectattr;
  ulong volume;
  volume =  (rectdata->maxX - rectdata->minX + 1)
    * (rectdata->maxY - rectdata->minY + 1)
    * (rectdata->maxZ - rectdata->minZ + 1);
  return (double) volume;
} /* encl_rect_area_attribute */

double encl_rect_diag_attribute(void *rectattr){/* Computes the square of the length of the diagonal */
  EnclRectData *rectdata = rectattr;
  ulong minx, miny, minz, maxx, maxy, maxz, l;

  minx = rectdata->minX;
  miny = rectdata->minY;
  minz = rectdata->minZ;
  maxx = rectdata->maxX;
  maxy = rectdata->maxY;
  maxz = rectdata->maxZ;
  l = (maxx-minx+1) * (maxx-minx+1)
    + (maxy-miny+1) * (maxy-miny+1)
    + (maxz-minz+1) * (maxz-minz+1);
  return(double) (l);
} /* encl_rect_diag_attribute */

/****** Typedefs and functions for moment of inertia attributes **************************/


void *new_inertia_data( AuxDataStore *store, void *data){
  InertiaData *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaData));
  check_alloc(inertiadata, 110);
  // init_inertia_data((void *) inertiadata, data);
  return(inertiadata);
} /* new_inertia_data */

void init_inertia_data(void *inertiaattr, bool init, ulong x, ulong y, ulong z, value gval){
  InertiaData *inertiadata = inertiaattr;
  //  double *init = (double *) data;
  inertiadata->area = (ulong) init;
  inertiadata->sumX = init ? (double) x : 0.;
  inertiadata->sumY = init ? (double) y : 0.;
  inertiadata->sumZ = init ? (double) z : 0.;
  inertiadata->sumR2 = init ? x*x + y*y+ z*z : 0.;
}

void *load_inertia_data( AuxDataStore *store, void *data){
  double *init = (double *) data;
  InertiaData *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaData));
  check_alloc(inertiadata, 111);
  
  inertiadata->area = init[0];
  inertiadata->sumX = init[1];
  inertiadata->sumY = init[2];
  inertiadata->sumZ = init[3];
  inertiadata->sumR2 = init[4];
  return(inertiadata);
} /* new_inertia_data */

void delete_inertia_data(void *inertiaattr)
{
  free(inertiaattr);
} /* delete_inertia_data */

void add_to_inertia_data(void *inertiaattr, void *data){
   InertiaData *inertiadata = inertiaattr;
   double *values = (double *) data;
   inertiadata->area++;
   inertiadata->sumX += values[0];
   inertiadata->sumY += values[1];
   inertiadata->sumZ += values[2];
   inertiadata->sumR2 += values[0]*values[0] + values[1]*values[1] + values[2]*values[2];;
} /* add_to_inertia_data */

void merge_inertia_data(void *inertiaattr, void *childattr)
{
   InertiaData *inertiadata = inertiaattr;
   InertiaData *childdata = childattr;

   inertiadata->area += childdata->area;
   inertiadata->sumX += childdata->sumX;
   inertiadata->sumY += childdata->sumY;
   inertiadata->sumZ += childdata->sumZ;
   inertiadata->sumR2 += childdata->sumR2;
} /* merge_inertia_data */

void merge_to_inertia_data( AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr){
   InertiaData *thisdata = *thisattr;
   InertiaData *inertiadata = inertiaattr;
   InertiaData *childdata = childattr;

   if (!thisdata) {
     thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaData));
     check_alloc(thisdata, 112);
     *thisattr = thisdata;
   }
   thisdata->area = inertiadata->area + childdata->area;
   thisdata->sumX = inertiadata->sumX + childdata->sumX;
   thisdata->sumY = inertiadata->sumY + childdata->sumY;
   thisdata->sumZ = inertiadata->sumZ + childdata->sumZ;
   thisdata->sumR2 = inertiadata->sumR2 + childdata->sumR2;
} /* merge_to_inertia_data */

void clone_inertia_data( AuxDataStore *store, void **thisattr, void *inertiaattr){
  InertiaData *thisdata = *thisattr;
  InertiaData *inertiadata = inertiaattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaData));
    check_alloc(thisdata, 113);
    *thisattr = thisdata;
  }
  thisdata->area = inertiadata->area;
  thisdata->sumX = inertiadata->sumX;
  thisdata->sumY = inertiadata->sumY;
  thisdata->sumZ = inertiadata->sumZ;
  thisdata->sumR2 = inertiadata->sumR2;
} /* clone_inertia_data */

double inertia_area_attribute(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;
   double area, inertia;
   area = inertiadata->area;
   return(area);
} /* inertia_area */

double inertia_attribute(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;
   double area, inertia;

   area = inertiadata->area;
   if(inertiadata->sumZ != 0){
     inertia = inertiadata->sumR2 -
       (inertiadata->sumX * inertiadata->sumX +
	inertiadata->sumY * inertiadata->sumY +
	inertiadata->sumZ * inertiadata->sumZ) / area
       + area / 4.0;  /* ??? */
   } else {
     inertia = inertiadata->sumR2 -
       (inertiadata->sumX * inertiadata->sumX +
	inertiadata->sumY * inertiadata->sumY ) / area
       + area / 6.0;  /* ??? */
   }
   return(inertia);

} /* inertia_attribute */

double inertia_div_a2_attribute(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;
   double inertia, area;

   area = (double)(inertiadata->area);

   if(inertiadata->sumZ != 0){
     inertia = inertiadata->sumR2 -
       (inertiadata->sumX * inertiadata->sumX +
	inertiadata->sumY * inertiadata->sumY +
	inertiadata->sumZ * inertiadata->sumZ) / area
       + area / 4.0;  /* ??? */
     return(inertia/pow(area,5./3.));
   } else {
     inertia = inertiadata->sumR2 -
       (inertiadata->sumX * inertiadata->sumX +
	inertiadata->sumY * inertiadata->sumY ) / area
       + area / 6.0;  /* ??? */
     return(inertia/pow(area,2));
   }
} /* inerti
} /* inertia_div_a2_attribute */

double mean_x_attribute(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;
   double area, sumx;

   area = inertiadata->area;
   sumx = inertiadata->sumX;
   return(sumx/area);
} /* mean_x_attribute */

double mean_y_attribute(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;
   double area, sumy;

   area = inertiadata->area;
   sumy = inertiadata->sumY;
   return(sumy/area);
} /* mean_y_attribute */

double mean_z_attribute(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;

   double area, sumz;

   area = inertiadata->area;
   sumz = inertiadata->sumZ;
   return(sumz/area);
} /* mean_z_attribute */



/* Inertia FULL analysis */



void *new_inertiafull_data( AuxDataStore *store, void *data){
  InertiaDataFull *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataFull));
  check_alloc(inertiadata, 114);
  /*
  inertiadata->area = 1;
  inertiadata->sumX = init[3]*init[0];
  inertiadata->sumY = init[3]*init[1];
  inertiadata->sumZ = init[3]*init[2];
  inertiadata->sumX2 = init[3]*(init[0]*init[0]);
  inertiadata->sumY2 = init[3]*(init[1]*init[1]);
  inertiadata->sumZ2 = init[3]*(init[2]*init[2]);
  inertiadata->sumXY = init[3]*init[0]*init[1];
  inertiadata->sumYZ = init[3]*init[1]*init[2];
  inertiadata->sumXZ = init[3]*init[0]*init[2];
  inertiadata->sumval = init[3];
  inertiadata->sumval2 = init[3]*init[3];
  */
  return(inertiadata);
} /* new_inertia_data */


void init_inertiafull_data(void *inertiaattr, bool init, ulong x, ulong y, ulong z, value gval){
  InertiaDataFull *inertiadata = inertiaattr;
  //  double *init = (double *) data;
  inertiadata->area = (ulong) init;
  inertiadata->sumX = init ? (double)  x : 0.;
  inertiadata->sumY = init ? (double)  y : 0.;
  inertiadata->sumZ = init ? (double)  z : 0.;
  inertiadata->sumX2 = init ? (double)  x * x: 0.;
  inertiadata->sumY2 = init ? (double)  y * y: 0.;
  inertiadata->sumZ2 = init ? (double)  z * z: 0.;
  inertiadata->sumXY = init ? (double)  x * y: 0.;
  inertiadata->sumYZ = init ? (double)  y * z: 0.;
  inertiadata->sumXZ = init ? (double)  z * x: 0.;
  inertiadata->sumval = init ? (double) 1: 0;
  inertiadata->sumval2 =  init ? (double)1 : 0;
}

void *load_inertiafull_data( AuxDataStore *store, void *data){
  InertiaDataFull *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataFull));
  check_alloc(inertiadata, 115);
    double *init = (double *) data;

  inertiadata->area = init[0];
  inertiadata->sumX = init[1];
  inertiadata->sumY = init[2];
  inertiadata->sumZ = init[3];
  //inertiadata->sumR2 = init[4];
  return(inertiadata);
} /* load_inertiafull_data */

void delete_inertiafull_data(void *inertiaattr)
{
  free(inertiaattr);
} /* delete_inertiafull_data */

void add_to_inertiafull_data(void *inertiaattr, void *data){
   InertiaDataFull *inertiadata = inertiaattr;
   double *values = (double *) data;
   inertiadata->area++;
   inertiadata->sumX += values[0];
   inertiadata->sumY += values[1];
   inertiadata->sumZ += values[2];
   inertiadata->sumX2 += values[0]*values[0];
   inertiadata->sumY2 += values[1]*values[1];
   inertiadata->sumZ2 += values[2]*values[2];
   inertiadata->sumXY += values[0]*values[1];
   inertiadata->sumYZ += values[1]*values[2];
   inertiadata->sumXZ += values[0]*values[2];
   inertiadata->sumval += 1;
   inertiadata->sumval2 += 1;
} /* add_to_inertiafull_data */

void merge_inertiafull_data(void *inertiaattr, void *childattr)
{
   InertiaDataFull *inertiadata = inertiaattr;
   InertiaDataFull *childdata = childattr;

   inertiadata->area += childdata->area;
   inertiadata->sumX += childdata->sumX;
   inertiadata->sumY += childdata->sumY;
   inertiadata->sumZ += childdata->sumZ;
   inertiadata->sumX2 += childdata->sumX2;
   inertiadata->sumY2 += childdata->sumY2;
   inertiadata->sumZ2 += childdata->sumZ2;
   inertiadata->sumXY += childdata->sumXY;
   inertiadata->sumYZ += childdata->sumYZ;
   inertiadata->sumXZ += childdata->sumXZ;
   inertiadata->sumval += childdata->sumval;
   inertiadata->sumval2 += childdata->sumval2;
} /* merge_inertia_data */


void merge_to_inertiafull_data( AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr){
   InertiaDataFull *thisdata = *thisattr;
   InertiaDataFull *inertiadata = inertiaattr;
   InertiaDataFull *childdata = childattr;

   if (!thisdata)  {
     thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataFull));
     check_alloc(thisdata, 116);
     *thisattr = thisdata;
   }
   thisdata->area = inertiadata->area + childdata->area;
   thisdata->sumX = inertiadata->sumX + childdata->sumX;
   thisdata->sumY = inertiadata->sumY + childdata->sumY;
   thisdata->sumZ = inertiadata->sumZ + childdata->sumZ;
   thisdata->sumX2 = inertiadata->sumX2 + childdata->sumX2;
   thisdata->sumY2 = inertiadata->sumY2 + childdata->sumY2;
   thisdata->sumZ2 = inertiadata->sumZ2 + childdata->sumZ2;
   thisdata->sumXY = inertiadata->sumXY + childdata->sumXY;
   thisdata->sumYZ = inertiadata->sumYZ + childdata->sumYZ;
   thisdata->sumXZ = inertiadata->sumXZ + childdata->sumXZ;
   thisdata->sumval = inertiadata->sumval + childdata->sumval;
   thisdata->sumval2 = inertiadata->sumval2 + childdata->sumval2;

   // thisdata->sumR2 = inertiadata->sumR2 + childdata->sumR2;
} /* merge_to_inertiafull_data */

void clone_inertiafull_data( AuxDataStore *store, void **thisattr, void *inertiaattr){
  InertiaDataFull *thisdata = *thisattr;
  InertiaDataFull *inertiadata = inertiaattr;

   if (!thisdata)  {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataFull));
    check_alloc(thisdata, 117);
    *thisattr = thisdata;
  }
  thisdata->area = inertiadata->area;
  thisdata->sumX = inertiadata->sumX;
  thisdata->sumY = inertiadata->sumY;
  thisdata->sumZ = inertiadata->sumZ;
  thisdata->sumX2 = inertiadata->sumX2;
  thisdata->sumY2 = inertiadata->sumY2;
  thisdata->sumZ2 = inertiadata->sumZ2;
  thisdata->sumXY = inertiadata->sumXY;
  thisdata->sumYZ = inertiadata->sumYZ;
  thisdata->sumXZ = inertiadata->sumXZ;
  thisdata->sumval = inertiadata->sumval;
  thisdata->sumval2 = inertiadata->sumval2;
  
} /* clone_inertiafull_data */

void *inertiafull_attribute_arr(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));

   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;
   return(inertia);
} /* inertia_attribute */

double inertiafull_area_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double inertia = inertiadata->area;

   return(inertia);
} /* inertia_attribute */



double  inertiafull_elon_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));
   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;
   double c_xyz[3]  = { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
   double corrVal    = inertia[10]/12;
   /* double corrValZ   = inertia[3] == 0? 0: inertia[10]/12; */

    double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[10] + corrVal,
			inertia[5] - inertia[2]*inertia[2]/inertia[10] + corrVal,
			inertia[6] - inertia[3]*inertia[3]/inertia[10] + corrVal,
			inertia[7] - inertia[1]*inertia[2]/inertia[10],
			inertia[8] - inertia[2]*inertia[3]/inertia[10],
			inertia[9] - inertia[1]*inertia[3]/inertia[10] };		
   
   double *eigval  = malloc(3  * sizeof(double));
   double *eigvec  = malloc(9  * sizeof(double));
  
   double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
   int ierr = rs (3, tens_mat, eigval, 1, eigvec);

   double elong = eigval[2]/eigval[1];
   /*  ncomp[i]    = (matrix[i*6]+matrix[i*6+1]+matrix[i*6+2])/pow(intens[i*2], 5.0/3.0);
   elong[i]    = eigval[i*3 + 2]/eigval[i*3 + 1];
   flat[i]     = eigval[i*3 + 1]/eigval[i*3 + 0];
   //  if ( eigval[i*3 + 0] < 0 || eigval[i*3 + 1] <0 || eigval[i*3 + 2] < 0) 
   // printf("%lf\n", eigval[i*3 + 0], eigval[i*3 + 1], eigval[i*3 + 2]);
   double ax_len[3] = {sqrt(20*eigval[i*3 + 2]/intens[i*2]), sqrt(20*eigval[i*3 + 1]/intens[i*2]), sqrt(20*eigval[i*3 + 0]/intens[i*2])};
   spars[i]     = intens[i*2]/(double) volume[i] * 3.14*ax_len[0]*ax_len[1]*ax_len[2]/(6*intens[i*2]);
   int pic = 2;
   if(eigval[i*3+2] == eigval[i*3+1])
    pic = eigval[i*3+1] == eigval[i*3] ? rand() % 3: rand() % 2 + 1;
   vec[3*i] = eigvec[i*9+pic*3];
   vec[3*i+1] = eigvec[i*9+pic*3+1];
   vec[3*i+2] = eigvec[i*9+pic*3+2]; */
   return(elong);
} /* mean_z_attribute */

double  inertiafull_flat_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));
   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;
   double c_xyz[3]  = { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
   double corrVal    = inertia[10]/12;

   if(inertia[3] == 0){
     return -1;
   }
   double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[10] + corrVal,
			inertia[5] - inertia[2]*inertia[2]/inertia[10] + corrVal,
			inertia[6] - inertia[3]*inertia[3]/inertia[10] + corrVal,
			inertia[7] - inertia[1]*inertia[2]/inertia[10],
			inertia[8] - inertia[2]*inertia[3]/inertia[10],
			inertia[9] - inertia[1]*inertia[3]/inertia[10] };		
   
   double *eigval  = malloc(3  * sizeof(double));
   double *eigvec  = malloc(9  * sizeof(double));
  
   double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
   int ierr = rs (3, tens_mat, eigval, 1, eigvec);
   // printf("%lf, %lf, %lf, %lf \n",inertia[3],eigval[0], eigval[1],eigval[2]);
   //printf("%lf, %lf, %lf, %lf \n",inertia[3], inertia[9], inertia[6], inertia[8]);
   double flat = eigval[1]/eigval[0];
   return(flat);
} /* mean_z_attribute */

double inertiafull_spar_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));
   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;
   double c_xyz[3]   = { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
   double corrVal    =  inertia[10]/12;

   double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[10] + corrVal,
			inertia[5] - inertia[2]*inertia[2]/inertia[10] + corrVal,
			inertia[6] - inertia[3]*inertia[3]/inertia[10] + corrVal,
			inertia[7] - inertia[1]*inertia[2]/inertia[10],
			inertia[8] - inertia[2]*inertia[3]/inertia[10],
			inertia[9] - inertia[1]*inertia[3]/inertia[10] };		
   
   double *eigval  = malloc(3  * sizeof(double));
   double *eigvec  = malloc(9  * sizeof(double));
  
   double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
   int ierr = rs (3, tens_mat, eigval, 1, eigvec);
   double ax_len[3] = {sqrt(20*eigval[2]/inertia[10]), sqrt(20*eigval[1]/inertia[10]), sqrt(20*eigval[0]/inertia[10])};
   double spars = inertia[10]/(double)inertia[0] * 3.14*ax_len[0]*ax_len[1]*ax_len[2]/(6*inertia[10]);
   
   //printf("%lf, %lf, %lf \n",ax_len[0], inertia[10], inertia[0]);

   return(spars);
} /* mean_z_attribute */

double  inertiafull_ncom_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));
   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;

   double c_xyz[3]  = { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
   double corrVal    = inertia[10]/12;
   double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[10] + corrVal,
			inertia[5] - inertia[2]*inertia[2]/inertia[10] + corrVal,
			inertia[6] - inertia[3]*inertia[3]/inertia[10] + corrVal,
			inertia[7] - inertia[1]*inertia[2]/inertia[10],
			inertia[8] - inertia[2]*inertia[3]/inertia[10],
			inertia[9] - inertia[1]*inertia[3]/inertia[10] };
   double ncomp    =  inertia[3] != 0? (matrix[0]+matrix[1]+matrix[2])/pow(inertia[10], 5.0/3.0):
     (matrix[0]+matrix[1])/pow(inertia[10], 2.0);

   // ncomp    =   (matrix[0]+matrix[1]+matrix[2])/pow(inertia[10], 5.0/3.0);
   return(ncomp);
} /* mean_z_attribute */



/* Inertia moments weighted by gvalues */


void *new_inertiaweight_data( AuxDataStore *store, void *data){
  InertiaDataFull *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataFull));
  check_alloc(inertiadata, 114);
  /*
  inertiadata->area = 1;
  inertiadata->sumX = init[3]*init[0];
  inertiadata->sumY = init[3]*init[1];
  inertiadata->sumZ = init[3]*init[2];
  inertiadata->sumX2 = init[3]*(init[0]*init[0]);
  inertiadata->sumY2 = init[3]*(init[1]*init[1]);
  inertiadata->sumZ2 = init[3]*(init[2]*init[2]);
  inertiadata->sumXY = init[3]*init[0]*init[1];
  inertiadata->sumYZ = init[3]*init[1]*init[2];
  inertiadata->sumXZ = init[3]*init[0]*init[2];
  inertiadata->sumval = init[3];
  inertiadata->sumval2 = init[3]*init[3];
  */
  return(inertiadata);
} /* new_inertia_data */


void init_inertiaweight_data(void *inertiaattr, bool init, ulong x, ulong y, ulong z, value gval){
  InertiaDataFull *inertiadata = inertiaattr;
  //  double *init = (double *) data;
  inertiadata->area = (ulong) init;
  inertiadata->sumX = init ? (double) gval * x : 0.;
  inertiadata->sumY = init ? (double) gval * y : 0.;
  inertiadata->sumZ = init ? (double) gval * z : 0.;
  inertiadata->sumX2 = init ? (double) gval * x * x: 0.;
  inertiadata->sumY2 = init ? (double) gval * y * y: 0.;
  inertiadata->sumZ2 = init ? (double) gval * z * z: 0.;
  inertiadata->sumXY = init ? (double) gval * x * y: 0.;
  inertiadata->sumYZ = init ? (double) gval * y * z: 0.;
  inertiadata->sumXZ = init ? (double) gval * z * x: 0.;
  inertiadata->sumval = init ? (double) gval: 0;
  inertiadata->sumval2 =  init ? (double)gval*gval : 0;
}

void *load_inertiaweight_data( AuxDataStore *store, void *data){
  InertiaDataFull *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataFull));
  check_alloc(inertiadata, 115);
    double *init = (double *) data;

  inertiadata->area = init[0];
  inertiadata->sumX = init[1];
  inertiadata->sumY = init[2];
  inertiadata->sumZ = init[3];
  //inertiadata->sumR2 = init[4];
  return(inertiadata);
} /* load_inertiaweight_data */

void delete_inertiaweight_data(void *inertiaattr)
{
  free(inertiaattr);
} /* delete_inertiafull_data */

void add_to_inertiaweight_data(void *inertiaattr, void *data){
   InertiaDataFull *inertiadata = inertiaattr;
   double *values = (double *) data;
   inertiadata->area++;
   inertiadata->sumX += values[3]*values[0];
   inertiadata->sumY += values[3]*values[1];
   inertiadata->sumZ += values[3]*values[2];
   inertiadata->sumX2 += values[3]*values[0]*values[0];
   inertiadata->sumY2 += values[3]*values[1]*values[1];
   inertiadata->sumZ2 += values[3]*values[2]*values[2];
   inertiadata->sumXY += values[3]*values[0]*values[1];
   inertiadata->sumYZ += values[3]*values[1]*values[2];
   inertiadata->sumXZ += values[3]*values[0]*values[2];
   inertiadata->sumval += values[3];
   inertiadata->sumval2 += values[3]*values[3];
} /* add_to_inertiaweight_data */

void merge_inertiaweight_data(void *inertiaattr, void *childattr)
{
   InertiaDataFull *inertiadata = inertiaattr;
   InertiaDataFull *childdata = childattr;

   inertiadata->area += childdata->area;
   inertiadata->sumX += childdata->sumX;
   inertiadata->sumY += childdata->sumY;
   inertiadata->sumZ += childdata->sumZ;
   inertiadata->sumX2 += childdata->sumX2;
   inertiadata->sumY2 += childdata->sumY2;
   inertiadata->sumZ2 += childdata->sumZ2;
   inertiadata->sumXY += childdata->sumXY;
   inertiadata->sumYZ += childdata->sumYZ;
   inertiadata->sumXZ += childdata->sumXZ;
   inertiadata->sumval += childdata->sumval;
   inertiadata->sumval2 += childdata->sumval2;
} /* merge_inertia_data */


void merge_to_inertiaweight_data( AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr){
   InertiaDataFull *thisdata = *thisattr;
   InertiaDataFull *inertiadata = inertiaattr;
   InertiaDataFull *childdata = childattr;

   if (!thisdata)  {
     thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataFull));
     check_alloc(thisdata, 116);
     *thisattr = thisdata;
   }
   thisdata->area = inertiadata->area + childdata->area;
   thisdata->sumX = inertiadata->sumX + childdata->sumX;
   thisdata->sumY = inertiadata->sumY + childdata->sumY;
   thisdata->sumZ = inertiadata->sumZ + childdata->sumZ;
   thisdata->sumX2 = inertiadata->sumX2 + childdata->sumX2;
   thisdata->sumY2 = inertiadata->sumY2 + childdata->sumY2;
   thisdata->sumZ2 = inertiadata->sumZ2 + childdata->sumZ2;
   thisdata->sumXY = inertiadata->sumXY + childdata->sumXY;
   thisdata->sumYZ = inertiadata->sumYZ + childdata->sumYZ;
   thisdata->sumXZ = inertiadata->sumXZ + childdata->sumXZ;
   thisdata->sumval = inertiadata->sumval + childdata->sumval;
   thisdata->sumval2 = inertiadata->sumval2 + childdata->sumval2;

   // thisdata->sumR2 = inertiadata->sumR2 + childdata->sumR2;
} /* merge_to_inertiaweight_data */

void clone_inertiaweight_data( AuxDataStore *store, void **thisattr, void *inertiaattr){
  InertiaDataFull *thisdata = *thisattr;
  InertiaDataFull *inertiadata = inertiaattr;

   if (!thisdata)  {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataFull));
    check_alloc(thisdata, 117);
    *thisattr = thisdata;
  }
  thisdata->area = inertiadata->area;
  thisdata->sumX = inertiadata->sumX;
  thisdata->sumY = inertiadata->sumY;
  thisdata->sumZ = inertiadata->sumZ;
  thisdata->sumX2 = inertiadata->sumX2;
  thisdata->sumY2 = inertiadata->sumY2;
  thisdata->sumZ2 = inertiadata->sumZ2;
  thisdata->sumXY = inertiadata->sumXY;
  thisdata->sumYZ = inertiadata->sumYZ;
  thisdata->sumXZ = inertiadata->sumXZ;
  thisdata->sumval = inertiadata->sumval;
  thisdata->sumval2 = inertiadata->sumval2;
  
} /* clone_inertiaweight_data */

void *inertiaweight_attribute_arr(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));

   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;
   return(inertia);
} /* inertia_attribute */

double inertiaweight_area_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double inertia = inertiadata->area;

   return(inertia);
} /* inertia_attribute */



double  inertiaweight_elon_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));
   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;
   double c_xyz[3]  = { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
   double corrVal    = inertia[10]/12;
   double corrValZ   = inertia[3] == 0? 0: inertia[10]/12;

   double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[10] + corrVal,
			inertia[5] - inertia[2]*inertia[2]/inertia[10] + corrVal,
			inertia[6] - inertia[3]*inertia[3]/inertia[10] + corrValZ,
			inertia[7] - inertia[1]*inertia[2]/inertia[10],
			inertia[8] - inertia[2]*inertia[3]/inertia[10],
			inertia[9] - inertia[1]*inertia[3]/inertia[10] };		
   
   double *eigval  = malloc(3  * sizeof(double));
   double *eigvec  = malloc(9  * sizeof(double));
  
   double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
   int ierr = rs (3, tens_mat, eigval, 1, eigvec);

   double elong = eigval[2]/eigval[1];
   /*  ncomp[i]    = (matrix[i*6]+matrix[i*6+1]+matrix[i*6+2])/pow(intens[i*2], 5.0/3.0);
   elong[i]    = eigval[i*3 + 2]/eigval[i*3 + 1];
   flat[i]     = eigval[i*3 + 1]/eigval[i*3 + 0];
   //  if ( eigval[i*3 + 0] < 0 || eigval[i*3 + 1] <0 || eigval[i*3 + 2] < 0) 
   // printf("%lf\n", eigval[i*3 + 0], eigval[i*3 + 1], eigval[i*3 + 2]);
   double ax_len[3] = {sqrt(20*eigval[i*3 + 2]/intens[i*2]), sqrt(20*eigval[i*3 + 1]/intens[i*2]), sqrt(20*eigval[i*3 + 0]/intens[i*2])};
   spars[i]     = intens[i*2]/(double) volume[i] * 3.14*ax_len[0]*ax_len[1]*ax_len[2]/(6*intens[i*2]);
   int pic = 2;
   if(eigval[i*3+2] == eigval[i*3+1])
    pic = eigval[i*3+1] == eigval[i*3] ? rand() % 3: rand() % 2 + 1;
   vec[3*i] = eigvec[i*9+pic*3];
   vec[3*i+1] = eigvec[i*9+pic*3+1];
   vec[3*i+2] = eigvec[i*9+pic*3+2]; */
   return(elong);
} /* mean_z_attribute */

double  inertiaweight_flat_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));
   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;
   double c_xyz[3]  = { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
   double corrVal    = inertia[10]/12;
   double corrValZ   = inertia[3] == 0? 0: inertia[10]/12;

   if(inertia[3] == 0){
     return -1;
   }
   double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[10] + corrVal,
			inertia[5] - inertia[2]*inertia[2]/inertia[10] + corrVal,
			inertia[6] - inertia[3]*inertia[3]/inertia[10] + corrValZ,
			inertia[7] - inertia[1]*inertia[2]/inertia[10],
			inertia[8] - inertia[2]*inertia[3]/inertia[10],
			inertia[9] - inertia[1]*inertia[3]/inertia[10] };		
   
   double *eigval  = malloc(3  * sizeof(double));
   double *eigvec  = malloc(9  * sizeof(double));
  
   double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
   int ierr = rs (3, tens_mat, eigval, 1, eigvec);
   // printf("%lf, %lf, %lf, %lf \n",inertia[3],eigval[0], eigval[1],eigval[2]);
   //printf("%lf, %lf, %lf, %lf \n",inertia[3], inertia[9], inertia[6], inertia[8]);
   double flat = eigval[1]/eigval[0];
   return(flat);
} /* mean_z_attribute */

double inertiaweight_spar_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));
   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;
   double c_xyz[3]  = { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
   double corrVal    = inertia[10]/12;
   double corrValZ   =  inertia[10]/12;

   double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[10] + corrVal,
			inertia[5] - inertia[2]*inertia[2]/inertia[10] + corrVal,
			inertia[6] - inertia[3]*inertia[3]/inertia[10] + corrValZ,
			inertia[7] - inertia[1]*inertia[2]/inertia[10],
			inertia[8] - inertia[2]*inertia[3]/inertia[10],
			inertia[9] - inertia[1]*inertia[3]/inertia[10] };		
   
   double *eigval  = malloc(3  * sizeof(double));
   double *eigvec  = malloc(9  * sizeof(double));
  
   double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
   int ierr = rs (3, tens_mat, eigval, 1, eigvec);
   double ax_len[3] = {sqrt(20*eigval[2]/inertia[10]), sqrt(20*eigval[1]/inertia[10]), sqrt(20*eigval[0]/inertia[10])};
   double spars = inertia[10]/(double)inertia[0] * 3.14*ax_len[0]*ax_len[1]*ax_len[2]/(6*inertia[10]);
   
   printf("%lf, %lf, %lf \n",ax_len[0], ax_len[1], ax_len[2]);

   return(spars);
} /* mean_z_attribute */

double  inertiaweight_ncom_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double *inertia = calloc(12, sizeof(double));
   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->sumX;
   inertia[2] = inertiadata->sumY;
   inertia[3] = inertiadata->sumZ;
   inertia[4] = inertiadata->sumX2;
   inertia[5] = inertiadata->sumY2;
   inertia[6] = inertiadata->sumZ2;
   inertia[7] = inertiadata->sumXY;
   inertia[8] = inertiadata->sumYZ;
   inertia[9] = inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;

   double c_xyz[3]  = { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
   double corrVal    = inertia[10]/12;
   double corrValZ   = inertia[3] == 0? 0: inertia[10]/12;

   double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[10] + corrVal,
			inertia[5] - inertia[2]*inertia[2]/inertia[10] + corrVal,
			inertia[6] - inertia[3]*inertia[3]/inertia[10] + corrValZ,
			inertia[7] - inertia[1]*inertia[2]/inertia[10],
			inertia[8] - inertia[2]*inertia[3]/inertia[10],
			inertia[9] - inertia[1]*inertia[3]/inertia[10] };
   double ncomp    = (matrix[0]+matrix[1]+matrix[2])/pow(inertia[10], 5.0/3.0);

   return(ncomp);
} /* mean_z_attribute */


/*    */
void init_attrib_array(void *attribute, value *gvals, bool border[6], ulong *dims, ulong attr_off[3], ulong lwb, ulong upb, ulong size_att){
  ulong width  = dims[0];
  ulong height = dims[1];
  ulong depth  = dims[2];
  ulong size2D = width * height;

  for (ulong i =lwb; i < upb; i++){
    ulong x, y, z;
    x = (i % size2D) % width;
    y = (i % size2D) / width;
    z = i / size2D;
    bool init = (((x == 0 && border[0]) || (x == width-1  && border[1])  ||
		  (y == 0 && border[2]) || (y == height-1 && border[3])  ||
		  (z == 0 && border[4]) || (z == depth-1  && border[5])) || (i >= size2D*depth));
    init_aux_data(attribute + i*size_att, !init, x + attr_off[0], y + attr_off[1], z+ attr_off[2], 1.);
  }
}


AttribStruct AttribsArray[NUMATTR] =
{
 {"Area", sizeof(AreaData), new_area_data, init_area_data, delete_area_data, add_to_area_data, merge_area_data, merge_to_area_data, clone_area_data, create_mpi_area_type, area_attribute, area_attribute}, //1
 {"Area of min. enclosing rectangle", sizeof(EnclRectData),  new_encl_rect_data, init_encl_rect_data, delete_encl_rect_data, add_to_encl_rect_data, merge_encl_rect_data, merge_to_encl_rect_data, clone_encl_rect_data, create_mpi_rect_type, encl_rect_area_attribute, encl_rect_area_attribute}, //2 TO CHANGE AREA!!!!
 {"Square of diagonal of min. enclosing rectangle", sizeof(EnclRectData), new_encl_rect_data, init_encl_rect_data, delete_encl_rect_data, add_to_encl_rect_data, merge_encl_rect_data, merge_to_encl_rect_data, clone_encl_rect_data, create_mpi_rect_type, encl_rect_diag_attribute, encl_rect_diag_attribute},//3
 {"Moment of Inertia", sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data, create_mpi_inertia_type, inertia_area_attribute, inertia_attribute}, //4
 {"(Moment of Inertia) / (area)^2", sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data,  create_mpi_inertia_type, inertia_area_attribute, inertia_div_a2_attribute},//5
 {"Mean X position", sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data,  create_mpi_inertia_type, inertia_area_attribute, mean_x_attribute},//6
 {"Mean Y position", sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data,  create_mpi_inertia_type, inertia_area_attribute, mean_y_attribute},//7
 {"Mean Z position", sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data,  create_mpi_inertia_type, inertia_area_attribute, mean_z_attribute},//8
 {"Full inertia: area", sizeof(InertiaDataFull), new_inertiafull_data, init_inertiafull_data, delete_inertiafull_data, add_to_inertiafull_data, merge_inertiafull_data, merge_to_inertiafull_data, clone_inertiafull_data,  create_mpi_inertiafull_type, inertiafull_area_attribute, inertiafull_area_attribute},//9
 {"Full inertia: elongation", sizeof(InertiaDataFull), new_inertiafull_data, init_inertiafull_data, delete_inertiafull_data, add_to_inertiafull_data, merge_inertiafull_data, merge_to_inertiafull_data, clone_inertiafull_data,  create_mpi_inertiafull_type, inertiafull_area_attribute, inertiafull_elon_attribute},//10
 {"Full inertia: flatness", sizeof(InertiaDataFull), new_inertiafull_data, init_inertiafull_data, delete_inertiafull_data, add_to_inertiafull_data, merge_inertiafull_data, merge_to_inertiafull_data, clone_inertiafull_data,  create_mpi_inertiafull_type, inertiafull_area_attribute,inertiafull_flat_attribute},//11
 {"Full inertia: sparseness", sizeof(InertiaDataFull), new_inertiafull_data, init_inertiafull_data, delete_inertiafull_data, add_to_inertiafull_data, merge_inertiafull_data, merge_to_inertiafull_data, clone_inertiafull_data,  create_mpi_inertiafull_type, inertiafull_area_attribute,inertiafull_spar_attribute}, //12
 {"Full inertia: non-compactness", sizeof(InertiaDataFull), new_inertiafull_data, init_inertiafull_data, delete_inertiafull_data, add_to_inertiafull_data, merge_inertiafull_data, merge_to_inertiafull_data, clone_inertiafull_data,  create_mpi_inertiafull_type, inertiafull_area_attribute, inertiafull_ncom_attribute} //13
};

