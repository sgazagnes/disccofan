#include "types.h"
#include "attributes.h"
#include "checks.h"
#include "tree_filt.h"
#include "writefile.h"

ulong DIMS[3];
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
  p = store->data + (ulong) store->size_item * store->item_curr;
  store->item_curr++;
  return((void *) p);    
}

void realloc_store(AuxDataStore *store){  
  char *p;
  if (store->item_curr < store->size_alloc){
    store->data       = realloc(store->data, store->item_curr * store->size_item);
    check_alloc(store->data, 101);
  }
}


/****** Typedefs and functions for area attributes ******************************/

void *new_area_data( AuxDataStore *store, double *init){
   AreaData *areadata;
   areadata = store ? get_new_aux_data(store) :  calloc(1, sizeof(AreaData));
   check_alloc(areadata, 102);
   areadata->area = 1;
   return(areadata);
} /* new_area_data */

void *load_area_data( AuxDataStore *store, double *init){
   AreaData *areadata;
   ulong area = *(ulong *)(init);
   areadata = store ? get_new_aux_data(store) :  calloc(1, sizeof(AreaData));
   check_alloc(areadata, 103);
   areadata->area = area;
   return(areadata);
} /* new_area_data */

void delete_area_data(void *areaattr){
  free(areaattr);
} /* delete_area_data */

void ResetAreaData(void *areaattr){
   AreaData *areadata = areaattr;
   areadata->area = 0;
} /* add_to_area_data */


void add_to_area_data(void *areaattr, double *init){
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

void *new_encl_rect_data( AuxDataStore *store, double *init){
   EnclRectData *rectdata;
   rectdata = store ? get_new_aux_data(store) :  calloc(1, sizeof(EnclRectData));
   check_alloc(rectdata, 106);
   rectdata->minX = rectdata->maxX = (ulong) init[0];
   rectdata->minY = rectdata->maxY = (ulong) init[1];
   rectdata->minZ = rectdata->maxZ = (ulong) init[2];
   return(rectdata);
} /* new_encl_rect_data */

void *load_encl_rect_data( AuxDataStore *store, double *init){
   EnclRectData *rectdata;
   rectdata = store ? get_new_aux_data(store) :  calloc(1, sizeof(EnclRectData));
   check_alloc(rectdata, 107);   
   rectdata->minX = (ulong) init[0];
   rectdata->maxX = (ulong) init[1];
   rectdata->minY = (ulong) init[2];
   rectdata->maxY = (ulong) init[3];
   rectdata->minZ = (ulong) init[4];
   rectdata->maxZ = (ulong) init[5];
   return(rectdata);
} /* new_encl_rect_data */


void delete_encl_rect_data(void *rectattr)
{
  free(rectattr);
} /* delete_encl_rect_data */

void add_to_encl_rect_data(void *rectattr, double *init){
   EnclRectData *rectdata = rectattr;

   rectdata->minX = MIN(rectdata->minX, (ulong) init[0]);
   rectdata->minY = MIN(rectdata->minY, (ulong) init[1]);
   rectdata->minZ = MIN(rectdata->minZ, (ulong) init[2]);
   rectdata->maxX = MAX(rectdata->maxX, (ulong) init[0]);
   rectdata->maxY = MAX(rectdata->maxY, (ulong) init[1]);
   rectdata->maxZ = MAX(rectdata->maxZ, (ulong) init[2]);
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
  double volume;

  volume = (rectdata->maxX - rectdata->minX + 1)
    * (rectdata->maxY - rectdata->minY + 1)
    * (rectdata->maxZ - rectdata->minZ + 1);
  return(volume);
} /* encl_rect_area_attribute */

double encl_rect_diag_attribute(void *rectattr){/* Computes the square of the length of the diagonal */
  EnclRectData *rectdata = rectattr;
  double minx, miny, minz, maxx, maxy, maxz, l;

  minx = rectdata->minX;
  miny = rectdata->minY;
  minz = rectdata->minZ;
  maxx = rectdata->maxX;
  maxy = rectdata->maxY;
  maxz = rectdata->maxZ;
  l = (maxx-minx+1) * (maxx-minx+1)
    + (maxy-miny+1) * (maxy-miny+1)
    + (maxz-minz+1) * (maxz-minz+1);
  return(l);
} /* encl_rect_diag_attribute */

/****** Typedefs and functions for moment of inertia attributes **************************/

void *new_inertia_data( AuxDataStore *store, double *init){
  InertiaData *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaData));
  check_alloc(inertiadata, 110);
  inertiadata->area = 1;
  inertiadata->sumX =  init[0];
  inertiadata->sumY =  init[1];
  inertiadata->sumZ =  init[2];
  inertiadata->sumR2 =  init[0]*init[0] + init[1]*init[1] + init[2]*init[2];
  return(inertiadata);
} /* new_inertia_data */

void *load_inertia_data( AuxDataStore *store, double *init){
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

void add_to_inertia_data(void *inertiaattr, double *init){
   InertiaData *inertiadata = inertiaattr;

   inertiadata->area ++;
   inertiadata->sumX += init[0];
   inertiadata->sumY += init[1];
   inertiadata->sumZ += init[2];
   inertiadata->sumR2 += init[0]*init[0] + init[1]*init[1] + init[2]*init[2];
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

double inertia_attribute(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;
   double area, inertia;

   area = inertiadata->area;
   inertia = inertiadata->sumR2 -
             (inertiadata->sumX * inertiadata->sumX +
              inertiadata->sumY * inertiadata->sumY +
              inertiadata->sumZ * inertiadata->sumZ) / area
             + area / 6.0;  /* ??? */
   return(inertia);
} /* inertia_attribute */

double inertia_div_a2_attribute(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;
   double inertia, area;

   area = (double)(inertiadata->area);
   inertia = inertiadata->sumR2 -
             (inertiadata->sumX * inertiadata->sumX +
              inertiadata->sumY * inertiadata->sumY +
              inertiadata->sumZ * inertiadata->sumZ) / area
             + area / 6.0;  /* ??? */
   return(inertia/pow(area,5.0/3.0));
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

/****** Typedefs and functions for moment of inertia attributes for EoR Application *****************/


void *new_inertiaeor_data( AuxDataStore *store, double *init){
  InertiaDataEoR *inertiadata;
  //debug("%lf, %lf, %lf", init[0], init[1], init[2]);

  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataEoR));
  check_alloc(inertiadata, 114);
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
  /*  inertiadata->sumXd = init[4]*init[0];
  inertiadata->sumYd = init[4]*init[1];
  inertiadata->sumZd = init[4]*init[2];
  inertiadata->sumX2d = init[4]*(init[0]*init[0]);
  inertiadata->sumY2d = init[4]*(init[1]*init[1]);
  inertiadata->sumZ2d = init[4]*(init[2]*init[2]);
  inertiadata->sumXYd = init[4]*init[0]*init[1];
  inertiadata->sumYZd = init[4]*init[1]*init[2];
  inertiadata->sumXZd = init[4]*init[0]*init[2];
  inertiadata->sumvald = init[4];
  inertiadata->sumval2d = init[4]*init[4];*/
  // inertiadata->sumR2 = x*x + y*y + z*z;
  return(inertiadata);
} /* new_inertia_data */

void *load_inertiaeor_data( AuxDataStore *store, double *init){
  InertiaDataEoR *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataEoR));
  check_alloc(inertiadata, 115);
  
  inertiadata->area = init[0];
  inertiadata->sumX = init[1];
  inertiadata->sumY = init[2];
  inertiadata->sumZ = init[3];
  //inertiadata->sumR2 = init[4];
  return(inertiadata);
} /* new_inertia_data */

void delete_inertiaeor_data(void *inertiaattr)
{
  free(inertiaattr);
} /* delete_inertia_data */

void add_to_inertiaeor_data(void *inertiaattr, double *init){
   InertiaDataEoR *inertiadata = inertiaattr;
   inertiadata->area ++;
   inertiadata->sumX += init[3]*init[0];
   inertiadata->sumY += init[3]*init[1];
   inertiadata->sumZ += init[3]*init[2];
   inertiadata->sumX2 += init[3]*init[0]*init[0];
   inertiadata->sumY2 += init[3]*init[1]*init[1];
   inertiadata->sumZ2 += init[3]*init[2]*init[2];
   inertiadata->sumXY += init[3]*init[0]*init[1];
   inertiadata->sumYZ += init[3]*init[1]*init[2];
   inertiadata->sumXZ += init[3]*init[0]*init[2];
   inertiadata->sumval += init[3];
   inertiadata->sumval2 += init[3]*init[3];
   /* inertiadata->sumXd += init[4]*init[0];
   inertiadata->sumYd += init[4]*init[1];
   inertiadata->sumZd += init[4]*init[2];
   inertiadata->sumX2d += init[4]*init[0]*init[0];
   inertiadata->sumY2d += init[4]*init[1]*init[1];
   inertiadata->sumZ2d += init[4]*init[2]*init[2];
   inertiadata->sumXYd += init[4]*init[0]*init[1];
   inertiadata->sumYZd += init[4]*init[1]*init[2];
   inertiadata->sumXZd += init[4]*init[0]*init[2];
   inertiadata->sumvald += init[4];
   inertiadata->sumval2d += init[4]*init[4];*/

   // inertiadata->sumR2 += x*x + y*y + z*z;
} /* add_to_inertia_data */

void merge_inertiaeor_data(void *inertiaattr, void *childattr)
{
   InertiaDataEoR *inertiadata = inertiaattr;
   InertiaDataEoR *childdata = childattr;

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
   /* inertiadata->sumXd += childdata->sumXd;
   inertiadata->sumYd += childdata->sumYd;
   inertiadata->sumZd += childdata->sumZd;
   inertiadata->sumX2d += childdata->sumX2d;
   inertiadata->sumY2d += childdata->sumY2d;
   inertiadata->sumZ2d += childdata->sumZ2d;
   inertiadata->sumXYd += childdata->sumXYd;
   inertiadata->sumYZd += childdata->sumYZd;
   inertiadata->sumXZd += childdata->sumXZd;
   inertiadata->sumvald += childdata->sumvald;
   inertiadata->sumval2d += childdata->sumval2d;*/

} /* merge_inertia_data */

void correct_inertiaeor_data(void *inertiaattr, void *childattr)
{
   InertiaDataEoR *inertiadata = inertiaattr;
   InertiaDataEoR *childdata = childattr;
   double offset[3][2] ={0};
   double cen[2];

   cen[0] = inertiadata->sumX/inertiadata->sumval;
   cen[1] = childdata->sumX/childdata->sumval;
   
    if(cen[0] - cen[1] > (double) DIMS[0]/2){
     offset[0][0] = + (double) DIMS[0];
     offset[0][1] = - (double) DIMS[0];
   } else if (cen[0] - cen[1] < - (double) DIMS[0]/2){
     offset[0][0] =  - (double) DIMS[0];
     offset[0][1] =  + (double) DIMS[0];
     }

    if((inertiadata->sumX + childdata->sumX + childdata->sumval*offset[0][0]) / (inertiadata->sumval + childdata->sumval) < 0 || (inertiadata->sumX + childdata->sumX + childdata->sumval*offset[0][0]) / (inertiadata->sumval + childdata->sumval)  >= (double) DIMS[0])
     offset[0][0] = 0;
   else
     offset[0][1] = 0;

   cen[0] = inertiadata->sumY/inertiadata->sumval;
   cen[1] = childdata->sumY/childdata->sumval;
   
   if(cen[0] - cen[1] > (double) DIMS[1]/2){
     offset[1][0] =   (double) DIMS[1];
     offset[1][1] = - (double) DIMS[1];
   } else if (cen[0] - cen[1] < - (double) DIMS[1]/2){
     offset[1][0] = - (double) DIMS[1];
     offset[1][1] =   (double) DIMS[1];
     }
   
   if((inertiadata->sumY + childdata->sumY + childdata->sumval*offset[1][0]) / (inertiadata->sumval + childdata->sumval) < 0 || (inertiadata->sumY + childdata->sumY + childdata->sumval*offset[1][0]) / (inertiadata->sumval + childdata->sumval)  >= (double) DIMS[1])
     offset[1][0] = 0;
   else
     offset[1][1] = 0;

    cen[0] = inertiadata->sumZ/inertiadata->sumval;
    cen[1] = childdata->sumZ/childdata->sumval;
   
     if(cen[0] - cen[1] > (double) DIMS[2]/2){
     offset[2][0] =   (double) DIMS[2];
     offset[2][1] = - (double) DIMS[2];
   } else if (cen[0] - cen[1] < - (double) DIMS[2]/2){
     offset[2][0] = - (double) DIMS[2];
     offset[2][1] =   (double) DIMS[2];
     }

     if((inertiadata->sumZ + childdata->sumZ + childdata->sumval*offset[2][0]) / (inertiadata->sumval + childdata->sumval) < 0 || (inertiadata->sumZ + childdata->sumZ + childdata->sumval*offset[2][0]) / (inertiadata->sumval + childdata->sumval)  >= (double) DIMS[2])
     offset[2][0] = 0;
   else
     offset[2][1] = 0;

     
    inertiadata->sumX2 = childdata->sumX2 + childdata->sumval*offset[0][0]*offset[0][0] + 2*offset[0][0]*childdata->sumX +  inertiadata->sumX2 + inertiadata->sumval*offset[0][1]*offset[0][1] + 2*offset[0][1]*inertiadata->sumX;
    inertiadata->sumY2 = childdata->sumY2 + childdata->sumval*offset[1][0]*offset[1][0] + 2*offset[1][0]*childdata->sumY +  inertiadata->sumY2 + inertiadata->sumval*offset[1][1]*offset[1][1] + 2*offset[1][1]*inertiadata->sumY;
    inertiadata->sumZ2 = childdata->sumZ2 + childdata->sumval*offset[2][0]*offset[2][0] + 2*offset[2][0]*childdata->sumZ +  inertiadata->sumZ2 + inertiadata->sumval*offset[2][1]*offset[2][1] + 2*offset[2][1]*inertiadata->sumZ;

    /* inertiadata->sumX2d = childdata->sumX2d + childdata->sumvald*offset[0][0]*offset[0][0] + 2*offset[0][0]*childdata->sumXd +  inertiadata->sumX2d + inertiadata->sumvald*offset[0][1]*offset[0][1] + 2*offset[0][1]*inertiadata->sumXd;
    inertiadata->sumY2d = childdata->sumY2d + childdata->sumvald*offset[1][0]*offset[1][0] + 2*offset[1][0]*childdata->sumYd +  inertiadata->sumY2d + inertiadata->sumvald*offset[1][1]*offset[1][1] + 2*offset[1][1]*inertiadata->sumYd;
    inertiadata->sumZ2d = childdata->sumZ2d + childdata->sumvald*offset[2][0]*offset[2][0] + 2*offset[2][0]*childdata->sumZd +  inertiadata->sumZ2d + inertiadata->sumvald*offset[2][1]*offset[2][1] + 2*offset[2][1]*inertiadata->sumZd;
    */
    inertiadata->sumXY = childdata->sumXY + offset[1][0]*childdata->sumX + offset[0][0]*childdata->sumY + childdata->sumval*offset[0][0]*offset[1][0] + inertiadata->sumXY + offset[1][1]*inertiadata->sumX + offset[0][1]*inertiadata->sumY + inertiadata->sumval*offset[0][1]*offset[1][1];
    inertiadata->sumYZ = childdata->sumYZ + offset[1][0]*childdata->sumZ + offset[2][0]*childdata->sumY + childdata->sumval*offset[1][0]*offset[2][0] + inertiadata->sumYZ + offset[1][1]*inertiadata->sumZ + offset[2][1]*inertiadata->sumY + inertiadata->sumval*offset[2][1]*offset[1][1];
    inertiadata->sumXZ = childdata->sumXZ + offset[2][0]*childdata->sumX + offset[0][0]*childdata->sumZ + childdata->sumval*offset[2][0]*offset[0][0] + inertiadata->sumXZ + offset[2][1]*inertiadata->sumX + offset[0][1]*inertiadata->sumZ + inertiadata->sumval*offset[0][1]*offset[2][1];
    /*
    inertiadata->sumXYd = childdata->sumXYd + offset[1][0]*childdata->sumXd + offset[0][0]*childdata->sumYd + childdata->sumvald*offset[0][0]*offset[1][0] + inertiadata->sumXY + offset[1][1]*inertiadata->sumX + offset[0][1]*inertiadata->sumYd + inertiadata->sumvald*offset[0][1]*offset[1][1];
    inertiadata->sumYZd = childdata->sumYZd + offset[1][0]*childdata->sumZd + offset[2][0]*childdata->sumYd + childdata->sumvald*offset[1][0]*offset[2][0] + inertiadata->sumYZ + offset[1][1]*inertiadata->sumZ + offset[2][1]*inertiadata->sumYd + inertiadata->sumvald*offset[2][1]*offset[1][1];
    inertiadata->sumXZd = childdata->sumXZd + offset[2][0]*childdata->sumXd + offset[0][0]*childdata->sumZd + childdata->sumvald*offset[2][0]*offset[0][0] + inertiadata->sumXZ + offset[2][1]*inertiadata->sumX + offset[0][1]*inertiadata->sumZd + inertiadata->sumvald*offset[0][1]*offset[2][1];
    */
    inertiadata->sumX = childdata->sumX + childdata->sumval*offset[0][0] + inertiadata->sumX +  inertiadata->sumval*offset[0][1];
    inertiadata->sumY = childdata->sumY + childdata->sumval*offset[1][0] + inertiadata->sumY +  inertiadata->sumval*offset[1][1];
    inertiadata->sumZ = childdata->sumZ + childdata->sumval*offset[2][0] + inertiadata->sumZ +  inertiadata->sumval*offset[2][1];
    /*
    inertiadata->sumXd = childdata->sumXd + childdata->sumvald*offset[0][0] + inertiadata->sumXd +  inertiadata->sumvald*offset[0][1];
    inertiadata->sumYd = childdata->sumYd + childdata->sumvald*offset[1][0] + inertiadata->sumYd +  inertiadata->sumvald*offset[1][1];
    inertiadata->sumZd = childdata->sumZd + childdata->sumvald*offset[2][0] + inertiadata->sumZd +  inertiadata->sumvald*offset[2][1];*/
       
    inertiadata->area += childdata->area;
    inertiadata->sumval += childdata->sumval;
    inertiadata->sumval2 += childdata->sumval2;

    //    inertiadata->sumvald += childdata->sumvald;
    //  inertiadata->sumval2d += childdata->sumval2d;

} /* merge_inertia_data */

void merge_to_inertiaeor_data( AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr){
   InertiaDataEoR *thisdata = *thisattr;
   InertiaDataEoR *inertiadata = inertiaattr;
   InertiaDataEoR *childdata = childattr;

   if (thisdata==NULL) {
     thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataEoR));
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
} /* merge_to_inertia_data */

void clone_inertiaeor_data( AuxDataStore *store, void **thisattr, void *inertiaattr){
  InertiaDataEoR *thisdata = *thisattr;
  InertiaDataEoR *inertiadata = inertiaattr;

  if (thisdata==NULL) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataEoR));
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
  
  /*thisdata->sumXd = inertiadata->sumXd;
  thisdata->sumYd = inertiadata->sumYd;
  thisdata->sumZd = inertiadata->sumZd;
  thisdata->sumX2d = inertiadata->sumX2d;
  thisdata->sumY2d = inertiadata->sumY2d;
  thisdata->sumZ2d = inertiadata->sumZ2d;
  thisdata->sumXYd = inertiadata->sumXYd;
  thisdata->sumYZd = inertiadata->sumYZd;
  thisdata->sumXZd = inertiadata->sumXZd;
  thisdata->sumvald = inertiadata->sumvald;
  thisdata->sumval2d = inertiadata->sumval2d;*/
  //thisdata->sumR2 = inertiadata->sumR2;
} /* clone_inertia_data */

void *inertiaeor_attribute_arr(void *inertiaattr){
   InertiaDataEoR *inertiadata = inertiaattr;
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
   /*   inertia[12] = inertiadata->sumXd;
   inertia[13] = inertiadata->sumYd;
   inertia[14] = inertiadata->sumZd;
   inertia[15] = inertiadata->sumX2d;
   inertia[16] = inertiadata->sumY2d;
   inertia[17] = inertiadata->sumZ2d;
   inertia[18] = inertiadata->sumXYd;
   inertia[19] = inertiadata->sumYZd;
   inertia[20] = inertiadata->sumXZd;
   inertia[21] = inertiadata->sumvald;
   inertia[22] = inertiadata->sumval2d;*/
   return(inertia);
} /* inertia_attribute */

double inertiaeor_attribute(void *inertiaattr){
   InertiaDataEoR *inertiadata = inertiaattr;
   double inertia = inertiadata->area;

   return(inertia);
} /* inertia_attribute */


/****** Typedefs and functions for area attributes ******************************/

void *new_mto_data( AuxDataStore *store, double *init){
   MtoData *mtodata;
   mtodata = store ? get_new_aux_data(store) :  calloc(1, sizeof(MtoData));
   check_alloc(mtodata, 102);
   mtodata->area = 1;
   mtodata->gval = (value) init[0];
   mtodata->volume =0;
   mtodata->power =0;
   return(mtodata);
} /* new_area_data */


void add_to_mto_data(void *areaattr, double *init){
   MtoData *areadata = areaattr;
   areadata->area++;
   
} /* add_to_area_data */

void merge_mto_data(void *areaattr, void *childattr){
   MtoData *areadata = areaattr;
   MtoData *childdata = childattr;

   float delta = (float) childdata->gval - areadata->gval;
   areadata->area += childdata->area;

   childdata->power += delta * (2 * childdata->volume + delta * childdata->area);
   areadata->power += areadata->power;

   childdata->volume += delta * childdata->area;
   areadata->volume += childdata->volume;
} /* merge_area_data */

void merge_to_mto_data( AuxDataStore *store, void **thisattr, void *areaattr, void *childattr){
  MtoData *thisdata = *thisattr;
  MtoData *areadata = areaattr;
  MtoData *childdata = childattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(MtoData));
    check_alloc(thisdata, 104);
    *thisattr = thisdata;
  }
  thisdata->area = areadata->area + childdata->area;
} /* merge_to_area_data */

void clone_mto_data( AuxDataStore *store, void **thisattr, void *areaattr){

  MtoData *thisdata = *thisattr;
  MtoData *areadata = areaattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(MtoData));
    check_alloc(thisdata, 105);
    *thisattr = thisdata;
  }
  thisdata->area = areadata->area;
  thisdata->gval = areadata->gval;
  thisdata->volume = areadata->volume;
  thisdata->power = areadata->power;

} /* clone_area_data */

void *mto_attribute_arr(void *inertiaattr){
   MtoData *inertiadata = inertiaattr;
   double *inertia = calloc(3, sizeof(double));

   inertia[0] = inertiadata->area;
   inertia[1] = inertiadata->power;
   inertia[2] = inertiadata->volume;
   return(inertia);
} /* inertia_attribute */

double mto_attribute(void *inertiaattr){
   MtoData *inertiadata = inertiaattr;
   double inertia = inertiadata->area;

   return(inertia);
} /* inertia_attribute */


/*   */

/****** Typedefs and functions for moment of inertia attributes for EoR Application 2 *****************


void *new_eor_data( AuxDataStore *store, double *init){
  DataEoR *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(DataEoR));
  check_alloc(inertiadata, 114);
  inertiadata->area = 1;
  inertiadata->sumX = 0;
  inertiadata->sumY = 0;
  inertiadata->sumZ = 0;
  inertiadata->sumX2 = 0;
  inertiadata->sumY2 = 0;
  inertiadata->sumZ2 = 0;
  inertiadata->sumXY = 0;
  inertiadata->sumYZ = 0;
  inertiadata->sumXZ = 0;
  inertiadata->sumX = 0;
  inertiadata->sumY = 0;
  inertiadata->sumZ = 0;
  inertiadata->sumX2 = 0;
  inertiadata->sumY2 = 0;
  inertiadata->sumZ2 = 0;
  inertiadata->sumXY = 0;
  inertiadata->sumYZ = 0;
  inertiadata->sumXZ = 0;
  inertiadata->sumval = init[3];
  inertiadata->oriX = init[0];
  inertiadata->oriY = init[1];
  inertiadata->oriZ = init[2];

  return(inertiadata);
} 

void *load_eor_data( AuxDataStore *store, double *init){
  DataEoR *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(DataEoR));
  check_alloc(inertiadata, 115);
  
  inertiadata->area = init[0];
  inertiadata->sumX = init[1];
  inertiadata->sumY = init[2];
  inertiadata->sumZ = init[3];
  //inertiadata->sumR2 = init[4];
  return(inertiadata);
} 

void delete_eor_data(void *inertiaattr)
{
  free(inertiaattr);
} 

void add_to_eor_data(void *inertiaattr, double *init){
   DataEoR *inertiadata = inertiaattr;
   inertiadata->area++;
   inertiadata->sumval += init[3];
   inertiadata->sumX += init[3]*(init[0]-inertiadata->oriX);
   inertiadata->sumY += init[3]*(init[1]-inertiadata->oriY);
   inertiadata->sumZ += init[3]*(init[2]-inertiadata->oriZ);
   inertiadata->sumX2 += init[3]*(init[0]-inertiadata->oriX)*(init[0]-inertiadata->oriX);
   inertiadata->sumY2 += init[3]*(init[1]-inertiadata->oriY)*(init[1]-inertiadata->oriY);
   inertiadata->sumZ2 += init[3]*(init[2]-inertiadata->oriZ)*(init[2]-inertiadata->oriZ);
   inertiadata->sumXY += init[3]*(init[0]-inertiadata->oriX)*(init[1]-inertiadata->oriY);
   inertiadata->sumYZ += init[3]*(init[1]-inertiadata->oriY)*(init[2]-inertiadata->oriZ);
   inertiadata->sumXZ += init[3]*(init[0]-inertiadata->oriX)*(init[2]-inertiadata->oriZ);
   // inertiadata->sumR2 += x*x + y*y + z*z;
} 

void merge_eor_data(void *inertiaattr, void *childattr)
{
   DataEoR *inertiadata = inertiaattr;
   DataEoR *childdata = childattr;
   double offset[3][2], delta[3], newori[3];

   offset[0][0] = (inertiadata->oriX - childdata->oriX);
   offset[0][1] = - offset[0][0];
  
   
   if((inertiadata->sumX + childdata->sumX - childdata->area*offset[0][0]) / (inertiadata->area + childdata->area) + inertiadata->oriX < 0 || (inertiadata->sumX + childdata->sumX - childdata->area*offset[0][0]) / (inertiadata->area + childdata->area) + inertiadata->oriX >= DIMS[0]){
     newori[0] = childdata->oriX;
     offset[0][0] = 0;
   } else{
     newori[0] = inertiadata->oriX;
     offset[0][1] = 0;
   }
   
   offset[1][0] = -(childdata->oriY - inertiadata->oriY);
   offset[1][1] = 0;
     newori[1] = inertiadata->oriY;

   offset[2][0] = -(childdata->oriZ - inertiadata->oriZ);
   offset[2][1] = 0;
     newori[2] = inertiadata->oriZ;


   inertiadata->sumX = childdata->sumX - childdata->area*offset[0][0] + inertiadata->sumX -  inertiadata->area*offset[0][1];
   inertiadata->sumY = childdata->sumY - childdata->area*offset[1][0] + inertiadata->sumY -  inertiadata->area*offset[1][1];
   inertiadata->sumZ = childdata->sumZ - childdata->area*offset[2][0] + inertiadata->sumZ -  inertiadata->area*offset[2][1];

   inertiadata->sumX2 = childdata->sumX2 + childdata->area*offset[0][0]*offset[0][0] - 2*offset[0][0]*childdata->sumX +  inertiadata->sumX2 + inertiadata->area*offset[0][1]*offset[0][1] - 2*offset[0][1]*inertiadata->sumX;
   inertiadata->sumY2 = childdata->sumY2 + childdata->area*offset[1][0]*offset[1][0] - 2*offset[1][0]*childdata->sumY +  inertiadata->sumY2 + inertiadata->area*offset[1][1]*offset[1][1] - 2*offset[1][1]*inertiadata->sumY;
   inertiadata->sumZ2 = childdata->sumZ2 + childdata->area*offset[2][0]*offset[2][0] - 2*offset[2][0]*childdata->sumZ +  inertiadata->sumZ2 + inertiadata->area*offset[2][1]*offset[2][1] - 2*offset[2][1]*inertiadata->sumZ;

   inertiadata->sumXY = childdata->sumXY - offset[1][0]*childdata->sumX - offset[0][0]*childdata->sumY + childdata->area*offset[0][0]*offset[1][0] + inertiadata->sumXY - offset[1][1]*inertiadata->sumX - offset[0][1]*inertiadata->sumY + inertiadata->area*offset[0][1]*offset[1][1];
   inertiadata->sumYZ = childdata->sumYZ - offset[1][0]*childdata->sumZ - offset[2][0]*childdata->sumY + childdata->area*offset[1][0]*offset[2][0] + inertiadata->sumYZ - offset[1][1]*inertiadata->sumZ - offset[2][1]*inertiadata->sumY + inertiadata->area*offset[2][1]*offset[1][1];
   inertiadata->sumXZ = childdata->sumXZ - offset[2][0]*childdata->sumX - offset[0][0]*childdata->sumZ + childdata->area*offset[2][0]*offset[0][0] + inertiadata->sumXZ - offset[2][1]*inertiadata->sumX - offset[0][1]*inertiadata->sumZ + inertiadata->area*offset[0][1]*offset[2][1];
       
   inertiadata->area += childdata->area;
   inertiadata->sumval += childdata->sumval;
   
   inertiadata->oriX = newori[0];
   inertiadata->oriY = newori[1];
   inertiadata->oriZ = newori[2];

} 

void merge_to_eor_data( AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr){
   DataEoR *thisdata = *thisattr;
   DataEoR *inertiadata = inertiaattr;
   DataEoR *childdata = childattr;
     double offset[3];
   offset[0] = ( -inertiadata->oriX + childdata->oriX);
   offset[1] = (childdata->oriY - inertiadata->oriY);
   offset[2] = (childdata->oriZ - inertiadata->oriZ);
   if (thisdata==NULL) {
     thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(DataEoR));
     check_alloc(thisdata, 116);
     *thisattr = thisdata;
   }
   thisdata->area = inertiadata->area + childdata->area;
   thisdata->sumval = inertiadata->sumval + childdata->sumval;

   thisdata->sumX = inertiadata->sumX + childdata->sumX + childdata->area*offset[0];
   thisdata->sumY = inertiadata->sumY + childdata->sumY + childdata->area*offset[1];
   thisdata->sumZ = inertiadata->sumZ + childdata->sumZ + childdata->area*offset[2];

   thisdata->sumX2 = inertiadata->sumX2 + childdata->sumX2 + childdata->area*offset[0]*offset[0] + 2*offset[0]*childdata->sumX;
   thisdata->sumY2 = inertiadata->sumY2 + childdata->sumY2 + childdata->area*offset[1]*offset[1] + 2*offset[1]*childdata->sumY;
   thisdata->sumZ2 = inertiadata->sumZ2 + childdata->sumZ2 + childdata->area*offset[2]*offset[2] + 2*offset[2]*childdata->sumZ;
   thisdata->sumXY = inertiadata->sumXY + childdata->sumXY + offset[1]*childdata->sumX + offset[0]*childdata->sumY + childdata->area*offset[0]*offset[1];
   thisdata->sumYZ = inertiadata->sumYZ + childdata->sumYZ + offset[1]*childdata->sumZ + offset[2]*childdata->sumY + childdata->area*offset[1]*offset[2];
   thisdata->sumXZ = inertiadata->sumXZ + childdata->sumXZ + offset[2]*childdata->sumX + offset[0]*childdata->sumZ + childdata->area*offset[0]*offset[2];
   thisdata->oriX = inertiadata->oriX;
   thisdata->oriY = inertiadata->oriY;
   thisdata->oriZ = inertiadata->oriZ;

   // thisdata->sumR2 = inertiadata->sumR2 + childdata->sumR2;
} 
void clone_eor_data( AuxDataStore *store, void **thisattr, void *inertiaattr){
  DataEoR *thisdata = *thisattr;
  DataEoR *inertiadata = inertiaattr;

  if (thisdata==NULL) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(DataEoR));
    check_alloc(thisdata, 117);
    *thisattr = thisdata;
  }
  thisdata->area = inertiadata->area;
  thisdata->sumval = inertiadata->sumval;
    
  thisdata->sumX = inertiadata->sumX;
  thisdata->sumY = inertiadata->sumY;
  thisdata->sumZ = inertiadata->sumZ;
  thisdata->sumX2 = inertiadata->sumX2;
  thisdata->sumY2 = inertiadata->sumY2;
  thisdata->sumZ2 = inertiadata->sumZ2;
  thisdata->sumXY = inertiadata->sumXY;
  thisdata->sumYZ = inertiadata->sumYZ;
  thisdata->sumXZ = inertiadata->sumXZ;
  thisdata->oriX = inertiadata->oriX;
  thisdata->oriY = inertiadata->oriY;
  thisdata->oriZ = inertiadata->oriZ;
  //thisdata->sumR2 = inertiadata->sumR2;
} 


void *eor_attribute(void *inertiaattr){
   DataEoR *inertiadata = inertiaattr;
   double *inertia = calloc(14, sizeof(double));

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
   inertia[11] = inertiadata->oriX;
   inertia[12] = inertiadata->oriY;
   inertia[13] = inertiadata->oriZ;

   return(inertia);
}


*/





AttribStruct AttribsArray[NUMATTR] =
{
  {"Area", sizeof(AreaData),  new_area_data, delete_area_data, add_to_area_data, merge_area_data, merge_to_area_data, clone_area_data, create_mpi_area_type, area_attribute, NULL},
  {"Area of min. enclosing rectangle", sizeof(EnclRectData),  new_encl_rect_data, delete_encl_rect_data, add_to_encl_rect_data, merge_encl_rect_data, merge_to_encl_rect_data, clone_encl_rect_data, create_mpi_rect_type, encl_rect_area_attribute,  NULL},
  {"Square of diagonal of min. enclosing rectangle", sizeof(EnclRectData),0, new_encl_rect_data, delete_encl_rect_data, add_to_encl_rect_data, merge_encl_rect_data, merge_to_encl_rect_data, clone_encl_rect_data, create_mpi_rect_type, encl_rect_diag_attribute, NULL},
  {"Moment of Inertia", sizeof(InertiaData), new_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data, create_mpi_inertia_type, inertia_attribute, NULL},
  {"Moment of Inertia for EoR", sizeof(InertiaDataEoR),  new_inertiaeor_data, delete_inertiaeor_data, add_to_inertiaeor_data, merge_inertiaeor_data, merge_to_inertiaeor_data, clone_inertiaeor_data, create_mpi_inertiaeor_type,inertiaeor_attribute, inertiaeor_attribute_arr},
  {"(Moment of Inertia) / (area)^2", sizeof(InertiaData),  new_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data,  create_mpi_inertia_type, inertia_div_a2_attribute, NULL},
  {"Mean X position", sizeof(InertiaData), new_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data,  create_mpi_inertia_type, mean_x_attribute, NULL},
  {"Mean Y position", sizeof(InertiaData),  new_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data,  create_mpi_inertia_type, mean_y_attribute, NULL},
  {"Mean Z position", sizeof(InertiaData),  new_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data,  create_mpi_inertia_type, mean_z_attribute, NULL},
  {"MTO", sizeof(MtoData),  new_mto_data, delete_area_data, add_to_mto_data, merge_mto_data, merge_to_mto_data, clone_mto_data, create_mpi_mto_type, mto_attribute, mto_attribute_arr},
};


DecisionStruct Decisions[NUMDECISIONS] =
{
  {"Direct", tree_filter_direct},
  {"Min", tree_filter_min},  
  {"Max", tree_filter_max},
  {"Subtractive", tree_filter_subtractive},
};
