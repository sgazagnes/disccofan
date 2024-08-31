#include "types.h"
#include "checks.h"
#include "filter.h"
#include "eispack.h"
#include "attributes.h"

//double pixdim[3];


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
    debug(" Reallocation of attribute storage (may cause errors), size item %ld",  store->size_item);
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
/* FUNCTIONS TO DEFINE ARE */

/****** Typedefs and functions for area attributes ******************************/

void *new_area_data(AuxDataStore *store){
   AreaData *areadata;
   areadata = store ? get_new_aux_data(store) :  calloc(1, sizeof(AreaData));
   check_alloc(areadata, 102);
   return(areadata);
} /* new_area_data */

void init_area_data(void *areaattr, bool init, ulong x, ulong y, ulong z, value gval){
  AreaData *areadata = areaattr;
  areadata->area = (ulong) init;
}


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
   area = pixdim[0]*pixdim[1]*pixdim[2]*areadata->area;
   return(area);
} /* area_attribute */




/****** Typedefs and functions for extent and  enclosing rectangle attributes *******/

void *new_extent_data( AuxDataStore *store){
   ExtentData *extentdata;
   extentdata = store ? get_new_aux_data(store) :  calloc(1, sizeof(ExtentData));
   check_alloc(extentdata, 106);
   return(extentdata);
} /* new_extent_data */


void init_extent_data(void *rectattr, bool init, ulong x, ulong y, ulong z, value gval){
  ExtentData *extentdata = rectattr;
  if(init){
    extentdata->minX = x;
    extentdata->maxX = x+1;
    extentdata->minY = y;
    extentdata->maxY = y+1;
    extentdata->minZ = z;
    extentdata->maxZ = z+1;
  } else {
    extentdata->minX = extentdata->minY = extentdata->minZ = ULONG_MAX;
    extentdata->maxX = extentdata->maxY = extentdata->maxZ = 0;
  }
}

void delete_extent_data(void *rectattr)
{
  free(rectattr);
} /* delete_extent_data */

void add_to_extent_data(void *rectattr, void *data){
   ExtentData *extentdata = rectattr;
   double *values = (double *) data; 

   extentdata->minX = MIN(extentdata->minX, (ulong) values[0]);
   extentdata->minY = MIN(extentdata->minY, (ulong) values[1]);
   extentdata->minZ = MIN(extentdata->minZ, (ulong) values[2]);
   extentdata->maxX = MAX(extentdata->maxX, (ulong) values[0]);
   extentdata->maxY = MAX(extentdata->maxY, (ulong) values[1]);
   extentdata->maxZ = MAX(extentdata->maxZ, (ulong) values[2]);
} /* add_to_extent_data */


void merge_extent_data(void *rectattr, void *childattr){
   ExtentData *extentdata = rectattr;
   ExtentData *childdata = childattr;

   extentdata->minX = MIN(extentdata->minX, childdata->minX);
   extentdata->minY = MIN(extentdata->minY, childdata->minY);
   extentdata->minZ = MIN(extentdata->minZ, childdata->minZ);
   extentdata->maxX = MAX(extentdata->maxX, childdata->maxX);
   extentdata->maxY = MAX(extentdata->maxY, childdata->maxY);
   extentdata->maxZ = MAX(extentdata->maxZ, childdata->maxZ);
} /* merge_extent_data */


void merge_to_extent_data( AuxDataStore *store, void **thisattr, void *rectattr, void *childattr){
  ExtentData *thisdata = *thisattr;
  ExtentData *extentdata = rectattr;
  ExtentData *childdata = childattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(ExtentData));
    check_alloc(thisdata, 108);
    *thisattr = thisdata;
  }
  thisdata->minX = MIN(extentdata->minX, childdata->minX);
  thisdata->minY = MIN(extentdata->minY, childdata->minY);
  thisdata->minZ = MIN(extentdata->minZ, childdata->minZ);
  thisdata->maxX = MAX(extentdata->maxX, childdata->maxX);
  thisdata->maxY = MAX(extentdata->maxY, childdata->maxY);
  thisdata->maxZ = MAX(extentdata->maxZ, childdata->maxZ);
} /* merge_to_extent_data */

void clone_extent_data( AuxDataStore *store, void **thisattr, void *rectattr){
  ExtentData *thisdata = *thisattr;
  ExtentData *extentdata = rectattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(ExtentData));
    check_alloc(thisdata, 109);
    *thisattr = thisdata;
  }
  thisdata->minX = extentdata->minX;
  thisdata->minY = extentdata->minY;
  thisdata->minZ = extentdata->minZ;
  thisdata->maxX = extentdata->maxX;
  thisdata->maxY = extentdata->maxY;
  thisdata->maxZ = extentdata->maxZ;
} /* clone_extent_data */

double extent_rectarea_attribute(void *rectattr){
  ExtentData *extentdata = rectattr;
  ulong volume;
  volume =  pixdim[0]*pixdim[1]*pixdim[2]*
    (extentdata->maxX - extentdata->minX) 
    * (extentdata->maxY - extentdata->minY)
    * (extentdata->maxZ - extentdata->minZ);
  return (double) volume;
} /* extent_area_attribute */



double extent_rectdiag_attribute(void *rectattr){
/* Computes the square of the length of the diagonal */
  ExtentData *extentdata = rectattr;
  ulong minx, miny, minz, maxx, maxy, maxz, l;

  minx = pixdim[0]*extentdata->minX;
  miny = pixdim[1]*extentdata->minY;
  minz = pixdim[2]*extentdata->minZ;
  maxx = pixdim[0]*extentdata->maxX;
  maxy = pixdim[1]*extentdata->maxY;
  maxz = pixdim[2]*extentdata->maxZ;
  l = (maxx-minx) * (maxx-minx)
    + (maxy-miny) * (maxy-miny)
    + (maxz-minz) * (maxz-minz);
  return(double) (l);
} /* extent_diag_attribute */

double extent_x_attribute(void *rectattr){
  /* Computes the X-extent */
  ExtentData *extentdata = rectattr;
  return (double) ( pixdim[0]*(extentdata->maxX - extentdata->minX));
} /* extent_xextent_attribute */

double extent_y_attribute(void *rectattr){
  /* Computes the Y-extent */
  ExtentData *extentdata = rectattr;
  return (double) ( pixdim[1]*(extentdata->maxY - extentdata->minY));
} /* extent_yextent_attribute */


double extent_z_attribute(void *rectattr){
  /* Computes the Z-extent */
  ExtentData *extentdata = rectattr;
  return (double) ( pixdim[2]*(extentdata->maxZ - extentdata->minZ));
} /* extent_zextent_attribute */



/****** Typedefs and functions for mean attributes **************************/


void *new_mean_data( AuxDataStore *store){
  MeanData *meandata;
  meandata = store ? get_new_aux_data(store) : calloc(1, sizeof(MeanData));
  check_alloc(meandata, 110);
  return(meandata);
} /* new_mean_data */

void init_mean_data(void *meanattr, bool init, ulong x, ulong y, ulong z, value gval){
  MeanData *meandata = meanattr;
  meandata->area = (ulong) init;
  meandata->sumX = init ? x : 0.;
  meandata->sumY = init ? y : 0.;
  meandata->sumZ = init ? z : 0.;
}/* init_mean_data*/

void delete_mean_data(void *meanattr)
{
  free(meanattr);
} /* delete_mean_data */

void add_to_mean_data(void *meanattr, void *data){
   MeanData *meandata = meanattr;
   double *values = (double *) data;
   meandata->area++;
   meandata->sumX += (ulong) values[0];
   meandata->sumY += (ulong) values[1];
   meandata->sumZ += (ulong) values[2];
} /* add_to_mean_data */

void merge_mean_data(void *meanattr, void *childattr)
{
   MeanData *meandata = meanattr;
   MeanData *childdata = childattr;
   meandata->area += childdata->area;
   meandata->sumX += childdata->sumX;
   meandata->sumY += childdata->sumY;
   meandata->sumZ += childdata->sumZ;
} /* merge_mean_data */

void merge_to_mean_data( AuxDataStore *store, void **thisattr, void *meanattr, void *childattr){
   MeanData *thisdata = *thisattr;
   MeanData *meandata = meanattr;
   MeanData *childdata = childattr;

   if (!thisdata) {
     thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(MeanData));
     check_alloc(thisdata, 112);
     *thisattr = thisdata;
   }
   thisdata->area = meandata->area + childdata->area;
   thisdata->sumX = meandata->sumX + childdata->sumX;
   thisdata->sumY = meandata->sumY + childdata->sumY;
   thisdata->sumZ = meandata->sumZ + childdata->sumZ;
} /* merge_to_mean_data */

void clone_mean_data( AuxDataStore *store, void **thisattr, void *meanattr){
  MeanData *thisdata = *thisattr;
  MeanData *meandata = meanattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(MeanData));
    check_alloc(thisdata, 113);
    *thisattr = thisdata;
  }
  thisdata->area = meandata->area;
  thisdata->sumX = meandata->sumX;
  thisdata->sumY = meandata->sumY;
  thisdata->sumZ = meandata->sumZ;
} /* clone_mean_data */

double mean_area_attribute(void *meanattr){
   MeanData *meandata = meanattr;
   return((double)pixdim[0]*pixdim[1]*pixdim[2]*meandata->area);
} /* mean_area */

double mean_x_attribute(void *meanattr){
   MeanData *meandata = meanattr;
   return((double) pixdim[0]*meandata->sumX / (double) meandata->area);
} /* mean_x_attribute */

double mean_y_attribute(void *meanattr){
   MeanData *meandata = meanattr;
   return((double) pixdim[1]*meandata->sumY / (double)  meandata->area);
} /* mean_y_attribute */

double mean_z_attribute(void *meanattr){
   MeanData *meandata = meanattr;
   return((double) pixdim[2]*meandata->sumZ / (double) meandata->area);
} /* mean_z_attribute */


/* Weighted mean data attribute */


void *new_wmean_data( AuxDataStore *store){
  WMeanData *wmeandata;
  wmeandata = store ? get_new_aux_data(store) : calloc(1, sizeof(WMeanData));
  check_alloc(wmeandata, 110);
  return(wmeandata);
} /* new_wmean_data */

void init_wmean_data(void *wmeanattr, bool init, ulong x, ulong y, ulong z, value gval){
  WMeanData *wmeandata = wmeanattr;
  wmeandata->sumGval = (double) init*gval;
  wmeandata->sumXd = init ? (double) gval*x : 0.;
  wmeandata->sumYd = init ? (double) gval*y : 0.;
  wmeandata->sumZd = init ? (double) gval*z : 0.;
}/* init_wmean_data*/

void delete_wmean_data(void *wmeanattr)
{
  free(wmeanattr);
} /* delete_wmean_data */

void add_to_wmean_data(void *wmeanattr, void *data){
   WMeanData *wmeandata = wmeanattr;
   double *values = (double *) data;
   wmeandata->sumGval+= (double) values[3];
   wmeandata->sumXd += (double) values[3]*values[0];
   wmeandata->sumYd += (double) values[3]*values[1];
   wmeandata->sumZd += (double) values[3]*values[2];
} /* add_to_wmean_data */

void merge_wmean_data(void *wmeanattr, void *childattr)
{
   WMeanData *wmeandata = wmeanattr;
   WMeanData *childdata = childattr;
   wmeandata->sumGval += childdata->sumGval;
   wmeandata->sumXd += childdata->sumXd;
   wmeandata->sumYd += childdata->sumYd;
   wmeandata->sumZd += childdata->sumZd;
} /* merge_wmean_data */

void merge_to_wmean_data( AuxDataStore *store, void **thisattr, void *wmeanattr, void *childattr){
   WMeanData *thisdata = *thisattr;
   WMeanData *wmeandata = wmeanattr;
   WMeanData *childdata = childattr;

   if (!thisdata) {
     thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(WMeanData));
     check_alloc(thisdata, 112);
     *thisattr = thisdata;
   }
   thisdata->sumGval = wmeandata->sumGval + childdata->sumGval;
   thisdata->sumXd = wmeandata->sumXd + childdata->sumXd;
   thisdata->sumYd = wmeandata->sumYd + childdata->sumYd;
   thisdata->sumZd = wmeandata->sumZd + childdata->sumZd;
} /* merge_to_wmean_data */

void clone_wmean_data( AuxDataStore *store, void **thisattr, void *wmeanattr){
  WMeanData *thisdata = *thisattr;
  WMeanData *wmeandata = wmeanattr;

  if (!thisdata) {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(WMeanData));
    check_alloc(thisdata, 113);
    *thisattr = thisdata;
  }
  thisdata->sumGval = wmeandata->sumGval;
  thisdata->sumXd = wmeandata->sumXd;
  thisdata->sumYd = wmeandata->sumYd;
  thisdata->sumZd = wmeandata->sumZd;
} /* clone_wmean_data */

double wmean_totalflux_attribute(void *wmeanattr){
   WMeanData *wmeandata = wmeanattr;
   return((double)wmeandata->sumGval);
} /* mean_area */


double wmean_x_attribute(void *wmeanattr){
   WMeanData *wmeandata = wmeanattr;
   return((double) pixdim[0]*wmeandata->sumXd / (double) wmeandata->sumGval);
} /* mean_x_attribute */

double wmean_y_attribute(void *wmeanattr){
   WMeanData *wmeandata = wmeanattr;
   return((double) pixdim[1]*wmeandata->sumYd / (double) wmeandata->sumGval);
} /* mean_y_attribute */

double wmean_z_attribute(void *wmeanattr){
   WMeanData *wmeandata = wmeanattr;
   return((double) pixdim[2]*wmeandata->sumZd / (double)  wmeandata->sumGval);
} /* mean_z_attribute */



/*------------------------ */
/* Moment of Inertia attributes + analysis */

void *new_inertia_data( AuxDataStore *store){
  InertiaData *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaData));
  check_alloc(inertiadata, 114);
  return(inertiadata);
} /* new_inertia_data */


void init_inertia_data(void *inertiaattr, bool init, ulong x, ulong y, ulong z, value gval){
  InertiaData *inertiadata = inertiaattr;
  //  double *init = (double *) data;
  inertiadata->area   = (ulong) init;
  inertiadata->sumX   = init ? (double)  x : 0.;
  inertiadata->sumY   = init ? (double)  y : 0.;
  inertiadata->sumZ   = init ? (double)  z : 0.;
  inertiadata->sumX2  = init ? (double)  x * x: 0.;
  inertiadata->sumY2  = init ? (double)  y * y: 0.;
  inertiadata->sumZ2  = init ? (double)  z * z: 0.;
  inertiadata->sumXY  = init ? (double)  x * y: 0.;
  inertiadata->sumYZ  = init ? (double)  y * z: 0.;
  inertiadata->sumXZ  = init ? (double)  z * x: 0.;
  inertiadata->sumval = init ? (double)  gval: 0;
  inertiadata->sumval2 = init ? (double) gval * gval: 0;
  inertiadata->sumXd  = init ? (double)  gval*x : 0.;
  inertiadata->sumYd  = init ? (double)  gval*y : 0.;
  inertiadata->sumZd  = init ? (double)  gval*z : 0.;
}

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
   inertiadata->sumX2 += values[0]*values[0];
   inertiadata->sumY2 += values[1]*values[1];
   inertiadata->sumZ2 += values[2]*values[2];
   inertiadata->sumXY += values[0]*values[1];
   inertiadata->sumYZ += values[1]*values[2];
   inertiadata->sumXZ += values[0]*values[2];
   inertiadata->sumval += values[3];
   inertiadata->sumval2 += values[3]*values[3];
   inertiadata->sumXd += values[3]*values[0];
   inertiadata->sumYd += values[3]*values[1];
   inertiadata->sumZd += values[3]*values[2];
} /* add_to_inertia_data */

void merge_inertia_data(void *inertiaattr, void *childattr)
{
   InertiaData *inertiadata = inertiaattr;
   InertiaData *childdata = childattr;
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
   inertiadata->sumXd += childdata->sumXd;
   inertiadata->sumYd += childdata->sumYd;
   inertiadata->sumZd += childdata->sumZd;
   //inertiadata->sumval2 += childdata->sumval2;
} /* merge_inertia_data */


void merge_to_inertia_data( AuxDataStore *store, void **thisattr, void *inertiaattr, void *childattr){
   InertiaData *thisdata = *thisattr;
   InertiaData *inertiadata = inertiaattr;
   InertiaData *childdata = childattr;

   if (!thisdata)  {
     thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaData));
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
   thisdata->sumXd = inertiadata->sumXd + childdata->sumXd;
   thisdata->sumYd = inertiadata->sumYd + childdata->sumYd;
   thisdata->sumZd = inertiadata->sumZd + childdata->sumZd;
} /* merge_to_inertia_data */

void clone_inertia_data( AuxDataStore *store, void **thisattr, void *inertiaattr){
  InertiaData *thisdata = *thisattr;
  InertiaData *inertiadata = inertiaattr;

   if (!thisdata)  {
    thisdata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaData));
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
  thisdata->sumXd = inertiadata->sumXd;
  thisdata->sumYd = inertiadata->sumYd;
  thisdata->sumZd = inertiadata->sumZd;
  
} /* clone_inertia_data */

void *inertia_arr(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;
   double *inertia = calloc(15, sizeof(double));
   inertia[0] = (double) pixdim[0]*pixdim[1]*pixdim[2]*inertiadata->area;
   inertia[1] = pixdim[0]*inertiadata->sumX;
   inertia[2] = pixdim[1]*inertiadata->sumY;
   inertia[3] = pixdim[2]*inertiadata->sumZ;
   inertia[4] = pixdim[0]*pixdim[0]*inertiadata->sumX2;
   inertia[5] = pixdim[1]*pixdim[1]*inertiadata->sumY2;
   inertia[6] = pixdim[2]*pixdim[2]*inertiadata->sumZ2;
   inertia[7] = pixdim[0]*pixdim[1]*inertiadata->sumXY;
   inertia[8] = pixdim[1]*pixdim[2]*inertiadata->sumYZ;
   inertia[9] = pixdim[0]*pixdim[2]*inertiadata->sumXZ;
   inertia[10] = inertiadata->sumval;
   inertia[11] = inertiadata->sumval2;
   inertia[12] = pixdim[0]*inertiadata->sumXd;
   inertia[13] = pixdim[1]*inertiadata->sumYd;
   inertia[14] = pixdim[2]*inertiadata->sumZd;
   return(inertia);
} /* inertia_array */


void *inertia_attribute_arr(void *inertiaattr){
  //  InertiaData *inertiadata = inertiaattr;
  double *inertia = inertia_arr(inertiaattr);
  double *attrArr = calloc(13, sizeof(double));
  
  attrArr[0] = inertia[0];  // Area
  attrArr[1] = inertia[10]; // Sum of gray values
  attrArr[2] = inertia[11]; // Sum of gray values squared
  attrArr[3] = inertia[1]/inertia[0]; // Mean X
  attrArr[4] = inertia[2]/inertia[0]; // Mean Y
  attrArr[5] = inertia[3]/inertia[0]; // Mean Z
  attrArr[6] = inertia[11]/inertia[10]; // Weighted Mean X
  attrArr[7] = inertia[12]/inertia[10]; // Weighted Mean Y
  attrArr[8] = inertia[13]/inertia[10]; // Weighted Mean Z
  
  double corrVal    = inertia[0]/12;
  
  if(inertia[3] == 0){		       
     double matrix[3] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] +  pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] +  pixdim[1]*corrVal,
       inertia[7] - inertia[1]*inertia[2]/inertia[0] };		
   
     double *eigval  = malloc(2  * sizeof(double));
     double *eigvec  = malloc(4  * sizeof(double));

     double tens_mat[4] = {matrix[0], matrix[2], matrix[2], matrix[1]};
     rs (2, tens_mat, eigval, 1, eigvec);
     
     attrArr[9] = eigval[1]/eigval[0]; // Elong
     attrArr[10] = -1;  // Flatness

          
     double ax_len[2] = {sqrt(4*eigval[1]/inertia[0]), sqrt(4*eigval[0]/inertia[0])};
     attrArr[11] =  M_PI*ax_len[0]*ax_len[1]/(inertia[0]); // Sparseness
     attrArr[12] =  (matrix[0]+matrix[1])/pow(inertia[0], 2.0);  // Non-compactness
     
     free(eigval);
     free(eigvec);

  } else {
     double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] +  pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] +  pixdim[1]*corrVal,
       inertia[6] - inertia[3]*inertia[3]/inertia[0] +  pixdim[2]*corrVal,
       inertia[7] - inertia[1]*inertia[2]/inertia[0],
       inertia[8] - inertia[2]*inertia[3]/inertia[0],
       inertia[9] - inertia[1]*inertia[3]/inertia[0] };		
   
     double *eigval  = malloc(3  * sizeof(double));
     double *eigvec  = malloc(9  * sizeof(double));
  
     double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
     rs (3, tens_mat, eigval, 1, eigvec);

     attrArr[9] = eigval[2]/eigval[1];  // Elong
     attrArr[10] = eigval[1]/eigval[0];  // Flat
     
     double ax_len[3] = {sqrt(20*eigval[2]/inertia[0]), sqrt(20*eigval[1]/inertia[0]), sqrt(20*eigval[0]/inertia[0])};
     attrArr[11] =  M_PI*ax_len[0]*ax_len[1]*ax_len[2]/(6*inertia[0]); // Sparseness
     attrArr[12] =  (matrix[0]+matrix[1]+matrix[2])/pow(inertia[0], 5.0/3.0); // Non-compactness
     
     free(eigval);
     free(eigvec);
   }
   
   free(inertia);
   return(attrArr);
} /* inertia_attribute_arr */



double inertia_area_attribute(void *inertiaattr){
   InertiaData *inertiadata = inertiaattr;
   double area = (double) pixdim[0]*pixdim[1]*pixdim[2]*inertiadata->area;
   return(area);
} /* inertia_area_attribute */


double  inertia_trace_attribute(void *inertiaattr){
  //   InertiaData *inertiadata = inertiaattr;
   double *inertia   = inertia_arr(inertiaattr);
   double corrVal    = inertia[0]/12;

   double trace = 0;
   if(inertia[3] == 0){
     double matrix[2] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] + pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] + pixdim[1]*corrVal };
     trace    =  (matrix[0]+matrix[1]);
   } else {
     double matrix[3] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] + pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] + pixdim[1]*corrVal,
       inertia[6] - inertia[3]*inertia[3]/inertia[0] + pixdim[2]*corrVal };
      trace    =  (matrix[0]+matrix[1]+matrix[2]);
   }
   free(inertia);
   return trace;
} /* inertia_trace_attribute */


double  inertia_elon_attribute(void *inertiaattr){
   double *inertia  = inertia_arr(inertiaattr); 
   double corrVal   = inertia[0]/12;
   double elong = 0;
   
   if(inertia[3] == 0){
		       
     double matrix[3] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] + pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] +  pixdim[1]*corrVal,
       inertia[7] - inertia[1]*inertia[2]/inertia[0] };		
   
     double *eigval  = malloc(2  * sizeof(double));
     double *eigvec  = malloc(4  * sizeof(double));

     double tens_mat[4] = {matrix[0], matrix[2], matrix[2], matrix[1]};
     rs (2, tens_mat, eigval, 1, eigvec);

     elong = eigval[1]/eigval[0];
     free(eigval);
     free(eigvec);
     
   } else {
     double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] +  pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] +  pixdim[1]*corrVal,
       inertia[6] - inertia[3]*inertia[3]/inertia[0] +  pixdim[2]*corrVal,
       inertia[7] - inertia[1]*inertia[2]/inertia[0],
       inertia[8] - inertia[2]*inertia[3]/inertia[0],
       inertia[9] - inertia[1]*inertia[3]/inertia[0] };		
     double *eigval  = malloc(3  * sizeof(double));
     double *eigvec  = malloc(9  * sizeof(double));
  
     double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
     rs (3, tens_mat, eigval, 1, eigvec);

     elong = eigval[2]/eigval[1];
     free(eigval);
     free(eigvec);
   }
   
   free(inertia);
   return elong;
} /* inertia_elon_attribute */


double  inertia_flat_attribute(void *inertiaattr){
  //InertiaData *inertiadata = inertiaattr;
   double *inertia  = inertia_arr(inertiaattr);
   
   //double c_xyz[3]  = { inertia[1]/inertia[0],inertia[2]/inertia[0],inertia[3]/inertia[0]} ;
   double corrVal    = inertia[0]/12;

   if(inertia[3] == 0){
     free(inertia);
     return -1;
   }
   
   double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] + pixdim[0]*corrVal,
     inertia[5] - inertia[2]*inertia[2]/inertia[0] + pixdim[1]*corrVal,
     inertia[6] - inertia[3]*inertia[3]/inertia[0] + pixdim[2]*corrVal,
     inertia[7] - inertia[1]*inertia[2]/inertia[0],
     inertia[8] - inertia[2]*inertia[3]/inertia[0],
     inertia[9] - inertia[1]*inertia[3]/inertia[0] };		
   
   double *eigval  = malloc(3  * sizeof(double));
   double *eigvec  = malloc(9  * sizeof(double));
  
   double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
    rs (3, tens_mat, eigval, 1, eigvec);
   double flat = eigval[1]/eigval[0];
   free(eigval);
   free(eigvec);
   free(inertia);
   
   return(flat);
} /* inertia_flat_attribute */



double inertia_spar_attribute(void *inertiaattr){
  //   InertiaData *inertiadata = inertiaattr;
   double *inertia  = inertia_arr(inertiaattr);
   double corrVal    =  inertia[0]/12;
   double spars = 0;
   
   if(inertia[3] == 0){
		       
     double matrix[3] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] + pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] + pixdim[1]*corrVal,
       inertia[7] - inertia[1]*inertia[2]/inertia[0] };		
   
     double *eigval  = malloc(2  * sizeof(double));
     double *eigvec  = malloc(4  * sizeof(double));
  
     double tens_mat[4] = {matrix[0], matrix[2], matrix[2], matrix[1]};
     rs (2, tens_mat, eigval, 1, eigvec);
     double ax_len[2] = {sqrt(4*eigval[1]/inertia[0]), sqrt(4*eigval[0]/inertia[0])};
     spars =  M_PI*ax_len[0]*ax_len[1]/(inertia[0]);
     free(eigval);
     free(eigvec);

   }else {
     
     double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] + pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] + pixdim[1]*corrVal,
       inertia[6] - inertia[3]*inertia[3]/inertia[0] + pixdim[2]*corrVal,
       inertia[7] - inertia[1]*inertia[2]/inertia[0],
       inertia[8] - inertia[2]*inertia[3]/inertia[0],
       inertia[9] - inertia[1]*inertia[3]/inertia[0] };		
   
     double *eigval  = malloc(3  * sizeof(double));
     double *eigvec  = malloc(9  * sizeof(double));
  
     double tens_mat[9] = {matrix[0], matrix[3], matrix[5], matrix[3], matrix[1],matrix[4], matrix[5], matrix[4],matrix[2]};
     rs (3, tens_mat, eigval, 1, eigvec);
     double ax_len[3] = {sqrt(20*eigval[2]/inertia[0]), sqrt(20*eigval[1]/inertia[0]), sqrt(20*eigval[0]/inertia[0])};
     spars =  M_PI*ax_len[0]*ax_len[1]*ax_len[2]/(6*inertia[0]);
     free(eigval);
     free(eigvec);

   }
   free(inertia);
   return spars;

} /* inertia_spar_attribute */

double  inertia_ncom_attribute(void *inertiaattr){
   double *inertia  = inertia_arr(inertiaattr);
   double ncomp = 0;
   
   if(inertia[3] == 0){
     double corrVal    = inertia[0]/12;
     double matrix[2] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] + pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] + pixdim[1]*corrVal };    
     ncomp    =  (matrix[0]+matrix[1])/pow(inertia[0], 2.0);
   } else {
     double corrVal    = inertia[0]/12;
     double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] + pixdim[0]*corrVal,
       inertia[5] - inertia[2]*inertia[2]/inertia[0] + pixdim[1]*corrVal,
       inertia[6] - inertia[3]*inertia[3]/inertia[0] + pixdim[2]*corrVal,
       inertia[7] - inertia[1]*inertia[2]/inertia[0],
       inertia[8] - inertia[2]*inertia[3]/inertia[0],
       inertia[9] - inertia[1]*inertia[3]/inertia[0] };
      ncomp    =  (matrix[0]+matrix[1]+matrix[2])/pow(inertia[0], 5.0/3.0);
   }
   
   free(inertia);
   return ncomp;
   
} /* inertia_ncom_attribute */



/* Inertia moments weighted by gvalues  INCORRECT !!


void *new_inertiaweight_data( AuxDataStore *store, void *data){
  InertiaDataFull *inertiadata;
  inertiadata = store ? get_new_aux_data(store) : calloc(1, sizeof(InertiaDataFull));
  check_alloc(inertiadata, 114);
  return(inertiadata);
} // new_inertia_data 


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
} /// load_inertiaweight_data 

void delete_inertiaweight_data(void *inertiaattr)
{
  free(inertiaattr);
} // delete_inertiafull_data 

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
} // add_to_inertiaweight_data 

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
} // merge_inertia_data 


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
} // merge_to_inertiaweight_data

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
  
} // clone_inertiaweight_data 

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
} // inertia_attribute 

double inertiaweight_area_attribute(void *inertiaattr){
   InertiaDataFull *inertiadata = inertiaattr;
   double inertia = inertiadata->area;

   return(inertia);
} // inertia_attribute 



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
   return(elong);
} // mean_z_attribute 

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
} // mean_z_attribute 

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
} // mean_z_attribute 

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
} //mean_z_attribute 

double inertiafull_meangrey_attribute(void *inertattr){
   InertiaDataFull *inertiadata = inertattr;

  return inertiadata->sumval/inertiadata->area;

}
*/
/* PET Segmentation */



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
    init_aux_data((attribute + i*size_att), !init, x + attr_off[0], y + attr_off[1], z+ attr_off[2], gvals[i]);
  }
}


AttribStruct AttribsArray[NUMATTR] =
{
 {"Area", sizeof(AreaData), new_area_data, init_area_data, delete_area_data, add_to_area_data, merge_area_data, merge_to_area_data, clone_area_data, create_mpi_area_type, area_attribute, area_attribute}, //0
 
 {"Area of minimal enclosing rectangle", sizeof(ExtentData),  new_extent_data, init_extent_data, delete_extent_data, add_to_extent_data, merge_extent_data, merge_to_extent_data, clone_extent_data, create_mpi_extent_type, extent_rectarea_attribute, extent_rectarea_attribute}, //1
 {"Square of diagonal of minimal enclosing rectangle", sizeof(ExtentData), new_extent_data, init_extent_data, delete_extent_data, add_to_extent_data, merge_extent_data, merge_to_extent_data, clone_extent_data, create_mpi_extent_type, extent_rectarea_attribute, extent_rectdiag_attribute},//2
 
 {"X-extent", sizeof(ExtentData), new_extent_data, init_extent_data, delete_extent_data, add_to_extent_data, merge_extent_data, merge_to_extent_data, clone_extent_data, create_mpi_extent_type,extent_rectarea_attribute, extent_x_attribute},//3
 {"Y-extent", sizeof(ExtentData), new_extent_data, init_extent_data, delete_extent_data, add_to_extent_data, merge_extent_data, merge_to_extent_data, clone_extent_data, create_mpi_extent_type,extent_rectarea_attribute, extent_y_attribute},//4
 {"Z-extent", sizeof(ExtentData), new_extent_data, init_extent_data, delete_extent_data, add_to_extent_data, merge_extent_data, merge_to_extent_data, clone_extent_data, create_mpi_extent_type,extent_rectarea_attribute, extent_z_attribute},//5
 
  {"Mean X position", sizeof(MeanData), new_mean_data, init_mean_data, delete_mean_data, add_to_mean_data, merge_mean_data, merge_to_mean_data, clone_mean_data,  create_mpi_mean_type, mean_area_attribute, mean_x_attribute},//6
 {"Mean Y position", sizeof(MeanData), new_mean_data, init_mean_data, delete_mean_data, add_to_mean_data, merge_mean_data, merge_to_mean_data, clone_mean_data,  create_mpi_mean_type, mean_area_attribute, mean_y_attribute},//7
 {"Mean Z position", sizeof(MeanData), new_mean_data, init_mean_data, delete_mean_data, add_to_mean_data, merge_mean_data, merge_to_mean_data, clone_mean_data,  create_mpi_mean_type, mean_area_attribute, mean_z_attribute},//8

 
 {"Weighted Mean X position", sizeof(WMeanData), new_wmean_data, init_wmean_data, delete_wmean_data, add_to_wmean_data, merge_wmean_data, merge_to_wmean_data, clone_wmean_data,  create_mpi_wmean_type, wmean_totalflux_attribute, wmean_x_attribute},//9
 {"Weighted Mean Y position", sizeof(WMeanData), new_wmean_data, init_wmean_data, delete_wmean_data, add_to_wmean_data, merge_wmean_data, merge_to_wmean_data, clone_wmean_data,  create_mpi_wmean_type, wmean_totalflux_attribute, wmean_y_attribute},//10
 {"Weighted Mean Z position", sizeof(WMeanData), new_wmean_data, init_wmean_data, delete_wmean_data, add_to_wmean_data, merge_wmean_data, merge_to_wmean_data, clone_wmean_data,  create_mpi_wmean_type, wmean_totalflux_attribute, wmean_z_attribute},//11
  {"Total flux", sizeof(WMeanData), new_wmean_data, init_wmean_data, delete_wmean_data, add_to_wmean_data, merge_wmean_data, merge_to_wmean_data, clone_wmean_data,  create_mpi_wmean_type, wmean_totalflux_attribute, wmean_totalflux_attribute},//12
 
 {"Trace of Inertia Matrice", sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data, create_mpi_inertia_type, inertia_area_attribute, inertia_trace_attribute}, //13

 {"Elongation", sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data, create_mpi_inertia_type, inertia_area_attribute, inertia_elon_attribute},//14
 {"Flatness", sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data, create_mpi_inertia_type, inertia_area_attribute, inertia_flat_attribute},//15
 {"Sparseness", sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data, create_mpi_inertia_type, inertia_area_attribute, inertia_spar_attribute},//16
 {"Non-compactnesss",sizeof(InertiaData), new_inertia_data, init_inertia_data, delete_inertia_data, add_to_inertia_data, merge_inertia_data, merge_to_inertia_data, clone_inertia_data, create_mpi_inertia_type, inertia_area_attribute, inertia_ncom_attribute}, //17
};

