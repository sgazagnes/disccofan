#include "types.h"
#include "attributes.h"
#include "lambdavec.h"
#include "flood.h"
#include "csl.h"


int find_scale_csl(LambdaVec *lvec, double attribute){
  int upper = lvec->num_lambdas-1, lower = 0, mid;

  if (attribute >=  (double) lvec->lambdas[upper])
    return upper+1;


  mid = (upper + lower) / 2;
  while (mid!=lower) {
    if(attribute >= (double) lvec->lambdas[mid])
      lower = mid;
    else
      upper = mid;

    mid = (upper + lower) / 2;
  }
  return lower;
} /* find_scale */

void node_differential_seg(Node *tree, LambdaVec *lvec, value *contrast, value *scale, value *luminance, value *temp_contrast, value *temp_scale, bool *temp_valid, ulong lwb, ulong upb, ulong current, value *maxDH, value *curDH, value *maxOrig, int *maxScale, int *curScale, double (*attribute)(void *))
{

  int myscale = -1, numscales = lvec->num_lambdas;
  idx parent;
  value DH = 0;

  if (is_levelroot(tree, current))
  {

    myscale = find_scale_csl(lvec, (*attribute)(tree->attribute + current * tree->size_attr));
    DH = ((tree->parent[current] != BOTTOM) && (myscale < numscales)) ? tree->gval[current] - tree->gval[tree->parent[current]] : 0;
  }
  if (is_levelroot(tree, current) && ((myscale == numscales) || (tree->parent[current] == BOTTOM)))
  {

    if (tree->parent[current] != BOTTOM)
    {
      *maxScale = numscales;
      *maxOrig = 0;
    }
    else
    {
      *maxScale = myscale;
      *maxOrig = tree->gval[current];
    }

    *maxDH = 0;
    *curDH = 0;
    *curScale = *maxScale;
  }
  else
  {

    parent = tree->parent[current];

    if (!temp_valid[parent])
    {
      // go into recursion to set parent values correctly

      node_differential_seg(tree, lvec, contrast, scale, luminance, temp_contrast, temp_scale, temp_valid, lwb, upb, parent, maxDH, curDH, maxOrig, maxScale, curScale, attribute);
    }
    else
    { // if the parent is valid, copy relevant values

      *maxScale = (int) scale[parent];
      *maxDH = contrast[parent];
      *maxOrig = luminance[parent];
      *curScale = (int) temp_scale[parent];
      *curDH = temp_contrast[parent];
    }

    if (is_levelroot(tree, current))
    {
      // if I have a level root, some things might change 

      if (myscale == *curScale)
      {
        // parent's area is in same scale class,  add current pixel's curDH 
        *curDH += DH;
      }
      else
      {
        // at scale class change, update current scale and DH 
        *curDH = DH;
        *curScale = myscale;
      }

      if (*curDH >= *maxDH)
      {
        // If updated curDH is higher than or equal to the maximum DH found     update maxDH, maxScale, and outOrig 
        *maxDH = *curDH;
        *maxScale = *curScale;
        *maxOrig = tree->gval[current];
      }
    }
  }

  if ((current >= lwb) && (current < upb))
  {
    scale[current] = (value)*maxScale;
    contrast[current] = *maxDH;
    luminance[current] = *maxOrig;
    temp_scale[current] = (value)*curScale;
    temp_contrast[current] = *curDH;
    temp_valid[current] = true;
  }
  return;
} /* node_differential_seg */


void tree_differential(Node *tree, ulong size_tile, LambdaVec *lvec, value *contrast, value *scale, value *luminance, value *temp_contrast, value *temp_scale, bool *temp_valid, double (*attribute)(void *))
{
  ;
#pragma omp parallel
  {
    int np_threads = omp_get_num_threads();
    int id = omp_get_thread_num();
    ulong lwb = id * size_tile / np_threads;
    ulong upb = (id + 1) * size_tile / np_threads;
    ulong v;
    value curDH, maxDH, maxOrig;
    int maxScale, curScale;

    for (v = lwb; v < upb; v++)
    {
      if (!temp_valid[v])
      {
        node_differential_seg(tree, lvec, contrast, scale, luminance, temp_contrast, temp_scale, temp_valid, lwb, upb, v, &maxDH, &curDH, &maxOrig, &maxScale, &curScale, attribute);
      }
    }
  }
}

void combine_results(Node *tree, ulong size, LambdaVec *lvec, value *out_dh, value *out_dh2, value *out_orig, value *out_orig2, value *out_scale, value *out_scale2)
{

#pragma omp parallel
  {
    int np_threads = omp_get_num_threads();
    int id = omp_get_thread_num();
    ulong lwb = id * size / np_threads;
    ulong upb = (id + 1) * size / np_threads;
    int numscales = lvec->num_lambdas;
    for (ulong i = lwb; i < upb; i++)
    {

      if (out_dh2[i] > out_dh[i])
      { /* maxDH of opening > maxDH of closing instance */
        out_scale[i] = (value)numscales + 1 + out_scale2[i];
        out_dh[i] = out_dh2[i];
        out_orig[i] = out_orig2[i];
      }
      else if (out_dh2[i] < out_dh[i])
      {
        out_orig[i] = data_properties.g_max_gval - out_orig[i];
        out_scale[i]++;
      }
      else
      { /* equality */
        out_scale[i] = 0;
        out_orig[i] = (out_dh[i] != 0) * (data_properties.g_max_gval - tree->gval[i]);
      }
    }
  }
}

