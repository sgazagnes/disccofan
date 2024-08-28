#include "types.h"
#include "attributes.h"
#include "lambdavec.h"
#include "flood.h"
#include "pattern.h"

int find_scale(LambdaVec *lvec, double attribute){
  int upper = lvec->num_lambdas-1, lower = 0, mid;

  if (attribute >=  (double) lvec->lambdas[upper])
    return upper;


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


void tree_pattern_spectrum(Node *tree, ulong size, LambdaVec *lvec, double *copy_area, value *gvals_par, double *spectrum, bool background, double (*area)(void *), double (*attribute)(void *))
{

  int numscales = lvec->num_lambdas;
  double *spectrum_th;

#pragma omp parallel
  {
    int np_threads = omp_get_num_threads();
    int id = omp_get_thread_num();
    int scale;
    ulong u, v;
    idx parent;
    value dh;
    double private;

#pragma omp single
    spectrum_th = calloc(numscales * np_threads, sizeof(double));
    check_alloc(spectrum_th, 601);
    if (np() > 1)
    {
#pragma omp for
      for (v = 0; v < size; v++)
      {
        if ((copy_area[v] != -10) && (tree->parent[v] != BOTTOM || background))
        {
          scale = find_scale(lvec, (*attribute)(tree->attribute + get_levelroot(tree, v) * tree->size_attr));
          if (scale < numscales)
          {
            parent = get_levelroot(tree, tree->parent[v]);
            private = copy_area[v];
            dh = tree->parent[v] != BOTTOM ? (tree->gval[v] - tree->gval[parent]) : tree->gval[v];
            spectrum_th[(numscales)*id + scale] += dh * private;
            u = v;
            if (tree->parent[v] != BOTTOM)
            {
              while ((tree->parent[parent] != BOTTOM || background) && (gvals_par[v] < tree->gval[parent]))
              {
                u = parent;
                scale = find_scale(lvec, (*attribute)(tree->attribute + u * tree->size_attr));
                if (scale >= numscales)
                  break;
                parent = get_levelroot(tree, tree->parent[u]);
                dh = parent != BOTTOM ? (tree->gval[u] - tree->gval[parent]) : tree->gval[u];
                spectrum_th[(numscales)*id + scale] += dh * private;
                if (parent == BOTTOM)
                  break;
              }
            }
          }
        }
      }
    }
    else
    {
#pragma omp for
      for (v = 0; v < size; v++)
      {
        if (is_levelroot(tree, v) && (tree->parent[v] != BOTTOM || background))
        {
          dh = tree->parent[v] != BOTTOM ? (tree->gval[v] - tree->gval[tree->parent[v]]) : tree->gval[v];
          if (dh)
          {
            scale = find_scale(lvec, (*attribute)(tree->attribute + v * tree->size_attr));
            if (scale < numscales)
              spectrum_th[(numscales)*id + scale] += dh * (*area)(tree->attribute + v * tree->size_attr);
          }
        }
      }
    }

#pragma omp for
    for (int i = 0; i < numscales; i++)
    {
      for (int t = 0; t < np_threads; t++)
      {
        spectrum[i] += spectrum_th[t * numscales + i];
      }
    }
  }
  free(spectrum_th);
} /* tree_pattern_spectrum */

void tree_pattern_spectrum2d(Node *tree,  ulong size, LambdaVec *lvec_attr1,LambdaVec *lvec_attr2, double *copy_attr, value *gvals_par, double* spectrum, int background,  double (*area)(void *), double (*attribute1)(void *), double (*attribute2)(void *)){

  int numscales_attr1 = lvec_attr1->num_lambdas;
  int numscales_attr2 = lvec_attr2->num_lambdas;

  double *spectrum_th;

  #pragma omp parallel
  {
    int np_threads = omp_get_num_threads();
    int   id	   = omp_get_thread_num();
    int scale_attr1, scale_attr2;
    ulong u, v;
    idx parent;
    value dh;
    double private;

#pragma omp single
    spectrum_th = calloc(numscales_attr1 * numscales_attr2 * np_threads, sizeof(double));
    check_alloc(spectrum_th, 601);

    if (np() > 1)
    {

#pragma omp for
      for (v = 0; v < size; v++)
      {
        if ((copy_attr[v] != -10) && (tree->parent[v] != BOTTOM || background))
        {
          // scale = find_scale(lvec, (*attribute)(tree->attribute + get_levelroot(tree, v) * tree->size_attr));
          scale_attr1 = find_scale(lvec_attr1, (*attribute1)(tree->attribute + get_levelroot(tree, v) * tree->size_attr));
          scale_attr2 = find_scale(lvec_attr2, (*attribute2)(tree->attribute + get_levelroot(tree, v) * tree->size_attr));
          if (scale_attr1 < numscales_attr1 && scale_attr2 < numscales_attr2)
          {
            parent = get_levelroot(tree, tree->parent[v]);
            private = copy_attr[v];
            dh = tree->parent[v] != BOTTOM ? (tree->gval[v] - tree->gval[parent]) : tree->gval[v];
            spectrum_th[(numscales_attr1 * numscales_attr2) * id + scale_attr1 + scale_attr2 * numscales_attr1] += dh * private;
            u = v;
            if (tree->parent[v] != BOTTOM)
            {
              while ((tree->parent[parent] != BOTTOM || background) && (gvals_par[v] < tree->gval[parent]))
              {
                u = parent;
                scale_attr1 = find_scale(lvec_attr1, (*attribute1)(tree->attribute + u* tree->size_attr));
                scale_attr2 = find_scale(lvec_attr2, (*attribute2)(tree->attribute + u* tree->size_attr));
                // scale = find_scale(lvec, (*attribute)(tree->attribute + u * tree->size_attr));
                if (scale_attr1 >= numscales_attr1 && scale_attr2 >= numscales_attr2)
                  break;
                parent = get_levelroot(tree, tree->parent[u]);
                dh = parent != BOTTOM ? (tree->gval[u] - tree->gval[parent]) : tree->gval[u];
                spectrum_th[(numscales_attr1 * numscales_attr2) * id + scale_attr1 + scale_attr2 * numscales_attr1] += dh * private;
                if (parent == BOTTOM)
                  break;
              }
            }
          }
        }
      }
    }
    else
    {
      // error("%d, %d", numscales_attr1,numscales_attr2);
#pragma omp for
      for (v = 0; v < size; v++)
      {
        if (is_levelroot(tree, v) && (tree->parent[v] != BOTTOM || background))
        {
          dh = tree->parent[v] != BOTTOM ? (tree->gval[v] - tree->gval[tree->parent[v]]) : tree->gval[v];
          if (dh)
          {
            scale_attr1 = find_scale(lvec_attr1, (*attribute1)(tree->attribute + v * tree->size_attr));
            scale_attr2 = find_scale(lvec_attr2, (*attribute2)(tree->attribute + v * tree->size_attr));

            if (scale_attr1 < numscales_attr1 || scale_attr2 < numscales_attr2){
              spectrum_th[(numscales_attr1 * numscales_attr2) * id + scale_attr1 + scale_attr2 * numscales_attr1] +=
                  dh * (*area)(tree->attribute + v * tree->size_attr);
            }
          }
        }
      }
    }

#pragma omp for
    for (int i = 0; i < numscales_attr1 * numscales_attr2; i++)
    {
      for (int t = 0; t < np_threads; t++)
      {
        spectrum[i] += spectrum_th[t * numscales_attr1 * numscales_attr2 + i];
      }
    }
  }
  free(spectrum_th);
} /* tree_pattern_spectrum_2D */
