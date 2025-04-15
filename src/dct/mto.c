#include "types.h"
#include "mto.h"
#include "tree_flood.h"

static const idx max_area = 4087;
//static const INT_TYPE max_area = 2;

static const float p1 = 1.683355084690155e-01;
static const float p2 = 3.770229379757511e+02;
static const float p3 = 1.176722049258011e+05;
static const float p4 = 6.239836661965291e+06;
static const float q1 = 1.354265276841128e+03;
static const float q2 = 2.091126298053044e+05;
static const float q3 = 1.424803575269314e+06;


void mt_objects_init(Node *mt, mt_object_data* mt_o)
{
  memset(mt_o, 0, sizeof(mt_object_data));
  
  mt_o->mt = mt;
  mt_o->flags = calloc(mt->size, sizeof(*mt_o->flags));
  
  mt_o->closest_significant_ancestors = malloc(mt->size *
    sizeof(*mt_o->closest_significant_ancestors));
    
  for (idx i = 0; i != mt->size; ++i)
  {
    mt_o->closest_significant_ancestors[i] = -3;
  }

  mt_o->main_branches = malloc(mt->size *
    sizeof(*mt_o->main_branches));
    
  mt_o->main_power_branches = malloc(mt->size *
    sizeof(*mt_o->main_power_branches));

  mt_o->bg_variance = 4;

  mt_o->gain = 5.5;
  mt_o->move_factor = 0.5;
  info("%f",mt_o->gain);

  info("DONE INIT MTO");
}



float mt_noise_variance(mt_object_data* mt_o,
  idx node_idx, float max_normalized_distance)
{
  Node* mt = mt_o->mt;

  assert(mt_o->relevant_indices_len > 0); 
  
  float variance = mt_o->bg_variance;
  if (mt_o->closest_significant_ancestors[node_idx] != MT_NO_PARENT)
  {
    variance +=
      mt->gval[mt_o->closest_significant_ancestors[node_idx]] /
      mt_o->gain;
  }
  
  if (max_normalized_distance >= 0)
  {
    float distance = (mt_o->closest_significant_ancestors[node_idx] != MT_NO_PARENT ? 
		      mt->gval[node_idx] - mt->gval[mt_o->closest_significant_ancestors[node_idx]] : 
		      mt->gval[node_idx]) ;
    
    float rms = sqrt(variance);
    
    if (distance / rms > max_normalized_distance)
    {
      float max_normalized_distance_sqr = max_normalized_distance *
        max_normalized_distance;
        
      float gain_sqr = mt_o->gain * mt_o->gain;
      
      float b = 2 * mt->gval[node_idx] * mt_o->gain;
        
      float f_a = b + max_normalized_distance_sqr - max_normalized_distance *
          sqrt(4 * mt_o->bg_variance * gain_sqr + 2 * b + max_normalized_distance_sqr);
          
      f_a /= 2 * mt_o->gain;
      
      variance = f_a / mt_o->gain + mt_o->bg_variance;
    }
  }
  
  return variance;
}

double mt_alternative_power_definition(mt_object_data* mt_o,
  idx node_idx, float max_normalized_distance)
{
  Node* mt = mt_o->mt;
  
  assert(mt_o->relevant_indices_len > 0); 
  

  double *attr = (*AttribsArray[9].attribute_arr)(mt->attribute[node_idx]); 
  
  idx parent_idx = mt->parent[node_idx];  
  
  // added distance for the power calculation.
  double delta;
  if ((mt_o->closest_significant_ancestors[node_idx] != MT_NO_PARENT))
  {
    delta = mt->gval[parent_idx] -
      mt->gval[mt_o->closest_significant_ancestors[node_idx]];      
  }
  else
  {
    delta = mt->gval[parent_idx];    
  }  
      
  if (max_normalized_distance >= 0)
  {
    double rms = sqrt(mt_noise_variance(mt_o,
      node_idx, max_normalized_distance));
    
    double distance = (mt_o->closest_significant_ancestors[node_idx] != MT_NO_PARENT ?  mt->gval[node_idx] -  mt->gval[mt_o->closest_significant_ancestors[node_idx]] : mt->gval[node_idx]);

    printf("%f\n", distance / rms);

    if (distance / rms > max_normalized_distance)
    {
      delta = max_normalized_distance * rms - mt->gval[parent_idx];
    }  
  }
    
  return attr[1] + delta * (2 * attr[2] + delta * attr[0]);
}

int mt_node_test_4(mt_object_data* mt_o, idx node_idx)
{
  // Point to max tree data
  Node* mt = mt_o->mt;

  float variance = mt_noise_variance(mt_o, node_idx, MT_NO_MAX_DISTANCE);
    
  float min_distance =
    *((float *)mt_o->node_significance_test_data);
  info("%f", min_distance);

  if (min_distance > 0 && ((mt_o->closest_significant_ancestors[node_idx] != MT_NO_PARENT ?  mt->gval[node_idx] -  mt->gval[mt_o->closest_significant_ancestors[node_idx]] : mt->gval[node_idx])
			   / sqrt(variance) < min_distance))
  {
    return 0;
  }
  
  float power = mt_alternative_power_definition(mt_o, node_idx, MT_NO_MAX_DISTANCE);
    
  ulong area = (*AttribsArray[9].attribute)(mt->attribute[node_idx]);
  float power_normalized = power / variance / area;
  //info("%f", power_normalized);
  if (area > max_area)
  {
    area = max_area;
  }  

  float area_to_2 = area * area;
  float area_to_3 = area_to_2 * area;
  
  float x = p1 * area_to_3 + p2 * area_to_2 + p3 * area + p4;
  x /= area_to_3 + q1 * area_to_2 + q2 * area + q3;

  return power_normalized > x;
}

void mt_node_test_4_data_free(mt_object_data* mt_o)
{
  free(mt_o->node_significance_test_data);
  
  mt_o->node_significance_test_data = NULL;
}


void node_significance_test_data_clear(mt_object_data* mt_o)
{
  if (mt_o->node_significance_test_data != NULL)
  {
    mt_o->node_significance_test_data_free(mt_o);
  }
  
  mt_o->node_significance_test_data_free = NULL;
  mt_o->node_significance_test_data = NULL;  
}


void mt_use_node_test_4(mt_object_data* mt_o,float significance_level_power,
  float min_distance)
{

  info("GHERE");
  node_significance_test_data_clear(mt_o);

  mt_o->node_significance_test_data = malloc(sizeof(float));

  *((float *)mt_o->node_significance_test_data) =
    min_distance;

  mt_o->node_significance_test = mt_node_test_4;

  mt_o->node_significance_test_data_free =  mt_node_test_4_data_free;

}


static void mt_relevant_nodes(mt_object_data* mt_o)
{
  Node *mt = mt_o->mt;

  idx *list_idx = malloc(mt->size * sizeof(idx));
  value *gval_arr = malloc(mt->size * sizeof(value));
  //PrioQueue *q    = create_prio_queue(mt->size);	
  ulong count = 0;
  info("%ld", mt->size);

  for (ulong i = 0; i < mt->size; i++){
    if(  is_levelroot(mt, i) && mt->parent[i] != BOTTOM){
      if(mt->gval[i] == mt->gval[mt->parent[i]])
	 error("EHER");
      list_idx[count]=i;
      gval_arr[count]=mt->gval[i];
      count++;
    }
  }
  //
  mt_o->relevant_indices_len = count;
  mt_o->relevant_indices = malloc(count *     sizeof(*mt_o->relevant_indices));
 
  info("Number of nodes to be tested: %d.\n",  mt_o->relevant_indices_len);

  ulong *ranks = calloc(count, sizeof(ulong));
  create_mappings(gval_arr, NULL, ranks, count, 0, count);

  for (int i = mt_o->relevant_indices_len; i--;)
  {
    mt_o->relevant_indices[i] = list_idx[ranks[i]];
    // info("%f", mt->gval[ list_idx[ranks[i]]]);
  }

  free(list_idx);
  free(gval_arr);
  free(ranks);
}


static void mt_update_parent_main_branch(
  mt_object_data* mt_o, idx node_idx)
{
  Node *mt = mt_o->mt;
  
  if (mt_o->closest_significant_ancestors[node_idx] == MT_NO_PARENT)
    return;
    
  idx ancestor_idx = mt_o-> closest_significant_ancestors[node_idx];
    
  if (mt_o->flags[ancestor_idx] & 4)
  {
    double *attr = (*AttribsArray[9].attribute_arr)(mt->attribute[node_idx]); 
    double *attranc = (*AttribsArray[9].attribute_arr)(mt->attribute[ancestor_idx]); 

    if (attranc[0] < attr[0])
    {
      mt_o->main_branches[ancestor_idx] = node_idx;
    }
  }
  else
  {
    mt_o->flags[ancestor_idx] |= 4;
    mt_o->main_branches[ancestor_idx] = node_idx;
  }
}

static void mt_significant_nodes(mt_object_data* mt_o)
{
  Node *mt = mt_o->mt;
  
  idx num_significant = 0;
  //idx i;
  for (idx i = (idx) mt_o->relevant_indices_len; i--; )
  {

    idx node_idx = mt_o->relevant_indices[i];
    idx parent_idx = mt->parent[node_idx];
    
    if (mt_o->flags[parent_idx] & 1)
    {
      mt_o->closest_significant_ancestors[node_idx] = parent_idx;
    }
    else if (mt_o->closest_significant_ancestors[parent_idx] != MT_NO_PARENT)
    {
      mt_o->closest_significant_ancestors[node_idx] =
        mt_o->closest_significant_ancestors[parent_idx];
    }
    if (mt_o->node_significance_test(mt_o, node_idx))
    {
      
      mt_o->flags[node_idx] |= 1;
      ++num_significant;
      
      mt_update_parent_main_branch(mt_o, node_idx);
    }
  }
  
 info("%d significant nodes.\n", num_significant);
  
  
  mt_o->num_significant_nodes = num_significant;
}

void mt_find_objects(mt_object_data* mt_o)
{
  // Count significant nodes and set object markers

  Node *mt = mt_o->mt;

  ulong num_objects = 0;
  ulong num_objects_nested = 0;
  
  ulong i;

  // Iterate over pixels in image
  for (i = 0; i != mt->size; ++i)  
  {
    // Skip if no significant flag
    if (!(mt_o->flags[i] & 1))
    {
      continue;
    }

    // Count and mark as object if no significant ancestor
    if (!(mt_o->closest_significant_ancestors[i] != MT_NO_PARENT))
    {
      ++num_objects;
      mt_o->flags[i] |= 8;

      continue;
    }

    // If it has a significant ancestor, and that ancestor's largest descendant is NOT this node,
        // mark it as a nested object
    idx parent = mt_o->closest_significant_ancestors[i];

    if (mt_o->main_branches[parent] != i)
    {
      ++num_objects_nested;
      mt_o->flags[i] |= 8;
      continue;
    }

    // i.e. significant nodes who ARE the significant descendant of their significant ancestor
    // are not marked as objects
  }

  // Get the total number of objects and print
  num_objects += num_objects_nested;


  info("Found %d objects (including %d nested).\n", num_objects,   num_objects_nested);
    
  
  mt_o->num_objects = num_objects;
}


void mt_main_power_branches(mt_object_data* mt_o)
{
  // Find the descendant of each node with the highest power

  Node *mt = mt_o->mt;
  
  idx i;
  // Iterate over image pixels
  for (i = 0; i != mt->size; ++i)  
  {
    // Skip the root
    if (mt->parent[i] != BOTTOM)
    {
      continue;
    }

    // Get the pixel's parent
    idx parent = mt->parent[i];

    // If the parent has a descendant, check if this node has a higher power and set accordingly
    if (mt_o->flags[parent] & 16)
    {
      double *attr = (*AttribsArray[9].attribute_arr)(mt->attribute[parent]); 
      double *attri = (*AttribsArray[9].attribute_arr)(mt->attribute[i]); 

      
      if (attr[1] < attri[1])
      {
        mt_o->main_power_branches[parent] = i;
      }
    }
    else
    // If the parent has no marked descendant, set this as the highest power descendant
    {
      mt_o->flags[parent] |= 16;
      mt_o->main_power_branches[parent] = i;
    }
  }
}

void mt_move_up(mt_object_data* mt_o, float move_factor, float gain, float bg_variance)
{
  // Move object markers up the tree

  Node *mt = mt_o->mt;

  // Skip this function if no move factor is specified
  /* if (mt_o->paras->move_factor == 0.0)
  {
    return;
    }*/

  // Iterate over image pixels
  idx i;  
  for (i = 0; i != mt->size; ++i)  
  {  
	// Skip if the node is not an object or is marked as 'don't move'
    if (!(mt_o->flags[i] & 8) || (mt_o->flags[i] & 32))
    {
      continue;
    }
    
    // Mark as not an object
    mt_o->flags[i] &= ~8;
    
    // Base = pixel value - distance from nearest significant ancestor
    float base = mt->gval[i] - (mt_o->closest_significant_ancestors[i] != MT_NO_PARENT ?  mt->gval[i] -  mt->gval[mt_o->closest_significant_ancestors[i]] : mt->gval[i]);
      
    // scale by gain, background variance, move_factor
    // seems a tad hacky
    base +=move_factor *
      sqrt(base / gain + bg_variance);
            
    idx next_idx = i;
    // Find next id - first pixel with value less than base or no descendants
    // Check pixels by descending through (significant) descendants
    while (mt->gval[next_idx] < base)
    {      
    
      if (mt_o->flags[next_idx] & 4)
      {
        next_idx = mt_o->main_branches[next_idx];        
      }
      else if (mt_o->flags[next_idx] & 16)
      {
        next_idx = mt_o->main_power_branches[next_idx];
      }
      else
      {
        break;
      }
    }
    
    // Mark the next ID as an object that cannot be moved
    mt_o->flags[next_idx] |= 8;
    mt_o->flags[next_idx] |= 32;
  }
}

void mt_objects(mt_object_data* mt_o)
{

  assert(mt_o->bg_variance > 0);
  assert(mt_o->gain > 0);
  assert(mt_o->move_factor >= 0);
  assert(mt_o->node_significance_test != NULL);
  
  mt_relevant_nodes(mt_o);
  
  mt_significant_nodes(mt_o);

  mt_find_objects(mt_o);
    
  mt_main_power_branches(mt_o);
    
  mt_move_up(mt_o,mt_o->move_factor,mt_o->gain ,mt_o->bg_variance );    
}
