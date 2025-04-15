#include "types.h"

#ifndef MT_OBJECTS_H
#define MT_OBJECTS_H

#define MT_UNASSIGNED -1
#define MT_IN_QUEUE -2
#define MT_NO_PARENT -3

#define MT_NO_OBJECT -1

#define MT_IS_ROOT(MT_PTR, IDX) ((MT_PTR)->nodes + IDX == \
  (MT_PTR)->root)

#define MT_CONN_12_WIDTH 5
#define MT_CONN_12_HEIGHT 5

#define MT_CONN_8_WIDTH 3
#define MT_CONN_8_HEIGHT 3

#define MT_CONN_4_WIDTH 3
#define MT_CONN_4_HEIGHT 3

extern const int mt_conn_12[MT_CONN_12_HEIGHT * MT_CONN_12_WIDTH];
extern const int mt_conn_8[MT_CONN_8_HEIGHT * MT_CONN_8_WIDTH];
extern const int mt_conn_4[MT_CONN_4_HEIGHT * MT_CONN_4_WIDTH];


struct mt_object_data;

void mt_objects_init(Node* mt, mt_object_data* mt_o);
void mt_objects_free(mt_object_data* mt_o);
void mt_objects(mt_object_data* mt_o);
void mt_object_ids(mt_object_data* mt_o);
  
void mt_use_node_test_1(mt_object_data* mt_o,
  float significance_level_power);
void mt_use_node_test_2(mt_object_data* mt_o,
  float significance_level_power);  
void mt_use_node_test_3(mt_object_data* mt_o,
  float significance_level_power);
void mt_use_node_test_4(
  mt_object_data* mt_o, float significance_level_power,
  float min_distance);
void mt_use_node_test_5(
  mt_object_data* mt_o, float significance_level_power,
  float min_distance);

void mt_set_bg_variance(mt_object_data* mt_o, float bg_variance);
void mt_set_gain(mt_object_data* mt_o, float gain);
void mt_set_move_factor(mt_object_data* mt_o, float move_factor);

void node_significance_test_data_clear(mt_object_data* mt_o);

#define MT_NO_MAX_DISTANCE -1.0

float mt_noise_variance(mt_object_data* mt_o,
  idx node_idx, float max_normalized_distance);
double mt_alternative_power_definition(mt_object_data* mt_o,
  idx node_idx, float max_normalized_distance);

#endif
