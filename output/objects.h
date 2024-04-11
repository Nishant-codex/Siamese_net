
#ifndef _BRIAN_OBJECTS_H
#define _BRIAN_OBJECTS_H

#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include "randomkit.h"
#include<vector>


namespace brian {

// In OpenMP we need one state per thread
extern std::vector< rk_state* > _mersenne_twister_states;

//////////////// clocks ///////////////////
extern Clock defaultclock;

//////////////// networks /////////////////
extern Network magicnetwork;

//////////////// dynamic arrays ///////////
extern std::vector<int32_t> _dynamic_array_spikemonitor_i;
extern std::vector<int32_t> _dynamic_array_spikemonitor_i_1;
extern std::vector<double> _dynamic_array_spikemonitor_t;
extern std::vector<double> _dynamic_array_spikemonitor_t_1;
extern std::vector<double> _dynamic_array_statemonitor_t;
extern std::vector<double> _dynamic_array_statemonitor_t_1;

//////////////// arrays ///////////////////
extern double *_array_defaultclock_dt;
extern const int _num__array_defaultclock_dt;
extern double *_array_defaultclock_t;
extern const int _num__array_defaultclock_t;
extern int64_t *_array_defaultclock_timestep;
extern const int _num__array_defaultclock_timestep;
extern int32_t *_array_neurongroup__spikespace;
extern const int _num__array_neurongroup__spikespace;
extern int32_t *_array_neurongroup__spikespace_1;
extern const int _num__array_neurongroup__spikespace_1;
extern int32_t *_array_neurongroup_i;
extern const int _num__array_neurongroup_i;
extern double *_array_neurongroup_I_0;
extern const int _num__array_neurongroup_I_0;
extern double *_array_neurongroup_I_0_1;
extern const int _num__array_neurongroup_I_0_1;
extern double *_array_neurongroup_I_1;
extern const int _num__array_neurongroup_I_1;
extern int32_t *_array_neurongroup_i_1;
extern const int _num__array_neurongroup_i_1;
extern double *_array_neurongroup_I_1_1;
extern const int _num__array_neurongroup_I_1_1;
extern double *_array_neurongroup_lastspike;
extern const int _num__array_neurongroup_lastspike;
extern double *_array_neurongroup_lastspike_1;
extern const int _num__array_neurongroup_lastspike_1;
extern char *_array_neurongroup_not_refractory;
extern const int _num__array_neurongroup_not_refractory;
extern char *_array_neurongroup_not_refractory_1;
extern const int _num__array_neurongroup_not_refractory_1;
extern double *_array_neurongroup_Th_s;
extern const int _num__array_neurongroup_Th_s;
extern double *_array_neurongroup_Th_s_1;
extern const int _num__array_neurongroup_Th_s_1;
extern double *_array_neurongroup_Th_v;
extern const int _num__array_neurongroup_Th_v;
extern double *_array_neurongroup_Th_v_1;
extern const int _num__array_neurongroup_Th_v_1;
extern double *_array_neurongroup_V;
extern const int _num__array_neurongroup_V;
extern double *_array_neurongroup_V_1;
extern const int _num__array_neurongroup_V_1;
extern int32_t *_array_spikemonitor__source_idx;
extern const int _num__array_spikemonitor__source_idx;
extern int32_t *_array_spikemonitor__source_idx_1;
extern const int _num__array_spikemonitor__source_idx_1;
extern int32_t *_array_spikemonitor_count;
extern const int _num__array_spikemonitor_count;
extern int32_t *_array_spikemonitor_count_1;
extern const int _num__array_spikemonitor_count_1;
extern int32_t *_array_spikemonitor_N;
extern const int _num__array_spikemonitor_N;
extern int32_t *_array_spikemonitor_N_1;
extern const int _num__array_spikemonitor_N_1;
extern int32_t *_array_statemonitor__indices;
extern const int _num__array_statemonitor__indices;
extern int32_t *_array_statemonitor__indices_1;
extern const int _num__array_statemonitor__indices_1;
extern double *_array_statemonitor_I_0;
extern const int _num__array_statemonitor_I_0;
extern double *_array_statemonitor_I_0_1;
extern const int _num__array_statemonitor_I_0_1;
extern double *_array_statemonitor_I_1;
extern const int _num__array_statemonitor_I_1;
extern double *_array_statemonitor_I_1_1;
extern const int _num__array_statemonitor_I_1_1;
extern int32_t *_array_statemonitor_N;
extern const int _num__array_statemonitor_N;
extern int32_t *_array_statemonitor_N_1;
extern const int _num__array_statemonitor_N_1;
extern double *_array_statemonitor_Th_s;
extern const int _num__array_statemonitor_Th_s;
extern double *_array_statemonitor_Th_s_1;
extern const int _num__array_statemonitor_Th_s_1;
extern double *_array_statemonitor_Th_v;
extern const int _num__array_statemonitor_Th_v;
extern double *_array_statemonitor_Th_v_1;
extern const int _num__array_statemonitor_Th_v_1;
extern double *_array_statemonitor_V;
extern const int _num__array_statemonitor_V;
extern double *_array_statemonitor_V_1;
extern const int _num__array_statemonitor_V_1;

//////////////// dynamic arrays 2d /////////
extern DynamicArray2D<double> _dynamic_array_statemonitor_I_0;
extern DynamicArray2D<double> _dynamic_array_statemonitor_I_0_1;
extern DynamicArray2D<double> _dynamic_array_statemonitor_I_1;
extern DynamicArray2D<double> _dynamic_array_statemonitor_I_1_1;
extern DynamicArray2D<double> _dynamic_array_statemonitor_Th_s;
extern DynamicArray2D<double> _dynamic_array_statemonitor_Th_s_1;
extern DynamicArray2D<double> _dynamic_array_statemonitor_Th_v;
extern DynamicArray2D<double> _dynamic_array_statemonitor_Th_v_1;
extern DynamicArray2D<double> _dynamic_array_statemonitor_V;
extern DynamicArray2D<double> _dynamic_array_statemonitor_V_1;

/////////////// static arrays /////////////
extern double *_timedarray_1_values;
extern const int _num__timedarray_1_values;
extern double *_timedarray_values;
extern const int _num__timedarray_values;

//////////////// synapses /////////////////

// Profiling information for each code object
}

void _init_arrays();
void _load_arrays();
void _write_arrays();
void _dealloc_arrays();

#endif


