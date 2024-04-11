
#include "objects.h"
#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include "randomkit.h"
#include<vector>
#include<iostream>
#include<fstream>

namespace brian {

std::vector< rk_state* > _mersenne_twister_states;

//////////////// networks /////////////////
Network magicnetwork;

//////////////// arrays ///////////////////
double * _array_defaultclock_dt;
const int _num__array_defaultclock_dt = 1;
double * _array_defaultclock_t;
const int _num__array_defaultclock_t = 1;
int64_t * _array_defaultclock_timestep;
const int _num__array_defaultclock_timestep = 1;
int32_t * _array_neurongroup__spikespace;
const int _num__array_neurongroup__spikespace = 2;
int32_t * _array_neurongroup__spikespace_1;
const int _num__array_neurongroup__spikespace_1 = 2;
int32_t * _array_neurongroup_i;
const int _num__array_neurongroup_i = 1;
double * _array_neurongroup_I_0;
const int _num__array_neurongroup_I_0 = 1;
double * _array_neurongroup_I_0_1;
const int _num__array_neurongroup_I_0_1 = 1;
double * _array_neurongroup_I_1;
const int _num__array_neurongroup_I_1 = 1;
int32_t * _array_neurongroup_i_1;
const int _num__array_neurongroup_i_1 = 1;
double * _array_neurongroup_I_1_1;
const int _num__array_neurongroup_I_1_1 = 1;
double * _array_neurongroup_lastspike;
const int _num__array_neurongroup_lastspike = 1;
double * _array_neurongroup_lastspike_1;
const int _num__array_neurongroup_lastspike_1 = 1;
char * _array_neurongroup_not_refractory;
const int _num__array_neurongroup_not_refractory = 1;
char * _array_neurongroup_not_refractory_1;
const int _num__array_neurongroup_not_refractory_1 = 1;
double * _array_neurongroup_Th_s;
const int _num__array_neurongroup_Th_s = 1;
double * _array_neurongroup_Th_s_1;
const int _num__array_neurongroup_Th_s_1 = 1;
double * _array_neurongroup_Th_v;
const int _num__array_neurongroup_Th_v = 1;
double * _array_neurongroup_Th_v_1;
const int _num__array_neurongroup_Th_v_1 = 1;
double * _array_neurongroup_V;
const int _num__array_neurongroup_V = 1;
double * _array_neurongroup_V_1;
const int _num__array_neurongroup_V_1 = 1;
int32_t * _array_spikemonitor__source_idx;
const int _num__array_spikemonitor__source_idx = 1;
int32_t * _array_spikemonitor__source_idx_1;
const int _num__array_spikemonitor__source_idx_1 = 1;
int32_t * _array_spikemonitor_count;
const int _num__array_spikemonitor_count = 1;
int32_t * _array_spikemonitor_count_1;
const int _num__array_spikemonitor_count_1 = 1;
int32_t * _array_spikemonitor_N;
const int _num__array_spikemonitor_N = 1;
int32_t * _array_spikemonitor_N_1;
const int _num__array_spikemonitor_N_1 = 1;
int32_t * _array_statemonitor__indices;
const int _num__array_statemonitor__indices = 1;
int32_t * _array_statemonitor__indices_1;
const int _num__array_statemonitor__indices_1 = 1;
double * _array_statemonitor_I_0;
const int _num__array_statemonitor_I_0 = (0, 1);
double * _array_statemonitor_I_0_1;
const int _num__array_statemonitor_I_0_1 = (0, 1);
double * _array_statemonitor_I_1;
const int _num__array_statemonitor_I_1 = (0, 1);
double * _array_statemonitor_I_1_1;
const int _num__array_statemonitor_I_1_1 = (0, 1);
int32_t * _array_statemonitor_N;
const int _num__array_statemonitor_N = 1;
int32_t * _array_statemonitor_N_1;
const int _num__array_statemonitor_N_1 = 1;
double * _array_statemonitor_Th_s;
const int _num__array_statemonitor_Th_s = (0, 1);
double * _array_statemonitor_Th_s_1;
const int _num__array_statemonitor_Th_s_1 = (0, 1);
double * _array_statemonitor_Th_v;
const int _num__array_statemonitor_Th_v = (0, 1);
double * _array_statemonitor_Th_v_1;
const int _num__array_statemonitor_Th_v_1 = (0, 1);
double * _array_statemonitor_V;
const int _num__array_statemonitor_V = (0, 1);
double * _array_statemonitor_V_1;
const int _num__array_statemonitor_V_1 = (0, 1);

//////////////// dynamic arrays 1d /////////
std::vector<int32_t> _dynamic_array_spikemonitor_i;
std::vector<int32_t> _dynamic_array_spikemonitor_i_1;
std::vector<double> _dynamic_array_spikemonitor_t;
std::vector<double> _dynamic_array_spikemonitor_t_1;
std::vector<double> _dynamic_array_statemonitor_t;
std::vector<double> _dynamic_array_statemonitor_t_1;

//////////////// dynamic arrays 2d /////////
DynamicArray2D<double> _dynamic_array_statemonitor_I_0;
DynamicArray2D<double> _dynamic_array_statemonitor_I_0_1;
DynamicArray2D<double> _dynamic_array_statemonitor_I_1;
DynamicArray2D<double> _dynamic_array_statemonitor_I_1_1;
DynamicArray2D<double> _dynamic_array_statemonitor_Th_s;
DynamicArray2D<double> _dynamic_array_statemonitor_Th_s_1;
DynamicArray2D<double> _dynamic_array_statemonitor_Th_v;
DynamicArray2D<double> _dynamic_array_statemonitor_Th_v_1;
DynamicArray2D<double> _dynamic_array_statemonitor_V;
DynamicArray2D<double> _dynamic_array_statemonitor_V_1;

/////////////// static arrays /////////////
double * _timedarray_1_values;
const int _num__timedarray_1_values = 20000;
double * _timedarray_values;
const int _num__timedarray_values = 20000;

//////////////// synapses /////////////////

//////////////// clocks ///////////////////
Clock defaultclock;  // attributes will be set in run.cpp

// Profiling information for each code object
}

void _init_arrays()
{
	using namespace brian;

    // Arrays initialized to 0
	_array_defaultclock_dt = new double[1];
    
	for(int i=0; i<1; i++) _array_defaultclock_dt[i] = 0;

	_array_defaultclock_t = new double[1];
    
	for(int i=0; i<1; i++) _array_defaultclock_t[i] = 0;

	_array_defaultclock_timestep = new int64_t[1];
    
	for(int i=0; i<1; i++) _array_defaultclock_timestep[i] = 0;

	_array_neurongroup__spikespace = new int32_t[2];
    
	for(int i=0; i<2; i++) _array_neurongroup__spikespace[i] = 0;

	_array_neurongroup__spikespace_1 = new int32_t[2];
    
	for(int i=0; i<2; i++) _array_neurongroup__spikespace_1[i] = 0;

	_array_neurongroup_i = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_i[i] = 0;

	_array_neurongroup_I_0 = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_I_0[i] = 0;

	_array_neurongroup_I_0_1 = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_I_0_1[i] = 0;

	_array_neurongroup_I_1 = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_I_1[i] = 0;

	_array_neurongroup_i_1 = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_i_1[i] = 0;

	_array_neurongroup_I_1_1 = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_I_1_1[i] = 0;

	_array_neurongroup_lastspike = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_lastspike[i] = 0;

	_array_neurongroup_lastspike_1 = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_lastspike_1[i] = 0;

	_array_neurongroup_not_refractory = new char[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_not_refractory[i] = 0;

	_array_neurongroup_not_refractory_1 = new char[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_not_refractory_1[i] = 0;

	_array_neurongroup_Th_s = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_Th_s[i] = 0;

	_array_neurongroup_Th_s_1 = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_Th_s_1[i] = 0;

	_array_neurongroup_Th_v = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_Th_v[i] = 0;

	_array_neurongroup_Th_v_1 = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_Th_v_1[i] = 0;

	_array_neurongroup_V = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_V[i] = 0;

	_array_neurongroup_V_1 = new double[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_V_1[i] = 0;

	_array_spikemonitor__source_idx = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_spikemonitor__source_idx[i] = 0;

	_array_spikemonitor__source_idx_1 = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_spikemonitor__source_idx_1[i] = 0;

	_array_spikemonitor_count = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_spikemonitor_count[i] = 0;

	_array_spikemonitor_count_1 = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_spikemonitor_count_1[i] = 0;

	_array_spikemonitor_N = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_spikemonitor_N[i] = 0;

	_array_spikemonitor_N_1 = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_spikemonitor_N_1[i] = 0;

	_array_statemonitor__indices = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_statemonitor__indices[i] = 0;

	_array_statemonitor__indices_1 = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_statemonitor__indices_1[i] = 0;

	_array_statemonitor_N = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_statemonitor_N[i] = 0;

	_array_statemonitor_N_1 = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_statemonitor_N_1[i] = 0;


	// Arrays initialized to an "arange"
	_array_neurongroup_i = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_i[i] = 0 + i;

	_array_neurongroup_i_1 = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_neurongroup_i_1[i] = 0 + i;

	_array_spikemonitor__source_idx = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_spikemonitor__source_idx[i] = 0 + i;

	_array_spikemonitor__source_idx_1 = new int32_t[1];
    
	for(int i=0; i<1; i++) _array_spikemonitor__source_idx_1[i] = 0 + i;


	// static arrays
	_timedarray_1_values = new double[20000];
	_timedarray_values = new double[20000];

	// Random number generator states
	for (int i=0; i<1; i++)
	    _mersenne_twister_states.push_back(new rk_state());
}

void _load_arrays()
{
	using namespace brian;

	ifstream f_timedarray_1_values;
	f_timedarray_1_values.open("static_arrays/_timedarray_1_values", ios::in | ios::binary);
	if(f_timedarray_1_values.is_open())
	{
		f_timedarray_1_values.read(reinterpret_cast<char*>(_timedarray_1_values), 20000*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _timedarray_1_values." << endl;
	}
	ifstream f_timedarray_values;
	f_timedarray_values.open("static_arrays/_timedarray_values", ios::in | ios::binary);
	if(f_timedarray_values.is_open())
	{
		f_timedarray_values.read(reinterpret_cast<char*>(_timedarray_values), 20000*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _timedarray_values." << endl;
	}
}

void _write_arrays()
{
	using namespace brian;

	ofstream outfile__array_defaultclock_dt;
	outfile__array_defaultclock_dt.open("results\\_array_defaultclock_dt_1978099143", ios::binary | ios::out);
	if(outfile__array_defaultclock_dt.is_open())
	{
		outfile__array_defaultclock_dt.write(reinterpret_cast<char*>(_array_defaultclock_dt), 1*sizeof(_array_defaultclock_dt[0]));
		outfile__array_defaultclock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_dt." << endl;
	}
	ofstream outfile__array_defaultclock_t;
	outfile__array_defaultclock_t.open("results\\_array_defaultclock_t_2669362164", ios::binary | ios::out);
	if(outfile__array_defaultclock_t.is_open())
	{
		outfile__array_defaultclock_t.write(reinterpret_cast<char*>(_array_defaultclock_t), 1*sizeof(_array_defaultclock_t[0]));
		outfile__array_defaultclock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_t." << endl;
	}
	ofstream outfile__array_defaultclock_timestep;
	outfile__array_defaultclock_timestep.open("results\\_array_defaultclock_timestep_144223508", ios::binary | ios::out);
	if(outfile__array_defaultclock_timestep.is_open())
	{
		outfile__array_defaultclock_timestep.write(reinterpret_cast<char*>(_array_defaultclock_timestep), 1*sizeof(_array_defaultclock_timestep[0]));
		outfile__array_defaultclock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_timestep." << endl;
	}
	ofstream outfile__array_neurongroup__spikespace;
	outfile__array_neurongroup__spikespace.open("results\\_array_neurongroup__spikespace_3522821529", ios::binary | ios::out);
	if(outfile__array_neurongroup__spikespace.is_open())
	{
		outfile__array_neurongroup__spikespace.write(reinterpret_cast<char*>(_array_neurongroup__spikespace), 2*sizeof(_array_neurongroup__spikespace[0]));
		outfile__array_neurongroup__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup__spikespace." << endl;
	}
	ofstream outfile__array_neurongroup__spikespace_1;
	outfile__array_neurongroup__spikespace_1.open("results\\_array_neurongroup__spikespace_1_1319845205", ios::binary | ios::out);
	if(outfile__array_neurongroup__spikespace_1.is_open())
	{
		outfile__array_neurongroup__spikespace_1.write(reinterpret_cast<char*>(_array_neurongroup__spikespace_1), 2*sizeof(_array_neurongroup__spikespace_1[0]));
		outfile__array_neurongroup__spikespace_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup__spikespace_1." << endl;
	}
	ofstream outfile__array_neurongroup_i;
	outfile__array_neurongroup_i.open("results\\_array_neurongroup_i_2649026944", ios::binary | ios::out);
	if(outfile__array_neurongroup_i.is_open())
	{
		outfile__array_neurongroup_i.write(reinterpret_cast<char*>(_array_neurongroup_i), 1*sizeof(_array_neurongroup_i[0]));
		outfile__array_neurongroup_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_i." << endl;
	}
	ofstream outfile__array_neurongroup_I_0;
	outfile__array_neurongroup_I_0.open("results\\_array_neurongroup_I_0_2472010253", ios::binary | ios::out);
	if(outfile__array_neurongroup_I_0.is_open())
	{
		outfile__array_neurongroup_I_0.write(reinterpret_cast<char*>(_array_neurongroup_I_0), 1*sizeof(_array_neurongroup_I_0[0]));
		outfile__array_neurongroup_I_0.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_I_0." << endl;
	}
	ofstream outfile__array_neurongroup_I_0_1;
	outfile__array_neurongroup_I_0_1.open("results\\_array_neurongroup_I_0_1_3599645968", ios::binary | ios::out);
	if(outfile__array_neurongroup_I_0_1.is_open())
	{
		outfile__array_neurongroup_I_0_1.write(reinterpret_cast<char*>(_array_neurongroup_I_0_1), 1*sizeof(_array_neurongroup_I_0_1[0]));
		outfile__array_neurongroup_I_0_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_I_0_1." << endl;
	}
	ofstream outfile__array_neurongroup_I_1;
	outfile__array_neurongroup_I_1.open("results\\_array_neurongroup_I_1_3830502043", ios::binary | ios::out);
	if(outfile__array_neurongroup_I_1.is_open())
	{
		outfile__array_neurongroup_I_1.write(reinterpret_cast<char*>(_array_neurongroup_I_1), 1*sizeof(_array_neurongroup_I_1[0]));
		outfile__array_neurongroup_I_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_I_1." << endl;
	}
	ofstream outfile__array_neurongroup_i_1;
	outfile__array_neurongroup_i_1.open("results\\_array_neurongroup_i_1_3692926075", ios::binary | ios::out);
	if(outfile__array_neurongroup_i_1.is_open())
	{
		outfile__array_neurongroup_i_1.write(reinterpret_cast<char*>(_array_neurongroup_i_1), 1*sizeof(_array_neurongroup_i_1[0]));
		outfile__array_neurongroup_i_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_i_1." << endl;
	}
	ofstream outfile__array_neurongroup_I_1_1;
	outfile__array_neurongroup_I_1_1.open("results\\_array_neurongroup_I_1_1_3612104487", ios::binary | ios::out);
	if(outfile__array_neurongroup_I_1_1.is_open())
	{
		outfile__array_neurongroup_I_1_1.write(reinterpret_cast<char*>(_array_neurongroup_I_1_1), 1*sizeof(_array_neurongroup_I_1_1[0]));
		outfile__array_neurongroup_I_1_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_I_1_1." << endl;
	}
	ofstream outfile__array_neurongroup_lastspike;
	outfile__array_neurongroup_lastspike.open("results\\_array_neurongroup_lastspike_1647074423", ios::binary | ios::out);
	if(outfile__array_neurongroup_lastspike.is_open())
	{
		outfile__array_neurongroup_lastspike.write(reinterpret_cast<char*>(_array_neurongroup_lastspike), 1*sizeof(_array_neurongroup_lastspike[0]));
		outfile__array_neurongroup_lastspike.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_lastspike." << endl;
	}
	ofstream outfile__array_neurongroup_lastspike_1;
	outfile__array_neurongroup_lastspike_1.open("results\\_array_neurongroup_lastspike_1_797426588", ios::binary | ios::out);
	if(outfile__array_neurongroup_lastspike_1.is_open())
	{
		outfile__array_neurongroup_lastspike_1.write(reinterpret_cast<char*>(_array_neurongroup_lastspike_1), 1*sizeof(_array_neurongroup_lastspike_1[0]));
		outfile__array_neurongroup_lastspike_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_lastspike_1." << endl;
	}
	ofstream outfile__array_neurongroup_not_refractory;
	outfile__array_neurongroup_not_refractory.open("results\\_array_neurongroup_not_refractory_1422681464", ios::binary | ios::out);
	if(outfile__array_neurongroup_not_refractory.is_open())
	{
		outfile__array_neurongroup_not_refractory.write(reinterpret_cast<char*>(_array_neurongroup_not_refractory), 1*sizeof(_array_neurongroup_not_refractory[0]));
		outfile__array_neurongroup_not_refractory.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_not_refractory." << endl;
	}
	ofstream outfile__array_neurongroup_not_refractory_1;
	outfile__array_neurongroup_not_refractory_1.open("results\\_array_neurongroup_not_refractory_1_4145701307", ios::binary | ios::out);
	if(outfile__array_neurongroup_not_refractory_1.is_open())
	{
		outfile__array_neurongroup_not_refractory_1.write(reinterpret_cast<char*>(_array_neurongroup_not_refractory_1), 1*sizeof(_array_neurongroup_not_refractory_1[0]));
		outfile__array_neurongroup_not_refractory_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_not_refractory_1." << endl;
	}
	ofstream outfile__array_neurongroup_Th_s;
	outfile__array_neurongroup_Th_s.open("results\\_array_neurongroup_Th_s_69136235", ios::binary | ios::out);
	if(outfile__array_neurongroup_Th_s.is_open())
	{
		outfile__array_neurongroup_Th_s.write(reinterpret_cast<char*>(_array_neurongroup_Th_s), 1*sizeof(_array_neurongroup_Th_s[0]));
		outfile__array_neurongroup_Th_s.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_Th_s." << endl;
	}
	ofstream outfile__array_neurongroup_Th_s_1;
	outfile__array_neurongroup_Th_s_1.open("results\\_array_neurongroup_Th_s_1_2610436805", ios::binary | ios::out);
	if(outfile__array_neurongroup_Th_s_1.is_open())
	{
		outfile__array_neurongroup_Th_s_1.write(reinterpret_cast<char*>(_array_neurongroup_Th_s_1), 1*sizeof(_array_neurongroup_Th_s_1[0]));
		outfile__array_neurongroup_Th_s_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_Th_s_1." << endl;
	}
	ofstream outfile__array_neurongroup_Th_v;
	outfile__array_neurongroup_Th_v.open("results\\_array_neurongroup_Th_v_1953766372", ios::binary | ios::out);
	if(outfile__array_neurongroup_Th_v.is_open())
	{
		outfile__array_neurongroup_Th_v.write(reinterpret_cast<char*>(_array_neurongroup_Th_v), 1*sizeof(_array_neurongroup_Th_v[0]));
		outfile__array_neurongroup_Th_v.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_Th_v." << endl;
	}
	ofstream outfile__array_neurongroup_Th_v_1;
	outfile__array_neurongroup_Th_v_1.open("results\\_array_neurongroup_Th_v_1_2639517742", ios::binary | ios::out);
	if(outfile__array_neurongroup_Th_v_1.is_open())
	{
		outfile__array_neurongroup_Th_v_1.write(reinterpret_cast<char*>(_array_neurongroup_Th_v_1), 1*sizeof(_array_neurongroup_Th_v_1[0]));
		outfile__array_neurongroup_Th_v_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_Th_v_1." << endl;
	}
	ofstream outfile__array_neurongroup_V;
	outfile__array_neurongroup_V.open("results\\_array_neurongroup_V_729996477", ios::binary | ios::out);
	if(outfile__array_neurongroup_V.is_open())
	{
		outfile__array_neurongroup_V.write(reinterpret_cast<char*>(_array_neurongroup_V), 1*sizeof(_array_neurongroup_V[0]));
		outfile__array_neurongroup_V.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_V." << endl;
	}
	ofstream outfile__array_neurongroup_V_1;
	outfile__array_neurongroup_V_1.open("results\\_array_neurongroup_V_1_4079630038", ios::binary | ios::out);
	if(outfile__array_neurongroup_V_1.is_open())
	{
		outfile__array_neurongroup_V_1.write(reinterpret_cast<char*>(_array_neurongroup_V_1), 1*sizeof(_array_neurongroup_V_1[0]));
		outfile__array_neurongroup_V_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_V_1." << endl;
	}
	ofstream outfile__array_spikemonitor__source_idx;
	outfile__array_spikemonitor__source_idx.open("results\\_array_spikemonitor__source_idx_1477951789", ios::binary | ios::out);
	if(outfile__array_spikemonitor__source_idx.is_open())
	{
		outfile__array_spikemonitor__source_idx.write(reinterpret_cast<char*>(_array_spikemonitor__source_idx), 1*sizeof(_array_spikemonitor__source_idx[0]));
		outfile__array_spikemonitor__source_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor__source_idx." << endl;
	}
	ofstream outfile__array_spikemonitor__source_idx_1;
	outfile__array_spikemonitor__source_idx_1.open("results\\_array_spikemonitor__source_idx_1_3782058880", ios::binary | ios::out);
	if(outfile__array_spikemonitor__source_idx_1.is_open())
	{
		outfile__array_spikemonitor__source_idx_1.write(reinterpret_cast<char*>(_array_spikemonitor__source_idx_1), 1*sizeof(_array_spikemonitor__source_idx_1[0]));
		outfile__array_spikemonitor__source_idx_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor__source_idx_1." << endl;
	}
	ofstream outfile__array_spikemonitor_count;
	outfile__array_spikemonitor_count.open("results\\_array_spikemonitor_count_598337445", ios::binary | ios::out);
	if(outfile__array_spikemonitor_count.is_open())
	{
		outfile__array_spikemonitor_count.write(reinterpret_cast<char*>(_array_spikemonitor_count), 1*sizeof(_array_spikemonitor_count[0]));
		outfile__array_spikemonitor_count.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_count." << endl;
	}
	ofstream outfile__array_spikemonitor_count_1;
	outfile__array_spikemonitor_count_1.open("results\\_array_spikemonitor_count_1_3225046912", ios::binary | ios::out);
	if(outfile__array_spikemonitor_count_1.is_open())
	{
		outfile__array_spikemonitor_count_1.write(reinterpret_cast<char*>(_array_spikemonitor_count_1), 1*sizeof(_array_spikemonitor_count_1[0]));
		outfile__array_spikemonitor_count_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_count_1." << endl;
	}
	ofstream outfile__array_spikemonitor_N;
	outfile__array_spikemonitor_N.open("results\\_array_spikemonitor_N_225734567", ios::binary | ios::out);
	if(outfile__array_spikemonitor_N.is_open())
	{
		outfile__array_spikemonitor_N.write(reinterpret_cast<char*>(_array_spikemonitor_N), 1*sizeof(_array_spikemonitor_N[0]));
		outfile__array_spikemonitor_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_N." << endl;
	}
	ofstream outfile__array_spikemonitor_N_1;
	outfile__array_spikemonitor_N_1.open("results\\_array_spikemonitor_N_1_292489421", ios::binary | ios::out);
	if(outfile__array_spikemonitor_N_1.is_open())
	{
		outfile__array_spikemonitor_N_1.write(reinterpret_cast<char*>(_array_spikemonitor_N_1), 1*sizeof(_array_spikemonitor_N_1[0]));
		outfile__array_spikemonitor_N_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_N_1." << endl;
	}
	ofstream outfile__array_statemonitor__indices;
	outfile__array_statemonitor__indices.open("results\\_array_statemonitor__indices_2854283999", ios::binary | ios::out);
	if(outfile__array_statemonitor__indices.is_open())
	{
		outfile__array_statemonitor__indices.write(reinterpret_cast<char*>(_array_statemonitor__indices), 1*sizeof(_array_statemonitor__indices[0]));
		outfile__array_statemonitor__indices.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor__indices." << endl;
	}
	ofstream outfile__array_statemonitor__indices_1;
	outfile__array_statemonitor__indices_1.open("results\\_array_statemonitor__indices_1_1801137096", ios::binary | ios::out);
	if(outfile__array_statemonitor__indices_1.is_open())
	{
		outfile__array_statemonitor__indices_1.write(reinterpret_cast<char*>(_array_statemonitor__indices_1), 1*sizeof(_array_statemonitor__indices_1[0]));
		outfile__array_statemonitor__indices_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor__indices_1." << endl;
	}
	ofstream outfile__array_statemonitor_N;
	outfile__array_statemonitor_N.open("results\\_array_statemonitor_N_4140778434", ios::binary | ios::out);
	if(outfile__array_statemonitor_N.is_open())
	{
		outfile__array_statemonitor_N.write(reinterpret_cast<char*>(_array_statemonitor_N), 1*sizeof(_array_statemonitor_N[0]));
		outfile__array_statemonitor_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_N." << endl;
	}
	ofstream outfile__array_statemonitor_N_1;
	outfile__array_statemonitor_N_1.open("results\\_array_statemonitor_N_1_993853559", ios::binary | ios::out);
	if(outfile__array_statemonitor_N_1.is_open())
	{
		outfile__array_statemonitor_N_1.write(reinterpret_cast<char*>(_array_statemonitor_N_1), 1*sizeof(_array_statemonitor_N_1[0]));
		outfile__array_statemonitor_N_1.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_N_1." << endl;
	}

	ofstream outfile__dynamic_array_spikemonitor_i;
	outfile__dynamic_array_spikemonitor_i.open("results\\_dynamic_array_spikemonitor_i_1976709050", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_i.is_open())
	{
        if (! _dynamic_array_spikemonitor_i.empty() )
        {
			outfile__dynamic_array_spikemonitor_i.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_i[0]), _dynamic_array_spikemonitor_i.size()*sizeof(_dynamic_array_spikemonitor_i[0]));
		    outfile__dynamic_array_spikemonitor_i.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_i." << endl;
	}
	ofstream outfile__dynamic_array_spikemonitor_i_1;
	outfile__dynamic_array_spikemonitor_i_1.open("results\\_dynamic_array_spikemonitor_i_1_2564775399", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_i_1.is_open())
	{
        if (! _dynamic_array_spikemonitor_i_1.empty() )
        {
			outfile__dynamic_array_spikemonitor_i_1.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_i_1[0]), _dynamic_array_spikemonitor_i_1.size()*sizeof(_dynamic_array_spikemonitor_i_1[0]));
		    outfile__dynamic_array_spikemonitor_i_1.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_i_1." << endl;
	}
	ofstream outfile__dynamic_array_spikemonitor_t;
	outfile__dynamic_array_spikemonitor_t.open("results\\_dynamic_array_spikemonitor_t_383009635", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_t.is_open())
	{
        if (! _dynamic_array_spikemonitor_t.empty() )
        {
			outfile__dynamic_array_spikemonitor_t.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_t[0]), _dynamic_array_spikemonitor_t.size()*sizeof(_dynamic_array_spikemonitor_t[0]));
		    outfile__dynamic_array_spikemonitor_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_t." << endl;
	}
	ofstream outfile__dynamic_array_spikemonitor_t_1;
	outfile__dynamic_array_spikemonitor_t_1.open("results\\_dynamic_array_spikemonitor_t_1_2351001028", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_t_1.is_open())
	{
        if (! _dynamic_array_spikemonitor_t_1.empty() )
        {
			outfile__dynamic_array_spikemonitor_t_1.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_t_1[0]), _dynamic_array_spikemonitor_t_1.size()*sizeof(_dynamic_array_spikemonitor_t_1[0]));
		    outfile__dynamic_array_spikemonitor_t_1.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_t_1." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_t;
	outfile__dynamic_array_statemonitor_t.open("results\\_dynamic_array_statemonitor_t_3983503110", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_t.is_open())
	{
        if (! _dynamic_array_statemonitor_t.empty() )
        {
			outfile__dynamic_array_statemonitor_t.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_t[0]), _dynamic_array_statemonitor_t.size()*sizeof(_dynamic_array_statemonitor_t[0]));
		    outfile__dynamic_array_statemonitor_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_t." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_t_1;
	outfile__dynamic_array_statemonitor_t_1.open("results\\_dynamic_array_statemonitor_t_1_2792580478", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_t_1.is_open())
	{
        if (! _dynamic_array_statemonitor_t_1.empty() )
        {
			outfile__dynamic_array_statemonitor_t_1.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_t_1[0]), _dynamic_array_statemonitor_t_1.size()*sizeof(_dynamic_array_statemonitor_t_1[0]));
		    outfile__dynamic_array_statemonitor_t_1.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_t_1." << endl;
	}

	ofstream outfile__dynamic_array_statemonitor_I_0;
	outfile__dynamic_array_statemonitor_I_0.open("results\\_dynamic_array_statemonitor_I_0_4257686315", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_I_0.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_I_0.n; n++)
        {
            if (! _dynamic_array_statemonitor_I_0(n).empty())
            {
                outfile__dynamic_array_statemonitor_I_0.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_I_0(n, 0)), _dynamic_array_statemonitor_I_0.m*sizeof(_dynamic_array_statemonitor_I_0(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_I_0.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_I_0." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_I_0_1;
	outfile__dynamic_array_statemonitor_I_0_1.open("results\\_dynamic_array_statemonitor_I_0_1_4180877954", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_I_0_1.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_I_0_1.n; n++)
        {
            if (! _dynamic_array_statemonitor_I_0_1(n).empty())
            {
                outfile__dynamic_array_statemonitor_I_0_1.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_I_0_1(n, 0)), _dynamic_array_statemonitor_I_0_1.m*sizeof(_dynamic_array_statemonitor_I_0_1(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_I_0_1.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_I_0_1." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_I_1;
	outfile__dynamic_array_statemonitor_I_1.open("results\\_dynamic_array_statemonitor_I_1_2327843773", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_I_1.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_I_1.n; n++)
        {
            if (! _dynamic_array_statemonitor_I_1(n).empty())
            {
                outfile__dynamic_array_statemonitor_I_1.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_I_1(n, 0)), _dynamic_array_statemonitor_I_1.m*sizeof(_dynamic_array_statemonitor_I_1(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_I_1.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_I_1." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_I_1_1;
	outfile__dynamic_array_statemonitor_I_1_1.open("results\\_dynamic_array_statemonitor_I_1_1_4176562357", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_I_1_1.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_I_1_1.n; n++)
        {
            if (! _dynamic_array_statemonitor_I_1_1(n).empty())
            {
                outfile__dynamic_array_statemonitor_I_1_1.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_I_1_1(n, 0)), _dynamic_array_statemonitor_I_1_1.m*sizeof(_dynamic_array_statemonitor_I_1_1(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_I_1_1.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_I_1_1." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_Th_s;
	outfile__dynamic_array_statemonitor_Th_s.open("results\\_dynamic_array_statemonitor_Th_s_3598580311", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_Th_s.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_Th_s.n; n++)
        {
            if (! _dynamic_array_statemonitor_Th_s(n).empty())
            {
                outfile__dynamic_array_statemonitor_Th_s.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_Th_s(n, 0)), _dynamic_array_statemonitor_Th_s.m*sizeof(_dynamic_array_statemonitor_Th_s(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_Th_s.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_Th_s." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_Th_s_1;
	outfile__dynamic_array_statemonitor_Th_s_1.open("results\\_dynamic_array_statemonitor_Th_s_1_2243319218", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_Th_s_1.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_Th_s_1.n; n++)
        {
            if (! _dynamic_array_statemonitor_Th_s_1(n).empty())
            {
                outfile__dynamic_array_statemonitor_Th_s_1.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_Th_s_1(n, 0)), _dynamic_array_statemonitor_Th_s_1.m*sizeof(_dynamic_array_statemonitor_Th_s_1(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_Th_s_1.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_Th_s_1." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_Th_v;
	outfile__dynamic_array_statemonitor_Th_v.open("results\\_dynamic_array_statemonitor_Th_v_2786528984", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_Th_v.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_Th_v.n; n++)
        {
            if (! _dynamic_array_statemonitor_Th_v(n).empty())
            {
                outfile__dynamic_array_statemonitor_Th_v.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_Th_v(n, 0)), _dynamic_array_statemonitor_Th_v.m*sizeof(_dynamic_array_statemonitor_Th_v(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_Th_v.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_Th_v." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_Th_v_1;
	outfile__dynamic_array_statemonitor_Th_v_1.open("results\\_dynamic_array_statemonitor_Th_v_1_2206046041", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_Th_v_1.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_Th_v_1.n; n++)
        {
            if (! _dynamic_array_statemonitor_Th_v_1(n).empty())
            {
                outfile__dynamic_array_statemonitor_Th_v_1.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_Th_v_1(n, 0)), _dynamic_array_statemonitor_Th_v_1.m*sizeof(_dynamic_array_statemonitor_Th_v_1(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_Th_v_1.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_Th_v_1." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_V;
	outfile__dynamic_array_statemonitor_V.open("results\\_dynamic_array_statemonitor_V_940519138", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_V.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_V.n; n++)
        {
            if (! _dynamic_array_statemonitor_V(n).empty())
            {
                outfile__dynamic_array_statemonitor_V.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_V(n, 0)), _dynamic_array_statemonitor_V.m*sizeof(_dynamic_array_statemonitor_V(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_V.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_V." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_V_1;
	outfile__dynamic_array_statemonitor_V_1.open("results\\_dynamic_array_statemonitor_V_1_2646276080", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_V_1.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_V_1.n; n++)
        {
            if (! _dynamic_array_statemonitor_V_1(n).empty())
            {
                outfile__dynamic_array_statemonitor_V_1.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_V_1(n, 0)), _dynamic_array_statemonitor_V_1.m*sizeof(_dynamic_array_statemonitor_V_1(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_V_1.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_V_1." << endl;
	}
	// Write last run info to disk
	ofstream outfile_last_run_info;
	outfile_last_run_info.open("results/last_run_info.txt", ios::out);
	if(outfile_last_run_info.is_open())
	{
		outfile_last_run_info << (Network::_last_run_time) << " " << (Network::_last_run_completed_fraction) << std::endl;
		outfile_last_run_info.close();
	} else
	{
	    std::cout << "Error writing last run info to file." << std::endl;
	}
}

void _dealloc_arrays()
{
	using namespace brian;


	// static arrays
	if(_timedarray_1_values!=0)
	{
		delete [] _timedarray_1_values;
		_timedarray_1_values = 0;
	}
	if(_timedarray_values!=0)
	{
		delete [] _timedarray_values;
		_timedarray_values = 0;
	}
}

