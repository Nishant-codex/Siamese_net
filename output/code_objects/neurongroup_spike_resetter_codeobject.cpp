#include "code_objects/neurongroup_spike_resetter_codeobject.h"
#include "objects.h"
#include "brianlib/common_math.h"
#include "brianlib/stdint_compat.h"
#include<cmath>
#include<ctime>
#include<iostream>
#include<fstream>
#include<climits>

////// SUPPORT CODE ///////
namespace {
        
    template < typename T1, typename T2 > struct _higher_type;
    template < > struct _higher_type<int,int> { typedef int type; };
    template < > struct _higher_type<int,long> { typedef long type; };
    template < > struct _higher_type<int,long long> { typedef long long type; };
    template < > struct _higher_type<int,float> { typedef float type; };
    template < > struct _higher_type<int,double> { typedef double type; };
    template < > struct _higher_type<int,long double> { typedef long double type; };
    template < > struct _higher_type<long,int> { typedef long type; };
    template < > struct _higher_type<long,long> { typedef long type; };
    template < > struct _higher_type<long,long long> { typedef long long type; };
    template < > struct _higher_type<long,float> { typedef float type; };
    template < > struct _higher_type<long,double> { typedef double type; };
    template < > struct _higher_type<long,long double> { typedef long double type; };
    template < > struct _higher_type<long long,int> { typedef long long type; };
    template < > struct _higher_type<long long,long> { typedef long long type; };
    template < > struct _higher_type<long long,long long> { typedef long long type; };
    template < > struct _higher_type<long long,float> { typedef float type; };
    template < > struct _higher_type<long long,double> { typedef double type; };
    template < > struct _higher_type<long long,long double> { typedef long double type; };
    template < > struct _higher_type<float,int> { typedef float type; };
    template < > struct _higher_type<float,long> { typedef float type; };
    template < > struct _higher_type<float,long long> { typedef float type; };
    template < > struct _higher_type<float,float> { typedef float type; };
    template < > struct _higher_type<float,double> { typedef double type; };
    template < > struct _higher_type<float,long double> { typedef long double type; };
    template < > struct _higher_type<double,int> { typedef double type; };
    template < > struct _higher_type<double,long> { typedef double type; };
    template < > struct _higher_type<double,long long> { typedef double type; };
    template < > struct _higher_type<double,float> { typedef double type; };
    template < > struct _higher_type<double,double> { typedef double type; };
    template < > struct _higher_type<double,long double> { typedef long double type; };
    template < > struct _higher_type<long double,int> { typedef long double type; };
    template < > struct _higher_type<long double,long> { typedef long double type; };
    template < > struct _higher_type<long double,long long> { typedef long double type; };
    template < > struct _higher_type<long double,float> { typedef long double type; };
    template < > struct _higher_type<long double,double> { typedef long double type; };
    template < > struct _higher_type<long double,long double> { typedef long double type; };
    template < typename T1, typename T2 >
    static inline typename _higher_type<T1,T2>::type
    _brian_mod(T1 x, T2 y)
    {{
        return x-y*floor(1.0*x/y);
    }}
    template < typename T1, typename T2 >
    static inline typename _higher_type<T1,T2>::type
    _brian_floordiv(T1 x, T2 y)
    {{
        return floor(1.0*x/y);
    }}
    #ifdef _MSC_VER
    #define _brian_pow(x, y) (pow((double)(x), (y)))
    #else
    #define _brian_pow(x, y) (pow((x), (y)))
    #endif

}

////// HASH DEFINES ///////



void _run_neurongroup_spike_resetter_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double A_0 = 4.89490361114548e-12;
const double A_1 = 5.008484746493212e-12;
const size_t _numI_0 = 1;
const size_t _numI_1 = 1;
const int32_t N = 1;
const double R_0 = 0.22058427457755816;
const double R_1 = 22.88393450483256;
const size_t _numTh_s = 1;
const size_t _numTh_v = 1;
const size_t _numV = 1;
const size_t _num_spikespace = 2;
const double a_r = 44.95461156822592;
const double a_s = - 0.029014898397803815;
const double b_r = - 0.04135663245495133;
const size_t _numdt = 1;
const double k_0 = 76.2517802375484;
const double k_1 = 72.18184923084418;
const double t_ref = 0.004243152630637912;
const size_t _numnot_refractory = 1;
    ///// POINTERS ////////////
        
    double*   _ptr_array_neurongroup_I_0 = _array_neurongroup_I_0;
    double*   _ptr_array_neurongroup_I_1 = _array_neurongroup_I_1;
    double*   _ptr_array_neurongroup_Th_s = _array_neurongroup_Th_s;
    double*   _ptr_array_neurongroup_Th_v = _array_neurongroup_Th_v;
    double*   _ptr_array_neurongroup_V = _array_neurongroup_V;
    int32_t* __restrict  _ptr_array_neurongroup__spikespace = _array_neurongroup__spikespace;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    char*   _ptr_array_neurongroup_not_refractory = _array_neurongroup_not_refractory;



	const int32_t *_events = _ptr_array_neurongroup__spikespace;
	const int32_t _num_events = _ptr_array_neurongroup__spikespace[N];

	//// MAIN CODE ////////////	
	// scalar code
	const size_t _vectorisation_idx = -1;
 	
 const double dt = _ptr_array_defaultclock_dt[0];
 const double _lio_1 = R_0 * exp((- k_0) * (t_ref - dt));
 const double _lio_2 = R_1 * exp((- k_1) * (t_ref - dt));


	
	for(int32_t _index_events=0; _index_events<_num_events; _index_events++)
	{
	    // vector code
		const size_t _idx = _events[_index_events];
		const size_t _vectorisation_idx = _idx;
                
        double I_0 = _ptr_array_neurongroup_I_0[_idx];
        double I_1 = _ptr_array_neurongroup_I_1[_idx];
        double Th_s = _ptr_array_neurongroup_Th_s[_idx];
        double Th_v = _ptr_array_neurongroup_Th_v[_idx];
        double V = _ptr_array_neurongroup_V[_idx];
        V = b_r + (a_r * V);
        Th_s += a_s;
        Th_v *= 1;
        I_0 = A_0 + (_lio_1 * I_0);
        I_1 = A_1 + (_lio_2 * I_1);
        _ptr_array_neurongroup_I_0[_idx] = I_0;
        _ptr_array_neurongroup_I_1[_idx] = I_1;
        _ptr_array_neurongroup_Th_s[_idx] = Th_s;
        _ptr_array_neurongroup_Th_v[_idx] = Th_v;
        _ptr_array_neurongroup_V[_idx] = V;

	}

}


