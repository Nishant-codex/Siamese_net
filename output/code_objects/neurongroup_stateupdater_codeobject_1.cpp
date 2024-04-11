#include "code_objects/neurongroup_stateupdater_codeobject_1.h"
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
        
    static double* _namespace_timedarray_1_values;
    static inline double _timedarray_1(const double t)
    {
        const double epsilon = 0.000050000000000000 / 8;
        int i = (int)((t/epsilon + 0.5)/8);
        if(i < 0)
           i = 0;
        if(i >= 20000)
            i = 20000-1;
        return _namespace_timedarray_1_values[i];
    }
    static inline int64_t _timestep(double t, double dt)
    {
        return (int64_t)((t + 1e-3*dt)/dt);
    }
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



void _run_neurongroup_stateupdater_codeobject_1()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double C = 2.732045518525941e-11;
const double El = - 0.05240582263619105;
const double G = 8.36170725325334e-09;
const size_t _numI_0 = 1;
const size_t _numI_1 = 1;
const int32_t N = 1;
const double Th_inf = - 0.04143648409446105;
const size_t _numTh_s = 1;
const size_t _numTh_v = 1;
const size_t _numV = 1;
const double a_v = 22.32734606686461;
const double b_s = 51.70829174324437;
const double b_v = 64.85415674574465;
const size_t _numdt = 1;
const double k_0 = 57.62701168252686;
const double k_1 = 5.9726331090284805;
const size_t _numlastspike = 1;
const size_t _numnot_refractory = 1;
const size_t _numt = 1;
    ///// POINTERS ////////////
        
    double*   _ptr_array_neurongroup_I_0_1 = _array_neurongroup_I_0_1;
    double*   _ptr_array_neurongroup_I_1_1 = _array_neurongroup_I_1_1;
    double*   _ptr_array_neurongroup_Th_s_1 = _array_neurongroup_Th_s_1;
    double*   _ptr_array_neurongroup_Th_v_1 = _array_neurongroup_Th_v_1;
    double*   _ptr_array_neurongroup_V_1 = _array_neurongroup_V_1;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double*   _ptr_array_neurongroup_lastspike_1 = _array_neurongroup_lastspike_1;
    char*   _ptr_array_neurongroup_not_refractory_1 = _array_neurongroup_not_refractory_1;
    double*   _ptr_array_defaultclock_t = _array_defaultclock_t;
    _namespace_timedarray_1_values = _timedarray_1_values;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double t = _ptr_array_defaultclock_t[0];
    const int32_t _lio_1 = _timestep(0.004317, dt);
    const double _lio_2 = (- dt) * k_0;
    const double _lio_3 = exp(_lio_2);
    const double _lio_4 = (- dt) * k_1;
    const double _lio_5 = exp(_lio_4);
    const double _lio_6 = (- b_s) * dt;
    const double _lio_7 = exp(_lio_6);
    const double _lio_8 = 0.0 - (1.0f*((((- G) * Th_inf) * b_v) - (a_v * _timedarray_1(t)))/(G * b_v));
    const double _lio_9 = 1.0f*(C * a_v)/((C * b_v) - G);
    const double _lio_10 = b_v * dt;
    const double _lio_11 = 1.0f*(G * dt)/C;
    const double _lio_12 = (- b_v) * dt;
    const double _lio_13 = 1.0f*((- G) * dt)/C;
    const double _lio_14 = 1.0f*((C * a_v) * (((- El) * G) - _timedarray_1(t)))/(G * ((C * b_v) - G));
    const double _lio_15 = 1.0f*a_v/((((((_brian_pow(C, 2)) * (_brian_pow(b_v, 2))) * k_0) + ((C * G) * (_brian_pow(k_0, 2)))) + ((_brian_pow(G, 2)) * b_v)) - (((((_brian_pow(C, 2)) * b_v) * (_brian_pow(k_0, 2))) + ((C * G) * (_brian_pow(b_v, 2)))) + ((_brian_pow(G, 2)) * k_0)));
    const double _lio_16 = C * b_v;
    const double _lio_17 = dt * k_0;
    const double _lio_18 = C * k_0;
    const double _lio_19 = 1.0f*a_v/((((((_brian_pow(C, 2)) * (_brian_pow(b_v, 2))) * k_1) + ((C * G) * (_brian_pow(k_1, 2)))) + ((_brian_pow(G, 2)) * b_v)) - (((((_brian_pow(C, 2)) * b_v) * (_brian_pow(k_1, 2))) + ((C * G) * (_brian_pow(b_v, 2)))) + ((_brian_pow(G, 2)) * k_1)));
    const double _lio_20 = dt * k_1;
    const double _lio_21 = C * k_1;
    const double _lio_22 = 1.0f*((((- G) * Th_inf) * b_v) - (a_v * _timedarray_1(t)))/(G * b_v);
    const double _lio_23 = _lio_8 + _lio_22;
    const double _lio_24 = _lio_15 * (((_lio_16 + _lio_18) + G) - ((_lio_16 + _lio_18) + G));
    const double _lio_25 = _lio_19 * (((_lio_16 + _lio_21) + G) - ((_lio_16 + _lio_21) + G));
    const double _lio_26 = (_lio_8 + (((_lio_14 * (exp(_lio_10) - exp(_lio_11))) * exp(_lio_12)) * exp(_lio_13))) + (_lio_22 * exp(_lio_12));
    const double _lio_27 = ((_lio_9 * (exp(_lio_10) - exp(_lio_11))) * exp(_lio_12)) * exp(_lio_13);
    const double _lio_28 = (((_lio_15 * (((((_lio_16 * exp(_lio_10)) * exp(_lio_17)) + ((_lio_18 * exp(_lio_17)) * exp(_lio_11))) + ((G * exp(_lio_10)) * exp(_lio_11))) - ((((_lio_16 * exp(_lio_10)) * exp(_lio_11)) + ((_lio_18 * exp(_lio_10)) * exp(_lio_17))) + ((G * exp(_lio_17)) * exp(_lio_11))))) * exp(_lio_12)) * exp(_lio_2)) * exp(_lio_13);
    const double _lio_29 = (((_lio_19 * (((((_lio_16 * exp(_lio_10)) * exp(_lio_20)) + ((_lio_21 * exp(_lio_20)) * exp(_lio_11))) + ((G * exp(_lio_10)) * exp(_lio_11))) - ((((_lio_16 * exp(_lio_10)) * exp(_lio_11)) + ((_lio_21 * exp(_lio_10)) * exp(_lio_20))) + ((G * exp(_lio_20)) * exp(_lio_11))))) * exp(_lio_12)) * exp(_lio_4)) * exp(_lio_13);
    const double _lio_30 = exp(_lio_12);
    const double _lio_31 = 0.0 - (1.0f*(((- El) * G) - _timedarray_1(t))/G);
    const double _lio_32 = 1.0f*1.0/((C * k_0) - G);
    const double _lio_33 = 1.0f*1.0/((C * k_1) - G);
    const double _lio_34 = 1.0f*(((- El) * G) - _timedarray_1(t))/G;
    const double _lio_35 = _lio_31 + _lio_34;
    const double _lio_36 = _lio_31 + (_lio_34 * exp(_lio_13));
    const double _lio_37 = ((_lio_32 * (exp(_lio_17) - exp(_lio_11))) * exp(_lio_2)) * exp(_lio_13);
    const double _lio_38 = ((_lio_33 * (exp(_lio_20) - exp(_lio_11))) * exp(_lio_4)) * exp(_lio_13);
    const double _lio_39 = exp(_lio_13);


    const int _N = N;
    
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        double I_0 = _ptr_array_neurongroup_I_0_1[_idx];
        double I_1 = _ptr_array_neurongroup_I_1_1[_idx];
        double Th_s = _ptr_array_neurongroup_Th_s_1[_idx];
        double Th_v = _ptr_array_neurongroup_Th_v_1[_idx];
        double V = _ptr_array_neurongroup_V_1[_idx];
        const double lastspike = _ptr_array_neurongroup_lastspike_1[_idx];
        char not_refractory = _ptr_array_neurongroup_not_refractory_1[_idx];
        not_refractory = _timestep(t - lastspike, dt) >= _lio_1;
        double _I_0;
        if(!not_refractory)
            _I_0 = I_0;
        else 
            _I_0 = _lio_3 * I_0;
        double _I_1;
        if(!not_refractory)
            _I_1 = I_1;
        else 
            _I_1 = _lio_5 * I_1;
        double _Th_s;
        if(!not_refractory)
            _Th_s = Th_s;
        else 
            _Th_s = _lio_7 * Th_s;
        double _Th_v;
        if(!not_refractory)
            _Th_v = _lio_23 + (((_lio_24 * I_0) + (_lio_25 * I_1)) + Th_v);
        else 
            _Th_v = _lio_26 + ((((_lio_27 * V) + (_lio_28 * I_0)) + (_lio_29 * I_1)) + (_lio_30 * Th_v));
        double _V;
        if(!not_refractory)
            _V = _lio_35 + V;
        else 
            _V = _lio_36 + (((_lio_37 * I_0) + (_lio_38 * I_1)) + (_lio_39 * V));
        if(not_refractory)
            I_0 = _I_0;
        if(not_refractory)
            I_1 = _I_1;
        if(not_refractory)
            Th_s = _Th_s;
        if(not_refractory)
            Th_v = _Th_v;
        if(not_refractory)
            V = _V;
        _ptr_array_neurongroup_I_0_1[_idx] = I_0;
        _ptr_array_neurongroup_I_1_1[_idx] = I_1;
        _ptr_array_neurongroup_Th_s_1[_idx] = Th_s;
        _ptr_array_neurongroup_Th_v_1[_idx] = Th_v;
        _ptr_array_neurongroup_V_1[_idx] = V;
        _ptr_array_neurongroup_not_refractory_1[_idx] = not_refractory;

    }

}


