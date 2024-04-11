#include <stdlib.h>
#include "objects.h"
#include <ctime>
#include <time.h>

#include "run.h"
#include "brianlib/common_math.h"
#include "randomkit.h"

#include "code_objects/neurongroup_spike_resetter_codeobject.h"
#include "code_objects/neurongroup_spike_resetter_codeobject_1.h"
#include "code_objects/neurongroup_spike_thresholder_codeobject.h"
#include "code_objects/after_run_neurongroup_spike_thresholder_codeobject.h"
#include "code_objects/neurongroup_spike_thresholder_codeobject_1.h"
#include "code_objects/after_run_neurongroup_spike_thresholder_codeobject_1.h"
#include "code_objects/neurongroup_stateupdater_codeobject.h"
#include "code_objects/neurongroup_stateupdater_codeobject_1.h"
#include "code_objects/spikemonitor_codeobject.h"
#include "code_objects/spikemonitor_codeobject_1.h"
#include "code_objects/statemonitor_codeobject.h"
#include "code_objects/statemonitor_codeobject_1.h"


#include <iostream>
#include <fstream>
#include <string>




int main(int argc, char **argv)
{
        

	brian_start();
        

	{
		using namespace brian;

		
                
        _array_defaultclock_dt[0] = 0.0001;
        _array_defaultclock_dt[0] = 0.0001;
        _array_defaultclock_dt[0] = 0.0001;
        _array_neurongroup_lastspike[0] = - 10000.0;
        _array_neurongroup_not_refractory[0] = true;
        _array_neurongroup_V[0] = - 0.052249;
        _array_neurongroup_Th_s[0] = - 0.03;
        _array_neurongroup_Th_v[0] = - 0.03;
        _array_neurongroup_I_0[0] = 1e-12;
        _array_neurongroup_I_1[0] = 1e-12;
        _array_statemonitor__indices[0] = 0;
        _array_defaultclock_dt[0] = 5e-05;
        _array_defaultclock_timestep[0] = 0;
        _array_defaultclock_t[0] = 0.0;
        magicnetwork.clear();
        magicnetwork.add(&defaultclock, _run_statemonitor_codeobject);
        magicnetwork.add(&defaultclock, _run_neurongroup_stateupdater_codeobject);
        magicnetwork.add(&defaultclock, _run_neurongroup_spike_thresholder_codeobject);
        magicnetwork.add(&defaultclock, _run_spikemonitor_codeobject);
        magicnetwork.add(&defaultclock, _run_neurongroup_spike_resetter_codeobject);
        magicnetwork.run(1.0, NULL, 10.0);
        _after_run_neurongroup_spike_thresholder_codeobject();
        _array_defaultclock_dt[0] = 5e-05;
        _array_neurongroup_lastspike_1[0] = - 10000.0;
        _array_neurongroup_not_refractory_1[0] = true;
        _array_neurongroup_V_1[0] = - 0.052249;
        _array_neurongroup_Th_s_1[0] = - 0.03;
        _array_neurongroup_Th_v_1[0] = - 0.03;
        _array_neurongroup_I_0_1[0] = 1e-12;
        _array_neurongroup_I_1_1[0] = 1e-12;
        _array_statemonitor__indices_1[0] = 0;
        _array_defaultclock_dt[0] = 5e-05;
        _array_defaultclock_timestep[0] = 0;
        _array_defaultclock_t[0] = 0.0;
        magicnetwork.clear();
        magicnetwork.add(&defaultclock, _run_statemonitor_codeobject_1);
        magicnetwork.add(&defaultclock, _run_neurongroup_stateupdater_codeobject_1);
        magicnetwork.add(&defaultclock, _run_neurongroup_spike_thresholder_codeobject_1);
        magicnetwork.add(&defaultclock, _run_spikemonitor_codeobject_1);
        magicnetwork.add(&defaultclock, _run_neurongroup_spike_resetter_codeobject_1);
        magicnetwork.run(1.0, NULL, 10.0);
        _after_run_neurongroup_spike_thresholder_codeobject_1();
        #ifdef DEBUG
        _debugmsg_spikemonitor_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_spikemonitor_codeobject_1();
        #endif

	}
        

	brian_end();
        

	return 0;
}