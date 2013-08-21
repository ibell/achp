
#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include "CPState.h"
class CompressorClass
{
public:
	double Tsat_s_K, Tsat_d_K, DT_sh_K, Tsat_s, Tsat_d, power_map, Vdot_ratio, P1, P2, F, T1_actual, v_map, v_actual, mdot, fp, eta_oi, Wdot, CycleEnergyIn, Vdot_pumped;
	std::vector<double> P, M;

	CoolPropStateClass inlet_state, outlet_state;

	void Calculate();
	void test();
};

#endif