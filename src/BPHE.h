#ifndef BPHE_H
#define BPHE_H

class BrazedPlateHeatExchanger
{
public:
	/// Inlet state for the hot stream
	CoolPropStateClass State_h_inlet;
		
	/// Inlet state for the cold stream
	CoolPropStateClass State_c_inlet;

	/// Mass flow rate of cold stream [kg/s]
	double mdot_c;

	/// Mass flow rate of hot stream [kg/s]
	double mdot_h;

	BrazedPlateHeatExchanger(){};
	void test();
	void calculate();
	double DetermineQmax(void);
};

#endif