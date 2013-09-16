#ifndef BPHE_H
#define BPHE_H

class BrazedPlateHeatExchanger
{
public:
	/// Inlet state for the hot stream
	CoolPropStateClassSI State_h_inlet;
		
	/// Inlet state for the cold stream
	CoolPropStateClassSI State_c_inlet;

	/// Saturation state for the cold stream at the inlet pressure
	CoolPropStateClassSI State_h_sat;

	/// Saturation state for the hot stream at the inlet pressure
	CoolPropStateClassSI State_c_sat;

	/// Mass flow rate of cold stream [kg/s]
	double mdot_c;

	/// Mass flow rate of hot stream [kg/s]
	double mdot_h;

	/// The list of enthalpies for the cold stream at cell boundaries
	std::vector<double> EnthalpyList_c;
	
	/// The list of enthalpies for the hot stream at cell boundaries
	std::vector<double> EnthalpyList_h;

	/// The list of temperatures of the cold stream at cell boundaries
	std::vector<double> TemperatureList_c;

	/// The list of temperatures of the hot stream at cell boundaries
	std::vector<double> TemperatureList_h;

	/// A vector of booleans for whether cell boundary index corresponds to phase change of the hot stream
	std::vector<bool> PhaseBoundary_h;
	
	/// A vector of booleans for whether cell boundary index corresponds to phase change of the hot stream
	std::vector<bool> PhaseBoundary_c;

	BrazedPlateHeatExchanger(){this->verbosity = 1;};

	/// Build the list of enthalpies for a given heat transfer rate
	void BuildEnthalpyLists(double Q);

	/// Calculate the saturation states for both streams
	void SaturationStates();

	/// How verbose the debugging should be [0: no output, 10: very annoying output]
	int verbosity;

	void _check();
	void test();
	void calculate();
	double DetermineQmax(void);
};

#endif