#ifndef ACHP_BPHE_H
#define ACHP_BPHE_H

#include "coolprop/CoolProp/CPState.h"
#include "thermalcorr/src/BrazedPlateHeatExchanger.h"
#include "thermalcorr/src/InternalFlow.h"

namespace TCBPHE = BrazedPlateHX;

class BPHECell
{
public:
	BPHECell()
	{
		// Initialize all values to huge
		charge_h = _HUGE; charge_c = _HUGE; HTC_h = _HUGE; HTC_c = _HUGE; T_i_h = _HUGE; T_i_c = _HUGE;
		T_o_h = _HUGE; T_o_c = _HUGE; Tmean_h = _HUGE; Tmean_c = _HUGE; UA_required = _HUGE; UA_available = _HUGE;
		w = _HUGE; DP_h = _HUGE; DP_c = _HUGE; Qdot = _HUGE; LMTD = _HUGE; 
		
		// Initialize all CoolPropStateClasses to NULL pointers
		CPS_h = NULL; CPS_c = NULL;
	};
	CoolPropStateClassSI *CPS_h, *CPS_c; ///< A pointer to a CoolPropStateClassSI instance
	double Qdot;     ///< The heat transfer rate [W]
	double charge_h; ///< The mass of hot fluid [kg]
	double charge_c; ///< The mass of cold fluid [kg]
	double HTC_h;    ///< The heat transfer coefficient of the hot fluid [W/m^2/K]
	double HTC_c;    ///< The heat transfer coefficient of the cold fluid [W/m^2/K]
	double T_i_h;  ///< The inlet temperature on the hot side [K]
	double T_i_c;  ///< The inlet temperature on the cold side [K]
	double T_o_h;  ///< The outlet temperature on the hot side [K]
	double T_o_c;  ///< The outlet temperature on the cold side [K]
	double Tmean_h;  ///< The mean temperature on the hot side [K]
	double Tmean_c;  ///< The mean temperature on the cold side [K]
	double LMTD;  ///< The log mean temperature difference [K]
	double UA_required;   ///< The required heat transfer conductance [W/K]
	double UA_available; ///< The available heat transfer conductance if all the BPHE is in this cell [W/K]
	double w;        ///< The fraction of the length of the HX in this cell [-]
	double DP_h;     ///< The pressure change on the hot side [Pa]
	double DP_c;     ///< The pressure change on the cold side [Pa]
};

class BrazedPlateHeatExchanger
{

public:
	enum {MORE_CHANNELS_HOT, MORE_CHANNELS_COLD};

	/// 
	BrazedPlateHeatExchanger(){this->verbosity = 1; plate_conductivity = _HUGE; more_channels = -1; Cooper_tune_factor = 1.0;};

	/// Inlet state for the hot stream
	CoolPropStateClassSI State_h_inlet;
		
	/// Inlet state for the cold stream
	CoolPropStateClassSI State_c_inlet;

	/// Saturation state for the cold stream at the inlet pressure
	CoolPropStateClassSI State_h_sat;

	/// Saturation state for the hot stream at the inlet pressure
	CoolPropStateClassSI State_c_sat;

	double mdot_c; ///< Mass flow rate of cold stream [kg/s]
	double mdot_h; ///< Mass flow rate of hot stream [kg/s]
	std::vector<double> EnthalpyList_c; ///< The list of enthalpies for the cold stream at cell boundaries
	std::vector<double> EnthalpyList_h; ///< The list of enthalpies for the hot stream at cell boundaries
	std::vector<double> TemperatureList_c; ///< The list of temperatures of the cold stream at cell boundaries
	std::vector<double> TemperatureList_h; ///< The list of temperatures of the hot stream at cell boundaries
	std::vector<double> DensityList_c; ///< The list of densities of the cold stream at cell boundaries
	std::vector<double> DensityList_h; ///< The list of densities of the hot stream at cell boundaries
	std::vector<bool> PhaseBoundary_h; ///< A vector of booleans for whether cell boundary index corresponds to phase change of the hot stream
	std::vector<bool> PhaseBoundary_c; ///< A vector of booleans for whether cell boundary index corresponds to phase change of the hot stream
	std::vector<int> CellPhaseList_c; ///< A vector of phase indices for the cells on the cold side
	std::vector<int> CellPhaseList_h; ///< A vector of phase indices for the cells on the hot side
	std::vector<CoolPropStateClassSI> CellStateList_h; ///< A vector of state class instances for the cells on the hot side
	std::vector<CoolPropStateClassSI> CellStateList_c; ///< A vector of state class instances for the cells on the cold side
	std::vector<BPHECell> CellList; ///< A vector of cell instances
	/// Build the list of enthalpies for a given heat transfer rate
	void BuildEnthalpyLists(double Q);

	/// Calculate the saturation states for both streams
	void SaturationStates();

	/// Calculate the actual heat transfer rate
	void CalculateQ(double Qmax);
	
	int verbosity; ///< How verbose the debugging should be [0: no output, 10: very annoying output]

	void _check();
	void test();
	void calculate();
	double DetermineQmax(void);
	double Cooper_tune_factor; ///< The tuning factor for Cooper pool boiling relationship

	double Qdot; ///< The heat transfer rate [kW]
	double Qmax; ///< The maximum heat transfer rate [kW]
	double DELTAP_h; ///< Total change in pressure on hot side [Pa]
	double DELTAP_c; ///< Total change in pressure on cold side [Pa]

	int Ngaps_h; ///< The number of gaps on the hot side [-]
	int Ngaps_c; ///< The number of gaps on the cold side [-]
	int Nplates; ///< The number of plates [-]
	int more_channels; ///< Which side gets the additional channel if uneven, one of MORE_CHANNELS_HOT or MORE_CHANNELS_COLD

	double A_wetted_h; ///< The wetted area on the hot side that contributes to heat transfer [m^2]
	double V_h; ///< The volume of the hot side [m^3]
	double A_flow_h; ///< The total flow area of all of the cold gaps [m^2]
	double Dh_h; ///< Hydraulic diameter of the hot side gaps [m]
	double A_wetted_c; ///< The wetted area on the cold side that contributes to heat transfer [m^2]
	double V_c; ///< The volume of the cold side [m^3]
	double A_flow_c; ///< The total flow area of all of the hot gaps [m^2]
	double Dh_c; ///< Hydraulic diameter of the cold side gaps [m]
	double R_plate; ///< The thermal resistance of the plate in [K/W]
	double plate_conductivity; ///< The thermal conductivity of the plate [W/m/K]
	
	TCBPHE::BPHEGeometry geo;

	void _OnePhaseH_OnePhaseC_Qimposed(BPHECell *);
	//void _TwoPhaseH_OnePhaseC_Qimposed(BPHECell *);
	void _OnePhaseH_TwoPhaseC_Qimposed(BPHECell *);
};

#endif