
#include "CoolProp/CoolProp.h"
#include "Compressor.h"
#include "time.h"
#include "ACHPcore.h"

/*!

Compressor Model based on 10-coefficient Model from `ANSI/AHRI standard 540 <http://www.ahrinet.org/App_Content/ahri/files/standards%20pdfs/ANSI%20standards%20pdfs/ANSI-ARI-540-2004%20latest.pdf>`_

Required Parameters:
    
===========   ==========  ========================================================================
Variable      Units       Description
===========   ==========  ========================================================================
M             varied      A numpy-like list of compressor map coefficients for mass flow
P             varied      A numpy-like list of compressor map coefficients for electrical power
inlet_state
outlet_state
fp            --          Fraction of electrical power lost as heat to ambient
Vdot_ratio    --          Displacement Scale factor
===========   ==========  ========================================================================

All variables are of double-type unless otherwise specified

*/

std::vector<OutputEntryClass> CompressorClass::OutputList()
{
	std::vector<OutputEntryClass> list;

	list.push_back(OutputEntryClass("M1","-",this->M[0]));
    list.push_back(OutputEntryClass("M2","-",this->M[1]));
    list.push_back(OutputEntryClass("M3","-",this->M[2]));
	list.push_back(OutputEntryClass("M4","-",this->M[3]));
	list.push_back(OutputEntryClass("M5","-",this->M[4]));
	list.push_back(OutputEntryClass("M6","-",this->M[5]));
	list.push_back(OutputEntryClass("M7","-",this->M[6]));
	list.push_back(OutputEntryClass("M8","-",this->M[7]));
	list.push_back(OutputEntryClass("M9","-",this->M[8]));
	list.push_back(OutputEntryClass("M10","-",this->M[9]));
	list.push_back(OutputEntryClass("P1","-",this->P[0]));
	list.push_back(OutputEntryClass("P2","-",this->P[1]));
	list.push_back(OutputEntryClass("P3","-",this->P[2]));
	list.push_back(OutputEntryClass("P4","-",this->P[3]));
	list.push_back(OutputEntryClass("P5","-",this->P[4]));
	list.push_back(OutputEntryClass("P6","-",this->P[5]));
	list.push_back(OutputEntryClass("P7","-",this->P[6]));
	list.push_back(OutputEntryClass("P8","-",this->P[7]));
	list.push_back(OutputEntryClass("P9","-",this->P[8]));
	list.push_back(OutputEntryClass("P10","-",this->P[9]));
	list.push_back(OutputEntryClass("Heat Loss Fraction","-",this->fp));
	list.push_back(OutputEntryClass("Displacement scale factor","-",this->Vdot_ratio));
	list.push_back(OutputEntryClass("Power","W",this->Wdot));
	list.push_back(OutputEntryClass("Mass flow rate","kg/s",this->mdot));
	list.push_back(OutputEntryClass("Inlet Temperature","K",this->inlet_state.T()));
	list.push_back(OutputEntryClass("Outlet Temperature","K",this->outlet_state.T()));
	list.push_back(OutputEntryClass("Inlet Enthalpy","J/kg",this->inlet_state.h()));
	list.push_back(OutputEntryClass("Outlet Enthalpy","J/kg",this->outlet_state.h()));
	list.push_back(OutputEntryClass("Overall isentropic efficiency","-",this->eta_oi));
	list.push_back(OutputEntryClass("Pumped flow rate","m^3/s",this->Vdot_pumped));
	return list;
}
void CompressorClass::calculate()
{	
	//_validate();
	long iFluid = get_Fluid_index("R134a");
    
    //  Calculate suction superheat and dew temperatures
	this->Tsat_s_K = IProps(iT, iP, this->inlet_state.p(), iQ, 1, iFluid);
	this->Tsat_d_K = IProps(iT, iP, this->outlet_state.p(), iQ, 1, iFluid);
	this->DT_sh_K = this->inlet_state.T()-this->Tsat_s_K;

	// Bug in CoolPropStateClass
    //this->Tsat_s_K = this->inlet_state.Tsat(1.0);
    //this->Tsat_d_K = this->outlet_state.Tsat(1.0);
    //this->DT_sh_K = this->inlet_state.superheat();
    
    //  Convert saturation temperatures in K to F
    Tsat_s = this->Tsat_s_K * 9/5 - 456.67;
    Tsat_d = this->Tsat_d_K * 9/5 - 456.67;

    //  Apply the 10 coefficient ARI map to saturation temps in F
    double power_map = this->P[0] + this->P[1] * Tsat_s + this->P[2] * Tsat_d + this->P[3] * Tsat_s*Tsat_s + this->P[4] * Tsat_s * Tsat_d + this->P[5] * Tsat_d*Tsat_d + this->P[6] * Tsat_s*Tsat_s*Tsat_s + this->P[7] * Tsat_d * Tsat_s * Tsat_s + this->P[8] * Tsat_d * Tsat_d * Tsat_s + this->P[9] * Tsat_d*Tsat_d*Tsat_d;
    double mdot_map = this->M[0] + this->M[1] * Tsat_s + this->M[2] * Tsat_d + this->M[3] * Tsat_s*Tsat_s + this->M[4] * Tsat_s * Tsat_d + this->M[5] * Tsat_d*Tsat_d + this->M[6] * Tsat_s*Tsat_s*Tsat_s + this->M[7] * Tsat_d * Tsat_s * Tsat_s + this->M[8] * Tsat_d * Tsat_d * Tsat_s + this->M[9] * Tsat_d*Tsat_d*Tsat_d;

    //  Convert mass flow rate to kg/s from lbm/h
    mdot_map *= 0.000125998 ;

    //  Add more mass flow rate to scale
    mdot_map *= this->Vdot_ratio;
    power_map *= this->Vdot_ratio;

    P1 = this->inlet_state.p();
    P2 = this->outlet_state.p();
    T1_actual = this->Tsat_s_K + this->DT_sh_K;

	double s1_actual = this->inlet_state.s();
    double h1_actual = this->inlet_state.h();
    double h2s_actual = IProps(iH, iS, s1_actual, iP, P2, iFluid);

	// Map correction
    v_map = 1 / IProps(iD, iT, this->Tsat_s_K + 20.0/9.0*5.0, iP, P1, iFluid);
    v_actual = 1 / this->inlet_state._rho;
    F = 0.75;
    mdot = (1 + F * (v_map / v_actual - 1)) * mdot_map;

    double T1_map = this->Tsat_s_K + 20 * 5 / 9;
    double s1_map = IProps(iS, iT, T1_map, iP, P1, iFluid);
    double h1_map = IProps(iH, iT, T1_map, iP, P1, iFluid);
    double h2s_map = IProps(iH, iS, s1_map, iP, P2, iFluid);

    //  Shaft power based on 20F superheat calculation from fit overall isentropic efficiency
    double power = power_map * (mdot / mdot_map) * (h2s_actual - h1_actual) / (h2s_map - h1_map);

    double h2 = power/1000 * (1 - this->fp) / mdot + h1_actual;
    this->eta_oi = mdot*(h2s_actual-h1_actual)/(power/1000);

	this->outlet_state.update(iH, h2, iP, P2);
    this->mdot = mdot;
    this->Wdot = power;
    this->CycleEnergyIn = power*(1 - this->fp);
	this->Vdot_pumped = mdot/this->inlet_state.rho();
}
void CompressorClass::speed_test(unsigned long long N)
{
	clock_t t1, t2;
	t1 = clock();
	for (unsigned long long i = 0; i < N; i++)
	{
		this->calculate();
	}
	t2 = clock();
	std::cout << format("Running the compressor took %g us/call",(double)(t2-t1)/CLOCKS_PER_SEC/(double)N*1e6).c_str(); 
}
void CompressorClass::test()
{
	double M[] = {217.3163128,5.094492028,-0.593170311,4.38E-02,-2.14E-02,1.04E-02,7.90E-05,-5.73E-05,1.79E-04,-8.08E-05};
	double P[] = {-561.3615705,-15.62601841,46.92506685,-0.217949552,0.435062616,-0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03};
	this->M = std::vector<double>(M,M+sizeof(M)/sizeof(double));
	this->P = std::vector<double>(P,P+sizeof(P)/sizeof(double));
	this->fp = 0.15;
	this->Vdot_ratio = 1.0;

	this->inlet_state = CoolPropStateClass("R134a");
	this->inlet_state.update(iT,300,iQ,1);
	std::string name = this->inlet_state.get_name();
	this->outlet_state = CoolPropStateClass("R134a");
	this->outlet_state.update(iT,350,iQ,1);
}