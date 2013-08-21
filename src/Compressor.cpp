
#include "CoolProp.h"
#include "Compressor.h"

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
        
    /*def OutputList(self):
        """
            Return a list of parameters for this component for further output
            
            It is a list of tuples, and each tuple is formed of items with indices:
                [0] Description of value
                
                [1] Units of value
                
                [2] The value itself
        """
        
        return [
            ('M1','-',self.M[0]),
            ('M2','-',self.M[1]),
            ('M3','-',self.M[2]),
            ('M4','-',self.M[3]),
            ('M5','-',self.M[4]),
            ('M6','-',self.M[5]),
            ('M7','-',self.M[6]),
            ('M8','-',self.M[7]),
            ('M9','-',self.M[8]),
            ('M10','-',self.M[9]),
            ('P1','-',self.P[0]),
            ('P2','-',self.P[1]),
            ('P3','-',self.P[2]),
            ('P4','-',self.P[3]),
            ('P5','-',self.P[4]),
            ('P6','-',self.P[5]),
            ('P7','-',self.P[6]),
            ('P8','-',self.P[7]),
            ('P9','-',self.P[8]),
            ('P10','-',self.P[9]),
            ('Heat Loss Fraction','-',self.fp),
            ('Displacement scale factor','-',self.Vdot_ratio),
            ('Power','W',self.W),
            ('Mass flow rate','kg/s',self.mdot_r),
            ('Inlet Temperature','K',self.Tin_r),
            ('Outlet Temperature','K',self.Tout_r),
            ('Inlet Enthalpy','J/kg',self.hin_r),
            ('Outlet Enthalpy','J/kg',self.hout_r),
            ('Overall isentropic efficiency','-',self.eta_oi),
            ('Pumped flow rate','m^3/s',self.Vdot_pumped)
         ]*/
        


void CompressorClass::Calculate()
{	
    long iFluid = get_Fluid_index(this->inlet_state.get_name());
    
    //  Calculate suction superheat and dew temperatures
    this->Tsat_s_K = this->inlet_state.Tsat(1.0);
    this->Tsat_d_K = this->outlet_state.Tsat(1.0);
    this->DT_sh_K = this->inlet_state.superheat();
    
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

    v_map = 1 / IProps(iD, iT, this->Tsat_s_K + 20.0/9.0*5.0, iP, P1, iFluid);
    v_actual = 1 / this->inlet_state._rho;
    F = 0.75;
    mdot = (1 + F * (v_map / v_actual - 1)) * mdot_map;

    double T1_map = this->Tsat_s_K + 20 * 5 / 9;
    double s1_map = IProps(iS, iT, T1_map, iP, P1, iFluid);
    double h1_map = IProps(iH, iT, T1_map, iP, P1, iFluid);
    double h2s_map = IProps(iH, iS, s1_map, iP, P2, iFluid);

    double s1_actual = this->inlet_state.s();
    double h1_actual = this->inlet_state.h();
    double h2s_actual = IProps(iH, iS, s1_actual, iP, P2, iFluid);

    //  Shaft power based on 20F superheat calculation from fit overall isentropic efficiency
    double power = power_map * (mdot / mdot_map) * (h2s_actual - h1_actual) / (h2s_map - h1_map);

    double h2 = power/1000 * (1 - this->fp) / mdot + h1_actual;
    this->eta_oi = mdot*(h2s_actual-h1_actual)/(power/1000);

	this->outlet_state = CoolPropStateClass(this->inlet_state.get_name());
	this->outlet_state.update(iH, h2, iP, P2);
    this->mdot = mdot;
    this->Wdot = power;
    this->CycleEnergyIn = power*(1 - this->fp);
    this->Vdot_pumped = mdot/IProps(iD, iT, this->inlet_state.T(), iP, P1, iFluid);
}
        
void CompressorClass::test()
{
	long iR134a = get_Fluid_index("R134a");
	double M[] = {217.3163128,5.094492028,-0.593170311,4.38E-02,-2.14E-02,1.04E-02,7.90E-05,-5.73E-05,1.79E-04,-8.08E-05};
	double P[] = {-561.3615705,-15.62601841,46.92506685,-0.217949552,0.435062616,-0.442400826,2.25E-04,2.37E-03,-3.32E-03,2.50E-03};
	this->M = std::vector<double>(M,M+sizeof(M)/sizeof(double));
	this->P = std::vector<double>(P,P+sizeof(P)/sizeof(double));
	this->fp = 0.15;
	this->Vdot_ratio = 1.0;

	this->inlet_state = CoolPropStateClass();
	this->outlet_state = CoolPropStateClass();

	Calculate();
}