
#include "Compressor.h"
#include "BPHE.h"
#include "time.h"
int main()
{
	

	//enable_TTSE_LUT("Nitrogen");
	set_standard_unit_system(UNIT_SYSTEM_SI);
	BrazedPlateHeatExchanger BPHE = BrazedPlateHeatExchanger();
	//BPHE.test();
	//for (double mdot = 0.0001; mdot < 0.1; mdot += 0.0001)
	//{
	//	BPHE.mdot_c = mdot;
	//	BPHE.SaturationStates();
	//	std::cout << mdot << " " << BPHE.calculate() << std::endl;
	//}
	
	BPHE = BrazedPlateHeatExchanger();

	BPHE.Nplates = 52;
	BPHE.plate_conductivity = 15.0; //[W/m-K]
	BPHE.more_channels = BPHE.MORE_CHANNELS_HOT;

	BPHE.geo.W = 0.200;
	BPHE.geo.Wp = 0.178;
	BPHE.geo.L = 0.4; 
	BPHE.geo.Lp = 0.455; // Center-to-center distance between ports
    BPHE.geo.PlateAmplitude = 0.002; //[m]
    BPHE.geo.PlateThickness = 0.0004; //[m]
    BPHE.geo.PlateWavelength =  0.007; //[m]
    BPHE.geo.InclinationAngle=  58.8/180.0*M_PI; //[rad]

	// Cold stream inlet
	BPHE.State_c_inlet = CoolPropStateClassSI("Water");
	//BPHE.State_c_inlet.update(iT,73.93+273.15-5.48,iQ,0.04854);
	BPHE.State_c_inlet.update(iP, 101325 ,iT, 300.316986084);
    BPHE.mdot_c = 0.262600004673; 
    
    // Hot stream inlet
    BPHE.State_h_inlet = CoolPropStateClassSI("Water");
	BPHE.State_h_inlet.update(iT, 311.979980469 ,iP, 101325);
	BPHE.mdot_h = 1.45700001717;

	BPHE.calculate();
	double rr = 0;

	///*
	//BPHE.mdot_h = BPHE.mdot_c*0.53940779782763737;
	//BPHE.calculate();
	//BPHE.verbosity = 10;
	//*/
	//
	//for (double mdot_ratio_h_c = 0.5; mdot_ratio_h_c < 0.9; mdot_ratio_h_c += 0.01)
	//{
	//	
	//	BPHE.mdot_h = BPHE.mdot_c*mdot_ratio_h_c;
	//	BPHE.calculate();
	//	std::cout << mdot_ratio_h_c << std::endl;
	//	//std::cout << mdot << " " << BPHE.calculate() << std::endl;
	//}

	CompressorClass C = CompressorClass();
	C.test();
	C.calculate();
	std::vector<OutputEntryClass> list = C.OutputList();
	//C.speed_test(1000);
	
	

	return 0;
}