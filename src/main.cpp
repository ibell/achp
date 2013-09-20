
#include "Compressor.h"
#include "BPHE.h"

int main()
{

	BrazedPlateHeatExchanger BPHE = BrazedPlateHeatExchanger();
	//BPHE.test();
	//for (double mdot = 0.0001; mdot < 0.1; mdot += 0.0001)
	//{
	//	BPHE.mdot_c = mdot;
	//	BPHE.SaturationStates();
	//	std::cout << mdot << " " << BPHE.calculate() << std::endl;
	//}

	BPHE = BrazedPlateHeatExchanger();

	// Cold stream inlet
	BPHE.State_c_inlet = CoolPropStateClassSI("Propane");
    BPHE.State_c_inlet.update(iT,290,iQ,0.3);
    BPHE.mdot_c = 0.01;
    
    // Hot stream inlet
    BPHE.State_h_inlet = CoolPropStateClassSI("Propane");
    BPHE.State_h_inlet.update(iT,320,iQ,1);
    BPHE.State_h_inlet.update(iT,360,iP,BPHE.State_h_inlet.p());

	BPHE.mdot_h = BPHE.mdot_c*0.53940779782763737;
	BPHE.calculate();
	
	//for (double mdot_ratio_h_c = 0.001; mdot_ratio_h_c < 1000; mdot_ratio_h_c *= 1.1)
	//{
	//	
	//	BPHE.mdot_h = BPHE.mdot_c*0.53940779782763737;
	//	BPHE.calculate();
	//	//std::cout << mdot << " " << BPHE.calculate() << std::endl;
	//}

	CompressorClass C = CompressorClass();
	C.test();
	C.calculate();
	std::vector<OutputEntryClass> list = C.OutputList();
	//C.speed_test(1000);
	
	

	return 0;
}