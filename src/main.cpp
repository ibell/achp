
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
	BPHE.State_h_inlet = CoolPropStateClassSI("Propane");
	BPHE.State_h_inlet.update(1,-1,2,101325);
	BPHE.test();
	for (double mdot = 0.0001; mdot < 0.1; mdot += 0.0001)
	{
		BPHE.mdot_c = mdot;
		BPHE.calculate();
		//std::cout << mdot << " " << BPHE.calculate() << std::endl;
	}

	CompressorClass C = CompressorClass();
	C.test();
	C.calculate();
	std::vector<OutputEntryClass> list = C.OutputList();
	//C.speed_test(1000);
	
	

	return 0;
}