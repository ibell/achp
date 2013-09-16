
#include "Compressor.h"
#include "BPHE.h"

int main()
{

	BrazedPlateHeatExchanger BPHE = BrazedPlateHeatExchanger();
	BPHE.test();
	for (double mdot = 0.001; mdot < 0.1; mdot += 0.001)
	{
		BPHE.mdot_c = mdot;
		BPHE.SaturationStates();
		std::cout << mdot << " " << BPHE.DetermineQmax() << std::endl;
	}

	CompressorClass C = CompressorClass();
	C.test();
	C.calculate();
	std::vector<OutputEntryClass> list = C.OutputList();
	//C.speed_test(1000);
	
	

	return 0;
}