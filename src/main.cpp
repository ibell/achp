
#include "Compressor.h"
#include "BPHE.h"

int main()
{

	BrazedPlateHeatExchanger BPHE = BrazedPlateHeatExchanger();
	BPHE.test();
	BPHE.calculate();


	CompressorClass C = CompressorClass();
	C.test();
	C.calculate();
	std::vector<OutputEntryClass> list = C.OutputList();
	C.speed_test(1000);
	
	

	return 0;
}