
#include "Compressor.h"

int main()
{
	CompressorClass C = CompressorClass();

	C.test();
	C.calculate();
	std::vector<OutputEntryClass> list = C.OutputList();
	C.speed_test(1000);
	
	return 0;
}