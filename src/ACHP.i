%module ACHP

//// *************************** EXCEPTION HANDLING ****************************
//// *************************** EXCEPTION HANDLING ****************************
//// *************************** EXCEPTION HANDLING ****************************
%include exception.i

// A generic exception handler.  Any exceptions thrown in C++ will be caught here
%exception {
	try {
		$action
	}
    catch(std::exception &e) {
		SWIG_exception(SWIG_RuntimeError,e.what());
	}
    catch(...) {
		SWIG_exception(SWIG_RuntimeError,"Unknown exception");
	}
}

%ignore CoolPropStateClassSI::CoolPropStateClassSI(Fluid*);
%ignore CoolPropStateClassSI::CoolPropStateClassSI();

// This allows for the use of STL vectors
%include "std_vector.i"
// This allows for the use of STL strings
%include "std_string.i"

namespace std {
   %template(vectord) vector<double>;
};

// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "CoolProp/GlobalConstants.h"
#include "ACHPcore.h"
#include "CoolProp/CPState.h"
#include "BPHE.h"
#include "Compressor.h"
%}

// This is where the parsing actually happens
%include "CoolProp/GlobalConstants.h"
%include "ACHPcore.h"
%include "CoolProp/CPState.h"
%include "BPHE.h"
%include "Compressor.h"