%module Compressor


// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "ACHPcore.h"
#include "Compressor.h"
%}

// This allows for the use of STL vectors
%include "std_vector.i"

namespace std {
   %template(vectord) vector<double>;
};

// This is where the parsing actually happens
%include "ACHPcore.h"
%include "Compressor.h"