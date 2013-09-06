#ifndef ACHPCORE_H
#define ACHPCORE_H

#include <iostream>
#include <string>
#include <vector>

class OutputEntryClass
{

private:
	std::string variable, units, string, value_string;
	double value_double;
	int type;
	enum types{TYPE_DOUBLE, TYPE_STRING};

public:
	OutputEntryClass(std::string variable, std::string units, std::string value)
	{
		this->variable = variable; this->units = units; this->value_string = value; this->type = TYPE_STRING;
	};
	OutputEntryClass(std::string variable, std::string units, double value)
	{
		this->variable = variable; this->units = units; this->value_double = value; this->type = TYPE_DOUBLE;
	}
};

class ACHPComponentClass
{
	virtual std::vector<OutputEntryClass> OutputList() = 0;
};

#endif