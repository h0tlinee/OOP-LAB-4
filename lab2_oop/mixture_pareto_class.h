#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include "std_pareto_class.h"


const double PI = 3.141592653589793;


class mixture {
private:
	double p;
	Pareto pareto1;
	Pareto pareto2;

public:
	mixture();
	mixture(Pareto& std1, Pareto& std2, double p);
};


