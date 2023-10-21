#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include "std_pareto_class.h"




class mixture {
private:
	double p;
	Pareto pareto1;
	Pareto pareto2;
	double randomizer();

public:
	mixture();
	mixture(Pareto& std1, Pareto& std2, double p);
	mixture(std::string file);
	void SetP(double p);
	double GetP();

	Pareto& component1();
	Pareto& component2();

	void save_atr(std::string file);
	double density(double x);
	double math_expect();
	double dispersion();
	double asymmetry();
	double excess();
	double modeling();
	std::vector<double> generate_x(int n);
	std::vector<std::pair<double,double>> generate_xy(int n);
	void xy_output(int n);
};


