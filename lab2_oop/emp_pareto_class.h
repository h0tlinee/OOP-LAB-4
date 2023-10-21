#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include "mixture_pareto_class.h"
#include "std_pareto_class.h"


class Empirical {
private:
	std::vector<double> selection;
	std::vector<double> empirical_density;
	int size = 0;
	int k = 0;
	double delta_calc();
	std::vector<double> create_intervals(const double delta);
	std::vector<double> create_empirical_density(const double delta, std::vector<double> intervals);
	double randomizer();
	std::vector<double> generate_x(int n);
	double qumulative_probability(int i);
	double top_bound();
public:
	Empirical(Pareto& std, int n, int _k);
	Empirical(mixture& mix, int n, int _k);
	Empirical(Empirical& ed, int n, int _k);
	Empirical(Empirical& emp);
	Empirical(std::vector<double>& selection);
	Empirical(std::string file);
	~Empirical();
	Empirical& operator = (const Empirical& emp);
	void SetK(int k);
	std::vector<double> GetSelection();
	std::vector<double> GetEmpDensity();
	int GetSize();
	int GetIntervalsNumber();
	double density(double x);
	double math_expect();
	double dispersion();
	double asymmetry();
	double excess();
	double rand_var();
	std::vector<std::pair<double, double>> generate_xy();
	void emp_output(std::ofstream& file);
};