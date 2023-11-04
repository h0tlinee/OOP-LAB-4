#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include "std_pareto_class.h"




class Empirical :public IDistribution, public IPersistend {
public:
	Empirical(const IDistribution& distr, int n, int k);
	Empirical(const Empirical& emp);
	Empirical(std::string file);
	Empirical(const std::vector<double>& selection);
	~Empirical();

	double density(double x) const override;
	double math_expect() const override;
	double dispersion() const override;
	double excess() const override;
	double asymmetry() const override;
	double modeling() const override;

	std::vector<double> generate_x(int n) const override;
	std::vector<std::pair<double, double>> generate_xy(int n) const override;
	std::vector<std::pair<double, double>> generate_graph_selection(const std::vector<double>& selection) const;

	Empirical& operator=(const Empirical& emp);

	void SetK(int k);

	std::vector<double> GetSelection() const;
	std::vector<double> GetEmpDensity() const;

	int GetSize() const;
	int GetIntervalsNumber() const;

	void save_attr(std::ofstream& output) override;
	void load_attr(std::ifstream& input) override;

	void emp_output(int n);

private:
	std::vector<double> selection;
	std::vector<double> empirical_density;
	int size;
	int k;

	double delta_calc() const;
	std::vector<double> create_intervals(const double delta) const;
	double randomizer() const;
	std::vector<double> create_empirical_density(const double delta, std::vector<double> intervals) const;
	double qumulative_probability(const int i) const;
	double top_bound() const;
};