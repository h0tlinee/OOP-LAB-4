#include "emp_pareto_class.h"
#include <fstream>


Empirical::Empirical(const IDistribution& distr, int n, int k) ://конструктор
	size(n > 1 ? n : throw 1), k(k > 2 ? k : (int)log2(size) + 1) {
	selection = distr.generate_x(size);
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

Empirical::Empirical(const Empirical& emp):
	size(emp.size > 1 ? emp.size : throw 1), k(emp.k > 2 ? emp.k : ((int)log2(size) + 1)), selection(emp.selection), empirical_density(emp.empirical_density){//выборка на основе выборки
}

Empirical::Empirical(std::string file) {
	std::ifstream input;
	input.open("emp_selection.txt");
	if (!input.is_open()) {
		throw 0;
	}
	double x;
	while (input >> x) {
		selection.push_back(x);
	}
	size = selection.size();
	k = (int)log2(size) + 1;
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

Empirical::Empirical(const std::vector<double>& selection) :
	size(selection.size()), k((int)log2(size) + 1), selection(selection) {
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

Empirical::~Empirical() {
	selection.clear();
	empirical_density.clear();
}

std::vector<double> Empirical::create_empirical_density(const double delta, std::vector<double> intervals) const {
		std::vector<double> result;
		int j = 0;
		for (int i = 0; i < intervals.size() - 1; ++i) {
			int count = 0;
			for (j; j < size; ++j) {
				if (i + 1 == intervals.size() - 1) {
					if (selection[j] <= intervals[i + 1]) {
						++count;
					}
					else {
						break;
					}
				}
				else {
					if (selection[j] < intervals[i + 1]) {
						++count;
					}
					else {
						break;
					}
				}
			}
			result.push_back(count / (size * delta));
		}
		return result;
}

double Empirical::randomizer() const {
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

double Empirical::qumulative_probability(const int i) const {
	double q = 0;
	for (int j = 0; j <= i; ++j) {
		q += empirical_density[j];
	}
	return q;
}

double Empirical::top_bound() const {
	double sum = 0;
	for (int i = 0; i < k; ++i) {
		sum += empirical_density[i];
	}
	return sum;
}

double Empirical::modeling() const {
	double r;
	do r = (double)rand() / RAND_MAX * top_bound();
	while (r == 0 || r == top_bound());
	for (int i = 0; i < k - 1; ++i) {
		if (r > qumulative_probability(i) and r < qumulative_probability(i + 1)) {
			do r = (double)rand() / RAND_MAX * ((selection[0] + delta_calc() * (i + 1)) - (selection[0] + delta_calc() * i)) + (selection[0] + delta_calc() * (i + 1));
			while (r == (selection[0] + delta_calc() * i) || r == (selection[0] + delta_calc() * (i + 1)));
			break;
		}
	}
	return r;
}

Empirical& Empirical::operator = (const Empirical& emp) {
	if (this == &emp) {
		return *this;
	}
	selection = emp.selection;
	empirical_density = emp.empirical_density;
	size = emp.size;
	k = emp.k;
	return *this;
}

void Empirical::SetK(const int k) {
	if (k >= 2) {
		this->k = k;
	}
	else {
		this->k = (int)log2(size) + 1;
	}
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

std::vector<double> Empirical::GetSelection() const {
	return selection;
}

std::vector<double> Empirical::GetEmpDensity() const {
	return empirical_density;
}

int Empirical::GetSize() const {
	return size;
}

int Empirical::GetIntervalsNumber() const {
	return k;
}

double Empirical::delta_calc() const {
	return (selection[size - 1] - selection[0]) / (k);
}

std::vector<double> Empirical::create_intervals(const double delta) const {
	std::vector<double> intervals;
	double slider = selection[0];
	for (int i = 0; i < k + 1; ++i) {
		intervals.push_back(slider);
		slider += delta;
	}
	return intervals;
}

double Empirical::density(const double x) const {
	double slider = selection[0];
	if (x < slider) {
		return 0;
	}
	for (int i = 0; i < k; ++i) {
		slider += delta_calc();
		if (i == k - 1) {
			if (x <= slider) {
				return empirical_density[i];
			}
		}
		else {
			if (x < slider) {
				return empirical_density[i];
			}
		}
	}
	return 0;
}

double Empirical::math_expect() const {
	double sum = 0;
	for (int i = 0; i < size; ++i) {
		sum += selection[i];
	}
	return sum / size;
}

double Empirical::dispersion() const {
	double sum = 0;
	double exp_val = math_expect();
	for (int i = 0; i < size; ++i) {
		sum += pow(selection[i] - exp_val, 2);
	}
	return sum / size;
}

double Empirical::asymmetry() const {
	double sum = 0;
	double exp_val = math_expect();
	for (int i = 0; i < size; ++i) {
		sum += pow(selection[i] - exp_val, 3);
	}
	return sum / (size * pow(dispersion(), 3 / 2));
}

double Empirical::excess() const {
	double sum = 0;
	double exp_val = math_expect();
	for (int i = 0; i < size; ++i) {
		sum += pow(selection[i] - exp_val, 4);
	}
	return (sum / (size * pow(dispersion(), 2))) - 3;
}

std::vector<double>  Empirical::generate_x(const int n) const {
	std::vector<double> result;
	for (int i = 0; i < n; ++i) {
		result.push_back(modeling());
	}
	sort(result.begin(), result.end());
	return result;
}

std::vector<std::pair<double, double>> Empirical::generate_xy(int n) const {
	std::vector<double> selection;
	for (int i = 0; i < n; i++) {
		selection.push_back(modeling());
	}
	std::vector <std::pair<double, double >> result;
	for (int j = 0; j < size; ++j) {
		result.push_back(std::make_pair(selection[j], density(selection[j])));
	}
	return result;
}

std::vector<std::pair<double, double>> Empirical::generate_graph_selection(const std::vector<double>& selection) const {
	std::vector <std::pair<double, double >> result;
	for (int j = 0; j < size; ++j) {
		result.push_back(std::make_pair(selection[j], density(selection[j])));
	}
	return result;
}

void Empirical::save_attr(std::ofstream& output) {
	if (!output.is_open()) {
		throw 0;
	}
	output.open("emp_selection.txt");
	for (int i = 0; i < size; ++i) {
		output << selection[i] <<std::endl;
	}
	output.close();
}

void Empirical::load_attr(std::ifstream& input) {
	input.open("emp_selection.txt");
	if (!input.is_open()) {
		throw 0;
	}
	double x;
	std::vector<double> result;
	while (input >> x) {
		result.push_back(x);
	}
	selection = result;
	size = result.size();
	k = (int)log2(size) + 1;
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
	
}

void Empirical::emp_output(int n){//вывод в файл
	std::ofstream file;
	auto pairs = generate_xy(n);
	file.open("emp_output.txt");
	for (int i = 0; i < size; ++i) {
		file << pairs[i].first << "\t" << pairs[i].second << std::endl;
	}
	file.close();
}



