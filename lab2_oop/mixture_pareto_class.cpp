
#include "std_pareto_class.h"

//открытые функции

template<class distrib1,class distrib2>
class mixture:public IDistribution,public IPersistend{
public:
	mixture(distrib1& d1, distrib2& d2, double p) ://конструктор с использованием инициализаторов
		d1(d1), d2(d2), p(p) {};
	mixture(std::ifstream& input) {
		this->load_attr(input);
	};
	double density(double x) const override;
	double math_expect() const override;
	double dispersion() const override;
	double excess() const override;
	double asymmetry() const override;
	double modeling() const override;
	void xy_output(int n);
	
	std::vector<double> generate_x(int n) const override;
	std::vector<std::pair<double, double>>generate_xy(int n) const override;

	void load_attr(std::ifstream& input) override;
	void save_attr(std::ofstream& output) override;

	distrib1& component1() {
		return d1;
	}
	distrib2& component2() {
		return d2;
	}

	void SetP(double p);
	double GetP() const;


private:
	distrib1 d1;
	distrib2 d2;
	double p;
	double randomizer() const;

};

template<class distr1,class distr2>
double mixture<distr1, distr2>::density(double x) const {
	return (1 - p) * d1.density(x) + p * d2.density(x);
}

template<class distr1, class distr2>
double mixture<distr1, distr2>::math_expect() const {
	return (1 - p) * d1.math_expect() + p * d2.math_expect();
}

template<class distr1, class distr2>
double mixture<distr1, distr2>::dispersion() const {
	return (1 - p) * (pow(d1.math_expect(), 2) + d1.dispersion()) +
		p * (pow(d2.math_expect(), 2) + d2.dispersion()) -
		pow(math_expect(), 2);
}

template<class distr1, class distr2>
double mixture<distr1, distr2>::excess() const {
	return ((1 - p) * (pow((d1.math_expect() - math_expect()), 4) + 6 * d1.dispersion() * pow((d1.math_expect() - math_expect()), 2) +
		4 * (d1.math_expect() - math_expect()) * pow(d1.dispersion(), 3 / 2) * d1.asymmetry() + pow(d1.dispersion(), 2) * d1.excess()) +
		p * (pow((d2.math_expect() - math_expect()), 4) + 6 * d2.dispersion() * pow((d2.math_expect() - math_expect()), 2) +
			4 * (d2.math_expect() - math_expect()) * pow(d2.dispersion(), 3 / 2) * d2.asymmetry() + pow(d2.dispersion(), 2) * d2.excess()) - 3) /
		pow(dispersion(), 2);
}

template<class distr1, class distr2>
double mixture<distr1, distr2>::asymmetry() const {
	return ((1 - p) * (pow((d1.math_expect() - math_expect()), 3) + 3 * (d1.math_expect() - math_expect()) * d1.dispersion() + pow(d1.dispersion(), 3 / 2) * d1.asymmetry()) +
		p * (pow((d2.math_expect() - math_expect()), 3) + 3 * (d2.math_expect() - math_expect()) * d2.dispersion() + pow(d2.dispersion(), 3 / 2) * d2.asymmetry())) /
		pow(dispersion(), 3 / 2);
}

template<class distr1, class distr2>
double mixture<distr1, distr2>::randomizer() const {
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

template<class distr1, class distr2>
double mixture<distr1, distr2>::modeling() const {
	double r = randomizer();
		if (r > p) {
			return d1.modeling();
		}
		else {
			return d2.modeling();
		}
}

template<class distr1, class distr2>
std::vector<double> mixture<distr1, distr2>::generate_x(int n) const {
	std::vector<double> x;
		for (int i = 0; i < n; ++i) {
			x.push_back(modeling());
		}
		sort(x.begin(), x.end());
		return x;
}

template<class distr1, class distr2>
std::vector<std::pair<double, double>> mixture<distr1,distr2>::generate_xy(int n) const {
	std::vector<std::pair<double, double>> result;
		auto x = this->generate_x(n);
		for (int i = 0; i < n; i++) {
			result.push_back(std::make_pair(x[i], this->density(x[i])));
		}
		return result;
}

template<class distr1, class distr2>
void mixture<distr1, distr2>::save_attr(std::ofstream& output)  {
	component1().save_attr(output);
	component2().save_attr(output);
	output << p<<std::endl;
}

template<class distr1,class distr2>
void mixture<distr1, distr2>::load_attr(std::ifstream& input)  {
	component1().load_attr(input);
	component2().load_attr(input);
	input >> p;
	
}


template<class distr1, class distr2>
double mixture<distr1, distr2>::GetP() const {
	return this->p;
}

template<class distr1,class distr2>
void mixture<distr1, distr2>::SetP(double _p) {
	if (_p < 0 || _p>1) {
		throw 1;
	}
	else {
		this->p = _p;
	}
}

template<class distr1,class distr2>
void mixture<distr1, distr2>::xy_output(int n) {
	auto result = this->generate_xy(n);
	std::ofstream file("mixture_xy_output.txt");
	for (int i = 0; i < n; i++) {
		file << result[i].first << '\t' << result[i].second << std::endl;
	}
	file.close();
}

