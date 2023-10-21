#include "mixture_pareto_class.h"

//открытые функции


void mixture::SetP(double _p) {
	if (_p < 0 || _p > 1)
	{
		throw 1;
	}
	else {
		this->p = _p;
	}
}

double mixture::GetP() {
	return this->p;
}


Pareto& mixture::component1() {
	return pareto1;
}

Pareto& mixture::component2() {
	return pareto2;
}

void mixture::save_atr(std::string file) {//сохранение атрибутов в файл
	std::ofstream output;
	output.open(file);
	std::vector<double> buffer;
	buffer.push_back(this->component1().GetForm());
	buffer.push_back(this->component1().GetMu());
	buffer.push_back(this->component1().GetLambda());
	buffer.push_back(this->component2().GetForm());
	buffer.push_back(this->component2().GetMu());
	buffer.push_back(this->component2().GetLambda());
	buffer.push_back(this->GetP());
	if (!output) {
		throw std::runtime_error("Невозможно найти файл");
	};
	for (int i = 0; i < 7; i++) {
		output << buffer[i] << " ";
	}
	output.close();

}

double mixture::density(double x) {
	return (1 - p) * pareto1.pareto_std_distribution(x) + p * pareto2.pareto_std_distribution(x);
}

double mixture::math_expect() {
	return (1 - p) * pareto1.math_expect() + p * pareto2.math_expect();
}

double mixture::dispersion() {
	return (1 - p) * (pow(pareto1.math_expect(), 2) + pareto1.dispersion()) +
		p * (pow(pareto2.math_expect(), 2) + pareto2.dispersion()) -
		pow(math_expect(), 2);
}

double mixture::asymmetry() {
	return ((1 - p) * (pow((pareto1.math_expect() - math_expect()), 3) + 3 * (pareto1.math_expect() - math_expect()) * pareto1.dispersion() + pow(pareto1.dispersion(), 3 / 2) * pareto1.asymmetry()) +
		p * (pow((pareto2.math_expect() - math_expect()), 3) + 3 * (pareto2.math_expect() - math_expect()) * pareto2.dispersion() + pow(pareto2.dispersion(), 3 / 2) * pareto2.asymmetry())) /
		pow(dispersion(), 3 / 2);
}

double mixture::excess() {
	return ((1 - p) * (pow((pareto1.math_expect() - math_expect()), 4) + 6 * pareto1.dispersion() * pow((pareto1.math_expect() - math_expect()), 2) +
		4 * (pareto1.math_expect() - math_expect()) * pow(pareto1.dispersion(), 3 / 2) * pareto1.asymmetry() + pow(pareto1.dispersion(), 2) * pareto1.excess_koef()) +
		p * (pow((pareto2.math_expect() - math_expect()), 4) + 6 * pareto2.dispersion() * pow((pareto2.math_expect() - math_expect()), 2) +
			4 * (pareto2.math_expect() - math_expect()) * pow(pareto2.dispersion(), 3 / 2) * pareto2.asymmetry() + pow(pareto2.dispersion(), 2) * pareto2.excess_koef()) - 3) /
		pow(dispersion(), 2);
}



double mixture::modeling() {
	double r = randomizer();
	if (r > p) {
		return pareto1.pareto_modeling();
	}
	else {
		return pareto2.pareto_modeling();
	}
}

std::vector<double> mixture::generate_x(int n) {
	std::vector<double> x;
	for (int i = 0; i < n; ++i) {
		x.push_back(modeling());
	}
	sort(x.begin(), x.end());
	return x;
}

std::vector<std::pair<double,double>> mixture::generate_xy(int n) {
	std::vector<std::pair<double, double>> result;
	auto x = this->generate_x(n);
	for (int i = 0; i < n; i++) {
		result.push_back(std::make_pair(x[i], this->density(x[i])));
	}
	return result;
}

void mixture::xy_output(int n) {
	auto result = this->generate_xy(n);
	std::ofstream file("mixture_xy_output.txt");
	for (int i = 0; i < n; i++) {
		file << result[i].first << '\t' << result[i].second << std::endl;
	}
	file.close();
}




//конструкторы


mixture::mixture():p(0.5) {//конструктор по умолчанию
	pareto1.SetForm(3.3);
	pareto2.SetForm(3.3);
	pareto1.SetLambda(1);
	pareto2.SetLambda(2);
	pareto1.SetMu(-2);
	pareto2.SetMu(2);
	

}

mixture::mixture(Pareto& std1, Pareto& std2, double _p)://конструктор из существующих двух стандратных распределений
	pareto1(std1), pareto2(std2),p(_p>0 and _p<1 ? _p:throw 1){}




mixture::mixture(std::string file) {//создание объекта с помощью чтения арибутов из файла
	std::ifstream input;
	input.open(file);
	if (!input) {
		throw std::runtime_error("Невозможно найти файл");
	};

	if (!file.empty()) {
		std::vector<double> buffer;

		for (int i = 0; i < 7; i++) {
			double buf;
			input >> buf;
			//std::cout << buf;
			if (!input) {
				throw std::runtime_error("Неправильные данные");
			}
			else {
				buffer.push_back(buf);
			}
		}
		pareto1.SetForm(buffer[0]);
		pareto1.SetMu(buffer[1]);
		pareto1.SetLambda(buffer[2]);
		pareto2.SetForm(buffer[3]);
		pareto2.SetMu(buffer[4]);
		pareto2.SetLambda(buffer[5]);
		this->SetP(buffer[6]);

		if (pareto1.GetForm() <= 1 || pareto1.GetLambda() <= 0 || pareto2.GetForm() <= 1 || pareto2.GetLambda() <= 0 ) {
			throw 1;
		}


	}
	else {
		throw std::runtime_error("Файл пустой");
	}

}


//закрытые функции

double mixture::randomizer(){
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}