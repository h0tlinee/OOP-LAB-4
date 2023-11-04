#include "std_pareto_class.h"


//работа с атрибутами

double Pareto::GetForm() const {
	return v;
}

double Pareto::GetMu() const {
	return mu;
}

double Pareto::GetLambda() const {
	return lambda;
}

void Pareto::SetForm(double form_param){
	if (form_param <= 1) {
		throw 1;
	}
	else {
		v = form_param;
	}
}

void Pareto::SetMu(double mu_param){
	mu = mu_param;
}

void Pareto::SetLambda(double lambda_param) {
	if (lambda_param <= 0) {
		throw 1;
	}
	else {
		lambda = lambda_param;
	}
}


//закрытые функции

double Pareto::randomizer() const {//гнератор случайных чисел
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

double Pareto::new_argument(const double x, const double lambda, const double mu) const {//пересчет аргумента для сдвига
	return (x - mu) / lambda;
}

double Pareto::normal_distribution_function(double form_param) const {//вычисление функции нормального распределения, в нашем варианте зависима только от v
	double distrib_f = (0.5 * (1 + erf(form_param / sqrt(2))));
	return distrib_f;
}

double Pareto::normal_distribution_plotnost(double x) const {//вычисление плотности нормального распределения
	double plotnost = ((1 / sqrt(2 * PI)) * (exp(-((pow(x, 2) / 2)))));
	return plotnost;
}

double Pareto::k_koef() const {//расчет коэффициента K
	return 2 * Pareto::normal_distribution_function(this->GetForm()) - 1 + (2 * this->GetForm() * Pareto::normal_distribution_plotnost(this->GetForm()) / (pow(this->GetForm(), 2) - 1));
}

double Pareto::modeling() const { //моделирование случайной величины
	double x1;
	double r1 = randomizer();
	double p = central_interval();
	if (r1 <= p) {//шаг1
		do {//шаг 2,3
			double r2 = randomizer();
			double r3 = randomizer();
			x1 = sqrt((-2) * log(r2)) * cos(2 * PI * r3);
		} while (x1 <= (-(this->GetForm())) || (x1 >= this->GetForm()));
		return x1 * this->GetLambda() + this->GetMu();
	}
	else {//шаг 4 
		double r4 = randomizer();
		double x2 = this->GetForm() / (pow(r4, (1 / (pow(this->GetForm(), 2) - 1))));
		if (r1 < ((1 + central_interval()) / 2)) {
			return x2 * this->GetLambda() + this->GetMu();
		}
		else {
			return ((-(x2)) * this->GetLambda() + this->GetMu());
		}
	}
}

std::vector<double> Pareto::generate_x(int n) const {//генерация выборки
	std::vector <double> x;
	for (int i = 0; i < n; i++) {
		x.push_back(modeling());
	}
	sort(x.begin(), x.end());
	return x;
}

std::vector<std::pair<double, double>> Pareto::generate_xy(int n) const {//выборка в формате x  y
	std::vector<std::pair<double, double>> result;
	auto x = this->generate_x(n);
	for (int i = 0; i < n; i++) {
		result.push_back(std::make_pair(x[i], this->density(x[i])));
	}
	return result;
}


//открытые функции
double Pareto::density(double x) const {//плотность парето-нормального распределения
	double new_arg = Pareto::new_argument(x, this->GetLambda(), this->GetMu());
	if (abs(new_arg) <= this->GetForm()) {
		return ((1 / k_koef()) * normal_distribution_plotnost(new_arg)) / this->GetLambda();
	}
	else {
		return ((1 / k_koef()) * normal_distribution_plotnost(this->GetForm()) * (pow((this->GetForm() / abs(new_arg)), pow(this->GetForm(), 2)))) / this->GetLambda();
	}
}

void Pareto::test() {
	std::cout << this->GetForm();
}

double Pareto::math_expect() const {
	return this->GetMu();
}

double Pareto::asymmetry() const {
	return 0;
}

double Pareto::dispersion() const { //увеличивается в lambda^2 при сдвиге
	return ((1 + (4 * pow(this->GetForm(), 3) * normal_distribution_plotnost(this->GetForm())) / (k_koef() * (pow(this->GetForm(), 2) - 1) * (pow(this->GetForm(), 2) - 3)))) * pow(this->GetLambda(), 2);
}

double Pareto::excess() const {
	return (1 / (k_koef() * (pow(dispersion(), 2)))) * ((3 * (2 * normal_distribution_function(this->GetForm()) - 1) -
		2 * this->GetForm() * normal_distribution_plotnost(this->GetForm()) * ((3 + pow(this->GetForm(), 2) - (pow(this->GetForm(), 4)) / (pow(this->GetForm(), 2) - 5))))) - 3;
}

double Pareto::central_interval() const {//вероятность попадания в центральный интервал
	return (2 * normal_distribution_function(this->GetForm()) - 1) / k_koef();
}

void Pareto::xy_output(int n) {//генерация случайных величин размера n
	auto result = this->generate_xy(n);
	std::ofstream file("std_xy_output.txt");
	for (int i = 0; i < n; i++) {
		file << result[i].first << '\t' << result[i].second << std::endl;
	}
	file.close();
}

Pareto::Pareto(double form_param, double mu_param, double lambda_param):
	v(form_param > 0 ? form_param : throw 1), mu(mu_param), lambda(lambda_param > 0 ? lambda_param : throw 1){
}

Pareto::Pareto(std::string file) {//создание объекта с помощью чтения арибутов из файла
	std::ifstream input;
	input.open(file);
	if (!input) {
		throw std::runtime_error("Ошибка при поиске файла");
	};

	if (!file.empty()) {
		Pareto::v;
		Pareto::mu;
		Pareto::lambda;
		std::vector<double> buffer;

		for (int i = 0; i < 3; i++) {
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
		v = buffer[0];
		mu = buffer[1];
		lambda = buffer[2];
		if (v <= 1 || lambda <= 0) {
			throw 1;
		}


	}
	else {
		throw std::runtime_error("Файл пустой");
	};

}

void Pareto::save_attr(std::ofstream& output){//сохранение атрибутов в файл
	if (!output) {
		throw 0;
	};
	output << v<<std::endl << mu<<std::endl << lambda<<std::endl;
}

void Pareto::load_attr(std::ifstream& input){
	double v, mu, lambda;
	
	if (!input.is_open()) {
		throw 0;
	}
	else {
		input >> v >> mu >> lambda;
		if (v <= 1 || lambda < 0) {
			throw 1;
		}
	}
	this->v = v;
	this->mu = mu;
	this->lambda = lambda;

}

Pareto::Pareto() :v(3.3), mu(0), lambda(1) {}//параметры выбранные по умолчанию