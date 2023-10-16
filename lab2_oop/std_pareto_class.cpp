#include "std_pareto_class.h"


//������ � ����������

double Pareto::GetForm() const {
	return v;
}

double Pareto::GetMu() const {
	return mu;
}

double Pareto::GetLambda() const {
	return lambda;
}

void Pareto::SetForm(double form_param) {
	if (form_param <= 1) {
		throw std::range_error("The new data is incorrect.Form parameter should be >1\n");
	}
	else {
		v = form_param;
	}
}

void Pareto::SetMu(double mu_param) {
	mu = mu_param;
}

void Pareto::SetLambda(double lambda_param) {
	if (lambda_param <= 0) {
		throw std::range_error("The new data is incorrect. Lambda parameter should be positive\n");
	}
	else {
		lambda = lambda_param;
	}
}


//�������� �������

double Pareto::randomizer() {//�������� ��������� �����
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

double Pareto::new_argument(const double x, const double lambda, const double mu) {//�������� ��������� ��� ������
	return (x - mu) / lambda;
}

double Pareto::normal_distribution_function(double form_param) {//���������� ������� ����������� �������������, � ����� �������� �������� ������ �� v
	double distrib_f = (0.5 * (1 + erf(form_param / sqrt(2))));
	return distrib_f;
}

double Pareto::normal_distribution_plotnost(double x) {//���������� ��������� ����������� �������������
	double plotnost = ((1 / sqrt(2 * PI)) * (exp(-((pow(x, 2) / 2)))));
	return plotnost;
}

double Pareto::k_koef() {//������ ������������ K
	return 2 * Pareto::normal_distribution_function(this->GetForm()) - 1 + (2 * this->GetForm() * Pareto::normal_distribution_plotnost(this->GetForm()) / (pow(this->GetForm(), 2) - 1));
}

double Pareto::pareto_modeling() { //������������� ��������� ��������
	double x1;
	double r1 = randomizer();
	double p = central_interval();
	if (r1 <= p) {//���1
		do {//��� 2,3
			double r2 = randomizer();
			double r3 = randomizer();
			x1 = sqrt((-2) * log(r2)) * cos(2 * PI * r3);
		} while (x1 <= (-(this->GetForm())) || (x1 >= this->GetForm()));
		return x1 * this->GetLambda() + this->GetMu();
	}
	else {//��� 4 
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

std::vector<double> Pareto::generate_x_std(int n) {//��������� �������
	std::vector <double> x;
	for (int i = 0; i < n; i++) {
		x.push_back(pareto_modeling());
	}
	sort(x.begin(), x.end());
	return x;
}

std::vector<std::pair<double, double>> Pareto::generate_xy_std(int n) {//������� � ������� x  y
	std::vector<std::pair<double, double>> result;
	auto x = this->generate_x_std(n);
	for (int i = 0; i < n; i++) {
		result.push_back(std::make_pair(x[i], this->pareto_std_distribution(x[i])));
	}
	return result;
}


//�������� �������
double Pareto::pareto_std_distribution(double x) {//��������� ������-����������� �������������
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

double Pareto::math_expect() {
	return this->GetMu();
}

double Pareto::asymmetry() {
	return 0;
}

double Pareto::dispersion() { //������������� � lambda^2 ��� ������
	return ((1 + (4 * pow(this->GetForm(), 3) * normal_distribution_plotnost(this->GetForm())) / (k_koef() * (pow(this->GetForm(), 2) - 1) * (pow(this->GetForm(), 2) - 3)))) * pow(this->GetLambda(), 2);
}

double Pareto::excess_koef() {
	return (1 / (k_koef() * (pow(dispersion(), 2)))) * ((3 * (2 * normal_distribution_function(this->GetForm()) - 1) -
		2 * this->GetForm() * normal_distribution_plotnost(this->GetForm()) * ((3 + pow(this->GetForm(), 2) - (pow(this->GetForm(), 4)) / (pow(this->GetForm(), 2) - 5))))) - 3;
}

double Pareto::central_interval() {//����������� ��������� � ����������� ��������
	return (2 * normal_distribution_function(this->GetForm()) - 1) / k_koef();
}

void Pareto::xy_output(int n) {//��������� ��������� ������� ������� n
	auto result = this->generate_xy_std(n);
	std::ofstream file("std_xy_output.txt");
	for (int i = 0; i < n; i++) {
		file << result[i].first << '\t' << result[i].second << std::endl;
	}
	file.close();
}

Pareto::Pareto(double form_param, double mu_param, double lambda_param) {
	//����������� ������ � �������� ����������
	v = form_param;
	mu = mu_param;
	lambda = lambda_param;
	if (v <= 1 || lambda <= 0) {
		throw std::range_error("Form parameter should be > 1 and lambda parameter should be positive");
	}
	else {
		std::cout << "Object has been generated succesfully" << std::endl;
		std::cout << "Params at moment: " << "form: " << v << " mu: " << mu << " lambda: " << lambda << std::endl;
	}
}

Pareto::Pareto(std::string file) {//�������� ������� � ������� ������ �������� �� �����
	std::ifstream input;
	input.open(file);
	if (!input) {
		throw std::runtime_error("Can't find file");
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
				throw std::runtime_error("Incorrect data");
			}
			else {
				buffer.push_back(buf);
			}
		}
		v = buffer[0];
		mu = buffer[1];
		lambda = buffer[2];
		if (v <= 1 || lambda <= 0) {
			throw std::range_error("Form parameter should be > 1 and lambda parameter should be positive");
		}


	}
	else {
		throw std::runtime_error("File is empty");
	}

	std::cout << "Object has been generated succesfully" << std::endl;
	std::cout << "Params at moment: " << "form: " << v << " mu: " << mu << " lambda: " << lambda << std::endl;

}

void Pareto::save_atr(std::string file) {//���������� ��������� � ����
	std::ofstream output;
	output.open(file);
	std::vector<double> buffer;
	buffer.push_back(this->GetForm());
	buffer.push_back(this->GetMu());
	buffer.push_back(this->GetLambda());
	if (!output) {
		throw std::runtime_error("Can't find file");
	};
	for (int i = 0; i < 3; i++) {
		output << buffer[i] << " ";
	}
	output.close();

}

Pareto::Pareto() {//��������� ��������� �� ���������
	v = 3.3;
	mu = 0;
	lambda = 1;
}