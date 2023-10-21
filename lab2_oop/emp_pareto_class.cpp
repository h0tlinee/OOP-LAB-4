#include "emp_pareto_class.h"



//3 конструктора с моделированием случайных величин

Empirical::Empirical(Pareto& std, int n, int _k)://эмпирическая выборка по стандратному распределению
	size(n > 1 ? n : throw 1), k(_k > 2 ? _k : (int)log2(size) + 1) {
		selection = std.generate_x_std(size);
		empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

Empirical::Empirical(mixture& mix, int n, int _k) ://эмпирическая выборка по смеси распределений
	size(n > 1 ? n : throw 1), k(_k > 2 ? _k : ((int)log2(size) + 1)) {
	selection = mix.generate_x(n);
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

//ВСТАВЬ КОНСТРУКТОР С ЭМП ВЫБОРКОЙ С НУЛЯ
Empirical::Empirical( Empirical& ed, int n, int _k) ://эмпирическая выборка
	size(n), k(_k > 2 ? k : (int)log2(size) + 1) {
	selection = ed.generate_x(size);
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}



Empirical::Empirical(Empirical& emp) ://конструктор копирования
	size(emp.size > 1 ? emp.size : throw 1), k(emp.k > 2 ? emp.k : ((int)log2(size) + 1)), selection(emp.selection), empirical_density(emp.empirical_density) {
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

Empirical::Empirical(std::vector<double>& selection) ://для абсолютно случайной любой выборки
	size(selection.size()), k((int)log2(size) + 1) {
	this->selection = selection;
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

Empirical& Empirical::operator = (const Empirical& emp) { //оператор присваивания
	if (this == &emp) {
		return *this;
	}
	selection = emp.selection;
	empirical_density = emp.empirical_density;
	size = emp.size;
	k = emp.k;
	return *this;
}

Empirical::Empirical(std::string file_name) {//конструктор чтением из файла
	std::ifstream file;
	file.open(file_name);
	if (!file.is_open()) {
		throw std::runtime_error("Невозможно найти файл");
	}
	double x;
	while (file >> x) {
		selection.push_back(x);
	}
	size = selection.size();
	k = (int)log2(size) + 1;
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}


//деструктор
Empirical::~Empirical() {
	selection.clear();
	empirical_density.clear();
}


//сеттер кол-ва интервалов
void Empirical::SetK( int k) {
	if (k >= 2) {
		this->k = k;
	}
	else {
		this->k = (int)log2(size) + 1;
	}
	empirical_density = create_empirical_density(delta_calc(), create_intervals(delta_calc()));
}

//рандомайзер

double Empirical::randomizer(){
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

//геттеры

std::vector<double> Empirical::GetSelection(){
	return selection;
}

std::vector<double> Empirical::GetEmpDensity(){
	return empirical_density;
}

int Empirical::GetSize(){
	return size;
}

int Empirical::GetIntervalsNumber(){
	return k;
}

//вычисления

double Empirical::delta_calc(){
	return (selection[size - 1] - selection[0]) / (k);
}

std::vector<double> Empirical::create_intervals(const double delta){
	std::vector<double> intervals;
	double slider = selection[0];
	for (int i = 0; i < k + 1; ++i) {
		intervals.push_back(slider);
		slider += delta;
	}
	return intervals;
}

std::vector<double> Empirical::create_empirical_density(const double delta, std::vector<double> intervals){
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



double Empirical::density(double x) {//плотность в точке
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

double Empirical::math_expect(){
	double sum = 0;
	for (int i = 0; i < size; ++i) {
		sum += selection[i];
	}
	return sum / size;
}

double Empirical::dispersion(){
	double sum = 0;
	double exp_val = math_expect();
	for (int i = 0; i < size; ++i) {
		sum += pow(selection[i] - exp_val, 2);
	}
	return sum / size;
}

double Empirical::asymmetry(){
	double sum = 0;
	double exp_val = math_expect();
	for (int i = 0; i < size; ++i) {
		sum += pow(selection[i] - exp_val, 3);
	}
	return sum / (size * pow(dispersion(), 3 / 2));
}

double Empirical::excess(){
	double sum = 0;
	double exp_val = math_expect();
	for (int i = 0; i < size; ++i) {
		sum += pow(selection[i] - exp_val, 4);
	}
	return (sum / (size * pow(dispersion(), 2))) - 3;
}



double Empirical::qumulative_probability(int i)  {//кумулятивка
	double q = 0;
	for (int j = 0; j <= i; ++j) {
		q += empirical_density[j];
	}
	return q;
}

double Empirical::top_bound(){
	double sum = 0;
	for (int i = 0; i < k; ++i) {
		sum += empirical_density[i];
	}
	return sum;
}

double Empirical::rand_var(){//моделирование случайной величины
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

std::vector<double>  Empirical::generate_x(int n) {//генерация выборки x
	std::vector<double> result;
	for (int i = 0; i < n; ++i) {
		result.push_back(rand_var());
	}
	sort(result.begin(), result.end());
	return result;
}

std::vector<std::pair<double, double>> Empirical::generate_xy(){//генерация выборки xy
	std::vector <std::pair<double, double >> result;
	for (int j = 0; j < size; ++j) {
		result.push_back(std::make_pair(selection[j], density(selection[j])));
	}
	return result;
}
/// 
void Empirical::emp_output(std::ofstream& file){//вывод в файл
	auto pairs = generate_xy();
	file.open("emp_output.txt");
	for (int i = 0; i < size; ++i) {
		file << pairs[i].first << "\t" << pairs[i].second << std::endl;
	}
	file.close();
}