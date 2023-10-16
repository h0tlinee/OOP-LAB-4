#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>



const double PI = 3.141592653589793;

class Pareto {
private:
	//атрибуты класса
	double v;
	double mu;
	double lambda;
	//функции, вызываваемые только в методах класса, пользователь не имеет к ним доступа
	double randomizer();
	double new_argument(double x, const double lambda, const double mu);
	double normal_distribution_function(double form_param);
	double normal_distribution_plotnost(double x);
	double k_koef();
	double pareto_modeling();
	std::vector<double> generate_x_std(int n);
	std::vector<std::pair<double, double>> generate_xy_std(int n);

public:
	Pareto(double form_param, double mu_param, double lambda_param);//конструктор класса с заданием параметров
	Pareto(std::string file);//конструктор с вводом из файла
	Pareto();//конструктор с параметрами по умолчанию
		


	//блок работающий с атрибутами
	double GetForm() const;
	double GetMu() const;
	double GetLambda() const;
	void SetForm(double form_param);
	void SetMu(double mu_param);
	void SetLambda(double lambda_param);
	//вычисление распределения в точке
	double pareto_std_distribution(double x);
	void test();
	double math_expect();
	double asymmetry();
	double dispersion();
	double excess_koef();
	double central_interval();
	//вывод в файл, сохранение
	void xy_output(int n);
	void save_atr(std::string file);
	


	

	
};