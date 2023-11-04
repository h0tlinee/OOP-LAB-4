#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>



const double PI = 3.141592653589793;

class IDistribution {
public:
	double virtual density(const double x) const = 0;
	double virtual math_expect() const = 0;
	double virtual dispersion() const = 0;
	double virtual excess() const = 0;
	double virtual asymmetry() const = 0;
	double virtual modeling() const = 0;
	std::vector<double> virtual generate_x(int n) const = 0;
	std::vector<std::pair<double, double>> virtual generate_xy(int n) const = 0;
};

class IPersistend {
	void virtual load_attr(std::ifstream& input) = 0;
	void virtual save_attr(std::ofstream& output) = 0;
};


class Pareto:public IDistribution,public IPersistend {
private:
	//атрибуты класса
	double v;
	double mu;
	double lambda;
	//функции, вызываваемые только в методах класса, пользователь не имеет к ним доступа
	double randomizer() const;
	double new_argument(double x, const double lambda, const double mu) const;
	double normal_distribution_function(double form_param) const;
	double normal_distribution_plotnost(double x) const;
	double k_koef() const;
	//double pareto_modeling();
	//std::vector<double> generate_x_std(int n);
	std::vector<std::pair<double, double>> generate_xy(int n) const override;

public:
	Pareto(double form_param, double mu_param, double lambda_param);//конструктор класса с заданием параметров
	Pareto(std::string file);//конструктор с вводом из файла
	Pareto();//конструктор с параметрами по умолчанию
	std::vector<double> generate_x(int n) const override;


	//блок работающий с атрибутами
	double GetForm() const;
	double GetMu() const;
	double GetLambda() const;
	void SetForm(double form_param);
	void SetMu(double mu_param);
	void SetLambda(double lambda_param);
	//вычисление распределения в точке
	double modeling() const override;
	double density(double x) const override;
	void test();
	double math_expect() const override;
	double asymmetry() const override;
	double dispersion() const override;
	double excess() const override;
	double central_interval() const;
	//вывод в файл, сохранение
	void xy_output(int n);
	void save_attr(std::ofstream& output) override;
	void load_attr(std::ifstream& input) override;
	


	

	
};