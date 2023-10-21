#include "std_pareto_class.h"
#include "mixture_pareto_class.h"
#include "emp_pareto_class.h"
#define CATCH_CONFIG_MAIN
#include "catch.hpp" 

using namespace std;
void calculate(Pareto& std);
void calculate(mixture& mix);
void calculate(Empirical& emp);

int run_unit_tests(int argc, char** argv) {
	int result = Catch::Session().run(argc, argv);
	return result;
}


void std_menu() {
	int input=0;
	double v, mu, lambda;
	Pareto std;
	fstream file;
	while (true) {
		cout << "Действия:" << endl << "1)Ввести параметры вручную" << endl << "2)Загрузить из файла" << endl << "3)Просмотреть текущие параметры"<<endl << "4)Продолжить..." << endl << "5)Назад"<<endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "Введите парамеры v, mu, lambda" << endl;
			cin >> v >> mu >> lambda;
			std = Pareto(v, mu, lambda);
			std.save_atr("atributes.txt");
			break;
		case 2:
			std = Pareto("atributes.txt");
			break;
		case 3:
			cout << "Текущие параметры Парето-нормального распределения:" << endl << "Параметр формы =  " << std.GetForm() << ", " << "параметр сдвига =  " << std.GetMu() << ", " << "параметр масштаба =  " << std.GetLambda() << endl;
			break;
		case 4:
			calculate(std);
			break;
		case 5:
			return;
		}
		

		
	
	}
}

void mix_menu() {
	int input = 0;
	double v1, mu1, lambda1,v2,mu2,lambda2,p;
	Pareto std1;
	Pareto std2;
	mixture mix;
	fstream file;
	while (true) {
		cout << "Действия:" << endl << "1)Ввести параметры вручную" << endl << "2)Загрузить из файла" << endl << "3)Просмотреть текущие параметры" << endl << "4)Продолжить..." << endl << "5)Назад" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "Введите парамеры v1, mu1, lambda1, v2, mu2, lambda2, p" << endl;
			cin >> v1 >> mu1 >> lambda1>>v2>>mu2>>lambda2>>p;
			std1 = Pareto(v1, mu1, lambda1);
			std2 = Pareto(v2, mu2, lambda2);
			mix = mixture(std1, std2, p);
			mix.save_atr("mixtr_attr.txt");
			break;
		case 2:
			mix = mixture("mixtr_attr.txt");
			break;
		case 3:
			cout << "Текущие параметры смеси распределений:" << endl <<
				"Параметр формы первого распределения: " << mix.component1().GetForm() << ", параметр сдвига: " << mix.component1().GetMu() << ",  параметр масштаба: " << mix.component1().GetLambda() << endl
				<< "Параметр формы второго распределения: " << mix.component2().GetForm() << ", параметр сдвига: " << mix.component2().GetMu() << ",  параметр масштаба: " << mix.component2().GetLambda() << endl;
			cout << "Параметр смеси: " << mix.GetP() << endl;
			break;
		case 4:
			calculate(mix);
		case 5:
			return;
		}




	}
}

void calculate(Pareto& std) {
	int input = 0, n=2000;
	double x = 0;
	
	while (true) {
		cout << "Дальнейшие действия: " << endl << "1)Посмотреть текущие параметры" << endl << "2)Вычислить плотность в точке" << endl <<
			"3)Вычислить характеристики" << endl << "4)Строим график" << endl << "5)Выход" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "Параметр формы = " << std.GetForm() << ", параметр сдвига = " << std.GetMu() << ", параметр масштаба = " << std.GetLambda() << endl;
			break;
		case 2:
			cout << "Введите точку x" << endl;
			cin >> x;
			cout<<"Плотность в выбранной точке = "<<std.pareto_std_distribution(x) << endl;
			break;
		case 3:
			cout << "Мат.ожидание = " << std.math_expect() << endl;
			cout << "Дисперсия = " << std.dispersion() << endl;
			cout << "Асиммметрия = " << std.asymmetry() << endl;
			cout << "Эксцесс = " << std.excess_koef() << endl;
			break;
		case 4:
			cout << "Введите размер выборки для графика:" << endl;
			cin >> n;
			std.xy_output(n);
			cout << "Данные успешно получены!" << endl;
			break;
		case 5:
			return;
		}
	}
}

void calculate(mixture& mix) {
	int input = 0, n = 2000;
	double x = 0;
	while (true) {
		cout << "Дальнейшие действия: " << endl << "1)Посмотреть текущие параметры" << endl << "2)Вычислить плотность в точке" << endl <<
			"3)Вычислить характеристики" << endl << "4)Строим график" << endl << "5)Выход" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "Первое распределение:" << endl << endl << "Параметр формы = " << mix.component1().GetForm() << ", параметр сдвига = " <<
				mix.component1().GetMu() << ", параметр масштаба = " << mix.component1().GetLambda() << endl<<endl;
			cout<<"Второе распределение: "<<endl<<endl<< "Параметр формы = " << mix.component2().GetForm() << ", параметр сдвига = " <<
				mix.component2().GetMu() << ", параметр масштаба = " << mix.component2().GetLambda() << endl;
			cout << "Параметр смеси = " << mix.GetP() << endl;
			break;
		case 2:
			cout << "Введите точку x" << endl;
			cin >> x;
			cout << "Плотность в выбранной точке = " << mix.density(x) << endl;
			break;
		case 3:
			cout << "Мат.ожидание = " << mix.math_expect() << endl;
			cout << "Дисперсия = " << mix.dispersion() << endl;
			cout << "Асиммметрия = " << mix.asymmetry() << endl;
			cout << "Эксцесс = " << mix.excess() << endl;
			break;
		case 4:
			cout << "Введите размер выборки для графика:" << endl;
			cin >> n;
			mix.xy_output(n);
			cout << "Данные успешно получены!" << endl;
			break;
		case 5:
			return;
		}
	}
}

void calculate(Empirical& emp) {
	int input = 0, n = 2000;
	double x = 0;
	vector<double> sample;
	vector<pair<double, double>> graph;
	ofstream output;
	while (true) {
		cout << "Дальнейшие действия: " << endl << "1)Вывести текущие параметры" << endl << "2)Вычислить плотность в точке" << endl <<
			"3)Вычислить характеристики" << endl << "4)Построить график" << endl << "5)Выход" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "Размер выборки = " << emp.GetSize() << endl;
			cout << "Кол-во интервалов = " << emp.GetIntervalsNumber() << endl;
			break;
		case 2:
			cout << "Введите точку x" << endl;
			cin >> x;
			cout << "Плотность в выбранной точке = " << emp.density(x) << endl;
			break;
		case 3:
			cout << "Мат.ожидание = " << emp.math_expect() << endl;
			cout << "Дисперсия = " << emp.dispersion() << endl;
			cout << "Асиммметрия = " << emp.asymmetry() << endl;
			cout << "Эксцесс = " << emp.excess() << endl;
			break;
		case 4:
			emp.emp_output(output);
			cout << "Успешно!" << endl;
			break;
		case 5:
			return;
		}
	}
}



void emp_menu() {
	int input = 0;
	double v1, mu1, lambda1,v2,mu2,lambda2,p;
	Pareto std1;
	Pareto std2;
	mixture mix;
	Empirical emp(mix,10,1);
	fstream file;
	while (true) {
		cout << "Действия(на основе чего строим эмпирическое распределение):" << endl << "1)Стандартное" << endl << "2)Смесь" << endl << "3)Назад" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "Введите парамеры v, mu, lambda" << endl;
			cin >> v1 >> mu1 >> lambda1;
			std1 = Pareto(v1, mu1, lambda1);
			emp = Empirical(std1, 2000, 1);
			calculate(emp);
			break;
		case 2:
			cout << "Введите парамеры v1, mu1, lambda1, v2, mu2, lambda2, p" << endl;
			cin >> v1 >> mu1 >> lambda1 >> v2 >> mu2 >> lambda2 >> p;
			std1 = Pareto(v1, mu1, lambda1);
			std2 = Pareto(v2, mu2, lambda2);
			mix = mixture(std1, std2, p);
			emp = Empirical(mix, 2000, 1);
			calculate(emp);
			break;
		case 3:
			return;

		}



	}
}







void menu() {
	int input;
	while (true) {
		cout << "Выберите распределение:" << endl << "1)Стандартное" << endl << "2)Смесь распределений" << endl << "3)Эмпирическое"<<endl<<"4)Выход" << endl;
		cin >> input;
		switch (input) {
		case 1:
			std_menu();
			break;
		case 2:
			mix_menu();
			break;
		case 3:
			emp_menu();
			break;
		case 4:
			return;
		}
	}
}






int main(int argc, char** argv) {
	setlocale(LC_ALL, "Russian");
	try {
		/*ofstream file;
		Pareto std1(3.3, 2, 2);
		Pareto std2(3.3, -2, 2);
		mixture mix(std1, std2, 0.5);
		mix.xy_output(2000);
		Empirical emp(mix, 2000, 1);
		file.open("emp_output.txt");
		for (int i = 0; i < emp.GetSize(); i++) {
			file << emp.GetSelection()[i] << "\t" << emp.density(emp.GetSelection()[i]) << endl;
		}
		file.close();*/
		//menu();
		/*ofstream file1,file2;

		Pareto std(3.3, 2, 1);
		Empirical emp(std,2000,1);
		file1.open("test1.txt");
		for (int i = 0; i < emp.GetSize(); i++) {
			file1 << emp.GetSelection()[i] << "\t" << emp.density(emp.GetSelection()[i]) << endl;
		}
		file1.close();

		Empirical emp1(emp, 2000, 1);
		file2.open("test2.txt");
		for (int i = 0; i < emp1.GetSize(); i++) {
			file2 << emp1.GetSelection()[i] << "\t" << emp1.density(emp1.GetSelection()[i]) << endl;
		}
		file2.close();*/

		
		

		
		


		
		
	}
	catch (runtime_error ex) {
		cout << "Ошибка при работе с файлами: " << ex.what();
	}
	catch (int error) {
		if (error == 1) {
			cout << "Параметры распределения указаны неккоректно" << endl;
		}
	}
	run_unit_tests(argc, argv);

}