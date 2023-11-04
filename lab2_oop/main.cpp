#include "std_pareto_class.h"
#include "mixture_pareto_class.cpp"
#include "emp_pareto_class.h"
#define CATCH_CONFIG_MAIN
#include "catch.hpp" 

using namespace std;
//void calculate(Pareto& std);
//void calculate(mixture& mix);
//void calculate(Empirical& emp);

int run_unit_tests(int argc, char** argv) {
	int result = Catch::Session().run(argc, argv);
	return result;
}


/*void std_menu() {
	int input=0;
	double v, mu, lambda;
	Pareto std;
	fstream file;
	while (true) {
		cout << "��������:" << endl << "1)������ ��������� �������" << endl << "2)��������� �� �����" << endl << "3)����������� ������� ���������"<<endl << "4)����������..." << endl << "5)�����"<<endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "������� �������� v, mu, lambda" << endl;
			cin >> v >> mu >> lambda;
			std = Pareto(v, mu, lambda);
			std.save_atr("atributes.txt");
			break;
		case 2:
			std = Pareto("atributes.txt");
			break;
		case 3:
			cout << "������� ��������� ������-����������� �������������:" << endl << "�������� ����� =  " << std.GetForm() << ", " << "�������� ������ =  " << std.GetMu() << ", " << "�������� �������� =  " << std.GetLambda() << endl;
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
		cout << "��������:" << endl << "1)������ ��������� �������" << endl << "2)��������� �� �����" << endl << "3)����������� ������� ���������" << endl << "4)����������..." << endl << "5)�����" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "������� �������� v1, mu1, lambda1, v2, mu2, lambda2, p" << endl;
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
			cout << "������� ��������� ����� �������������:" << endl <<
				"�������� ����� ������� �������������: " << mix.component1().GetForm() << ", �������� ������: " << mix.component1().GetMu() << ",  �������� ��������: " << mix.component1().GetLambda() << endl
				<< "�������� ����� ������� �������������: " << mix.component2().GetForm() << ", �������� ������: " << mix.component2().GetMu() << ",  �������� ��������: " << mix.component2().GetLambda() << endl;
			cout << "�������� �����: " << mix.GetP() << endl;
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
		cout << "���������� ��������: " << endl << "1)���������� ������� ���������" << endl << "2)��������� ��������� � �����" << endl <<
			"3)��������� ��������������" << endl << "4)������ ������" << endl << "5)�����" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "�������� ����� = " << std.GetForm() << ", �������� ������ = " << std.GetMu() << ", �������� �������� = " << std.GetLambda() << endl;
			break;
		case 2:
			cout << "������� ����� x" << endl;
			cin >> x;
			cout<<"��������� � ��������� ����� = "<<std.pareto_std_distribution(x) << endl;
			break;
		case 3:
			cout << "���.�������� = " << std.math_expect() << endl;
			cout << "��������� = " << std.dispersion() << endl;
			cout << "����������� = " << std.asymmetry() << endl;
			cout << "������� = " << std.excess_koef() << endl;
			break;
		case 4:
			cout << "������� ������ ������� ��� �������:" << endl;
			cin >> n;
			std.xy_output(n);
			cout << "������ ������� ��������!" << endl;
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
		cout << "���������� ��������: " << endl << "1)���������� ������� ���������" << endl << "2)��������� ��������� � �����" << endl <<
			"3)��������� ��������������" << endl << "4)������ ������" << endl << "5)�����" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "������ �������������:" << endl << endl << "�������� ����� = " << mix.component1().GetForm() << ", �������� ������ = " <<
				mix.component1().GetMu() << ", �������� �������� = " << mix.component1().GetLambda() << endl<<endl;
			cout<<"������ �������������: "<<endl<<endl<< "�������� ����� = " << mix.component2().GetForm() << ", �������� ������ = " <<
				mix.component2().GetMu() << ", �������� �������� = " << mix.component2().GetLambda() << endl;
			cout << "�������� ����� = " << mix.GetP() << endl;
			break;
		case 2:
			cout << "������� ����� x" << endl;
			cin >> x;
			cout << "��������� � ��������� ����� = " << mix.density(x) << endl;
			break;
		case 3:
			cout << "���.�������� = " << mix.math_expect() << endl;
			cout << "��������� = " << mix.dispersion() << endl;
			cout << "����������� = " << mix.asymmetry() << endl;
			cout << "������� = " << mix.excess() << endl;
			break;
		case 4:
			cout << "������� ������ ������� ��� �������:" << endl;
			cin >> n;
			mix.xy_output(n);
			cout << "������ ������� ��������!" << endl;
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
		cout << "���������� ��������: " << endl << "1)������� ������� ���������" << endl << "2)��������� ��������� � �����" << endl <<
			"3)��������� ��������������" << endl << "4)��������� ������" << endl << "5)�����" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "������ ������� = " << emp.GetSize() << endl;
			cout << "���-�� ���������� = " << emp.GetIntervalsNumber() << endl;
			break;
		case 2:
			cout << "������� ����� x" << endl;
			cin >> x;
			cout << "��������� � ��������� ����� = " << emp.density(x) << endl;
			break;
		case 3:
			cout << "���.�������� = " << emp.math_expect() << endl;
			cout << "��������� = " << emp.dispersion() << endl;
			cout << "����������� = " << emp.asymmetry() << endl;
			cout << "������� = " << emp.excess() << endl;
			break;
		case 4:
			emp.emp_output(output);
			cout << "�������!" << endl;
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
		cout << "��������(�� ������ ���� ������ ������������ �������������):" << endl << "1)�����������" << endl << "2)�����" << endl << "3)�����" << endl;
		cin >> input;
		switch (input) {
		case 1:
			cout << "������� �������� v, mu, lambda" << endl;
			cin >> v1 >> mu1 >> lambda1;
			std1 = Pareto(v1, mu1, lambda1);
			emp = Empirical(std1, 2000, 1);
			calculate(emp);
			break;
		case 2:
			cout << "������� �������� v1, mu1, lambda1, v2, mu2, lambda2, p" << endl;
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
		cout << "�������� �������������:" << endl << "1)�����������" << endl << "2)����� �������������" << endl << "3)������������"<<endl<<"4)�����" << endl;
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
	}*/
//}






int main(int argc, char** argv) {
	setlocale(LC_ALL, "Russian");
	try {
		
		
		
		


		
		
	}
	catch (int error){
		if (error == 0){
			cout << "������ ��� ������ � �������: "<<endl;
		}
		if (error == 1) {
			cout << "��������� ������������� ������� �����������" << endl;
		}
		
	}

	//run_unit_tests(argc, argv);

}