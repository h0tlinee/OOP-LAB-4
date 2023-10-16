#include "std_pareto_class.h"
#define CATCH_CONFIG_MAIN
#include "catch.hpp" 

using namespace std;

int run_unit_tests(int argc, char** argv) {
	int result = Catch::Session().run(argc, argv);
	return result;
}


int main(int argc, char** argv) {
	setlocale(LC_ALL, "Russian");
	try {
		
		Pareto pareto("atributes.txt");
		cout << "���������: " << pareto.pareto_std_distribution(0) << endl;
		cout << "�����������: " << pareto.math_expect() << endl;
		cout << "����������: " << pareto.asymmetry() << endl;
		cout << "�������: " << pareto.excess_koef() << endl;
		cout << "���������: " << pareto.dispersion() << endl;
		cout << "����������� ��������: " << pareto.central_interval() << endl;
		pareto.xy_output(2000);
		
		
	}
	catch (range_error ex) {
		cout << "Error:  " << ex.what();
	}
	catch (runtime_error ex) {
		cout << "Error while working with files: " << ex.what();
	}
	//run_unit_tests(argc, argv);

}