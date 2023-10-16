#include "mixture_pareto_class.h"


mixture::mixture() {//конструктор по умолчанию
	p = 0.5;
	Pareto std1;
	Pareto std2;
	std1.SetForm(3.3);
	std2.SetForm(3.3);
	std1.SetLambda(1);
	std2.SetLambda(2);
	std1.SetMu(-2);
	std2.SetMu(2);
	

}

mixture::mixture(Pareto& std1, Pareto& std2, double _p)://конструктор из существующих двух стандратных распределений
	p(_p<0 and _p>1 ? _p:throw 1), pareto1(std1), pareto2(std2){}

