#include "catch.hpp"
#include "std_pareto_class.h"

bool equal(const double& x, const double& y) {
    if (abs(x - y) <= 0.005) {
        return true;
    }
    else return false;
}

TEST_CASE("Standart pareto distribution") {
    Pareto pareto(3.3,0,1);
    CHECK(equal(pareto.dispersion(),1)==true);
    CHECK(equal(pareto.excess_koef(), 0.05) == true);
    CHECK(equal(pareto.central_interval(), 0.999) == true);
    CHECK(equal(pareto.pareto_std_distribution(0), 0.399) == true);
    CHECK(equal(pareto.math_expect(), 0) == true);
    CHECK(equal(pareto.asymmetry(), 0) == true);
}

TEST_CASE("Shift-scale transformation") {
    Pareto pareto(3.3, 2, 2);
    CHECK(equal(pareto.dispersion(), 4.01) == true);
    CHECK(equal(abs(pareto.excess_koef()), 2.809) == true);
    CHECK(equal(pareto.central_interval(), 0.999) == true);
    CHECK(equal(pareto.pareto_std_distribution(0), 0.12) == true);
    CHECK(equal(pareto.math_expect(), 2) == true);
    CHECK(equal(pareto.asymmetry(), 0) == true);
}