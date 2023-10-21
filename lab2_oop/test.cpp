#include "catch.hpp"
#include "std_pareto_class.h"
#include "mixture_pareto_class.h"
#include "emp_pareto_class.h"

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

TEST_CASE("Trivial mixture case") {
    Pareto std1(3.3, 2, 2);
    Pareto std2(3.3, -2, 2);
    mixture mix(std1, std2, 0.5);
    CHECK(equal(mix.math_expect(), 0) == true);
    CHECK(equal(mix.dispersion(), 8.01) == true);
    CHECK(equal(mix.asymmetry(), 0) == true);
    CHECK(equal(mix.excess(), 0.998) == true);

}

TEST_CASE("Math expect case") {
    Pareto std1(2.3,2,1);
    Pareto std2(3.3,-2,2);
    mixture mix(std1, std2, 0.5);
    CHECK(equal(mix.math_expect(), 0) == true);
    CHECK(equal(mix.dispersion(), 6.58) == true);
    CHECK(equal(mix.asymmetry(), -1.31) == true);
    CHECK(equal(mix.excess(), 1.327) == true);
}

TEST_CASE("Dispersion test") {
    Pareto std1(3.3, 0, 1);
    Pareto std2(3.3, 0, 3);
    mixture mix(std1, std2, 0.5);
    CHECK(equal(mix.math_expect(), 0) == true);
    CHECK(equal(mix.dispersion(), 5.015) == true);
    CHECK(equal(mix.asymmetry(), 0) == true);
    CHECK(equal(mix.excess(), -4.917) == true);
}
