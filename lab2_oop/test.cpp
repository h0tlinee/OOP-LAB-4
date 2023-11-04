#include "catch.hpp"
#include "std_pareto_class.h"
#include "mixture_pareto_class.cpp"
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
    CHECK(equal(pareto.excess(), 0.05) == true);
    CHECK(equal(pareto.central_interval(), 0.999) == true);
    CHECK(equal(pareto.density(0), 0.399) == true);
    CHECK(equal(pareto.math_expect(), 0) == true);
    CHECK(equal(pareto.asymmetry(), 0) == true);
}

TEST_CASE("Shift-scale transformation") {
    Pareto pareto(3.3, 2, 2);
    CHECK(equal(pareto.dispersion(), 4.01) == true);
    CHECK(equal(abs(pareto.excess()), 2.809) == true);
    CHECK(equal(pareto.central_interval(), 0.999) == true);
    CHECK(equal(pareto.density(0), 0.12) == true);
    CHECK(equal(pareto.math_expect(), 2) == true);
    CHECK(equal(pareto.asymmetry(), 0) == true);
}

TEST_CASE("Trivial mixture case") {
    Pareto std1(3.3, 2, 2);
    Pareto std2(3.3, -2, 2);
    mixture<Pareto,Pareto> mix(std1, std2, 0.5);
    CHECK(equal(mix.math_expect(), 0) == true);
    CHECK(equal(mix.dispersion(), 8.01) == true);
    CHECK(equal(mix.asymmetry(), 0) == true);
    CHECK(equal(mix.excess(), 0.998) == true);

}

TEST_CASE("Math expect case") {
    Pareto std1(2.3,2,1);
    Pareto std2(3.3,-2,2);
    mixture<Pareto,Pareto> mix(std1, std2, 0.5);
    CHECK(equal(mix.math_expect(), 0) == true);
    CHECK(equal(mix.dispersion(), 6.58) == true);
    CHECK(equal(mix.asymmetry(), -1.31) == true);
    CHECK(equal(mix.excess(), 1.327) == true);
}

TEST_CASE("Dispersion test") {
    Pareto std1(3.3, 0, 1);
    Pareto std2(3.3, 0, 3);
    mixture<Pareto,Pareto> mix(std1, std2, 0.5);
    CHECK(equal(mix.math_expect(), 0) == true);
    CHECK(equal(mix.dispersion(), 5.015) == true);
    CHECK(equal(mix.asymmetry(), 0) == true);
    CHECK(equal(mix.excess(), -4.917) == true);
}

TEST_CASE("Late binding") {
    Pareto s1(3.3, 6, 2);
    Pareto s2(3.3, -6, 3);
    mixture<Pareto, Pareto> mx1(s1, s2, 0.5);
    Pareto s3(3.3, -10, 1);
    mixture<Pareto, mixture<Pareto, Pareto>>mx2(s3, mx1, 0.5);
    CHECK(equal(mx2.density(0), 0.00505) == true);
    CHECK(equal(mx2.math_expect(), -5) == true);
    CHECK(equal(mx2.asymmetry(), 6.176) == true);
    CHECK(equal(mx2.excess(), 2.15575) == true);
    CHECK(equal(mx2.dispersion(), 46.76) == true);
}

TEST_CASE("Empirical distribution") {
    Pareto s1(3.3, 0, 1);
    Empirical emp(s1,2000,1);
    CHECK(emp.GetSize() == 2000);
    CHECK(emp.GetIntervalsNumber() == 11);
    CHECK(equal(emp.math_expect(),0.000)==true);
}