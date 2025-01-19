#include "quadpack_cpp/quadpack.hpp"

#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <iomanip>

const double EPSABS = 0.0;
const double EPSREL = std::pow(10, std::log10(std::numeric_limits<double>::epsilon())/2.0 + 1.0);

// *********************************************** //
// prototype
// *********************************************** //

void test_dqag();
void test_dqagi();
void test_dqagp();
void test_dqags();
void test_dqawc();
void test_dqawf();
void test_dqawo();
void test_dqaws();
void test_dqng();

// *********************************************** //
// main
// *********************************************** //

int main()
{
    test_dqag();
    test_dqagi();
    test_dqagp();
    test_dqags();
    test_dqawc();
    test_dqawf();
    test_dqawo();
    test_dqaws();
    test_dqng();
    return 0;
}

// *********************************************** //
// check_result
// *********************************************** //

void check_result(std::string routine, const double result, const double answer, const int neval)
{
    // double relerr = std::abs((result - answer)/answer);
    double relerr = std::abs(result - answer);
    std::cout << std::scientific;
    if (relerr < EPSREL) {
        std::cout << std::setw(10) << routine << " "
                  << std::setw(15) << result << " "
                  << std::setw(15) << answer << " "
                  << std::setw(15) << relerr << " "
                  << std::setw(8)  << neval  << std::endl;
                  
    } else {
        std::cout << std::setw(10) << routine << " "
                  << std::setw(15) << result << " "
                  << std::setw(15) << answer << " "
                  << std::setw(15) << relerr << " "
                  << std::setw(8)  << neval  << " "
                  << "failed" << std::endl;
    }
}

// *********************************************** //
// test_dqag
// *********************************************** //

double func_dqag(double x, void *data)
{
    return 2.0 /(2.0 + sin(10.0*M_PI*x));
}


void test_dqag()
{
    const double a = 0.0;
    const double b = 1.0;
    const int limit = 200;
    const int lenw = limit*4;
    const double answer = 2.0 / std::sqrt(3.0);

    double abserr, result, work[lenw];
    int ier, iwork[limit], last, neval;

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    for (int key = 1; key <=5; ++key) {
        quadpack_cpp::dqag(func_dqag, a, b, epsabs, epsrel, key, result, abserr, neval, 
            ier, limit, lenw, last, iwork, work, nullptr);

        std::string routine = "dqag-" + std::to_string(key);
        check_result(routine, result, answer, neval);
    }

    return;
}

// *********************************************** //
// test_dqagi
// *********************************************** //

double func_dqagi(double x, void *data)
{
    if (x > 0.0) {
        return std::sqrt(x)*std::log(x)/((x + 1.0)*(x + 2.0));
    } else {
        return 0.0;
    }
}


void test_dqagi()
{
    const double boun  = 0.0;
    const int    inf   = 1;
    const int    limit = 100;
    const int    lenw  = limit*4;
    const double answer = std::sqrt(2.0)*M_PI*std::log(2.0);

    double abserr, result, work[lenw];
    int ier, iwork[limit], last, neval;

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    quadpack_cpp::dqagi(func_dqagi, boun, inf, epsabs, epsrel, result, abserr, neval,
        ier, limit, lenw, last, iwork, work, nullptr);

    std::string routine = "dqagi";
    check_result(routine, result, answer, neval);

    return;
}

// *********************************************** //
// test_dqagp
// *********************************************** //

double func_dqagp(double x, void *data)
{
    if (x != 1.0/7.0 && x != 2.0/3.0) {
        return std::pow(std::abs(x - 1.0/7.0), -0.25) 
             * std::pow(std::abs(x - 2.0/3.0), -0.55);
    } else {
        return 0.0;
    }
}


void test_dqagp()
{
    const int npts2 = 4;
    const int limit = 100;
    const int leniw = limit*2 + npts2;
    const int lenw  = limit*4 + npts2;
    const double a  = 0.0;
    const double b  = 1.0;
    const double points[npts2] = {1.0/7.0, 2.0/3.0, 0.0, 0.0};
    const double answer = 4.25368768812224946110743394858422;

    double abserr, result, work[lenw];
    int ier, iwork[leniw], last, neval;

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    quadpack_cpp::dqagp(func_dqagp, a, b, npts2, points, epsabs, epsrel, result, abserr, 
        neval, ier, leniw, lenw, last, iwork, work, nullptr);

    std::string routine = "dqagp";
    check_result(routine, result, answer, neval);

    return;
}

// *********************************************** //
// test_dqags
// *********************************************** //

double func_dqags(double x, void *data) 
{
    if (x > 0.0) {
        return 1.0 / std::sqrt(x);
    } else {
        return 0.0;
    }
}

void test_dqags()
{
    const double a = 0.0;
    const double b = 1.0;
    const int limit = 100;
    const int lenw = limit*4;
    const double answer = 2.0;

    double abserr, result, work[lenw];
    int ier, iwork[limit], last, neval;

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    quadpack_cpp::dqags(func_dqags, a, b, epsabs, epsrel, result, abserr, neval, 
        ier, limit, lenw, last, iwork, work, nullptr);

    check_result("dqags", result, answer, neval);

    return;
}

// *********************************************** //
// test_dqawc
// *********************************************** //

double func_dqawc(double x, void *data)
{
    return 1.0 /(x*x + 1.0e-4);
}

void test_dqawc()
{
    const double a = -1.0;
    const double b =  1.0;
    const double c =  0.5;
    const int limit = 200;
    const int lenw  = limit*4;
    const double answer = -628.461728506562366229080921522473;

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    double abserr, result, work[lenw];
    int ier, iwork[limit], last, neval;

    quadpack_cpp::dqawc(func_dqawc, a, b, c, epsabs, epsrel, result, abserr, neval,
        ier, limit, lenw, last, iwork, work, nullptr);

    check_result("dqawc", result, answer, neval);

    return;
}

// *********************************************** //
// test_dqawf
// *********************************************** //

double func_dqawf(double x, void *data)
{
    if (x > 0.0) {
        return std::sin(50.0*x)/(x*std::sqrt(x));
    } else {
        return 0.0;
    }
}

void test_dqawf()
{
    const double a      = 0.0;
    const double omega  = 8.0;
    const int integr    = 2;
    const int limlst    = 100;
    const int limit     = 500;
    const int leniw     = limit*2 + limlst;
    const int maxp1     = 100;
    const int lenw      = leniw*2 + maxp1*25;
    const double answer = std::sqrt(29.0*M_PI) - std::sqrt(21.0*M_PI);

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    double abserr, result, work[lenw];
    int ier, iwork[leniw], last, lst, neval;

    quadpack_cpp::dqawf(func_dqawf, a, omega, integr, epsrel, result, abserr, neval,
                   ier, limlst, lst, leniw, maxp1, lenw, iwork, work, nullptr);

    check_result("dqawf", result, answer, neval);

    return;
}

// *********************************************** //
// test_dqawo
// *********************************************** //

double func_dqawo(double x, void *data)
{
    if (x > 0.0) {
        return std::exp(-x)*std::log(x);
    } else {
        return 0.0;
    }
}

void test_dqawo()
{
    const double a      = 0.0;
    const double b      = 1.0;
    const double omega  = 10.0;
    const int integr    = 1;
    const int limit     = 100;
    const int leniw     = limit*2;
    const int maxp1     = 21;
    const int lenw      = limit*4 + maxp1*25;
    const double answer = -0.177639206511388980501003222731069;

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    double abserr, result, work[lenw];
    int ier, iwork[leniw], last, neval;

    quadpack_cpp::dqawo(func_dqawo, a, b, omega, integr, epsabs, epsrel, result, abserr,
                   neval, ier, leniw, maxp1, lenw, last, iwork, work, nullptr);

    check_result("dqawo", result, answer, neval);

    return;
}

// *********************************************** //
// test_dqaws
// *********************************************** //

double func_dqaws(double x, void *data)
{
    return std::sin(10.0*x);
}

void test_dqaws()
{
    const double a = 0.0;
    const double b = 1.0;
    const double alfa = -0.5;
    const double beta = -0.5;
    const int integr = 1;
    const int limit  = 100;
    const int lenw   = limit*4;
    const double answer = 0.535019056922365344126359;

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    double abserr, result, work[lenw];
    int ier, iwork[limit], last, neval;

    quadpack_cpp::dqaws(func_dqaws, a, b, alfa, beta, integr, epsabs, epsrel, result,
        abserr, neval, ier, limit, lenw, last, iwork, work, nullptr);

    check_result("dqaws", result, answer, neval);

    return;
}

// *********************************************** //
// test_dqng
// *********************************************** //

double func_dqng(double x, void *data)
{
    return std::exp(x) / (x*x + 1.0);
}

void test_dqng()
{
    const double a = 0.0;
    const double b = 1.0;
    const double answer = 1.27072413983362022013785374440150;

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    double abserr, result;
    int ier, neval;

    quadpack_cpp::dqng(func_dqng, a, b, epsabs, epsrel, result, abserr, neval, ier, nullptr);

    check_result("dqng", result, answer, neval);

    return;
}