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

// *********************************************** //
// main
// *********************************************** //

int main()
{
    test_dqag();
    test_dqagi();
    test_dqagp();
    test_dqags();
    return 0;
}

// *********************************************** //
// check_result
// *********************************************** //

void check_result(std::string routine, const double result, const double answer, const int neval)
{
    double relerr = std::abs((result - answer)/answer);
    std::cout << std::scientific;
    if (relerr < EPSREL) {
        std::cout << std::setw(10) << routine << " "
                  << result << " "
                  << answer << " "
                  << relerr << " "
                  << neval  << std::endl;
                  
    } else {
        std::cout << std::setw(10) << routine << " "
                  << result << " "
                  << answer << " "
                  << relerr << " "
                  << neval  << " "
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

// double func_dqags(double x, void *data) 
// {
//     // double f = 1.4;
//     double f = *(double*)(data);
//     return 1.0 / std::sqrt(std::log(x)/f + 1.0/x - 1.0);
// }

void test_dqags()
{
    const double a = 0.0;
    const double b = 1.0;
    const int limit = 100;
    const int lenw = limit*4;
    const double answer = 2.0;
    double f = 1.4;
    // const double answer = 2.567172178745525; 
    // double f = 2.0;
    // const double answer = 2.080463121721123;

    double abserr, result, work[lenw];
    int ier, iwork[limit], last, neval;

    const double epsabs = EPSABS;
    const double epsrel = EPSREL;

    quadpack_cpp::dqags(func_dqags, a, b, epsabs, epsrel, result, abserr, neval, 
        ier, limit, lenw, last, iwork, work, &f);

    check_result("dqags", result, answer, neval);

    return;
}