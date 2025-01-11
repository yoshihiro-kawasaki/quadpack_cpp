#ifndef QUADPACK_HPP_
#define QUADPACK_HPP_

#include <iostream>
#include <string>
#include <limits>
#include <cmath>

namespace quadpack_cpp
{

typedef double (*QUADPACK_CPP_FUNCTION)(const double x, void *data);

typedef double (*QUADPACK_CPP_WEIGHT_FUNCTION)(const double x, 
                                               const double a, 
                                               const double b, 
                                               const double c, 
                                               const double d, 
                                               const int i);

void dqag(QUADPACK_CPP_FUNCTION f, 
          const double a, 
          const double b, 
          const double epsabs, 
          const double epsrel,
          const int key, 
          double &result,
          double &abserr, 
          int &neval, 
          int &ier,
          const int limit, 
          const int lenw, 
          int &last, 
          int *iwork,
          double *work,
          void *data);

void dqage(QUADPACK_CPP_FUNCTION f, 
           const double a, 
           const double b, 
           const double epsabs, 
           const double epsrel,
           const int key, 
           const int limit, 
           double &result, 
           double &abserr, 
           int &neval,
           int &ier, 
           double *alist, 
           double *blist,
           double *rlist, 
           double *elist, 
           int *iord, 
           int &last,
           void *data);

void dqagi(QUADPACK_CPP_FUNCTION f,
           const double bound,
           const int inf,
           const double epsabs,
           const double epsrel,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int limit,
           const int lenw,
           int &last,
           int *iwork,
           double *work,
           void *data);

void dqagie(QUADPACK_CPP_FUNCTION f,
            const double bound,
            const int inf,
            const double epsabs,
            const double epsrel,
            const int limit,
            double &result,
            double &abserr,
            int &neval,
            int &ier,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            int *iord,
            int &last,
            void *data);

void dqagp(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           const int npts2,
           const double *points,
           const double epsabs,
           const double epsrel,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int leniw,
           const int lenw,
           int &last,
           int *iwork,
           double *work,
           void *data);

void dqagpe(QUADPACK_CPP_FUNCTION f,
            double a,
            double b,
            const int npts2,
            const double *points,
            const double epsabs,
            const double epsrel,
            const int limit,
            double &result,
            double &abserr,
            int &neval,
            int &ier,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            double *pts,
            int *iord,
            int *level,
            int *ndin,
            int &last,
            void *data);

void dqags(QUADPACK_CPP_FUNCTION f, 
           const double a, 
           const double b, 
           const double epsabs, 
           const double epsrel, 
           double &result,
           double &abserr, 
           int &neval, 
           int &ier,
           const int limit, 
           const int lenw, 
           int &last, 
           int *iwork,
           double *work,
           void *data);


void dqagse(QUADPACK_CPP_FUNCTION f, 
            const double a, 
            const double b, 
            const double epsabs, 
            const double epsrel, 
            const int limit, 
            double &result, 
            double &abserr, 
            int &neval,
            int &ier, 
            double *alist, 
            double *blist,
            double *rlist, 
            double *elist, 
            int *iord, 
            int &last,
            void *data);


void dqawc(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           const double c,
           const double epsabs,
           const double epsrel,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int limit,
           const int lenw,
           int &last,
           int *iwork,
           double *work,
           void *data);

void dqawce(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double c,
            const double epsabs,
            const double epsrel,
            const int limit,
            double &result,
            double &abserr,
            int &neval,
            int &ier,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            int *iord,
            int &last,
            void *data);

void dqc25c(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double c,
            double &result,
            double &abserr,
            int &krul,
            int &neval,
            void *data);

void dqcheb(const double *x,
            double *fval,
            double *cheb12,
            double *cheb24);

void dqelg(int &n,
           double *epstab,
           double &result,
           double &abserr,
           double *res3la,
           int &nres);

void dqk15(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk15i(QUADPACK_CPP_FUNCTION f,
            const double boun,
            const int inf,
            const double a,
            const double b,
            double &result,
            double &abserr,
            double &resabs,
            double &resasc,
            void *data);

void dqk15w(QUADPACK_CPP_FUNCTION f,
            QUADPACK_CPP_WEIGHT_FUNCTION w,
            const double p1,
            const double p2,
            const double p3,
            const double p4,
            const int kp,
            const double a,
            const double b,
            double &result,
            double &abserr,
            double &resabs,
            double &resasc,
            void *data);

void dqk21(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk31(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk41(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk51(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk61(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqpsrt(const int limit,
            const int last,
            int &maxerr,
            double &ermax,
            const double *elist,
            int *iord,
            int &nrmax);

double dqwgtc(const double x,
              const double c,
              const double p2,
              const double p3,
              const double p4,
              const int kp);

void xerror(const std::string messg, 
           const int nerr, 
           const int level);

};

#endif /* QUADPACK_HPP_ */