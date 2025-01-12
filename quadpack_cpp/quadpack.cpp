#include "quadpack.hpp"

namespace quadpack_cpp
{

namespace {
    const double d1mach[] = {
        std::numeric_limits<double>::min(),
        std::numeric_limits<double>::max(),
        std::pow(static_cast<double>(std::numeric_limits<double>::radix), -std::numeric_limits<double>::digits),
        std::numeric_limits<double>::epsilon(),
        std::log10(static_cast<double>(std::numeric_limits<double>::radix))
    };
    const int limexp = 50;
}

// C++ array index --> Fortran array index
#define ARRAYF(a, i)         (a[(i)-1])
// C++ array index --> Fortran 2d-array index
#define MATF(a, rows,  i, j) (a[((i) - 1) * (rows) + ((j) - 1)]) 

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
          void *data)
{
    int lvl, l1, l2, l3;
    //
    // check validity of lenw.
    //
    // first executable statement  dqag
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (limit < 1 || lenw < limit*4) goto label_10;
    //
    // prepare call for dqage.
    //
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    //
    dqage(f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, ier,
        &ARRAYF(work, 1), &ARRAYF(work, l1), &ARRAYF(work, l2), &ARRAYF(work, l3), iwork, last, data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
    label_10:
    if (ier == 6) lvl = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqag";
        xerror(messg, ier,lvl);
    }
    return;
}

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
           void *data)
{
    double area1, a1, b1, defab1, error1;
    double area2, a2, b2, defab2, error2;
    double area;
    double area12;
    double erro12;
    double errsum;
    double errmax;
    double errbnd;
    int maxerr;
    double resabs, defabs;
    int iroff1, iroff2, k, keyf, nrmax;
    double epmach, uflow;
    // first executable statement  dqage
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    // test on validity of parameters
    // ------------------------------
    //
    ier    = 0;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    if (epsabs <= 0.0 && epsrel < std::max(50.0*epmach, 5.0e-29)) ier = 6;
    if (ier == 6) goto label_999;
    //
    // first approximation to the integral
    // -----------------------------------
    //
    keyf = key;
    if (key <= 0) keyf = 1;
    if (key >= 7) keyf = 6;
    neval = 0;
    if (keyf == 1) dqk15(f,a,b,result,abserr,defabs,resabs, data);
    if (keyf == 2) dqk21(f,a,b,result,abserr,defabs,resabs, data);
    if (keyf == 3) dqk31(f,a,b,result,abserr,defabs,resabs, data);
    if (keyf == 4) dqk41(f,a,b,result,abserr,defabs,resabs, data);
    if (keyf == 5) dqk51(f,a,b,result,abserr,defabs,resabs, data);
    if (keyf == 6) dqk61(f,a,b,result,abserr,defabs,resabs, data);
    last = 1;
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord, 1)  = 1;
    //
    // test on accuracy.
    //
    errbnd = std::max(epsabs, epsrel*std::abs(result));
    if (abserr <= 50.0*epmach*defabs && abserr > errbnd) ier = 2;
    if(limit == 1) ier = 1;
    if(ier != 0 || (abserr <= errbnd && abserr != resabs) || abserr == 0.0) goto label_60;
    //
    // initialization
    // --------------
    //
    errmax = abserr;
    maxerr = 1;
    area   = result;
    errsum = abserr;
    nrmax  = 1;
    iroff1 = 0;
    iroff2 = 0;
    //
    // main for-loop
    // -------------
    //
    for (last = 2; last <= limit; ++last) {
        //
        // bisect the subinterval with the largest error estimate.
        //
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2 = b1;
        b2 = ARRAYF(blist, maxerr);
        if (keyf == 1) dqk15(f, a1, b1, area1, error1, resabs, defab1, data);
        if (keyf == 2) dqk21(f, a1, b1, area1, error1, resabs, defab1, data);
        if (keyf == 3) dqk31(f, a1, b1, area1, error1, resabs, defab1, data);
        if (keyf == 4) dqk41(f, a1, b1, area1, error1, resabs, defab1, data);
        if (keyf == 5) dqk51(f, a1, b1, area1, error1, resabs, defab1, data);
        if (keyf == 6) dqk61(f, a1, b1, area1, error1, resabs, defab1, data);
        if (keyf == 1) dqk15(f, a2, b2, area2, error2, resabs, defab2, data);
        if (keyf == 2) dqk21(f, a2, b2, area2, error2, resabs, defab2, data);
        if (keyf == 3) dqk31(f, a2, b2, area2, error2, resabs, defab2, data);
        if (keyf == 4) dqk41(f, a2, b2, area2, error2, resabs, defab2, data);
        if (keyf == 5) dqk51(f, a2, b2, area2, error2, resabs, defab2, data);
        if (keyf == 6) dqk61(f, a2, b2, area2, error2, resabs, defab2, data);
        //
        // improve previous approximations to integral
        // and error and test for accuracy.
        neval += 1;
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto label_5;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) <= 1.0e-5*std::abs(area12) && erro12 >= 0.99*errmax) iroff1 += 1;
        if (last > 10 && erro12 > errmax) iroff2 += 1;
        label_5:
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist, last)   = area2;
        errbnd = std::max(epsabs, epsrel*std::abs(area));
        if (errsum <= errbnd) goto label_8;
        //
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 >= 6 || iroff2 >= 20) ier = 2;
        //
        // set error flag in the case that the number of subintervals
        // equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour
        // at a point of the integration range.
        //
        if(std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2*epmach)*(std::abs(a2) + 1.0e3*uflow)) ier = 3;
        //
        // append the newly-created intervals to the list.
        //
        label_8:
        if (error2 > error1) goto label_10;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto label_20;
        label_10:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with the largest error estimate (to be bisected next).
        //
        label_20:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if(ier != 0.0 || errsum <= errbnd) goto label_40;
    }
    //
    // compute final result.
    // ---------------------
    //
    label_40:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
    label_60:
    if (keyf != 1) neval = (10*keyf + 1)*(2*neval + 1);
    if (keyf == 1) neval = 30*neval + 15;
    label_999:
    return;
}


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
           void *data)
{
    int lvl, l1, l2, l3;
    // first executable statement  dqagi
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if(limit < 1 || lenw < limit*4) goto label_10;
    //
    // prepare call for dqagie.
    //
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    //
    dqagie(f, bound, inf, epsabs, epsrel, limit, result, abserr, neval, ier,
        &ARRAYF(work, 1), &ARRAYF(work, l1), &ARRAYF(work, l2), &ARRAYF(work, l3), iwork, last, data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
    label_10:
    if (ier == 6) lvl = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqagi";
        xerror(messg, ier, lvl);   
    }
    return;
}

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
            void *data)
{
    double area1, a1, b1, defab1, error1;
    double area2, a2, b2, defab2, error2;
    double area12;
    double erro12;
    double errmax;
    double erlast;
    double area;
    double errsum;
    double errbnd;
    double small;
    double erlarg;
    int maxerr;
    int nres;
    int numrl2;
    bool extrap;
    bool noext;
    double rlist2[52];
    double abseps, boun, correc, defabs, dres, ertest, resabs, reseps, res3la[3];
    int id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, nrmax;
    double epmach, uflow, oflow;
    // first executable statement  dqagie
    epmach = ARRAYF(d1mach, 4);
    //
    // test on validity of parameters
    // -----------------------------
    //
    ier    = 0;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    ARRAYF(alist, 1) = 0.0;
    ARRAYF(blist, 1) = 1.0;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    ARRAYF(iord,  1) = 0;
    if (epsabs <= 0.0 && epsrel < std::max(50.0*epmach, 5.0e-29)) ier = 6;
    if (ier == 6) goto label_999;
    //
    //
    // first approximation to the integral
    // -----------------------------------
    //
    // determine the interval to be mapped onto (0,1).
    // if inf = 2 the integral is computed as i = i1+i2, where
    // i1 = integral of f over (-infinity,0),
    // i2 = integral of f over (0,+infinity).
    //
    boun = bound;
    if (inf == 2) boun = 0.0;
    dqk15i(f, boun, inf, 0.0, 1.0, result, abserr, defabs, resabs, data);
    //
    // test on accuracy
    //
    last = 1;
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord,  1) = 1;
    dres = std::abs(result);
    errbnd = std::max(epsabs, epsrel*dres);
    if (abserr <= 1.0e2*epmach*defabs && abserr > errbnd) ier = 2;
    if (limit == 1) ier = 1;
    if (ier != 0 || (abserr <= errbnd && abserr != resabs) || abserr == 0.0) goto label_130;
    //
    // initialization
    // --------------
    //
    uflow = ARRAYF(d1mach, 1);
    oflow = ARRAYF(d1mach, 2);
    ARRAYF(rlist2, 1) = result;
    errmax = abserr;
    maxerr = 1;
    area = result;
    errsum = abserr;
    abserr = oflow;
    nrmax = 1;
    nres = 0;
    ktmin = 0;
    numrl2 = 2;
    extrap = false;
    noext = false;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if(dres >= (1.0 - 50.0*epmach)*defabs) ksgn = 1;
    //
    // main for-loop
    // -------------
    //
    for (last = 2; last <= limit; ++last) {
        //
        // bisect the subinterval with the nrmax-th largest error estimate.
        //
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5*(ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2 = b1;
        b2 = ARRAYF(blist, maxerr);
        erlast = errmax;
        dqk15i(f, boun, inf, a1, b1, area1, error1, resabs, defab1, data);
        dqk15i(f, boun, inf, a2, b2, area2, error2, resabs, defab2, data);
        // 
        // improve previous approximations to integral 
        // and error and test for accuracy.
        //
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto label_15;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) > 1.0e-5*std::abs(area12) || erro12 < 0.99*errmax) goto label_10;
        if (extrap)  iroff2 += 1;
        if (!extrap) iroff1 += 1;
        label_10:
        if (last > 10 && erro12 > errmax) iroff3 += 1;
        label_15:
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist, last)   = area2;
        errbnd = std::max(epsabs, epsrel*std::abs(area));
        // 
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
        //
        // set error flag in the case that the number of subintervals 
        // equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour 
        // at a point of the integration range.
        //
        if (std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2*epmach)*(std::abs(a2) + 1.0e3*uflow)) ier = 4;
        //
        // append the newly-created intervals to the list.
        //
        if (error2 > error1) goto label_20;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto label_30;
        label_20:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to be bisected next).
        //
        label_30:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        if (errsum <= errbnd) goto label_115;
        if (ier != 0)  goto label_100;
        if (last == 2) goto label_80;
        if (noext) goto label_90;
        erlarg -= erlast;
        if (std::abs(b1 - a1) > small) erlarg += erro12;
        if (extrap) goto label_40;
        // 
        // test whether the interval to be bisected next is the
        // smallest interval.
        // 
        if(std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto label_90;
        extrap = true;
        nrmax  = 2;
        label_40:
        if(ierro == 3 || erlarg <= ertest) goto label_60;
        //
        // the smallest interval has the largest error.
        // before bisecting decrease the sum of the errors over the
        // larger intervals (erlarg) and perform extrapolation.
        // 
        id     = nrmax;
        jupbnd = last;
        if (last > (2 + limit/2)) jupbnd = limit + 3 - last;
        for (k = id; k <= jupbnd; ++k) {
            maxerr = ARRAYF(iord, nrmax);
            errmax = ARRAYF(elist, maxerr);
            if (std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto label_90;
            nrmax += 1;
        }
        //
        // perform extrapolation.
        // 
        label_60:
        numrl2 += 1;
        ARRAYF(rlist2, numrl2) = area;
        dqelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin += 1;
        if(ktmin > 5 && abserr < 1.0e-3*errsum) ier = 5;
        if(abseps >= abserr) goto label_70;
        ktmin  = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel*std::abs(reseps));
        if (abserr <= ertest) goto label_100;
        //
        // prepare bisection of the smallest interval.
        // 
        label_70:
        if (numrl2 == 1) noext = true;
        if (ier == 5) goto label_100;
        maxerr = ARRAYF(iord, 1);
        errmax = ARRAYF(elist, maxerr);
        nrmax  = 1;
        extrap = false;
        small  = small*0.5;
        erlarg = errsum;
        goto label_90;
        label_80:
        small  = 0.375;
        erlarg = errsum;
        ertest = errbnd;
        ARRAYF(rlist2, 2) = area;
        label_90:
        continue;
    }
    //
    // set final result and error estimate.
    //
    label_100:
    if (abserr == oflow) goto label_115;
    if ((ier + ierro) == 0) goto label_110;
    if (ierro == 3) abserr += correc;
    if (ier == 0) ier = 3;
    if (result != 0.0 && area != 0.0) goto label_105;
    if (abserr > errsum) goto label_115;
    if (area == 0.0) goto label_130;
    goto label_110;
    label_105:
    if (abserr/std::abs(result) > errsum/std::abs(area)) goto label_115;
    //
    // test on divergence.
    //
    label_110:
    if(ksgn == -1 && std::max(std::abs(result), std::abs(area)) <= defabs*1.0e-2) goto label_130;
    if(1.0e-2 > (result/area) || (result/area) > 1.0e2 || errsum > std::abs(area)) ier = 6;
    goto label_130;
    // 
    // compute global integral sum.
    //
    label_115:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
    label_130:
    neval = 30*last - 15;
    if (inf == 2)  neval = 2*neval;
    if (ier > 2) ier -= 1;
    label_999:
    return;
}

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
           void *data)
{
    int limit, lvl, l1, l2, l3, l4;
    // first executable statement  dqagp
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (leniw < (3*npts2 - 2) || lenw < (leniw*2 - npts2) || npts2 < 2) goto label_10;
    //
    // prepare call for dqagpe.
    //
    limit = (leniw - npts2) / 2;
    l1    = limit + 1;
    l2    = limit + l1;
    l3    = limit + l2;
    l4    = limit + l3;
    //
    dqagpe(f , a, b, npts2, points, epsabs, epsrel, limit, result, abserr, 
        neval, ier, &ARRAYF(work, 1), &ARRAYF(work, l1), &ARRAYF(work, l2), &ARRAYF(work, l3), &ARRAYF(work, l4),
        &ARRAYF(iwork, 1), &ARRAYF(iwork, l1), &ARRAYF(iwork, l2), last, data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
    label_10:
    if (ier == 6) lvl = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqagp";
        xerror(messg, ier, lvl);   
    }
    return;
}

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
            void *data)
{
    double abseps, correc, defabs, dres, ertest, resa, reseps, Result,
           res3la[3], sign, temp, resabs;
    int i, id, ierro, ind1, ind2, ip1, iroff1, iroff2, iroff3, j, jlow, 
        jupbnd, k, ksgn, ktmin, levcur, levmax, nint, nintp1, npts, nrmax;
    double area1, a1, b1, defab1, error1;
    double area2, a2, b2, defab2, error2;
    double area12;
    double erro12;
    double rlist2[52];
    double erlast;
    double errsum;
    double errbnd;
    double area;
    double erlarg;
    double errmax;
    bool extrap;
    bool noext;
    int maxerr;
    int nres;
    int numrl2;
    double epmach, uflow, oflow;
    // first executable statement  dqagpe
    epmach = ARRAYF(d1mach, 4);
    //
    // test on validity of parameters
    // -----------------------------
    //
    ier    = 0;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    ARRAYF(iord,  1) = 0;
    ARRAYF(level, 1) = 0;
    npts = npts2 - 2;
    if (npts2 < 2 || limit <= npts || (epsabs <= 0.0 && epsrel < std::max(50.0*epmach, 5.0e-29))) ier = 6;
    if (ier == 6) goto label_999;
    //
    // if any break points are provided, sort them into an
    // ascending sequence.
    //
    sign = 1.0;
    if (a > b) sign = -1.0;
    ARRAYF(pts, 1) = std::min(a, b);
    if (npts == 0) goto label_15;
    for (i = 1; i <= npts; ++i) {
        ARRAYF(pts, i+1) = ARRAYF(points, i);
    }
    label_15:
    ARRAYF(pts, npts + 2) = std::max(a, b);
    nint = npts + 1;
    a1 = ARRAYF(pts, 1);
    if (npts == 0) goto label_40;
    nintp1 = nint + 1;
    for (i = 1; i <= nint; ++i) {
        ip1 = i + 1;
        for (j = ip1; j <= nintp1; ++j) {
          if (ARRAYF(pts, i) <= ARRAYF(pts, j)) continue;
          temp           = ARRAYF(pts, i);
          ARRAYF(pts, i) = ARRAYF(pts, j);
          ARRAYF(pts, j) = temp;
        }
    }
    if (ARRAYF(pts, 1) != std::min(a, b) || ARRAYF(pts, nintp1) != std::max(a, b)) ier = 6;
    if (ier == 6) goto label_999;
    //
    // compute first integral and error approximations.
    // ------------------------------------------------
    //
    label_40:
    resabs = 0.0;
    for (i = 1; i <= nint; ++i) {
        b1 = ARRAYF(pts, i+1);
        dqk21(f, a1, b1, area1, error1, defabs, resa, data);
        abserr = abserr + error1;
        result = result + area1;
        ARRAYF(ndin, i) = 0;
        if(error1 == resa && error1 != 0.0) ARRAYF(ndin, i) = 1;
        resabs = resabs + defabs;
        ARRAYF(level, i) = 0;
        ARRAYF(elist, i) = error1;
        ARRAYF(alist, i) = a1;
        ARRAYF(blist, i) = b1;
        ARRAYF(rlist, i) = area1;
        ARRAYF(iord,  i) = i;
        a1 = b1;
    }
    errsum = 0.0;
    for (i = 1; i <= nint; ++i) {
        if (ARRAYF(ndin, i) == 1) ARRAYF(elist, i) = abserr;
        errsum = errsum + ARRAYF(elist, i);
    }
    //
    // test on accuracy.
    //
    last   = nint;
    neval  = 21 * nint;
    dres   = std::abs(result);
    errbnd = std::max(epsabs, epsrel*dres);
    if (abserr <= 1.0e2*epmach*resabs && abserr > errbnd) ier = 2;
    if (nint == 1) goto label_80;
    for (i = 1; i <= npts; ++i) {
        jlow = i + 1;
        ind1 = ARRAYF(iord, i);
        for (j = jlow; j <= nint; ++j) {
            ind2 = ARRAYF(iord, j);
            if (ARRAYF(elist, ind1) > ARRAYF(elist, ind2)) continue;
            ind1 = ind2;
            k = j;
        }
        if (ind1 == ARRAYF(iord, i)) continue;
        ARRAYF(iord, k) = ARRAYF(iord, i);
        ARRAYF(iord, i) = ind1;
    }
    if (limit < npts2) ier = 1;
    label_80:
    if (ier != 0 || abserr <= errbnd) goto label_210;
    //
    // initialization
    // --------------
    //
    ARRAYF(rlist2, 1) = result;
    maxerr = ARRAYF(iord, 1);
    errmax = ARRAYF(elist, maxerr);
    area = result;
    nrmax = 1;
    nres = 0;
    numrl2 = 1;
    ktmin = 0;
    extrap = false;
    noext = false;
    erlarg = errsum;
    ertest = errbnd;
    levmax = 1;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ierro = 0;
    uflow = ARRAYF(d1mach, 1);
    oflow = ARRAYF(d1mach, 2);
    abserr = oflow;
    ksgn = -1;
    if (dres >= (1.0 - 50.0*epmach)*resabs) ksgn = 1;
    //
    // main for-loop
    // ------------
    //
    for (last = npts2; last <= limit; ++last) {
        //
        // bisect the subinterval with the nrmax-th largest error
        // estimate.
        //
        levcur = ARRAYF(level, maxerr) + 1;
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2 = b1;
        b2 = ARRAYF(blist, maxerr);
        erlast = errmax;
        dqk21(f, a1, b1, area1, error1, resa, defab1, data);
        dqk21(f, a2, b2, area2, error2, resa, defab2, data);
        //
        // improve previous approximations to integral
        // and error and test for accuracy.
        //
        neval  = neval + 42;
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto label_95;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) > 1.0e-5*std::abs(area12) || erro12 < 0.99*errmax) goto label_90;
        if (extrap)  iroff2 = iroff2 + 1;
        if (!extrap) iroff1 = iroff1 + 1;
        label_90:
        if (last > 10 && erro12 > errmax) iroff3 = iroff3 + 1;
        label_95:
        ARRAYF(level, maxerr) = levcur;
        ARRAYF(level, last)   = levcur;
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist, last)   = area2;
        errbnd = std::max(epsabs, epsrel*std::abs(area));
        //
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
        //
        // set error flag in the case that the number of
        // subintervals equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour
        // at a point of the integration range
        //
        if (std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2*epmach)*(std::abs(a2) + 1.03*uflow)) ier = 4;
        //
        // append the newly-created intervals to the list.
        //
        if(error2 > error1) goto label_100;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto label_110;
        label_100:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to be bisected next).
        //
        label_110:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if (errsum <= errbnd) goto label_190;
        // jump out of for-loop
        if (ier != 0) goto label_170;
        if (noext) goto label_160;
        erlarg = erlarg - erlast;
        if (levcur + 1 <= levmax) erlarg = erlarg + erro12;
        if (extrap) goto label_120;
        //
        // test whether the interval to be bisected next is the
        // smallest interval.
        //
        if (ARRAYF(level, maxerr) + 1 <= levmax) goto label_160;
        extrap = true;
        nrmax = 2;
        label_120:
        if (ierro == 3 || erlarg <= ertest) goto label_140;
        //
        // the smallest interval has the largest error.
        // before bisecting decrease the sum of the errors over
        // the larger intervals (erlarg) and perform extrapolation.
        //
        id     = nrmax;
        jupbnd = last;
        if (last > (2 + limit/2)) jupbnd = limit + 3 - last;
        for (k = id; k <= jupbnd; ++k) {
            maxerr = ARRAYF(iord, nrmax);
            errmax = ARRAYF(elist, maxerr);
            // jump out of do-loop
            if (ARRAYF(level, maxerr) + 1 <= levmax) goto label_160;
            nrmax = nrmax+1;
        }
        //
        // perform extrapolation.
        //
        label_140:
        numrl2 = numrl2 + 1;
        ARRAYF(rlist2, numrl2) = area;
        if (numrl2 <= 2) goto label_155;
        dqelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin = ktmin + 1;
        if (ktmin > 5 && abserr < 1.0e-3*errsum) ier = 5;
        if (abseps >= abserr) goto label_150;
        ktmin = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel*std::abs(reseps));
        // jump out of do-loop
        if (abserr < ertest) goto label_170;
        //
        // prepare bisection of the smallest interval.
        //
        label_150:
        if (numrl2 == 1) noext = true;
        if(ier >= 5) goto label_170;
        label_155:
        maxerr = ARRAYF(iord, 1);
        errmax = ARRAYF(elist, maxerr);
        nrmax = 1;
        extrap = false;
        levmax = levmax + 1;
        erlarg = errsum;
        label_160:
        continue;
    }
    //
    // set the final result.
    // ---------------------
    //
    //
    label_170:
    if (abserr == oflow) goto label_190;
    if ((ier + ierro) == 0) goto label_180;
    if (ierro == 3) abserr = abserr + correc;
    if (ier == 0) ier = 3;
    if (result != 0.0 && area != 0.0) goto label_175;
    if (abserr > errsum) goto label_190;
    if (area == 0.0) goto label_210;
    goto label_180;
    label_175:
    if (abserr/std::abs(result) > errsum/std::abs(area)) goto label_190;
    //
    // test on divergence.
    //
    label_180:
    if(ksgn == -1 && std::max(std::abs(result), std::abs(area)) <= resabs*1.0e-2) goto label_210;
    if (1.0e-2 > result/area || result/area > 1.0e2 || errsum > std::abs(area)) ier = 6;
    goto label_210;
    //
    // compute global integral sum.
    //
    label_190:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
    label_210:
    if (ier > 2) ier -= 1;
    result *= sign;
    label_999:
    return;
}

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
           void *data)
{
    int lvl, l1, l2, l3;
    //
    // check validity of limit and lenw.
    //
    // first executable statement  dqags
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (limit < 1 || lenw < 4*limit) goto label_10;
    //
    // prepare call for dqagse.
    //
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    //
    dqagse(f, a, b, epsabs, epsrel, limit, result, abserr, neval, ier,
        &ARRAYF(work, 1), &ARRAYF(work, l1), &ARRAYF(work, l2), &ARRAYF(work, l3), iwork, last, data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
    label_10:
    if (ier == 6) lvl = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqags";
        xerror(messg, ier, lvl);
    }
    return;
}


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
            void *data) 
{
    double abseps, correc, defabs, dres, ertest, resabs, reseps, res3la[3];
    double uflow, oflow, epmach;
    int id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, nrmax;
    double area12;
    double erro12;
    double area1, a1, b1, defab1, error1;
    double area2, a2, b2, defab2, error2;
    double rlist2[limexp + 2];
    int maxerr;
    int nres;
    int numrl2;
    double errmax;
    double erlast;
    double area;
    double errsum;
    double errbnd;
    double small;
    double erlarg;
    bool extrap;
    bool noext;
    // first executable statement  dqagse
    epmach = ARRAYF(d1mach, 4);
    //
    // test on validity of parameters
    // ------------------------------
    //
    ier    = 0;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    if (epsabs <= 0.0 && epsrel < std::max(50.0*epmach, 5.0e-29)) ier = 6;
    if (ier == 6) goto label_999;
    // 
    // first approximation to the integral
    // -----------------------------------
    //
    uflow = ARRAYF(d1mach, 1);
    oflow = ARRAYF(d1mach, 2);
    ierro = 0;
    dqk21(f, a, b, result, abserr, defabs, resabs, data);
    // 
    // test on accuracy
    // 
    dres   = std::abs(result);
    errbnd = std::max(epsabs, epsrel*dres);
    last   = 1;
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord,  1) = 1;
    if (abserr <= 1.0e2*epmach*defabs && abserr > errbnd) ier = 2;
    if (limit == 1) ier = 1;
    if (ier != 0 || (abserr <= errbnd && abserr != resabs) || abserr == 0.0) goto label_140;
    //
    // initialization
    // --------------
    //
    ARRAYF(rlist2, 1) = result;
    errmax = abserr;
    maxerr = 1;
    area   = result;
    errsum = abserr;
    abserr = oflow;
    nrmax  = 1;
    nres   = 0;
    numrl2 = 2;
    ktmin  = 0;
    extrap = false;
    noext  = false;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn   = -1;
    if (dres >= (1.0 - 50.0*epmach)*defabs) ksgn = 1;
    //
    // main for-loop
    // ---------
    //
    for (last = 2; last <= limit; ++last) { // label_90
        //
        // bisect the subinterval with the nrmax-th largest error 
        // estimate.
        //
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5*(ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2 = b1;
        b2 = ARRAYF(blist, maxerr);
        erlast = errmax;
        dqk21(f, a1, b1, area1, error1, resabs, defab1, data);
        dqk21(f, a2, b2, area2, error2, resabs, defab2, data);
        // 
        // improve previous approximations to integral 
        // and error and test for accuracy.
        //
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto label_15;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) > 1.0e-5*std::abs(area12) || erro12 < 0.99*errmax) goto label_10;
        if (extrap)  iroff2 += 1;
        if (!extrap) iroff1 += 1;
        label_10:
        if (last > 10 && erro12 > errmax) iroff3 += 1;
        label_15:
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist, last)   = area2;
        errbnd = std::max(epsabs, epsrel*std::abs(area));
        // 
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
        //
        // set error flag in the case that the number of subintervals 
        // equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour 
        // at a point of the integration range.
        //
        if (std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2*epmach)*(std::abs(a2) + 1.0e3*uflow)) ier = 4;
        //
        // append the newly-created intervals to the list.
        //
        if (error2 > error1) goto label_20;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto label_30;
        label_20:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to be bisected next).
        //
        label_30:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if (errsum <= errbnd) goto label_115;
        // jump out of for-loop
        if (ier != 0)  goto label_100;
        if (last == 2) goto label_80;
        if (noext) goto label_90;
        erlarg -= erlast;
        if (std::abs(b1-a1) > small) erlarg += erro12;
        if (extrap) goto label_40;
        // 
        // test whether the interval to be bisected next is the
        // smallest interval.
        // 
        if(std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto label_90;
        extrap = true;
        nrmax  = 2;
        label_40:
        if(ierro == 3 || erlarg <= ertest) goto label_60;
        //
        // the smallest interval has the largest error.
        // before bisecting decrease the sum of the errors over the
        // larger intervals (erlarg) and perform extrapolation.
        // 
        id     = nrmax;
        jupbnd = last;
        if (last > (2 + limit/2)) jupbnd = limit + 3 - last;
        for (k = id; k <= jupbnd; ++k) {
            maxerr = ARRAYF(iord, nrmax);
            errmax = ARRAYF(elist, maxerr);
            // jump out of for-loop
            if (std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto label_90;
            nrmax += 1;
        }
        //
        // perform extrapolation.
        // 
        label_60:
        numrl2 += 1;
        ARRAYF(rlist2, numrl2) = area;
        dqelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin += 1;
        if(ktmin > 5 && abserr < 1.0e-3*errsum) ier = 5;
        if(abseps >= abserr) goto label_70;
        ktmin  = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel*std::abs(reseps));
        // jump out of for-loop
        if (abserr <= ertest) goto label_100;
        //
        // prepare bisection of the smallest interval.
        // 
        label_70:
        if (numrl2 == 1) noext = true;
        if (ier == 5) goto label_100;
        maxerr = ARRAYF(iord, 1);
        errmax = ARRAYF(elist, maxerr);
        nrmax  = 1;
        extrap = false;
        small  = small*0.5;
        erlarg = errsum;
        goto label_90;
        label_80:
        small  = std::abs(b - a)*0.375;
        erlarg = errsum;
        ertest = errbnd;
        ARRAYF(rlist2, 2) = area;
        label_90:
        continue;
    }
    //
    // set final result and error estimate.
    //
    label_100:
    if (abserr == oflow) goto label_115;
    if ((ier + ierro) == 0) goto label_110;
    if (ierro == 3) abserr += correc;
    if (ier == 0) ier = 3;
    if (result != 0.0 && area != 0.0) goto label_105;
    if (abserr > errsum) goto label_115;
    if (area == 0.0) goto label_130;
    goto label_110;
    label_105:
    if (abserr/std::abs(result) > errsum/std::abs(area)) goto label_115;
    //
    // test on divergence.
    //
    label_110:
    if(ksgn == -1 && std::max(std::abs(result), std::abs(area)) <= defabs*1.0e-2) goto label_130;
    if(1.0e-2 > (result/area) || (result/area) > 1.0e2 || errsum > std::abs(area)) ier = 6;
    goto label_130;
    // 
    // compute global integral sum.
    //
    label_115:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
    label_130:
    if(ier > 2) ier -= 1;
    label_140:
    neval = 42*last - 21;
    label_999:
    return;
}

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
           void *data)
{
    int lvl, l1, l2, l3;
    // first executable statement  dqawc
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if(limit < 1 || lenw < limit*4) goto label_10;
    //
    // prepare call for dqawce.
    //
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    dqawce(f, a, b, c, epsabs, epsrel, limit, result, abserr, neval, ier,
        &ARRAYF(work, 1), &ARRAYF(work, l1), &ARRAYF(work, l2), &ARRAYF(work, l3), iwork, last, data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
    label_10:
    if (ier == 6) ier = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqawc";
        xerror(messg, ier, lvl);
    }
    return;
}

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
            void *data)
{
    double aa, bb, epmach, uflow;
    int iroff1, iroff2, k, krule, nev, nrmax;
    double area1, a1, b1, error1;
    double area2, a2, b2, error2;
    double area12;
    double erro12;
    double errmax;
    double area;
    double errsum;
    double errbnd;
    int maxerr;
    // first executable statement  dqawce
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    //
    // test on validity of parameters
    // ------------------------------
    //
    ier = 6;
    neval = 0;
    last = 0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    ARRAYF(iord,  1) = 0;
    result = 0.0;
    abserr = 0.0;
    if(c == a || c== b || (epsabs <= 0.0 && epsrel < std::max(50.0*epmach, 5.0e-29))) goto label_999;
    //
    // first approximation to the integral
    // -----------------------------------
    //
    aa = a;
    bb = b;
    if (a <= b) goto label_10;
    aa = b;
    bb = a;
    label_10:
    ier = 0;
    krule = 1;
    dqc25c(f, aa, bb, c, result, abserr, krule, neval, data);
    last = 1;
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord,  1) = 1;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    //
    // test on accuracy
    //
    errbnd = std::max(epsabs, epsrel*std::abs(result));
    if(limit == 1) ier = 1;
    if(abserr < std::min(1.0e-2*std::abs(result), errbnd) || ier == 1) goto label_70;
    //
    // initialization
    // --------------
    //
    ARRAYF(alist, 1) = aa;
    ARRAYF(blist, 1) = bb;
    ARRAYF(rlist, 1) = result;
    errmax = abserr;
    maxerr = 1;
    area = result;
    errsum = abserr;
    nrmax = 1;
    iroff1 = 0;
    iroff2 = 0;
    //
    // main for-loop
    // ------------
    //
    for (last = 2; last <= limit; ++last) {
        //
        // bisect the subinterval with nrmax-th largest
        // error estimate.
        //
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        b2 = ARRAYF(blist, maxerr);
        if (c <= b1 && c > a1) b1 = 0.5 * (c + b2);
        if (c > b1  && c < b2) b1 = 0.5 * (a1 + c);
        a2 = b1;
        krule = 2;
        dqc25c(f, a1, b1, c, area1, error1, krule, nev, data);
        neval = neval + nev;
        dqc25c(f, a2, b2, c, area2, error2, krule, nev, data);
        neval = neval + nev;
        //
        // improve previous approximations to integral
        // and error and test for accuracy.
        //
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area + area12 - ARRAYF(rlist, maxerr);
        if (std::abs(ARRAYF(rlist, maxerr) - area12) < 1.0e-5*std::abs(area12) 
           && erro12 >= 0.99*errmax && krule == 0) iroff1 = iroff1 + 1;
        if (last > 10 && erro12 > errmax && krule == 0) iroff2 = iroff2 + 1;
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist, last)   = area2;
        errbnd = std::max(epsabs, epsrel*std::abs(area));
        if (errsum <= errbnd) goto label_15;
        //
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 >= 6 && iroff2 > 20) ier = 2;
        //
        // set error flag in the case of bad integrand behaviour
        // at a point of the integration range.
        //
        if (std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2*epmach)*(std::abs(a2) + 1.0e3*uflow)) ier = 3;
        //
        // append the newly-created intervals to the list.
        //
        label_15:
        if(error2 > error1) goto label_20;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto label_30;
        label_20:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to be bisected next).
        //
        label_30:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if(ier != 0 || errsum <= errbnd) goto label_50;
    }
    //
    // compute final result.
    // ---------------------
    //
    label_50:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
    label_70:
    if (aa == b) result = -result;
    label_999:
    return;
}

void dqawf(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double omega,
           const int integr,
           const double epsabs,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int limlst,
           int &lst,
           const int leniw,
           const int maxp1,
           const int lenw,
           int *iwork,
           double *work,
           void *data)
{
    int last, limit, ll2, lvl, l1, l2, l3, l4, l5, l6;
    // first executable statement  dqawf
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (limlst < 3 || leniw < limlst + 2 || maxp1 < 1 || lenw < (leniw*2 + maxp1*25)) goto label_10;
    //
    // prepare call for dqawfe
    //
    limit = (leniw - limlst) / 2;
    l1    = limlst + 1;
    l2    = limlst + l1;
    l3    = limit  + l2;
    l4    = limit  + l3;
    l5    = limit  + l4;
    l6    = limit  + l5;
    ll2   = limit  + l1;
    dqawfe(f, a, omega, integr, epsabs, limlst, limit, maxp1, result, abserr, neval, ier,
        &ARRAYF(work,  1), &ARRAYF(work, l1), &ARRAYF(iwork, 1), lst, &ARRAYF(work, l2),
        &ARRAYF(work, l3), &ARRAYF(work, l4), &ARRAYF(work, l5), &ARRAYF(iwork, l1), &ARRAYF(iwork, ll2), &ARRAYF(work, l6),
        data);
    //
    // call error handler if necessary
    //
    lvl = 0;
    label_10:
    if (ier == 6) ier = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqawf";
        xerror(messg, ier, lvl);
    }
    return;
}

void dqawfe(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double omega,
            const int integr,
            const double epsabs,
            const int limlst,
            const int limit,
            const int maxp1,
            double &result,
            double &abserr,
            int &neval,
            int &ier,
            double *rslst,
            double *erlst,
            int *ierlst,
            int &lst,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            int *iord,
            int *nnlog,
            double *chebmo, // **chebmo
            void *data)
{
    double abseps, correc, dl, dla, drl, ep, eps, fact, p1, reseps, res3la[3];
    int ktmin, l, last, ll, momcom, nev, nres, numrl2;
    double psum[52];
    double c1, c2;
    double cycle;
    double errsum;
    double epsa;
    double uflow;
    const double p = 0.9;
    //
    // test on validity of parameters
    // ------------------------------
    //
    // first executable statement  dqawfe
    result = 0.0;
    abserr = 0.0;
    neval  = 0;
    lst    = 0;
    ier    = 0;
    if ((integr != 1 && integr != 2) || epsabs <= 0.0 || limlst < 3) ier = 6;
    if (ier == 6) goto label_999;
    if (omega != 0.0) goto label_10;
    //
    // integration by dqagie if omega is zero
    // --------------------------------------
    //
    if(integr == 1) dqagie(f, 0.0, 1, epsabs, 0.0,limit, result, abserr, neval, ier,
                        alist, blist, rlist, elist, iord, last, data);
    ARRAYF(rslst, 1)  = result;
    ARRAYF(erlst, 1)  = abserr;
    ARRAYF(ierlst, 1) = ier;
    lst = 1;
    goto label_999;
    //
    // initializations
    // ---------------
    //
    label_10:
    l = std::abs(omega);
    dl = 2*l + 1;
    cycle = dl*M_PI/std::abs(omega);
    ier = 0;
    ktmin = 0;
    neval = 0;
    numrl2 = 0;
    nres = 0;
    c1 = a;
    c2 = cycle + a;
    p1 = 1.0 - p;
    uflow = ARRAYF(d1mach, 1);
    eps = epsabs;
    if(epsabs > uflow/p1) eps = epsabs*p1;
    ep = eps;
    fact = 1.0;
    correc = 0.0;
    abserr = 0.0;
    errsum = 0.0;
    //
    // main for-loop
    // -------------
    //
    for (lst = 1; lst <= limlst; ++lst) {
        //
        // integrate over current subinterval.
        //
        dla = lst;
        epsa = eps*fact;
        dqawoe(f, c1, c2, omega, integr, epsa, 0.0, limit, lst, maxp1,
            ARRAYF(rslst, lst), ARRAYF(erlst, lst), nev, ARRAYF(ierlst, lst), 
            last, alist, blist, rlist, elist, iord, nnlog, momcom, chebmo, data);
        neval = neval + nev;
        fact = fact*p;
        errsum = errsum + ARRAYF(erlst, lst);
        drl = 50.0*std::abs(ARRAYF(rslst, lst));
        //
        // test on accuracy with partial sum
        //
        if (errsum + drl <= epsabs && lst >= 6) goto label_80;
        correc = std::max(correc, ARRAYF(erlst, lst));
        if (ARRAYF(ierlst, lst) != 0) eps = std::max(ep, correc*p1);
        if (ARRAYF(ierlst, lst) != 0) ier = 7;
        if (ier == 7 && errsum + drl <= correc*10.0 && lst > 5) goto label_80;
        numrl2 = numrl2 + 1;
        if (lst > 1) goto label_20;
        ARRAYF(psum, 1) = ARRAYF(rslst, 1);
        goto label_40;
        label_20:
        ARRAYF(psum, numrl2) = ARRAYF(psum, ll) + ARRAYF(rslst, lst);
        if(lst == 2) goto label_40;
        //
        // test on maximum number of subintervals
        //
        if (lst == limlst) ier = 1;
        //
        // perform new extrapolation
        //
        dqelg(numrl2, psum, reseps, abseps, res3la, nres);
        //
        // test whether extrapolated result is influenced by roundoff
        //
        ktmin = ktmin + 1;
        if (ktmin >= 15 && abserr <= 1.0e-3*(errsum + drl)) ier = 4;
        if (abseps > abserr && lst != 3) goto label_30;
        abserr = abseps;
        result = reseps;
        ktmin = 0;
        //
        // if ier is not 0, check whether direct result (partial sum)
        // or extrapolated result yields the best integral
        // approximation
        // 
        if ((abserr + 10.0*correc) <= epsabs || (abserr <= epsabs && 10.0*correc >= epsabs)) goto label_60;
        label_30:
        if (ier != 0 && ier != 7) goto label_60;
        label_40:
        ll = numrl2;
        c1 = c2;
        c2 = c2 + cycle;
    }
    //
    // set final result and error estimate
    // -----------------------------------
    //
    label_60:
    if (ier == 0) goto label_999;
    if (result != 0.0 && ARRAYF(psum, numrl2) != 0.0) goto label_70;
    if (abserr > errsum) goto label_80;
    if (ARRAYF(psum, numrl2) == 0.0) goto label_999;
    label_70:
    if (abserr/std::abs(result) > (errsum + drl)/std::abs(ARRAYF(psum, numrl2))) goto label_80;
    if (ier >= 1 && ier != 7) abserr = abserr + drl;
    goto label_999;
    label_80:
    result = ARRAYF(psum, numrl2);
    abserr = errsum + drl;
    label_999:
    return;
}

void dqawoe(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double omega,
            const int integr,
            const double epsabs,
            const double epsrel,
            const int limit,
            const int icall,
            const int maxp1,
            double &result,
            double &abserr,
            int  &neval,
            int &ier,
            int &last,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            int *iord,
            int *nnlog,
            int &momcom,
            double *chebmo, // ** chebmo
            void * data)
{
    double abseps, correc, defab1, defab2, defabs, domega, dres, ertest, resabs, reseps, res3la[3], width;
    int id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, nev, nrmax, nrmom;
    bool extall, done, test;
    double rlist2[52];
    int maxerr;
    double errmax;
    double erlast;
    double area;
    double errsum;
    double errbnd;
    double a1, area1, b1, error1;
    double a2, area2, b2, error2;
    int nres;
    int numrl2;
    double small;
    double erlarg;
    double area12;
    double erro12;
    bool extrap;
    bool noext;
    double epmach, uflow, oflow;
    // first executable statement  dqawoe
    epmach = ARRAYF(d1mach, 4);
    //
    // test on validity of parameters
    // ------------------------------
    //
    ier = 0;
    neval = 0;
    last = 0;
    result = 0.0;
    abserr = 0.0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    ARRAYF(iord,  1) = 0;
    ARRAYF(nnlog, 1) = 0;
    if((integr != 1 && integr != 2) || 
       (epsabs <= 0.0 && epsrel < std::max(50.0*epmach, 5.0e-29)) || 
       icall < 1 || maxp1 < 1) ier = 6;
    if(ier == 6) goto label_999;
    //
    // first approximation to the integral
    // -----------------------------------
    //
    domega = std::abs(omega);
    nrmom = 0;
    if (icall > 1) goto label_5;
    momcom = 0;
    label_5:
    dqc25f(f, a, b, domega, integr, nrmom, maxp1, 0, result, abserr,
        neval, defabs, resabs, momcom, chebmo, data);
    //
    // test on accuracy.
    //
    dres = std::abs(result);
    errbnd = std::max(epsabs, epsrel*dres);
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord,  1) = 1;
    if (abserr <= 1.0e2*epmach*defabs && abserr > errbnd) ier = 2;
    if (limit == 1) ier = 1;
    if (ier != 0 || abserr <= errbnd) goto label_200;
    //
    // initializations
    // ---------------
    //
    uflow = ARRAYF(d1mach, 1);
    oflow = ARRAYF(d1mach, 2);
    errmax = abserr;
    maxerr = 1;
    area = result;
    errsum = abserr;
    abserr = oflow;
    nrmax = 1;
    extrap = false;
    noext = false;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ktmin = 0;
    small = std::abs(b - a)*0.75;
    nres = 0;
    numrl2 = 0;
    extall = false;
    if (0.5*std::abs(b - a)*domega > 2.0) goto label_10;
    numrl2 = 1;
    extall = true;
    ARRAYF(rlist2, 1) = result;
    label_10:
    if (0.25*std::abs(b - a)*domega <= 2.0) extall = true;
    ksgn = -1;
    if (dres >= (1.0 - 50.0*epmach)*defabs) ksgn = 1;
    //
    // main for-loop
    // -------------
    //
    for (last = 2; last <= limit; ++last) {
        //
        // bisect the subinterval with the nrmax-th largest
        // error estimate.
        //
        nrmom = ARRAYF(nnlog, maxerr) + 1;
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2 = b1;
        b2 = ARRAYF(blist, maxerr);
        erlast = errmax;
        dqc25f(f, a1, b1, domega, integr, nrmom, maxp1, 0, 
            area1, error1, nev, resabs, defab1, momcom, chebmo, data);
        neval = neval + nev;
        dqc25f(f, a2, b2, domega, integr, nrmom, maxp1, 1,
            area2, error2, nev, resabs, defab2, momcom, chebmo, data);
        neval = neval + nev;
        //
        // improve previous approximations to integral
        // and error and test for accuracy.
        //
        area12 = area1  + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto label_25;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) > 1.0e-5*std::abs(area12) || erro12 < 0.99*errmax) goto label_20;
        if (extrap) iroff2 = iroff2 + 1;
        if (!extrap) iroff1 = iroff1 + 1;
        label_20:
        if (last > 10 && erro12 > errmax) iroff3 = iroff3 + 1;
        label_25:
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist,   last) = area2;
        ARRAYF(nnlog, maxerr) = nrmom;
        ARRAYF(nnlog,   last) = nrmom;
        errbnd = std::max(epsabs, epsrel*std::abs(area));
        //
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
        //
        // set error flag in the case that the number of
        // subintervals equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour
        // at a point of the integration range.
        //
        if(std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2*epmach)*(std::abs(a2) + 1.0e3*uflow)) ier = 4;
        //
        // append the newly-created intervals to the list.
        //
        if (error2 > error1) goto label_30;
        ARRAYF(alist,   last) = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist,   last) = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist,   last) = error2;
        goto label_40;
        label_30:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist,   last) = a1;
        ARRAYF(blist,   last) = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist,   last) = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist,   last) = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to bisected next).
        //
        label_40:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if (errsum <= errbnd) goto label_170;
        if (ier != 0) goto label_150;
        if (last == 2 && extall) goto label_120;
        if (noext) goto label_140;
        if (!extall) goto label_50;
        erlarg = erlarg - erlast;
        if (std::abs(b1 - a1) > small) erlarg = erlarg + erro12;
        if (extrap) goto label_70;
        //
        // test whether the interval to be bisected next is the
        // smallest interval.
        //
        label_50:
        width = std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr));
        if (width > small) goto label_140;
        if (extall) goto label_60;
        //
        // test whether we can start with the extrapolation procedure
        // (we do this if we integrate over the next interval with
        // use of a gauss-kronrod rule - see subroutine dqc25f).
        //
        small = small*0.5;
        if (0.25*width*domega > 2.0) goto label_140;
        extall = true;
        goto label_130;
        label_60:
        extrap = true;
        nrmax = 2;
        label_70:
        if (ierro == 3 || erlarg <= ertest) goto label_90;
        //
        // the smallest interval has the largest error.
        // before bisecting decrease the sum of the errors over
        // the larger intervals (erlarg) and perform extrapolation.
        //
        jupbnd = last;
        if (last > (limit/2 + 2)) jupbnd = limit + 3 - last;
        id = nrmax;
        for (k = id; k <= jupbnd; ++k) {
            maxerr = ARRAYF(iord,   nrmax);
            errmax = ARRAYF(elist, maxerr);
            if (std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto label_140;
            nrmax = nrmax + 1;
        }
        // 
        // perform extrapolation.
        //
        label_90:
        numrl2 = numrl2 + 1;
        ARRAYF(rlist2, numrl2) = area;
        if (numrl2 < 3) goto label_110;
        dqelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin = ktmin + 1;
        if (ktmin > 5 && abserr < 1.0e-3*errsum) ier = 5;
        if (abseps >= abserr) goto label_100;
        ktmin = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel*std::abs(reseps));
        // jump out of for-loop
        if (abserr <= ertest) goto label_150;
        //
        // prepare bisection of the smallest interval.
        //
        label_100:
        if (numrl2 == 1) noext = true;
        if (ier == 5) goto label_150;
        label_110:
        maxerr = ARRAYF(iord, 1);
        errmax = ARRAYF(elist, maxerr);
        nrmax = 1;
        extrap = false;
        small = small*0.5;
        erlarg = errsum;
        goto label_140;
        label_120:
        small = small*0.5;
        numrl2 = numrl2 + 1;
        ARRAYF(rlist2, numrl2) = area;
        label_130:
        ertest = errbnd;
        erlarg = errsum;
        label_140:
        continue;
    }
    //
    // set the final result.
    // ---------------------
    //
    label_150:
    if (abserr == oflow || nres == 0) goto label_170;
    if (ier + ierro == 0) goto label_165;
    if (ierro == 3) abserr = abserr + correc;
    if (ier == 0) ier = 3;
    if (result != 0.0 && area != 0.0) goto label_160;
    if (abserr > errsum) goto label_170;
    if (area == 0.0) goto label_190;
    goto label_165;
    label_160:
    if (abserr/std::abs(result) > errsum/std::abs(area)) goto label_170;
    //
    // test on divergence.
    //
    label_165:
    if (ksgn == -1 && std::max(std::abs(result), std::abs(area)) <= defabs*1.0e-2) goto label_190;
    if(1.0e-2 > result/area || result/area > 1.0e2 || errsum >= std::abs(area)) ier = 6;
    goto label_190;
    //
    // compute global integral sum.
    //
    label_170:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
    label_190:
    if (ier > 2) ier -= 1;
    label_200:
    if (integr == 2 && omega < 0.0) result = -result;
    label_999:
    return;
}

void dqc25c(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double c,
            double &result,
            double &abserr,
            int &krul,
            int &neval,
            void *data)
{
    double ak22, amom0, amom1, amom2, cc, p2, p3, p4, resabs, resasc, u;
    int i, isym, k;
    double fval[25];
    double cheb12[13];
    double cheb24[25];
    double res12;
    double res24;
    double hlgth;
    double centr;
    const int kp = 0;
    double x[11];
    for (k = 1; k <= 11; ++k) ARRAYF(x, k) = std::cos(k*M_PI/24.0);
    // first executable statement  dqc25c
    cc = (2.0*c - b - a) / (b - a);
    if (std::abs(cc) < 1.1) goto label_10;
    //
    // apply the 15-point gauss-kronrod scheme.
    //
    krul = krul - 1;
    dqk15w(f, dqwgtc, c, p2, p3, p4, kp, a, b, result, abserr, resabs, resasc, data);
    neval = 15;
    if (resasc == abserr) krul = krul + 1;
    goto label_50;
    //
    // use the generalized clenshaw-curtis method.
    //
    label_10:
    hlgth = 0.5 * (b - a);
    centr = 0.5 * (b + a);
    neval = 25;
    ARRAYF(fval,  1) = 0.5 * f(hlgth + centr, data);
    ARRAYF(fval, 13) = f(centr, data);
    ARRAYF(fval, 25) = 0.5 * f(centr - hlgth, data);
    for (i = 2; i <= 12; ++i) {
        u = hlgth * ARRAYF(x, i-1);
        isym = 26 - i;
        ARRAYF(fval, i)    = f(u + centr, data);
        ARRAYF(fval, isym) = f(centr - u, data);
    }
    //
    // compute the chebyshev series expansion.
    //
    dqcheb(x, fval, cheb12, cheb24);
    //
    // the modified chebyshev moments are computed by forward
    // recursion, using amom0 and amom1 as starting values.
    //
    amom0 = std::log(std::abs((1.0 - cc) / (1.0 + cc)));
    amom1 = 2.0 + cc*amom0;
    res12 = ARRAYF(cheb12, 1)*amom0 + ARRAYF(cheb12, 2)*amom1;
    res24 = ARRAYF(cheb24, 1)*amom0 + ARRAYF(cheb24, 2)*amom1;
    for (k = 3; k <= 13; ++k) {
        amom2 = 2.0*cc*amom1 - amom0;
        ak22 = (k - 2) * (k - 2);
        if ((k/2)*2 == k) amom2 = amom2 - 4.0/(ak22 - 1.0);
        res12 = res12 + ARRAYF(cheb12, k)*amom2;
        res24 = res24 + ARRAYF(cheb24, k)*amom2;
        amom0 = amom1;
        amom1 = amom2;
    }
    for (k = 14; k <= 25; ++k) {
        amom2 = 2.0*cc*amom1 - amom0;
        ak22 = (k - 2) * (k - 2);
        if ((k/2)*2 == k) amom2 = amom2 - 4.0/(ak22 - 1.0);
        res24 = res24 +  ARRAYF(cheb24, k)*amom2;
        amom0 = amom1;
        amom1 = amom2;
    }
    result = res24;
    abserr = std::abs(res24 - res12);
    label_50:
    return;
}

void dqc25f(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double omega,
            const int integr,
            const int nrmom,
            const int maxp1,
            const int ksave,
            double &result,
            double &abserr,
            int &neval,
            double &resabs,
            double &resasc,
            int &momcom,
            double *chebmo, // chebmo
            void *data)
{

}

void dqcheb(const double *x,
            double *fval,
            double *cheb12,
            double *cheb24)
{
    double alam, alam1, alam2, part1, part2, part3, v[12];
    int i, j;
    // first executable statement  dqcheb
    for (i = 1; i <= 12; ++i) {
        j = 26 - i;
        ARRAYF(v, i)    = ARRAYF(fval, i) - ARRAYF(fval, j);
        ARRAYF(fval, i) = ARRAYF(fval, i) + ARRAYF(fval, j);
    }
    alam1 = ARRAYF(v, 1) - ARRAYF(v, 9);
    alam2 = ARRAYF(x, 6) * (ARRAYF(v, 3) - ARRAYF(v, 7) - ARRAYF(v, 11));
    ARRAYF(cheb12,  4) = alam1 + alam2;
    ARRAYF(cheb12, 10) = alam1 - alam2;
    alam1 = ARRAYF(v, 2) - ARRAYF(v, 8) - ARRAYF(v, 10);
    alam2 = ARRAYF(v, 4) - ARRAYF(v, 6) - ARRAYF(v, 12);
    alam  = ARRAYF(x, 3)*alam1 + ARRAYF(x, 9)*alam2;
    ARRAYF(cheb24,  4) = ARRAYF(cheb12, 4) + alam;
    ARRAYF(cheb24, 22) = ARRAYF(cheb12, 4) - alam;
    alam = ARRAYF(x, 9)*alam1 - ARRAYF(x, 3)*alam2;
    ARRAYF(cheb24, 10) = ARRAYF(cheb12, 10) + alam;
    ARRAYF(cheb24, 16) = ARRAYF(cheb12, 10) - alam;
    part1 = ARRAYF(x, 4) * ARRAYF(v, 5);
    part2 = ARRAYF(x, 8) * ARRAYF(v, 9);
    part3 = ARRAYF(x, 6) * ARRAYF(v, 7);
    alam1 = ARRAYF(v, 1) + part1 + part2;
    alam2 = ARRAYF(x, 2)*ARRAYF(v, 3) + part3 + ARRAYF(x, 10)*ARRAYF(v, 11);
    ARRAYF(cheb12,  2) = alam1 + alam2;
    ARRAYF(cheb12, 12) = alam1 - alam2;
    alam = ARRAYF(x, 1)*ARRAYF(v, 2) + ARRAYF(x, 3)*ARRAYF(v,  4) + ARRAYF(x,  5)*ARRAYF(v,  6) 
         + ARRAYF(x, 7)*ARRAYF(v, 8) + ARRAYF(x, 9)*ARRAYF(v, 10) + ARRAYF(x, 11)*ARRAYF(v, 12);
    ARRAYF(cheb24,  2) = ARRAYF(cheb12, 2) + alam;
    ARRAYF(cheb24, 24) = ARRAYF(cheb12, 2) - alam;
    alam = ARRAYF(x, 11)*ARRAYF(v, 2) - ARRAYF(x, 9)*ARRAYF(v,  4) + ARRAYF(x, 7)*ARRAYF(v,  6)
         - ARRAYF(x,  5)*ARRAYF(v, 8) + ARRAYF(x, 3)*ARRAYF(v, 10) - ARRAYF(x, 1)*ARRAYF(v, 12);
    ARRAYF(cheb24, 12) = ARRAYF(cheb12, 12) + alam;
    ARRAYF(cheb24, 14) = ARRAYF(cheb12, 12) - alam;
    alam1 = ARRAYF(v, 1) - part1 + part2;
    alam2 = ARRAYF(x, 10)*ARRAYF(v, 3) - part3 + ARRAYF(x, 2)*ARRAYF(v, 11);
    ARRAYF(cheb12, 6) = alam1 + alam2;
    ARRAYF(cheb12, 8) = alam1 - alam2;
    alam = ARRAYF(x,  5)*ARRAYF(v, 2) - ARRAYF(x, 9)*ARRAYF(v,  4) - ARRAYF(x, 1)*ARRAYF(v,  6)
         - ARRAYF(x, 11)*ARRAYF(v, 8) + ARRAYF(x, 3)*ARRAYF(v, 10) + ARRAYF(x, 7)*ARRAYF(v, 12);
    ARRAYF(cheb24,  6) = ARRAYF(cheb12, 6) + alam;
    ARRAYF(cheb24, 20) = ARRAYF(cheb12, 6) - alam;
    alam = ARRAYF(x, 7)*ARRAYF(v, 2) - ARRAYF(x, 3)*ARRAYF(v,  4) - ARRAYF(x, 11)*ARRAYF(v,  6)
         + ARRAYF(x, 1)*ARRAYF(v, 8) - ARRAYF(x, 9)*ARRAYF(v, 10) - ARRAYF(x,  5)*ARRAYF(v, 12);
    ARRAYF(cheb24,  8) = ARRAYF(cheb12, 8) + alam;
    ARRAYF(cheb24, 18) = ARRAYF(cheb12, 8) - alam;
    for (i = 1; i <= 6; ++i) {
        j = 14 - i;
        ARRAYF(v, i)    = ARRAYF(fval, i) - ARRAYF(fval, j);
        ARRAYF(fval, i) = ARRAYF(fval, i) + ARRAYF(fval, j);
    }
    alam1 = ARRAYF(v, 1) + ARRAYF(x, 8)*ARRAYF(v, 5);
    alam2 = ARRAYF(x, 4)*ARRAYF(v, 3);
    ARRAYF(cheb12,  3) = alam1 + alam2;
    ARRAYF(cheb12, 11) = alam1 - alam2;
    ARRAYF(cheb12,  7) = ARRAYF(v, 1) - ARRAYF(v, 5);
    alam = ARRAYF(x, 2)*ARRAYF(v, 2) + ARRAYF(x, 6)*ARRAYF(v, 4) + ARRAYF(x, 10)*ARRAYF(v, 6);
    ARRAYF(cheb24,  3) = ARRAYF(cheb12, 3) + alam;
    ARRAYF(cheb24, 23) = ARRAYF(cheb12, 3) - alam;
    alam = ARRAYF(x, 6) * (ARRAYF(v, 2) - ARRAYF(v, 4) - ARRAYF(v, 6));
    ARRAYF(cheb24,  7) = ARRAYF(cheb12, 7) + alam;
    ARRAYF(cheb24, 19) = ARRAYF(cheb12, 7) - alam;
    alam = ARRAYF(x, 10)*ARRAYF(v, 2) - ARRAYF(x, 6)*ARRAYF(v, 4) + ARRAYF(x, 2)*ARRAYF(v, 6);
    ARRAYF(cheb24, 11) = ARRAYF(cheb12, 11) + alam;
    ARRAYF(cheb24, 15) = ARRAYF(cheb12, 11) - alam;
    for (i = 1; i <= 3; ++i) {
        j = 8 - i;
        ARRAYF(v, i)    = ARRAYF(fval, i) - ARRAYF(fval, j);
        ARRAYF(fval, i) = ARRAYF(fval, i) + ARRAYF(fval, j);
    }
    ARRAYF(cheb12, 5) = ARRAYF(v, 1) + ARRAYF(x, 8)*ARRAYF(v, 3);
    ARRAYF(cheb12, 9) = ARRAYF(fval, 1) - ARRAYF(x, 8)*ARRAYF(fval, 3);
    alam = ARRAYF(x, 4) * ARRAYF(v, 2);
    ARRAYF(cheb24,  5) = ARRAYF(cheb12, 5) + alam;
    ARRAYF(cheb24, 21) = ARRAYF(cheb12, 5) - alam;
    alam = ARRAYF(x, 8)*ARRAYF(fval, 2) - ARRAYF(fval, 4);
    ARRAYF(cheb24,  9) = ARRAYF(cheb12, 9) + alam;
    ARRAYF(cheb24, 17) = ARRAYF(cheb12, 9) - alam;
    ARRAYF(cheb12,  1) = ARRAYF(fval, 1) + ARRAYF(fval, 3);
    alam =ARRAYF(fval, 2) + ARRAYF(fval, 4);
    ARRAYF(cheb24,  1) = ARRAYF(cheb12, 1) + alam;
    ARRAYF(cheb24, 25) = ARRAYF(cheb12, 1) - alam;
    ARRAYF(cheb12, 13) = ARRAYF(v, 1) - ARRAYF(v, 3);
    ARRAYF(cheb24, 13) = ARRAYF(cheb12, 13);
    alam = 1.0/6.0;
    for (i = 2; i <= 12; ++i) {
        ARRAYF(cheb12, i) = ARRAYF(cheb12, i)*alam;
    }
    alam = 0.5 * alam;
    ARRAYF(cheb12,  1) = ARRAYF(cheb12,  1)*alam;
    ARRAYF(cheb12, 13) = ARRAYF(cheb12, 13)*alam;
    for (i = 2; i <= 24; ++i) {
        ARRAYF(cheb24, i) = ARRAYF(cheb24, i)*alam;
    }
    ARRAYF(cheb24,  1) = 0.5 * alam * ARRAYF(cheb24,  1);
    ARRAYF(cheb24, 25) = 0.5 * alam * ARRAYF(cheb24, 25);
    return;
}

void dqelg(int &n,
           double *epstab,
           double &result,
           double &abserr,
           double *res3la,
           int &nres)
{
    double delta1, delta2, delta3, epsinf,
           err1, err2, err3, e0, e1, e1abs,
           e2, e3, res, ss, tol1, tol2, tol3;
    double epmach, oflow;
    int i, ib, ib2, ie, indx, k1, k2, k3, num;
    int newelm;
    double error;

    // first executable statement  dqelg
    epmach = ARRAYF(d1mach, 4);
    oflow  = ARRAYF(d1mach, 2);
    nres  += 1;
    abserr = oflow;
    result = ARRAYF(epstab, n);
    if(n < 3) goto label_100;
    // limexp = 50
    ARRAYF(epstab, n+2) = ARRAYF(epstab, n);
    newelm = (n - 1)/2;
    ARRAYF(epstab, n) = oflow;
    num    = n;
    k1     = n;
    for (i = 1; i <= newelm; ++i) {
        k2     = k1 - 1;
        k3     = k1 - 2;
        res    = ARRAYF(epstab, k1+2);
        e0     = ARRAYF(epstab, k3);
        e1     = ARRAYF(epstab, k2);
        e2     = res;
        e1abs  = std::abs(e1);
        delta2 = e2 - e1;
        err2   = std::abs(delta2);
        tol2   = std::max(std::abs(e2), e1abs)*epmach;
        delta3 = e1 - e0;
        err3   = std::abs(delta3);
        tol3   = std::max(e1abs, std::abs(e0))*epmach;
        if (err2 > tol2 || err3 > tol3) goto label_10;
        //
        // if e0, e1 and e2 are equal to within machine
        // accuracy, convergence is assumed.
        // result = e2
        // abserr = abs(e1-e0)+abs(e2-e1)
        //
        result = res;
        abserr = err2 + err3;
        // jump out of for-loop
        goto label_100;
        label_10:
        e3           = ARRAYF(epstab, k1);
        ARRAYF(epstab, k1) = e1;
        delta1       = e1 - e3;
        err1         = std::abs(delta1);
        tol1         = std::max(e1abs, std::abs(e3))*epmach;
        //
        // if two elements are very close to each other, omit
        // a part of the table by adjusting the value of n
        //
        if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) goto label_20;
        ss = 1.0/delta1 +1.0/delta2 - 1.0/delta3;
        epsinf = std::abs(ss*e1);
        //
        // test to detect irregular behaviour in the table, and
        // eventually omit a part of the table adjusting the value
        // of n.
        // 
        if (epsinf > 1.0e-4) goto label_30;
        label_20:
        n = i + i - 1;
        // jump out of for-loop
        goto label_50;
        //
        // compute a new element and eventually adjust
        // the value of result.
        //
        label_30:
        res          = e1 + 1.0/ss;
        ARRAYF(epstab, k1) = res;
        k1           = k1 - 2;
        error        = err2 + std::abs(res - e2) + err3;
        if (error > abserr) continue; // goto label_40;
        abserr = error;
        result = res;
        // label_40:
        // continue;
    }
    //
    // shift the table.
    //
    label_50:
    if (n == limexp) n = 2*(limexp/2) - 1;
    ib = 1;
    if ((num/2)*2 == num) ib = 2;
    ie = newelm + 1;
    for (i = 1; i <= ie; ++i) {
        ib2 = ib + 2;
        ARRAYF(epstab, ib) = ARRAYF(epstab, ib2);
        ib = ib2;
    }
    if(num == n) goto label_80;
    indx = num - n + 1;
    for (i = 1; i <= n; ++i) {
        ARRAYF(epstab, i) = ARRAYF(epstab, indx);
        indx += 1;
    }
    label_80:
    if(nres >= 4) goto label_90;
    ARRAYF(res3la, nres) = result;
    abserr = oflow;
    goto label_100;
    //
    // compute error estimate
    // 
    label_90:
    abserr = std::abs(result - ARRAYF(res3la, 3)) + std::abs(result - ARRAYF(res3la, 2)) + std::abs(result - ARRAYF(res3la, 1));
    ARRAYF(res3la, 1) = ARRAYF(res3la, 2);
    ARRAYF(res3la, 2) = ARRAYF(res3la, 3);
    ARRAYF(res3la, 3) = result;
    label_100:
    abserr = std::max(abserr, 5.0*epmach*std::abs(result));
    return;
}


void dqk15(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data)
{
    double dhlgth, fc, fsum, fv1[7], fv2[7];
    double uflow, epmach;
    int j, jtw, jtwm1;
    double centr;
    double hlgth;
    double absc;
    double fval1;
    double fval2;
    double resg;
    double resk;
    double reskh;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 15-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 7-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 7-point gauss rule
    // 
    // wgk    - weights of the 15-point kronrod rule
    //
    // wg     - weights of the 7-point gauss rule
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    // 
    const double wg[4] = {
        1.29484966168869693270611432679082018329e-1,
        2.79705391489276667901467771423779582487e-1,
        3.81830050505118944950369775488975133878e-1,
        4.17959183673469387755102040816326530612e-1
    };
    const double xgk[8] = {
        9.91455371120812639206854697526328516642e-1,
        9.49107912342758524526189684047851262401e-1,
        8.64864423359769072789712788640926201211e-1,
        7.41531185599394439863864773280788407074e-1,
        5.86087235467691130294144838258729598437e-1,
        4.05845151377397166906606412076961463347e-1,
        2.07784955007898467600689403773244913480e-1,
        0.00000000000000000000000000000000000000
    };
    const double wgk[8] = {
        2.29353220105292249637320080589695919936e-2,
        6.30920926299785532907006631892042866651e-2,
        1.04790010322250183839876322541518017444e-1,
        1.40653259715525918745189590510237920400e-1,
        1.69004726639267902826583426598550284106e-1,
        1.90350578064785409913256402421013682826e-1,
        2.04432940075298892414161999234649084717e-1,
        2.09482141084727828012999174891714263698e-1
    };
    //
    // first executable statement  dqk15
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 15-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    fc     = f(centr, data);
    resg   = fc * ARRAYF(wg, 4);
    resk   = fc * ARRAYF(wgk, 8);
    resabs = std::abs(resk);
    for (j = 1; j <= 3; ++j) {
        jtw = j * 2;
        absc = hlgth * ARRAYF(xgk, jtw);
        fval1 = f(centr - absc, data);
        fval2 = f(centr + absc, data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum = fval1 + fval2;
        resg = resg + ARRAYF(wg, j) * fsum;
        resk = resk + ARRAYF(wgk, jtw) * fsum;
        resabs = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 4; ++j) {
        jtwm1      = 2*j - 1;
        absc       = hlgth * ARRAYF(xgk, jtwm1);
        fval1      = f(centr - absc, data);
        fval2      = f(centr + absc, data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum       = fval1 + fval2;
        resk       = resk + ARRAYF(wgk, jtwm1) * fsum;
        resabs     = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 8) * std::abs(fc - reskh);
    for (j = 1; j <= 7; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2*abserr/resasc, 1.5));
    if (resabs > uflow/(50.0*epmach)) abserr = std::max((epmach*50.0)*resabs, abserr);
    return;
}

void dqk15i(QUADPACK_CPP_FUNCTION f,
            const double boun,
            const int inf,
            const double a,
            const double b,
            double &result,
            double &abserr,
            double &resabs,
            double &resasc,
            void *data)
{
    double absc, dinf, fc, fsum, fv1[7], fv2[7];
    double uflow, epmach;
    int j;
    double centr;
    double hlgth;
    double absc1;
    double absc2;
    double tabsc1;
    double tabsc2;
    double fval1;
    double fval2;
    double resg;
    double resk;
    double reskh;
    //
    // the abscissae and weights are supplied for the interval
    // (-1,1).  because of symmetry only the positive abscissae and
    // their corresponding weights are given.
    //
    // xgk    - abscissae of the 15-point kronrod rule
    //          xgk(2), xgk(4), ... abscissae of the 7-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 7-point gauss rule
    //
    // wgk    - weights of the 15-point kronrod rule
    //
    // wg     - weights of the 7-point gauss rule, corresponding
    //          to the abscissae xgk(2), xgk(4), ...
    //          wg(1), wg(3), ... are set to zero.
    //
    const double wg[8] = {
        0.00000000000000000000000000000000000000,
        1.29484966168869693270611432679082018329e-1,
        0.00000000000000000000000000000000000000,
        2.79705391489276667901467771423779582487e-1,
        0.00000000000000000000000000000000000000,
        3.81830050505118944950369775488975133878e-1,
        0.00000000000000000000000000000000000000,
        4.17959183673469387755102040816326530612e-1
    };
    const double xgk[8] = {
        9.91455371120812639206854697526328516642e-1,
        9.49107912342758524526189684047851262401e-1,
        8.64864423359769072789712788640926201211e-1,
        7.41531185599394439863864773280788407074e-1,
        5.86087235467691130294144838258729598437e-1,
        4.05845151377397166906606412076961463347e-1,
        2.07784955007898467600689403773244913480e-1,
        0.00000000000000000000000000000000000000
    };
    const double wgk[8] = {
        2.29353220105292249637320080589695919936e-2,
        6.30920926299785532907006631892042866651e-2,
        1.04790010322250183839876322541518017444e-1,
        1.40653259715525918745189590510237920400e-1,
        1.69004726639267902826583426598550284106e-1,
        1.90350578064785409913256402421013682826e-1,
        2.04432940075298892414161999234649084717e-1,
        2.09482141084727828012999174891714263698e-1
    };
    // first executable statement  dqk15i
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    dinf   = std::max(1.0, double(inf));
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    tabsc1 = boun + dinf * (1.0 - centr) / centr;
    fval1  = f(tabsc1, data);
    if(inf == 2) fval1 = fval1 + f(-tabsc1, data);
    fc     = (fval1 / centr) / centr;
    //
    // compute the 15-point kronrod approximation to
    // the integral, and estimate the error.
    //
    resg   = ARRAYF(wg, 8)  * fc;
    resk   = ARRAYF(wgk, 8) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 7; ++j) {
        absc   = hlgth * ARRAYF(xgk, j);
        absc1  = centr - absc;
        absc2  = centr + absc;
        tabsc1 = boun + dinf * (1.0 - absc1) / absc1;
        tabsc2 = boun + dinf * (1.0 - absc2) / absc2;
        fval1  = f(tabsc1, data);
        fval2  = f(tabsc2, data);
        if(inf == 2) fval1 = fval1 + f(-tabsc1, data);
        if(inf == 2) fval2 = fval2 + f(-tabsc2, data);
        fval1  = (fval1 / absc1) / absc1;
        fval2  = (fval2 / absc2) / absc2;
        ARRAYF(fv1, j) = fval1;
        ARRAYF(fv2, j) = fval2;
        fsum   = fval1 + fval2;
        resg   = resg + ARRAYF(wg, j)  * fsum;
        resk   = resk + ARRAYF(wgk, j) * fsum;
        resabs = resabs + ARRAYF(wgk, j) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 8) * std::abs(fc - reskh);
    for (j = 1; j <= 7; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * hlgth;
    resasc = resasc * hlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2*abserr/resasc, 1.5));
    if (resabs > uflow/(50.0*epmach)) abserr = std::max((epmach*50.0)*resabs, abserr);
    return;
}

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
            void *data)
{
    double absc1, absc2, dhlgth, fc, fsum, fv1[7], fv2[7];
    int j, jtw, jtwm1;
    double centr;
    double hlgth;
    double absc;
    double fval1;
    double fval2;
    double resg;
    double resk;
    double reskh;
    double epmach, uflow;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.

    // xgk    - abscissae of the 15-point gauss-kronrod rule
    //          xgk(2), xgk(4), ... abscissae of the 7-point
    //          gauss rule
    //          xgk(1), xgk(3), ... abscissae which are optimally
    //          added to the 7-point gauss rule

    // wgk    - weights of the 15-point gauss-kronrod rule

    // wg     - weights of the 7-point gauss rule
    //
    const double xgk[8] = {
        9.91455371120812639206854697526328516642e-1,
        9.49107912342758524526189684047851262401e-1,
        8.64864423359769072789712788640926201211e-1,
        7.41531185599394439863864773280788407074e-1,
        5.86087235467691130294144838258729598437e-1,
        4.05845151377397166906606412076961463347e-1,
        2.07784955007898467600689403773244913480e-1,
        0.00000000000000000000000000000000000000
    };
    const double wgk[8] = {
        2.29353220105292249637320080589695919936e-2,
        6.30920926299785532907006631892042866651e-2,
        1.04790010322250183839876322541518017444e-1,
        1.40653259715525918745189590510237920400e-1,
        1.69004726639267902826583426598550284106e-1,
        1.90350578064785409913256402421013682826e-1,
        2.04432940075298892414161999234649084717e-1,
        2.09482141084727828012999174891714263698e-1
    };
    const double wg[4] = {
        1.29484966168869693270611432679082018329e-1,
        2.79705391489276667901467771423779582487e-1,
        3.81830050505118944950369775488975133878e-1,
        4.17959183673469387755102040816326530612e-1
    };
    // first executable statement  dqk15w
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 15-point kronrod approximation to the
    // integral, and estimate the error.
    //
    fc = f(centr, data) * w(centr, p1, p2, p3, p4, kp);
    resg = ARRAYF(wg,  4) * fc;
    resk = ARRAYF(wgk, 8) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 3; ++j) {
        jtw = j*2;
        absc = hlgth * ARRAYF(xgk, jtw); 
        absc1 = centr - absc;
        absc2 = centr + absc;
        fval1 = f(absc1, data) * w(absc1, p1, p2, p3, p4, kp);
        fval2 = f(absc2, data) * w(absc2, p1, p2, p3, p4, kp);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum = fval1 + fval2;
        resg = resg + ARRAYF(wg,    j)*fsum;
        resk = resk + ARRAYF(wgk, jtw)*fsum;
        resabs = resabs + ARRAYF(wgk, jtw)*(std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 4; ++j) {
        jtwm1 = j*2 - 1;
        absc  = hlgth * ARRAYF(xgk, jtwm1);
        absc1 = centr - absc;
        absc2 = centr + absc;
        fval1 = f(absc1, data) * w(absc1, p1, p2, p3, p4, kp);
        fval2 = f(absc2, data) * w(absc2, p1, p2, p3, p4, kp);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum = fval1 + fval2;
        resk = resk + ARRAYF(wgk, jtwm1)*fsum;
        resabs = resabs + ARRAYF(wgk, jtwm1)*(std::abs(fval1) + std::abs(fval2));
    }
    reskh = resk*0.5;
    resasc = ARRAYF(wgk, 8) * std::abs(fc - reskh);
    for (j = 1; j <= 7; ++j) {
        resasc = resasc + ARRAYF(wgk, j)*(std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk*hlgth;
    resabs = resabs*dhlgth;
    resasc = resasc*dhlgth;
    abserr = std::abs((resk - resg)*hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc*std::min(1.0, std::pow(2.0e2*abserr/resasc, 1.5));
    if (resabs > uflow/(50.0*epmach)) abserr = std::max((epmach*50.0)*resabs, abserr);
    return;
}

void dqk21(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data)
{
    double dhlgth, fc, fsum, fv1[10], fv2[10];
    double uflow, epmach;
    int j, jtw, jtwm1;
    double centr;
    double hlgth;
    double absc;
    double fval1;
    double fval2;
    double resg;
    double resk;
    double reskh;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 21-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 10-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 10-point gauss rule
    //
    // wgk    - weights of the 21-point kronrod rule
    //
    // wg     - weights of the 10-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // c bell labs, nov. 1981.
    //
    const double wg[5] = {
        6.66713443086881375935688098933317928579e-2,
        1.49451349150580593145776339657697332403e-1,
        2.19086362515982043995534934228163192459e-1,
        2.69266719309996355091226921569469352860e-1,
        2.95524224714752870173892994651338329421e-1
    };
    const double xgk[11] = {
        9.95657163025808080735527280689002847921e-1,
        9.73906528517171720077964012084452053428e-1,
        9.30157491355708226001207180059508346225e-1,
        8.65063366688984510732096688423493048528e-1,
        7.80817726586416897063717578345042377163e-1,
        6.79409568299024406234327365114873575769e-1,
        5.62757134668604683339000099272694140843e-1,
        4.33395394129247190799265943165784162200e-1,
        2.94392862701460198131126603103865566163e-1,
        1.48874338981631210884826001129719984618e-1,
        0.00000000000000000000000000000000000000
    };
    const double wgk[11] = {
        1.16946388673718742780643960621920483962e-2,
        3.25581623079647274788189724593897606174e-2,
        5.47558965743519960313813002445801763737e-2,
        7.50396748109199527670431409161900093952e-2,
        9.31254545836976055350654650833663443900e-2,
        1.09387158802297641899210590325804960272e-1,
        1.23491976262065851077958109831074159512e-1,
        1.34709217311473325928054001771706832761e-1,
        1.42775938577060080797094273138717060886e-1,
        1.47739104901338491374841515972068045524e-1,
        1.49445554002916905664936468389821203745e-1
    };
    // first executable statement  dqk21
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 21-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    resg   = 0.0;
    fc     = f(centr, data);
    resk   = ARRAYF(wgk, 11) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 5; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, data);
        fval2            = f(centr + absc, data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg, j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 5; ++j) {
        jtwm1      = 2*j - 1;
        absc       = hlgth * ARRAYF(xgk, jtwm1);
        fval1      = f(centr - absc, data);
        fval2      = f(centr + absc, data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum       = fval1 + fval2;
        resk       = resk + ARRAYF(wgk, jtwm1) * fsum;
        resabs     = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 11) * std::abs(fc - reskh);
    for (j = 1; j <= 10; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2*abserr/resasc, 1.5));
    if (resabs > uflow/(50.0*epmach)) abserr = std::max((epmach*50.0)*resabs, abserr);
    return;
}


void dqk31(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data)
{
    double dhlgth, fc, fsum, fv1[15], fv2[15];
    double uflow, epmach;
    int j, jtw, jtwm1;
    double centr;
    double hlgth;
    double absc;
    double fval1;
    double fval2;
    double resg;
    double resk;
    double reskh;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 31-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 15-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 15-point gauss rule
    //
    // wgk    - weights of the 31-point kronrod rule
    //
    // wg     - weights of the 15-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    //
    const double wg[8] = {
        3.07532419961172683546283935772044177217e-2,
        7.03660474881081247092674164506673384667e-2,
        1.07159220467171935011869546685869303416e-1,
        1.39570677926154314447804794511028322521e-1,
        1.66269205816993933553200860481208811131e-1,
        1.86161000015562211026800561866422824506e-1,
        1.98431485327111576456118326443839324819e-1,
        2.02578241925561272880620199967519314839e-1
    };
    const double xgk[16] = {
        9.98002298693397060285172840152271209073e-1,
        9.87992518020485428489565718586612581147e-1,
        9.67739075679139134257347978784337225283e-1,
        9.37273392400705904307758947710209471244e-1,
        8.97264532344081900882509656454495882832e-1,
        8.48206583410427216200648320774216851366e-1,
        7.90418501442465932967649294817947346862e-1,
        7.24417731360170047416186054613938009631e-1,
        6.50996741297416970533735895313274692547e-1,
        5.70972172608538847537226737253910641238e-1,
        4.85081863640239680693655740232350612866e-1,
        3.94151347077563369897207370981045468363e-1,
        2.99180007153168812166780024266388962662e-1,
        2.01194093997434522300628303394596207813e-1,
        1.01142066918717499027074231447392338787e-1,
        0.00000000000000000000000000000000000000
    };
    const double wgk[16] = {
        5.37747987292334898779205143012764981831e-3,
        1.50079473293161225383747630758072680946e-2,
        2.54608473267153201868740010196533593973e-2,
        3.53463607913758462220379484783600481226e-2,
        4.45897513247648766082272993732796902233e-2,
        5.34815246909280872653431472394302967716e-2,
        6.20095678006706402851392309608029321904e-2,
        6.98541213187282587095200770991474757860e-2,
        7.68496807577203788944327774826590067221e-2,
        8.30805028231330210382892472861037896016e-2,
        8.85644430562117706472754436937743032123e-2,
        9.31265981708253212254868727473457185619e-2,
        9.66427269836236785051799076275893351367e-2,
        9.91735987217919593323931734846031310596e-2,
        1.00769845523875595044946662617569721916e-1,
        1.01330007014791549017374792767492546771e-1
    };
    // first executable statement  dqk31
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 31-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    fc     = f(centr, data);
    resg   = ARRAYF(wg, 8) * fc;
    resk   = ARRAYF(wgk, 16) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 7; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, data);
        fval2            = f(centr + absc, data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg, j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 8; ++j) {
        jtwm1      = 2*j - 1;
        absc       = hlgth * ARRAYF(xgk, jtwm1);
        fval1      = f(centr - absc, data);
        fval2      = f(centr + absc, data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum       = fval1 + fval2;
        resk       = resk + ARRAYF(wgk, jtwm1) * fsum;
        resabs     = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 16) * std::abs(fc - reskh);
    for (j = 1; j <= 15; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2*abserr/resasc, 1.5));
    if (resabs > uflow/(50.0*epmach)) abserr = std::max((epmach*50.0)*resabs, abserr);
    return;
}

void dqk41(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data)
{
    double dhlgth, fc, fsum, fv1[20], fv2[20];
    double uflow, epmach;
    int j, jtw, jtwm1;
    double centr;
    double hlgth;
    double absc;
    double fval1;
    double fval2;
    double resg;
    double resk;
    double reskh;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 41-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 20-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 20-point gauss rule
    //
    // wgk    - weights of the 41-point kronrod rule
    //
    // wg     - weights of the 20-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    //
    const double wg[10] = {
        1.76140071391521183118619623518528163621e-2,
        4.06014298003869413310399522749321098791e-2,
        6.26720483341090635695065351870416063516e-2,
        8.32767415767047487247581432220462061002e-2,
        1.01930119817240435036750135480349876167e-1,
        1.18194531961518417312377377711382287005e-1,
        1.31688638449176626898494499748163134916e-1,
        1.42096109318382051329298325067164933035e-1,
        1.49172986472603746787828737001969436693e-1,
        1.52753387130725850698084331955097593492e-1
    };
    const double xgk[21] = {
        9.98859031588277663838315576545863010000e-1,
        9.93128599185094924786122388471320278223e-1,
        9.81507877450250259193342994720216944567e-1,
        9.63971927277913791267666131197277221912e-1,
        9.40822633831754753519982722212443380274e-1,
        9.12234428251325905867752441203298113049e-1,
        8.78276811252281976077442995113078466711e-1,
        8.39116971822218823394529061701520685330e-1,
        7.95041428837551198350638833272787942959e-1,
        7.46331906460150792614305070355641590311e-1,
        6.93237656334751384805490711845931533386e-1,
        6.36053680726515025452836696226285936743e-1,
        5.75140446819710315342946036586425132814e-1,
        5.10867001950827098004364050955250998425e-1,
        4.43593175238725103199992213492640107840e-1,
        3.73706088715419560672548177024927237396e-1,
        3.01627868114913004320555356858592260615e-1,
        2.27785851141645078080496195368574624743e-1,
        1.52605465240922675505220241022677527912e-1,
        7.65265211334973337546404093988382110048e-2,
        0.00000000000000000000000000000000000000e0
    };
    const double wgk[21] = {
        3.07358371852053150121829324603098748803e-3,
        8.60026985564294219866178795010234725213e-3,
        1.46261692569712529837879603088683561639e-2,
        2.03883734612665235980102314327547051228e-2,
        2.58821336049511588345050670961531429995e-2,
        3.12873067770327989585431193238007378878e-2,
        3.66001697582007980305572407072110084875e-2,
        4.16688733279736862637883059368947380440e-2,
        4.64348218674976747202318809261075168421e-2,
        5.09445739237286919327076700503449486648e-2,
        5.51951053482859947448323724197773291948e-2,
        5.91114008806395723749672206485942171364e-2,
        6.26532375547811680258701221742549805858e-2,
        6.58345971336184221115635569693979431472e-2,
        6.86486729285216193456234118853678017155e-2,
        7.10544235534440683057903617232101674129e-2,
        7.30306903327866674951894176589131127606e-2,
        7.45828754004991889865814183624875286161e-2,
        7.57044976845566746595427753766165582634e-2,
        7.63778676720807367055028350380610018008e-2,
        7.66007119179996564450499015301017408279e-2
    };
    // first executable statement  dqk31
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 41-point gauss-kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    resg   = 0.0;
    fc     = f(centr, data);
    resk   = ARRAYF(wgk, 21) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 10; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, data);
        fval2            = f(centr + absc, data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg, j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 10; ++j) {
        jtwm1              = 2*j - 1;
        absc               = hlgth * ARRAYF(xgk, jtwm1);
        fval1              = f(centr - absc, data);
        fval2              = f(centr + absc, data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum               = fval1 + fval2;
        resk               = resk + ARRAYF(wgk, jtwm1) * fsum;
        resabs             = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 21) * std::abs(fc - reskh);
    for (j = 1; j <= 20; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2*abserr/resasc, 1.5));
    if (resabs > uflow/(50.0*epmach)) abserr = std::max((epmach*50.0)*resabs, abserr);
    return;
}

void dqk51(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data)
{
    double dhlgth, fc, fsum, fv1[25], fv2[25];
    double uflow, epmach;
    int j, jtw, jtwm1;
    double centr;
    double hlgth;
    double absc;
    double fval1;
    double fval2;
    double resg;
    double resk;
    double reskh;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 51-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 25-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 25-point gauss rule
    //
    // wgk    - weights of the 51-point kronrod rule
    //
    // wg     - weights of the 25-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    //
    const double wg[13] = {
        1.13937985010262879479029641132347736033e-2,
        2.63549866150321372619018152952991449360e-2,
        4.09391567013063126556234877116459536608e-2,
        5.49046959758351919259368915404733241601e-2,
        6.80383338123569172071871856567079685547e-2,
        8.01407003350010180132349596691113022902e-2,
        9.10282619829636498114972207028916533810e-2,
        1.00535949067050644202206890392685826988e-1,
        1.08519624474263653116093957050116619340e-1,
        1.14858259145711648339325545869555808641e-1,
        1.19455763535784772228178126512901047390e-1,
        1.22242442990310041688959518945851505835e-1,
        1.23176053726715451203902873079050142438e-1
    };
    const double xgk[26] = {
        9.99262104992609834193457486540340593705e-1,
        9.95556969790498097908784946893901617258e-1,
        9.88035794534077247637331014577406227072e-1,
        9.76663921459517511498315386479594067745e-1,
        9.61614986425842512418130033660167241692e-1,
        9.42974571228974339414011169658470531905e-1,
        9.20747115281701561746346084546330631575e-1,
        8.94991997878275368851042006782804954175e-1,
        8.65847065293275595448996969588340088203e-1,
        8.33442628760834001421021108693569569461e-1,
        7.97873797998500059410410904994306569409e-1,
        7.59259263037357630577282865204360976388e-1,
        7.17766406813084388186654079773297780598e-1,
        6.73566368473468364485120633247622175883e-1,
        6.26810099010317412788122681624517881020e-1,
        5.77662930241222967723689841612654067396e-1,
        5.26325284334719182599623778158010178037e-1,
        4.73002731445714960522182115009192041332e-1,
        4.17885382193037748851814394594572487093e-1,
        3.61172305809387837735821730127640667422e-1,
        3.03089538931107830167478909980339329200e-1,
        2.43866883720988432045190362797451586406e-1,
        1.83718939421048892015969888759528415785e-1,
        1.22864692610710396387359818808036805532e-1,
        6.15444830056850788865463923667966312817e-2,
        0.00000000000000000000000000000000000000e0 
    };
    const double wgk[26] = {
        1.98738389233031592650785188284340988943e-3,
        5.56193213535671375804023690106552207018e-3,
        9.47397338617415160720771052365532387165e-3,
        1.32362291955716748136564058469762380776e-2,
        1.68478177091282982315166675363363158404e-2,
        2.04353711458828354565682922359389736788e-2,
        2.40099456069532162200924891648810813929e-2,
        2.74753175878517378029484555178110786148e-2,
        3.07923001673874888911090202152285856009e-2,
        3.40021302743293378367487952295512032257e-2,
        3.71162714834155435603306253676198759960e-2,
        4.00838255040323820748392844670756464014e-2,
        4.28728450201700494768957924394951611020e-2,
        4.55029130499217889098705847526603930437e-2,
        4.79825371388367139063922557569147549836e-2,
        5.02776790807156719633252594334400844406e-2,
        5.23628858064074758643667121378727148874e-2,
        5.42511298885454901445433704598756068261e-2,
        5.59508112204123173082406863827473468203e-2,
        5.74371163615678328535826939395064719948e-2,
        5.86896800223942079619741758567877641398e-2,
        5.97203403241740599790992919325618538354e-2,
        6.05394553760458629453602675175654271623e-2,
        6.11285097170530483058590304162927119227e-2,
        6.14711898714253166615441319652641775865e-2,
        6.15808180678329350787598242400645531904e-2
    };
    // first executable statement  dqk31
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 51-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    fc     = f(centr, data);
    resg   = ARRAYF(wg, 13) * fc;
    resk   = ARRAYF(wgk, 26) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 12; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, data);
        fval2            = f(centr + absc, data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg, j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 13; ++j) {
        jtwm1              = 2*j - 1;
        absc               = hlgth * ARRAYF(xgk, jtwm1);
        fval1              = f(centr - absc, data);
        fval2              = f(centr + absc, data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum               = fval1 + fval2;
        resk               = resk + ARRAYF(wgk, jtwm1) * fsum;
        resabs             = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 26) * std::abs(fc - reskh);
    for (j = 1; j <= 25; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2*abserr/resasc, 1.5));
    if (resabs > uflow/(50.0*epmach)) abserr = std::max((epmach*50.0)*resabs, abserr);
    return;
}

void dqk61(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data)
{
    double dhlgth, fc, fsum, fv1[30], fv2[30];
    double uflow, epmach;
    int j, jtw, jtwm1;
    double centr;
    double hlgth;
    double absc;
    double fval1;
    double fval2;
    double resg;
    double resk;
    double reskh;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 61-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 30-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 30-point gauss rule
    //
    // wgk    - weights of the 61-point kronrod rule
    //
    // wg     - weights of the 30-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    //
    const double wg[15] = {
        7.96819249616660561546588347467362245048e-3,
        1.84664683110909591423021319120472690962e-2,
        2.87847078833233693497191796112920436396e-2,
        3.87991925696270495968019364463476920332e-2,
        4.84026728305940529029381404228075178153e-2,
        5.74931562176190664817216894020561287971e-2,
        6.59742298821804951281285151159623612374e-2,
        7.37559747377052062682438500221907341538e-2,
        8.07558952294202153546949384605297308759e-2,
        8.68997872010829798023875307151257025768e-2,
        9.21225222377861287176327070876187671969e-2,
        9.63687371746442596394686263518098650964e-2,
        9.95934205867952670627802821035694765299e-2,
        1.01762389748405504596428952168554044633e-1,
        1.02852652893558840341285636705415043868e-1
        };
    const double xgk[31] = {
        9.99484410050490637571325895705810819469e-1,
        9.96893484074649540271630050918695283341e-1,
        9.91630996870404594858628366109485724851e-1,
        9.83668123279747209970032581605662801940e-1,
        9.73116322501126268374693868423706884888e-1,
        9.60021864968307512216871025581797662930e-1,
        9.44374444748559979415831324037439121586e-1,
        9.26200047429274325879324277080474004086e-1,
        9.05573307699907798546522558925958319569e-1,
        8.82560535792052681543116462530225590057e-1,
        8.57205233546061098958658510658943856821e-1,
        8.29565762382768397442898119732501916439e-1,
        7.99727835821839083013668942322683240736e-1,
        7.67777432104826194917977340974503131695e-1,
        7.33790062453226804726171131369527645669e-1,
        6.97850494793315796932292388026640068382e-1,
        6.60061064126626961370053668149270753038e-1,
        6.20526182989242861140477556431189299207e-1,
        5.79345235826361691756024932172540495907e-1,
        5.36624148142019899264169793311072794164e-1,
        4.92480467861778574993693061207708795644e-1,
        4.47033769538089176780609900322854000162e-1,
        4.00401254830394392535476211542660633611e-1,
        3.52704725530878113471037207089373860654e-1,
        3.04073202273625077372677107199256553531e-1,
        2.54636926167889846439805129817805107883e-1,
        2.04525116682309891438957671002024709524e-1,
        1.53869913608583546963794672743255920419e-1,
        1.02806937966737030147096751318000592472e-1,
        5.14718425553176958330252131667225737491e-2,
        0.00000000000000000000000000000000000000
        };
    const double wgk[31] = {
        1.38901369867700762455159122675969968105e-3,
        3.89046112709988405126720184451550327852e-3,
        6.63070391593129217331982636975016813363e-3,
        9.27327965951776342844114689202436042127e-3,
        1.18230152534963417422328988532505928963e-2,
        1.43697295070458048124514324435800101958e-2,
        1.69208891890532726275722894203220923686e-2,
        1.94141411939423811734089510501284558514e-2,
        2.18280358216091922971674857383389934015e-2,
        2.41911620780806013656863707252320267604e-2,
        2.65099548823331016106017093350754143665e-2,
        2.87540487650412928439787853543342111447e-2,
        3.09072575623877624728842529430922726353e-2,
        3.29814470574837260318141910168539275106e-2,
        3.49793380280600241374996707314678750972e-2,
        3.68823646518212292239110656171359677370e-2,
        3.86789456247275929503486515322810502509e-2,
        4.03745389515359591119952797524681142161e-2,
        4.19698102151642461471475412859697577901e-2,
        4.34525397013560693168317281170732580746e-2,
        4.48148001331626631923555516167232437574e-2,
        4.60592382710069881162717355593735805947e-2,
        4.71855465692991539452614781810994864829e-2,
        4.81858617570871291407794922983045926058e-2,
        4.90554345550297788875281653672381736059e-2,
        4.97956834270742063578115693799423285392e-2,
        5.04059214027823468408930856535850289022e-2,
        5.08817958987496064922974730498046918534e-2,
        5.12215478492587721706562826049442082511e-2,
        5.14261285374590259338628792157812598296e-2,
        5.14947294294515675583404336470993075327e-2
        };
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 51-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    resg   = 0.0;
    fc     = f(centr, data);
    resk   = ARRAYF(wgk, 31) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 15; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, data);
        fval2            = f(centr + absc, data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg, j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 15; ++j) {
        jtwm1              = 2*j - 1;
        absc               = hlgth * ARRAYF(xgk, jtwm1);
        fval1              = f(centr - absc, data);
        fval2              = f(centr + absc, data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum               = fval1 + fval2;
        resk               = resk + ARRAYF(wgk, jtwm1) * fsum;
        resabs             = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 31) * std::abs(fc - reskh);
    for (j = 1; j <= 30; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2*abserr/resasc, 1.5));
    if (resabs > uflow/(50.0*epmach)) abserr = std::max((epmach*50.0)*resabs, abserr);
    return;
}




void dqpsrt(const int limit,
            const int last,
            int &maxerr,
            double &ermax,
            const double *elist,
            int *iord,
            int &nrmax)
{
    double errmax, errmin;
    int i, ibeg, ido, isucc, j, jbnd, jupbn, k;
    //
    // check whether the list contains more than
    // two error estimates.
    //
    // first executable statement  dqpsrt
    if (last > 2) goto label_10;
    ARRAYF(iord, 1) = 1;
    ARRAYF(iord, 2) = 2;
    goto label_90;
    //
    // this part of the routine is only executed if, due to a
    // difficult integrand, subdivision increased the error
    // estimate. in the normal case the insert procedure should
    // start after the nrmax-th largest error estimate.
    //
    label_10:
    errmax = ARRAYF(elist, maxerr);
    if (nrmax == 1) goto label_30;
    ido = nrmax - 1;
    for (i = 1; i <= ido; ++i) {
        isucc = ARRAYF(iord, nrmax-1);
        // jump out of do-loop
        if (errmax <= ARRAYF(elist, isucc)) goto label_30;
        ARRAYF(iord, nrmax) = isucc;
        nrmax -= 1;
    }
    //
    // compute the number of elements in the list to be maintained
    // in descending order. this number depends on the number of
    // subdivisions still allowed.
    //
    label_30:
    jupbn = last;
    if (last > (limit/2 + 2)) jupbn = limit + 3 - last;
    errmin = ARRAYF(elist, last);
    //
    // insert errmax by traversing the list top-down,
    // starting comparison from the element elist(iord(nrmax+1)).
    // 
    jbnd = jupbn - 1;
    ibeg = nrmax + 1;
    if (ibeg > jbnd) goto label_50;
    for (i = ibeg; i <= jbnd; ++i) {
        isucc = ARRAYF(iord, i);
        // jump out of for-loop
        if(errmax >= ARRAYF(elist, isucc)) goto label_60;
        ARRAYF(iord, i-1) = isucc;
    }
    label_50:
    ARRAYF(iord, jbnd)  = maxerr;
    ARRAYF(iord, jupbn) = last;
    goto label_90;
    //
    // insert errmin by traversing the list bottom-up.
    //
    label_60:
    ARRAYF(iord, i-1) = maxerr;
    k = jbnd;
    for (j = i; j <= jbnd; ++j) {
        isucc = ARRAYF(iord, k);
        // jump out of for-loop
        if (errmin < ARRAYF(elist, isucc)) goto label_80;
        ARRAYF(iord, k+1) = isucc;
        k -= 1;
    }
    ARRAYF(iord, i) = last;
    goto label_90;
    label_80:
    ARRAYF(iord, k+1) = last;
    //
    // set maxerr and ermax.
    //
    label_90:
    maxerr = ARRAYF(iord, nrmax);
    ermax  = ARRAYF(elist, maxerr);
    return;
}

double dqwgtc(const double x,
              const double c,
              const double p2,
              const double p3,
              const double p4,
              const int kp)
{
    return 1.0 / (x - c);
}

double dqwgtf(const double x,
              const double omega,
              const double p2,
              const double p3,
              const double p4,
              const int integr)
{
    double omx;
    // first executable statement  dqwgtf
    omx = omega * x;
    if (integr == 1) {
        return std::cos(omx);
    } else if (integr == 2) {
        return std::sin(omx);
    } else {
        return 0.0;
    }
}

void xerror(const std::string messg, 
            const int nerr, 
            const int level)
{
    std::cerr << nerr << " " << messg << std::endl;
    if (level == 2) std::exit(1);
    return;
}

}; // end namespace quadpack