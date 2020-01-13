/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpecialFunctions.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

double SpecialFunctions::lngamma(double a)
{
    int j;
    double xx, y, tmp, ser;
    static const double cof[6] = {76.18009172947146,  -86.50532032941677,    24.01409824083091,
                                  -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
    y = xx = a;
    tmp = xx + 5.5;
    tmp -= (xx + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j < 6; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / xx);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::gamma(double a)
{
    return exp(lngamma(a));
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::incompleteGamma(double a, double x)
{
    const int ITMAX = 100;
    const double EPS = 3.0e-7;
    if (x < 0.0 || a <= 0.0) throw FATALERROR("Invalid arguments in routine SpecialFunctions::incompletegamma");
    double gln = lngamma(a);
    if (x == 0) return 0.0;
    if (x < (a + 1.0))
    {
        double sum;
        double ap = a;
        double del = sum = 1.0 / a;
        for (int n = 1; n <= ITMAX; n++)
        {
            ++ap;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum) * EPS) return sum * exp(-x + a * log(x) - gln);
        }
        throw FATALERROR("a too large, ITMAX too small in routine SpecialFunctions::incompletegamma");
    }
    else
    {
        const double FPMIN = 1.0e-30;
        double an, del;
        double b = x + 1.0 - a;
        double c = 1.0 / FPMIN;
        double d = 1.0 / b;
        double h = d;
        int i;
        for (i = 1; i <= ITMAX; i++)
        {
            an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (fabs(d) < FPMIN) d = FPMIN;
            c = b + an / c;
            if (fabs(c) < FPMIN) c = FPMIN;
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if (fabs(del - 1.0) < EPS) break;
        }
        if (i > ITMAX) throw FATALERROR("a too large, ITMAX too small in SpecialFunctions::incompletegamma");
        return 1.0 - exp(-x + a * log(x) - gln) * h;
    }
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::beta(double a, double b)
{
    return exp(lngamma(a) + lngamma(b) - lngamma(a + b));
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::beta(double x, double a, double b)
{
    return betaRegularized(x, a, b) * beta(a, b);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::betaRegularized(double x, double a, double b)
{
    double bt;
    if (x < 0.0 || x > 1.0) throw FATALERROR("x should be between 0 and 1 (x = " + std::to_string(x) + ")");
    if (x == 0.0 || x == 1.0)
        bt = 0.0;
    else
        bt = exp(lngamma(a + b) - lngamma(a) - lngamma(b) + a * log(x) + b * log(1.0 - x));
    double anew, bnew, xnew, voorfac, cterm;
    if (x < (a + 1.0) / (a + b + 2.0))
    {
        voorfac = bt / a;
        cterm = 0.0;
        anew = a;
        bnew = b;
        xnew = x;
    }
    else
    {
        voorfac = -bt / b;
        cterm = 1.0;
        anew = b;
        bnew = a;
        xnew = 1.0 - x;
    }
    const int MAXIT = 100;
    const double EPS = std::numeric_limits<double>::epsilon();
    const double FPMIN = std::numeric_limits<double>::min() / EPS;
    double qab = anew + bnew;
    double qap = anew + 1.0;
    double qam = anew - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * xnew / qap;
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;
    int m;
    for (m = 1; m <= MAXIT; m++)
    {
        int m2 = m * 2;
        double aa = m * (bnew - m) * xnew / ((qam + m2) * (anew + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -(anew + m) * (qab + m) * xnew / ((anew + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (fabs(del - 1.0) <= EPS) break;
    }
    if (m > MAXIT) throw FATALERROR("Maximum number of iterations reached");
    return cterm + voorfac * h;
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::bessi1(double x)
{
    double ax, ans, y;
    if ((ax = fabs(x)) < 3.75)
    {
        y = x / 3.75;
        y *= y;
        ans = ax
              * (0.5
                 + y
                       * (0.87890594
                          + y
                                * (0.51498869
                                   + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
    }
    else
    {
        y = 3.75 / x;
        ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
        ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
        ans *= (exp(ax) / sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::bessk1(double x)
{
    double ans, y;
    if (x <= 0.0) throw FATALERROR("x should be positive (x = " + std::to_string(x) + ")");
    if (x <= 2.0)
    {
        y = x * x / 4.0;
        ans = (log(x / 2.0) * bessi1(x))
              + (1.0 / x)
                    * (1.0
                       + y
                             * (0.15443144
                                + y
                                      * (-0.67278579
                                         + y
                                               * (-0.18156897
                                                  + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
    }
    else
    {
        y = 2.0 / x;
        ans = (exp(-x) / sqrt(x))
              * (1.25331414
                 + y
                       * (0.23498619
                          + y
                                * (-0.3655620e-1
                                   + y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3)))))));
    }
    return ans;
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::dawson(double x)
{
    const int NMAX = 6;
    const double H = 0.4, a1 = 2.0 / 3.0, a2 = 0.4, a3 = 2.0 / 7.0;
    int i, n0;
    static bool init = true;
    double d1, d2, e1, e2, sum, x2, xp, xx, ans;
    static double c[NMAX];

    if (init)
    {
        init = false;
        for (i = 0; i < NMAX; i++)
        {
            double qq = (2.0 * i + 1.0) * H;
            c[i] = exp(-qq * qq);
        }
    }
    if (fabs(x) < 0.2)
    {
        x2 = x * x;
        ans = x * (1.0 - a1 * x2 * (1.0 - a2 * x2 * (1.0 - a3 * x2)));
    }
    else
    {
        xx = fabs(x);
        n0 = 2 * int(0.5 * xx / H + 0.5);
        xp = xx - n0 * H;
        e1 = exp(2.0 * xp * H);
        e2 = e1 * e1;
        d1 = n0 + 1.0;
        d2 = d1 - 2.0;
        sum = 0.0;
        for (i = 0; i < NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2) sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));
        ans = 0.5641895835 * exp(-xp * xp) * sum;
        if (x <= 0) ans *= -1;
    }
    return ans;
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::rf(double x, double y, double z)
{
    const double ERRTOL = 0.0025, TINY = 1.5e-38, BIG = 3.0e37, THIRD = 1.0 / 3.0;
    const double C1 = 1.0 / 24.0, C2 = 0.1, C3 = 3.0 / 44.0, C4 = 1.0 / 14.0;
    double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt, zt;

    if (x < 0.0 || y < 0.0 || z < 0.0 || x + y < TINY || x + z < TINY || y + z < TINY || x > BIG || y > BIG || z > BIG)
        throw FATALERROR("Invalid argument(s)");
    xt = x;
    yt = y;
    zt = z;
    do
    {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        ave = THIRD * (xt + yt + zt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
    } while (fabs(delx) > ERRTOL || fabs(dely) > ERRTOL || fabs(delz) > ERRTOL);
    e2 = delx * dely - delz * delz;
    e3 = delx * dely * delz;
    return (1.0 + (C1 * e2 - C2 - C3 * e3) * e2 + C4 * e3) / sqrt(ave);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::rd(double x, double y, double z)
{
    const double ERRTOL = 0.0015, TINY = 1.0e-25, BIG = 4.5e21;
    const double C1 = 3.0 / 14.0, C2 = 1.0 / 6.0, C3 = 9.0 / 22.0;
    const double C4 = 3.0 / 26.0, C5 = 0.25 * C3, C6 = 1.5 * C4;
    double alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx, sqrty, sqrtz, sum, xt, yt, zt;

    if (x < 0.0 || y < 0.0 || x + y < TINY || z < TINY || x > BIG || y > BIG || z > BIG)
        throw FATALERROR("Invalid argument(s)");
    xt = x;
    yt = y;
    zt = z;
    sum = 0.0;
    fac = 1.0;
    do
    {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        sum += fac / (sqrtz * (zt + alamb));
        fac *= 0.25;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        ave = 0.2 * (xt + yt + 3.0 * zt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
    } while (fabs(delx) > ERRTOL || fabs(dely) > ERRTOL || fabs(delz) > ERRTOL);
    ea = delx * dely;
    eb = delz * delz;
    ec = ea - eb;
    ed = ea - 6.0 * eb;
    ee = ed + ec + ec;
    return 3.0 * sum
           + fac * (1.0 + ed * (-C1 + C5 * ed - C6 * delz * ee) + delz * (C2 * ee + delz * (-C3 * ec + delz * C4 * ea)))
                 / (ave * sqrt(ave));
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::rj(double x, double y, double z, double p)
{
    const double ERRTOL = 0.0015, TINY = 2.5e-13, BIG = 9.0e11;
    const double C1 = 3.0 / 14.0, C2 = 1.0 / 3.0, C3 = 3.0 / 22.0, C4 = 3.0 / 26.0;
    const double C5 = 0.75 * C3, C6 = 1.5 * C4, C7 = 0.5 * C2, C8 = C3 + C3;
    double a = 0.0, alamb, alpha, ans, ave, b = 0.0, beta, delp, delx, dely, delz, ea, eb, ec, ed, ee, fac, pt,
           rcx = 0.0, rho, sqrtx, sqrty, sqrtz, sum, tau, xt, yt, zt;

    if (x < 0.0 || y < 0.0 || z < 0.0 || x + y < TINY || x + z < TINY || y + z < TINY || fabs(p) < TINY || x > BIG
        || y > BIG || z > BIG || fabs(p) > BIG)
        throw FATALERROR("Invalid argument(s)");
    sum = 0.0;
    fac = 1.0;
    if (p > 0.0)
    {
        xt = x;
        yt = y;
        zt = z;
        pt = p;
    }
    else
    {
        if (x <= y && x <= z)
            xt = x;
        else if (y <= z)
            xt = y;
        else
            xt = z;
        if (x >= y && x >= z)
            zt = x;
        else if (y >= z)
            zt = y;
        else
            zt = z;
        yt = x + y + z - xt - zt;
        a = 1.0 / (yt - p);
        b = a * (zt - yt) * (yt - xt);
        pt = yt + b;
        rho = xt * zt / yt;
        tau = p * pt / yt;
        rcx = rc(rho, tau);
    }
    do
    {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        alpha = pow(pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz, 2);
        beta = pt * pow(pt + alamb, 2);
        sum += fac * rc(alpha, beta);
        fac *= 0.25;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        pt = 0.25 * (pt + alamb);
        ave = 0.2 * (xt + yt + zt + pt + pt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
        delp = (ave - pt) / ave;
    } while (fabs(delx) > ERRTOL || fabs(dely) > ERRTOL || fabs(delz) > ERRTOL || fabs(delp) > ERRTOL);
    ea = delx * (dely + delz) + dely * delz;
    eb = delx * dely * delz;
    ec = delp * delp;
    ed = ea - 3.0 * ec;
    ee = eb + 2.0 * delp * (ea - ec);
    ans = 3.0 * sum
          + fac
                * (1.0 + ed * (-C1 + C5 * ed - C6 * ee) + eb * (C7 + delp * (-C8 + delp * C4))
                   + delp * ea * (C2 - delp * C3) - C2 * delp * ec)
                / (ave * sqrt(ave));
    if (p <= 0.0) ans = a * (b * ans + 3.0 * (rcx - rf(xt, yt, zt)));
    return ans;
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::rc(double x, double y)
{
    const double ERRTOL = 0.0012, TINY = 1.69e-38, SQRTNY = 1.3e-19, BIG = 3.0e37;
    const double TNBG = TINY * BIG, COMP1 = 2.236 / SQRTNY, COMP2 = TNBG * TNBG / 25.0;
    const double THIRD = 1.0 / 3.0, C1 = 0.3, C2 = 1.0 / 7.0, C3 = 0.375, C4 = 9.0 / 22.0;
    double alamb, ave, s, w, xt, yt;

    if (x < 0.0 || y == 0.0 || (x + fabs(y)) < TINY || (x + fabs(y)) > BIG || (y < -COMP1 && x > 0.0 && x < COMP2))
        throw FATALERROR("Invalid argument(s)");
    if (y > 0.0)
    {
        xt = x;
        yt = y;
        w = 1.0;
    }
    else
    {
        xt = x - y;
        yt = -y;
        w = sqrt(x) / sqrt(xt);
    }
    do
    {
        alamb = 2.0 * sqrt(xt) * sqrt(yt) + yt;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        ave = THIRD * (xt + yt + yt);
        s = (yt - ave) / ave;
    } while (fabs(s) > ERRTOL);
    return w * (1.0 + s * s * (C1 + s * (C2 + s * (C3 + s * C4)))) / sqrt(ave);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::EllipticF(double z, double k)
{
    double a = (1.0 - z) * (1.0 + z);
    double b = (1.0 - z * k) * (1.0 + z * k);
    double c = 1.0;
    return z * rf(a, b, c);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::EllipticE(double z, double k)
{
    double a = (1.0 - z) * (1.0 + z);
    double b = (1.0 - z * k) * (1.0 + z * k);
    double c = 1.0;
    return z * (rf(a, b, c) - z * z * k * k * rd(a, b, c) / 3.0);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::EllipticPi(double z, double nu, double k)
{
    double a = (1.0 - z) * (1.0 + z);
    double b = (1.0 - z * k) * (1.0 + z * k);
    double c = 1.0;
    double nuz2 = -nu * z * z;
    double d = 1.0 + nuz2;
    return z * (rf(a, b, c) - nuz2 * rj(a, b, c, d) / 3.0);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::EllipticK(double k)
{
    return EllipticF(1.0, k);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::EllipticE(double k)
{
    return EllipticE(1.0, k);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::EllipticPi(double nu, double k)
{
    return EllipticPi(1.0, nu, k);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::functionX(double s)
{
    const double s0 = 0.01;
    double ans;

    double as = fabs(s);
    if (as == 0.0) throw FATALERROR("This function diverges for s = 0");
    if (as < s0)
    {
        double z = log(as / 2.0);
        double s2 = as * as;
        ans = -z - (1.0 + 2.0 * z) * s2 / 4.0 - (7.0 + 12 * z) * s2 * s2 / 32.0;
    }
    else if (as < 1.0 - s0)
    {
        double t = sqrt((1.0 - as) * (1.0 + as));
        ans = log((1.0 + t) / as) / t;
    }
    else if (as < 1.0 + s0)
    {
        double z = as - 1.0;
        ans = 1.0 + z * (-2.0 / 3.0 + z * (7.0 / 15.0 + z * (-12.0 / 35.0 + z * 83.0 / 315.0)));
    }
    else
    {
        ans = acos(1.0 / as) / sqrt((as - 1.0) * (as + 1.0));
    }
    return ans;
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::LambertW(double z)
{
    int i;
    const double eps = 4.0e-16;
    const double em1 = 0.3678794411714423215955237701614608;
    double p, e, t, w;
    if (z < -em1 || std::isinf(z) || std::isnan(z)) throw FATALERROR("Bad argument (z = " + std::to_string(z) + ")");
    if (z == 0.0) return 0.0;
    if (z < -em1 + 1e-4)  // series near -em1 in sqrt(q)
    {
        double q = z + em1, r = sqrt(q), q2 = q * q, q3 = q2 * q;
        return -1.0 + 2.331643981597124203363536062168 * r - 1.812187885639363490240191647568 * q
               + 1.936631114492359755363277457668 * r * q - 2.353551201881614516821543561516 * q2
               + 3.066858901050631912893148922704 * r * q2 - 4.175335600258177138854984177460 * q3
               + 5.858023729874774148815053846119 * r * q3
               - 8.401032217523977370984161688514 * q3 * q;  // error approx 1e-16
    }
    // initial approx for iteration...
    if (z < 1.0)  // series near 0
    {
        p = sqrt(2.0 * (2.7182818284590452353602874713526625 * z + 1.0));
        w = -1.0 + p * (1.0 + p * (-0.333333333333333333333 + p * 0.152777777777777777777777));
    }
    else  // asymptotic
    {
        w = log(z);
    }
    if (z > 3.0) w -= log(w);  // useful?
    // Halley iteration
    for (i = 0; i < 10; i++)
    {
        e = exp(w);
        t = w * e - z;
        p = w + 1.0;
        t /= e * p - 0.5 * (p + 1.0) * t / p;
        w -= t;
        if (fabs(t) < eps * (1.0 + fabs(w))) return w;  // rel-abs error
    }
    throw FATALERROR("No convergence at z = " + std::to_string(z));
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::LambertW1(double z)
{
    const double eps = 1.0e-12;
    const double em1 = 0.3678794411714423215955237701614608;
    static const double c[12] = {-1.0,
                                 2.331643981597124203363536062168,
                                 -1.812187885639363490240191647568,
                                 1.936631114492359755363277457668,
                                 -2.353551201881614516821543561516,
                                 3.066858901050631912893148922704,
                                 -4.175335600258177138854984177460,
                                 5.858023729874774148815053846119,
                                 -8.401032217523977370984161688514,
                                 12.250753501314460424,
                                 -18.100697012472442755,
                                 27.029044799010561650};
    if (z < -em1 || z > 0.0 || std::isinf(z) || std::isnan(z))
        throw FATALERROR("Bad argument (z = " + std::to_string(z) + ")");
    if (z == 0.0) return -DBL_MAX;
    double q = z + em1;
    double r = -sqrt(q);
    double t8 = c[8] + r * (c[9] + r * (c[10] + r * c[11]));
    double t5 = c[5] + r * (c[6] + r * (c[7] + r * t8));
    double t1 = c[1] + r * (c[2] + r * (c[3] + r * (c[4] + r * t5)));
    double w0 = c[0] + r * t1;
    if (q < 3.0e-3)
        return w0;
    else
    {
        double w, e, p, t;
        if (z < -1e-6)
            w = w0;
        else
        {
            double l1 = log(-z);
            double l2 = log(-l1);
            w = l1 - l2 + l2 / l1;
        }
        for (int i = 0; i < 10; i++)
        {
            e = exp(w);
            t = w * e - z;
            p = w + 1.0;
            t /= e * p - 0.5 * (p + 1.0) * t / p;
            w -= t;
            if (fabs(t) < eps * (1.0 + fabs(w))) return w;
        }
    }
    throw FATALERROR("No convergence at z = " + std::to_string(z));
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::DebyeD(int n, double x)
{
    if (x < 0.0 || n < 1 || n > 20)
        throw FATALERROR("Bad arguments (x = " + std::to_string(x) + ", n = " + std::to_string(n) + ")");

    // 1/(2pi)
    const double M_1_2PI = .159154943091895335768883763373;

    // for values up to 4.80 the list of zeta functions and the sum up to k < K are huge enough to gain
    // numeric stability in the sum
    if (x >= 3.)
    {
        // list of n! zeta(n+1) for n =0 up to the maximum n implemented.
        // Limited by the cancellation of digits encountered at smaller x and larger n.
        // Digits := 30 :
        // for n from 1 to 30 do
        //   printf("%.29e, ", evalf(n!*Zeta(n+1))) ;
        // od:
        static const double nzetan[] = {0.,
                                        1.64493406684822643647241516665e+00,
                                        2.40411380631918857079947632302e+00,
                                        6.49393940226682914909602217926e+00,
                                        2.48862661234408782319527716750e+01,
                                        1.22081167438133896765742151575e+02,
                                        7.26011479714984435324654235892e+02,
                                        5.06054987523763947046857360209e+03,
                                        4.04009783987476348853278236554e+04,
                                        3.63240911422382626807143525567e+05,
                                        3.63059331160662871299061884284e+06,
                                        3.99266229877310867023270732405e+07,
                                        4.79060379889831452426876764501e+08,
                                        6.22740219341097176419285340896e+09,
                                        8.71809578301720678451912203103e+10,
                                        1.30769435221891382089009990749e+12,
                                        2.09229496794815109066316556880e+13,
                                        3.55688785859223715975612396717e+14,
                                        6.40238592281892140073564945334e+15,
                                        1.21645216453639396669876696274e+17,
                                        2.43290316850786132173725681824e+18,
                                        5.10909543543702856776502748606e+19,
                                        1.12400086178089123060215294900e+21};

        // n!*zeta(n) is the integral for x=infinity, 27.1.3
        double sum = nzetan[n];

        // the number of terms needed in the k-sum for x=0,1,2,3...
        // reflects the n=1 case, because higher n need less terms.
        static int kLim[] = {0, 0, 0, 13, 10, 8, 7, 6, 5, 5, 4, 4, 4, 3};
        const int kLimSize = sizeof(kLim) / sizeof(kLim[0]);
        const int kmax = (static_cast<int>(x) < kLimSize) ? kLim[static_cast<int>(x)] : 3;
        // Abramowitz Stegun 27.1.2
        for (int k = 1; k <= kmax; k++)
        {
            // do not use x(k+1)=xk+x to avoid loss of precision
            const double xk = x * k;
            double ksum = 1. / xk;
            double tmp = n * ksum / xk;  // n/(xk)^2
            for (int s = 1; s <= n; s++)
            {
                ksum += tmp;
                tmp *= (n - s) / xk;
            }
            sum -= exp(-xk) * ksum * pow(x, n + 1.);
        }
        return sum * n / pow(x, n);
    }
    else
    {
        // list of absolute values of Bernoulli numbers of index 2*k,
        // multiplied by (2*pi)^k/(2k)!, and 2 subtracted, k=0,1,2,3,4
        // Digits := 60 :
        // interface(prettyprint=0) :
        // for k from 1 to 70 do
        //   printf("%.30e,\n",evalf( abs((2*Pi)^(2*k)*bernoulli(2*k)/(2*k)!)-2 )) ;
        // od;
        static const double koeff[] = {0.,
                                       1.289868133696452872944830333292e+00,
                                       1.646464674222763830320073930823e-01,
                                       3.468612396889827942903585958184e-02,
                                       8.154712395888678757370477017305e-03,
                                       1.989150255636170674291917800638e-03,
                                       4.921731066160965972759960954793e-04,
                                       1.224962701174096585170902102707e-04,
                                       3.056451881730374346514297527344e-05,
                                       7.634586529999679712923289243879e-06,
                                       1.907924067745592226304077366899e-06,
                                       4.769010054554659800072963735060e-07,
                                       1.192163781025189592248804158716e-07,
                                       2.980310965673008246931701326140e-08,
                                       7.450668049576914109638408036805e-09,
                                       1.862654864839336365743529470042e-09,
                                       4.656623667353010984002911951881e-10,
                                       1.164154417580540177848737197821e-10,
                                       2.910384378208396847185926449064e-11,
                                       7.275959094757302380474472711747e-12,
                                       1.818989568052777856506623677390e-12,
                                       4.547473691649305030453643155957e-13,
                                       1.136868397525517121855436593505e-13,
                                       2.842170965606321353966861428348e-14,
                                       7.105427382674227346596939068119e-15,
                                       1.776356842186163180619218277278e-15,
                                       4.440892101596083967998640188409e-16,
                                       1.110223024969096248744747318102e-16,
                                       2.775557561945046552567818981300e-17,
                                       6.938893904331845249488542992219e-18,
                                       1.734723476023986745668411013469e-18,
                                       4.336808689994439570027820336642e-19,
                                       1.084202172491329082183740080878e-19,
                                       2.710505431220232916297046799365e-20,
                                       6.776263578041593636171406200902e-21,
                                       1.694065894509399669649398521836e-21,
                                       4.235164736272389463688418879636e-22,
                                       1.058791184067974064762782460584e-22,
                                       2.646977960169798160618902050189e-23,
                                       6.617444900424343177893912768629e-24,
                                       1.654361225106068880734221123349e-24,
                                       4.135903062765153408791935838694e-25,
                                       1.033975765691286264082026643327e-25,
                                       2.584939414228213340076225223666e-26,
                                       6.462348535570530772269628236053e-27,
                                       1.615587133892632406631747637268e-27,
                                       4.038967834731580698317525293132e-28,
                                       1.009741958682895139216954234507e-28,
                                       2.524354896707237808750799932127e-29,
                                       6.310887241768094478219682436680e-30,
                                       1.577721810442023614704107565240e-30,
                                       3.944304526105059031370476640000e-31,
                                       9.860761315262647572437533499000e-32,
                                       2.465190328815661892443976898000e-32,
                                       6.162975822039154730370601500000e-33,
                                       1.540743955509788682510501190000e-33,
                                       3.851859888774471706184973900000e-34,
                                       9.629649721936179265360991000000e-35,
                                       2.407412430484044816328953000000e-35,
                                       6.018531076210112040809600000000e-36,
                                       1.504632769052528010200750000000e-36,
                                       3.761581922631320025497600000000e-37,
                                       9.403954806578300063715000000000e-38,
                                       2.350988701644575015901000000000e-38,
                                       5.877471754111437539470000000000e-39,
                                       1.469367938527859384580000000000e-39,
                                       3.673419846319648458500000000000e-40,
                                       9.183549615799121117000000000000e-41,
                                       2.295887403949780249000000000000e-41,
                                       5.739718509874450320000000000000e-42,
                                       1.434929627468612270000000000000e-42};

        double sum = 0.;
        // Abramowitz-Stegun 27.1.1
        const double x2pi = x * M_1_2PI;
        int koeffSize = sizeof(koeff) / sizeof(koeff[0]);
        for (int k = 1; k < koeffSize - 1; k++)
        {
            const double sumold = sum;
            // do not precompute x2pi^2 to avoid loss of precision
            sum += (2. + koeff[k]) * pow(x2pi, 2. * k) / (2 * k + n);
            k++;
            sum -= (2. + koeff[k]) * pow(x2pi, 2. * k) / (2 * k + n);
            if (sum == sumold) break;
        }
        sum += 1. / n - x / (2 * (1 + n));
        return sum * n;
    }
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::gln(double p, double x)
{
    const double q = 1.0 - p;
    if (q == 0.0)
        return log(x);
    else if (fabs(q) < 1e-3)
    {
        double lnx = log(x);
        double s = q * lnx;
        return lnx * (1.0 + 0.5 * s + 1.0 / 6.0 * s * s + 1.0 / 24.0 * s * s * s);
    }
    else
        return (pow(x, q) - 1.0) / q;
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::gln2(double p, double x1, double x2)
{
    return pow(x2, 1.0 - p) * gln(p, x1 / x2);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::gexp(double p, double x)
{
    const double q = 1.0 - p;
    if (q == 0.0)
        return exp(x);
    else if (fabs(q) < 1e-3)
    {
        double x2 = x * x;
        return exp(x)
               * (1.0 - 0.5 * x2 * q + 1.0 / 24.0 * x * x2 * (8.0 + 3.0 * x) * q * q
                  - 1.0 / 48.0 * x2 * x2 * (12.0 + 8.0 * x + x2) * q * q * q);
    }
    else
        return pow(1.0 + q * x, 1.0 / q);
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::lnmean(double x1, double x2)
{
    if (x1 > x2) std::swap(x1, x2);
    if (x1 <= 0) return 0.;

    double x = x2 / x1 - 1.;
    if (x < 1e-3)
    {
        return x1
               / (1. - 1. / 2. * x + 1. / 3. * x * x - 1. / 4. * x * x * x + 1. / 5. * x * x * x * x
                  - 1. / 6. * x * x * x * x * x);
    }
    else
    {
        return (x2 - x1) / log(x2 / x1);
    }
}

////////////////////////////////////////////////////////////////////

double SpecialFunctions::lnmean(double x1, double x2, double lnx1, double lnx2)
{
    if (x1 > x2)
    {
        std::swap(x1, x2);
        std::swap(lnx1, lnx2);
    }
    if (x1 <= 0) return 0.;

    double x = x2 / x1 - 1.;
    if (x < 1e-3)
    {
        return x1
               / (1. - 1. / 2. * x + 1. / 3. * x * x - 1. / 4. * x * x * x + 1. / 5. * x * x * x * x
                  - 1. / 6. * x * x * x * x * x);
    }
    else
    {
        return (x2 - x1) / (lnx2 - lnx1);
    }
}

////////////////////////////////////////////////////////////////////
