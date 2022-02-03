#include"Lambert.h"
#include"Vec.h"

#include<cmath>

double hypergeometricF(double z, double tol)
{
    double Sj = 1.0;
    double Cj = 1.0;
    double err = 1.0;
    double Cj1 = 0.0;
    double Sj1 = 0.0;
    double j = 0;
    while (err > tol)
    {
        Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1);
        Sj1 = Sj + Cj1;
        err = fabs(Cj1);
        Sj = Sj1;
        Cj = Cj1;
        j += 1.0;
    }
    return Sj;
}

void dTdx(double &DT, double &DDT, double &DDDT, double x, double T, double l2, double l3)
{
    double umx2 = 1.0 - x * x;
    double y = sqrt(1.0 - l2 * umx2);
    double y2 = y * y;
    double y3 = y2 * y;
    DT = 1.0 / umx2 * (3.0 * T * x - 2.0 + 2.0 * l3 * x / y);
    DDT = 1.0 / umx2 * (3.0 * T + 5.0 * x * DT + 2.0 * (1.0 - l2) * l3 / y3);
    DDDT = 1.0 / umx2 * (7.0 * x * DDT + 8.0 * DT - 6.0 * (1.0 - l2) * l2 * l3 * x / y3 / y2);
}

void x2tof2(double &tof, double x, double lambda, double l2)
{
    double a = 1.0 / (1.0 - x * x);
    if (a > 0) // ellipse
    {
        double alfa = 2.0 * acos(x);
        double beta = 2.0 * asin(sqrt(l2 / a));
        if (lambda < 0.0) beta = -beta;
        tof = ((a * sqrt(a) * ((alfa - sin(alfa)) - (beta - sin(beta)))) / 2.0);
    }
    else
    {
        double alfa = 2.0 * acosh(x);
        double beta = 2.0 * asinh(sqrt(-l2 / a));
        if (lambda < 0.0) beta = -beta;
        tof = (-a * sqrt(-a) * ((beta - sinh(beta)) - (alfa - sinh(alfa))) / 2.0);
    }
}

void x2tof(double &tof, const double x, double lambda, double l2)
{
    double battin = 0.01;
    double lagrange = 0.2;
    double dist = fabs(x - 1);
    if (dist < lagrange && dist > battin)
    { // We use Lagrange tof expression
        x2tof2(tof, x, lambda, l2);
        return;
    }
    double K = l2;
    double E = x * x - 1.0;
    double rho = fabs(E);
    double z = sqrt(1 + K * E);
    if (dist < battin)
    { // We use Battin series tof expression
        double eta = z - lambda * x;
        double S1 = 0.5 * (1.0 - lambda - x * eta);
        double Q = hypergeometricF(S1, 1e-11);
        Q = 4.0 / 3.0 * Q;
        tof = (eta * eta * eta * Q + 4.0 * lambda * eta) / 2.0;
        return;
    }
    else
    { // We use Lancaster tof expresion
        double y = sqrt(rho);
        double g = x * z - lambda * E;
        double d = 0.0;
        if (E < 0)
        {
            double l = acos(g);
            d = l;
        }
        else
        {
            double f = y * (z - lambda * x);
            d = log(f + g);
        }
        tof = (x - lambda * z - d / y) / E;
        return;
    }
}

static void Householder(double T, double &x0, double eps, unsigned int itrs, double lambda, double l2, double l3)
{
    unsigned int it = 0;
    double err = 1.0;
    double xnew = 0.0;
    double tof = 0.0, delta = 0.0, DT = 0.0, DDT = 0.0, DDDT = 0.0;

    while ((err > eps) && (it < itrs))
    {
        x2tof(tof, x0, lambda, l2);
        dTdx(DT, DDT, DDDT, x0, tof, l2, l3);
        delta = tof - T;
        double DT2 = DT * DT;
        xnew = x0 - delta * (DT2 - delta * DDT / 2.0) / (DT * (DT2 - delta * DDT) + DDDT * delta * delta / 6.0);
        err = fabs(x0 - xnew);
        x0 = xnew;
        ++it;
    }
}

void LambertDIzzo(const Vec3 &r1, const Vec3 &r2, double tof, double mu, Vec3 &v1, Vec3 &v2)
{
    Vec3 c = r2 - r1;
    double C = Vec3::Length(c);
    double R1 = Vec3::Length(r1);
    double R2 = Vec3::Length(r2);
    double s = 0.5 * (R1 + R2 + C);
    Vec3 ir1 = r1 / R1;
    Vec3 ir2 = r2 / R2;
    Vec3 ih = Vec3::Cross(ir1, ir2);
    ih = ih / Vec3::Length(ih);
    double ll = 1 - C / s;
    double l = sqrt(ll);
    double lll = ll * l;

    Vec3 it1, it2;

    if (ih.z < 0)
    {
        l = -l;
        it1 = Vec3::Cross(ir1, ih);
        it2 = Vec3::Cross(ir2, ih);
    }
    else
    {
        it1 = Vec3::Cross(ih, ir1);
        it2 = Vec3::Cross(ih, ir2);
    }
        
    it1 = it1 / Vec3::Length(it1);
    it2 = it2 / Vec3::Length(it2);
        
    double T = sqrt(2 * mu / (s * s * s)) * tof;

    double T00 = acos(l) + l * sqrt(1 - ll);
    double T1 = 2.0 / 3.0 * (1 - lll);

    double x = 0.0;

    if (T >= T00)
        x = -(T - T00) / (T - T00 + 4);
    else if (T <= T1)
        x = 2.5 * T1 * (T1 - T) / ((1 - ll * lll) * T) + 1;
    else
        x = pow(T / T00, 0.69314718055994529 / log(T1 / T00)) - 1;

    Householder(T, x, 1e-5, 15, l, ll, lll);

    double gamma = sqrt(0.5 * mu * s);
    double rho = (R1 - R2) / C;
    double sigma = sqrt(1 - rho * rho);

    double y = sqrt(1 - ll + ll * x * x);
    double arg1 = l * y - x;
    double arg2 = l * y + x;
    double arg3 = y + l * x;
    double Vr1 = gamma * (arg1 - rho * arg2) / R1;
    double Vr2 = -gamma * (arg1 + rho * arg2) / R2;
    double Vt = gamma * sigma * arg3;
    double Vt1 = Vt / R1;
    double Vt2 = Vt / R2;
    v1 = ir1 * Vr1 + it1 * Vt1;
    v2 = ir2 * Vr2 + it2 * Vt2;
}

static double C2(double phi)
{
    return (1 - cos(sqrt(phi))) / phi;
}

static double C3(double phi)
{
    double sphi = sqrt(phi);
    return (sphi - sin(sphi)) / (sphi * phi);
}

void LambertUV(const Vec3 &r1, const Vec3 &r2, double tof, double mu, Vec3 &v1, Vec3 &v2)
{
    double smu = sqrt(mu);

    double R1 = Vec3::Length(r1);
    double R2 = Vec3::Length(r2);

    double gamma = Vec3::Dot(r1, r2) / R1 / R2;

    double A = sqrt(R1 * R2 * (1 + gamma));

    if (A == 0)
    {
        v1 = Vec3();
        v2 = Vec3();
        return;
    }

    double c2 = 0.5;
    double c3 = 1.0 / 6.0;

    double psi = 0.0;
    double psi_u = 4 * M_PI * M_PI;
    double psi_l = -psi_u;

    unsigned int itrs = 200;
    double tol = 1e-4;
    for (unsigned int i = 0; i < itrs; ++i)
    {
        double B = R1 + R2 + A * (psi * c3 - 1) / sqrt(c2);

        if (A > 0.0 && B < 0.0)
        {
            psi_l += M_PI;
            B = -B;
        }

        double chi = sqrt(B / c2);
        double chi3 = chi * chi * chi;
        double dt_ = (chi3 * c3 + A * sqrt(B)) / smu;

        if (fabs(tof - dt_) < tol)
        {
            double f = 1 - B / R1;
            double g = A * sqrt(B / mu);
            double g_d = 1 - B / R2;

            v1 = (r2 - r1 * f) / g;
            v2 = (r2 * g_d - r1) / g;
            return;
        }

        if (dt_ < tof) psi_l = psi;
        else psi_u = psi;

        psi = (psi_l + psi_u) * 0.5;
        c2 = C2(psi);
        c3 = C3(psi);
    }

    v1 = Vec3();
    v2 = Vec3();
}