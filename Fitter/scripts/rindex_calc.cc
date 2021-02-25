#include <iostream>
#include <cmath>

using namespace std;

double rindex_our(double *par)
{
    double fa = par[0];
    double fb = par[1];
    double fc = par[2];
    double A = 1.2055e-2*2/3.;
    double rho_ratio = 34.49/(44.66e-3);
    double l  = 0.128;
    double l1 = 91.012;
    double l2 = 89.892;
    double l3 = 214.02;

    return sqrt((3/(1-(A*rho_ratio*(fa/(l1-1/l/l)+fb/(l2-1/l/l)+fc/(l3-1/l/l)))))-2);
}

double rindex_ba(double *par)
{
    double a0  = par[0];
    double aUV = par[1];
    double aIR = par[2];
    double l = par[3];
    double lUV = 0.1066;
    double lIR = 0.9083;

    double A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) ;

    return sqrt(1+3*A/(3-A));
}

double raylength(double l, double rindex, double delta)
{
    double kT = 2.24442E-9;
    double kB = 1.380649E-23;
    double T = 90; // K
    double f = 1e22;

    double pi = 3.141592653;
    double rayL = 1 / (8*pow(pi, 3)/3/pow(l, 4)
                * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f * (6+3*delta)/(6-7*delta));

    return rayL;
}

int main(int argc, char *argv[])
{
    if (argc != 7) {
        cout << "Not Correct Parameter Number !" << endl;
        return 0;
    }

    double par[4] = {atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[6]) };
    double rindex;
    if (atoi(argv[1]) == 1) {
        cout << "rindex @128nm is " << rindex_our(par) << endl; rindex = rindex_our(par);
    } else if (atoi(argv[1]) == 0) {
        cout << "rindex @128nm is " << rindex_ba(par) << endl; rindex = rindex_ba(par); }

    cout << "raylength @128nm is " << raylength(atof(argv[6]), rindex, atof(argv[5])) << endl;

    return 1;
}
