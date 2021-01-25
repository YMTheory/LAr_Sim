#include <iostream>
#include <cmath>

using namespace std;

double raylength(double l, double rindex, double delta)
{
    double kT = 2.24442E-9;
    double kB = 1.380649E-23;
    double T = 90; // K
    double f = 1e22;

    double pi = 3.141592653;
    double rayL = 1 / (8*pow(pi, 3)/3/pow(l, 4)
                * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f) ;
                //* ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f * (6+3*delta)/(6-7*delta));

    return rayL;
}


int main(int argc, char *argv[])
{
    if (argc != 3) {
        cout << "Not Correct Parameter Number !" << endl;
        return 0;
    }

    double rindex = atof(argv[1]);
    double delta  = atof(argv[2]);
    double rayL   = raylength(0.128, rindex, delta); 

    cout << "raylength @128nm " << rayL << endl;

    return 1;
}
