#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "LArMathMinimizer.hh"

using namespace std;

LArMathMinimizer::LArMathMinimizer()
{}

LArMathMinimizer::~LArMathMinimizer()
{}

void LArMathMinimizer::Initialize()
{
    LArRindex::Initialize();
    LArTrans::Initialize();
    LArGroupVelocity::Initialize();

}


double LArMathMinimizer::GetChi2(const double* xx)
{
    LArRindex::setp0(xx[0]);
    LArRindex::setp1(xx[1]);
    LArRindex::setp2(xx[2]);
    LArRindex::setrhoratio(xx[12]);

    LArTrans::setdelta(xx[3]);
    LArTrans::setA1(xx[4]);
    LArTrans::setmu1(xx[5]);
    LArTrans::setsigma1(xx[6]);
    LArTrans::setA2(xx[7]);
    LArTrans::setmu2(xx[8]);
    LArTrans::setsigma2(xx[9]);

    // pull term
    //LArGroupVelocity::setnulambda(xx[10]);
    LArTrans::setnuf(xx[10]);

    double chi2 = 0;
    chi2 += LArRindex::GetChi2();
    //LArRindex::GetChi2();
    chi2 += LArTrans::GetChi2();
    
    chi2 += LArGroupVelocity::GetChi2();

    return chi2;
}


int LArMathMinimizer::Minimization()
{

    ROOT::Minuit2::Minuit2Minimizer minimum (ROOT::Minuit2::kMigrad);

    // set tolerance , etc...
    minimum.SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum.SetMaxIterations(10000);  // for GSL
    minimum.SetTolerance(0.001);
    minimum.SetPrintLevel(1);

    // function wrapper for minimizer
    ROOT::Math::Functor f(&GetChi2, 11);
    double step[11];
    for (int i=0; i<11; i++) {
        step[i] = 0.001;
    }

    // start point
    double variable[11];
    variable[0] = 0.335;
    variable[1] = 0.099;
    variable[2] = 0.008;
    variable[3] = 0.2;
    variable[4] = 0.94;
    variable[5] = 126;
    variable[6] = 1;
    variable[7] = 0.4;
    variable[8] = 140;
    variable[9] = 1.5;
    variable[10] = 0;

    minimum.SetFunction(f);

    // Set the free variables to be minimized !
    minimum.SetVariable(0, "p0", variable[0], step[0]);
    minimum.SetVariable(1, "p1", variable[1], step[1]);
    minimum.SetVariable(2, "p2", variable[2], step[2]);
    minimum.SetVariable(3, "delta", variable[3], step[3]);
    minimum.SetVariable(4, "peakratio", variable[4], step[4]);
    minimum.SetVariable(5, "mu1", variable[5], step[5]);
    minimum.SetVariable(6, "sigma1", variable[6], step[6]);
    minimum.SetVariable(7, "A2", variable[7], step[7]);
    minimum.SetVariable(8, "mu2", variable[8], step[8]);
    minimum.SetVariable(9, "sigma2", variable[9], step[9]);
    minimum.SetVariable(10, "nu_f", variable[10], step[10]);

    // do the minimization
    minimum.Minimize();

    /*
    // Cov Matrix
    cout << "\n";
    cout << "Cov Matrix : " << std::endl; 
    double hess[11][11];
    for (int i=0; i<11; i++) {
        for (int j=0; j<11; j++) {
            hess[i][j] = minimum.CovMatrix(i, j) ;
            cout << hess[i][j] << " ";
        }
        cout << "\n" ;
    }
    */

    cout << "\n";
    cout << "Hessian Matrix : " << std::endl; 
    const int num = 11*11;
    double  h[num];
    minimum.GetHessianMatrix(h);
    for (int i=0; i<num; i++) {
        cout << h[i] << " ";
        if ( (i+1)%11 == 0)
            cout << "\n";
    }
    
    //const double *xs = minimum.X();
    //std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
    //          << minimum.MinValue() << std::endl;

    //// expected minimum is 0
    //if (minimum.MinValue() < 1.E-4 && f(xs) < 1.E-4)
    //    std::cout << "Minimizer " << minName << " - " << algoName
    //              << "   converged to the right minimum" << std::endl;
    //else
    //{
    //    std::cout << "Minimizer " << minName << " - " << algoName
    //              << "   failed to converge !!!" << std::endl;
    //    Error("NumericalMinimization", "fail to converge");
    //}

   return 0;
}
