#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "LArMathMinimizer_new.hh"

using namespace std;

bool LArMathMinimizer_new::m_fit_purified = LArConfiguration::fit_purified;


LArMathMinimizer_new::LArMathMinimizer_new()
{}

LArMathMinimizer_new::~LArMathMinimizer_new()
{}

void LArMathMinimizer_new::Initialize()
{
    LArRindex_new::Initialize();
    LArTrans_new::Initialize();

}


double LArMathMinimizer_new::GetChi2(const double* xx)
{
    LArRindex_new::setrho(xx[0]);
    //LArRindex::setrhoratio(xx[12]);

    LArTrans_new::setdelta(xx[1]);
    LArTrans_new::setA1(xx[2]);
    LArTrans_new::setmu1(xx[3]);
    LArTrans_new::setsigma1(xx[4]);
    LArTrans_new::setA2(xx[5]);
    LArTrans_new::setmu2(xx[6]);
    LArTrans_new::setsigma2(xx[7]);

    // pull term
    //LArGroupVelocity::setnulambda(xx[10]);
    LArTrans_new::setnuf(xx[8]);
    LArTrans_new::settemp(xx[10]);
    LArRindex_new::seta0(xx[11]) ;
    LArRindex_new::setaUV(xx[12]);
    LArRindex_new::setaIR(xx[13]);
    LArTrans_new::setp0(xx[14]);
    LArTrans_new::setp1(xx[15]);


    // scale for purified spectra
    LArTrans_new::setscale(xx[9]);

    double chi2 = 0;
    chi2 += LArRindex_new::GetChi2();

    chi2 += LArTrans_new::GetChi2();
    
    return chi2;
}


int LArMathMinimizer_new::Minimization()
{

    //ROOT::Minuit2::Minuit2Minimizer_new minimum (ROOT::Minuit2::kMigrad);
    ROOT::Math::Minimizer* minimum = 
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(1);

    // function wrapper for Minimizer_new
    ROOT::Math::Functor f(&GetChi2, 16);
    double step[16];
    for (int i=0; i<16; i++) {
        step[i] = 0.001;
    }
    step[14] = 1e-12;
    step[15] = 1e-12;

    // start point
    double variable[16];
    variable[0] = 1;
    variable[1] = 0.2;
    variable[2] = 0.937;
    variable[3] = 126.51;
    variable[4] = 1;
    variable[5] = 0.4;
    variable[6] = 140.121;
    variable[7] = 1.537;
    variable[8] = 0;
    variable[9] = 1.00;
    variable[10] = 85;
    variable[11] = 0.335;
    variable[12] = 0.099;
    variable[13] = 0.008;
    variable[14] = p0;
    variable[15] = p1;

    minimum->SetFunction(f);

    // Set the free variables to be minimized !
    minimum->SetVariable(0, "rho", variable[0], step[0]);
    minimum->SetVariable(1, "delta", variable[1], step[1]);
    minimum->SetVariable(2, "peakratio", variable[2], step[2]);
    minimum->SetVariable(3, "mu1", variable[3], step[3]);
    minimum->SetVariable(4, "sigma1", variable[4], step[4]);
    minimum->SetVariable(5, "A2", variable[5], step[5]);
    minimum->SetVariable(6, "mu2", variable[6], step[6]);
    minimum->SetVariable(7, "sigma2", variable[7], step[7]);
    minimum->SetVariable(8, "nu_f", variable[8], step[8]);
    minimum->SetVariable(9, "scale", variable[9], step[9]);
    minimum->SetVariable(10, "temperature", variable[10], step[10]);
    minimum->SetVariable(11, "vara0", variable[11], step[11]);
    minimum->SetVariable(12, "varaUV", variable[12], step[12]);
    minimum->SetVariable(13, "varaIR", variable[13], step[13]);
    minimum->SetVariable(14, "varp0", variable[14], step[14]);
    minimum->SetVariable(15, "varp1", variable[15], step[15]);

    if (m_fit_purified) {
        minimum->FixVariable(2);
        minimum->FixVariable(3);
        minimum->FixVariable(4);
        minimum->FixVariable(6);
        minimum->FixVariable(7);
    } else {
        minimum->FixVariable(9);
    }
    
    //minimum->FixVariable(14);
    //minimum->FixVariable(15);

    // do the minimization
    minimum->Minimize();

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

    //cout << "\n";
    //cout << "Hessian Matrix : " << std::endl; 
    //const int num = 11*11;
    //double  h[num];
    //minimum.GetHessianMatrix(h);
    //for (int i=0; i<num; i++) {
    //    cout << h[i] << " ";
    //    if ( (i+1)%11 == 0)
    //        cout << "\n";
    //}
    
    //const double *xs = minimum.X();
    //std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
    //          << minimum.MinValue() << std::endl;

    //// expected minimum is 0
    //if (minimum.MinValue() < 1.E-4 && f(xs) < 1.E-4)
    //    std::cout << "Minimizer_new " << minName << " - " << algoName
    //              << "   converged to the right minimum" << std::endl;
    //else
    //{
    //    std::cout << "Minimizer_new " << minName << " - " << algoName
    //              << "   failed to converge !!!" << std::endl;
    //    Error("NumericalMinimization", "fail to converge");
    //}

   return 0;
}


bool LArMathMinimizer_new::Plot()
{
    LArRindex_new::Plot();
    LArTrans_new::Plot();
}



