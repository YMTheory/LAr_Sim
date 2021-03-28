#include "LArChiFunction_new.hh"
#include "LArConfiguration.hh"

#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2D.h>

using namespace std;

int LArChiFunction_new::usePull = LArConfiguration::use_pullterm;

double LArChiFunction_new::m_chi2;
double LArChiFunction_new::m_chi2Min;
bool LArChiFunction_new::m_DoFit;
int LArChiFunction_new::m_nParameters = 13;
double LArChiFunction_new::m_bestFit[20];
double LArChiFunction_new::m_fitError[20];
double LArChiFunction_new::m_ratio;

bool LArChiFunction_new::m_fit_purified = LArConfiguration::fit_purified;

// Lagrange Multipliers
double LArChiFunction_new::m_factor1 = LArConfiguration::factor1;
double LArChiFunction_new::m_factor2 = LArConfiguration::factor2;
double LArChiFunction_new::m_factor3 = LArConfiguration::factor3;

LArChiFunction_new::LArChiFunction_new()
{}

LArChiFunction_new::~LArChiFunction_new()
{}

void LArChiFunction_new::Initialize()
{
    LArRindex_new::Initialize();
    LArTrans_new::Initialize();

}

double LArChiFunction_new::GetChi2()
{
    double chi2 = 0;
    chi2 += LArRindex_new::GetChi2();
    //LArRindex_new::GetChi2();
    chi2 += LArTrans_new::GetChi2();
    
    // Lagrange Multipliers:
    chi2 += m_factor1 * m_ratio;
    double rindex = LArRindex_new::CalcRindex(0.128);
    LArTrans_new::setrindex(rindex);
    chi2 += m_factor2 * LArTrans_new::CalcRayLength(0.128);
    chi2 += m_factor3 * rindex;

    return chi2;
}

void LArChiFunction_new::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


void LArChiFunction_new::SetParameters(double *xx)
{
    LArRindex_new::settemp(xx[0]);
    LArTrans_new::settemp(xx[0]);

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
    LArRindex_new::seta0(xx[10]) ;
    LArRindex_new::setaUV(xx[11]);
    LArRindex_new::setaIR(xx[12]);
    LArRindex_new::setp0(xx[13]);
    LArRindex_new::setp1(xx[14]);
    LArTrans_new::setp0(xx[15]);
    LArTrans_new::setp1(xx[16]);

    // scale for purified spectra
    LArTrans_new::setscale(xx[9]);
}


double LArChiFunction_new::GetChiSquare(double maxChi2)
{
    LArMinuit = new TMinuit();
    LArMinuit->SetFCN(ChisqFCN);
    LArMinuit->SetPrintLevel(1);

    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    LArMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    LArMinuit->mnparm(iPar, "temp", 86, 0.01, 80, 90, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "delta", 0.2, 0., 1.0, 0.001, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "peakratio", 0.94, 0.9, 1.0, 0.001, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "mu1", 126.51, 0.001, 123, 129, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "sigma1", 1, 0.01, 0, 5, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "A2", 0.4, 0.01, 0, 1, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "mu2", 140.121, 0.001, 137, 143, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "sigma2", 1.537, 0.001, 0, 5, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "nu_f", 0, 0.0001, 0, 1, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "scale", 1, 0.001, 0.9, 1.1, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "vara0", 0.335, 0.0001, 0, 1, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "varaUV", 0.099, 0.0001, 0, 0.1, ierrflag); iPar++;
    LArMinuit->mnparm(iPar, "varaIR", 0.008, 0.0001, 0, 0.1, ierrflag); iPar++;


    // Minimization strategy
    LArMinuit->SetErrorDef(1);
    arglist[0] = 2;
    LArMinuit->mnexcm("SET STR", arglist, 1, ierrflag);

    arglist[0] = 50000; // maxCalls
    arglist[1] = 0.01; // tolerance
    LArMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    LArMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    for (Int_t i=0; i<m_nParameters; i++) {
        LArMinuit->GetParameter(i, m_bestFit[i], m_fitError[i]);
    }

    m_chi2Min = min ;
    //m_chi2Min = min;

    //LArMinuit->mnexcm("SHOw COVariance", arglist, 0, ierrflag);

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;

    cout << "=======> 128nm refractive index = " << LArRindex_new::CalcRindex(0.128) << endl;
    //cout << "=======> 128nm refractive index = " << LArRindex_new::CalcRindex(0.128)*(1+LArRindex_new::getnulambda()) << endl;
    cout << "=======> Absorption ratio = " << m_bestFit[4] << endl;
    LArTrans_new::setrindex(LArRindex_new::CalcRindex(0.128));
    cout << "=======> 128nm Rayleigh scattering length = " << LArTrans_new::CalcRayLength(0.128) << endl;
    cout << "=======> chi2min = " << m_chi2Min - m_factor1*m_ratio - m_factor2*LArTrans_new::CalcRayLength(0.128) << endl;

    //cout << "scanning-outputs" << m_bestFit[4] << " " << LArTrans_new::CalcRayLength(0.128) << " " << m_chi2Min - m_factor1*m_ratio - m_factor2*LArTrans_new::CalcRayLength(0.128) << endl;
    //cout << "scanning" << LArRindex_new::CalcRindex(0.128) << " " << m_chi2Min-m_factor3*LArRindex_new::CalcRindex(0.128) << endl;

    delete LArMinuit;
    return min;

}

void LArChiFunction_new::Plot()
{
    LArRindex_new::Plot();
    LArTrans_new::Plot();
}


