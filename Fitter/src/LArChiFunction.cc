#include "LArChiFunction.hh"

using namespace std;

double LArChiFunction::m_chi2;
double LArChiFunction::m_chi2Min;
bool LArChiFunction::m_DoFit;

LArChiFunction::LArChiFunction()
{}

LArChiFunction::~LArChiFunction()
{}

void LArChiFunction::Initialize()
{
    LArRindex::Initialize();
    LArTrans::Initialize();
}

double LArChiFunction::GetChi2()
{
    double chi2 = 0;
    chi2 += LArRindex::GetChi2();
    chi2 += LArTrans::GetChi2();

    return chi2;
}

void LArChiFunction::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


void LArChiFunction::SetParameters(double *par)
{
    LArRindex::setp0(par[0]);
    LArRindex::setp1(par[1]);
    LArRindex::setp2(par[2]);

    LArTrans::setdelta(par[3]);
    LArTrans::setA1(par[4]);
    LArTrans::setmu1(par[5]);
    LArTrans::setsigma1(par[6]);
    LArTrans::setA2(par[7]);
    LArTrans::setmu2(par[8]);
    LArTrans::setsigma2(par[9]);

    // pull term
    LArRindex::setnulambda(par[10]);
    LArTrans::setnuf(par[11]);
}


double LArChiFunction::GetChiSquare(double maxChi2)
{
    LArMinuit = new TMinuit();
    LArMinuit->SetFCN(ChisqFCN);
    LArMinuit->SetPrintLevel(1);

    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    LArMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    LArMinuit->mnparm(iPar, "p0", 0.335, 0.01, 0., 1., ierrflag);     iPar++;
    LArMinuit->mnparm(iPar, "p1", 0.099, 0.001, 0., 1., ierrflag);     iPar++;
    LArMinuit->mnparm(iPar, "p2", 0.008, 0.001, 0., 1., ierrflag);    iPar++;
    LArMinuit->mnparm(iPar, "delta", 0.0, 0.01, 0., 1., ierrflag);   iPar++;
    LArMinuit->mnparm(iPar, "A1", 0.3, 0.01, 0., 0., ierrflag);      iPar++;
    LArMinuit->mnparm(iPar, "mu1",126, 0.1, 123, 129, ierrflag);     iPar++;
    LArMinuit->mnparm(iPar, "sigma1", 1, 0.01, 0.5, 2, ierrflag);    iPar++;
    LArMinuit->mnparm(iPar, "A2", 0.3, 0.01, 0., 0., ierrflag);      iPar++;
    LArMinuit->mnparm(iPar, "mu2",140, 0.1, 135, 145, ierrflag);     iPar++;
    LArMinuit->mnparm(iPar, "sigma2", 1.5, 0.01, 0.5, 2, ierrflag);  iPar++;
    LArMinuit->mnparm(iPar, "nu_lambda", 0, 0.01, 0, 1, ierrflag);   iPar++;
    LArMinuit->mnparm(iPar, "nu_f", 0, 0.01, 0, 1, ierrflag);   iPar++;

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

    //LArMinuit->mnexcm("SHOw COVariance", arglist, 0, ierrflag);

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;

    delete LArMinuit;
    return min;

}

void LArChiFunction::Plot()
{
    LArRindex::Plot();
    LArTrans::Plot();
}





