#include "LArChiFunction.hh"

using namespace std;

LArChiFunction::LArChiFunction(LArRindex* rdx, LArTrans* trans)
{
    gRdx = rdx;
    gTrans = trans;
}

LArChiFunction::~LArChiFunction()
{
    delete gRdx;
    delete gTrans;
}

void LArChiFunction::Initialize()
{
    gRdx->LoadData();
    gTrans->LoadData();
}

double LArChiFunction::GetChi2()
{
    double chi2 = 0;
    chi2 += gRdx->GetChi2();
    chi2 += gTrans->GetChi2();

    return chi2;
}

void LArChiFunction::SetParameters(double *par)
{
    gRdx->SetParameters:
}







