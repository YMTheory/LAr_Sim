#include "LArTrans.hh"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

double gRayLength(double* x, double* p)
{
    double l = x[0];  // wavelength um
    double rindex = p[0];  // rindex

    double kT = 2.24442E-9;
    double kB = 1.380649E-23;
    double T = 90; // K
    double f = 1e22;
    
    double pi = TMath::Pi();

    double rayL = 1 / (8*TMath::Power(pi, 3)/3/TMath::Power(l, 4)
                * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f );

    return rayL;
}

double gRayLength_delta(double* x, double* p)
{
    double l = x[0];  // wavelength um
    double rindex = p[0];  // rindex
    double delta = p[1];   // depolarization

    double kT = 2.24442E-9;
    double kB = 1.380649E-23;
    double T = 90; // K
    double f = 1e22;
    
    double pi = TMath::Pi();

    double rayL = 1 / (8*TMath::Power(pi, 3)/3/TMath::Power(l, 4)
                * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f * (6+3*delta)/(6-7*delta));

    return rayL;

}


double gAbs(Double_t* x, Double_t* p)
{
    Double_t l = x[0]*1000;
    Double_t A2 = p[0];
    Double_t mu1 = p[1];
    Double_t sigma1 = p[2];
    Double_t mu2 = p[3];
    Double_t sigma2 = p[4];
    Double_t A1 = p[5];

    Double_t A_abs = A1*TMath::Exp(-(l-mu1)*(l-mu1)/2/sigma1/sigma1) + A2*TMath::Exp(-(l-mu2)*(l-mu2)/2/sigma2/sigma2);
    Double_t T_abs = TMath::Exp( -A_abs*TMath::Log(10.) );
    return T_abs;
}


double gCorr(Double_t* x, Double_t* p)
{
    Double_t wl = x[0]*1000.;
    Double_t E = 1240./wl;

    Double_t a = 38.;
    Double_t E0 = 12.8234;
    Double_t gamma = 0.42357;
    Double_t b = 247.188;
    Double_t E1 = 18.8448;

    Double_t n_MgF2 = TMath::Sqrt(1+a*(E0*E0-E*E)/((E0*E0-E*E)*(E0*E0-E*E)+gamma*E*E) + b/(E1*E1-E*E));

    Double_t n_vac = 1.;
    Double_t n_LAr = p[0];

    return 1/TMath::Power( (1-TMath::Power((n_vac-n_MgF2)/(n_vac+n_MgF2),2))/(1-TMath::Power((n_MgF2-n_LAr)/(n_MgF2+n_LAr) ,2))  ,2);
}


LArTrans::LArTrans(LArRindex* rdx)
{
    depolarization = 1; // 0 for no, 1 for yes
    if (depolarization==0)
        fRayLength = new TF1("fRayLength", gRayLength, 0.11, 0.15, 2);
    else if (depolarization==1)
        fRayLength = new TF1("fRayLength", gRayLength_delta, 0.11, 0.15, 2);

    fAbs = new TF1("fAbs", gAbs, 0.11, 0.15, 6);
    fCorr = new TF1("fCorr", gCorr, 0.1, 0.2, 1);

    gData = new TGraphErrors();
    gCalc = new TGraphErrors();

    gRdx = rdx;
}


void LArTrans::LoadData()
{
    ifstream in; in.open("./data/G140ppb.txt");
    string line;
    double m_wl, m_tran, m_tran_err;
    Int_t idx = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> m_wl >> m_tran >> m_tran_err;
        if( m_wl>125 ) {
            gData->SetPoint(idx, m_wl/1000., m_tran);
            gData->SetPointError(idx, 0, m_tran_err);
            idx++;
        }
    }
}

double LArTrans::GetChi2()
{
    double chi2 = 0;

    double* datay = gData->GetY();
    double* datae = gData->GetEY();
    double* predy = gCalc->GetY();
    for (int i=0; i<gData->GetN(); i++) {
        chi2 += (predy[i]-datay[i]) * (predy[i]-datay[i]) / datae[i]/datae[i]; 
    }

    return chi2;
}

void LArTrans::Calculate()
{
    double* datax = gData->GetX();
    for (int i=0; i<gData->GetN(); i++) {
        double wl = gData->GetX()[i];
        double rindex = gRdx->CalcRindex(datax[i]);
        fRayLength->SetParameter(0, rindex);
        double d = 5.8;
        double rayL = fRayLength->Eval(wl);
        double T_Ray = TMath::Exp(-d/rayL);
        fCorr->SetParameter(0, rindex);
        double corr = fCorr->Eval(wl);
        //fAbs->SetParameters(0.3985, 126.5, 1.04, 0.4122, 140.1, 1.566);
        double T_abs = fAbs->Eval(wl);
        double trans_pred = T_Ray * T_abs * corr;
        //trans_pred *= (1+nu_R) * (1+nu_f);   // nuisance parameters
        
        gCalc->SetPoint(i, datax[i], trans_pred);

    }
}

void LArTrans::SetParameters(double* par)
{
    fRayLength->SetParameter(0, par[0]);
    fRayLength->SetParameter(1, par[1]);
    
    fAbs->SetParameter(0, par[2]);
    fAbs->SetParameter(1, par[3]);
    fAbs->SetParameter(2, par[4]);
    fAbs->SetParameter(3, par[5]);
    fAbs->SetParameter(4, par[6]);
    fAbs->SetParameter(5, par[7]);
}

void LArTrans::Plot()
{}




