#include "LArTrans.hh"
#include "LArConfiguration.hh"

#include <iostream>
#include <fstream>
#include <sstream>

#include "TCanvas.h"

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
    //Double_t A1 = p[0];
    Double_t mu1 = p[1];
    Double_t sigma1 = p[2];
    Double_t A2 = p[3];
    Double_t mu2 = p[4];
    Double_t sigma2 = p[5];
    Double_t A1 = A2 * p[0]; // p0 is the peak amp ratio;

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

double LArTrans::m_delta;
double LArTrans::m_A1;
double LArTrans::m_mu1;
double LArTrans::m_sigma1;
double LArTrans::m_A2;
double LArTrans::m_mu2;
double LArTrans::m_sigma2;
int LArTrans::depolarization = LArConfiguration::use_depolarization;
int LArTrans::fixratio = LArConfiguration::fix_absratio;

double LArTrans::sigma_R = 0.04;
double LArTrans::sigma_f = 1;
double LArTrans::m_nuR = 0;
double LArTrans::m_nuf = 0;

TGraphErrors* LArTrans::gData;
TGraphErrors* LArTrans::gCalc;
TF1* LArTrans::fRayLength;
TF1* LArTrans::fAbs;
TF1* LArTrans::fCorr;


LArTrans::LArTrans()
{}

void LArTrans::Initialize()
{
    if (depolarization==0) {
        fRayLength = new TF1("fRayLength", gRayLength, 0.11, 0.15, 2);
        cout << "+++++++++++++++++++++++ Did not consider depolarization" << endl;
    }
    else if (depolarization==1) {
        fRayLength = new TF1("fRayLength", gRayLength_delta, 0.11, 0.15, 2);
        cout << "+++++++++++++++++++++++ Did consider depolarization" << endl;
    }

    if (fixratio) {
        cout << "+++++++++++++++++++++++ fix absorption ratio" << endl;
    }

    fAbs = new TF1("fAbs", gAbs, 0.11, 0.15, 6);
    fCorr = new TF1("fCorr", gCorr, 0.1, 0.2, 1);

    gData = new TGraphErrors();
    gCalc = new TGraphErrors();

    LoadData();
}

void LArTrans::LoadData()
{
    cout << "+++++++++++++++++++ Load Transmission Data ..." << endl;
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

void LArTrans::SetParameters()
{
    fRayLength->SetParameter(1, m_delta);
    fAbs->SetParameters(m_A1, m_mu1, m_sigma1, m_A2, m_mu2, m_sigma2);
}


double LArTrans::GetChi2()
{
    Calculate();

    double chi2 = 0;

    double* datay = gData->GetY();
    double* datae = gData->GetEY();
    double* predy = gCalc->GetY();
    for (int i=0; i<gData->GetN(); i++) {
        predy[i] *= (1+m_nuf);
        chi2 += (predy[i]-datay[i]) * (predy[i]-datay[i]) / datae[i]/datae[i]; 
    }
    chi2 += (m_nuf)*(m_nuf)/(sigma_f)/(sigma_f);
        
    return chi2;
}

void LArTrans::Calculate()
{
    SetParameters();
    double* datax = gData->GetX();
    for (int i=0; i<gData->GetN(); i++) {
        double wl = gData->GetX()[i];
        double rindex = LArRindex::CalcRindex(datax[i]);
        fRayLength->SetParameter(0, rindex);
        double d = 5.8;
        double rayL = fRayLength->Eval(wl);
        double T_Ray = TMath::Exp(-d/rayL);
        fCorr->SetParameter(0, rindex);
        double corr = fCorr->Eval(wl);
        double T_abs = fAbs->Eval(wl);
        double trans_pred = T_Ray * T_abs * corr;
        //trans_pred *= (1+nu_R) * (1+nu_f);   // nuisance parameters
        
        gCalc->SetPoint(i, datax[i], trans_pred);

    }
}

void LArTrans::Plot()
{
    TGraph* pred_graph = new TGraph();
    for (Int_t i=0; i<1000; i++) {
        Double_t wl = ((150.-125.)/1000*i + 125. )/1000.;
        Double_t rindex = LArRindex::CalcRindex(wl);
        fRayLength->SetParameter(0, rindex);
        fCorr->SetParameter(0, rindex);
        Double_t corr = fCorr->Eval(wl);
        Double_t rayL = fRayLength->Eval(wl);
        Double_t T_abs = fAbs->Eval(wl);
        Double_t T_Ray = TMath::Exp(-5.8/(rayL)) ;
        Double_t trans_pred = T_Ray * T_abs * corr;
        pred_graph->SetPoint(i, wl, trans_pred);

    }
    pred_graph->SetLineColor(kGreen+1);
    pred_graph->SetLineWidth(3);

    gData->SetMarkerColor(kBlue+1);
    gData->SetLineColor(kBlue+1);
    gData->SetMarkerSize(0.5);
    gData->SetMarkerStyle(20);
    gData->SetLineWidth(2);

    TCanvas* cg = new TCanvas("cg", "", 800, 600);
    gData->SetTitle("Transmission; wavelength/um; transmission");
    //gTransData->GetXaxis()->SetLimits(0.11, 0.7);
    //gTransData->GetYaxis()->SetRangeUser(1.2, 1.45);
    gData->Draw("AP");
    pred_graph->Draw("L SAME");
    cg->SaveAs("Trans.pdf");
}




