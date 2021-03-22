#include "LArTrans_new.hh"
#include "LArConfiguration.hh"

#include <iostream>
#include <fstream>
#include <sstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraphErrors.h"

using namespace std;

double gRayLength_new(double* x, double* p)
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

double gRayLength_delta_new(double* x, double* p)
{
    double l = x[0];  // wavelength um
    double rindex = p[0];  // rindex
    double delta = p[1];   // depolarization
    double T = p[2];
    
    double n_p0 = p[3];
    double n_p1 = p[4];
    double kT = n_p0 * T + n_p1;

    //double kT = 2.25282E-9;
    //double kT = 2.24442E-9;
    double kB = 1.380649E-23;
    //double T = 90; // K
    double f = 1e22;
    
    double pi = TMath::Pi();

    double rayL = 1 / (8*TMath::Power(pi, 3)/3/TMath::Power(l, 4)
                * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f * (6+3*delta)/(6-7*delta));

    //cout << rindex << " " << delta << " " << rayL << endl;

    return rayL;

}


double gAbs_new(Double_t* x, Double_t* p)
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


double gCorr_new(Double_t* x, Double_t* p)
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

double LArTrans_new::m_delta;
double LArTrans_new::m_A1;
double LArTrans_new::m_mu1;
double LArTrans_new::m_sigma1;
double LArTrans_new::m_A2;
double LArTrans_new::m_mu2;
double LArTrans_new::m_sigma2;
double LArTrans_new::m_scale;
double LArTrans_new::m_kappaT;
double LArTrans_new::m_meankappaT;
double LArTrans_new::m_sigmakappaT;
double LArTrans_new::m_temp;
double LArTrans_new::temp = 86;
double LArTrans_new::sigma_temp = 3./TMath::Sqrt(12);
double LArTrans_new::m_p0;
double LArTrans_new::m_p1;
bool LArTrans_new::m_loadData = false;
bool LArTrans_new::m_toyMC = LArConfiguration::m_toyMC;
// from fitting 
double LArTrans_new::p0 =  6.0712e-11;
double LArTrans_new::p1 = -3.1699e-09;
double LArTrans_new::sigmap0 = TMath::Sqrt(1.49189964e-24);
double LArTrans_new::sigmap1 = TMath::Sqrt(1.11432969e-20);

int LArTrans_new::depolarization = LArConfiguration::use_depolarization;
int LArTrans_new::fixratio = LArConfiguration::fix_absratio;
bool LArTrans_new::m_fit_purified = LArConfiguration::fit_purified;

double LArTrans_new::sigma_R = 0.04;
double LArTrans_new::sigma_f = 1;
double LArTrans_new::m_nuR = 0;
double LArTrans_new::m_nuf = 0;

TGraphErrors* LArTrans_new::gData;
TGraphErrors* LArTrans_new::gCalc;
TGraphErrors* LArTrans_new::gtoyMC;
TF1* LArTrans_new::fRayLength;
TF1* LArTrans_new::fAbs;
TF1* LArTrans_new::fCorr;


LArTrans_new::LArTrans_new()
{}

void LArTrans_new::Initialize()
{
    if (depolarization==0) {
        fRayLength = new TF1("fRayLength", gRayLength_new, 0.11, 0.15, 2);
        cout << "+++++++++++++++++++++++ Did not consider depolarization" << endl;
    }
    else if (depolarization==1) {
        fRayLength = new TF1("fRayLength", gRayLength_delta_new, 0.11, 0.15, 5);
        cout << "+++++++++++++++++++++++ Did consider depolarization" << endl;
    }

    if (fixratio) {
        cout << "+++++++++++++++++++++++ fix absorption ratio" << endl;
    }

    if (m_fit_purified) {
        cout << "+++++++++++++++++++++++ fit purified LAr data" << endl;
    } else {
        cout << "+++++++++++++++++++++++ fit Xe-doped LAr data" << endl;
    }

    fAbs = new TF1("fAbs", gAbs_new, 0.11, 0.15, 6);
    fCorr = new TF1("fCorr", gCorr_new, 0.1, 0.2, 1);

    gData = new TGraphErrors();
    gCalc = new TGraphErrors();

    LoadData();
    
    if(m_toyMC) toyMC();
}

void LArTrans_new::LoadData()
{
    cout << "+++++++++++++++++++ Load Transmission Data ..." << endl;
    ifstream in; 
    if (m_fit_purified)
        in.open("./data/data2012.txt");
    else
        in.open("./data/G140ppb.txt");
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

    double kappaT0 = 1.9355e-9; // 84K
    double kappaT1 = 1.9903e-9; // 85K
    double kappaT2 = 2.0473e-9; // 86K
    double kappaT3 = 2.1064e-9; // 87K

    m_meankappaT = (kappaT0 + kappaT1 + kappaT2 + kappaT3) / 4.; // consider kappaT values at 4 conditions ...
    m_sigmakappaT = max(kappaT3-m_meankappaT, m_meankappaT-kappaT0);

    m_loadData = true;
}

void LArTrans_new::SetParameters()
{
    fRayLength->SetParameter(1, m_delta);
    fRayLength->SetParameter(2, m_temp);
    fRayLength->SetParameter(3, m_p0);
    fRayLength->SetParameter(4, m_p1); 

    fAbs->SetParameters(m_A1, m_mu1, m_sigma1, m_A2, m_mu2, m_sigma2);
}


double LArTrans_new::GetChi2()
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

    // add pull term
    chi2 += (m_nuf)*(m_nuf)/(sigma_f)/(sigma_f);
    //chi2 += ((m_kappaT - m_meankappaT) / m_sigmakappaT ) * ((m_kappaT - m_meankappaT) / m_sigmakappaT);
    chi2 += TMath::Power((m_p0-p0)/sigmap0, 2);  
    chi2 += TMath::Power((m_p1-p1)/sigmap1, 2);  

    // temperature pull term
    chi2 += TMath::Power((m_temp - temp)/sigma_temp, 2);

    return chi2;
}

void LArTrans_new::Calculate()
{
    SetParameters();
    double* datax = gData->GetX();
    for (int i=0; i<gData->GetN(); i++) {
        double wl = gData->GetX()[i];
        double rindex = LArRindex_new::CalcRindex(datax[i]);
        fRayLength->SetParameter(0, rindex);
        double d = 5.8;
        double rayL = fRayLength->Eval(wl);
        double T_Ray = TMath::Exp(-d/rayL);
        fCorr->SetParameter(0, rindex);
        double corr = fCorr->Eval(wl);
        double T_abs = fAbs->Eval(wl);
        double trans_pred = T_Ray * T_abs * corr * m_scale;
        //trans_pred *= (1+nu_R) * (1+nu_f);   // nuisance parameters
        
        //cout << wl << " " << gData->GetY()[i] << " " << trans_pred << endl;
        gCalc->SetPoint(i, datax[i], trans_pred);

    }
}

void LArTrans_new::Plot()
{
    TGraph* pred_graph = new TGraph();
    TGraph* abs_graph  = new TGraph();
    TGraph* ray_graph  = new TGraph();
    TGraphErrors* data_graph = new TGraphErrors();
    for(Int_t i=0; i<gData->GetN(); i++) {
        Double_t corr = fCorr->Eval(gData->GetPointX(i));
        data_graph->SetPoint(i, gData->GetPointX(i), gData->GetPointY(i)/corr);
        data_graph->SetPointError(i, 0, gData->GetEY()[i]);
    }
    for (Int_t i=0; i<1000; i++) {
        Double_t wl = ((150.-125.)/1000*i + 125. )/1000.;
        Double_t rindex = LArRindex_new::CalcRindex(wl);
        fRayLength->SetParameter(0, rindex);
        fCorr->SetParameter(0, rindex);
        Double_t corr = fCorr->Eval(wl);
        Double_t rayL = fRayLength->Eval(wl);
        Double_t T_abs = fAbs->Eval(wl);
        Double_t T_Ray = TMath::Exp(-5.8/(rayL)) ;
        Double_t trans_pred = T_Ray * T_abs;
        pred_graph->SetPoint(i, wl, trans_pred);
        abs_graph->SetPoint(i, wl, T_abs);
        ray_graph->SetPoint(i, wl, T_Ray);
    }
    pred_graph->SetLineColor(kGreen+1);
    pred_graph->SetLineWidth(3);
    abs_graph->SetLineColor(kOrange+1);
    abs_graph->SetLineWidth(3);
    ray_graph->SetLineColor(kViolet+1);
    ray_graph->SetLineWidth(3);

    data_graph->SetMarkerColor(kBlue+1);
    data_graph->SetLineColor(kBlue+1);
    data_graph->SetMarkerSize(0.5);
    data_graph->SetMarkerStyle(20);
    data_graph->SetLineWidth(2);

    TCanvas* cg = new TCanvas("cg", "", 800, 600);
    data_graph->SetTitle("Transmission; wavelength/um; transmission");
    //gTransData->GetXaxis()->SetLimits(0.11, 0.7);
    //gTransData->GetYaxis()->SetRangeUser(1.2, 1.45);
    data_graph->Draw("AP");
    pred_graph->Draw("L SAME");
    ray_graph->Draw("L SAME");
    abs_graph->Draw("L SAME");

    TLegend* ll = new TLegend();
    ll->AddEntry(pred_graph, "Transmission");
    ll->AddEntry(abs_graph, "absorption");
    ll->AddEntry(ray_graph, "rayleigh");
    ll->Draw("SAME");

    cg->SaveAs("Trans-new.root");
}



double LArTrans_new::CalcRayLength(double wl)
{
    double rindex = LArRindex_new::CalcRindex(wl);
    fRayLength->SetParameter(0, rindex);
    return fRayLength->Eval(wl);
}


void LArTrans_new::toyMC()
{
    if(!m_loadData)
        LoadData();

    for(int i=0; i<gData->GetN(); i++ )
    {
        double tmp = gRandom->Gaus(gData->GetY()[i], gData->GetEY()[i]);
        gtoyMC->SetPoint(i, gData->GetX()[i], tmp);
        gtoyMC->SetPointError(i, 0, gData->GetEY()[i]);
    }
}
