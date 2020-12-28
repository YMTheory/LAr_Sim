#include <iostream>
#include <fstream>
#include <sstream>

#include "TGraphErrors.h"
#include "TF1.h"
#include "TGraph.h"

using namespace std;

TGraphErrors* LoadData();
Double_t gRindex(Double_t* x, Double_t* p) 
{
    double A = 1.2055e-2*2/3.;
    double rho_ratio = 34.49/(44.66e-3);
    double l1 = 91.012;
    double l2 = 89.892;
    double l3 = 214.02;
    //double c = 4.333;
    double l = x[0];
    double a = p[0];
    double b = p[1];
    double c = p[2];

    return TMath::Sqrt((3/(1-(A*rho_ratio*(a/(l1-1/l/l)+b/(l2-1/l/l)+c/(l3-1/l/l)))))-2);

}

TF1* fRindex = new TF1("fRindex", gRindex, 0.1, 0.7, 3);

Double_t gRindex20(Double_t* x, Double_t* p)
{
    Double_t lUV = 0.1066; //um
    Double_t lIR = 0.9083; //um
    
    Double_t l = x[0];
    Double_t a0 = p[0];
    Double_t aUV = p[1];
    Double_t aIR = p[2];
    
    Double_t A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) ;

    return TMath::Sqrt(1+3*A/(3-A));
}

TF1* fRindex20 = new TF1("fRindex20", gRindex20, 0.1, 0.7, 3);

Double_t gRayLength(Double_t* x, Double_t* p)
{
    Double_t l = x[0];  // wavelength um
    Double_t rindex = p[0];  // rindex
    Double_t delta = p[1];   // depolarization

    Double_t kT = 2.24442E-9;
    Double_t kB = 1.380649E-23;
    Double_t T = 90; // K
    Double_t f = 1e22;
    
    Double_t pi = TMath::Pi();

    Double_t rayL = 1 / (8*TMath::Power(pi, 3)/3/TMath::Power(l, 4)
                * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f );

    //Double_t rayL = 1 / (8*TMath::Power(pi, 3)/3/TMath::Power(l, 4)
    //            * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f * (6+3*delta)/(6-7*delta));

    return rayL;
}

TF1* fRayLength = new TF1("fRayLength", gRayLength, 0.118, 0.150, 2);

Double_t gCorr(Double_t* x, Double_t* p)
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

    return TMath::Power( (1-TMath::Power((n_vac-n_MgF2)/(n_vac+n_MgF2),2))/(1-TMath::Power((n_MgF2-n_LAr)/(n_MgF2+n_LAr) ,2))  ,2);
}

TF1* fCorr = new TF1("fCorr", gCorr, 0.1, 0.2, 1);
TGraphErrors* calcAbs(TGraphErrors* graph);
TGraph* gaussFit(TGraphErrors* graph, Double_t low, Double_t high);

void absorb_fit()
{
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    TGraphErrors* trans_graph = LoadData();
    TGraphErrors* abs_graph = calcAbs(trans_graph);
    
    // Plotting configurations :
    TCanvas* cc = new TCanvas();

    abs_graph->SetTitle(";wavelength/nm;absorbance");
    abs_graph->SetMarkerColor(kBlue+1);
    abs_graph->SetMarkerSize(0.5);
    abs_graph->SetMarkerStyle(20);
    abs_graph->SetLineColor(kBlue+1);
    abs_graph->SetLineWidth(2);

    abs_graph->GetYaxis()->SetRangeUser(-0.02, 0.55);

    // Fitting  configurations :
    TGraph* fit_graph = gaussFit(abs_graph, 125., 150.);

    abs_graph->Draw("AP");
    fit_graph->Draw("L SAME");

    TGraph* zero_graph = new TGraph();
    zero_graph->SetPoint(0, 118, 0);
    zero_graph->SetPoint(1, 150, 0);
    zero_graph->SetLineColor(kViolet+1);
    zero_graph->SetLineWidth(3);
    zero_graph->SetLineStyle(kDashed);
    zero_graph->Draw("L SAME");

    cc->SaveAs("absFitXe.pdf");    
}

TGraphErrors* LoadData()
{
    ifstream in; in.open("./G140ppb.txt");
    string line;
    Double_t m_wl, m_trans, m_err;
    TGraphErrors* graph = new TGraphErrors();
    Int_t idx = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> m_wl >> m_trans >> m_err;
        graph->SetPoint(idx, m_wl, m_trans);
        graph->SetPointError(idx, 0, m_err);
        idx++;
    }

    return graph;
}

TGraphErrors* calcAbs(TGraphErrors* graph)
{
    TGraphErrors* out = new TGraphErrors();

    Int_t N = graph->GetN();
    for(Int_t i=0; i<N; i++) {
        Double_t wl = graph->GetPointX(i)/1000.;
        Double_t rindex = fRindex->Eval(wl);
        Double_t T_att = graph->GetPointY(i);
        fCorr->SetParameter(0, rindex);
        Double_t corr = fCorr->Eval(wl); 
        T_att *= corr;
        Double_t A_att = -TMath::Log(T_att)/TMath::Log(10) ;
        Double_t err = graph->GetEY()[i];
        Double_t A_att_err = 1/TMath::Log(10)/T_att * corr * err;

        Double_t d = 5.8;  // cm
        //fRindex->SetParameters(3.45227e-01, 1.89568e-02, 4.01598e+00);   // Zhou's model
        fRindex->SetParameters(0.2075, 0.0415, 4.333);
        fRayLength->SetParameters(rindex, 0.3);
        Double_t L_Ray = fRayLength->Eval(wl, rindex);
        Double_t T_Ray = TMath::Exp(-d/L_Ray);
        Double_t A_Ray = -TMath::Log(T_Ray) / TMath::Log(10);

        Double_t A_abs = A_att - A_Ray;
        Double_t A_abs_err = A_att_err;

        out->SetPoint(i, wl*1000, A_abs);
        out->SetPointError(i, 0, A_abs_err);
    }

    return out;
}

TGraph* gaussFit(TGraphErrors* graph, Double_t low, Double_t high)
{
    Double_t abs1 = 126.;
    Double_t abs2 = 141.;
    Int_t N = 100;

    TF1* func = new TF1("func", "[0]*TMath::Exp(-(x-[1])*(x-[1])/2/[2]/[2]) + [3]*TMath::Exp(-(x-[4])*(x-[4])/2/[5]/[5])", low, high);
    func->SetParameter(1, abs1);
    func->SetParameter(2, 2);
    func->SetParameter(4, abs2);
    func->SetParameter(5, 4);
    graph->Fit(func, "RE");

    TGraph* out = new TGraph();
    for(Int_t i=0; i<N; i++) {
        Double_t delta = (high-low)/N;
        out->SetPoint(i, low+i*delta, func->Eval(low+i*delta));
    }

    out->SetLineColor(kGreen+1);
    out->SetLineWidth(2);

    return out;
}