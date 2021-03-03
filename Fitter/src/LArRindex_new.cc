#include "LArRindex_new.hh"

#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"

using namespace std;


double gRindex_new(Double_t* x, Double_t* p) 
{
    double l = x[0];
    double lUV = 0.1066;
    double lIR = 0.9083;

    double a0  = 0.335;
    double aUV = 0.099; 
    double aIR = 0.008;

    double rho = p[0];

    double A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) * rho;
    double n = TMath::Sqrt(1 + 3*A/(3-A));
    return n;
}

double LArRindex_new::m_rho = 1;
TGraphErrors* LArRindex_new::gData = new TGraphErrors();
TGraphErrors* LArRindex_new::gCalc = new TGraphErrors();
TF1* LArRindex_new::fRindex;

LArRindex_new::LArRindex_new()
{;}

LArRindex_new::~LArRindex_new()
{;}

void LArRindex_new::Initialize()
{

    fRindex = new TF1("fRindex", gRindex_new, 0.1, 0.7, 1);
    LoadData();
}

void LArRindex_new::LoadData()
{
    double m_wavelength[9] = {0.3612, 0.3650, 0.4063, 0.4358, 0.4753, 0.5086, 0.5461, 0.5780, 0.6439};  // um
    double m_rindex0[9] = {1.2395, 1.2392, 1.2372, 1.2361, 1.2349, 1.2341, 1.2334, 1.2328, 1.2321}; // 83.81 K
    double m_rindex1[9] = {1.2370, 1.2367, 1.2347, 1.2336, 1.2324, 1.2316, 1.2308, 1.2303, 1.2296}; // 86 K
    double m_rindex2[9] = {1.2349, 1.2346, 1.2326, 1.2315, 1.2303, 1.2295, 1.2287, 1.2282, 1.2274}; // 88 K

    for(int i=0; i<9; i++) {
        double mean = (m_rindex0[i] + m_rindex1[i] + m_rindex2[i]) / 3.;
        double err  = max(abs(m_rindex0[i]-mean), abs(m_rindex2[i]-mean));
        gData->SetPoint(i, m_wavelength[i], mean);
        gData->SetPointError(i, 0, err);
    }
}

double LArRindex_new::GetChi2()
{
    Calculate();

    double chi2 = 0;

    double *datay = gData->GetY();
    double *datae = gData->GetEY();
    double *calcy = gCalc->GetY();

    for(int i=0; i<gData->GetN(); i++) {
        double pred = calcy[i];
        double y    = datay[i];
        double e    = datae[i];
        chi2 += (pred - y) * (pred - y) /e /e;
    }
    
    return chi2;
}


void LArRindex_new::Calculate()
{
    SetParameters();
    double *datax = gData->GetX();
    for(int i=0; i<gData->GetN(); i++) {
        double calc = fRindex->Eval(datax[i]);
        gCalc->SetPoint(i, datax[i], calc);
    }
}


void LArRindex_new::SetParameters()
{
    fRindex->SetParameter(0, m_rho);
}

void LArRindex_new::Plot()
{
    SetParameters();

    gData->SetMarkerStyle(20);
    gData->SetMarkerColor(kBlue+1);
    gData->SetLineColor(kBlue+1);
    gData->SetLineWidth(2);

    TGraphErrors* gDraw = new TGraphErrors();
    const Int_t N = 1000;
    for(int i=0; i<N; i++) {
        Double_t tmpx = 0.11+0.6/1000*i;
        gDraw->SetPoint(i, tmpx, fRindex->Eval(tmpx));     
    }
    gDraw->SetLineColor(kGreen+1);
    gDraw->SetLineWidth(2);

    TCanvas* cc = new TCanvas("cc", "", 800, 600); cc->cd();
    gData->GetXaxis()->SetLimits(0.11, 0.7);
    gData->GetYaxis()->SetRangeUser(1.2, 1.45);
    gData->SetTitle("rindex fitting; wavelength/um; rindex");

    gData->Draw("AP");
    gDraw->Draw("L SAME");

    cc->SaveAs("rindex-new.root");
}
