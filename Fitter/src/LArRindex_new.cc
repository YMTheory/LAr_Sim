#include "LArRindex_new.hh"
#include "LArConfiguration.hh"

#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TRandom.h"

using namespace std;


double m_rho90K = 0.03449;
double gRindex_new(Double_t* x, Double_t* p) 
{
    double l = x[0];
    double lUV = 0.1066;
    double lIR = 0.9083;

    double a0  = p[1]; //0.335;
    double aUV = p[2]; //0.099; 
    double aIR = p[3]; //0.008;

    double T = p[0];
    double p0 = p[4];
    double p1 = p[5];

    double rho = ( p0*T + p1 ) / m_rho90K; // density scaling

    //cout << a0 << " " << aUV << " " << aIR << endl;

    double A = ( a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) ) * rho;
    double n = TMath::Sqrt(1 + 3*A/(3-A));
    return n;
}

double LArRindex_new::p0 = -0.000158;
double LArRindex_new::p1 = 0.0487;
double LArRindex_new::sigma_p0 = 0.000004;
double LArRindex_new::sigma_p1 = 0.0003;
double LArRindex_new::m_p0;
double LArRindex_new::m_p1;
double LArRindex_new::m_temp = 86;
double LArRindex_new::m_rho = 1;
double LArRindex_new::m_a0;
double LArRindex_new::m_aUV;
double LArRindex_new::m_aIR;
bool LArRindex_new::m_loadData = false;
bool LArRindex_new::m_toyMC = LArConfiguration::m_toyMC;

TGraphErrors* LArRindex_new::gData = new TGraphErrors();
TGraphErrors* LArRindex_new::gCalc = new TGraphErrors();
TGraphErrors* LArRindex_new::gtoyMC = new TGraphErrors();
TF1* LArRindex_new::fRindex;

LArRindex_new::LArRindex_new()
{;}

LArRindex_new::~LArRindex_new()
{;}

void LArRindex_new::Initialize()
{

    std::cout << "===========> LArRindex_new Initialization" << std::endl;
    fRindex = new TF1("fRindex", gRindex_new, 0.1, 0.7, 6);
    LoadData();
    if (m_toyMC) toyMC();
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

    m_loadData = true;
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

    // pull term due to fitting error :
    chi2 += ((m_a0 - 0.335)/0.003 )  * ((m_a0-0.335)/0.003 ); 
    chi2 += ((m_aUV - 0.099)/0.003 ) * ((m_aUV-0.099)/0.003 ); 
    chi2 += ((m_aIR - 0.008)/0.003 ) * ((m_aIR-0.008)/0.003 ); 

    chi2 += ((m_p0 - p0)/sigma_p0) * ((m_p0 - p0)/sigma_p0);
    chi2 += ((m_p1 - p1)/sigma_p1) * ((m_p1 - p1)/sigma_p1);
    
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
    fRindex->SetParameter(0, m_temp);
    fRindex->SetParameter(1, m_a0);
    fRindex->SetParameter(2, m_aUV);
    fRindex->SetParameter(3, m_aIR);
    fRindex->SetParameter(4, m_p0);
    fRindex->SetParameter(5, m_p1);
}

void LArRindex_new::Plot()
{
    SetParameters();

    std::cout << "Rindex @128nm = " << fRindex->Eval(0.128) << std::endl;

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



void LArRindex_new::toyMC()
{
    if (!m_loadData)
        LoadData();
    
    for(int i=0; i<gData->GetN(); i++ ) {
        double tmp = gRandom->Gaus(gData->GetY()[i], gData->GetEY()[i]);
        gtoyMC->SetPoint(i, gData->GetX()[i], tmp);
        gtoyMC->SetPointError(i, 0, gData->GetEY()[i]);
    }
}
