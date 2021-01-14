#include "LArChiFunction.hh"
#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2D.h>

using namespace std;

int LArChiFunction::usePull = 1;

double LArChiFunction::m_chi2;
double LArChiFunction::m_chi2Min;
bool LArChiFunction::m_DoFit;
int LArChiFunction::m_nParameters = 12;
double LArChiFunction::m_bestFit[20];
double LArChiFunction::m_fitError[20];

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

    if(LArRindex::getoption() ==1) {
        LArMinuit->mnparm(iPar, "p0", 0.207, 0.01, 0., 1., ierrflag);     iPar++;
        LArMinuit->mnparm(iPar, "p1", 0.0415, 0.001, 0., 1., ierrflag);     iPar++;
        LArMinuit->mnparm(iPar, "p2", 4.333, 0.001, 0., 5., ierrflag);    iPar++;
    } else if (LArRindex::getoption() == 0) {
        LArMinuit->mnparm(iPar, "p0", 0.335, 0.001, 0., 1., ierrflag);     iPar++;
        LArMinuit->mnparm(iPar, "p1", 0.099, 0.001, 0., 1., ierrflag);     iPar++;
        LArMinuit->mnparm(iPar, "p2", 0.008, 0.001, 0., 1., ierrflag);    iPar++;
    }
    LArMinuit->mnparm(iPar, "delta", 0.0, 0.01, 0., 1., ierrflag);   iPar++;
    LArMinuit->mnparm(iPar, "ratio", 0.947, 0.001, 0.90, 1.0, ierrflag);      iPar++;
    //LArMinuit->mnparm(iPar, "A1", 0.3, 0.01, 0., 0., ierrflag);      iPar++;
    LArMinuit->mnparm(iPar, "mu1",126, 0.1, 123, 129, ierrflag);     iPar++;
    LArMinuit->mnparm(iPar, "sigma1", 1, 0.01, 0.5, 2, ierrflag);    iPar++;
    LArMinuit->mnparm(iPar, "A2", 0.4, 0.01, 0., 1., ierrflag);      iPar++;
    LArMinuit->mnparm(iPar, "mu2",140, 0.1, 135, 145, ierrflag);     iPar++;
    LArMinuit->mnparm(iPar, "sigma2", 1.5, 0.01, 0.5, 2, ierrflag);  iPar++;
    LArMinuit->mnparm(iPar, "nu_lambda", 0, 0.01, 0, 1, ierrflag);   iPar++;
    LArMinuit->mnparm(iPar, "nu_f", 0, 0.01, 0, 1, ierrflag);   iPar++;
    
    if(LArTrans::getdepolarization() == 0)
        LArMinuit->FixParameter(3);

    if(LArTrans::getfixratio() == 1)
        LArMinuit->FixParameter(4);  // fixed peak ratio;

    if (usePull == 0) {
        LArMinuit->FixParameter(10);
        LArMinuit->FixParameter(11);
    }

    // use values from Zhou's note
    //LArMinuit->FixParameter(0);
    //LArMinuit->FixParameter(1);
    //LArMinuit->FixParameter(2);

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

    m_chi2Min = min;

    //LArMinuit->mnexcm("SHOw COVariance", arglist, 0, ierrflag);

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;

    cout << "=======> 128nm refractive index = " << LArRindex::CalcRindex(0.128) << endl;
    LArTrans::setrindex(LArRindex::CalcRindex(0.128));
    cout << "=======> 128nm Rayleigh scattering length = " << LArTrans::CalcRayLength(0.128) << endl;

    delete LArMinuit;
    return min;

}

void LArChiFunction::Plot()
{
    LArRindex::Plot();
    LArTrans::Plot();
}

void LArChiFunction::Profile1D(int index, double min, double max, double step, double CI)
{
    TGraph* graph = new TGraph(); graph->SetName("graph");
    double pars[20];
    for(int i=0; i<m_nParameters; i++) pars[i] = m_bestFit[i];
    double left_margin, right_margin;
    double delta_left = 100; double delta_right = 100;
    double maxChi2 = 0;
    double bestpoint = m_bestFit[index];
    double maxing = 50;

    for(int i=0; i<(max-min)/step; i++) {
        double value = min+i*step;
        pars[index] = value;
        SetParameters(pars);

        double tmp_chi2 = GetChi2()-m_chi2Min;   // subtract minimum chi2;
        graph->SetPoint(i, value, tmp_chi2);

        if(value<bestpoint and TMath::Abs(tmp_chi2-CI*CI)<delta_left ) {
            delta_left = tmp_chi2-CI*CI;
            left_margin = value;
        }
        else if(value>bestpoint and TMath::Abs(tmp_chi2-CI*CI)<delta_right ) {
            delta_right = TMath::Abs(tmp_chi2-CI*CI);
            right_margin = value;
        }

        if(tmp_chi2>maxChi2) maxChi2 = tmp_chi2;
    }
    cout << "68% CI is " << "[ " << left_margin << ", " << right_margin << " ], with best fit value = " <<  bestpoint << endl;

    TGraph* zone = new TGraph(); zone->SetName("zone");
    zone->SetPoint(0, left_margin, 0);
    zone->SetPoint(1, left_margin, maxChi2);
    zone->SetPoint(2, right_margin, maxChi2);
    zone->SetPoint(3, right_margin, 0);

    TGraph* line = new TGraph(); line->SetName("line");
    line->SetPoint(0, bestpoint, 0);
    line->SetPoint(1, bestpoint, maxChi2);

    TCanvas* cc = new TCanvas();
    graph->GetYaxis()->SetRangeUser(0, maxChi2);
    graph->SetLineColor(kBlue+1);
    graph->Draw("AL");
    zone->SetFillStyle(3013);
    zone->SetFillColor(29);
    zone->Draw("F SAME");
    line->SetLineColor(kBlue+2);
    line->SetLineWidth(3);
    line->SetLineStyle(kDashed);
    line->Draw("L SAME");

    cc->SaveAs("delta.root");

    return;
}


void LArChiFunction::Profile2D(double* min, double* max, double* step)
{
    double pars[20];
    for(int i=0; i<m_nParameters; i++) pars[i] = m_bestFit[i];
    int i0 = 4; int i1 = 3;

    //TH2D* hist = new TH2D("hist", "", 100, min[0], max[0], 100, min[1], max[1]); 
    TGraph* g1 = new TGraph(); g1->SetName("sigma1"); int id1 = 0;
    TGraph* g5 = new TGraph(); g5->SetName("sigma5"); int id5 = 0;

    float ibin = 1000;
    float jbin = 1000;
    for(int i=0; i<int(ibin); i++) {
        for(int j=0; j<int(jbin); j++) {

            // processing output
            if ( j==0 ) cout << "processing " << ( i/ibin ) * 100  << "%..." << endl;


            double Xe_ratio = min[0] + step[0]*i;
            double delta = min[1] + step[1]*j;
            pars[i0] = Xe_ratio; pars[i1] = delta;
            SetParameters(pars);

            //double LRay = LArTrans::CalcRayLength(0.128);
            double tmp_chi2 = GetChi2()-m_chi2Min;

            //int fillx = int( (Xe_ratio - min[0]) /step[0] )+1;
            //int filly = int( (delta-min[1]) / step[1] ) +1;

            //hist->SetBinContent(fillx, filly, tmp_chi2);
            if (TMath::Abs(tmp_chi2 -1) < 0.01) { g1->SetPoint(id1, Xe_ratio, delta); id1++;}
            if (TMath::Abs(tmp_chi2 -25) < 0.05) { g5->SetPoint(id5, Xe_ratio, delta); id5++;}
        }
    }
    
    TGraph* bestpoint = new TGraph(); bestpoint->SetName("best");
    bestpoint->SetPoint(0, m_bestFit[i0], m_bestFit[i1]);
    bestpoint->SetMarkerStyle(20);
    bestpoint->SetMarkerColor(kBlue+1);

    g1->SetMarkerColor(kPink+2);
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(0.5);

    g5->SetMarkerColor(kViolet+3);
    g5->SetMarkerStyle(20);
    g5->SetMarkerSize(0.5);


    TCanvas* cc = new TCanvas();
    //hist->Draw("COLZ");
    g5->Draw("AP");
    g1->Draw("P SAME");
    bestpoint->Draw("P SAME");
    
    cc->SaveAs("test.root");
}



