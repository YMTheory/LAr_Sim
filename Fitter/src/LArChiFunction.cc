#include "LArChiFunction.hh"
#include "LArConfiguration.hh"

#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2D.h>

using namespace std;

int LArChiFunction::usePull = LArConfiguration::use_pullterm;

double LArChiFunction::m_chi2;
double LArChiFunction::m_chi2Min;
bool LArChiFunction::m_DoFit;
int LArChiFunction::m_nParameters = 13;
double LArChiFunction::m_bestFit[20];
double LArChiFunction::m_fitError[20];
double LArChiFunction::m_ratio;

bool LArChiFunction::m_fit_purified = LArConfiguration::fit_purified;

// Lagrange Multipliers
double LArChiFunction::m_factor1 = LArConfiguration::factor1;
double LArChiFunction::m_factor2 = LArConfiguration::factor2;
double LArChiFunction::m_factor3 = LArConfiguration::factor3;

LArChiFunction::LArChiFunction()
{}

LArChiFunction::~LArChiFunction()
{}

void LArChiFunction::Initialize()
{
    LArRindex::Initialize();
    LArTrans::Initialize();
    LArGroupVelocity::Initialize();

}

double LArChiFunction::GetChi2()
{
    double chi2 = 0;
    chi2 += LArRindex::GetChi2();
    //LArRindex::GetChi2();
    chi2 += LArTrans::GetChi2();
    
    chi2 += LArGroupVelocity::GetChi2();

    // Lagrange Multipliers:
    chi2 += m_factor1 * m_ratio;
    double rindex = LArRindex::CalcRindex(0.128);
    LArTrans::setrindex(rindex);
    chi2 += m_factor2 * LArTrans::CalcRayLength(0.128);
    chi2 += m_factor3 * rindex;

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
    LArRindex::setrhoratio(par[12]);

    LArTrans::setdelta(par[3]);
    LArTrans::setA1(par[4]);
    LArTrans::setmu1(par[5]);
    LArTrans::setsigma1(par[6]);
    LArTrans::setA2(par[7]);
    LArTrans::setmu2(par[8]);
    LArTrans::setsigma2(par[9]);

    // pull term
    //LArRindex::setnulambda(par[10]);
    LArTrans::setnuf(par[10]);

    // scale only for purified spectra:
    LArTrans::setscale(par[11]);
    // Lagrange Multipliers:
    m_ratio = par[4];
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

    if(LArRindex::getoption() ==1 or LArRindex::getoption()==2 ) {
        LArMinuit->mnparm(iPar, "p0", 0.2075, 0.01, 0.1, 0.3, ierrflag);     iPar++;
        LArMinuit->mnparm(iPar, "p1", 0.0415, 0.001, 0.01, 0.1, ierrflag);     iPar++;
        LArMinuit->mnparm(iPar, "p2", 4.333, 0.001, 0., 5., ierrflag);    iPar++;
    } else if (LArRindex::getoption() == 0) {
        LArMinuit->mnparm(iPar, "p0", 0.335, 0.001, 0., 1., ierrflag);     iPar++;
        LArMinuit->mnparm(iPar, "p1", 0.099, 0.001, 0., 1., ierrflag);     iPar++;
        LArMinuit->mnparm(iPar, "p2", 0.008, 0.001, 0., 1., ierrflag);    iPar++;
    }

    LArMinuit->mnparm(iPar, "delta", 0.0, 0.01, 0., 1., ierrflag);   iPar++;
    LArMinuit->mnparm(iPar, "peakratio", 0.947, 0.001, 0.90, 1.0, ierrflag);      iPar++;
    //LArMinuit->mnparm(iPar, "A1", 0.3, 0.01, 0., 0., ierrflag);      iPar++;
    LArMinuit->mnparm(iPar, "mu1",126.51, 0.1, 123, 129, ierrflag);     iPar++;
    LArMinuit->mnparm(iPar, "sigma1", 1, 0.01, 0.5, 2, ierrflag);    iPar++;
    LArMinuit->mnparm(iPar, "A2", 0.4, 0.01, 0., 1., ierrflag);      iPar++;
    LArMinuit->mnparm(iPar, "mu2",140.121, 0.1, 135, 145, ierrflag);     iPar++;
    LArMinuit->mnparm(iPar, "sigma2", 1.537, 0.01, 0.5, 2, ierrflag);  iPar++;
    //LArMinuit->mnparm(iPar, "nu_lambda", 0, 0.01, -1, 1, ierrflag);   iPar++;
    LArMinuit->mnparm(iPar, "nu_f", 0, 0.01, -1, 1, ierrflag);   iPar++;
    //LArMinuit->mnparm(iPar, "rhoratio", 34.49/(44.66e-3), 1, 200, 1200, ierrflag); iPar++;
    
    LArMinuit->mnparm(iPar, "scale", 1, 0.0001, 0.9, 1.1, ierrflag); iPar++;


    //LArMinuit->mnparm(iPar, "delta", 0.2, 0.01, 0., 1., ierrflag);   iPar++;
    //LArMinuit->mnparm(iPar, "peakratio", 1.2, 0.001, 0.0, 2.0, ierrflag);      iPar++;
    ////LArMinuit->mnparm(iPar, "A1", 0.3, 0.01, 0., 0., ierrflag);      iPar++;
    //LArMinuit->mnparm(iPar, "mu1",126, 0.1, 123, 129, ierrflag);     iPar++;
    //LArMinuit->mnparm(iPar, "sigma1", 0.5, 0.01, 0.0, 20, ierrflag);    iPar++;
    //LArMinuit->mnparm(iPar, "A2", 0.1, 0.01, 0., 1., ierrflag);      iPar++;
    //LArMinuit->mnparm(iPar, "mu2",140, 0.1, 135, 145, ierrflag);     iPar++;
    //LArMinuit->mnparm(iPar, "sigma2",11.5, 0.01, 0.5, 10, ierrflag);  iPar++;
    //LArMinuit->mnparm(iPar, "nu_lambda", 0, 0.01, 0, 1, ierrflag);   iPar++;
    //LArMinuit->mnparm(iPar, "nu_f", 0, 0.01, 0, 1, ierrflag);   iPar++;
    //LArMinuit->mnparm(iPar, "rhoratio", 34.49/(44.66e-3), 1, 200, 1200, ierrflag); iPar++;

    //LArMinuit->FixParameter(0);
    //LArMinuit->FixParameter(1);
    //LArMinuit->FixParameter(2);

    if(LArRindex::getoption() == 2) {
        LArMinuit->FixParameter(0);
        LArMinuit->FixParameter(1);
        LArMinuit->FixParameter(2);
    }

    if(LArTrans::getdepolarization() == 0)
        LArMinuit->FixParameter(3);

    if(LArTrans::getfixratio() == 1)
        LArMinuit->FixParameter(4);  // fixed peak ratio;

    if (usePull == 0) {
        LArMinuit->FixParameter(10);
    }

    if (m_fit_purified) {
        LArMinuit->FixParameter(5);
        LArMinuit->FixParameter(6);
        LArMinuit->FixParameter(8);
        LArMinuit->FixParameter(9);
    } else {
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

    m_chi2Min = min ;
    //m_chi2Min = min;

    //LArMinuit->mnexcm("SHOw COVariance", arglist, 0, ierrflag);

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;

    cout << "=======> 128nm refractive index = " << LArRindex::CalcRindex(0.128) << endl;
    //cout << "=======> 128nm refractive index = " << LArRindex::CalcRindex(0.128)*(1+LArRindex::getnulambda()) << endl;
    cout << "=======> Absorption ratio = " << m_bestFit[4] << endl;
    LArTrans::setrindex(LArRindex::CalcRindex(0.128));
    cout << "=======> 128nm Rayleigh scattering length = " << LArTrans::CalcRayLength(0.128) << endl;
    cout << "=======> chi2min = " << m_chi2Min - m_factor1*m_ratio - m_factor2*LArTrans::CalcRayLength(0.128) << endl;

    //cout << "scanning-outputs" << m_bestFit[4] << " " << LArTrans::CalcRayLength(0.128) << " " << m_chi2Min - m_factor1*m_ratio - m_factor2*LArTrans::CalcRayLength(0.128) << endl;
    //cout << "scanning" << LArRindex::CalcRindex(0.128) << " " << m_chi2Min-m_factor3*LArRindex::CalcRindex(0.128) << endl;

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

        if (value==0) cout << "if delta=0, chi2 = " << tmp_chi2 << endl;

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
    //int i0 = 4; int i1 = 3;
    int i0 = 5; int i1 = 8;

    //TH2D* hist = new TH2D("hist", "", 100, min[0], max[0], 100, min[1], max[1]); 
    TGraph* g1 = new TGraph(); g1->SetName("sigma1"); int id1 = 0;
    TGraph* g5 = new TGraph(); g5->SetName("sigma5"); int id5 = 0;

    float ibin = 1000;
    float jbin = 1000;
    for(int i=0; i<int(ibin); i++) {
        for(int j=0; j<int(jbin); j++) {

            // processing output
            if ( j==0 ) cout << "processing " << ( i/ibin ) * 100  << "%..." << endl;


            double p0 = min[0] + step[0]*i;
            double p1 = min[1] + step[1]*j;
            pars[i0] = p0; pars[i1] = p1;
            SetParameters(pars);

            //double LRay = LArTrans::CalcRayLength(0.128);
            double tmp_chi2 = GetChi2()-m_chi2Min;

            //cout << "1st peak: " << pars[i0] << " ,2nd peak: " << pars[i1] << " ,chi2=" << tmp_chi2 << endl;
            if (TMath::Abs(tmp_chi2 -2.30) < 0.01) { g1->SetPoint(id1, p0, p1); id1++;}
            if (TMath::Abs(tmp_chi2 -11.83) < 0.05) { g5->SetPoint(id5, p0, p1); id5++;}
        }
    }
    
    TGraph* bestpoint = new TGraph(); bestpoint->SetName("best");
    bestpoint->SetPoint(0, m_bestFit[i0], m_bestFit[i1]);
    bestpoint->SetMarkerStyle(20);
    bestpoint->SetMarkerColor(kBlue+1);

    g1->SetMarkerColor(kPink+2);
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(0.5);

    g5->SetMarkerColor(kGreen+2);
    g5->SetMarkerStyle(20);
    g5->SetMarkerSize(0.5);


    TCanvas* cc = new TCanvas();
    //hist->Draw("COLZ");
    g5->Draw("AP");
    g1->Draw("P SAME");
    bestpoint->Draw("P SAME");
    
    cc->SaveAs("test.root");
}

void LArChiFunction::Lray_profile()
{
    TGraph* graph = new TGraph();

    double par[20];
    for (int i=0; i<20; i++) par[i] = m_bestFit[i];  // best fit values
    int iteration = 20;
    int n =1;
    double step0 = 2*n*m_fitError[0]/iteration;
    double step1 = 2*n*m_fitError[1]/iteration;
    double step2 = 2*n*m_fitError[2]/iteration;
    double step3 = 2*n*m_fitError[3]/iteration;

    int idx = 0;
    for (int i0=0; i0<iteration; i0++){
        for (int i1=0; i1<iteration; i1++){
            for (int i2=0; i2<iteration; i2++){
                for (int i3=0; i3<iteration; i3++){

                    // processing output
                    if ( i1==0 and i2==0 and i3==0 ) cout << "processing " << ( i0/100.) * 100  << "%..." << endl;

                    par[0] = m_bestFit[0] - n*m_fitError[0] + step0*i0; 
                    par[1] = m_bestFit[1] - n*m_fitError[1] + step1*i1; 
                    par[2] = m_bestFit[2] - n*m_fitError[2] + step2*i2; 
                    par[3] = m_bestFit[3] - n*m_fitError[3] + step3*i3; 
                    SetParameters(par);

                    double tmp_chi2 = GetChi2();

                    double rindex = LArRindex::CalcRindex(0.128);
                    LArTrans::setrindex(rindex);
                    
                    //cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << LArTrans::CalcRayLength(0.128) << " " << tmp_chi2 << endl;
                    graph->SetPoint(idx, LArTrans::CalcRayLength(0.128), tmp_chi2 ); idx++;
                }
            }
        }
    }

    graph->SetName("LRay");
    graph->SaveAs("LRay.root");
}



