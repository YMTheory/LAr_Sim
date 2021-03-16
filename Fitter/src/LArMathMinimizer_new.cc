#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2D.h>

#include "LArMathMinimizer_new.hh"

using namespace std;

bool LArMathMinimizer_new::m_fit_purified = LArConfiguration::fit_purified;
double LArMathMinimizer_new::m_chi2Min;
double LArMathMinimizer_new::m_bestFit[20];
int LArMathMinimizer_new::m_npar;


LArMathMinimizer_new::LArMathMinimizer_new()
{}

LArMathMinimizer_new::~LArMathMinimizer_new()
{}

void LArMathMinimizer_new::Initialize()
{
    LArRindex_new::Initialize();
    LArTrans_new::Initialize();

}


double LArMathMinimizer_new::GetChi2(const double* xx)
{
    LArRindex_new::settemp(xx[0]);
    LArTrans_new::settemp(xx[0]);
    //LArRindex::setrhoratio(xx[12]);

    LArTrans_new::setdelta(xx[1]);
    LArTrans_new::setA1(xx[2]);
    LArTrans_new::setmu1(xx[3]);
    LArTrans_new::setsigma1(xx[4]);
    LArTrans_new::setA2(xx[5]);
    LArTrans_new::setmu2(xx[6]);
    LArTrans_new::setsigma2(xx[7]);

    // pull term
    //LArGroupVelocity::setnulambda(xx[10]);
    LArTrans_new::setnuf(xx[8]);
    LArRindex_new::seta0(xx[10]) ;
    LArRindex_new::setaUV(xx[11]);
    LArRindex_new::setaIR(xx[12]);
    LArRindex_new::setp0(xx[13]);
    LArRindex_new::setp1(xx[14]);
    LArTrans_new::setp0(xx[15]);
    LArTrans_new::setp1(xx[16]);


    // scale for purified spectra
    LArTrans_new::setscale(xx[9]);

    double chi2 = 0;
    chi2 += LArRindex_new::GetChi2();

    chi2 += LArTrans_new::GetChi2();
    
    return chi2;
}


int LArMathMinimizer_new::Minimization()
{

    //ROOT::Minuit2::Minuit2Minimizer_new minimum (ROOT::Minuit2::kMigrad);
    ROOT::Math::Minimizer* minimum = 
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(1);

    // function wrapper for Minimizer_new
    ROOT::Math::Functor f(&GetChi2, 17);
    double step[17];
    for (int i=0; i<17; i++) {
        step[i] = 0.001;
    }
    step[13] = 1e-8;
    step[14] = 1e-5;
    step[15] = 1e-12;
    step[16] = 1e-12;

    // start point
    double variable[17];
    variable[0] = 86;
    variable[1] = 0.2;
    variable[2] = 0.937;
    variable[3] = 126.51;
    variable[4] = 1;
    variable[5] = 0.4;
    variable[6] = 140.121;
    variable[7] = 1.537;
    variable[8] = 0;
    variable[9] = 1.00;
    variable[10] = 0.335;
    variable[11] = 0.099;
    variable[12] = 0.008;
    variable[13] = LArRindex_new::getp0_init();
    variable[14] = LArRindex_new::getp1_init();
    variable[15] = LArTrans_new::getp0_init();
    variable[16] = LArTrans_new::getp1_init();

    minimum->SetFunction(f);

    // Set the free variables to be minimized !
    minimum->SetVariable(0, "temperature", variable[0], step[0]);
    minimum->SetVariable(1, "delta", variable[1], step[1]);
    minimum->SetVariable(2, "peakratio", variable[2], step[2]);
    minimum->SetVariable(3, "mu1", variable[3], step[3]);
    minimum->SetVariable(4, "sigma1", variable[4], step[4]);
    minimum->SetVariable(5, "A2", variable[5], step[5]);
    minimum->SetVariable(6, "mu2", variable[6], step[6]);
    minimum->SetVariable(7, "sigma2", variable[7], step[7]);
    minimum->SetVariable(8, "nu_f", variable[8], step[8]);
    minimum->SetVariable(9, "scale", variable[9], step[9]);
    minimum->SetVariable(10, "vara0", variable[10], step[10]);
    minimum->SetVariable(11, "varaUV", variable[11], step[11]);
    minimum->SetVariable(12, "varaIR", variable[12], step[12]);
    minimum->SetVariable(13, "densityp0", variable[13], step[13]);
    minimum->SetVariable(14, "densityp1", variable[14], step[14]);
    minimum->SetVariable(15, "kappap0", variable[15], step[15]);
    minimum->SetVariable(16, "kappap1", variable[16], step[16]);

    if (m_fit_purified) {
        minimum->FixVariable(2);
        minimum->FixVariable(3);
        minimum->FixVariable(4);
        minimum->FixVariable(6);
        minimum->FixVariable(7);
    } else {
        minimum->FixVariable(9);
    }
    
    //minimum->FixVariable(14);
    //minimum->FixVariable(15);

    // do the minimization
    minimum->Minimize();

    m_chi2Min = minimum->MinValue();

    m_npar = minimum->NDim(); 
    const double *xx = minimum->X();
    for (int i=0; i<m_npar; i++) {
        m_bestFit[i] = xx[i];
    }


    /*
    // Cov Matrix
    cout << "\n";
    cout << "Cov Matrix : " << std::endl; 
    double hess[11][11];
    for (int i=0; i<11; i++) {
        for (int j=0; j<11; j++) {
            hess[i][j] = minimum.CovMatrix(i, j) ;
            cout << hess[i][j] << " ";
        }
        cout << "\n" ;
    }
    */

    //cout << "\n";
    //cout << "Hessian Matrix : " << std::endl; 
    //const int num = 11*11;
    //double  h[num];
    //minimum.GetHessianMatrix(h);
    //for (int i=0; i<num; i++) {
    //    cout << h[i] << " ";
    //    if ( (i+1)%11 == 0)
    //        cout << "\n";
    //}
    
    //const double *xs = minimum.X();
    //std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
    //          << minimum.MinValue() << std::endl;

    //// expected minimum is 0
    //if (minimum.MinValue() < 1.E-4 && f(xs) < 1.E-4)
    //    std::cout << "Minimizer_new " << minName << " - " << algoName
    //              << "   converged to the right minimum" << std::endl;
    //else
    //{
    //    std::cout << "Minimizer_new " << minName << " - " << algoName
    //              << "   failed to converge !!!" << std::endl;
    //    Error("NumericalMinimization", "fail to converge");
    //}

   return 0;
}


bool LArMathMinimizer_new::Plot()
{
    LArRindex_new::Plot();
    LArTrans_new::Plot();
}

void LArMathMinimizer_new::Profile1D(int index, double min, double max, double step, double CI)
{
    TGraph* graph = new TGraph(); graph->SetName("graph");
    double pars[20];
    for(int i=0; i<m_npar; i++) pars[i] = m_bestFit[i];
    double left_margin, right_margin;
    double delta_left = 100; double delta_right = 100;
    double maxChi2 = 0;
    double bestpoint = m_bestFit[index];
    double maxing = 50;

    for(int i=0; i<(max-min)/step; i++) {
        double value = min+i*step;
        pars[index] = value;

        double tmp_chi2 = GetChi2(pars)-m_chi2Min;   // subtract minimum chi2;
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

    cc->SaveAs("profile1d.root");

    return;

}

void LArMathMinimizer_new::Profile2D(int *index, double* min, double* max, int *num)
{
    double pars[20];
    for(int i=0; i<m_npar; i++) pars[i] = m_bestFit[i];
    int i0 = index[0]; int i1 = index[1];

    TGraph* g1 = new TGraph(); g1->SetName("sigma1"); int id1 = 0;
    TGraph* g5 = new TGraph(); g5->SetName("sigma5"); int id5 = 0;

    int ibin = num[0]; int jbin = num[1];
    for(int i=0; i<ibin; i++) {
        for(int j=0; j<jbin; j++) {
            // processing output
            if ( j==0 ) cout << "processing " << ( i/float(ibin)) * 100  << "%..." << endl;

            double p0 = min[0] + (max[0]-min[0]) / ibin * i;
            double p1 = min[1] + (max[1]-min[1]) / jbin * j;
            pars[i0] = p0;
            pars[i1] = p1;

            double tmp_chi2 = GetChi2(pars) - m_chi2Min;

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
    
    cc->SaveAs("Profile2D.root");
}

