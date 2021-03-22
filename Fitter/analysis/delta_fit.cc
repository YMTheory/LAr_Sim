#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include <iostream>
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"

using namespace RooFit;

void delta_fit()
{
    
    // Histogram Loading
    
    TH1F* hist = new TH1F("hist", "", 100, 0, 0.5);

    ifstream in; in.open("../delta_toyMC.txt");
    string line; double tmp;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmp;
        hist->Fill(tmp);
    }

    gStyle->SetOptFit(1101);

    RooRealVar x("x", "x", 0, 0.5);
    RooRealVar mean("mean", "mean", 0.32, 0, 0.5);
    RooRealVar sigma("sigma", "sigma", 0.05, 0, 0.2);
    RooRealVar alpha("alpha", "alpha", 80, 0, 400);
    RooRealVar n("n", "n", 80, 0, 200);

    RooCBShape cball("cball", "crystal ball", x, mean, sigma, alpha, n);

    RooDataHist dh("dh", "dh", x, Import(*hist));

    cball.fitTo(dh);

    RooPlot* xframe = x.frame(Title("crystal ball pdf"));
    dh.plotOn(xframe);
    cball.plotOn(xframe);

    TPaveText* pt2 = new TPaveText(100,200, 110,200);
    pt2->SetBorderSize(0);
    pt2->SetFillColor(0);
    pt2->SetTextAlign(12);
    pt2->SetTextSize(0.04);
    pt2->SetTextFont(42);
    TString Par5V2 = Form("%2.1f", xframe->chiSquare());
    TString Par52 = "#chi^{2}#/n.d.f = " + Par5V2;
    TText *text2 = pt2->AddText(Par52);

    TCanvas* cc = new TCanvas("cball fitting", "cball fitting", 600, 600);
    xframe->SetTitle("CB Fitting; delta");
    xframe->Draw();
    //pt2->Draw("SAME");

    
}

