#include "LArRindex.hh"

#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"

using namespace std;


double gRindex(Double_t* x, Double_t* p) 
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


double gRindex20(double *x, double *p) {
    double l = x[0];
    double lUV = 0.1066;
    double lIR = 0.9083;

    double a0 = p[0];
    double aUV = p[1];
    double aIR = p[2];

    double A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR);
    double n = TMath::Sqrt(1 + 3*A/(3-A));

    return n;

}

LArRindex::LArRindex()
{
    option = 0;
    if (option==1)
        fRindex = new TF1("fRindex", gRindex, 0.1, 0.7, 3);
    else if (option==0)
        fRindex = new TF1("fRindex", gRindex20, 0.1, 0.7, 3);
    gData = new TGraphErrors();
    gCalc = new TGraphErrors();
    gData128nm = new TGraphErrors();
}

double LArRindex::GetChi2()
{
    double chi2 = 0;

    double* datay = gData->GetY();
    double* datae = gData->GetEY();
    double* calce = gCalc->GetY();

    for (int i=0; i<gData->GetN(); i++ )
    {
        double pred = calce[i];
        double y = datay[i];
        double e = datae[i];
        chi2 += (pred-y) * (pred-y) / e / e;
    }

    //128nm group velocity
    if (option == 0) {

        // delta chi2 at 128 nm for Babicz
        Double_t a0  = fRindex->GetParameter(0);
        Double_t aUV = fRindex->GetParameter(1);
        Double_t aIR = fRindex->GetParameter(2);
        Double_t part1 = -3*(0.323001*aIR + 115.417*aUV)*(a0 - 0.0202616*aIR + 3.26346*aUV)/(3 -a0+ 0.0202616*aIR - 3.26346*aUV)/(3 -a0+ 0.0202616*aIR - 3.26346*aUV);
        Double_t part2 = 3*(-0.323001*aIR - 115.417*aUV)/(3-a0 +0.0202616*aIR-3.26346*aUV);
        Double_t part3 = 2*TMath::Sqrt(1+(3*(a0-0.0202616*aIR+3.26346*aUV))/(3-a0+0.0202616*aIR-3.26346*aUV));
        Double_t dndl = (part1+part2)/part3;
        Double_t dataY = 2.238 + dndl * 0.128;
        Double_t predY = fRindex->Eval(0.128);
        //predY *= (1 + nu_lambda);   // nuisance parameter for resonance peak wavelength

        Double_t dataE = 0.03*0.3;

        chi2 += (dataY-predY)*(dataY-predY)/dataE/dataE;  // statistic part

        //chi2 += (nu_lambda/sigma_lambda) * (nu_lambda/sigma_lambda);   // systematic part
    }
    else if (option==1) {
        // delta chi2 at 128 for Zhou
        Double_t fa = fRindex->GetParameter(0);
        Double_t fb = fRindex->GetParameter(1);
        Double_t fc = fRindex->GetParameter(2);
        Double_t part1 = 9.57976*(-1.06128*fa - 1.14526*fb - 0.0407477*fc);
        Double_t part2 = TMath::Sqrt(-2+(3/(1-6.3865*(0.0333591*fa + 0.0346538*fb + 0.0065366*fc))));
        Double_t part3 = (1 - 6.3865*(0.0333591*fa + 0.0346538*fb + 0.0065366*fc)) * (1 - 6.3865*(0.0333591*fa + 0.0346538*fb + 0.0065366*fc) );
        Double_t dndl = part1/(part2*part3);
        Double_t dataY = 2.238 + dndl*0.128;
        Double_t predY = fRindex->Eval(0.128);
        //predY *= (1 + nu_lambda);   // nuisance parameter for resonance peak wavelength

        Double_t dataE = 0.03*0.3;

        chi2 += (dataY-predY)*(dataY-predY)/dataE/dataE;

        //chi2 += (nu_lambda/sigma_lambda) * (nu_lambda/sigma_lambda);   // systematic part
    }

    return chi2;
}

void LArRindex::LoadData()
{
    Int_t i = 0;
    gData->SetPoint(i, 0.3612, 1.2326); i++;
    gData->SetPoint(i, 0.3650, 1.2331); i++;
    gData->SetPoint(i, 0.4063, 1.2308); i++;
    gData->SetPoint(i, 0.4358, 1.2297); i++;
    gData->SetPoint(i, 0.4753, 1.2285); i++;
    gData->SetPoint(i, 0.5086, 1.2277); i++;
    gData->SetPoint(i, 0.5461, 1.2269); i++;
    gData->SetPoint(i, 0.5780, 1.2264); i++;
    gData->SetPoint(i, 0.6439, 1.2256); i++;

    for(int j=0; j<i; j++) {
        gData->SetPointError(j, 0, gData->GetY()[j]*0.001);
    }

}

void LArRindex::SetParameters(double* par)
{
    fRindex->SetParameter(0, par[0]);
    fRindex->SetParameter(1, par[1]);
    fRindex->SetParameter(2, par[2]);
}

void LArRindex::Calculate()
{
    double* datax = gData->GetX();
    for (int i=0; i<gData->GetN(); i++) {
        double calc = fRindex->Eval(datax[i]);
        gCalc->SetPoint(i, datax[i], calc);
    }
}


void LArRindex::Plot()
{

    gData128nm->SetPoint(0, 0.128, fRindex->Eval(0.128));
    gData128nm->SetPointError(0, 0, 0.03*0.3);

    gData->SetMarkerStyle(20);
    gData->SetMarkerColor(kBlue+1);
    gData->SetLineColor(kBlue+1);
    gData->SetLineWidth(2);
    gData128nm->SetMarkerStyle(20);
    gData128nm->SetMarkerColor(kBlue+1);
    gData128nm->SetLineColor(kBlue+1);
    gData128nm->SetLineWidth(2);

    TGraphErrors* gDraw = new TGraphErrors();
    const Int_t N = 1000;
    for(int i=0; i<N; i++) {
        Double_t tmpx = 0.11+0.6/1000*i;
        gDraw->SetPoint(i, tmpx, fRindex->Eval(tmpx));     
    }
    gDraw->SetLineColor(kGreen+1);
    gDraw->SetLineWidth(2);

    TCanvas* cc = new TCanvas();
    gData->GetXaxis()->SetLimits(0.11, 0.7);
    gData->GetYaxis()->SetRangeUser(1.2, 1.45);

    gData->Draw("AP");
    gData128nm->Draw("P SAME");
    gDraw->Draw("L SAME");

    cc->SaveAs("test.pdf");
}
