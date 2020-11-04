#include <iostream>
#include <TF1.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TMinuit.h>

using namespace std;

Int_t option = 1; // 0 for Babicz, 1 for Zhou

Double_t rindex128_data;

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


TGraphErrors* gRindexData = new TGraphErrors();
void SetData()
{
    Int_t i = 0;
    //gRindexData->SetPoint(i, 0.1280, 1.3570); i++;
    // triple point: 83.81 K
    //gRindexData->SetPoint(i, 0.3612, 1.2395); i++;
    //gRindexData->SetPoint(i, 0.3650, 1.2392); i++;
    //gRindexData->SetPoint(i, 0.4063, 1.2372); i++;
    //gRindexData->SetPoint(i, 0.4358, 1.2361); i++;
    //gRindexData->SetPoint(i, 0.4753, 1.2349); i++;
    //gRindexData->SetPoint(i, 0.5086, 1.2341); i++;
    //gRindexData->SetPoint(i, 0.5461, 1.2334); i++;
    //gRindexData->SetPoint(i, 0.5780, 1.2328); i++;
    //gRindexData->SetPoint(i, 0.6439, 1.2321); i++;
    
    // 90 K
    gRindexData->SetPoint(i, 0.3612, 1.2326); i++;
    gRindexData->SetPoint(i, 0.3650, 1.2331); i++;
    gRindexData->SetPoint(i, 0.4063, 1.2308); i++;
    gRindexData->SetPoint(i, 0.4358, 1.2297); i++;
    gRindexData->SetPoint(i, 0.4753, 1.2285); i++;
    gRindexData->SetPoint(i, 0.5086, 1.2277); i++;
    gRindexData->SetPoint(i, 0.5461, 1.2269); i++;
    gRindexData->SetPoint(i, 0.5780, 1.2264); i++;
    gRindexData->SetPoint(i, 0.6439, 1.2256); i++;

    for(int j=0; j<i; j++) {
        gRindexData->SetPointError(j, 0, gRindexData->GetPointY(j)*0.001);
    }

}

void  SetParameters(double *par)
{
    if(option == 1) {
        fRindex->SetParameter(0, par[0]);
        fRindex->SetParameter(1, par[1]);
        fRindex->SetParameter(2, par[2]);
    } else if (option == 0) {
        fRindex20->SetParameter(0, par[0]);
        fRindex20->SetParameter(1, par[1]);
        fRindex20->SetParameter(2, par[2]);
    }
}

Double_t GetChi2()
{
    Double_t chi2 = 0;
    for (int i=0; i<gRindexData->GetN(); i++) {
        Double_t dataX = gRindexData->GetPointX(i);
        Double_t dataY = gRindexData->GetPointY(i);
        Double_t predY;
        if(option==1)
            predY = fRindex->Eval(dataX);
        else if (option==0)
            predY = fRindex20->Eval(dataX);
        Double_t dataE = gRindexData->GetEY()[i];
        chi2 += (predY-dataY)*(predY-dataY)/dataE/dataE;
    }

    if (option == 0) {
        
        // delta chi2 at 128 nm for Babicz
        Double_t a0  = fRindex20->GetParameter(0);
        Double_t aUV = fRindex20->GetParameter(1);
        Double_t aIR = fRindex20->GetParameter(2);
        Double_t part1 = -3*(0.323001*aIR + 115.417*aUV)*(a0 - 0.0202616*aIR + 3.26346*aUV)/(3 -a0+ 0.0202616*aIR - 3.26346*aUV)/(3 -a0+ 0.0202616*aIR - 3.26346*aUV);
        Double_t part2 = 3*(-0.323001*aIR - 115.417*aUV)/(3-a0 +0.0202616*aIR-3.26346*aUV);
        Double_t part3 = 2*TMath::Sqrt(1+(3*(a0-0.0202616*aIR+3.26346*aUV))/(3-a0+0.0202616*aIR-3.26346*aUV));
        Double_t dndl = (part1+part2)/part3;
        Double_t dataY = 2.238 + dndl * 0.128;
        Double_t predY = fRindex20->Eval(0.128);
        rindex128_data = dataY;

        Double_t dataE = 0.03*0.3;

        chi2 += (dataY-predY)*(dataY-predY)/dataE/dataE;
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
        rindex128_data = dataY;

        Double_t dataE = 0.03*0.3;

        chi2 += (dataY-predY)*(dataY-predY)/dataE/dataE;
        
    }

    //cout << "Delta_chi2 = " << chi2 << endl;
    return chi2;
}


void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


double GetChiSquare(double maxChi2) 
{
    TMinuit* minuit = new TMinuit();
    minuit->SetFCN(ChisqFCN);
    minuit->SetPrintLevel(1);

    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    minuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    // configure parameters
    if(option == 1) {
        minuit->mnparm(iPar, "fa", 0.2, 0.001, 0., 1., ierrflag); iPar++;
        minuit->mnparm(iPar, "fb", 0.04, 0.001, 0., 1., ierrflag); iPar++;
        minuit->mnparm(iPar, "fc", 4.3, 0.01, 0., 10., ierrflag); iPar++;
    } 
    else if (option==0) {
        minuit->mnparm(iPar, "a0", 0.3, 0.001, 0., 1., ierrflag); iPar++;
        minuit->mnparm(iPar, "a1", 0.1, 0.001, 0., 1., ierrflag); iPar++;
        minuit->mnparm(iPar, "a2", 0.01, 0.001, 0., 1., ierrflag); iPar++;
    }

    // Minimization strategy
    minuit->SetErrorDef(1);
    arglist[0] = 2;
    minuit->mnexcm("SET STR", arglist, 1, ierrflag);

    arglist[0] = 5000; // maxCalls
    arglist[1] = 0.01; // tolerance
    minuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    minuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    //minuit->mnexcm("SHOw COVariance", arglist, 0, ierrflag);

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;
    delete minuit;
    return min;

}

void Draw()
{
    gRindexData->SetMarkerColor(kViolet+1);
    gRindexData->SetMarkerStyle(20);
    gRindexData->SetLineColor(kViolet+1);

    TGraphErrors* gRindex128 = new TGraphErrors();
    cout << "rindex @128 nm " << rindex128_data << endl;
    gRindex128->SetPoint(0, 0.1280, rindex128_data);
    gRindex128->SetPointError(0, 0, 0.03*0.3);
    gRindex128->SetMarkerColor(kViolet+1);
    gRindex128->SetMarkerStyle(20);
    gRindex128->SetLineColor(kViolet+1);

    TGraph* gCalc = new TGraph();
    const Int_t N = 1000;
    for(int i=0; i<N; i++) {
        Double_t tmpx = 0.11+0.6/1000*i;
        if (option == 0)
            gCalc->SetPoint(i, tmpx, fRindex20->Eval(tmpx));     
        else if (option == 1)
            gCalc->SetPoint(i, tmpx, fRindex->Eval(tmpx));     
    }
    gCalc->SetMarkerColor(kGreen+2);
    gCalc->SetMarkerStyle(21);
    gCalc->SetLineColor(kGreen+2);
    gCalc->SetLineWidth(2);

    TCanvas* cc = new TCanvas();
    gRindexData->SetTitle("Chi2; wavelength/um; rindex");
    gRindexData->GetXaxis()->SetLimits(0.11, 0.7);
    gRindexData->GetYaxis()->SetRangeUser(1.2, 1.45);
    gRindexData->Draw("AP");
    gCalc->Draw("L SAME");
    gRindex128->Draw("P SAME");
    //cc->SaveAs("rindexZhou_GV.pdf");

}

void chi2fit()
{
    if (option == 0) 
        cout << "Use Babicz's Model Fitting" << endl;
    else if (option == 1)
        cout << "Use Zhou's Model Fitting" << endl;
    SetData();
    cout << "Total Measurement data number: " << gRindexData->GetN() << endl;

    GetChiSquare(1);

    Draw();
}

