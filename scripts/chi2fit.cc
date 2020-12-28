#include <iostream>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMinuit.h>

using namespace std;

Int_t option = 0; // 0 for Babicz, 1 for Zhou
Int_t constrain = 1; // whether use amplitude ratio constraints

Double_t abs_ratio = 0.947;  // fitting from Raz's data

Double_t rindex128_data;
Double_t ray128_data;

Int_t n_par = 4;
Double_t m_parameters[10]; 
Double_t m_parerror[10];

// systematic uncertainties :
Double_t sigma_lambda = 0.01;   // Argon resonance peak1
Double_t sigma_R = 0.04;        // Xe Absorption peak ratio
Double_t sigma_f = 1;           // Fresnel Correction 
// nuisance parameters:
Double_t nu_lambda = 0;
Double_t nu_R = 0;
Double_t nu_f = 0;

// Lagrange multipler method:
Double_t global_l = 0;

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

    //Double_t rayL = 1 / (8*TMath::Power(pi, 3)/3/TMath::Power(l, 4)
    //            * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f );

    Double_t rayL = 1 / (8*TMath::Power(pi, 3)/3/TMath::Power(l, 4)
                * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f * (6+3*delta)/(6-7*delta));

    return rayL;
}

TF1* fRayLength = new TF1("fRayLength", gRayLength, 0.118, 0.150, 2);


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

TGraphErrors* gTransData = new TGraphErrors();
void LoadTrans()
{
    ifstream in; in.open("./G140ppb.txt");
    string line;
    Double_t m_wl, m_tran, m_tran_err;
    Int_t idx = 0;
    while(getline(in, line)) {
        istringstream ss(line);
        ss >> m_wl >> m_tran >> m_tran_err;
        if( m_wl>125 ) {
            gTransData->SetPoint(idx, m_wl/1000., m_tran);
            gTransData->SetPointError(idx, 0, m_tran_err);
            idx++;
        }
    }
}


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

    return 1/TMath::Power( (1-TMath::Power((n_vac-n_MgF2)/(n_vac+n_MgF2),2))/(1-TMath::Power((n_MgF2-n_LAr)/(n_MgF2+n_LAr) ,2))  ,2);
}

TF1* fCorr = new TF1("fCorr", gCorr, 0.1, 0.2, 1);

Double_t gAbs(Double_t* x, Double_t* p)
{
    Double_t l = x[0]*1000;
    Double_t A2 = p[0];
    Double_t mu1 = p[1];
    Double_t sigma1 = p[2];
    Double_t mu2 = p[3];
    Double_t sigma2 = p[4];
    Double_t A1;
    if (constrain == 0) // two independent peaks
        A1 = p[5];
    else 
        A1 = abs_ratio * A2;

    Double_t A_abs = A1*TMath::Exp(-(l-mu1)*(l-mu1)/2/sigma1/sigma1) + A2*TMath::Exp(-(l-mu2)*(l-mu2)/2/sigma2/sigma2);
    Double_t T_abs = TMath::Exp( -A_abs*TMath::Log(10.) );
    return T_abs;
}
TF1* fAbs = new TF1("fAbs", gAbs, 0.125, 0.15, 6);
    
void  SetParameters(double *par)
{
    if(option == 1) {
        fRindex->SetParameter(0, par[0]);
        fRindex->SetParameter(1, par[1]);
        fRindex->SetParameter(2, par[2]);
        fRayLength->SetParameter(1, par[3]);
        fAbs->SetParameter(0, par[4]);
        fAbs->SetParameter(1, par[5]);
        fAbs->SetParameter(2, par[6]);
        fAbs->SetParameter(3, par[7]);
        fAbs->SetParameter(4, par[8]);
        nu_lambda = par[9];
        nu_R = par[10];
        nu_f = par[11];
        if (constrain == 0)
            fAbs->SetParameter(5, par[12]);
    } else if (option == 0) {
        fRindex20->SetParameter(0, par[0]);
        fRindex20->SetParameter(1, par[1]);
        fRindex20->SetParameter(2, par[2]);
        fRayLength->SetParameter(1, par[3]);
        fAbs->SetParameter(0, par[4]);
        fAbs->SetParameter(1, par[5]);
        fAbs->SetParameter(2, par[6]);
        fAbs->SetParameter(3, par[7]);
        fAbs->SetParameter(4, par[8]);
        nu_lambda = par[9];
        nu_R = par[10];
        nu_f = par[11];
        if (constrain == 0)
            fAbs->SetParameter(5, par[12]);
    }
}

Double_t GetChi2()
{
    Double_t chi2 = 0;

    // rindex data part
    for (int i=0; i<gRindexData->GetN(); i++) {
        Double_t dataX = gRindexData->GetPointX(i);
        Double_t dataY = gRindexData->GetPointY(i);
        Double_t predY;
        if(option==1)
            predY = fRindex->Eval(dataX);
        else if (option==0)
            predY = fRindex20->Eval(dataX);
        Double_t dataE = gRindexData->GetEY()[i];
        //std::cout <<  "rindex fitting: " << dataX << " " << dataY << " " << predY << std::endl;
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
        rindex128_data = predY;
        predY *= (1 + nu_lambda);   // nuisance parameter for resonance peak wavelength

        Double_t dataE = 0.03*0.3;

        chi2 += (dataY-predY)*(dataY-predY)/dataE/dataE;  // statistic part

        chi2 += (nu_lambda/sigma_lambda) * (nu_lambda/sigma_lambda);   // systematic part
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
        rindex128_data = predY;
        predY *= (1 + nu_lambda);   // nuisance parameter for resonance peak wavelength

        Double_t dataE = 0.03*0.3;

        chi2 += (dataY-predY)*(dataY-predY)/dataE/dataE;

        chi2 += (nu_lambda/sigma_lambda) * (nu_lambda/sigma_lambda);   // systematic part
    }

    // transmission part:
    for (Int_t ibin=0; ibin<gTransData->GetN(); ibin++) {
        Double_t wl = gTransData->GetPointX(ibin);
        Double_t rindex;
        if (option == 0)
            rindex = fRindex20->Eval(wl);
        else if (option == 1)
            rindex = fRindex->Eval(wl);
        fRayLength->SetParameter(0, rindex);
        Double_t d = 5.8;
        Double_t rayL = fRayLength->Eval(wl);
        Double_t T_Ray = TMath::Exp(-d/rayL);
        fCorr->SetParameter(0, rindex);
        Double_t corr = fCorr->Eval(wl);
        //fAbs->SetParameters(0.3985, 126.5, 1.04, 0.4122, 140.1, 1.566);
        Double_t T_abs = fAbs->Eval(wl);
        Double_t trans_pred = T_Ray * T_abs * corr;
        trans_pred *= (1+nu_R) * (1+nu_f);   // nuisance parameters
        Double_t trans_data = gTransData->GetPointY(ibin);
        Double_t trans_err  =gTransData->GetEY()[ibin];
        
        chi2 += (trans_pred - trans_data)*(trans_pred-trans_data) / trans_err / trans_err;   // statistic part

        chi2 += (nu_R/sigma_R)*(nu_R/sigma_R);
        chi2 += (nu_f/sigma_f)*(nu_f/sigma_f);
    }

    //cout << "Delta_chi2 = " << chi2 << endl;

    // use Lagrange multiple mothed:
    Double_t n128 = fRindex20->Eval(0.128);
    fRayLength->SetParameter(0, n128);
    chi2 += global_l * fRayLength->Eval(0.128);

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
        minuit->mnparm(iPar, "fa", 0.2, 0.001, 0., 1., ierrflag);     iPar++;
        minuit->mnparm(iPar, "fb", 0.04, 0.001, 0., 1., ierrflag);    iPar++;
        minuit->mnparm(iPar, "fc", 4.3, 0.01, 0., 10., ierrflag);     iPar++;
        minuit->mnparm(iPar, "delta", 0.0, 0.01, 0., 5., ierrflag);   iPar++;
        minuit->mnparm(iPar, "A1", 0.3, 0.01, 0., 1., ierrflag);      iPar++;
        minuit->mnparm(iPar, "mu1",126, 0.1, 123, 129, ierrflag);     iPar++;
        minuit->mnparm(iPar, "sigma1", 1, 0.01, 0.5, 2, ierrflag);    iPar++;
        minuit->mnparm(iPar, "mu2",140, 0.1, 135, 145, ierrflag);     iPar++;
        minuit->mnparm(iPar, "sigma2", 1.5, 0.01, 0.5, 2, ierrflag);  iPar++;
        minuit->mnparm(iPar, "nu_lambda", 0, 0.01, 0., 1., ierrflag); iPar++;
        minuit->mnparm(iPar, "nu_R", 0, 0.01, 0., 1., ierrflag); iPar++;
        minuit->mnparm(iPar, "nu_f", 0, 0.01, 0., 1., ierrflag); iPar++;
        if (constrain == 0)
            minuit->mnparm(iPar, "A1", 0.2, 0.01, 0., 1., ierrflag);  iPar++;

    } 
    else if (option==0) {
        minuit->mnparm(iPar, "a0", 0.3, 0.001, 0., 1., ierrflag);     iPar++;
        minuit->mnparm(iPar, "aUV", 0.1, 0.001, 0., 1., ierrflag);    iPar++;
        minuit->mnparm(iPar, "aIR", 0.01, 0.001, 0., 1., ierrflag);   iPar++;
        minuit->mnparm(iPar, "delta", 0.0, 0.01, 0., 5., ierrflag);   iPar++;
        minuit->mnparm(iPar, "A2", 0.3, 0.01, 0., 0., ierrflag);      iPar++;
        minuit->mnparm(iPar, "mu1",126, 0.1, 123, 129, ierrflag);     iPar++;
        minuit->mnparm(iPar, "sigma1", 1, 0.01, 0.5, 2, ierrflag);    iPar++;
        minuit->mnparm(iPar, "mu2",140, 0.1, 135, 145, ierrflag);     iPar++;
        minuit->mnparm(iPar, "sigma2", 1.5, 0.01, 0.5, 2, ierrflag);  iPar++;
        minuit->mnparm(iPar, "nu_lambda", 0, 0.01, 0., 1., ierrflag); iPar++;
        minuit->mnparm(iPar, "nu_R", 0, 0.01, 0., 1., ierrflag); iPar++;
        minuit->mnparm(iPar, "nu_f", 0, 0.01, 0., 1., ierrflag); iPar++;
        if (constrain == 0)
            minuit->mnparm(iPar, "A1", 0.2, 0.01, 0., 1., ierrflag);  iPar++;
    }

    //minuit->FixParameter(4);

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

    for (Int_t i=0; i<n_par; i++) {
        minuit->GetParameter(i, m_parameters[i], m_parerror[i]);
    }

    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << endl;
    cout << " ====================== " << endl;

    Double_t n128 = fRindex20->Eval(0.128);
    fRayLength->SetParameter(0, n128);
    delete minuit;
    return min;

}

void DrawRindex()
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
    //cc->SaveAs("rindexBa_allerr.png");

}


void DrawTrans()
{
    TGraph* pred_graph = new TGraph();
    for (Int_t i=0; i<1000; i++) {
        Double_t wl = ((150.-125.)/1000*i + 125. )/1000.;
        Double_t rindex;
        if (option == 0)
            rindex = fRindex20->Eval(wl);
        else if (option == 1)
            rindex = fRindex->Eval(wl);
        fRayLength->SetParameter(0, rindex);
        fRayLength->SetParameter(1, m_parameters[3]);
        fCorr->SetParameter(0, rindex);
        Double_t corr = fCorr->Eval(wl);
        Double_t rayL = fRayLength->Eval(wl);
        Double_t T_abs = fAbs->Eval(wl);
        Double_t T_Ray = TMath::Exp(-5.8/(rayL)) ;
        Double_t trans_pred = T_Ray * T_abs * corr;
        pred_graph->SetPoint(i, wl, trans_pred);

        if (i == 128)
            cout << "Rayleigh scattering length @128nm : " << rayL << endl;
    }

    gTransData->SetMarkerColor(kBlue+1);
    gTransData->SetLineColor(kBlue+1);
    gTransData->SetMarkerSize(0.5);
    gTransData->SetMarkerStyle(20);
    gTransData->SetLineWidth(2);

    pred_graph->SetLineColor(kGreen+1);
    pred_graph->SetLineWidth(3);

    TCanvas* cg = new TCanvas();
    gTransData->SetTitle("Transmission; wavelength/um; transmission");
    //gTransData->GetXaxis()->SetLimits(0.11, 0.7);
    //gTransData->GetYaxis()->SetRangeUser(1.2, 1.45);
    gTransData->Draw("AP");
    pred_graph->Draw("L SAME");

    //cg->SaveAs("transBa_allerr.png");
}


void chi2fit()
{
    SetData();
    LoadTrans();

    GetChiSquare(1);

    DrawRindex();
    DrawTrans();

    if (option == 0) 
        cout << "Use Babicz's Model Fitting" << endl;
    else if (option == 1)
        cout << "Use Zhou's Model Fitting" << endl;
    cout << "Total Measurement data number: " << gRindexData->GetN()+gTransData->GetN()+1 << endl;

}

