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


TGraphErrors* gRindexData = new TGraphErrors();
void SetData()
{

    Int_t i = 0;
    
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

void SetParameters(double *par)
{
    fRindex->SetParameter(0, par[0]);
    fRindex->SetParameter(1, par[1]);
    fRindex->SetParameter(2, par[2]);
}

double GetChi2()
{
    double chi2 = 0;

    for (int i=0; i<gRindexData->GetN(); i++) {

        Double_t dataX = gRindexData->GetPointX(i);
        Double_t dataY = gRindexData->GetPointY(i);
        Double_t predY = fRindex->Eval(dataX);
        Double_t dataE = gRindexData->GetEY()[i];
        chi2 += (predY-dataY)*(predY-dataY)/dataE/dataE;
    }

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

    minuit->mnparm(iPar, "fa", 0.2, 0.001, 0., 1., ierrflag);     iPar++;
    minuit->mnparm(iPar, "fb", 0.04, 0.001, 0., 1., ierrflag);    iPar++;
    minuit->mnparm(iPar, "fc", 4.3, 0.01, 0., 10., ierrflag);     iPar++;


    // Minimization strategy
    minuit->SetErrorDef(1);
    arglist[0] = 2;
    minuit->mnexcm("SET STR", arglist, 1, ierrflag);

    arglist[0] = 50000; // maxCalls
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


void fit_test()
{
    SetData();

    GetChiSquare(1);
}
