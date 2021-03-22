double gRindex20(Double_t* x, Double_t* p)
{
    double lUV = 0.1066; //um
    double lIR = 0.9083; //um
    
    double l = x[0];
    double a0 = p[0];
    double aUV = p[1];
    double aIR = p[2];
    
    double_t A = a0 + aUV*l*l/(l*l-lUV*lUV) + aIR*l*l/(l*l-lIR*lIR) ;

    return TMath::Sqrt(1+3*A/(3-A));
}
TF1* fRindex20 = new TF1("fRindex20", gRindex20, 0.1, 0.7, 3);


double gRayLength(double* x, double* p)
{
    double l = x[0];  // wavelength um
    double rindex = p[0];  // rindex
    double delta = p[1];

    double kT = 2.24442E-9;
    double kB = 1.380649E-23;
    double T = 90; // K
    double f = 1e22;
    
    double pi = TMath::Pi();

    double rayL = 1 / (8*TMath::Power(pi, 3)/3/TMath::Power(l, 4)
                * ((rindex*rindex-1)*(rindex*rindex+2)/3)*((rindex*rindex-1)*(rindex*rindex+2)/3) * kT * kB * T * f * (6+3*delta)/(6-7*delta));

    return rayL;
}

TF1* fRayLength = new TF1("fRayLength", gRayLength, 0.118, 0.150, 2);

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
    Double_t A1 = p[0];
    Double_t mu1 = p[1];
    Double_t sigma1 = p[2];
    double A2 = p[3];
    Double_t mu2 = p[4];
    Double_t sigma2 = p[5];

    Double_t A_abs = A1*TMath::Exp(-(l-mu1)*(l-mu1)/2/sigma1/sigma1) + A2*TMath::Exp(-(l-mu2)*(l-mu2)/2/sigma2/sigma2);
    Double_t T_abs = TMath::Exp( -A_abs*TMath::Log(10.) );
    return T_abs;
}
TF1* fAbs = new TF1("fAbs", gAbs, 0.125, 0.15, 6);

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
    gRindexData->SetPoint(i, 0.1280, 1.3570); i++;
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



TGraphErrors* gRindexSample = new TGraphErrors();
TGraphErrors* gTransSample = new TGraphErrors();
void sample(int seed)
{
    TRandom3 r; r.SetSeed(seed);
    for (int i=0; i<gRindexData->GetN(); i++) {
        double sample = r.Gaus(gRindexData->GetY()[i], gRindexData->GetEY()[i]);
        gRindexSample->SetPoint(i, gRindexData->GetX()[i], sample);
    }
    for (int i=0; i<gTransData->GetN(); i++) {
        double sample = r.Gaus(gTransData->GetY()[i], gTransData->GetEY()[i]);
        gTransSample->SetPoint(i, gTransData->GetX()[i], sample);
    }
}

double GetChi2(double *par)
{
    double chi2 = 0;
    fRindex20->SetParameters(par[0], par[1], par[2]);
    for(int i=0; i<gRindexSample->GetN(); i++) {
        double pred = fRindex20->Eval(gRindexSample->GetX()[i]);
        chi2 += (pred-gRindexSample->GetY()[i]) * (pred-gRindexSample->GetY()[i]) / gRindexData->GetEY()[i] / gRindexData->GetEY()[i];
    }

    
    for (int i=0; i<gTransSample->GetN(); i++) {
        double wl = gTransSample->GetX()[i];
        double rindex = fRindex20->Eval(gTransSample->GetX()[i]);
        fRayLength->SetParameters(rindex, par[3]);
        fCorr->SetParameter(0, rindex);
        fAbs->SetParameters(par[4], par[5], par[6], par[7], par[8], par[9]);
        
        Double_t d = 5.8;
        Double_t rayL = fRayLength->Eval(wl);
        Double_t T_Ray = TMath::Exp(-d/rayL);
        Double_t corr = fCorr->Eval(wl);
        Double_t T_abs = fAbs->Eval(wl);
        Double_t pred = T_Ray * T_abs * corr;

        chi2 += (pred-gTransSample->GetY()[i]) * (pred-gTransSample->GetY()[i]) / gTransData->GetEY()[i] / gTransData->GetEY()[i];
    }
    

    return chi2;
}


double LHR(double val1, double val2)
{
    return -2*(val1 - val2);
}


void toyMC()
{
    SetData();
    LoadTrans();

    double par0[10] = {3.35312e-01, 9.87354e-02, 7.99851e-03, 2.71931e-01, 3.85373e-01, 1.26513e+02, 9.99938e-01, 4.10933e-01, 1.40121e+02, 1.53676e+00};
    double par1[10] = {3.35649e-01,  9.90729e-02, 9.14486e-03, 0, 4.00371e-01, 1.26507e+02, 1.05158e+00, 4.13184e-01, 1.40109e+02, 1.58052e+00};

    int counter = 0;
    TH1F* hist = new TH1F("hist", "", 100, 0, 400);
    while (counter < 10000) {
        sample(counter+888);
        double m_lhr = LHR( GetChi2(par0), GetChi2(par1) );
        hist->Fill(m_lhr);
        counter++;
    }

    hist->Draw();
}
