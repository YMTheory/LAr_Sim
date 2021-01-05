#ifndef LArRindex_h
#define LArRindex_h

#include "TGraphErrors.h"
#include "TF1.h"

class LArRindex
{
    public:
        LArRindex() ;

        void LoadData();
        void Calculate();
        double GetChi2();
        void SetParameters(double* par);
        double GetParameter(int i) {return fRindex->GetParameter(i);}
        void Plot();

        double CalcRindex(double wl) {return fRindex->Eval(wl);}

    private:
        TGraphErrors* gData;
        TGraphErrors* gData128nm;
        TGraphErrors* gCalc;

        TF1* fRindex;

    private:
        Int_t option; // 0 for Babicz, 1 for ours

};

extern double gRindex20(double* x, double* p);
extern double gRindex(double* x, double* p);

#endif
