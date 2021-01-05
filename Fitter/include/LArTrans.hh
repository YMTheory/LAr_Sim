#ifndef LArTrans_h
#define LArTrans_h

#include "LArRindex.hh"

#include "TGraphErrors.h"
#include "TF1.h"

class LArTrans
{
    public:
        LArTrans(LArRindex* rdx);

        void LoadData();
        void SetParameters(double *par);
        void Calculate();
        double GetChi2();
        void Plot();

    private:
        TGraphErrors* gData;
        TGraphErrors* gCalc;

        TF1* fRayLength;
        TF1* fAbs;
        TF1* fCorr;

        int depolarization;

    private:
        LArRindex* gRdx;
};

extern double gRayLength(double* x, double* p);
extern double gRayLength_delta(double* x, double* p);
extern double gAbs(double* x, double* p);
extern double gCorr(double* x, double* p);

#endif
