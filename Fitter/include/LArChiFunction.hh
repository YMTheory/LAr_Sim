#ifndef LArChiFunction_h
#define LArChiFunction_h

#include "LArRindex.hh"
#include "LArTrans.hh"
#include <TMinuit.h>

class LArChiFunction
{
    public:
        LArChiFunction(LArRindex* rdx, LArTrans* trans);
        ~LArChiFunction();

        void Initialize();

        void GetChiSquare();
        void SetParameters(double *par);
        double GetChi2();

        void Plot();

    private:
        void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);
        TMinuit* LArMinuit;
    
    private:
        LArRindex* gRdx;
        LArTrans* gTrans;

};

#endif
