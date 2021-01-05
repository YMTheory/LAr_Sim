#ifndef LArChiFunction_h
#define LArChiFunction_h

#include "LArRindex.hh"
#include "LArTrans.hh"
#include <TMinuit.h>

class LArChiFunction
{
    public:
        LArChiFunction();
        ~LArChiFunction();

        void Initialize();

        double GetChiSquare(double maxChi2);
        static void SetParameters(double *par);
        static double GetChi2();

        static void Plot();

    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);
        TMinuit* LArMinuit;
    
    private:
        static double m_chi2;
        static double m_chi2Min;
        static bool m_DoFit;


};

#endif
