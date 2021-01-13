#ifndef LArChiFunction_h
#define LArChiFunction_h

#include "LArRindex.hh"
#include "LArTrans.hh"
#include <TMinuit.h>
#include <TGraph.h>

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

        static void Profile1D(int index, double min, double max, double step, double CI);
        static void Profile2D(double* min, double* max, double* step);

    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);
        TMinuit* LArMinuit;
    
    private:
        static double m_chi2;
        static double m_chi2Min;
        static bool   m_DoFit;
        static int    m_nParameters;
        static double m_bestFit[20];
        static double m_fitError[20];


};

#endif
