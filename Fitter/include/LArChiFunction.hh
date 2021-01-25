#ifndef LArChiFunction_h
#define LArChiFunction_h

#include "LArRindex.hh"
#include "LArTrans.hh"
#include "LArGroupVelocity.hh"
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
        static void Lray_profile();

        static void setfactor1(double factor1) {m_factor1 = factor1;}
        static void setfactor2(double factor2) {m_factor2 = factor2;}
        static void setfactor3(double factor3) {m_factor3 = factor3;}

    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);
        TMinuit* LArMinuit;
    
    private:
        static int usePull;

        static double m_chi2;
        static double m_chi2Min;
        static bool   m_DoFit;
        static int    m_nParameters;
        static double m_bestFit[20];
        static double m_fitError[20];

        static double m_factor1;
        static double m_factor2;
        static double m_factor3;

        static double m_ratio;

};

#endif
