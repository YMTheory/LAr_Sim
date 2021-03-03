#ifndef LArTrans_new_h
#define LArTrans_new_h

#include "LArRindex_new.hh"

#include "TGraphErrors.h"
#include "TF1.h"

class LArTrans_new
{
    public:
        LArTrans_new();
        ~LArTrans_new();

        static void Initialize();

        static void LoadData();
        static void Calculate();
        static double GetChi2();
        static void Plot();

        static void setrindex(double rindex) {fRayLength->SetParameter(0, rindex);}
        static void setdelta(double delta)   {m_delta = delta;}
        static double getdelta()             {return m_delta;}
        static void setA1(double A1)         {m_A1 = A1;}
        static double getA1()                {return m_A1;}
        static void setmu1(double mu1)       {m_mu1 = mu1;}
        static double getmu1()               {return m_mu1;}
        static void setsigma1(double sigma1) {m_sigma1 = sigma1;}
        static double getsigma1()            {return m_sigma1;}
        static void setA2(double A2)         {m_A2 = A2;}
        static double getA2()                {return m_A2;}
        static void setmu2(double mu2)       {m_mu2 = mu2;}
        static double getmu2()               {return m_mu2;}
        static void setsigma2(double sigma2) {m_sigma2 = sigma2;}
        static double getsigma2()            {return m_sigma2;}
        static void setnuf(double nu_f)      {m_nuf = nu_f;}
        static double getnuf()               {return m_nuf;}
        static int getdepolarization()       {return depolarization;}
        static int getfixratio()             {return fixratio;}
        static double getscale()             {return m_scale;}
        static void setscale(double scale)   {m_scale = scale;}
        static void SetParameters();

        static double CalcRayLength(double wl) {return fRayLength->Eval(wl);}

    private:
        static TGraphErrors* gData;
        static TGraphErrors* gCalc;

        static TF1* fRayLength;
        static TF1* fAbs;
        static TF1* fCorr;

    private:
        static int depolarization;
        static int fixratio;
        static bool m_fit_purified;

        static double m_delta;
        static double m_A1;
        static double m_mu1;
        static double m_sigma1;
        static double m_A2;
        static double m_mu2;
        static double m_sigma2;

        static double m_nuR;
        static double sigma_R;
        static double m_nuf;
        static double sigma_f;

        static double m_scale;

};

extern double gRayLength_new(double* x, double* p);
extern double gRayLength_delta_new(double* x, double* p);
extern double gAbs_new(double* x, double* p);
extern double gCorr_new(double* x, double* p);

#endif
