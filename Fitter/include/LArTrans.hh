#ifndef LArTrans_h
#define LArTrans_h

#include "LArRindex.hh"

#include "TGraphErrors.h"
#include "TF1.h"

class LArTrans
{
    public:
        LArTrans();
        ~LArTrans();

        static void Initialize();

        static void LoadData();
        static void Calculate();
        static double GetChi2();
        static void Plot();

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
        static void SetParameters();

    private:
        static TGraphErrors* gData;
        static TGraphErrors* gCalc;

        static TF1* fRayLength;
        static TF1* fAbs;
        static TF1* fCorr;

    private:
        static int depolarization;
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

};

extern double gRayLength(double* x, double* p);
extern double gRayLength_delta(double* x, double* p);
extern double gAbs(double* x, double* p);
extern double gCorr(double* x, double* p);

#endif
