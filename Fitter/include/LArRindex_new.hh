#ifndef LArRindex_new_h
#define LArRindex_new_h

#include <TGraphErrors.h>
#include "TF1.h"
#include "TRandom.h"

class LArRindex_new
{
    public:
        LArRindex_new();
        ~LArRindex_new();

        static void Initialize();

        static void LoadData();
        static void Calculate();
        static double GetChi2();


        static void settemp(double temp)    {m_temp = temp;}
        static double gettemp()             {return m_temp;}
        static void setrho(double rho)      {m_rho = rho;}
        static double getrho()              {return m_rho;}
        static void seta0(double a0)        {m_a0 = a0;}
        static double geta0()               {return m_a0;}
        static void setaUV(double aUV)      {m_aUV = aUV;}
        static double getaUV()              {return m_aUV;}
        static void setaIR(double aIR)      {m_aIR = aIR;}
        static double getaIR()              {return m_aIR;}
        static void setp0(double p0)        {m_p0 = p0;}
        static double getp0()               {return m_p0;}
        static void setp1(double p1)        {m_p1 = p1;}
        static double getp1()               {return m_p1;}
        static double getp0_init()          {return p0;}
        static void setp0_init(double p0_init)               {p0 = p0_init;}
        static double getp1_init()                      {return p1;}
        static void setp1_init(double p1_init)               {p1 = p1_init;}
        static double CalcRindex(double wl) {return fRindex->Eval(wl);}

        static void setseed(double seed)  {gRandom->SetSeed(seed);}
        static void setMC(bool MC)        {m_toyMC = MC;}

        static void SetParameters();    

        static void toyMC();

        static void Plot();

    private:
        static double m_p0;
        static double m_p1;
        static double p0;
        static double p1;
        static double sigma_p0;
        static double sigma_p1;

        static double m_temp;
        static double temp;
        static double sigma_temp;

        static double m_rho;
        static double m_a0;
        static double m_aUV;
        static double m_aIR;

        static bool   m_loadData;
        static bool   m_toyMC;


    private:
        static TGraphErrors* gData;
        static TGraphErrors* gCalc;
        static TGraphErrors* gtoyMC;

        static TF1* fRindex;

};

extern double m_rho90K;
extern double gRindex_new(double *x, double *p);

#endif
