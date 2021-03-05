#ifndef LArRindex_new_h
#define LArRindex_new_h

#include <TGraphErrors.h>
#include "TF1.h"

class LArRindex_new
{
    public:
        LArRindex_new();
        ~LArRindex_new();

        static void Initialize();

        static void LoadData();
        static void Calculate();
        static double GetChi2();


        static void setrho(double rho) {m_rho = rho;}
        static double getrho()         {return m_rho;}
        static void seta0(double a0)   {m_a0 = a0;}
        static double geta0()          {return m_a0;}
        static void setaUV(double aUV)   {m_aUV = aUV;}
        static double getaUV()          {return m_aUV;}
        static void setaIR(double aIR)   {m_aIR = aIR;}
        static double getaIR()          {return m_aIR;}
        static double CalcRindex(double wl) {return fRindex->Eval(wl);}

        static void SetParameters();    

        static void Plot();

    private:
        static double m_rho;
        static double m_a0;
        static double m_aUV;
        static double m_aIR;

    private:
        static TGraphErrors* gData;
        static TGraphErrors* gCalc;

        static TF1* fRindex;

};

extern double gRindex_new(double *x, double *p);

#endif
