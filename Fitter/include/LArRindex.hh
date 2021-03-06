#ifndef LArRindex_h
#define LArRindex_h

#include "TGraphErrors.h"
#include "TF1.h"

class LArRindex
{
    public:
        LArRindex() ;
        ~LArRindex();

        static void Initialize();

        static void LoadData();
        static void Calculate();
        static double GetChi2();
        
        static void setp0(double p0) {m_p0 = p0;}
        static double getp0()        {return m_p0;}
        static void setp1(double p1) {m_p1 = p1;}
        static double getp1()        {return m_p1;}
        static void setp2(double p2) {m_p2 = p2;}
        static double getp2()        {return m_p2;}
        static void setrhoratio(double rhoratio) {m_rhoratio = rhoratio;}
        static double getrhoratio()  {return m_rhoratio;}
        static int getoption()       {return option;}
        static void SetParameters();

        static void Plot();

        static double CalcRindex(double wl);

    private:
        static TGraphErrors* gData;
        static TGraphErrors* gData128nm;
        static TGraphErrors* gCalc;

        static TF1* fRindex;

    private:
        static int option; // 0 for Babicz, 1 for ours
        
        static double m_p0;
        static double m_p1;
        static double m_p2;
        static double m_rhoratio;
        
    private:
        static double gRindex128();
        static double gRindex128_20();

};

extern double gRindex20(double* x, double* p);
extern double gRindex(double* x, double* p);
extern double gRindex_scale(double* x, double* p);

#endif
