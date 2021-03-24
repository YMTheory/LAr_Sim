#ifndef LArMathMinimizer_new_h
#define LArMathMinimizer_new_h

#include "LArRindex_new.hh"
#include "LArTrans_new.hh"
#include "LArConfiguration.hh"

class LArMathMinimizer_new
{
    public:
        LArMathMinimizer_new();
        ~LArMathMinimizer_new();

        void Initialize();
        
        static double GetChi2(const double *par);

        //static int Minimization(char * minName ,
        //         const char *algoName ,
        //         int randomSeed );

        static int Minimization(int num, int* ivar, double *init);

        static bool Plot();

        static void Profile1D(int index, double min, double max, double step, double CI);
        static void Profile2D(int* index, double *min, double *max, int *num);
    
    public:
        static void setlagRindex(double lagRindex) {m_lag_rindex = lagRindex;}
        static double getlagRindex()               {return m_lag_rindex;}
        static void setlagRayL(double lagrayL)     {m_lag_rayL = lagrayL;}
        static double getlagRayL()                 {return m_lag_rayL;}

        static void setDelta(bool delta)           {m_delta = delta;}
        static bool getDelta()                     {return m_delta;}

        static double getChiMin()                  {return m_chi2Min;}

    private:
        static bool m_fit_purified;

        static double m_chi2Min;
        static int m_npar;
        static double m_bestFit[20];

    private:
        static double m_lag_rindex;
        static double m_lag_rayL;

        static bool m_delta;

};

#endif
