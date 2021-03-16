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

        static int Minimization();

        static bool Plot();

        static void Profile1D(int index, double min, double max, double step, double CI);
        static void Profile2D(int* index, double *min, double *max, int *num);

    private:
        static bool m_fit_purified;

        static double m_chi2Min;
        static int m_npar;
        static double m_bestFit[20];
};

#endif
