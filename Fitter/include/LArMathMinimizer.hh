#ifndef LArMathMinimizer_h
#define LArMathMinimizer_h

#include "LArRindex.hh"
#include "LArTrans.hh"
#include "LArGroupVelocity.hh"
#include "LArConfiguration.hh"

class LArMathMinimizer
{
    public:
        LArMathMinimizer();
        ~LArMathMinimizer();

        void Initialize();
        
        static double GetChi2(const double *par);

        //static int Minimization(char * minName ,
        //         const char *algoName ,
        //         int randomSeed );

        static int Minimization();

        static bool Plot();

    private:
        static bool m_fit_purified;
};

#endif
