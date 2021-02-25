#ifndef LArConfiguration_h
#define LArConfiguration_h

class LArConfiguration
{
    public:
        static int use_pullterm;
        
        static int rindex_model;
        
        static int use_depolarization;
        static int fix_absratio;
        static bool fit_purified;

        // Lagrange multiplier
        static double factor1;  // Xe peak ratio
        static double factor2;  // Rayleigh scattering length
        static double factor3;  // refractive index
};

#endif
