#include "LArRindex.hh"
#include "LArGroupVelocity.hh"

double LArGroupVelocity::m_data;
double LArGroupVelocity::m_data_err;
double LArGroupVelocity::m_calc;
double LArGroupVelocity::m_p0;
double LArGroupVelocity::m_p1;
double LArGroupVelocity::m_p2;

LArGroupVelocity::LArGroupVelocity()
{}

LArGroupVelocity::~LArGroupVelocity()
{}

void LArGroupVelocity::Initialize()
{
    setdata(7.46);      // group velocity 
    setdataerr(0.273);   // group velocity error
}

void LArGroupVelocity::Calculate()
{
    double rindex = LArRindex::CalcRindex(0.128);
    setp0(LArRindex::getp0());
    setp1(LArRindex::getp1());
    setp2(LArRindex::getp2());

    if (LArRindex::getoption() == 0) {
        double a0 = m_p0;
        double aUV = m_p1;
        double aIR = m_p2;
        Double_t part1 = -3*(0.323001*aIR + 115.417*aUV)*(a0 - 0.0202616*aIR + 3.26346*aUV)/(3 -a0+ 0.0202616*aIR - 3.26346*aUV)/(3 -a0+ 0.0202616*aIR - 3.26346*aUV);
        Double_t part2 = 3*(-0.323001*aIR - 115.417*aUV)/(3-a0 +0.0202616*aIR-3.26346*aUV);
        Double_t part3 = 2*TMath::Sqrt(1+(3*(a0-0.0202616*aIR+3.26346*aUV))/(3-a0+0.0202616*aIR-3.26346*aUV));
        Double_t dndl = (part1+part2)/part3;
        
        double vlight = 299792458 * 1e-9; // m/ns 

        m_calc = (rindex - 0.128 * dndl) / vlight;
    } else if (LArRindex::getoption() == 1) {
        double fa = m_p0;
        double fb = m_p1;
        double fc = m_p2;
        Double_t part1 = 9.57976*(-1.06128*fa - 1.14526*fb - 0.0407477*fc);
        Double_t part2 = TMath::Sqrt(-2+(3/(1-6.3865*(0.0333591*fa + 0.0346538*fb + 0.0065366*fc))));
        Double_t part3 = (1 - 6.3865*(0.0333591*fa + 0.0346538*fb + 0.0065366*fc)) * (1 - 6.3865*(0.0333591*fa + 0.0346538*fb + 0.0065366*fc) );
        Double_t dndl = part1/(part2*part3);

        double vlight = 299792458 * 1e-9; // m/ns 

        m_calc = (rindex - 0.128 * dndl) / vlight;
    } else if (LArRindex::getoption() == 2) {
        Double_t m_rhoratio = LArRindex::getrhoratio();
        Double_t dndl = -0.005356*m_rhoratio / ( TMath::Sqrt(-2+3/(1-0.0002948*m_rhoratio)) * (1-0.000294811*m_rhoratio)*(1-0.000294811*m_rhoratio) );
        double vlight = 299792458 * 1e-9; // m/ns 

        m_calc = (rindex - 0.128 * dndl) / vlight;

    }
}

double LArGroupVelocity::GetChi2()
{
    Calculate();

    double chi2 = (m_data - m_calc) * (m_data - m_calc) / m_data_err / m_data_err;
    return chi2;
}





