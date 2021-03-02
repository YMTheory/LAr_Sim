#include "LArRindex_new.hh"

using namespace std;

double LArRindex_new::m_rho = 34;
TGraphErrors* LArRindex_new::gData = new TGraphErrors();
TGraphErrors* LArRindex_new::gCalc = new TGraphErrors();

LArRindex_new::LArRindex_new()
{;}

LArRindex_new::~LArRindex_new()
{;}

void LArRindex_new::Initialize()
{
    //std::cout << "==========> Rindex Initialization " << std::endl;

    LoadData();
}

void LArRindex_new::LoadData()
{
    double m_wavelength[9] = {0.3612, 0.3650, 0.4063, 0.4358, 0.4753, 0.5086, 0.5461, 0.5780, 0.6439};  // um
    double m_rindex0[9] = {1.2395, 1.2392, 1.2372, 1.2361, 1.2349, 1.2341, 1.2334, 1.2328, 1.2321}; // 83.81 K
    double m_rindex1[9] = {1.2370, 1.2367, 1.2347, 1.2336, 1.2324, 1.2316, 1.2308, 1.2303, 1.2296}; // 86 K
    double m_rindex2[9] = {1.2349, 1.2346, 1.2326, 1.2315, 1.2303, 1.2295, 1.2287, 1.2282, 1.2274}; // 88 K

    for(int i=0; i<9; i++) {
        double mean = (m_rindex0[i] + m_rindex1[i] + m_rindex2[i]) / 3.;
        double err  = max(abs(m_rindex0[i]-mean), abs(m_rindex2[i]-mean));
        gData->SetPoint(i, m_wavelength[i], mean);
        gData->SetPointError(i, 0, err);
    }
}

double LArRindex_new::GetChi2()
{;}

void LArRindex_new::Calculate()
{;}
