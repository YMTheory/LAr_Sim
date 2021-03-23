#include "LArConfiguration.hh"

typedef LArConfiguration LArConfig;

int LArConfig::use_pullterm = 1;
int LArConfig::rindex_model = 0; // 0 for model 1, 1 for model 2, 2 for model 2 scaling
int LArConfig::use_depolarization = 1;
int LArConfig::fix_absratio = 0;

double LArConfig::factor1 = 0.;
double LArConfig::factor2 = 0.;
double LArConfig::factor3 = 0.;

bool LArConfig::fit_purified = false;

bool LArConfig::m_toyMC = false;
