#include "LAr.hh"

int main()
{
    LArChiFunction* LArFCN = new LArChiFunction();
    LArFCN->Initialize();
    LArFCN->GetChiSquare(10000);
    LArFCN->Plot();

    //LArFCN->Profile1D(3, 0.0, 0.45, 0.0001, 5);
    
    double min[2]  = {0.90, 0.0};
    double max[2]  = {1.0, 0.4};
    double step[2] = {0.001, 0.0004};
    LArFCN->Profile2D(min, max, step);

}
