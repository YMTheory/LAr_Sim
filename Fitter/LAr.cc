#include "LAr.hh"

int main()
{
    LArChiFunction* LArFCN = new LArChiFunction();
    LArFCN->Initialize();
    LArFCN->setfactor1(0);
    LArFCN->setfactor2(0);
    LArFCN->setfactor3(0);
    LArFCN->GetChiSquare(10000);

    //for (int i=0; i<100; i++) {
    //        LArFCN->setfactor3(-1000+20*i);
    //        LArFCN->GetChiSquare(10000);
    //}
    
    LArFCN->Plot();

    //LArFCN->Profile1D(3, 0.0, 0.5, 0.0001, 5);
    
    //double min[2]  = {126.013, 139.621};
    //double max[2]  = {127.013, 140.621};
    //double step[2] = {0.001, 0.001}; // 1000 bins for x and y
    //LArFCN->Profile2D(min, max, step);
    
    //LArFCN->Lray_profile();

}
