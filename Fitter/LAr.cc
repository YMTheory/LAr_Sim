#include "LAr.hh"

int main()
{
    //LArChiFunction* LArFCN = new LArChiFunction();
    //LArFCN->Initialize();
    //LArFCN->setfactor1(0);
    //LArFCN->setfactor2(0);
    //LArFCN->setfactor3(0);
    //LArFCN->GetChiSquare(10000);

    //for (int i=0; i<100; i++) {
    //        LArFCN->setfactor3(-1000+20*i);
    //        LArFCN->GetChiSquare(10000);
    //}
    
    //LArFCN->Plot();

    //LArFCN->Profile1D(3, 0.0, 0.5, 0.0001, 5);
    
    //double min[2]  = {126.013, 139.621};
    //double max[2]  = {127.013, 140.621};
    //double step[2] = {0.001, 0.001}; // 1000 bins for x and y
    //LArFCN->Profile2D(min, max, step);
    
    //LArFCN->Lray_profile();

    //LArMathMinimizer* minimizer;
    //minimizer->Initialize();
    //minimizer->Minimization();
    //minimizer->Plot();

    LArMathMinimizer_new * minimizer = new LArMathMinimizer_new();
    minimizer->Initialize();
    for(int i=10000; i<20000; i++) {
        LArRindex_new::setseed(i);
        LArRindex_new::toyMC();
        LArTrans_new::setseed(i);
        LArTrans_new::toyMC();
        minimizer->Minimization();
    }

    //for (int j=0; j<200; j++) {
    //    std::cout << "Processing " << j/2. << "% ...." << std::endl;
    //    for(int i=0; i<200; i++) {
    //        minimizer->setlagRindex(1000/200*j - 500);
    //        minimizer->setlagRayL(2/200.*i-1);
    //        minimizer->Minimization();
    //    } 
    //}

    //minimizer->Profile1D(1, 0.0, 0.6, 0.001, 1);
    //minimizer->Profile1D(2, 0.90, 0.99, 0.0001, 1);
    //int index[2] = {1, 2};
    //double min[2] = {0.1, 0.88};
    //double max[2] = {0.4, 0.99};
    //int num[2] = {1000, 1000};
    //minimizer->Profile2D(index, min, max, num);

    //minimizer->Plot();
}
