#include "LAr.hh"

using namespace std;

int main(int argc, char* argv[])
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

    //LArMathMinimizer_new * minimizer = new LArMathMinimizer_new();
    //minimizer->Initialize();
    //minimizer->Minimization();
    //for(int i=10000; i<20000; i++) {
    //    LArRindex_new::setseed(i);
    //    LArRindex_new::toyMC();
    //    LArTrans_new::setseed(i);
    //    LArTrans_new::toyMC();
    //    minimizer->Minimization();
    //}

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


    //double chi2[10];
    //double chi2_nodelta[10];
    //double delta_chi2[10];

    //if (argc != 3) { 
    //    cout << "Wrong Argument Number ..." << endl;
    //    return 0; 
    //}

    //int begin = int(atof(argv[1]));
    //int end   = int(atof(argv[2]));

    //// profile
    LArMathMinimizer_new * minimizer = new LArMathMinimizer_new();
    minimizer->Initialize();
    int num = 0;
    int ivar[2] = {1, 2};
    double init[2] = {0.0, 0};
    minimizer->Minimization(num, ivar, init);
    minimizer->Plot();

    //for(int i=begin; i<end; i++) {
    //    std::cout << "Processing " << i / 500. *100 << "% ..." << std::endl;
    //    for(int j=0; j<600; j++) {
    //        init[0] = 0.5/500 *i ;
    //        init[1] = 0.12/600*j + 0.88;
    //        minimizer->Minimization(num, ivar, init);
    //    }
    //}

    



    // toyMC generation on profile distribution:
    //LArRindex_new::toyMC();
    //LArTrans_new::toyMC();

    //LArRindex_new::setMC(true);
    //LArTrans_new::setMC(true);
    //minimizer->setlagRindex(10);
    //minimizer->Minimization();
    //for (int i=0; i<1; i++) {
    //    minimizer->setDelta(true);
    //    minimizer->Minimization();
    //    chi2[i] = minimizer->getChiMin();

    //    //minimizer->setDelta(false);
    //    //minimizer->Minimization();
    //    //chi2_nodelta[i] = minimizer->getChiMin();

    //    //delta_chi2[i] = chi2_nodelta[i] - chi2[i];
    //}

}
