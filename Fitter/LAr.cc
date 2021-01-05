#include "LAr.hh"

int main()
{
    LArChiFunction* LArFCN = new LArChiFunction();
    LArFCN->Initialize();
    LArFCN->GetChiSquare(10000);
    LArFCN->Plot();
}
