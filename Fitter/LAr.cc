#include "LAr.hh"

int main()
{
    LArRindex* exp1 = new LArRindex(); 

    LArTrans* exp2 = new LArTrans(exp1);
    

    delete exp1;
    delete exp2;
}
