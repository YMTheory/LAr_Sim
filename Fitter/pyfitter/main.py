import numpy as np

from LArRindex import LArRindex



if __name__ == "__main__" :
    
    LArRindex.LoadData()
    LArRindex.Calculate()
    print(LArRindex.GetChi2())



