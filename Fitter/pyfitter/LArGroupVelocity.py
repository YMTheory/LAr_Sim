import numpy as np
from LArRindex import LArRindex


class LArGroupVelocity(object):

    v_data = 7.46
    v_data_err = 0.273

    v_calc = 0
    chi2 = 0

    @staticmethod
    def Calculate():
        rindex128 = LArRindex.rindex_func(0.128)
        a0  = LArRindex.geta0()
        aUV = LArRindex.getaUV()
        aIR = LArRindex.getaIR()

        part1 = -3*(0.323001*aIR + 115.417*aUV)*(a0 - 0.0202616*aIR + 3.26346*aUV)/(3 -a0+ 0.0202616*aIR - 3.26346*aUV)/(3 -a0+ 0.0202616*aIR - 3.26346*aUV);
        part2 = 3*(-0.323001*aIR - 115.417*aUV)/(3-a0 +0.0202616*aIR-3.26346*aUV);
        part3 = 2*np.sqrt(1+(3*(a0-0.0202616*aIR+3.26346*aUV))/(3-a0+0.0202616*aIR-3.26346*aUV));
        dndl = (part1+part2)/part3;
        
        vlight = 299792458 * 1e-9; # m/ns 

        LArGroupVelocity.v_calc = (rindex128 - 0.128 * dndl) / vlight;
        return LArGroupVelocity.v_calc


    @staticmethod
    def GetChi2():
        LArGroupVelocity.chi2 += (LArGroupVelocity.v_data - LArGroupVelocity.v_calc)**2 / LArGroupVelocity.v_data_err**2
        return LArGroupVelocity.chi2


    @staticmethod
    def getDataX():
        return 0.128

    @staticmethod
    def getDataY():
        return LArGroupVelocity.v_data

    @staticmethod
    def getDataYerr():
        return LArGroupVelocity.v_data_err
    

