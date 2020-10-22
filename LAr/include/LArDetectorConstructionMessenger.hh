/*************************************************************************
 @Author: MiaoYu ---> miaoyu@ihep.ac.cn
 @Created Time : Tue Sep 29 14:10:10 2020
 @File Name: LArDetectorConstructionMessenger.hh
 ************************************************************************/

#ifndef LArDetectorConstructionMessenger_h
#define LArDetectorConstructionMessenger_h 1

class LArDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class LArDetectorConstructionMessenger : public G4UImessenger
{
    public:
        LArDetectorConstructionMessenger( LArDetectorConstruction* mpga );
        ~LArDetectorConstructionMessenger();


        virtual void SetNewValue( G4UIcommand *cmd, G4String newValues );
        virtual G4String GetCurrentValue( G4UIcommand* cmd );

    private:
        LArDetectorConstruction*         fTarget;

        G4UIdirectory*                  fDirectory;

        G4UIcmdWithADoubleAndUnit*      fRadiusCmd;
};

#endif
