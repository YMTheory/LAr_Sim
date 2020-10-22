/*************************************************************************
 @Author: MiaoYu ---> miaoyu@ihep.ac.cn
 @Created Time : Wed Sep 30 09:13:41 2020
 @File Name: LArAnalysisMessenger.cc
 ************************************************************************/

#include "LArAnalysisMessenger.hh"

LArAnalysisMessenger::LArAnalysisMessenger(
        LArAnalysisManager* analysisManager)
    : LArAnalysis(analysisManager)
{
    LArAnalysisDir = new G4UIdirectory("/LAr/analy/");
    LArAnalysisDir->SetGuidance("analysis control");

    outputFileCmd = new G4UIcmdWithAString("/LAr/analy/outputFile", this);
    outputFileCmd->SetGuidance("specify output file name");
    outputFileCmd->SetParameterName("choice", true);
    outputFileCmd->SetDefaultValue("test");
    //outputFileCmd->AvailableForStates(G4State_Idle);

}

LArAnalysisMessenger::~LArAnalysisMessenger()
{
    delete LArAnalysisDir;
    delete outputFileCmd;
}


void LArAnalysisMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{
    if ( cmd == outputFileCmd ) {
        LArAnalysis->SetOutputFileName( newValue );
    }
    
    else
        G4cout << "Error: Unknown Command !" << G4endl;
}
