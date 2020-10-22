#ifndef LArAnalysisMessenger_h
#define LArAnalysisMessenger_h 1

#include "globals.hh"
#include "LArAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImessenger.hh"

class LArAnalysisManager;

class LArAnalysisMessenger : public G4UImessenger
{

    public:
        LArAnalysisMessenger(LArAnalysisManager*);
        ~LArAnalysisMessenger();
            
        void SetNewValue(G4UIcommand*, G4String);

    private:
        // pointer to LArAnalysisManager
        LArAnalysisManager* LArAnalysis;
        G4UIdirectory* LArAnalysisDir;
        G4UIcmdWithAString* outputFileCmd;

};
#endif
