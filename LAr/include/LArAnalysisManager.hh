#ifndef LArAnalysisManager_h
#define LArAnalysisManager_h 1

#include "globals.hh"
#include "LArAnalysisMessenger.hh"

class LArAnalysisMessenger;

class LArAnalysisManager  {
    
    public:

        //LArAnalysisManager();
        virtual ~LArAnalysisManager();

        void book();

        void finish();
        
        //method to call to create an instance of this class
        static LArAnalysisManager* getInstance();
        
        void SetOutputFileName(G4String);

        void analysePhotonNumber(G4int number);
        void analyseEventID(G4int evtid);
        void analyseInitPhotonNumber(G4int number);
        void analyseAddNtupleRow();


    private:
        LArAnalysisManager();

        G4String outputFileName;

        static LArAnalysisManager* instance;
        
        LArAnalysisMessenger* analysisMessenger;

};

#endif
