#ifndef LArAnalysisManager_h
#define LArAnalysisManager_h 1

#include "globals.hh"
#include "LArAnalysisMessenger.hh"
#include <vector>
#include "G4ThreeVector.hh"

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
        void AddInitPos(G4ThreeVector);
        void ClearInitPos();
        void analyseAddNtupleRow();

        int getPosVecLength() {return vecInitPosX.size();}


    private:
        LArAnalysisManager();

        G4String outputFileName;

        static LArAnalysisManager* instance;
        
        LArAnalysisMessenger* analysisMessenger;

    private:
        std::vector<double> vecInitPosX;
        std::vector<double> vecInitPosY;
        std::vector<double> vecInitPosZ;

};

#endif
