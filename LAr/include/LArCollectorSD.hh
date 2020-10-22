#ifndef LArCollectorSD_h
#define LArCollectorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "LArCollectorHit.hh"

class G4Step;
class G4HCofThisEvent;

class LArCollectorSD : public G4VSensitiveDetector
{
    public:
        LArCollectorSD(const G4String& name,
                const G4String& hitsCollectionName,
                G4int nofPmts);
        virtual ~LArCollectorSD();

        // methods from base class
        virtual void   Initialize  (G4HCofThisEvent* hitCollection);
        virtual G4bool ProcessHits (G4Step* step, G4TouchableHistory* history);
        virtual void   EndOfEvent  (G4HCofThisEvent* hitCollection);

    private:
        LArCollectorHitsCollection* fHitsCollection;
        G4int fNofPmts;

};

#endif
