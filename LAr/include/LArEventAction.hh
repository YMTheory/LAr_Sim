#ifndef LArEventAction_h
#define LArEventAction_h 1

#include "G4UserEventAction.hh"
#include "LArCollectorHit.hh"
#include "globals.hh"

#include <vector>
#include <array>

/// Event action class
///

class LArEventAction : public G4UserEventAction
{
  public:
    LArEventAction();
    virtual ~LArEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);


  private:
    // methods
    LArCollectorHitsCollection* GetHitsCollection(G4int hcID,
                                           const G4Event* event) const;
    //void PrintEventStatistics(G4int pmtTrackID) const;

    // data members
    G4int fCollectorID;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
