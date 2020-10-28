#ifndef LArSteppingAction_h
#define LArSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class LArEventAction;

class G4LogicalVolume;

/// Stepping action class
/// 

class LArSteppingAction : public G4UserSteppingAction
{
  public:
    LArSteppingAction();
    virtual ~LArSteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
