#include "LArSteppingAction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArSteppingAction::LArSteppingAction()
: G4UserSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArSteppingAction::~LArSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LArSteppingAction::UserSteppingAction(const G4Step* step)
{
    //G4cout << step->GetTrack()->GetTrackID() << " " 
    //    << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << " "
    //    << step->GetPreStepPoint() ->GetPhysicalVolume()->GetName() << " " << G4endl;
    //    << step->GetPostStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

