#include "LArActionInitialization.hh"
#include "LArPrimaryGeneratorAction.hh"
#include "LArRunAction.hh"
#include "LArEventAction.hh"
#include "LArSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArActionInitialization::LArActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArActionInitialization::~LArActionInitialization()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LArActionInitialization::Build() const
{
    SetUserAction(new LArPrimaryGeneratorAction);


    LArEventAction* eventAction = new LArEventAction();
    SetUserAction(eventAction);

    LArRunAction* runAction = new LArRunAction();
    SetUserAction(runAction);

    LArSteppingAction* stepAction = new LArSteppingAction();
    SetUserAction(stepAction);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
