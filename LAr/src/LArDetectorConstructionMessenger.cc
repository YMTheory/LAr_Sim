#include "LArDetectorConstructionMessenger.hh"
#include "LArDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

LArDetectorConstructionMessenger::
LArDetectorConstructionMessenger( LArDetectorConstruction* mpga )
    :fTarget(mpga) {

    fDirectory = new G4UIdirectory("/LAr/det/");
    fDirectory->SetGuidance("detector setup commands");

    fRadiusCmd = new G4UIcmdWithADoubleAndUnit("/LAr/det/setRadius", this);
    fRadiusCmd->SetGuidance("Set CD Radius");
    fRadiusCmd->SetParameterName("Radius", true, true);
    fRadiusCmd->SetDefaultValue(85);
    fRadiusCmd->SetDefaultUnit("cm");

}   

LArDetectorConstructionMessenger::
~LArDetectorConstructionMessenger() {
    delete fRadiusCmd;
}


void LArDetectorConstructionMessenger::SetNewValue(
        G4UIcommand* cmd, G4String newValue) {

    if ( cmd == fRadiusCmd ) {
        fTarget->setRadius( fRadiusCmd->GetNewDoubleValue(newValue) );
    }

    else 
        G4cout << "Error: Unknown Command " << G4endl;

}


