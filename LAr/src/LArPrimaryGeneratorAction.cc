#include "LArPrimaryGeneratorAction.hh"
#include "LArParticleSource.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArPrimaryGeneratorAction::LArPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
    fParticleSource   = new LArParticleSource();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArPrimaryGeneratorAction::~LArPrimaryGeneratorAction()
{
  delete fParticleSource;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LArPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //this function is called at the begining of ecah event
    //

    // In order to avoid dependence of PrimaryGeneratorAction
    // on DetectorConstruction class we get Envelope volume
    // from G4LogicalVolumeStore.

    
    //G4cout << "Genarate Primary Particles ..." << G4endl; 

    
    fParticleSource -> GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

