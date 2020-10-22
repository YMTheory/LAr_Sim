#include "LArPrimaryGeneratorAction.hh"
//#include "LArParticleSource.hh"

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
    fParticleGun      = new G4ParticleGun();
    //fParticleSource   = new LArParticleSource();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArPrimaryGeneratorAction::~LArPrimaryGeneratorAction()
{
  delete fParticleGun;
  //delete fParticleSource;
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

    
    // Generate more than one particle each time
    NumberOfParticlesToBeGenerated = 10000;
    fParticleGun = new G4ParticleGun(NumberOfParticlesToBeGenerated);

    const G4double pi = 3.141592653;

    // particle definition
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName="opticalphoton");

    for( G4int iP = 0; iP<NumberOfParticlesToBeGenerated; iP++ ) {

        // set particle type
        fParticleGun->SetParticleDefinition(particle);

        // set optical photon energy/wavelength
        fParticleGun->SetParticleEnergy(2*eV);

        // set momentum direction
        G4double mom_x, mom_y, mom_z, mom_theta, mom_phi;
        while(1)  {
            mom_theta = G4UniformRand() * pi;
            mom_phi   = G4UniformRand() * pi * 2;
            mom_x     = sin(mom_theta) * cos(mom_phi);
            if ( mom_x >0 ) break;
        }
        mom_y  =  sin(mom_theta) * sin(mom_phi);
        mom_z  =  cos(mom_theta);

        fParticleGun->SetParticleMomentumDirection( G4ThreeVector(mom_x, mom_y, mom_z));

        fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));

        fParticleGun->GeneratePrimaryVertex( anEvent );
    }
     

    fParticleGun->GeneratePrimaryVertex(anEvent);
    //fParticleSource -> GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

