#ifndef LArPrimaryGeneratorAction_h
#define LArPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "LArParticleSource.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued 
/// in front of the phantom across 80% of the (X,Y) phantom size.

class LArPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    LArPrimaryGeneratorAction();    
    virtual ~LArPrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         
  
  private:
    LArParticleSource* fParticleSource;
    G4int NumberOfParticlesToBeGenerated;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
