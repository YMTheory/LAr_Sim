#include "LArPhysicsList.hh"

#include "G4OpticalPhysics.hh"
#include "G4SystemOfUnits.hh"

// particles
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"

LArPhysicsList::LArPhysicsList() : G4VModularPhysicsList()
{
    defaultCutValue = 1.0*mm;

}

LArPhysicsList::~LArPhysicsList() {}

void LArPhysicsList::SetCuts() {
    SetCutsWithDefault();
}


void LArPhysicsList::ConstructParticle()
{
    G4OpticalPhoton::OpticalPhotonDefinition();
}

void LArPhysicsList::ConstructProcess()
{
    AddTransportation();
    ConstructOpticalProcess();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"

void LArPhysicsList::ConstructOpticalProcess()
{
    G4OpAbsorption* theAbsProcess     = new G4OpAbsorption();
    G4OpRayleigh* theRayProcess       = new G4OpRayleigh();
    G4OpBoundaryProcess* theBdProcess = new G4OpBoundaryProcess();
    theBdProcess->SetInvokeSD(false);
    auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    while( (*particleIterator)() ){

        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        if (theAbsProcess->IsApplicable(*particle)) {
            pmanager -> AddDiscreteProcess(theAbsProcess);
            G4cout << " ===> Registered Absorption Process " << G4endl;
        }

        if(theRayProcess->IsApplicable(*particle)) {
            pmanager -> AddDiscreteProcess(theRayProcess);
            G4cout << " ===> Registered Rayleigh Scatterinng Process " << G4endl;
        }

        if(theBdProcess->IsApplicable(*particle)) {
            pmanager -> AddDiscreteProcess(theBdProcess);
            G4cout << " ===> Registered Boundary Process " << G4endl;
        }
    }
}
