#include "LArRunAction.hh"
#include "LArEventAction.hh"
#include "LArAnalysisManager.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArRunAction::LArRunAction()
: G4UserRunAction()
{
    // set printing event number per each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);
    LArAnalysisManager* analysis = LArAnalysisManager::getInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArRunAction::~LArRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LArRunAction::BeginOfRunAction(const G4Run*)
{
    G4cout << "Begin of One Run" << G4endl;

    LArAnalysisManager* analysis = LArAnalysisManager::getInstance();
    analysis->book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LArRunAction::EndOfRunAction(const G4Run* )
{
    LArAnalysisManager* analysis = LArAnalysisManager::getInstance();
    analysis->finish();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

