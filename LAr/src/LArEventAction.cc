#include "LArEventAction.hh"
#include "LArAnalysisManager.hh"
#include "LArCollectorHit.hh"
#include "LArCollectorSD.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4VHitsCollection.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArEventAction::LArEventAction()
: G4UserEventAction(),
  fCollectorID(-1)
{
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArEventAction::~LArEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LArCollectorHitsCollection* 
LArEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<LArCollectorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("LArEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LArEventAction::BeginOfEventAction(const G4Event* )
{
    //G4cout << "Begin of Event " << evt->GetEventID() << G4endl;
    LArAnalysisManager* analysis = LArAnalysisManager::getInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LArEventAction::EndOfEventAction(const G4Event* event)
{   

    LArAnalysisManager* analysis = LArAnalysisManager::getInstance();
    G4int evtid = event->GetEventID();
    analysis->analyseEventID( evtid );

    //// Get hits collections IDs (only once)
    if(fCollectorID == -1)  {
        fCollectorID
          = G4SDManager::GetSDMpointer()->GetCollectionID("CollectorHitsCollection");
    }

    //// Get hits collections
    auto CollectorHC = GetHitsCollection(fCollectorID, event);

    if ( CollectorHC->entries()>0 ) {
        // Get hit with total values
        //auto CollectorHit = (*CollectorHC)[CollectorHC->entries()-1];

        // Print per event
    //    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    //    if( ( printModulo > 1 ) && ( evtid % printModulo == 0 ) ) {
    //        G4cout << "---> End of event: " << evtid << G4endl;     
    //        PrintEventStatistics( CollectorHit->GetTrackID() );
    //    }
        analysis->analysePhotonNumber(CollectorHC->entries());
    }  else  {  // no hits in SD
        analysis->analysePhotonNumber(0);
    }

    analysis->analyseAddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
