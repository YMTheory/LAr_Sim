/*************************************************************************
 @Author: MiaoYu ---> miaoyu@ihep.ac.cn
 @Created Time : Fri Sep 25 16:13:44 2020
 @File Name: LArCollectorSD.cc
 ************************************************************************/

#include "LArCollectorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

LArCollectorSD::LArCollectorSD( const G4String& name, 
                  const G4String& hitsCollectionName,
                  G4int nofPmts)
    : G4VSensitiveDetector(name),
    fHitsCollection(NULL),
    fNofPmts(nofPmts)
{
    collectionName.insert(hitsCollectionName);
}

LArCollectorSD::~LArCollectorSD()
{;}


void LArCollectorSD::Initialize(G4HCofThisEvent* hce)
{
    // Create hits collection
    fHitsCollection
        = new LArCollectorHitsCollection( SensitiveDetectorName, collectionName[0]);
    
    // Add this collection in hce
    G4int hcID 
        = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection( hcID, fHitsCollection );

    // Create hits
    //for( G4int i=0; i<fNofPmts+1; i++ )  {
    //    fHitsCollection->insert(new LArCollectorHit());
    //}
}

G4bool LArCollectorSD::ProcessHits( G4Step* aStep, G4TouchableHistory*)
{
    G4double edep = aStep->GetTotalEnergyDeposit();
    G4double stepLength = aStep->GetStepLength();
    if(edep == 0. && stepLength == 0. ) return false;

    auto touchable = (aStep->GetPreStepPoint()->GetTouchable());

    // Get pmt id
    auto pmtNumber = touchable->GetReplicaNumber(1);

    //auto hit = (*fHitsCollection)[pmtNumber];
    //if ( ! hit ) {
    //    G4ExceptionDescription msg;
    //    msg << "Cannot access hit " << pmtNumber; 
    //    G4Exception("B4cCalorimeterSD::ProcessHits()",
    //            "MyCode0004", FatalException, msg);
    //}         

    LArCollectorHit* hit = new LArCollectorHit();
    G4int trackId = aStep->GetTrack()->GetTrackID();
    hit->SetTrackID(trackId);

    fHitsCollection->insert(hit);

    return true;
}



void LArCollectorSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}
