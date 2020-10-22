/*************************************************************************
 @Author: MiaoYu ---> miaoyu@ihep.ac.cn
 @Created Time : Tue Sep 29 16:11:41 2020
 @File Name: LArAnalysisManager.cc
 ************************************************************************/

#include "g4root.hh"

#include "LArAnalysisManager.hh"

LArAnalysisManager* LArAnalysisManager::instance = 0;

LArAnalysisManager::LArAnalysisManager()
    : outputFileName("user")
{
    G4AnalysisManager::Instance();

    // creating the messenger
    analysisMessenger = new LArAnalysisMessenger(this);

    G4cout << " +++++ LArAnalysisManager created" << G4endl;
}

LArAnalysisManager::~LArAnalysisManager()
{
    delete instance;
    instance = 0;

    delete G4AnalysisManager::Instance();
}

void LArAnalysisManager::book()
{
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    //Open an output file
    man->OpenFile(outputFileName);
    man->SetVerboseLevel(1);
    man->SetFirstHistoId(1);
    man->SetFirstNtupleId(1);

    G4cout << "Open output file: " << outputFileName << G4endl;
    man->CreateNtuple("photon", "Hits Info on SD");
    man->CreateNtupleIColumn("EventID");
    man->CreateNtupleIColumn("nPhoton");
    man->CreateNtupleIColumn("nInitPhoton");
    man->CreateNtupleDColumn("InitPosY");
    man->FinishNtuple();
    G4cout << "Created ntuple for photon counting" << G4endl;
}

void LArAnalysisManager::finish()
{
    G4cout << "Going to save ntuples" << G4endl;
    // Save histograms
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile();
}

LArAnalysisManager* LArAnalysisManager::getInstance()
{
    if (instance==0) { instance = new LArAnalysisManager(); }
    return instance;
}


void LArAnalysisManager::SetOutputFileName(G4String newName)
{
  
  outputFileName = newName;
}

void LArAnalysisManager::analyseEventID( G4int evtid )
{
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->FillNtupleIColumn( 0, evtid );
}


void LArAnalysisManager::analysePhotonNumber(G4int number)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->FillNtupleIColumn( 1, number );
}

void LArAnalysisManager::analyseInitPhotonNumber(G4int number)
{
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->FillNtupleIColumn( 2, number );
}


void LArAnalysisManager::analyseAddNtupleRow()
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man->AddNtupleRow();
}

