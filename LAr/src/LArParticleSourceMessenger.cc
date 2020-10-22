/*************************************************************************
 @Author: MiaoYu ---> miaoyu@ihep.ac.cn
 @Created Time : Mon Sep 28 10:48:06 2020
 @File Name: LArParticleSourceMessenger.cc
 ************************************************************************/

#include "LArParticleSourceMessenger.hh"
#include "LArParticleSource.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalPhoton.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"
#include "G4Tokenizer.hh"

LArParticleSourceMessenger::LArParticleSourceMessenger
(LArParticleSource* fPtclGun) 
    : fParticleGun(fPtclGun)
{
    gunDirectory = new G4UIdirectory("/LAr/gun/");
    gunDirectory->SetGuidance("Particle Source control commands");

    listCmd = new G4UIcmdWithoutParameter("/LAr/gun/List",this);
    listCmd->SetGuidance("List available particles");
    listCmd->SetGuidance("Invoke G4ParticleTable");

    // particle definition
    particleCmd = new G4UIcmdWithAString("/LAr/gun/particle", this);
    particleCmd->SetGuidance("Set particle to be generated.");
    particleCmd->SetGuidance(" (opticalphoton is default)");
    particleCmd->SetGuidance(" (ion can be specified for shooting ions)");
    particleCmd->SetParameterName("particleName",true);
    particleCmd->SetDefaultValue("opticalphoton");
    particleCmd->SetCandidates("opticalphoton");

    // particle number
    numberCmd = new G4UIcmdWithAnInteger("/LAr/gun/number", this);
    numberCmd->SetGuidance("Set initial particle number");
    numberCmd->SetParameterName("particleNumber", true);
    numberCmd->SetDefaultValue(10000);

   
    // particles energy
    energyCmd = new G4UIcmdWithADoubleAndUnit("/LAr/gun/energy",this);
    energyCmd->SetGuidance("Set kinetic energy.");
    energyCmd->SetParameterName("Energy",true,true);
    energyCmd->SetDefaultUnit("eV");
    //energyCmd->SetUnitCategory("Energy");
    //energyCmd->SetUnitCandidates("eV keV MeV GeV TeV");
    

    positionCmd = new G4UIcmdWith3VectorAndUnit("/LAr/gun/position",this);
    positionCmd->SetGuidance("Set starting position of the particle.");
    positionCmd->SetParameterName("X","Y","Z",true,true);
    positionCmd->SetDefaultUnit("cm");
    //positionCmd->SetUnitCategory("Length");
    //positionCmd->SetUnitCandidates("microm mm cm m km");

    // source position distribution type
    typeCmd = new G4UIcmdWithAString("/LAr/gun/type", this);
    typeCmd->SetGuidance("Sets source distribution type.");
    typeCmd->SetParameterName("DisType", true, true);
    typeCmd->SetDefaultValue("Point");
    typeCmd->SetCandidates("Point Uniform");


    // centre coordinates
    centreCmd = new G4UIcmdWith3VectorAndUnit("/LAr/gun/center", this);
    centreCmd->SetGuidance("Set centre coordinates of source.");
    centreCmd->SetParameterName("X","Y","Z",true,true);
    centreCmd->SetDefaultUnit("cm");
    centreCmd->SetUnitCandidates("nm um mm cm m km");

    // energy distribution
    energytypeCmd = new G4UIcmdWithAString("/LAr/gun/energytype", this);
    energytypeCmd->SetGuidance("Sets energy distribution type");
    energytypeCmd->SetGuidance("Possible variables are: Mono");
    energytypeCmd->SetParameterName("EnergyDis",true,true);
    energytypeCmd->SetDefaultValue("Mono");
    energytypeCmd->SetCandidates("Mono");

    
    // verbosity
    verbosityCmd = new G4UIcmdWithAnInteger("/LAr/gun/verbose",this);
    verbosityCmd->SetGuidance("Set Verbose level for gun");
    verbosityCmd->SetGuidance(" 0 : Silent");
    verbosityCmd->SetGuidance(" 1 : Limited information");
    verbosityCmd->SetGuidance(" 2 : Detailed information");
    verbosityCmd->SetParameterName("level",false);
    verbosityCmd->SetRange("level>=0 && level <=2");

}

LArParticleSourceMessenger::~LArParticleSourceMessenger()
{
    delete listCmd;
    delete typeCmd;
    delete centreCmd;
    delete energytypeCmd;
    delete particleCmd;
    delete positionCmd;
    delete energyCmd;
    delete numberCmd;

    delete gunDirectory;
}

void LArParticleSourceMessenger::SetNewValue
(G4UIcommand* cmd, G4String newValues) {
    
    
    if ( cmd == typeCmd )
        fParticleGun->SetPosDisType(newValues);

    else if ( cmd == centreCmd ) 
        fParticleGun->SetCentreCoords( centreCmd->GetNew3VectorValue(newValues) );

    else if ( cmd == energytypeCmd )
        fParticleGun->SetEnergyDisType(newValues);

    else if ( cmd == particleCmd ) {
        G4ParticleDefinition *pd = particleTable->FindParticle(newValues);
        if( pd != NULL )
            fParticleGun->SetParticleDefinition(pd);
    }

    else if ( cmd == listCmd )
        particleTable->DumpTable();

    else if ( cmd == energyCmd ) {
        fParticleGun->SetEnergyDisType("Mono");
        fParticleGun->SetMonoEnergy(energyCmd->GetNewDoubleValue(newValues));
    }

    else if ( cmd == positionCmd ) {
        fParticleGun->SetPosDisType("Point");
        fParticleGun->SetCentreCoords(positionCmd->GetNew3VectorValue(newValues));
    }

    else if ( cmd == numberCmd ) {
        fParticleGun->SetParticleNumber(numberCmd->GetNewIntValue(newValues));
    }

    else if ( cmd == verbosityCmd )
        fParticleGun->SetVerbosity(verbosityCmd->GetNewIntValue(newValues));

    else
        G4cout << "Error: Unknow Command" << G4endl;
        
}
