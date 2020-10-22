#include "LArParticleSource.hh"
//#include "LArAnalysisManager.hh"

#include "G4Event.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

LArParticleSource::LArParticleSource ()  {

    NumberOfParticlesToBeGenerated = 10000; 
    particle_definition            = NULL;
    G4ThreeVector zero (0., 0., 0.);
    particle_momentum_direction    = G4ParticleMomentum(1., 0., 0.);
    particle_time                  = 0.0;
    particle_energy                = 1.0*MeV;
    particle_polarization          = G4ThreeVector(1., 0., 0.);

    SourcePosType                  = "Point";
    CentreCoords                   = zero;

    EnergyDisType                  = "Mono";
    MonoEnergy                     = 1*MeV;

    FV_radius                      = 85*cm;
    //theMessenger                   = new LArParticleSourceMessenger(this);
}

LArParticleSource::~LArParticleSource () {
    //delete theMessenger;
}


// set particle number in one event
void LArParticleSource::SetParticleNumber(G4int number) 
{
    NumberOfParticlesToBeGenerated = number;
}

// set particle Position distribution type
void LArParticleSource::SetPosDisType(G4String PosType)
{
    // candidates: Point/Uniform
    SourcePosType = PosType;
}

// set particle PointType position
void LArParticleSource::SetCentreCoords(G4ThreeVector coordsOfCentre)
{
    CentreCoords = coordsOfCentre;
}

// generate point source
void LArParticleSource::GeneratePointSource()
{
    if(SourcePosType == "Point" ) 
        particle_position = CentreCoords;
    else{
        if(verbosityLevel >= 1)
            G4cout << "Error SourcePosType is not set to Point" << G4endl;
    }
}

// generate uniform source
void LArParticleSource::GenerateUniformSource()
{
    if (SourcePosType == "Uniform" ) {
        G4double phi = G4UniformRand()*2*pi;
        G4double costheta = -1 * 2*G4UniformRand();
        G4double u = G4UniformRand();

        G4double posr = FV_radius * std::cbrt(u);
        G4double posx = posr * std::sin(phi) * costheta;
        G4double posy = posr * std::sin(phi) * (std::sqrt(1-costheta*costheta));
        G4double posz = posr * std::cos(phi);

        particle_position = G4ThreeVector(posx, posy, posz);
    } else {
        if(verbosityLevel >= 1)
            G4cout << "Error SourcePosType is not set to Uniform" << G4endl;
    }

}

void LArParticleSource::SetParticleMomentumDirection
(G4ParticleMomentum aDirection)  {
    particle_momentum_direction = aDirection.unit();
}

void LArParticleSource::GenerateIsotropicFlux()
{
    G4double rndm, rndm2;
    G4double px, py, pz;

    G4double sintheta, sinphi, costheta, cosphi;
    rndm = G4UniformRand();
    costheta = -1 + rndm * 2;
    sintheta = std::sqrt(1. - costheta*costheta);

    rndm2 = G4UniformRand();
    G4double Phi = 2 * pi * rndm2; 
    sinphi = std::sin(Phi);
    cosphi = std::cos(Phi);

    px = -sintheta * cosphi;
    py = -sintheta * sinphi;
    pz = -costheta;

    G4double ResMag = std::sqrt((px*px) + (py*py) + (pz*pz));
    px = px/ResMag;
    py = py/ResMag;
    pz = pz/ResMag;

    particle_momentum_direction.setX(px);
    particle_momentum_direction.setY(py);
    particle_momentum_direction.setZ(pz);

    // particle_momentum_direction now holds unit momentum vector.
    if(verbosityLevel >= 2)
        G4cout << "Generating isotropic vector: " << particle_momentum_direction << G4endl;
}


// set particle energy type
void LArParticleSource::SetEnergyDisType(G4String DisType)
{
    // candidates: Mono, Spec
    EnergyDisType = DisType;
}

void LArParticleSource::SetMonoEnergy(G4double menergy) 
{
    MonoEnergy = menergy;
}

// mono energy
void LArParticleSource::GenerateMonoEnergetic()
{
    particle_energy = MonoEnergy;
}

void LArParticleSource::SetVerbosity(G4int vL)
{
    verbosityLevel = vL;
    G4cout << " Verbosity Level Set to : " << verbosityLevel << G4endl;
}


void LArParticleSource::SetParticleDefinition(G4ParticleDefinition* aParticleDefinition)
{
    particle_definition = aParticleDefinition;
}


void LArParticleSource::GeneratePrimaryVertex(G4Event* event)  {
    if( particle_definition == NULL ) {
        G4cout << "No particle has been defined !" << G4endl;
        return;
    }


    //LArAnalysisManager* analysis = LArAnalysisManager::getInstance();
    //analysis->analyseInitPhotonNumber(NumberOfParticlesToBeGenerated);

    for( G4int i=0; i<NumberOfParticlesToBeGenerated; i++ ) {
        // Position
        if( SourcePosType == "Point" ) {
            GeneratePointSource();
        } else if( SourcePosType == "Uniform" ) {  
            GenerateUniformSource();
        } else {
            G4cout << "Error:: SourcePosType undefined" << G4endl;
            G4cout << "Generating point source" << G4endl;
            GeneratePointSource();
        }

        GenerateIsotropicFlux();

        // Energy Stuff
        if(EnergyDisType == "Mono" )
            GenerateMonoEnergetic();
        else
            G4cout << "Error: EnergyDisType has unusual value" << G4endl;


        // create a new verbosityLevel
        G4PrimaryVertex* vertex = new G4PrimaryVertex(particle_position, particle_time);
        if(verbosityLevel >= 2)
            G4cout << "Creating primaries and assigning to vertex" << G4endl;
        G4double mass   = 0.;   // optical photon
        G4double energy = particle_energy;
        G4double pmom   = std::sqrt(energy*energy - mass*mass);
        G4double px     = pmom*particle_momentum_direction.x();
        G4double py     = pmom*particle_momentum_direction.y();
        G4double pz     = pmom*particle_momentum_direction.z();
        if(verbosityLevel >= 2){
            G4cout << "Particle name: " 
                << particle_definition->GetParticleName() << G4endl; 
            G4cout << "       Energy: "<<particle_energy << G4endl;
            G4cout << "     Position: "<<particle_position<< G4endl; 
            G4cout << "    Direction: "<<particle_momentum_direction << G4endl;
            G4cout << " NumberOfParticlesToBeGenerated: "
                << NumberOfParticlesToBeGenerated << G4endl;
        }

        G4PrimaryParticle* particle =
            new G4PrimaryParticle(particle_definition,px,py,pz);
        particle->SetMass( mass );
        //particle->SetCharge( particle_charge );
        particle->SetPolarization(particle_polarization.x(),
                particle_polarization.y(),
                particle_polarization.z());
        vertex->SetPrimary( particle );
        event->AddPrimaryVertex( vertex );
    }
    if(verbosityLevel > 1)
        G4cout << " Primary Vetex generated "<< G4endl;   


}
