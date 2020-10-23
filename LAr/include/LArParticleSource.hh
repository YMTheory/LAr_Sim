#ifndef LArParticleSource_h
#define LArParticleSource_h 1


#include "G4VPrimaryGenerator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"

#include "LArParticleSourceMessenger.hh"

class LArParticleSource : public G4VPrimaryGenerator {

    public:
        LArParticleSource   ();
        ~LArParticleSource  ();
    
        void GeneratePrimaryVertex(G4Event* event);

    public:
        // particle number
        void SetParticleNumber(G4int);

        // position distribution
        void SetPosDisType(G4String);
        void SetCentreCoords(G4ThreeVector);
        void SetFVRadius(G4double);
        void GeneratePointSource();
        void GenerateUniformSource();

        // angular distribution
        void SetParticleMomentumDirection(G4ParticleMomentum);
        void GenerateIsotropicFlux();

        // energy distribution
        void SetEnergyDisType(G4String);
        void SetMonoEnergy(G4double);
        inline G4double GetParticleEnergy() { return particle_energy; }
        void GenerateMonoEnergetic();
        void GenerateEnergySpectrum();

        // verbosity
        void SetVerbosity(G4int);

        // particle properties
        void SetParticleDefinition(G4ParticleDefinition* aParticleDefinition);

        
    private:
        // position distribution
        G4String                SourcePosType;
        G4ThreeVector           CentreCoords;
        G4double                FV_radius;
        
        // energy distribution
        G4String EnergyDisType;
        G4double MonoEnergy;

        // particle properties
        G4int                   NumberOfParticlesToBeGenerated;
        G4ParticleDefinition*   particle_definition;
        G4ParticleMomentum      particle_momentum_direction;
        G4double                particle_energy;
        G4ThreeVector           particle_position;
        G4double                particle_time;
        G4ThreeVector           particle_polarization;
    
        //std::vector<G4double>   vec_particle_y;

        // Verbose
        G4int verbosityLevel;

    private:
        LArParticleSourceMessenger* theMessenger;
};

#endif
