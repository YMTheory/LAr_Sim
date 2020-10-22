#ifndef LArParticleSourceMessenger_h
#define LArParticleSourceMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class LArParticleSource;

class G4ParticleTable;
class G4UIcommand;
class G4UImessenger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

class LArParticleSourceMessenger: public G4UImessenger {
    public:
        LArParticleSourceMessenger(LArParticleSource* fPtclGun);
        ~LArParticleSourceMessenger();

        void SetNewValue(G4UIcommand* cmd, G4String newValues);

    private:
        LArParticleSource* fParticleGun;
        G4ParticleTable* particleTable;

    private:
        G4UIdirectory                *gunDirectory;

        G4UIcmdWithAnInteger         *numberCmd;
        G4UIcmdWithAString           *particleCmd;
        G4UIcmdWithAString           *typeCmd;
        G4UIcmdWith3VectorAndUnit    *centreCmd;
        G4UIcmdWith3VectorAndUnit    *positionCmd;
        G4UIcmdWithADoubleAndUnit    *energyCmd;
        G4UIcmdWithAString           *energytypeCmd;
        G4UIcmdWithAnInteger         *verbosityCmd;
        G4UIcmdWithoutParameter      *listCmd;

};


#endif
