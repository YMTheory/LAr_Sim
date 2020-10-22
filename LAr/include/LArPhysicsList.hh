#ifndef LArPhysicsList_h
#define LArPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;

class LArPhysicsList : public G4VModularPhysicsList
{
    public:
        LArPhysicsList();
        virtual ~LArPhysicsList();

    public:
        virtual void SetCuts();

        virtual void ConstructParticle();
        virtual void ConstructProcess();

        void ConstructOpticalProcess();
};

#endif
