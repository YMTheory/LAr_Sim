#ifndef LArDetectorConstruction_h
#define LArDetectorConstruction_h

#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class  G4Material;

class LArDetectorConstruction : public G4VUserDetectorConstruction
{
    public:
        LArDetectorConstruction();
        virtual ~LArDetectorConstruction();

        virtual G4VPhysicalVolume* Construct();
        //virtual void ConstructSDandField();

        void DefineMaterials();
        G4VPhysicalVolume* DefineVolumes();

    private:
        G4bool fCheckOverlaps;

        G4Material* fAir;
        G4Material* fLAr;
};

#endif
