#ifndef LArDetectorConstruction_h
#define LArDetectorConstruction_h

#include "G4VUserDetectorConstruction.hh"
#include "LArDetectorConstructionMessenger.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class  G4Material;

class LArDetectorConstruction : public G4VUserDetectorConstruction
{
    public:
        LArDetectorConstruction();
        virtual ~LArDetectorConstruction();

        virtual G4VPhysicalVolume* Construct();
        virtual void ConstructSDandField();

        void DefineMaterials();
        G4VPhysicalVolume* DefineVolumes();
    
    public:
        void setRadius(G4double fR) { radius = fR; }

    private:
        LArDetectorConstructionMessenger* fMessenger;

    private:
        G4bool fCheckOverlaps;

        G4Material* fAir;
        G4Material* fLAr;

        G4double radius;
};

#endif
