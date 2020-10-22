#include "LArDetectorConstruction.hh"
#include "LArDetectorConstructionMessenger.hh"

#include "G4NistManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

LArDetectorConstruction::LArDetectorConstruction()
    : G4VUserDetectorConstruction(),
    fCheckOverlaps(true),
    fAir(NULL),
    fLAr(NULL)
{
    radius = 85*cm;    

    fMessenger  = new LArDetectorConstructionMessenger(this);
}


LArDetectorConstruction::~LArDetectorConstruction()
{}

G4VPhysicalVolume* LArDetectorConstruction::Construct()
{
    DefineMaterials();
    
    return DefineVolumes();
}

void LArDetectorConstruction::DefineMaterials()
{
    // material: air
    G4NistManager* nist = G4NistManager::Instance();
    fAir = nist->FindOrBuildMaterial("G4_AIR");
    G4MaterialPropertiesTable* air_mpt = new G4MaterialPropertiesTable();
    air_mpt -> AddConstProperty("RINDEX", 1.0);
    fAir -> SetMaterialPropertiesTable(air_mpt);

    // material: LAr
    G4double z, a, density;
    G4String name;
    density = 1.390*g/cm3;
    a = 39.95*g/mole;
    fLAr = new G4Material(name="liquidArgon", z=18., a, density);

}


G4VPhysicalVolume* LArDetectorConstruction::DefineVolumes()
{
    G4cout << "Initializa CD : radius = " << radius << G4endl;

    G4double worldLength = 3 * radius;

    ///// World Construction /////
    G4Box* worldSD = 
        new G4Box("worldSD", 
                 worldLength/2., worldLength/2., worldLength/2.);
    G4LogicalVolume* worldLV = 
    new G4LogicalVolume(
            worldSD,
            fAir,
            "worldLV"
        );
    G4VPhysicalVolume* worldPV =
    new G4PVPlacement(
            0,
            G4ThreeVector(),
            worldLV,
            "worldPV",
            0,
            false,
            0,
            fCheckOverlaps
        );

    
    ///// Spheric CD /////
    G4Sphere* cdSD = 
        new G4Sphere("cdSD", 
                    0, radius,
                    0, twopi,
                    0, pi
                    );
    G4LogicalVolume* cdLV =
    new G4LogicalVolume(
            cdSD,
            fLAr,
            "cdLV"
        );

    new G4PVPlacement(
            0,
            G4ThreeVector(),
            cdLV,
            "cdPV",
            worldLV,
            false,
            0,
            fCheckOverlaps
        );

    return worldPV;
}
