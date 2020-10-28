#include "LArDetectorConstruction.hh"
#include "LArDetectorConstructionMessenger.hh"
#include "LArCollectorSD.hh"

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
#include "G4SDManager.hh"

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
    const G4int nEntries = 2;
    G4double PhotonEnergy[nEntries] = {1.*eV, 10.*eV};

    // material: air
    G4NistManager* nist = G4NistManager::Instance();
    fAir = nist->FindOrBuildMaterial("G4_AIR");
    G4MaterialPropertiesTable* air_mpt = new G4MaterialPropertiesTable();

    G4double airRindex[nEntries] = {1.0, 1.0};
    G4double airAbsLength[nEntries] = {10000*m, 10000*m};
    G4double airRayleigh[nEntries] = {10000*m, 10000*m};
    air_mpt -> AddProperty("RINDEX", PhotonEnergy, airRindex, nEntries);
    air_mpt -> AddProperty("ABSLENGTH", PhotonEnergy, airAbsLength, nEntries);
    air_mpt -> AddProperty("RAYLEIGH", PhotonEnergy, airRayleigh, nEntries);
    
    fAir -> SetMaterialPropertiesTable(air_mpt);

    // material: LAr
    G4double z, a, density;
    G4String name;
    density = 1.390*g/cm3;
    a = 39.95*g/mole;
    fLAr = new G4Material(name="liquidArgon", z=18., a, density);
    G4MaterialPropertiesTable* lar_mpt = new G4MaterialPropertiesTable();

    G4double larRindex[nEntries] = {1.4, 1.4};
    G4double larAbsLength[nEntries] = {1*m, 1*m};
    G4double larRayleigh[nEntries] = {1*m, 1*m};
    lar_mpt -> AddProperty("RINDEX", PhotonEnergy, larRindex, nEntries);
    lar_mpt -> AddProperty("ABSLENGTH", PhotonEnergy, larAbsLength, nEntries);
    lar_mpt -> AddProperty("RAYLEIGH", PhotonEnergy, larRayleigh, nEntries);
    fLAr -> SetMaterialPropertiesTable(lar_mpt);
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

    G4VisAttributes* cdVisAtt = new G4VisAttributes(G4Colour(0, 1, 1));
    cdVisAtt->SetForceSolid();
    cdLV->SetVisAttributes(cdVisAtt);

    /// photon collector
    G4double inner = radius;
    G4double outer = radius+1*mm;
    G4Sphere* colSD = 
    new G4Sphere(
            "colSD",
            inner, outer,
            0, twopi,
            0, pi
        );

    
    G4LogicalVolume* colLV = 
    new G4LogicalVolume(
            colSD,
            //fLAr,
            fAir,
            "colLV"
        );

    new G4PVPlacement(
            0,
            G4ThreeVector(),
            colLV,
            "colPV",
            worldLV,
            false,
            0,
            fCheckOverlaps
        );
    G4VisAttributes* colVisAtt = new G4VisAttributes(G4Colour(1, 1, 0));
    colLV->SetVisAttributes(colVisAtt);

    return worldPV;
}

void LArDetectorConstruction::ConstructSDandField()
{
    G4cout << " ----> Add Sensitive Detector "  << G4endl;

    auto collectorSD =
        new LArCollectorSD("colSD", "CollectorHitsCollection", 1);
    G4SDManager::GetSDMpointer()->AddNewDetector(collectorSD);
    SetSensitiveDetector("colLV", collectorSD);
}
