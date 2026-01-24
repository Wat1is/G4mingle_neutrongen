//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4dDetectorConstruction.cc
/// \brief Implementation of the B4dDetectorConstruction class

#include "ModelDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4AutoDelete.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4SDChargedFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDKineticEnergyFilter.hh"
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSPopulation.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4PSPassageCellCurrent.hh"
#include "G4PSCellFlux.hh"
#include "G4PSTrackLength.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSFlatSurfaceFlux.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSNofStep.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "F05Field.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4EqMagElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4DormandPrince745.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4MagIntegratorDriver.hh"

#include "GB01BOptrMultiParticleChangeCrossSection.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModelDetectorConstruction::ModelDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModelDetectorConstruction::~ModelDetectorConstruction()
{ 
	if (fEMfield) delete fEMfield;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ModelDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModelDetectorConstruction::DefineMaterials()
{ 
	// Lead material defined using NIST Manager
	G4NistManager* man = G4NistManager::Instance();

	// define NIST materials

	G4Isotope* H2 = new G4Isotope("H2", 1, 2, 2.01410177812*g / mole);
	G4Isotope* H3 = new G4Isotope("H3", 1, 3, 3.0160492779*g / mole);

	G4Element* H = man->FindOrBuildElement("H");
	G4Element* O = man->FindOrBuildElement("O");
	G4Element* C = man->FindOrBuildElement("C");
	G4Element* Ti = man->FindOrBuildElement("Ti");
	G4Element* Al = man->FindOrBuildElement("Al");
	

	G4Element* D = new G4Element("Deuterium", "D", 1);
	D->AddIsotope(H2, 100.*perCent);
	G4Element* T = new G4Element("Tritium", "T", 1);
	T->AddIsotope(H3, 100.*perCent);

	G4Material* air = man->FindOrBuildMaterial("G4_AIR");
	G4Material* water = man->FindOrBuildMaterial("G4_WATER");
	G4Material* vacuum = man->FindOrBuildMaterial("G4_Galactic");
	G4Material* lead = man->FindOrBuildMaterial("G4_Pb");
	G4Material* aluminium = man->FindOrBuildMaterial("G4_Al");

	/*G4Material* tubeair = man->BuildMaterialWithNewDensity("TubeAir", "G4_AIR", 0.000000012*g/cm3);*/

	G4Material* tar = new G4Material("TargetMaterial", 0.001*g/cm3, 2);
	tar->AddElement(T, 0.6);
	tar->AddElement(D, 0.4);

	G4Material* DTgas = new G4Material("DeuteriumTritiumGas", 0.00000000096*g / cm3, 2);
	DTgas->AddElement(T, 0.6);
	DTgas->AddElement(D, 0.4);

	G4Material* ethylenglycol = new G4Material("EthylenGlycol", 1.114000*g / cm3, 3);
	ethylenglycol->AddElement(H, 6);
	ethylenglycol->AddElement(C, 2);
	ethylenglycol->AddElement(O, 2);

	G4Material* coolingwater = new G4Material("CoolingWater", 1.13*g / cm3, 2);
	coolingwater->AddMaterial(water, 0.7);
	coolingwater->AddMaterial(ethylenglycol, 0.3);

	G4Material* coolingmat = new G4Material("CoolingMat", 1.62782*g / cm3, 2);
	coolingmat->AddMaterial(coolingwater, 0.4740675);
	coolingmat->AddMaterial(aluminium, 0.5259325);

	G4Material* stilbene = new G4Material("Stilbene", 1.16*g / cm3, 2);
	stilbene->AddElement(C, 14);
	stilbene->AddElement(H, 12);

	G4Material* al2o3 = new G4Material("Al2O3", 3.85*g / cm3, 2);
	al2o3->AddElement(Al, 2);
	al2o3->AddElement(O, 3);

	G4Material* gettermat = new G4Material("GetterMat", 3.76*g / cm3, 3);
	gettermat->AddElement(T, 1);
	gettermat->AddElement(D, 1);
	gettermat->AddElement(Ti, 1);

	/*G4Material* steel = new G4Material("Steel", 7.87*g / cm3, 4);
	steel->AddElement(Fe, 99.6 * perCent);
	steel->AddElement(C, 0.31* perCent);
	steel->AddElement(S, 0.050* perCent);
	steel->AddElement(P, 0.04* perCent);*/

	// Print materials
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ModelDetectorConstruction::DefineVolumes()
{
	// Get materials
	/*auto water = G4Material::GetMaterial("G4_WATER");*/

	auto air = G4Material::GetMaterial("G4_AIR");
	auto vacuum = G4Material::GetMaterial("G4_Galactic");
	auto tar = G4Material::GetMaterial("TargetMaterial");
	auto lead = G4Material::GetMaterial("G4_Pb");
	auto aluminium = G4Material::GetMaterial("G4_Al");
	auto coolingmat = G4Material::GetMaterial("CoolingMat");
	auto dtgas = G4Material::GetMaterial("DeuteriumTritiumGas");
	auto gettermat = G4Material::GetMaterial("GetterMat");
	auto al2o3 = G4Material::GetMaterial("Al2O3");
	auto stilbene = G4Material::GetMaterial("Stilbene");
	



	//     
	// World
	//
	auto worldS
		= new G4Box("World",           // its name
			200.*cm, 200.*cm, 200.*cm); // its size
	auto worldLV
		= new G4LogicalVolume(
			worldS,           // its solid
			air,  // its material
			"World");         // its name
	auto worldPV
		=new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(),  // at (0,0,0)
			worldLV,          // its logical volume                         
			"World",          // its name
			0,                // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps 

	auto shieldingS
		= new G4Tubs("Shielding",     // its name
			78.*mm, 89. * mm, 346. * mm, 0., 360.); // its size
	auto shieldingLV
		= new G4LogicalVolume(
			shieldingS,    // its solid
			lead, // its material
			"Shielding");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, -48.*mm),
			shieldingLV,          // its logical volume                         
			"Shielding",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto extercasingS
		= new G4Tubs("ExterCasing",     // its name
			75.*mm, 78. * mm, 346. * mm, 0., 360.); // its size
	auto extercasingLV
		= new G4LogicalVolume(
			extercasingS,    // its solid
			aluminium, // its material
			"ExterCasing");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, -48.*mm),
			extercasingLV,          // its logical volume                         
			"ExterCasing",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto topcasingS
		= new G4Tubs("TopCasing",     // its name
			0.*mm, 75. * mm, 7.5 * mm, 0., 360.); // its size
	auto topcasingLV
		= new G4LogicalVolume(
			topcasingS,    // its solid
			aluminium, // its material
			"TopCasing");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, 290.5*mm),
			topcasingLV,          // its logical volume                         
			"TopCasing",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto bottomcasingS
		= new G4Tubs("BottomCasing",     // its name
			58.*mm, 75. * mm, 7.5 * mm, 0., 360.); // its size
	auto bottomcasingLV
		= new G4LogicalVolume(
			bottomcasingS,    // its solid
			aluminium, // its material
			"BottomCasing");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, -386.5*mm),
			bottomcasingLV,          // its logical volume                         
			"BottomCasing",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto getterS
		= new G4Tubs("Getter",     // its name
			0.*mm, 58. * mm, 10 * mm, 0., 360.); // its size
	auto getterLV
		= new G4LogicalVolume(
			getterS,    // its solid
			gettermat, // its material
			"Getter");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, -384.*mm),
			getterLV,          // its logical volume                         
			"Getter",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto coolingS
		= new G4Tubs("Cooling",     // its name
			58.*mm, 75. * mm, 331 * mm, 0., 360.); // its size
	auto coolingLV
		= new G4LogicalVolume(
			coolingS,    // its solid
			coolingmat, // its material
			"Cooling");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, -48.*mm),
			coolingLV,          // its logical volume                         
			"Cooling",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto innercasingS
		= new G4Tubs("InnerCasing",     // its name
			48.*mm, 58. * mm, 328.5 * mm, 0., 360.); // its size
	auto innercasingLV
		= new G4LogicalVolume(
			innercasingS,    // its solid
			aluminium, // its material
			"InnerCasing");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, -45.5*mm),
			innercasingLV,          // its logical volume                         
			"InnerCasing",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto reactchambertopS
		= new G4Tubs("ReactionChamberTop",     // its name
			0., 48. * mm, 54. * mm, 0., 360.); // its size
	auto reactchambertopLV
		= new G4LogicalVolume(
			reactchambertopS,    // its solid
			dtgas, // its material
			"ReactionChamberTop");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, 229.*mm),
			reactchambertopLV,          // its logical volume                         
			"ReactionChamberTop",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto insulatortopS
		= new G4Tubs("InsulatorTop",     // its name
			0., 26.5 * mm, 54. * mm, 0., 360.); // its size
	auto insulatortopLV
		= new G4LogicalVolume(
			insulatortopS,    // its solid
			al2o3, // its material
			"InsulatorTop");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, 0.*mm),
			insulatortopLV,          // its logical volume                         
			"InsulatorTop",    // its name
			reactchambertopLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto reactchamberbottomS
		= new G4Tubs("ReactionChamberBottom",     // its name
			0., 48. * mm, 99.5 * mm, 0., 360.); // its size
	auto reactchamberbottomLV
		= new G4LogicalVolume(
			reactchamberbottomS,    // its solid
			dtgas, // its material
			"ReactionChamberBottom");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, -274.5*mm),
			reactchamberbottomLV,          // its logical volume                         
			"ReactionChamberBottom",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto insulatorbottomS
		= new G4Tubs("InsulatorBottom",     // its name
			0., 26.5 * mm, 88. * mm, 0., 360.); // its size
	auto insulatorbottomLV
		= new G4LogicalVolume(
			insulatorbottomS,    // its solid
			al2o3, // its material
			"InsulatorBottom");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, 11.5*mm),
			insulatorbottomLV,          // its logical volume                         
			"InsulatorBottom",    // its name
			reactchamberbottomLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto getterbottomS
		= new G4Tubs("GetterBottom",     // its name
			0., 29.5 * mm, 11.5 * mm, 0., 360.); // its size
	auto getterbottomLV
		= new G4LogicalVolume(
			getterbottomS,    // its solid
			gettermat, // its material
			"GetterBottom");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, -88*mm),
			getterbottomLV,          // its logical volume                         
			"GetterBottom",    // its name
			reactchamberbottomLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps


	auto reactchamberS
		= new G4Tubs("ReactionChamber",     // its name
			0., 48. * mm, 175. * mm, 0., 360.); // its size
	auto reactchamberLV
		= new G4LogicalVolume(
			reactchamberS,    // its solid
			dtgas, // its material
			"ReactionChamber");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, 0.*cm),
			reactchamberLV,          // its logical volume                         
			"ReactionChamber",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto targetS
		= new G4Tubs("Target",     // its name
			0., 10. * mm, 175. * mm, 0., 360.); // its size
	auto targetLV
		= new G4LogicalVolume(
			targetS,    // its solid
			tar, // its material
			"Target");  // its name
	new G4PVPlacement(
			0,                // rotation
			G4ThreeVector(0.*cm, 0.*cm, 0.*cm),
			targetLV,          // its logical volume                         
			"Target",    // its name
			reactchamberLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto stilbenedetector1S=new G4Tubs("StilbeneDetector1",0, 25.*mm, 25.*mm, 0.*deg, 360.*deg); 
	auto stilbenedetector1LV=new G4LogicalVolume(stilbenedetector1S,vacuum,"StilbeneDetector1");
	G4RotationMatrix* xRot = new G4RotationMatrix;
	xRot->rotateX(90.*deg);
	new G4PVPlacement(xRot,G4ThreeVector(0.*mm, 1026.*mm, 0.*mm),stilbenedetector1LV,
			"StilbeneDetector1",worldLV,false,0,fCheckOverlaps);

	auto stilbenedetector2S=new G4Tubs("StilbeneDetector2",0, 25.*mm, 25.*mm, 0.*deg, 360.*deg); 
	auto stilbenedetector2LV=new G4LogicalVolume(stilbenedetector2S,stilbene,"StilbeneDetector2");
	new G4PVPlacement(xRot,G4ThreeVector(0.*mm, -1026.*mm, 0.*mm),stilbenedetector2LV,
			"StilbeneDetector2",worldLV,false,0,fCheckOverlaps);


	G4double R = 1000*mm;
	G4double Ro = R + 0.001;

	auto detectorR0S
		= new G4Tubs("DetectorR0",     // its name
			R, Ro, 25.*mm, 0.*deg, 360.*deg); // its size
	auto detectorR0LV
		= new G4LogicalVolume(
			detectorR0S,    // its solid
			vacuum, // its material
			"DetectorR0");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, 0.*mm),  // at (0,0,0)
			detectorR0LV,          // its logical volume                         
			"DetectorR0",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorR1upS
		= new G4Tubs("DetectorR1Up",     // its name
			R, Ro, 25.*mm, 0.*deg, 360.*deg); // its size
	auto detectorR1upLV
		= new G4LogicalVolume(
			detectorR1upS,    // its solid
			vacuum, // its material
			"DetectorR1Up");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, 50.*mm),  // at (0,0,0)
			detectorR1upLV,          // its logical volume                         
			"DetectorR1Up",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorR2upS
		= new G4Tubs("DetectorR2Up",     // its name
			R, Ro, 25.*mm, 0.*deg, 360.*deg); // its size
	auto detectorR2upLV
		= new G4LogicalVolume(
			detectorR2upS,    // its solid
			vacuum, // its material
			"DetectorR2Up");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, 100.*mm),  // at (0,0,0)
			detectorR2upLV,          // its logical volume                         
			"DetectorR2Up",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorR3upS
		= new G4Tubs("DetectorR3Up",     // its name
			R, Ro, 25.*mm, 0.*deg, 360.*deg); // its size
	auto detectorR3upLV
		= new G4LogicalVolume(
			detectorR3upS,    // its solid
			vacuum, // its material
			"DetectorR3Up");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, 150.*mm),  // at (0,0,0)
			detectorR3upLV,          // its logical volume                         
			"DetectorR3Up",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorR1lowS
		= new G4Tubs("DetectorR1Low",     // its name
			R, Ro, 25.*mm, 0.*deg, 360.*deg); // its size
	auto detectorR1lowLV
		= new G4LogicalVolume(
			detectorR1lowS,    // its solid
			vacuum, // its material
			"DetectorR1Low");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, -50.*mm),  // at (0,0,0)
			detectorR1lowLV,          // its logical volume                         
			"DetectorR1Low",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorR2lowS
		= new G4Tubs("DetectorR2Low",     // its name
			R, Ro, 25.*mm, 0.*deg, 360.*deg); // its size
	auto detectorR2lowLV
		= new G4LogicalVolume(
			detectorR2lowS,    // its solid
			vacuum, // its material
			"DetectorR2Low");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, -100.*mm),  // at (0,0,0)
			detectorR2lowLV,          // its logical volume                         
			"DetectorR2Low",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorR3lowS
		= new G4Tubs("DetectorR3Low",     // its name
			R, Ro, 25.*mm, 0.*deg, 360.*deg); // its size
	auto detectorR3lowLV
		= new G4LogicalVolume(
			detectorR3lowS,    // its solid
			vacuum, // its material
			"DetectorR3Low");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, -150.*mm),  // at (0,0,0)
			detectorR3lowLV,          // its logical volume                         
			"DetectorR3Low",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps


	auto detectorA0upS
		= new G4Tubs("DetectorA0Up",     // its name
			0.*mm, 25.*mm, 0.0005*mm, 0.*deg, 360.*deg); // its size
	auto detectorA0upLV
		= new G4LogicalVolume(
			detectorA0upS,    // its solid
			vacuum, // its material
			"DetectorA0Up");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, 798.0005*mm),  // at (0,0,0)
			detectorA0upLV,          // its logical volume                         
			"DetectorA0Up",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorA1upS
		= new G4Tubs("DetectorA1Up",     // its name
			25.*mm, 50.*mm, 0.0005*mm, 0.*deg, 360.*deg); // its size
	auto detectorA1upLV
		= new G4LogicalVolume(
			detectorA1upS,    // its solid
			vacuum, // its material
			"DetectorA1Up");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, 798.0005*mm),  // at (0,0,0)
			detectorA1upLV,          // its logical volume                         
			"DetectorA1Up",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorA2upS
		= new G4Tubs("DetectorA2Up",     // its name
			50.*mm, 75.*mm, 0.0005*mm, 0.*deg, 360.*deg); // its size
	auto detectorA2upLV
		= new G4LogicalVolume(
			detectorA2upS,    // its solid
			vacuum, // its material
			"DetectorA2Up");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, 798.0005*mm),  // at (0,0,0)
			detectorA2upLV,          // its logical volume                         
			"DetectorA2Up",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorA3upS
		= new G4Tubs("DetectorA3Up",     // its name
			75.*mm, 100.*mm, 0.0005*mm, 0.*deg, 360.*deg); // its size
	auto detectorA3upLV
		= new G4LogicalVolume(
			detectorA3upS,    // its solid
			vacuum, // its material
			"DetectorA3Up");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, 798.0005*mm),  // at (0,0,0)
			detectorA3upLV,          // its logical volume                         
			"DetectorA3Up",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps


	auto detectorA0lowS
		= new G4Tubs("DetectorA0Low",     // its name
			0.*mm, 25.*mm, 0.0005*mm, 0.*deg, 360.*deg); // its size
	auto detectorA0lowLV
		= new G4LogicalVolume(
			detectorA0lowS,    // its solid
			vacuum, // its material
			"DetectorA0Low");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, -894.0005*mm),  // at (0,0,0)
			detectorA0lowLV,          // its logical volume                         
			"DetectorA0Low",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorA1lowS
		= new G4Tubs("DetectorA1Low",     // its name
			25.*mm, 50.*mm, 0.0005*mm, 0.*deg, 360.*deg); // its size
	auto detectorA1lowLV
		= new G4LogicalVolume(
			detectorA1lowS,    // its solid
			vacuum, // its material
			"DetectorA1Low");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, -894.0005*mm),  // at (0,0,0)
			detectorA1lowLV,          // its logical volume                         
			"DetectorA1Low",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorA2lowS
		= new G4Tubs("DetectorA2Low",     // its name
			50.*mm, 75.*mm, 0.0005*mm, 0.*deg, 360.*deg); // its size
	auto detectorA2lowLV
		= new G4LogicalVolume(
			detectorA2lowS,    // its solid
			vacuum, // its material
			"DetectorA2Low");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, -894.0005*mm),  // at (0,0,0)
			detectorA2lowLV,          // its logical volume                         
			"DetectorA2Low",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

	auto detectorA3lowS
		= new G4Tubs("DetectorA3Low",     // its name
			75.*mm, 100.*mm, 0.0005*mm, 0.*deg, 360.*deg); // its size
	auto detectorA3lowLV
		= new G4LogicalVolume(
			detectorA3lowS,    // its solid
			vacuum, // its material
			"DetectorA3Low");  // its name
	new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0.*mm, 0.*mm, -894.0005*mm),  // at (0,0,0)
			detectorA3lowLV,          // its logical volume                         
			"DetectorA3Low",    // its name
			worldLV,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			fCheckOverlaps);  // checking overlaps

  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  /*auto yellow = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
  yellow->SetVisibility(true);
  yellow->SetForceWireframe(true);
  auto white = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  white->SetVisibility(true);
  white->SetForceWireframe(true);
  auto red = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  red->SetVisibility(true);
  red->SetForceWireframe(true);
  auto blue = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  blue->SetVisibility(true);
  blue->SetForceWireframe(true);
  auto gray = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  gray->SetVisibility(true);
  gray->SetForceWireframe(true);
  auto green = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
  green->SetVisibility(true);
  green->SetForceWireframe(true);
  auto magenta = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
  magenta->SetVisibility(true);
  magenta->SetForceWireframe(true);
  auto black = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0));
  black->SetVisibility(true);
  black->SetForceSolid(true);

	  plasticfrontLV->SetVisAttributes(white);*/


  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreadLocal F05Field* ModelDetectorConstruction::fEMfield = 0;

void ModelDetectorConstruction::ConstructSDandField()
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SDman->SetVerboseLevel(1);

  auto detector1 = new G4MultiFunctionalDetector("Detector01");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector1);
  auto detector2 = new G4MultiFunctionalDetector("Detector02");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector2);
  auto detector3 = new G4MultiFunctionalDetector("Detector03");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector3);
  auto detector4 = new G4MultiFunctionalDetector("Detector04");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector4);
  auto detector5 = new G4MultiFunctionalDetector("Detector05");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector5);
  auto detector6 = new G4MultiFunctionalDetector("Detector06");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector6);
  auto detector7 = new G4MultiFunctionalDetector("Detector07");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector7);
  auto detector8 = new G4MultiFunctionalDetector("Detector08");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector8);
  auto detector9 = new G4MultiFunctionalDetector("Detector09");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector9);
  auto detector10 = new G4MultiFunctionalDetector("Detector10");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector10);
  auto detector11 = new G4MultiFunctionalDetector("Detector11");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector11);
  auto detector12 = new G4MultiFunctionalDetector("Detector12");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector12);
  auto detector13 = new G4MultiFunctionalDetector("Detector13");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector13);
  auto detector14 = new G4MultiFunctionalDetector("Detector14");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector14);
  auto detector15 = new G4MultiFunctionalDetector("Detector15");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector15);
  auto detector16 = new G4MultiFunctionalDetector("Detector16");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector16);
  auto detector17 = new G4MultiFunctionalDetector("Detector17");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector17);
  auto detector18 = new G4MultiFunctionalDetector("Detector18");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector18);
  auto detector19 = new G4MultiFunctionalDetector("Detector19");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector19);
  auto detector20 = new G4MultiFunctionalDetector("Detector20");
  G4SDManager::GetSDMpointer()->AddNewDetector(detector20);

	auto neutronfilter = new G4SDParticleFilter("NeutronFilter", "neutron");
	auto protonfilter = new G4SDParticleFilter("ProtonFilter", "proton");
	auto deuteronfilter = new G4SDParticleFilter("DeuteronFilter", "deuteron");
	auto tritonfilter = new G4SDParticleFilter("TritonFilter", "triton");
	auto he3filter = new G4SDParticleFilter("He3Filter", "He3");
	auto alphafilter = new G4SDParticleFilter("AlphaFilter", "alpha");

  G4VPrimitiveScorer* scorer;


  char detectorname[12];
  char scorername[10];


  for (G4int idet = 1; idet < 4; idet++) {  // Loop for all MFD.
	  std::sprintf(detectorname, "Detector%02d", idet);
	  G4MultiFunctionalDetector* mfd = (G4MultiFunctionalDetector*)(SDman->FindSensitiveDetector(detectorname));

		  /*std::sprintf(scorername, "DTrack%02d", idet);
		  G4String PSnameT1(scorername);
		  scorer = new G4PSTrackLength(PSnameT1);
		  scorer->SetFilter(deuteronfilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "TTrack%02d", idet);
		  G4String PSnameT2(scorername);
		  scorer = new G4PSTrackLength(PSnameT2);
		  scorer->SetFilter(tritonfilter);    
		  mfd->RegisterPrimitive(scorer);*/

		  std::sprintf(scorername, "DEdep%02d", idet);
		  G4String PSnameE1(scorername);
		  scorer = new G4PSEnergyDeposit(PSnameE1);
		  scorer->SetFilter(deuteronfilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "TEdep%02d", idet);
		  G4String PSnameE2(scorername);
		  scorer = new G4PSEnergyDeposit(PSnameE2);
		  scorer->SetFilter(tritonfilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "PEdep%02d", idet);
		  G4String PSnameE3(scorername);
		  scorer = new G4PSEnergyDeposit(PSnameE3);
		  scorer->SetFilter(protonfilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "He3Edep%02d", idet);
		  G4String PSnameE4(scorername);
		  scorer = new G4PSEnergyDeposit(PSnameE4);
		  scorer->SetFilter(he3filter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "AEdep%02d", idet);
		  G4String PSnameE5(scorername);
		  scorer = new G4PSEnergyDeposit(PSnameE5);
		  scorer->SetFilter(alphafilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "NEdep%02d", idet);
		  G4String PSnameE6(scorername);
		  scorer = new G4PSEnergyDeposit(PSnameE6);
		  scorer->SetFilter(neutronfilter);    
		  mfd->RegisterPrimitive(scorer);

		  /*std::sprintf(scorername, "DFlSCur%02d", idet);
		  G4String PSnameFC1(scorername);
		  scorer = new G4PSTrackLength(PSnameFC1,fFlux_In);
		  scorer->SetFilter(deuteronfilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "TFlSCur%02d", idet);
		  G4String PSnameFC2(scorername);
		  scorer = new G4PSTrackLength(PSnameFC2,fFlux_In);
		  scorer->SetFilter(tritonfilter);    
		  mfd->RegisterPrimitive(scorer);*/

		  /*std::sprintf(scorername, "DNSec%02d", idet);
		  G4String PSnameNS1(scorername);
		  scorer = new G4PSNofSecondary(PSnameNS1);
		  scorer->SetFilter(deuteronfilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "TNSec%02d", idet);
		  G4String PSnameNS2(scorername);
		  scorer = new G4PSNofSecondary(PSnameNS2);
		  scorer->SetFilter(tritonfilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "PNSec%02d", idet);
		  G4String PSnameNS3(scorername);
		  scorer = new G4PSNofSecondary(PSnameNS3);
		  scorer->SetFilter(protonfilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "He3NSec%02d", idet);
		  G4String PSnameNS4(scorername);
		  scorer = new G4PSNofSecondary(PSnameNS4);
		  scorer->SetFilter(he3filter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "ANSec%02d", idet);
		  G4String PSnameNS5(scorername);
		  scorer = new G4PSNofSecondary(PSnameNS5);
		  scorer->SetFilter(alphafilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "NNSec%02d", idet);
		  G4String PSnameNS6(scorername);
		  scorer = new G4PSNofSecondary(PSnameNS6);
		  scorer->SetFilter(neutronfilter);    
		  mfd->RegisterPrimitive(scorer);*/

		  /*std::sprintf(scorername, "DSteps%02d", idet);
		  G4String PSnameNSt1(scorername);
		  scorer = new G4PSNofStep(PSnameNSt1);
		  scorer->SetFilter(deuteronfilter);    
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "TSteps%02d", idet);
		  G4String PSnameNSt2(scorername);
		  scorer = new G4PSNofStep(PSnameNSt2);
		  scorer->SetFilter(tritonfilter);    
		  mfd->RegisterPrimitive(scorer);*/

		  std::sprintf(scorername, "DAbs%02d", idet);
		  G4String PSname3(scorername);
		  scorer = new G4PSPopulation(PSname3);
		  scorer->SetFilter(deuteronfilter);
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "TAbs%02d", idet);
		  G4String PSname4(scorername);
		  scorer = new G4PSPopulation(PSname4);
		  scorer->SetFilter(tritonfilter);
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "PAbs%02d", idet);
		  G4String PSname2(scorername);
		  scorer = new G4PSPopulation(PSname2);
		  scorer->SetFilter(protonfilter);
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "He3Abs%02d", idet);
		  G4String PSname5(scorername);
		  scorer = new G4PSPopulation(PSname5);
		  scorer->SetFilter(he3filter);
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "AAbs%02d", idet);
		  G4String PSname6(scorername);
		  scorer = new G4PSPopulation(PSname6);
		  scorer->SetFilter(alphafilter);
		  mfd->RegisterPrimitive(scorer);
		  std::sprintf(scorername, "NAbs%02d", idet);
		  G4String PSnameA0(scorername);
		  scorer = new G4PSPopulation(PSnameA0);
		  scorer->SetFilter(neutronfilter);    
		  mfd->RegisterPrimitive(scorer);


	  for (G4int i = 0; i < 390; i++) {
		  std::sprintf(scorername, "N%02d%03d", idet, (i+1));
		  G4String PSname(scorername);
		  G4double kmin = (0 + 0.05*i)*MeV;
		  G4double kmax = (0.05 + 0.05*i)*MeV;
		  //-- Particle with kinetic energy filter.
		  G4SDParticleWithEnergyFilter* pkinEFilter = new G4SDParticleWithEnergyFilter("NeutronEnergyFilter", kmin, kmax);
		  pkinEFilter->add("neutron");  // Accept only neutrons
		  scorer = new G4PSPopulation(PSname);
		  scorer->SetFilter(pkinEFilter);    // Assign filter.
		  mfd->RegisterPrimitive(scorer);  // Register it to MultiFunctionalDetector.
	  }
  }

  for (G4int idet = 4; idet < 21; idet++) {  // Loop for all MFD.
	  std::sprintf(detectorname, "Detector%02d", idet);
	  G4MultiFunctionalDetector* mfd = (G4MultiFunctionalDetector*)(SDman->FindSensitiveDetector(detectorname));

		  std::sprintf(scorername, "NAbs%02d", idet);
		  G4String PSname1(scorername);
		  scorer = new G4PSCellFlux(PSname1);
		  scorer->SetFilter(neutronfilter);    
		  mfd->RegisterPrimitive(scorer);

	  for (G4int i = 0; i < 400; i++) {
		  std::sprintf(scorername, "N%02d%03d", idet, (i+1));
		  G4String PSname(scorername);
		  G4double kmin = (0 + 0.05*i)*MeV;
		  G4double kmax = (0.05 + 0.05*i)*MeV;
		  //-- Particle with kinetic energy filter.
		  G4SDParticleWithEnergyFilter* pkinEFilter = new G4SDParticleWithEnergyFilter("NeutronEnergyFilter", kmin, kmax);
		  pkinEFilter->add("neutron");  // Accept only neutrons
		  scorer = new G4PSCellFlux(PSname);
		  scorer->SetFilter(pkinEFilter);    // Assign filter.
		  mfd->RegisterPrimitive(scorer);  // Register it to MultiFunctionalDetector.
	  }

		  std::sprintf(scorername, "NPopAbs%02d", idet);
		  G4String PSname2(scorername);
		  scorer = new G4PSPopulation(PSname2);
		  scorer->SetFilter(neutronfilter);
		  mfd->RegisterPrimitive(scorer);
  }

	SetSensitiveDetector("Target", detector1);
	SetSensitiveDetector("ReactionChamber", detector2);
	SetSensitiveDetector("World", detector3);
	SetSensitiveDetector("DetectorR0", detector4);
	SetSensitiveDetector("DetectorR1Up", detector5);
	SetSensitiveDetector("DetectorR2Up", detector6);
	SetSensitiveDetector("DetectorR3Up", detector7);
	SetSensitiveDetector("DetectorR1Low", detector8);
	SetSensitiveDetector("DetectorR2Low", detector9);
	SetSensitiveDetector("DetectorR3Low", detector10);
	SetSensitiveDetector("DetectorA0Up", detector11);
	SetSensitiveDetector("DetectorA1Up", detector12);
	SetSensitiveDetector("DetectorA2Up", detector13);
	SetSensitiveDetector("DetectorA3Up", detector14);
	SetSensitiveDetector("DetectorA0Low", detector15);
	SetSensitiveDetector("DetectorA1Low", detector16);
	SetSensitiveDetector("DetectorA2Low", detector17);
	SetSensitiveDetector("DetectorA3Low", detector18);
	SetSensitiveDetector("StilbeneDetector1", detector19);
	SetSensitiveDetector("StilbeneDetector2", detector20);


	// -- Fetch volume for biasing:
	G4LogicalVolume* logicTest = G4LogicalVolumeStore::GetInstance()->GetVolume("Target");

	// ----------------------------------------------
	// -- operator creation and attachment to volume:
	// ----------------------------------------------
	GB01BOptrMultiParticleChangeCrossSection* testMany =
		new GB01BOptrMultiParticleChangeCrossSection();
	testMany->AddParticle("deuteron");
	testMany->AddParticle("triton");
	testMany->AttachTo(logicTest);
	G4cout << " Attaching biasing operator " << testMany->GetName()
		<< " to logical volume " << logicTest->GetName()
		<< G4endl;

	if (!fEMfield) {
		fEMfield = new F05Field();

		// Create an equation of motion for this field
		G4EqMagElectricField* fEquation = new G4EqMagElectricField(fEMfield);

		G4int nvar = 8;

		// Create the Runge-Kutta 'stepper' using the efficient 'DoPri5' method
		G4MagIntegratorStepper* fStepper = new G4DormandPrince745(fEquation, nvar);

		// Get the global field manager
		G4FieldManager* fFieldMgr = G4TransportationManager::GetTransportationManager()
			->GetFieldManager();
		// Set this field to the global field manager
		fFieldMgr->SetDetectorField(fEMfield);


		G4double fMinStep = 0.01*mm; // minimal step of 10 microns

									 // The driver will ensure that integration is control to give
									 //   acceptable integration error
		G4MagInt_Driver* fIntgrDriver = new G4MagInt_Driver(fMinStep,
			fStepper,
			fStepper->GetNumberOfVariables());

		G4ChordFinder* fChordFinder = new G4ChordFinder(fIntgrDriver);
		fFieldMgr->SetChordFinder(fChordFinder);
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
