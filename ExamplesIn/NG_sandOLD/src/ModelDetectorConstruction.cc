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
	G4Element* Cl = man->FindOrBuildElement("Cl");
	G4Element* Si = man->FindOrBuildElement("Si");


	G4Element* D = new G4Element("Deuterium", "D", 1);
	D->AddIsotope(H2, 100.*perCent);
	G4Element* T = new G4Element("Tritium", "T", 1);
	T->AddIsotope(H3, 100.*perCent);

	G4Material* air = man->FindOrBuildMaterial("G4_AIR");
	G4Material* water = man->FindOrBuildMaterial("G4_WATER");
	G4Material* vacuum = man->FindOrBuildMaterial("G4_Galactic");
	G4Material* lead = man->FindOrBuildMaterial("G4_Pb");
	G4Material* aluminium = man->FindOrBuildMaterial("G4_Al");
	G4Material* copper = man->FindOrBuildMaterial("G4_Cu");

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

	G4Material* pvc = new G4Material("PVC", 1.406 * g / cm3, 3);
	pvc->AddElement(H, 0.048382);
	pvc->AddElement(C, 0.384361);
	pvc->AddElement(Cl, 0.567257);

	G4Material* insoil = new G4Material("InsulationOil", 0.81 * g / cm3, 2);
	insoil->AddElement(H, 41);
	insoil->AddElement(C, 20);

	G4Material* sand = new G4Material("Sand", 1.45 * g / cm3, 2);
	sand->AddElement(Si, 1);
	sand->AddElement(O, 2);

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
	auto copper = G4Material::GetMaterial("G4_Cu");
	auto insulationoil = G4Material::GetMaterial("InsulationOil");
	auto pvc = G4Material::GetMaterial("PVC");
	auto sand = G4Material::GetMaterial("Sand");




	//
	// World
	//
	auto worldS = new G4Box("World", 101. * cm, 101. * cm, 60. * cm);
	auto worldLV = new G4LogicalVolume(worldS, air, "World");
	auto worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "World",
		0, false, 0, fCheckOverlaps);

	/* auto sandbarrelS = new G4Tubs("SandBarrel", 78. * mm, 287. * mm, 300. * mm, 0. * deg, 360. * deg);
	 auto sandbarrelLV = new G4LogicalVolume(sandbarrelS, sand, "SandBarrel");
	 new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), sandbarrelLV,
		 "SandBarrel", worldLV, false, 0, fCheckOverlaps);
	*/
	/*
	auto shieldingS = new G4Tubs("Shielding", 0. * mm, 89. * mm, 400. * mm, 0. * deg, 360. * deg);
	auto shieldingLV = new G4LogicalVolume(shieldingS, lead, "Shielding");
	new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), shieldingLV,
		"Shielding", worldLV, false, 0, fCheckOverlaps);
	*/

	auto casingS = new G4Tubs("Casing",0.*mm,78.*mm, 400.*mm,0.*deg, 360. * deg);
	auto casingLV = new G4LogicalVolume(casingS,aluminium,"Casing");
	new G4PVPlacement(0,G4ThreeVector(0.*mm, 0.*mm,0.*mm),casingLV,
                   "Casing",worldLV,false,0,fCheckOverlaps);

    auto coolingS = new G4Tubs("Cooling",58.*mm,75.*mm,331*mm,0. * deg,360. * deg);
	auto coolingLV = new G4LogicalVolume(coolingS,coolingmat,"Cooling");
	new G4PVPlacement(0,G4ThreeVector(0.*mm, 0.*mm,-48.*mm),coolingLV,
                   "Cooling",casingLV,false,0,fCheckOverlaps);
	auto coolingupS = new G4Tubs("CoolingUp", 0. * mm, 75. * mm, 51 * mm, 0. * deg, 360. * deg);
	auto coolingupLV = new G4LogicalVolume(coolingupS, coolingmat, "CoolingUp");
	new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 349. * mm), coolingupLV,
		"CoolingUp", casingLV, false, 0, fCheckOverlaps);
	auto coolingdownS = new G4Tubs("CoolingDown", 58. * mm, 75. * mm, 3 * mm, 0. * deg, 360. * deg);
	auto coolingdownLV = new G4LogicalVolume(coolingdownS, coolingmat, "CoolingDown");
	new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, -397. * mm), coolingdownLV,
		"CoolingDown", casingLV, false, 0, fCheckOverlaps);

	auto getterS = new G4Tubs("Getter",0.*mm,58.*mm,13*mm,0. * deg,360. * deg);
	auto getterLV = new G4LogicalVolume(getterS,gettermat,"Getter");
	new G4PVPlacement(0,G4ThreeVector(0.*mm, 0.*mm, -387.*mm),getterLV,
                   "Getter",casingLV,false,0,fCheckOverlaps);

    auto reactchamberS = new G4Tubs("ReactionChamber",0.*mm,48.*mm,328.5*mm,0. * deg,360. * deg);
	auto reactchamberLV = new G4LogicalVolume(reactchamberS,dtgas,"ReactionChamber");
	new G4PVPlacement(0,G4ThreeVector(0.*mm, 0.*mm, -45.5*mm),reactchamberLV,
                   "ReactionChamber",casingLV,false,0,fCheckOverlaps);

    auto insulatortopS = new G4Tubs("InsulatorTop",0.,26.5*mm,54.*mm,0. * deg,360. * deg);
	auto insulatortopLV = new G4LogicalVolume(insulatortopS,al2o3,"InsulatorTop");
	new G4PVPlacement(0,G4ThreeVector(0.*mm,0.*mm,274.5*mm),insulatortopLV,
                   "InsulatorTop",reactchamberLV,false,0,fCheckOverlaps);
	auto insulatorbottomS = new G4Tubs("InsulatorBottom",0.,26.5*mm,88.*mm, 0. * deg, 360. * deg);
	auto insulatorbottomLV = new G4LogicalVolume(insulatorbottomS,al2o3,"InsulatorBottom");
	new G4PVPlacement(0,G4ThreeVector(0.*mm, 0.*mm,-217.5*mm),insulatorbottomLV,
                   "InsulatorBottom",reactchamberLV,false,0,fCheckOverlaps);
	auto insulatortopoutsideS = new G4Tubs("InsulatorTopOutside", 0., 38 * mm, 47.5 * mm, 0. * deg, 360. * deg);
	auto insulatortopoutsideLV = new G4LogicalVolume(insulatortopoutsideS, al2o3, "InsulatorTopOutside");
	new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, -3.5 * mm), insulatortopoutsideLV,
		"InsulatorTopOutside", coolingupLV, false, 0, fCheckOverlaps);

	auto oilinsulatorS = new G4Tubs("OilInsulator", 0., 46.5 * mm, 3.5 * mm, 0. * deg, 360. * deg);
	auto oilinsulatorLV = new G4LogicalVolume(oilinsulatorS, insulationoil, "OilInsulator");
	new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 47.5 * mm), oilinsulatorLV,
		"OilInsulator", coolingupLV, false, 0, fCheckOverlaps);
	auto cableinsulatorS = new G4Tubs("CableInsulator", 0., 11. * mm, 3.5 * mm, 0. * deg, 360. * deg);
	auto cableinsulatorLV = new G4LogicalVolume(cableinsulatorS, pvc, "CableInsulator");
	new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), cableinsulatorLV,
		"CableInsulator", oilinsulatorLV, false, 0, fCheckOverlaps);
	auto cableS = new G4Tubs("Cable", 0., 9. * mm, 3.5 * mm, 0. * deg, 360. * deg);
	auto cableLV = new G4LogicalVolume(cableS, copper, "Cable");
	new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), cableLV,
		"Cable", cableinsulatorLV, false, 0, fCheckOverlaps);

	auto getterbottomS = new G4Tubs("GetterBottom",0.,29.5*mm,11.5*mm,0. * deg,360. * deg);
	auto getterbottomLV = new G4LogicalVolume(getterbottomS,gettermat,"GetterBottom");
	new G4PVPlacement(0,G4ThreeVector(0.*mm, 0.*mm,-317*mm),getterbottomLV,
                   "GetterBottom",reactchamberLV,false,0,fCheckOverlaps);

	auto targetS = new G4Tubs("Target",0.,10.*mm,175.*mm,0. * deg,360. * deg);
	auto targetLV = new G4LogicalVolume(targetS,tar,"Target");
	new G4PVPlacement(0,G4ThreeVector(0.*mm, 0.*mm, 45.5*mm),targetLV,
                   "Target",reactchamberLV,false,0,fCheckOverlaps);

	/*auto stilbenedetector1S=new G4Tubs("StilbeneDetector1",0, 25.*mm, 25.*mm, 0.*deg, 360.*deg);
	auto stilbenedetector1LV=new G4LogicalVolume(stilbenedetector1S,vacuum,"StilbeneDetector1");
	G4RotationMatrix* xRot = new G4RotationMatrix;
	xRot->rotateX(90.*deg);
	new G4PVPlacement(xRot,G4ThreeVector(0.*mm, 1000.*mm, 0.*mm),stilbenedetector1LV,
			"StilbeneDetector1",worldLV,false,0,fCheckOverlaps);*/

	auto detectorS = new G4Tubs("BeltDetector", 1000 * mm, 1000.1 * mm, 22.5 * mm, 0. * deg, 360. * deg);
	auto detectorLV = new G4LogicalVolume(detectorS, vacuum, "BeltDetector");
	fdetectorPV = new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), detectorLV,
		"BeltDetector", worldLV, false, 0, fCheckOverlaps);


  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  /*auto yellow = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
  auto white = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  auto red = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  auto blue = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  auto gray = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  auto green = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
  auto magenta = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
  auto black = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0));

  auto brasscol = new G4VisAttributes(G4Colour(0.71, 0.651, 0.259));
  auto alcol = new G4VisAttributes(G4Colour(0.647, 0.649, 0.653));
  auto nbcol = new G4VisAttributes(G4Colour(0.498, 0.498, 0.498));
  auto steelgrey = new G4VisAttributes(G4Colour(0.796, 0.8, 0.808));
  auto watercol = new G4VisAttributes(G4Colour(0.247, 0.278, 0.8));
  auto vacuumcol = new G4VisAttributes(G4Colour(0.9609, 0.9844, 0.9844));
  auto hecol = new G4VisAttributes(G4Colour(0.953, 0.769, 0.914));
  auto pvccolor = new G4VisAttributes(G4Colour(0.992, 0.373, 0.0));
  auto cucolor = new G4VisAttributes(G4Colour(0.72, 0.45, 0.2));

  auto oilcol = new G4VisAttributes(G4Colour(0.67, 1.0, 0.67));
  auto insulatorcol = new G4VisAttributes(G4Colour(1.0, 0.5, 0.16));
  auto casingcol = new G4VisAttributes(G4Colour(0.33, 0.6, 1.0));
  auto shieldingcol = new G4VisAttributes(G4Colour(0.44, 0.54, 0.57));
  auto chambercol = new G4VisAttributes(G4Colour(1.0, 0.5, 0.9));
  auto targetcol = new G4VisAttributes(G4Colour(1.0, 0.33, 0.6));
  auto coolingcol = new G4VisAttributes(G4Colour(1.0, 0.84, 0.84));
  auto gettercol = new G4VisAttributes(G4Colour(0.22, 0.67, 0.78));

      shieldingLV->SetVisAttributes(shieldingcol);
	  casingLV->SetVisAttributes(casingcol);
	  coolingLV->SetVisAttributes(coolingcol);
	  coolingupLV->SetVisAttributes(coolingcol);
	  coolingdownLV->SetVisAttributes(coolingcol);
	  getterLV->SetVisAttributes(gettercol);
	  reactchamberLV->SetVisAttributes(chambercol);
	  insulatortopLV->SetVisAttributes(insulatorcol);
	  insulatorbottomLV->SetVisAttributes(insulatorcol);
	  insulatortopoutsideLV->SetVisAttributes(insulatorcol);
	  oilinsulatorLV->SetVisAttributes(oilcol);
	  cableinsulatorLV->SetVisAttributes(black);
	  cableLV->SetVisAttributes(cucolor);
	  getterbottomLV->SetVisAttributes(gettercol);
	  targetLV->SetVisAttributes(targetcol);
	  teflon1LV->SetVisAttributes(white);
	  teflon2LV->SetVisAttributes(white);
	  teflon3LV->SetVisAttributes(white);
	  teflon4LV->SetVisAttributes(white);
	  teflon5LV->SetVisAttributes(white);
	  teflon6LV->SetVisAttributes(white);
	  teflon7LV->SetVisAttributes(white);
	  teflon8LV->SetVisAttributes(white);
	  teflon9LV->SetVisAttributes(white);
	  teflon10LV->SetVisAttributes(white);
	  teflon11LV->SetVisAttributes(white);
	  teflon12LV->SetVisAttributes(white);
	  teflon13LV->SetVisAttributes(white);
	  teflon14LV->SetVisAttributes(white);
	  teflon15LV->SetVisAttributes(white);
	  teflon16LV->SetVisAttributes(white);
	  teflon17LV->SetVisAttributes(white);
	  teflon18LV->SetVisAttributes(white);
	  teflon19LV->SetVisAttributes(white);
	  teflon20LV->SetVisAttributes(white);*/

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreadLocal F05Field* ModelDetectorConstruction::fEMfield = 0;

void ModelDetectorConstruction::ConstructSDandField()
{

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


		G4double fMinStep = 0.01 * mm; // minimal step of 10 microns

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
