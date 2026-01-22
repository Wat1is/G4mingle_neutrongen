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
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class

#include "ModelRunAction.hh"
#include "ModelAnalysis.hh"

#include "G4Run.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModelRunAction::ModelRunAction()
 : G4UserRunAction()
{ 
	char detectorname[12];
	for (G4int idet = 1; idet < 21; idet++) {  // Loop for all MFD.
		std::sprintf(detectorname, "Detector%02d", idet);
		G4String mfdName(detectorname);
		MFDname.push_back(G4String(mfdName));
	}

  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);   

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModelRunAction::~ModelRunAction()
{
  delete G4AnalysisManager::Instance(); 
  MFDname.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModelRunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysismanager = G4AnalysisManager::Instance();
  analysismanager->SetVerboseLevel(1);

  // Open an output file
  //
 
  analysismanager->OpenFile("Result");

  G4cout << "Using " << analysismanager->GetType() << G4endl;

  // Book histograms, ntuple


  // Creating ntuple
  //

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  analysismanager->CreateNtuple("Results", "Detector_results");
 
  G4int mfdnumber = MFDname.size();
  for (G4int idet = 0; idet < mfdnumber; idet++) {  // Loop for all MFD.
	  G4String detName = MFDname[idet];
	  //--- Seek and Obtain MFD objects from SDmanager.
	  G4MultiFunctionalDetector* mfd =(G4MultiFunctionalDetector*)(SDman->FindSensitiveDetector(detName));
	  //
	  if (mfd) {
		  //--- Loop over the registered primitive scorers.
		  for (G4int iscor = 0; iscor < mfd->GetNumberOfPrimitives(); iscor++) {
			  // Get Primitive Scorer object.
			  G4VPrimitiveScorer* scorer = mfd->GetPrimitive(iscor);
			  G4String scorername = scorer->GetName();
			  analysismanager->CreateNtupleDColumn(scorername);
		  }
	  }
  }
  analysismanager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModelRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  auto analysismanager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  //
  analysismanager->Write();
  analysismanager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
