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
/// \file B4dEventAction.cc
/// \brief Implementation of the B4dEventAction class

#include "ModelEventAction.hh"
#include "ModelAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModelEventAction::ModelEventAction()
 : G4UserEventAction()
{
	char detectorname[12];
	for (G4int idet = 1; idet < 21; idet++) {  // Loop for all MFD.
		std::sprintf(detectorname, "Detector%02d", idet);
		G4String mfdName(detectorname);
		MFDname.push_back(G4String(mfdName));
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModelEventAction::~ModelEventAction()
{
  MFDname.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4THitsMap<G4double>* 
ModelEventAction::GetHitsCollection(G4int hcID, const G4Event* event) const
{
  auto hitsCollection = static_cast<G4THitsMap<G4double>*>(event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("B4dEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ModelEventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
  G4double sumValue = 0.;
  for ( auto it : *hitsMap->GetMap() ) {
    // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>
    sumValue += *(it.second);
  }
  return sumValue;  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModelEventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModelEventAction::EndOfEventAction(const G4Event* event)
{  

  auto analysisManager = G4AnalysisManager::Instance();
   // Get hist collections IDs
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	//=================================================
	//  Initalize RunMaps for accumulation.
	//  Get CollectionIDs for HitCollections.
	//=================================================

	G4int mfdNum=MFDname.size();
	//G4int mfdNum=18;
	for (G4int idet = 0; idet < mfdNum; idet++) {  // Loop for all MFD.
		G4String detectorname = MFDname[idet];
		//--- Seek and Obtain MFD objects from SDmanager.
		G4MultiFunctionalDetector* mfd =(G4MultiFunctionalDetector*)(SDman->FindSensitiveDetector(detectorname));
		if (mfd) {
			//--- Loop over the registered primitive scorers.
			for (G4int iscor = 0; iscor < mfd->GetNumberOfPrimitives(); iscor++) {
				// Get Primitive Scorer object.
				G4VPrimitiveScorer* scorer = mfd->GetPrimitive(iscor);
				// collection name and collectionID for HitsCollection,
				G4String scorername = scorer->GetName();
				G4String fullCollectionName = detectorname + "/" + scorername;
				G4int    hitcollectionID = SDman->GetCollectionID(fullCollectionName);

				if (hitcollectionID >= 0) {
					//G4cout << "++ " << fullCollectionName << " id " << hitcollectionID << G4endl;
					// gain score quantity and store it into 
			        auto value = GetSum(GetHitsCollection(hitcollectionID, event));
			        analysisManager->FillNtupleDColumn(hitcollectionID, value);
				}
				else {
					G4cout << "** collection " << fullCollectionName << " not found. "
						<< G4endl;
				}
			}
		}
	} 


  analysisManager->AddNtupleRow();  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
