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
/// \file B4dEventAction.hh
/// \brief Definition of the B4dEventAction class

#ifndef ModelEventAction_h
#define ModelEventAction_h 1

#include "G4UserEventAction.hh"

#include "G4THitsMap.hh"
#include "globals.hh"
#include <vector>

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in the hits collections.

class ModelEventAction : public G4UserEventAction
{
public:
  ModelEventAction();
  virtual ~ModelEventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);
    
private:
  // methods
  G4THitsMap<G4double>* GetHitsCollection(G4int hcID, const G4Event* event) const;
  G4double GetSum(G4THitsMap<G4double>* hitsMap) const;
  
  //data members 
  std::vector<G4String> MFDname;                  
  /* G4String MFDname[18] = {"Detector01","Detector02","Detector03","Detector04","Detector05","Detector06","Detector07",
	"Detector08","Detector09","Detector10","Detector11","Detector12","Detector13","Detector14","Detector15","Detector16",
	"Detector17","Detector18"}; */

};
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
