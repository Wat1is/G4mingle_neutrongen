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
/// \file field/field05/src/F05Field.cc
/// \brief Implementation of the F05Field class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05Field.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05Field::F05Field() : G4ElectroMagneticField()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05Field::~F05Field()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05Field::GetFieldValue(const G4double Point[4], G4double* Bfield) const
{
	// Point[0],Point[1],Point[2] are x-, y-, z-cordinates, Point[3] is time

	const G4double R1 = 20.0 * mm;
	const G4double R2 = 96.0 * mm;
	const G4double U = 8.0e7 * volt / m;
	const G4double Rmin = 10.0 * mm;
	const G4double Rmax = 48.0 * mm;
	const G4double Zmax = 180.0 * mm;

	const G4double x = Point[0];
	const G4double y = Point[1];
	const G4double z = Point[2];
	const G4double posR = std::sqrt(x * x + y * y);

	G4double Ex = 0.0;
	G4double Ey = 0.0;

	if (posR > 0.0 && posR >= Rmin && posR <= Rmax && std::fabs(z) <= Zmax && R2 > R1)
	{
		const G4double Er = U / (std::log(R2 / R1) * posR);
		const G4double cosTheta = x / posR;
		const G4double sinTheta = y / posR;
		Ex = -Er * cosTheta;
		Ey = -Er * sinTheta;
	}

	Bfield[0] = 0.0;
	Bfield[1] = 0.0;
	Bfield[2] = 0.0;
	Bfield[3] = Ex;
	Bfield[4] = Ey;
	Bfield[5] = 0.0;
}
