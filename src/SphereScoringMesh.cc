#include "SphereScoringMesh.hh"

#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cmath>

SphereScoringMesh::SphereScoringMesh(const G4String& name)
    : G4VScoringMesh(name)
{
    fShape = MeshShape::sphere;
    fDivisionAxisNames[0] = "theta";
    fDivisionAxisNames[1] = "phi";
    fDivisionAxisNames[2] = "";
    fNSegment[0] = 1;
    fNSegment[1] = 1;
    fNSegment[2] = 1;
    fSize[0] = 0.;
    fSize[1] = 0.;
    fSize[2] = 0.;
}

void SphereScoringMesh::SetupGeometry(G4VPhysicalVolume* world)
{
    AttachToWorld(world);
}

void SphereScoringMesh::WorkerConstruct(G4VPhysicalVolume* world)
{
    AttachToWorld(world);
}

void SphereScoringMesh::AttachToWorld(G4VPhysicalVolume* world)
{
    if (!world) return;

    auto* logical = world->GetLogicalVolume();
    if (!logical) return;

    logical->SetSensitiveDetector(fMFD);
    SetMeshElementLogical(logical);
}

void SphereScoringMesh::List() const
{
    G4cout << "Spherical scoring mesh <" << fWorldName << ">" << G4endl;
    G4cout << "  center: "
           << fCenterPosition.x() / cm << " "
           << fCenterPosition.y() / cm << " "
           << fCenterPosition.z() / cm << " cm" << G4endl;
    G4cout << "  radius: " << GetSphereRadius() / cm << " cm" << G4endl;
    G4cout << "  bins: theta=" << GetThetaBins()
           << " phi=" << GetPhiBins() << G4endl;
}

void SphereScoringMesh::Draw(RunScore*, G4VScoreColorMap*, G4int)
{
    G4cout << "[SphereScore] Drawing is not implemented for spherical scoring meshes." << G4endl;
}

void SphereScoringMesh::DrawColumn(RunScore*, G4VScoreColorMap*, G4int, G4int)
{
    G4cout << "[SphereScore] Column drawing is not implemented for spherical scoring meshes." << G4endl;
}

void SphereScoringMesh::SetSphereCenter(const G4ThreeVector& center)
{
    fCenterPosition = center;
}

void SphereScoringMesh::SetSphereRadius(G4double radius)
{
    fSize[0] = radius;
    fSize[1] = radius;
    fSize[2] = radius;
    sizeIsSet = (radius > 0.);
}

void SphereScoringMesh::SetSphereBins(G4int nTheta, G4int nPhi)
{
    fNSegment[0] = std::max(1, nTheta);
    fNSegment[1] = std::max(1, nPhi);
    fNSegment[2] = 1;
    nMeshIsSet = true;
}

G4ThreeVector SphereScoringMesh::GetSphereCenter() const
{
    return fCenterPosition;
}

G4double SphereScoringMesh::GetSphereRadius() const
{
    return fSize[0];
}

G4int SphereScoringMesh::GetThetaBins() const
{
    return fNSegment[0];
}

G4int SphereScoringMesh::GetPhiBins() const
{
    return fNSegment[1];
}

G4int SphereScoringMesh::GetSphereIndex(const G4ThreeVector& globalPoint) const
{
    const auto rel = globalPoint - fCenterPosition;
    const auto r = rel.mag();

    G4double theta = 0.;
    G4double phi = 0.;
    if (r > 0.) {
        const auto cosTheta = std::clamp(rel.z() / r, -1.0, 1.0);
        theta = std::acos(cosTheta);
        phi = std::atan2(rel.y(), rel.x());
        if (phi < 0.) phi += twopi;
    }

    auto thetaBin = static_cast<G4int>(std::floor(theta / pi * GetThetaBins()));
    auto phiBin = static_cast<G4int>(std::floor(phi / twopi * GetPhiBins()));

    thetaBin = std::clamp(thetaBin, 0, GetThetaBins() - 1);
    phiBin = std::clamp(phiBin, 0, GetPhiBins() - 1);

    return GetSphereIndex(thetaBin, phiBin);
}

G4int SphereScoringMesh::GetSphereIndex(G4int thetaBin, G4int phiBin) const
{
    thetaBin = std::clamp(thetaBin, 0, GetThetaBins() - 1);
    phiBin = std::clamp(phiBin, 0, GetPhiBins() - 1);
    return thetaBin * GetPhiBins() + phiBin;
}

void SphereScoringMesh::GetThetaPhiBins(G4int index, G4int& thetaBin, G4int& phiBin) const
{
    const auto nPhi = GetPhiBins();
    thetaBin = index / nPhi;
    phiBin = index % nPhi;
}

G4double SphereScoringMesh::GetThetaMin(G4int thetaBin) const
{
    return pi * thetaBin / GetThetaBins();
}

G4double SphereScoringMesh::GetThetaMax(G4int thetaBin) const
{
    return pi * (thetaBin + 1) / GetThetaBins();
}

G4double SphereScoringMesh::GetPhiMin(G4int phiBin) const
{
    return twopi * phiBin / GetPhiBins();
}

G4double SphereScoringMesh::GetPhiMax(G4int phiBin) const
{
    return twopi * (phiBin + 1) / GetPhiBins();
}
