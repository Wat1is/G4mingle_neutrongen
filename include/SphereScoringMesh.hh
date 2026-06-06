#ifndef SphereScoringMesh_hh
#define SphereScoringMesh_hh 1

#include "G4ThreeVector.hh"
#include "G4VScoringMesh.hh"

class G4VPhysicalVolume;

class SphereScoringMesh : public G4VScoringMesh
{
public:
    explicit SphereScoringMesh(const G4String& name);
    ~SphereScoringMesh() override = default;

    void WorkerConstruct(G4VPhysicalVolume* world) override;
    void List() const override;
    void Draw(RunScore* map, G4VScoreColorMap* colorMap, G4int axflg = 111) override;
    void DrawColumn(RunScore* map, G4VScoreColorMap* colorMap, G4int idxProj, G4int idxColumn) override;

    void SetSphereCenter(const G4ThreeVector& center);
    void SetSphereRadius(G4double radius);
    void SetSphereBins(G4int nTheta, G4int nPhi);

    G4ThreeVector GetSphereCenter() const;
    G4double GetSphereRadius() const;
    G4int GetThetaBins() const;
    G4int GetPhiBins() const;

    G4int GetSphereIndex(const G4ThreeVector& globalPoint) const;
    G4int GetSphereIndex(G4int thetaBin, G4int phiBin) const;
    void GetThetaPhiBins(G4int index, G4int& thetaBin, G4int& phiBin) const;

    G4double GetThetaMin(G4int thetaBin) const;
    G4double GetThetaMax(G4int thetaBin) const;
    G4double GetPhiMin(G4int phiBin) const;
    G4double GetPhiMax(G4int phiBin) const;

protected:
    void SetupGeometry(G4VPhysicalVolume* world) override;

private:
    void AttachToWorld(G4VPhysicalVolume* world);
};

#endif
