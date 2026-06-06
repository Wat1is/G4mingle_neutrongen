#ifndef PSSphereSurfaceCurrent_hh
#define PSSphereSurfaceCurrent_hh 1

#include "G4THitsMap.hh"
#include "G4VPrimitiveScorer.hh"

class SphereScoringMesh;

class PSSphereSurfaceCurrent : public G4VPrimitiveScorer
{
public:
    PSSphereSurfaceCurrent(const G4String& name, SphereScoringMesh* mesh, G4int direction);
    ~PSSphereSurfaceCurrent() override = default;

    void Initialize(G4HCofThisEvent* hce) override;
    void clear() override;
    void PrintAll() override;

protected:
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

private:
    G4bool ComputeCrossingPoint(const G4ThreeVector& pre,
                                const G4ThreeVector& post,
                                G4ThreeVector& crossing) const;

    SphereScoringMesh* fMesh = nullptr;
    G4int fDirection = 0;
    G4int fHCID = -1;
    G4THitsMap<G4double>* fEvtMap = nullptr;
};

#endif
