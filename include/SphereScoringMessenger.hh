#ifndef SphereScoringMessenger_hh
#define SphereScoringMessenger_hh 1

#include "G4UImessenger.hh"

class G4ScoringManager;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcommand;
class G4UIdirectory;
class SphereScoringMesh;

class SphereScoringMessenger : public G4UImessenger
{
public:
    explicit SphereScoringMessenger(G4ScoringManager* scoringManager);
    ~SphereScoringMessenger() override;

    void SetNewValue(G4UIcommand* command, G4String newValue) override;

private:
    SphereScoringMesh* GetCurrentSphereMesh(const G4String& commandName) const;

    G4ScoringManager* fScoringManager = nullptr;
    G4UIcmdWithAString* fCreateSphereCmd = nullptr;
    G4UIdirectory* fSphereDir = nullptr;
    G4UIcmdWith3VectorAndUnit* fCenterCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fRadiusCmd = nullptr;
    G4UIcommand* fBinCmd = nullptr;
    G4UIcommand* fSurfaceCurrentCmd = nullptr;
};

#endif
