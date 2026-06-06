#include "SphereScoringMessenger.hh"

#include "PSSphereSurfaceCurrent.hh"
#include "SphereScoringMesh.hh"

#include "G4ScoringManager.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "G4ios.hh"

SphereScoringMessenger::SphereScoringMessenger(G4ScoringManager* scoringManager)
    : fScoringManager(scoringManager)
{
    fCreateSphereCmd = new G4UIcmdWithAString("/score/create/sphere", this);
    fCreateSphereCmd->SetGuidance("Create spherical scoring mesh.");
    fCreateSphereCmd->SetParameterName("MeshName", false);

    fSphereDir = new G4UIdirectory("/score/sphere/");
    fSphereDir->SetGuidance("Spherical scoring mesh commands.");

    fCenterCmd = new G4UIcmdWith3VectorAndUnit("/score/sphere/center", this);
    fCenterCmd->SetGuidance("Set center of the current spherical scoring mesh.");
    fCenterCmd->SetParameterName("X", "Y", "Z", false, false);
    fCenterCmd->SetDefaultUnit("mm");

    fRadiusCmd = new G4UIcmdWithADoubleAndUnit("/score/sphere/radius", this);
    fRadiusCmd->SetGuidance("Set radius of the current spherical scoring mesh.");
    fRadiusCmd->SetParameterName("R", false);
    fRadiusCmd->SetRange("R>0.");
    fRadiusCmd->SetDefaultUnit("mm");

    fBinCmd = new G4UIcommand("/score/sphere/nSphereBin", this);
    fBinCmd->SetGuidance("Set theta/phi binning of the current spherical scoring mesh.");
    auto* param = new G4UIparameter("nTheta", 'i', false);
    param->SetParameterRange("nTheta>0");
    fBinCmd->SetParameter(param);
    param = new G4UIparameter("nPhi", 'i', false);
    param->SetParameterRange("nPhi>0");
    fBinCmd->SetParameter(param);

    fSurfaceCurrentCmd = new G4UIcommand("/score/quantity/sphereSurfaceCurrent", this);
    fSurfaceCurrentCmd->SetGuidance("Create a theta/phi-binned sphere surface-current primitive scorer.");
    param = new G4UIparameter("quantityName", 's', false);
    fSurfaceCurrentCmd->SetParameter(param);
    param = new G4UIparameter("direction", 'i', false);
    param->SetParameterRange("direction == -1 || direction == 0 || direction == 1");
    fSurfaceCurrentCmd->SetParameter(param);
}

SphereScoringMessenger::~SphereScoringMessenger()
{
    delete fSurfaceCurrentCmd;
    delete fBinCmd;
    delete fRadiusCmd;
    delete fCenterCmd;
    delete fSphereDir;
    delete fCreateSphereCmd;
}

void SphereScoringMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == fCreateSphereCmd) {
        if (fScoringManager->FindMesh(newValue)) {
            G4cerr << "ERROR : Scoring mesh <" << newValue << "> already exists." << G4endl;
            return;
        }

        fScoringManager->RegisterScoringMesh(new SphereScoringMesh(newValue));
        return;
    }

    if (command == fCenterCmd) {
        auto* mesh = GetCurrentSphereMesh(command->GetCommandPath());
        if (!mesh) return;
        mesh->SetSphereCenter(fCenterCmd->GetNew3VectorValue(newValue));
        return;
    }

    if (command == fRadiusCmd) {
        auto* mesh = GetCurrentSphereMesh(command->GetCommandPath());
        if (!mesh) return;
        mesh->SetSphereRadius(fRadiusCmd->GetNewDoubleValue(newValue));
        return;
    }

    if (command == fBinCmd) {
        auto* mesh = GetCurrentSphereMesh(command->GetCommandPath());
        if (!mesh) return;

        G4Tokenizer next(newValue);
        const auto nTheta = StoI(next());
        const auto nPhi = StoI(next());
        mesh->SetSphereBins(nTheta, nPhi);
        return;
    }

    if (command == fSurfaceCurrentCmd) {
        auto* mesh = GetCurrentSphereMesh(command->GetCommandPath());
        if (!mesh) return;

        G4Tokenizer next(newValue);
        const auto psName = next();
        const auto direction = StoI(next());

        if (direction != -1 && direction != 0 && direction != 1) {
            G4cerr << "ERROR : /score/quantity/sphereSurfaceCurrent direction must be -1, 0, or 1." << G4endl;
            return;
        }

        if (!mesh->ReadyForQuantity()) {
            G4cerr << "ERROR : Spherical scoring mesh is not ready. Set radius and nSphereBin first." << G4endl;
            return;
        }

        if (mesh->FindPrimitiveScorer(psName)) {
            G4cerr << "ERROR : Primitive scorer <" << psName
                   << "> already exists in mesh <" << mesh->GetWorldName() << ">." << G4endl;
            return;
        }

        auto* scorer = new PSSphereSurfaceCurrent(psName, mesh, direction);
        scorer->SetNijk(mesh->GetThetaBins(), mesh->GetPhiBins(), 1);
        mesh->SetPrimitiveScorer(scorer);
        return;
    }
}

SphereScoringMesh* SphereScoringMessenger::GetCurrentSphereMesh(const G4String& commandName) const
{
    auto* mesh = dynamic_cast<SphereScoringMesh*>(fScoringManager->GetCurrentMesh());
    if (!mesh) {
        G4cerr << "ERROR : " << commandName
               << " requires the current scoring mesh to be a sphere mesh." << G4endl;
    }
    return mesh;
}
