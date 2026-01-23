#include <G4VUserDetectorConstruction.hh>
#include <G4tgbVolumeMgr.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4Neutron.hh>
#include <G4Proton.hh>
#include <G4Deuteron.hh>
#include <G4Triton.hh>
#include <atomic>
#include "GB01BOptrMultiParticleChangeCrossSection.hh"
#include <G4LogicalVolumeStore.hh>

namespace {
    std::atomic<bool> gEnableDTXSBoost{ false };
}

class Detector : public G4VUserDetectorConstruction
{
	public:
		G4VPhysicalVolume* Construct() override {
			G4tgbVolumeMgr::GetInstance()->AddTextFile("detector.tg");
			return G4tgbVolumeMgr::GetInstance()->ReadAndConstructDetector();
		} ///< load detector definition from a text file "detector.tg"
        void ConstructSDandField() override
        {
            // One operator per thread (MT-safe)
            static thread_local bool attached = false;
            static thread_local GB01BOptrMultiParticleChangeCrossSection* op = nullptr;

            if (attached) return;
            if (!gEnableDTXSBoost.load(std::memory_order_relaxed)) return;

            op = new GB01BOptrMultiParticleChangeCrossSection();
            op->AddParticle("deuteron");
            op->AddParticle("triton");

            // Attach to ALL logical volumes => “everywhere”
            for (auto* lv : *G4LogicalVolumeStore::GetInstance())
            {
                if (lv) op->AttachTo(lv);
            }

            attached = true;
        }
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <G4LogicalVolumeStore.hh>
#include <G4GenericBiasingPhysics.hh>
#include <G4VModularPhysicsList.hh>

#include <G4UIdirectory.hh>



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4GeneralParticleSource.hh>

class Generator : public G4VUserPrimaryGeneratorAction
{
	private:
		G4GeneralParticleSource* fGPS;
	public:
		Generator() : G4VUserPrimaryGeneratorAction() {
			fGPS = new G4GeneralParticleSource; }
		~Generator() { delete fGPS; }
		void GeneratePrimaries(G4Event *evt) { fGPS->GeneratePrimaryVertex(evt); }
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <G4UserSteppingAction.hh>

namespace {
    std::atomic<bool> gKeepOnlyDTpnEnabled{ false };



    inline void RegisterDTBiasingPhysics(G4VModularPhysicsList* pl)
    {
        // Required so GB01 operators can see wrapped biasing processes for these particles.
        auto* biasPhys = new G4GenericBiasingPhysics();
        biasPhys->Bias("deuteron");
        biasPhys->Bias("triton");
        pl->RegisterPhysics(biasPhys);
    }
}


class KillboxSteppingAction : public G4UserSteppingAction
{
public:
    void UserSteppingAction(const G4Step* step) override
    {
        auto* track = step->GetTrack();

        if (gKeepOnlyDTpnEnabled.load(std::memory_order_relaxed))
        {
            auto* def = track->GetParticleDefinition();

            const bool keep =
                (def == G4Neutron::NeutronDefinition()) ||
                (def == G4Proton::ProtonDefinition()) ||
                (def == G4Deuteron::DeuteronDefinition()) ||
                (def == G4Triton::TritonDefinition());

            if (!keep) {
                track->SetTrackStatus(fStopAndKill);
                return;
            }
        }

        const auto* prePV = step->GetPreStepPoint() ? step->GetPreStepPoint()->GetPhysicalVolume() : nullptr;
        const auto* postPV = step->GetPostStepPoint() ? step->GetPostStepPoint()->GetPhysicalVolume() : nullptr;

        // Only do anything if we actually entered a different volume
        if (!postPV || postPV == prePV)
            return;

        bool enteredKillbox = false;


        // Material naming rule: "Killbox"
        if (!enteredKillbox)
        {
            const auto* postMat = step->GetPostStepPoint()->GetMaterial();
            if (postMat && postMat->GetName() == "Killbox")
                enteredKillbox = true;
        }

        if (enteredKillbox)
        {
            track->SetTrackStatus(fStopAndKill);
            return;
        }
    }
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <G4VUserActionInitialization.hh>
#include <G4PhysListFactory.hh>
#include <G4PhysListFactoryMessenger.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithABool.hh>
#include <G4RunManagerFactory.hh>

class Action : public G4VUserActionInitialization, public G4UImessenger
{
private:
    G4UIcmdWithAString* fCmdPhys; ///< macro cmd to select a physics list
    G4UIcmdWithABool* fCmdKeepOnly; ///< enable n/p/d/t filter
    G4UIdirectory* fDirBias = nullptr; ///< enable isotropic biasing (WIP)
    G4UIcmdWithABool* fCmdDTXSBoost = nullptr; ///< enable dt biasing
public:
    Action() : G4VUserActionInitialization(), G4UImessenger() {
        fCmdPhys = new G4UIcmdWithAString("/physics_lists/select", this);
        fCmdPhys->SetGuidance("Select a physics list");
        fCmdPhys->SetGuidance("Candidates are specified in G4PhysListFactory.cc");
        fCmdPhys->SetParameterName("name of a physics list", false);
        fCmdPhys->AvailableForStates(G4State_PreInit);
        fCmdKeepOnly = new G4UIcmdWithABool("/particle/KeepOnlyNG", this);
        fCmdKeepOnly->SetGuidance("true: kill everything except neutron/proton/deuteron/triton. false: disable.");
        fCmdDTXSBoost = new G4UIcmdWithABool("/particle/enableDTXSBoost", this);
        fCmdDTXSBoost->SetGuidance("Enable GB01 D/T cross-section boost everywhere. Must be set before /run/initialize.");
        fCmdDTXSBoost->AvailableForStates(G4State_PreInit);
    }
    ~Action() {
        delete fCmdDTXSBoost;
        delete fDirBias;
        delete fCmdPhys;
        delete fCmdKeepOnly;
    }
    void Build() const {
        SetUserAction(new Generator);
        SetUserAction(new KillboxSteppingAction);
    }
    void SetNewValue(G4UIcommand* cmd, G4String value) {
        if (cmd == fCmdKeepOnly) {
            gKeepOnlyDTpnEnabled.store(fCmdKeepOnly->GetNewBoolValue(value), std::memory_order_relaxed);
            return;
        }
        if (cmd == fCmdDTXSBoost) {
            gEnableDTXSBoost.store(fCmdDTXSBoost->GetNewBoolValue(value), std::memory_order_relaxed);
            return;
        }
        if (cmd != fCmdPhys) return;
        auto run = G4RunManager::GetRunManager();
        G4PhysListFactory factory;
        G4VModularPhysicsList* physicsList1 = factory.GetReferencePhysList(value);
        RegisterDTBiasingPhysics(physicsList1);
        run->SetUserInitialization(physicsList1);
    } ///< for UI
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <G4RunManagerFactory.hh>
#include <G4AnalysisManager.hh>
#include <G4ScoringManager.hh>
#include <G4VisExecutive.hh>
#include <G4UIExecutive.hh>
#include <G4VisManager.hh>
#include <G4UImanager.hh>

int main(int argc,char** argv)
{
	auto *run = G4RunManagerFactory::CreateRunManager();

	G4ScoringManager::GetScoringManager(); // enable macro commands in /score/


    G4PhysListFactory factory;
    auto* pl0 = factory.ReferencePhysList();
    if (auto* modular = dynamic_cast<G4VModularPhysicsList*>(pl0)) {
        RegisterDTBiasingPhysics(modular);
    }
    run->SetUserInitialization(pl0);

	run->SetUserInitialization(new Detector);

	run->SetUserInitialization(new Action);



	G4VisManager *vis = new G4VisExecutive("quiet"); vis->Initialize();

	if (argc==1) { // interactive mode
		//G4UIExecutive ui(argc, argv);  TO FIX AV FALSE POSITIVE ON KEYBOARD READS ON UI MODE
		//ui.SessionStart();
	} else { // batch mode
		G4String command = "/control/execute ";
		G4UImanager::GetUIpointer()->ApplyCommand(command+argv[1]);
	}

	delete vis; delete run;
	return 0;
}
