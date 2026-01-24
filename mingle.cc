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
#include <G4PhysicalVolumeStore.hh>
#include <cmath>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>
#include <G4ChordFinder.hh>
#include <G4TransportationManager.hh>
#include <G4FieldManager.hh>
#include <G4ElectroMagneticField.hh>
#include <G4EqMagElectricField.hh>
#include <G4DormandPrince745.hh>
#include <G4MagIntegratorDriver.hh>

namespace {
    // ---- D/T XS boost (GB01) ----
    std::atomic<bool> gEnableDTXSBoost{ false };
    // Empty => global
    G4String gDTXSBoostPVName = "";

    // ---- EM field (F05Field) ----
    std::atomic<bool> gEnableEMField{ false };
    // Empty => global (TransportationManager FieldManager)
    G4String gEMFieldTargetPVName = "";
    // Min integration step for the field driver
    G4double gEMFieldMinStep = 0.01 * mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F05Field : public G4ElectroMagneticField
{
public:
    F05Field() : G4ElectroMagneticField() {}
    ~F05Field() override = default;

    G4bool DoesFieldChangeEnergy() const override { return true; }

    // field[6]: Bx,By,Bz,Ex,Ey,Ez
    // Point[4]: x,y,z,t
    void GetFieldValue(const G4double Point[4], G4double* Bfield) const override
    {
        // Radial electric field like in your ExamplesIn F05Field.cc
        G4double Ex = 0.0, Ey = 0.0;
        G4double R2 = 96.0;
        G4double R1 = 20.0;
        G4double U = 80000000 * volt / m; // same as example

        G4double posR = std::sqrt(Point[0] * Point[0] + Point[1] * Point[1]);
        G4double posZ = std::sqrt(Point[2] * Point[2]);

        if (posR > 0.0) {
            G4double Er = U / ((std::log(R2 / R1)) * posR);

            if (posR >= 10.0 && posR <= 48.0 && posZ <= 180.0) {
                G4double cos_theta = Point[0] / posR;
                G4double sin_theta = Point[1] / posR;
                Ex = -Er * cos_theta;
                Ey = -Er * sin_theta;
            }
        }

        // No magnetic field in this example
        Bfield[0] = 0.0;
        Bfield[1] = 0.0;
        Bfield[2] = 0.0;

        // Electric field
        Bfield[3] = Ex;
        Bfield[4] = Ey;
        Bfield[5] = 0.0;
    }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Detector : public G4VUserDetectorConstruction
{
	public:
		G4VPhysicalVolume* Construct() override {
			G4tgbVolumeMgr::GetInstance()->AddTextFile("detector.tg");
			return G4tgbVolumeMgr::GetInstance()->ReadAndConstructDetector();
		} ///< load detector definition from a text file "detector.tg"
        void ConstructSDandField() override
        {
            // ============================================================
            // (A) D/T XS boost operator (GB01) — attach once per thread
            // ============================================================
            if (gEnableDTXSBoost.load(std::memory_order_relaxed))
            {
                static thread_local bool biasAttached = false;
                static thread_local GB01BOptrMultiParticleChangeCrossSection* op = nullptr;

                if (!biasAttached)
                {
                    op = new GB01BOptrMultiParticleChangeCrossSection();
                    op->AddParticle("deuteron");
                    op->AddParticle("triton");

                    if (gDTXSBoostPVName.empty())
                    {
                        // Global: attach to ALL logical volumes
                        for (auto* lv : *G4LogicalVolumeStore::GetInstance()) {
                            if (lv) op->AttachTo(lv);
                        }
                    }
                    else
                    {
                        // Targeted: match ALL physical volumes with this name
                        bool foundAny = false;
                        for (auto* pv : *G4PhysicalVolumeStore::GetInstance()) {
                            if (!pv) continue;
                            if (pv->GetName() != gDTXSBoostPVName) continue;
                            foundAny = true;
                            op->AttachTo(pv->GetLogicalVolume());
                        }
                        if (!foundAny) {
                            G4cout << "[Biasing] WARNING: no physical volume found with name '"
                                << gDTXSBoostPVName << "'; no targeted D/T XS boost attached." << G4endl;
                        }
                    }

                    biasAttached = true;
                }
            }

            // ============================================================
            // (B) EM field (F05Field) — attach once per thread
            // ============================================================
            if (gEnableEMField.load(std::memory_order_relaxed))
            {
                static thread_local bool fieldAttached = false;

                static thread_local F05Field* field = nullptr;
                static thread_local G4EqMagElectricField* equation = nullptr;
                static thread_local G4DormandPrince745* stepper = nullptr;
                static thread_local G4MagInt_Driver* driver = nullptr;
                static thread_local G4ChordFinder* chordFinder = nullptr;

                // For targeted PV attachment
                static thread_local G4FieldManager* localFieldMgr = nullptr;

                if (!fieldAttached)
                {
                    field = new F05Field();

                    equation = new G4EqMagElectricField(field);

                    // ExamplesIn uses DormandPrince745 with nvar=8
                    stepper = new G4DormandPrince745(equation, 8);

                    driver = new G4MagInt_Driver(gEMFieldMinStep, stepper, stepper->GetNumberOfVariables());
                    chordFinder = new G4ChordFinder(driver);

                    if (gEMFieldTargetPVName.empty())
                    {
                        // Global field manager (same pattern as the ExamplesIn setup)
                        auto* globalFM = G4TransportationManager::GetTransportationManager()->GetFieldManager();
                        globalFM->SetDetectorField(field);
                        globalFM->SetChordFinder(chordFinder);

                        G4cout << "[Field] Enabled global EM field (F05Field). minStep="
                            << gEMFieldMinStep / mm << " mm" << G4endl;
                    }
                    else
                    {
                        // Targeted: attach to ALL physical volumes matching name
                        bool foundAny = false;

                        localFieldMgr = new G4FieldManager();
                        localFieldMgr->SetDetectorField(field);
                        localFieldMgr->SetChordFinder(chordFinder);

                        for (auto* pv : *G4PhysicalVolumeStore::GetInstance()) {
                            if (!pv) continue;
                            if (pv->GetName() != gEMFieldTargetPVName) continue;
                            foundAny = true;
                            pv->GetLogicalVolume()->SetFieldManager(localFieldMgr, true);
                        }

                        if (!foundAny) {
                            G4cout << "[Field] WARNING: no physical volume found with name '"
                                << gEMFieldTargetPVName << "'; no targeted EM field attached." << G4endl;
                        }
                        else {
                            G4cout << "[Field] Enabled targeted EM field (F05Field) for PV name '"
                                << gEMFieldTargetPVName << "'. minStep="
                                << gEMFieldMinStep / mm << " mm" << G4endl;
                        }
                    }

                    fieldAttached = true;
                }
            }
        }

};

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
#include <G4LogicalVolumeStore.hh>
#include <G4GenericBiasingPhysics.hh>
#include <G4VModularPhysicsList.hh>

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
#include <G4UIcmdWithADoubleAndUnit.hh>

class Action : public G4VUserActionInitialization, public G4UImessenger
{
private:
    G4UIcmdWithAString* fCmdPhys; ///< macro cmd to select a physics list
    G4UIcmdWithABool* fCmdKeepOnly; ///< enable n/p/d/t filter
    G4UIdirectory* fDirBias = nullptr; ///< enable isotropic biasing (WIP)
    G4UIcmdWithABool* fCmdDTXSBoost = nullptr; ///< enable dt biasing
    G4UIcmdWithAString* fCmdDTXSBoostVolume = nullptr; //< dt biasing target volume (defaults to world if none given)
    G4UIcmdWithABool* fCmdEnableEMField = nullptr; ///< enable EM field
    G4UIcmdWithAString* fCmdEMFieldTargetPV = nullptr; ///< EM field target (physical volume name,defaults to world if none given)
    G4UIcmdWithADoubleAndUnit* fCmdEMFieldMinStep = nullptr; ///<EM field min step
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

        fCmdDTXSBoostVolume = new G4UIcmdWithAString("/particle/DTXSBoostVolume", this);
        fCmdDTXSBoostVolume->SetGuidance("Set physical volume name to apply D/T XS boost only there. Empty = global.");
        fCmdDTXSBoostVolume->SetParameterName("pvName", false);
        fCmdDTXSBoostVolume->AvailableForStates(G4State_PreInit);

        fCmdEnableEMField = new G4UIcmdWithABool("/particle/EnableEMField", this);
        fCmdEnableEMField->SetGuidance("Enable EM field (F05Field). Must be set before /run/initialize.");
        fCmdEnableEMField->AvailableForStates(G4State_PreInit);

        fCmdEMFieldTargetPV = new G4UIcmdWithAString("/particle/EMFieldTargetPV", this);
        fCmdEMFieldTargetPV->SetGuidance("Physical volume name to attach EM field only there. Empty = global.");
        fCmdEMFieldTargetPV->SetParameterName("pvName", false);
        fCmdEMFieldTargetPV->AvailableForStates(G4State_PreInit);

        fCmdEMFieldMinStep = new G4UIcmdWithADoubleAndUnit("/particle/EMFieldMinStep", this);
        fCmdEMFieldMinStep->SetGuidance("Minimum integration step for EM field tracking (default 0.01 mm).");
        fCmdEMFieldMinStep->SetUnitCategory("Length");
        fCmdEMFieldMinStep->SetDefaultUnit("mm");
        fCmdEMFieldMinStep->AvailableForStates(G4State_PreInit);
    }
    ~Action() {
        delete fCmdDTXSBoost;
        delete fDirBias;
        delete fCmdPhys;
        delete fCmdKeepOnly;
        delete fCmdDTXSBoostVolume;
        delete fCmdEnableEMField;
        delete fCmdEMFieldTargetPV;
        delete fCmdEMFieldMinStep;
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
        if (cmd == fCmdDTXSBoostVolume) {
            gDTXSBoostPVName = value;   // empty string => global
            return;
        }
        if (cmd == fCmdEnableEMField) {
            gEnableEMField.store(fCmdEnableEMField->GetNewBoolValue(value), std::memory_order_relaxed);
            return;
        }
        if (cmd == fCmdEMFieldTargetPV) {
            gEMFieldTargetPVName = value; // empty => global
            return;
        }
        if (cmd == fCmdEMFieldMinStep) {
            gEMFieldMinStep = fCmdEMFieldMinStep->GetNewDoubleValue(value);
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
