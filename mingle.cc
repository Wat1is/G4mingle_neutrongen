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

    enum class EMFieldMode {
        UniformOnly = 0,   // use uniform E only
        RadialE = 1    // add radial E contribution (plus uniform E offset)
    };
    G4ThreeVector gFieldB = G4ThreeVector(0., 0., 0.) * tesla;
    G4ThreeVector gFieldE = G4ThreeVector(0., 0., 0.) * (volt / m);
    EMFieldMode gEMFieldMode = EMFieldMode::UniformOnly;
    // Radial E parameters:
    // Er = U / (ln(R2/R1) * r) within (r in [Rmin,Rmax]) and |z|<=Zmax, directed inward in xy-plane
    G4double gRadialU = 8.0e7 * (volt / m);  // previously "80000000 * volt/m"
    G4double gRadialR1 = 20.0 * mm;
    G4double gRadialR2 = 96.0 * mm;
    G4double gRadialRmin = 10.0 * mm;
    G4double gRadialRmax = 48.0 * mm;
    G4double gRadialZmax = 180.0 * mm;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ParamEMField : public G4ElectroMagneticField
{
public:
    ParamEMField() : G4ElectroMagneticField() {}
    ~ParamEMField() override = default;

    G4bool DoesFieldChangeEnergy() const override { return true; }

    // field[6] = {Bx,By,Bz,Ex,Ey,Ez}
    void GetFieldValue(const G4double Point[4], G4double* field) const override
    {
        // Uniform B
        field[0] = gFieldB.x();
        field[1] = gFieldB.y();
        field[2] = gFieldB.z();

        // Start with uniform E
        G4double Ex = gFieldE.x();
        G4double Ey = gFieldE.y();
        G4double Ez = gFieldE.z();

        if (gEMFieldMode == EMFieldMode::RadialE)
        {
            const G4double x = Point[0];
            const G4double y = Point[1];
            const G4double z = Point[2];

            const G4double r = std::sqrt(x * x + y * y);
            const G4double az = std::fabs(z);

            // avoid r=0 and invalid log
            if (r > 0.0 && gRadialR2 > gRadialR1 && gRadialR1 > 0.0)
            {
                if (r >= gRadialRmin && r <= gRadialRmax && az <= gRadialZmax)
                {
                    const G4double denom = std::log(gRadialR2 / gRadialR1) * r;
                    if (denom != 0.0)
                    {
                        const G4double Er = gRadialU / denom;  // units: (V/m)/(1* m) -> V/m

                        // inward radial direction in xy-plane
                        const G4double cx = x / r;
                        const G4double cy = y / r;

                        Ex += -Er * cx;
                        Ey += -Er * cy;
                    }
                }
            }
        }

        field[3] = Ex;
        field[4] = Ey;
        field[5] = Ez;
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

                static thread_local ParamEMField* field = nullptr;
                static thread_local G4EqMagElectricField* equation = nullptr;
                static thread_local G4DormandPrince745* stepper = nullptr;
                static thread_local G4MagInt_Driver* driver = nullptr;
                static thread_local G4ChordFinder* chordFinder = nullptr;

                // Local manager for targeted PV case
                static thread_local G4FieldManager* localFM = nullptr;

                if (!fieldAttached)
                {
                    field = new ParamEMField();
                    equation = new G4EqMagElectricField(field);
                    stepper = new G4DormandPrince745(equation, 8);

                    driver = new G4MagInt_Driver(gEMFieldMinStep, stepper, stepper->GetNumberOfVariables());
                    chordFinder = new G4ChordFinder(driver);

                    if (gEMFieldTargetPVName.empty())
                    {
                        // Global field manager (affects whole geometry)
                        auto* globalFM = G4TransportationManager::GetTransportationManager()->GetFieldManager();
                        globalFM->SetDetectorField(field);
                        globalFM->SetChordFinder(chordFinder);

                        G4cout << "[Field] Global EM field enabled. "
                            << "B(T)=(" << gFieldB.x() / tesla << "," << gFieldB.y() / tesla << "," << gFieldB.z() / tesla << ") "
                            << "E(V/m)=(" << gFieldE.x() / (volt / m) << "," << gFieldE.y() / (volt / m) << "," << gFieldE.z() / (volt / m) << ") "
                            << "minStep=" << gEMFieldMinStep / mm << " mm"
                            << G4endl;
                    }
                    else
                    {
                        // Targeted field manager (only inside volumes with matching PHYSICAL name)
                        localFM = new G4FieldManager();
                        localFM->SetDetectorField(field);
                        localFM->SetChordFinder(chordFinder);

                        bool foundAny = false;
                        for (auto* pv : *G4PhysicalVolumeStore::GetInstance())
                        {
                            if (!pv) continue;
                            if (pv->GetName() != gEMFieldTargetPVName) continue;

                            foundAny = true;
                            pv->GetLogicalVolume()->SetFieldManager(localFM, true); // force to daughters
                        }

                        if (!foundAny)
                        {
                            G4cout << "[Field] WARNING: no PV found with name '" << gEMFieldTargetPVName
                                << "'. No targeted EM field attached." << G4endl;
                        }
                        else
                        {
                            G4cout << "[Field] Targeted EM field enabled for PV name '" << gEMFieldTargetPVName
                                << "'. minStep=" << gEMFieldMinStep / mm << " mm" << G4endl;
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
#include <G4UIcmdWith3VectorAndUnit.hh>


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


    G4UIcmdWith3VectorAndUnit* fCmdEMFieldB = nullptr; // EM Field Params start
    G4UIcmdWith3VectorAndUnit* fCmdEMFieldE = nullptr;

    G4UIcmdWithAString* fCmdEMFieldMode = nullptr;

    G4UIcmdWithADoubleAndUnit* fCmdRadialU = nullptr;
    G4UIcmdWithADoubleAndUnit* fCmdRadialR1 = nullptr;
    G4UIcmdWithADoubleAndUnit* fCmdRadialR2 = nullptr;
    G4UIcmdWithADoubleAndUnit* fCmdRadialRmin = nullptr;
    G4UIcmdWithADoubleAndUnit* fCmdRadialRmax = nullptr;
    G4UIcmdWithADoubleAndUnit* fCmdRadialZmax = nullptr; // EM Field params finish
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

        fCmdEMFieldB = new G4UIcmdWith3VectorAndUnit("/particle/EMFieldB", this);
        fCmdEMFieldB->SetGuidance("Uniform magnetic field vector. Example: /particle/EMFieldB 0 0 1 tesla");
        fCmdEMFieldB->SetUnitCategory("Magnetic flux density");
        fCmdEMFieldB->SetDefaultUnit("tesla");
        fCmdEMFieldB->AvailableForStates(G4State_PreInit);

        fCmdEMFieldE = new G4UIcmdWith3VectorAndUnit("/particle/EMFieldE", this);
        fCmdEMFieldE->SetGuidance("Uniform electric field vector. Example: /particle/EMFieldE 0 0 1 volt/m");
        fCmdEMFieldE->SetUnitCategory("Electric field");
        fCmdEMFieldE->SetDefaultUnit("volt/m");
        fCmdEMFieldE->AvailableForStates(G4State_PreInit);

        fCmdEMFieldMode = new G4UIcmdWithAString("/particle/EMFieldMode", this);
        fCmdEMFieldMode->SetGuidance("Field mode: 'uniform' or 'radialE'. PreInit only.");
        fCmdEMFieldMode->SetParameterName("mode", false);
        fCmdEMFieldMode->AvailableForStates(G4State_PreInit);

        // Radial-E parameters (all PreInit only)
        fCmdRadialU = new G4UIcmdWithADoubleAndUnit("/particle/EMFieldRadialU", this);
        fCmdRadialU->SetGuidance("Radial E scale U (units Electric field).");
        fCmdRadialU->SetUnitCategory("Electric field");
        fCmdRadialU->SetDefaultUnit("volt/m");
        fCmdRadialU->AvailableForStates(G4State_PreInit);

        fCmdRadialR1 = new G4UIcmdWithADoubleAndUnit("/particle/EMFieldRadialR1", this);
        fCmdRadialR1->SetGuidance("Radial model R1.");
        fCmdRadialR1->SetUnitCategory("Length");
        fCmdRadialR1->SetDefaultUnit("mm");
        fCmdRadialR1->AvailableForStates(G4State_PreInit);

        fCmdRadialR2 = new G4UIcmdWithADoubleAndUnit("/particle/EMFieldRadialR2", this);
        fCmdRadialR2->SetGuidance("Radial model R2.");
        fCmdRadialR2->SetUnitCategory("Length");
        fCmdRadialR2->SetDefaultUnit("mm");
        fCmdRadialR2->AvailableForStates(G4State_PreInit);

        fCmdRadialRmin = new G4UIcmdWithADoubleAndUnit("/particle/EMFieldRadialRmin", this);
        fCmdRadialRmin->SetGuidance("Radial model active region minimum r.");
        fCmdRadialRmin->SetUnitCategory("Length");
        fCmdRadialRmin->SetDefaultUnit("mm");
        fCmdRadialRmin->AvailableForStates(G4State_PreInit);

        fCmdRadialRmax = new G4UIcmdWithADoubleAndUnit("/particle/EMFieldRadialRmax", this);
        fCmdRadialRmax->SetGuidance("Radial model active region maximum r.");
        fCmdRadialRmax->SetUnitCategory("Length");
        fCmdRadialRmax->SetDefaultUnit("mm");
        fCmdRadialRmax->AvailableForStates(G4State_PreInit);

        fCmdRadialZmax = new G4UIcmdWithADoubleAndUnit("/particle/EMFieldRadialZmax", this);
        fCmdRadialZmax->SetGuidance("Radial model active region |z| <= Zmax.");
        fCmdRadialZmax->SetUnitCategory("Length");
        fCmdRadialZmax->SetDefaultUnit("mm");
        fCmdRadialZmax->AvailableForStates(G4State_PreInit);
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
        delete fCmdEMFieldB;
        delete fCmdEMFieldE;
        delete fCmdEMFieldMode;
        delete fCmdRadialU;
        delete fCmdRadialR1;
        delete fCmdRadialR2;
        delete fCmdRadialRmin;
        delete fCmdRadialRmax;
        delete fCmdRadialZmax;
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
        if (cmd == fCmdEMFieldB) {
            gFieldB = fCmdEMFieldB->GetNew3VectorValue(value);
            return;
        }
        if (cmd == fCmdEMFieldE) {
            gFieldE = fCmdEMFieldE->GetNew3VectorValue(value);
            return;
        }
        if (cmd == fCmdEMFieldMode) {
            if (value == "radialE" || value == "radial") gEMFieldMode = EMFieldMode::RadialE;
            else gEMFieldMode = EMFieldMode::UniformOnly;
            return;
        }
        if (cmd == fCmdRadialU) { gRadialU = fCmdRadialU->GetNewDoubleValue(value);   return; }
        if (cmd == fCmdRadialR1) { gRadialR1 = fCmdRadialR1->GetNewDoubleValue(value);  return; }
        if (cmd == fCmdRadialR2) { gRadialR2 = fCmdRadialR2->GetNewDoubleValue(value);  return; }
        if (cmd == fCmdRadialRmin) { gRadialRmin = fCmdRadialRmin->GetNewDoubleValue(value);return; }
        if (cmd == fCmdRadialRmax) { gRadialRmax = fCmdRadialRmax->GetNewDoubleValue(value);return; }
        if (cmd == fCmdRadialZmax) { gRadialZmax = fCmdRadialZmax->GetNewDoubleValue(value);return; }
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
