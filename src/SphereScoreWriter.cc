#include "SphereScoreWriter.hh"

#include "SphereScoringMesh.hh"

#include "G4StatDouble.hh"
#include "G4SystemOfUnits.hh"
#include "G4VScoringMesh.hh"
#include "G4ios.hh"

#include <algorithm>
#include <fstream>
#include <iomanip>

void SphereScoreWriter::DumpQuantityToFile(const G4String& psName,
                                           const G4String& fileName,
                                           const G4String& option)
{
    auto* sphereMesh = dynamic_cast<SphereScoringMesh*>(fScoringMesh);
    if (!sphereMesh) {
        G4VScoreWriter::DumpQuantityToFile(psName, fileName, option);
        return;
    }

    G4String opt = option;
    std::transform(opt.begin(), opt.end(), opt.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    if (!opt.empty() && opt.find("csv") == std::string::npos) {
        G4cerr << "ERROR : DumpToFile : Unknown option for sphere mesh -> "
               << option << G4endl;
        return;
    }

    std::ofstream out(fileName);
    if (!out) {
        G4cerr << "ERROR : DumpToFile : File open error -> "
               << fileName << G4endl;
        return;
    }

    DumpSphereQuantityToStream(sphereMesh, psName, out);
}

void SphereScoreWriter::DumpAllQuantitiesToFile(const G4String& fileName,
                                                const G4String& option)
{
    auto* sphereMesh = dynamic_cast<SphereScoringMesh*>(fScoringMesh);
    if (!sphereMesh) {
        G4VScoreWriter::DumpAllQuantitiesToFile(fileName, option);
        return;
    }

    G4String opt = option;
    std::transform(opt.begin(), opt.end(), opt.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    if (!opt.empty() && opt.find("csv") == std::string::npos) {
        G4cerr << "ERROR : DumpToFile : Unknown option for sphere mesh -> "
               << option << G4endl;
        return;
    }

    std::ofstream out(fileName);
    if (!out) {
        G4cerr << "ERROR : DumpToFile : File open error -> "
               << fileName << G4endl;
        return;
    }

    const auto scoreMap = sphereMesh->GetScoreMap();
    for (const auto& entry : scoreMap) {
        out << "# primitive scorer name: " << entry.first << G4endl;
        DumpSphereQuantityToStream(sphereMesh, entry.first, out);
    }
}

void SphereScoreWriter::DumpSphereQuantityToStream(SphereScoringMesh* mesh,
                                                   const G4String& psName,
                                                   std::ostream& out)
{
    const auto scoreMap = mesh->GetScoreMap();
    const auto found = scoreMap.find(psName);
    if (found == scoreMap.end()) {
        G4cerr << "ERROR : DumpToFile : Unknown quantity, \""
               << psName << "\"." << G4endl;
        return;
    }

    auto* score = found->second->GetMap();
    const auto unitValue = mesh->GetPSUnitValue(psName);

    out << "theta_bin,phi_bin,theta_min_deg,theta_max_deg,phi_min_deg,phi_max_deg,value" << G4endl;
    out << std::setprecision(16);

    for (G4int thetaBin = 0; thetaBin < mesh->GetThetaBins(); ++thetaBin) {
        for (G4int phiBin = 0; phiBin < mesh->GetPhiBins(); ++phiBin) {
            const auto index = mesh->GetSphereIndex(thetaBin, phiBin);
            G4double value = 0.;
            const auto valueIt = score->find(index);
            if (valueIt != score->end() && valueIt->second != nullptr) {
                value = valueIt->second->sum_wx() / unitValue * fact;
            }

            out << thetaBin << ","
                << phiBin << ","
                << mesh->GetThetaMin(thetaBin) / deg << ","
                << mesh->GetThetaMax(thetaBin) / deg << ","
                << mesh->GetPhiMin(phiBin) / deg << ","
                << mesh->GetPhiMax(phiBin) / deg << ","
                << value << G4endl;
        }
    }

    out << std::setprecision(6);
}
