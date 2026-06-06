#ifndef SphereScoreWriter_hh
#define SphereScoreWriter_hh 1

#include "G4VScoreWriter.hh"

class SphereScoringMesh;

class SphereScoreWriter : public G4VScoreWriter
{
public:
    SphereScoreWriter() = default;
    ~SphereScoreWriter() override = default;

    void DumpQuantityToFile(const G4String& psName,
                            const G4String& fileName,
                            const G4String& option) override;
    void DumpAllQuantitiesToFile(const G4String& fileName,
                                 const G4String& option) override;

private:
    void DumpSphereQuantityToStream(SphereScoringMesh* mesh,
                                    const G4String& psName,
                                    std::ostream& out);
};

#endif
