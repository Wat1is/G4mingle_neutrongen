#include "PSSphereSurfaceCurrent.hh"

#include "SphereScoringMesh.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cmath>

PSSphereSurfaceCurrent::PSSphereSurfaceCurrent(const G4String& name,
                                               SphereScoringMesh* mesh,
                                               G4int direction)
    : G4VPrimitiveScorer(name), fMesh(mesh), fDirection(direction)
{
}

void PSSphereSurfaceCurrent::Initialize(G4HCofThisEvent* hce)
{
    fEvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
    if (fHCID < 0) {
        fHCID = GetCollectionID(0);
    }
    hce->AddHitsCollection(fHCID, fEvtMap);
}

void PSSphereSurfaceCurrent::clear()
{
    if (fEvtMap) fEvtMap->clear();
}

void PSSphereSurfaceCurrent::PrintAll()
{
    if (fEvtMap) fEvtMap->PrintAllHits();
}

G4bool PSSphereSurfaceCurrent::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    if (!fMesh || !step || !fEvtMap) return false;

    const auto* preStep = step->GetPreStepPoint();
    const auto* postStep = step->GetPostStepPoint();
    if (!preStep || !postStep) return false;

    const auto center = fMesh->GetSphereCenter();
    const auto radius = fMesh->GetSphereRadius();
    if (radius <= 0.) return false;

    const auto pre = preStep->GetPosition() - center;
    const auto post = postStep->GetPosition() - center;
    const auto preR = pre.mag();
    const auto postR = post.mag();

    const auto crossedOut = (preR < radius && postR >= radius);
    const auto crossedIn = (preR > radius && postR <= radius);

    if (fDirection == 1 && !crossedOut) return false;
    if (fDirection == -1 && !crossedIn) return false;
    if (fDirection == 0 && !crossedOut && !crossedIn) return false;

    G4ThreeVector crossing;
    if (!ComputeCrossingPoint(pre, post, crossing)) {
        crossing = (crossedOut ? post : pre);
        if (crossing.mag2() > 0.) {
            crossing = crossing.unit() * radius;
        }
    }

    const auto globalCrossing = crossing + center;
    const auto index = fMesh->GetSphereIndex(globalCrossing);
    const auto weight = step->GetTrack() ? step->GetTrack()->GetWeight() : 1.0;

    fEvtMap->add(index, weight);
    return true;
}

G4bool PSSphereSurfaceCurrent::ComputeCrossingPoint(const G4ThreeVector& pre,
                                                    const G4ThreeVector& post,
                                                    G4ThreeVector& crossing) const
{
    if (!fMesh) return false;

    const auto radius = fMesh->GetSphereRadius();
    const auto direction = post - pre;
    const auto a = direction.mag2();
    if (a <= 0.) return false;

    const auto b = 2.0 * pre.dot(direction);
    const auto c = pre.mag2() - radius * radius;
    const auto discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.) return false;

    const auto root = std::sqrt(discriminant);
    const auto invDenom = 1.0 / (2.0 * a);
    const auto t0 = (-b - root) * invDenom;
    const auto t1 = (-b + root) * invDenom;

    auto pick = -1.0;
    if (t0 >= 0.0 && t0 <= 1.0) pick = t0;
    if (t1 >= 0.0 && t1 <= 1.0 && (pick < 0.0 || t1 < pick)) pick = t1;
    if (pick < 0.0) return false;

    crossing = pre + pick * direction;
    return true;
}
