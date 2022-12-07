/*
attempt to add track information - snl
*/

#include "e17011_simTrackInformation.hh"
#include "G4ios.hh"

G4Allocator<e17011_simTrackInformation> aTrackInformationAllocator;

e17011_simTrackInformation::e17011_simTrackInformation()
{
  //  G4cout<<"e17011_simTrackInformation"<<G4endl;
    originalTrackID = 0;
    particleDefinition = 0;
    originalPosition = G4ThreeVector(0.,0.,0.);
    originalMomentum = G4ThreeVector(0.,0.,0.);
    originalEnergy = 0.;
    originalTime = 0.;
}

e17011_simTrackInformation::e17011_simTrackInformation(const G4Track* aTrack)
{
  // G4cout<<"e17011_simTrackInformation2"<<G4endl;
    originalTrackID = aTrack->GetTrackID();
    particleDefinition = aTrack->GetDefinition();
    originalPosition = aTrack->GetPosition();
    originalMomentum = aTrack->GetMomentum();
    originalEnergy = aTrack->GetTotalEnergy();
    originalTime = aTrack->GetGlobalTime();
}

e17011_simTrackInformation::e17011_simTrackInformation(const e17011_simTrackInformation* aTrackInfo)
{
  //G4cout<<"e17011_simTrackInformation3"<<G4endl;
    originalTrackID = aTrackInfo->originalTrackID;
    particleDefinition = aTrackInfo->particleDefinition;
    originalPosition = aTrackInfo->originalPosition;
    originalMomentum = aTrackInfo->originalMomentum;
    originalEnergy = aTrackInfo->originalEnergy;
    originalTime = aTrackInfo->originalTime;
}

e17011_simTrackInformation::~e17011_simTrackInformation(){;
  // G4cout<<"~e17011_simTrackInformation"<<G4endl;
}

void e17011_simTrackInformation::Print() const
{
    G4cout 
     << "Original track ID " << originalTrackID 
     << " at " << originalPosition << G4endl;
}
