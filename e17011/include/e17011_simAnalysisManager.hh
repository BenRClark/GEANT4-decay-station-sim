//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//

#ifndef e17011_simAnalysisManager_h
#define e17011_simAnalysisManager_h 1


class G4VPhysicalVolume;
class G4Event;
class G4Run;
class G4Track;
class G4Step;
class e17011_simDetectorConstruction;

#include "G4ClassificationOfNewTrack.hh"
#include "G4TrackStatus.hh"
#include "G4Types.hh"
#include "globals.hh"


class e17011_simAnalysisManager;
extern e17011_simAnalysisManager *ge17011_simAnalysisManager; // global e17011_simAnalysisManager

class e17011_simAnalysisManager {
  
public:
  e17011_simAnalysisManager() {
    if (ge17011_simAnalysisManager)
      delete ge17011_simAnalysisManager;
    ge17011_simAnalysisManager = this;
  }
  
  virtual ~e17011_simAnalysisManager() {
    if (ge17011_simAnalysisManager == this)
      ge17011_simAnalysisManager = (e17011_simAnalysisManager *)0;
  }
  
  static e17011_simAnalysisManager *GetAnalysisManager() {
    return ge17011_simAnalysisManager;
  }
  
public:


  // G4VUserDetectorConstruction
  virtual void Construct(const G4VPhysicalVolume */*theWorldWolume*/) {;}
  
  // G4VUserPhysicsList
  virtual void ConstructParticle() {;}
  virtual void ConstructProcess() {;}
  virtual void SetCuts() {;}
  
  // G4VUserPrimaryGeneratorAction
  virtual void GeneratePrimaries(const G4Event */*anEvent*/, const G4double /*beam*/) {;}
  virtual void ResetStartLocation(G4double &, G4double &) {;}
  
  virtual void DetectorInfo(e17011_simDetectorConstruction*) {;}

  // G4UserRunAction
  virtual void BeginOfRunAction(const G4Run */*aRun*/, G4int) {;}
  virtual void EndOfRunAction(const G4Run */*aRun*/, G4int) {;}
  
  // G4UserEventAction
  virtual void BeginOfEventAction(const G4Event */*anEvent*/) {;}
  virtual void EndOfEventAction(const G4Event */*anEvent*/) {;}
  
  // G4UserStackingAction
  virtual void ClassifyNewTrack(
		   const G4Track */*aTrack*/,
		   G4ClassificationOfNewTrack */*classification*/) {;}
  virtual void NewStage() {;}
  virtual void PrepareNewEvent() {;}
  virtual void PrepareNewRun(G4double, G4double, G4double, G4double, G4int, G4int) {;}
  
  // G4UserTrackingAction
  virtual void PreUserTrackingAction(const G4Track */*aTrack*/) {;}
  virtual void PostUserTrackingAction(const G4Track */*aTrack*/,
				      G4TrackStatus */*status*/) {;}
  
  // G4UserSteppingAction
  virtual void UserSteppingAction(const G4Step */*aStep*/) {;}
  
  virtual void OpticalPhotons(G4int, G4double){;}
  virtual void PSF(G4int, G4double, G4double){;}
};

#endif

