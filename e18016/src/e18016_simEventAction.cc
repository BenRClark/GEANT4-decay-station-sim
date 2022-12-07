//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// CHANGE HISTORY
// --------------
#include "G4ios.hh"
#include "e18016_simEventActionMessenger.hh"
#include "e18016_simEventAction.hh"
#include "e18016_simDetectorHit.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "e18016_simAnalysisManager.hh"
#include "G4SDManager.hh"


extern G4bool drawEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

e18016_simEventAction::e18016_simEventAction()
  : drawFlag("all")
{
  eventMessenger = new e18016_simEventActionMessenger(this);
  //  G4cout<<"gedssdEventAction"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

e18016_simEventAction::~e18016_simEventAction()
{
  //  G4cout<<"~e18016_simEventAction"<<G4endl;
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e18016_simEventAction::BeginOfEventAction(const G4Event* evt )
{
  //   G4cout<<"BeginOfEventAction"<<G4endl;
  if(ge18016_simAnalysisManager) ge18016_simAnalysisManager->BeginOfEventAction(evt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e18016_simEventAction::EndOfEventAction(const G4Event* evt)
{
  //  G4cout<<"EndOfEventAction"<<G4endl;

//   G4SDManager *SDman = G4SDManager::GetSDMpointer();
//   G4String colNam;
//   G4int DetectorColID = SDman->GetCollectionID(colNam="DetectorCollection");
//   G4HCofThisEvent *HCE = evt->GetHCofThisEvent();
//   e18016_simDetectorHitsCollection *DHC = 0;

//   if(HCE){
//     DHC = (e18016_simDetectorHitsCollection*)(HCE->GetHC(DetectorColID));
//   }

//   if(DHC){
//     int n_hit = DHC->entries();

//     for(int i=0; i<n_hit; i++){
//       G4String name = (*DHC)[i]->GetParticlename();
//       G4String parentname = (*DHC)[i]->GetParentName();

//       G4cout << "This particle - " << name << " - came from " << parentname << G4endl;
//     }

//   }
  if(ge18016_simAnalysisManager) ge18016_simAnalysisManager->EndOfEventAction(evt);

  G4int event_id       = evt->GetEventID();
  //analysis
  //e18016_simAnalysisManager::getInstance()->EndOfEvent();
  // visualisation
  if (event_id < 100 && G4VVisManager::GetConcreteInstance()) {
    G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
    for (G4int i=0; i<n_trajectories; i++) {
      G4Trajectory* trj = (G4Trajectory *)
        ((*(evt->GetTrajectoryContainer()))[i]);
      if (drawFlag == "all") {
        trj->DrawTrajectory();
      } else if (drawFlag == "charged" && trj->GetCharge() != 0.) {
        trj->DrawTrajectory();
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

































