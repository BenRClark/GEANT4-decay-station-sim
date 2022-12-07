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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "e17011_simRunAction.hh"
#include "e17011_simAnalysisManager.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

//#include "G4MPImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

e17011_simRunAction::e17011_simRunAction()
{
  // G4cout<<"e17011_simRunAction"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

e17011_simRunAction::~e17011_simRunAction()
{
  // G4cout<<"~e17011_simRunAction"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simRunAction::BeginOfRunAction(const G4Run* aRun)
{
  // Creation of the analysis manager
  //e17011_simAnalysisManager* analysis = e17011_simAnalysisManager::getInstance();
  //analysis->BeginOfRun();
  //G4int rank = G4MPImanager::GetManager()->GetRank();
  // G4cout<<"e17011_simRunAction BeginOfRunAction"<<G4endl;
  G4int rank = 0;
  if(ge17011_simAnalysisManager) ge17011_simAnalysisManager->BeginOfRunAction(aRun,rank);
  G4int RunN = aRun->GetRunID();
  if ( RunN % 1000 == 0 ) 
    G4cout << "### Run : " << RunN << G4endl;

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/clear/view");
      UI->ApplyCommand("/vis/draw/current");
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simRunAction::EndOfRunAction(const G4Run* aRun)
{
  //  G4cout<<"e17011_simRunAction EndOfRunAction"<<G4endl;
  // Get the analysis manager
  //e17011_simAnalysisManager* analysis = e17011_simAnalysisManager::getInstance();
  //analysis->EndOfRun();
  //G4int rank = G4MPImanager::GetManager()->GetRank();
  G4int rank = 0;
  if(ge17011_simAnalysisManager) ge17011_simAnalysisManager->EndOfRunAction(aRun,rank);
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/show/view");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




