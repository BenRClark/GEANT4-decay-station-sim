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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "gedssd_bcPhysicsList.hh"
#include "gedssd_bcPhysicsListMessenger.hh"

//#include "gedssd_bcPhysListParticles.hh"
#include "G4EmStandardPhysics.hh"
//#include "gedssd_bcPhysListEmLowEnergy.hh"
//#include "gedssd_bcPhysListHadron.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BIC_HP.hh"

#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4QStoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
//the next three were added, as they were also included in Sean's file
#include "G4RadioactiveDecay.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4DecayPhysics.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

gedssd_bcPhysicsList::gedssd_bcPhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  //defaultCutValue = 0.1*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  DetectorCuts = 0;
  TargetCuts   = 0;

  pMessenger = new gedssd_bcPhysicsListMessenger(this);

  SetVerboseLevel(1);

  //default physics
  particleList = new G4DecayPhysics();

  //default physics
  raddecayList = new G4RadioactiveDecayPhysics();

  // EM physics
  emPhysicsList = new G4EmStandardPhysics();
  
  // Had physics 
  hadPhysicsList = 0;
  nhadcomp = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

gedssd_bcPhysicsList::~gedssd_bcPhysicsList()
{
  delete pMessenger;
  delete raddecayList;
  delete emPhysicsList;
  if (hadPhysicsList) delete hadPhysicsList;
  if (nhadcomp > 0) {
    for(G4int i=0; i<nhadcomp; i++) {
      delete hadronPhys[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::ConstructParticle()
{
  particleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::ConstructProcess()
{
  AddTransportation();
  // em
  emPhysicsList->ConstructProcess();
  // decays
  particleList->ConstructProcess();
  raddecayList->ConstructProcess();//this was commented out in Sean's version
  
// had
  if (nhadcomp > 0) {
    for(G4int i=0; i<nhadcomp; i++) {
      (hadronPhys[i])->ConstructProcess();
    }
  }

  if (hadPhysicsList) hadPhysicsList->ConstructProcess();
  G4cout << "### gedssd_bcPhysicsList::ConstructProcess is done" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::SelectPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "gedssd_bcPhysicsList::SelectPhysicsList: <" << name << ">" << G4endl;
  }
  // default  Had physics
  if (name == "Hadron" && !hadPhysicsList) {
    // hadPhysicsList = new gedssd_bcPhysListHadron("hadron");
  } else if (name == "QGSP_BERT") {
    AddExtraBuilders(false);
    hadPhysicsList = new HadronPhysicsQGSP_BERT("std-hadron");
  } else if (name == "QGSP_BIC" && !hadPhysicsList) {
    AddExtraBuilders(false);
    hadPhysicsList = new HadronPhysicsQGSP_BIC("std-hadron");
  } else if (name == "QGSP_BERT_HP"  && !hadPhysicsList) {
    AddExtraBuilders(true);
    hadPhysicsList = new HadronPhysicsQGSP_BERT_HP("std-hadron");
  } else if (name == "QGSP_BIC_HP"  && !hadPhysicsList) {
    AddExtraBuilders(true);
    hadPhysicsList = new HadronPhysicsQGSP_BIC_HP("std-hadron");
    //} else if (name == "LowEnergy_EM") {
  } else if (name == "LowEnergy_EM_Livermore") {
      delete emPhysicsList;
      // emPhysicsList = new gedssd_bcPhysListEmLowEnergy("lowe-em");
      emPhysicsList = new G4EmLivermorePhysics();
  } else if (name == "Standard_EM") {
      delete emPhysicsList;
      emPhysicsList = new G4EmStandardPhysics();
  } else if(name == "LowEnergy_EM_Penelope"){
      delete emPhysicsList;
      emPhysicsList = new G4EmPenelopePhysics();
  } else {
      G4cout << "gedssd_bcPhysicsList WARNING wrong or unkonwn <" 
	     << name << "> Physics " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::AddExtraBuilders(G4bool flagHP)
{
  nhadcomp = 5;
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronElasticPhysics("elastic",verboseLevel,
						   flagHP));
  hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("ionBIC"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::SetCuts()
{
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  G4cout << "world cuts are set" << G4endl;

//   if( !TargetCuts ) SetTargetCut(cutForElectron);
//   G4Region* region = (G4RegionStore::GetInstance())->GetRegion("Target");
//   region->SetProductionCuts(TargetCuts);
//   G4cout << "Target cuts are set" << G4endl;

 G4Region* region = (G4RegionStore::GetInstance())->GetRegion("GeThickDetector"
);
  if( !DetectorCuts ) SetDetectorCut(cutForElectron);
  //region = (G4RegionStore::GetInstance())->GetRegion("Detector");
  region->SetProductionCuts(DetectorCuts);
  G4cout << "Detector cuts are set" << G4endl;

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::SetTargetCut(G4double cut)
{
  if( !TargetCuts ) TargetCuts = new G4ProductionCuts();

  TargetCuts->SetProductionCut(cut, idxG4GammaCut);
  TargetCuts->SetProductionCut(cut, idxG4ElectronCut);
  TargetCuts->SetProductionCut(cut, idxG4PositronCut);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void gedssd_bcPhysicsList::SetDetectorCut(G4double cut)
{
  if( !DetectorCuts ) DetectorCuts = new G4ProductionCuts();

  DetectorCuts->SetProductionCut(cut, idxG4GammaCut);
  DetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
  DetectorCuts->SetProductionCut(cut, idxG4PositronCut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
