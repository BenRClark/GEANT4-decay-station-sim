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
#include "G4RunManager.hh"
#include "G4UImanager.hh"
// #ifdef G4UI_USE_GAG
// #include "G4UIGAG.hh"
// #endif
// #include "G4UIterminal.hh"
// #include "G4UItcsh.hh"
// #ifdef G4UI_USE_XM
// #include "G4UIXm.hh"
// #endif 


#include "e17011_simDetectorConstruction.hh"
#include "e17011_simPhysicsList.hh"
#include "e17011_simEventAction.hh"
#include "e17011_simRunAction.hh"
#include "e17011_simTrackingAction.hh"
#include "e17011_simSteppingAction.hh"
#include "e17011_simPrimaryGeneratorAction.hh"
#include "e17011_simAnalysisManager.hh"
#include "e17011_simAnalysisExample.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{
  // random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // Creation of the analysis 
  e17011_simAnalysisExample* analysis = new e17011_simAnalysisExample();
  //e17011_simAnalysisManager* analysis = e17011_simAnalysisManager::getInstance();

  // set mandatory initialization classes

  e17011_simDetectorConstruction* Detector = new e17011_simDetectorConstruction;
  runManager->SetUserInitialization(Detector);
  runManager->SetUserInitialization(new e17011_simPhysicsList);

  G4cout << "driver done 1"  << G4endl;

  // set mandatory user action class
  runManager->SetUserAction(new e17011_simPrimaryGeneratorAction);
  G4cout << "driver done 2"  << G4endl;
  runManager->SetUserAction(new e17011_simRunAction);
  G4cout << "driver done 3"  << G4endl;
  runManager->SetUserAction(new e17011_simEventAction);
  G4cout << "driver done 4"  << G4endl;
  runManager->SetUserAction(new e17011_simTrackingAction); //added SNL
  G4cout << "driver done 5"  << G4endl;
  runManager->SetUserAction(new e17011_simSteppingAction);
  G4cout << "driver done 6"  << G4endl;
  //
  //
   G4UIsession* session=0;



  // Initialize G4 kernel
  //
  //********* NOTE: this happens in the macro files. If you initialize the run here, 
  //********* you are unable to add any additional physics lists items in a macro.
  //********* So, if you are planning on using macros, don't initialize here!!!
  //runManager->Initialize();

  

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac");
#else
    UImanager->ApplyCommand("/control/execute init.mac");
#endif
    ui->SessionStart();
    delete ui;
#endif
  }



















   
//    if (argc==1)   // Define UI session for interactive mode.
//      {
// #ifdef G4UI_USE_XM
//       session = new G4UIXm(argc,argv);
// #else           
// #ifdef G4UI_USE_TCSH
//       session = new G4UIterminal(new G4UItcsh);      
// #else
//       session = new G4UIterminal();
// #endif
// #endif
//     }
 
// #ifdef G4VIS_USE
//   // visualization manager
//   G4VisManager* visManager = new G4VisExecutive;
//   visManager->Initialize();
// #endif

//   // Initialize G4 kernel
//   // do this at run time so the geometry/physics can be changed
//   //  runManager->Initialize();

//   // get the pointer to the User Interface manager 
//   G4UImanager* UI = G4UImanager::GetUIpointer();  
//   if (session)   // Define UI session for interactive mode.
//     {
//       session->SessionStart();
//       delete session;
//     }
//   else           // Batch mode
//     { 
//       // G4String command = "/control/execute/vis.mac ";//does not let gamma.mac etc run
//       G4String command = "/control/execute ";
//       G4String fileName = argv[1];
//       UI->ApplyCommand(command+fileName);
//     }
  
//   // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete analysis;
  delete runManager;

  return 0;
}








