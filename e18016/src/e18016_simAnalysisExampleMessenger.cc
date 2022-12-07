//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// AnalysisExampleMessenger program
// --------------------------------------------------------------

#include "e18016_simAnalysisExampleMessenger.hh"

#include "e18016_simAnalysisExample.hh"
#include "e18016_simEventActionMessenger.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

e18016_simAnalysisExampleMessenger::e18016_simAnalysisExampleMessenger
   (e18016_simAnalysisExample* SA):analysisexample(SA) {
  // cout<<"e18016_simAnalysisExampleMessenger"<<endl;

     analysisDir = new G4UIdirectory("/e18016_sim/analysis/");
     analysisDir->SetGuidance("event control");

     trackparentCmd = new G4UIcmdWithAnInteger("/e18016_sim/analysis/trackparent",this);
     trackparentCmd->SetGuidance("Z of parent ion to be tracked");
     trackparentCmd->SetParameterName("parentZ",false);
     trackparentCmd->SetRange("parentZ>=0");
     trackparentCmd->AvailableForStates(G4State_Idle);

     trackdaughterCmd = new G4UIcmdWithAnInteger("/e18016_sim/analysis/trackdaughter",this);
     trackdaughterCmd->SetGuidance("Z of daughter ion to be tracked");
     trackdaughterCmd->SetParameterName("daughterZ",false);
     trackdaughterCmd->SetRange("daughterZ>=0");
     trackdaughterCmd->AvailableForStates(G4State_Idle);

     filenameCmd = new G4UIcmdWithAString("/e18016_sim/analysis/filename",this);
     filenameCmd->SetGuidance("Analysis file name");
     filenameCmd->SetParameterName("choice",false);
     filenameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

     dirnameCmd = new G4UIcmdWithAString("/e18016_sim/analysis/dirname",this);
     dirnameCmd->SetGuidance("Directory for output files");
     dirnameCmd->SetParameterName("dirname",false);
     dirnameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}


//ooooooooooooooooooooooooooooooooooooooooo
e18016_simAnalysisExampleMessenger::~e18016_simAnalysisExampleMessenger() {
  //  cout<<"~e18016_simAnalysisExampleMessenger"<<endl;

  delete trackparentCmd;
  delete trackdaughterCmd;
  delete filenameCmd;
  delete dirnameCmd;
}


//ooooooooooooooooooooooooooooooooooooooooo
void e18016_simAnalysisExampleMessenger::SetNewValue(G4UIcommand* command, 
  G4String newValue) {

//   if(command == trackparentCmd)
//     analysisexample->SetparentZ(trackparentCmd->GetNewIntValue(newValue));

//    if(command == trackdaughterCmd)
//      analysisexample->SetdaughterZ(trackdaughterCmd->GetNewIntValue(newValue));

   if(command == filenameCmd)
     analysisexample->Setfilename(newValue);

   if(command == dirnameCmd)
     analysisexample->Setdirname(newValue);
}



