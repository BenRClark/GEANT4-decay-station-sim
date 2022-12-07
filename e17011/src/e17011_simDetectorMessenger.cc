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

#include "e17011_simDetectorMessenger.hh"

#include "e17011_simDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

e17011_simDetectorMessenger::e17011_simDetectorMessenger(e17011_simDetectorConstruction* myDet)
:myDetector(myDet)
{ 
  // G4cout<<"gedssdnrDetectorMessenger"<<G4endl;
  e17011_simDir = new G4UIdirectory("/e17011_sim/");
  e17011_simDir->SetGuidance("UI commands specific to this example.");
  
  detDir = new G4UIdirectory("/e17011_sim/det/");
  detDir->SetGuidance("detector control.");
  
  GeThickDetectorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/e17011_sim/det/setGeThickDetectorThickness",this);
  GeThickDetectorThicknessCmd->SetGuidance("Set the Thick Ge Detector Thickness.");
  GeThickDetectorThicknessCmd->SetUnitCategory("Length");
  GeThickDetectorThicknessCmd->SetParameterName("choice",false);
  GeThickDetectorThicknessCmd->AvailableForStates(G4State_PreInit);

  GeThickDetectorRadiusCmd = new G4UIcmdWithADoubleAndUnit("/e17011_sim/det/setGeThickDetectorRadius",this);
  GeThickDetectorRadiusCmd->SetGuidance("Set the Thick Ge Detector Radius.");
  GeThickDetectorRadiusCmd->SetUnitCategory("Length");
  GeThickDetectorRadiusCmd->SetParameterName("choice",false);
  GeThickDetectorRadiusCmd->AvailableForStates(G4State_PreInit);
  
  GeThickDetectorMatCmd = new G4UIcmdWithAString("/e17011_sim/det/setGeThickDetectorMate",this);
  GeThickDetectorMatCmd->SetGuidance("Select Material of the Thick Ge Detector.");
  GeThickDetectorMatCmd->SetParameterName("choice",false);
  GeThickDetectorMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  GeThinDetectorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/e17011_sim/det/setGeThinDetectorThickness",this);
  GeThinDetectorThicknessCmd->SetGuidance("Set the Thin Ge Detector Thickness.");
  GeThinDetectorThicknessCmd->SetUnitCategory("Length");
  GeThinDetectorThicknessCmd->SetParameterName("choice",false);
  GeThinDetectorThicknessCmd->AvailableForStates(G4State_PreInit);

  GeThinDetectorRadiusCmd = new G4UIcmdWithADoubleAndUnit("/e17011_sim/det/setGeThinDetectorRadius",this);
  GeThinDetectorRadiusCmd->SetGuidance("Set the Thin Ge Detector Radius.");
  GeThinDetectorRadiusCmd->SetUnitCategory("Length");
  GeThinDetectorRadiusCmd->SetParameterName("choice",false);
  GeThinDetectorRadiusCmd->AvailableForStates(G4State_PreInit);
  
  GeThinDetectorMatCmd = new G4UIcmdWithAString("/e17011_sim/det/setGeThinDetectorMate",this);
  GeThinDetectorMatCmd->SetGuidance("Select Material of the Thin Ge Detector.");
  GeThinDetectorMatCmd->SetParameterName("choice",false);
  GeThinDetectorMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  DetectorSpacingCmd = new G4UIcmdWithADoubleAndUnit("/e17011_sim/det/setDetectororSpacing",this);
  DetectorSpacingCmd->SetGuidance("Set the inter-Detector Spacing.");
  DetectorSpacingCmd->SetUnitCategory("Length");
  DetectorSpacingCmd->SetParameterName("choice",false);
  DetectorSpacingCmd->AvailableForStates(G4State_PreInit);

  UseGeThickDetectorCmd = new G4UIcmdWithABool("/e17011_sim/det/UseGeThickDetector",this);
  UseGeThickDetectorCmd->SetGuidance("activate the thick detector");
  UseGeThickDetectorCmd->SetParameterName("choice",false);
  UseGeThickDetectorCmd->AvailableForStates(G4State_PreInit);

  UseGeThinDetectorCmd = new G4UIcmdWithABool("/e17011_sim/det/UseGeThinDetector",this);
  UseGeThinDetectorCmd->SetGuidance("activate the thin detector");
  UseGeThinDetectorCmd->SetParameterName("choice",false);
  UseGeThinDetectorCmd->AvailableForStates(G4State_PreInit);

  UseGeThickDetectorCryoCmd = new G4UIcmdWithABool("/e17011_sim/det/UseGeThickDetectorCryo",this);
  UseGeThickDetectorCryoCmd->SetGuidance("activate the thick detector Cryo");
  UseGeThickDetectorCryoCmd->SetParameterName("choice",false);
  UseGeThickDetectorCryoCmd->AvailableForStates(G4State_PreInit);

  UseGeThickDetectorCryoEndCapCmd = new G4UIcmdWithABool("/e17011_sim/det/UseGeThickDetectorCryoEndCap",this);
  UseGeThickDetectorCryoEndCapCmd->SetGuidance("activate the thick detectorCryo EndCap");
  UseGeThickDetectorCryoEndCapCmd->SetParameterName("choice",false);
  UseGeThickDetectorCryoEndCapCmd->AvailableForStates(G4State_PreInit);

  UseSegaCmd = new G4UIcmdWithABool("/e17011_sim/det/UseSega",this);
  UseSegaCmd->SetGuidance("activate Sega");
  UseSegaCmd->SetParameterName("choice",false);
  UseSegaCmd->AvailableForStates(G4State_PreInit);

  UseCloverCmd = new G4UIcmdWithABool("/e17011_sim/det/UseClover",this);
  UseCloverCmd->SetGuidance("activate Clover");
  UseCloverCmd->SetParameterName("choice",false);
  UseCloverCmd->AvailableForStates(G4State_PreInit);

  // BC - use LaBr3 

  UseLaBr3Cmd = new G4UIcmdWithABool("/e17011_sim/det/UseLaBr3",this);
  UseLaBr3Cmd->SetGuidance("activate LaBr3");
  UseLaBr3Cmd->SetParameterName("choice",false);
  UseLaBr3Cmd->AvailableForStates(G4State_PreInit);
  
  // BC - use LaBr3 Frame

  UseLaBr3FrameCmd = new G4UIcmdWithABool("/e17011_sim/det/UseLaBr3Frame",this);
  UseLaBr3FrameCmd->SetGuidance("activate LaBr3 Frame");
  UseLaBr3FrameCmd->SetParameterName("choice",false);
  UseLaBr3FrameCmd->AvailableForStates(G4State_PreInit);
  
  //BC - use 3Hen
  
  Use3HenCmd = new G4UIcmdWithABool("/e17011_sim/det/Use3Hen",this);
  Use3HenCmd->SetGuidance("activate 3Hen");
  Use3HenCmd->SetParameterName("choice",false);
  Use3HenCmd->AvailableForStates(G4State_PreInit);
  
  //BC - use MTC
  UseMTCCmd = new G4UIcmdWithABool("/e17011_sim/det/UseMTC",this);
  UseMTCCmd->SetGuidance("activate MTC");
  UseMTCCmd->SetParameterName("choice",false);
  UseMTCCmd->AvailableForStates(G4State_PreInit);
  
  //BC - use LENDA
  UseLENDACmd = new G4UIcmdWithABool("/e17011_sim/det/UseLENDA",this);
  UseLENDACmd->SetGuidance("activate LENDA");
  UseLENDACmd->SetParameterName("choice",false);
  UseLENDACmd->AvailableForStates(G4State_PreInit);
  
  //BC - use Getest
  UseGetestCmd = new G4UIcmdWithABool("/e17011_sim/det/UseGetest",this);
  UseGetestCmd->SetGuidance("activate Getest");
  UseGetestCmd->SetParameterName("choice",false);
  UseGetestCmd->AvailableForStates(G4State_PreInit);
  
  //BC - use Xe test for EXO
  UseEXOtestCmd = new G4UIcmdWithABool("/e17011_sim/det/UseEXOtest",this);
  UseEXOtestCmd ->SetGuidance("activate EXOtest");
  UseEXOtestCmd ->SetParameterName("choice",false);
  UseEXOtestCmd ->AvailableForStates(G4State_PreInit);
  
  //BC - use EJ204 Scintillator
  UseEJ204ScintCmd = new G4UIcmdWithABool("/e17011_sim/det/UseEJ204Scint",this);
  UseEJ204ScintCmd ->SetGuidance("activate EJ204 Scintillator");
  UseEJ204ScintCmd ->SetParameterName("choice",false);
  UseEJ204ScintCmd ->AvailableForStates(G4State_PreInit);
  
  //BC - use PSPMT
  UsePSPMTCmd = new G4UIcmdWithABool("/e17011_sim/det/UsePSPMT",this);
  UsePSPMTCmd ->SetGuidance("activate PSPMT");
  UsePSPMTCmd ->SetParameterName("choice",false);
  UsePSPMTCmd ->AvailableForStates(G4State_PreInit);
  
  //BC - use SiDSSD
  UseSiDSSDCmd = new G4UIcmdWithABool("/e17011_sim/det/UseSiDSSD",this);
  UseSiDSSDCmd ->SetGuidance("activate SiDSSD");
  UseSiDSSDCmd ->SetParameterName("choice",false);
  UseSiDSSDCmd ->AvailableForStates(G4State_PreInit);

  //BC - use beam pipe
  UsePipeCmd = new G4UIcmdWithABool("/e17011_sim/det/UsePipe",this);
  UsePipeCmd ->SetGuidance("activate Beam Pipe");
  UsePipeCmd ->SetParameterName("choice",false);
  UsePipeCmd ->AvailableForStates(G4State_PreInit);

  //BC - use CeBr3 Scintillator
  UseCeBr3ScintCmd = new G4UIcmdWithABool("/e17011_sim/det/UseCeBr3Scint",this);
  UseCeBr3ScintCmd ->SetGuidance("activate CeBr3 Scintillator");
  UseCeBr3ScintCmd ->SetParameterName("choice",false);
  UseCeBr3ScintCmd ->AvailableForStates(G4State_PreInit);
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

e17011_simDetectorMessenger::~e17011_simDetectorMessenger()
{
  //G4cout<<"~gedssdnrDetectorMessenger"<<G4endl;
  delete GeThickDetectorMatCmd;
  delete GeThickDetectorRadiusCmd;
  delete GeThickDetectorThicknessCmd;
  delete GeThinDetectorMatCmd;
  delete GeThinDetectorRadiusCmd;
  delete GeThinDetectorThicknessCmd;
  delete DetectorSpacingCmd;
  delete UseGeThickDetectorCmd;
  delete UseGeThinDetectorCmd;
  delete UseGeThickDetectorCryoCmd;
  delete UseGeThickDetectorCryoEndCapCmd;
  delete UseSegaCmd;
  delete UseCloverCmd;
  delete detDir;
  delete e17011_simDir;
  delete UseLaBr3Cmd;
  delete UseLaBr3FrameCmd;
  delete Use3HenCmd;
  delete UseMTCCmd;
  delete UseLENDACmd;
  delete UseGetestCmd;
  delete UseEXOtestCmd;
  delete UseEJ204ScintCmd;
  delete UsePSPMTCmd;
  delete UseSiDSSDCmd;
  delete UsePipeCmd;
  delete UseCeBr3ScintCmd;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void e17011_simDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == GeThickDetectorRadiusCmd ) 
    { myDetector->setGeThickDetectorRadius(GeThickDetectorRadiusCmd->GetNewDoubleValue(newValue));}
  else if( command == GeThickDetectorMatCmd )
    { myDetector->setGeThickDetectorMaterial(newValue);} 
  else if (command == GeThickDetectorThicknessCmd ) 
    { myDetector->setGeThickDetectorThickness(GeThickDetectorThicknessCmd->GetNewDoubleValue(newValue));}
  else if ( command == GeThinDetectorRadiusCmd ) 
    { myDetector->setGeThinDetectorRadius(GeThinDetectorRadiusCmd->GetNewDoubleValue(newValue));}
  else if( command == GeThinDetectorMatCmd )
    { myDetector->setGeThinDetectorMaterial(newValue);} 
  else if (command == GeThinDetectorThicknessCmd ) 
    { myDetector->setGeThinDetectorThickness(GeThickDetectorThicknessCmd->GetNewDoubleValue(newValue));}
  else if (command == DetectorSpacingCmd ) 
    { myDetector->setDetectorSpacing(DetectorSpacingCmd->GetNewDoubleValue(newValue));}
  else if (command == UseGeThickDetectorCmd ) 
    { myDetector->setUseGeThickDetector(UseGeThickDetectorCmd->GetNewBoolValue(newValue));}
  else if (command == UseGeThickDetectorCryoCmd ) 
    { myDetector->setUseGeThickDetectorCryo(UseGeThickDetectorCryoCmd->GetNewBoolValue(newValue));}
  else if (command == UseGeThickDetectorCryoEndCapCmd ) 
    { myDetector->setUseGeThickDetectorCryoEndCap(UseGeThickDetectorCryoEndCapCmd->GetNewBoolValue(newValue));}
  else if (command == UseGeThinDetectorCmd ) 
    { myDetector->setUseGeThinDetector(UseGeThinDetectorCmd->GetNewBoolValue(newValue));}
  else if (command == UseSegaCmd ) 
    { myDetector->setUseSega(UseSegaCmd->GetNewBoolValue(newValue));}
  else if (command == UseCloverCmd ) 
    { myDetector->setUseClover(UseCloverCmd->GetNewBoolValue(newValue));}
  else if (command == UseLaBr3Cmd )
    { myDetector->setUseLaBr3(UseLaBr3Cmd->GetNewBoolValue(newValue));}
  else if (command == UseLaBr3FrameCmd )
    { myDetector->setUseLaBr3Frame(UseLaBr3FrameCmd->GetNewBoolValue(newValue));}
  else if (command == Use3HenCmd )
    { myDetector->setUse3Hen(Use3HenCmd->GetNewBoolValue(newValue));}
  else if (command == UseMTCCmd )
    { myDetector->setUseMTC(UseMTCCmd->GetNewBoolValue(newValue));}
  else if (command == UseLENDACmd )
    { myDetector->setUseLENDA(UseLENDACmd->GetNewBoolValue(newValue));}
  else if (command == UseGetestCmd )
    { myDetector->setUseGetest(UseGetestCmd->GetNewBoolValue(newValue));}
  else if (command == UseEXOtestCmd )
    { myDetector->setUseEXOtest(UseEXOtestCmd->GetNewBoolValue(newValue));}
  else if (command == UseEJ204ScintCmd )
    { myDetector->setUseEJ204Scint(UseEJ204ScintCmd->GetNewBoolValue(newValue));}
  else if (command == UsePSPMTCmd )
    { myDetector->setUsePSPMT(UsePSPMTCmd->GetNewBoolValue(newValue));}
  else if (command == UseSiDSSDCmd )
    { myDetector->setUseSiDSSD(UseSiDSSDCmd->GetNewBoolValue(newValue));}
  else if (command == UsePipeCmd )
    { myDetector->setUsePipe(UsePipeCmd->GetNewBoolValue(newValue));}
  else if (command == UseCeBr3ScintCmd )
    { myDetector->setUseCeBr3Scint(UseCeBr3ScintCmd->GetNewBoolValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
