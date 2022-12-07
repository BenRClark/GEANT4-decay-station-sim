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

#ifndef e17011_simDetectorMessenger_h
#define e17011_simDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class e17011_simDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class e17011_simDetectorMessenger: public G4UImessenger
{
  public:
    e17011_simDetectorMessenger(e17011_simDetectorConstruction*);
   ~e17011_simDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    e17011_simDetectorConstruction* myDetector;
    
    G4UIdirectory*             e17011_simDir;
    G4UIdirectory*             detDir;

    G4UIcmdWithAString*        GeThickDetectorMatCmd;
    G4UIcmdWithADoubleAndUnit* GeThickDetectorRadiusCmd;
    G4UIcmdWithADoubleAndUnit* GeThickDetectorThicknessCmd;

    G4UIcmdWithAString*        GeThinDetectorMatCmd;
    G4UIcmdWithADoubleAndUnit* GeThinDetectorRadiusCmd;
    G4UIcmdWithADoubleAndUnit* GeThinDetectorThicknessCmd;

    G4UIcmdWithADoubleAndUnit* DetectorSpacingCmd;
   
    G4UIcmdWithABool*          UseGeThickDetectorCmd;
    G4UIcmdWithABool*          UseGeThinDetectorCmd;
    G4UIcmdWithABool*          UseGeThickDetectorCryoCmd;
    G4UIcmdWithABool*          UseGeThickDetectorCryoEndCapCmd;
    G4UIcmdWithABool*          UseSegaCmd;
    G4UIcmdWithABool*          UseCloverCmd;
    G4UIcmdWithABool*          UseCloverBracketCmd;  // for Clover Bracket

    G4UIcmdWithABool*          UseLaBr3Cmd;  //BC - LaBr3
    G4UIcmdWithABool*          UseLaBr3FrameCmd;  //BC - LaBr3 Frame
    G4UIcmdWithABool*          Use3HenCmd;   //BC - 3Hen
    G4UIcmdWithABool*          UseMTCCmd;    //BC - Moving tape collector   
    G4UIcmdWithABool*          UseLENDACmd;  //BC - LENDA   
    G4UIcmdWithABool*          UseGetestCmd; //BC - Getest
    G4UIcmdWithABool*          UseEXOtestCmd;//BC - Xe test for EXO
    G4UIcmdWithABool*          UseEJ204ScintCmd;//BC - EJ204 Scint e14057
    G4UIcmdWithABool*          UsePSPMTCmd;//BC - PSPMT e14057 + e16032
    G4UIcmdWithABool*          UseSiDSSDCmd;//BC - SiDSSD e14057
    G4UIcmdWithABool*          UsePipeCmd;//BC - Beam Pipe e14057 + e16032
    G4UIcmdWithABool*          UseCeBr3ScintCmd;//BC - CeBr3 Scint e16032
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

