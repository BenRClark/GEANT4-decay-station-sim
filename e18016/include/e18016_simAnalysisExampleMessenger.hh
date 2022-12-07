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
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// AnalysisExampleMessenger header
// --------------------------------------------------------------

#ifndef e18016_simAnalysisExampleMessenger_h
#define e18016_simAnalysisExampleMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class e18016_simAnalysisExample;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;


class e18016_simAnalysisExampleMessenger: public G4UImessenger {

  public:
    e18016_simAnalysisExampleMessenger(e18016_simAnalysisExample*);
   ~e18016_simAnalysisExampleMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);

  private:
    e18016_simAnalysisExample*  analysisexample;   
  
    G4UIdirectory*        analysisDir;

    G4UIcmdWithAnInteger* trackparentCmd;
    G4UIcmdWithAnInteger* trackdaughterCmd;
    G4UIcmdWithAString*   filenameCmd;
    G4UIcmdWithAString*   dirnameCmd;
};

#endif

