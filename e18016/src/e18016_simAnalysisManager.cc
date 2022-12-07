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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "e18016_simAnalysisManager.hh"

#include "G4ios.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Track.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4TrackStatus.hh"
#include "G4Step.hh"
#include "e18016_simDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


e18016_simAnalysisManager *ge18016_simAnalysisManager = (e18016_simAnalysisManager *)0;
//e18016_simAnalysisManager *ge18016_simAnalysisManager = (e18016_simAnalysisManager *);
//G4cout<<"mad it into AnalysisManage"<<G4endl;-it can't compile this-needs a constructor etc

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

