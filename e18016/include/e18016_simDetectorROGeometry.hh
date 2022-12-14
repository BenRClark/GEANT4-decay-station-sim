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
#ifndef e18016_simDetectorROGeometry_h
#define e18016_simDetectorROGeometry_h 1

#include "G4VReadOutGeometry.hh"

class e18016_simDetectorROGeometry : public G4VReadOutGeometry
{
public:
  e18016_simDetectorROGeometry();
  e18016_simDetectorROGeometry(G4String);
  ~e18016_simDetectorROGeometry();

private:
  G4VPhysicalVolume* Build();

  //#include "e18016_simDetectorParameterDef.hh"

};

#endif
