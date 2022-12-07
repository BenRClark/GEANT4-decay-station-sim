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
#include "e17011_simDetectorROGeometry.hh"
#include "e17011_simDummySD.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

 

e17011_simDetectorROGeometry::e17011_simDetectorROGeometry()
  : G4VReadOutGeometry()
{
  // G4cout<<"e17011_simDetectorROGeometry"<<G4endl;
//   G4double fDetectorThickness = 1.3* cm;
//   G4double fDetectorRadius = 0.5 * cm;
//   //fDetectorLength = 5.0 * cm;      
//   //fTargetLength  = 1.0 * cm;         
// //--------- Sizes of the principal geometrical components (solids)  ---------
//   G4double fWorldLength = fDetectorThickness + 1*cm;
//   G4double fWorldRadius = fDetectorRadius + 1.*cm;

  //#include "e17011_simDetectorParameterDef.icc"
}


e17011_simDetectorROGeometry::e17011_simDetectorROGeometry(G4String aString)
  : G4VReadOutGeometry(aString)
{
  //#include "e17011_simDetectorParameterDef.icc"
 
}

e17011_simDetectorROGeometry::~e17011_simDetectorROGeometry()
{
  // G4cout<<"~e17011_simDetectorROGeometry"<<G4endl;
}

G4VPhysicalVolume* e17011_simDetectorROGeometry::Build()
{
  // G4cout<<"e17011_simDetectorROGeometry Build"<<G4endl;
  G4double fDetectorThickness = 1.3* cm;
  G4double fDetectorRadius = 0.5 * cm;
  //fDetectorLength = 5.0 * cm;      
  //fTargetLength  = 1.0 * cm;         
//--------- Sizes of the principal geometrical components (solids)  ---------
  G4double fWorldLength = fDetectorThickness + 1*cm;
  G4double fWorldRadius = fDetectorRadius + 10.*cm;

  // A dummy material is used to fill the volumes of the readout geometry.
  // ( It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)
  G4Material* dummyMat  = new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);

  //Builds the ReadOut World:
  G4Tubs *ROWorldTubs = new G4Tubs("ROWorldTubs", 0,fWorldRadius,fWorldLength/2.,0.,twopi);
  G4LogicalVolume *ROWorldLog = new G4LogicalVolume(ROWorldTubs, dummyMat,
						    "ROWorldLogical", 0, 0, 0);
  G4PVPlacement *ROWorldPhys = new G4PVPlacement(0,G4ThreeVector(),
						 "ROWorldPhysical",
						 ROWorldLog,
						 0,false,0);
  // DSSD volume:
  G4VSolid * e17011_simROsol
    = new G4Tubs("e17011_simRO",0,fDetectorRadius,fDetectorThickness/2.,0.,twopi);
  G4LogicalVolume * e17011_simROlog
    = new G4LogicalVolume(e17011_simROsol,dummyMat,"e17011_simROlogical",0,0,0);
  G4VPhysicalVolume * e17011_simROphys
    = new G4PVPlacement(0,G4ThreeVector(),"e17011_simROphysical",e17011_simROlog,
			ROWorldPhys,false,0);



  // DSSD readout pixels

//   //Flags the cells as sensitive .The pointer here serves
//   // as a flag only to check for sensitivity.
//   // (Could we make it by a simple cast of a non-NULL value ?)
  e17011_simDummySD * dummySensi = new e17011_simDummySD;
  e17011_simROlog->SetSensitiveDetector(dummySensi);
  
  return ROWorldPhys;
}
