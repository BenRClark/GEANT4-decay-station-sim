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

#include "e17011_simDetectorHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4AttValue.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"

G4Allocator<e17011_simDetectorHit> e17011_simDetectorHitAllocator;

e17011_simDetectorHit::e17011_simDetectorHit()
{;
  // G4cout<<"e17011_simDetectorHit Constructor"<<G4endl;
}

e17011_simDetectorHit::~e17011_simDetectorHit()
{;
  // G4cout<<"~e17011_simDetectorHit Destructor"<<G4endl;
}

e17011_simDetectorHit::e17011_simDetectorHit(const e17011_simDetectorHit &right)
  : G4VHit()
{
  //G4cout<<"e17011_simDetectorHit()"<<G4endl;

  edep = right.edep;
  pos = right.pos;
  gtime = right.gtime;
  stepno = right.stepno;
  parentno = right.parentno;
  particletype = right.particletype;
  processname = right.processname;
  particlename = right.particlename;
  parentname = right.parentname;
  parentenergy = right.parentenergy;
  parentmodir = right.parentmodir;
  kineticenergy = right.kineticenergy;

  volname = right.volname;
  volcopyno = right.volcopyno;
  trackno = right.trackno;
  prepos = right.prepos;
  postpos = right.postpos;
  deltapos = right.deltapos;

  fulldis = right.fulldis;

  stepleng = right.stepleng;


}

const e17011_simDetectorHit& e17011_simDetectorHit::operator=(const e17011_simDetectorHit &right)
{
  // G4cout<<"e17011_simDetectorHit operator="<<G4endl;
  edep = right.edep;
  pos = right.pos;
  gtime = right.gtime;
  stepno = right.stepno;
  parentno = right.parentno;
  particletype = right.particletype;
  processname = right.processname;
  particlename = right.particlename;
  parentname = right.parentname;
  parentenergy = right.parentenergy;
  parentmodir = right.parentmodir;
  kineticenergy = right.kineticenergy;

  volname = right.volname;
  volcopyno = right.volcopyno;
  trackno = right.trackno;
  prepos = right.prepos;
  postpos = right.postpos;
  deltapos = right.deltapos;

  fulldis = right.fulldis;
  stepleng = right.stepleng;
  return *this;
}

G4int e17011_simDetectorHit::operator==(const e17011_simDetectorHit &right) const
{
  return (this==&right) ? 1 : 0;
}

std::map<G4String,G4AttDef> e17011_simDetectorHit::fAttDefs;

void e17011_simDetectorHit::Draw()
{
  // G4cout<<"e17011_simDetectorHit Draw"<<G4endl;
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

const std::map<G4String,G4AttDef>* e17011_simDetectorHit::GetAttDefs() const
{
  //  G4cout<<"e17011_simDetectorHit GeAttDefs"<<G4endl;
  // G4AttDefs have to have long life.  Use static member...
  if (fAttDefs.empty()) {
    fAttDefs["HitType"] =
      G4AttDef("HitType","Type of hit","Physics","","G4String");
  }
  return &fAttDefs;
}

std::vector<G4AttValue>* e17011_simDetectorHit::CreateAttValues() const
{
  //G4cout<<"e17011_simDetectorHit CreateAttValues"<<G4endl;
  // Create expendable G4AttsValues for picking...
  std::vector<G4AttValue>* attValues = new std::vector<G4AttValue>;
  attValues->push_back
    (G4AttValue("HitType","e17011_simDetectorHit",""));
  G4cout << "Checking...\n" << G4AttCheck(attValues, GetAttDefs());
  return attValues;
}

void e17011_simDetectorHit::Print()
{;}


