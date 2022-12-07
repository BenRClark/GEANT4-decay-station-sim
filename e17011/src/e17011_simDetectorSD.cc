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

#include "e17011_simDetectorSD.hh"
#include "e17011_simDetectorHit.hh"
#include "e17011_simTrackInformation.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

e17011_simDetectorSD::e17011_simDetectorSD(G4String name)
  :G4VSensitiveDetector(name), HCID(-1)
{  
  // G4cout<<"e17011_simDetectorSD"<<G4endl;
  G4String HCname;
  //collectionName.insert(HCname="DetectorCollection");
  collectionName.insert(HCname=name);
  positionResolution = 5*cm;
}

e17011_simDetectorSD::~e17011_simDetectorSD(){;
  // G4cout<<"~e17011_simDetectorSD"<<G4endl;
}

void e17011_simDetectorSD::Initialize(G4HCofThisEvent*HCE)
{
  // G4cout<<"e17011_simDetectorSD Initialize"<<G4endl;
  //static int HCID = -1;
  scintCollection = new e17011_simDetectorHitsCollection
    (GetName()/*SensitiveDetectorName*/,collectionName[0]); 
  if(HCID<0)
    { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,scintCollection);
}

G4bool e17011_simDetectorSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  // G4cout<<"e17011_simDetectorSD ProcessHits"<<G4endl;
  G4double edep = aStep->GetTotalEnergyDeposit();
  // G4cout<<"SD ProcessHits    Edep: "<<edep<<G4endl;
  //G4VPhysicalVolume *physVol = ROhist->GetVolume();
  
  if(edep==0.) return true; //snl changed

  e17011_simDetectorHit* aHit;
  int nHit = scintCollection->entries();
  //  G4cout << " processing hit with enetries " << nHit << G4endl;
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  //  G4cout << " prestep point " << aStep->GetPreStepPoint()->GetPosition() << G4endl;
  
  
  //   for(int i=0;i<nHit;i++)
  //   {
  //    aHit = (*scintCollection)[i];
  //    G4ThreeVector pos = aHit->GetPos();
  //    G4double dist2 = sqr(pos.x()-hitpos.x())
  //                    +sqr(pos.y()-hitpos.y())+sqr(pos.z()-hitpos.z());
  //    if(dist2<=sqr(positionResolution))
  //    aHit->AddEdep(edep);
  //    return true;
  //   }
  
  // G4cout<<"particle type: "<< aStep->GetTrack()->GetDefinition()->GetParticleType()<<"   particle name: "<< aStep->GetTrack()->GetDefinition()->GetParticleName()<<G4endl;
  
  // G4cout<<" Before '' aHit = new e17011_simDetectorHit();'' "<<G4endl; 
  
  aHit = new e17011_simDetectorHit();
  
  // G4cout << " Looking at Parent ID " << aStep->GetTrack()->GetParentID() << G4endl;

  // if(aStep->GetTrack()->GetParentID()==1){
  //   // if(aStep->GetTrack()->GetParentID()!=0){
  //   // G4cout<<" In loop"<<G4endl;
    
  //   e17011_simTrackInformation *info = (e17011_simTrackInformation*)(aStep->GetTrack()->GetUserInformation());
  //   // G4cout<<" In loop have info"<<G4endl;
  //   G4String pnam = info->GetOriginalParticle()->GetParticleName();
  //   // G4cout<<" In loop pnam  "<<pnam<<G4endl;
  //   aHit->SetParentName(pnam);
  //   G4double paren = info->GetOriginalEnergy();
  //   // G4cout<<" In loop paren "<<paren<<G4endl;
  //   aHit->SetParentEnergy(paren);
  //   G4ThreeVector parmodir = info->GetOriginalMomentum();
  //   // G4cout<<" In loop parmodir "<<parmodir<<G4endl;
  //   aHit->SetParentMoDir(parmodir);
  //   // G4cout<<" Leaving loop "<<G4endl;
  // }
  // else{
  //   aHit->SetParentName("decayingnucleus");
  // }

  // G4cout<<"SD ProcessHits 2"<<G4endl;

  G4double gt = aStep->GetPreStepPoint()->GetGlobalTime();
  aHit->SetGtime( gt );
  G4int st = aStep->GetTrack()->GetCurrentStepNumber();
  aHit->SetStepno( st );
  G4int tr = aStep->GetTrack()->GetTrackID();
  aHit->SetTrackno( tr );
  G4double ke = aStep->GetTrack()->GetKineticEnergy();
  aHit->SetKineticEnergy(ke);
  G4int pa = aStep->GetTrack()->GetParentID();
  aHit->SetParentno( pa );
  G4String ty = aStep->GetTrack()->GetDefinition()->GetParticleType();
  aHit->SetParticletype( ty );
  G4String nam = aStep->GetTrack()->GetDefinition()->GetParticleName();
  aHit->SetParticlename( nam );

  G4String volname = aStep->GetTrack()->GetVolume()->GetName();
  aHit->SetVolName( volname );
  //G4int volcopyno = aStep->GetTrack()->GetVolume()->GetCopyNo();

  // G4cout<<"SD ProcessHits 3 "<<gt<<"  "<<st<<"  "<<tr<<"  "<<ke<<"  "<<pa<<"  "<<ty<<"  "<<nam<<"  "<<volname<<G4endl;
  if(volname == "SegaActiveArea"){
    // G4cout<<"SD ProcessHits SeGA"<<G4endl;
    //only the prestep point is certain to be in appropriate SD volume
    G4int volcopyno = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(5)->GetCopyNo();
    aHit->SetVolCopyNo( volcopyno );
  }
  else if(volname == "CloverCrystal"){
    //   G4cout<<"SD ProcessHits Clover"<<G4endl;
    G4int crystalcopyno = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(0)->GetCopyNo();
    G4int Cloverno = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetCopyNo();
    G4int volcopyno = Cloverno * 4 + crystalcopyno;
    aHit->SetVolCopyNo( volcopyno );

  }

  else if(volname == "CloverCrystal2"){
    //   G4cout<<"SD ProcessHits Clover"<<G4endl;
    G4int crystalcopyno = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(0)->GetCopyNo();
    G4int Cloverno = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetCopyNo();
    G4int volcopyno = Cloverno * 4 + crystalcopyno;
    aHit->SetVolCopyNo( volcopyno );

  }
  //BC - handle interaction with LaBr3crystal (logical volume)
  else if(volname == "LaBr3Crys"){
    //   G4cout<<"SD ProcessHits LaBr3"<<G4endl;
    //The GetVolume(#) refers to the volume depth.  This must be correct and indicates how 
    //many far nested your sensitive volume is.  Here, it is 1... 
    //from http://geant4.web.cern.ch/geant4/G4UsersDocuments/UsersGuides/ForApplicationDeveloper/html/Detector/geomTouch.html  ::   
    //by depth it is meant the number of levels up in the tree to be 
    //considered (the default and current one is 0)
    G4int volcopyno = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetCopyNo();
    aHit->SetVolCopyNo( volcopyno );
  
  }

  //BC - handle interaction with EJ204 scintillator (logical volume)
  else if(volname == "EJ204Scint"){
    //   G4cout<<"SD ProcessHits EJ204 Scint"<<G4endl;
    //The GetVolume(#) refers to the volume depth.  This must be correct and indicates how 
    //many far nested your sensitive volume is.  Here, it is 1... 
    //from http://geant4.web.cern.ch/geant4/G4UsersDocuments/UsersGuides/ForApplicationDeveloper/html/Detector/geomTouch.html  ::   
    //by depth it is meant the number of levels up in the tree to be 
    //considered (the default and current one is 0)
    G4int volcopyno = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetCopyNo();
    aHit->SetVolCopyNo( volcopyno );
  
  }

  G4ThreeVector pre = aStep->GetPreStepPoint()->GetPosition();
  aHit->SetPrePosition( pre );
  G4ThreeVector post =aStep->GetPostStepPoint()->GetPosition();
  aHit->SetPostPosition( post );
  G4ThreeVector del = aStep->GetDeltaPosition();
  aHit->SetDeltaPosition( del );

  G4double dlen = aStep->GetStepLength();
  aHit->SetDeltaLength( dlen );
    
  G4double slen = aStep->GetStepLength();
  aHit->SetStepLength( slen );
    
  //  G4cout<<"SD ProcessHits "<<pre<<" "<<post<<"  "<<del<<"  "<<dlen<<"  "<<slen<<G4endl;


  //     G4Track *fTrack = aStep->GetTrack();
  //    const G4VProcess *proc = aStep->GetTrack()->GetCreatorProcess();
  //    G4String pr = proc->GetProcessName();
  //    G4String pr = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
  //     aHit->SetProcessname( pr );

  aHit->SetEdep( edep );
  
  // G4cout<<"SD ProcessHits About to insert hit"<<G4endl;

  //aHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  scintCollection->insert( aHit );

  // G4cout<<"Made it to the end of Hit Processing"<<G4endl;

  return true;
}

void e17011_simDetectorSD::EndOfEvent(G4HCofThisEvent*)
{;}

void e17011_simDetectorSD::clear()
{
} 

void e17011_simDetectorSD::DrawAll()
{
} 

void e17011_simDetectorSD::PrintAll()
{
} 

