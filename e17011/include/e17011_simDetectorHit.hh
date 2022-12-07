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

#ifndef e17011_simDetectorHit_h
#define e17011_simDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class G4AttDef;

class e17011_simDetectorHit : public G4VHit
{
  public:

      e17011_simDetectorHit();
      ~e17011_simDetectorHit();
      e17011_simDetectorHit(const e17011_simDetectorHit &right);
      const e17011_simDetectorHit& operator=(const e17011_simDetectorHit &right);
      G4int operator==(const e17011_simDetectorHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      const std::map<G4String,G4AttDef>* GetAttDefs() const;
      std::vector<G4AttValue>* CreateAttValues() const;
      void Print();

  private:
      G4double edep;
      G4ThreeVector pos;
      G4double gtime;
      G4int stepno;
      G4int parentno;
      G4int trackno;
      G4String particletype;
      G4String particlename;
      G4String parentname;
      G4String processname;
      G4double kineticenergy;
      G4String volname;
      G4int volcopyno;
      G4ThreeVector prepos;
      G4ThreeVector postpos;
      G4ThreeVector deltapos;
      G4double fulldis;
      G4double parentenergy;
      G4ThreeVector parentmodir;
      G4double stepleng;
  
      static std::map<G4String,G4AttDef> fAttDefs;

  public:
      inline void SetEdep(G4double de)
      { edep = de; }
      inline void AddEdep(G4double de)
      { edep += de; }
      inline G4double GetEdep()
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos()
      { return pos; }
      inline void SetGtime(G4double gt)
      {gtime = gt;}
      inline void SetKineticEnergy(G4double ke)
      {kineticenergy = ke;}
       inline G4double GetGtime()
      {return gtime;}
    inline void SetStepno(G4int st)
        {stepno = st;}
    inline G4int GetStepno()
        {return stepno;}
  inline void SetParentno(G4int pa)
  {parentno = pa;}
  inline G4int GetParentno()
  {return parentno;}
  inline void SetParticletype(G4String ty)
  {particletype = ty;}
  inline G4String GetParticletype()
  {return particletype;}
  inline void SetProcessname(G4String pr)
  {processname = pr;}
  inline G4String GetProcessname()
  {return processname;}
  inline void SetParticlename(G4String nam)
  {particlename = nam;}
  inline G4String GetParticlename()
  {return particlename;}
  inline void SetParentName(G4String nam)
  {parentname = nam;}
  inline G4String GetParentName()
  {return parentname;}
  inline void SetParentEnergy(G4double paren)
  {parentenergy = paren;}
  inline G4double GetParentEnergy()
  {return parentenergy;}
  inline G4double GetKineticEnergy()
  {return kineticenergy;}
  inline void SetParentMoDir(G4ThreeVector parmodir)
  {parentmodir = parmodir;}
  inline G4ThreeVector GetParentMoDir()
  {return parentmodir;}
   
  inline void SetVolName(G4String nam)
  {volname = nam;}
  inline G4String GetVolName()
  {return volname;}
  inline void SetVolCopyNo(G4int num)
  {volcopyno = num;}
  inline G4int GetVolCopyNo()
  {return volcopyno;}
  inline void SetTrackno(G4int tr)
  {trackno = tr;}
  inline G4int GetTrackno()
  {return trackno;}
  inline void SetPrePosition(G4ThreeVector pre)
  {prepos = pre;};
  inline G4ThreeVector GetPrePosition()
  {return prepos;};
  inline void SetPostPosition(G4ThreeVector post)
  {postpos = post;};
  inline G4ThreeVector GetPostPosition()
  {return postpos;};
  inline void SetDeltaPosition(G4ThreeVector del)
  {deltapos = del;};
  inline G4ThreeVector GetDeltaPosition()
  {return deltapos;};

  inline void SetStepLength(G4double slen)
  {stepleng = slen;};
  inline G4double GetStepLength()
  {return stepleng;};

  inline void SetDeltaLength(G4double dlen)
  {fulldis = dlen;};
  inline G4double GetDeltaLength()
  {return fulldis;};

};

typedef G4THitsCollection<e17011_simDetectorHit> e17011_simDetectorHitsCollection;

extern G4Allocator<e17011_simDetectorHit> e17011_simDetectorHitAllocator;

inline void* e17011_simDetectorHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) e17011_simDetectorHitAllocator.MallocSingle();
  return aHit;
}

inline void e17011_simDetectorHit::operator delete(void *aHit)
{
  e17011_simDetectorHitAllocator.FreeSingle((e17011_simDetectorHit*) aHit);
}

#endif


