//This is the original
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


//
// **********************************************************************



#ifndef e17011_simPrimaryGeneratorAction_h
#define e17011_simPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4GeneralParticleSource;
class G4Event;

class e17011_simPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  e17011_simPrimaryGeneratorAction();
  ~e17011_simPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);

private:
  G4GeneralParticleSource* particleGun;
  
};

#endif






//This is the RD Gamma Stuff


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

//
// **********************************************************************
/*
//BC - use this one for no 60Co level scheme, prior to 11/5/2015
#ifndef e17011_simPrimaryGeneratorAction_h
#define e17011_simPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TH1.h"

//class G4GeneralParticleSource;
class G4ParticleGun;
class G4Event;

class e17011_simPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  e17011_simPrimaryGeneratorAction();
  ~e17011_simPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);
  
private:
  G4ParticleGun* particleGun;
  
  G4double position_X;
  G4double position_Y;
  G4double position_Z;
  
  G4double Theta;
  G4double Phi;

  G4int counter;

};

#endif
*/

//BC modification to include 60Co level scheme
/*
#ifndef e17011_simPrimaryGeneratorAction_h
#define e17011_simPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TH1.h"

//class G4GeneralParticleSource;
class G4ParticleGun;
class G4Event;

class e17011_simPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  e17011_simPrimaryGeneratorAction();
  ~e17011_simPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);
  
private:
  G4ParticleGun* particleGun;
  
  G4double position_X;
  G4double position_Y;
  G4double position_Z;
  
  G4double Theta;
  G4double Phi;

  G4int counter;

  //BC - add 60Co
  G4bool Full_60Co;

};

#endif
*/


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


//
// **********************************************************************

/*
#ifndef e17011_simPrimaryGeneratorAction_h
#define e17011_simPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TH1.h"

//class G4GeneralParticleSource;
class G4ParticleGun;
class G4Event;

class e17011_simPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  e17011_simPrimaryGeneratorAction();
  ~e17011_simPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);

private:
  //G4GeneralParticleSource* particleGun;
  G4ParticleGun* particleGun;

  double Q_Value_Silver_3223;
  
  double Branch_Silver_1938;
  double Branch_Silver_2788;
  double Branch_Silver_3224;

  double Branch_Silver_EO; //percentage of E0 transitions from the 1285

  double Branch_Silver_pp;
  double Branch_Silver_ce;

  TH1D *BDecay_Energy_Dist_Silver_3223;



  double Branch_7minus;
  double Energy_Offset_1plus;
  
  //Decay from the 7- state in Co68
  TH1D *BDecay_Energy_Dist_3934p5;
  TH1D *BDecay_Energy_Dist_3557p1;
  TH1D *BDecay_Energy_Dist_3443p6;
  TH1D *BDecay_Energy_Dist_3120p1;
  TH1D *BDecay_Energy_Dist_2848p2;
  
  double Q_Value_3934p5;
  double Q_Value_3557p1;
  double Q_Value_3443p6;
  double Q_Value_3120p1;
  double Q_Value_2848p2;

  double Branching_Ratio_3934p5;
  double Branching_Ratio_3557p1;
  double Branching_Ratio_3443p6;
  double Branching_Ratio_3120p1;
  double Branching_Ratio_2848p2;

  //Decay from the 1+ state in Co68

  TH1D *BDecay_Energy_Dist_5774;
  TH1D *BDecay_Energy_Dist_5549;
  TH1D *BDecay_Energy_Dist_5512; 
  TH1D *BDecay_Energy_Dist_4167;
  TH1D *BDecay_Energy_Dist_4026;
  TH1D *BDecay_Energy_Dist_3987;
  TH1D *BDecay_Energy_Dist_3542;
  TH1D *BDecay_Energy_Dist_3148;
  TH1D *BDecay_Energy_Dist_2745;
  TH1D *BDecay_Energy_Dist_2511;
  TH1D *BDecay_Energy_Dist_2033;
  
  double Q_Value_5774;
  double Q_Value_5549;
  double Q_Value_5512;
  double Q_Value_4167;
  double Q_Value_4026;
  double Q_Value_3987;
  double Q_Value_3542;
  double Q_Value_3148;
  double Q_Value_2745;
  double Q_Value_2511;
  double Q_Value_2033;

  double Branching_Ratio_5774;
  double Branching_Ratio_5549;
  double Branching_Ratio_5512;
  double Branching_Ratio_4167;
  double Branching_Ratio_4026;
  double Branching_Ratio_3987;
  double Branching_Ratio_3542;
  double Branching_Ratio_3148;
  double Branching_Ratio_2745;
  double Branching_Ratio_2511;
  double Branching_Ratio_2033;

  double Pair_Prod_Branch;
  double Intern_Conv_Branch;
  double Theta;
  double Phi;

  G4double position_X;
  G4double position_Y;
  G4double position_Z;

  double E0_counter;
  double E0_counter2;
  
  G4bool Full_118Ag;
  G4bool Full_68Co;
  G4bool Just_2511E0_ic;
  G4bool Just_2511E0_pp;
  G4bool Both_2511E0_pp_ic;
  G4bool Background;
  
  //BC - add 60Co
  G4bool Full_60Co;
};

#endif
*/


