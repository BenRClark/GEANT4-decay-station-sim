//This is the original stuff

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

// **********************************************************************

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

#include "e17011_simPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

e17011_simPrimaryGeneratorAction::e17011_simPrimaryGeneratorAction()
{
  // G4cout<<"e17011_simPrimaryGeneratorAction"<<G4endl;
  particleGun = new G4GeneralParticleSource ();
}

e17011_simPrimaryGeneratorAction::~e17011_simPrimaryGeneratorAction()
{
  // G4cout<<"~e17011_simPrimaryGeneratorAction"<<G4endl;
  delete particleGun;
}

void e17011_simPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // G4cout<<"GeneratePrimaries"<<G4endl;
  particleGun->GeneratePrimaryVertex(anEvent);
}




/**************************
 //this is the RD gamma stuff
/*************************

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

// **********************************************************************/

/*

#include "G4Event.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4PrimaryVertex.hh"
#include "TMath.h"
#include "TRandom.h"

#include "e17011_simPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

e17011_simPrimaryGeneratorAction::e17011_simPrimaryGeneratorAction()
{
  particleGun = new G4ParticleGun(); 
  counter =0;
}


e17011_simPrimaryGeneratorAction::~e17011_simPrimaryGeneratorAction()
{
  delete particleGun;
  G4cout<<counter<<" Gammas in the detector"<<G4endl;
}

//BC - use this for no 60Co decay, before 11/5/2015
/*
void e17011_simPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // G4cout<<"GeneratePrimaries"<<G4endl;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle= particleTable->FindParticle("e-");
  G4ParticleDefinition* particle2= particleTable->FindParticle("e+");
  G4ParticleDefinition* particle3= particleTable->FindParticle("gamma");
  
  position_Y = gRandom->Gaus(0.0,0.95)*cm;
  // position_X = 0.25*cm;
  position_X = gRandom->Gaus(0.0,1.13)*cm;
  // position_Y = 0.25*cm;
  // position_Z = gRandom->Uniform(-0.500,-0.490)*cm;  //0 to 100 microns deep
  // position_Z = gRandom->Uniform(-0.495,-0.485)*cm;  //50 to 150 microns deep
  // position_Z = gRandom->Uniform(-0.490,-0.480)*cm;  //100 to 200 microns deep
  // position_Z = gRandom->Uniform(-0.485,-0.475)*cm;  //150 to 250 microns deep
  // position_Z = gRandom->Uniform(-0.480,-0.470)*cm;  //200 to 300 microns deep
  // position_Z = gRandom->Uniform(-0.475,-0.465)*cm;  //250 to 350 microns deep
  // position_Z = gRandom->Uniform(-0.470,-0.460)*cm;  //300 to 400 microns deep
  // position_Z = gRandom->Uniform(-0.465,-0.455)*cm;  //350 to 450 microns deep
  // position_Z = gRandom->Uniform(-0.460,-0.450)*cm;  //400 to 500 microns deep
  // position_Z = gRandom->Uniform(-0.428,-0.416)*cm;  //720 to 840 microns deep  
  // position_Z = gRandom->Uniform(-0.410,-0.400)*cm;  //900 to 1000 microns deep
  // position_Z = gRandom->Uniform(-0.310,-0.300)*cm;     //1900 to 2000 microns deep
  position_Z = gRandom->Uniform(-2.505,-2.465)*cm;  //2450 to 2850 microns deep

  
  if(position_X >= -26 && position_X <= 26 && position_Y > -27 && position_Y < 25) {
    
    
    //emit gamma
    
    particleGun->SetParticleDefinition(particle3);
    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
    //particleGun->SetParticleEnergy(114.0*keV);
    
    Theta = gRandom->Uniform(0.0,3.141592653);
    Phi = gRandom->Uniform(0.0,6.2831853072);
    
    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
    
    particleGun->GeneratePrimaryVertex(anEvent);
    
    counter++;
  }
  
  }*/
 /*
 //BC - trying to include 60Co... should just use G4DecayManager I guess
void e17011_simPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //hardcoded for 60Co test.  Set to false for e14057 simulations
  Full_60Co = false;
  //Full_60Co = false;
  // G4cout<<"GeneratePrimaries"<<G4endl;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle= particleTable->FindParticle("e-");
  G4ParticleDefinition* particle2= particleTable->FindParticle("e+");
  G4ParticleDefinition* particle3= particleTable->FindParticle("gamma");
  
  //position_Y = gRandom->Gaus(0.0,0.95)*cm;
  // position_X = 0.25*cm;
  //position_X = gRandom->Gaus(0.0,1.13)*cm;
  // position_Y = 0.25*cm;
  // position_Z = gRandom->Uniform(-0.500,-0.490)*cm;  //0 to 100 microns deep
  // position_Z = gRandom->Uniform(-0.495,-0.485)*cm;  //50 to 150 microns deep
  // position_Z = gRandom->Uniform(-0.490,-0.480)*cm;  //100 to 200 microns deep
  // position_Z = gRandom->Uniform(-0.485,-0.475)*cm;  //150 to 250 microns deep
  // position_Z = gRandom->Uniform(-0.480,-0.470)*cm;  //200 to 300 microns deep
  // position_Z = gRandom->Uniform(-0.475,-0.465)*cm;  //250 to 350 microns deep
  // position_Z = gRandom->Uniform(-0.470,-0.460)*cm;  //300 to 400 microns deep
  // position_Z = gRandom->Uniform(-0.465,-0.455)*cm;  //350 to 450 microns deep
  // position_Z = gRandom->Uniform(-0.460,-0.450)*cm;  //400 to 500 microns deep
  // position_Z = gRandom->Uniform(-0.428,-0.416)*cm;  //720 to 840 microns deep  
  // position_Z = gRandom->Uniform(-0.410,-0.400)*cm;  //900 to 1000 microns deep
  // position_Z = gRandom->Uniform(-0.310,-0.300)*cm;     //1900 to 2000 microns deep
  //position_Z = gRandom->Uniform(-2.505,-2.465)*cm;  //2450 to 2850 microns deep

  //use for e14057 implant depth
  // position_Y = gRandom->Gaus(0.0,0.95)*cm;
  // position_X = gRandom->Gaus(0.0,1.13)*cm;
  // position_Z = gRandom->Uniform(-2.505,-2.465)*cm;  //2450 to 2850 microns deep

  //60Co hardcoded positions
  position_X = 0.15*cm;
  position_Y = 0.0*cm;
  position_Z = -2.9*cm;

  //uncomment for e14057 implant depth
  // if(position_X >= -26 && position_X <= 26 && position_Y > -27 && position_Y < 25) {
    
    
  

  if(Full_60Co){
    //numbers from NNDC 60Co decay (nuclear data chart)
    double random = gRandom->Uniform(0,100.0);

    //decay from 60Co to 4+ in 60Ni
    if(random <= 99.88){
      //do random for branches from 2506 keV level
      double random2 = gRandom->Uniform(0,1.0);
      //emit 1173: 99.85/99.857502 = 0.99992487
      //also, emit 1332
      if(random2 <= 0.99992487){
	//1173
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(1173.2*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//1332
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(1332.5*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      }
      //emit 347.1: 0.0075/99.857502 = 0.00007510702
      else if(random2 > 0.99992487 && random2 <= 0.99999997997){
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(347.1*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
	
	//do random for branches from 2159 keV level
	double random3 = gRandom->Uniform(0,1.0);
	//emit 826 and a 1332
	if(random3 <= 0.86363636){
	  //826
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(826.1*keV);
      
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
      
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	  particleGun->GeneratePrimaryVertex(anEvent);

	  //1332
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(1332.5*keV);
      
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
      
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	  particleGun->GeneratePrimaryVertex(anEvent);
	}else{
	  //emit 2159
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2158.6*keV);
      
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
      
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	  particleGun->GeneratePrimaryVertex(anEvent);
	}
	
      }else{
	//emit 2506
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2505.7*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
	
      }
    }
    //decay from 60Co to 2+ in 60Ni
    else{
      //1332
      particleGun->SetParticleDefinition(particle3);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(1332.5*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
    }

  }else{
    //emit gamma
    
    particleGun->SetParticleDefinition(particle3);
    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
    //particleGun->SetParticleEnergy(114.0*keV);
    
    Theta = gRandom->Uniform(0.0,3.141592653);
    Phi = gRandom->Uniform(0.0,6.2831853072);
    
    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
    
    particleGun->GeneratePrimaryVertex(anEvent);
    
    counter++;
    
  }
  
}
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

// **********************************************************************

/*
#include "G4Event.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4PrimaryVertex.hh"
#include "TMath.h"
#include "TRandom.h"

#include "e17011_simPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

e17011_simPrimaryGeneratorAction::e17011_simPrimaryGeneratorAction()
{
  // G4cout<<"e17011_simPrimaryGeneratorAction"<<G4endl;
  // particleGun = new G4GeneralParticleSource ();
  particleGun = new G4ParticleGun(); 
  
  
  //Set up the variables

  Full_118Ag = false;
  Full_68Co = false;
  Just_2511E0_ic = false;
  Just_2511E0_pp = false;
  Both_2511E0_pp_ic = false;
  Background = true;

  //for the silver
  Q_Value_Silver_3223 = 3924.4; // keV
  

  //For the silver simulation all that currently matters are the relative intensities
  //3.16 % abs. intensity for the 1938 keV transition feeding the 0+
  //10.13 % abs. intensity for the 2788 keV transition
  //9.07 % abs. intenstity for the 3224 keV transition

  Branch_Silver_1938 = 0.0516;
  Branch_Silver_2788 = 0.0813;
  Branch_Silver_3224 = 0.0907;

  //Total is 22.36 % of decays taken into account.  (Need to divide number of 118Ag ions by 0.2236 for accurate stats)
  
  Branch_Silver_EO = 1; //percentage of E0 transitions from the 1285

  Branch_Silver_pp = 0.85;
  Branch_Silver_ce = 0.15;


  BDecay_Energy_Dist_Silver_3223 = new TH1D("BDecay_Silver_3223","BDecay_Silver_3223",1000000,0,10000);

  for(int ii = 0; ii<1000000; ii++) {
    BDecay_Energy_Dist_Silver_3223->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_Silver_3223-(ii/100.0)),2.0)*((ii/100.0)+511.0));
  }

  

  //  position_X = 0.25*cm;
  //  position_Y = 0.25*cm;
  //  position_Z = -0.49*cm;

  //using the data from e11503 there were 306824 implants and an absolute 478 keV gamma ray intensity of 0.0218
  //from the mueller paper the relative intensity of the 478 was 15.5 and the relative intensity of the 1515 was 10.2
  //therefore the absolute efficiency of the 1515 is 0.01433.  The beta decay feeding of the 4026 keV state in 68Ni is
  //6.8 %.  Therefore x+y = 306824 and 0.068x = 0.01433*306824 so x = 64658.7 and therefore the 1+ branch ratio equals 
  //64658.7/306824.0 = 0.21.  The branch to the 7- is then 79 percent.  Not very good. 
 
  Branch_7minus = 11.0; // %  Based on comissioning run data.  THe raw SeGA spectrum was gated on the 2033 and the intensity of the 815 was looked at. 
  
  //7 minus levels
  Q_Value_3934p5 = 7615.5;  //keV
  Q_Value_3557p1 = 7992.9;  //keV
  Q_Value_3443p6 = 8106.4;  //keV
  Q_Value_3120p1 = 8429.9;  //keV
  Q_Value_2848p2 = 8701.8;  //keV

  Branching_Ratio_3934p5 = 1.50;  // %
  Branching_Ratio_3557p1 = 28.0;  // %
  Branching_Ratio_3443p6 = 62.0;  // % 
  Branching_Ratio_3120p1 = 5.00;  // %
  Branching_Ratio_2848p2 = 3.50;  // %

  BDecay_Energy_Dist_3934p5 = new TH1D("BDecay_3934p5","BDecay_3934p5",1000000,0,10000);
  BDecay_Energy_Dist_3557p1 = new TH1D("BDecay_3557p1","BDecay_3557p1",1000000,0,10000);
  BDecay_Energy_Dist_3443p6 = new TH1D("BDecay_3443p6","BDecay_3443p6",1000000,0,10000);
  BDecay_Energy_Dist_3120p1 = new TH1D("BDecay_3120p1","BDecay_3120p1",1000000,0,10000);
  BDecay_Energy_Dist_2848p2 = new TH1D("BDecay_2848p2","BDecay_2848p2",1000000,0,10000);

  for(int ii = 0; ii<1000000; ii++) {
    BDecay_Energy_Dist_3934p5->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_3934p5-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_3557p1->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_3557p1-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_3443p6->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_3443p6-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_3120p1->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_3120p1-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_2848p2->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_2848p2-(ii/100.0)),2.0)*((ii/100.0)+511.0));
  }


  //1 plus levels

  Energy_Offset_1plus = 100.0; //keV
  
  Q_Value_5774 = Energy_Offset_1plus + 5776.0;  //keV
  Q_Value_5549 = Energy_Offset_1plus + 6001.0;  //keV
  Q_Value_5512 = Energy_Offset_1plus + 6038.0;  //keV
  Q_Value_4167 = Energy_Offset_1plus + 7383.0;  //keV
  Q_Value_4026 = Energy_Offset_1plus + 7524.0;  //keV
  Q_Value_3987 = Energy_Offset_1plus + 7563.0;  //keV
  Q_Value_3542 = Energy_Offset_1plus + 8008.0;  //keV
  Q_Value_3148 = Energy_Offset_1plus + 8402.0;  //keV
  Q_Value_2745 = Energy_Offset_1plus + 8805.0;  //keV
  Q_Value_2511 = Energy_Offset_1plus + 9039.0;  //keV
  Q_Value_2033 = Energy_Offset_1plus + 9517.0;  //keV

  Branching_Ratio_5774 = 7.4;  // %
  Branching_Ratio_5549 = 4.2;  // %
  Branching_Ratio_5512 = 5.5;  // %
  Branching_Ratio_4167 = 7.5;  // %
  Branching_Ratio_4026 = 7.2;  // %
  Branching_Ratio_3987 = 4.5;  // %
  Branching_Ratio_3542 = 5.3;  // %
  Branching_Ratio_3148 = 6.8;  // %
  Branching_Ratio_2745 = 7.5;  // %
  Branching_Ratio_2511 = 3.8;  // %
  Branching_Ratio_2033 = 40.3;  // %

  BDecay_Energy_Dist_5774 = new TH1D("BDecay_5774","BDecay_5774",1000000,0,10000);
  BDecay_Energy_Dist_5549 = new TH1D("BDecay_5549","BDecay_5549",1000000,0,10000);
  BDecay_Energy_Dist_5512 = new TH1D("BDecay_5512","BDecay_5512",1000000,0,10000);
  BDecay_Energy_Dist_4167 = new TH1D("BDecay_4167","BDecay_4167",1000000,0,10000);
  BDecay_Energy_Dist_4026 = new TH1D("BDecay_4026","BDecay_4026",1000000,0,10000);
  BDecay_Energy_Dist_3987 = new TH1D("BDecay_3987","BDecay_3987",1000000,0,10000);
  BDecay_Energy_Dist_3542 = new TH1D("BDecay_3542","BDecay_3542",1000000,0,10000);
  BDecay_Energy_Dist_3148 = new TH1D("BDecay_3148","BDecay_3148",1000000,0,10000);
  BDecay_Energy_Dist_2745 = new TH1D("BDecay_2745","BDecay_2745",1000000,0,10000);
  BDecay_Energy_Dist_2511 = new TH1D("BDecay_2511","BDecay_2511",1000000,0,10000);
  BDecay_Energy_Dist_2033 = new TH1D("BDecay_2033","BDecay_2033",1000000,0,10000);

  for(int ii = 0; ii<1000000; ii++) {
    BDecay_Energy_Dist_5774->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_5774-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_5549->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_5549-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_5512->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_5512-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_4167->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_4167-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_4026->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_4026-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_3987->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_3987-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_3542->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_3542-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_3148->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_3148-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_2745->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_2745-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_2511->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_2511-(ii/100.0)),2.0)*((ii/100.0)+511.0));
    BDecay_Energy_Dist_2033->Fill(ii, TMath::Power((TMath::Power((ii/100.0),2.0)+2*(ii/100.0)*511),0.5)*TMath::Power((Q_Value_2033-(ii/100.0)),2.0)*((ii/100.0)+511.0));
  }

  Pair_Prod_Branch = 92.0;  //% from the NNDC for the 2511 keV transition
  E0_counter=0;
  E0_counter2=0;

  
}

e17011_simPrimaryGeneratorAction::~e17011_simPrimaryGeneratorAction()
{
  // G4cout<<"~e17011_simPrimaryGeneratorAction"<<G4endl;
  delete particleGun;

  delete BDecay_Energy_Dist_3934p5; 
  delete BDecay_Energy_Dist_3557p1; 
  delete BDecay_Energy_Dist_3443p6; 
  delete BDecay_Energy_Dist_3120p1; 
  delete BDecay_Energy_Dist_2848p2; 
  
  delete BDecay_Energy_Dist_5774;
  delete BDecay_Energy_Dist_5549; 
  delete BDecay_Energy_Dist_5512; 
  delete BDecay_Energy_Dist_4167; 
  delete BDecay_Energy_Dist_4026; 
  delete BDecay_Energy_Dist_3987; 
  delete BDecay_Energy_Dist_3542; 
  delete BDecay_Energy_Dist_3148; 
  delete BDecay_Energy_Dist_2745; 
  delete BDecay_Energy_Dist_2511;
  delete BDecay_Energy_Dist_2033;

  G4cout<<"Number of E0 Transitions following 1515 keV gamma ray: "<<E0_counter<<G4endl;
  G4cout<<"Number of E0 Transitions following beta decay to the 2511 keV state: "<<E0_counter2<<G4endl;
}

void e17011_simPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // G4cout<<"GeneratePrimaries"<<G4endl;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle= particleTable->FindParticle("e-");
  G4ParticleDefinition* particle2= particleTable->FindParticle("e+");
  G4ParticleDefinition* particle3= particleTable->FindParticle("gamma");

  position_X = gRandom->Gaus(0.0,0.75)*cm;
  // position_X = 0.25*cm;
  position_Y = gRandom->Gaus(0.0,0.75)*cm;
  // position_Y = 0.25*cm;
  // position_Z = gRandom->Uniform(-0.500,-0.490)*cm;  //0 to 100 microns deep
  // position_Z = gRandom->Uniform(-0.495,-0.485)*cm;  //50 to 150 microns deep
  // position_Z = gRandom->Uniform(-0.490,-0.480)*cm;  //100 to 200 microns deep
  // position_Z = gRandom->Uniform(-0.485,-0.475)*cm;  //150 to 250 microns deep
  // position_Z = gRandom->Uniform(-0.480,-0.470)*cm;  //200 to 300 microns deep
  // position_Z = gRandom->Uniform(-0.475,-0.465)*cm;  //250 to 350 microns deep
  // position_Z = gRandom->Uniform(-0.470,-0.460)*cm;  //300 to 400 microns deep
  // position_Z = gRandom->Uniform(-0.465,-0.455)*cm;  //350 to 450 microns deep
  // position_Z = gRandom->Uniform(-0.460,-0.450)*cm;  //400 to 500 microns deep
   position_Z = gRandom->Uniform(-0.410,-0.400)*cm;  //900 to 1000 microns deep
    
  

  if(Full_118Ag) {

    double random = gRandom->Uniform(0,0.2236);

    if(random<=Branch_Silver_2788) {
      //emit beta and 2788 gamma
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_Silver_3223->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
      //emit 2788 keV gamma
      particleGun->SetParticleDefinition(particle3);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(2788.0*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
    }

    if(random> Branch_Silver_2788 && random<= (Branch_Silver_2788 +Branch_Silver_3224)) {
      //emit beta and 3224 gamma
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_Silver_3223->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
      //emit 3224 keV gamma
      particleGun->SetParticleDefinition(particle3);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(3224.0*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
    }

    if(random>(Branch_Silver_2788 +Branch_Silver_3224)) {
      //emit beta and 1938 gamma
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_Silver_3223->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
      //emit 1938 keV gamma
      particleGun->SetParticleDefinition(particle3);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(1938.0*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);

      double random2 = gRandom->Uniform(0,1.0);

      if(random2<=Branch_Silver_EO) {
	//emit E0 
	double random3 = gRandom->Uniform(0,1.0);

	if(random3<= Branch_Silver_ce) {
	  //emit conversion electron
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(230.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	  
	  // G4cout<<" Conv Elec Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
	  
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	}
	else {
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(250.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	  
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);
	  

	}
	
      }
      
      
      else{
	//emit 798 and 487 
	//emit 798 keV gamma
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(798.0*keV);
	
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
	
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	//	particleGun->GeneratePrimaryVertex(anEvent);
	
	//emit 487 keV gamma
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(487.0*keV);
	
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
	
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	//	particleGun->GeneratePrimaryVertex(anEvent);
	
      }
    }
  }




  
  if(Full_68Co) {
  
    double random = gRandom->Uniform(0,100);
  
    if(random <= Branch_7minus) {
      //Beta decay from the 7 minus state

      double random2 = gRandom->Uniform(0,100);
      double value = 0;
    
      if(random2 < Branching_Ratio_3934p5) {
	//Decay to the 3934.5 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3934p5->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 377 keV gamma
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(377.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//follow cascade
	double random3 = gRandom->Uniform(0,100);
          
	if(random3 <= 69.9) {
	  //emit 709 keV gamma followed by 815 followed by 2033
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(709.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(815.0*keV);
	
	  //	  Theta = gRandom->Uniform(0.0,3.141592653);
	  //	  Phi = gRandom->Uniform(0.0,6.2831853072);
	  
	  //temporarily set opposite to 709 to mitigate fake summing in SeGA
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(-(TMath::Sin(Theta)*TMath::Cos(Phi)),-(TMath::Sin(Theta)*TMath::Sin(Phi)),-(TMath::Cos(Theta))));

	  //  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  // particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);	
	
	}
	else {
	  //emit 114 keV gamma
	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(114.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  double random4 = gRandom->Uniform(0,100);
	
	  if(random4 <= 54.4) {
	    //emit 324 keV then 272 then 815 then 2033
	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(324.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);


	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(272.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(815.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    // particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2033.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	  }

	  else{
	    //emit 595 kev then 815 then 2033
	  
	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(595.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(815.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    // particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2033.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);
	  }
	
	}
      
      }
      value = value + Branching_Ratio_3934p5;
    


      if(random2 >= value && random2 < value + Branching_Ratio_3557p1){
	//Decay to the 3557.1 keV state
      
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3557p1->GetRandom()*keV);

	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
       
	particleGun->GeneratePrimaryVertex(anEvent);

      
	//gamma decay 
	double random3 = gRandom->Uniform(0,100);
          
	if(random3 <= 69.9) {
	  //emit 709 keV gamma followed by 815 followed by 2033
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(709.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(815.0*keV);
	
	  //	  Theta = gRandom->Uniform(0.0,3.141592653);
	  //  Phi = gRandom->Uniform(0.0,6.2831853072);

	  //temporarily set it to be opposite the 709 (mitigate artificial summing in SeGA)
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(-(TMath::Sin(Theta)*TMath::Cos(Phi)),-(TMath::Sin(Theta)*TMath::Sin(Phi)),-(TMath::Cos(Theta))));

	  //  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  // particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);	
	
	}
	else {
	  //emit 114 keV gamma
	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(114.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  double random4 = gRandom->Uniform(0,100);
	
	  if(random4 <= 54.4) {
	    //emit 324 keV then 272 then 815 then 2033
	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(324.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);


	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(272.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(815.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    //  particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2033.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	  }

	  else{
	    //emit 595 kev then 815 then 2033
	  
	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(595.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(815.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    // particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2033.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);
	  }	
	}
      }
      value = value + Branching_Ratio_3557p1;
    
      if(random2 >= value && random2 < value + Branching_Ratio_3443p6){
	//Decay to the 3443.6 keV state

	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3443p6->GetRandom()*keV);

	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	double random4 = gRandom->Uniform(0,100);
	
	if(random4 <= 54.4) {
	  //emit 324 keV then 272 then 815 then 2033
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(324.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);


	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(272.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(815.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  //	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);

	}

	else{
	  //emit 595 kev then 815 then 2033
	  
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(595.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(815.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  //  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);
	}	
      }
    
      value = value + Branching_Ratio_3443p6;
    
    
    
      if(random2 >= value && random2 < value + Branching_Ratio_3120p1){
	//Decay to the 3120.1 keV state
      
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3120p1->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
	//emit 272 gamma then 815 then 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(272.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(815.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	//	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      }
    
      value = value + Branching_Ratio_3120p1;
    
    
      if(random2 >= value && random2 < value + Branching_Ratio_2848p2){
	//Decay to the 2848.2 keV state
      
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_2848p2->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(815.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	//	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
    }
  
    else {
      //Beta decay from the 1 plus state
    
      double random2 = gRandom->Uniform(0,100);
      double value = 0;
    
      if(random2 < Branching_Ratio_5774) {
	//Decay to the 5774 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_5774->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 3741 keV gamma and 2033 keV gamma
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(3741.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_5774;
    
      if(random2 >= value && random2 < value + Branching_Ratio_5549){
	//decay to the 5549 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_5549->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 3516 keV gamma and 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(3516.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_5549; 


      if(random2 >= value && random2 < value + Branching_Ratio_5512){
	//decay to the 5512 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_5512->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 3479 keV gamma and 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(3479.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_5512; 





      if(random2 >= value && random2 < value + Branching_Ratio_4167){
	//decay to 4167 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_4167->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 1422 keV gamma 

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(1422.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//branch 709 and 2745

	double random3 = gRandom->Uniform();
      
	if(random3 < 29.0) {
	  //emit 709 and 2033
	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(709.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	}
	else {
	  //emit 2745 keV
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2745.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	}
      
      }
      value = value +  Branching_Ratio_4167; 

      if(random2 >= value && random2 < value + Branching_Ratio_4026){
	//decay to 4026 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_4026->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	//	G4cout<<" bdecay Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);


	//assume all population yields 1515 keV gamma as from Mueller Paper
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(1515*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//Once the 2511 keV state is populated the decay can proceed by the 478 keV gamma or by the E0 decay to probably the ground state
	//In the Mueller paper the 511 keV peak (treated as a positron-electron anhilation) was 20 % relative intensity and the
	//478 keV gamma ray was 15.5 relative.  The only means of depopulating the 2511 are by 478 keV gamma or E0 which should be 99% pair production
	//and 1% conversion electron.

	double random3 = gRandom->Uniform(0.0,100.0);
      
	//	if(random3 <= 43.66) {
	if(random3 <= 95.0) {
	  //emit gamma rays 478 and 2033
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(478.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	}

	else {
	  E0_counter = E0_counter+1;
	  //	G4cout<<"Number of E0 Transitions: "<<E0_counter<<G4endl;
	
	  //E0 transition
	  double random4 = gRandom->Uniform(0.0, 100.0);
	
	  if(random4 <=  Pair_Prod_Branch) {
	    //pair produce
	    //electron
	    particleGun->SetParticleDefinition(particle);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(744.5*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);
	  
	    // G4cout<<" Pair Prod Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
	  
	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    //positron
	    particleGun->SetParticleDefinition(particle2);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(744.5*keV);
         
	    particleGun->SetParticleMomentumDirection(G4ThreeVector(-(TMath::Sin(Theta)*TMath::Cos(Phi)),-(TMath::Sin(Theta)*TMath::Sin(Phi)),-(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);
	  

	  }
	  else {
	    //conversion electron
	    particleGun->SetParticleDefinition(particle);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2511*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);
	  
	    // G4cout<<" Conv Elec Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
	  
	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);
	  }
	}
      
      }
      value = value + Branching_Ratio_4026;

      if(random2 >= value && random2 < value + Branching_Ratio_3987){
	//decay to 3987 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3987->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//gamma cascade 1139, 815, 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(1139.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
          
	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(815.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	//	particleGun->GeneratePrimaryVertex(anEvent);
        
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_3987; 

      if(random2 >= value && random2 < value + Branching_Ratio_3542){
	//decay to the 3542 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3542->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 694, 815, and 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(694.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
          
	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(815.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	//	particleGun->GeneratePrimaryVertex(anEvent);
        
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value +  Branching_Ratio_3542; 

      if(random2 >= value && random2 < value + Branching_Ratio_3148){
	//decay to 3148 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3148->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 1115 and 2033 keV
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(1115.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
        
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      
      }
      value = value + Branching_Ratio_3148; 

      if(random2 >= value && random2 < value + Branching_Ratio_2745){
	//decay to 2745 keV state
      
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_2745->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);
      
	double random3 = gRandom->Uniform();
      
	if(random3 < 29.0) {
	  //emit 709 and 2033
	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(709.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	}
	else {
	  //emit 2745 keV
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2745.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	}      
      
      }
      value = value + Branching_Ratio_2745; 


      if(random2 >= value && random2 < value + Branching_Ratio_2511){
	//decay to 2511 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_2511->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//Once the 2511 keV state is populated the decay can proceed by the 478 keV gamma or by the E0 decay to probably the ground state
	//In the Mueller paper the 511 keV peak (treated as a positron-electron anhilation) was 20 % relative intensity and the
	//478 keV gamma ray was 15.5 relative.  The only means of depopulating the 2511 are by 478 keV gamma or E0 which should be 99% pair production
	//and 1% conversion electron.

	double random3 = gRandom->Uniform(0.0,100.0);
      
	//	if(random3 <= 43.66) {
	if(random3 <= 95.00) {
	  //emit gamma rays
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(478*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	}

	else {
	  E0_counter2 = E0_counter2+1;
	  //	G4cout<<"Number of E0 Transitions: "<<E0_counter<<G4endl;
	  //E0 transition
	  double random4 = gRandom->Uniform(0.0, 100.0);
	
	  if(random4 <=  Pair_Prod_Branch) {
	    //pair produce
	    //electron
	    particleGun->SetParticleDefinition(particle);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(744.5*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);
	  
	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    //positron
	    particleGun->SetParticleDefinition(particle2);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(744.5*keV);
         
	    particleGun->SetParticleMomentumDirection(G4ThreeVector(-(TMath::Sin(Theta)*TMath::Cos(Phi)),-(TMath::Sin(Theta)*TMath::Sin(Phi)),-(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);
	  

	  }
	  else {
	    //conversion electron
	    particleGun->SetParticleDefinition(particle);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2511*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);
	  
	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);
	  }
	}
      
      }
      value = value +  Branching_Ratio_2511; 

      if(random2 >= value && random2 < value + Branching_Ratio_2033){
	//decay to 2033 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_2033->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
	
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      
      }
    }
  }
  

  if(Just_2511E0_ic) {
    
    particleGun->SetParticleDefinition(particle);
    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
    particleGun->SetParticleEnergy(BDecay_Energy_Dist_4026->GetRandom()*keV);
    
    Theta = gRandom->Uniform(0.0,3.141592653);
    Phi = gRandom->Uniform(0.0,6.2831853072);
    
    //	G4cout<<" bdecay Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
    
    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
    
    particleGun->GeneratePrimaryVertex(anEvent);
    
    // particleGun->SetParticleDefinition(particle3);
    // particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
    // particleGun->SetParticleEnergy(1515*keV);
    
    // Theta = gRandom->Uniform(0.0,3.141592653);
    // Phi = gRandom->Uniform(0.0,6.2831853072);
    
    // particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
    
    // particleGun->GeneratePrimaryVertex(anEvent);
    
    double random3 = gRandom->Uniform(0,100.0);
    
    if(random3 > 92.0) {
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(2511*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      // G4cout<<" Conv Elec Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
    }
  }
    

  if(Just_2511E0_pp) {
    particleGun->SetParticleDefinition(particle);
    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
    particleGun->SetParticleEnergy(BDecay_Energy_Dist_4026->GetRandom()*keV);
    
    Theta = gRandom->Uniform(0.0,3.141592653);
    Phi = gRandom->Uniform(0.0,6.2831853072);
    
    //	G4cout<<" bdecay Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
    
    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
    
    particleGun->GeneratePrimaryVertex(anEvent);

    double random3 = gRandom->Uniform(0,100.0);
    
    if(random3 <= 92.0) {

      //electron
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(744.5*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      // G4cout<<" Pair Prod Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
      //positron
      particleGun->SetParticleDefinition(particle2);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(744.5*keV);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector(-(TMath::Sin(Theta)*TMath::Cos(Phi)),-(TMath::Sin(Theta)*TMath::Sin(Phi)),-(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
    }
  }

  if(Both_2511E0_pp_ic) {
    
    particleGun->SetParticleDefinition(particle);
    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
    particleGun->SetParticleEnergy(BDecay_Energy_Dist_4026->GetRandom()*keV);
    
    Theta = gRandom->Uniform(0.0,3.141592653);
    Phi = gRandom->Uniform(0.0,6.2831853072);
    
    //	G4cout<<" bdecay Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
    
    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
    
    particleGun->GeneratePrimaryVertex(anEvent);

    double random3 = gRandom->Uniform(0,100.0);
    
    if(random3 <= 92.0) {

      //electron
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(744.5*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      // G4cout<<" Pair Prod Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
      //positron
      particleGun->SetParticleDefinition(particle2);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(744.5*keV);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector(-(TMath::Sin(Theta)*TMath::Cos(Phi)),-(TMath::Sin(Theta)*TMath::Sin(Phi)),-(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
    }

    else {
      
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(2511*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      // G4cout<<" Conv Elec Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
    }

  }
  
  if(Background) {
    //this removes only the E0 transitions

    double random = gRandom->Uniform(0,100);
  
    if(random <= Branch_7minus) {
      //Beta decay from the 7 minus state

      double random2 = gRandom->Uniform(0,100);
      double value = 0;
    
      if(random2 < Branching_Ratio_3934p5) {
	//Decay to the 3934.5 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3934p5->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 377 keV gamma
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(377.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//follow cascade
	double random3 = gRandom->Uniform(0,100);
          
	if(random3 <= 69.9) {
	  //emit 709 keV gamma followed by 815 followed by 2033
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(709.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(815.0*keV);
	
	  //	  Theta = gRandom->Uniform(0.0,3.141592653);
	  //	  Phi = gRandom->Uniform(0.0,6.2831853072);
	  
	  //temporarily set opposite to 709 to mitigate fake summing in SeGA
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(-(TMath::Sin(Theta)*TMath::Cos(Phi)),-(TMath::Sin(Theta)*TMath::Sin(Phi)),-(TMath::Cos(Theta))));

	  //  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  // particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);	
	
	}
	else {
	  //emit 114 keV gamma
	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(114.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  double random4 = gRandom->Uniform(0,100);
	
	  if(random4 <= 54.4) {
	    //emit 324 keV then 272 then 815 then 2033
	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(324.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);


	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(272.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(815.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    // particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2033.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	  }

	  else{
	    //emit 595 kev then 815 then 2033
	  
	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(595.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(815.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    // particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2033.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);
	  }
	
	}
      
      }
      value = value + Branching_Ratio_3934p5;
    


      if(random2 >= value && random2 < value + Branching_Ratio_3557p1){
	//Decay to the 3557.1 keV state
      
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3557p1->GetRandom()*keV);

	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
       
	particleGun->GeneratePrimaryVertex(anEvent);

      
	//gamma decay 
	double random3 = gRandom->Uniform(0,100);
          
	if(random3 <= 69.9) {
	  //emit 709 keV gamma followed by 815 followed by 2033
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(709.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(815.0*keV);
	
	  //	  Theta = gRandom->Uniform(0.0,3.141592653);
	  //  Phi = gRandom->Uniform(0.0,6.2831853072);

	  //temporarily set it to be opposite the 709 (mitigate artificial summing in SeGA)
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(-(TMath::Sin(Theta)*TMath::Cos(Phi)),-(TMath::Sin(Theta)*TMath::Sin(Phi)),-(TMath::Cos(Theta))));

	  //  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  // particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);	
	
	}
	else {
	  //emit 114 keV gamma
	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(114.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  double random4 = gRandom->Uniform(0,100);
	
	  if(random4 <= 54.4) {
	    //emit 324 keV then 272 then 815 then 2033
	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(324.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);


	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(272.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(815.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    //  particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2033.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	  }

	  else{
	    //emit 595 kev then 815 then 2033
	  
	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(595.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(815.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    // particleGun->GeneratePrimaryVertex(anEvent);

	    particleGun->SetParticleDefinition(particle3);
	    particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	    particleGun->SetParticleEnergy(2033.0*keV);
	  
	    Theta = gRandom->Uniform(0.0,3.141592653);
	    Phi = gRandom->Uniform(0.0,6.2831853072);

	    particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	    particleGun->GeneratePrimaryVertex(anEvent);
	  }	
	}
      }
      value = value + Branching_Ratio_3557p1;
    
      if(random2 >= value && random2 < value + Branching_Ratio_3443p6){
	//Decay to the 3443.6 keV state

	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3443p6->GetRandom()*keV);

	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	double random4 = gRandom->Uniform(0,100);
	
	if(random4 <= 54.4) {
	  //emit 324 keV then 272 then 815 then 2033
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(324.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);


	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(272.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(815.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  //	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);

	}

	else{
	  //emit 595 kev then 815 then 2033
	  
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(595.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(815.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  //  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);
	}	
      }
    
      value = value + Branching_Ratio_3443p6;
    
    
    
      if(random2 >= value && random2 < value + Branching_Ratio_3120p1){
	//Decay to the 3120.1 keV state
      
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3120p1->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
	//emit 272 gamma then 815 then 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(272.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(815.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	//	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      }
    
      value = value + Branching_Ratio_3120p1;
    
    
      if(random2 >= value && random2 < value + Branching_Ratio_2848p2){
	//Decay to the 2848.2 keV state
      
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_2848p2->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(815.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	//	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
    }
  
    else {
      //Beta decay from the 1 plus state
    
      double random2 = gRandom->Uniform(0,100);
      double value = 0;
    
      if(random2 < Branching_Ratio_5774) {
	//Decay to the 5774 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_5774->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 3741 keV gamma and 2033 keV gamma
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(3741.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_5774;
    
      if(random2 >= value && random2 < value + Branching_Ratio_5549){
	//decay to the 5549 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_5549->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 3516 keV gamma and 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(3516.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_5549; 


      if(random2 >= value && random2 < value + Branching_Ratio_5512){
	//decay to the 5512 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_5512->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 3479 keV gamma and 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(3479.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_5512; 





      if(random2 >= value && random2 < value + Branching_Ratio_4167){
	//decay to 4167 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_4167->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 1422 keV gamma 

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(1422.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//branch 709 and 2745

	double random3 = gRandom->Uniform();
      
	if(random3 < 29.0) {
	  //emit 709 and 2033
	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(709.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	}
	else {
	  //emit 2745 keV
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2745.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	}
      
      }
      value = value +  Branching_Ratio_4167; 

      if(random2 >= value && random2 < value + Branching_Ratio_4026){
	double random3 = gRandom->Uniform(0.0,100.0);
      
	//if(random3 <= 43.66) {
	if(random3 <= 95.00) {
	  //decay to 4026 keV state
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(BDecay_Energy_Dist_4026->GetRandom()*keV);
      
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);

	  //	G4cout<<" bdecay Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;

	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	  particleGun->GeneratePrimaryVertex(anEvent);


	  //assume all population yields 1515 keV gamma as from Mueller Paper
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(1515*keV);
      
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
      
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	  particleGun->GeneratePrimaryVertex(anEvent);

	  //emit gamma rays 478 and 2033
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(478.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	}      
      }
      value = value + Branching_Ratio_4026;

      if(random2 >= value && random2 < value + Branching_Ratio_3987){
	//decay to 3987 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3987->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//gamma cascade 1139, 815, 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(1139.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
          
	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(815.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	//	particleGun->GeneratePrimaryVertex(anEvent);
        
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_3987; 

      if(random2 >= value && random2 < value + Branching_Ratio_3542){
	//decay to the 3542 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3542->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 694, 815, and 2033
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(694.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
          
	particleGun->GeneratePrimaryVertex(anEvent);
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(815.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	//	particleGun->GeneratePrimaryVertex(anEvent);
        
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value +  Branching_Ratio_3542; 

      if(random2 >= value && random2 < value + Branching_Ratio_3148){
	//decay to 3148 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_3148->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	//emit 1115 and 2033 keV
      
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(1115.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
        
	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
      
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      
      }
      value = value + Branching_Ratio_3148; 

      if(random2 >= value && random2 < value + Branching_Ratio_2745){
	//decay to 2745 keV state
      
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_2745->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

	particleGun->GeneratePrimaryVertex(anEvent);
      
	double random3 = gRandom->Uniform();
      
	if(random3 < 29.0) {
	  //emit 709 and 2033
	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(709.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	}
	else {
	  //emit 2745 keV
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2745.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	}      
      
      }
      value = value + Branching_Ratio_2745; 


      if(random2 >= value && random2 < value + Branching_Ratio_2511){
	//decay to 2511 keV state

	double random3 = gRandom->Uniform(0.0,100.0);
	
	//	if(random3 <= 43.66) {
	if(random3 <= 95.00) {
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(BDecay_Energy_Dist_2511->GetRandom()*keV);
	  
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	  
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	  
	  particleGun->GeneratePrimaryVertex(anEvent);
	  
	  //Once the 2511 keV state is populated the decay can proceed by the 478 keV gamma or by the E0 decay to probably the ground state
	  //In the Mueller paper the 511 keV peak (treated as a positron-electron anhilation) was 20 % relative intensity and the
	  //478 keV gamma ray was 15.5 relative.  The only means of depopulating the 2511 are by 478 keV gamma or E0 which should be 99% pair production
	  //and 1% conversion electron.
	  
	  
	  //emit gamma rays
	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(478*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);

	  particleGun->SetParticleDefinition(particle3);
	  particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	  particleGun->SetParticleEnergy(2033.0*keV);
	
	  Theta = gRandom->Uniform(0.0,3.141592653);
	  Phi = gRandom->Uniform(0.0,6.2831853072);
	
	  particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
	
	  particleGun->GeneratePrimaryVertex(anEvent);
	
	}

      }
      value = value +  Branching_Ratio_2511; 

      if(random2 >= value && random2 < value + Branching_Ratio_2033){
	//decay to 2033 keV state
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(BDecay_Energy_Dist_2033->GetRandom()*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);

	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);

	particleGun->SetParticleDefinition(particle3);
	particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
	particleGun->SetParticleEnergy(2033.0*keV);
      
	Theta = gRandom->Uniform(0.0,3.141592653);
	Phi = gRandom->Uniform(0.0,6.2831853072);
	
	particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
	particleGun->GeneratePrimaryVertex(anEvent);
      
      
      }
    }
  }
}


    /*
    
      double random = gRandom->Uniform(0,100);
    
      if(random <= Branch_7minus) {
      //Beta decay from the 7 minus state
      
      double random2 = gRandom->Uniform(0,100);
      double value = 0;
    
      if(random2 < Branching_Ratio_3934p5) {
      //Decay to the 3934.5 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_3934p5->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

      particleGun->GeneratePrimaryVertex(anEvent);      
      }
      value = value + Branching_Ratio_3934p5;
    


      if(random2 >= value && random2 < value + Branching_Ratio_3557p1){
      //Decay to the 3557.1 keV state
      
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_3557p1->GetRandom()*keV);

      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
       
      particleGun->GeneratePrimaryVertex(anEvent);

      }
      value = value + Branching_Ratio_3557p1;
    
      if(random2 >= value && random2 < value + Branching_Ratio_3443p6){
      //Decay to the 3443.6 keV state

      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_3443p6->GetRandom()*keV);

      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

      particleGun->GeneratePrimaryVertex(anEvent);

		
      }
    
      value = value + Branching_Ratio_3443p6;
    
    
    
      if(random2 >= value && random2 < value + Branching_Ratio_3120p1){
      //Decay to the 3120.1 keV state
      
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_3120p1->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
     
      }
    
      value = value + Branching_Ratio_3120p1;
    
    
      if(random2 >= value && random2 < value + Branching_Ratio_2848p2){
      //Decay to the 2848.2 keV state
      
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_2848p2->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
          
      }
      }
  
      else {
      //Beta decay from the 1 plus state
    
      double random2 = gRandom->Uniform(0,100);
      double value = 0;
    
      if(random2 < Branching_Ratio_5774) {
      //Decay to the 5774 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_5774->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

      particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_5774;
    
      if(random2 >= value && random2 < value + Branching_Ratio_5549){
      //decay to the 5549 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_5549->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_5549; 


      if(random2 >= value && random2 < value + Branching_Ratio_5512){
      //decay to the 5512 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_5512->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_5512; 


      if(random2 >= value && random2 < value + Branching_Ratio_4167){
      //decay to 4167 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_4167->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value +  Branching_Ratio_4167; 

      if(random2 >= value && random2 < value + Branching_Ratio_4026){
      //decay to 4026 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_4026->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      //	G4cout<<" bdecay Theta: "<<Theta<<" Phi: "<<Phi<<G4endl;

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_4026;

      if(random2 >= value && random2 < value + Branching_Ratio_3987){
      //decay to 3987 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_3987->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

      particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value + Branching_Ratio_3987; 

      if(random2 >= value && random2 < value + Branching_Ratio_3542){
      //decay to the 3542 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_3542->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);
      
      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

      particleGun->GeneratePrimaryVertex(anEvent);
      
      }
      value = value +  Branching_Ratio_3542; 

      if(random2 >= value && random2 < value + Branching_Ratio_3148){
      //decay to 3148 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_3148->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);    
      
      }
      value = value + Branching_Ratio_3148; 

      if(random2 >= value && random2 < value + Branching_Ratio_2745){
      //decay to 2745 keV state
      
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_2745->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

      particleGun->GeneratePrimaryVertex(anEvent);
      }

      if(random2 >= value && random2 < value + Branching_Ratio_2511){
      //decay to 2511 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_2511->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));

      particleGun->GeneratePrimaryVertex(anEvent);      
      }
      value = value +  Branching_Ratio_2511; 

      if(random2 >= value && random2 < value + Branching_Ratio_2033){
      //decay to 2033 keV state
      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(position_X,position_Y,position_Z));
      particleGun->SetParticleEnergy(BDecay_Energy_Dist_2033->GetRandom()*keV);
      
      Theta = gRandom->Uniform(0.0,3.141592653);
      Phi = gRandom->Uniform(0.0,6.2831853072);

      particleGun->SetParticleMomentumDirection(G4ThreeVector((TMath::Sin(Theta)*TMath::Cos(Phi)),(TMath::Sin(Theta)*TMath::Sin(Phi)),(TMath::Cos(Theta))));
      
      particleGun->GeneratePrimaryVertex(anEvent);     
      
      }
      }
      }

 
  
  }
   */




