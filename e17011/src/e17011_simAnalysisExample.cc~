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

#include "e17011_simDetectorConstruction.hh"
#include "e17011_simAnalysisExample.hh"
#include "e17011_simAnalysisExampleMessenger.hh"

#include "e17011_simDetectorHit.hh"
//#include "e17011_simscintHit.hh"
//#include "e17011_simPMTHit.hh"

#include "G4ios.hh"

#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Run.hh"
#include "G4Track.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4TrackStatus.hh"
#include "G4Step.hh"
#include "G4Types.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <sstream>
#include <iostream>
#include <time.h>
#include <cstring>
#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include <math.h>

static const G4double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;

#if defined (G4ANALYSIS_USE_ROOT)

#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TRandom.h"

#endif /* defined (G4ANALYSIS_USE_ROOT) */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

e17011_simAnalysisExample::e17011_simAnalysisExample(){

  //cout<<"e17011_simAnalysisExample"<<endl;
  analysisexample = new e17011_simAnalysisExampleMessenger(this);

  filecount = 0;  

  LastDoItTime = (time_t)0;
  
#if defined (G4ANALYSIS_USE_ROOT)
  
  if (gSystem) gSystem->ProcessEvents();
  
  e17011_simtree = (TTree *)0; // e17011_sim out tree

  hEnergyDeposit = (TH1D*)0;
  hEnergyDepositPixel = (TH1D*)0;
  hEnergyDepositElectron = (TH1D*)0;
  hEnergyDepositGamma = (TH1D*)0;
  hEnergyDepositGamma_nocenter = (TH1D*)0;
  hEnergyDepositGamma2D = (TH2D*)0;

  hBeam_energy = (TH1D *)0; // beam energy
  hPrimaryPos  = (TH2D *)0; // primary position  

  //hGammaSingle = (TH1D*)0;
  //hGammaSingleAdd = (TH1D*)0;
  //for(int i=0; i<68;i++){
  //    hCrystalEdep[i] = (TH1D*)0;
  //}

  //hScintEdep_xy = (TH2D*)0; 
  //hScintEdep_yz = (TH2D*)0; 


#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

e17011_simAnalysisExample::~e17011_simAnalysisExample()
{
  // cout<<"~e17011_simAnalysisExample"<<endl;
  //OnceAWhileDoIt();
  
  delete analysisexample;

#if defined (G4ANALYSIS_USE_ROOT)
  
  if (gSystem) gSystem->ProcessEvents();
  
  
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::Construct(const G4VPhysicalVolume */*theWorldWolume*/)
{
  // cout<<"e17011_simAnalysisExample Construct"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  if (gSystem) gSystem->ProcessEvents();
  
  
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::ConstructParticle()
{
  //  cout<<"e17011_simAnalysisExample ConstructParticle"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  
  if (gSystem) gSystem->ProcessEvents();
  
  
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::ConstructProcess()
{
  //  cout<<"e17011_simAnalysisExample ConstructProcess"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  
  if (gSystem) gSystem->ProcessEvents();
  
  
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::SetCuts()
{
  //  G4cout<<"e17011_simAnalysisExample SetCuts"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  
  if (gSystem) gSystem->ProcessEvents();
  
  
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::GeneratePrimaries(const G4Event *anEvent,G4double beam)
{
  // G4cout<<"e17011_simAnalysisExample GeneratePrimaries"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  if (gSystem) gSystem->ProcessEvents();
  
  beam_energy = 0;
  beam_energy = beam;
  if (hBeam_energy) hBeam_energy->Fill(beam_energy / keV, 1.); // Step Length

  G4ThreeVector startpos = anEvent->GetPrimaryVertex()->GetPosition();

  //G4G4cout << "prim pos " << startpos.x()/micrometer << " " <<startpos.y()/micrometer << " " << startpos.z()/micrometer << G4endl;

  hPrimaryPos->Fill(startpos.x()/micrometer,startpos.y()/micrometer);
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//void e17011_simAnalysisExample::ResetStartLocation(G4double &newx, G4double &newy)
//{
  // #if defined (G4ANALYSIS_USE_ROOT)
  //   if (gSystem) gSystem->ProcessEvents();
  
  //   hPrimaryPosin->GetRandom2(newx,newy);
  //   G4cout << " random start value " << newx << " " << newy << G4endl;

  //   //hPrimaryPos->Fill(startpos.x()/micrometer,startpos.y()/micrometer);
  // #endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  //   OnceAWhileDoIt();
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::DetectorInfo(e17011_simDetectorConstruction* myDet){

  cout << "go the detector and the Ge detectors are " << myDet->GetUseGeThickDetector() << " and "
       << myDet->GetUseGeThinDetector() << G4endl;


  UseGeThickDetector = Int_t(myDet->GetUseGeThickDetector());
  UseGeThinDetector = Int_t(myDet->GetUseGeThinDetector());

  GeThickDetectorThickness = myDet->GetGeThickDetectorThickness();
  GeThinDetectorThickness = myDet->GetGeThinDetectorThickness();

  GeThickDetectorRadius = myDet->GetGeThickDetectorRadius();
  GeThinDetectorRadius = myDet->GetGeThinDetectorRadius();

  GeThickx = myDet->GetGeThickDetectorPos().x();
  GeThicky = myDet->GetGeThickDetectorPos().y();
  GeThickz = myDet->GetGeThickDetectorPos().z();
  GeThinx = myDet->GetGeThinDetectorPos().x();
  GeThiny = myDet->GetGeThinDetectorPos().y();
  GeThinz = myDet->GetGeThinDetectorPos().z();

  GeSpacing = myDet->GetGeSpacing();

  e17011_siminfo = (TTree *)gROOT->FindObject("e17011_siminfo");
  if(e17011_siminfo);
  else{
    e17011_siminfo = new TTree("e17011_siminfo","e17011_siminfo variables");
    e17011_siminfo->Branch("UseGeThickDetector",&UseGeThickDetector,"UseGeThickDetector/I");
    e17011_siminfo->Branch("UseGeThinDetector",&UseGeThinDetector,"UseGeThinDetector/I");
    e17011_siminfo->Branch("GeThickDetectorThickness",&GeThickDetectorThickness,"GeThickDetectorThickness/D");
    e17011_siminfo->Branch("GeThinDetectorThickness",&GeThinDetectorThickness,"GeThinDetectorThickness/D");
    e17011_siminfo->Branch("GeThickDetectorRadius",&GeThickDetectorRadius,"GeThickDetectorRadius/D");
    e17011_siminfo->Branch("GeThinDetectorRadius",&GeThinDetectorRadius,"GeThinDetectorRadius/D");
    e17011_siminfo->Branch("GeThickx",&GeThickx,"GeThickx/D");
    e17011_siminfo->Branch("GeThicky",&GeThicky,"GeThicky/D");
    e17011_siminfo->Branch("GeThickz",&GeThickz,"GeThickz/D");
    e17011_siminfo->Branch("GeThinx",&GeThinx,"GeThinx/D");
    e17011_siminfo->Branch("GeThiny",&GeThiny,"GeThiny/D");
    e17011_siminfo->Branch("GeThinz",&GeThinz,"GeThinz/D");
    e17011_siminfo->Branch("GeSpacing",&GeSpacing,"GeSpacing/D");
  }
  e17011_siminfo->Fill();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::BeginOfRunAction(const G4Run */*aRun*/, G4int rank)
{
  //  cout<<"e17011_simAnalysisExample BeginOfRunAction"<<endl;
  eventcounter = 0;
  // stripwidthx = 5; //units of mm //3 mm for 60 mm diameter
  // stripwidthy = 5; //untis of mm //3 mm for 60 mm diameter
  stripwidthx = 2; //units of mm //3 mm for 60 mm diameter
  stripwidthy = 2; //untis of mm //3 mm for 60 mm diameter
  stripwidthz = (GeThickDetectorThickness/mm)/2.;
  decaycount = 0;
  decaycount2 = 0;
  number0=0;
  numberother=0;
  gamma189count =0;
  clover_hit =0;
  electroncount =0;
  gammacount=0;
  implantcorr=0;
  goode = 0;
  goodestrip = 0;
  decay_in_same_imp_strip=0;
  decay_in_1_strip_away=0;
  decay_1_away_xy=0;
  avgimpx=0;
  avgimpy=0;
  avgdecx=0;
  avgdecy=0;
  decayelectrons=0;
  avgoffbygreaterthan1x=0;
  avgoffbygreaterthan1y=0;
  offcount=0;
  Mult1both=0;
  Mult1front=0;
  Mult1back=0;
  Mult2both=0;
  Mult2front=0;
  Mult2back=0;
  Multgreaterthan2=0;
  AlgorithmMiss=0;
  Mult1bothTotal=0;
  Mult1frontTotal=0;
  Mult1backTotal=0;
  Mult2bothTotal=0;
  Mult2frontTotal=0;
  Mult2backTotal=0;
  Multgreaterthan2Total=0;
  Multgreaterthan2Total=0;
  Mult1bothNoLeft=0;
  Mult1frontNoLeft=0;
  Mult1backNoLeft=0;
  Mult2bothNoLeft=0;
  Mult2frontNoLeft=0;
  Mult2backNoLeft=0;
  Multgreaterthan2NoLeft=0;
  Mult1both189=0;
  Mult1front189=0;
  Mult1back189=0;
  Mult2both189=0;
  Mult2front189=0;
  Mult2back189=0;
  Multgreaterthan2189=0;
  Mult1bothFound189=0;
  Mult1frontFound189=0;
  Mult1backFound189=0;
  Mult2bothFound189=0;
  Mult2frontFound189=0;
  Mult2backFound189=0;
  Multgreaterthan2Found189=0;
  twobytwoA=0;
  twobytwoB=0;
  twobytwoAcorrect=0;
  twobytwoBcorrect=0;
  twobythree=0;
  twobythreecorrect=0;

#if defined (G4ANALYSIS_USE_ROOT)

  stringstream fil;
  fil << dirname << filename << "_" << rank << ".root";
  G4String filenam = fil.str();
  file1 = new TFile(filenam.c_str(),"RECREATE");

  filecount++;
  G4cout << "starting run " << filenam << " " << filecount << G4endl;

  //scint_zpixel = 25.; //length of a pixel in the z direction in mm.  Currently a 25 mm pixel length

  if (gSystem) gSystem->ProcessEvents();


  timer = new TStopwatch();


  // G4cout << "starting run 2" << G4endl;

  e17011_simtree = (TTree *)gROOT->FindObject("e17011_simtree");
  if(e17011_simtree);
  else
    {
      e17011_simtree = new TTree("e17011_simtree","DSSD hit variables");
      e17011_simtree->Branch("e17011_sim_energy",&e17011_sim_energy,"e17011_sim_energy/D");
      e17011_simtree->Branch("totE",&totE,"totE/D");
      e17011_simtree->Branch("beam_energy",&beam_energy,"beam_energy/D");
      e17011_simtree->Branch("pathlength",&pathlength,"pathlength/D");
      e17011_simtree->Branch("totaledep",&totaledep,"totaledep/D");
      e17011_simtree->Branch("tedepx",&tedepx,"tedepx[20]/D");
      e17011_simtree->Branch("tedepy",&tedepy,"tedepy[20]/D");

    }

  //   // Energy deopisitions into each crystal
  //   for(int i = 0; i<68; i++){
  //       stringstream temptitle;
  //       temptitle << "Energy" << i;
  //       string title = temptitle.str();
  //       hCrystalEdep[i] = (TH1D*)gROOT->FindObject(title.c_str());
  //       if(hCrystalEdep[i]) hCrystalEdep[i]->Reset();
  //       else{
  // 	  hCrystalEdep[i] = new TH1D(title.c_str(),title.c_str(),4000,0,4000);
  //       }
  //   }
  //  if(hEvsStripx) hEvsStripx->Reset();
  //   else {
  //     hEvsStripx = new TH2D("hEvsStripx","hEvsStripx", 17,0,17,10000,0,10000);
  //   }

  //   if(hPositionOfGamma) hPositionOfGamma->Reset();
  //   else {
  //     hPositionOfGamma = new TH2D("hPositionOfGamma","hPositionOfGamma", 17,0,17,17,0,17);
  //   }

  //   if(hFinalPosition) hFinalPosition->Reset();
  //   else {
  //     hFinalPosition = new TH2D("hFinalPosition","hFinalPosition", 17,0,17,17,0,17);
  //   }

  //   if(hEvsStripy) hEvsStripy->Reset();
  //   else {
  //     hEvsStripy = new TH2D("hEvsStripy","hEvsStripy", 17,0,17,10000,0,10000);
  //   } 
  
  //   if(hEnergyvsMulf) hEnergyvsMulf->Reset();
  //   else {
  //     hEnergyvsMulf = new TH2D("hEnergyvsMulf","hEnergyvsMulf", 10,0,10,10000,0,10000);
  //   }

  //   if(hEnergyvsMulb) hEnergyvsMulb->Reset();
  //   else {
  //     hEnergyvsMulb = new TH2D("hEnergyvsMulb","hEnergyvsMulb", 10,0,10,10000,0,10000);
  //   }

  //   if(hMultvsE) hMultvsE->Reset();
  //   else {
  //     hMultvsE = new TH2D("hMultvsE","hMultvsE", 10000,0,10000,30,0,30);
  //   }

  //   if(hMultvsEcountf2) hMultvsEcountf2->Reset();
  //   else {
  //     hMultvsEcountf2 = new TH2D("hMultvsEcountf2","hMultvsEcountf2", 10000,0,10000,30,0,30);
  //   }
  //   if(hTracklengthvsE) hTracklengthvsE->Reset();
  //   else {
  //     hTracklengthvsE = new TH2D("hTracklengthvsE","hTracklengthvsE", 10000,0,10000,200,0,20);
  //   }
  
  //   if(hTrackLengLowE) hTrackLengLowE->Reset();
  //   else {
  //     hTrackLengLowE = new TH1D("hTrackLengLowE","hTrackLengLowE",200,0,20);
  //   }

  //   if(hTrackLengHighestE) hTrackLengHighestE->Reset();
  //   else {
  //     hTrackLengHighestE = new TH1D("hTrackLengHighestE","hTrackLengHighestE",200,0,20);
  //   }

  //   if(hTrackLengHighE) hTrackLengHighE->Reset();
  //   else {
  //     hTrackLengHighE = new TH1D("hTrackLengHighE","hTrackLengHighE",200,0,20);
  //   }
  
  //   if(hHighETotalE) hHighETotalE->Reset();
  //   else {
  //     hHighETotalE = new TH1D("hHighETotalE","hHighETotalE", 10000,0,10000);
  //   }

  //   if(hHighestETotalE) hHighestETotalE->Reset();
  //   else {
  //     hHighestETotalE = new TH1D("hHighestETotalE","hHighestETotalE", 10000,0,10000);
  //   }

  //   if(hLowETotalE) hLowETotalE->Reset();
  //   else {
  //     hLowETotalE = new TH1D("hLowETotalE","hLowETotalE", 10000,0,10000);
  //   }
  
  //   if(hEnergyImplantPixelMult1BothSides) hEnergyImplantPixelMult1BothSides->Reset();
  //   else {
  //     hEnergyImplantPixelMult1BothSides = new TH1D("hEnergyImplantPixelMult1BothSides","hEnergyImplantPixelMult1BothSides", 10000,0,10000);
  //   }
  
  //   if(hEnergyImplantPixel) hEnergyImplantPixel->Reset();
  //   else {
  //     hEnergyImplantPixel = new TH1D("hEnergyImplantPixel","hEnergyImplantPixel", 10000,0,10000);
  //   }

  //   if(hChangeinZlow) hChangeinZlow->Reset();
  //   else {
  //     hChangeinZlow = new TH1D("hChangeinZlow","hChangeinZlow",100,-5,5);
  //   }

  //   if(hChangeinZhigh) hChangeinZhigh->Reset();
  //   else {
  //     hChangeinZhigh = new TH1D("hChangeinZhigh","hChangeinZhigh",100,-5,5);
  //   }
  
  //   if(hHighEFractionInImplant) hHighEFractionInImplant->Reset();
  //   else {
  //     hHighEFractionInImplant = new TH1D("hHighEFractionInImplant","hHighEFractionInImplant", 100,0,1);
  //   }

  // if(hHighestEFractionInImplant) hHighestEFractionInImplant->Reset();
  //   else {
  //     hHighestEFractionInImplant = new TH1D("hHighestEFractionInImplant","hHighestEFractionInImplant", 100,0,1);
  //   }

  //  if(hLowEFractionInImplant) hLowEFractionInImplant->Reset();
  //   else {
  //     hLowEFractionInImplant = new TH1D("hLowEFractionInImplant","hLowEFractionInImplant", 100,0,1);
  //   }

  if(hmultvsmult) hmultvsmult->Reset();
  else {
    hmultvsmult = new TH2D("hmultvsmult","hmultvsmult", 40,0,40,40,0,40);
  }

  // if(hElectronMult) hElectronMult->Reset();
  // else {
  //   hElectronMult = new TH1D("hElectronMult","hElectronMult", 17,0,17);
  // }

  // if(hGammaMult) hGammaMult->Reset();
  // else {
  //   hGammaMult = new TH1D("hGammaMult","hGammaMult", 17,0,17);
  // }
  
  // if(hTotalMult) hTotalMult->Reset();
  // else {
  //   hTotalMult = new TH1D("hTotalMult","hTotalMult", 17,0,17);
  // }

  // if(haverageEvsmult) haverageEvsmult->Reset();
  // else {
  //   haverageEvsmult = new TH2D("haverageEvsmult","haverageEvsmult", 17,0,17,5000,0,5000);
  // }

  if(hStripMultx) hStripMultx->Reset();
  else {
    hStripMultx = new TH1D("hStripMultx","hStripMultx",40,0,40);
  }

  if(hStripMulty) hStripMulty->Reset();
  else {
    hStripMulty = new TH1D("hStripMulty","hStripMulty",40,0,40);
  }

  // if(hEnergyDeposit) hEnergyDeposit->Reset();
  // else {
  //   hEnergyDeposit = new TH1D("hEnergyDeposit","hEnergyDeposit",5000,0,/*15000*/5000);
  // }

  // if(hPixelAndAlgorithmDifference) hPixelAndAlgorithmDifference->Reset();
  //   else {
  //     hPixelAndAlgorithmDifference = new TH1D("hPixelAndAlgorithmDifference","hPixelAndAlgorithmDifference",5000,0,/*15000*/5000);
  //   }

  // if(hEnergyDepositmult1both) hEnergyDepositmult1both->Reset();
  //   else {
  //     hEnergyDepositmult1both = new TH1D("hEnergyDepositmult1both","hEnergyDepositmult1both",5000,0,/*15000*/5000);
  //   }

  // if(hEnergyDepositmult1front) hEnergyDepositmult1front->Reset();
  //   else {
  //     hEnergyDepositmult1front = new TH1D("hEnergyDepositmult1front","hEnergyDepositmult1front",5000,0,/*15000*/5000);
  //   }

  // if(hEnergyDepositmult1back) hEnergyDepositmult1back->Reset();
  //   else {
  //     hEnergyDepositmult1back = new TH1D("hEnergyDepositmult1back","hEnergyDepositmult1back",5000,0,/*15000*/5000);
  //   }

  // if(hEnergyDepositmult2both) hEnergyDepositmult2both->Reset();
  //   else {
  //     hEnergyDepositmult2both = new TH1D("hEnergyDepositmult2both","hEnergyDepositmult2both",5000,0,/*15000*/5000);
  //   }

  // if(hEnergyDepositmult2front) hEnergyDepositmult2front->Reset();
  //   else {
  //     hEnergyDepositmult2front = new TH1D("hEnergyDepositmult2front","hEnergyDepositmult2front",5000,0,/*15000*/5000);
  //   }

  // if(hEnergyDepositmult2back) hEnergyDepositmult2back->Reset();
  //   else {
  //     hEnergyDepositmult2back = new TH1D("hEnergyDepositmult2back","hEnergyDepositmult2back",5000,0,/*15000*/5000);
  //   }

  // if(hEnergyDepositmultgreaterthan2) hEnergyDepositmultgreaterthan2->Reset();
  //   else {
  //     hEnergyDepositmultgreaterthan2 = new TH1D("hEnergyDepositmultgreaterthan2","hEnergyDepositmultgreaterthan2",5000,0,/*15000*/5000);
  //   }

  // if(hEnergyDepositAllPlusLeftovers) hEnergyDepositAllPlusLeftovers->Reset();
  //   else {
  //     hEnergyDepositAllPlusLeftovers = new TH1D("hEnergyDepositAllPlusLeftovers","hEnergyDepositAllPlusLeftovers",5000,0,/*15000*/5000);
  //   }

  // if(hPercentMatch) hPercentMatch->Reset();
  //   else {
  //     hPercentMatch = new TH1D("hPercentMatch","hPercentMatch",200,0,1);
  //   }

  if(hEnergyDepositcal) hEnergyDepositcal->Reset();
  else {
    hEnergyDepositcal = new TH1D("hEnergyDepositcal","hEnergyDepositcal",5000,0,5000);
  }
  
  if(hEnergyDepositcal_gammagated) hEnergyDepositcal_gammagated->Reset();
  else {
    hEnergyDepositcal_gammagated = new TH1D("hEnergyDepositcal_gg","hEnergyDepositcal_gg",5000,0,5000);
  }

  if(hEnergyDepositcal_gammagated_antigated) hEnergyDepositcal_gammagated_antigated->Reset();
   else {
    hEnergyDepositcal_gammagated_antigated = new TH1D("hEnergyDepositcal_gg_ag","hEnergyDepositcal_gg_ag",5000,0,5000);
  }

  //max strip energies
  if(hStripMaxEnergyx) hStripMaxEnergyx->Reset();
  else {
    hStripMaxEnergyx = new TH1D("hStripMaxEnergyx","hStripMaxEnergyx",5000,0,5000);
  }

  if(hStripMaxEnergyy) hStripMaxEnergyy->Reset();
  else {
    hStripMaxEnergyy = new TH1D("hStripMaxEnergyy","hStripMaxEnergyy",5000,0,5000);
  }

  
  //   if(hEnergyDepositTest) hEnergyDepositTest->Reset();
  //   else {
  //     hEnergyDepositTest = new TH1D("hEnergyDepositTest","hEnergyDepositTest",5000,0,/*15000*/5000);
  //   }
  //   if(hEnergyDepositAlgorithmTest) hEnergyDepositAlgorithmTest->Reset();
  //   else {
  //     hEnergyDepositAlgorithmTest = new TH1D("hEnergyDepositAlgorithmTest","hEnergyDepositAlgorithmTest",5000,0,/*15000*/5000);
  //   }
  // if(hEnergyDepositx) hEnergyDepositx->Reset();
  //   else {
  //     hEnergyDepositx = new TH1D("hEnergyDepositx","hEnergyDepositx",5000,0,/*15000*/5000);
  //   }
  // if(hEnergyDeposity) hEnergyDeposity->Reset();
  //   else {
  //     hEnergyDeposity = new TH1D("hEnergyDeposity","hEnergyDeposity",5000,0,/*15000*/5000);
  //   }

  //  if(hNumberof188) hNumberof188->Reset();
  //   else {
  //     hNumberof188 = new TH1D("hNumberof188","hNumberof188",10,0,10);
  //   }
 
  //  if(hEnergyDepositLeftOverx) hEnergyDepositLeftOverx->Reset();
  //   else {
  //     hEnergyDepositLeftOverx = new TH1D("hEnergyDepositLeftOverx","hEnergyDepositLeftOverx",5000,0,/*15000*/5000);
  //   }

  //  if(hEnergyDepositLeftOvery) hEnergyDepositLeftOvery->Reset();
  //   else {
  //     hEnergyDepositLeftOvery = new TH1D("hEnergyDepositLeftOvery","hEnergyDepositLeftOvery",5000,0,/*15000*/5000);
  //   }

  //  if(hEnergyDepositLowz) hEnergyDepositLowz->Reset();
  //  else {
  //      hEnergyDepositLowz = new TH1D("hEnergyDepositLowz","hEnergyDepositLowz",5000,0,15000);
  //  }
  //  if(hEnergyDepositHiz) hEnergyDepositHiz->Reset();
  //  else {
  //      hEnergyDepositHiz = new TH1D("hEnergyDepositHiz","hEnergyDepositHiz",5000,0,15000);
  //  }
  
  //  if(hEnergyDepositMaxStrip) hEnergyDepositMaxStrip->Reset();
  //  else {
  //    hEnergyDepositMaxStrip = new TH1D("hEnergyDepositMaxStrip","hEnergyDepositMaxStrip",5000,0,/*15000*/5000);
  //  }

  //  if(hEnergyDepositAllStrips) hEnergyDepositAllStrips->Reset();
  //  else {
  //    hEnergyDepositAllStrips = new TH1D("hEnergyDepositAllStrips","hEnergyDepositAllStrips",5000,0,/*15000*/5000);
  //  }
  
  //  if(helectronx) helectronx->Reset();
  //  else {
  //    helectronx = new TH1D("helectronx","helectronx",5000,0,/*15000*/5000);
  //  }
  
  //  if(helectrony) helectrony->Reset();
  //  else {
  //    helectrony = new TH1D("helectrony","helectrony",5000,0,/*15000*/5000);
  //  }
  
  //  if(hnotelectronx) hnotelectronx->Reset();
  //  else {
  //    hnotelectronx = new TH1D("hnotelectronx","hnotelectronx",5000,0,/*15000*/5000);
  //  }

  // if(hnotelectrony) hnotelectrony->Reset();
  //  else {
  //    hnotelectrony = new TH1D("hnotelectrony","hnotelectrony",5000,0,/*15000*/5000);
  //  }


  if(hBetaDecay) hBetaDecay->Reset();
  else {
    hBetaDecay = new TH1D("hBetaDecay","hBetaDecay",1500,0,1500);
  }

  // if(hEnergyDepositElectron) hEnergyDepositElectron->Reset();
  // else {
  //     hEnergyDepositElectron = new TH1D("hEnergyDepositElectron","hEnergyDepositElectron",5000,0,15000);
  // }
  // if(hEnergyDepositElectron2D) hEnergyDepositElectron2D->Reset();
  // else {
  //   hEnergyDepositElectron2D = new TH2D("hEnergyDepositElectron2D","hEnergyDepositElectron2D",100,-50,50,100,-50,50);
  // }

  // if(hEnergyDepositGamma) hEnergyDepositGamma->Reset();
  // else {
  //     hEnergyDepositGamma = new TH1D("hEnergyDepositGamma","hEnergyDepositGamma",5000,0,15000);
  // }
  // if(hEnergyDepositGamma_nocenter) hEnergyDepositGamma_nocenter->Reset();
  // else {
  //     hEnergyDepositGamma_nocenter = new TH1D("hEnergyDepositGamma_nocenter","hEnergyDepositGamma_nocenter",8192,0,8192);
  // }

  // if(hEnergyDepositPixel) hEnergyDepositPixel->Reset();
  // else {
  //     hEnergyDepositPixel = new TH1D("hEnergyDepositPixel","hEnergyDepositPixel",5000,0,15000);
  // }

  /********/
  //   for(G4int x=0;x<100;x++){
  //     for(G4int y=0;y<100;y++){
  //       stringstream temptitle;
  //       temptitle << "hEnergyDepositPixelGamma" << x << "_"<<y;
  //       string title = temptitle.str();
  //       if(hEnergyDepositPixelGamma[x][y]) hEnergyDepositPixelGamma[x][y]->Reset();
  //       else {
  // 	hEnergyDepositPixelGamma[x][y] = new TH1D(title.c_str(),title.c_str(),2000,0,2000);
  //       }
  //     }
  //   }


  // if(hEnergyDepositOther) hEnergyDepositOther->Reset();
  // else {
  //     hEnergyDepositOther = new TH1D("hEnergyDepositOther","hEnergyDepositOther",8192,0,8192);
  // }

  // if(hEnergyDepositGamma2D) hEnergyDepositGamma2D->Reset();
  // else {
  //   hEnergyDepositGamma2D = new TH2D("hEnergyDepositGamma2D","hEnergyDepositGamma2D",100,0,100,100,0,100);
  // }

  // if(hEnergyDepositLocMax2D) hEnergyDepositLocMax2D->Reset();
  // else {
  //   hEnergyDepositLocMax2D = new TH2D("hEnergyDepositLocMax2D","hEnergyDepositLocMax2D",100,0,100,100,0,100);
  // }

  // if(hEnergyDepositStripx) hEnergyDepositStripx->Reset();
  // else {
  //     hEnergyDepositStripx = new TH1D("hEnergyDepositStripx","hEnergyDepositStripx",8192,0,8192);
  // }
  // if(hEnergyDepositStripy) hEnergyDepositStripy->Reset();
  // else {
  //     hEnergyDepositStripy = new TH1D("hEnergyDepositStripy","hEnergyDepositStripy",8192,0,8192);
  // }


  /*if(hStripMultx) hStripMultx->Reset();
    else {
    hStripMultx = new TH1D("hStripMultx","hStripMultx",100,0,100);
    }*/

  // beam energy 
  // hBeam_energy = (TH1D *)gROOT->FindObject("hBeam_energy");
  // if (hBeam_energy) hBeam_energy->Reset();
  // else
  //   {
  //     hBeam_energy = new TH1D("hBeam_energy","Beam energy", 4000, 0.0, 4000.0); // in keV
  //     if (hBeam_energy) hBeam_energy->SetXTitle("[keV]");
  //   }

  // hPrimaryPosin = (TH2D *)gROOT->FindObject("hPrimaryPosin");
  // if (hPrimaryPosin) hPrimaryPosin->Reset();
  // else{
  //     hPrimaryPosin = new TH2D("hPrimaryPosin","hPrimaryPosin",250,-125,125,250,-125,125);
  // }

  // hPrimaryPos = (TH2D *)gROOT->FindObject("hPrimaryPos");
  // if (hPrimaryPos) hPrimaryPos->Reset();
  // else{
  //     hPrimaryPos = new TH2D("hPrimaryPos","hPrimaryPos",250,-125,125,250,-125,125);
  // }

  // hGammaGamma = (TH2D *)gROOT->FindObject("hGammaGamma");
  // if (hGammaGamma) hGammaGamma->Reset();
  // else{
  //     hGammaGamma = new TH2D("hGammaGamma","hGammaGamma",2000,0,2000,2000,0,2000);
  // }
  
  // hPixelGamma = (TH2D *)gROOT->FindObject("hPixelGamma");
  // if (hPixelGamma) hPixelGamma->Reset();
  // else{
  //     hPixelGamma = new TH2D("hPixelGamma","hPixelGamma",20,0,19,20,0,19);
  // }

  hPixele = (TH2D *)gROOT->FindObject("hPixele");
  if (hPixele) hPixele->Reset();
  else{
    hPixele = new TH2D("hPixele","hPixele",40,0,40,40,0,40);
  }

  hPixelFe = (TH2D *)gROOT->FindObject("hPixelFe");
  if (hPixelFe) hPixelFe->Reset();
  else{
    hPixelFe = new TH2D("hPixelFe","hPixelFe",40,0,40,40,0,40);
  }

  // hgammaelectrontime= (TH1D *)gROOT->FindObject("hgammaelectrontime");
  // if (hgammaelectrontime) hgammaelectrontime->Reset();
  // else{
  //   hgammaelectrontime = new TH1D("hgammaelectrontime","hgammaelectrontime", 200, 0, 200); // in degrees
  // }
 
  hBetaDist= (TH1D *)gROOT->FindObject("hBetaDist");
  if (hBetaDist) hBetaDist->Reset();
  else{
    hBetaDist = new TH1D("hBetaDist","hBetaDist", 10000, 0, 10000); 
  }

  hBetaDistAl= (TH1D *)gROOT->FindObject("hBetaDistAl");
  if (hBetaDistAl) hBetaDistAl->Reset();
  else{
    hBetaDistAl = new TH1D("hBetaDistAl","hBetaDistAl", 60, 0, 15); 
  }

  // // gamma angular correlation from simulation
  //   hgangcortheory = (TH1D *)gROOT->FindObject("hgangcortheory");
  //   if (hgangcortheory) hgangcortheory->Reset();
  //   else{
  //     hgangcortheory = new TH1D("hgangcortheory","hgangcortheory", 360, 0.0, 360.0); // in degrees
  //   }
  //   // gamma angular correlation from "experiment"
  //   hgangcorexper = (TH1D *)gROOT->FindObject("hgangcorexper");
  //   if (hgangcorexper) hgangcorexper->Reset();
  //   else{
  //     hgangcorexper = new TH1D("hgangcorexper","hgangcorexper", 360, 0.0, 360.0); // in degrees
  //   }

  for(G4int x=0;x<16;x++){
    stringstream temptitle;
    temptitle << "hEnergyDepositSega_" << x;
    string title = temptitle.str();
    if(hEnergyDepositSega[x]) hEnergyDepositSega[x]->Reset();
    else {
      hEnergyDepositSega[x] = new TH1D(title.c_str(),title.c_str(),12000,0,12000);
    }
  }
  // gamma angular correlation from "experiment"
  hSega = (TH1D *)gROOT->FindObject("hSega");
  if (hSega) hSega->Reset();
  else{
    hSega = new TH1D("hSega","hSega", 12000, 0.0, 12000);
  }

  // sega addback
  hSegaab = (TH1D *)gROOT->FindObject("hSegaab");
  if (hSegaab) hSegaab->Reset();
  else{
    hSegaab = new TH1D("hSegaab","hSegaab", 12000, 0.0, 12000);
  }
  
  //BC - add loop for filling LaBr3 histograms
  numLaBr3 = 22;
  for(G4int x = 0; x < numLaBr3; x++){
    stringstream temptitle;
    temptitle << "hEnergyDepositLaBr3_" << x;
    string title = temptitle.str();
    if(hEnergyDepositLaBr3[x]) hEnergyDepositLaBr3[x]->Reset();
    else{
      //from http://root.cern.ch/root/html/TH1D.html
      //TH1D(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup)
      hEnergyDepositLaBr3[x] = new TH1D(title.c_str(),title.c_str(),12000,0,12000);
    }
  }
  // gamma angular correlation from "experiment"
  hLaBr3 = (TH1D *)gROOT->FindObject("hLaBr3");
  if (hLaBr3) hLaBr3->Reset();
  else{
    hLaBr3 = new TH1D("hLaBr3","hLaBr3", 12000, 0.0, 12000); // in degrees
  }
  
  //LaBr3 gamma-gamma
  hLaBr3gg = (TH2D*)gROOT->FindObject("hLaBr3gg");
  if (hLaBr3gg) hLaBr3gg->Reset();
  else{
    hLaBr3gg = new TH2D("hLaBr3gg", "hLaBr3gg", 2500, 0, 2500, 2500, 0, 2500);
  }
  
  //BC - add Getest histogram
  hGetest = (TH1D *)gROOT->FindObject("hGetest");
  if (hGetest) hGetest->Reset();
  else{
    hGetest = new TH1D("hGetest","hGetest", 6000, 0, 3000);
  }
  
  
  //BC - add EXOtest histogram
  hEXOtest = (TH1D *)gROOT->FindObject("hEXOtest");
  if (hEXOtest) hEXOtest->Reset();
  else{
    hEXOtest = new TH1D("hEXOtest","hEXOtest", 6000, 0, 3000);
  }
  
  //BC - EJ204 Scintillator histograms
  //just a single one for now, which represents dynode signal
  hEJ204Scint = (TH1D *)gROOT->FindObject("hEJ204Scint");
  if (hEJ204Scint) hEJ204Scint->Reset();
  else{
    hEJ204Scint = new TH1D("hEJ204Scint","hEJ204Scint", 3000, 0, 3000);
  }


  //BC - CeBr3 Scintillator histograms
  //just a single one for now, which represents dynode signal
  hCeBr3Scint = (TH1D *)gROOT->FindObject("hCeBr3Scint");
  if (hCeBr3Scint) hCeBr3Scint->Reset();
  else{
    hCeBr3Scint = new TH1D("hCeBr3Scint","hCeBr3Scint", 3000, 0, 3000);
  }

  
  

  for(G4int x=0;x<72;x++){
    stringstream temptitle;
    temptitle << "hEnergyDepositClover_" << x;
    string title = temptitle.str();
    if(hEnergyDepositClover[x]) hEnergyDepositClover[x]->Reset();
    else {
      hEnergyDepositClover[x] = new TH1D(title.c_str(),title.c_str(),5500,0,5500);
    }
  }

  for(G4int x=0;x<40;x++){
    stringstream temptitle;
    temptitle << "hEnergyDepositGeDSSDx_" << x;
    string title = temptitle.str();
    if(hEnergyDepositGeDSSDx[x]) hEnergyDepositGeDSSDx[x]->Reset();
    else {
      hEnergyDepositGeDSSDx[x] = new TH1D(title.c_str(),title.c_str(),5500,0,5500);
    }
  }

  for(G4int x=0;x<40;x++){
     stringstream temptitle;
     temptitle << "hEnergyDepositGeDSSDy_" << x;
     string title = temptitle.str();
     if(hEnergyDepositGeDSSDy[x]) hEnergyDepositGeDSSDy[x]->Reset();
     else {
       hEnergyDepositGeDSSDy[x] = new TH1D(title.c_str(),title.c_str(),5500,0,5500);
     }
   }

  for(G4int x=0;x<40;x++){
     stringstream temptitle;
     temptitle << "hEnergyDepositGeDSSDy_gg_" << x;
     string title = temptitle.str();
     if(hEnergyDepositGeDSSDy_gg[x]) hEnergyDepositGeDSSDy_gg[x]->Reset();
     else {
       hEnergyDepositGeDSSDy_gg[x] = new TH1D(title.c_str(),title.c_str(),5500,0,5500);
     }
   }
  

  // for(G4int x=0;x<16;x++){
  //   for(G4int y=0;y<16;y++){
  //    stringstream temptitle;
  //    temptitle << "hEnergyDepositGeDSSDxy_" << x<<"_"<<y;
  //    string title = temptitle.str();
  //    if(hEnergyDepositGeDSSDxy[x][y]) hEnergyDepositGeDSSDxy[x][y]->Reset();
  //    else {
  //      hEnergyDepositGeDSSDxy[x][y] = new TH1D(title.c_str(),title.c_str(),10000,0,10000);
  //    }
  //   }
  //  }

  // for(G4int x=0;x<16;x++){
  //     stringstream temptitle;
  //     temptitle << "hEnergyDepositNoAlgorithmGeDSSDx_" << x;
  //     string title = temptitle.str();
  //     if(hEnergyDepositNoAlgorithmGeDSSDx[x]) hEnergyDepositNoAlgorithmGeDSSDx[x]->Reset();
  //     else {
  //       hEnergyDepositNoAlgorithmGeDSSDx[x] = new TH1D(title.c_str(),title.c_str(),5500,0,5500);
  //     }
  //   }

  // for(G4int x=0;x<16;x++){
  //   stringstream temptitle;
  //   temptitle << "hEnergyDepositNoAlgorithmGeDSSDy_" << x;
  //   string title = temptitle.str();
  //   if(hEnergyDepositNoAlgorithmGeDSSDy[x]) hEnergyDepositNoAlgorithmGeDSSDy[x]->Reset();
  //   else {
  //     hEnergyDepositNoAlgorithmGeDSSDy[x] = new TH1D(title.c_str(),title.c_str(),5500,0,5500);
  //   }
  // }

  // for(G4int x=0;x<16;x++){
  //   stringstream temptitle;
  //   temptitle << "hEdepMultFront" << x;
  //   string title = temptitle.str();
  //   if(hEdepMultFront[x]) hEdepMultFront[x]->Reset();
  //   else {
  //     hEdepMultFront[x] = new TH1D(title.c_str(),title.c_str(),10000,0,10000);
  //   }
  // }
  // for(G4int x=0;x<16;x++){
  //   stringstream temptitle;
  //   temptitle << "hEdepMultBack" << x;
  //   string title = temptitle.str();
  //   if(hEdepMultBack[x]) hEdepMultBack[x]->Reset();
  //   else {
  //     hEdepMultBack[x] = new TH1D(title.c_str(),title.c_str(),10000,0,10000);
  //   }
  // }

  //  hGeDSSD = (TH1D *)gROOT->FindObject("hGeDSSD");
  //  if (hGeDSSD) hGeDSSD->Reset();
  //  else{
  //    hGeDSSD = new TH1D("hGeDSSD","hGeDSSD", 5500, 0.0, 5500); 
  //  }

  //  hmaxtype = (TH1D *)gROOT->FindObject("hmaxtype");
  //  if (hmaxtype) hmaxtype->Reset();
  //  else{
  //    hmaxtype = new TH1D("hmaxtype","hmaxtype", 150, 0, 150); 
  //  }

  //  hEventsSeparateFromElectron = (TH1D *)gROOT->FindObject("hEventsSeparateFromElectron");
  //  if (hEventsSeparateFromElectron) hEventsSeparateFromElectron->Reset();
  //  else{
  //    hEventsSeparateFromElectron = new TH1D("hEventsSeparateFromElectron","hEventsSeparateFromElectron", 5500, 0.0, 5500); 
  //  }

  // hGeDSSDsecondmax = (TH1D *)gROOT->FindObject("hGeDSSDsecondmax");
  //  if (hGeDSSDsecondmax) hGeDSSDsecondmax->Reset();
  //  else{
  //    hGeDSSDsecondmax = new TH1D("hGeDSSDsecondmax","hGeDSSDsecondmax", 5500, 0.0, 5500); 
  //  }

  //  hMult1x = (TH1D *)gROOT->FindObject("hMult1x");
  //  if (hMult1x) hMult1x->Reset();
  //  else{
  //    hMult1x = new TH1D("hMult1x","hMult1x", 5500, 0.0, 5500); 
  //  }

  // hMultmanyx = (TH1D *)gROOT->FindObject("hMultmanyx");
  //  if (hMultmanyx) hMultmanyx->Reset();
  //  else{
  //    hMultmanyx = new TH1D("hMultmanyx","hMultmanyx", 5500, 0.0, 5500); 
  //  }

  //  hBackMultforFrontMult1 = (TH1D *)gROOT->FindObject("hBackMultforFrontMult1");
  //  if (hBackMultforFrontMult1) hBackMultforFrontMult1->Reset();
  //  else{
  //    hBackMultforFrontMult1 = new TH1D("hBackMultforFrontMult1","hBackMultforFrontMult1", 17, 0.0, 17); 
  //  }

  //  hBackMultforFrontMult2 = (TH1D *)gROOT->FindObject("hBackMultforFrontMult2");
  //  if (hBackMultforFrontMult2) hBackMultforFrontMult2->Reset();
  //  else{
  //    hBackMultforFrontMult2 = new TH1D("hBackMultforFrontMult2","hBackMultforFrontMult2", 17, 0.0, 17); 
  //  }



  //look at the energy spectrum for strips that are assigned a strip number greater than the 16 that actually exist
  /* hGeDSSDImaginaryStrips = (TH1D *)gROOT->FindObject("hGeDSSDImaginaryStrips");
     if (hGeDSSDImaginaryStrips) hGeDSSDImaginaryStrips->Reset();
     else{
     hGeDSSDImaginaryStrips = new TH1D("hGeDSSDImaginaryStrips","hGeDSSDImaginaryStrips", 5500, 0.0, 5500); // in degrees
     }*/

  // gamma angular correlation from "experiment"
  hClover = (TH1D *)gROOT->FindObject("hClover");
  if (hClover) hClover->Reset();
  else{
    hClover = new TH1D("hClover","hClover", 5500, 0.0, 5500); // in degrees
  }


  for(G4int x=0;x<18;x++){
    stringstream temptitle;
    temptitle << "hEnergyDepositClover_addback_" << x;
    string title = temptitle.str();
    if(hEnergyDepositClover_addback[x]) hEnergyDepositClover_addback[x]->Reset();
    else {
      hEnergyDepositClover_addback[x] = new TH1D(title.c_str(),title.c_str(),5500,0,5500);
    }
  }
  // gamma angular correlation from "experiment"
  hClover_addback = (TH1D *)gROOT->FindObject("hClover_addback");
  if (hClover_addback) hClover_addback->Reset();
  else{
    hClover_addback = new TH1D("hClover_addback","hClover_addback", 5500, 0.0, 5500); // in degrees
  }

  hClover_arrayadd = (TH1D *)gROOT->FindObject("hClover_arrayadd");
  if (hClover_arrayadd) hClover_arrayadd->Reset();
  else{
    hClover_arrayadd = new TH1D("hClover_arrayadd","hClover_arrayadd", 5500, 0.0, 5500); // in degrees
  }

  // sega tigress singles add
  /* hsega_tigress_s = (TH1D *)gROOT->FindObject("hsega_tigress_s");
     if (hsega_tigress_s) hsega_tigress_s->Reset();
     else{
     hsega_tigress_s = new TH1D("hsega_tigress_s","hsega_tigress_s", 5500, 0.0, 5500); // in degrees
     }

     hsega_tigress = (TH1D *)gROOT->FindObject("hsega_tigress");
     if (hsega_tigress) hsega_tigress->Reset();
     else{
     hsega_tigress = new TH1D("hsega_tigress","hsega_tigress", 5500, 0.0, 5500); // in degrees
     }*/

  //for each file wait 5 minutes before continuing to let other
  //files be created
  //sleep(300);


  //   hgAngCorEffCorrection = (TH1D*)gROOT->FindObject("hgAngCorEffCorrection");
  //   if (hgAngCorEffCorrection) hgAngCorEffCorrection->Reset();
  //   else{
  //     hgAngCorEffCorrection = new TH1D("hgAngCorEffCorrection","hgAngCorEffCorrection",180,0,180);
  //   }

  //   for(G4int x1=2;x1<17;x1++){
  //     for(G4int y1=2;y2<17;y2++){
  //       G4ThreeVector gamma1(x1-9,y1-9,0);

  //       for(G4int x2=2;x2<17;x2++){
  // 	for(G4int y2=2;y2<17;y2++){
  // 	  G4ThreeVector gamma2(x-9,y-9,0);
  // 	  G4double dotprod;
  // 	  dotprod= (gamma1/gamma1.mag()) * (gamma2/gamma2.mag());
  // 	  G4double angle;
  // 	  angle = acos(dotprod)/CLHEP::degree;
  // 	  G4cout << "dotprod " << dotprod << " angle " << angle << G4endl;
	  
  // 	  hgangcorexper->Fill(angle);
  // 	}
  //       }
  //     }
  //   }
  
  eventtimer = new TStopwatch();
 
  OnceAWhileDoIt(true); // do it now

#endif /* defined (G4ANALYSIS_USE_ROOT) */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::EndOfRunAction(const G4Run */*aRun*/, G4int rank)
{
  //cout<<"EndOfRunAction"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  if (gSystem) gSystem->ProcessEvents();

  stringstream fil;
  fil << dirname << filename << "_" << rank << ".root";
  G4String filenam = fil.str();
  file1=gROOT->GetFile(filenam.c_str());
  //file1=(TFile*)gROOT->FindObject(filenam.c_str());
  //file1 = new TFile(filenam.c_str(),"RECREATE");
  //e17011_siminfo->SetDirectory(file1);
  file1->cd();
  gROOT->SetStyle("Pub");
  e17011_siminfo->Write();

  //e17011_simtree->SetDirectory(file1);
  e17011_simtree->Write();
  
  //hEnergyDeposit->Write();
  // hBeam_energy->Write();
  //hPrimaryPos->Write();

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->Divide(3,5);
  for(int i =1;i<41;i++){
    c1->cd(i);
    hEnergyDepositGeDSSDx[i-1]->GetXaxis()->SetLabelSize(0.06);
    hEnergyDepositGeDSSDx[i-1]->GetYaxis()->SetLabelSize(0.06);
    hEnergyDepositGeDSSDx[i-1]->Draw();
  }
  c1->Write();

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->Divide(2,5);
  for(int i =1;i<10;i++){
    c2->cd(i);
    hEnergyDepositClover_addback[i-1]->GetXaxis()->SetLabelSize(0.06);
    hEnergyDepositClover_addback[i-1]->GetYaxis()->SetLabelSize(0.06);
    hEnergyDepositClover_addback[i-1]->Draw();
  }
  c2->cd(10);
  hClover_addback->GetXaxis()->SetLabelSize(0.06);
  hClover_addback->GetYaxis()->SetLabelSize(0.06);
  hClover_addback->Draw();
  c2->Write();

  file1->Write();
  file1->Close();
  timer->Stop();
  timer->Print();

#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt(true); // do it now
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::BeginOfEventAction(const G4Event */*anEvent*/)
{
  // cout<<"e17011_simAnalysisExample BeginOfEventAction"<<endl;
  // G4cout << "***************************** starting event ************************" << G4endl;


#if defined (G4ANALYSIS_USE_ROOT)
  if (gSystem) gSystem->ProcessEvents();
  //G4cout << "here at root" << G4endl;

  if(eventcounter%1000==0){
    eventtimer->Stop();
    //eventtimer->Print();
    eventtimer->Start();
    G4cout << eventcounter << " done " << G4endl;
  }
  eventcounter++;
  

  //   pixelhits = 0;
  //   pmthits = 0;
  gmult = 0;

  e17011_sim_energy           = 0;

  totE          = 0;
  edep          = 0;
  ID            = 0;
  Z             = 0;
  particletype  = "";
  particleName  = "";
  process       = "";
  creatorID     = 0;
  StepNo        = 0;
  gtime         = 0;

  volname       = "";
  volcopyno     = -1;
  e17011_simcopyno   = -1;
  e17011_sim_edep = 0;

  //   psfx = 0;
  //   psfy = 0;
  e17011_simap_len = 0;

  pathlength = 0;

  particletypebit = 0;
  empty = 0;

#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::EndOfEventAction(const G4Event *evt)
{    

  // At the end of an event take the collection of hits
  // inside the DSSD and track individual particles.

#if defined (G4ANALYSIS_USE_ROOT)
  // G4cout << "USE ROOT defined " << G4endl;
  if (gSystem) gSystem->ProcessEvents();
 
  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  //G4int DetectorColID = SDman->GetCollectionID(colNam="DetectorCollection");
  G4int DetectorColID = SDman->GetCollectionID(colNam="GeThickDetectorCollection");
  G4int GeThinDetectorColID = SDman->GetCollectionID(colNam="GeThinDetectorCollection");
  G4int SegaDetectorID = SDman->GetCollectionID(colNam="SegaCollection");
  G4int CloverColID = SDman->GetCollectionID(colNam="CloverCollection");
  //BC - add int that identifies the LaBr3 array
  G4int LaBr3ColID = SDman->GetCollectionID(colNam="LaBr3Collection");
  //BC - add int that identifies the Getest
  G4int GetestColID = SDman->GetCollectionID(colNam="GetestCollection");
  //BC - add int that identifies the EXOtest
  G4int EXOtestColID = SDman->GetCollectionID(colNam="EXOtestCollection");
  //BC - add int that identifies the EJ204
  G4int EJ204ScintColID = SDman->GetCollectionID(colNam="EJ204ScintCollection");
  //BC - add int that identifies the CeBr3
  G4int CeBr3ScintColID = SDman->GetCollectionID(colNam="CeBr3ScintCollection");

  G4HCofThisEvent *HCE = evt->GetHCofThisEvent();
  e17011_simDetectorHitsCollection *DHC = 0;
  e17011_simDetectorHitsCollection *DHCthin = 0;
  e17011_simDetectorHitsCollection *DHCsega = 0;
  e17011_simDetectorHitsCollection *DHCClover = 0;
  //BC - add LaBr3 to detector hits collection
  e17011_simDetectorHitsCollection *DHClabr3 = 0;
  //BC - add Getest to detector hits collection
  e17011_simDetectorHitsCollection *DHCGetest = 0;
  //BC - add EXOtest to detector hits collection
  e17011_simDetectorHitsCollection *DHCEXOtest = 0;
  //BC - add EXOtest to detector hits collection
  e17011_simDetectorHitsCollection *DHCEJ204Scint = 0;
  e17011_simDetectorHitsCollection *DHCCeBr3Scint = 0;
  

  //G4int gammaincenter;
  //gammaincenter = 0;
  if(HCE){
    DHC = (e17011_simDetectorHitsCollection*)(HCE->GetHC(DetectorColID));
    DHCthin = (e17011_simDetectorHitsCollection*)(HCE->GetHC(GeThinDetectorColID));
    DHCsega = (e17011_simDetectorHitsCollection*)(HCE->GetHC(SegaDetectorID));
    DHCClover = (e17011_simDetectorHitsCollection*)(HCE->GetHC(CloverColID));
    //BC - handle HC (hit collection) for LaBr3
    DHClabr3 = (e17011_simDetectorHitsCollection*)(HCE->GetHC(LaBr3ColID));
    //BC - handle hit collection for Getest
    DHCGetest = (e17011_simDetectorHitsCollection*)(HCE->GetHC(GetestColID));
    DHCEXOtest = (e17011_simDetectorHitsCollection*)(HCE->GetHC(EXOtestColID));
    DHCEJ204Scint = (e17011_simDetectorHitsCollection*)(HCE->GetHC(EJ204ScintColID));
    DHCCeBr3Scint = (e17011_simDetectorHitsCollection*)(HCE->GetHC(CeBr3ScintColID));
  }

  if(DHCthin){
    int n_hit = DHCthin->entries();
    for(int i=0; i<n_hit; i++){
      G4String name = (*DHCthin)[i]->GetParticlename();
      G4String parentname = (*DHCthin)[i]->GetParentName();
      //G4double parentenergy = (*DHCthin)[i]->GetParentEnergy();
      //G4double edep = (*DHCthin)[i]->GetEdep()/keV;
      G4ThreeVector pos = (*DHCthin)[i]->GetPostPosition();
      //       G4cout << "particle - " << name << " - came from " 
      // 	     << parentname << " " << parentenergy  
      // 	     << " and deposited " << edep/keV << " keV at " 
      // 	     << G4int(pos.x()/mm + 50.5) << " " << G4int(pos.y()/mm + 50.5) << G4endl;
    }    
  }





  G4double totaledep_stripx[100];
  G4double totaledep_stripy[100];
  G4double totaledep_xy[100][100];
  //G4double GeDSSDedep = 0;
  G4int stripmaxx=0;
  G4int stripmaxy=0;
  G4int stripimpx=0;
  G4int stripimpy=0;
  G4int stripdecx=0;
  G4int stripdecy=0;
  G4double exmax=0;
  G4double eymax=0;

  //G4double decposx=0;
  //G4double decposy=0;
  //G4double impposx=0;
  //G4double impposy=0;

  G4double gammatime=0;
  G4double electrontime=0;
  G4double e=0;
  G4double eAl=0;

  for(int i=0;i<100;i++){
    totaledep_stripx[i] = totaledep_stripy[i] = 0;
    for(int j=0;j<100;j++){
      totaledep_xy[i][j] = 0;
    }
  }
  
  
  for(G4int type=0;type<3;type++){
    for(G4int x=0;x<100;x++){
      for(G4int y=0;y<100;y++){
	for(G4int z=0;z<2;z++){
	  totaledep[type][x][y][z]=0;
	}
      }
    }
  }


  G4ThreeVector gamma1modir;
  G4ThreeVector gamma2modir;

  int parentID;
  double tot=0;
  bool implanted=false;
  double time1=0;
  double time2=0;
  int tvar=0;

  if(DHC){

    //  G4cout<<"Processing GeDSSD"<<G4endl;
    
    int n_hit = DHC->entries(); //get the entries
    //G4String namemax[n_hit];  //array of names
    G4String *namemax = new G4String[n_hit];  //array of names

    int prevstep = 0; 
    
    //loop over hits
    // G4cout<<"Processing GeDSSD, n_hit: "<<n_hit<<G4endl;

    for(int i=0; i<n_hit; i++){
      // for every particle that deposited energy in the GeDSSD determine
      // its name, what beta decay product it came from and the beta decay
      // product's energy.
      G4String name = (*DHC)[i]->GetParticlename();
      G4String parentname = (*DHC)[i]->GetParentName();
      G4double parentenergy = (*DHC)[i]->GetParentEnergy();
      G4ThreeVector parmodir = (*DHC)[i]->GetParentMoDir();
      G4int parentno = (*DHC)[i]->GetParentno();
      parentID=parentno;
      namemax[i]=name;

      // retrieve the energy deposited in this step, the position and the step number
      G4double edep = (*DHC)[i]->GetEdep()/ keV;
      G4ThreeVector pos = (*DHC)[i]->GetPostPosition();
      G4int stepno = (*DHC)[i]->GetStepno();
      G4double globaltime = (*DHC)[i]->GetGtime()/second;
    

      // determine the strip number
      // center strip is centered around .67,0,0 cm
      // detector is centered around x,y = 0,0.  position must 
      // be offset by one -4.5 corrects strip number for this offset
      // radius to ensure that the stripnumx and stripnumy variables are non-negative

      // G4int stripx = G4int((pos.x()/mm + (GeThickDetectorRadius/mm)- (5.0))/(stripwidthx));   5mm strips
      G4int stripx = G4int((pos.x()/mm + (GeThickDetectorRadius/mm)- (2.0))/(stripwidthx));
      G4int stripy = G4int((pos.y()/mm + (GeThickDetectorRadius/mm))/(stripwidthy)); 
                                                          //y is not offset at all

      
      //look for implants
      if(parentID==0){
       	stripimpx=stripx-1;//0.25,0.25 and 0,0 0.25,0 0.25,0.25 in "real" strip numbers
	      stripimpy=stripy-1;
	      hPixelFe->Fill(stripimpx,stripimpy);
	
	      avgimpx+=stripimpx;
	      avgimpy+=stripimpy;
	      implanted=true;
      }

      //Fill the strip arrays
      if(edep<=5000){
	      totaledep_stripx[stripx]+=edep;
	      totaledep_stripy[stripy]+=edep;
      }

      if(globaltime<=0.5 && edep>6 && edep<3000 && globaltime>=.001){//mimic a correlation window of 500m
       	//added energy thresholds-should be in keV, but if not, then I'll have to check again
	      //also added the min implant time

	      if(edep>exmax){
	        exmax=edep;
	        stripmaxx=stripx;
	      }
	
	      if(edep>eymax){
	        eymax=edep;
	        stripmaxy=stripy;
	      }
	
	      //re-create the beta-time spectrum
	      if(tvar==0){
	        time1=globaltime;
	        tvar=1;
	        decayelectrons++;
	      }	 
      }
      
      if(tvar==1){
	      if(globaltime>time1){
	        time2=globaltime;
	      }
      }

      // technically there are no physical strips in the z direction but
      // this is included to see if there is an advantage of trying to determine
      // z position
      //G4int stripz = G4int((pos.z()/mm + ((GeThickDetectorThickness/mm)/2.))/(stripwidthz));

      // This definition of stripz will be used for 2 cm detectors an will devide the GeDSSD into
      // the first 2 mm and then everything else.  Applicable for radioxenon project.
      G4int stripz = -1;
      if( (pos.z()/mm + ((GeThickDetectorThickness/mm)/2.)) < 0.5){
	      stripz = 0;
      }
      else stripz = 1;


      //      G4cout << "zpos calc "<<pos.z()/mm<<" "<<(GeThickDetectorThickness/mm)/2.<<" "<< stripwidthz
      //	     << " " << stripz << G4endl;


      if(name == "gamma" && parentname == "gamma" && stepno == 1){
	      //G4cout << "what is happening " << parentenergy/keV << G4endl;

	      if((parentenergy/keV) > 201.99 && (parentenergy/keV) < 202.01){
	        gamma1modir = parmodir;
	        //G4cout << "found 202 " <<" " << parentenergy << " " << parmodir.mag() << G4endl;
	      }
	      else if((parentenergy/keV) > 197.99 && (parentenergy/keV) < 198.01){
	        gamma2modir = parmodir;
	      //G4cout << " found 198 " << G4endl;
	      }
      }

      //do gamma multiplicities
      if(name == "gamma" && parentname == "gamma"){
	      gmult++;
      }

      if(name == "Co54[0.0]" && stepno == 1){
	      decaycount++;
      }

      if(name=="gamma" && parentno==1){
	      //	hPixelGamma->Fill(stripx,stripy);	
      }
           
      if(name=="gamma" && parentno==0){
	      //	cout<<"gamma"<<endl;
	      gammatime= globaltime;
      }

      if(name=="e-" && parentno==0){
	      //	cout<<"electron"<<endl;
	      electrontime=globaltime;
      }
      
      if(name=="e-" && parentno==1 && e==0){
	      e= (*DHC)[i]->GetKineticEnergy()/ MeV;
	      e = e+ (*DHC)[i]->GetEdep()/ MeV;
	      prevstep=stepno;
	      //this won't quite be the full kinetic energy is 
        //some was lost rather than being deposited in step 0, 
	      //but should make for a good approximation
	      // cout<<e<<endl;
	      // hBetaDist->Fill(e);
      }

      if(name=="e-" && parentno==1 && stepno==n_hit-1 && pos.z()>0){ 
        //pos.z to check which side it's coming out of. <0 1mm side
      	eAl = (*DHC)[i]->GetKineticEnergy()/ MeV;
      }
      /*if(edep>0){
	        hEnergyDeposit->Fill(edep);
	      }*/
      //G4cout << "edep and location " << stepno<< " " << name << " " << parentname << " " << parentenergy
      // << " " << edep << " " << stripx << " " << stripy << " " << stripz << G4endl;
      if(parentname == "e-" ){
	      totaledep[0][stripx][stripy][stripz] += edep;
	      //	hPixele->Fill(stripx,stripy);
      }
      else if(parentname == "gamma"){
	      totaledep[1][stripx][stripy][stripz] += edep;
	      //	hPixelGamma->Fill(stripx,stripy);
      }
      else{
	      totaledep[2][stripx][stripy][stripz] += edep;
	      //	hPixelFe->Fill(stripx,stripy);
	      /*	if(name!="Fe66[0.0]"){
		    stripdecx=stripx;
		    stripdecy=stripy;
		    avgdecx+=stripx;
		    avgdecy+=stripy;
		    }*/
      }
      
    }//end loop over all hits in the event

    //  G4cout<<"Done processing hits"<<G4endl;

    hBetaDist->Fill(e);
    if(eAl>0) hBetaDistAl->Fill(eAl);
    tot=0;

    GeDSSDcal=0;  //to use as calorimeter
    countf =0; //front mult counter reset
    countb =0; //back mult counter reset
    max_energy_x =0;
    max_strip_x =0;
    max_energy_y =0;
    max_strip_y =0;
    
    
    //G4cout<<"Looking at strips"<<G4endl;
    
    //loop over all strips
    for(int x=1;x<41;x++){
      int y=x-1;

      // GeDSSDedep+=totaledep_stripx[x];

      //x strips
      if(totaledep_stripx[x]>0){
	
	      countf = countf + 1; //increment front mult counter
	      //G4cout<<"Front: "<<countf<<G4endl;
	
	      //add a resolution
	      if(totaledep_stripx[x]>300) {
	        // totaledep_stripx[x]= gRandom->Gaus(totaledep_stripx[x],(0.0047*totaledep_stripx[x]/2.35));
          //from 662 before
	        totaledep_stripx[x]= gRandom->Gaus(totaledep_stripx[x],(0.01*totaledep_stripx[x]/2.35));
                               //from 662 after
	      }else{
	      //totaledep_stripx[x] = gRandom->Gaus(totaledep_stripx[x],(0.009*totaledep_stripx[x]/2.35));
                                //from 247.7 before
	        totaledep_stripx[x] = gRandom->Gaus(totaledep_stripx[x],(0.05*totaledep_stripx[x]/2.35));
                                //from123 after
	      }
	
	      //	G4cout<<"Resolutions added to front"<<G4endl;
	
	      //add the energy from the front strips together to obtain calorimeter
	      GeDSSDcal+=totaledep_stripx[x];

	      //fill x strip energy histograms
	      hEnergyDepositGeDSSDx[y]->Fill(totaledep_stripx[x]);
	      //	G4cout<<"Done with front"<<G4endl;

	      //Look for max energy strip number
	      if(totaledep_stripx[x]>=max_energy_x) {
	        max_energy_x = totaledep_stripx[x];
	        max_strip_x = x;
	      }
	       
      }

      // y strips
      if(totaledep_stripy[x]>0){
	
	      countb = countb + 1; //increment back mult counter
	      //	G4cout<<"Back: "<<countb<<G4endl;
		
	      if(totaledep_stripy[x]>300) {//add a resolution
	        totaledep_stripy[x]= gRandom->Gaus(totaledep_stripy[x],(0.0050*totaledep_stripy[x]/2.35));
          //from 662 before
	      }else{
	        totaledep_stripy[x] = gRandom->Gaus(totaledep_stripy[x],(0.03*totaledep_stripy[x]/2.35));
          //from 247.7keV SRM source before
	      }
	      //	G4cout<<"Resoultions done"<<G4endl;

	      //fill y strip energy histograms
	      hEnergyDepositGeDSSDy[y]->Fill(totaledep_stripy[x]);

	      //Look for max energy strip number
	      if(totaledep_stripy[x]>=max_energy_y) {
	        max_energy_y = totaledep_stripy[x];
	        max_strip_y = x;
	      }
	
      }
      //  G4cout<<"Done with back"<<G4endl;
    }
    
    //  G4cout<<"Done Looking at strips"<<G4endl;

    //fill mult histograms
    hStripMultx->Fill(countf);
    hStripMulty->Fill(countb);
    hmultvsmult->Fill(countf,countb);

    //fill max energy histograms
    hStripMaxEnergyx->Fill(max_energy_x);
    hStripMaxEnergyy->Fill(max_energy_y);
    
    //fill calorimeter histogram
    if(GeDSSDcal>0) hEnergyDepositcal->Fill(GeDSSDcal);
    
 
    stripdecx=stripmaxx-1;//to get back to "real" strip numbers
    stripdecy=stripmaxy-1;
    avgdecx+=stripdecx;
    avgdecy+=stripdecy;
    //decayelectrons++;
    //	cout<<"max stripx "<<stripdecx+1<<" max strip y "<<stripdecy+1<<endl;
    if(implanted){
      //corr1
      if(stripdecx==stripimpx && stripdecy==stripimpy){
	      //corr3
	      // if(((stripdecx-stripimpx)>=-1) && ((stripdecx-stripimpx)<=1)){
	      // if(((stripdecy-stripimpy)>=-1) && ((stripdecy-stripimpy)<=1)){
	      //  cout<<"here!"<<endl;
	      hBetaDecay->Fill(time1*1000);
	      hBetaDecay->Fill(time2*1000);
	      hPixele->Fill(stripdecx,stripdecy);
	      //	}
      }
    }
      
    //if((stripmaxx==/*stripimpx*/9)&&(stripmaxy==/*stripimpy*/9)) implantcorr++;
    if((stripimpx==stripdecx)&&(stripimpy==stripdecy)){
      decay_in_same_imp_strip++;
    }else if((stripimpx-stripdecx==-1)||(stripimpx-stripdecx==1)){
      if(stripimpy==stripdecy) decay_in_1_strip_away++;
      if((stripimpy-stripdecy==-1)||(stripimpy-stripdecy==1)) decay_1_away_xy++;
    }else if((stripimpy-stripdecy==-1)||(stripimpy-stripdecy==1)){
      if(stripimpx==stripdecx) decay_in_1_strip_away++;
    }
    
    else{
      if(stripdecx>-1 && stripdecy>-1){//within corr window and in a real strip
	      avgoffbygreaterthan1x+=abs(stripdecx-stripimpx);
	      avgoffbygreaterthan1y+=abs(stripdecy-stripimpy);
	      offcount++;
      }
    }
  }
  else{
    G4cout << "no decay empty event " << empty++ << G4endl;
  }
 
 
  // hits in the Sega detectors
  if(DHCsega){

    G4double segaE[16];
    for(int x=0;x<16;x++){
      segaE[x]=0;
    }
    
    int n_hit = DHCsega->entries();
    for(int i=0; i<n_hit; i++){
      G4String name = (*DHCsega)[i]->GetParticlename();
      G4String parentname = (*DHCsega)[i]->GetParentName();
      //G4double parentenergy = (*DHCsega)[i]->GetParentEnergy();
      G4double edep = (*DHCsega)[i]->GetEdep()/keV;
      G4ThreeVector pos = (*DHCsega)[i]->GetPostPosition();
      G4int copyno = (*DHCsega)[i]->GetVolCopyNo();

      segaE[copyno] += edep;
    }

    //Look for the 478 keV gamma in SeGA
    G4bool SeGA_478 = false;
    for(int i=0; i<16; i++) {
      if(segaE[i]>=475 && segaE[i]<=481) {
	      SeGA_478 = true;
      }
    }

    G4bool Back_478 = false;
    for(int i=0; i<40; i++) {
      if(totaledep_stripy[i]>=475 && totaledep_stripy[i]<=481) {
	      Back_478 = true;
      }
    }

    G4bool Front_478 = false;
    for(int i=0; i<40; i++) {
      if(totaledep_stripx[i]>=474 && totaledep_stripx[i]<=482) {
	      Front_478 = true;
      }
    }

    //look for the pair production from 2033 keV Gamma
    G4bool Back_1011 = false;
    for(int i=0; i<40; i++) {
      if(totaledep_stripy[i]>=1000 && totaledep_stripy[i]<=1020) {
	      Back_1011 = true;
      }
    }

    G4bool Front_1011 = false;
    for(int i=0; i<40; i++) {
      if(totaledep_stripx[i]>=1000 && totaledep_stripx[i]<=1020) {
	      Front_1011 = true;
      }
    }
    
    //Look for 2033 in SeGA
    G4bool SeGA_2033 = false;
    for(int i=0; i<16; i++) {
      if(segaE[i]>=2020 && segaE[i]<=2045) {
	      SeGA_2033 = true;
      }
    }

    //Look for 511 in SeGA
    G4bool SeGA_511 = false;
    for(int i=0; i<16; i++) {
      if(segaE[i]>=508  && segaE[i]<=513) {
	      SeGA_511 = true;
      }
    }
    
    G4double segaadd = 0;
    for(int x=0;x<16;x++){
      if(segaE[x] > 0){

	      //Add resolutions
	      //G4double rando = ((G4UniformRand() * 0.002) -0.001);
	      if(segaE[x]>300) {
	        segaE[x] = gRandom->Gaus(segaE[x],(0.0037*segaE[x]/2.35));
	      }else if(segaE[x]>200){
	        segaE[x] = gRandom->Gaus(segaE[x],(0.012*segaE[x]/2.35));
	      }else{
	        segaE[x] = gRandom->Gaus(segaE[x],(0.022*segaE[x]/2.35));//.022
	      }
          
	      hEnergyDepositSega[x]->Fill(segaE[x]);
	
	      //Gated on 511 in sega
	      if(/*(segaE[x] > 1510.0 && segaE[x] < 1520.0) && */SeGA_511 == true  &&  Front_1011 == false  &&    Back_1011 == false) {
	        if(GeDSSDcal>0) {
	          hEnergyDepositcal_gammagated->Fill(GeDSSDcal);
	        }
	        for(int yy=0; yy<40; yy++) {
	          if(totaledep_stripy[yy]>0) {
	            hEnergyDepositGeDSSDy_gg[yy]->Fill(totaledep_stripy[yy]);
	          }
	        }
	      }

	      //Gated on the 1515 and antigated on the 478
	      if((segaE[x] > 1510.0 && segaE[x] < 1520.0) && (SeGA_478 == false) && (Front_478 == false) && (Back_478 == false) && (SeGA_2033 == false)) {
	        if(GeDSSDcal>0) {
	          hEnergyDepositcal_gammagated_antigated->Fill(GeDSSDcal);
	        }      
	      }
	

	      // //Gated on 1938 keV transition
	      // if(segaE[x] > 1933.0 && segaE[x] < 1943.0) {
	      //   if(GeDSSDcal>0) {
	      //     hEnergyDepositcal_gammagated->Fill(GeDSSDcal);
	      //   }
	      //   for(int yy=0; yy<16; yy++) {
	      //     if(totaledep_stripy[yy]>0) {
	      //       hEnergyDepositGeDSSDy_gg[yy]->Fill(totaledep_stripy[yy]);
	      //     }
	      //   }
	      // }



      }
      segaadd = segaadd+segaE[x];
    }
    
    if(segaadd>0) hSegaab->Fill(segaadd);    
  } // end sega
  
  //BC - Hits in LaBr3 detectors
  if(DHClabr3){
  
    //G4cout << "in LaBr3" << endl;
    
    //define number of LaBr3 detectors on these 2 lines...
    //separated b/c of ISO C++ warning of variable length array if use numLaBr3
    G4double labr3E[numLaBr3];
    
    for(int x=0;x<numLaBr3;x++){
      labr3E[x]=0;
    }
    
    int n_hit = DHClabr3->entries();
    
    //G4cout << "n_hit = " << n_hit << endl;
    
    if(n_hit > 0){
      for(int i=0; i<n_hit; i++){
        G4String name = (*DHClabr3)[i]->GetParticlename();
        G4String parentname = (*DHClabr3)[i]->GetParentName();
        //G4double parentenergy = (*DHClabr3)[i]->GetParentEnergy();
        G4double edep = (*DHClabr3)[i]->GetEdep()/keV;
        G4ThreeVector pos = (*DHClabr3)[i]->GetPostPosition();
        G4int copyno = (*DHClabr3)[i]->GetVolCopyNo();
  
        labr3E[copyno] += edep;
      }
    
      //assuming ~3% resolution when looking for gammas
      //EDIT: I fit the Saing Gobain BrilLanCe spec sheet values to 
      //a 3rd order polynomial on a logarithm scale.  Using this equation
      //for the resolution
      
      // add a resolution to the LaBr3 energy
      //  int count=0;
      for(int x=0;x<numLaBr3;x++){
  	    if(labr3E[x]>0){
  	    
  	      //G4cout << "labr3E[" << x << "] = " << labr3E[x] << endl;
  	      //  count++;
  	      //resolution of 3% across all energy ranges
  	      // cout<<x<<",  ";
  	      double log_E = log10(labr3E[x]);
  	      double log_resolution = -0.0159*pow(log_E,3) + 0.1816*pow(log_E,2) - 1.0939*log_E + 2.4554;
  	      double resolution = pow(10.0, log_resolution)/100; //turn into number b/t 0 and 1

	      //set resolution to equation fit from e14057 data
	      resolution = 2.0085*pow(labr3E[x],-0.629);
  	      
  	      
  	      // if(labr3E[x] < 122){
  	       
  	      //  resolution = 0.06; 
  	       
  	      //   //labr3E[x] = gRandom->Gaus(labr3E[x],(0.03*labr3E[x]/2.355));//3% is for LaBr3
  	      // }
  	      //G4cout << "labr3E[" << x << "] = " << labr3E[x] << ", log_E = " << log_E << ", resolution = " << resolution << G4endl;
  	      
  	      labr3E[x] = gRandom->Gaus(labr3E[x],(resolution*labr3E[x]/2.355));
  	      
  	    }
  	      
  	    //put into histograms
  	    hEnergyDepositLaBr3[x]->Fill(labr3E[x]);//energy of each ind. crystal
	    hLaBr3->Fill(labr3E[x]);//energy from all crystals
	    
	    //Fill 2D LaBr3 gamma-gamma histo
	    for(int y = x+1; y < numLaBr3; y++){
	      if(labr3E[y]>0){
		hLaBr3gg->Fill(labr3E[x], labr3E[y]);
		hLaBr3gg->Fill(labr3E[y], labr3E[x]);
	      }
	    }
  	    
      }
      
      //G4double labr3add = 0;
      //for(int x=0;x<numLaBr3;x++){
        //if(labr3E[x] > 0){
  
  	      //Add resolutions
  	      //G4double rando = ((G4UniformRand() * 0.002) -0.001);
  	      //if(labr3E[x]>300) {
  	        //using 3% resolution for now
  	        //labr3E[x] = gRandom->Gaus(labr3E[x],(0.03*labr3E[x]/2.355)); //FWHM of Gauss = 2.355*sigma
  	      //}else{
  	        //labr3E[x] = gRandom->Gaus(labr3E[x],(0.02*labr3E[x]/2.355)); 
  	      //}
            
  	      //hEnergyDepositLaBr3[x]->Fill(labr3E[x]);
      
    } 
    
  }  //BC - end LaBr3 array
  
  
  //BC - Hits in Getest
  if(DHCGetest){
  
    int n_hit = DHCGetest->entries();
    
    //G4cout << "n_hit = " << n_hit << endl;
    G4double GetestE = 0;
    
    if(n_hit > 0){
      for(int i=0; i<n_hit; i++){
        G4String name = (*DHCGetest)[i]->GetParticlename();
        G4String parentname = (*DHCGetest)[i]->GetParentName();
        //G4double parentenergy = (*DHCGetest)[i]->GetParentEnergy();
        G4double edep = (*DHCGetest)[i]->GetEdep()/keV;
        G4ThreeVector pos = (*DHCGetest)[i]->GetPostPosition();
        G4int copyno = (*DHCGetest)[i]->GetVolCopyNo();
  
        GetestE += edep;
        
        
      }
      
      double resolution = 0.0015;
      
      if(GetestE > 0){
        
        GetestE = gRandom->Gaus(GetestE,(resolution*GetestE/2.355));
        
      }
      
      //fill histogram
      hGetest->Fill(GetestE);
    }
  
  } //BC - end Getest hits
  
  
  //BC - Hits in EXOtest
  if(DHCEXOtest){
  
    int n_hit = DHCEXOtest->entries();
    
    //G4cout << "n_hit = " << n_hit << endl;
    G4double EXOtestE = 0;
    
    if(n_hit > 0){
      for(int i=0; i<n_hit; i++){
        G4String name = (*DHCEXOtest)[i]->GetParticlename();
        G4String parentname = (*DHCEXOtest)[i]->GetParentName();
        //G4double parentenergy = (*DHCEXOtest)[i]->GetParentEnergy();
        G4double edep = (*DHCEXOtest)[i]->GetEdep()/keV;
        G4ThreeVector pos = (*DHCEXOtest)[i]->GetPostPosition();
        G4int copyno = (*DHCEXOtest)[i]->GetVolCopyNo();
  
        EXOtestE += edep;
        
        
      }
      
      if(EXOtestE > 0){
      
        double resolution = 0.04;  //100 keV at 2.5 MeV
        //double resolution = 0.0004;  //1 keV at 2.5 MeV
        //double resolution = 0.004;   //10 keV at 2.5 MeV
        
        EXOtestE = gRandom->Gaus(EXOtestE,(resolution*EXOtestE/2.355));
        
      }
      
      //fill histogram
      hEXOtest->Fill(EXOtestE);
    }
  
  } //BC - end EXOtest hits


  //BC - Hits in EJ204 Scintillator
  if(DHCEJ204Scint){
    
    int n_hit = DHCEXOtest->entries();
    
    //G4cout << "n_hit = " << n_hit << endl;
    G4double EJ204ScintE = 0;
    
    if(n_hit > 0){
      for(int i=0; i<n_hit; i++){
        G4String name = (*DHCEJ204Scint)[i]->GetParticlename();
        G4String parentname = (*DHCEJ204Scint)[i]->GetParentName();
        //G4double parentenergy = (*DHCEJ204Scint)[i]->GetParentEnergy();
        G4double edep = (*DHCEJ204Scint)[i]->GetEdep()/keV;
        G4ThreeVector pos = (*DHCEJ204Scint)[i]->GetPostPosition();
        G4int copyno = (*DHCEJ204Scint)[i]->GetVolCopyNo();
  
        EJ204ScintE += edep;
        
        
      }
      
      // if(EJ204ScintE > 0){
      
      //   double resolution = 0.04;  //100 keV at 2.5 MeV
      //   //double resolution = 0.0004;  //1 keV at 2.5 MeV
      //   //double resolution = 0.004;   //10 keV at 2.5 MeV
        
      //   EJ204ScintE = gRandom->Gaus(EJ204ScintE,(resolution*EJ204ScintE/2.355));
        
      // }
      
      //fill histogram
      hEJ204Scint->Fill(EJ204ScintE);
    }

  }


  //BC - hits in CeBr3 scintillator
  {
    
    int n_hit = DHCEXOtest->entries();
    
    //G4cout << "n_hit = " << n_hit << endl;
    G4double CeBr3ScintE = 0;
    
    if(n_hit > 0){
      for(int i=0; i<n_hit; i++){
        G4String name = (*DHCCeBr3Scint)[i]->GetParticlename();
        G4String parentname = (*DHCCeBr3Scint)[i]->GetParentName();
        //G4double parentenergy = (*DHCCeBr3Scint)[i]->GetParentEnergy();
        G4double edep = (*DHCCeBr3Scint)[i]->GetEdep()/keV;
        G4ThreeVector pos = (*DHCCeBr3Scint)[i]->GetPostPosition();
        G4int copyno = (*DHCCeBr3Scint)[i]->GetVolCopyNo();
  
        CeBr3ScintE += edep;
        
        
      }
      
      // if(CeBr3ScintE > 0){
      
      //   double resolution = 0.04;  //100 keV at 2.5 MeV
      //   //double resolution = 0.0004;  //1 keV at 2.5 MeV
      //   //double resolution = 0.004;   //10 keV at 2.5 MeV
        
      //   CeBr3ScintE = gRandom->Gaus(CeBr3ScintE,(resolution*CeBr3ScintE/2.355));
        
      // }
      
      //fill histogram
      hCeBr3Scint->Fill(CeBr3ScintE);
    }

  }
  
  

  if(DHCClover){
    
    //cout << "in clover " << endl;

    G4double CloverE[72];
    G4double CloverE_addback[18];
    G4double totE = 0;
    G4double CloverAdd =0;
    G4int crysno;
    for(int x=0;x<72;x++){
      CloverE[x]=0;
    }
    for(int x=0;x<9;x++){
      CloverE_addback[x]=0;
    }
    
    int n_hit = DHCClover->entries();
    double nh = DHCClover->entries();
    clover_hit=clover_hit+nh;
    if((n_hit>0)/*&&(n_hit<25)*/){     /*cout << "found clover event " << n_hit << endl;*/
      // cout <<n_hit<<" crystal numbers: ";
      // cout <<n_hit<<endl;
      for(int i=0; i<n_hit; i++){
	      G4String name = (*DHCClover)[i]->GetParticlename();
	      G4double edep = (*DHCClover)[i]->GetEdep()/keV;
	      G4ThreeVector pos = (*DHCClover)[i]->GetPostPosition();
	      G4int copyno = (*DHCClover)[i]->GetVolCopyNo();
	      G4int parentno = (*DHCClover)[i]->GetParentno();
	      G4int trackno = (*DHCClover)[i]->GetTrackno();
	      //cout<<"name+track"<<name <<trackno<<" parent "<<parentno<<" crys "<<copyno<<" edep "<<edep<<endl;
	      //cout<<name<<" "<<edep<<endl;
	      /*	if(name=="e-" && parentno==1){//makes no difference
	      //cout <<"here"<<endl;
	      }else{
	      //cout<<name<<" "<<edep<<endl;*/
	      CloverE[copyno] += edep;
	      totE += edep;
	      crysno=copyno;
	      // cout<<crysno<<" ";
	      //	}
      }
      
      // add a resolution to the clover energy
      //  int count=0;
      for(int x=0;x<72;x++){
	      if(CloverE[x]>0){
	        //  count++;
	        //resolution of 0.3% across all energy ranges
	        // cout<<x<<",  ";
	        if(CloverE[x]>300){
	          CloverE[x] = gRandom->Gaus(CloverE[x],(0.0037*CloverE[x]/2.35));//.3% is ok based on SeGA 723.3
	        }else if(CloverE[x]>200){
	          CloverE[x] = gRandom->Gaus(CloverE[x],(0.012*CloverE[x]/2.35));//based on SeGA res at 247
	        }else{
	          CloverE[x] = gRandom->Gaus(CloverE[x],(0.022*CloverE[x]/2.35));
	        }
	      }
      }

      for(int x=0;x<72;x++){
	      if(CloverE[x]>0){
	        CloverE_addback[x/4] += CloverE[x]; 
	      }
      }
      //  double max =0; 
      for(int x=0;x<72;x++){
	      //cout << "energy " << x << " " << CloverE[x] << endl;
	      if(CloverE[x]>0){
	        //cout << "filling " << endl;
	        //  if(count<3){
	        // if(CloverE[x]>max) max = CloverE[x];
	        //	  if(x<20){
	        hEnergyDepositClover[x]->Fill(CloverE[x]);//energy of each ind. crystal
	        hClover->Fill(CloverE[x]);//energy from all crystals
	        //	  }
	        //increment only for specific energies so that they may be used as efficiencies
	        /* if((CloverE[x]>3984)&&(CloverE[x]<4015)){
	          hEnergyDepositClover[x]->Fill(CloverE[x]);
	          hClover->Fill(x);
	          }*/
	        //  hEnergyDeposit->Fill(CloverE[x]);//trying to make spectrum with both clover and gedssd events
	        //	hsega_tigress_s->Fill(CloverE[x]);
	        // }
	  
	      }//end if CloverE>0
	      //cout << "end loop " << endl;
      }//end Clover for
      // hClover->Fill(max);
      //  double max2;
      for(int x=0;x<18;x++){
	      if(CloverE_addback[x]>0){
	        // if(CloverE_addback[x]>max2) max2 =CloverE_addback[x];
	        // if(count<3){
	        // if(x<5) {
	        CloverAdd +=CloverE_addback[x];
	        hEnergyDepositClover_addback[x]->Fill(CloverE_addback[x]);
              //energy from 4 crystals added to make a clover
	        hClover_addback->Fill(CloverE_addback[x]);//energy from all clovers
	        //  }
	        //increment only for specific energies so that they may be used as efficiencies
	        /*if((CloverE_addback[x]>3984)&&(CloverE_addback[x]<4015)){
	            hEnergyDepositClover_addback[0]->Fill(x);
	            //hClover_addback->Fill(CloverE_addback[x]);
	        }*/
	        //hsega_tigress->Fill(CloverE_addback[x]);
	        //}
	      }//end Clover_addback>0
      }//end addback for
      // hClover_addback->Fill(max2);
      //if((CloverAdd>3985)&&(CloverAdd<4015)){//increment only for specific energies so that they may be used as efficiencies
      // if(count<3){
      hClover_arrayadd->Fill(CloverAdd);
      //hEnergyDeposit->Fill(CloverAdd);
      // }
      //   hClover_arrayadd->Fill(totE);
      //  if(parentID==1) hPixelGamma->Fill(crysno,crysno/4);
    }//end n_hit>0
    // cout<<" "<<endl;
  }//end Clover
 







  //   //calculate angle between vectors
  //   //check if vectors are unit vectors
  //   if(gamma1modir.mag() != 0 && gamma2modir.mag() != 0){
  //     //G4cout << "unit vector ? " << gamma1modir.mag() << " " << gamma2modir.mag() << G4endl;
  //     //G4cout << "gamma1 " << gamma1modir.x() << " " <<gamma1modir.y() << " " << gamma1modir.z() << G4endl;
  //     //G4cout << "gamma2 " << gamma2modir.x() << " " <<gamma2modir.y() << " " << gamma2modir.z() << G4endl;
  //     G4double dotprod;
  //     dotprod= (gamma1modir/gamma1modir.mag()) * (gamma2modir/gamma2modir.mag());
  //     G4double angle;
  //     angle = acos(dotprod)/CLHEP::degree;
  //     //G4cout << " angle " << angle << G4endl;

  //     hgangcortheory->Fill(angle);
  //   }

  //   //hEnergyDepositElectron2D->Fill(pos.x()/mm + 0.5,pos.y()/mm + 0.5,edep/keV);
  //   //hEnergyDepositGamma2D->Fill(pos.x()/mm + 0.5,pos.y()/mm + 0.5,edep/keV);


  //   // define variables
  //   G4double tedep, tedepe, tedepg, tedepo, tedeplowz, tedephiz;
  //   tedep = tedepe = tedepg = tedepo = tedeplowz = tedephiz = 0;

  //   // strips in x and z
  //   G4int tedepx_size = 20;
  //   for(G4int i=0;i<tedepx_size;i++) tedepx[i] = tedepy[i] = 0;

  //   //fill spectra

  //   // loop over the totaledep array and fill the appropriate spectra
  //   for(G4int type=0; type<3; type++){
 
  //     //for(G4int x=0;x<tedepx_size;x++){
  //     for(G4int x=1;x<17;x++){

  //       // for(G4int y=0;y<tedepx_size;y++){
  //       for(G4int y=1;y<17;y++){
  // 	G4double tempsum = 0;

  // 	for(G4int z=0;z<2;z++){

  // 	  tedepx[x] += totaledep[type][x][y][z];
  // 	  tedepy[y] += totaledep[type][x][y][z];

  // 	  tempsum +=totaledep[type][x][y][z];

  // 	  // accuulate total energy deposited, treat GeDSSD as one crystal calorimeter
  // 	  tedep += totaledep[type][x][y][z];
  // 	  if(z==0) tedeplowz += totaledep[type][x][y][z];
  // 	  else tedephiz += totaledep[type][x][y][z];


  // 	  // accumulate total energy deposited for electron
  // 	  if(type == 0){
  // 	    tedepe += totaledep[type][x][y][z];
  // 	  }
  // 	  // accumulate total energy deposited for gamma
  // 	  if(type == 1){
  // 	    tedepg += totaledep[type][x][y][z];
  // 	  }
  // 	  // accumulate total energy deposited for other
  // 	  if(type == 2){
  // 	    tedepo += totaledep[type][x][y][z];
  // 	  }
  // 	}

  // 	// 	if(tempsum > 0){
  // // 	  G4cout << "pixel E " << tempsum << " " << x << " " << y << " " << G4endl;
  // // 	}
  //       }
  //     }
  //   }

  //   //G4cout << "total E" << tedep << " " << tedepe << " " <<tedepg << " " << tedepo << G4endl;
  //   // if(tedep!=0){
  //   //add a resolution
  //   if(tedep>0){
  //     // tedep = gRandom->Gaus(tedep,(0.03*tedep/2.35));
  //   // }
  //     if(tedep>300) {
  //       tedep = gRandom->Gaus(tedep,(0.003*tedep/2.35));
  //     }else{
  //       tedep = gRandom->Gaus(tedep,(0.03*tedep/2.35));
  //     }
  //     if(tedepg>300) {
  //       tedepg = gRandom->Gaus(tedepg,(0.003*tedepg/2.35));
  //     }else{
  //       tedepg = gRandom->Gaus(tedepg,(0.03*tedepg/2.35));
  //     }
  //     // hEnergyDeposit->Fill(tedep);
  //   //}//there is a large energy deposit at 0 that overwhelems the rest of the spectrum
  //   hEnergyDepositLowz->Fill(tedeplowz);
  //   hEnergyDepositElectron->Fill(tedepe);
  //   if(tedepg>0)  hEnergyDepositGamma->Fill(tedepg);
  //   hEnergyDepositOther->Fill(tedepo);
  //   }


  // //   for(G4int x=0;x<tedepx_size;x++){
  // //     for(G4int y=0;y<tedepx_size;y++){
  // //       //G4double e;
  // //       //e = totaledep[1][x][y][0] + totaledep[1][x][y][1];
  // //       hEnergyDepositPixelGamma[x][y]->Fill(totaledep[1][x][y][0]+totaledep[1][x][y][1]);
  // //     }
  // //   }

  
  //   if(tedepg < 201 && tedepg < 203) goode++;

  //   G4int stripmultx, stripmulty;
  //   stripmultx = stripmulty = 0;

  //   G4double maxEstripx = 0;
  //   G4int maxEstriplocx = 0;
  //   G4double maxEstripy = 0;
  //   G4int maxEstriplocy = 0;

  //   G4bool xval =false;
  //   G4bool yval =false;
  //   for(G4int i=0;i<tedepx_size;i++){

  //     if(tedepx[i] > maxEstripx){
  //       maxEstripx = tedepx[i];
  //       maxEstriplocx = i;
  //     }

  //     if(tedepx[i] > 0){
  //       hEnergyDepositStripx->Fill(tedepx[i]);
  //       stripmultx++;

  //       //perform simple gamma-gamma
  //       for(G4int j=i+1;j<tedepx_size;j++){
  // 	if(tedepx[j] > 0){
  // 	  hGammaGamma->Fill(tedepx[i],tedepx[j]);
  // 	  hGammaGamma->Fill(tedepx[j],tedepx[i]);
  // 	}

  //       }

  //       //G4cout << "stripx " << i << " " << tedepx[i] << G4endl;
  //       if(tedepx[i]>201 && tedepx[i]<203) xval = true;
  //     }
    
  //     if(tedepy[i] > maxEstripy){
  //       maxEstripy = tedepy[i];
  //       maxEstriplocy = i;
  //     }
  //     if(tedepy[i] > 0){
  //       hEnergyDepositStripy->Fill(tedepy[i]);
  //       stripmulty++;
  //       //G4cout << "stripy " << i << " " << tedepy[i] << G4endl;
  //       if(tedepy[i]>201 && tedepy[i]<203) yval = true;
  //     }

  //   }

  //   if((stripmultx>1) && (stripmulty>1))   hEnergyDepositHiz->Fill(tedephiz);
  //   hEnergyDepositLocMax2D->Fill(maxEstriplocx,maxEstriplocy);
  
  //   if(xval || yval) goodestrip++;

  //   //G4cout << " ** end event ** "<< G4endl;
  //   G4int gamma1x, gamma1y;
  //   G4int gamma2x, gamma2y;

  //   gamma1x = gamma1y = gamma2x = gamma2y = 0;

  //   for(G4int x=0;x<tedepx_size;x++){
  //     //find an x strip with some energy
  //     if(tedepx[x] > 0){
  //       //match this x strip to a corresponding y strip within a 1% energy difference
  //       for(G4int y=0;y<tedepx_size;y++){
  // 	if(tedepy[y] > 0){
  // 	  G4double diff = abs(tedepx[x] - tedepy[y]);
  // 	  G4double perdiff = diff / tedepx[x];
  // 	  G4double avg = (tedepx[x] + tedepy[y])/2.;
  // 	  if(perdiff < 0.01){
  // 	    //match found zero values and put into 2D array
  // 	    //G4cout << "energy in DSSD " << x << " " << y << " " << avg/keV <<G4endl;
  // 	    totaledep_xy[x][y] = avg;
  // 	    hEnergyDepositPixel->Fill(avg);
  // 	    if(x!=9 && y!=9) hEnergyDepositGamma_nocenter->Fill(avg);

  // 	    if(avg>201 && avg < 203){
  // 	      gamma1x = x; gamma1y = y;
  // 	    }
  // 	    if(avg > 197 && avg < 199){
  // 	      gamma2x = x; gamma2y = y;
  // 	    }
  // 	  }
  // 	}
  //       }
  //     }
  //   }

  //   if(gamma1x != 0 && gamma2x != 0){
  //     G4ThreeVector gamma1(gamma1x-9,gamma1y-9,0);
  //     G4ThreeVector gamma2(gamma2x-9,gamma2y-9,0);
  //     G4double dotprod;
  //     dotprod= (gamma1/gamma1.mag()) * (gamma2/gamma2.mag());
  //     G4double angle;
  //     angle = acos(dotprod)/CLHEP::degree;
  //     //G4cout << "dotprod " << dotprod << " angle " << angle << G4endl;

  //     hgangcorexper->Fill(angle);
  //   }

  //hStripMultx->Fill(stripmultx);


  if (e17011_simtree){

    e17011_simtree->Fill();
  }



#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::ClassifyNewTrack(
						const G4Track *aTrack,
						G4ClassificationOfNewTrack */*classification_ptr*/)
{
  // cout<<"ClassifyNewTrack"<<endl;
  // G4ClassificationOfNewTrack &classification = (*classification_ptr);
  
#if defined (G4ANALYSIS_USE_ROOT)
  
  if (gSystem) gSystem->ProcessEvents();
  
  G4String particleName = aTrack->GetDefinition()->GetParticleName();
  
  if (particleName == "opticalphoton")
    {
      Double_t aWaveLength = 0.0; // will be in [nanometer]
      aWaveLength = (LambdaE / aTrack->GetTotalEnergy()) / nanometer; // in [nanometer]
      // aWaveLength = (LambdaE / aTrack->GetKineticEnergy()) / nanometer; // in [nanometer]
      if (hOPWaveLength) hOPWaveLength->Fill(aWaveLength, aTrack->GetWeight()); // Wavelength of the produced optical photon
    }
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::NewStage()
{
  // cout<<"e17011_simAnalysisExample NewStage"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  
  if (gSystem) gSystem->ProcessEvents();
  
  
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::PrepareNewEvent()
{
  //  cout<<"e17011_simAnalysisExample PreparNewEvent"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  
  if (gSystem) gSystem->ProcessEvents();
  
  
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

void e17011_simAnalysisExample::PrepareNewRun(G4double x, G4double y, G4double z, G4double zpos,
					     G4int r, G4int c)
{
  // cout<<"e17011_simAnalysisExample PrepareNewRun"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  
  if (gSystem) gSystem->ProcessEvents();
  
  //set up detector parameters

  //position of the scint fiber, the scint fiber length
  //   scint_fiber_x = x;
  //   scint_fiber_y = y;
  //   scint_fiber_z = z;
  //   scint_fiber_zpos = zpos;
  //   scint_rows = r;
  //   scint_columns = c;
  
  //set up the dimensions of the scint energy deposition array
  // intialize all elements to zero
  //   for(G4int i=0; i<scint_rows; i++){
  //     for(G4int j=0; j<scint_columns; j++){
  //       for(G4int k=0; k<int((scint_fiber_z/mm)/scint_zpixel)+1; k++){
  // 	temp1d.push_back(0);
  //       }
  //       temp2d.push_back(temp1d);
  //       temp2dzero.push_back(temp1d);
  //       temp1d.clear();
  //     }
  //     scint_fiber_edep.push_back(temp2d);
  //     scint_fiber_edep_zero.push_back(temp2d);
  //     temp2d.clear();
  //   }


#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::PreUserTrackingAction(const G4Track */*aTrack*/)
{
  //   cout<<"e17011_simAnalysisExample PreUserTrackingAction"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  if (gSystem) gSystem->ProcessEvents();
  
  
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::PostUserTrackingAction(const G4Track */*aTrack*/,
						      G4TrackStatus */*status_ptr*/)
{
  // cout<<"e17011_simAnalysisExample PostUserTrackingAction"<<endl;
  // G4TrackStatus &status = (*status_ptr);
  
#if defined (G4ANALYSIS_USE_ROOT)
  if (gSystem) gSystem->ProcessEvents();
  
  
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::UserSteppingAction(const G4Step *aStep)
{
  //  cout<<"e17011_simAnalysisExample UserSteppingAction"<<endl;
#if defined (G4ANALYSIS_USE_ROOT)
  

  //Double_t aStepLength = (aStep->GetStepLength() / cm); // in [cm]
  // Double_t aTotalEnergyLoss = (aStep->GetTotalEnergyDeposit() / MeV); // in [MeV]
  //Double_t aTotalEnergyLoss = -(aStep->GetDeltaEnergy() / MeV); // in [MeV]
  //Double_t aTotalEnergyLossPerLengthUnit = -1.0; // will be in [MeV/cm]
  //Double_t aWeight = aStep->GetTrack()->GetWeight();
  
  //if (hStepLength) hStepLength->Fill(aStepLength, aWeight); // Step Length
  
  //if (hStepTotELoss) hStepTotELoss->Fill(aTotalEnergyLoss, aWeight); // Step Total Energy Loss
  
  //if (aStepLength != 0.0) aTotalEnergyLossPerLengthUnit = aTotalEnergyLoss / aStepLength;
  //else aTotalEnergyLossPerLengthUnit = aTotalEnergyLoss / 1.0e-7; // in [MeV/cm]
  //if (hTotELossNorm) hTotELossNorm->Fill(aTotalEnergyLossPerLengthUnit, aWeight); // Total Energy Loss per Step Length Unit
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
  
  OnceAWhileDoIt();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void e17011_simAnalysisExample::OnceAWhileDoIt(const G4bool DoItNow)
{
  //  cout<<"e17011_simAnalysisExample OnceAWhileDoIt"<<endl;
  time_t Now = time(0); // get the current time (measured in seconds)
  if ( (!DoItNow) && (LastDoItTime > (Now - 10)) ) return; // every 10 seconds
  LastDoItTime = Now;
  
#if defined (G4ANALYSIS_USE_ROOT)
  
  TVirtualPad *CurrentPad = gPad;
  
  const char OnceAWhileCanvas[] = "e17011_simAnalysisExampleCanvas";
  
  TCanvas *c = (TCanvas *)gROOT->GetListOfCanvases()->FindObject(OnceAWhileCanvas);
  
  if ( DoItNow && (!c) && ((G4RunManager::GetRunManager())->GetVerboseLevel() > 1) )
    {
      // first, try to make sure we do not run in batch mode
      if (!gApplication)
	new TApplication("Application", ((int *)0), ((char **)0));
      // then, create the TCanvas
      c = new TCanvas(OnceAWhileCanvas, OnceAWhileCanvas);
      if (c) {
	c->Divide(2,2);
	c->cd(1); if (hStepLength) hStepLength->Draw();
	c->cd(2); if (hStepTotELoss) hStepTotELoss->Draw();
	c->cd(3); if (hOPWaveLength) hOPWaveLength->Draw();
	c->cd(4); if (hTotELossNorm) hTotELossNorm->Draw();
	c->cd(0);
      }
    }
  
  if (c) {
    //c->cd(0);
    c->Draw();
    c->Update();
  }
  
  if (CurrentPad) CurrentPad->cd();
  
  if (gSystem) gSystem->ProcessEvents();
  
#endif /* defined (G4ANALYSIS_USE_ROOT) */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// end of file e17011_simAnalysisExample.cc by Jacek M. Holeczek

void e17011_simAnalysisExample::Setfilename(G4String nam)
{
  filename = nam;
}

void e17011_simAnalysisExample::Setdirname(G4String nam)
{
  dirname = nam;
}

//void e17011_simAnalysisExample::OpticalPhotons(G4int status, G4double val)
//{
  //   if(status == 1){
  //     opticalphoton++;
  //   }
  //   else if(status == 9){
  //     opticalkill++;
  //   }
  //   else{
  //     G4cout << "wandering gamma" << G4endl;
  //   }

//}

//void e17011_simAnalysisExample::PSF(G4int status, G4double val1, G4double val2)
//{
  
  //   if(status == 0){
  //     e17011_simap_len = e17011_simap_len + val1;
  //   }
  //   else if(status == 1){
  //     psfx = val1;
  //     psfy = val2;
  //     htally->Fill(psfx,psfy);
  //     hpsf->Fill(psfx,psfy,exp(-e17011_simap_len/5));
  //   }

//}
