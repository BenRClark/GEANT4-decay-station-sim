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

#ifndef e17011_simAnalysisExample_h
#define e17011_simAnalysisExample_h 1

//#include "e17011_simAnalysisManager.hh"
#include "e17011_simAnalysisManager.hh"


class G4VPhysicalVolume;
class G4Event;
class G4Run;
class G4Track;
class G4Step;
class e17011_simAnalysisExampleMessenger;
//class ExCPAnalysisExampleMessenger;

class e17011_simDetectorConstruction;
//class ExCPDetectorConstruction;

#include "G4ClassificationOfNewTrack.hh"
#include "G4ThreeVector.hh"
#include "G4TrackStatus.hh"
#include "G4Types.hh"
#include <time.h>
//#include <vector>
#include <iostream>
#include <cctype>
using namespace std;


#if defined (G4UI_USE_ROOT) || defined (G4UI_BUILD_ROOT_SESSION)
#if !defined (G4ANALYSIS_USE_ROOT)
#define G4ANALYSIS_USE_ROOT 1
#endif /* !defined (G4ANALYSIS_USE_ROOT) */
#endif /* defined (G4UI_USE_ROOT) || defined (G4UI_BUILD_ROOT_SESSION) */

#if defined (G4ANALYSIS_USE_ROOT)

class TH1D;
class TH2D;
class TFile;
class TNtuple;
class TStopwatch;
class TTree;

//   struct declare_event //event information
//   {
//     double parent_alpha_E;
//     double alpha1_E;
//   };
//   extern declare_event user_event;

//G4int ex_check;

#endif /* defined (G4ANALYSIS_USE_ROOT) */


//class e17011_simAnalysisExample : public e17011_simAnalysisManager {
class e17011_simAnalysisExample : public e17011_simAnalysisManager {
  
public:
  e17011_simAnalysisExample();
  ~e17011_simAnalysisExample();
  
public:
  // G4VUserDetectorConstruction
  void Construct(const G4VPhysicalVolume*);
  
  // G4VUserPhysicsList
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();
  
  // G4VUserPrimaryGeneratorAction
  void GeneratePrimaries(const G4Event*, const G4double beam);
  // void ResetStartLocation(G4double &, G4double &);

  void DetectorInfo(e17011_simDetectorConstruction*);
  // void DetectorInfo(ExCPDetectorConstruction*);
  
  // G4UserRunAction
  void BeginOfRunAction(const G4Run*, G4int);
  void EndOfRunAction(const G4Run*, G4int);
  
  // G4UserEventAction
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
  // G4UserStackingAction
  void ClassifyNewTrack(const G4Track*, G4ClassificationOfNewTrack*);
  void NewStage();
  void PrepareNewEvent();
  void PrepareNewRun(G4double, G4double, G4double, G4double, G4int, G4int);
  
  // G4UserTrackingAction
  void PreUserTrackingAction(const G4Track*);
  void PostUserTrackingAction(const G4Track*, G4TrackStatus*);
  
  // G4UserSteppingAction
  void UserSteppingAction(const G4Step*);
  
  // once a while do "something"
  void OnceAWhileDoIt(const G4bool DoItNow = false);

  //void SetparentZ(G4int val) {parentZ = val;};
  //void SetdaughterZ(G4int val) {daughterZ = val;};
  void Setfilename(G4String nam);
  void Setdirname(G4String nam);

  // void OpticalPhotons(G4int, G4double);
  // void PSF(G4int, G4double, G4double);

private:

  time_t LastDoItTime; // used in OnceAWhileDoIt method
  
  G4int UseGeThickDetector;
  G4int UseGeThinDetector;

  G4double GeThickDetectorThickness;
  G4double GeThinDetectorThickness;

  G4double GeThickDetectorRadius;
  G4double GeThinDetectorRadius;

  G4double GeThickx;
  G4double GeThicky;
  G4double GeThickz;
  G4double GeThinx;
  G4double GeThiny;
  G4double GeThinz;

  G4double GeSpacing;

  G4int pixelhits;
  G4int gmult;
  G4double stripwidthx;
  G4double stripwidthy;
  G4double stripwidthz;

  G4double daughter_energy, e17011_sim_energy;

  G4int filecount;
  G4double totE;
  G4double edep;
  G4int ID;
  G4double Z;
  G4String particletype;
  G4String particleName;
  G4String process;
  G4int creatorID;
  G4int StepNo;
  G4double gtime;


  

  G4double avgimpx;
  G4double avgimpy;
  
  G4double countf;
  G4double countb;

  G4double max_energy_x;
  G4double max_energy_y;
  G4double max_strip_x;
  G4double max_strip_y;
  
    G4String volname;
    G4int volcopyno, e17011_simcopyno;

    G4double e17011_sim_edep, crystal_edep[68];

  G4int decaycount;
  G4int decaycount2;
  G4int gamma189count;
  G4int gammacount;
  G4int electroncount;
  G4int Mult1both;
  G4int Mult1front;
  G4int Mult1back;
  G4int Mult2both;
  G4int Mult2front;
  G4int Mult2back;
  G4int Multgreaterthan2;
  G4int Mult1bothTotal;
  G4int Mult1frontTotal;
  G4int Mult1backTotal;
  G4int Mult2bothTotal;
  G4int Mult2frontTotal;
  G4int Mult2backTotal;
  G4int Multgreaterthan2Total;
  G4int Mult1bothNoLeft;
  G4int Mult1frontNoLeft;
  G4int Mult1backNoLeft;
  G4int Mult2bothNoLeft;
  G4int Mult2frontNoLeft;
  G4int Mult2backNoLeft;
  G4int Multgreaterthan2NoLeft;
  G4int Mult1both189;
  G4int Mult1front189;
  G4int Mult1back189;
  G4int Mult2both189;
  G4int Mult2front189;
  G4int Mult2back189;
  G4int Multgreaterthan2189;
  G4int Mult1bothFound189;
  G4int Mult1frontFound189;
  G4int Mult1backFound189;
  G4int Mult2bothFound189;
  G4int Mult2frontFound189;
  G4int Mult2backFound189;
  G4int Multgreaterthan2Found189;
  G4int AlgorithmMiss;
  G4double clover_hit;
  G4double implantcorr;
  G4double decay_in_same_imp_strip;
  G4double decay_in_1_strip_away;
  G4double decay_1_away_xy;
 
  G4double avgdecx;
  G4double avgdecy;
  G4double decayelectrons;
  G4double avgoffbygreaterthan1x;
  G4double avgoffbygreaterthan1y;
  G4double offcount;
  G4int number0;
  G4int numberother;
  G4int twobytwoB;
  G4int twobytwoA;
  G4int twobytwoBcorrect;
  G4int twobytwoAcorrect;
  G4int twobythree;
  G4int twobythreecorrect;

  G4ThreeVector prepos;
  G4ThreeVector postpos;
  G4ThreeVector deltapos;

  G4double beam_energy;

  G4double GeDSSDcal;

  G4double tedepx[40];
  G4double tedepy[40];
  //for the full scint array


  G4double pathlength;
  G4double totaledep[3][100][100][2];
  G4double totalelectronedep;
  G4double totalgammaedep;
  G4double e17011_simap_len;
  G4double psfx;
  G4double psfy;

  G4double goode;
  G4double goodestrip;

  G4int eventcounter;
  G4int particletypebit;
  G4int empty;
  e17011_simAnalysisExampleMessenger* analysisexample;

  //char ch[] = "";
  //const char *cp = ch;
  //char filename[100];
  //const char *filepointer;

  G4String filename;
  G4String dirname;
  //C4Char *filepointer = filename;

#if defined (G4ANALYSIS_USE_ROOT)
  
  TFile *file1; // test file
  TFile *rfile1; //retrieve file
  TNtuple *e17011_simntup; // e17011_sim ntuple output
  TTree *e17011_simtree;
  TTree *e17011_siminfo;

  TStopwatch *timer;
  TStopwatch *eventtimer;

  TH1D *hEnergyDeposit;
  TH1D *hPixelAndAlgorithmDifference;
  TH1D *hEnergyDepositmult1both;
  TH1D *hEnergyDepositmult1front;
  TH1D *hEnergyDepositmult1back;
  TH1D *hEnergyDepositmult2both;
  TH1D *hEnergyDepositmult2front;
  TH1D *hEnergyDepositmult2back;
  TH1D *hEnergyDepositmultgreaterthan2;
  TH1D *hEnergyDepositcal;
  TH1D *hEnergyDepositTest;
  TH1D *hEnergyDepositLeftOverx;
  TH1D *hEnergyDepositLeftOvery;
  TH1D *hEnergyDepositAlgorithmTest;
  TH1D *hEnergyDepositAllPlusLeftovers;
  TH1D *hPercentMatch;
  TH1D *hElectronMult;
  TH1D *hTotalMult;
  TH1D *hGammaMult;
  TH1D *hNumberof188;
  TH1D *hEnergyDepositx;
  TH1D *hEnergyDeposity;
  TH1D *hEnergyDepositMaxStrip;
  TH1D *hEnergyDepositAllStrips;
  TH1D *hEnergyDepositLowz;
  TH1D *hEnergyDepositHiz;
  TH1D *hEnergyDepositPixel;
  TH1D *hEnergyDepositElectron;
  TH2D *hEnergyDepositElectron2D;
  TH2D *hTracklengthvsE;
  TH2D *hMultvsE;
  TH2D *hMultvsEcountf2;
  TH1D *hEnergyDepositGamma;
  TH1D *hEnergyDepositGamma_nocenter;
  TH1D *hEnergyImplantPixelMult1BothSides;
  TH1D *hEnergyImplantPixel;
  TH2D *hEnergyDepositGamma2D;
  TH2D *hEnergyDepositLocMax2D;
  TH2D *hEvsStripx;
  TH2D *hEvsStripy;
  TH2D *hFinalPosition;
  TH2D *hEnergyvsMulf;
  TH2D *hEnergyvsMulb;
  TH2D *hPositionOfGamma;
  TH1D *hEdepMultFront[40];
  TH1D *hEdepMultBack[40];
  //TH1D *hEnergyDepositPixelGamma[100][100];
  TH1D *hEnergyDepositOther;
  TH1D *helectronx;
  TH1D *helectrony;
  TH1D *hnotelectronx;
  TH1D *hnotelectrony;
  TH1D *hmaxtype;
  TH1D *hTrackLengLowE;
  TH1D *hTrackLengHighE;
  TH1D *hTrackLengHighestE;
  TH1D *hBackMultforFrontMult1;
  TH1D *hBackMultforFrontMult2;
  TH1D *hHighETotalE;
  TH1D *hHighestETotalE;
  TH1D *hLowETotalE;
  TH1D *hHighEFractionInImplant;
  TH1D *hLowEFractionInImplant;
  TH1D *hHighestEFractionInImplant;
  TH1D *hChangeinZlow;
  TH1D *hChangeinZhigh;

  TH2D *hGammaGamma;
  TH2D *hmultvsmult;
  TH2D *haverageEvsmult;

  TH1D *hEnergyDepositStripx;
  TH1D *hEnergyDepositStripy;

  TH1D *hStripMaxEnergyx;
  TH1D *hStripMaxEnergyy;
 
  TH1D *hStripMultx;
  TH1D *hStripMulty;
  TH1D *hMult1x;
  TH1D *hMultmanyx;

  TH1D *hBetaDecay;

  // TH1D *hStripMultx;

  TH1D *hStepLength; // Step Length
  TH1D *hStepTotELoss; // Step Total Energy Loss
  TH1D *hTotELossNorm; // Total Energy Loss per Step Length Unit
  
  TH1D *hOPWaveLength; // Wavelength of the produced optical photon
  
  TH1D *hAlpha1_energy; // alpha 1 energy deposited in detector
  TH1D *hAlpha2_energy; // alpha 2 energy deposited in detector
  TH1D *hParent_energy; // energy deposited from parent ion (including implantation)
  TH1D *hDaughter_energy; // energy deposited from daughter ion

  TH1D *hTot_energy; // energy deposited in DSSD
  TH1D *hBeam_energy; // energy deposited in DSSD
  TH2D *hPrimaryPosin; //primary map input from external source
  TH2D *hPrimaryPos; //primary position of particle

  TH1D *hCrystalEdep[68]; //energy deposited in each crystal
  TH1D *hGammaSingle;
  TH1D *hGammaSingleAdd;
  TH1D *hParent_tracklen; //
  
  TH2D *hScintEdep_xy;
  TH2D *hScintEdep_yz;
  
  TH2D *hScintPhoton_xy;
  TH2D *htally;
  TH2D *hpsf;

  TH1D *hgangcortheory;
  TH1D *hgangcorexper;

  TH1D *hEnergyDepositSega[16];
  TH1D *hSega;
  TH1D *hSegaab;  //Sega addback
  TH1D *hEnergyDepositClover[72];
  TH1D *hClover;
  TH1D *hEnergyDepositClover_addback[18];
  TH1D *hClover_addback;
  TH1D *hClover_arrayadd;
  TH2D *hClover_gg;
  TH2D *hClover_gg_esum;
  TH1D *hsega_tigress;
  TH1D *hsega_tigress_s;
  TH1D *hGeDSSD;
  TH1D *hGeDSSDsecondmax;
  TH1D *hEventsSeparateFromElectron;
  TH1D *hEnergyDepositGeDSSDx[40];
  TH1D *hEnergyDepositGeDSSDy[40];

  TH1D *hEnergyDepositGeDSSDx_gg[40];
  TH1D *hEnergyDepositGeDSSDy_gg[40];
  //  TH1D *hEnergyDepositGeDSSDy_gg_ag[40];

  TH1D *hEnergyDepositNoAlgorithmGeDSSDx[40];
  TH1D *hEnergyDepositNoAlgorithmGeDSSDy[40];
  TH1D *hEnergyDepositGeDSSDxy[40][40];
  TH2D *hPixelGamma;
  TH2D *hPixele;
  TH2D *hPixelFe;
  TH1D *hgammaelectrontime;
  TH1D *hBetaDist;
  TH1D *hBetaDistAl;

  TH1D *hEnergyDepositcal_gammagated;
  TH1D *hEnergyDepositcal_gammagated_antigated;


  //TH1D *hGeDSSDImaginaryStrips;
  
  
  //BC - add histograms for LaBr3 array
  TH1D *hEnergyDepositLaBr3[22];
  TH1D *hLaBr3;
  TH2D *hLaBr3gg;
  int numLaBr3;
  
  //BC - add histograms for Getest
  TH1D *hEnergyDepositGetest;
  TH1D *hGetest;
  
  //BC - add histograms for EXOtest
  TH1D *hEnergyDepositEXOtest;
  TH1D *hEXOtest;
  
  //BC - add histograms for EJ204 Scintillator
  TH1D *hEnergyDepositEJ204Scint[65];
  TH1D *hEJ204Scint;
  int numScint;

  //BC - add histograms for CeBR3 Scintillator
  TH1D *hEnergyDepositCeBr3Scint[257];
  TH1D *hCeBr3Scint;
  int numScintCeBr3;

#endif /* defined (G4ANALYSIS_USE_ROOT) */
  



};

#endif
