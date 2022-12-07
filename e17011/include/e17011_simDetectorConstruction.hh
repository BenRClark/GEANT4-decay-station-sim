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
#ifndef e17011_simDetectorConstruction_h
#define e17011_simDetectorConstruction_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
//#include "e17011_simMagneticField.hh"
#include "DetectionSystemTigress.hh"

class G4Tubs;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Region;
class e17011_simDetectorMessenger;
class e17011_simMaterial;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class e17011_simDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     e17011_simDetectorConstruction();
    ~e17011_simDetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
     const 
     G4VPhysicalVolume* GetDetector() {return physiGeThickDetector;};

  //   G4double GetDetectoFullLength() {return fDetectorLength;};
  //G4double GetTargetFullLength()  {return fTargetLength;};

  // get world paramters
     G4double GetWorldFullLength()   {return fWorldLength;}; 
     G4double GetWorldRadius()   {return fWorldRadius;}; 

  // get thick ge detector parameters
     G4double GetGeThickDetectorThickness() {return fGeThickDetectorThickness;};
     G4double GetGeThickDetectorRadius()  {return fGeThickDetectorRadius;};

  // get thin ge detector parameters
     G4double GetGeThinDetectorThickness() {return fGeThinDetectorThickness;};
     G4double GetGeThinDetectorRadius()  {return fGeThinDetectorRadius;};
     
  // get the boolean values to determine which detector(s) is being used
     G4bool GetUseGeThickDetector() {return fUseGeThickDetector;};
     G4bool GetUseGeThinDetector() {return fUseGeThinDetector;};

  // get the inter detector spacing
     G4double GetGeSpacing() {return fSpacing;};

  // get the position of the thick and thin detector
     G4ThreeVector GetGeThickDetectorPos() {return fGeThickDetectorPos;};
     G4ThreeVector GetGeThinDetectorPos() {return fGeThinDetectorPos;};

  // set the detector material for the thick detector
     void setGeThickDetectorMaterial(G4String);
  // set the detector material for the thin detector
     void setGeThinDetectorMaterial(G4String);

  // set the thick ge detector parameters
     void setGeThickDetectorRadius (G4double value) { fGeThickDetectorRadius = value; };
     void setGeThickDetectorThickness(G4double value){ fGeThickDetectorThickness = value;};  

  // set the thin ge detector parameters
     void setGeThinDetectorRadius (G4double value) { fGeThinDetectorRadius = value; };
     void setGeThinDetectorThickness(G4double value){ fGeThinDetectorThickness = value;};  

  // set the inter detector spacing
     void setDetectorSpacing( G4double value) {fSpacing = value;};

  // set which Ge detector(s) will be used
     void setUseGeThickDetector( G4bool value) {fUseGeThickDetector = value;};
     void setUseGeThinDetector( G4bool value) {fUseGeThinDetector = value;};

  // set which aspects of the planar Ge crystat will be used
     void setUseGeThickDetectorCryo( G4bool value) {fUseGeThickDetectorCryo = value;};
     void setUseGeThickDetectorCryoEndCap( G4bool value) {fUseGeThickDetectorCryoEndCap = value;};

  // set which external Ge array will be used
     void setUseSega( G4bool value) {fUseSega = value;};
     void setUseClover( G4bool value) {fUseClover = value;};
  void setUseCloverBracket( G4bool value) {fUseCloverBracket = value;};

  // BC - set whether or not the LaBr3 array will be used
     void setUseLaBr3( G4bool value) {fUseLaBr3 = value;}; 
     
  // BC - set whether or not LaBr3 Frame will be used
     void setUseLaBr3Frame( G4bool value) {fUseLaBr3Frame = value;};
     
  // BC - set whether or not 3Hen will be used
     void setUse3Hen( G4bool value) {fUse3Hen = value;};
     
  // BC - set whether or not Moving Tape Collector will be used
     void setUseMTC( G4bool value) {fUseMTC = value;};
     
  // BC - set whether or not LENDA will be used
     void setUseLENDA( G4bool value) {fUseLENDA = value;};
     
  // BC - set whether or not Ge test will be used
     void setUseGetest( G4bool value) {fUseGetest = value;};  
     
  // BC - set whether or not EXO test will be used
     void setUseEXOtest( G4bool value) {fUseEXOtest = value;};
     
  // BC - set whether or not EJ204 scintillator will be used (e14057)
     void setUseEJ204Scint( G4bool value) {fUseEJ204Scint = value;};
     
  // BC - set whether or not PSPMT will be used (e14057 + e16032)
     void setUsePSPMT( G4bool value) {fUsePSPMT = value;};
     
  // BC - set whether or not SiDSSD will be used (e14057)
     void setUseSiDSSD( G4bool value) {fUseSiDSSD = value;};

  // BC - set whether or not beam pipe will be used (e14057 + e16032)
     void setUsePipe( G4bool value) {fUsePipe = value;};

  //TO - face plate
     void setUsePipeFacePlate( G4bool value) {fUsePipeFacePlate = value;};
  
  // BC - set whether or not CeBr3 scintillator will be used (e16032)
     void setUseCeBr3Scint( G4bool value) {fUseCeBr3Scint = value;};
  

  //void setDetectorLength(G4double value){ fDetectorLength = value;};  
     
  private:
     void DefineMaterials();
     
  private: 

     e17011_simDetectorMessenger* detectorMessenger;  // pointer to the Messenger
     e17011_simMaterial* materialsManager;            // material manager

     G4Tubs*             solidWorld;              // pointer to the solid envelope 
     G4LogicalVolume*   logicWorld;               // pointer to the logical envelope
     G4VPhysicalVolume* physiWorld;               // pointer to the physical envelope
                  
     G4Tubs*             solidGeThickDetector;    // pointer to the solid Detector
     G4LogicalVolume*   logicGeThickDetector;     // pointer to the logical Detector
     G4VPhysicalVolume* physiGeThickDetector;     // pointer to the physical Detector

     G4Tubs*            solidGeThinDetector;      // pointer to the solid Detector
     G4LogicalVolume*   logicGeThinDetector;      // pointer to the logical Detector
     G4VPhysicalVolume* physiGeThinDetector;      // pointer to the physical Detector
     
     G4Tubs*            solidGeThickDetectorCryo; // pointer to the solid Detector cryostat
     G4LogicalVolume*   logicGeThickDetectorCryo; // pointer to the logical Detector cryostat
     G4VPhysicalVolume* physiGeThickDetectorCryo; // pointer to the physical Detector cryostat

     G4Tubs*            solidGeCryoEndCap;        // pointer to the solid Detector end cap
     G4LogicalVolume*   logicGeCryoEndCap;        // pointer to the logical Detector end cap
     G4VPhysicalVolume* physiGeCryoEndCap;        // pointer to the physical Detector end cap

     G4Tubs*            solidSegaCryo;            // pointer to the solid Sega cryostat
     G4LogicalVolume*   logicSegaCryo;            // pointer to the logical Sega cryostat
     G4VPhysicalVolume* physiSegaCryo;            // pointer to the physical Sega cryostat

     G4Tubs*            solidSegaVac2;            // pointer to the solid Sega vacuum inside the cryo
     G4LogicalVolume*   logicSegaVac2;            // pointer to the logical Sega vacuum inside the cryo
     G4VPhysicalVolume* physiSegaVac2;            // pointer to the physical Sega vacuum inside the cryo
 
     G4Tubs*            solidSegaDetCup;          // pointer to the solid Sega detector cup
     G4LogicalVolume*   logicSegaDetCup;          // pointer to the logical Sega detector cup
     G4VPhysicalVolume* physiSegaDetCup;          // pointer to the physical Sega detector cup

     G4Tubs*            solidSegaVac1;            // pointer to the solid Sega vacuum inside the det cup
     G4LogicalVolume*   logicSegaVac1;            // pointer to the logical Sega vacuum inside the det cup
     G4VPhysicalVolume* physiSegaVac1;            // pointer to the physical Sega vacuum inside the det cup

     G4Tubs*            solidSegaCrystal;         // pointer to the solid Sega crystal
     G4LogicalVolume*   logicSegaCrystal;         // pointer to the logical Sega crystal
     G4VPhysicalVolume* physiSegaCrystal;         // pointer to the physical Sega crystal

     G4Tubs*            solidSegaActiveArea;      // pointer to the solid Sega active area
     G4LogicalVolume*   logicSegaActiveArea;      // pointer to the logical Sega active area
     G4VPhysicalVolume* physiSegaActiveArea;      // pointer to the physical Sega active area

     G4Tubs*            solidSegaDLinner;         // pointer to the solid inner dead layer of SeGA
     G4LogicalVolume*   logicSegaDLinner;         // pointer to the logical Sega crystal
     G4VPhysicalVolume* physiSegaDLinner;         // pointer to the physical Sega crystal

     G4Tubs*            solidSegaCentralCon;      // pointer to the solid Sega central contact
     G4LogicalVolume*   logicSegaCentralCon;      // pointer to the logical Sega central contact
     G4VPhysicalVolume* physiSegaCentralCon;      // pointer to the physical Sega central contact
      
//      G4Box*             solidTIGRESSCrystal;      // pointer to the solid TIGRESS Crystal     
//      G4LogicalVolume*   logicTIGRESSCrystal;      // pointer to the logical TIGRESS Crystal     
//      G4VPhysicalVolume* physiTIGRESSCrystal;      // pointer to the physical TIGRESS Crystal     

//      G4Box*             solidTIGRESSVacuum;       // pointer to the solid TIGRESS Vacuum     
//      G4LogicalVolume*   logicTIGRESSVacuum;       // pointer to the logical TIGRESS Vacuum     
//      G4VPhysicalVolume* physiTIGRESSVacuum;       // pointer to the physical TIGRESS Vacuum     

//      G4Box*             solidTIGRESSCryo;         // pointer to the solid TIGRESS Cryo
//      G4LogicalVolume*   logicTIGRESSCryo;         // pointer to the logical TIGRESS Cryo     
//      G4VPhysicalVolume* physiTIGRESSCryo;         // pointer to the physical TIGRESS Cryo     

     G4Tubs* solidCloverBool1;                    // pointer to the solid boolean clover
     G4Box* solidCloverBool2;                     // pointer to the solid boolean clover

     G4Tubs* solidCloverBool3;                    // pointer to the solid boolean clover
     G4Box* solidCloverBool4; 

  //     G4IntersectionSolid* solidCloverCrystal;     // pointer to the solid clover crystal
     G4Tubs* solidCloverCrystal;                  // pointer to the solid clover crystal
  G4Tubs* solidCloverCrystal2;
  G4Tubs* solidCloverCrystalDeadLayer;
  G4Tubs* solidCloverCrystalDeadLayerWasher;
    G4Tubs* solidCloverCrystalDeadLayerWasher2;
     G4LogicalVolume* logicCloverCrystal;         // pointer to the logical clover crystal
     G4LogicalVolume* logicCloverCrystal2;
  G4LogicalVolume* logicCloverCrystalDeadLayer;
  G4LogicalVolume* logicCloverCrystalDeadLayerWasher;
    G4LogicalVolume* logicCloverCrystalDeadLayerWasher2;
     G4VPhysicalVolume* physiCloverCrystal;       // pointer to the physical clover crystal
     G4VPhysicalVolume* physiCloverCrystal2;
  G4VPhysicalVolume* physiCloverCrystalDead;
  G4VPhysicalVolume* physiCloverCrystalDead2;
  G4VPhysicalVolume* physiCloverCrystalDeadWasher;
  G4VPhysicalVolume* physiCloverCrystalDeadWasher2;

     G4Tubs* solidCloverCore;                     // pointer to the solid clover core
     G4Tubs* solidCloverCore2;
     G4LogicalVolume* logicCloverCore;            // pointer to the logical clover core
     G4LogicalVolume* logicCloverCore2; 
     G4VPhysicalVolume* physiCloverCore;          // pointer to the physical clover core
     G4VPhysicalVolume* physiCloverCore2; 

     G4Box* solidCloverCryo;                      // pointer to the solid clover cryostat
      G4Box* solidCloverCryo2;  
     G4LogicalVolume* logicCloverCryo;            // pointer to the logical clover cryostat
     G4VPhysicalVolume* physiCloverCryo;          // pointer to the physical clover cryostat
  G4LogicalVolume* logicCloverCryo2;

     G4Box *solidCloverCryoVac;                   // pointer to the solid clover cryostat vacuum
     G4Box *solidCloverCryoVac2;
     G4LogicalVolume *logicCloverCryoVac;         // pointer to the logical cryostat vacuum
     G4VPhysicalVolume* physiCloverCryoVac;       // pointer to the physical cryostat vacuum
     G4VPhysicalVolume* physiCloverCryoVac2; 
     G4LogicalVolume *logicCloverCryoVac2; 

  G4Box *outbracket;             //pointer to the outerbracket (for use in G4SubtractionSolid)
  G4Box *innerbracket;           // pointer to inner bracket (for use in G4SubtractionSolid)
     G4SubtractionSolid *solidCloverBracket;                   // pointer to the solid clover cryostat vacuum
     G4LogicalVolume *logicCloverBracket;         // pointer to the logical cryostat vacuum
     G4VPhysicalVolume* physiCloverBracket;       // pointer to the physical cryostat vacuum


     G4Tubs*            solidBCSpipe;            // pointer to the solid BCSpipe
     G4LogicalVolume*   logicBCSpipe;            // pointer to the logical BCS pipe
     G4VPhysicalVolume* physiBCSpipe;            // pointer to the physical BCS pipe

     G4Tubs*            solidBCSpipefaceplate;            // pointer to the solid BCSpipe
     G4LogicalVolume*   logicBCSpipefaceplate;            // pointer to the logical BCS pipe
     G4VPhysicalVolume* physiBCSpipefaceplate;            // pointer to the physical BCS pipe

     G4Tubs*            solidBCSvac;            // pointer to the solid BCSvac
     G4LogicalVolume*   logicBCSvac;            // pointer to the logical BCS vac
     G4VPhysicalVolume* physiBCSvac;            // pointer to the physical BCS vac

     G4Box*            solidBCSdssd;            // pointer to the solid BCS dssd
     G4LogicalVolume*   logicBCSdssd;            // pointer to the logical BCS dssd
     G4VPhysicalVolume* physiBCSdssd;            // pointer to the physical BCS dssd

     G4Box*            solidPIN;            // pointer to the solid BCS dssd
     G4LogicalVolume*   logicPIN;            // pointer to the logical BCS dssd
     G4VPhysicalVolume* physiPIN;            // pointer to the physical BCS dssd

     G4Tubs*            solidBeamEndCap;            // pointer to the solid end of the beamline
     G4LogicalVolume*   logicBeamEndCap;            // pointer to the logical end of the beamline
     G4VPhysicalVolume* physiBeamEndCap;            // pointer to the physical end of the beamline

     G4Box*            solidBCSPlastic;            // pointer to the solid BCS plastic material
     G4LogicalVolume*   logicBCSPlastic;            // pointer to the logical plastic material
     G4VPhysicalVolume* physiBCSPlastic; 

     //BC - LaBr3 solids, logicals, and physical volumes
     G4Tubs*            solidLaBr3Crys;           // pointer to solid LaBr3 scintillator crystal
     G4LogicalVolume*   logicLaBr3Crys;           // pointer to the logical LaBr3 scintillator crystal
     G4VPhysicalVolume* physiLaBr3Crys;           // pointer to the physical LaBr3 scintillator crystal
     G4Tubs*            solidLaBr3AlCover;        // pointer to the solid LaBr3 Al cover
     G4LogicalVolume*   logicLaBr3AlCover;        // pointer to the logical LaBr3 Al cover
     G4VPhysicalVolume* physiLaBr3AlCover;        // pointer to the physical LaBr3 Al cover
     G4Tubs*            solidLaBr3PMT;            // pointer to the solid LaBr3 PMT
     G4LogicalVolume*   logicLaBr3PMT;            // pointer to the logical LaBr3 PMT
     G4VPhysicalVolume* physiLaBr3PMT;            // pointer to the physical LaBr3 PMT
     G4Tubs*            solidLaBr3PMTvac;         // pointer to the solid LaBr3 PMT vac
     G4LogicalVolume*   logicLaBr3PMTvac;         // pointer to the logical LaBr3 PMT vac
     G4VPhysicalVolume* physiLaBr3PMTvac;         // pointer to the physical LaBr3 PMT vac
     
     //BC - LaBr3 frame solid, logical, and physical volume
     G4Tubs*            solidLaBr3Frame;          // pointer to solid LaBr3 frame
     G4LogicalVolume*   logicLaBr3Frame;          // pointer to the logical LaBr3 frame
     G4VPhysicalVolume* physiLaBr3Frame;          // pointer to the physical LaBr3 frame
     
     //BC - 3Hen solids, logicals, and physical volumes
     G4Tubs*            solid3HenHDPE;            // pointer to solid 3Hen HDPE (high-density polyethylene)
     G4LogicalVolume*   logic3HenHDPE;            // pointer to the logical 3Hen HDPE
     G4VPhysicalVolume* physi3HenHDPE;            // pointer to physical 3Hen HDPE
     G4Tubs*            solid3HenTube1in;         // solid 3Hen 10 atm 3He 1 inch tube
     G4LogicalVolume*   logic3HenTube1in;         // logical 3Hen 10 atm 3He 1 inch tube
     G4VPhysicalVolume* physi3HenTube1in;         // physical 3Hen 10 atm 3He 1 inch tube
     G4Tubs*            solid3HenTubeCover1in;    // solid stainless steel 1 inch tube cover
     G4LogicalVolume*   logic3HenTubeCover1in;    // logical stainless steel 1 inch tube cover
     G4VPhysicalVolume* physi3HenTubeCover1in;    // physical stainless steel 1 inch tube cover
     G4Tubs*            solid3HenTube2in;         // solid 3Hen 10 atm 3He 2 inch tube
     G4LogicalVolume*   logic3HenTube2in;         // logical 3Hen 10 atm 3He 2 inch tube
     G4VPhysicalVolume* physi3HenTube2in;         // physical 3Hen 10 atm 3He 2 inch tube
     G4Tubs*            solid3HenTubeCover2in;    // solid stainless steel 2 inch tube cover
     G4LogicalVolume*   logic3HenTubeCover2in;    // logical stainless steel 2 inch tube cover
     G4VPhysicalVolume* physi3HenTubeCover2in;    // physical stainless steel 2 inch tube cover
     
     //BC - MTC solids, logicals, and physical volumes
     G4Box*             solidMTCbox;            // pointer to solid MTC large box
     G4LogicalVolume*   logicMTCbox;            // pointer to the logical MTC box
     G4VPhysicalVolume* physiMTCbox;            // pointer to physical MTC box
     G4Box*             solidMTCtapeS1_1;       // pointer to solid top (1) tape section 1
     G4LogicalVolume*   logicMTCtapeS1_1;       // pointer to the logical top (1) tape section 1
     G4VPhysicalVolume* physiMTCtapeS1_1;       // pointer to physical top (1) tape section 1
     G4Box*             solidMTCtapeS1_2;       // pointer to solid bottom (2) tape section 1
     G4LogicalVolume*   logicMTCtapeS1_2;       // pointer to the logical bottom (2) tape section 1
     G4VPhysicalVolume* physiMTCtapeS1_2;       // pointer to physical bottom (2) tape section 1
     G4Box*             solidMTCtapeS2_1;       // pointer to solid top (1) tape section 2
     G4LogicalVolume*   logicMTCtapeS2_1;       // pointer to the logical top (1) tape section 2
     G4VPhysicalVolume* physiMTCtapeS2_1;       // pointer to physical top (1) tape section 2
     G4Box*             solidMTCtapeS2_2;       // pointer to solid bottom (2) tape section 2
     G4LogicalVolume*   logicMTCtapeS2_2;       // pointer to the logical bottom (2) tape section 2
     G4VPhysicalVolume* physiMTCtapeS2_2;       // pointer to physical bottom (2) tape section 2
     G4Box*             solidMTCtapeS3;         // pointer to solid bottom (2) tape section 3
     G4LogicalVolume*   logicMTCtapeS3;         // pointer to the logical bottom (2) tape section 3
     G4VPhysicalVolume* physiMTCtapeS3;         // pointer to physical bottom (2) tape section 3
     
     
     //BC - LENDA solids, logicals, and physical volumes
     G4Box*             solidLENDAscintbar;     // pointer to the solid LENDA scintillator bar
     G4LogicalVolume*   logicLENDAscintbar;     // pointer to the logical LENDA scintillator bar
     G4VPhysicalVolume* physiLENDAscintbar;     // pointer to the physical LENDA scintillator bar
     G4Tubs*            solidH6410_1;           // pointer to the solid LENDA H6410 PMT assembly
     G4LogicalVolume*   logicH6410_1;           // pointer to the logical LENDA H6410 PMT assembly
     G4VPhysicalVolume* physiH6410_1;           // pointer to the physical LENDA H6410 PMT assembly
     G4Tubs*            solidH6410_2;           // pointer to the solid LENDA H6410 PMT assembly
     G4LogicalVolume*   logicH6410_2;           // pointer to the logical LENDA H6410 PMT assembly
     G4VPhysicalVolume* physiH6410_2;           // pointer to the physical LENDA H6410 PMT assembly
     
     //BC - Getest solid, logical, and physical volume
     G4Box*             solidGetest;            // pointer to the solid Ge cube for effic test
     G4LogicalVolume*   logicGetest;            // pointer to the logical Ge cube for effic test
     G4VPhysicalVolume* physiGetest;            // pointer to the physical Ge cube for effic test
     
     //BC - EXOtest solid, logical, and physical volume
     G4Box*             solidEXOtest;           // pointer to the solid Xe cube for effic test
     G4LogicalVolume*   logicEXOtest;           // pointer to the logical Xe cube for effic test
     G4VPhysicalVolume* physiEXOtest;           // pointer to the physical Xe cube for effic test
     
     //BC - EJ204 Scint solid, logical, and physical volume
     //treating as one for now... will make pixels later if needed
     G4Box*             solidEJ204Scint;        // pointer to the solid EJ204 scint
     G4LogicalVolume*   logicEJ204Scint;        // pointer to the logical EJ204 scint
     G4VPhysicalVolume* physiEJ204Scint;        // pointer to the physical EJ204 scint
     
     //BC - PSPMT Scint solid, logical, and physical volume
     //treating as one for now... will make pixels later if needed
     G4Box*             solidPSPMT;        // pointer to the solid PSPMT  
     G4LogicalVolume*   logicPSPMT;        // pointer to the logical PSPMT  
     G4VPhysicalVolume* physiPSPMT;        // pointer to the physical PSPMT  

     G4Box*             solidPSPMTVac;            // pointer to the solid PSPMT vacuum
     G4LogicalVolume*   logicPSPMTVac;            // pointer to the logical PSPMT vacuum
     G4VPhysicalVolume* physiPSPMTVac;            // pointer to the physical PSPMT vacuum

     G4Box*             solidPSPMTWindow;        // pointer to the solid PSPMT window
     G4LogicalVolume*   logicPSPMTWindow;        // pointer to the logical PSPMT window
     G4VPhysicalVolume* physiPSPMTWindow;        // pointer to the physical PSPMT window

     G4Box*             solidPSPMTCathode;        // pointer to the solid PSPMT cathode
     G4LogicalVolume*   logicPSPMTCathode;        // pointer to the logical PSPMT cathode
     G4VPhysicalVolume* physiPSPMTCathode;        // pointer to the physical PSPMT cathode
     
     //BC - SiDSSD solid, logical, and physical volume
     G4Box*             solidSiDSSD;        // pointer to the solid SiDSSD
     G4LogicalVolume*   logicSiDSSD;        // pointer to the logical SiDSSD
     G4VPhysicalVolume* physiSiDSSD;        // pointer to the physical SiDSSD

     //BC - CeBr3 Scint solid, logical, and physical volume
     //treating as one for now... will make pixels later if needed
     G4Box*             solidCeBr3Scint;        // pointer to the solid CeBr3 scint
     G4LogicalVolume*   logicCeBr3Scint;        // pointer to the logical CeBr3 scint
     G4VPhysicalVolume* physiCeBr3Scint;        // pointer to the physical CeBr3 scint

     
     
     G4Material* DefaultMater;                    // Default material
     G4Material* GeThickDetectorMater;            // Thick Ge Detector material
     G4Material* GeThickDetectorCryoMater;        // Cryostat material around Ge detector
     G4Material* GeThinDetectorMater;             // Thick Ge Detector material
     G4Material* SegaCryoMater;                   // SeGA cryostat material
     G4Material* SegaVacuumMater;                 // SeGA vacuum material
     G4Material* SegaMater;                       // SeGA material
     G4Material* CloverCryoMater;                 // Clover cryostat material
     G4Material* CloverBracketMater;              // Clover bracker material
     G4Material* VacuumMater;                     // vacuum material
     G4Material* CloverCrystalMater;              // CLover crystal material
     G4Material* BCSdssdMater;                    // DSSD material
     G4Material* BCSPlasticMater;                 //plastic for the stands behind the DSSD
     G4Material* LaBr3ScintMater;                 //LaBr3 scintillator
     G4Material* LaBr3CoverMater;                 //LaBr3 Aluminum covering
     G4Material* LaBr3FrameMater;                 //LaBr3 Aluminum frame
     G4Material* HDPE_3Hen_Mater;                 //3Hen HDPE
     G4Material* Stainless_3Hen_Mater;            //3Hen Stainless steel tube
     G4Material* He3_10atm_3Hen_Mater;            //3Hen 10 atm 3He
     G4Material* tape_35mm_Mater;                 //MTC 35mm tape
     G4Material* BC408_Mater;                     //LENDA BC408 scintillator
     G4Material* H6410_Mater;                     //Using Al covering for H6410
     G4Material* Getest_Mater;                    //Ge for test cube
     G4Material* EXOtest_Mater;                   //Xe for EXO test cube
     G4Material* EJ204Scint_Mater;                //EJ204 Scintillator material
     G4Material* PSPMT_Mater;                     //PSPMT material
     G4Material* PSPMTWindow_Mater;               //PSPMT window material
     G4Material* PSPMTCathode_Mater;              //PSPMT cathode material
     G4Material* SiDSSD_Mater;                    //SiDSSD material
     G4Material* CeBr3Scint_Mater;                //CeBr3 Scintillator material
     

     G4double fWorldLength;                       // Full length the world volume
     G4double fWorldRadius;                       // Radius of  the world volume

     G4double fGeThickDetectorThickness;          // Thickness of thick Ge detector
     G4double fGeThinDetectorThickness;           // Thickness of thin Ge detector
     G4double fGeThickDetectorCryoThickness;      // Thickness of the cryostat
     G4double fGeCryoEndCapThickness;             // Thickness of the end caps

     G4double fGeThickDetectorRadius;             // Radius of the thick Ge detector
     G4double fGeThinDetectorRadius;              // Radius of the thin Ge detector
     G4double fGeThickDetectorCryoOutRad;         // Outer radius of the detector cryostat
     G4double fGeThickDetectorCryoInnRad;         // Inner radius of the detector cryostat
     G4double fGeCryoEndCapOutRad;                // Thickness of the cryostat endcap

     G4double fPipeOutRad;                        // Outer radius of the e14057 detector pipe
     G4double fPipeInRad;                         // Inner radius of the e14057 detector pipe
     G4double fPipeLength;                        // Length of the e14057 detector pipe

     G4double fPipeFacePlateOutRad;                        // Outer radius of the e14057 detector pipe
     G4double fPipeFacePlateInRad;                        // Outer radius of the e14057 detector pipe
    G4double fPipeFacePlateDepth;                        // Length of the e14057 detector pipe
    G4double fPipe_zpos;
    G4double fPipeFacePlate_zpos;  

     G4double fSpacing;                           //Interdetector spacing

     G4bool fUseGeThinDetector;                   //flag to determine if thin detector is used
     G4bool fUseGeThickDetector;                  //flag to determine if thick detector is used
     G4bool fUseGeThickDetectorCryo;              //flag to determine if detector cryostat is used
     G4bool fUseGeThickDetectorCryoEndCap;        //flag to determine if detector endcapis used
     G4bool fUseSega;                             //flag to determine if SeGA is used
     G4bool fUseClover;                           //flag to determine if Clover is used
  G4bool fUseCloverBracket;                       // flag to determine if Clover Bracket is used
     G4bool fUseLaBr3;                            //flag to determine if LaBr3 is used
     G4bool fUseLaBr3Frame;                       //flag to determine if LaBr3 frame is used
     G4bool fUse3Hen;                             //flag to determine if 3Hen is used
     G4bool fUseMTC;                              //flag to determine if MTC is used
     G4bool fUseLENDA;                            //flag to determine if LENDA is used
     G4bool fUseGetest;                           //flag to determine if Getest is used
     G4bool fUseEXOtest;                          //flag to determine if EXOtest is used
     G4bool fUseEJ204Scint;                       //flag to determine if EJ204 Scint is used
     G4bool fUsePSPMT;                            //flag to determine if PSPMT is used
     G4bool fUseSiDSSD;                           //flag to determine if SiDSSD is used
     G4bool fUsePipe;                             //flag to determine if beam pipe is used
     G4bool fUsePipeFacePlate;                             //flag to determine if beam pipe is used
  G4bool fUseCeBr3Scint;                       //flag to determine if CeBr3 Scint is used
  


     G4ThreeVector fGeThickDetectorPos;           //location in cryostat of thick detector
     G4ThreeVector fGeThinDetectorPos;            //location in cryostat of thin detector
     G4ThreeVector fGeThickDetectorCryoPos;       // location of the cryostat in the world
     G4ThreeVector fGeCryoEndCapPos1;             // location of the endcap in the cryostat
     G4ThreeVector fGeCryoEndCapPos2;             // location of the other endcap in the cryostat

     G4ThreeVector fSegaCryoPos;                  //location of the sega cryo
     G4ThreeVector fSegaVac2Pos;                  // location of the vacuum region inside the sega cryo
     G4ThreeVector fSegaDetCupPos;                // location of the detector cup inside the vacuum region
     G4ThreeVector fSegaVac1Pos;                  // location of the vacuum region inside the detector cup
     G4ThreeVector fSegaCrystalPos;               // location of the crystal inside the vacuum region
     G4ThreeVector fSegaActivePos;                // location of the active area inside the crystal
     G4ThreeVector fSegaDLinnerPos;                  //location of the inner dead layer of sega
     G4ThreeVector fSegaCentralConPos;            // location of the central contact inside the active area
     G4ThreeVector fPSPMTVacPos;                  // location of the vacuum for the PSPMT (so not a solid block) 

     G4double fsegacryolength;                    // The length of the SeGA cryostat
     G4double fsegacryothick;                     // The thickness of the cryostat wall
     G4double fsegacryooutrad;                    // The outer radius of the cryostat
     G4double fsegavac2length;                    // The length of the vacuum inside the SeGA cryostat
     G4double fsegavac2thick;                     // The thickness of the vacuum inside the SeGA cryostat
     G4double fsegadetcuplength;                  // The length of the detector cup surounding the crystal
     G4double fsegadetcupthick;                   // The thickness of the detector cup
     G4double fsegavac1length;                    // The length of the vacuum inside the SeGA detector cup
     G4double fsegavac1thick;                     // The thickness of the vacuum inside the SeGA detector cup
     G4double fsegacrystallength;                 // The length of the detector crystal
     G4double fsegacrystalthick;                  // The thickness of the crystal
     G4double fsegaactivelength;                  // The length of the active area
     G4double fsegaactivethick;                   // The thickness of the active area
     G4double fsegaDLinnerthick;                  //thickness of the inner dead layer 
     G4double fsegaDLinnerlength;                 //length of the inner dead layer
     G4double fsegacentralconlength;              // The length of the central contact
     G4double fsegacentralconthick;               // The thickness of the central contact
     G4double fBCSPlasticlength;                  //width of the BCSplatic bits
     G4double fBCSPlasticSpacing;                 //spacing between BCS blocks
     G4double fPSPMTVacHeight;                    //height of PSPMT vacuum
     G4double fPSPMTVacWidth;                     //width of PSPMT vacuum
     G4double fPSPMTVacDepth;                     //depth of PSPMT vacuum
  G4double CloverCrystalLength2;
 
     G4ThreeVector CloverCrystalPos;              // location of the Clover crystal inside the vacuum
  G4ThreeVector CloverCrystalDeadPos;
  G4ThreeVector CloverCrystalWasherPos;
     G4ThreeVector CloverCryoPos;                 // location of the Clover cry in the world
     G4ThreeVector CloverCryoVacuumPos;           // location of the CloverSS vacuum inside the cryostat

     G4ThreeVector CloverBracketPos;              //location of the brackets for Clovers

     //BC - LaBr3 array constants needed
     G4double flabr3cryslength;                   // The length of the LaBr3 scintillator crystal
     G4double flabr3crysthick;                    // The thickness of the LaBr3 scintillator crystal
     G4double flabr3crysoutrad;                   // The outer radius of the LaBr3 scitillator crystal
     G4double flabr3AlCoverlength;                // The length of the LaBr3 Al cover
     G4double flabr3AlCoverthick;                 // The thickness of the LaBr3 Al cover
     G4double flabr3AlCoveroutrad;                // The outer radius of the LaBr3 Al cover
     G4double flabr3AlCoverinrad;                 // The inner radius of the LaBr3 Al cover
     G4double flabr3PMToutrad;                    // The outer radius of the LaBr3 PM tube
     G4double flabr3PMTlength;                    // The length of the LaBr3 PM tube
     
     //BC - LaBr3 frame constants needed
     G4double flabr3framelength;                  // The length of the LaBr3 frame
     G4double flabr3framethick;                   // The thickness of the LaBr3 frame
     G4double flabr3frameinrad;                   // The inner radius of the LaBr3 frame
     G4double flabr3frameoutrad;                  // The outer radius of the LaBr3 frame


  G4double CloverBracketLength;
  G4double CloverBracketThickness;
  G4double CloverBracketOutx;
  G4double CloverBracketOuty;
  G4double CloverBracketOutz;
  G4double CloverBracketInx;
  G4double CloverBracketIny;
  G4double CloverBracketInz;
     
     
     //BC - 3Hen constants as defined in http://www.phy.ornl.gov/hribf/equipment/leribss/3hen-prop.pdf
     G4double f3HenRmin;                          // inner hole radius
     G4double f3HenRmax;                          // outer radius of HDPE
     G4double f3HenR1;                            // distance from center to first ring of 1 in tubes
     G4double f3HenR2;                            // distance from center to first ring of 2 in tubes
     G4double f3HenR3;                            // distance from center to second ring of 2 in tubes
     G4double f3HenR4;                            // distance from center to third ring of 3 in tubes
     G4double f3HenTubeoutrad1in;                 // outer radius of 3He 1 inch tube
     G4double f3HenTubeCoveroutrad1in;            // outer radius of 3He 1 inch tube cover
     G4double f3HenTubeCoverinnerrad1in;          // inner radius of 3He 1 inch tube cover
     G4double f3HenTubeoutrad2in;                 // outer radius of 3He 2 inch tube
     G4double f3HenTubeCoveroutrad2in;            // outer radius of 3He 2 inch tube cover
     G4double f3HenTubeCoverinnerrad2in;          // inner radius of 3He 2 inch tube cover
     G4double f3HenLength;                        // length of 3Hen
     
     //BC - MTC constants
     G4double fMTCtapeL1;                         // length of first section of tape from box to bend
     G4double fMTCtapeL2;                         // length of second section of tape from bend to target
     G4double fMTCtapeL3;                         // length of tape at target position connecting tape sections
     G4double fMTCtapeWidth;                      // 35mm tape is the width
     G4double fMTCtapeThickness;                  // thickness of 35mm tape (13 mm to 75 mm)
     G4double fMTCboxX;                           // MTC box X coordinate
     G4double fMTCboxY;                           // MTC box Y coordinate
     G4double fMTCboxZ;                           // MTC box Z coordinate
     G4double fMTCboxxpos;                        // MTC box X position
     G4double fMTCboxypos;                        // MTC box Y position
     G4double fMTCboxzpos;                        // MTC box Z position
     G4double fMTCtapeSep;                        // distance between top and bottom tape
     
     //BC - LENDA constants
     G4double fscintHeight;                       // height of scintillator bar
     G4double fscintWidth;                        // width of scintillator bar
     G4double fscintLength;                       // length of scintillator bar
     G4double fH6410Diamater;                     // diameter of H6410 PMT assembly
     G4double fH6410Height;                       // height of H6410 PMT assembly
     G4double fH6410PMTDiamter;                   // diamter of H6410 PMT
     G4double fH6410PMTHeight;                    // height of H6410 PMT 
     G4double fLENDAradius;                       // radial distance from target position to LENDA
     
     //BC - Getest constants
     G4double fGetestHeight;                      // height of Ge cube for effic test
     G4double fGetestWidth;                       // width of Ge cube for effic test
     G4double fGetestDepth;                       // depth of Ge cube for effic test
     
     //BC - EXOtest constants
     G4double fEXOtestHeight;                     // height of Xe cube for effic test
     G4double fEXOtestWidth;                      // width of Xe cube for effic test
     G4double fEXOtestDepth;                      // depth of Xe cube for effic test
     
     //BC - EJ204 Scintillator constants
     G4double fEJ204ScintHeight;                  // height of EJ204 Scint
     G4double fEJ204ScintWidth;                   // width of EJ204 Scint
     G4double fEJ204ScintDepth;                   // depth of EJ204 Scint
     G4double fEJ204Scint_zpos;
     
     //BC - PSPMT constants
     G4double fPSPMTHeight;                  // height of PSPMT
     G4double fPSPMTWidth;                   // width of PSPMT 
     G4double fPSPMTDepth;                   // depth of PSPMT 
     G4double fPSPMT_zpos;                   // z position of PSPMT 
     G4double fPSPMTWindowHeight;            // height of PSPMT window
     G4double fPSPMTWindowWidth;             // width of PSPMT window
     G4double fPSPMTWindowDepth;             // depth of PSPMT window
     G4double fPSPMTWindow_zpos;             // z position of PSPMT window
     G4double fPSPMTCathodeHeight;           // height of PSPMT cathode
     G4double fPSPMTCathodeWidth;            // width of PSPMT cathode
     G4double fPSPMTCathodeDepth;            // depth of PSPMT cathode
     G4double fPSPMTCathode_zpos;            // z position of PSPMT cathode
     
     //BC - SiDSSD constants
     G4double fSiDSSDHeight;                      // height of SiDSSD
     G4double fSiDSSDWidth;                       // width of SiDSSD
     G4double fSiDSSDDepth;                       // depth of SiDSSD
     G4double fSiDSSD_zpos;

     //BC - CeBr3 Scintillator constants
     G4double fCeBr3ScintHeight;                  // height of CeBr3 Scint
     G4double fCeBr3ScintWidth;                   // width of CeBr3 Scint
     G4double fCeBr3ScintDepth;                   // depth of CeBr3 Scint
     G4double fCeBr3Scint_zpos;
  G4double zOffset;
  G4double yOffset;
  G4double downOffset;
  G4double upOffset;

  //G4Region*   targetRegion;
     G4Region*   GeThickDetectorRegion;           

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
