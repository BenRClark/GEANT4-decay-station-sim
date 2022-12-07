#ifndef Materials_H
#define Materials_H 1

#include "globals.hh"
#include "G4Material.hh"

class Materials
{
  public:

  // Elements
  G4Element* elementH;
  G4Element* elementC;
  G4Element* elementN;
  G4Element* elementO;
  G4Element* elementMg;
  G4Element* elementAl;
  G4Element* elementSi;
  G4Element* elementTi;
  G4Element* elementV;
  G4Element* elementGe;
  G4Element* elementFe;
  G4Element* elementMo;
  G4Element* elementPt;
  G4Element* elementAu;

  // Materials

  G4Material* vacuum;
  G4Material* polyethylene;
  G4Material* Al;
  G4Material* Nb;
  G4Material* C;
  G4Material* Au;
  G4Material* Ir;
  G4Material* Si;
  G4Material* Ge;
  G4Material* Be;



//     G4Material* Air;
//     G4Material* Vacuum;
//     G4Material* Pb;
//     G4Material* Ar;
//     G4Material* Xenon;
//     G4Material* Ge;
//     G4Material* Al;
//     G4Material* Cu;
//     G4Material* LN2;
//     G4Material* BGO;
//     G4Material* CsI;
//     G4Material* Delrin;
//     G4Material* Polyethylene;
//     G4Material* Hevimetal;
//     G4Material* Be;
//     G4Material* Fe;
//     G4Material* Pyrex;
//     G4Material* silicon_dioxide;
//     G4Material* boron_trioxide;
//     G4Material* sodium_oxide;
//     G4Material* aluminium_oxide;
//     G4Material* potassium_oxide;
//     G4Material* calcium_oxide;
//     G4Material* dichlorine;
//     G4Material* magnesium_oxide;
//     G4Material* iron3_oxide;
//     G4Material* lanthanum_bromide;
//     G4Material* cerium_doped_lanthanum_bromide;    
//     G4Material* cadmium_zinc_telluride;    
//     G4Material* silicon;
    
//     // Collimator material
//     G4Material* WCu;
    
//     // Building materials
//     G4Material* sand;
//     G4Material* limestone;
//     G4Material* tricalcium_silicate;
//     G4Material* dicalcium_silicate;
//     G4Material* tricalcium_aluminate;
//     G4Material* tetracalcium_aluminoferrite;
//     G4Material* gypsum;
//     G4Material* K2O;
//     G4Material* Na2O;
//     G4Material* water;
//     G4Material* portland_cement;
//     G4Material* concrete;

    Materials();
    ~Materials();

//     void defineAllMaterials();

  private:

    //Made Elements
    G4Element* elH;
    G4Element* elC;
    G4Element* elN;
    G4Element* elO;
    G4Element* elAl;
    G4Element* elSi;
    G4Element* elS;
    G4Element* elCa;
    G4Element* elFe;
    G4Element* elNi;
    G4Element* elCu;
    G4Element* elGe;
    G4Element* elI;
    G4Element* elCs;
    G4Element* elTa;
    G4Element* elW;
    G4Element* elBi;
    G4Element* elB;
    G4Element* elK;
    G4Element* elNa;
    G4Element* elCl;
    G4Element* elMg;
    G4Element* elLa;
    G4Element* elBr;
    G4Element* elCe;
    G4Element* elCd;
    G4Element* elZn;
    G4Element* elTe;

//     void makeAir();
//     void makeVacuum();
//     void makeLead();
//     void makeXenon();
//     void makeArgon();
//     void makeGermanium();
//     void makeAluminum();
//     void makeCopper();
//     void makeLiquidNitrogen();
//     void makeBismuthGermanate();
//     void makeCesiumIodide();
//     void makeDelrin();
//     void makePolyethylene();
//     void makeHevimetal();
//     void makeWCu();
//     void makeBeryllium();
//     void makeIron();
//     void makeSiliconDioxide();
//     void makeBoronTrioxide();
//     void makeSodiumOxide();
//     void makeAluminiumOxide();
//     void makePotassiumOxide();
//     void makeCalciumOxide();
//     void makeDichlorine();
//     void makeMagnesiumOxide();
//     void makeIron3Oxide();
//     void makePyrex();   

//     // Building materials
//     void makeSand();
//     void makeLimestone();
//     void makeTricalciumSilicate();
//     void makeDicalciumSilicate();
//     void makeTricalciumAluminate();
//     void makeTetracalciumAluminoferrite();
//     void makeGypsum();
//     void makeK2O();
//     void makeNa2O();
//     void makeWater();
//     void makePortlandCement();
//     void makeConcrete();
    
//     void defineAllElements();
};

#endif

