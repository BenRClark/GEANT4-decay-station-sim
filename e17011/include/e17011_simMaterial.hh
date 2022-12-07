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
#ifndef e17011_simMaterial_HH
#define e17011_simMaterial_HH
////////////////////////////////////////////////////////////////////////////////
#include "G4Material.hh"
#include "G4NistManager.hh"
#include <vector>

class e17011_simMaterialMessenger;
////////////////////////////////////////////////////////////////////////////////
//
class e17011_simMaterial
{
public:

  e17011_simMaterial ();
  ~e17011_simMaterial ();

public:

  //  void  AddMaterial (G4String, G4String, G4double, G4String, 
  //	     G4double tem = STP_Temperature, 
  //	     G4double pres = STP_Pressure);
  //  G4Material* GetMaterial (G4int i)  {return Material[i];};
  G4Material* GetMaterial (G4String name) {return G4Material::GetMaterial(name);};
  //  G4int GetMaterialIndex (G4String);
  //  G4int GetNbOfMaterial () {return Material.size();};
  //  void  DeleteMaterial (G4int);
  //  void  DeleteMaterial (G4String);

  void  ListMaterial();

private:

  e17011_simMaterialMessenger         *materialMessenger;


  // Elements
  
  G4Element* elementH;
  G4Element* elementHe3;
  G4Element* elementC;
  G4Element* elementN;
  G4Element* elementO;
  G4Element* elementMg;
  G4Element* elementAl;
  G4Element* elementSi;
  G4Element* elementTi;
  G4Element* elementV;
  G4Element* elementCr;
  G4Element* elementMn;
  G4Element* elementNi;
  G4Element* elementGe;
  G4Element* elementFe;
  G4Element* elementMo;
  G4Element* elementPt;
  G4Element* elementAu;
  //BC - add for La and Br
  G4Element* elementLa;
  G4Element* elementBr;
  //BC - add for Xe
  G4Element* elementXe;
  //BC - add for PSPMT materials
  G4Element* elementB;
  G4Element* elementSb;
  G4Element* elementCs;
  G4Element* elementNa;
  G4Element* elementK;
  G4Element* elementCe;


  // Materials

  G4Material* vacuum;
  G4Material* polyethylene;
  G4Material* Al;
  G4Material* Nb;
  G4Material* C;
  G4Material* Au;
  G4Material* Ir;
  G4Material* Si;
  G4Material* Cr;
  G4Material* Mn;
  G4Material* Ni;
  G4Material* Ge;
  G4Material* Be;
  G4Material* Cs;
  G4Material* He3;
  G4Material* He3_10atm;
  //BC - add for La, Br, and LaBr3
  G4Material* La;
  G4Material* Br;
  G4Material* LaBr3;
  //BC - making cellulous acetate for 35 mm film (moving tape collector)
  G4Material* cellulous_acetate;
  //BC - making stainless steel
  G4Material* StainlessSteel;
  //BC - making BC408 for LENDA scintillators
  G4Material* BC408;  
  //BC - making liquid Xe
  //G4Material* Xe;
  G4Material* LiquidXeEXO;
  //BC - making EJ204 for e14057 pixelated scintillator
  G4Material* EJ204;
  //BC - making CeBr3 for e16032 pixelated scintillator
  G4Material* CeBr3;
  //BC - making materials for e14057 PSPMT
  G4Material* Borosilicate;
  G4Material* Bialkali;


//   std::vector<G4Material*>   Material;
//   std::vector<G4Element*>    Element;
//   std::vector<G4Isotope*>    Isotope;

// private:
//   static const G4String        ELU[110];
//   static const G4String        ELL[110];
//   static const G4String        EUU[110];
//   static const G4double        A[110];
       
};
////////////////////////////////////////////////////////////////////////////////
#endif
