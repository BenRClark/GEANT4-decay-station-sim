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
////////////////////////////////////////////////////////////////////////////////
//
#include "e17011_simMaterial.hh"
//#include "e17011_simMaterialData.hh"
#include "e17011_simMaterialMessenger.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include <vector>
#include <iomanip>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"  
////////////////////////////////////////////////////////////////////////////////
//
e17011_simMaterial::e17011_simMaterial ()
{

  //BC - using G4NistManager
  G4NistManager* man = G4NistManager::Instance();

  // Elements
  // G4cout<<"e17011_simMaterial"<<G4endl;

  elementH  = new G4Element("Hydrogen",  "H",  1.,  1.0079*g/mole);
  elementHe3 = new G4Element("Helium-3", "He3", 2., 3.0160293191*g/mole);
  elementC  = new G4Element("Carbon",    "C",  6.,  12.011*g/mole);
  elementN  = new G4Element("Nitrogen",  "N",  7.,  14.007*g/mole);
  elementO  = new G4Element("Oxygen",    "O",  8., 15.9994*g/mole);
  elementMg = new G4Element("Magnesium", "Mg",12., 24.3050*g/mole);
  elementAl = new G4Element("Aluminium", "Al",13., 26.9815*g/mole);
  elementSi = new G4Element("Silicon",   "Si",14., 28.0855*g/mole);
  elementTi = new G4Element("Titanium",  "Ti",22.,   47.90*g/mole);
  elementV  = new G4Element("Vanadium",  "V", 23., 50.9415*g/mole);
  elementCr = new G4Element("Chromium",   "Cr",24., 51.9961*g/mole);
  elementMn = new G4Element("Manganese", "Mn",25., 54.938045*g/mole);
  elementNi = new G4Element("Nickel",    "Ni",28., 58.6934*g/mole);
  elementGe = new G4Element("Germanium", "Ge",32.,   72.61*g/mole);
  elementFe = new G4Element("Iron",      "Fe",26.,  55.845*g/mole);
  elementMo = new G4Element("Molybdenum","Mo",42.,   95.94*g/mole);
  //elementCs = new G4Element("Cesium",    "Cs",55.,  132.90545*g/mole);
  elementPt = new G4Element("Platinum",  "Pt",78.,  195.08*g/mole);
  elementAu = new G4Element("Gold",      "Au",79.,  196.97*g/mole);
  //BC - add for La and Br
  elementLa = new G4Element("Lanthanum", "La",57.,  138.9055*g/mole);
  elementBr = new G4Element("Bromine",   "Br",35.,   79.904*g/mole);
  //BC - add for Xe
  elementXe = new G4Element("Xenon",     "Xe",54.,  131.293*g/mole); 
  //BC - add for PSPMT materials
  elementB = man->FindOrBuildElement("B");
  elementCs = man->FindOrBuildElement("Cs");
  elementSb = man->FindOrBuildElement("Sb");
  elementNa = man->FindOrBuildElement("Na");
  elementK = man->FindOrBuildElement("K");
  elementCe = new G4Element("Cerium",   "Ce",58.,  140.42*g/mole);
  
  
  // Materials

  vacuum = new G4Material("vacuum", 1, 1.00794*g/mole, 
			  1.0E-25*g/cm3, kStateGas, 0.1*kelvin, 1.0E-19*pascal);

  polyethylene = new G4Material("polyethylene", 1.36*g/cm3, 2);
  polyethylene->AddElement(elementC,1);
  polyethylene->AddElement(elementH,2);

  Al = new G4Material("Al", 13, 26.98153*g/mole, 2.70*g/cm3);
  Nb = new G4Material("Nb", 41, 92.90638*g/mole, 8.57*g/cm3);
  C  = new G4Material("C",   6, 12.011*g/mole,   2.15*g/cm3);
  Au = new G4Material("Au", 79., 196.9*g/mole, 19.32*g/cm3);
  Be = new G4Material("Be", 4, 9.012182*g/mole, 1.84*g/cm3);
  Si = new G4Material("Si", 14., 28.0855*g/mole, 2.33*g/cm3);
  Cr = new G4Material("Cr", 24., 51.9961*g/mole, 7.19*g/cm3);
  Mn = new G4Material("Mn", 25., 54.938045*g/mole, 7.21*g/cm3);
  Ni = new G4Material("Ni", 28., 58.6934*g/mole, 8.91*g/cm3);
  Ge = new G4Material("Ge", 32., 76.61*g/mole, 5.325*g/cm3);
  Cs = new G4Material("Cs", 55., 132.90545*g/mole, 1.93*g/cm3);
  Ir = new G4Material("Ir", 77., 192.217*g/mole, 22.65*g/cm3);
  
  
  He3 = new G4Material("He3", 2., 3.0160293191*g/mole, 0.0001786*g/cm3);  //0.1786 g/L, 1000cm3 in L
  He3_10atm = new G4Material("He3_10atm", 2., 3.0160293191*g/mole, 1.3456*kg/m3);  //10 atm: http://yeroc.us/calculators/gas-density.php

  //BC - add for La, Br, and LaBr3
  G4int natoms;
  La = new G4Material("La", 57., 138.9055*g/mole, 6.16*g/cm3);
  Br = new G4Material("Br", 35., 79.904*g/mole, 3.10*g/cm3);
  LaBr3 = new G4Material("LaBr3", 5.06*g/cm3, 2);
  LaBr3->AddElement(elementLa, natoms = 1);
  LaBr3->AddElement(elementBr, natoms = 3);
  
  // BC - define stainless steel
  G4double densityNIST;
  G4int ncomponents;
  G4double fractionmass;
  StainlessSteel = new G4Material("StainlessSteel", densityNIST= 8.06*g/cm3, ncomponents=6);
  StainlessSteel->AddElement(elementC, fractionmass=0.001);
  StainlessSteel->AddElement(elementSi, fractionmass=0.007);
  StainlessSteel->AddElement(elementCr, fractionmass=0.18);
  StainlessSteel->AddElement(elementMn, fractionmass=0.01);
  StainlessSteel->AddElement(elementFe, fractionmass=0.712);
  StainlessSteel->AddElement(elementNi, fractionmass=0.09);

  //BC - define PSPMT materials
  Borosilicate = new G4Material("Borosilicate glass", densityNIST= 2.23*g/cm3, ncomponents=6);
  Borosilicate->AddElement(elementB, fractionmass=0.040064);
  Borosilicate->AddElement(elementO, fractionmass=0.539562); 
  Borosilicate->AddElement(elementNa, fractionmass=0.028191);
  Borosilicate->AddElement(elementAl, fractionmass=0.011644);
  Borosilicate->AddElement(elementSi, fractionmass=0.377220);
  Borosilicate->AddElement(elementK, fractionmass=0.003321);

  Bialkali = new G4Material("Bialkali", densityNIST= 4.28*g/cm3, ncomponents=3);
  Bialkali->AddElement(elementK, fractionmass=0.133);
  Bialkali->AddElement(elementCs, fractionmass=0.452);
  Bialkali->AddElement(elementSb, fractionmass=0.415); 

	
  //BC - define cellulous acetate for 35mm film
  //formula: http://chemistry.about.com/od/factsstructures/ig/Chemical-Structures---C/Cellulose-Acetate.htm
  //density:  http://www.azom.com/article.aspx?ArticleID=383
  //molar mass = 6*C + 10*H + 5*O, which will be used for formula when separating hydroxides
	
  cellulous_acetate = new G4Material("CellulousAcetate", 1.28*g/cm3, 3);
  cellulous_acetate->AddElement(elementC, natoms = 6);
  cellulous_acetate->AddElement(elementH, natoms = 10);
  cellulous_acetate->AddElement(elementO, natoms = 5);
	
  //BC - define BC408 for LENDA
  //formula:  http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=16213094
  //density: http://www.phys.ufl.edu/courses/phy4803L/group_I/muon/bicron_bc400-416.pdf
  
  BC408 = new G4Material("BC408", 1.032*g/cm3, 2);
  BC408->AddElement(elementC, natoms = 27);
  BC408->AddElement(elementH, natoms = 30);
  
  //BC - define Liquid Xe for xenon efficiency measurements (EXO)
  LiquidXeEXO = new G4Material("LiquidXeEXO", 54., 135.5069*g/mole, 3.1*g/cm3); //80% 136Xe (108.7257752), 20% 134Xe (26.7810788) for EXO: 
  
  //BC - define EJ204 for e14057 scintillator
  //info: http://www.eljentechnology.com/index.php/products/plastic-scintillators/60-product-template
  //polymer base is Polyvinyltoluene, which is C-27, H-30 just like BC408
	
  EJ204 = new G4Material("EJ204", 1.023*g/cm3, 2);
  EJ204->AddElement(elementC, natoms = 27);
  EJ204->AddElement(elementH, natoms = 30);
	
	
  //BC - define CeBr3 for e16032 scintillator
  CeBr3 = new G4Material("CeBr3",5.1*g/cm3,2);
  CeBr3->AddElement(elementCe,natoms = 1);
  CeBr3->AddElement(elementBr,natoms = 3);

  //   Material.clear();
  //   Element.clear();
  //   Isotope.clear();
  //   // some default materials vacuum (0), air (1)  and aluminium (2)  defined here
  //   // examples of vacuum
  //   //
  //   //  G4double a,z;

  //   static G4bool bmat = false ;

  //   if (!bmat) {

  //     // vacuum
  //     G4double  density    = universe_mean_density;    //from PhysicalConstants.h
  //     G4double pressure    = 3.e-18*pascal;
  //     G4double temperature = 2.73*kelvin;
  //     AddMaterial("Vacuum", "H", density,"gas",temperature,pressure);
  //     // air
  //     density   = 1.290*mg/cm3;
  //     AddMaterial("Air", "N0.78-O0.22", density, "gas");

  //     //realistic vacuum, scaled from definition of air above 1atm (760 torr) = 1.29*mg/cm3
  //     density = 1.7e-9*mg/cm3;
  //     pressure = 1.32e-9*atmosphere;
  //     temperature = 80*kelvin;
  //     AddMaterial("LowPressAir", "N0.78-O0.22", density, "gas",temperature,pressure);
  
  //     // aluminium
  //     density=2.700*g/cm3 ;
  //     AddMaterial ("Aluminium", "Al", density,"");
    
  //      //silicon
  //     density=2.3290*g/cm3 ;
  //     AddMaterial ("Silicon", "Si", density,"");

  //     bmat              = true;
  //   }




  // create commands for interactive definition of the geometry 
  materialMessenger = new e17011_simMaterialMessenger(this);
}
////////////////////////////////////////////////////////////////////////////////
//
e17011_simMaterial::~e17011_simMaterial ()
{
  delete materialMessenger;
  // G4cout<<"~e17011_simMaterial"<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
/*
  void e17011_simMaterial::AddMaterial (G4String name, G4String formula, G4double density,
  G4String state, G4double tem, G4double pres)
  {
  G4int isotope, Z;
  size_t i;
  for (i = 0; i<Material.size(); i++) {
  if (Material[i]->GetName() == name) {
  G4cerr <<" AddMaterial : material " <<name
  <<" already exists." <<G4endl;
  G4cerr <<"--> Command rejected." <<G4endl;
  return;
  }
  }
  
  char *tokenPtr1 = NULL;
  char *sname     = NULL;
  G4String s, s1("0123456789");
  G4String element, isotopename;
  G4int ncomponents, natoms;
  G4double fatoms = 0.;
  size_t ls, id=0, ll, lr;
  ncomponents = 0;

  sname       = new char[strlen(formula)+1];
  strcpy(sname,formula);
  tokenPtr1 = strtok(sname,"-");

  while (tokenPtr1 != NULL) {
  ncomponents++;
  tokenPtr1 = strtok( NULL, "-");
  }
  delete[] sname;

  G4Material* aMaterial = 0;
  //  G4cout << name <<" "<< formula << " " << density/(g/cm3) << " " << tem <<" " <<pres << G4endl;
 
  if (state == "") {
  aMaterial = new G4Material(name, density, ncomponents);
  } else if (state == "solid" && tem > 0.) {
  aMaterial = new G4Material(name, density, ncomponents, 
  kStateSolid, tem );
  } else if (state == "gas" && pres > 0.) {
  aMaterial = new G4Material(name, density, ncomponents, 
  kStateGas, tem, pres );
  }
  if (aMaterial == 0) {
  G4cerr <<" AddMaterial : Name " <<name <<"." <<G4endl;
  G4cerr <<"--> Command failed." <<G4endl;
  return;
  }

  sname=new char[strlen(formula)+1];
  strcpy(sname,formula);
  tokenPtr1 = strtok(sname,"-");

  while (tokenPtr1 != NULL) {
  isotope = 0;
  //      G4cout << tokenPtr1 << G4endl;
  s       = G4String(tokenPtr1);
  ls      = s.length();
  ll      = s.find("(");
  lr      = s.find(")");
  if (ll == lr) {
  id = s.find_first_of(s1);
  element = s.substr(0,id);
      
  if (element.length() == 1) element.insert(0," ");
  for (i = 0; i<110; i++) {
  if (element == ELU[i]) break;
  }
  if (i == 110) {
  for (i = 0; i<110; i++) {
  if (element == ELL[i]) break;
  }
  if (i == 110) {
  for (i = 0; i<110; i++) {
  if (element == EUU[i]) break;
  }
  }
  }
      
  if (i == 110) {
  G4cerr <<"AddMaterial : Invalid element in material formula."
  <<element <<G4endl;
  G4cerr <<"--> Command rejected." <<G4endl;
  //        delete aMaterial;
  //	Material[NbMat] = NULL;
  return;
  }

  Z       = i+1;
  element = ELU[i];
  if (id == std::string::npos) {
  natoms = 1;
  } else {
  natoms = atoi((s.substr(id,ls-id)).c_str());
  }
  if (natoms < 1) fatoms = atof((s.substr(id,ls-id)).c_str());
  //	G4cout << "   Elements = " << element << G4endl;
  //G4cout << "   Nb of atoms = " << natoms << G4endl;
  } else {
  element = s.substr(0,ll);
  isotope = atoi((s.substr(ll+1,lr-ll)).c_str());
  if (element.length() == 1) element.insert(0," ");
  for (i = 0; i<110; i++) {
  if (element == ELU[i]) break;
  }
  if (i == 110) {
  for (i = 0; i<110; i++) {
  if (element == ELL[i]) break;
  }
  if (i == 110) {
  for (i = 0; i<110; i++) {
  if (element == EUU[i]) break;
  }
  }
  }
  if (i == 110) {
  G4cerr <<"AddMaterial : Invalid element in material formula."
  <<element <<G4endl;
  G4cerr <<"--> Command rejected." <<G4endl;
  //        delete aMaterial;
  //	Material[NbMat] = NULL;
  return;
  }

  Z           = i+1;
  element     = ELU[i];
  isotopename = element+s.substr(ll+1,lr-ll-1);
  if (lr == std::string::npos ) {
  natoms = 1;
  } else {
  natoms = atoi((s.substr(lr+1,ls-lr)).c_str());
  }  
  if (natoms < 1)  fatoms = atof((s.substr(id,ls-id)).c_str());
  if (fatoms == 0.) natoms = 1;
  //
  //	G4cout << "   Elements = " << element << G4endl;
  //   G4cout << "   Isotope Nb = " << isotope << G4endl;
  //	G4cout << "   Nb of atoms = " << natoms << G4endl;
  }
  if (isotope != 0) {
  if (G4Isotope::GetIsotope(isotopename) == NULL) {
  //        G4Isotope* aIsotope = new G4Isotope(isotopename, Z, isotope, A[Z-1]*g/mole);
  G4Isotope* aIsotope = new G4Isotope(isotopename, Z, isotope, isotope*g/mole);
  G4Element* aElement = new G4Element(isotopename, element, 1);
  aElement->AddIsotope(aIsotope, 100.*perCent);
  Isotope.push_back(aIsotope);
  if (natoms>0) { 
  aMaterial->AddElement(aElement, natoms);
  } else {
  aMaterial->AddElement(aElement, fatoms);
  }
  Element.push_back(aElement);
  } else {
  if (natoms>0) { 
  aMaterial->AddElement( G4Element::GetElement(isotopename) , natoms);
  } else {
  aMaterial->AddElement( G4Element::GetElement(isotopename) , fatoms);
  }
  }      
  } else {
  if ( G4Element::GetElement(element) == NULL) {
  G4Element* aElement = new G4Element(element, element, Z, A[Z-1]*g/mole);
  if (natoms>0) { 
  aMaterial->AddElement(aElement, natoms);
  } else {
  aMaterial->AddElement(aElement, fatoms);
  }
  Element.push_back(aElement);
  } else {
  if (natoms>0) { 
  aMaterial->AddElement( G4Element::GetElement(element) , natoms);
  } else {
  aMaterial->AddElement( G4Element::GetElement(element) , fatoms);
  } 
  }
  }
  tokenPtr1 = strtok( NULL, "-");
  //      s.empty();
  //element.erase();
  //
  }

  delete[] sname;
  Material.push_back(aMaterial);
  G4cout <<" Material:" <<name <<" with formula: " <<formula <<" added! "
  <<G4endl;
  G4cout <<"     Nb of Material = " <<Material.size() <<G4endl;
  G4cout <<"     Nb of Isotope =  " <<Isotope.size() <<G4endl;
  G4cout <<"     Nb of Element =  " <<Element.size() <<G4endl;
  }
  ////////////////////////////////////////////////////////////////////////////////
  //
  void e17011_simMaterial::DeleteMaterial (G4int j)
  {
  size_t i(j-1);
  if (i > Material.size()) {
  G4cerr <<"DeleteMaterial : Invalid material index " <<j <<"." <<G4endl;
  G4cerr <<"--> Command rejected." <<G4endl;
  } else {
  G4cerr <<"It seems there is no mechanism in G4 for deleting a material yet!"
  <<G4endl;
  G4cerr <<"--> Command rejected." <<G4endl;
  }
  }
  ////////////////////////////////////////////////////////////////////////////////
  //
  void e17011_simMaterial::DeleteMaterial (G4String )
  {
  G4cerr <<"It seems there is no mechanism in G4 for deleting a material yet!"
  <<G4endl;
  G4cerr <<"--> Command rejected." <<G4endl;
  }
  ////////////////////////////////////////////////////////////////////////////////
  //
  G4int e17011_simMaterial::GetMaterialIndex (G4String name)
  {
  size_t i ;
  for (i = 0; i < Material.size(); i++) {
  if (Material[i]->GetName() == name) break;
  }
  G4int k = G4int(i);
  if (i == Material.size()) k = -1;
  return k;
  }
*/
/*
////////////////////////////////////////////////////////////////////////////////
//
void e17011_simMaterial::ListMaterial ()
{
G4cout <<" There are" <<std::setw(3) <<Material.size()
<<" materials defined."  <<G4endl;
for (size_t i = 0; i< Material.size(); i++) 
G4cout <<"     Material Index " <<std::setw(3) <<i+1 <<" "
<<std::setw(14) <<Material[i]->GetName()
<<"  density: " <<std::setw(6) <<std::setprecision(3)
<<G4BestUnit(Material[i]->GetDensity(),"Volumic Mass") <<G4endl;
G4cout <<G4endl;
  
}
////////////////////////////////////////////////////////////////////////////////
*/
