//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: e17011_simVisManager.cc,v 1.12 2003/11/06 14:37:10 johna Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//
// John Allison 24th January 1998.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4VIS_USE

#include "e17011_simVisManager.hh"

// Supported drivers...

// Not needing external packages or libraries...
#include "G4ASCIITree.hh"
#include "G4DAWNFILE.hh"
//#include "G4GAGTree.hh"
#include "G4HepRepFile.hh"
#include "G4HepRep.hh"
#include "G4RayTracer.hh"
#include "G4VRML1File.hh"
#include "G4VRML2File.hh"

// Needing external packages or libraries...

#ifdef G4VIS_USE_DAWN
#include "G4FukuiRenderer.hh"
#endif

#ifdef G4VIS_USE_OPENGLX
#include "G4OpenGLImmediateX.hh"
#include "G4OpenGLStoredX.hh"
#endif

#ifdef G4VIS_USE_OPENGLWIN32
#include "G4OpenGLImmediateWin32.hh"
#include "G4OpenGLStoredWin32.hh"
#endif

#ifdef G4VIS_USE_OPENGLXM
#include "G4OpenGLImmediateXm.hh"
#include "G4OpenGLStoredXm.hh"
#endif

#ifdef G4VIS_USE_OIX
#include "G4OpenInventorX.hh"
#endif

#ifdef G4VIS_USE_OIWIN32
#include "G4OpenInventorWin32.hh"
#endif

#ifdef G4VIS_USE_VRML
#include "G4VRML1.hh"
#include "G4VRML2.hh"
#endif

//BC - add OPENGLQT
#ifdef G4VIS_USE_OPENGLQT
#include "G4OpenGLImmediateQt.hh"
#include "G4OpenGLStoredQt.hh"
#endif



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

e17011_simVisManager::e17011_simVisManager () {
  // G4cout<<"e17011_simVisManager"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void e17011_simVisManager::RegisterGraphicsSystems () {
  // G4cout<<"e17011_simVisManager RegisterGraphicsSystems"<<G4endl;
  // Graphics Systems not needing external packages or libraries...
  RegisterGraphicsSystem (new G4ASCIITree);
  RegisterGraphicsSystem (new G4DAWNFILE);
//  RegisterGraphicsSystem (new G4GAGTree);
  RegisterGraphicsSystem (new G4HepRepFile);
  RegisterGraphicsSystem (new G4HepRep);
  RegisterGraphicsSystem (new G4RayTracer);
  RegisterGraphicsSystem (new G4VRML1File);
  RegisterGraphicsSystem (new G4VRML2File);

  // Graphics systems needing external packages or libraries...

#ifdef G4VIS_USE_DAWN
  RegisterGraphicsSystem (new G4FukuiRenderer);
#endif

//BC - add OPENGL stuff
// Register OGL graphics system with generic nicknames.
// Note: Any graphics system can be used with any UI session except
// OPENGLQT - this must have a UI Qt session unless in batch mode.  This
// is handled in /vis/sceneHandler/create but a potential fallback is
// defined here.
#ifdef G4VIS_USE_OPENGL
  G4VGraphicsSystem* ogl  = 0;
  G4VGraphicsSystem* ogli = 0;
  G4VGraphicsSystem* ogls = 0;
  G4VGraphicsSystem* ogl_fallback = 0;
  G4VGraphicsSystem* ogli_fallback = 0;
  G4VGraphicsSystem* ogls_fallback = 0;
#ifdef G4VIS_USE_OPENGLQT
  ogl  = new G4OpenGLStoredQt;
  ogli = new G4OpenGLImmediateQt;
  ogls = new G4OpenGLStoredQt;
#if defined G4VIS_USE_OPENGLXM
  ogl_fallback  = new G4OpenGLStoredXm;
  ogli_fallback = new G4OpenGLImmediateXm;
  ogls_fallback = new G4OpenGLStoredXm;
#elif defined G4VIS_USE_OPENGLX
  ogl_fallback  = new G4OpenGLStoredX;
  ogli_fallback = new G4OpenGLImmediateX;
  ogls_fallback = new G4OpenGLStoredX;
#elif defined G4VIS_USE_OPENGLWIN32
  ogl_fallback  = new G4OpenGLStoredWin32;
  ogli_fallback = new G4OpenGLImmediateWin32;
  ogls_fallback = new G4OpenGLStoredWin32;
#endif
#elif defined G4VIS_USE_OPENGLWT
  ogl  = new G4OpenGLImmediateWt;
  ogli = new G4OpenGLImmediateWt;
  ogls = new G4OpenGLImmediateWt;
#elif defined G4VIS_USE_OPENGLXM
  ogl  = new G4OpenGLStoredXm;
  ogli = new G4OpenGLImmediateXm;
  ogls = new G4OpenGLStoredXm;
#elif defined G4VIS_USE_OPENGLWIN32
  ogl  = new G4OpenGLStoredWin32;
  ogli = new G4OpenGLImmediateWin32;
  ogls = new G4OpenGLStoredWin32;
#elif defined G4VIS_USE_OPENGLX
  ogl  = new G4OpenGLStoredX;
  ogli = new G4OpenGLImmediateX;
  ogls = new G4OpenGLStoredX;
#endif

  if (ogl) {
    //ogl->SetNickname("OGL");
    RegisterGraphicsSystem (ogl);
  }
  if (ogli) {
    //ogli->SetNickname("OGLI");
    RegisterGraphicsSystem (ogli);
  }
  if (ogls) {
    //ogls->SetNickname("OGLS");
    RegisterGraphicsSystem (ogls);
  }
  if (ogl_fallback) {
    //ogl_fallback->SetNickname("OGL_FALLBACK");
    RegisterGraphicsSystem (ogl_fallback);
  }
  if (ogli_fallback) {
    //ogli_fallback->SetNickname("OGLI_FALLBACK");
    RegisterGraphicsSystem (ogli_fallback);
  }
  if (ogls_fallback) {
    //ogls_fallback->SetNickname("OGLS_FALLBACK");
    RegisterGraphicsSystem (ogls_fallback);
  }
#endif


#ifdef G4VIS_USE_OPENGLQT
  G4VGraphicsSystem* gsi = new G4OpenGLImmediateQt;
  G4VGraphicsSystem* gss = new G4OpenGLStoredQt;
  RegisterGraphicsSystem (gsi);
  RegisterGraphicsSystem (gss);
  G4VGraphicsSystem* gsi_fallback = 0;
  G4VGraphicsSystem* gss_fallback = 0;
#if defined G4VIS_USE_OPENGLXM
  gsi_fallback = new G4OpenGLImmediateXm;
  gss_fallback = new G4OpenGLStoredXm;
#elif defined G4VIS_USE_OPENGLX
  gsi_fallback = new G4OpenGLImmediateX;
  gss_fallback = new G4OpenGLStoredX;
#elif defined G4VIS_USE_OPENGLWIN32
  gsi_fallback = new G4OpenGLImmediateWin32;
  gss_fallback = new G4OpenGLStoredWin32;
#endif
  if (gsi_fallback) {
    //gsi_fallback->SetNickname(gsi->GetNickname()+"_FALLBACK");
    RegisterGraphicsSystem (gsi_fallback);
  }
  if (gss_fallback) {
    //gss_fallback->SetNickname(gss->GetNickname()+"_FALLBACK");
    RegisterGraphicsSystem (gss_fallback);
  }
#endif

#ifdef G4VIS_USE_OPENGLWIN32
  RegisterGraphicsSystem (new G4OpenGLImmediateWin32);
  RegisterGraphicsSystem (new G4OpenGLStoredWin32);
#endif

#ifdef G4VIS_USE_OIX
  RegisterGraphicsSystem (new G4OpenInventorX);
  RegisterGraphicsSystem (new G4OpenInventorXtExtended);
#endif

#ifdef G4VIS_USE_OIWIN32
  RegisterGraphicsSystem (new G4OpenInventorWin32);
#endif

#ifdef G4VIS_USE_VRML
  RegisterGraphicsSystem (new G4VRML1);
  RegisterGraphicsSystem (new G4VRML2);
#endif




  if (fVerbose > 0) {
    G4cout <<
      "\nYou have successfully chosen to use the following graphics systems."
	 << G4endl;
    PrintAvailableGraphicsSystems(GetVerbosity());
  }
}

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
