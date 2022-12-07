# $Id: GNUmakefile,v 1.4 2006/12/21 17:52:39 flei Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := e17011_sim
G4TARGET := $(name)
G4EXLIB := true
G4WORKDIR := /home/brc339/simulations/geant4_workdir

ifndef G4INSTALL
  G4INSTALL = /opt/geant4.10.04.p02
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk

ifdef ROOTSYS
#G4ANALYSIS_USE_ROOT

CPPFLAGS += -DG4ANALYSIS_USE_ROOT
CPPFLAGS += -D_REENTRANT -I$(ROOTSYS)/include/root -m64 -pthread


EXTRALIBS += -L$(ROOTSYS)/lib/root \
-lCore -lCint -lRIO -lNet\
-lHist -lGraf -lGraf3d -lGpad \
-lTree -lRint -lPostscript \
-lMatrix -lPhysics -lMathCore -lThread -lGui -pthread -lMinuit -lHtml\
-lm -ldl -lpthread -rdynamic 

#the following lines about the root flas and libs were commented out before the update
ROOTCPPFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
#CPPFLAGS +=  -DG4ANALYSIS_USE_ROOT $(ROOTCPPFLAGS)
#new in update
CPPFLAGS +=  $(ROOTCPPFLAGS)

ROOTLIBS = $(shell $(ROOTSYS)/bin/root-config --nonew --glibs)
ROOTLIBS := $(filter-out -lNew,$(ROOTLIBS))
ROOTLIBS += -lMinuit -lHtml
LDLIBS += $(ROOTLIBS)


endif

ifdef G4ANALYSIS_USE_AIDA
   CPPFLAGS += -DG4ANALYSIS_USE_AIDA
   CPPFLAGS += `aida-config --include`
   LDLIBS += `aida-config --lib`
endif

