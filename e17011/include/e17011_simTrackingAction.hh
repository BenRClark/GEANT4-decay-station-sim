#ifndef e17011_simTrackingAction_h
#define e17011_simTrackingAction_h 1

#include "G4UserTrackingAction.hh"


class e17011_simTrackingAction : public G4UserTrackingAction {

  public:
    e17011_simTrackingAction(){};
    virtual ~e17011_simTrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

};

#endif
