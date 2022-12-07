#ifndef e18016_simTrackingAction_h
#define e18016_simTrackingAction_h 1

#include "G4UserTrackingAction.hh"


class e18016_simTrackingAction : public G4UserTrackingAction {

  public:
    e18016_simTrackingAction(){};
    virtual ~e18016_simTrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

};

#endif
