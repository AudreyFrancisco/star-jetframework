#ifndef STCHARGEDPARTICLES_H
#define STCHARGEDPARTICLES_H
// $Id$

// C++ includes
#include <set>

// STAR includes
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

// ROOT classes
class TClonesArray;
class TObjArray;
class TList;
class TH1;
class TH1D;
class TH2;
class TH2D;
class TProfile;
class TVector3;

// STAR classes
class StMaker;
class StChain;
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;


// star jet-frameworks classes
class StJetFrameworkPicoBase;
class StCentMaker;

// centrality class
class StRefMultCorr;


class StChargedParticles : public StMaker {
 public:

  // debug flags for specifics
  enum fDebugFlagEnum {
    kDebugNothing, // don't want lowest elements to be used
    kDebugEmcTrigger,
    kDebugGeneralEvt,
    kDebugCentrality,
  };

  StChargedParticles();
  StChargedParticles(const char *name, bool dohistos, const char* outName);
  virtual ~StChargedParticles();

  // needed class functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t *opt="");
  virtual Int_t Finish();

  // booking of histograms (optional)
  void    DeclareHistograms();
  void    WriteHistograms();


 protected:

  StChargedParticles(const StChargedParticles&);            // not implemented
  StChargedParticles &operator=(const StChargedParticles&); // not implemented

  ClassDef(StChargedParticles, 0) // track/cluster QA task
};
#endif
