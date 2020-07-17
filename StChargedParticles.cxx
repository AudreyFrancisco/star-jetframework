// ################################################################
// Author:  Audrey Francisco for the STAR Collaboration
// Affiliation: Yale University
//
// ################################################################
// $Id$

#include "StChargedParticles.h"
#include "StMemStat.h"

// C++ includes
#include <sstream>
#include <fstream>

// ROOT includes
#include <TChain.h>
#include <TClonesArray.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include "TFile.h"
#include "TVector3.h"

// StRoot includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoArrays.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

// jet-framework includes
#include "StPicoConstants.h"
#include "runlistP12id.h" // Run12 pp
#include "runlistP16ij.h"
#include "runlistP17id.h" // SL17i - Run14, now SL18b (March20)
#include "runlistRun14AuAu_P18ih.h" // new Run14 AuAu
#include "runlistIso.h"
#include "StJetFrameworkPicoBase.h"
#include "StCentMaker.h"

// towers/clusters related includes:
#include "StMuDSTMaker/COMMON/StMuTrack.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"


class StJetFrameworkPicoBase;

ClassImp(StChargedParticles)
//
//________________________________________________________________________
StChargedParticles::StChargedParticles() :
  StMaker()
{
  // Default constructor.
  cout << "StChargedParticles::DefaultConstructor()\n";
}
//
//________________________________________________________________________
StChargedParticles::StChargedParticles(const char *name, bool doHistos = kFALSE, const char* outName = "") :
//  StJetFrameworkPicoBase(name),
  StMaker(name)
{
  cout << "StChargedParticles::StandardConstructor()\n";
  // Standard constructor.
  if (!name) return;
  SetName(name);
}
//
//________________________________________________________________________
StChargedParticles::~StChargedParticles()
{
  // free up histogram objects if they exist
}
//
//_____________________________________________________________________________
Int_t StChargedParticles::Init() {
cout << "StChargedParticles::Init()\n";
  return kStOK;
}
//
// finish running - write to output file and close
//_____________________________________________________________________________
Int_t StChargedParticles::Finish() {
  cout << "StChargedParticles::Finish()\n";
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}
//
//________________________________________________________________________
void StChargedParticles::Clear(Option_t *opt) {

  cout << "StChargedParticles::Clear()\n";
}
//
// Main loop, called for each event.
//________________________________________________________________________
int StChargedParticles::Make()
{
  cout << "StChargedParticles::Make()\n";
  // get PicoDstMaker
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer
  mBaseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!mBaseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  return kStOK;
}
//________________________________________________________________________
