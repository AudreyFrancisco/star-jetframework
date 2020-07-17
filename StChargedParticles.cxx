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
  StMaker(),
  doWriteHistos(kFALSE),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),
  doppAnalysis(kFALSE),
  fCorrPileUp(kFALSE),
  fDoEffCorr(kFALSE),
  fMaxEventTrackPt(30.0),
  doRejectBadRuns(kFALSE),
  fEventZVtxMinCut(-40.0),
  fEventZVtxMaxCut(40.0),
  fEventVzDiffCut(5.),
  fEventVrCut(2.),
  fCentralitySelectionCut(-99),
  fRequireCentSelection(kFALSE),
  mOutName(""),
  fAnalysisMakerName(""),
  fTracksName(""),
  fTrackPtMinCut(0.2),
  fTrackPtMaxCut(30.0),
  fTrackPhiMinCut(0.0),
  fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0),
  fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15),
  fTracknHitsRatio(0.2),
  fTracknHitsRatioMax(1.05),
  fTrackChargePos(-1),
  fGoodTrackCounter(0),
  fCentralityScaled(0.),
  fmycentral(0.),
  ref16(-99), ref9(-99),
  Bfield(0.0),
  mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fTriggerToUse(0), // kTriggerANY
  fMBEventType(2),  // kVPDMB
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mCentMaker(0x0),
  mBaseMaker(0x0)
  {
  // Default constructor.
  cout << "StChargedParticles::DefaultConstructor()\n";
}
//
//________________________________________________________________________
StChargedParticles::StChargedParticles(const char *name, bool doHistos = kFALSE, const char* outName = "") :
//  StJetFrameworkPicoBase(name),
  StMaker(name),
  doWriteHistos(doHistos),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  doppAnalysis(kFALSE),
  fCorrPileUp(kFALSE),
  fDoEffCorr(kFALSE),
  fMaxEventTrackPt(30.0),
  doRejectBadRuns(kFALSE),
  fEventZVtxMinCut(-40.0),
  fEventZVtxMaxCut(40.0),
  fEventVzDiffCut(5.),
  fEventVrCut(2.),
  fCentralitySelectionCut(-99),
  fRequireCentSelection(kFALSE),
  mOutName(outName),
  fAnalysisMakerName(name),
  fTracksName("Tracks"),
  fTrackPtMinCut(0.2),
  fTrackPtMaxCut(30.0),
  fTrackPhiMinCut(0.0),
  fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0),
  fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15),
  fTracknHitsRatio(0.2),
  fTracknHitsRatioMax(1.05),
  fTrackChargePos(-1),
  fGoodTrackCounter(0),
  fCentralityScaled(0.),
  fmycentral(0.),
  ref16(-99), ref9(-99),
  Bfield(0.0),
  mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fTriggerToUse(0),  // kTriggerANY
  fMBEventType(2),   // kVPDMB
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mCentMaker(0x0),
  mBaseMaker(0x0)
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

  // Destructor
  if(VzHist)                    delete VzHist;
  if(ZDCHist)                   delete ZDCHist;
  if(refMultHist)               delete refMultHist;
  if(refMultPileupHist)         delete refMultPileupHist;
  if(refMultNoPileupHist)       delete refMultNoPileupHist;
  if(RawrefMultHist)            delete RawrefMultHist;
  if(EventStat)                 delete EventStat;
  if(TrackStat)                 delete TrackStat;
  if(ZDCCoincidence)            delete ZDCCoincidence;

  if(DcaHist)                   delete DcaHist;
  if(DcaHistBTOFMatched)        delete DcaHistBTOFMatched;

  if(TOF_ZDCCoincidence)        delete TOF_ZDCCoincidence;
  if(refMult_ZDCCoincidence)    delete refMult_ZDCCoincidence;
  if(TOFMult_refMultHist)       delete TOFMult_refMultHist;
  if(TOF_BEMC)                  delete TOF_BEMC;
  if(BEMC_refMultHist)          delete BEMC_refMultHist;
  if(Vz_rankVzHist)             delete Vz_rankVzHist;
  if(TOF_VzHist)                delete TOF_VzHist;
  if(refMult_VzHist)            delete refMult_VzHist;
  if(TOF_rankVzHist)            delete TOF_rankVzHist;
  if(refMult_rankVzHist)        delete refMult_rankVzHist;


  if(TOF_refMultHist)           delete TOF_refMultHist;
  if(Vz_vpdVzHist)              delete Vz_vpdVzHist;
  if(refMult_ZDCHist)           delete refMult_ZDCHist;
  if(ZDCEastWestHist)           delete ZDCEastWestHist;

  // Histogramming
  // Event
  if(hVtxXvsY)                  delete hVtxXvsY;
  // Track
  if(hGlobalPtot)               delete hGlobalPtot;
  if(hGlobalPtotCut)            delete hGlobalPtotCut;
  if(hPrimaryPtot)              delete hPrimaryPtot;
  if(hPrimaryPtotCut)           delete hPrimaryPtotCut;
  if(hTransvMomentum)           delete hTransvMomentum;
  if(hGlobalPhiVsPt)            delete hGlobalPhiVsPt;
  if(hNSigmaProton)             delete hNSigmaProton;
  if(hNSigmaPion)               delete hNSigmaPion;
  if(hNSigmaElectron)           delete hNSigmaElectron;
  if(hNSigmaKaon)               delete hNSigmaKaon;
  if(hTofBeta)                  delete hTofBeta;

  if(Ptdist)                    delete Ptdist;
  if(EventCent)                 delete EventCent;

  if(runidvsrefmult)            delete runidvsrefmult;
  if(runidvszdcand)             delete runidvszdcand;
  if(runidvstofmult)            delete runidvstofmult;
  if(runidvstofmatched)         delete runidvstofmatched;
  if(runidvsbemcmatched)        delete runidvsbemcmatched;
  if(VzvsrefMult)               delete VzvsrefMult;
  if(DeltaVzvsrefMult)          delete DeltaVzvsrefMult;


  if(fHistCentrality)           delete fHistCentrality;
  if(fHistCentralityAfterCuts)  delete fHistCentralityAfterCuts;
  if(fHistMultiplicity)         delete fHistMultiplicity;
  if(fHistMultiplicityCorr)     delete fHistMultiplicityCorr;
  if(fHistEventPileUp)          delete fHistEventPileUp;
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

  // open output file
  if(doWriteHistos && mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE"); //"RECREATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    cout<<GetName()<<endl;

    // write histograms and output file before closing
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StChargedParticles::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");
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
