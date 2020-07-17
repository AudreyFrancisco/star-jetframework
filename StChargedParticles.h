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
class TH1F;
class TH2;
class TH2F;
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
  // switches
  Bool_t               doWriteHistos;           // write QA histos
  Bool_t               doUsePrimTracks;         // primary track switch
  Int_t                fDebugLevel;             // debug printout level
  Int_t                fRunFlag;                // Run Flag numerator value
  Bool_t               doppAnalysis;            // use pp analysis data
  Bool_t               fDoEffCorr;              // efficiency correction to tracks
  Bool_t               fCorrPileUp;             // correct for Pile-Up using the CME method (isobar parameters)

  // event cuts
  Double_t             fMaxEventTrackPt;        // max track pt in the event (to cut on)
  Bool_t               doRejectBadRuns;         // switch to reject bad runs and thus skip from analysis
  Double_t             fEventZVtxMinCut;        // min event z-vertex cut
  Double_t             fEventZVtxMaxCut;        // max event z-vertex cut
  Double_t             fEventVzDiffCut;        // max event z-vertex - z-vpd cut
  Double_t             fEventVrCut;             // max event sqrt(vx2 + vy2)vertex cut
  Int_t                fCentralitySelectionCut; // centrality selection cut
  Bool_t               fRequireCentSelection;   // require particular centrality bin

  // names
  TString              mOutName;                // name of output file
  TString              fAnalysisMakerName;      // name of this analysis maker
  TString              fTracksName;             // name of track collection


  Double_t             fTrackPtMinCut;          // min track pt cut
  Double_t             fTrackPtMaxCut;          // max track pt cut
  Double_t             fTrackPhiMinCut;         // min track phi cut
  Double_t             fTrackPhiMaxCut;         // max track phi cut
  Double_t             fTrackEtaMinCut;         // min track eta cut
  Double_t             fTrackEtaMaxCut;         // max track eta cut
  Double_t             fTrackDCAcut;            // max track dca cut
  Int_t                fTracknHitsFit;          // requirement for track hits
  Double_t             fTracknHitsRatio;        // requirement for nHitsFit / nHitsMax
  Double_t             fTracknHitsRatioMax;        // requirement for nHitsFit / nHitsMax
  Int_t                fTrackChargePos;   // track charge sign
  Int_t                fGoodTrackCounter;       // good tracks - passed quality cuts
  // centrality
  Double_t             fCentralityScaled;       // scaled by 5% centrality
  Int_t                fmycentral;       // scaled by 5% centrality
  Int_t                ref16;                   // multiplicity bin (16)
  Int_t                ref9;                    // multiplicity bin (9)

  // event
  Float_t              Bfield;                  // event Bfield
  TVector3             mVertex;                 // event vertex 3-vector
  Double_t             zVtx;                    // z-vertex component
  Int_t                fRunNumber;              // Run number

  //event
  UInt_t               fTriggerToUse;
  UInt_t               fMBEventType;

 private:
  StPicoDstMaker      *mPicoDstMaker; // PicoDstMaker object
  StPicoDst           *mPicoDst;      // PicoDst object
  StPicoEvent         *mPicoEvent;    // PicoEvent object
  StCentMaker         *mCentMaker;    // Centrality maker object
  StJetFrameworkPicoBase *mBaseMaker; // Base maker object


  // histograms
  //QA histograms
  TH1F *VzHist;//!
  TH1F *ZDCHist;//!
  TH1F *refMultHist;//!
  TH1F *refMultPileupHist;//!
  TH1F *refMultNoPileupHist;//!
  TH1F *RawrefMultHist;//!
  TH1F *EventStat;//!
  TH1F *TrackStat;//!
  TH1F *ZDCCoincidence;//!

  TH1F *DcaHist;//!
  TH1F *DcaHistBTOFMatched;//!

  TH2F *TOF_ZDCCoincidence;//!
  TH2F *refMult_ZDCCoincidence;//!
  TH2F *TOFMult_refMultHist;//!
  TH2F *TOF_BEMC;//!
  TH2F *BEMC_refMultHist;//!
  TH2F *Vz_rankVzHist;//!
  TH2F *TOF_VzHist;//!
  TH2F *refMult_VzHist;//!
  TH2F *TOF_rankVzHist;//!
  TH2F *refMult_rankVzHist;//!


  TH2F *TOF_refMultHist;//!
  TH2F *Vz_vpdVzHist;//!
  TH2F *refMult_ZDCHist;//!
  TH2F *ZDCEastWestHist;//!

  // Histogramming
  // Event
  TH2F *hVtxXvsY;//!
  // Track
  TH1F *hGlobalPtot;//!
  TH1F *hGlobalPtotCut;//!
  TH1F *hPrimaryPtot;//!
  TH1F *hPrimaryPtotCut;//!
  TH1F *hTransvMomentum;//!
  TH2F *hGlobalPhiVsPt[2];//!
  TH1F *hNSigmaProton;//!
  TH1F *hNSigmaPion;//!
  TH1F *hNSigmaElectron;//!
  TH1F *hNSigmaKaon;//!
  TH1F *hTofBeta;//!

  TH1F *Ptdist[10];//!
  TH1F *EventCent;//!

  TProfile * runidvsrefmult;//!
  TProfile * runidvszdcand;//!
  TProfile * runidvstofmult;//!
  TProfile * runidvstofmatched;//!
  TProfile * runidvsbemcmatched;//!
  TProfile * VzvsrefMult;//!
  TProfile * DeltaVzvsrefMult;//!

  TH1F     *fHistCentrality;//!
  TH1F     *fHistCentralityAfterCuts;//!
  TH1F     *fHistMultiplicity;//!
  TH1F     *fHistMultiplicityCorr;//!
  TH1F     *fHistEventPileUp;//!

  StChargedParticles(const StChargedParticles&);            // not implemented
  StChargedParticles &operator=(const StChargedParticles&); // not implemented

  ClassDef(StChargedParticles, 0) // track/cluster QA task
};
#endif
