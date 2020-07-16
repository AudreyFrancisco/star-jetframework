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
#include "StEmcUtil/geometry/StEmcGeom.h"

// ROOT classes
class TClonesArray;
class TObjArray;
class TList;
class TH1;
class TH1F;
class TH2;
class TH2F;
class THnSparse;
class TProfile;
class TVector3;

// STAR classes
class StMaker;
class StChain;
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoBTowHit;

// EMC and tower related classes
class StEmcGeom;
class StBemcTables; //v3.14
class StEmcCluster;
class StEmcCollection;

// star jet-frameworks classes
class StJetFrameworkPicoBase;
class StEmcPosition2;
class StCentMaker;

// centrality class
class StRefMultCorr;


class StPicoTrackClusterQA : public StMaker {
//class StPicoTrackClusterQA : public StJetFrameworkPicoBase {
 public:

  // debug flags for specifics
  enum fDebugFlagEnum {
    kDebugNothing, // don't want lowest elements to be used
    kDebugEmcTrigger,
    kDebugGeneralEvt,
    kDebugCentrality,
  };

  StPicoTrackClusterQA();
  StPicoTrackClusterQA(const char *name, bool dohistos, const char* outName);
  virtual ~StPicoTrackClusterQA();

  // needed class functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t *opt="");
  virtual Int_t Finish();

  // booking of histograms (optional)
  void    DeclareHistograms();
  void    WriteHistograms();

  // THnSparse Setup
  virtual THnSparse*      NewTHnSparseFTracks(const char* name, UInt_t entries);
  virtual void GetDimParamsTracks(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  virtual THnSparse*      NewTHnSparseFTowers(const char* name, UInt_t entries);
  virtual void GetDimParamsTowers(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

  // switches
  virtual void         SetUsePrimaryTracks(Bool_t P)    { doUsePrimTracks       = P; }
  virtual void         SetDebugLevel(Int_t l)           { fDebugLevel           = l; }
  virtual void         SetRunFlag(Int_t f)              { fRunFlag              = f; }
  virtual void         SetdoppAnalysis(Bool_t pp)       { doppAnalysis          = pp;}
  virtual void         SetTurnOnCentSelection(Bool_t o) { fRequireCentSelection = o; }
  virtual void         SetCentralityBinCut(Int_t c)     { fCentralitySelectionCut = c; }

  // event setters
  virtual void         SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
  virtual void         SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
  virtual void         SetMaxEventTowerEt(Double_t mxEt)  { fMaxEventTowerEt = mxEt; }
  virtual void         SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }

  // track / cluster setters
  virtual void         SetTrackPtRange(Double_t ptmi, Double_t ptma) { fTrackPtMinCut = ptmi; fTrackPtMaxCut = ptma; }
  virtual void         SetTrackPhiRange(Double_t phimi, Double_t phima) { fTrackPhiMinCut = phimi; fTrackPhiMaxCut = phima; }
  virtual void         SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
  virtual void         SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
  virtual void         SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
  virtual void         SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }
  virtual void         SetTracknHitsRatioMax(Double_t r)     { fTracknHitsRatioMax = r   ; }
  virtual void         SetTrackSign(Bool_t pos)     { fTrackChargePos = pos   ; cout <<"selecting " << fTrackChargePos << "only "; }//0 for negative tracks, 1 for positive, do not call for all
  virtual void         SetClusterPtRange(Double_t mi, Double_t ma) { fClusterPtMinCut = mi; fClusterPtMaxCut = ma; }

  // tower setters
  virtual void         SetTowerERange(Double_t enmi, Double_t enmx) { fTowerEMinCut = enmi; fTowerEMaxCut = enmx; }
  virtual void         SetTowerEtaRange(Double_t temi, Double_t temx) { fTowerEtaMinCut = temi; fTowerEtaMaxCut = temx; }
  virtual void         SetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fTowerPhiMinCut = tpmi; fTowerPhiMaxCut = tpmx; }

  // event selection
  virtual void         SetTriggerToUse(UInt_t ttu)        { fTriggerToUse = ttu; }
  virtual void         SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
  virtual void         SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }
  virtual void         SetDoTowerQAforHT(Bool_t m)        { fDoTowerQAforHT = m; }
  virtual void         SetPileUpCorrection(Bool_t m)      { fCorrPileUp = m; }

  // efficiency correction setter
  virtual void         SetDoEffCorr(Bool_t effcorr)     { fDoEffCorr = effcorr; }

  // common setters
  void                 SetClusName(const char *n)       { fCaloName      = n;  }
  void                 SetTracksName(const char *n)     { fTracksName    = n;  }

  /* define if tower status should be used to reject towers, or if all
   * towers should be accepted - default is to accept all towers, then
   * generate a bad tower list for the entire data set.
  */

 protected:
  Double_t             GetMaxTrackPt();               // find max track pt in event
  void                 FillTriggerIDs(TH1* h);
  Int_t                GetRunNo(int runid);

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
  TString              fCaloName;               // name of calo cluster collection

  Double_t             fTrackPtMinCut;          // min track pt cut
  Double_t             fTrackPtMaxCut;          // max track pt cut
  Double_t             fClusterPtMinCut;        // min cluster pt cut
  Double_t             fClusterPtMaxCut;        // max cluster pt cut
  Double_t             fTrackPhiMinCut;         // min track phi cut
  Double_t             fTrackPhiMaxCut;         // max track phi cut
  Double_t             fTrackEtaMinCut;         // min track eta cut
  Double_t             fTrackEtaMaxCut;         // max track eta cut
  Double_t             fTrackDCAcut;            // max track dca cut
  Int_t                fTracknHitsFit;          // requirement for track hits
  Double_t             fTracknHitsRatio;        // requirement for nHitsFit / nHitsMax
  Double_t             fTracknHitsRatioMax;        // requirement for nHitsFit / nHitsMax
  Int_t		             fTrackChargePos;		// track charge sign
  Int_t                fGoodTrackCounter;       // good tracks - passed quality cuts
  // centrality
  Double_t             fCentralityScaled;       // scaled by 5% centrality
  Int_t                ref16;                   // multiplicity bin (16)
  Int_t                ref9;                    // multiplicity bin (9)

  // event
  Float_t              Bfield;                  // event Bfield
  TVector3             mVertex;                 // event vertex 3-vector
  Double_t             zVtx;                    // z-vertex component
  Int_t                fRunNumber;              // Run number

 private:
  StMuDstMaker        *mMuDstMaker;   // MuDstMaker object
  StMuDst             *mMuDst;        // muDst object
  StMuEvent           *mMuInputEvent; // muDst event object
  StPicoDstMaker      *mPicoDstMaker; // PicoDstMaker object
  StPicoDst           *mPicoDst;      // PicoDst object
  StPicoEvent         *mPicoEvent;    // PicoEvent object
  StCentMaker         *mCentMaker;    // Centrality maker object
  StJetFrameworkPicoBase *mBaseMaker; // Base maker object


  //bool              *mTowerStatusArr; // tower status array

  // centrality objects
  StRefMultCorr       *grefmultCorr;

  // histograms
  //QA histograms
  TH1D *VzHist;//!
  TH1D *ZDCHist;//!
  TH1D *refMultHist;//!
  TH1D *refMultPileupHist;//!
  TH1D *refMultNoPileupHist;//!
  TH1D *RawrefMultHist;//!
  TH1D *EventStat;//!
  TH1D *TrackStat;//!
  TH1D *ZDCCoincidence;//!

  TH1D *DcaHist;//!
  TH1D *DcaHistBTOFMatched;//!

  TH2D *TOF_ZDCCoincidence;//!
  TH2D *refMult_ZDCCoincidence;//!
  TH2D *TOFMult_refMultHist;//!
  TH2D *TOF_BEMC;//!
  TH2D *BEMC_refMultHist;//!
  TH2D *Vz_rankVzHist;//!
  TH2D *TOF_VzHist;//!
  TH2D *refMult_VzHist;//!
  TH2D *TOF_rankVzHist;//!
  TH2D *refMult_rankVzHist;//!


  TH2D *TOF_refMultHist;//!
  TH2D *Vz_vpdVzHist;//!
  TH2D *refMult_ZDCHist;//!
  TH2D *ZDCEastWestHist;//!

  // Histogramming
  // Event
  TH2D *hVtxXvsY;//!
  // Track
  TH1D *hGlobalPtot;//!
  TH1D *hGlobalPtotCut;//!
  TH1D *hPrimaryPtot;//!
  TH1D *hPrimaryPtotCut;//!
  TH1D *hTransvMomentum;//!
  TH2D *hGlobalPhiVsPt[2];//!
  TH1D *hNSigmaProton;//!
  TH1D *hNSigmaPion;//!
  TH1D *hNSigmaElectron;//!
  TH1D *hNSigmaKaon;//!
  TH1D *hTofBeta;//!

  TH1D *Ptdist[10];//!
  TH1D *EventCent;//!

  TProfile * runidvsrefmult;//!
  TProfile * runidvszdcand;//!
  TProfile * runidvstofmult;//!
  TProfile * runidvstofmatched;//!
  TProfile * runidvsbemcmatched;//!
  TProfile * VzvsrefMult;//!
  TProfile * DeltaVzvsrefMult;//!

  TH1F           *fHistCentrality;//!
  TH1F           *fHistCentralityAfterCuts;//!
  TH1F           *fHistMultiplicity;//!
  TH1F           *fHistMultiplicityCorr;//!
  // bad and dead tower list
  std::set<Int_t>        badTowers;
  std::set<Int_t>        deadTowers;

  // bad run list
  void                   ResetBadRunList( );
  Bool_t                 AddBadRuns(TString csvfile);
  Bool_t                 IsRunOK( Int_t mRunId );
  std::set<Int_t>        badRuns;

  StChargedParticles(const StChargedParticles&);            // not implemented
  StChargedParticles &operator=(const StChargedParticles&); // not implemented

  ClassDef(StChargedParticles, 1) // track/cluster QA task
};
#endif
