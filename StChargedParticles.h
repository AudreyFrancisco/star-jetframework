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
  virtual void         SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }

  // track / cluster setters
  virtual void         SetTrackPtRange(Double_t ptmi, Double_t ptma) { fTrackPtMinCut = ptmi; fTrackPtMaxCut = ptma; }
  virtual void         SetTrackPhiRange(Double_t phimi, Double_t phima) { fTrackPhiMinCut = phimi; fTrackPhiMaxCut = phima; }
  virtual void         SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
  virtual void         SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
  virtual void         SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
  virtual void         SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }
  virtual void         SetTracknHitsRatioMax(Double_t r)     { fTracknHitsRatioMax = r   ; }
  virtual void         SetTrackSign(int pos)     { fTrackChargePos = pos   ; cout <<"selecting " << fTrackChargePos << "only "; }//0 for negative tracks, 1 for positive, do not call for all

  //pt histos
  virtual void	       SetHistoPtMax(Double_t pt) 	  { ptmax = pt; }
  virtual void         SetHistoNPtBins(Int_t n)		  { nptbins = n; }

  virtual void	       SetVarPtHistBinWidths(Double_t w0, Double_t w1, Double_t w2, Double_t w3, Double_t w4, Double_t w5)	{ ptW0 = w0, ptW1 = w1, ptW2 = w2, ptW3 = w3, ptW4 = w4, ptW5 = w5; }
  virtual void	       SetVarPtHistUpperLims(Double_t lim0, Double_t lim1, Double_t lim2, Double_t lim3, Double_t lim4)	{maxpt0 = lim0, maxpt1 = lim1, maxpt2 = lim2, maxpt3 = lim3, maxpt4 = lim4; }

  // event selection
  virtual void         SetTriggerToUse(UInt_t ttu)        { fTriggerToUse = ttu; }
  virtual void         SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }
  virtual void         SetPileUpCorrection(Bool_t m)      { fCorrPileUp = m; }

  // efficiency correction setter
  virtual void         SetDoEffCorr(Bool_t effcorr)     { fDoEffCorr = effcorr; }

  // common setters
  void                 SetTracksName(const char *n)     { fTracksName    = n;  }


 protected:
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

  // pt histos
  Double_t 	       ptmax;			//max pt 
  Int_t		       nptbins;			//number of pt bins
  Double_t 	       ptW0,ptW1,ptW2,ptW3,ptW4,ptW5;	//pt bin widths
  Double_t	       maxpt0,maxpt1,maxpt2,maxpt3,maxpt4; //pt bin upper limits
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

  // centrality objects
  StRefMultCorr       *grefmultCorr;






  // histograms
  //QA histograms
  TH1F *fVzHist;//!
  TH1F *fZDCHist;//!
  TH1F *frefMultHist;//!
  TH1F *frefMultPileupHist;//!
  TH1F *frefMultNoPileupHist;//!
  TH1F *fRawrefMultHist;//!
  TH1F *fEventStat;//!
  TH1F *fTrackStat;//!
  TH1F *fZDCCoincidence;//!

  TH1F *fDcaHist;//!
  TH1F *fDcaHistBTOFMatched;//!

  TH2F *fTOF_ZDCCoincidence;//!
  TH2F *frefMult_ZDCCoincidence;//!
  TH2F *fTOFMult_refMultHist;//!
  TH2F *fTOF_BEMC;//!
  TH2F *fBEMC_refMultHist;//!
  TH2F *fVz_rankVzHist;//!
  TH2F *fTOF_VzHist;//!
  TH2F *frefMult_VzHist;//!
  TH2F *fTOF_rankVzHist;//!
  TH2F *frefMult_rankVzHist;//!


  TH2F *fTOF_refMultHist;//!
  TH2F *fVz_vpdVzHist;//!
  TH2F *frefMult_ZDCHist;//!
  TH2F *fZDCEastWestHist;//!

  // Histogramming
  // Event
  TH2F *fhVtxXvsY;//!
  // Track
  TH1F *fhGlobalPtot;//!
  TH1F *fhGlobalPtotCut;//!
  TH1F *fhPrimaryPtot;//!
  TH1F *fhPrimaryPtotCut;//!
  TH1F *fhTransvMomentum;//!
  TH2F *fhGlobalPhiVsPt[2];//!
  TH1F *fhNSigmaProton;//!
  TH1F *fhNSigmaPion;//!
  TH1F *fhNSigmaElectron;//!
  TH1F *fhNSigmaKaon;//!
  TH1F *fhTofBeta;//!

  TH1F *fIntPtdist;//!
  TH1F *fIntVarPtdist;//!
  
  TH1F *fPtdist[10];//!
  TH1F *fVarPtdist[10];//!
  TH1F *fHistNTrackvsEta;//!
  TH1F *fHistNTrackvsPhi;//!
  TH1F *fEventCent;//!

  TProfile *frunidvsrefmult;//!
  TProfile *frunidvszdcand;//!
  TProfile *frunidvstofmult;//!
  TProfile *frunidvstofmatched;//!
  TProfile *frunidvsbemcmatched;//!
  TProfile *fVzvsrefMult;//!
  TProfile *fDeltaVzvsrefMult;//!

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
