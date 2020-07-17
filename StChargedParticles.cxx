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
  fDoEffCorr(kFALSE),
  fCorrPileUp(kFALSE),
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
  // mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fTriggerToUse(0), // kTriggerANY
  fMBEventType(2),  // kVPDMB
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mCentMaker(0x0),
  mBaseMaker(0x0),
  grefmultCorr(0x0)
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
  fDoEffCorr(kFALSE),
  fCorrPileUp(kFALSE),
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
  // mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fTriggerToUse(0),  // kTriggerANY
  fMBEventType(2),   // kVPDMB
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mCentMaker(0x0),
  mBaseMaker(0x0),
  grefmultCorr(0x0)
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
  if(fVzHist)                    delete fVzHist;
  if(fZDCHist)                   delete fZDCHist;
  if(frefMultHist)               delete frefMultHist;
  if(frefMultPileupHist)         delete frefMultPileupHist;
  if(frefMultNoPileupHist)       delete frefMultNoPileupHist;
  if(fRawrefMultHist)            delete fRawrefMultHist;
  if(fEventStat)                 delete fEventStat;
  if(fTrackStat)                 delete fTrackStat;
  if(fZDCCoincidence)            delete fZDCCoincidence;

  if(fDcaHist)                   delete fDcaHist;
  if(fDcaHistBTOFMatched)        delete fDcaHistBTOFMatched;

  if(fTOF_ZDCCoincidence)        delete fTOF_ZDCCoincidence;
  if(frefMult_ZDCCoincidence)    delete frefMult_ZDCCoincidence;
  if(fTOFMult_refMultHist)       delete fTOFMult_refMultHist;
  if(fTOF_BEMC)                  delete fTOF_BEMC;
  if(fBEMC_refMultHist)          delete fBEMC_refMultHist;
  if(fVz_rankVzHist)             delete fVz_rankVzHist;
  if(fTOF_VzHist)                delete fTOF_VzHist;
  if(frefMult_VzHist)            delete frefMult_VzHist;
  if(fTOF_rankVzHist)            delete fTOF_rankVzHist;
  if(frefMult_rankVzHist)        delete frefMult_rankVzHist;


  if(fTOF_refMultHist)           delete fTOF_refMultHist;
  if(fVz_vpdVzHist)              delete fVz_vpdVzHist;
  if(frefMult_ZDCHist)           delete frefMult_ZDCHist;
  if(fZDCEastWestHist)           delete fZDCEastWestHist;

  // Histogramming
  // Event
  if(fhVtxXvsY)                  delete fhVtxXvsY;
  // Track
  if(fhGlobalPtot)               delete fhGlobalPtot;
  if(fhGlobalPtotCut)            delete fhGlobalPtotCut;
  if(fhPrimaryPtot)              delete fhPrimaryPtot;
  if(fhPrimaryPtotCut)           delete fhPrimaryPtotCut;
  if(fhTransvMomentum)           delete fhTransvMomentum;
  if(fhGlobalPhiVsPt)            delete fhGlobalPhiVsPt;
  if(fhNSigmaProton)             delete fhNSigmaProton;
  if(fhNSigmaPion)               delete fhNSigmaPion;
  if(fhNSigmaElectron)           delete fhNSigmaElectron;
  if(fhNSigmaKaon)               delete fhNSigmaKaon;
  if(fhTofBeta)                  delete fhTofBeta;

  if(fPtdist)                    delete fPtdist;
  if(fEventCent)                 delete fEventCent;

  if(frunidvsrefmult)            delete frunidvsrefmult;
  if(frunidvszdcand)             delete frunidvszdcand;
  if(frunidvstofmult)            delete frunidvstofmult;
  if(frunidvstofmatched)         delete frunidvstofmatched;
  if(frunidvsbemcmatched)        delete frunidvsbemcmatched;
  if(fVzvsrefMult)               delete fVzvsrefMult;
  if(fDeltaVzvsrefMult)          delete fDeltaVzvsrefMult;


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

  return kStOK;
}
//
// Declare histograms and global objects for running
//________________________________________________________________________
void StChargedParticles::DeclareHistograms() {
    // declare histograms
    double pi = 1.0*TMath::Pi();

    // set binning for run based corrections - run dependent
    Int_t nRunBins = 1; // - just a default
    if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200)   nRunBins = 857 + 43;
    if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) nRunBins = 830; //1654;
    if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) nRunBins = 1359;
    if(fRunFlag == StJetFrameworkPicoBase::RunIsobar) nRunBins = 1428; //FIXME
    Double_t nRunBinsMax = (Double_t)nRunBins + 0.5;

    // tweak refmult plot binnings for pp datasets
    int nFactor = (doppAnalysis) ? 5 : 1;

    // binning for cent histograms
    int nHistCentBins = 20;

    // binning for mult histograms - pp : AuAu
    double kHistMultMax = (doppAnalysis) ? 100. : 800.;
    int kHistMultBins = (doppAnalysis) ? 100. : 400.;

    // basic event QA'
    //Event
    fhVtxXvsY = new TH2F("hVtxXvsY", "hVtxXvsY",200,-10.,10.,200,-10.,10.);
    // Track
    fhGlobalPtot = new TH1F("hGlobalPtot", "Global track momentum;p (GeV/c)", 100, 0., 1. );
    hGlobalPtotCut = new TH1F("hGlobalPtotCut", "Global track momentum after cut;p (GeV/c)", 100, 0., 1. );
    hPrimaryPtot = new TH1F("hPrimaryPtot", "Primary track momentum;p (GeV/c)", 100, 0., 1. );
    hPrimaryPtotCut = new TH1F("hPrimaryPtotCut", "Primary track momentum after cut;p (GeV/c)", 100, 0., 1. );
    hTransvMomentum = new TH1F("hTransvMomentum", "Track transverse momentum;p_{T} (GeV/c)", 200, 0., 2.);

    for(int i=0; i<2; i++) {
      fhGlobalPhiVsPt[i] = new TH2F(Form("hGlobalPhiVsPt_%d",i), Form("#phi vs. p_{T} for charge: %d;p_{T} (GeV/c);#phi (rad)", (i==0) ? 1 : -1),
          300, 0., 3., 630, -3.15, 3.15);
    }

    for(int i=0; i<10; i++) {
      fPtdist[i] = new TH1F(Form("Ptdist_%d",i), Form("p_{T} for centrality bin %d", i),
               100, 0, 10);
    }

    fEventCent = new TH1F("EventCent","EventCent",10, -0.5, 9.5);

    fhNSigmaPion = new TH1F("hNSigmaPion", "n#sigma(#pi);n#sigma(#pi)",400, -10., 10.);
    fhNSigmaElectron = new TH1F("hNSigmaElectron","n#sigma(e);n#sigma(e)",400,-10.,10.);
    fhNSigmaKaon = new TH1F("hNSigmaKaon","n#sigma(K);n#sigma(K)",400, -10., 10.);
    fhNSigmaProton = new TH1F("hNSigmaProton","n#sigma(p);n#sigma(p)",400, -10., 10.);

    // TofPidTrait
    fhTofBeta = new TH1F("hTofBeta","BTofPidTraits #beta;#beta",2000, 0., 2.);

    //QA histograms

    fDcaHist = new TH1F("Dca","Dca",200,0,3.1);
    fDcaHist->GetXaxis()->SetTitle("DCA");
    fDcaHist->GetYaxis()->SetTitle("Counts");


    fDcaHistBTOFMatched = new TH1F("DcaBTOFMatched","DcaBTOFMatched",200,0,3.1);
    fDcaHistBTOFMatched->GetXaxis()->SetTitle("DCA(BTOFMatched)");
    fDcaHistBTOFMatched->GetYaxis()->SetTitle("Counts");

    fVzHist = new TH1F("Vz","Vz",350,-35,35);
    fVzHist->GetXaxis()->SetTitle("Vz");
    fVzHist->GetYaxis()->SetTitle("Counts");

    fZDCHist = new TH1F("ZDC","ZDC",1000,0,5000);
    fZDCHist->GetXaxis()->SetTitle("ZDC");
    fZDCHist->GetYaxis()->SetTitle("Counts");

    frefMultHist = new TH1F("refMult","refMult",800,-0.5,799.5);
    frefMultHist->GetXaxis()->SetTitle("refMult");
    frefMultHist->GetYaxis()->SetTitle("Counts");


    frefMultPileupHist = new TH1F("refMultPileup","refMultPileup",800,-0.5,799.5);
    frefMultPileupHist->GetXaxis()->SetTitle("refMultPileup");
    frefMultPileupHist->GetYaxis()->SetTitle("Counts");
    frefMultNoPileupHist = new TH1F("refMultNoPileup","refMultNoPileup",800,-0.5,799.5);
    frefMultNoPileupHist->GetXaxis()->SetTitle("refMultNoPileup");
    frefMultNoPileupHist->GetYaxis()->SetTitle("Counts");


    fRawrefMultHist = new TH1F("RawrefMult","refMult",800,-0.5,799.5);
    fRawrefMultHist->GetXaxis()->SetTitle("refMult");
    fRawrefMultHist->GetYaxis()->SetTitle("Counts");

    fEventStat = new TH1F("EventStat","Event Stat.",20,-0.5,19.5);
    fEventStat->GetXaxis()->SetTitle("CutsId");
    fEventStat->GetYaxis()->SetTitle("Counts");

    fTrackStat = new TH1F("TrackStat","Track Stat.",20,-0.5,19.5);
    fTrackStat->GetXaxis()->SetTitle("CutsId");
    fTrackStat->GetYaxis()->SetTitle("Counts");

    fTOFMult_refMultHist = new TH2F("TOFMult_refMult","TOFMult v.s. refMult",700,-0.5,3499.5,700,-0.5,699.5);
    fTOFMult_refMultHist->GetXaxis()->SetTitle("TOFMult");
    fTOFMult_refMultHist->GetYaxis()->SetTitle("refMult");

    fTOF_ZDCCoincidence = new TH2F("TOF_ZDCCoincidence","TOF v.s. ZDCCoincidence",700,-0.5,699.5,1000,0,60000);
    fTOF_ZDCCoincidence->GetXaxis()->SetTitle("TOF Matched");
    fTOF_ZDCCoincidence->GetYaxis()->SetTitle("ZDCCoincidence");

    frefMult_ZDCCoincidence = new TH2F("refMult_ZDCCoincidence","refMult v.s. ZDCCoincidence",700,-0.5,699.5,1000,0,60000);
    frefMult_ZDCCoincidence->GetXaxis()->SetTitle("refMult");
    frefMult_ZDCCoincidence->GetYaxis()->SetTitle("ZDCCoincidence");


    fTOF_refMultHist = new TH2F("TOF_refMult","TOF v.s. refMult",700,-0.5,699.5,700,-0.5,699.5);
    fTOF_refMultHist->GetXaxis()->SetTitle("TOF Matched");
    fTOF_refMultHist->GetYaxis()->SetTitle("refMult");

    //BEMC_refMultHist = new TH2F("BEMC_refMult","BEMC v.s. refMult",1000,0,3000,700,-0.5,699.5);
    fBEMC_refMultHist = new TH2F("BEMC_refMult","BEMC v.s. refMult",700,-0.5,699.5,700,-0.5,699.5);
    fBEMC_refMultHist->GetXaxis()->SetTitle("BEMC Matched");
    fBEMC_refMultHist->GetYaxis()->SetTitle("refMult");

    fTOF_BEMC = new TH2F("TOF_BEMC","TOF v.s. BEMC Matched",700,-0.5,699.5,700,-0.5,699.5);
    fTOF_BEMC->GetXaxis()->SetTitle("TOF Matched");
    fTOF_BEMC->GetYaxis()->SetTitle("BEMC Matched");

    frefMult_VzHist= new TH2F("refMult_Vz","refMult v.s Vz",700,-0.5,699.5,350,-35,35);
    frefMult_VzHist->GetXaxis()->SetTitle("refMult");
    frefMult_VzHist->GetYaxis()->SetTitle("Vz");

    fTOF_VzHist= new TH2F("TOF_Vz","TOF v.s Vz",700,-0.5,699.5,350,-35,35);
    fTOF_VzHist->GetXaxis()->SetTitle("TOF");
    fTOF_VzHist->GetYaxis()->SetTitle("Vz");

    fTOF_rankVzHist= new TH2F("TOF_rankVz","TOF v.s rankVz",700,-0.5,699.5,350,-35,35);
    fTOF_rankVzHist->GetXaxis()->SetTitle("TOF");
    fTOF_rankVzHist->GetYaxis()->SetTitle("rankVz");

    frefMult_rankVzHist= new TH2F("refMult_rankVz","refMult v.s rankVz",700,-0.5,699.5,350,-35,35);
    frefMult_rankVzHist->GetXaxis()->SetTitle("refMult");
    frefMult_rankVzHist->GetYaxis()->SetTitle("rankVz");

    fVz_rankVzHist= new TH2F("Vz_rankVz","Vz v.s rankVz",350,-35,35,350,-35,35);
    fVz_rankVzHist->GetXaxis()->SetTitle("Vz");
    fVz_rankVzHist->GetYaxis()->SetTitle("rankVz");

    fVz_vpdVzHist= new TH2F("Vz_vpdVz","Vz v.s vpdVz",350,-35,35,350,-35,35);
    fVz_vpdVzHist->GetXaxis()->SetTitle("Vz");
    fVz_vpdVzHist->GetYaxis()->SetTitle("vpdVz");

    frefMult_ZDCHist= new TH2F("refMult_ZDC","refMult v.s. ZDC",700,-0.5,699.5,1000,0,5000);
    frefMult_ZDCHist->GetXaxis()->SetTitle("refMult");
    frefMult_ZDCHist->GetYaxis()->SetTitle("ZDC");

    fZDCCoincidence = new TH1F("ZDCCoincidence","ZDCCoincidence",1000,0,60000);
    fZDCCoincidence->GetXaxis()->SetTitle("ZDC X Rate");
    fZDCCoincidence->GetYaxis()->SetTitle("Counts");

    fZDCEastWestHist= new TH2F("ZDCEastWest","zdc East v.s. West",500,0,2500,500,0,2500);
    fZDCEastWestHist->GetXaxis()->SetTitle("zdcE");
    fZDCEastWestHist->GetYaxis()->SetTitle("zdcW");

    int runbins = 1428;
    int runmin = 19084005;
    int runmax = 19130031;
    //Run-by-run
    frunidvsrefmult= new TProfile("runidvsrefmult","Run Id-RefMult",runbins,runmin,runmax);//,5000,0,5000);
    frunidvsrefmult->GetXaxis()->SetTitle("RunId");
    frunidvsrefmult->GetYaxis()->SetTitle("RefMult");
    frunidvsrefmult->Sumw2();

    frunidvszdcand= new TProfile("runidvszdcand","Run Id-ZDC coincidence",runbins,runmin,runmax);//,5000,0,5000);
    frunidvszdcand->GetXaxis()->SetTitle("RunId");
    frunidvszdcand->GetYaxis()->SetTitle("ZDCAnd");
    frunidvszdcand->Sumw2();

    frunidvstofmatched= new TProfile("runidvstofmatched","Run Id-TOFMatched",runbins,runmin,runmax);//,5000,0,5000);
    frunidvstofmatched->GetXaxis()->SetTitle("RunId");
    frunidvstofmatched->GetYaxis()->SetTitle("TOFMatched");
    frunidvstofmatched->Sumw2();

    frunidvsbemcmatched= new TProfile("runidvsbemcmatched","Run Id-BEMCMatched",runbins,runmin,runmax);//,5000,0,5000);
    frunidvsbemcmatched->GetXaxis()->SetTitle("RunId");
    frunidvsbemcmatched->GetYaxis()->SetTitle("BEMCMatched");
    frunidvsbemcmatched->Sumw2();

    frunidvstofmult= new TProfile("runidvstofmult","Run Id-TOFMult",runbins,runmin,runmax);//,5000,0,5000);
    frunidvstofmult->GetXaxis()->SetTitle("RunId");
    frunidvstofmult->GetYaxis()->SetTitle("TOFMult");
    frunidvstofmult->Sumw2();

    fVzvsrefMult = new TProfile("VzvsrefMult","Vz-refMult",350,-35,35);
    fVzvsrefMult->GetXaxis()->SetTitle("Vz");
    fVzvsrefMult->GetYaxis()->SetTitle("refMult");
    fVzvsrefMult->Sumw2();

    fDeltaVzvsrefMult = new TProfile("DeltaVzvsrefMult","Delta(Vz)-refMult",350,-35,35);
    fDeltaVzvsrefMult->GetXaxis()->SetTitle("Delta(Vz)");
    fDeltaVzvsrefMult->GetYaxis()->SetTitle("refMult");
    fDeltaVzvsrefMult->Sumw2();

    ffHistCentrality = new TH1F("fHistCentrality", "No. events vs centrality", nHistCentBins, 0, 100);
    ffHistCentrality->Sumw2();
    ffHistCentralityAfterCuts = new TH1F("fHistCentralityAfterCuts", "No. events vs centrality after cuts", nHistCentBins, 0, 100);
    ffHistCentralityAfterCuts->Sumw2();
    ffHistMultiplicity = new TH1F("fHistMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);
    ffHistMultiplicity->Sumw2();
    ffHistMultiplicityCorr = new TH1F("fHistMultiplicityCorr", "No. events vs corr. multiplicity", kHistMultBins, 0, kHistMultMax);
    ffHistMultiplicityCorr->Sumw2();

    // run range for runID histogram
    int nRunBinSize = 200;
    double runMin = 0, runMax = 0;
    if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200)   { runMin = 13000000.; runMax = 13100000.; }
    if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) { runMin = 15050000.; runMax = 15200000.; }
    if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) { runMin = 17050000.; runMax = 17150000.; }
    if(fRunFlag == StJetFrameworkPicoBase::RunIsobar) { runMin = 19084005.; runMax = 19130031.; nRunBinSize = 1; } //FIXME

    // Switch on Sumw2 for all histos - (except profiles)
    //SetSumw2();
}
//
// write histograms
//________________________________________________________________________
void StChargedParticles::WriteHistograms() {

  fVzHist->Write();
  fZDCHist->Write();
  frefMultHist->Write();
  frefMultPileupHist->Write();
  frefMultNoPileupHist->Write();
  fRawrefMultHist->Write();
  fEventStat->Write();
  fTrackStat->Write();
  fZDCCoincidence->Write();

  fDcaHist->Write();
  fDcaHistBTOFMatched->Write();

  fTOF_ZDCCoincidence->Write();
  frefMult_ZDCCoincidence->Write();
  fTOFMult_refMultHist->Write();
  fTOF_BEMC->Write();
  fBEMC_refMultHist->Write();
  fVz_rankVzHist->Write();
  fTOF_VzHist->Write();
  frefMult_VzHist->Write();
  fTOF_rankVzHist->Write();
  frefMult_rankVzHist->Write();


  fTOF_refMultHist->Write();
  fVz_vpdVzHist->Write();
  frefMult_ZDCHist->Write();
  fZDCEastWestHist->Write();

  // Histogramming
  // Event
  fhVtxXvsY->Write();
  // Track
  fhGlobalPtot->Write();
  fhGlobalPtotCut->Write();
  fhPrimaryPtot->Write();
  fhPrimaryPtotCut->Write();
  fhTransvMomentum->Write();
  fhGlobalPhiVsPt[0]->Write();
  fhGlobalPhiVsPt[1]->Write();
  fhNSigmaProton->Write();
  fhNSigmaPion->Write();
  fhNSigmaElectron->Write();
  fhNSigmaKaon->Write();
  fhTofBeta->Write();

  for(int i=0; i<10; i++){ fPtdist[i]->Write();}
  fEventCent->Write();

  frunidvsrefmult->Write();
  frunidvszdcand->Write();
  frunidvstofmult->Write();
  frunidvstofmatched->Write();
  frunidvsbemcmatched->Write();
  fVzvsrefMult->Write();
  fDeltaVzvsrefMult->Write();

  ffHistCentrality->Write();
  ffHistCentralityAfterCuts->Write();
  ffHistMultiplicity->Write();
  ffHistMultiplicityCorr->Write();
  ffHistEventPileUp->Write();
}
//
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
  fEventStat->Fill(1);
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

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
