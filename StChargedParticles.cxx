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
    hVtxXvsY = new TH2D("hVtxXvsY",
    "hVtxXvsY",
    200,-10.,10.,200,-10.,10.);
    // Track
    hGlobalPtot = new TH1D("hGlobalPtot",
        "Global track momentum;p (GeV/c)",
        100, 0., 1. );
    hGlobalPtotCut = new TH1D("hGlobalPtotCut",
        "Global track momentum after cut;p (GeV/c)",
        100, 0., 1. );
    hPrimaryPtot = new TH1D("hPrimaryPtot",
        "Primary track momentum;p (GeV/c)",
        100, 0., 1. );
    hPrimaryPtotCut = new TH1D("hPrimaryPtotCut",
        "Primary track momentum after cut;p (GeV/c)",
        100, 0., 1. );
    hTransvMomentum = new TH1D("hTransvMomentum",
        "Track transverse momentum;p_{T} (GeV/c)",
        200, 0., 2.);

    for(int i=0; i<2; i++) {
      hGlobalPhiVsPt[i] = new TH2D(Form("hGlobalPhiVsPt_%d",i),
          Form("#phi vs. p_{T} for charge: %d;p_{T} (GeV/c);#phi (rad)", (i==0) ? 1 : -1),
          300, 0., 3.,
          630, -3.15, 3.15);
    }

    for(int i=0; i<10; i++) {
      Ptdist[i] = new TH1D(Form("Ptdist_%d",i),
          Form("p_{T} for centrality bin %d", i),
               100, 0, 10);
    }

    EventCent = new TH1D("EventCent","EventCent",10, -0.5, 9.5);

    hNSigmaPion = new TH1D("hNSigmaPion",
        "n#sigma(#pi);n#sigma(#pi)",
        400, -10., 10.);
    hNSigmaElectron = new TH1D("hNSigmaElectron",
        "n#sigma(e);n#sigma(e)",
        400,-10.,10.);
    hNSigmaKaon = new TH1D("hNSigmaKaon",
        "n#sigma(K);n#sigma(K)",
        400, -10., 10.);
    hNSigmaProton = new TH1D("hNSigmaProton",
        "n#sigma(p);n#sigma(p)",
        400, -10., 10.);

    // TofPidTrait
    hTofBeta = new TH1D("hTofBeta",
        "BTofPidTraits #beta;#beta",
        2000, 0., 2.);

    //QA histograms

    DcaHist = new TH1D("Dca","Dca",200,0,3.1);
    DcaHist->GetXaxis()->SetTitle("DCA");
    DcaHist->GetYaxis()->SetTitle("Counts");


    DcaHistBTOFMatched = new TH1D("DcaBTOFMatched","DcaBTOFMatched",200,0,3.1);
    DcaHistBTOFMatched->GetXaxis()->SetTitle("DCA(BTOFMatched)");
    DcaHistBTOFMatched->GetYaxis()->SetTitle("Counts");

    VzHist = new TH1D("Vz","Vz",350,-35,35);
    VzHist->GetXaxis()->SetTitle("Vz");
    VzHist->GetYaxis()->SetTitle("Counts");

    ZDCHist = new TH1D("ZDC","ZDC",1000,0,5000);
    ZDCHist->GetXaxis()->SetTitle("ZDC");
    ZDCHist->GetYaxis()->SetTitle("Counts");

    refMultHist = new TH1D("refMult","refMult",800,-0.5,799.5);
    refMultHist->GetXaxis()->SetTitle("refMult");
    refMultHist->GetYaxis()->SetTitle("Counts");


    refMultPileupHist = new TH1D("refMultPileup","refMultPileup",800,-0.5,799.5);
    refMultPileupHist->GetXaxis()->SetTitle("refMultPileup");
    refMultPileupHist->GetYaxis()->SetTitle("Counts");
    refMultNoPileupHist = new TH1D("refMultNoPileup","refMultNoPileup",800,-0.5,799.5);
    refMultNoPileupHist->GetXaxis()->SetTitle("refMultNoPileup");
    refMultNoPileupHist->GetYaxis()->SetTitle("Counts");


    RawrefMultHist = new TH1D("RawrefMult","refMult",800,-0.5,799.5);
    RawrefMultHist->GetXaxis()->SetTitle("refMult");
    RawrefMultHist->GetYaxis()->SetTitle("Counts");

    EventStat = new TH1D("EventStat","Event Stat.",20,-0.5,19.5);
    EventStat->GetXaxis()->SetTitle("CutsId");
    EventStat->GetYaxis()->SetTitle("Counts");

    TrackStat = new TH1D("TrackStat","Track Stat.",20,-0.5,19.5);
    TrackStat->GetXaxis()->SetTitle("CutsId");
    TrackStat->GetYaxis()->SetTitle("Counts");

    TOFMult_refMultHist = new TH2D("TOFMult_refMult","TOFMult v.s. refMult",700,-0.5,3499.5,700,-0.5,699.5);
    TOFMult_refMultHist->GetXaxis()->SetTitle("TOFMult");
    TOFMult_refMultHist->GetYaxis()->SetTitle("refMult");

    TOF_ZDCCoincidence = new TH2D("TOF_ZDCCoincidence","TOF v.s. ZDCCoincidence",700,-0.5,699.5,1000,0,60000);
    TOF_ZDCCoincidence->GetXaxis()->SetTitle("TOF Matched");
    TOF_ZDCCoincidence->GetYaxis()->SetTitle("ZDCCoincidence");

    refMult_ZDCCoincidence = new TH2D("refMult_ZDCCoincidence","refMult v.s. ZDCCoincidence",700,-0.5,699.5,1000,0,60000);
    refMult_ZDCCoincidence->GetXaxis()->SetTitle("refMult");
    refMult_ZDCCoincidence->GetYaxis()->SetTitle("ZDCCoincidence");


    TOF_refMultHist = new TH2D("TOF_refMult","TOF v.s. refMult",700,-0.5,699.5,700,-0.5,699.5);
    TOF_refMultHist->GetXaxis()->SetTitle("TOF Matched");
    TOF_refMultHist->GetYaxis()->SetTitle("refMult");

    //BEMC_refMultHist = new TH2D("BEMC_refMult","BEMC v.s. refMult",1000,0,3000,700,-0.5,699.5);
    BEMC_refMultHist = new TH2D("BEMC_refMult","BEMC v.s. refMult",700,-0.5,699.5,700,-0.5,699.5);
    BEMC_refMultHist->GetXaxis()->SetTitle("BEMC Matched");
    BEMC_refMultHist->GetYaxis()->SetTitle("refMult");

    TOF_BEMC = new TH2D("TOF_BEMC","TOF v.s. BEMC Matched",700,-0.5,699.5,700,-0.5,699.5);
    TOF_BEMC->GetXaxis()->SetTitle("TOF Matched");
    TOF_BEMC->GetYaxis()->SetTitle("BEMC Matched");

    refMult_VzHist= new TH2D("refMult_Vz","refMult v.s Vz",700,-0.5,699.5,350,-35,35);
    refMult_VzHist->GetXaxis()->SetTitle("refMult");
    refMult_VzHist->GetYaxis()->SetTitle("Vz");

    TOF_VzHist= new TH2D("TOF_Vz","TOF v.s Vz",700,-0.5,699.5,350,-35,35);
    TOF_VzHist->GetXaxis()->SetTitle("TOF");
    TOF_VzHist->GetYaxis()->SetTitle("Vz");

    TOF_rankVzHist= new TH2D("TOF_rankVz","TOF v.s rankVz",700,-0.5,699.5,350,-35,35);
    TOF_rankVzHist->GetXaxis()->SetTitle("TOF");
    TOF_rankVzHist->GetYaxis()->SetTitle("rankVz");

    refMult_rankVzHist= new TH2D("refMult_rankVz","refMult v.s rankVz",700,-0.5,699.5,350,-35,35);
    refMult_rankVzHist->GetXaxis()->SetTitle("refMult");
    refMult_rankVzHist->GetYaxis()->SetTitle("rankVz");

    Vz_rankVzHist= new TH2D("Vz_rankVz","Vz v.s rankVz",350,-35,35,350,-35,35);
    Vz_rankVzHist->GetXaxis()->SetTitle("Vz");
    Vz_rankVzHist->GetYaxis()->SetTitle("rankVz");

    Vz_vpdVzHist= new TH2D("Vz_vpdVz","Vz v.s vpdVz",350,-35,35,350,-35,35);
    Vz_vpdVzHist->GetXaxis()->SetTitle("Vz");
    Vz_vpdVzHist->GetYaxis()->SetTitle("vpdVz");

    refMult_ZDCHist= new TH2D("refMult_ZDC","refMult v.s. ZDC",700,-0.5,699.5,1000,0,5000);
    refMult_ZDCHist->GetXaxis()->SetTitle("refMult");
    refMult_ZDCHist->GetYaxis()->SetTitle("ZDC");

    ZDCCoincidence = new TH1D("ZDCCoincidence","ZDCCoincidence",1000,0,60000);
    ZDCCoincidence->GetXaxis()->SetTitle("ZDC X Rate");
    ZDCCoincidence->GetYaxis()->SetTitle("Counts");

    ZDCEastWestHist= new TH2D("ZDCEastWest","zdc East v.s. West",500,0,2500,500,0,2500);
    ZDCEastWestHist->GetXaxis()->SetTitle("zdcE");
    ZDCEastWestHist->GetYaxis()->SetTitle("zdcW");

    int runbins = 1428;
    int runmin = 19084005;
    int runmax = 19130031;
    //Run-by-run
    runidvsrefmult= new TProfile("runidvsrefmult","Run Id-RefMult",runbins,runmin,runmax);//,5000,0,5000);
    runidvsrefmult->GetXaxis()->SetTitle("RunId");
    runidvsrefmult->GetYaxis()->SetTitle("RefMult");
    runidvsrefmult->Sumw2();

    runidvszdcand= new TProfile("runidvszdcand","Run Id-ZDC coincidence",runbins,runmin,runmax);//,5000,0,5000);
    runidvszdcand->GetXaxis()->SetTitle("RunId");
    runidvszdcand->GetYaxis()->SetTitle("ZDCAnd");
    runidvszdcand->Sumw2();

    runidvstofmatched= new TProfile("runidvstofmatched","Run Id-TOFMatched",runbins,runmin,runmax);//,5000,0,5000);
    runidvstofmatched->GetXaxis()->SetTitle("RunId");
    runidvstofmatched->GetYaxis()->SetTitle("TOFMatched");
    runidvstofmatched->Sumw2();

    runidvsbemcmatched= new TProfile("runidvsbemcmatched","Run Id-BEMCMatched",runbins,runmin,runmax);//,5000,0,5000);
    runidvsbemcmatched->GetXaxis()->SetTitle("RunId");
    runidvsbemcmatched->GetYaxis()->SetTitle("BEMCMatched");
    runidvsbemcmatched->Sumw2();

    runidvstofmult= new TProfile("runidvstofmult","Run Id-TOFMult",runbins,runmin,runmax);//,5000,0,5000);
    runidvstofmult->GetXaxis()->SetTitle("RunId");
    runidvstofmult->GetYaxis()->SetTitle("TOFMult");
    runidvstofmult->Sumw2();

    VzvsrefMult = new TProfile("VzvsrefMult","Vz-refMult",350,-35,35);
    VzvsrefMult->GetXaxis()->SetTitle("Vz");
    VzvsrefMult->GetYaxis()->SetTitle("refMult");
    VzvsrefMult->Sumw2();

    DeltaVzvsrefMult = new TProfile("DeltaVzvsrefMult","Delta(Vz)-refMult",350,-35,35);
    DeltaVzvsrefMult->GetXaxis()->SetTitle("Delta(Vz)");
    DeltaVzvsrefMult->GetYaxis()->SetTitle("refMult");
    DeltaVzvsrefMult->Sumw2();


    fHistCentrality = new TH1F("fHistCentrality", "No. events vs centrality", nHistCentBins, 0, 100);
    fHistCentrality->Sumw2();
    fHistCentralityAfterCuts = new TH1F("fHistCentralityAfterCuts", "No. events vs centrality after cuts", nHistCentBins, 0, 100);
    fHistCentralityAfterCuts->Sumw2();
    fHistMultiplicity = new TH1F("fHistMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);
    fHistMultiplicity->Sumw2();
    fHistMultiplicityCorr = new TH1F("fHistMultiplicityCorr", "No. events vs corr. multiplicity", kHistMultBins, 0, kHistMultMax);
    fHistMultiplicityCorr->Sumw2();

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

  VzHist->Write();
  ZDCHist->Write();
  refMultHist->Write();
  refMultPileupHist->Write();
  refMultNoPileupHist->Write();
  RawrefMultHist->Write();
  EventStat->Write();
  TrackStat->Write();
  ZDCCoincidence->Write();

  DcaHist->Write();
  DcaHistBTOFMatched->Write();

  TOF_ZDCCoincidence->Write();
  refMult_ZDCCoincidence->Write();
  TOFMult_refMultHist->Write();
  TOF_BEMC->Write();
  BEMC_refMultHist->Write();
  Vz_rankVzHist->Write();
  TOF_VzHist->Write();
  refMult_VzHist->Write();
  TOF_rankVzHist->Write();
  refMult_rankVzHist->Write();


  TOF_refMultHist->Write();
  Vz_vpdVzHist->Write();
  refMult_ZDCHist->Write();
  ZDCEastWestHist->Write();

  // Histogramming
  // Event
  hVtxXvsY->Write();
  // Track
  hGlobalPtot->Write();
  hGlobalPtotCut->Write();
  hPrimaryPtot->Write();
  hPrimaryPtotCut->Write();
  hTransvMomentum->Write();
  hGlobalPhiVsPt[0]->Write();
  hGlobalPhiVsPt[1]->Write();
  hNSigmaProton->Write();
  hNSigmaPion->Write();
  hNSigmaElectron->Write();
  hNSigmaKaon->Write();
  hTofBeta->Write();

  for(int i=0; i<10; i++){ Ptdist[i]->Write();}
  EventCent->Write();

  runidvsrefmult->Write();
  runidvszdcand->Write();
  runidvstofmult->Write();
  runidvstofmatched->Write();
  runidvsbemcmatched->Write();
  VzvsrefMult->Write();
  DeltaVzvsrefMult->Write();

  fHistCentrality->Write();
  fHistCentralityAfterCuts->Write();
  fHistMultiplicity->Write();
  fHistMultiplicityCorr->Write();
  fHistEventPileUp->Write();
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
  EventStat->Fill(1);
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
