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
#include <THnSparse.h>
#include "TVector3.h"

// StRoot includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoArrays.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"

// jet-framework includes
#include "StPicoConstants.h"
#include "runlistP12id.h" // Run12 pp
#include "runlistP16ij.h"
#include "runlistP17id.h" // SL17i - Run14, now SL18b (March20)
#include "runlistRun14AuAu_P18ih.h" // new Run14 AuAu
#include "runlistIso.h"
#include "StEmcPosition2.h"
#include "StJetFrameworkPicoBase.h"
#include "StCentMaker.h"

// towers/clusters related includes:
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcClusterCollection.h"
#include "StEmcCollection.h"
#include "StEmcCluster.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
//#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcRawHit.h"
#include "StEmcModule.h"
#include "StEmcDetector.h"

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
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  doppAnalysis(kFALSE),
  fCorrPileUp(kFALSE),
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
  //mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fMBEventType(2),  // kVPDMB
  fTriggerToUse(0), // kTriggerANY
  mMuDstMaker(0x0),
  mMuDst(0x0),
  mMuInputEvent(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mCentMaker(0x0),
  mBaseMaker(0x0),
  grefmultCorr(0x0)
{
  // Default constructor.
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
  //mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fMBEventType(2),   // kVPDMB
  fTriggerToUse(0),  // kTriggerANY
  mMuDstMaker(0x0),
  mMuDst(0x0),
  mMuInputEvent(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mCentMaker(0x0),
  mBaseMaker(0x0),
  grefmultCorr(0x0)
{
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
  if(VzHist)         delete VzHist;
  if(ZDCHist)         delete ZDCHist;
  if(refMultHist)         delete refMultHist;
  if(refMultPileupHist)         delete refMultPileupHist;
  if(refMultNoPileupHist)         delete refMultNoPileupHist;
  if(RawrefMultHist)         delete RawrefMultHist;
  if(EventStat)         delete EventStat;
  if(TrackStat)         delete TrackStat;
  if(ZDCCoincidence)         delete ZDCCoincidence;

  if(DcaHist)         delete DcaHist;
  if(DcaHistBTOFMatched)         delete DcaHistBTOFMatched;

  if(TOF_ZDCCoincidence)         delete TOF_ZDCCoincidence;
  if(refMult_ZDCCoincidence)         delete refMult_ZDCCoincidence;
  if(TOFMult_refMultHist)         delete TOFMult_refMultHist;
  if(TOF_BEMC)         delete TOF_BEMC;
  if(BEMC_refMultHist)         delete BEMC_refMultHist;
  if(Vz_rankVzHist)         delete Vz_rankVzHist;
  if(TOF_VzHist)         delete TOF_VzHist;
  if(refMult_VzHist)         delete refMult_VzHist;
  if(TOF_rankVzHist)         delete TOF_rankVzHist;
  if(refMult_rankVzHist)         delete refMult_rankVzHist;


  if(TOF_refMultHist)         delete TOF_refMultHist;
  if(Vz_vpdVzHist)         delete Vz_vpdVzHist;
  if(refMult_ZDCHist)         delete refMult_ZDCHist;
  if(ZDCEastWestHist)         delete ZDCEastWestHist;

  // Histogramming
  // Event
  if(hVtxXvsY)         delete hVtxXvsY;
  // Track
  if(hGlobalPtot)         delete hGlobalPtot;
  if(hGlobalPtotCut)         delete hGlobalPtotCut;
  if(hPrimaryPtot)         delete hPrimaryPtot;
  if(hPrimaryPtotCut)         delete hPrimaryPtotCut;
  if(hTransvMomentum)         delete hTransvMomentum;
  if(hGlobalPhiVsPt)         delete hGlobalPhiVsPt;
  if(hNSigmaProton)         delete hNSigmaProton;
  if(hNSigmaPion)         delete hNSigmaPion;
  if(hNSigmaElectron)         delete hNSigmaElectron;
  if(hNSigmaKaon)         delete hNSigmaKaon;
  if(hTofBeta)         delete hTofBeta;

  if(Ptdist)         delete Ptdist;
  if(EventCent)         delete EventCent;

  if(runidvsrefmult)         delete runidvsrefmult;
  if(runidvszdcand)         delete runidvszdcand;
  if(runidvstofmult)         delete runidvstofmult;
  if(runidvstofmatched)         delete runidvstofmatched;
  if(runidvsbemcmatched)         delete runidvsbemcmatched;
  if(VzvsrefMult)         delete VzvsrefMult;
  if(DeltaVzvsrefMult)         delete DeltaVzvsrefMult;


  if(fHistCentrality)         delete fHistCentrality;
  if(fHistCentralityAfterCuts)         delete fHistCentralityAfterCuts;
  if(fHistMultiplicity)         delete fHistMultiplicity;
  if(fHistMultiplicityCorr)         delete fHistMultiplicityCorr;
  if(fHistEventPileUp)		delete fHistEventPileUp;
}
//
//_____________________________________________________________________________
Int_t StChargedParticles::Init() {
  // declare histograms
  DeclareHistograms();

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

}
//
// Main loop, called for each event.
//________________________________________________________________________
int StChargedParticles::Make()
{
  EventStat->Fill(1);
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;


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

  // get vertex 3-vector and declare variables
  mVertex = mPicoEvent->primaryVertex(); //pRcVx
  zVtx = mVertex.z(); //zTpc

  float yVtx = mVertex.y();
  float xVtx = mVertex.x();
  float vzVPD = mPicoEvent->vzVpd();
  float vDiff = TMath::Abs(zVtx - vzVPD);
  float vrVtx = TMath::Sqrt(xVtx*xVtx+yVtx*yVtx);

  // track variables
  int ntracks = mPicoDst->numberOfTracks();
  int ntracksCovM = mPicoDst->numberOfTrackCovMatrices();

  EventStat->Fill(2);

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  //if(fDebugLevel == kDebugEmcTrigger)
  //cout<<"EventTriggers: ";
  double triggerids [4] = {600001, 600011, 600021, 600031};
  bool IsTrigger=false;
  for(int j=0; j<4; j++){
    IsTrigger = (IsTrigger || (mPicoEvent->isTrigger(triggerids[j]))==1);
  }
  if(!IsTrigger){ return kStOk;}
  EventStat->Fill(3);


  int nBtofMatch =  mPicoEvent->nBTOFMatch();
  if (nBtofMatch < 1 ) { return kStOk; }


// ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    LOG_WARN << " No CenttMaker! Skip! " << endm;
    return kStWarn;
  }

  // centrality variables
  int grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  int refMult =  mCentMaker->GetrefMult();  // see StPicoEvent
  ref9 = mCentMaker->GetRef9();   // binning from central -> peripheral
  ref16 = mCentMaker->GetRef16(); // binning from central -> peripheral
  int cent16 = mCentMaker->GetCent16(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int cent9 = mCentMaker->GetCent9(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%
  //if(cent16 == -1 && fDebugLevel != 99) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them CHECK TODO
  fmycentral = 8-cent9; // WARNING!! RefMultCorr convention 0->peripheral, 8-> central; so we do 8-centrlaity
  // fill histograms
  fHistCentrality->Fill(fCentralityScaled);
  fHistMultiplicity->Fill(grefMult);


  if(fabs(zVtx) < fEventZVtxMaxCut) fHistMultiplicityCorr->Fill(refCorr2);

  if( fmycentral<0 || fmycentral>8 ) return kStOK;

  EventStat->Fill(4);


  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
  }

  EventStat->Fill(5);
  if(fabs(zVtx)<70) EventStat->Fill(6);
  if(fabs(zVtx)<50) EventStat->Fill(7);


  // cut event on max track pt > 30.0 GeV
  // if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField();

  VzvsrefMult->Fill(zVtx,grefMult);
  DeltaVzvsrefMult->Fill(vDiff,grefMult);

  // Z-vertex cut - per the Aj analysis (-40, 40)
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;
  EventStat->Fill(8);
  if(vrVtx > fEventVrCut) return kStOk;
  EventStat->Fill(9);
  if(vDiff > fEventVzDiffCut) return kStOk;
  EventStat->Fill(10);

  int RunId_Order =  GetRunNo(fRunNumber);
  bool doTrackQA = kTRUE;
  if(fCorrPileUp){
    refMultHist->Fill(grefMult);
    if(mCentMaker->Refmult_check(nBtofMatch,refCorr2,3,4)) {
      refMultNoPileupHist->Fill(grefMult);
      doTrackQA = kTRUE;
    }
    else{
      refMultPileupHist->Fill(grefMult);
      fHistEventPileUp->Fill(RunId_Order + 1., 1);
      doTrackQA = kFALSE;
      return kStOk;
   }
  }
  EventStat->Fill(11);

  runidvsrefmult->Fill(mPicoEvent->runId(),grefMult);
  runidvstofmult->Fill(mPicoEvent->runId(),mPicoEvent->btofTrayMultiplicity());
  runidvstofmatched->Fill(mPicoEvent->runId(),nBtofMatch);
  runidvsbemcmatched->Fill(mPicoEvent->runId(),mPicoEvent->nBEMCMatch());
  runidvszdcand->Fill(mPicoEvent->runId(),mPicoEvent->ZDCx());
  VzHist->Fill(zVtx);
  hVtxXvsY->Fill(xVtx,yVtx);

  ZDCHist->Fill(mPicoEvent->ZdcSumAdcEast()+mPicoEvent->ZdcSumAdcWest());
  ZDCEastWestHist->Fill(mPicoEvent->ZdcSumAdcEast(),mPicoEvent->ZdcSumAdcWest());
  TOFMult_refMultHist->Fill(mPicoEvent->btofTrayMultiplicity(),grefMult);
  TOF_refMultHist->Fill(nBtofMatch,grefMult);

  TOF_ZDCCoincidence->Fill(nBtofMatch,mPicoEvent->ZDCx());
  refMult_ZDCCoincidence->Fill(grefMult,mPicoEvent->ZDCx());

  TOF_BEMC->Fill(nBtofMatch,mPicoEvent->nBEMCMatch());
  BEMC_refMultHist->Fill(mPicoEvent->nBEMCMatch(),grefMult);
  TOF_VzHist->Fill(nBtofMatch,zVtx);
  TOF_rankVzHist->Fill(nBtofMatch,mPicoEvent->ranking());

  refMult_VzHist->Fill(grefMult,zVtx);
  refMult_rankVzHist->Fill(grefMult,mPicoEvent->ranking());

  Vz_vpdVzHist->Fill(zVtx,vzVPD);
  Vz_rankVzHist->Fill(zVtx,mPicoEvent->ranking());
  refMult_ZDCHist->Fill(grefMult,mPicoEvent->ZdcSumAdcEast()+mPicoEvent->ZdcSumAdcWest());
  ZDCCoincidence->Fill(mPicoEvent->ZDCx());


  ///Add by YU, to check the Event per CENT, Jun.14
  EventCent->Fill(fmycentral);
     // Track loop
  for(Int_t iTrk=0; iTrk<ntracks; iTrk++) {

    // Retrieve i-th pico track
    StPicoTrack *gTrack = mPicoDst->track(iTrk);

    if(!gTrack) continue;
    //std::cout << "Track #[" <</ (iTrk+1) << "/" << NoGlobalTracks << "]"  << std::endl;
   //------------------------------------------------------------------------
      // Grigory's default histograms QA block for TPC tracks
    //________________________________________________________________________
    hGlobalPtot->Fill( gTrack->gMom().Mag() );
      if( gTrack->isPrimary() ) {
      hPrimaryPtot->Fill( gTrack->pMom().Mag() );
    }
     // Simple single-track cut
      if( gTrack->gMom().Mag() < 0.1 ||
          //          gTrack->gDCA(pRcVx).Mag()>50. ) {
        TMath::Abs(gTrack->gDCAxy(xVtx, yVtx))>50. ){

          continue;
        }

      hGlobalPtotCut->Fill( gTrack->gMom().Mag() );
      if( gTrack->isPrimary() ) {
        hPrimaryPtotCut->Fill( gTrack->pMom().Mag() );
      }
      if( gTrack->charge() > 0 ) {
        hGlobalPhiVsPt[0]->Fill( gTrack->gMom().Pt(),
            gTrack->gMom().Phi() );
      }
      else {
        hGlobalPhiVsPt[1]->Fill( gTrack->gMom().Pt(),
            gTrack->gMom().Phi() );
      }
      hNSigmaElectron->Fill( gTrack->nSigmaElectron() );
      hNSigmaPion->Fill( gTrack->nSigmaPion() );
      hNSigmaKaon->Fill( gTrack->nSigmaKaon() );
      hNSigmaProton->Fill( gTrack->nSigmaProton() );

      hTransvMomentum->Fill( gTrack->gMom().Pt() );

      // Check if track has TOF signal
      /*if( gTrack->isTofTrack() ) {
        // Retrieve corresponding trait
        StPicoBTofPidTraits *trait = mPicoDst->btofPidTraits( gTrack->bTofPidTraitsIndex() );
        if( !trait ) {
          std::cout << "O-oh... No BTofPidTrait # " << gTrack->bTofPidTraitsIndex()
            << " for track # " << iTrk << std::endl;
          std::cout << "Check that you turned on the branch!" << std::endl;
          continue;
        }
        // Fill beta
        hTofBeta->Fill( trait->btofBeta() );
      } //if( isTofTrack() )
      */
      //________________________________________________________________________
      //
      //Grigory's block for TPC tracks ends here
      //________________________________________________________________________
      //


      DcaHist->Fill(gTrack->gDCA(mVertex).Mag());
      if(gTrack->isTofTrack())DcaHistBTOFMatched->Fill(gTrack->gDCA(mVertex).Mag());

     //------------------------------------------------------------------------
      // Prithwish's analysis cut block
      //________________________________________________________________________
      //if(fabs(zTpc)>20.) continue;
      //if(sqrt(pRcVx.x()*pRcVx.x()+yVtx*yVtx) >2.0) continue;
      //if(fabs(zTpc-zVpd)>2) continue;
      //if(nBTOFMatch<10) continue;
      TrackStat->Fill(1);
      if(!gTrack->isPrimary()) continue;
      TrackStat->Fill(2);
      if(gTrack->nHitsFit()<=fTracknHitsFit) continue;
      TrackStat->Fill(3);
      //if(TMath::Abs(gTrack->gDCA(pRcVx).Mag())>3) continue; //We can also use this cut //PT Jan15, 2018
      //if(TMath::Abs(gTrack->gDCAxy(pRcVx.x(), pRcVx.y())) > MAX_gDCAxy) continue;
      if(TMath::Abs(gTrack->gDCA(mVertex).Mag())>fTrackDCAcut) continue; //We can also use this cut //PT Jan15, 2018
      TrackStat->Fill(4);

      ///Add by Yu, just to check the Pt SPECTRUM
      if(fabs(gTrack->pMom().PseudoRapidity())<=fTrackEtaMaxCut){
        Ptdist[fmycentral]->Fill(gTrack->gMom().Pt()); ///------------ AUDREY
      }
    }

  return kStOK;
}
//
//
//
// Function: Track Quality Cuts

//
// Returns pt of hardest track in the event
//______________________________________________________________________________________________
Double_t StChargedParticles::GetMaxTrackPt()
{
  // get # of tracks
  int nTrack = mPicoDst->numberOfTracks();
  double fMaxTrackPt = -99;

  // loop over all tracks
  for(int i = 0; i < nTrack; i++) {
    // get track pointer
    StPicoTrack *track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) continue;

    // apply standard track cuts - (can apply more restrictive cuts below)
    //if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = track->pMom();
    } else {
      // get global track vector
      mTrkMom = track->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();

    // get max track
    if(pt > fMaxTrackPt) { fMaxTrackPt = pt; }
  }

  return fMaxTrackPt;
}
//
//________________________________________________________________________
void StChargedParticles::FillTriggerIDs(TH1 *h) {
  // All non-test triggers for Run12 Run14

  // Run14 AuAu (200 GeV) - 51, 0-50
  unsigned int triggersRun14[] = {440001, 440004, 440005, 440006, 440007, 440015, 440016, 440017, 440050, 440061, 440064, 450005, 450008, 450009, 450010, 450011, 450012, 450013, 450014, 450015, 450018, 450020, 450021, 450023, 450024, 450025, 450050, 450060, 450103, 450201, 450202, 450203, 450211, 450212, 450213, 450600, 450601, 460001, 460002, 460003, 460005, 460007, 460012, 460101, 460102, 460111, 460201, 460202, 460203, 460212, 490016};

  // Run12 pp (200 GeV) - 27, 0-26
  unsigned int triggersRun12[] = {370001, 370011, 370021, 370022, 370031, 370032, 370301, 370341, 370361, 370501, 370511, 370521, 370522, 370531, 370541, 370542, 370546, 370601, 370611, 370621, 370641, 370701, 370801, 370980, 370981, 370982, 370983};
  unsigned int triggersRunIso[] = {600001, 600011, 600021, 600031, 600201, 600202, 600203, 600211, 600212, 600213, 600221, 600222, 600231, 600232};

  // get size of trigger ID arrays:
  size_t nRun12IDs = sizeof(triggersRun12)/sizeof(triggersRun12[0]);
  size_t nRun14IDs = sizeof(triggersRun14)/sizeof(triggersRun14[0]);
  size_t nRunIsoIDs = sizeof(triggersRunIso)/sizeof(triggersRunIso[0]);
  int nLoopMax = 0;
  if(StJetFrameworkPicoBase::Run12_pp200)   nLoopMax = nRun12IDs;
  if(StJetFrameworkPicoBase::Run14_AuAu200) nLoopMax = nRun14IDs;
  if(StJetFrameworkPicoBase::RunIsobar) nLoopMax = nRunIsoIDs;

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  for(unsigned int i = 0; i < mytriggers.size(); i++) {

    // check for valid, non-test trigger ID
    if(mytriggers[i] > 1000) {
      for(int j = 0; j < nLoopMax; j++) {
        if(mytriggers[i] == triggersRun12[j] && fRunFlag == StJetFrameworkPicoBase::Run12_pp200)   h->Fill(j + 1);
        if(mytriggers[i] == triggersRun14[j] && fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) h->Fill(j + 1);
        if(mytriggers[i] == triggersRunIso[j] && fRunFlag == StJetFrameworkPicoBase::RunIsobar) h->Fill(j + 1);

      } // loops over ID's
    }   // non-test trigger
  }     // loop over triggers

  // label bins of the analysis trigger selection summary
  for(int i = 0; i < nLoopMax; i++) {
    if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200)   h->GetXaxis()->SetBinLabel(i+1, Form("%i", triggersRun12[i]));
    if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) h->GetXaxis()->SetBinLabel(i+1, Form("%i", triggersRun14[i]));
    if(fRunFlag == StJetFrameworkPicoBase::RunIsobar) h->GetXaxis()->SetBinLabel(i+1, Form("%i", triggersRun14[i]));
  }

  // set x-axis labels vertically
  h->LabelsOption("v");

}
//
// function to do some event QA
//________________________________________________________________________________________
// FIXME TODO - this function needs to be updated when ADDING runs
// _________________________________________________________________________________
Int_t StChargedParticles::GetRunNo(int runid){
  // Run12 pp (200 GeV) - 857
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200) {
    for(int i = 0; i < 857; i++) {
      if(Run12pp_IdNo[i] == runid) {
        return i;
      }
    }
  }

  // Run14 AuAu (200 GeV) - new picoDst production is 830
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
    for(int i = 0; i < 830; i++) {
      if(Run14AuAu_P18ih_IdNo[i] == runid) {
        return i;
      }
    }
  }

  // Run16 AuAu (200 GeV) -  1359 for Run16 AuAu
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
    for(int i = 0; i < 1359; i++){
      if(Run16AuAu_IdNo[i] == runid) {
        return i;
      }
    }
  }

  // Run18 Isobar (200 GeV)
  if(fRunFlag == StJetFrameworkPicoBase::RunIsobar) {
    // 405 runs on May 23
    for(int i = 0; i < 1428; i++){ //FIXME
      if(RunIso_IdNo[i] == runid) {
        return i;
      }
    }
  }
  cout<<" *********** RunID not matched with list ************!!!! "<<endl;
  return -999;
}
