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
  if(fhGlobalPhiVsPt0)            delete fhGlobalPhiVsPt0;
  if(fhGlobalPhiVsPt1)            delete fhGlobalPhiVsPt1;
  if(fhNSigmaProton)             delete fhNSigmaProton;
  if(fhNSigmaPion)               delete fhNSigmaPion;
  if(fhNSigmaElectron)           delete fhNSigmaElectron;
  if(fhNSigmaKaon)               delete fhNSigmaKaon;
  if(fhTofBeta)                  delete fhTofBeta;

  if(fPtdist[0])                    delete fPtdist[0];
  if(fPtdist[1])                    delete fPtdist[1];
  if(fPtdist[2])                    delete fPtdist[2];
  if(fPtdist[3])                    delete fPtdist[3];
  if(fPtdist[4])                    delete fPtdist[4];
  if(fPtdist[5])                    delete fPtdist[5];
  if(fPtdist[6])                    delete fPtdist[6];
  if(fPtdist[7])                    delete fPtdist[7];
  if(fPtdist[8])                    delete fPtdist[8];
  if(fPtdist[9])                    delete fPtdist[9];
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
  // declare histograms
  DeclareHistograms();
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
    fhGlobalPtotCut = new TH1F("hGlobalPtotCut", "Global track momentum after cut;p (GeV/c)", 100, 0., 1. );
    fhPrimaryPtot = new TH1F("hPrimaryPtot", "Primary track momentum;p (GeV/c)", 100, 0., 1. );
    fhPrimaryPtotCut = new TH1F("hPrimaryPtotCut", "Primary track momentum after cut;p (GeV/c)", 100, 0., 1. );
    fhTransvMomentum = new TH1F("hTransvMomentum", "Track transverse momentum;p_{T} (GeV/c)", 200, 0., 2.);


    fhGlobalPhiVsPt[0] = new TH2F("hGlobalPhiVsPt0", "#phi vs. p_{T} for charge: 1;p_{T} (GeV/c);#phi (rad)",
          300, 0., 3., 630, -3.15, 3.15);
    fhGlobalPhiVsPt[1] = new TH2F("hGlobalPhiVsPt1", "#phi vs. p_{T} for charge: -1;p_{T} (GeV/c);#phi (rad)",
          300, 0., 3., 630, -3.15, 3.15);


    fPtdist[0]= new TH1F("Ptdist0", "p_{T} for centrality bin 0",100, 0, 10);
    fPtdist[1]= new TH1F("Ptdist1", "p_{T} for centrality bin 1",100, 0, 10);
    fPtdist[2]= new TH1F("Ptdist2", "p_{T} for centrality bin 2",100, 0, 10);
    fPtdist[3]= new TH1F("Ptdist3", "p_{T} for centrality bin 3",100, 0, 10);
    fPtdist[4]= new TH1F("Ptdist4", "p_{T} for centrality bin 4",100, 0, 10);
    fPtdist[5]= new TH1F("Ptdist5", "p_{T} for centrality bin 5",100, 0, 10);
    fPtdist[6]= new TH1F("Ptdist6", "p_{T} for centrality bin 6",100, 0, 10);
    fPtdist[7]= new TH1F("Ptdist7", "p_{T} for centrality bin 7",100, 0, 10);
    fPtdist[8]= new TH1F("Ptdist8", "p_{T} for centrality bin 8",100, 0, 10);
    fPtdist[9]= new TH1F("Ptdist9", "p_{T} for centrality bin 9",100, 0, 10);

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

    fHistCentrality = new TH1F("fHistCentrality", "No. events vs centrality", nHistCentBins, 0, 100);
    fHistCentrality->Sumw2();
    fHistCentralityAfterCuts = new TH1F("fHistCentralityAfterCuts", "No. events vs centrality after cuts", nHistCentBins, 0, 100);
    fHistCentralityAfterCuts->Sumw2();
    fHistMultiplicity = new TH1F("fHistMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);
    fHistMultiplicity->Sumw2();
    fHistMultiplicityCorr = new TH1F("fHistMultiplicityCorr", "No. events vs corr. multiplicity", kHistMultBins, 0, kHistMultMax);
    fHistMultiplicityCorr->Sumw2();

    fHistEventPileUp = new TH1F("fHistEventPileUp", "PileUp events", nRunBins, 0.5, nRunBinsMax);
    fHistEventPileUp->Sumw2();

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

  fPtdist[0]->Write();
  fPtdist[1]->Write();
  fPtdist[2]->Write();
  fPtdist[3]->Write();
  fPtdist[4]->Write();
  fPtdist[5]->Write();
  fPtdist[6]->Write();
  fPtdist[7]->Write();
  fPtdist[8]->Write();
  fPtdist[9]->Write();
  fEventCent->Write();

  frunidvsrefmult->Write();
  frunidvszdcand->Write();
  frunidvstofmult->Write();
  frunidvstofmatched->Write();
  frunidvsbemcmatched->Write();
  fVzvsrefMult->Write();
  fDeltaVzvsrefMult->Write();

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
//________________________________________________________________________
