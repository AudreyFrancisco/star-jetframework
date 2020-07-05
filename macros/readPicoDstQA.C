// ************************************** //
// Sept28, 2017
// some current notes:
// Run16 AuAu 200 GeV:  P16ij on RCAS has no EmcalTriggers and no clusters
//     ::P16ij has triggers that need to be accessed a different way
//
//

#include <TSystem>

// basic STAR classes
class StMemStat;
class StMaker;
class StChain;
class StPicoDstMaker;
class StRefMultCorr;
// jet-framework STAR classes
class StJetMakerTask;
class StJetFrameworkPicoBase;
class StMyAnalysisMaker;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// constants
const double pi = 1.0*TMath::Pi();
Double_t ZVtxMin = -30.0;
Double_t ZVtxMax = 30.0;

bool doCentSelection = kFALSE; //kTRUE;
bool dopp = kFALSE; // FIXME kTRUE for pp data
int RunYear = 99;   // FIXME
// kTRUE for local tests, kFALSE for job submission
bool doTEST = kFALSE;  //FIXME FIXME!!!! be aware before submission
bool doPileUP = kFALSE;
StChain *chain;

//__________________________________________________________________________________________________________
void readPicoDstQA(const Char_t *inputFile="Run_15164046_files.list", const Char_t *outputFile="towerQA.root", Int_t nEv = 10, const Char_t *fEPoutJobappend="_this_is_a_test")
{
        Int_t nEvents = 100000;
//        Int_t nEvents = 10000;
//        Int_t nEvents = 50000;
//        Int_t nEvents = 100000;        
//        if(nEv > 100) nEvents = 100000000;

        // run enumerators
        enum RunType_t {
          mJobSubmission = 0, // 0th spot - default unless set
          mRun07 =  7, mRun08 = 8,  mRun09 = 9, mRun10 = 10,
          mRun11 = 11, mRun12 = 12, mRun13 = 13,
          mRun14 = 14,
          mRun15 = 15, mRun16 = 16, mRun17 = 17, mRun18 = 18,
          mRunIso = 99
	};

        // Load necessary libraries and macros
        LoadLibs();
        LoadMacros();

        // =============================================================================== //
        // =============================================================================== //
        // =============================================================================== //
        // over-ride functions
        if(dopp) {
          doCentSelection = kFALSE;  // can't ask for a particular centrality if requesting pp collisions
        }

        // input file for tests (based on Run) - updated for new Runs as needed
        if((RunYear == mRun12) && doTEST && dopp) inputFile = "testLIST_Run12pp.list";
        if((RunYear == mRun14) && doTEST)         inputFile = "Run14_P18ih_SL20a_15110029.list";//"Run_15164046_files.list"; //"Run_15151042_files.list"; //"testLIST_Run14.list";
        if((RunYear == mRun16) && doTEST)         inputFile = "test_run17124003_files.list";
        if((RunYear == mRun17) && doTEST && dopp) inputFile = "Run17pp_510GeV.list"; // "filelist_pp2017.list";
        if((RunYear == mRunIso) && doTEST) inputFile = "runlistP2.list";//"testIsoP2.list";//RunIso.list";
	cout<<"inputFileName = "<<inputFile<<endl;

        // centrality global flags - no centrality for pp collisions
        Int_t CentralitySelection = StJetFrameworkPicoBase::kCent005;
        Int_t CentralitySelection2 = StJetFrameworkPicoBase::kCent6080;
        //Int_t CentralityDefinition = StJetFrameworkPicoBase::kgrefmult;
        Int_t CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi;
        if(RunYear == mRun12) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30;  // no centrality defintion for Run 12 pp 
        //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30; // Run14 P17id (NEW - from Nick Oct 23)
        if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi; // Run14 P18ih (NEW - from Nick June10, 2019)
        if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id;           // Run16
        if(RunYear == mRunIso) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult;           // Run16
        //if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_VpdMBnoVtx;
        cout<<"Centrality definition: "<<CentralityDefinition<<endl;

        // Run/Event Flag
        Int_t RunFlag;
        if(RunYear == mRun14) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;
        if(RunYear == mRun16) RunFlag = StJetFrameworkPicoBase::Run16_AuAu200;
        if(RunYear == mRunIso) RunFlag = StJetFrameworkPicoBase::RunIsobar;
        if(RunYear == mRun11 && dopp) RunFlag = StJetFrameworkPicoBase::Run11_pp500;
        if(RunYear == mRun12 && dopp) RunFlag = StJetFrameworkPicoBase::Run12_pp200;
        if(RunYear == mRun13 && dopp) RunFlag = StJetFrameworkPicoBase::Run13_pp510;
        if(RunYear == mRun17 && dopp) RunFlag = StJetFrameworkPicoBase::Run17_pp510;
        Bool_t RejectBadRuns = kFALSE; // switch to load and than omit bad runs
        Int_t fBadRunListVers = StJetFrameworkPicoBase::fBadRuns_w_missing_HT; // fBadRuns_w_missing_HT, fBadRuns_wo_missing_HT,

        // trigger flags - update default
        Int_t EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT1; // kIsHT1 or kIsHT2 or kIsHT3
        if(RunYear == mRun12) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2;
        if(RunYear == mRun14) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // kIsHT2 Run14
        if(RunYear == mRun16) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT1; // kIsHT1 Run16
        if(RunYear == mRun17) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT3;
        Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB; //kVPDMB5 is default
        if(RunYear == mRun12) MBEventType = StJetFrameworkPicoBase::kRun12main; // default for Run12 pp
        if(RunYear == mRun17) MBEventType = StJetFrameworkPicoBase::kVPDMB; // default for Run17 pp
        //Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerHT;     // kTriggerANY, kTriggerMB, kTriggerHT
        Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerANY;    // kTriggerANY, kTriggerMB, kTriggerHT     FIXME
        Int_t TowerListToUse = 136; // Run14 136-122: jet-hadron, early jet shape // been using 51,  3, 79   - doesn't matter for charged jets
        if(dopp) TowerListToUse = 169;
        // Run12: 1 - Raghav's list, 102 - my initial list, 169 - new list
        // Run14: 136 - main list (updated version for AuAu 200 GeV Run14), 122 - past used list
        // Run14 P18ih: 999 (initial) 

        // track flags
        bool usePrimaryTracks;
        if(RunYear == mRun12) usePrimaryTracks = kTRUE;
        if(RunYear == mRun14) usePrimaryTracks = kTRUE;  // = kTRUE for Run14, kFALSE for Run16
        if(RunYear == mRun16) usePrimaryTracks = kFALSE;
        if(RunYear == mRun17) usePrimaryTracks = kTRUE;
        if(RunYear == mRunIso) usePrimaryTracks = kTRUE;

        // update settings for new centrality definitions
        if(CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30 ||
           CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30 ||
           CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi
        ) { ZVtxMin = -28.0; ZVtxMax = 28.0; }

        // =============================================================================== //
        // =============================================================================== //
        // =============================================================================== //

        // open and close output .root file (so it exist and can be updated by Analysis Tasks)
        TFile *fout = new TFile(outputFile, "RECREATE");
        //fout->cd();
        fout->Close();
	
        // create chain
        StChain* chain = new StChain();

	// create the picoMaker maker
	//StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
        StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst"); // updated Aug6th
        picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));

        // create base class maker pointer
        StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase("baseClassMaker");
        baseMaker->SetRunFlag(RunFlag);                  // run flag (year)
        baseMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        baseMaker->SetBadRunListVers(fBadRunListVers);          // switch to select specific bad run version file
        baseMaker->SetBadTowerListVers(TowerListToUse);
        cout<<baseMaker->GetName()<<endl;  // print name of class instance

        // update the below when running analysis - minjet pt and bias requirement
        StCentMaker *CentMaker = new StCentMaker("CentMaker", picoMaker, outputFile, kFALSE);
        CentMaker->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
        CentMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
        CentMaker->SetRunFlag(RunFlag);                         // Run Flag
        CentMaker->SetdoppAnalysis(dopp);                       // pp-analysis switch
        CentMaker->SetCentralityDef(CentralityDefinition);      // centrality definition
        CentMaker->SetUseBBCCoincidenceRate(kFALSE);            // BBC or ZDC (default) rate used?
        CentMaker->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        CentMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
	cout<<CentMaker->GetName()<<endl;  // print name of class instance

        // =======================================================================================================
	StPicoTrackClusterQA *Task = new StPicoTrackClusterQA("TrackClusterQAtest", kTRUE, outputFile);
        Task->SetTrackPtRange(0., 100.0);
        Task->SetTrackPhiRange(0., 2.0*pi);
        Task->SetTrackEtaRange(-0.5, 0.5);
        Task->SetEventZVtxRange(-30., 30.);      // can be tighter for Run16 (-20,20)
        Task->SetClusterPtRange(0., 500.0);
        Task->SetTowerERange(0., 500.0);
        Task->SetUsePrimaryTracks(usePrimaryTracks);
        Task->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        Task->SetRunFlag(RunFlag);                      // RunFlag
        //Task->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task->SetTurnOnCentSelection(kFALSE);  // run analysis for specific centrality
        Task->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task->SetDebugLevel(8);
        //Task->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task->SetHadronicCorrFrac(1.0);
        Task->SetDoTowerQAforHT(kFALSE);
        Task->SetdoppAnalysis(dopp);
        Task->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        //Task->SetPileUpCorrection(doPileUP);
	//Task->SetMaxEventTrackPt(100.0); // TEST
        
        StPicoTrackClusterQA *Task0 = new StPicoTrackClusterQA("TrackClusterQA010", kTRUE, outputFile);
        Task0->SetTrackPtRange(0.2, 30.0);
        Task0->SetTrackPhiRange(0., 2.0*pi);
        Task0->SetTrackEtaRange(-0.5,0.5);
        Task0->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task0->SetClusterPtRange(0.2, 100.0);
        Task0->SetTowerERange(0.2, 100.0);
        Task0->SetUsePrimaryTracks(usePrimaryTracks);
        Task0->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        Task0->SetRunFlag(RunFlag);                      // RunFlag
        //Task0->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task0->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        Task0->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task->SetDebugLevel(8);
        //Task->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task0->SetHadronicCorrFrac(1.0);
        Task0->SetDoTowerQAforHT(kFALSE);
        Task0->SetdoppAnalysis(dopp);
        Task0->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        //Task0->SetMaxEventTrackPt(100.0); // TEST
        //Task0->SetPileUpCorrection(doPileUP);

        StPicoTrackClusterQA *Task6 = new StPicoTrackClusterQA("TrackClusterQA6080", kTRUE, outputFile);
        Task6->SetTrackPtRange(0.2, 30.0);
        Task6->SetTrackPhiRange(0., 2.0*pi);
        Task6->SetTrackEtaRange(-0.5, 0.5);
        Task6->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task6->SetClusterPtRange(0.2, 100.0);
        Task6->SetTowerERange(0.2, 100.0);
        Task6->SetUsePrimaryTracks(usePrimaryTracks);
        Task6->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        Task6->SetRunFlag(RunFlag);                      // RunFlag
        //Task6->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task6->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        Task6->SetCentralityBinCut(CentralitySelection2); // specific centrality range to run
        //Task->SetDebugLevel(8);
        //Task->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task6->SetHadronicCorrFrac(1.0);
        Task6->SetDoTowerQAforHT(kFALSE);
        Task6->SetdoppAnalysis(dopp);
        Task6->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        //Task6->SetMaxEventTrackPt(100.0); // TEST
        //Task6->SetPileUpCorrection(doPileUP);

	/*StPicoTrackClusterQA *Taskp = new StPicoTrackClusterQA("TrackClusterQAp", kTRUE, outputFile);
	Taskp->SetTrackSign(1);
        Taskp->SetTrackPtRange(0.2, 30.0);
        Taskp->SetTrackPhiRange(0., 2.0*pi);
        Taskp->SetTrackEtaRange(-1.0, 1.0);
        Taskp->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Taskp->SetClusterPtRange(0.2, 100.0);
        Taskp->SetTowerERange(0.2, 100.0);
        Taskp->SetUsePrimaryTracks(usePrimaryTracks);
        Taskp->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        Taskp->SetRunFlag(RunFlag);                      // RunFlag
        Taskp->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Taskp->SetTurnOnCentSelection(kFALSE);  // run analysis for specific centrality
        Taskp->SetCentralityBinCut(CentralitySelection2); // specific centrality range to run
        //Task->SetDebugLevel(8);
        //Task->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Taskp->SetHadronicCorrFrac(1.0);
        Taskp->SetDoTowerQAforHT(kFALSE);
        Taskp->SetdoppAnalysis(dopp);
        Taskp->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Taskp->SetMaxEventTrackPt(100.0); // TEST

        
	StPicoTrackClusterQA *Task2p = new StPicoTrackClusterQA("TrackClusterQA010p", kTRUE, outputFile);
        Task2p->SetTrackPtRange(0.2, 30.0);
	Task2p->SetTrackSign(1);
        Task2p->SetTrackPhiRange(0., 2.0*pi);
        Task2p->SetTrackEtaRange(-1.0, 1.0);
        Task2p->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task2p->SetClusterPtRange(0.2, 100.0);
        Task2p->SetTowerERange(0.2, 100.0);
        Task2p->SetUsePrimaryTracks(usePrimaryTracks);
        Task2p->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        Task2p->SetRunFlag(RunFlag);                      // RunFlag
        Task2p->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task2p->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        Task2p->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task->SetDebugLevel(8);
        //Task->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task2p->SetHadronicCorrFrac(1.0);
        Task2p->SetDoTowerQAforHT(kFALSE);
        Task2p->SetdoppAnalysis(dopp);
        Task2p->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Task2p->SetMaxEventTrackPt(100.0); // TEST

        
        StPicoTrackClusterQA *TaskPp = new StPicoTrackClusterQA("TrackClusterQA6080p", kTRUE, outputFile);
	TaskPp->SetTrackSign(1);
        TaskPp->SetTrackPtRange(0.2, 30.0);
        TaskPp->SetTrackPhiRange(0., 2.0*pi);
        TaskPp->SetTrackEtaRange(-1.0, 1.0);
        TaskPp->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        TaskPp->SetClusterPtRange(0.2, 100.0);
        TaskPp->SetTowerERange(0.2, 100.0);
        TaskPp->SetUsePrimaryTracks(usePrimaryTracks);
        TaskPp->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        TaskPp->SetRunFlag(RunFlag);                      // RunFlag
        TaskPp->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        TaskPp->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        TaskPp->SetCentralityBinCut(CentralitySelection2); // specific centrality range to run
        //Task->SetDebugLevel(8);
        //Task->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        TaskPp->SetHadronicCorrFrac(1.0);
        TaskPp->SetDoTowerQAforHT(kFALSE);
        TaskPp->SetdoppAnalysis(dopp);
        TaskPp->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        TaskPp->SetMaxEventTrackPt(100.0); // TEST


        StPicoTrackClusterQA *Taskn = new StPicoTrackClusterQA("TrackClusterQAn", kTRUE, outputFile);
	Taskn->SetTrackSign(0);	
        Taskn->SetTrackPtRange(0.2, 30.0);
        Taskn->SetTrackPhiRange(0., 2.0*pi);
        Taskn->SetTrackEtaRange(-1.0, 1.0);
        Taskn->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Taskn->SetClusterPtRange(0.2, 100.0);
        Taskn->SetTowerERange(0.2, 100.0);
        Taskn->SetUsePrimaryTracks(usePrimaryTracks);
        Taskn->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        Taskn->SetRunFlag(RunFlag);                      // RunFlag
        Taskn->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Taskn->SetTurnOnCentSelection(kFALSE);  // run analysis for specific centrality
        Taskn->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task->SetDebugLevel(8);
        //Task->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Taskn->SetHadronicCorrFrac(1.0);
        Taskn->SetDoTowerQAforHT(kFALSE);
        Taskn->SetdoppAnalysis(dopp);
        Taskn->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Taskn->SetMaxEventTrackPt(100.0); // TEST
        
	StPicoTrackClusterQA *TaskCn = new StPicoTrackClusterQA("TrackClusterQA010n", kTRUE, outputFile);
	TaskCn->SetTrackSign(0);	
        TaskCn->SetTrackPtRange(0.2, 30.0);
        TaskCn->SetTrackPhiRange(0., 2.0*pi);
        TaskCn->SetTrackEtaRange(-1.0, 1.0);
        TaskCn->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        TaskCn->SetClusterPtRange(0.2, 100.0);
        TaskCn->SetTowerERange(0.2, 100.0);
        TaskCn->SetUsePrimaryTracks(usePrimaryTracks);
        TaskCn->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        TaskCn->SetRunFlag(RunFlag);                      // RunFlag
        TaskCn->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        TaskCn->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        TaskCn->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //TaCsk->SetDebugLevel(8);
        //TaCsk->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        TaskCn->SetHadronicCorrFrac(1.0);
        TaskCn->SetDoTowerQAforHT(kFALSE);
        TaskCn->SetdoppAnalysis(dopp);
        TaskCn->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        TaskCn->SetMaxEventTrackPt(100.0); // TEST
        
        
	StPicoTrackClusterQA *Task3 = new StPicoTrackClusterQA("TrackClusterQA6080n", kTRUE, outputFile);
	Task3->SetTrackSign(0);
        Task3->SetTrackPtRange(0.2, 30.0);
        Task3->SetTrackPhiRange(0., 2.0*pi);
        Task3->SetTrackEtaRange(-1.0, 1.0);
        Task3->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task3->SetClusterPtRange(0.2, 100.0);
        Task3->SetTowerERange(0.2, 100.0);
        Task3->SetUsePrimaryTracks(usePrimaryTracks);
        Task3->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        Task3->SetRunFlag(RunFlag);                      // RunFlag
        Task3->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task3->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        Task3->SetCentralityBinCut(CentralitySelection2); // specific centrality range to run
        //Task->SetDebugLevel(8);
        //Task->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task3->SetHadronicCorrFrac(1.0);
        Task3->SetDoTowerQAforHT(kFALSE);
        Task3->SetdoppAnalysis(dopp);
        Task3->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Task3->SetMaxEventTrackPt(100.0); // TEST
*/
        // =======================================================================================================
        // QA task - MB events, little to no cuts
        StPicoTrackClusterQA *TaskA = new StPicoTrackClusterQA("TrackClusterQAnocuts", kTRUE, outputFile);
        TaskA->SetTrackPtRange(0., 100.0);
        TaskA->SetTrackPhiRange(0., 2.0*pi);
        TaskA->SetTrackEtaRange(-10.0, 10.0);
        TaskA->SetEventZVtxRange(-40., 40.);  // TEST
        TaskA->SetClusterPtRange(0., 500.0);
        TaskA->SetTowerERange(0., 500.0);
        TaskA->SetUsePrimaryTracks(usePrimaryTracks);
        TaskA->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        TaskA->SetRunFlag(RunFlag);                      // RunFlag
        //TaskA->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        TaskA->SetTurnOnCentSelection(kFALSE);  // run analysis for specific centrality
        TaskA->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        TaskA->SetHadronicCorrFrac(1.0);
        TaskA->SetDoTowerQAforHT(kFALSE);
        TaskA->SetdoppAnalysis(dopp);
        TaskA->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        //TaskA->SetDebugLevel(99); // TEST
        TaskA->SetMaxEventTrackPt(100.0); // TEST

        // ======================================================================================================= 
 /*       StPicoTrackClusterQA *Task1 = new StPicoTrackClusterQA("TrackClusterQAHT1", kTRUE, outputFile);
        Task1->SetTrackPtRange(0.2, 30.0);
        Task1->SetTrackPhiRange(0., 2.0*pi);
        Task1->SetTrackEtaRange(-1.0, 1.0);
        Task1->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task1->SetClusterPtRange(0.2, 100.0);
        Task1->SetTowerERange(0.2, 100.0);
        Task1->SetUsePrimaryTracks(usePrimaryTracks);
        Task1->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT1);    // kIsHT1 or kIsHT2 or kIsHT3
        Task1->SetRunFlag(RunFlag);                      // RunFlag
        Task1->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task1->SetTurnOnCentSelection(kFALSE);  // run analysis for specific centrality
        Task1->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task1->SetDebugLevel(8);
        //Task1->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task1->SetHadronicCorrFrac(1.0);
        Task1->SetDoTowerQAforHT(kTRUE);
        Task1->SetdoppAnalysis(dopp);
        Task1->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Task1->SetMaxEventTrackPt(100.0); // TEST

        // ======================================================================================================= 
        StPicoTrackClusterQA *Task1a = new StPicoTrackClusterQA("TrackClusterQAHT1_010", kTRUE, outputFile);
        Task1a->SetTrackPtRange(0.2, 30.0);
        Task1a->SetTrackPhiRange(0., 2.0*pi);
        Task1a->SetTrackEtaRange(-1.0, 1.0);
        Task1a->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task1a->SetClusterPtRange(0.2, 100.0);
        Task1a->SetTowerERange(0.2, 100.0);
        Task1a->SetUsePrimaryTracks(usePrimaryTracks);
        Task1a->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT1);    // kIsHT1 or kIsHT2 or kIsHT3
        Task1a->SetRunFlag(RunFlag);                      // RunFlag
        Task1a->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task1a->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        Task1a->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task1->SetDebugLevel(8);
        //Task1->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task1a->SetHadronicCorrFrac(1.0);
        Task1a->SetDoTowerQAforHT(kTRUE);
        Task1a->SetdoppAnalysis(dopp);
        Task1a->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Task1a->SetMaxEventTrackPt(100.0); // TEST

        // ======================================================================================================= 
        // QA task - HT events
        StPicoTrackClusterQA *Task1b = new StPicoTrackClusterQA("TrackClusterQAHT1_6080", kTRUE, outputFile);
        Task1b->SetTrackPtRange(0.2, 30.0);
        Task1b->SetTrackPhiRange(0., 2.0*pi);
        Task1b->SetTrackEtaRange(-1.0, 1.0);
        Task1b->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task1b->SetClusterPtRange(0.2, 100.0);
        Task1b->SetTowerERange(0.2, 100.0);
        Task1b->SetUsePrimaryTracks(usePrimaryTracks);
        Task1b->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT1);    // kIsHT1 or kIsHT2 or kIsHT3
        Task1b->SetRunFlag(RunFlag);                      // RunFlag
        Task1b->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task1b->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        Task1b->SetCentralityBinCut(CentralitySelection2); // specific centrality range to run
        //Task1->SetDebugLevel(8);
        //Task1->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task1b->SetHadronicCorrFrac(1.0);
        Task1b->SetDoTowerQAforHT(kTRUE);
        Task1b->SetdoppAnalysis(dopp);
        Task1b->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Task1b->SetMaxEventTrackPt(100.0); // TEST

        // ======================================================================================================= 
        // QA task - HT events
        StPicoTrackClusterQA *Task2 = new StPicoTrackClusterQA("TrackClusterQAHT2", kTRUE, outputFile);
        Task2->SetTrackPtRange(0.2, 30.0);
        Task2->SetTrackPhiRange(0., 2.0*pi);
        Task2->SetTrackEtaRange(-1.0, 1.0);
        Task2->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task2->SetClusterPtRange(0.2, 100.0);
        Task2->SetTowerERange(0.2, 100.0);
        Task2->SetUsePrimaryTracks(usePrimaryTracks);
        Task2->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT2);    // kIsHT1 or kIsHT2 or kIsHT3
        Task2->SetRunFlag(RunFlag);                      // RunFlag
        Task2->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task2->SetTurnOnCentSelection(kFALSE);  // run analysis for specific centrality
        Task2->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task2->SetDebugLevel(8);
        //Task2->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task2->SetHadronicCorrFrac(1.0);
        Task2->SetDoTowerQAforHT(kTRUE);
        Task2->SetdoppAnalysis(dopp);
        Task2->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Task2->SetMaxEventTrackPt(100.0); // TEST

        // QA task - HT events
        StPicoTrackClusterQA *Task2a = new StPicoTrackClusterQA("TrackClusterQAHT2_010", kTRUE, outputFile);
        Task2a->SetTrackPtRange(0.2, 30.0);
        Task2a->SetTrackPhiRange(0., 2.0*pi);
        Task2a->SetTrackEtaRange(-1.0, 1.0);
        Task2a->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task2a->SetClusterPtRange(0.2, 100.0);
        Task2a->SetTowerERange(0.2, 100.0);
        Task2a->SetUsePrimaryTracks(usePrimaryTracks);
        Task2a->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT2);    // kIsHT1 or kIsHT2 or kIsHT3
        Task2a->SetRunFlag(RunFlag);                      // RunFlag
        Task2a->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task2a->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        Task2a->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task2->SetDebugLevel(8);
        //Task2->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task2a->SetHadronicCorrFrac(1.0);
        Task2a->SetDoTowerQAforHT(kTRUE);
        Task2a->SetdoppAnalysis(dopp);
        Task2a->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Task2a->SetMaxEventTrackPt(100.0); // TEST

        StPicoTrackClusterQA *Task2b = new StPicoTrackClusterQA("TrackClusterQAHT2_6080", kTRUE, outputFile);
        Task2b->SetTrackPtRange(0.2, 30.0);
        Task2b->SetTrackPhiRange(0., 2.0*pi);
        Task2b->SetTrackEtaRange(-1.0, 1.0);
        Task2b->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task2b->SetClusterPtRange(0.2, 100.0);
        Task2b->SetTowerERange(0.2, 100.0);
        Task2b->SetUsePrimaryTracks(usePrimaryTracks);
        Task2b->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT2);    // kIsHT1 or kIsHT2 or kIsHT3
        Task2b->SetRunFlag(RunFlag);                      // RunFlag
        Task2b->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task2b->SetTurnOnCentSelection(kTRUE);  // run analysis for specific centrality
        Task2b->SetCentralityBinCut(CentralitySelection2); // specific centrality range to run
        //Task2->SetDebugLevel(8);
        //Task2->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task2b->SetHadronicCorrFrac(1.0);
        Task2b->SetDoTowerQAforHT(kTRUE);
        Task2b->SetdoppAnalysis(dopp);
        Task2b->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
        Task2b->SetMaxEventTrackPt(100.0); // TEST
*/
        // =======================================================================================================
/*
        // QA task - HT1
        StPicoTrackClusterQA *Task1 = new StPicoTrackClusterQA("TrackClusterQAHT1", kTRUE, outputFile);
        Task1->SetTrackPtRange(0.2, 30.0);
        Task1->SetTrackPhiRange(0., 2.0*pi);
        Task1->SetTrackEtaRange(-1.0, 1.0);
        Task1->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task1->SetClusterPtRange(0.2, 100.0);
        Task1->SetTowerERange(0.2, 100.0);
        Task1->SetUsePrimaryTracks(usePrimaryTracks);
        Task1->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT1);    // kIsHT1 or kIsHT2 or kIsHT3
        Task1->SetRunFlag(RunFlag);                      // RunFlag
        Task1->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task1->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        Task1->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        Task1->SetHadronicCorrFrac(1.0);
        Task1->SetDoTowerQAforHT(kTRUE);
        Task1->SetdoppAnalysis(dopp);
        Task1->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs

        // QA task - HT2
        StPicoTrackClusterQA *Task2 = new StPicoTrackClusterQA("TrackClusterQAHT2", kTRUE, outputFile);
        Task2->SetTrackPtRange(0.2, 30.0);
        Task2->SetTrackPhiRange(0., 2.0*pi);
        Task2->SetTrackEtaRange(-1.0, 1.0);
        Task2->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task2->SetClusterPtRange(0.2, 100.0);
        Task2->SetTowerERange(0.2, 100.0);
        Task2->SetUsePrimaryTracks(usePrimaryTracks);
        Task2->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT2);    // kIsHT1 or kIsHT2 or kIsHT3
        Task2->SetRunFlag(RunFlag);                      // RunFlag
        Task2->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task2->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        Task2->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        Task2->SetHadronicCorrFrac(1.0);
        Task2->SetDoTowerQAforHT(kTRUE);
        Task2->SetdoppAnalysis(dopp);
        Task2->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
*/

        // initialize chain
        chain->Init();
        cout<<"chain->Init();"<<endl;
        int total = picoMaker->chain()->GetEntries();
        cout << " Total entries = " << total << endl;
        if(nEvents>total) nEvents = total;
  
        for (Int_t i=0; i<nEvents; i++){
          if(i%100==0) cout << "Working on eventNumber " << i << endl;

          chain->Clear();
          int iret = chain->Make(i);	
          if (iret) { cout << "Bad return code!" << iret << endl; break;}

          total++;		
	}
	
	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;
	
	delete chain;	

        // close output file if open
        if(fout->IsOpen()) fout->Close();
}

void LoadLibs()
{
  // load fastjet libraries 3.x
  //gSystem->Load("libCGAL"); - not installed 
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjettools");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");

  // add include path to use its functionality - FIXME update with your own path
  gSystem->AddIncludePath("-I/star/u/jmazer19/Y2017/STAR/FastJet/fastjet-install/include");

  // load the system libraries - these were defaults
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  // these are needed for new / additional classes
  gSystem->Load("libStPicoEvent");
  gSystem->Load("libStPicoDstMaker");

  // my libraries
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StMyAnalysisMaker");

  gSystem->ListLibraries();
} 

void LoadMacros()
{
}
