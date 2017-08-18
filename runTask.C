void runTask ( const char* runMode, const char* analysisMode,
               const char* inputName,
               const char* inputOptions = "",
               const char* softVersions = "",
               const char* analysisOptions = "",
               TString taskOptions = "" )
{

  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  gROOT->LoadMacro(gSystem->ExpandPathName("$TASKDIR/AliTaskSubmitter.cxx+"));
  AliTaskSubmitter sub;

  if ( ! sub.SetupAnalysis(runMode,analysisMode,inputName,inputOptions,softVersions,analysisOptions, "libPWGmuon.so AliAnalysisTaskDimu.cxx AddTaskDimuonAnalysis.C",". $ALICE_ROOT/include $ALICE_PHYSICS/include","TaskDimu") ) return;

//  sub.SetAliPhysicsBuildDir("$ALICE_WORK_DIR/BUILD/AliPhysics-latest-ali-master/AliPhysics");

//  if ( ! sub.SetupAnalysis(runMode,analysisMode,inputName,inputOptions,softVersions,analysisOptions, "PWGmuon.par AliUtilityDimuonSource.cxx AliAnalysisTaskDimu.cxx AddTaskDimuonAnalysis.C",". $ALICE_ROOT/include $ALICE_PHYSICS/include","TaskDimu") ) return;

//  sub.SetProofNworkers(1);

  Bool_t isMC = sub.IsMC();

  AliAnalysisAlien* plugin = (AliAnalysisAlien*)AliAnalysisManager::GetAnalysisManager()->GetGridHandler();

  if ( plugin ) plugin->SetGridWorkingDir("analysis"); // REMEMBER TO CHANGE

  gROOT->LoadMacro("./AddTaskDimuonAnalysis.C");
  AliAnalysisTaskDimu* task = AddTaskDimuonAnalysis(isMC);
  task->GetMuonPairCuts()->GetMuonTrackCuts().SetAllowDefaultParams(kTRUE);
  if ( isMC ) {
    if ( ! taskOptions.IsNull() ) task->SelectPairTypes(taskOptions.Data());
    // if ( taskOptions.Contains("Z0",TString::kIgnoreCase) ) task->SelectPairTypes("Z0");
  }

  // AliLog::SetClassDebugLevel("AliAnalysisTaskDimu",1);
//  AliLog::SetClassDebugLevel("AliMuonTrackSmearing",1);

  AliMuonEventCuts* eventCuts = BuildMuonEventCuts(sub.GetMap());
  if ( isMC ) eventCuts->SetTrigClassPatterns("ANY,MULU:Lpt2","");
  else eventCuts->SetTrigClassPatterns("kMUU7"); // REMEMBER TO CUT
  task->SetMuonEventCuts(eventCuts);

  Double_t trackletDistCuts[] = {0.1, 0.5};
  Int_t nTrackletDistCuts = sizeof(trackletDistCuts)/sizeof(trackletDistCuts[0]);
  task->SetTrackletDistCuts ( trackletDistCuts, nTrackletDistCuts );

  // if ( 0 ) {
  //   // task->GetMuonPairCuts()->GetMuonTrackCuts().SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca );
  //   task->GetMuonPairCuts()->GetMuonTrackCuts().SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca | AliMuonTrackCuts::kMuMatchLpt );
  //   eventCuts->SetTrigClassPatterns("ANY");
  //   task->SetMuonEventCuts(eventCuts);
  //   AliAnalysisTaskDimu* taskNoCut = AddTaskDimuonAnalysis(isMC,"NoTrigCut");
  //   taskNoCut->GetMuonPairCuts()->GetMuonTrackCuts().SetAllowDefaultParams(kTRUE);
  //   taskNoCut->GetMuonPairCuts()->GetMuonTrackCuts().SetFilterMask(AliMuonTrackCuts::kMuEta | AliMuonTrackCuts::kMuThetaAbs | AliMuonTrackCuts::kMuPdca );
  //   if ( isMC ) {
  //     if ( ! taskOptions.IsNull() ) taskNoCut->SelectPairTypes(taskOptions.Data());
  //     // if ( taskOptions.Contains("Z0",TString::kIgnoreCase) ) taskNoCut->SelectPairTypes("Z0");
  //   }
  //   taskNoCut->SetMuonEventCuts(eventCuts);
  // }

  sub.StartAnalysis();
}
