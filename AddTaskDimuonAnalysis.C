#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliMuonEventCuts.h"
#include "AliMuonPairCuts.h"
#include "AliAnalysisTaskDimu.h"
#endif

AliAnalysisTaskDimu* AddTaskDimuonAnalysis(Bool_t isMC = kFALSE, TString changeName = "")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddtaskDimu", "No analysis manager to connect to.");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD") && !type.Contains("AOD")) {
    ::Error("AddtaskDimu", "Dimu task needs the manager to have an ESD or AOD input handler.");
    return NULL;
  }

  // Create container
  TString outputfile = mgr->GetCommonFileName();
  if ( ! outputfile.IsNull() ) outputfile += ":PWG_Dimu" + changeName;
  else outputfile = "DimuAnalysis" + changeName + ".root";

  TString containerName = "DimuOut" + changeName;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName.Data(),AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // Create cuts
  TString cutsName = "StdMuonEventCuts" + changeName;
  AliMuonEventCuts* muonEventCuts = new AliMuonEventCuts(cutsName.Data(), cutsName.Data());
  if ( isMC ) muonEventCuts->SetTrigClassPatterns("ANY");

  cutsName = "StdMuonPairCuts" + changeName;
  AliMuonPairCuts* muonPairCuts = new AliMuonPairCuts(cutsName.Data(), cutsName.Data());
  muonPairCuts->SetIsMC(isMC);

  // Create task
  TString taskName = "DimuTask" + changeName;
  AliAnalysisTaskDimu *dimuAnalysisTask = new AliAnalysisTaskDimu(taskName.Data());
  dimuAnalysisTask->SetMuonEventCuts(muonEventCuts);
  dimuAnalysisTask->SetMuonPairCuts(muonPairCuts);
  mgr->AddTask(dimuAnalysisTask);

   // Connect containers
   mgr->ConnectInput  (dimuAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (dimuAnalysisTask,  1, coutput1);

   return dimuAnalysisTask;
}
