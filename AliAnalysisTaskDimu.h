#ifndef ALIANALYSISTASKDIMU_H
#define ALIANALYSISTASKDIMU_H

/* $Id$ */

//
// AliAnalysisTaskDimu
// Analysis task for muon pairs in the spectrometer
//
//  Author: Diego Stocco
//

#include <vector>
#include "TString.h"
#include "AliAnalysisTaskSE.h"
#include "AliMuonEventCuts.h"
#include "AliMuonPairCuts.h"
#include "AliUtilityDimuonSource.h"

class TObjArray;
class THnSparse;
class AliMergeableCollection;

class AliAnalysisTaskDimu : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDimu();
  AliAnalysisTaskDimu ( const char *name );
  virtual ~AliAnalysisTaskDimu();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
  virtual void Terminate(Option_t *option);

  /// Get muon event cuts
  AliMuonEventCuts* GetMuonEventCuts() { return &fMuonEventCuts; }
  /// Get muon pair cuts
  AliMuonPairCuts* GetMuonPairCuts() { return &fMuonPairCuts; }

  /// Set muon event cuts
  void SetMuonEventCuts ( AliMuonEventCuts* muonEventCuts ) { fMuonEventCuts = *muonEventCuts; }
  /// Set muon pair cuts
  void SetMuonPairCuts ( AliMuonPairCuts* muonPairCuts ) { fMuonPairCuts = *muonPairCuts; }

  /// Comma separated list of pair types to be kept
  /// e.g. J/psi,Upsilon,Z0,Uncorrelated
  /// See AliUtilityMuonSource and TDatabasePDF for naming conventions
  void SelectPairTypes ( TString selectedPairTypes ) { fSelectedPairTypes = selectedPairTypes; }

  void SetTrackletDistCuts ( Double_t* cuts, Int_t nCuts );

  enum {
    kStepReconstructed,  ///< Reconstructed tracks
    kStepGeneratedMC,    ///< Generated tracks (MC)
    kNsteps              ///< Number of steps
  };

  enum {
    kHvarPt,         ///< Pt
    kHvarY,          ///< Rapidity
    kHvarPhi,        ///< Phi
    kHvarInvMass,    ///< Invariant mass
    kHcentrality,    ///< event centrality
    kHtracklets,     ///< Number of SPD tracklets
    kNvars           ///< THnSparse dimensions
  };

 private:
  TObject* GetMergeableObject ( TString identifier, TString objectName );

  AliAnalysisTaskDimu(const AliAnalysisTaskDimu&);
  AliAnalysisTaskDimu& operator=(const AliAnalysisTaskDimu&);

  AliMuonEventCuts fMuonEventCuts;  ///< Muon event cuts
  AliMuonPairCuts fMuonPairCuts;  ///< Muon track cuts
  AliUtilityDimuonSource fUtilityDimuonSource; //!<! Utility to get the dimuon sources
  TString fSelectedPairTypes; ///< Selected pair types
  AliMergeableCollection* fMergeableCollection; //!<! collection of mergeable objects
  THnSparse* fSparse; ///< CF container
  std::vector<Double_t> fTrackletDistCuts; // Number of tracklet distance cuts

  ClassDef(AliAnalysisTaskDimu, 1); // Muon pair analysis
};

class AliVParticle;

class AliTrackMore : public TObject
{
public:
  AliTrackMore(AliVParticle* track);
  virtual ~AliTrackMore();

  // void SetPassTrigClassCut ( Int_t itrig ) { fTrigClassCut |= (1<<itrig); }
  // Bool_t MatchTrigPtForTrigClass ( Int_t itrig ) const { return (fTrigClassCut >> itrig) & 0x1; }

  void SetParticleType ( Int_t itype ) { fParticleType = itype; }
  Int_t GetParticleType () const { return fParticleType; }

  void SetAncestor ( Int_t ancestor ) { fAncestor = ancestor; }
  Int_t GetAncestor () const { return fAncestor; }

  void SetLabel ( Int_t label ) { fLabel = label; }
  Int_t GetLabel () const { return fLabel; }

  void SetHistory ( TString history ) { fHistory = history; }
  TString GetHistory () const { return fHistory; }

  AliVParticle* GetTrack() { return fTrack; }

private:
  AliVParticle* fTrack; // AliVParticle NOT OWNER
  // UInt_t fTrigClassCut; // Trigger class cut
  Int_t fParticleType; // Particle type
  Int_t fAncestor; // Ancestor index
  Int_t fLabel; // Position in MCEvent of MC particle
  TString fHistory; // Track history

  ClassDef(AliTrackMore,0); // AliTrackMore
};

#endif
