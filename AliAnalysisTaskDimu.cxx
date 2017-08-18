/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskDimu
/// Analysis task for muon paits in the spectrometer.
/// The output is a list of CF containers.
/// The macro class can run on AODs or ESDs.
/// If Monte Carlo information is present, some basics checks are performed.
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

#define AliAnalysisTaskDimu_cxx

#include "AliAnalysisTaskDimu.h"

#include <algorithm>

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
//#include "TMCProcess.h"
#include "TDatabasePDG.h"
#include "TList.h"
#include "TPaveStats.h"
#include "TPRegexp.h"
#include "THashList.h"
#include "AliMultiplicity.h"

// STEER includes
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

// CORRFW includes
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#include "AliCFEffGrid.h"

// PWG includes
#include "AliMergeableCollection.h"
#include "AliCounterCollection.h"
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"
#include "AliAnalysisMuonUtility.h"
#include "AliUtilityMuonAncestor.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDimu) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskDimu::AliAnalysisTaskDimu() :
AliAnalysisTaskSE(),
fSelectedPairTypes(""),
fMergeableCollection(0x0),
fSparse(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskDimu::AliAnalysisTaskDimu ( const char *name ) :
AliAnalysisTaskSE(name),
fSelectedPairTypes(""),
fMergeableCollection(0x0),
fSparse(0x0)
{
  //
  /// Constructor.
  //
  DefineOutput(1, AliMergeableCollection::Class());
}


//________________________________________________________________________
AliAnalysisTaskDimu::~AliAnalysisTaskDimu()
{
  //
  /// Destructor
  //
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    delete fMergeableCollection;
  }
  delete fSparse;
}

//________________________________________________________________________
void AliAnalysisTaskDimu::NotifyRun()
{
  /// Set run number for cuts
  fMuonPairCuts.SetRun(fInputHandler);
}

//________________________________________________________________________
void AliAnalysisTaskDimu::SetTrackletDistCuts ( Double_t* cuts, Int_t nCuts )
{
  /// Set cuts on tracklet distance
  /// defined as in AliMultiplicity::CalcDist();
  if ( ! fTrackletDistCuts.empty() ) {
    fTrackletDistCuts.clear();
  }
  for ( Int_t icut=0; icut<nCuts; ++icut ) {
    fTrackletDistCuts.push_back(cuts[icut]);
  }
  std::sort(fTrackletDistCuts.begin(),fTrackletDistCuts.end(),std::greater<Double_t>());
}

//________________________________________________________________________
TObject* AliAnalysisTaskDimu::GetMergeableObject ( TString identifier, TString objectName )
{
  /// Get mergeable object

  TObject* obj = fMergeableCollection->GetObject(identifier.Data(), objectName.Data());
  if ( obj ) return obj;


  if ( objectName == "DimuSparse" ) obj = fSparse->Clone(objectName.Data());
  else if ( objectName == "nevents" ) {
    TH1* histo = new TH1D(objectName.Data(),objectName.Data(),1,0.5,1.5);
    obj = histo;
  }
  else {
    AliError(Form("Unknown object %s\n",objectName.Data()));
  }

  fMergeableCollection->Adopt(identifier, obj);
  AliInfo(Form("Mergeable object collection size %g MB", fMergeableCollection->EstimateSize()/1024.0/1024.0));
  return obj;
}

//___________________________________________________________________________
void AliAnalysisTaskDimu::UserCreateOutputObjects()
{

  Int_t nPtBins = 100;
  Double_t ptMin = 0., ptMax = 100.;
  TString ptTitle("p_{T}"), ptUnits("GeV/c");

  Int_t nEtaBins = 25;
  Double_t etaMin = -4.5, etaMax = -2.;
  TString etaTitle("y"), etaUnits("");

  Int_t nPhiBins = 36;
  Double_t phiMin = 0., phiMax = 2.*TMath::Pi();
  TString phiTitle("#phi"), phiUnits("rad");

  // Int_t nInvMassBins = 1500;
  // Double_t invMassMin = 0., invMassMax = 150.;
  Int_t nInvMassBins = 750;
  Double_t invMassMin = 0., invMassMax = 15.;
  TString invMassTitle("M_{#mu#mu}"), invMassUnits("GeV/c^{2}");

  Int_t nMultBins = 10;
  Double_t multBinMin = 0., multBinMax = 100.;
  TString multBinTitle(Form("Centrality (%s)",fMuonEventCuts.GetCentralityEstimator().Data())), multBinUnits("");

  Int_t nTrackletBins = 150;
  Double_t trackletMin = -0.5, trackletMax = 150.-0.5;
  TString trackletTitle("SPD tracklets"), trackletUnits("");

//  Int_t nMotherTypeBins = kNpairTypes;
//  Double_t motherTypeMin = -0.5, motherTypeMax = (Double_t)kNpairTypes - 0.5;
//  TString motherType("MotherType"), motherTypeTitle("motherType"), motherTypeUnits("");

  Int_t nbins[kNvars] = {nPtBins, nEtaBins, nPhiBins, nInvMassBins, nMultBins, nTrackletBins};
  Double_t xmin[kNvars] = {ptMin, etaMin, phiMin, invMassMin, multBinMin, trackletMin};
  Double_t xmax[kNvars] = {ptMax, etaMax, phiMax, invMassMax, multBinMax, trackletMax};
  TString axisTitle[kNvars] = {ptTitle, etaTitle, phiTitle, invMassTitle, multBinTitle, trackletTitle};
  TString axisUnits[kNvars] = {ptUnits, etaUnits, phiUnits, invMassUnits, multBinUnits, trackletUnits};

  fSparse = new THnSparseF("BaseDimuSparse","Sparse for tracks",kNvars,nbins);

  TString histoTitle = "";
  for ( Int_t idim = 0; idim<kNvars; idim++ ) {
    histoTitle = Form("%s (%s)", axisTitle[idim].Data(), axisUnits[idim].Data());
    histoTitle.ReplaceAll("()","");
    fSparse->GetAxis(idim)->SetTitle(histoTitle.Data());

    Double_t array[nbins[idim]+1];
    for ( Int_t ibin=0; ibin<=nbins[idim]; ibin++ ) array[ibin] = xmin[idim] + ibin * (xmax[idim]-xmin[idim])/nbins[idim] ;
    fSparse->SetBinEdges(idim, array);
  }

  fMergeableCollection = new AliMergeableCollection(GetOutputSlot(1)->GetContainer()->GetName());
  fMuonEventCuts.Print("mask");
  fMuonPairCuts.Print("mask");

  AliInfo(Form("The task will store the results for %s",fSelectedPairTypes.IsNull()?"all particles":fSelectedPairTypes.Data()));

  TString trackletDistCuts = "";
  for ( auto& val : fTrackletDistCuts ) {
    trackletDistCuts += Form("  %g",val);
  }
  if ( trackletDistCuts.IsNull() ) trackletDistCuts = "none";
  AliInfo(Form("Cuts on tracklet distance: %s",trackletDistCuts.Data()));

  PostData(1,fMergeableCollection);
}

//________________________________________________________________________
void AliAnalysisTaskDimu::UserExec ( Option_t * /*option*/ )
{
  //
  /// Fill output objects
  //

  if ( ! fMuonEventCuts.IsSelected(fInputHandler) ) return;

  AliMultiplicity* mult = dynamic_cast<AliMultiplicity*>(InputEvent()->GetMultiplicity());
  int nTrackletDistCuts = fTrackletDistCuts.size();
  std::vector<Int_t> nTrackletsPerCut(nTrackletDistCuts+1,0);
  std::vector<TString> trackletDistCutsName;
  for ( auto & val : fTrackletDistCuts ) {
    trackletDistCutsName.push_back(Form("trackletDistCuts_%g",val));
  }
  trackletDistCutsName.push_back(Form("trackletDistCuts_none"));

  const TObjArray* selectTrigClasses = fMuonEventCuts.GetSelectedTrigClassesInEvent(fInputHandler);

  Double_t containerInput[kNvars];
  containerInput[kHcentrality] = fMuonEventCuts.GetCentrality(InputEvent());
  AliTrackMore* trackMore = 0x0, *trackMore2 = 0x0;
  AliVParticle* track = 0x0, *track2 = 0x0;

  Int_t nSteps = MCEvent() ? 2 : 1;
  for ( Int_t istep = 0; istep<nSteps; ++istep ) {
    std::vector<TString> selTrigClasses;
    if ( istep == kStepReconstructed ) {
      TIter nextTrig(selectTrigClasses);
      TObject* obj;
      while ( (obj = nextTrig()) ) selTrigClasses.push_back(obj->GetName());
    }
    else selTrigClasses.push_back("generated");

    for ( auto& trigClass : selTrigClasses ) {
      TString identifier = Form("/%s",trigClass.Data());
      static_cast<TH1*>(GetMergeableObject(identifier, "nevents"))->Fill(1.);
    }

    Int_t nTracks = ( istep == kStepReconstructed ) ? AliAnalysisMuonUtility::GetNTracks(InputEvent()) : MCEvent()->GetNumberOfTracks();


    // First select tracks
    TObjArray selectedTracks(nTracks);
    selectedTracks.SetOwner();
    Int_t nSelected = 0;
    for (Int_t itrack = 0; itrack < nTracks; itrack++) {
      track = ( istep == kStepReconstructed ) ? AliAnalysisMuonUtility::GetTrack(itrack,InputEvent()) : MCEvent()->GetTrack(itrack);

      // In case of MC we usually ask that the particle is a muon
      // However, in W or Z simulations, Pythia stores both the initial muon
      // (before ISR, FSR and kt kick) and the final state one.
      // The first muon is of course there only for information and should be rejected.
      // The Pythia code for initial state particles is 21
      // When running with POWHEG, Pythia puts the hard process input of POWHEG in the stack
      // with state 21, and then re-add it to stack before applying ISR, FSR and kt kick.
      // This muon produces the final state muon, and its status code is 11
      // To avoid all problems, keep only final state muon (status code <10)
      // FIXME: is the convention valid for other generators as well?
      Bool_t isSelected = ( istep == kStepReconstructed ) ? fMuonPairCuts.GetMuonTrackCuts().IsSelected(track) : ( TMath::Abs(track->PdgCode()) == 13 && AliAnalysisMuonUtility::GetStatusCode(track) < 10 ) && track->Eta() < -2.5 && track->Eta() > -4.;
      if ( ! isSelected ) continue;

      // Add per trigger information
      trackMore = new AliTrackMore(track);
      trackMore->SetParticleType(fUtilityDimuonSource.GetParticleType(track,MCEvent()));
      trackMore->SetHistory(AliAnalysisMuonUtility::GetTrackHistory(track,MCEvent()));
      trackMore->SetLabel((istep==kStepReconstructed)?track->GetLabel():itrack);
      // if ( istep == kStepReconstructed ) {
      //   for ( auto& trigClass : selTrigClasses ) {
      //     if ( fMuonPairCuts.GetMuonTrackCuts().TrackPtCutMatchTrigClass(track,fMuonEventCuts.GetTrigClassPtCutLevel(trigClass)) ) trackMore->SetPassTrigClassCut(itrig);
      //   }
      // }

      selectedTracks[nSelected++] = trackMore;
    } // loop on tracks

    if ( nSelected < 2 ) continue;

    // Loop on selected tracks
    for ( Int_t itrack=0; itrack<nSelected; itrack++) {
      trackMore = static_cast<AliTrackMore*>(selectedTracks.UncheckedAt(itrack));
      track = trackMore->GetTrack();

      // Check dimuons
      for ( Int_t jtrack=itrack+1; jtrack<nSelected; jtrack++ ) {
        trackMore2 = static_cast<AliTrackMore*>(selectedTracks.UncheckedAt(jtrack));
        track2 = trackMore2->GetTrack();
        // if ( track->Charge() * track2->Charge() >= 0 ) continue;
        TString chargeType = ( track->Charge() * track2->Charge() >= 0 ) ? "SS" : "OS";

        Int_t commonAncestor = fUtilityDimuonSource.GetCommonAncestor(track,track2,MCEvent());
        TString pairType = fUtilityDimuonSource.GetPairType(trackMore->GetParticleType(), trackMore2->GetParticleType(), commonAncestor, MCEvent());

        if ( ! fSelectedPairTypes.IsNull() ) {
          TPRegexp re(Form("(^|,)%s(,|$)",pairType.Data()));
          if ( ! fSelectedPairTypes.Contains(re) ) continue;
        }

        TLorentzVector dimuPair = AliAnalysisMuonUtility::GetTrackPair(track,track2);

        Double_t phi = dimuPair.Phi();

        if ( phi < 0. ) phi += 2.*TMath::Pi(); // phi in [0,2pi]

        containerInput[kHvarPt]         = dimuPair.Pt();
        containerInput[kHvarY]          = dimuPair.Rapidity();
        containerInput[kHvarPhi]        = phi;
        containerInput[kHvarInvMass]    = dimuPair.M();

        if ( mult ) {
          std::fill(nTrackletsPerCut.begin(), nTrackletsPerCut.end(), 0);
          Int_t nTracklets = mult->GetNumberOfTracklets();
          for ( Int_t itrk=0; itrk<nTracklets; ++itrk ) {
            if ( TMath::Abs(phi - mult->GetPhi(itrk)) > TMath::Pi()/2. ) continue; // MODIFY ME!
            Double_t dist = mult->CalcDist(itrk);
            // printf ("Track %i  dist %g  pass:",itrk,dist); // REMEMBER TO CUT
            for ( Int_t icut=0; icut<nTrackletDistCuts; ++icut ) {
              // Cuts are ordered, so if it does not pass this cut
              // it will not pass the following either
              if ( dist > fTrackletDistCuts[icut] ) break;
              // printf (" %g (%i)",fTrackletDistCuts[icut],nTrackletsPerCut[icut]); // REMEMBER TO CUT
              ++nTrackletsPerCut[icut];
            }
            // printf("\n"); // REMEMBER TO CUT
            ++nTrackletsPerCut[nTrackletDistCuts];
          }
        }


        AliDebug(1,Form("Srcs: %i %i  ancestor %i Type %s\n%s\n%s\n",trackMore->GetParticleType(), trackMore2->GetParticleType(), commonAncestor, pairType.Data(), trackMore->GetHistory().Data(), trackMore2->GetHistory().Data()));

        for ( auto& trigClass : selTrigClasses ) {
          if ( istep == kStepReconstructed ) {
            if ( ! fMuonPairCuts.TrackPtCutMatchTrigClass(track,track2,fMuonEventCuts.GetTrigClassPtCutLevel(trigClass)) ) continue;
          }
          for ( Int_t icut=0; icut<nTrackletDistCuts+1; ++icut ) {
            containerInput[kHtracklets] = nTrackletsPerCut[icut];
            TString identifier = Form("/%s/%s/%s/%s",trigClass.Data(),trackletDistCutsName[icut].Data(),pairType.Data(),chargeType.Data());
            static_cast<THnSparse*>(GetMergeableObject(identifier, "DimuSparse"))->Fill(containerInput,1.);
          } // loop on tracklets cuts
        } // loop on selected trigger classes
      } // loop on second track
    } // loop on tracks
  } // loop on container steps

  PostData(1,fMergeableCollection);
}


//________________________________________________________________________
void AliAnalysisTaskDimu::Terminate(Option_t *)
{
  //
  /// Draw some histograms at the end.
  //

  fMergeableCollection = static_cast<AliMergeableCollection*>(GetOutputData(1));

  if ( ! fMergeableCollection ) return;

  Int_t srcColors[] = {kBlack, kRed, kSpring, kTeal, kBlue, kViolet, kMagenta, kOrange, kGray};
  Int_t nColors = sizeof(srcColors)/sizeof(srcColors[0]);

  TList* trigClasses = fMergeableCollection->CreateListOfKeys(0);
  TList* trackletDistCuts = fMergeableCollection->CreateListOfKeys(1);
  TList* srcs = fMergeableCollection->CreateListOfKeys(2);
  TList* chargeTypes = fMergeableCollection->CreateListOfKeys(3);
  TIter nextClass(trigClasses);
  TIter nextTrackletDistCut(trackletDistCuts);
  TIter nextSrc(srcs);
  TIter nextCharge(chargeTypes);
  TObjString *trigClass = 0x0, *trackletDistCut = 0x0, *src = 0x0, *chargeType = 0x0;

  TString genName = "generated";

  THashList histoList;

  while ( (trigClass = static_cast<TObjString*>(nextClass())) ) {
    nextTrackletDistCut.Reset();
    while ( (trackletDistCut = static_cast<TObjString*>(nextTrackletDistCut())) ) {
      nextCharge.Reset();
      while ( (chargeType = static_cast<TObjString*>(nextCharge())) ) {
        nextSrc.Reset();
        while ( (src = static_cast<TObjString*>(nextSrc()) ) ) {
          TString identifier =  Form("/%s/%s/%s/%s",trigClass->GetName(),trackletDistCut->GetName(),src->GetName(),chargeType->GetName());
          THnSparse* sparse = static_cast<THnSparse*>(fMergeableCollection->GetObject(Form("%s/DimuSparse",identifier.Data()))); ;
          if ( ! sparse ) continue;
          AliCFGridSparse gridSparse;
          gridSparse.SetGrid(static_cast<THnSparse*>(sparse->Clone()));
          AliAnalysisMuonUtility::SetSparseRange(&gridSparse, kHvarY, "", -3.999, -2.501);
          for ( Int_t iproj=0; iproj<kNvars; ++iproj ) {
            TH1* histo = gridSparse.Project(iproj);
            if ( histo->GetEntries() == 0 ) {
              delete histo;
              continue;
            }
            TString histoName = Form("%s_%s_%s_%s_proj%i",trigClass->GetName(),trackletDistCut->GetName(),chargeType->GetName(),src->GetName(),iproj);
            histo->SetName(histoName.Data());
            histo->SetDirectory(0);
            histo->Sumw2();
            histoList.Add(histo);
          } // loop on projections
        } // loop on sources
      } // loop on OS/SS
    } // loop on tracklet dist cuts
  } // loop on trigger classes

  nextClass.Reset();
  while ( (trigClass = static_cast<TObjString*>(nextClass())) ) {
    nextTrackletDistCut.Reset();
    while ( (trackletDistCut = static_cast<TObjString*>(nextTrackletDistCut())) ) {
      nextCharge.Reset();
      while ( (chargeType = static_cast<TObjString*>(nextCharge())) ) {
        for ( Int_t ieff=0; ieff<2; ieff++ ) {
          if ( ieff == 1 && trigClass->String() == genName ) continue;
          TCanvas* can = NULL;
          TLegend* leg = NULL;
          nextSrc.Reset();
          Int_t isrc = -1;
          while ( (src = static_cast<TObjString*>(nextSrc()) ) ) {
            ++isrc;
            for ( Int_t iproj=0; iproj<kNvars; ++iproj ) {
              TString histoName = Form("%s_%s_%s_%s_proj%i",trigClass->GetName(),trackletDistCut->GetName(),chargeType->GetName(),src->GetName(),iproj);
              TH1* histo = static_cast<TH1*>(histoList.FindObject(histoName.Data()));
              if ( ! histo ) continue;
              if ( ieff == 1 ) {
                TString genHistoName = histoName;
                genHistoName.ReplaceAll(trigClass->GetName(),genName.Data());
                TH1* genHisto = static_cast<TH1*>(histoList.FindObject(genHistoName.Data()));
                if ( ! genHisto ) continue;
                if ( iproj == kHvarInvMass ) {
                  TAxis* axis = histo->GetXaxis();
                  Int_t minBin = axis->FindBin(60.001);
                  Int_t maxBin = axis->FindBin(119.999);
                  Double_t num = histo->Integral(minBin,maxBin);
                  Double_t den = genHisto->Integral(minBin,maxBin);
                  printf("\nEff for %s in (%g<%s<%g): %g / %g = %g\n", histoName.Data(),axis->GetBinLowEdge(minBin),axis->GetTitle(),axis->GetBinUpEdge(maxBin),num,den,den==0.?0.:num/den);
                }
                histoName.Append("_Efficiency");
                histo = static_cast<TH1*>(histo->Clone(histoName));
                // Reset maximum or the "beutify" later on will not work properly
                histo->SetMaximum(-1111);
                histo->Divide(genHisto);
              }
              if ( ! can ) {
                TString canName = Form("%s_%s_%s_%s", GetName(), trigClass->GetName(), trackletDistCut->GetName(), chargeType->GetName());
                if ( ieff == 1 ) canName.Append("_Efficiency");
                can = new TCanvas(canName.Data(),canName.Data(),200+50*ieff,100+50*ieff,800,600);
                can->Divide(3,2);
                leg = new TLegend(0.5,0.5,0.9,0.9);
              }
              can->cd(iproj+1);
              if ( ( iproj == kHvarPt || iproj == kHvarInvMass ) && ieff == 0 ) {
                gPad->SetLogy();
              }
              Int_t icolor = ( isrc < nColors ) ? srcColors[isrc] : isrc+2;
              histo->SetLineColor(icolor);
              histo->SetMarkerColor(icolor);
              histo->SetMarkerStyle(20+isrc);

              //          histo->GetYaxis()->SetRangeUser(minY,maxY);
              TString drawOpt = ( gPad->GetListOfPrimitives() == 0 ) ? "e" : "esames";
              histo->Draw(drawOpt.Data());
              gPad->Modified();
              gPad->Update();
              TPaveStats* paveStats = static_cast<TPaveStats*>(histo->FindObject("stats"));
              if ( paveStats ) paveStats->SetTextColor(icolor);
              if ( iproj == 0 ) leg->AddEntry(histo,src->GetName(),"lp");
            } // loop on projections
          } // loop on srcs

          // Change scale
          if ( ! can ) continue;
          for ( Int_t ipad=1; ipad<=4; ipad++ ) {
            // Draw legend
            can->cd(ipad);
            if ( ipad == 1 && leg->GetNRows() > 0 ) leg->Draw();
            // Beautify canvases
            TIter nextObj(gPad->GetListOfPrimitives());
            TObject* obj = 0x0;
            Double_t maxY = 0.;
            std::vector<TH1*> hList;
            while ( (obj = nextObj()) ) {
              if ( ! obj->InheritsFrom(TH1::Class()) ) continue;
              TH1* histo = static_cast<TH1*>(obj);
              maxY = TMath::Max(histo->GetMaximum(),maxY);
              hList.push_back(histo);
            }
            for ( TH1* histo : hList ) {
              Double_t minY = histo->GetYaxis()->GetXmin();
              maxY *= 1.1;
              if ( gPad->GetLogy() ) {
                minY = 0.1;
                maxY *= 2.;
              }
              histo->GetYaxis()->SetRangeUser(minY,maxY);
            }
            gPad->Modified();
            gPad->Update();
          } // loop on pad
        } // loop on yields/efficiency
      } // loop on OS/SS
    } // loop on tracklet dist cuts
  } // loop on event type
  delete trigClasses;
  delete srcs;
  delete chargeTypes;
}


///////////////////////////////////////////////////////////////////////////////
//
// AliTrackMore
//
///////////////////////////////////////////////////////////////////////////////

class AliTrackMore;

//_____________________________________________________________________________
AliTrackMore::AliTrackMore ( AliVParticle* track ):
TObject(),
fTrack(track),
// fTrigClassCut(0),
fParticleType(-1),
fAncestor(-1),
fLabel(-1),
fHistory("")
{
  /// Ctr
}


//_____________________________________________________________________________
AliTrackMore::~AliTrackMore()
{
  /// Destructor (does nothing since fTrack is not owner)
}
