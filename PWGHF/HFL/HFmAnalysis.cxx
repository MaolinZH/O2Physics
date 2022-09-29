// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HFmAnalysis.cxx
/// \brief Task derived from the DQ framework and used to extract the observables on single muons needed for the HF-muon analysis.
/// \author Maolin Zhang <maolin.zhang@cern.ch>, CCNU

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;

namespace o2::aod
{
namespace dq_analysis_flags
{
DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
}

DECLARE_SOA_TABLE(EventCuts, "AOD", "EVENTCUTS", dq_analysis_flags::IsEventSelected);
} // namespace o2::aod

using MyCollisions = o2::soa::Join<aod::Collisions, aod::EvSels>;
using MyMcCollisions = o2::soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
using MyEventsSelected = o2::soa::Join<aod::Collisions, aod::EvSels, o2::aod::EventCuts>;
using MyMcEventsSelected = o2::soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, o2::aod::EventCuts>;
using MyMuons = o2::soa::Join<aod::FullFwdTracks, aod::McFwdTrackLabels>;

constexpr static uint32_t gMuonFillMap(VarManager::ObjTypes::Muon | VarManager::ObjTypes::MuonCov);
constexpr static uint32_t gEventFillMap(VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision);
constexpr static uint32_t gTrackMCFillMap(VarManager::ObjTypes::ParticleMC);

void defineHistograms(HistogramManager*, TString);

struct HfmEventSelection {
  Configurable<bool> applySoftwareTrigger{"appllySoftwareTrigger", false, "whether to apply the software trigger"};
  Configurable<float> softwareTrigger{"softwareTrigger", VarManager::kIsMuonSingleLowPt7, "software tigger flag"};
  Configurable<bool> applyCutZVtx{"applyzVtx", false, "whether to apply the VtxZ cut"};
  Configurable<float> zVtxMin{"zVtxMin", -10., "min. z of primary vertex [cm]"};
  Configurable<float> zVtxMax{"zVtxMax", 10., "max. z of primary vertex [cm]"};

  HistogramManager* fHistMan;
  OutputObj<THashList> fOutputList{"output"};

  float* fValues;
  AnalysisCompositeCut* fEventCut;

  void init(o2::framework::InitContext&)
  {
    VarManager::SetDefaultVarNames();
    fValues = new float[VarManager::kNVars];
    fEventCut = new AnalysisCompositeCut("event selection", "Event Selection", true);

    AnalysisCut fCut;
    if (applyZVtx)
      fCut.AddCut(VarManager::kVtxZ, zVtxMin, zVtxMax, false);
    if (applySoftwareTrigger)
      fCut.AddCut(softwareTrigger, 0.5, 1.5, false);
    fEventCut->AddCut(&fCut);

    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    DefineHistograms(fHistMan, "EventBeforeCuts;EventAfterCuts;");
    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  Produces<aod::EventCuts> eventSel;

  template <uint32_t TEventFillMap, typename TEvent>
  void runEventSel(TEvent const& event, aod::BCs const& bcs)
  {
    // select events
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillEvent<TEventFillMap>(event, fValues);
    fHistMan->FillHistClass("EventBeforeCuts", fValues);

    if (fEventCut->IsSelected(fValues)) {
      fHistMan->FillHistClass("EventAfterCuts", fValues);
      eventSel(1);
    } else {
      eventSel(0);
    }
  }

  void processEventDataAO2D(MyCollisions::iterator const& event, aod::BCs const& bcs)
  {
    runEventSel<gEventFillMap>(event, bcs);
  }

  void processEventMCAO2D(MyMcCollisions::iterator const& event, aod::BCs const& bcs)
  {
    runEventSel<gEventFillMap>(event, bcs);
  }

  PROCESS_SWITCH(HfmEventSelection, processEventDataAO2D, "run event selection with Data AO2D file", true);
  PROCESS_SWITCH(HfmEventSelection, processEventMCAO2D, "run event selection with MC AO2D file", false);
};

struct MuonSelection {
  Configurable<std::string> fConfigCuts{"cfgMuonCuts", "muonQualityCuts", "muon selection"};

  float* fValues;
  AnalysisCompositeCut* fTrackCut;

  o2::framework::HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  void init(o2::framework::InitContext&)
  {
    AxisSpec axispT{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{500, -5., 5., "#eta"};
    AxisSpec axisDCA{500, 0., 5., "DCA (cm)"};
    AxisSpec axisSign{5, -2.5, 2.5, "Charge"};
    AxisSpec axisP{500, 0., 500., "p (GeV/#it{c})"};
    AxisSpec axisVtxZ{400, -20., 20., "Vertex Z (cm)"};

    AxisSpec axisTrkType{5, -0.5, 4.5, "TrackType"};
    AxisSpec axisMCmask{200, -0.5, 199.5, "mcMask"};
    // kinematics for MC
    AxisSpec axispTGen{200, 0., 200., "#it{p}_{T} Truth (GeV/#it{c})"};
    AxisSpec axisEtaGen{500, -5., 5., "#eta_{Truth}"};
    AxisSpec axisPGen{500, 0., 500., "p_{Truth} (GeV/#it{c})"};
    AxisSpec axispTDif{200, -2., 2., "#it{p}_{T} diff (GeV/#it{c})"};
    AxisSpec axisEtaDif{200, -2., 2., "#eta diff"};
    AxisSpec axisPDif{200, -2., 2., "p diff (GeV/#it{c})"};

    HistogramConfigSpec cfgTHnMu{HistType::kTHnSparseD, {axispT, axisEta, axisDCA, axisSign, axisP, axisVtxZ, axisTrkType, axisMCmask}, 8};
    HistogramConfigSpec cfgTHnPt{HistType::kTHnSparseD, {axispT, axispTGen, axispTDif, axisTrkType}, 4};
    HistogramConfigSpec cfgTHnEta{HistType::kTHnSparseD, {axisEta, axisEtaGen, axisEtaDif, axisTrkType}, 4};
    HistogramConfigSpec cfgTHnP{HistType::kTHnSparseD, {axisP, axisPGen, axisPDif, axisTrkType}, 4};

    registry.add("hMuBcuts", "", cfgTHnMu);
    registry.add("hMuAcuts", "", cfgTHnMu);
    registry.add("hPTBcuts", "", cfgTHnPt);
    registry.add("hPTAcuts", "", cfgTHnPt);
    registry.add("hEtaBcuts", "", cfgTHnEta);
    registry.add("hEtaAcuts", "", cfgTHnEta);
    registry.add("hPBcuts", "", cfgTHnP);
    registry.add("hPAcuts", "", cfgTHnP);

    VarManager::SetDefaultVarNames();
    fValues = new float[VarManager::kNVars];

    fTrackCut = new AnalysisCompositeCut(true);
    TString selectStr = fConfigCuts.value;
    fTrackCut->AddCut(dqcuts::GetAnalysisCut(selectStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, typename TEvent, typename TMuons>
  void runDataMuonSel(TEvent const& event, aod::BCs const& bcs, TMuons const& tracks)
  {
    // select muons in data
    if (event.isEventSelected() == 0)
      return;

    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables, fValues);
    VarManager::FillEvent<gEventFillMap>(event, fValues);

    // loop over muon tracks
    for (auto const& track : tracks) {
      VarManager::FillTrack<TMuonFillMap>(track, fValues);

      // compute DCAXY
      std::vector<double> vecCovs{
        track.cXX(), track.cXY(), track.cYY(),
        track.cPhiX(), track.cPhiY(), track.cPhiPhi(),
        track.cTglX(), track.cTglY(), track.cTglPhi(), track.cTglTgl(),
        track.c1PtX(), track.c1PtY(), track.c1PtPhi(), track.c1PtTgl(), track.c1Pt21Pt2()};

      SMatrix55 trackCovs(vecCovs.begin(), vecCovs.end());
      SMatrix5 trackPars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());

      o2::track::TrackParCovFwd trackParcov{track.z(), trackPars, trackCovs, track.chi2()};
      trackParcov.propagateToZlinear(event.posZ());

      const auto dcaX(trackParcov.getX() - event.posX());
      const auto dcaY(trackParcov.getY() - event.posY());
      const auto dcaXY(std::sqrt(dcaX * dcaX + dcaY * dcaY));
      // Before Muon Cuts
      registry.fill(HIST("hMuBcuts"),
                    fValues[VarManager::kPt],
                    fValues[VarManager::kEta], dcaXY,
                    fValues[VarManager::kCharge], track.p(),
                    fValues[VarManager::kVtxZ],
                    fValues[VarManager::kMuonTrackType], 0);
      // After Muon Cuts
      if (fTrackCut->IsSelected(fValues))
        registry.fill(HIST("hMuAcuts"),
                      fValues[VarManager::kPt],
                      fValues[VarManager::kEta], dcaXY,
                      fValues[VarManager::kCharge], track.p(),
                      fValues[VarManager::kVtxZ],
                      fValues[VarManager::kMuonTrackType], 0);
    } // end loop over muon tracks
  }

  template <uint32_t TEventFillMap, uint32_t TMuonFillMap, uint32_t TTrackMCFillMap, typename TEvent, typename TMuons, typename TMC>
  void runMCMuonSel(TEvent const& event, aod::BCs const& bcs, TMuons const& tracks, TMC const& mc)
  {
    // select muons in MC
    if (event.isEventSelected() == 0)
      return;

    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables, fValues);
    VarManager::FillEvent<gEventFillMap>(event, fValues);

    // loop over muon tracks
    for (auto const& track : tracks) {
      VarManager::FillTrack<TMuonFillMap>(track, fValues);

      if (!track.has_mcParticle()) {
        continue;
      }
      auto mcParticle = track.mcParticle();

      VarManager::FillTrack<TTrackMCFillMap>(mcParticle, fValues);

      // compute DCAXY
      std::vector<double> vecCovs{
        track.cXX(), track.cXY(), track.cYY(),
        track.cPhiX(), track.cPhiY(), track.cPhiPhi(),
        track.cTglX(), track.cTglY(), track.cTglPhi(), track.cTglTgl(),
        track.c1PtX(), track.c1PtY(), track.c1PtPhi(), track.c1PtTgl(), track.c1Pt21Pt2()};

      SMatrix55 trackCovs(vecCovs.begin(), vecCovs.end());
      SMatrix5 trackPars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());

      o2::track::TrackParCovFwd trackParcov{track.z(), trackPars, trackCovs, track.chi2()};
      trackParcov.propagateToZlinear(event.posZ());

      const auto dcaX(trackParcov.getX() - event.posX());
      const auto dcaY(trackParcov.getY() - event.posY());
      const auto dcaXY(std::sqrt(dcaX * dcaX + dcaY * dcaY));
      // Before Muon Cuts
      registry.fill(HIST("hMuBcuts"),
                    fValues[VarManager::kPt],
                    fValues[VarManager::kEta], dcaXY,
                    fValues[VarManager::kCharge], track.p(),
                    fValues[VarManager::kVtxZ],
                    fValues[VarManager::kMuonTrackType],
                    track.mcMask());

      registry.fill(HIST("hPTBcuts"),
                    fValues[VarManager::kPt], fValues[VarManager::kMCPt], fValues[VarManager::kMCPt] - fValues[VarManager::kPt], fValues[VarManager::kMuonTrackType]);
      registry.fill(HIST("hEtaBcuts"),
                    fValues[VarManager::kEta], fValues[VarManager::kMCEta], fValues[VarManager::kMCEta] - fValues[VarManager::kEta], fValues[VarManager::kMuonTrackType]);
      registry.fill(HIST("hPBcuts"),
                    fValues[VarManager::kP], fValues[VarManager::kMCPt] * std::cosh(fValues[VarManager::kMCEta]), fValues[VarManager::kMCPt] * std::cosh(fValues[VarManager::kMCEta]) - fValues[VarManager::kP], fValues[VarManager::kMuonTrackType]);
      // After Muon Cuts
      if (fTrackCut->IsSelected(fValues)) {
        registry.fill(HIST("hMuAcuts"),
                      fValues[VarManager::kPt],
                      fValues[VarManager::kEta], dcaXY,
                      fValues[VarManager::kCharge], track.p(),
                      fValues[VarManager::kVtxZ],
                      fValues[VarManager::kMuonTrackType],
                      track.mcMask());

        registry.fill(HIST("hPTAcuts"),
                      fValues[VarManager::kPt], fValues[VarManager::kMCPt], fValues[VarManager::kMCPt] - fValues[VarManager::kPt], fValues[VarManager::kMuonTrackType]);
        registry.fill(HIST("hEtaAcuts"),
                      fValues[VarManager::kEta], fValues[VarManager::kMCEta], fValues[VarManager::kMCEta] - fValues[VarManager::kEta], fValues[VarManager::kMuonTrackType]);
        registry.fill(HIST("hPAcuts"),
                      fValues[VarManager::kP], fValues[VarManager::kMCPt] * std::cosh(fValues[VarManager::kMCEta]), fValues[VarManager::kMCPt] * std::cosh(fValues[VarManager::kMCEta]) - fValues[VarManager::kP], fValues[VarManager::kMuonTrackType]);
      }
    } // end loop over muon tracks
  }

  void processMuonDataAO2D(MyEventsSelected::iterator const& event, aod::BCs const& bcs,
                           aod::FullFwdTracks const& tracks)
  {
    runDataMuonSel<gEventFillMap, gMuonFillMap>(event, bcs, tracks);
  }

  void processMuonMCAO2D(MyMcEventsSelected::iterator const& event, aod::BCs const& bcs,
                         MyMuons const& tracks, aod::McParticles const& mc)
  {
    runMCMuonSel<gEventFillMap, gMuonFillMap, gTrackMCFillMap>(event, bcs, tracks, mc);
  }

  PROCESS_SWITCH(MuonSelection, processMuonDataAO2D, "run muonselection with Data AO2D file", true);
  PROCESS_SWITCH(MuonSelection, processMuonMCAO2D, "run muonselection with MC AO2D file", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfmEventSelection>(cfgc),
    adaptAnalysisTask<MuonSelection>(cfgc),
  };
}

void defineHistograms(HistogramManager* histMan, TString histClasses)
{
  //
  // Define here the histograms for all the classes required in analysis.
  //  The histogram classes are provided in the histClasses string, separated by semicolon ";"
  //  The histogram classes and their components histograms are defined below depending on the name of the histogram class
  //
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));

  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    if (classStr.Contains("Event"))
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", "trigger,all");
  }
}
