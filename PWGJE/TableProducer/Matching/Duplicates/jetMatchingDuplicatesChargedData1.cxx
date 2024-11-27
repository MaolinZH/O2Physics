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

// jet matching duplicates charged data task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/TableProducer/Matching/jetMatchingDuplicates.cxx"

using Charged1JetDataMatchingDupliacates = JetMatchingDuplicates<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>,
                                                                 soa::Join<aod::Charged1Jets, aod::Charged1JetConstituents>,
                                                                 aod::ChargedJetsMatchedToCharged1Jets,
                                                                 aod::Charged1JetsMatchedToChargedJets,
                                                                 aod::JTracks,
                                                                 aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<Charged1JetDataMatchingDupliacates>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-data-ch-1"}));
  return WorkflowSpec{tasks};
}