#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "TDatabasePDG.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using CC = CCs::iterator;
using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCMu, aod::McTrackLabels>;
using TC = TCs::iterator;


struct strangeness_tutorial {
  
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaReco{"LambdaReco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaTruth{"LambdaTruth", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.97, "V0 CosPA"}; 
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};
 
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};

  TDatabasePDG* pdg = nullptr;
  Preslice<aod::McParticles> partPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<CCs> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<TCs> trackPerMcParticle = aod::mctracklabel::mcParticleId;

  void init(InitContext const&)
  {
    
    pdg = TDatabasePDG::Instance();
    
    
    AxisSpec LambdaMassAxis = {100, 0.9f, 1.3f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};

    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    rLambdaReco.add("hMassLambda", "Histogram of Minv of lambda particles from MC reconstructed", {HistType::kTH1F, {LambdaMassAxis}});
    rLambdaReco.add("hPtLambda", "Histogram of pT of lambda particles from MC reconstructed", {HistType::kTH1F, {ptAxis}});
    rLambdaTruth.add("hMassLambda", "Histogram of Minv of lambda particles from MC truth", {HistType::kTH1F, {LambdaMassAxis}});
    rLambdaTruth.add("hPtLambda", "Histogram of pT of lambda particles from MC truth", {HistType::kTH1F, {ptAxis}});
  }

  
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
 
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                          nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                          aod::v0data::dcaV0daughters < v0setting_dcav0dau);
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr,aod::McTrackLabels>;

  void processReco(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                    soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>> const& V0s, aod::McParticles const& mcParticles,
                    DaughterTracks const& // no need to define a variable for tracks, if we don't access them directly
                    )
  {
    
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    for (const auto& v0 : V0s) {
    const auto& posDaughterTrack = v0.posTrack_as<DaughterTracks>();
    const auto& negDaughterTrack = v0.negTrack_as<DaughterTracks>();
    
    if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0setting_cospa)
      continue;
    if (v0.v0radius() < v0setting_radius)
      continue;
    if (TMath::Abs(posDaughterTrack.tpcNSigmaPr()) > NSigmaTPCProton) {
    continue;
    }
    if (TMath::Abs(negDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
    continue;
    }

    if (1.1< v0.mLambda() && v0.mLambda() < 1.15) {       
    if (0.3<posDaughterTrack.tpcInnerParam() && posDaughterTrack.tpcInnerParam()<4){
      if(0.16<negDaughterTrack.tpcInnerParam() && negDaughterTrack.tpcInnerParam()<4){
        
          if(v0.pt()>0.5 && v0.pt()<2.0){
            if(v0.eta()>-0.8 && v0.eta()<0.8){
              rLambdaReco.fill(HIST("hMassLambda"), v0.mLambda());
              std::cout<<"processReco";
              rLambdaReco.fill(HIST("hPtLambda"), v0.pt());
            }}}}}}};

    PROCESS_SWITCH(strangeness_tutorial, processReco, "Process reconstructed data", true);

    
  void processTruth(aod::McCollisions const& mccollisions, CCs const& collisions,
    aod::McParticles const& McParts, TCs const& tracks){
    
    // loop over all genererated collisions
    for (auto mccollision : mccollisions) {
     
    // get McParticles which belong to mccollision
      auto partSlice = McParts.sliceBy(partPerMcCollision, mccollision.globalIndex());
      for (auto McPart : partSlice) {
   
          TParticlePDG* lambda = pdg->GetParticle(3122);
          if(McPart.has_daughters()){
            if (McPart.pdgCode() == 3122) {
              
              if(McPart.isPhysicalPrimary()){
                if(McPart.pt()>0.5 && McPart.pt()<2.0){
                  if(McPart.eta()>-0.8 && McPart.eta()<0.8){
                    rLambdaTruth.fill(HIST("hPtLambda"), McPart.pt());
                    rLambdaTruth.fill(HIST("hMassLambda"), lambda->Mass());}}}}}
    }}};
    PROCESS_SWITCH(strangeness_tutorial, processTruth, "Process MC truth data", true);
  };


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeness_tutorial>(cfgc)};
}

