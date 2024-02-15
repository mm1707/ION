#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


struct strangeness_tutorial {

  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda{"Lambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.97, "V0 CosPA"}; 
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};
  
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};
  void init(InitContext const&)
  {
    AxisSpec LambdaMassAxis = {200, 1.05f, 1.5f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};

    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    rLambda.add("hMassLambda", "Histogram of Minv of lambda particles", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hNSigmaPosProtonFromLambda", "hNSigmaPosProtonFromLambda", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
    rLambda.add("hNSigmaNegPionFromLambda", "hNSigmaNegPionFromLambda", {HistType::kTH2F, {{ptAxis},{100, -5.f, 5.f}}});
    rLambda.add("hPtLambda", "Histogram of pT of lambda particles", {HistType::kTH1F, {ptAxis}});
    rLambda.add("hMassPtLambda", "2D Histogram of Minv vs pT", {HistType::kTH2F, {{ptAxis}, {LambdaMassAxis}}});

    rLambda.add("hMassLambdaPt1", "hMassLambdaPt1", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt2", "hMassLambdaPt2", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt3", "hMassLambdaPt3", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt4", "hMassLambdaPt4", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt5", "hMassLambdaPt5", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt6", "hMassLambdaPt6", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt7", "hMassLambdaPt7", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt8", "hMassLambdaPt8", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt9", "hMassLambdaPt9", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt10", "hMassLambdaPt10", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt11", "hMassLambdaPt11", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt12", "hMassLambdaPt12", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt13", "hMassLambdaPt13", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt14", "hMassLambdaPt14", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt15", "hMassLambdaPt15", {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaPt16", "hMassLambdaPt16", {HistType::kTH1F, {LambdaMassAxis}});
  }

 
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                          nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                          aod::v0data::dcaV0daughters < v0setting_dcav0dau);
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               soa::Filtered<aod::V0Datas> const& V0s, DaughterTracks const&)
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

    rLambda.fill(HIST("hMassLambda"), v0.mLambda());

    
    rLambda.fill(HIST("hNSigmaPosProtonFromLambda"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPr());
    rLambda.fill(HIST("hNSigmaNegPionFromLambda"),  negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi());
       
    if (0.3<posDaughterTrack.tpcInnerParam() && posDaughterTrack.tpcInnerParam()<4){
      if(0.16<negDaughterTrack.tpcInnerParam() && negDaughterTrack.tpcInnerParam()<4){
        
        rLambda.fill(HIST("hPtLambda"), v0.pt());
        rLambda.fill(HIST("hMassPtLambda"), v0.pt(), v0.mLambda());
        float range=0.5;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt1"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt2"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt3"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt4"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt5"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt6"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt7"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt8"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt9"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt10"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt11"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt12"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt13"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt14"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt15"), v0.mLambda());
        range+=0.125;
        if(v0.pt()>range && v0.pt()<range+0.125)
          rLambda.fill(HIST("hMassLambdaPt16"), v0.mLambda());
        
   
      }}}}};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeness_tutorial>(cfgc)};
}
