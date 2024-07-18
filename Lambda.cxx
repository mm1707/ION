#include <iostream>
#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoUniverse;
using namespace o2::aod::pidutils;

struct lambda 
{
    SliceCache cache;
    using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
    Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;
    
    Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
    Configurable<float> ConfZVertexCut{"ConfZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};
    Configurable<float> ConfEta{"ConfEta", 0.8, "Eta cut for the global track"};
    
    Filter collisionFilter = (nabs(aod::collision::posZ) < ConfZVertexCut);

    using FilteredFDCollisions = soa::Filtered<o2::aod::FDCollisions>;
    using FilteredFDCollision = FilteredFDCollisions::iterator;
    
    HistogramRegistry rLambda{"Lambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
    HistogramRegistry rProton{"Proton", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
    HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
    
    Configurable<int> ConfChargePart1{"ConfChargePart1", 1, "sign of particle 1"};
    Configurable<float> ConfHPtPart1{"ConfHPtPart1", 4.05f, "higher limit for pt of particle 1"};
    Configurable<float> ConfLPtPart1{"ConfLPtPart1", 0.5f, "lower limit for pt of particle 1"};
    Configurable<float> ConfDCAV0{"ConfDCAV0", 0.5f, "limit for DCA of V0"};
    Configurable<float> Confmom{"Confmom", 0.75, "momentum threshold for particle identification using TOF"};
    Configurable<float> ConfNsigmaTPCParticle{"ConfNsigmaTPCParticle", 3.0, "TPC Sigma for particle momentum < Confmom"};
    Configurable<float> ConfNsigmaCombinedParticle{"ConfNsigmaCombinedParticle", 3.0, "TPC and TOF Sigma (combined) for particle momentum > Confmom"};
    Configurable<float> Confv0radiusMin{"Confv0radiusMin", 0.2, "v0 transverse radius of the decay vertex min"};
    Configurable<float> Confv0radiusMax{"Confv0radiusMax", 100, "v0 transverse radius of the decay vertex max"};
    Configurable<float> ConfDecayVertexMax{"ConfDecayVertexMax", 100, "v0 decay vertex max in x, y and z"};

    Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 3, "NSigmaTPCPion"};
    Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 3, "NSigmaTPCProton"};

    Configurable<float> ConfHPtPart2{"ConfHPtPart2", 4.0f, "higher limit for pt of particle 2"};
    Configurable<float> ConfLPtPart2{"ConfLPtPart2", 0.3f, "lower limit for pt of particle 2"};
    Configurable<float> ConfV0InvMassLowLimit{"ConfV0InvV0MassLowLimit", 1.111, "Lower limit of the V0 invariant mass"};
    Configurable<float> ConfV0InvMassUpLimit{"ConfV0InvV0MassUpLimit", 1.119, "Upper limit of the V0 invariant mass"};
  
    Configurable<bool> ConfisLambda{"ConfisLambda", true, "Is V0 lambda"};
    Configurable<bool> ConfisAntilambda{"ConfisAntilambda", false, "Is V0 Antilambda"};

    Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == ConfChargePart1) && (nabs(aod::femtouniverseparticle::eta) < ConfEta) && (aod::femtouniverseparticle::pt < ConfHPtPart1) && (aod::femtouniverseparticle::pt > ConfLPtPart1);
    Partition<FemtoFullParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && (aod::femtouniverseparticle::pt < ConfHPtPart2) && (aod::femtouniverseparticle::pt > ConfLPtPart2);

   
    void init(InitContext const&)
        {
        AxisSpec LambdaMassAxis = {100, 1.11f, 1.12f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
        AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
        AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
        AxisSpec MultAxis = {100, 0.0f, 4000.0f, "multiplicity"};
        
        rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
        rEventSelection.add("hMultNtr", "hMultNtr", {HistType::kTH1F, {MultAxis}});
        
        rLambda.add("posDaughter/hNSigmaProtonpTandTPC", "n#sigma_{TPC} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
        rLambda.add("posDaughter/hNSigmaProtonpTandTOF", "n#sigma_{TOF} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f}}});
        rLambda.add("posDaughter/hNSigmaPionpTandTPC", "n#sigma_{TPC} Pion vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis},{100, -5.f, 5.f, "#sigma_{TPC}"}}});
        rLambda.add("posDaughter/hNSigmaPionpTandTOF", "n#sigma_{TOF} Pion vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis},{100, -5.f, 5.f, "#sigma_{TOF}"}}});
        rLambda.add("posDaughter/hNSigmaProtonTPC", "n#sigma_{TPC} Proton", {HistType::kTH1F, {{100, -5.f, 5.f, "#sigma_{TPC}"}}});
        rLambda.add("posDaughter/hNSigmaProtonTOF", "n#sigma_{TOF} Proton", {HistType::kTH1F, {{100, -5.f, 5.f, "#sigma_{TOF}"}}});
        rLambda.add("posDaughter/hNSigmaPionTPC", "n#sigma_{TPC} Pion", {HistType::kTH1F, {{100, -5.f, 5.f, "#sigma_{TPC}"}}});
        rLambda.add("posDaughter/hNSigmaPionTOF", "n#sigma_{TOF} Pion", {HistType::kTH1F, {{100, -5.f, 5.f, "#sigma_{TOF}"}}});
        rLambda.add("posDaughter/TPC_dEdx", "TPC dE/dx;#it{p}_{T} (GeV/#it{c});#it{dE/dx};", {HistType::kTH2F, {{100, 0, 10}, {100, -50, 200}}});
        //rLambda.add("posDaughter/TOF_Beta", "TOF Signal;#it{p}_{T} (GeV/#it{c});#TOF Beta;", {HistType::kTH2F, {{100, 0, 10}, {100, 0, 5}}});
        rLambda.add("posDaughter/hPtProton", "#it{p}_{T} of Proton", {HistType::kTH1F, {ptAxis}});
        rLambda.add("posDaughter/hPtPion", "#it{p}_{T} of Pion", {HistType::kTH1F, {ptAxis}});
        rLambda.add("posDaughter/hSign", "Sign of positive daughter particles", {HistType::kTH1F, {{3, -1.f, 1.0001}}});

        rLambda.add("negDaughter/hNSigmaProtonpTandTPC", "n#sigma_{TPC} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f, "n#sigma_{TPC}"}}});
        rLambda.add("negDaughter/hNSigmaProtonpTandTOF", "n#sigma_{TOF} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f, "n#sigma_{TOF}"}}});
        rLambda.add("negDaughter/hNSigmaPionpTandTPC", "n#sigma_{TPC} Pion vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis},{100, -5.f, 5.f, "n#sigma_{TPC}"}}});
        rLambda.add("negDaughter/hNSigmaPionpTandTOF", "n#sigma_{TOF} Pion vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis},{100, -5.f, 5.f, "n#sigma_{TOF}"}}});     
        rLambda.add("negDaughter/hNSigmaProtonTPC", "n#sigma_{TPC} Proton", {HistType::kTH1F, {{100, -5.f, 5.f, "n#sigma_{TPC}"}}});
        rLambda.add("negDaughter/hNSigmaProtonTOF", "n#sigma_{TOF} Proton", {HistType::kTH1F, {{100, -5.f, 5.f, "n#sigma_{TOF}"}}});
        rLambda.add("negDaughter/hNSigmaPionTPC", "n#sigma_{TPC} Pion", {HistType::kTH1F, {{100, -5.f, 5.f, "n#sigma_{TPC}"}}});
        rLambda.add("negDaughter/hNSigmaPionTOF", "n#sigma_{TOF} Pion", {HistType::kTH1F, {{100, -5.f, 5.f, "n#sigma_{TOF}"}}});
        rLambda.add("negDaughter/TPC_dEdx", "TPC dE/dx;#it{p}_{T} (GeV/#it{c});#it{dE/dx};", {HistType::kTH2F, {{100, 0, 10}, {100, -50, 200}}});
        //rLambda.add("negDaughter/TOF_Beta", "TOF Signal;#it{p}_{T} (GeV/#it{c});#TOF Beta;", {HistType::kTH2F, {{100, 0, 10}, {100, 0, 5}}});
        rLambda.add("negDaughter/hPtPion", "#it{p}_{T} of Pion", {HistType::kTH1F, {ptAxis}});
        rLambda.add("negDaughter/hPtProton", "#it{p}_{T} of Proton", {HistType::kTH1F, {ptAxis}});
        rLambda.add("negDaughter/hSign", "Sign of negative daughter particles", {HistType::kTH1F, {{3, -1.f, 1.0001}}});
        
        rLambda.add("hMassLambda", "#it{M}_{inv} of V0", {HistType::kTH1F, {LambdaMassAxis}});
        rLambda.add("hPtLambda", "#it{p}_{T} of V0", {HistType::kTH1F, {ptAxis}});
        rLambda.add("hMassPtLambda", "#it{M}_{inv} vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {LambdaMassAxis}}});
        rLambda.add("hEtaLambda", "#eta of V0", kTH1F, {{100, -1., 1., "#eta"}});
        rLambda.add("hPhiLambda", "#varphi of V0", kTH1F, {{360, 0, 6.28, "#varphi"}});
        rLambda.add("hDCALambda", "DCA of V0", kTH1F, {{100, 0.f, 0.1, "DCA (cm)"}});
        rLambda.add("hDCAdaughterLambda", "DCA of daughter particles", kTH1F, {{100, 0.f, 2.0f}});
        rLambda.add("hTransRadiusLambda", "#it{r}_{xy} of the decay vertex", kTH1F, {{100, 0.f, 100.f, "#it{r}_{xy} (cm)"}});
    
        rProton.add("hNSigmaProtonTPC", "n#sigma_{TPC} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f, "n#sigma_{TPC}"}}});
        rProton.add("hNSigmaProtonTOF", "n#sigma_{TOF} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {100, -5.f, 5.f, "n#sigma_{TOF}"}}});
        rProton.add("hpT", "#it{p}_{T} of Proton", kTH1F, {{ptAxis}});
        rProton.add("hEta", "#eta of Proton", kTH1F, {{100, -1., 1., "#eta"}});
        rProton.add("hTOF", "n#sigma_{TOF} of Proton ", kTH1F, {{100, -5, 5, "n#sigma_{TOF}"}});
        rProton.add("hTPC", "n#sigma_{TPC} of Proton ", kTH1F, {{100, -5, 5, "n#sigma_{TOF}"}});
        rProton.add("hDCAZ", "DCA to z vertex of Proton", kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
        rProton.add("hDCAXY", "DCA to xy vertex of Proton", kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
        rProton.add("hDCAXYandpT", "DCA to xy vertex of Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {100, -0.1, 0.1, "DCA_{xy} (cm)"}}});
        rProton.add("hPhi", "#varphi of Proton", kTH1F, {{360, 0, 6.28, "#varphi"}});
        }
        
        bool IsNSigmaCombined(float nsigmaTPCParticle, float nsigmaTOFParticle, float mom)
        {
          if (mom <= Confmom) {
            if (TMath::Abs(nsigmaTPCParticle) < ConfNsigmaTPCParticle) {
              return true;
            } else {
              return false;
            }
          } else {
            if (TMath::Hypot(nsigmaTOFParticle, nsigmaTPCParticle) < ConfNsigmaCombinedParticle) {
              return true;
            } else {
              return false;
            }
          }
        }
        
    void analysis(FilteredFDCollision& col, FemtoFullParticles const& parts)
        {
        rEventSelection.fill(HIST("hVertexZRec"), col.posZ());
        rEventSelection.fill(HIST("hMultNtr"), col.multNtr());
        
        auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

        
        for (auto& part : groupPartsTwo){
              if (part.mLambda() > ConfV0InvMassLowLimit && part.mLambda() < ConfV0InvMassUpLimit){
                
                const auto& posChild = parts.iteratorAt(part.index() - 2);
                const auto& negChild = parts.iteratorAt(part.index() - 1);
                bool posChildPi=false, posChildPr=false, negChildPi=false, negChildPr=false;
                if (IsNSigmaCombined(std::abs(posChild.tpcNSigmaPi()), std::abs(posChild.tofNSigmaPi()), posChild.p()) && ConfisAntilambda) {
                  rLambda.fill(HIST("posDaughter/hNSigmaPionpTandTPC"), posChild.tpcInnerParam(), posChild.tpcNSigmaPi());
                  rLambda.fill(HIST("posDaughter/hNSigmaPionpTandTOF"), posChild.pt(), posChild.tofNSigmaPi());
                  rLambda.fill(HIST("posDaughter/hNSigmaPionTPC"), posChild.tpcNSigmaPi());
                  rLambda.fill(HIST("posDaughter/hNSigmaPionTOF"), posChild.tofNSigmaPi());
                  rLambda.fill(HIST("posDaughter/TPC_dEdx"), posChild.pt(), posChild.tpcSignal());
                  //rLambda.fill(HIST("posDaughter/TOF_Beta"), posChild.pt(), posChild.beta());
                  rLambda.fill(HIST("posDaughter/hPtPion"), posChild.pt());
                  rLambda.fill(HIST("posDaughter/hSign"), posChild.sign());
                  posChildPi=true;
                }
                if (IsNSigmaCombined(std::abs(posChild.tpcNSigmaPr()), std::abs(posChild.tofNSigmaPr()), posChild.p()) && ConfisLambda) {
                  rLambda.fill(HIST("posDaughter/hNSigmaProtonpTandTPC"), posChild.tpcInnerParam(), posChild.tpcNSigmaPr());
                  rLambda.fill(HIST("posDaughter/hNSigmaProtonpTandTOF"), posChild.pt(), posChild.tofNSigmaPr());
                  rLambda.fill(HIST("posDaughter/hNSigmaProtonTPC"), posChild.tpcNSigmaPr());
                  rLambda.fill(HIST("posDaughter/hNSigmaProtonTOF"), posChild.tofNSigmaPr());
                  rLambda.fill(HIST("posDaughter/TPC_dEdx"), posChild.pt(), posChild.tpcSignal());
                  //rLambda.fill(HIST("posDaughter/TOF_Beta"), posChild.pt(), posChild.beta());
                  rLambda.fill(HIST("posDaughter/hPtProton"), posChild.pt());
                  rLambda.fill(HIST("posDaughter/hSign"), posChild.sign());
                  posChildPr=true;
                }
                if (IsNSigmaCombined(std::abs(negChild.tpcNSigmaPi()), std::abs(negChild.tofNSigmaPi()), negChild.p())&& ConfisLambda) {
                  rLambda.fill(HIST("negDaughter/hNSigmaPionpTandTPC"), negChild.tpcInnerParam(), negChild.tpcNSigmaPi());
                  rLambda.fill(HIST("negDaughter/hNSigmaPionpTandTOF"), negChild.pt(), negChild.tofNSigmaPi());
                  rLambda.fill(HIST("negDaughter/hNSigmaPionTPC"), negChild.tpcNSigmaPi());
                  rLambda.fill(HIST("negDaughter/hNSigmaPionTOF"), negChild.tofNSigmaPi());
                  rLambda.fill(HIST("negDaughter/TPC_dEdx"), negChild.pt(), negChild.tpcSignal());
                  //rLambda.fill(HIST("negDaughter/TOF_Beta"), negChild.pt(), negChild.beta());
                  rLambda.fill(HIST("negDaughter/hPtPion"), negChild.pt());
                  rLambda.fill(HIST("negDaughter/hSign"), negChild.sign());
                  negChildPi=true;
                }
                if (IsNSigmaCombined(std::abs(negChild.tpcNSigmaPr()), std::abs(negChild.tofNSigmaPr()), negChild.p()) && ConfisAntilambda) {
                  rLambda.fill(HIST("negDaughter/hNSigmaProtonpTandTPC"), negChild.tpcInnerParam(), negChild.tpcNSigmaPr());
                  rLambda.fill(HIST("negDaughter/hNSigmaProtonpTandTOF"), negChild.pt(), negChild.tofNSigmaPr());
                  rLambda.fill(HIST("negDaughter/hNSigmaProtonTPC"), negChild.tpcNSigmaPr());
                  rLambda.fill(HIST("negDaughter/hNSigmaProtonTOF"), negChild.tofNSigmaPr());
                  rLambda.fill(HIST("negDaughter/TPC_dEdx"), negChild.pt(), negChild.tpcSignal());
                  //rLambda.fill(HIST("negDaughter/TOF_Beta"), negChild.pt(), negChild.beta());
                  rLambda.fill(HIST("negDaughter/hPtProton"), negChild.pt());
                  rLambda.fill(HIST("negDaughter/hSign"), negChild.sign());
                  negChildPr=true;
                  }
                  if((posChildPi && negChildPr && ConfisAntilambda) || (negChildPi && posChildPr && ConfisLambda)){
                    rLambda.fill(HIST("hMassLambda"), part.mLambda());
                    rLambda.fill(HIST("hPtLambda"), part.pt());
                    rLambda.fill(HIST("hMassPtLambda"), part.pt(), part.mLambda());
                    rLambda.fill(HIST("hEtaLambda"), part.eta());
                    rLambda.fill(HIST("hPhiLambda"), part.phi());
                    rLambda.fill(HIST("hDCALambda"), part.dcaXY()); 
                    rLambda.fill(HIST("hDCAdaughterLambda"), part.daughDCA()); 
                    rLambda.fill(HIST("hTransRadiusLambda"), part.transRadius()); 
                    }
                }
              }
                
        for (auto& part : groupPartsOne){
            if(TMath::Abs(part.dcaZ())<0.2 && TMath::Abs(part.dcaXY())<0.1){
              if (IsNSigmaCombined(TMath::Abs(part.tpcNSigmaPr()), TMath::Abs(part.tofNSigmaPr()), part.p())){
                rProton.fill(HIST("hNSigmaProtonTPC"), part.tpcInnerParam(), part.tpcNSigmaPr());
                rProton.fill(HIST("hNSigmaProtonTOF"), part.pt(), part.tofNSigmaPr());
                rProton.fill(HIST("hTOF"), part.tofNSigmaPr());
                rProton.fill(HIST("hpT"), part.pt());
                rProton.fill(HIST("hTPC"), part.tpcNSigmaPr());
                rProton.fill(HIST("hDCAZ"), part.dcaZ());
                rProton.fill(HIST("hDCAXY"), part.dcaXY());
                rProton.fill(HIST("hEta"), part.eta());
                rProton.fill(HIST("hPhi"), part.phi());
                rProton.fill(HIST("hDCAXYandpT"), part.pt(), part.dcaXY());
              }
            }
          }
        }PROCESS_SWITCH(lambda, analysis, "Enable analysis", true);
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<lambda>(cfgc),
  };
  return workflow;
}