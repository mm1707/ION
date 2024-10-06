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

struct qaPlots {
  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  Configurable<float> ConfZVertexCut{"ConfZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};
  Configurable<float> ConfEta{"ConfEta", 0.8, "Eta cut for the global track"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < ConfZVertexCut);

  using FilteredFDCollisions = soa::Filtered<o2::aod::FDCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  HistogramRegistry rLambdaReco{"Lambda Reco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rProtonReco{"Proton Reco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaTruth{"Lambda Truth", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rProtonTruth{"Proton Truth", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryPDG{"PDGHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry originRegistry{"OriginRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  Configurable<int32_t> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<int> ConfChargePart1{"ConfChargePart1", 1, "sign of particle 1"};
  Configurable<float> ConfHPtPart1{"ConfHPtPart1", 4.05f, "higher limit for pt of particle 1"};
  Configurable<float> ConfLPtPart1{"ConfLPtPart1", 0.5f, "lower limit for pt of particle 1"};
  Configurable<float> ConfDCAV0{"ConfDCAV0", 0.5f, "limit for DCA of V0"};
  Configurable<float> Confmom{"Confmom", 0.75, "momentum threshold for particle identification using TOF"};
  Configurable<float> ConfNsigmaTPCParticle{"ConfNsigmaTPCParticle", 3.0, "TPC Sigma for particle momentum < Confmom"};
  Configurable<float> ConfNsigmaCombinedParticle{"ConfNsigmaCombinedParticle", 3.0, "TPC and TOF Sigma (combined) for particle momentum > Confmom"};
  Configurable<bool> ConfNoPDGPartOne{"ConfNoPDGPartOne", false, "0: selecting part by PDG, 1: no PID selection"};

  Configurable<bool> ConfNoPDGPartTwo{"ConfNoPDGPartTwo", false, "0: selecting part by PDG, 1: no PID selection"};
  Configurable<int> ConfV0Type1{"ConfV0Type1", 0, "select one of the V0s (lambda = 0, anti-lambda = 1, k0 = 2) for v0-v0 and Track-v0 combination"};
  Configurable<float> ConfHPtPart2{"ConfHPtPart2", 5.f, "higher limit for pt of particle 2"};
  Configurable<float> ConfLPtPart2{"ConfLPtPart2", 0.3f, "lower limit for pt of particle 2"};
  Configurable<float> ConfV0InvMassLowLimit{"ConfV0InvV0MassLowLimit", 1.111, "Lower limit of the V0 invariant mass"};
  Configurable<float> ConfV0InvMassUpLimit{"ConfV0InvV0MassUpLimit", 1.119, "Upper limit of the V0 invariant mass"};
  Configurable<int32_t> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 3122, "Particle 2 - PDG code"};

  Configurable<bool> ConfisLambda{"ConfisLambda", true, "Is V0 lambda"};
  Configurable<bool> ConfisAntilambda{"ConfisAntilambda", false, "Is V0 Antilambda"};
  static constexpr UInt_t V0ChildTable[][2] = {{0, 1}, {1, 0}, {1, 1}}; // Table to select the V0 children

  Partition<FemtoFullParticles> partsOneReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == ConfChargePart1) && (nabs(aod::femtouniverseparticle::eta) < ConfEta) && (aod::femtouniverseparticle::pt < ConfHPtPart1) && (aod::femtouniverseparticle::pt > ConfLPtPart1);
  Partition<FemtoFullParticles> partsTwoReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && (aod::femtouniverseparticle::pt < ConfHPtPart2) && (aod::femtouniverseparticle::pt > ConfLPtPart2);

  Partition<FemtoFullParticles> partsOneGen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (ConfNoPDGPartOne || aod::femtouniverseparticle::pidcut == uint32_t(ConfPDGCodePartOne)) &&
                                              aod::femtouniverseparticle::pt < ConfHPtPart1 && aod::femtouniverseparticle::pt > ConfLPtPart1&& nabs(aod::femtouniverseparticle::eta) < ConfEta;
  Partition<FemtoFullParticles> partsTwoGen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (ConfNoPDGPartTwo || aod::femtouniverseparticle::pidcut == uint32_t(ConfPDGCodePartTwo)) &&
                                              aod::femtouniverseparticle::pt < ConfHPtPart2 && aod::femtouniverseparticle::pt > ConfLPtPart2&& nabs(aod::femtouniverseparticle::eta) < ConfEta;

  void init(InitContext const&)
  {
    AxisSpec LambdaMassAxis = {100, 1.11f, 1.12f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec MultAxis = {100, 0.0f, 4000.0f, "multiplicity"};

    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hMultNtr", "hMultNtr", {HistType::kTH1F, {MultAxis}});

    rLambdaReco.add("posDaughter/hNSigmaProtonpTandTPC", "n#sigma_{TPC} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025}}});
    rLambdaReco.add("posDaughter/hNSigmaProtonpTandTOF", "n#sigma_{TOF} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025}}});
    rLambdaReco.add("posDaughter/hNSigmaPionpTandTPC", "n#sigma_{TPC} Pion vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "#sigma_{TPC}"}}});
    rLambdaReco.add("posDaughter/hNSigmaPionpTandTOF", "n#sigma_{TOF} Pion vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "#sigma_{TOF}"}}});
    rLambdaReco.add("posDaughter/hNSigmaProtonTPC", "n#sigma_{TPC} Proton", {HistType::kTH1F, {{200, -4.975, 5.025, "#sigma_{TPC}"}}});
    rLambdaReco.add("posDaughter/hNSigmaProtonTOF", "n#sigma_{TOF} Proton", {HistType::kTH1F, {{200, -4.975, 5.025, "#sigma_{TOF}"}}});
    rLambdaReco.add("posDaughter/hNSigmaPionTPC", "n#sigma_{TPC} Pion", {HistType::kTH1F, {{200, -4.975, 5.025, "#sigma_{TPC}"}}});
    rLambdaReco.add("posDaughter/hNSigmaPionTOF", "n#sigma_{TOF} Pion", {HistType::kTH1F, {{200, -4.975, 5.025, "#sigma_{TOF}"}}});
    rLambdaReco.add("posDaughter/TPC_dEdx", "TPC dE/dx;#it{p}_{T} (GeV/#it{c});#it{dE/dx};", {HistType::kTH2F, {{100, 0, 10}, {100, -50, 200}}});
    // rLambdaReco.add("posDaughter/TOF_Beta", "TOF Signal;#it{p}_{T} (GeV/#it{c});#TOF Beta;", {HistType::kTH2F, {{100, 0, 10}, {100, 0, 5}}});
    rLambdaReco.add("posDaughter/hPtProton", "#it{p}_{T} of Proton", {HistType::kTH1F, {ptAxis}});
    rLambdaReco.add("posDaughter/hPtPion", "#it{p}_{T} of Pion", {HistType::kTH1F, {ptAxis}});
    rLambdaReco.add("posDaughter/hSign", "Sign of positive daughter particles", {HistType::kTH1F, {{3, -1.f, 1.0001}}});

    rLambdaReco.add("negDaughter/hNSigmaProtonpTandTPC", "n#sigma_{TPC} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "n#sigma_{TPC}"}}});
    rLambdaReco.add("negDaughter/hNSigmaProtonpTandTOF", "n#sigma_{TOF} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "n#sigma_{TOF}"}}});
    rLambdaReco.add("negDaughter/hNSigmaPionpTandTPC", "n#sigma_{TPC} Pion vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "n#sigma_{TPC}"}}});
    rLambdaReco.add("negDaughter/hNSigmaPionpTandTOF", "n#sigma_{TOF} Pion vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "n#sigma_{TOF}"}}});
    rLambdaReco.add("negDaughter/hNSigmaProtonTPC", "n#sigma_{TPC} Proton", {HistType::kTH1F, {{200, -4.975, 5.025, "n#sigma_{TPC}"}}});
    rLambdaReco.add("negDaughter/hNSigmaProtonTOF", "n#sigma_{TOF} Proton", {HistType::kTH1F, {{200, -4.975, 5.025, "n#sigma_{TOF}"}}});
    rLambdaReco.add("negDaughter/hNSigmaPionTPC", "n#sigma_{TPC} Pion", {HistType::kTH1F, {{200, -4.975, 5.025, "n#sigma_{TPC}"}}});
    rLambdaReco.add("negDaughter/hNSigmaPionTOF", "n#sigma_{TOF} Pion", {HistType::kTH1F, {{200, -4.975, 5.025, "n#sigma_{TOF}"}}});
    rLambdaReco.add("negDaughter/TPC_dEdx", "TPC dE/dx;#it{p}_{T} (GeV/#it{c});#it{dE/dx};", {HistType::kTH2F, {{100, 0, 10}, {100, -50, 200}}});
    // rLambdaReco.add("negDaughter/TOF_Beta", "TOF Signal;#it{p}_{T} (GeV/#it{c});#TOF Beta;", {HistType::kTH2F, {{100, 0, 10}, {100, 0, 5}}});
    rLambdaReco.add("negDaughter/hPtPion", "#it{p}_{T} of Pion", {HistType::kTH1F, {ptAxis}});
    rLambdaReco.add("negDaughter/hPtProton", "#it{p}_{T} of Proton", {HistType::kTH1F, {ptAxis}});
    rLambdaReco.add("negDaughter/hSign", "Sign of negative daughter particles", {HistType::kTH1F, {{3, -1.f, 1.0001}}});

    rLambdaReco.add("hMassLambda", "#it{M}_{inv} of V0", {HistType::kTH1F, {LambdaMassAxis}});
    rLambdaReco.add("hPtLambda", "#it{p}_{T} of V0", {HistType::kTH1F, {ptAxis}});
    rLambdaReco.add("hMassPtLambda", "#it{M}_{inv} vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {LambdaMassAxis}}});
    rLambdaReco.add("hEtaLambda", "#eta of V0", kTH1F, {{100, -1., 1., "#eta"}});
    rLambdaReco.add("hPhiLambda", "#varphi of V0", kTH1F, {{360, 0, 6.28, "#varphi"}});
    rLambdaReco.add("hDCALambda", "DCA of V0", kTH1F, {{100, 0.f, 0.1, "DCA (cm)"}});
    rLambdaReco.add("hDCAdaughterLambdaReco", "DCA of daughter particles", kTH1F, {{100, 0.f, 2.0f}});
    rLambdaReco.add("hTransRadiusLambda", "#it{r}_{xy} of the decay vertex", kTH1F, {{100, 0.f, 100.f, "#it{r}_{xy} (cm)"}});

    rProtonReco.add("hNSigmaProtonTPC", "n#sigma_{TPC} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "n#sigma_{TPC}"}}});
    rProtonReco.add("hNSigmaProtonTOF", "n#sigma_{TOF} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "n#sigma_{TOF}"}}});
    rProtonReco.add("hpT", "#it{p}_{T} of Proton", kTH1F, {{ptAxis}});
    rProtonReco.add("hEta", "#eta of Proton", kTH1F, {{100, -1., 1., "#eta"}});
    rProtonReco.add("hTOF", "n#sigma_{TOF} of Proton ", kTH1F, {{200, -4.975, 5.025, "n#sigma_{TOF}"}});
    rProtonReco.add("hTPC", "n#sigma_{TPC} of Proton ", kTH1F, {{200, -4.975, 5.025, "n#sigma_{TOF}"}});
    rProtonReco.add("hDCAZ", "DCA to z vertex of Proton", kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    rProtonReco.add("hDCAXY", "DCA to xy vertex of Proton", kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    rProtonReco.add("hDCAXYandpT", "DCA to xy vertex of Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {100, -0.1, 0.1, "DCA_{xy} (cm)"}}});
    rProtonReco.add("hPhi", "#varphi of Proton", kTH1F, {{360, 0, 6.28, "#varphi"}});

    rProtonTruth.add("hNSigmaProtonTPC", "n#sigma_{TPC} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "n#sigma_{TPC}"}}});
    rProtonTruth.add("hNSigmaProtonTOF", "n#sigma_{TOF} Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {200, -4.975, 5.025, "n#sigma_{TOF}"}}});
    rProtonTruth.add("hpT", "#it{p}_{T} of Proton", kTH1F, {{ptAxis}});
    rProtonTruth.add("hEta", "#eta of Proton", kTH1F, {{100, -1., 1., "#eta"}});
    rProtonTruth.add("hTOF", "n#sigma_{TOF} of Proton ", kTH1F, {{200, -4.975, 5.025, "n#sigma_{TOF}"}});
    rProtonTruth.add("hTPC", "n#sigma_{TPC} of Proton ", kTH1F, {{200, -4.975, 5.025, "n#sigma_{TOF}"}});
    rProtonTruth.add("hDCAZ", "DCA to z vertex of Proton", kTH1F, {{100, -0.1, 0.1, "DCA_{z} (cm)"}});
    rProtonTruth.add("hDCAXY", "DCA to xy vertex of Proton", kTH1F, {{100, -0.1, 0.1, "DCA_{xy} (cm)"}});
    rProtonTruth.add("hDCAXYandpT", "DCA to xy vertex of Proton vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {100, -0.1, 0.1, "DCA_{xy} (cm)"}}});
    rProtonTruth.add("hPhi", "#varphi of Proton", kTH1F, {{360, 0, 6.28, "#varphi"}});

    rLambdaTruth.add("hMassLambda", "#it{M}_{inv} of V0", {HistType::kTH1F, {LambdaMassAxis}});
    rLambdaTruth.add("hPtLambda", "#it{p}_{T} of V0", {HistType::kTH1F, {ptAxis}});
    rLambdaTruth.add("hMassPtLambda", "#it{M}_{inv} vs #it{p}_{T}", {HistType::kTH2F, {{ptAxis}, {LambdaMassAxis}}});
    rLambdaTruth.add("hEtaLambda", "#eta of V0", kTH1F, {{100, -1., 1., "#eta"}});
    rLambdaTruth.add("hPhiLambda", "#varphi of V0", kTH1F, {{360, 0, 6.28, "#varphi"}});
    rLambdaTruth.add("hDCALambda", "DCA of V0", kTH1F, {{100, 0.f, 0.1, "DCA (cm)"}});
    rLambdaTruth.add("hDCAdaughterLambdaTruth", "DCA of daughter particles", kTH1F, {{100, 0.f, 2.0f}});
    rLambdaTruth.add("hTransRadiusLambda", "#it{r}_{xy} of the decay vertex", kTH1F, {{100, 0.f, 100.f, "#it{r}_{xy} (cm)"}});

    originRegistry.add("hOrigin/Particle1", "Origin", {HistType::kTH1F, {{10, 0, 10}}});
    originRegistry.add("hOrigin/Particle2", "Origin", {HistType::kTH1F, {{10, 0, 10}}});

    registryPDG.add("PDG/Particle1", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {8001, -4000.5, 4000.5}}});
    registryPDG.add("PDG/Particle2", "PDG;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {8001, -4000.5, 4000.5}}});
  }

  bool IsNSigmaTPC(float nsigmaTPCParticle)
  {
    if (TMath::Abs(nsigmaTPCParticle) < ConfNsigmaTPCParticle) {
      return true;
    } else {
      return false;
    }
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

  bool invMLambda(float invMassLambda, float invMassAntiLambda)
  {
    if ((invMassLambda < ConfV0InvMassLowLimit || invMassLambda > ConfV0InvMassUpLimit) && (invMassAntiLambda < ConfV0InvMassLowLimit || invMassAntiLambda > ConfV0InvMassUpLimit)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool IsParticleTPC(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};

    return IsNSigmaTPC(tpcNSigmas[id]);
  }

  void analysisReco(FilteredFDCollision& col, FemtoFullParticles const& parts, aod::FDMCParticles const&)
  {
    rEventSelection.fill(HIST("hVertexZRec"), col.posZ());
    rEventSelection.fill(HIST("hMultNtr"), col.multNtr());
    auto groupPartsOne = partsOneReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwoReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    for (auto& part : groupPartsTwo) {
      if (!invMLambda(part.mLambda(), part.mAntiLambda()))
        continue;
      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);
      if (!IsParticleTPC(posChild, V0ChildTable[ConfV0Type1][0]) || !IsParticleTPC(negChild, V0ChildTable[ConfV0Type1][1]))
        continue;

      if (ConfisAntilambda) {
        rLambdaReco.fill(HIST("posDaughter/hNSigmaPionpTandTPC"), posChild.tpcInnerParam(), posChild.tpcNSigmaPi());
        rLambdaReco.fill(HIST("posDaughter/hNSigmaPionpTandTOF"), posChild.pt(), posChild.tofNSigmaPi());
        rLambdaReco.fill(HIST("posDaughter/hNSigmaPionTPC"), posChild.tpcNSigmaPi());
        rLambdaReco.fill(HIST("posDaughter/hNSigmaPionTOF"), posChild.tofNSigmaPi());
        rLambdaReco.fill(HIST("posDaughter/TPC_dEdx"), posChild.pt(), posChild.tpcSignal());
        // rLambdaReco.fill(HIST("posDaughter/TOF_Beta"), posChild.pt(), posChild.beta());
        rLambdaReco.fill(HIST("posDaughter/hPtPion"), posChild.pt());
        rLambdaReco.fill(HIST("posDaughter/hSign"), posChild.sign());
        rLambdaReco.fill(HIST("negDaughter/hNSigmaProtonpTandTPC"), negChild.tpcInnerParam(), negChild.tpcNSigmaPr());
        rLambdaReco.fill(HIST("negDaughter/hNSigmaProtonpTandTOF"), negChild.pt(), negChild.tofNSigmaPr());
        rLambdaReco.fill(HIST("negDaughter/hNSigmaProtonTPC"), negChild.tpcNSigmaPr());
        rLambdaReco.fill(HIST("negDaughter/hNSigmaProtonTOF"), negChild.tofNSigmaPr());
        rLambdaReco.fill(HIST("negDaughter/TPC_dEdx"), negChild.pt(), negChild.tpcSignal());
        // rLambdaReco.fill(HIST("negDaughter/TOF_Beta"), negChild.pt(), negChild.beta());
        rLambdaReco.fill(HIST("negDaughter/hPtProton"), negChild.pt());
        rLambdaReco.fill(HIST("negDaughter/hSign"), negChild.sign());
        rLambdaReco.fill(HIST("hMassLambda"), part.mLambda());
        rLambdaReco.fill(HIST("hPtLambda"), part.pt());
        rLambdaReco.fill(HIST("hMassPtLambda"), part.pt(), part.mLambda());
        rLambdaReco.fill(HIST("hEtaLambda"), part.eta());
        rLambdaReco.fill(HIST("hPhiLambda"), part.phi());
        rLambdaReco.fill(HIST("hDCALambda"), part.dcaXY());
        rLambdaReco.fill(HIST("hDCAdaughterLambdaReco"), part.daughDCA());
        rLambdaReco.fill(HIST("hTransRadiusLambda"), part.transRadius());
      }

      if (ConfisLambda) {
        if (!part.has_fdMCParticle()) {
          continue;
        }
        const auto mcParticle = part.fdMCParticle();
        registryPDG.fill(HIST("PDG/Particle2"), mcParticle.pt(), mcParticle.pdgMCTruth());
        originRegistry.fill(HIST("hOrigin/Particle2"), mcParticle.partOriginMCTruth());

        rLambdaReco.fill(HIST("posDaughter/hNSigmaProtonpTandTPC"), posChild.tpcInnerParam(), posChild.tpcNSigmaPr());
        rLambdaReco.fill(HIST("posDaughter/hNSigmaProtonpTandTOF"), posChild.pt(), posChild.tofNSigmaPr());
        rLambdaReco.fill(HIST("posDaughter/hNSigmaProtonTPC"), posChild.tpcNSigmaPr());
        rLambdaReco.fill(HIST("posDaughter/hNSigmaProtonTOF"), posChild.tofNSigmaPr());
        rLambdaReco.fill(HIST("posDaughter/TPC_dEdx"), posChild.pt(), posChild.tpcSignal());
        // rLambdaReco.fill(HIST("posDaughter/TOF_Beta"), posChild.pt(), posChild.beta());
        rLambdaReco.fill(HIST("posDaughter/hPtProton"), posChild.pt());
        rLambdaReco.fill(HIST("posDaughter/hSign"), posChild.sign());
        rLambdaReco.fill(HIST("negDaughter/hNSigmaPionpTandTPC"), negChild.tpcInnerParam(), negChild.tpcNSigmaPi());
        rLambdaReco.fill(HIST("negDaughter/hNSigmaPionpTandTOF"), negChild.pt(), negChild.tofNSigmaPi());
        rLambdaReco.fill(HIST("negDaughter/hNSigmaPionTPC"), negChild.tpcNSigmaPi());
        rLambdaReco.fill(HIST("negDaughter/hNSigmaPionTOF"), negChild.tofNSigmaPi());
        rLambdaReco.fill(HIST("negDaughter/TPC_dEdx"), negChild.pt(), negChild.tpcSignal());
        // rLambdaReco.fill(HIST("negDaughter/TOF_Beta"), negChild.pt(), negChild.beta());
        rLambdaReco.fill(HIST("negDaughter/hPtPion"), negChild.pt());
        rLambdaReco.fill(HIST("negDaughter/hSign"), negChild.sign());
        rLambdaReco.fill(HIST("hMassLambda"), part.mLambda());
        rLambdaReco.fill(HIST("hPtLambda"), part.pt());
        rLambdaReco.fill(HIST("hMassPtLambda"), part.pt(), part.mLambda());
        rLambdaReco.fill(HIST("hEtaLambda"), part.eta());
        rLambdaReco.fill(HIST("hPhiLambda"), part.phi());
        rLambdaReco.fill(HIST("hDCALambda"), part.dcaXY());
        rLambdaReco.fill(HIST("hDCAdaughterLambdaReco"), part.daughDCA());
        rLambdaReco.fill(HIST("hTransRadiusLambda"), part.transRadius());
      }
    }

    for (auto& part : groupPartsOne) {
      if (part.sign() != ConfChargePart1)
        continue;
      if (TMath::Abs(part.eta()) > ConfEta)
        continue;
      if ((part.pt() > ConfHPtPart1) && (part.pt() < ConfLPtPart1))
        continue;
      if (IsNSigmaCombined(TMath::Abs(part.tpcNSigmaPr()), TMath::Abs(part.tofNSigmaPr()), part.p())) {
        rProtonReco.fill(HIST("hNSigmaProtonTPC"), part.tpcInnerParam(), part.tpcNSigmaPr());
        rProtonReco.fill(HIST("hNSigmaProtonTOF"), part.pt(), part.tofNSigmaPr());
        rProtonReco.fill(HIST("hTOF"), part.tofNSigmaPr());
        rProtonReco.fill(HIST("hpT"), part.pt());
        rProtonReco.fill(HIST("hTPC"), part.tpcNSigmaPr());
        rProtonReco.fill(HIST("hDCAZ"), part.dcaZ());
        rProtonReco.fill(HIST("hDCAXY"), part.dcaXY());
        rProtonReco.fill(HIST("hEta"), part.eta());
        rProtonReco.fill(HIST("hPhi"), part.phi());
        rProtonReco.fill(HIST("hDCAXYandpT"), part.pt(), part.dcaXY());
        if (!part.has_fdMCParticle()) {
          continue;
        }
        const auto mcParticle = part.fdMCParticle();
        registryPDG.fill(HIST("PDG/Particle1"), mcParticle.pt(), mcParticle.pdgMCTruth());
        originRegistry.fill(HIST("hOrigin/Particle1"), mcParticle.partOriginMCTruth());
      }
    }
  }
  PROCESS_SWITCH(qaPlots, analysisReco, "Enable analysis of MC Reconstructed", true);

  void analysisTruth(FilteredFDCollision& col, FemtoFullParticles const&)
  {
    auto groupPartsOne = partsOneGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwoGen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    for (auto& part : groupPartsOne) {
      rProtonTruth.fill(HIST("hNSigmaProtonTPC"), part.tpcInnerParam(), part.tpcNSigmaPr());
      rProtonTruth.fill(HIST("hNSigmaProtonTOF"), part.pt(), part.tofNSigmaPr());
      rProtonTruth.fill(HIST("hTOF"), part.tofNSigmaPr());
      rProtonTruth.fill(HIST("hpT"), part.pt());
      rProtonTruth.fill(HIST("hTPC"), part.tpcNSigmaPr());
      rProtonTruth.fill(HIST("hDCAZ"), part.dcaZ());
      rProtonTruth.fill(HIST("hDCAXY"), part.dcaXY());
      rProtonTruth.fill(HIST("hEta"), part.eta());
      rProtonTruth.fill(HIST("hPhi"), part.phi());
      rProtonTruth.fill(HIST("hDCAXYandpT"), part.pt(), part.dcaXY());
    }

    for (auto& part : groupPartsTwo) {
      rLambdaTruth.fill(HIST("hMassLambda"), part.mLambda());
      rLambdaTruth.fill(HIST("hPtLambda"), part.pt());
      rLambdaTruth.fill(HIST("hMassPtLambda"), part.pt(), part.mLambda());
      rLambdaTruth.fill(HIST("hEtaLambda"), part.eta());
      rLambdaTruth.fill(HIST("hPhiLambda"), part.phi());
      rLambdaTruth.fill(HIST("hDCALambda"), part.dcaXY());
      rLambdaTruth.fill(HIST("hDCAdaughterLambdaTruth"), part.daughDCA());
      rLambdaTruth.fill(HIST("hTransRadiusLambda"), part.transRadius());
    }
  }
  PROCESS_SWITCH(qaPlots, analysisTruth, "Enable analysis of MC Truth", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<qaPlots>(cfgc),
  };
  return workflow;
}
