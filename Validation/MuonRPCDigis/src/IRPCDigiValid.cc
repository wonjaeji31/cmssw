#include "Validation/MuonRPCDigis/interface/IRPCDigiValid.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCDigi/interface/IRPCDigiTime.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <algorithm>
#include <cmath>
#include <vector>
#include <fmt/format.h>

using namespace std;
using namespace edm;

namespace {
  // iRPC TDC time = 25 ns * BX + 2.5 ns * SBX + 0.2 ns * fine (see IRPCDigiTime).
  // 2.5 ns is not a multiple of 0.2 ns, so the LR/HR times lie on a 0.1 ns grid,
  // their average on a 0.05 ns grid.
  constexpr double kTdcGrid = 0.1;
  // Signal speed used in IRPCDigiTime::coordinateY (0.66 c), in cm/ns
  const double kSignalSpeed = 0.66 * CLHEP::c_light * CLHEP::ns / CLHEP::cm;
  // Y = v * (tLR - tHR) / 2 lies on a grid of v * 0.1 ns / 2 (about 0.99 cm)
  const double kCoordinateYGrid = kSignalSpeed * kTdcGrid / 2.;
  constexpr int kMaxMultiplicity = 20;
  constexpr int kMaxBx = 5;  // digitizer bxRange
}  // namespace

IRPCDigiValid::IRPCDigiValid(const ParameterSet &ps) {
  simHitToken_ = consumes<PSimHitContainer>(
      ps.getUntrackedParameter<edm::InputTag>("simHitTag", edm::InputTag("g4SimHits:MuonRPCHits")));
  irpcDigiToken_ = consumes<IRPCDigiCollection>(
      ps.getUntrackedParameter<edm::InputTag>("irpcDigiTag", edm::InputTag("simMuonIRPCDigis")));

  rpcGeomToken_ = esConsumes();
}

void IRPCDigiValid::analyze(const Event &event, const EventSetup &eventSetup) {
  // Get the RPC Geometry
  auto rpcGeom = eventSetup.getHandle(rpcGeomToken_);

  edm::Handle<PSimHitContainer> simHitHandle;
  edm::Handle<IRPCDigiCollection> irpcDigiHandle;
  event.getByToken(simHitToken_, simHitHandle);
  event.getByToken(irpcDigiToken_, irpcDigiHandle);

  // Loop over SimHits - only iRPC rolls (RE3/1, RE4/1)
  std::map<const RPCRoll *, std::vector<LocalPoint>> detToSimHitsMap;
  for (auto simIt = simHitHandle->begin(); simIt != simHitHandle->end(); ++simIt) {
    const RPCDetId rsid = simIt->detUnitId();
    const RPCRoll *roll = dynamic_cast<const RPCRoll *>(rpcGeom->roll(rsid));
    if (!roll)
      continue;

    // Only process iRPC rolls
    if (!roll->isIRPC())
      continue;

    detToSimHitsMap[roll].push_back(simIt->localPosition());

    const int region = rsid.region();
    const int disk = region * rsid.station();
    const GlobalPoint gp = roll->toGlobal(simIt->localPosition());

    auto match = hXY_Endcap_.find(disk);
    if (match != hXY_Endcap_.end())
      match->second->Fill(gp.x(), gp.y());
  }
  for (const auto &detToSimHits : detToSimHitsMap) {
    hNSimHitPerRoll_->Fill(std::min<int>(detToSimHits.second.size(), kMaxMultiplicity));
  }

  // Loop over IRPCDigis
  std::map<const RPCRoll *, std::vector<IRPCDigi>> detToDigisMap;
  for (auto detUnitIt = irpcDigiHandle->begin(); detUnitIt != irpcDigiHandle->end(); ++detUnitIt) {
    const RPCDetId rsid = (*detUnitIt).first;
    const RPCRoll *roll = dynamic_cast<const RPCRoll *>(rpcGeom->roll(rsid));
    if (!roll)
      continue;

    const IRPCDigiCollection::Range &range = (*detUnitIt).second;

    for (auto digiIt = range.first; digiIt != range.second; ++digiIt) {
      // Strip profile
      hStripProf_->Fill(digiIt->strip());

      // Bunch crossing (low radius FEB)
      hBxDist_->Fill(digiIt->bx());

      // iRPC timing
      IRPCDigiTime digiTime(*digiIt);
      hDigiTime_->Fill(digiTime.time());
      hDigiTimeLR_->Fill(digiTime.timeLR());
      hDigiTimeHR_->Fill(digiTime.timeHR());
      hDigiTimeDiff_->Fill(digiTime.timeLR() - digiTime.timeHR());

      // Coordinate Y from signal propagation
      hCoordinateY_->Fill(digiTime.coordinateY());

      detToDigisMap[roll].push_back(*digiIt);
    }
  }
  for (const auto &detToDigis : detToDigisMap) {
    const auto &digis = detToDigis.second;
    const int nDigi = digis.size();
    hNDigiPerRoll_->Fill(std::min(nDigi, kMaxMultiplicity));

    // SimHit comparisons, only for nSimHit==1
    const auto roll = detToDigis.first;
    const auto detId = roll->id();
    auto simMatch = detToSimHitsMap.find(roll);
    if (simMatch == detToSimHitsMap.end())
      continue;

    const auto &simHits = simMatch->second;
    if (simHits.size() != 1)
      continue;
    const LocalPoint &simLp = simHits[0];

    // Y from timing: all strips of one cluster share the same TDC values,
    // so one entry per roll
    IRPCDigiTime firstDigiTime(digis.front());
    hCoordinateYRes_->Fill(firstDigiTime.coordinateY() - simLp.y());

    // Residual of x-position, only for nDigi==1.
    // Compare at the SimHit local y: iRPC strips are long and trapezoidal.
    if (nDigi != 1)
      continue;
    const double dStrip = (digis.front().strip() - 0.5) - roll->strip(simLp);
    const double dx = dStrip * roll->localPitch(simLp);
    hResStrip_->Fill(dStrip);
    hRes_->Fill(dx);

    const int disk = detId.region() * detId.station();
    auto matchDisk = hResEndcapDisks_.find(disk);
    if (matchDisk != hResEndcapDisks_.end())
      matchDisk->second->Fill(dx);
  }
}

void IRPCDigiValid::bookHistograms(DQMStore::IBooker &booker, edm::Run const &run, edm::EventSetup const &eSetup) {
  booker.setCurrentFolder("RPCDigisV/IRPCDigis");

  // iRPC chambers are at R < 350 cm
  const double maxXY = 400;
  const int nbinsXY = 160;  // bin width: 5cm

  // XY plots - iRPC is in RE3 and RE4 only
  for (int disk : {-4, -3, 3, 4}) {
    const std::string meName = fmt::format("XY_Endcap_{}{:1d}", disk > 0 ? "p" : "m", std::abs(disk));
    const std::string meTitle = fmt::format("X-Y view of iRPC Endcap{:+1d};X (cm);Y (cm)", disk);
    hXY_Endcap_[disk] = booker.book2D(meName, meTitle, nbinsXY, -maxXY, maxXY, nbinsXY, -maxXY, maxXY);
    hXY_Endcap_[disk]->setOption("colz");
  }

  // Strip profile: one bin per strip, centred on the strip number
  hStripProf_ = booker.book1D("Strip_Profile", "iRPC Strip Profile;Strip Number", 100, 0.5, 100.5);

  // Bunch crossing: one bin per BX
  const int nbinsBx = 2 * kMaxBx + 1;
  hBxDist_ =
      booker.book1D("Bunch_Crossing", "iRPC Bunch Crossing;Bunch crossing", nbinsBx, -kMaxBx - 0.5, kMaxBx + 0.5);

  // iRPC timing: one bin per possible value, centred on it (same rule as BX,
  // strip and Phase2 time), so that no digitized value falls on a bin edge.
  const int nTimeSteps = 150;  // LR/HR: +-15 ns on the 0.1 ns grid
  const double maxTime = (nTimeSteps + 0.5) * kTdcGrid;
  hDigiTimeLR_ =
      booker.book1D("DigiTimeLR", "iRPC Digi time Low Radius;Digi time (ns)", 2 * nTimeSteps + 1, -maxTime, maxTime);
  hDigiTimeHR_ =
      booker.book1D("DigiTimeHR", "iRPC Digi time High Radius;Digi time (ns)", 2 * nTimeSteps + 1, -maxTime, maxTime);
  const int nAvgSteps = 200;  // LR+HR average: +-10 ns on the 0.05 ns grid
  const double maxAvgTime = (nAvgSteps + 0.5) * kTdcGrid / 2.;
  hDigiTime_ = booker.book1D(
      "DigiTime", "iRPC Digi time (average LR+HR);Digi time (ns)", 2 * nAvgSteps + 1, -maxAvgTime, maxAvgTime);
  const int nDiffSteps = 100;  // LR-HR: +-10 ns on the 0.1 ns grid
  const double maxDiff = (nDiffSteps + 0.5) * kTdcGrid;
  hDigiTimeDiff_ = booker.book1D(
      "DigiTimeDiff", "iRPC time difference (LR-HR);#Delta t (ns)", 2 * nDiffSteps + 1, -maxDiff, maxDiff);

  // Coordinate Y: one bin per possible value (about 0.99 cm), centred on it
  const int nYSteps = 100;  // about +-99 cm
  const double maxY = (nYSteps + 0.5) * kCoordinateYGrid;
  hCoordinateY_ = booker.book1D("CoordinateY", "iRPC Y coordinate from timing;Y (cm)", 2 * nYSteps + 1, -maxY, maxY);
  hCoordinateYRes_ = booker.book1D(
      "CoordinateY_SimHit_Difference", "iRPC Y from timing - SimHit local y;#Delta y (cm)", 80, -20., 20.);

  // Multiplicity (last bin includes overflow)
  const int nbinsMult = kMaxMultiplicity + 1;
  hNSimHitPerRoll_ = booker.book1D(
      "NSimHitPerRoll", "SimHit multiplicity per Roll;Multiplicity", nbinsMult, -0.5, kMaxMultiplicity + 0.5);
  hNDigiPerRoll_ =
      booker.book1D("NDigiPerRoll", "Digi multiplicity per Roll;Multiplicity", nbinsMult, -0.5, kMaxMultiplicity + 0.5);

  // Residual of x-position, evaluated at the SimHit local y
  hRes_ = booker.book1D("Digi_SimHit_Difference", "Digi-SimHit difference;dx (cm)", 100, -8, 8);
  hResStrip_ =
      booker.book1D("Digi_SimHit_StripDifference", "Digi-SimHit difference in strip units;#Delta strip", 100, -5, 5);

  // iRPC only in RE3 and RE4
  for (int disk : {-4, -3, 3, 4}) {
    const std::string meName = fmt::format("Residual_Endcap_Disk_{}{:1d}", disk > 0 ? "p" : "m", std::abs(disk));
    const std::string meTitle = fmt::format("Residual of iRPC Endcap Disk{:+1d};dx (cm)", disk);
    hResEndcapDisks_[disk] = booker.book1D(meName, meTitle, 100, -8, 8);
  }
}

DEFINE_FWK_MODULE(IRPCDigiValid);
