#include "Validation/MuonRPCDigis/interface/RPCDigiPhase2Valid.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCDigi/interface/RPCDigiPhase2Time.h"

#include <algorithm>
#include <cmath>
#include <vector>
#include <fmt/format.h>

using namespace std;
using namespace edm;

namespace {
  // Phase2 digi time = 25 ns * bx + (25/16) ns * sbx (see RPCDigiPhase2Time)
  constexpr double kSubBxWidth = 25. / 16.;
  constexpr int kTimeHalfRange = 16;  // in units of sub-BX: +-25 ns
  constexpr int kMaxMultiplicity = 20;
  constexpr int kMaxBx = 5;  // digitizer bxRange
}  // namespace

RPCDigiPhase2Valid::RPCDigiPhase2Valid(const ParameterSet &ps) {
  simHitToken_ = consumes<PSimHitContainer>(
      ps.getUntrackedParameter<edm::InputTag>("simHitTag", edm::InputTag("g4SimHits:MuonRPCHits")));
  rpcDigiPhase2Token_ = consumes<RPCDigiPhase2Collection>(
      ps.getUntrackedParameter<edm::InputTag>("rpcDigiPhase2Tag", edm::InputTag("simMuonRPCDigisPhase2")));

  rpcGeomToken_ = esConsumes();
}

void RPCDigiPhase2Valid::analyze(const Event &event, const EventSetup &eventSetup) {
  // Get the RPC Geometry
  auto rpcGeom = eventSetup.getHandle(rpcGeomToken_);

  edm::Handle<PSimHitContainer> simHitHandle;
  edm::Handle<RPCDigiPhase2Collection> rpcDigiPhase2Handle;
  event.getByToken(simHitToken_, simHitHandle);
  event.getByToken(rpcDigiPhase2Token_, rpcDigiPhase2Handle);

  // Loop over SimHits
  std::map<const RPCRoll *, std::vector<LocalPoint>> detToSimHitsMap;
  for (auto simIt = simHitHandle->begin(); simIt != simHitHandle->end(); ++simIt) {
    const RPCDetId rsid = simIt->detUnitId();
    const RPCRoll *roll = dynamic_cast<const RPCRoll *>(rpcGeom->roll(rsid));
    if (!roll)
      continue;

    // Phase2 digis are not produced for iRPC rolls: keep them out of the
    // SimHit-digi comparison, but keep them in the SimHit maps below.
    if (!roll->isIRPC())
      detToSimHitsMap[roll].push_back(simIt->localPosition());

    const int region = rsid.region();
    const GlobalPoint gp = roll->toGlobal(simIt->localPosition());

    hRZ_->Fill(gp.z(), gp.perp());

    if (region == 0) {
      // Barrel
      hXY_Barrel_->Fill(gp.x(), gp.y());

      const int station = rsid.station();
      const int layer = rsid.layer();
      const int stla = (station <= 2) ? (2 * (station - 1) + layer) : (station + 2);

      auto match = hZPhi_.find(stla);
      if (match != hZPhi_.end()) {
        const double phiInDeg = 180. * gp.barePhi() / TMath::Pi();
        match->second->Fill(gp.z(), phiInDeg);
      }
    } else {
      // Endcap
      const int disk = region * rsid.station();
      auto match = hXY_Endcap_.find(disk);
      if (match != hXY_Endcap_.end())
        match->second->Fill(gp.x(), gp.y());
    }
  }
  for (const auto &detToSimHits : detToSimHitsMap) {
    hNSimHitPerRoll_->Fill(std::min<int>(detToSimHits.second.size(), kMaxMultiplicity));
  }

  // Loop over RPCDigiPhase2
  std::map<const RPCRoll *, std::vector<int>> detToDigiStripsMap;
  for (auto detUnitIt = rpcDigiPhase2Handle->begin(); detUnitIt != rpcDigiPhase2Handle->end(); ++detUnitIt) {
    const RPCDetId rsid = (*detUnitIt).first;
    const RPCRoll *roll = dynamic_cast<const RPCRoll *>(rpcGeom->roll(rsid));
    if (!roll)
      continue;
    const int region = rsid.region();

    const RPCDigiPhase2Collection::Range &range = (*detUnitIt).second;

    for (auto digiIt = range.first; digiIt != range.second; ++digiIt) {
      // Strip profile
      const int strip = digiIt->strip();
      hStripProf_->Fill(strip);

      if (region == 0) {
        const int station = rsid.station();
        if (station == 1 or station == 2)
          hStripProf_RB12_->Fill(strip);
        else if (station == 3 or station == 4)
          hStripProf_RB34_->Fill(strip);
      } else {
        hStripProf_Endcap_->Fill(strip);
      }

      // Bunch crossing
      const int bx = digiIt->bx();
      hBxDist_->Fill(bx);
      if (rsid.station() == 4) {
        if (region == 1)
          hBxDisc_4Plus_->Fill(bx);
        else if (region == -1)
          hBxDisc_4Min_->Fill(bx);
      }

      // Sub-BX distribution (Phase2 specific)
      hSubBxDist_->Fill(digiIt->sbx());

      // Phase2 timing
      RPCDigiPhase2Time digiTime(*digiIt);
      const double time = digiTime.time();
      hDigiTime_->Fill(time);
      if (region == 0)
        hDigiTimeBarrel_->Fill(time);
      else
        hDigiTimeEndcap_->Fill(time);

      detToDigiStripsMap[roll].push_back(strip);
    }
  }
  for (const auto &detToDigiStrips : detToDigiStripsMap) {
    const auto &digiStrips = detToDigiStrips.second;
    const int nDigi = digiStrips.size();
    hNDigiPerRoll_->Fill(std::min(nDigi, kMaxMultiplicity));

    // Fill residual plots, only for nDigi==1 and nSimHit==1
    const auto roll = detToDigiStrips.first;
    const auto detId = roll->id();
    if (nDigi != 1)
      continue;
    auto simMatch = detToSimHitsMap.find(roll);
    if (simMatch == detToSimHitsMap.end())
      continue;

    const auto &simHits = simMatch->second;
    if (simHits.size() != 1)
      continue;

    // Compare the strip centre and the SimHit at the SimHit local y.
    // Endcap strips are trapezoidal, so the strip x depends on y.
    const LocalPoint &simLp = simHits[0];
    const double dStrip = (digiStrips[0] - 0.5) - roll->strip(simLp);
    const double dx = dStrip * roll->localPitch(simLp);
    hResStrip_->Fill(dStrip);
    hRes_->Fill(dx);
    if (roll->isBarrel()) {
      const int wheel = detId.ring();
      const int station = detId.station();
      const int layer = detId.layer();
      const int stla = (station <= 2) ? (2 * (station - 1) + layer) : (station + 2);

      auto matchLayer = hResBarrelLayers_.find(stla);
      if (matchLayer != hResBarrelLayers_.end())
        matchLayer->second->Fill(dx);

      auto matchWheel = hResBarrelWheels_.find(wheel);
      if (matchWheel != hResBarrelWheels_.end())
        matchWheel->second->Fill(dx);
    } else {
      const int disk = detId.region() * detId.station();
      auto matchDisk = hResEndcapDisks_.find(disk);
      if (matchDisk != hResEndcapDisks_.end())
        matchDisk->second->Fill(dx);

      auto matchRing = hResEndcapRings_.find(detId.ring());
      if (matchRing != hResEndcapRings_.end())
        matchRing->second->Fill(dx);
    }
  }
}

void RPCDigiPhase2Valid::bookHistograms(DQMStore::IBooker &booker, edm::Run const &run, edm::EventSetup const &eSetup) {
  booker.setCurrentFolder("RPCDigisV/RPCDigiPhase2");

  // Define binnings of 2D-histograms
  const double maxZ = 1100;
  const int nbinsZ = 220;  // bin width: 10cm
  const double maxXY = 800;
  const int nbinsXY = 160;  // bin width: 10cm
  const double minR = 100, maxR = 800;
  const int nbinsR = 70;    // bin width: 10cm
  const int nbinsPhi = 72;  // bin width: 5 degree
  const double maxBarrelZ = 700;
  const int nbinsBarrelZ = 140;  // bin width: 10cm

  // RZ plot
  hRZ_ = booker.book2D("RZ", "R-Z view;Z (cm);R (cm)", nbinsZ, -maxZ, maxZ, nbinsR, minR, maxR);
  hRZ_->setOption("colz");

  // XY plots
  hXY_Barrel_ = booker.book2D("XY_Barrel", "X-Y view of Barrel", nbinsXY, -maxXY, maxXY, nbinsXY, -maxXY, maxXY);
  hXY_Barrel_->setOption("colz");
  for (int disk = 1; disk <= 4; ++disk) {
    const std::string meNameP = fmt::format("XY_Endcap_p{:1d}", disk);
    const std::string meNameN = fmt::format("XY_Endcap_m{:1d}", disk);
    const std::string meTitleP = fmt::format("X-Y view of Endcap{:+1d};X (cm);Y (cm)", disk);
    const std::string meTitleN = fmt::format("X-Y view of Endcap{:+1d};X (cm);Y (cm)", -disk);
    hXY_Endcap_[disk] = booker.book2D(meNameP, meTitleP, nbinsXY, -maxXY, maxXY, nbinsXY, -maxXY, maxXY);
    hXY_Endcap_[-disk] = booker.book2D(meNameN, meTitleN, nbinsXY, -maxXY, maxXY, nbinsXY, -maxXY, maxXY);
    hXY_Endcap_[disk]->setOption("colz");
    hXY_Endcap_[-disk]->setOption("colz");
  }

  // Z-phi plots
  for (int layer = 1; layer <= 6; ++layer) {
    const std::string meName = fmt::format("ZPhi_Layer{:1d}", layer);
    const std::string meTitle = fmt::format("Z-#phi view of Layer{:1d};Z (cm);#phi (degree)", layer);
    hZPhi_[layer] = booker.book2D(meName, meTitle, nbinsBarrelZ, -maxBarrelZ, maxBarrelZ, nbinsPhi, -180, 180);
    hZPhi_[layer]->setOption("colz");
  }

  // Strip profile: one bin per strip, centred on the strip number
  hStripProf_ = booker.book1D("Strip_Profile", "Strip Profile;Strip Number", 100, 0.5, 100.5);
  hStripProf_RB12_ = booker.book1D("Strip_Profile_RB12", "Strip Profile RB1 and RB2;Strip Number", 92, 0.5, 92.5);
  hStripProf_RB34_ = booker.book1D("Strip_Profile_RB34", "Strip Profile RB3 and RB4;Strip Number", 62, 0.5, 62.5);
  hStripProf_Endcap_ = booker.book1D("Strip_Profile_Endcap", "Strip Profile Endcap;Strip Number", 40, 0.5, 40.5);

  // Bunch crossing: one bin per BX
  const int nbinsBx = 2 * kMaxBx + 1;
  hBxDist_ = booker.book1D("Bunch_Crossing", "Bunch Crossing;Bunch crossing", nbinsBx, -kMaxBx - 0.5, kMaxBx + 0.5);
  hBxDisc_4Plus_ =
      booker.book1D("BxDisc_4Plus", "BX of Endcap Disk+4;Bunch crossing", nbinsBx, -kMaxBx - 0.5, kMaxBx + 0.5);
  hBxDisc_4Min_ =
      booker.book1D("BxDisc_4Min", "BX of Endcap Disk-4;Bunch crossing", nbinsBx, -kMaxBx - 0.5, kMaxBx + 0.5);

  // Sub-BX (Phase2 specific): 16 sub-BX per BX
  hSubBxDist_ = booker.book1D("SubBx", "Sub-BX distribution;Sub-BX", 16, -0.5, 15.5);

  // Phase2 timing: one bin per sub-BX, centred on the digitized time values
  const int nbinsTime = 2 * kTimeHalfRange + 1;
  const double maxTime = (kTimeHalfRange + 0.5) * kSubBxWidth;
  hDigiTime_ = booker.book1D("DigiTime", "RPCDigiPhase2 time;Digi time (ns)", nbinsTime, -maxTime, maxTime);
  hDigiTimeBarrel_ =
      booker.book1D("DigiTimeBarrel", "RPCDigiPhase2 Barrel time;Digi time (ns)", nbinsTime, -maxTime, maxTime);
  hDigiTimeEndcap_ =
      booker.book1D("DigiTimeEndcap", "RPCDigiPhase2 Endcap time;Digi time (ns)", nbinsTime, -maxTime, maxTime);

  // SimHit and Digi multiplicity per roll (last bin includes overflow)
  const int nbinsMult = kMaxMultiplicity + 1;
  hNSimHitPerRoll_ = booker.book1D(
      "NSimHitPerRoll", "SimHit multiplicity per Roll;Multiplicity", nbinsMult, -0.5, kMaxMultiplicity + 0.5);
  hNDigiPerRoll_ =
      booker.book1D("NDigiPerRoll", "Digi multiplicity per Roll;Multiplicity", nbinsMult, -0.5, kMaxMultiplicity + 0.5);

  // Residual of SimHit-Digi x-position, evaluated at the SimHit local y
  hRes_ = booker.book1D("Digi_SimHit_Difference", "Digi-SimHit difference;dx (cm)", 100, -8, 8);
  hResStrip_ =
      booker.book1D("Digi_SimHit_StripDifference", "Digi-SimHit difference in strip units;#Delta strip", 100, -5, 5);

  for (int layer = 1; layer <= 6; ++layer) {
    const std::string meName = fmt::format("Residual_Barrel_Layer{:1d}", layer);
    const std::string meTitle = fmt::format("Residual of Barrel Layer{:1d};dx (cm)", layer);
    hResBarrelLayers_[layer] = booker.book1D(meName, meTitle, 100, -8, 8);
  }

  hResBarrelWheels_[-2] = booker.book1D("Residual_Barrel_Wheel_m2", "Residual of Barrel Wheel-2;dx (cm)", 100, -8, 8);
  hResBarrelWheels_[-1] = booker.book1D("Residual_Barrel_Wheel_m1", "Residual of Barrel Wheel-1;dx (cm)", 100, -8, 8);
  hResBarrelWheels_[+0] = booker.book1D("Residual_Barrel_Wheel_00", "Residual of Barrel Wheel 0;dx (cm)", 100, -8, 8);
  hResBarrelWheels_[+1] = booker.book1D("Residual_Barrel_Wheel_p1", "Residual of Barrel Wheel+1;dx (cm)", 100, -8, 8);
  hResBarrelWheels_[+2] = booker.book1D("Residual_Barrel_Wheel_p2", "Residual of Barrel Wheel+2;dx (cm)", 100, -8, 8);

  for (int disk = 1; disk <= 4; ++disk) {
    const std::string meNameP = fmt::format("Residual_Endcap_Disk_p{:1d}", disk);
    const std::string meNameN = fmt::format("Residual_Endcap_Disk_m{:1d}", disk);
    const std::string meTitleP = fmt::format("Residual of Endcap Disk{:+1d};dx (cm)", disk);
    const std::string meTitleN = fmt::format("Residual of Endcap Disk{:+1d};dx (cm)", -disk);
    hResEndcapDisks_[+disk] = booker.book1D(meNameP, meTitleP, 100, -8, 8);
    hResEndcapDisks_[-disk] = booker.book1D(meNameN, meTitleN, 100, -8, 8);
  }

  // Ring 1 of RE3/RE4 is iRPC (no Phase2 digis), RE1/1 and RE2/1 have no RPC
  hResEndcapRings_[2] = booker.book1D("Residual_Endcap_Ring2", "Residual of Endcap Ring2;dx (cm)", 100, -8, 8);
  hResEndcapRings_[3] = booker.book1D("Residual_Endcap_Ring3", "Residual of Endcap Ring3;dx (cm)", 100, -8, 8);
}

DEFINE_FWK_MODULE(RPCDigiPhase2Valid);
