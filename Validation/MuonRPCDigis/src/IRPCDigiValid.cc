#include "Validation/MuonRPCDigis/interface/IRPCDigiValid.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCDigi/interface/IRPCDigiTime.h"

#include <cmath>
#include <fmt/format.h>

using namespace std;
using namespace edm;

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
  std::map<const RPCRoll *, std::vector<double>> detToSimHitXsMap;
  for (auto simIt = simHitHandle->begin(); simIt != simHitHandle->end(); ++simIt) {
    const RPCDetId rsid = simIt->detUnitId();
    const RPCRoll *roll = dynamic_cast<const RPCRoll *>(rpcGeom->roll(rsid));
    if (!roll)
      continue;

    // Only process iRPC rolls
    if (!roll->isIRPC())
      continue;

    if (detToSimHitXsMap.find(roll) == detToSimHitXsMap.end())
      detToSimHitXsMap[roll] = std::vector<double>();
    detToSimHitXsMap[roll].push_back(simIt->localPosition().x());

    const int region = rsid.region();
    const int disk = region * rsid.station();
    const GlobalPoint gp = roll->toGlobal(simIt->localPosition());

    auto match = hXY_Endcap_.find(disk);
    if (match != hXY_Endcap_.end())
      match->second->Fill(gp.x(), gp.y());
  }
  for (const auto &detToSimHitXs : detToSimHitXsMap) {
    hNSimHitPerRoll_->Fill(detToSimHitXs.second.size());
  }

  // Loop over IRPCDigis
  std::map<const RPCRoll *, std::vector<double>> detToDigiXsMap;
  for (auto detUnitIt = irpcDigiHandle->begin(); detUnitIt != irpcDigiHandle->end(); ++detUnitIt) {
    const RPCDetId rsid = (*detUnitIt).first;
    const RPCRoll *roll = dynamic_cast<const RPCRoll *>(rpcGeom->roll(rsid));
    if (!roll)
      continue;

    const IRPCDigiCollection::Range &range = (*detUnitIt).second;

    for (auto digiIt = range.first; digiIt != range.second; ++digiIt) {
      // Strip profile
      hStripProf_->Fill(digiIt->strip());

      // Bunch crossing
      hBxDist_->Fill(digiIt->bx());

      // iRPC timing
      IRPCDigiTime digiTime(*digiIt);
      hDigiTime_->Fill(digiTime.time());
      hDigiTimeLR_->Fill(digiTime.timeLR());
      hDigiTimeHR_->Fill(digiTime.timeHR());
      hDigiTimeDiff_->Fill(digiTime.timeLR() - digiTime.timeHR());

      // Coordinate Y from signal propagation
      hCoordinateY_->Fill(digiTime.coordinateY());

      // Keep digi position
      const double digiX = roll->centreOfStrip(digiIt->strip()).x();
      if (detToDigiXsMap.find(roll) == detToDigiXsMap.end())
        detToDigiXsMap[roll] = std::vector<double>();
      detToDigiXsMap[roll].push_back(digiX);
    }
  }
  for (const auto &detToDigiXs : detToDigiXsMap) {
    const auto digiXs = detToDigiXs.second;
    const int nDigi = digiXs.size();
    hNDigiPerRoll_->Fill(nDigi);

    // Fill residual plots, only for nDigi==1 and nSimHit==1
    const auto roll = detToDigiXs.first;
    const auto detId = roll->id();
    if (nDigi != 1)
      continue;
    if (detToSimHitXsMap.find(roll) == detToSimHitXsMap.end())
      continue;

    const auto simHitXs = detToSimHitXsMap[roll];
    const int nSimHit = simHitXs.size();
    if (nSimHit != 1)
      continue;

    const double dx = digiXs[0] - simHitXs[0];
    hRes_->Fill(dx);

    const int disk = detId.region() * detId.station();
    auto matchDisk = hResEndcapDisks_.find(disk);
    if (matchDisk != hResEndcapDisks_.end())
      matchDisk->second->Fill(dx);
  }
}

void IRPCDigiValid::bookHistograms(DQMStore::IBooker &booker, edm::Run const &run, edm::EventSetup const &eSetup) {
  booker.setCurrentFolder("RPCDigisV/IRPCDigis");

  const double maxXY = 800;
  const int nbinsXY = 160;

  // XY plots - iRPC is in RE3 and RE4 only
  for (int disk : {-4, -3, 3, 4}) {
    const std::string meName = fmt::format("XY_Endcap_{:+1d}", disk);
    const std::string meTitle = fmt::format("X-Y view of iRPC Endcap{:+1d};X (cm);Y (cm)", disk);
    hXY_Endcap_[disk] = booker.book2D(meName, meTitle, nbinsXY, -maxXY, maxXY, nbinsXY, -maxXY, maxXY);
    hXY_Endcap_[disk]->setOption("colz");
  }

  // Strip profile
  hStripProf_ = booker.book1D("Strip_Profile", "iRPC Strip Profile;Strip Number", 100, 0, 100);

  // Bunch crossing
  hBxDist_ = booker.book1D("Bunch_Crossing", "iRPC Bunch Crossing;Bunch crossing", 20, -10., 10.);

  // iRPC timing
  hDigiTime_ = booker.book1D("DigiTime", "iRPC Digi time (average LR+HR);Digi time (ns)", 100, -50., 50.);
  hDigiTimeLR_ = booker.book1D("DigiTimeLR", "iRPC Digi time Low Radius;Digi time (ns)", 100, -50., 50.);
  hDigiTimeHR_ = booker.book1D("DigiTimeHR", "iRPC Digi time High Radius;Digi time (ns)", 100, -50., 50.);
  hDigiTimeDiff_ = booker.book1D("DigiTimeDiff", "iRPC time difference (LR-HR);#Delta t (ns)", 100, -25., 25.);

  // Coordinate Y
  hCoordinateY_ = booker.book1D("CoordinateY", "iRPC Y coordinate from timing;Y (cm)", 100, -100., 100.);

  // Multiplicity
  hNSimHitPerRoll_ = booker.book1D("NSimHitPerRoll", "SimHit multiplicity per Roll;Multiplicity", 10, 0, 10);
  hNDigiPerRoll_ = booker.book1D("NDigiPerRoll", "Digi multiplicity per Roll;Multiplicity", 10, 0, 10);

  // Residual
  hRes_ = booker.book1D("Digi_SimHit_Difference", "Digi-SimHit difference;dx (cm)", 100, -8, 8);

  // iRPC only in RE3 and RE4
  for (int disk : {-4, -3, 3, 4}) {
    const std::string meName = fmt::format("Residual_Endcap_Disk_{:+1d}", disk);
    const std::string meTitle = fmt::format("Residual of iRPC Endcap Disk{:+1d};dx (cm)", disk);
    hResEndcapDisks_[disk] = booker.book1D(meName, meTitle, 100, -8, 8);
  }
}

DEFINE_FWK_MODULE(IRPCDigiValid);
