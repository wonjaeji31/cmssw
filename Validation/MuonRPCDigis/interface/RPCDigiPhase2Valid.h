#ifndef Validation_MuonRPCDigis_RPCDigiPhase2Valid_h
#define Validation_MuonRPCDigis_RPCDigiPhase2Valid_h

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include <map>

#include "DataFormats/RPCDigi/interface/RPCDigiPhase2Collection.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

class RPCDigiPhase2Valid : public DQMEDAnalyzer {
public:
  RPCDigiPhase2Valid(const edm::ParameterSet &ps);
  ~RPCDigiPhase2Valid() override = default;

protected:
  void analyze(const edm::Event &e, const edm::EventSetup &c) override;
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;

private:
  // RZ and XY views (SimHit based)
  MonitorElement *hRZ_;

  MonitorElement *hXY_Barrel_;
  std::map<int, MonitorElement *> hXY_Endcap_;
  std::map<int, MonitorElement *> hZPhi_;

  // Strip profile
  MonitorElement *hStripProf_;
  MonitorElement *hStripProf_RB12_, *hStripProf_RB34_;
  MonitorElement *hStripProf_Endcap_;

  // Bunch crossing distributions
  MonitorElement *hBxDist_;
  MonitorElement *hBxDisc_4Plus_;
  MonitorElement *hBxDisc_4Min_;

  // Sub-BX distribution (Phase2 specific)
  MonitorElement *hSubBxDist_;

  // Timing information (Phase2 specific)
  MonitorElement *hDigiTime_;
  MonitorElement *hDigiTimeBarrel_;
  MonitorElement *hDigiTimeEndcap_;

  // Multiplicity plots
  MonitorElement *hNSimHitPerRoll_, *hNDigiPerRoll_;

  // Residual plots
  MonitorElement *hRes_;
  MonitorElement *hResStrip_;
  std::map<int, MonitorElement *> hResBarrelLayers_;
  std::map<int, MonitorElement *> hResBarrelWheels_;
  std::map<int, MonitorElement *> hResEndcapDisks_;
  std::map<int, MonitorElement *> hResEndcapRings_;

  // Tokens
  edm::EDGetTokenT<edm::PSimHitContainer> simHitToken_;
  edm::EDGetTokenT<RPCDigiPhase2Collection> rpcDigiPhase2Token_;

  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeomToken_;
};

#endif
