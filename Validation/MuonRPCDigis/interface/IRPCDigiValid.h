#ifndef IRPCDigiValid_h
#define IRPCDigiValid_h

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include <map>

#include "DataFormats/RPCDigi/interface/IRPCDigiCollection.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

class IRPCDigiValid : public DQMEDAnalyzer {
public:
  IRPCDigiValid(const edm::ParameterSet &ps);
  ~IRPCDigiValid() override = default;

protected:
  void analyze(const edm::Event &e, const edm::EventSetup &c) override;
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;

private:
  // XY views (iRPC is endcap only: RE3/1, RE4/1)
  std::map<int, MonitorElement *> hXY_Endcap_;

  // Strip profile
  MonitorElement *hStripProf_;

  // Bunch crossing distributions
  MonitorElement *hBxDist_;

  // Timing information (iRPC specific)
  MonitorElement *hDigiTime_;
  MonitorElement *hDigiTimeLR_;
  MonitorElement *hDigiTimeHR_;
  MonitorElement *hDigiTimeDiff_;  // LR-HR time difference

  // Coordinate Y from timing (iRPC specific)
  MonitorElement *hCoordinateY_;

  // Multiplicity plots
  MonitorElement *hNSimHitPerRoll_, *hNDigiPerRoll_;

  // Residual plots
  MonitorElement *hRes_;
  std::map<int, MonitorElement *> hResEndcapDisks_;

  // Tokens
  edm::EDGetTokenT<edm::PSimHitContainer> simHitToken_;
  edm::EDGetTokenT<IRPCDigiCollection> irpcDigiToken_;

  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeomToken_;
};

#endif
