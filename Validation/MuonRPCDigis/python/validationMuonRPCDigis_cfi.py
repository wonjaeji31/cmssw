import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
validationMuonRPCDigis = DQMEDAnalyzer('RPCDigiValid',
    # Tag for Digis event data retrieval
    rpcDigiTag = cms.untracked.InputTag("simMuonRPCDigis"),
    # Tag for simulated hits event data retrieval
    simHitTag = cms.untracked.InputTag("g4SimHits", "MuonRPCHits"),
    # Flag to turn on/off timing plots
    digiTime = cms.untracked.bool(False)
)
from Configuration.Eras.Modifier_fastSim_cff import fastSim
fastSim.toModify(validationMuonRPCDigis, simHitTag = "MuonSimHits:MuonRPCHits")
from Configuration.Eras.Modifier_phase2_common_cff import phase2_common
phase2_common.toModify(validationMuonRPCDigis, digiTime = True)

# Phase2 RPC Digi validation (requires phase2_rpc_devel era)
validationMuonRPCDigisPhase2 = DQMEDAnalyzer('RPCDigiPhase2Valid',
    rpcDigiPhase2Tag = cms.untracked.InputTag("simMuonRPCDigisPhase2"),
    simHitTag = cms.untracked.InputTag("g4SimHits", "MuonRPCHits"),
)
fastSim.toModify(validationMuonRPCDigisPhase2, simHitTag = "MuonSimHits:MuonRPCHits")

# iRPC Digi validation (requires phase2_rpc_devel era)
validationMuonIRPCDigis = DQMEDAnalyzer('IRPCDigiValid',
    irpcDigiTag = cms.untracked.InputTag("simMuonIRPCDigis"),
    simHitTag = cms.untracked.InputTag("g4SimHits", "MuonRPCHits"),
)
fastSim.toModify(validationMuonIRPCDigis, simHitTag = "MuonSimHits:MuonRPCHits")
