import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

SiPixelLorentzAnglePCLHarvester = DQMEDHarvester(
    "SiPixelLorentzAnglePCLHarvester",
    newmodulelist = cms.vstring("BPix_BmI_SEC8_LYR2_LDR14_MOD1", "BPix_BmO_SEC3_LYR2_LDR5_MOD1", "BPix_BmO_SEC3_LYR2_LDR5_MOD2", "BPix_BmO_SEC3_LYR2_LDR5_MOD3", "BPix_BpO_SEC1_LYR2_LDR1_MOD1", "BPix_BpO_SEC1_LYR2_LDR1_MOD2", "BPix_BpO_SEC1_LYR2_LDR1_MOD3"),
    dqmDir = cms.string('AlCaReco/SiPixelLorentzAngle'),
    record = cms.string("SiPixelLorentzAngleRcd"),
    fitProbCut = cms.double(0.5)
)
