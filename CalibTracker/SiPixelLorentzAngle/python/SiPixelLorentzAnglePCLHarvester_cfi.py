import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

SiPixelLorentzAnglePCLHarvester = DQMEDHarvester(
    "SiPixelLorentzAnglePCLHarvester",
    newmodulelist = cms.string('newmodule.txt'),
    dqmDir = cms.string('AlCaReco/SiPixelLorentzAngle'),
    record = cms.string("SiPixelLorentzAngleRcd"),
    fitProbCut = cms.double(0.5)
)
