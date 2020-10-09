
import FWCore.ParameterSet.Config as cms

#
from CondTools.SiPixel.SiPixelGainCalibrationService_cfi import *
from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import siPixelClusters as _siPixelClusters
siPixelClustersPreSplitting = _siPixelClusters.clone()

from Configuration.ProcessModifiers.siPixelDigiMorphing_cff import *
siPixelDigiMorphing.toModify(
    siPixelClustersPreSplitting,
    src = cms.InputTag('siPixelDigisMorphed')
    )
