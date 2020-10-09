import FWCore.ParameterSet.Config as cms

siPixelDigisMorphed = cms.EDProducer(
    "SiPixelDigiMorphing",
    src = cms.InputTag('siPixelDigis'),
    )


