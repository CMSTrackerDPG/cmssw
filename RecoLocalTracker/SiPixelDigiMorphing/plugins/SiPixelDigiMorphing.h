#ifndef RecoLocalTracker_SiPixelDigiMorphing_SiPixelDigiMorphing_h
#define RecoLocalTracker_SiPixelDigiMorphing_SiPixelDigiMorphing_h

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

class SiPixelDigiMorphing : public edm::stream::EDProducer<> {
  
 public:
  
  explicit SiPixelDigiMorphing(const edm::ParameterSet& conf);
  ~SiPixelDigiMorphing() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event& e, const edm::EventSetup& c) override;

 private:
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> tPixelDigi;
  edm::EDPutTokenT<edm::DetSetVector<PixelDigi>> tPutPixelDigi;
};

#endif
