#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "SiPixelDigiMorphing.h"

SiPixelDigiMorphing::SiPixelDigiMorphing(edm::ParameterSet const& conf)
  :tPutPixelDigi(produces<edm::DetSetVector<PixelDigi>>())
{
  tPixelDigi = consumes<edm::DetSetVector<PixelDigi>>(conf.getParameter<edm::InputTag>("src"));     
}

SiPixelDigiMorphing::~SiPixelDigiMorphing() = default;

void SiPixelDigiMorphing::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("src", edm::InputTag("siPixelDigis"));
  descriptions.add("SiPixelDigiMorphing", desc);
}

void SiPixelDigiMorphing::produce(edm::Event& e, const edm::EventSetup& es) {
  
  edm::Handle<edm::DetSetVector<PixelDigi>> inputDigi;
  e.getByToken(tPixelDigi, inputDigi);
  
  auto outputDigis = std::make_unique<edm::DetSetVector<PixelDigi>>();
  edm::DetSetVector<PixelDigi>::const_iterator DSViter = (*inputDigi).begin();
  for (; DSViter != (*inputDigi).end(); DSViter++) {
    DetId detIdObject(DSViter->detId());
    edm::DetSet<PixelDigi>* detDigis = nullptr;
    auto rawId = DSViter->detId();
    edm::DetSet<PixelDigi>::const_iterator  di;
    detDigis = &(outputDigis->find_or_insert(rawId));
    for(di = DSViter->data.begin(); di != DSViter->data.end(); di++) {
      int adc = di->adc();    
      int col = di->column();  
      int row = di->row();    
      // int time = di->time();
      (*detDigis).data.emplace_back(row, col, adc, 1);      
    }
  }  
  e.put(tPutPixelDigi, std::move(outputDigis));
}

DEFINE_FWK_MODULE(SiPixelDigiMorphing);
