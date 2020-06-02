#include <cuda_runtime.h>
#include "CUDADataFormats/Common/interface/Product.h"
#include "CUDADataFormats/Common/interface/HostProduct.h"
#include "CUDADataFormats/SiPixelDigi/interface/SiPixelDigisCUDA.h"
#include "DataFormats/SiPixelDigi/interface/SiPixelDigisSoA.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/GeomDetEnumerators.h"
#include "HeterogeneousCore/CUDACore/interface/ScopedContext.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/pixelCPEforGPU.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "CUDADataFormats/SiPixelCluster/interface/SiPixelClustersCUDA.h" 
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h" 



class SiPixelValidateDigiClusterGPU : public DQMEDAnalyzer {

public:
  explicit SiPixelValidateDigiClusterGPU(const edm::ParameterSet& iConfig);
  ~SiPixelValidateDigiClusterGPU() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  void analyze(const edm::Event& iEvent, edm::EventSetup const& iSetup) override;
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  edm::EDGetTokenT<cms::cuda::Product<SiPixelDigisCUDA>> gpuDigiToken_;
  edm::EDGetTokenT<cms::cuda::Product<SiPixelClustersCUDA>> gpuClusterToken_;
  
  
  cms::cuda::host::unique_ptr<uint32_t[]> pdigi_;
  cms::cuda::host::unique_ptr<uint32_t[]> rawIdArr_;
  cms::cuda::host::unique_ptr<uint16_t[]> adc_;
  cms::cuda::host::unique_ptr<int32_t[]> clus_;


 

  std::string topFolderName_;

  MonitorElement* hnDigis;
  MonitorElement* hnClusters;
  MonitorElement* hadc;
};

SiPixelValidateDigiClusterGPU::SiPixelValidateDigiClusterGPU(const edm::ParameterSet& iConfig)
    : gpuDigiToken_(consumes<cms::cuda::Product<SiPixelDigisCUDA>>(iConfig.getParameter<edm::InputTag>("src"))),
      gpuClusterToken_(consumes<cms::cuda::Product<SiPixelClustersCUDA>>(iConfig.getParameter<edm::InputTag>("clussrc")))

{
  topFolderName_ = "SiPixelHeterogeneousV/PixelDigisGPU";
}

void SiPixelValidateDigiClusterGPU::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& run , edm::EventSetup const& es) {
  ibooker.cd();
  ibooker.setCurrentFolder(topFolderName_);
  hnDigis = ibooker.book1D("nDigis", ";nDigis;#entries", 1000, 0, 1000000);
  hadc = ibooker.book1D("adc", "Digi ADC values", 500, 0, 50000);
 hnClusters = ibooker.book1D("nClusters", "Number of Clusters", 1000, 0, 1000000);
}

void SiPixelValidateDigiClusterGPU::analyze(const edm::Event& iEvent, edm::EventSetup const& es) {
  cms::cuda::Product<SiPixelDigisCUDA>  const& inputDataWrapped = iEvent.get(gpuDigiToken_);
  cms::cuda::ScopedContextAnalyze ctx{inputDataWrapped};
  auto const&  cudaDigis = ctx.get(inputDataWrapped);
  auto const&  cudaClusters = ctx.get(iEvent, gpuClusterToken_);
  

  auto nDigis_ = cudaDigis.nDigis();
  auto nClusters_ = cudaClusters.nClusters();

  hnDigis->Fill(nDigis_);
  hnClusters->Fill(nClusters_);

  pdigi_ = cudaDigis.pdigiToHostAsync(ctx.stream());
  rawIdArr_ = cudaDigis.rawIdArrToHostAsync(ctx.stream());
  adc_ = cudaDigis.adcToHostAsync(ctx.stream());
  clus_ = cudaDigis.clusToHostAsync(ctx.stream());  


 

  for(uint32_t idigi = 0; idigi < nDigis_; idigi++) {
    DetId detId(rawIdArr_[idigi]);
    if(detId.null())  continue; 
    if(detId.det() != DetId::Detector::Tracker)   continue;
    auto adc_idigi = adc_[idigi];
    hadc->Fill(adc_idigi);
  }

}

void SiPixelValidateDigiClusterGPU::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("siPixelClustersCUDAPreSplitting"));
 desc.add<edm::InputTag>("clussrc", edm::InputTag("siPixelClustersCUDAPreSplitting")); 

  descriptions.add("SiPixelValidateDigiClusterGPU", desc);
}
DEFINE_FWK_MODULE(SiPixelValidateDigiClusterGPU);
