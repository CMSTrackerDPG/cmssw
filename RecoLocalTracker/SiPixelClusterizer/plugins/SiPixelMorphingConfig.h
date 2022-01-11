#ifndef RecoLocalTracker_SiPixelClusterizer_plugins_SiPixelMorphingConfig_h
#define RecoLocalTracker_SiPixelClusterizer_plugins_SiPixelMorphingConfig_h

#include <cstdint>
#include <vector>

struct SiPixelMorphingConfig {
  int32_t nrows_;
  int32_t ncols_;
  int32_t nrocs_;
  int32_t iters_;
  std::vector<int32_t> kernel1_;
  std::vector<int32_t> kernel2_;
  uint32_t fakeAdc_;
};

#endif  // RecoLocalTracker_SiPixelClusterizer_plugins_SiPixelMorphingConfig_h
