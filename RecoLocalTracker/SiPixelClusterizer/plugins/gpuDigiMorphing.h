#ifndef RecoLocalTracker_SiPixelClusterizer_plugins_gpuDigiMorhping_h
#define RecoLocalTracker_SiPixelClusterizer_plugins_gpuDigiMorhping_h

#include <cstdint>
#include <cstdio>

#include "CUDADataFormats/SiPixelCluster/interface/gpuClusteringConstants.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cuda_assert.h"

//local include(s)
#include "SiPixelMorphingConfig.h"

int countBitsSet(unsigned int v) {
  /**
    count the number of bits set in v
    **/
  int c;  // c accumulates the total bits set in v
  for (c = 0; v; c++) {
    v &= v - 1;  // clear the least significant bit set
  }
  return c;
}

int countKernelOverlap(SiPixelMorphingConfig const& c) {
  /**
    to determine an upper bound for extra digis
    count the overlapping bits in the two kernels dilate and erode
     - e.g. dilate + erode ->
     x  x  x  x  x             x                   x
     x  x  x  x  x          x  x  x             x  x  x
     x  x  O  x  x   +   x  x  O  x  x  ->   x  x  O  x  x
     x  x  x  x  x          x  x  x             x  x  x
     x  x  x  x  x             x                   x
     - overlapping bits (Xs) 12
    **/
  assert(c.kernel1_.size() == c.kernel2_.size());
  int bitsSet = 0;
  for (auto i = 0u; i < c.kernel1_.size(); ++i) {
    bitsSet += countBitsSet(c.kernel1_[i] & c.kernel2_[i]);
  }
  assert(bitsSet > 0);
  // extract 1 for original digi in the middle
  return bitsSet - 1;
}

int getUpperBoundForFakeDigis(uint32_t wordCounter, SiPixelMorphingConfig const& c) {
  return wordCounter * countKernelOverlap(c);
}

int getKernelSizeFromConfig(SiPixelMorphingConfig const& digiMorphingConfig)
{
  return digiMorphingConfig.kernel1_.size();
}

std::vector<std::vector<int>> constructMorphingKernelsFromConfig(SiPixelMorphingConfig const& digiMorphingConfig)
{
  int kernelSize = getKernelSizeFromConfig(digiMorphingConfig);
  std::vector<int> kernel_dilate_h(kernelSize*kernelSize, 0);
  std::vector<int> kernel_erode_h(kernelSize*kernelSize, 0);
  for (int i = 0; i*i < kernelSize; ++i) {
    int row_d = digiMorphingConfig.kernel1_[i];
    int row_e = digiMorphingConfig.kernel2_[i];
    for (int j = 0; j*j < kernelSize; ++j) {
      kernel_dilate_h[i * kernelSize + (kernelSize - j - 1)] = (row_d % 2);
      kernel_erode_h[i * kernelSize + (kernelSize - j - 1)] = (row_e % 2);
      row_d /= 2;
      row_e /= 2;
    }
  }
  return {kernel_dilate_h, kernel_erode_h};
}

#endif  // RecoLocalTracker_SiPixelClusterizer_plugins_gpuDigiMorhping_h
