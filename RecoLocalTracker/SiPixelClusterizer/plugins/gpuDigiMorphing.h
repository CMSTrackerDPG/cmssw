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

int getUpperBoundForFlaggedDigis(uint32_t wordCounter, SiPixelMorphingConfig const& c) {
  return wordCounter * countKernelOverlap(c);
}

#endif  // RecoLocalTracker_SiPixelClusterizer_plugins_gpuDigiMorhping_h
