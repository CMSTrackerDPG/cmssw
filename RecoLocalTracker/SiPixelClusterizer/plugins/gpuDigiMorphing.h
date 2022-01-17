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
    divide the bits in kernel dilate by 2
     - e.g. 24 / 2 = 12
     x  x  x  x  x
     x  x  x  x  x
     x  x  O  x  x
     x  x  x  x  x
     x  x  x  x  x
    **/
  int bitsSet = 0;
  for (auto i = 0u; i < c.kernel1_.size(); ++i) {
    bitsSet += countBitsSet(c.kernel1_[i]);
  }
  assert(bitsSet > 0);
  return (bitsSet - 1)/ 2;
}

int getUpperBoundForFakeDigis(uint32_t wordCounter, SiPixelMorphingConfig const& c) {
  return wordCounter * countKernelOverlap(c);
}

int getKernelSizeFromConfig(SiPixelMorphingConfig const& digiMorphingConfig) {
  return digiMorphingConfig.kernel1_.size();
}

auto constructMorphingKernelsFromConfig(SiPixelMorphingConfig const& digiMorphingConfig, cudaStream_t stream) {
  auto kernelSize = getKernelSizeFromConfig(digiMorphingConfig);
  auto twoKernelSize = 2 * kernelSize * kernelSize;
  auto kernels_d = cms::cuda::make_device_unique<int[]>(twoKernelSize, stream);
  std::vector<int> kernels_h(twoKernelSize, 0);

  for (int i = 0; i < kernelSize; ++i) {
    int row_d = digiMorphingConfig.kernel1_[i];
    int row_e = digiMorphingConfig.kernel2_[i];
    for (int j = 0; j < kernelSize; ++j) {
      kernels_h[i * kernelSize + (kernelSize - j - 1)] = (row_d % 2);
      kernels_h[kernelSize * kernelSize + i * kernelSize + (kernelSize - j - 1)] = (row_e % 2);
      row_d /= 2;
      row_e /= 2;
    }
  }

  cudaCheck(cudaMemcpyAsync(kernels_d.get(), kernels_h.data(), twoKernelSize * sizeof(int), cudaMemcpyDefault, stream));

  return kernels_d;
}

namespace gpudigimorphing {
  using FLAG_KERNEL_TYPE = unsigned long long;
  const int FLAG_TYPE_BITS = sizeof(FLAG_KERNEL_TYPE) * 8;

  const int divideModuleCols = 8;
  const int divideModuleRows = 2;
  const int moduleConvolutions = divideModuleCols * divideModuleRows;

  __device__ void binprintf(long v)
{
    uint64_t mask=(long)(1)<<((sizeof(long)<<3)-1);
    while(mask) {
        printf("%d", (v&mask ? 1 : 0));
        mask >>= 1;
    }
    printf("\n");
}

__device__ void printModule(FLAG_KERNEL_TYPE* p){
  for(int i=0;i<82;++i){
    binprintf(*(p+i));
  }
}

  inline __device__ int getIndex(int row, int col, int maxRows, int iters) { return row % maxRows + iters; }

  inline __device__ void setBit(FLAG_KERNEL_TYPE& num, int bit, int rocWidth) {
    int b = bit % rocWidth;
    b = (FLAG_TYPE_BITS - 1 - b);
    num |= 1UL << b;
  }

  __device__ void convolutionByBitManipulation(int const* kernel,
                                               SiPixelMorphingConfig const& morphingConfig,
                                               FLAG_KERNEL_TYPE* nonConvolutedPixels,
                                               FLAG_KERNEL_TYPE* convolutedPixels,
                                               int heightMin,
                                               int heightMax,
                                               int rocHeight,
                                               bool isDilate) {
    int kernelRadius = morphingConfig.iters_;
    int kernelSize = 2 * kernelRadius + 1;
    for (int i = threadIdx.x + heightMin; i < heightMax; i += blockDim.x) {
      int index = getIndex(i, 0, rocHeight, morphingConfig.iters_);
      FLAG_KERNEL_TYPE hits = 0;
      if(!isDilate) hits-=1; // get maximum value storable, full ones
      for (int r1 = -kernelRadius; r1 <= kernelRadius and (isDilate or hits > 0); ++r1) {
        for (int r2 = -kernelRadius; r2 <= kernelRadius and (isDilate or hits > 0); ++r2) {
          if (kernel[(1 + r1) * kernelSize + (1 + r2)] == 1) {
            auto newIndex = index + r1;
            // assert(newIndex >= 0 && newIndex < p.convolutionHeight);
            auto newRow = nonConvolutedPixels[newIndex];
            if (r2 < 0)
              newRow >>= std::abs(r2);
            else
              newRow <<= r2;
            if (isDilate)
              hits |= newRow;
            else
              hits &= newRow;
          }
        }
      }
      // if (hits)
      convolutedPixels[index] = hits;
    }
  }

  struct ROCParameters {
    int rocWidth;
    int rocHeight;
    int convolutionHeight;
  };

  __device__ ROCParameters computeROCParameters(SiPixelMorphingConfig const& morphingConfig) {
    return ROCParameters{.rocWidth = morphingConfig.ncols_ / divideModuleCols,
                         .rocHeight = morphingConfig.nrows_ / divideModuleRows,
                         .convolutionHeight = morphingConfig.nrows_ / divideModuleRows + 2 * morphingConfig.iters_};
  }

  __global__ void clusterHealingWithDigiMorphing_kernel(
      SiPixelDigisCUDASOAView const digisView,
      SiPixelDigisCUDASOAView fakeDigisView,
      SiPixelMorphingConfig const morphingConfig,
      int const* kernels,
      uint32_t* moduleStart,  // index of the first pixel of each module)
      int numDigis,
      int* fakeCounter) {
    assert(morphingConfig.ncols_ / divideModuleCols + 2 * morphingConfig.iters_ < FLAG_TYPE_BITS);

    auto p = computeROCParameters(morphingConfig);
    assert(p.rocWidth == 52);
    assert(p.rocHeight==80);
    assert(p.convolutionHeight==82);

#ifdef GPU_DEBUG
    const int width = FLAG_TYPE_BITS;  // 64
    assert(width % FLAG_TYPE_BITS == 0);
    if (blockIdx.x * blockDim.x + threadIdx.x == 0)
      printf("Start kernel clusterHealingWithDigiMorphing_kernel, width is %d\n", width);
#endif

    extern __shared__ FLAG_KERNEL_TYPE s[];  // shared memory has size width*height*2/8 bytes
    FLAG_KERNEL_TYPE* modulePixels = s;
    FLAG_KERNEL_TYPE* dilatedPixels = s + p.convolutionHeight;
    FLAG_KERNEL_TYPE* erodedPixels = s + 2 * p.convolutionHeight;
    int kernelRadius = morphingConfig.iters_;
    int kernelSize = 2 * kernelRadius + 1;
    auto kernelDilate = kernels;
    auto kernelErode = kernels + kernelSize * kernelSize;
    #ifdef GPU_DEBUG
    if(threadIdx.x + blockIdx.x*blockDim.x == 0){
      printf("%d%d%d\n", *kernelDilate, *(kernelDilate+1), *(kernelDilate+2));
      printf("%d%d%d\n", *(kernelDilate+3), *(kernelDilate+4), *(kernelDilate+5));
      printf("%d%d%d\n", *(kernelDilate+6), *(kernelDilate+7), *(kernelDilate+8));
    }
    #endif
    // set to zero
    for (int i = threadIdx.x; i < p.convolutionHeight; i += blockDim.x) {
      modulePixels[i] = dilatedPixels[i] = erodedPixels[i] = 0;
    }

    __syncthreads();

    auto firstModule = blockIdx.x;
    auto endModule = moduleStart[0];

    // iterate modules
    for (auto module = firstModule; module < endModule * moduleConvolutions; module += gridDim.x) {
      auto firstPixel = moduleStart[1 + (module / moduleConvolutions)];
      auto first = firstPixel + threadIdx.x;
      auto thisModuleId = digisView.moduleInd(firstPixel);

      int rocId = (module % moduleConvolutions);
      int widthMin = (rocId % divideModuleCols) * p.rocWidth;
      int widthMax = widthMin + p.rocWidth;
      int heightMin = (rocId / divideModuleCols) * p.rocHeight;
      int heightMax = heightMin + p.rocHeight;

      FLAG_KERNEL_TYPE num = 0;
      for (int i = first; i < numDigis; i += blockDim.x) {
        if (digisView.moduleInd(i) == gpuClustering::invalidModuleId)
          continue;
        if (digisView.moduleInd(i) != thisModuleId)
          break;
        assert(digisView.moduleInd(i) == thisModuleId);
        if (widthMin <= digisView.yy(i) && digisView.yy(i) < widthMax && heightMin <= digisView.xx(i) &&
            digisView.xx(i) < heightMax) {
          num = 0;
          #ifdef GPU_DEBUG
          if (thisModuleId % 2000 ==1701) {
        if ((blockIdx.x % moduleConvolutions) == 8) {
          printf("Hit at %d row%d col%d\n", thisModuleId, digisView.xx(i), digisView.yy(i));
        }}
          #endif

          int index = getIndex(digisView.xx(i), digisView.yy(i), p.rocHeight, morphingConfig.iters_);
          setBit(num, digisView.yy(i) + morphingConfig.iters_, p.rocWidth);
          atomicAdd(modulePixels + index, num);
        }
      }
      __syncthreads();

#ifdef GPU_DEBUG
      // print one pixel ROC of a module
      if (thisModuleId % 2000 == 1701) {
        if (threadIdx.x == 0 && (blockIdx.x % moduleConvolutions) == 8) {
          printf("ROC BEFORE\n");
          printModule(modulePixels);
        }
      }
      __syncthreads();
#endif

      // dilate
      convolutionByBitManipulation(
          kernelDilate, morphingConfig, modulePixels, dilatedPixels, heightMin, heightMax, p.rocHeight, true);
      __syncthreads();

      #ifdef GPU_DEBUG
      // print one pixel ROC of a module
      if (thisModuleId % 2000 ==1701) {
        if (threadIdx.x == 0 && (blockIdx.x % moduleConvolutions) == 8) {
          printf("ROC AFTER DILATE\n");
          printModule(dilatedPixels);
        }
      }
      __syncthreads();
#endif

      // erode
      convolutionByBitManipulation(
          kernelErode, morphingConfig, dilatedPixels, erodedPixels, heightMin, heightMax, p.rocHeight, false);
      __syncthreads();

#ifdef GPU_DEBUG
      // print one pixel ROC of a module
      if (thisModuleId % 2000 ==1701) {
        if (threadIdx.x == 0 && (blockIdx.x % moduleConvolutions) == 8) {
          printf("ROC AFTER ERODE\n");
          printModule(erodedPixels);
        }
      }
      __syncthreads();
#endif

      // compare morphed (erodedPixels) pixels and originals (modulePixels)
      for (int row = threadIdx.x + heightMin; row < heightMax; row += blockDim.x) {
        int col = widthMax - 1;
        int index = getIndex(row, col, p.rocHeight, morphingConfig.iters_);
        FLAG_KERNEL_TYPE hits = (erodedPixels[index] & (~modulePixels[index]));
        hits >>= (FLAG_TYPE_BITS - (morphingConfig.ncols_ / divideModuleCols + 2 * morphingConfig.iters_) +
                  morphingConfig.iters_);
        while (hits) {
          if (hits & 1)  // check last bit
          {
            assert(col >= 0 && col < morphingConfig.ncols_);
            assert(row >= 0 && row < morphingConfig.nrows_);
            int old = atomicAdd(fakeCounter, 1);
            fakeDigisView.moduleInd()[old] = thisModuleId;
            fakeDigisView.xx()[old] = row;
            fakeDigisView.yy()[old] = col;
            fakeDigisView.adc()[old] = morphingConfig.fakeAdc_ * 10;  // calibrate fake digis
            fakeDigisView.clus()[old] = old;
          }
          hits >>= 1;
          col--;
        }
      }
      __syncthreads();

    }  // for modules
  }
}  // namespace gpudigimorphing

#endif  // RecoLocalTracker_SiPixelClusterizer_plugins_gpuDigiMorhping_h
