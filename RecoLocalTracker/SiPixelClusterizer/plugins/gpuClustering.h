#ifndef RecoLocalTracker_SiPixelClusterizer_plugins_gpuClustering_h
#define RecoLocalTracker_SiPixelClusterizer_plugins_gpuClustering_h

#include <cstdint>
#include <cstdio>

#include "CUDADataFormats/SiPixelCluster/interface/gpuClusteringConstants.h"
#include "CUDADataFormats/SiPixelDigi/interface/SiPixelDigisCUDASOAView.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "HeterogeneousCore/CUDAUtilities/interface/HistoContainer.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cuda_assert.h"

namespace gpuClustering {

  constexpr uint32_t maxPixInModule = 6000;

  template <bool IsPhase2>
  constexpr int32_t nbins = 1024;

  template <>
  constexpr int32_t nbins<true> = 1024;

  template <>
  constexpr int32_t nbins<false> = phase1PixelTopology::numColsInModule + 2;

  template <bool IsPhase2>
  constexpr int32_t nbits = 10;

  template <>
  constexpr int32_t nbits<true> = 10;

  template <>
  constexpr int32_t nbits<false> = 9;

  template <bool IsPhase2, bool IsFakeHist = false>
  using Hist = cms::cuda::HistoContainer<uint16_t,
                                         nbins<IsPhase2>,
                                         maxPixInModule,
                                         nbits<IsPhase2>,
                                         typename std::conditional<IsFakeHist, uint32_t, uint16_t>::type>;

#ifdef GPU_DEBUG
  __device__ uint32_t gMaxHit = 0;
#endif

  template <bool isPhase2>
  __global__ void countModules(uint16_t const* __restrict__ id,
                               uint32_t* __restrict__ moduleStart,
                               int32_t* __restrict__ clusterId,
                               int numElements) {
    int first = blockDim.x * blockIdx.x + threadIdx.x;
    constexpr int nMaxModules = isPhase2 ? phase2PixelTopology::numberOfModules : phase1PixelTopology::numberOfModules;
    assert(nMaxModules < maxNumModules);
    for (int i = first; i < numElements; i += gridDim.x * blockDim.x) {
      clusterId[i] = i;
      if (invalidModuleId == id[i])
        continue;
      auto j = i - 1;
      while (j >= 0 and id[j] == invalidModuleId)
        --j;
      if (j < 0 or id[j] != id[i]) {
        // boundary...
        auto loc = atomicInc(moduleStart, nMaxModules);
        moduleStart[loc + 1] = i;
      }
    }
  }

  template <bool isPhase2, bool isFakeHist>
  __device__ void prepareHist(uint16_t const* __restrict__ id,
                              uint16_t const* __restrict__ y,
                              Hist<isPhase2, isFakeHist>* hist,
                              typename Hist<isPhase2, isFakeHist>::Counter* ws,
                              int numElements,
                              int idOffset,
                              uint16_t thisModuleId,
                              uint32_t firstPixel) {
    auto first = isFakeHist ? threadIdx.x : (firstPixel + threadIdx.x);

#ifdef GPU_DEBUG
    __shared__ uint32_t totGood;
    totGood = 0;
    __syncthreads();
#endif

    for (auto j = threadIdx.x; j < Hist<isPhase2, isFakeHist>::totbins(); j += blockDim.x) {
      hist->off[j] = 0;
    }
    __syncthreads();

    // fill histo
    for (int i = first; i < numElements; i += blockDim.x) {
      if (id[i] == invalidModuleId)  // skip invalid pixels
        continue;
      if (isFakeHist && id[i] != thisModuleId)  // for fake histograms we need further filtering
        continue;
      hist->count(y[i]);
#ifdef GPU_DEBUG
      atomicAdd(&totGood, 1);
#endif
    }
    __syncthreads();
    if (threadIdx.x < 32)
      ws[threadIdx.x] = 0;  // used by prefix scan...
    __syncthreads();
    hist->finalize(ws);
    __syncthreads();
#ifdef GPU_DEBUG
    assert(hist->size() == totGood);
    if (thisModuleId % 100 == 1)
      if (threadIdx.x == 0)
        printf("histo size %d\n", hist->size());
#endif
    for (int i = first; i < numElements; i += blockDim.x) {
      if (id[i] == invalidModuleId)  // skip invalid pixels
        continue;
      if (isFakeHist && id[i] != thisModuleId)  // for fake histograms we need further filtering
        continue;
      assert(id[i] == thisModuleId);
      hist->fill(y[i], i + idOffset);
    }
    __syncthreads();  // for hit filling!
  }

  template <bool isPhase2, bool fromFake, bool toFake, class NNType>
  __device__ void fillNearestNeighbours(uint16_t const* __restrict__ fromid,
                                        uint16_t const* __restrict__ fromx,
                                        uint16_t const* __restrict__ fromy,
                                        uint16_t const* __restrict__ tox,
                                        uint16_t const* __restrict__ toy,
                                        Hist<isPhase2, fromFake>* fromHist,
                                        Hist<isPhase2, toFake>* toHist,
                                        int maxiter,
                                        int maxNeighbours,
                                        int idFromOffset,
                                        int idToOffset,
                                        uint16_t thisModuleId,
                                        NNType nf[16][10],
                                        uint8_t nnf[16],
                                        unsigned int& k) {
    for (auto j = threadIdx.x; j < fromHist->size(); j += blockDim.x, ++k) {
      assert(int(k) < maxiter);
      auto p = fromHist->begin() + j;
      auto i = *p + idFromOffset;
      assert(fromid[i] == thisModuleId);  // same module
      auto toHistYBeigin = toHist->begin(Hist<isPhase2, toFake>::bin(fromy[i]));

      using typePtr = typename std::conditional<fromFake == toFake, decltype(p), decltype(toHistYBeigin)>::type;
      typePtr pp;
      if constexpr (fromFake == toFake)
        pp = ++p;
      else
        pp = toHistYBeigin;

      int be = Hist<isPhase2, toFake>::bin(fromy[i] + 1);
      auto e = toHist->end(be);
      assert(0 == nnf[k]);
      for (; pp < e; ++pp) {
        auto m = (*pp) + idToOffset;
        // assert(m != i);
        assert(int(toy[m]) - int(fromy[i]) >= 0);
        assert(int(toy[m]) - int(fromy[i]) <= 1);
        if (std::abs(int(tox[m]) - int(fromx[i])) > 1)
          continue;
        auto l = nnf[k]++;
        assert(int(l) < maxNeighbours);
        nf[k][l] = *pp;
      }
    }
  }

  template <bool isPhase2>
  __global__ void findClus(uint16_t const* __restrict__ id,           // module id of each pixel
                           uint16_t const* __restrict__ x,            // local coordinates of each pixel
                           uint16_t const* __restrict__ y,            //
                           uint32_t const* __restrict__ moduleStart,  // index of the first pixel of each module
                           uint32_t* __restrict__ nClustersInModule,  // output: number of clusters found in each module
                           uint32_t* __restrict__ moduleId,           // output: module id of each module
                           int32_t* __restrict__ clusterId,           // output: cluster id of each pixel
                           int numElements,
                           SiPixelDigisCUDASOAView const fakeDigisView = SiPixelDigisCUDASOAView{},
                           int* numFakes = nullptr) {
    __shared__ int msize;

    auto firstModule = blockIdx.x;
    auto endModule = moduleStart[0];

    constexpr int nMaxModules = isPhase2 ? phase2PixelTopology::numberOfModules : phase1PixelTopology::numberOfModules;
    assert(nMaxModules < maxNumModules);

    for (auto module = firstModule; module < endModule; module += gridDim.x) {
      auto firstPixel = moduleStart[1 + module];
      auto thisModuleId = id[firstPixel];
      assert(thisModuleId < nMaxModules);

#ifdef GPU_DEBUG
      if (thisModuleId % 100 == 1)
        if (threadIdx.x == 0)
          printf("start clusterizer for module %d in block %d\n", thisModuleId, blockIdx.x);
#endif

      auto first = firstPixel + threadIdx.x;

      // find the index of the first pixel not belonging to this module (or invalid)
      msize = numElements;
      __syncthreads();

      // skip threads not associated to an existing pixel
      for (int i = first; i < numElements; i += blockDim.x) {
        if (id[i] == invalidModuleId)  // skip invalid pixels
          continue;
        if (id[i] != thisModuleId) {  // find the first pixel in a different module
          atomicMin(&msize, i);
          break;
        }
      }

      //init hist  (ymax=416 < 512 : 9bits)
      __shared__ Hist<isPhase2, false> hist;
      __shared__ Hist<isPhase2, true> fakeHist;
      __shared__ typename Hist<isPhase2, false>::Counter ws[32];
      __shared__ typename Hist<isPhase2, true>::Counter ws2[32];

      bool useFakeDigisForClustering = (numFakes != nullptr && *numFakes != 0);

      if (useFakeDigisForClustering) {
        prepareHist<isPhase2, true>(
            fakeDigisView.moduleInd(), fakeDigisView.yy(), &fakeHist, ws2, *numFakes, 0, thisModuleId, 0);
      }

      assert((msize == numElements) or ((msize < numElements) and (id[msize] != thisModuleId)));

      // limit to maxPixInModule  (FIXME if recurrent (and not limited to simulation with low threshold) one will need to implement something cleverer)
      if (0 == threadIdx.x) {
        if (msize - firstPixel > maxPixInModule) {
          printf("too many pixels in module %d: %d > %d\n", thisModuleId, msize - firstPixel, maxPixInModule);
          msize = maxPixInModule + firstPixel;
        }
      }

      __syncthreads();
      assert(msize - firstPixel <= maxPixInModule);

      prepareHist<isPhase2, false>(id, y, &hist, ws, msize, -int(firstPixel), thisModuleId, firstPixel);

#ifdef __CUDA_ARCH__
      // assume that we can cover the whole module with up to 16 blockDim.x-wide iterations
      constexpr int maxiter = 16;
      if (threadIdx.x == 0 && (hist.size() / blockDim.x) >= maxiter)
        printf("THIS IS NOT SUPPOSED TO HAPPEN too many hits in module %d: %d for block size %d\n",
               thisModuleId,
               hist.size(),
               blockDim.x);
#else
      auto maxiter = hist.size();
#endif
      // allocate space for duplicate pixels: a pixel can appear more than once with different charge in the same event
      constexpr int maxNeighbours = 10;
      assert((hist.size() / blockDim.x) <= maxiter);
      // nearest neighbour
      uint16_t nn[maxiter][maxNeighbours];
      uint8_t nnn[maxiter];  // number of nn
      // nearest fake neighbours
      uint32_t nf[maxiter][maxNeighbours];
      uint8_t nnf[maxiter];  // number of nf
      for (uint32_t k = 0; k < maxiter; ++k) {
        nnn[k] = 0;
        nnf[k] = 0;
      }

      __syncthreads();  // for hit filling!

#ifdef GPU_DEBUG
      // look for anomalous high occupancy
      __shared__ uint32_t n40, n60;
      n40 = n60 = 0;
      __syncthreads();
      for (auto j = threadIdx.x; j < Hist<isPhase2, false>::nbins(); j += blockDim.x) {
        if (hist.size(j) > 60)
          atomicAdd(&n60, 1);
        if (hist.size(j) > 40)
          atomicAdd(&n40, 1);
      }
      __syncthreads();
      if (0 == threadIdx.x) {
        if (n60 > 0)
          printf("columns with more than 60 px %d in %d\n", n60, thisModuleId);
        else if (n40 > 0)
          printf("columns with more than 40 px %d in %d\n", n40, thisModuleId);
      }
      __syncthreads();
#endif

      // fill NN
      auto k = 0U;
      for (auto j = threadIdx.x; j < hist.size(); j += blockDim.x, ++k) {
        assert(k < maxiter);
        auto p = hist.begin() + j;
        auto i = *p + firstPixel;
        assert(id[i] != invalidModuleId);
        assert(id[i] == thisModuleId);  // same module
        int be = Hist<isPhase2, false>::bin(y[i] + 1);
        auto e = hist.end(be);
        ++p;
        assert(0 == nnn[k]);
        for (; p < e; ++p) {
          auto m = (*p) + firstPixel;
          assert(m != i);
          assert(int(y[m]) - int(y[i]) >= 0);
          assert(int(y[m]) - int(y[i]) <= 1);
          if (std::abs(int(x[m]) - int(x[i])) > 1)
            continue;
          auto l = nnn[k]++;
          assert(l < maxNeighbours);
          nn[k][l] = *p;
        }
      }

      if (useFakeDigisForClustering) {
        fillNearestNeighbours<isPhase2, true, false>(fakeDigisView.moduleInd(),
                                                     fakeDigisView.xx(),
                                                     fakeDigisView.yy(),
                                                     x,
                                                     y,
                                                     &fakeHist,
                                                     &hist,
                                                     maxiter,
                                                     maxNeighbours,
                                                     0,
                                                     -int(firstPixel),
                                                     thisModuleId,
                                                     nn,
                                                     nnn,
                                                     k);
        k = 0U;
        fillNearestNeighbours<isPhase2, false, true>(id,
                                                     x,
                                                     y,
                                                     fakeDigisView.xx(),
                                                     fakeDigisView.yy(),
                                                     &hist,
                                                     &fakeHist,
                                                     maxiter,
                                                     maxNeighbours,
                                                     -int(firstPixel),
                                                     0,
                                                     thisModuleId,
                                                     nf,
                                                     nnf,
                                                     k);
        fillNearestNeighbours<isPhase2, true, true>(fakeDigisView.moduleInd(),
                                                    fakeDigisView.xx(),
                                                    fakeDigisView.yy(),
                                                    fakeDigisView.xx(),
                                                    fakeDigisView.yy(),
                                                    &fakeHist,
                                                    &fakeHist,
                                                    maxiter,
                                                    maxNeighbours,
                                                    0,
                                                    0,
                                                    thisModuleId,
                                                    nf,
                                                    nnf,
                                                    k);
      }

      // for each pixel, look at all the pixels until the end of the module;
      // when two valid pixels within +/- 1 in x or y are found, set their id to the minimum;
      // after the loop, all the pixel in each cluster should have the id equeal to the lowest
      // pixel in the cluster ( clus[i] == i ).
      bool more = true;
      int nloops = 0;
      if (!useFakeDigisForClustering) {
        while (__syncthreads_or(more)) {
          if (1 == nloops % 2) {
            for (auto j = threadIdx.x, k = 0U; j < hist.size(); j += blockDim.x, ++k) {
              auto p = hist.begin() + j;
              auto i = *p + firstPixel;
              auto m = clusterId[i];
              while (m != clusterId[m])
                m = clusterId[m];
              clusterId[i] = m;
            }
          } else {
            more = false;
            for (auto j = threadIdx.x, k = 0U; j < hist.size(); j += blockDim.x, ++k) {
              auto p = hist.begin() + j;
              auto i = *p + firstPixel;
              for (int kk = 0; kk < nnn[k]; ++kk) {
                auto l = nn[k][kk];
                auto m = l + firstPixel;
                assert(m != i);
                auto old = atomicMin_block(&clusterId[m], clusterId[i]);
                // do we need memory fence?
                if (old != clusterId[i]) {
                  // end the loop only if no changes were applied
                  more = true;
                }
                atomicMin_block(&clusterId[i], old);
              }  // nnloop
            }    // pixel loop
          }
          ++nloops;
        }  // end while
      } else {
        while (__syncthreads_or(more)) {
          if (1 == nloops % 2) {
            for (auto j = threadIdx.x, k = 0U; j < hist.size(); j += blockDim.x, ++k) {
              auto p = hist.begin() + j;
              auto i = *p + firstPixel;
              auto m = clusterId[i];
              while (m != clusterId[m])
                m = clusterId[m];
              clusterId[i] = m;
            }
          } else {
            more = false;
            for (auto j = threadIdx.x, k = 0U; j < hist.size(); j += blockDim.x, ++k) {
              auto p = hist.begin() + j;
              auto i = *p + firstPixel;
              for (int kk = 0; kk < nnn[k]; ++kk) {
                auto l = nn[k][kk];
                auto m = l + firstPixel;
                assert(m != i);
                auto old = atomicMin_block(&clusterId[m], clusterId[i]);
                // do we need memory fence?
                if (old != clusterId[i]) {
                  // end the loop only if no changes were applied
                  more = true;
                }
                atomicMin_block(&clusterId[i], old);
              }  // nnloop
            }    // pixel loop
          }
          ++nloops;
        }  // end while
      }

#ifdef GPU_DEBUG
      {
        __shared__ int n0;
        if (threadIdx.x == 0)
          n0 = nloops;
        __syncthreads();
        auto ok = n0 == nloops;
        assert(__syncthreads_and(ok));
        if (thisModuleId % 100 == 1)
          if (threadIdx.x == 0)
            printf("# loops %d\n", nloops);
      }
#endif

      __shared__ unsigned int foundClusters;
      foundClusters = 0;
      __syncthreads();

      // find the number of different clusters, identified by a pixels with clus[i] == i;
      // mark these pixels with a negative id.
      for (int i = first; i < msize; i += blockDim.x) {
        if (id[i] == invalidModuleId)  // skip invalid pixels
          continue;
        if (clusterId[i] == i) {
          auto old = atomicInc(&foundClusters, 0xffffffff);
          clusterId[i] = -(old + 1);
        }
      }
      __syncthreads();

      // propagate the negative id to all the pixels in the cluster.
      for (int i = first; i < msize; i += blockDim.x) {
        if (id[i] == invalidModuleId)  // skip invalid pixels
          continue;
        if (clusterId[i] >= 0) {
          // mark each pixel in a cluster with the same id as the first one
          clusterId[i] = clusterId[clusterId[i]];
        }
      }
      __syncthreads();

      // adjust the cluster id to be a positive value starting from 0
      for (int i = first; i < msize; i += blockDim.x) {
        if (id[i] == invalidModuleId) {  // skip invalid pixels
          clusterId[i] = invalidClusterId;
          continue;
        }
        clusterId[i] = -clusterId[i] - 1;
      }
      __syncthreads();

      if (threadIdx.x == 0) {
        nClustersInModule[thisModuleId] = foundClusters;
        moduleId[module] = thisModuleId;
#ifdef GPU_DEBUG
        if (foundClusters > gMaxHit) {
          gMaxHit = foundClusters;
          if (foundClusters > 8)
            printf("max hit %d in %d\n", foundClusters, thisModuleId);
        }
#endif
#ifdef GPU_DEBUG
        if (thisModuleId % 100 == 1)
          printf("%d clusters in module %d\n", foundClusters, thisModuleId);
#endif
      }
    }  // module loop
  }
}  // namespace gpuClustering

#endif  // RecoLocalTracker_SiPixelClusterizer_plugins_gpuClustering_h
