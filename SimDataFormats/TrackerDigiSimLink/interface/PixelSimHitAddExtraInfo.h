#ifndef PixelSimHitAddExtraInfo_h
#define PixelSimHitAddExtraInfo_h

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include <vector>
#include <cstdint>

class PixelSimHitAddExtraInfo {
public:
  PixelSimHitAddExtraInfo(size_t Hindex, Local3DPoint entryP, Local3DPoint exitP, unsigned int ch) {
    index_ = Hindex;
    theEntryPoint_ = entryP;
    theExitPoint_ = exitP;
    chan_.push_back(ch);
  };
  PixelSimHitAddExtraInfo(){};
  ~PixelSimHitAddExtraInfo(){};
  size_t hitIndex() const { return index_; };
  Local3DPoint entryPoint() const { return theEntryPoint_; };
  Local3DPoint exitPoint() const { return theExitPoint_; }
  std::vector<unsigned int> channel() const { return chan_; };

  inline bool operator<(const PixelSimHitAddExtraInfo& other) const { return hitIndex() < other.hitIndex(); }

  void addDigiInfo(unsigned int theDigiChannel) { chan_.push_back(theDigiChannel); }
  bool isInTheList(unsigned int channelToCheck) {
    bool result_in_the_list = false;
    for (unsigned int icheck = 0; icheck < chan_.size(); icheck++) {
      if (channelToCheck == chan_[icheck])
        result_in_the_list = true;
    }
    return result_in_the_list;
  }

private:
  size_t index_;
  Local3DPoint theEntryPoint_;
  Local3DPoint theExitPoint_;
  std::vector<unsigned int> chan_;
};
#endif
