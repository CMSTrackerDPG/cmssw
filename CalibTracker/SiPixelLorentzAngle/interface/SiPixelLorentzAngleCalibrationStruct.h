#ifndef CalibTracker_SiPixelLorentzAngle_SiPixelLorentzAngleCalibrationStruct_h
#define CalibTracker_SiPixelLorentzAngle_SiPixelLorentzAngleCalibrationStruct_h

#include "DQMServices/Core/interface/DQMStore.h"
#include <unordered_map>

struct SiPixelLorentzAngleCalibrationHistograms {
public:
  SiPixelLorentzAngleCalibrationHistograms() = default;

  using MonitorMap = std::unordered_map<uint32_t, dqm::reco::MonitorElement*>;

  int nlay;
  std::vector<int> nModules_;
  std::vector<int> newDetIds_;
  std::vector<int> newModule_;
  std::vector<int> newLayer_;
  std::map<uint32_t, std::vector<uint32_t> > detIdsList;

  MonitorMap h_drift_depth_adc_;
  MonitorMap h_drift_depth_adc2_;
  MonitorMap h_drift_depth_noadc_;
  MonitorMap h_drift_depth_;
  MonitorMap h_mean_;

  dqm::reco::MonitorElement* h_tracks_;
};

#endif
