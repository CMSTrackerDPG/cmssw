#include "CalibTracker/Records/interface/SiStripGainRcd.h"  

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "CalibTracker/SiStripCommon/interface/ShallowTools.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
#include "CondFormats/DataRecord/interface/SiStripLorentzAngleRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"

#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"

#include "CalibTracker/SiStripLorentzAngle/interface/LorentzAngleRunInfo.h"

using namespace edm;
using namespace reco;
using namespace std;

LorentzAngleRunInfo::LorentzAngleRunInfo(const edm::ParameterSet& iConfig)
{
   usesResource("TFileService");
   
   prawid                 = &rawid               ;
   pBdotY                 = &BdotY               ;
   plocalB                = &localB              ;
   pglobalZofunitlocalY   = &globalZofunitlocalY ;
   plorentzAngle          = &lorentzAngle        ;
   pdriftx                = &driftx              ;
   pdrifty                = &drifty              ;
   pdriftz                = &driftz              ;
  
}

LorentzAngleRunInfo::~LorentzAngleRunInfo()
{
}


void LorentzAngleRunInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}

// ------------ method called when starting to processes a run  ------------
void LorentzAngleRunInfo::beginRun(const edm::Run & irun, const edm::EventSetup & iSetup)
{
   run = irun.run();
   edm::ESHandle<TrackerGeometry> theTrackerGeometry;         iSetup.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry );  
   edm::ESHandle<MagneticField> magfield;                     iSetup.get<IdealMagneticFieldRecord>().get(magfield);                      
   edm::ESHandle<SiStripLorentzAngle> SiStripLorentzAngle;    iSetup.get<SiStripLorentzAngleDepRcd>().get(SiStripLorentzAngle);      

   auto dets = theTrackerGeometry -> detsTIB();
   dets.insert(dets.end(),(theTrackerGeometry -> detsTOB()).begin(),(theTrackerGeometry -> detsTOB()).end());
   for ( auto det : dets )
   {
      float    bdoty = -999;
      float    locb  = -999;
      float    globalzy = -999;
      float    larcd = -999;
      float    drfx = -999;
      float    drfy = -999;
      float    drfz = -999;
   
      const StripGeomDetUnit* theStripDet = dynamic_cast<const StripGeomDetUnit*>( theTrackerGeometry->idToDet( det->geographicalId () ) );
      int detid = det->geographicalId().rawId();
      if ( theStripDet ) 
      {
         bdoty      = (theStripDet->surface()).toLocal( magfield->inTesla(theStripDet->surface().position())).y();
         locb       = magfield->inTesla(theStripDet->surface().position()).mag() ;
         globalzy   = (theStripDet->toGlobal(LocalVector(0,1,0))).z();
         larcd      = SiStripLorentzAngle -> getLorentzAngle(detid);
         LocalVector drift = shallow::drift( theStripDet, *magfield, *SiStripLorentzAngle);
         drfx   = drift.x();
         drfy   = drift.y();
         drfz   = drift.z();
         // fill the vectors
         prawid ->push_back( detid );         
         pBdotY->push_back( bdoty ); 
         plocalB->push_back( locb );
         pdriftx->push_back( drfx );
         pdrifty->push_back( drfy );
         pdriftz->push_back( drfz );
         pglobalZofunitlocalY->push_back( globalzy );
         plorentzAngle->push_back( larcd );
      }
   }
   
   runinfo->Fill();


}
// ------------ method called when ending the processing of a irun  ------------
void LorentzAngleRunInfo::endRun(const edm::Run & irun, const edm::EventSetup & iSetup)
{
}

void LorentzAngleRunInfo::beginJob()
{

  edm::Service<TFileService> fs;

  runinfo = fs->make<TTree>("tree","tree with run info for LA");
  
  runinfo -> Branch("run",&run,"run/i");
  // vectors
  runinfo -> Branch("rawid"              ,"std::vector<unsigned int>",&prawid              );
  runinfo -> Branch("BdotY"              ,"std::vector<float>"       ,&pBdotY              );
  runinfo -> Branch("localB"             ,"std::vector<float>"       ,&plocalB             );
  runinfo -> Branch("globalZofunitlocalY","std::vector<float>"       ,&pglobalZofunitlocalY);
  runinfo -> Branch("lorentzAngle"       ,"std::vector<float>"       ,&plorentzAngle       );
  runinfo -> Branch("driftx"             ,"std::vector<float>"       ,&pdriftx             );
  runinfo -> Branch("drifty"             ,"std::vector<float>"       ,&pdrifty             );
  runinfo -> Branch("driftz"             ,"std::vector<float>"       ,&pdriftz             );
  
}

void LorentzAngleRunInfo::endJob()
{
}
