#ifndef LORENTZANGLE_RUNINFO_ANALYZER
#define LORENTZANGLE_RUNINFO_ANALYZER

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "TTree.h"

class LorentzAngleRunInfo : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>
{
   public:
      explicit LorentzAngleRunInfo(const edm::ParameterSet&);
      ~LorentzAngleRunInfo();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void beginRun(const edm::Run &, const edm::EventSetup &);
      void endRun(const edm::Run &, const edm::EventSetup &);
      
      TTree * runinfo;
      unsigned int                 run;
      std::vector<unsigned int>    rawid               ;
      std::vector<float>           BdotY               ;
      std::vector<float>           localB              ;
      std::vector<float>           globalZofunitlocalY ;
      std::vector<float>           lorentzAngle        ;
      std::vector<float>           driftx              ;
      std::vector<float>           drifty              ;
      std::vector<float>           driftz              ;
      // pointer
      std::vector<unsigned int> * prawid               ;
      std::vector<float>        * pBdotY               ;
      std::vector<float>        * plocalB              ;
      std::vector<float>        * pglobalZofunitlocalY ;
      std::vector<float>        * plorentzAngle        ;
      std::vector<float>        * pdriftx              ;
      std::vector<float>        * pdrifty              ;
      std::vector<float>        * pdriftz              ;

};
#endif



