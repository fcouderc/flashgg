#ifndef FLASHgg_VertexSelectorBase_h
#define FLASHgg_VertexSelectorBase_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "flashgg/MicroAODFormats/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "flashgg/MicroAODFormats/interface/VertexCandidateMap.h"

// Shamelessly patterned on https://github.com/cms-sw/cmssw/blob/CMSSW_7_2_X/RecoParticleFlow/PFProducer/interface/BlockElementLinkerBase.h

namespace flashgg {

  class VertexSelectorBase {

  public:
    VertexSelectorBase(const edm::ParameterSet& conf):
      _selectorName( conf.getParameter<std::string>("VertexSelectorName") ) { }

    VertexSelectorBase(const VertexSelectorBase& ) = delete;
    VertexSelectorBase& operator=(const VertexSelectorBase&) = delete;

    virtual edm::Ptr<reco::Vertex> select(const edm::Ptr<flashgg::Photon>&,const edm::Ptr<flashgg::Photon>&,const edm::PtrVector<reco::Vertex>&,
					  const VertexCandidateMap&,const edm::PtrVector<reco::Conversion>&) const = 0;
    
    const std::string& name() const { return _selectorName; }

  private:
    const std::string _selectorName;

      
  };
}

#include "FWCore/PluginManager/interface/PluginFactory.h"
typedef edmplugin::PluginFactory< flashgg::VertexSelectorBase* (const edm::ParameterSet&) > FlashggVertexSelectorFactory;

#endif
