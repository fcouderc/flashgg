#ifndef FLASHgg_MCTruthUtils_h
#define FLASHgg_MCTruthUtils_h

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "flashgg/MicroAODFormats/interface/Photon.h"
#include "flashgg/MicroAODFormats/interface/VertexCandidateMap.h"



namespace flashgg {
  class MCTruthUtils {
  public:
    MCTruthUtils(void);
    ~MCTruthUtils(void) {};
    
    /// return index of closest to the truth recoVertex. If no vertex is matched then return -1
    int getMCTruthVertexIndex( const edm::Event & event, const edm::PtrVector<reco::Vertex>&, double dzMatch = 0.1);

  private:
    edm::InputTag genParticlesTag_;

  };



}


#endif
