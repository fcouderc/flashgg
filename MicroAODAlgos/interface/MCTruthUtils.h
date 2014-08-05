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

  struct phoMCTruthInfo {
    float  eMC;
    int    mcMatching;  /// -1 - not matched; 0 - matched; 2 - isolated photon; 3 - photon from resonance (Higgs only for now)
    float  eIso;
    float  eRelIso;
  };

  class MCTruthUtils {
  public:
    MCTruthUtils(void);
    ~MCTruthUtils(void) {};
    
    /// return index of closest to the truth recoVertex. If no vertex is matched then return -1
    int getMCTruthVertexIndex( const edm::Event &event, const edm::PtrVector<reco::Vertex>&, double dzMatch = 0.1);
    phoMCTruthInfo getMCMatchingPhoton(   const edm::Event &event, const pat::Photon &p );

  private:
    edm::InputTag genParticlesTag_;
    edm::InputTag stableGenParticleTag_;

  };



}


#endif
