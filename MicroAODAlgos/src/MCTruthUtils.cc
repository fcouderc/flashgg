#include "flashgg/MicroAODAlgos/interface/MCTruthUtils.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace std;
using namespace flashgg;
using namespace edm;

MCTruthUtils::MCTruthUtils(void) : 
  genParticlesTag_( edm::InputTag("prunedGenParticles") )
 { }


int MCTruthUtils::getMCTruthVertexIndex( const Event & event, const PtrVector<reco::Vertex>& vertices, double dzMatch ) {
  Handle<std::vector<reco::GenParticle> > genParticlesHandle; event.getByLabel(genParticlesTag_, genParticlesHandle);
  
  reco::Vertex::Point hardVertex(0,0,0);
  if( genParticlesHandle.isValid() )
    for( vector<reco::GenParticle>::const_iterator imc = genParticlesHandle->begin(); imc != genParticlesHandle->end(); ++imc )  {
      if( fabs( imc->pdgId() ) < 10 || fabs(imc->pdgId() ) == 25 ) hardVertex.SetCoordinates(imc->vx(),imc->vy(),imc->vz());
    }
  
  
  int  ivMatch = 0;
  double dzMin = 999;
  for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {
    double dz = fabs( vertices[iv]->z() - hardVertex.z() );
    if( dz < dzMin ) {
      ivMatch = iv;
      dzMin   = dz;
    }
  }
  

  if( dzMin < dzMatch ) return ivMatch;
  
  return -1;
}

