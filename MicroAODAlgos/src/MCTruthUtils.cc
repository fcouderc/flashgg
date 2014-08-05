#include "flashgg/MicroAODAlgos/interface/MCTruthUtils.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

using namespace std;
using namespace flashgg;
using namespace edm;

MCTruthUtils::MCTruthUtils(void) : 
  genParticlesTag_( edm::InputTag("prunedGenParticles") ),
  stableGenParticleTag_( edm::InputTag("packedGenParticles"))
 { }


int MCTruthUtils::getMCTruthVertexIndex( const Event & event, const PtrVector<reco::Vertex>& vertices, double dzMatch ) {
  Handle<std::vector<reco::GenParticle> > genParticlesHandle; event.getByLabel(genParticlesTag_, genParticlesHandle);
  
  reco::Vertex::Point hardVertex(0,0,0);
  if( genParticlesHandle.isValid() )
  for( vector<reco::GenParticle>::const_iterator imc = genParticlesHandle->begin(); imc != genParticlesHandle->end(); ++imc )  {
      if( fabs( imc->pdgId() ) < 10 || fabs(imc->pdgId() ) == 25 ) {
	hardVertex.SetCoordinates(imc->vx(),imc->vy(),imc->vz());
	break;
      }
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


phoMCTruthInfo MCTruthUtils::getMCMatchingPhoton( const Event &event, const pat::Photon &p ) {
  
  
  ///--- hard vertex  
  reco::Vertex::Point hardVertex(0,0,0);
  Handle<std::vector<reco::GenParticle> > genParticlesHandle1; event.getByLabel(genParticlesTag_, genParticlesHandle1);
  if( genParticlesHandle1.isValid() )
  for( vector<reco::GenParticle>::const_iterator imc = genParticlesHandle1->begin(); imc != genParticlesHandle1->end(); ++imc )  {
      if( fabs( imc->pdgId() ) < 10 || fabs(imc->pdgId() ) == 25 ) {
	hardVertex.SetCoordinates(imc->vx(),imc->vy(),imc->vz());
	break;
      }
  }
  
  
  math::XYZVector pPho( p.superCluster()->x() - hardVertex.x(),
			p.superCluster()->y() - hardVertex.y(),
			p.superCluster()->z() - hardVertex.z() 
			);
  
  
  ///--- match the photon
  Handle<std::vector<pat::PackedGenParticle> > genParticlesHandle2; event.getByLabel(stableGenParticleTag_, genParticlesHandle2);
  double matchedPhotonE(-1),matchedPhotonEta(-1),matchedPhotonPhi(-1);
  bool hasHiggsAncestor = false; 
  double drMin = 9999;

  phoMCTruthInfo result;
  result.mcMatching = -1;
  result.eMC        = -1;
  result.eIso       = -1;
  result.eRelIso    = -1;

  cout << "================================================ " << endl;
  cout << " Photon Pt: " << p.pt() << " eta: " << pPho.Eta() << " phi: " << pPho.Phi() << endl;
  if(  genParticlesHandle2.isValid() ) {
    ///--- first do matching
    for( vector<reco::GenParticle>::const_iterator imc = genParticlesHandle1->begin(); imc != genParticlesHandle1->end(); ++imc )
     if(  imc->pdgId() == 22 && imc->status() == 1 && imc->pt() > 5 ) {
       /// stable photon coming from hard vertex
       double dR = deltaR( pPho.Eta(), pPho.Phi(), imc->eta(), imc->phi() );
       if( dR < drMin ) {
	 drMin = dR; 
	 matchedPhotonE   = imc->energy();	 
	 matchedPhotonEta = imc->eta();	 
	 matchedPhotonPhi = imc->phi();	 
	 hasHiggsAncestor = false;
	 for( unsigned imoth = 0 ; imoth < imc->numberOfMothers(); imoth++ ) {
	   if     ( imc->mother(imoth)->pdgId() == 25 ) hasHiggsAncestor = true; 
	   else if( imc->mother(imoth)->pdgId() == 22 ) {
	     for( unsigned igmoth = 0 ; igmoth < imc->mother(imoth)->numberOfMothers(); igmoth++ )
	       if( imc->mother(imoth)->mother(igmoth)->pdgId() == 25 ) { hasHiggsAncestor = true; break; }		
	     }
	   if( hasHiggsAncestor ) break;
	 }
       }
     }
    
	   
    ///--- dR matching
    if( drMin > 0.1 ) return result;

    ///--- then compute generated isolation    
    double eIso = 0;
    for( vector<pat::PackedGenParticle>::const_iterator imc = genParticlesHandle2->begin(); imc != genParticlesHandle2->end(); ++imc ) {
      if( imc->status() == 1 ) {
	double dR = deltaR( matchedPhotonEta, matchedPhotonPhi, imc->eta(), imc->phi() );
	if( dR < 0.3 ) eIso += imc->energy();	
      }
    }

    eIso -=  matchedPhotonE;
    double eRelIso = eIso / matchedPhotonE;
    result.eMC  = matchedPhotonE;
    result.eIso    = eIso    > 0 ? eIso    : 0; 
    result.eRelIso = eRelIso > 0 ? eRelIso : 0; 
    result.mcMatching = 0;
    if( eRelIso < 0.15   ) result.mcMatching = 1;
    if( hasHiggsAncestor ) result.mcMatching = 2;

    //  cout << "      ===> matched photon PT: " << matchedPhotonE / cosh(matchedPhotonEta) << "  iso: " << eIso << " relIso: " << eRelIso << " : matching: " << result.mcMatching << endl;    
  }
  
  return result;
}






