
// -*- C++ -*-
//
// Package:    MicroAODAnalyzers/flashggCommissioning
// Class:      flashggCommissioning
// 
/**\class flashggCommissioning flashggCommissioning.cc MicroAODAnalyzers/flashggCommissioning/plugins/flashggCommissioning.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Carlotta Favaro
//         Created:  Wed, 18 Jun 2014 13:44:16 GMT
//
//


#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "flashgg/MicroAODAlgos/interface/VertexSelectorBase.h"
#include "flashgg/MicroAODFormats/interface/Photon.h"
#include "flashgg/MicroAODFormats/interface/DiPhotonCandidate.h"

#include "flashgg/MicroAODAlgos/interface/PhotonIdUtils.h"

#include "TTree.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees 

// per-event tree:
struct eventInfo {
  
  float nphoton;

};

// per-photon tree:
struct photonInfo {

  float pt;
  float eta;
  float phi;
  float energy;
  float e1x5;
  float e2x5;
  float e3x3;
  float e5x5;
  float sigmaIetaIeta;
  float sigmaEtaEta;
  float maxEnergyXtal;
  

};

// **********************************************************************

using namespace std;
using namespace edm;

// **********************************************************************

class flashggCommissioning : public edm::EDAnalyzer {
   public:
      explicit flashggCommissioning(const edm::ParameterSet&);
      ~flashggCommissioning();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:

      edm::Service<TFileService> fs_;

     

      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void initEventStructure();
  
      edm::EDGetTokenT<edm::View<flashgg::Photon> >            photonToken_;
      edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
      edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexToken_; 
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> >       pfcandidateToken_;

      TTree* photonTree; 
      photonInfo phoInfo;
      // add all variables as private members
      
      flashgg::PhotonIdUtils phou;
      


};

// ******************************************************************************************
// ******************************************************************************************

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
flashggCommissioning::flashggCommissioning(const edm::ParameterSet& iConfig):
  photonToken_(consumes<View<flashgg::Photon> >(iConfig.getUntrackedParameter<InputTag> ("PhotonTag", InputTag("flashggPhotons")))),
  diphotonToken_(consumes<View<flashgg::DiPhotonCandidate> >(iConfig.getUntrackedParameter<InputTag> ("DiPhotonTag", InputTag("flashggDiPhotons")))),
  vertexToken_(consumes<View<reco::Vertex> >(iConfig.getUntrackedParameter<InputTag> ("VertexTag", InputTag("offlineSlimmedPrimaryVertices")))),
  pfcandidateToken_(consumes<View<pat::PackedCandidate> >(iConfig.getUntrackedParameter<InputTag> ("PFCandidatesTag", InputTag("packedPFCandidates"))))
{
 
}


flashggCommissioning::~flashggCommissioning()
{
 
  
}



void
flashggCommissioning::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // ********************************************************************************

  // access edm objects

  Handle<View<flashgg::Photon> > photons;
  iEvent.getByToken(photonToken_,photons);
  const PtrVector<flashgg::Photon>& photonPointers = photons->ptrVector();

  Handle<View<flashgg::DiPhotonCandidate> > diphotons;
  iEvent.getByToken(diphotonToken_,diphotons);
  const PtrVector<flashgg::DiPhotonCandidate>& diphotonPointers = diphotons->ptrVector();  
   
  /*
  Handle<View<reco::Vertex> > primaryVertices;
  iEvent.getByToken(vertexToken_,primaryVertices);
  const PtrVector<reco::Vertex>& pvPointers = primaryVertices->ptrVector();
  */

  Handle<View<pat::PackedCandidate> > pfcandidates;
  iEvent.getByToken(pfcandidateToken_,pfcandidates);
  //const PtrVector<pat::PackedCandidate>& pfCandPointers = pfcandidates->ptrVector();

  // cout << "size = " << pvPointers.size() << " " << pfCandPointers.size() << endl;
  
  // ********************************************************************************

  for( size_t ipho = 0; ipho < photonPointers.size(); ipho++ ) {

    initEventStructure();

    Ptr<flashgg::Photon> phoPtr = photonPointers[ipho];
    // cout << " photon pt = " << phoPtr->pt() << endl;

    phoInfo.pt =  phoPtr->pt();
    phoInfo.eta    = phoPtr->eta();
    phoInfo.phi    = phoPtr->phi();
    phoInfo.energy = phoPtr->energy();
    phoInfo.e1x5   = phoPtr->e1x5();
    phoInfo.e2x5   = phoPtr->e2x5();
    phoInfo.e3x3   = phoPtr->e3x3();
    phoInfo.e5x5   = phoPtr->e5x5();
    phoInfo.sigmaIetaIeta = phoPtr->sigmaIetaIeta();
    phoInfo.sigmaEtaEta   = phoPtr->sigmaEtaEta();
    phoInfo.maxEnergyXtal = phoPtr->maxEnergyXtal();


    // cout << " isolation = " << phou.pfIsoChgWrtVtx( phoPtr, pvPointers[0], pfCandPointers, 0.3, 0.01, 0.1, 0.01 ) << endl;
    //cout << " isolation = " << phou.pfIsoGamma( phoPtr, pfCandPointers, 0.2, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0) << endl;
    
    photonTree->Fill();

  
    

  }   // end photon loop

  // ********************************************************************************

  for( size_t idipho = 0; idipho < diphotonPointers.size(); idipho++ ) {

    //cout << "dummy loop to shut up warnings....leadingpho_pt =  " << diphotonPointers[idipho]->leadingPhoton()->pt() << endl;

  }  // end diphoton candidate loop
  

}


void 
flashggCommissioning::beginJob()
{
  photonTree = fs_->make<TTree>("photonTree","per-photon tree");
  photonTree->Branch("photonBranch",&phoInfo.pt,"phoPt/F:phoEta/F:phoPhi/F:phoE/F:phoE1x5/F:phoE2x5/F:phoE3x3/F:phoE5x5/F:phoSigmaIEtaIEta/F:phoSigmaEtaEta/F:phoEmax/F");
}

void 
flashggCommissioning::endJob() 
{
}

void
flashggCommissioning::initEventStructure() 
{
  
  phoInfo.pt  = -999;
  phoInfo.eta = -999;
  phoInfo.phi = -999;
  phoInfo.energy = -999;
  phoInfo.e1x5 = -999;
  phoInfo.e2x5 = -999;
  phoInfo.e3x3 = -999;
  phoInfo.e5x5 = -999;
  phoInfo.sigmaIetaIeta = -999;
  phoInfo.sigmaEtaEta = -999;
  phoInfo.maxEnergyXtal = -999;

}


/*
void 
flashggCommissioning::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

/*
void 
flashggCommissioning::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

/*
void 
flashggCommissioning::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

/*
void 
flashggCommissioning::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

void
flashggCommissioning::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(flashggCommissioning);
