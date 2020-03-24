// system include files
#include <iostream>
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/JetReco/interface/Jet.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "EDBRChannels.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include <TFormula.h>

#define Pi 3.141593
using namespace std;
//
// class declaration
//


#include "VVV/VVVTreeMaker.h"
#include "VVV/VVV_METInfo.h"
#include "VVV/VVV_JetInfo.h"
#include "VVV/VVV_setDummyValues.h"
#include "VVV/VVV_HLT.h"
#include "VVV/VVV_GenInfo.h"
#include "VVV/VVV_PuppiAK8.h"
#include "VVV/VVV_weight_filter.h"
#include "VVV/VVV_leptonic_W.h"
#include "VVV/VVV_lepton.h"
#include "VVV/VVV_AK4chs_info.h"
#include "VVV/VVV_inv_massInfo.h"







float
VVVTreeMaker::dEtaInSeed( const pat::Electron*  ele ){
    return ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull() ? ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();

}







// ------------ method called for each event  ------------
void
VVVTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
    
    setDummyValues(); //Initalize variables with dummy values

    nevent = iEvent.eventAuxiliary().event();
    run    = iEvent.eventAuxiliary().run();
    ls     = iEvent.eventAuxiliary().luminosityBlock();

    //std::cout<< "num of run:" << run << "  lumi:" << ls << "  n event:" << nevent << std::endl;

    HLTStore(iEvent);


    GENStore(iEvent);

    if (RunOnSig_||RunOnMC_){
    Weight_Store(iEvent);
    }

    Filter_Store(iEvent);


    edm::Handle<edm::View<reco::Candidate> > metHandle;
    iEvent.getByToken(metSrc_, metHandle);

    edm::Handle<edm::View<reco::Candidate> > leptonicVs;
    iEvent.getByToken(leptonicVSrc_, leptonicVs);

    bool doPuppi  = iEvent.getByToken(puppijetInputToken_, puppijets_ );

    if((puppijets_->size()!= 0 )  && (leptonicVs->size()!= 0) ){
              
       edm::Handle<reco::VertexCollection> vertices;
       iEvent.getByToken(vtxToken_, vertices);
        if (vertices->empty()) return; // skip the event if no PV found
        nVtx = vertices->size();
        reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
        for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
            // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
            // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
            if (  /// !vtx->isFake() &&
  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
                && fabs(vtx->position().Z())<=24.0) {
                    firstGoodVertex = vtx;
                    break;
                }
        }
        if ( firstGoodVertex==vertices->end() ) return; // skip event if there are no good PVs
        // ***************************************************************** //


        Lepton_num_Store( iEvent );

        Lepton_fromLeptonicW_Store( leptonicVs, vertices );

        MET_Store(iEvent, metHandle);        

        LeptonicW_Store(iEvent, leptonicVs, metHandle );


        if( doPuppi ){//1

            PuppiAK8_Store( iEvent );
     
            VVV_AK4chs_Store( iEvent );

            VVV_inv_mass_Store();

        }//1

        if((nLooseEle==1||nLooseMu==1)&&jetAK8puppi_ptJEC>200&&IDLoose==1&& fabs(jetAK8puppi_eta)<2.4 && IDLoose>0&&((jetAK8puppi_ptJEC_2>100)?(fabs(jetAK8puppi_eta_2)<2.4 && IDLoose_2>0):1)&&((jetAK8puppi_ptJEC_3>100)?(fabs(jetAK8puppi_eta_3)<2.4 && IDLoose_3>0):1)) outTree_->Fill();
    outTreew_->Fill();
	}
    
    else {
        outTreew_->Fill();
    }

if (CoutOrNot_){
EDBR_Debug();
}

}
//-------------------------------------------------------------------------------------------------------------------------------------//


//define this as a plug-in
DEFINE_FWK_MODULE(VVVTreeMaker);
