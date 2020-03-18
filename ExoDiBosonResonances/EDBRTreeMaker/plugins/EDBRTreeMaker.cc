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


#include "userFun/EDBRTreeMaker.h"
#include "userFun/findLastCopy.h"
#include "userFun/EDBR_METInfo.h"
#include "userFun/EDBR_JetInfo.h"
#include "userFun/EDBR_setDummyValues.h"
#include "userFun/EDBR_HLT.h"
#include "userFun/EDBR_GenInfo.h"
#include "userFun/EDBR_PuppiAK8.h"







float
EDBRTreeMaker::dEtaInSeed( const pat::Electron*  ele ){
    return ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull() ? ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();

}







// ------------ method called for each event  ------------
void
EDBRTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
    
    setDummyValues(); //Initalize variables with dummy values

    nevent = iEvent.eventAuxiliary().event();
    run    = iEvent.eventAuxiliary().run();
    ls     = iEvent.eventAuxiliary().luminosityBlock();

    //std::cout<< "num of run:" << run << "  lumi:" << ls << "  n event:" << nevent << std::endl;
    Handle<TriggerResults> trigRes;
    iEvent.getByToken(hltToken_, trigRes);
    HLTStore(trigRes);


    edm::Handle<edm::View<reco::Candidate> > leptonicVs;
    iEvent.getByToken(leptonicVSrc_, leptonicVs);


    edm::Handle<double> rho;
    iEvent.getByToken(rhoToken_      , rho     );
    double fastJetRho = *(rho.product());
    useless = fastJetRho;


    edm::Handle<edm::View<pat::Jet> > ak4jets;
    iEvent.getByToken(ak4jetsSrc_, ak4jets);
    


    edm::Handle<edm::View<reco::Candidate> > metHandle;
    iEvent.getByToken(metSrc_, metHandle);

  
    edm::Handle<edm::View<pat::Muon>> loosemus;
    iEvent.getByToken(loosemuonToken_,loosemus);
    edm::Handle<edm::View<pat::Muon>> goodmus;
    iEvent.getByToken(goodMuSrc_, goodmus);
    edm::Handle<edm::View<pat::Muon>> mus;
    iEvent.getByToken(MuSrc_, mus);
    edm::Handle<edm::View<pat::Electron>> eles;
    iEvent.getByToken(EleSrc_, eles);



    edm::Handle<edm::View<pat::Electron>> looseels;
    iEvent.getByToken(looseelectronToken_, looseels);

    edm::Handle<edm::View<reco::GenParticle> > genParticles;//define genParticle
    iEvent.getByToken(genSrc_, genParticles);
    GENStore(genParticles,RunOnMC_,mus,eles);


    if (RunOnSig_||RunOnMC_){
        //  L1 prefiring
        edm::Handle< double > theprefweight;
        iEvent.getByToken(prefweight_token, theprefweight ) ;
        L1prefiring =(*theprefweight);
        edm::Handle< double > theprefweightup;
        iEvent.getByToken(prefweightup_token, theprefweightup ) ;
        L1prefiringup =(*theprefweightup);        
        edm::Handle< double > theprefweightdown;
        iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
        L1prefiringdown =(*theprefweightdown);
        
 

        edm::Handle<GenEventInfoProduct> genEvtInfo;
        iEvent.getByToken(GenToken_,genEvtInfo);
        theWeight = genEvtInfo->weight();
        if(theWeight>0) nump = nump+1;
        if(theWeight<0) numm = numm+1;


        edm::Handle<std::vector<PileupSummaryInfo>>  PupInfo;
        iEvent.getByToken(PUToken_, PupInfo);
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            nBX = PVI->getBunchCrossing();
            if(nBX == 0) { // "0" is the in-time crossing, negative values are the early crossings, positive are late
                npT = PVI->getTrueNumInteractions();
                npIT = PVI->getPU_NumInteractions();
            }
        }
    }



    //filter
    iEvent.getByToken(noiseFilterToken_, noiseFilterBits_);
    const edm::TriggerNames &names = iEvent.triggerNames(*noiseFilterBits_);
    for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i) {
        if (names.triggerName(i) == HBHENoiseFilter_Selector_)
            passFilter_HBHE_ = noiseFilterBits_->accept(i); // TO BE USED
        if (names.triggerName(i) == HBHENoiseIsoFilter_Selector_)
            passFilter_HBHEIso_ = noiseFilterBits_->accept(i); // TO BE USED
        if (names.triggerName(i) == GlobalHaloNoiseFilter_Selector_)
            passFilter_GlobalHalo_ = noiseFilterBits_->accept(i); // TO BE USED
        if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
            passFilter_ECALDeadCell_ = noiseFilterBits_->accept(i); // under scrutiny
        if (names.triggerName(i) == GoodVtxNoiseFilter_Selector_)
            passFilter_GoodVtx_ = noiseFilterBits_->accept(i); // TO BE USED
        if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
            passFilter_EEBadSc_ = noiseFilterBits_->accept(i); // under scrutiny
    }

    edm::Handle<bool> badMuonResultHandle;
    edm::Handle<bool> badChargedHadronResultHandle;
    iEvent.getByToken(badMuon_Selector_, badMuonResultHandle);
    iEvent.getByToken(badChargedHadron_Selector_, badChargedHadronResultHandle);
    passFilter_badMuon_ = *badMuonResultHandle;
    passFilter_badChargedHadron_ = *badChargedHadronResultHandle;
    edm::Handle< bool > passecalBadCalibFilterUpdate ;
    iEvent.getByToken(ecalBadCalibFilterUpdate_token,passecalBadCalibFilterUpdate);
    passecalBadCalibFilterUpdate_ =  (*passecalBadCalibFilterUpdate );




        bool doPuppi  = iEvent.getByToken(puppijetInputToken_, puppijets_ );

    if((puppijets_->size()!= 0 )  && (leptonicVs->size()!= 0) ){

       const reco::Candidate& leptonicV = leptonicVs->at(0);


       const reco::Candidate& metCand = metHandle->at(0);
       const reco::Candidate& lepton = (*leptonicV.daughter(0));
       nLooseMu = loosemus->size();
       nLooseEle = looseels->size();

       edm::Handle<reco::VertexCollection> vertices;
       iEvent.getByToken(vtxToken_, vertices);
        // edm::Handle<reco::VertexCollection> vertices;
        // iEvent.getByLabel("offlineSlimmedPrimaryVertices", vertices);
        // iEvent.getByToken(vtxToken_, vertices);
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

        MET_Store(iEvent);        

        /// For the time being, set these to 1


        ptlep1       = leptonicV.daughter(0)->pt();
        ptlep2       = leptonicV.daughter(1)->pt();
        etalep1      = leptonicV.daughter(0)->eta();
        etalep2      = leptonicV.daughter(1)->eta();
        philep1      = leptonicV.daughter(0)->phi();
        philep2      = leptonicV.daughter(1)->phi();
        lep          = std::max(abs(leptonicV.daughter(0)->pdgId()), abs(leptonicV.daughter(1)->pdgId()));
        double energylep1     = leptonicV.daughter(0)->energy();
        met          = metCand.pt();
        metPhi         = metCand.phi();
        //cout<<met<<" met candidate "<<metPhi<<endl;

        ptVlep       = leptonicV.pt();
        yVlep        = leptonicV.eta();
        phiVlep      = leptonicV.phi();
        massVlep     = leptonicV.mass();
        mtVlep       = leptonicV.mt();


        ////////////////////////lep ID  ////////////////////////////////////
        if( leptonicV.daughter(0)->isMuon()||leptonicV.daughter(1)->isMuon()){

            const pat::Muon *mu1 = abs(leptonicV.daughter(0)->pdgId())==13 ?
                                                  (pat::Muon*)leptonicV.daughter(0):
                                                  (pat::Muon*)leptonicV.daughter(1);
            isHighPt = mu1->isHighPtMuon(vertices->at(0));
            trackIso = mu1->trackIso();
            muchaiso=mu1->pfIsolationR04().sumChargedHadronPt;
            muneuiso=mu1->pfIsolationR04().sumNeutralHadronEt;
            muphoiso=mu1->pfIsolationR04().sumPhotonEt;
            muPU=mu1->pfIsolationR04().sumPUPt;
            muisolation = (muchaiso+ std::max(0.0,muneuiso+muphoiso-0.5*muPU))/mu1->pt();

        }
        if( leptonicV.daughter(0)->isElectron()||leptonicV.daughter(1)->isElectron() ) {
            const pat::Electron *el1 = leptonicV.daughter(0)->isElectron() ?
                                                  (pat::Electron*)leptonicV.daughter(0):
                                                  (pat::Electron*)leptonicV.daughter(1);
            double etaSC1         = el1->superCluster()->eta();
            double d01            = (-1)*el1->gsfTrack()->dxy(firstGoodVertex->position());
            isHEEP = false;
            et = el1->energy()!=0. ? el1->et()/el1->energy()*el1->caloEnergy() : 0.;
            if( et > 35. ) {
                if( fabs(etaSC1) < 1.4442 ){
                    iso = el1->dr03EcalRecHitSumEt() + el1->dr03HcalDepth1TowerSumEt();
                    isoCut = 2 + 0.03*et + 0.28*fastJetRho;
                    if( el1->ecalDriven() == 1 && dEtaInSeed( el1 ) < 0.004 && el1->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
                        el1->hadronicOverEm() < (1./el1->superCluster()->energy()+0.05) &&
                        (el1->full5x5_e2x5Max()/el1->full5x5_e5x5() > 0.94 || el1->full5x5_e1x5()/el1->full5x5_e5x5() > 0.83) &&
                        el1->dr03TkSumPt() < 5. && el1->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&//numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
                        iso < isoCut && fabs(d01) < 0.02 ) isHEEP = true;
                    }
                    if( fabs(etaSC1) > 1.566 && fabs(etaSC1) < 2.5 ){
                        iso = el1->dr03EcalRecHitSumEt() + el1->dr03HcalDepth1TowerSumEt();
                        if( et <= 50 )
                            isoCut = 2.5 + 0.28*fastJetRho;
                        else
                            isoCut = 2.5+0.03*(et-50.) + 0.28*fastJetRho;
                        if( el1->ecalDriven() == 1 && dEtaInSeed( el1 ) < 0.006 && el1->deltaPhiSuperClusterTrackAtVtx() < 0.06 &&
                            el1->hadronicOverEm() < (5./el1->superCluster()->energy()+0.05) && el1->full5x5_sigmaIetaIeta() < 0.03 &&
                            el1->dr03TkSumPt() < 5. && el1->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&//numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
                            iso < isoCut && fabs(d01) < 0.05 ) isHEEP = true;
                    }
            }
        }

        /////////////////////////Leptonic Part///////////////////////////////
        TLorentzVector  glepton, gleptonicV,gleptonicV_new,gleptonicV_JEC_up,gleptonicV_JEC_down,gleptonicV_JER_up,gleptonicV_JER_down;
        glepton.SetPtEtaPhiE(ptlep1, etalep1, philep1, energylep1);
        math::XYZTLorentzVector neutrinoP4 = getNeutrinoP4(MET_et, MET_phi, glepton, 1);
        reco::CandidateBaseRef METBaseRef = metHandle->refAt(0);
        reco::ShallowCloneCandidate neutrino(METBaseRef, 0 , neutrinoP4);
        reco::CompositeCandidate WLeptonic;
        WLeptonic.addDaughter(lepton);
        WLeptonic.addDaughter(neutrino);
        AddFourMomenta addP4;
        addP4.set(WLeptonic);
        gleptonicV.SetPtEtaPhiM(WLeptonic.pt(),WLeptonic.eta(),WLeptonic.phi(),WLeptonic.mass());
        ptVlepJEC       = WLeptonic.pt();
        yVlepJEC        = WLeptonic.eta();
        phiVlepJEC      = WLeptonic.phi();
        massVlepJEC     = WLeptonic.mass();
        if (RunOnMC_){ 
        math::XYZTLorentzVector     neutrinoP4_new = getNeutrinoP4(MET_et_new, MET_phi_new, glepton, 1);
        math::XYZTLorentzVector     neutrinoP4_JEC_up = getNeutrinoP4(MET_et_JEC_up, MET_phi_JEC_up, glepton, 1);
        math::XYZTLorentzVector     neutrinoP4_JEC_down = getNeutrinoP4(MET_et_JEC_down, MET_phi_JEC_down, glepton, 1);
        math::XYZTLorentzVector     neutrinoP4_JER_up = getNeutrinoP4(MET_et_JER_up, MET_phi_JER_up, glepton, 1);
        math::XYZTLorentzVector     neutrinoP4_JER_down = getNeutrinoP4(MET_et_JER_down, MET_phi_JER_down, glepton, 1);
        reco::ShallowCloneCandidate neutrino_new(METBaseRef, 0, neutrinoP4_new);
        reco::ShallowCloneCandidate neutrino_JEC_up(METBaseRef, 0, neutrinoP4_JEC_up);
        reco::ShallowCloneCandidate neutrino_JEC_down(METBaseRef, 0, neutrinoP4_JEC_down);
        reco::ShallowCloneCandidate neutrino_JER_up(METBaseRef, 0, neutrinoP4_JER_up);
        reco::ShallowCloneCandidate neutrino_JER_down(METBaseRef, 0, neutrinoP4_JER_down);
        reco::CompositeCandidate    WLeptonic_new;
        reco::CompositeCandidate    WLeptonic_JEC_up;
        reco::CompositeCandidate    WLeptonic_JEC_down;
        reco::CompositeCandidate    WLeptonic_JER_up;
        reco::CompositeCandidate    WLeptonic_JER_down;
        WLeptonic_new.addDaughter(lepton);
        WLeptonic_new.addDaughter(neutrino_new);
        WLeptonic_JEC_up.addDaughter(lepton);
        WLeptonic_JEC_up.addDaughter(neutrino_JEC_up);
        WLeptonic_JEC_down.addDaughter(lepton);
        WLeptonic_JEC_down.addDaughter(neutrino_JEC_down);
        WLeptonic_JER_up.addDaughter(lepton);
        WLeptonic_JER_up.addDaughter(neutrino_JER_up);
        WLeptonic_JER_down.addDaughter(lepton);
        WLeptonic_JER_down.addDaughter(neutrino_JER_down);
        AddFourMomenta addP4_new;
        addP4_new.set(WLeptonic_new);
        AddFourMomenta addP4_JEC_up;
        addP4_JEC_up.set(WLeptonic_JEC_up);
        AddFourMomenta addP4_JEC_down;
        addP4_JEC_down.set(WLeptonic_JEC_down);
        AddFourMomenta addP4_JER_up;
        addP4_JER_up.set(WLeptonic_JER_up);
        AddFourMomenta addP4_JER_down;
        addP4_JER_down.set(WLeptonic_JER_down);
        gleptonicV_new.SetPtEtaPhiM(WLeptonic_new.pt(),WLeptonic_new.eta(),WLeptonic_new.phi(),WLeptonic_new.mass());
        gleptonicV_JEC_up.SetPtEtaPhiM(WLeptonic_JEC_up.pt(),WLeptonic_JEC_up.eta(),WLeptonic_JEC_up.phi(),WLeptonic_JEC_up.mass());
        gleptonicV_JEC_down.SetPtEtaPhiM(WLeptonic_JEC_down.pt(),WLeptonic_JEC_down.eta(),WLeptonic_JEC_down.phi(),WLeptonic_JEC_down.mass());
        gleptonicV_JER_down.SetPtEtaPhiM(WLeptonic_JER_down.pt(),WLeptonic_JER_down.eta(),WLeptonic_JER_down.phi(),WLeptonic_JER_down.mass());
        gleptonicV_JER_up.SetPtEtaPhiM(WLeptonic_JER_up.pt(),WLeptonic_JER_up.eta(),WLeptonic_JER_up.phi(),WLeptonic_JER_up.mass());
        
        ptVlepJEC_new    = WLeptonic_new.pt();
        yVlepJEC_new     = WLeptonic_new.eta();
        phiVlepJEC_new   = WLeptonic_new.phi();
        massVlepJEC_new  = WLeptonic_new.mass();
        mtVlepJEC_new    = WLeptonic_new.mt();
        //cout<<ptVlep<<" lep W "<<ptVlepJEC<<"   "<<yVlep<<" lep W "<<yVlepJEC<<"   "<<phiVlep<<" lep W "<<phiVlepJEC<<"   "<<massVlep<<" lep W "<<massVlepJEC<<"   "<<endl;
        //cout<<ptVlep<<" lep Wnew "<<ptVlepJEC_new<<"   "<<yVlep<<" lep W "<<yVlepJEC_new<<"   "<<phiVlep<<" lep W "<<phiVlepJEC_new<<"   "<<massVlep<<" lep W "<<massVlepJEC_new<<"   "<<endl;
        
        ptVlepJEC_JEC_up    = WLeptonic_JEC_up.pt();
        yVlepJEC_JEC_up     = WLeptonic_JEC_up.eta();
        phiVlepJEC_JEC_up   = WLeptonic_JEC_up.phi();
        massVlepJEC_JEC_up  = WLeptonic_JEC_up.mass();
        mtVlepJEC_JEC_up    = WLeptonic_JEC_up.mt();
        
        ptVlepJEC_JEC_down    = WLeptonic_JEC_down.pt();
        yVlepJEC_JEC_down     = WLeptonic_JEC_down.eta();
        phiVlepJEC_JEC_down   = WLeptonic_JEC_down.phi();
        massVlepJEC_JEC_down  = WLeptonic_JEC_down.mass();
        mtVlepJEC_JEC_down    = WLeptonic_JEC_down.mt();
        
        ptVlepJEC_JER_up    = WLeptonic_JER_up.pt();
        yVlepJEC_JER_up     = WLeptonic_JER_up.eta();
        phiVlepJEC_JER_up   = WLeptonic_JER_up.phi();
        massVlepJEC_JER_up  = WLeptonic_JER_up.mass();
        mtVlepJEC_JER_up    = WLeptonic_JER_up.mt();
        
        ptVlepJEC_JER_down    = WLeptonic_JER_down.pt();
        yVlepJEC_JER_down     = WLeptonic_JER_down.eta();
        phiVlepJEC_JER_down   = WLeptonic_JER_down.phi();
        massVlepJEC_JER_down  = WLeptonic_JER_down.mass();
        mtVlepJEC_JER_down    = WLeptonic_JER_down.mt();}
        //cout<<"mtVlepJEC"<<mtVlepJEC<<endl;
        ////////////////////////JEC for AK8/////////////////////////////////

        reco::Candidate::LorentzVector uncorrPrunedJet;



        if( doPuppi ){//1


            PuppiAK8_Store( fastJetRho );

            int nak4 = 0;
            double tj1=-10.0, tj2=-10.0;
 
            for (size_t ik=0; ik<ak4jets->size();ik++)
            {//3
                double corr = 1;
                reco::Candidate::LorentzVector uncorrJet;
                if( doCorrOnTheFly_ ){
                    uncorrJet = (*ak4jets)[ik].correctedP4(0);
                    jecAK4_->setJetEta( uncorrJet.eta() );
                    jecAK4_->setJetPt ( uncorrJet.pt() );
                    jecAK4_->setJetE ( uncorrJet.energy() );
                    jecAK4_->setRho ( fastJetRho );
                    jecAK4_->setNPV ( vertices->size() );
                    jecAK4_->setJetA ( (*ak4jets)[ik].jetArea() );
                    corr = jecAK4_->getCorrection();
                } else {uncorrJet = (*ak4jets)[ik].p4();}
    
                //if( (corr*uncorrJet.pt())>20 && (fabs((*ak4jets)[ik].eta()) < 5.0) && looseJetID((*ak4jets)[ik])>0 && dtemp>0.8 && nak4<8){
                if( (corr*uncorrJet.pt())>20 && (fabs((*ak4jets)[ik].eta()) < 5.0) && tightJetID((*ak4jets)[ik])>0 && nak4<8){
                    ak4jet_hf[nak4]=(*ak4jets)[ik].hadronFlavour();
                    ak4jet_pf[nak4]=(*ak4jets)[ik].partonFlavour();
                    ak4jet_pt[nak4] =  corr*uncorrJet.pt();
                    ak4jet_pt_uncorr[nak4] =  uncorrJet.pt();
                    ak4jet_eta[nak4] = (*ak4jets)[ik].eta();
                    ak4jet_phi[nak4] = (*ak4jets)[ik].phi();
                    ak4jet_e[nak4] =   corr*uncorrJet.energy();
                    ak4jet_csv[nak4] = (*ak4jets)[ik].bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
                    ak4jet_icsv[nak4] = (*ak4jets)[ik].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
                        ak4jet_deepcsvudsg[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probudsg");
                        ak4jet_deepcsvb[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probb");
                        ak4jet_deepcsvc[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probc");
                        ak4jet_deepcsvbb[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probbb");
                        ak4jet_deepcsvcc[nak4] = (*ak4jets)[ik].bDiscriminator("pfDeepCSVJetTags:probcc");
                    ak4jet_IDLoose[nak4] = tightJetID((*ak4jets)[ik]);
                    ak4jet_IDTight[nak4] = tightJetID((*ak4jets)[ik]);
                    if(ak4jet_pt[nak4]>tj1 ) {
                        if(tj1>tj2) {tj2=tj1; nj2=nj1;}
                        tj1=ak4jet_pt[nak4]; nj1=nak4;
                    }
                    else if(ak4jet_pt[nak4]>tj2){
                        tj2=ak4jet_pt[nak4]; nj2=nak4;}
                    nak4 = nak4 + 1;
                }
            
            }//3
         

            if(nj1>-1 && nj2>-1 && ak4jet_pt[nj1]>30. && ak4jet_pt[nj2]>30.) {
                vbfeta=fabs(ak4jet_eta[nj1]-ak4jet_eta[nj2]);
                TLorentzVector vbfj1, vbfj2;
                vbfj1.SetPtEtaPhiE(ak4jet_pt[nj1], ak4jet_eta[nj1], ak4jet_phi[nj1], ak4jet_e[nj1]);
                vbfj2.SetPtEtaPhiE(ak4jet_pt[nj2], ak4jet_eta[nj2], ak4jet_phi[nj2], ak4jet_e[nj2]);
                vbfmjj=(vbfj1+vbfj2).Mag();
            }

            if(vbfeta>4.0 && vbfmjj>400) {vbftag=1;}
	    

            
            deltaRlepjet = deltaR(etalep1,philep1,jetAK8puppi_eta,jetAK8puppi_phi);
            deltaRlepjet_2 = deltaR(etalep1,philep1,jetAK8puppi_eta_2,jetAK8puppi_phi_2);
            TLorentzVector ghadronicVpuppi, gravitonpuppiJEC,ghadronicVpuppi_2, gravitonpuppiJEC_2;
            ghadronicVpuppi.SetPtEtaPhiM(jetAK8puppi_ptJEC, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2.SetPtEtaPhiM(jetAK8puppi_ptJEC_2, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC = gleptonicV + ghadronicVpuppi+ ghadronicVpuppi_2;
            candMasspuppiJEC     = gravitonpuppiJEC.Mag();
            m_jlv     = (gleptonicV + ghadronicVpuppi).Mag();

            TLorentzVector lvw[3];
            if (RunOnMC_){ 
            TLorentzVector ghadronicVpuppi_new, gravitonpuppiJEC_new,ghadronicVpuppi_2_new;
            ghadronicVpuppi_new.SetPtEtaPhiM(jetAK8puppi_ptJEC_new, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_new.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_new, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_new = gleptonicV_new + ghadronicVpuppi_new+ ghadronicVpuppi_2_new;
            candMasspuppiJEC_new     = gravitonpuppiJEC_new.Mag();
            m_jlv_new     = (gleptonicV_new + ghadronicVpuppi_new).Mag();
	    //cout<<WLeptonic.pt()<<" old "<<WLeptonic.eta()<<"   "<<WLeptonic.phi()<<"   "<<WLeptonic.mass()<<endl;
	    //cout<<WLeptonic_new.pt()<<" new "<<WLeptonic_new.eta()<<"   "<<WLeptonic_new.phi()<<"   "<<WLeptonic_new.mass()<<endl;
            //cout<<jetAK8puppi_sd<<"  jetAK8puppi_sd  "<<ghadronicVpuppi.Mag()<<"   "<<ghadronicVpuppi_new.Mag()<<"   "<<ghadronicVpuppi.E()<<"   "<<ghadronicVpuppi_new.E()<<endl;    
            TLorentzVector ghadronicVpuppi_JEC_up, gravitonpuppiJEC_JEC_up,ghadronicVpuppi_2_JEC_up;
            ghadronicVpuppi_JEC_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_JEC_up, jetAK8puppi_eta, jetAK8puppi_phi,jetAK8puppi_sd);
            ghadronicVpuppi_2_JEC_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JEC_up, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_JEC_up = gleptonicV_JEC_up + ghadronicVpuppi_JEC_up+ ghadronicVpuppi_2_JEC_up;
            candMasspuppiJEC_JEC_up     = gravitonpuppiJEC_JEC_up.Mag();
            m_jlv_JEC_up     = (gleptonicV_JEC_up + ghadronicVpuppi_JEC_up).Mag();
	    //cout<<ghadronicVpuppi_2_JEC_up.Pt()<<"  "<<candMasspuppiJEC_JEC_up<<endl;
            //cout<<jetAK8puppi_ptJEC_JEC_up<<"  "<<gleptonicV_JEC_up.Pt()<<"  "<<m_jlv_JEC_up<<endl;
            
            TLorentzVector ghadronicVpuppi_JEC_down, gravitonpuppiJEC_JEC_down,ghadronicVpuppi_2_JEC_down;
            ghadronicVpuppi_JEC_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_JEC_down, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_JEC_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JEC_down, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_JEC_down = gleptonicV_JEC_down + ghadronicVpuppi_JEC_down+ ghadronicVpuppi_2_JEC_down;
            candMasspuppiJEC_JEC_down     = gravitonpuppiJEC_JEC_down.Mag();
            m_jlv_JEC_down     = (gleptonicV_JEC_down + ghadronicVpuppi_JEC_down).Mag();
            
            TLorentzVector ghadronicVpuppi_JER_up, gravitonpuppiJEC_JER_up,ghadronicVpuppi_2_JER_up;
            ghadronicVpuppi_JER_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_JER_up, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_JER_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JER_up, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_JER_up = gleptonicV_JER_up + ghadronicVpuppi_JER_up+ ghadronicVpuppi_2_JER_up;
            candMasspuppiJEC_JER_up     = gravitonpuppiJEC_JER_up.Mag();
            m_jlv_JER_up     = (gleptonicV_JER_up + ghadronicVpuppi_JER_up).Mag();
            
            TLorentzVector ghadronicVpuppi_JER_down, gravitonpuppiJEC_JER_down,ghadronicVpuppi_2_JER_down;
            ghadronicVpuppi_JER_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_JER_down, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_JER_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JER_down, jetAK8puppi_eta_2, jetAK8puppi_phi_2,jetAK8puppi_sd_2);
            gravitonpuppiJEC_JER_down = gleptonicV_JER_down + ghadronicVpuppi_JER_down+ ghadronicVpuppi_2_JER_down;
            candMasspuppiJEC_JER_down     = gravitonpuppiJEC_JER_down.Mag();
            m_jlv_JER_down     = (gleptonicV_JER_down + ghadronicVpuppi_JER_down).Mag();

	    //swap var and var_new
            double jetAK8puppi_ptJEC_tmp = jetAK8puppi_ptJEC ; jetAK8puppi_ptJEC = jetAK8puppi_ptJEC_new ; jetAK8puppi_ptJEC_new = jetAK8puppi_ptJEC_tmp;
            double jetAK8puppi_ptJEC_2_tmp = jetAK8puppi_ptJEC_2 ; jetAK8puppi_ptJEC_2 = jetAK8puppi_ptJEC_2_new ; jetAK8puppi_ptJEC_2_new = jetAK8puppi_ptJEC_2_tmp;
            double jetAK8puppi_ptJEC_3_tmp = jetAK8puppi_ptJEC_3 ; jetAK8puppi_ptJEC_3 = jetAK8puppi_ptJEC_3_new ; jetAK8puppi_ptJEC_3_new = jetAK8puppi_ptJEC_3_tmp;
            double jetAK8puppi_e_tmp = jetAK8puppi_e ; jetAK8puppi_e = jetAK8puppi_e_new ; jetAK8puppi_e_new = jetAK8puppi_e_tmp;
            double jetAK8puppi_e_2_tmp = jetAK8puppi_e_2 ; jetAK8puppi_e_2 = jetAK8puppi_e_2_new ; jetAK8puppi_e_2_new = jetAK8puppi_e_2_tmp;
            double jetAK8puppi_e_3_tmp = jetAK8puppi_e_3 ; jetAK8puppi_e_3 = jetAK8puppi_e_3_new ; jetAK8puppi_e_3_new = jetAK8puppi_e_3_tmp;
            double ptVlepJEC_tmp = ptVlepJEC ; ptVlepJEC = ptVlepJEC_new ; ptVlepJEC_new = ptVlepJEC_tmp;
            double yVlepJEC_tmp = yVlepJEC ; yVlepJEC = yVlepJEC_new ; yVlepJEC_new = yVlepJEC_tmp;
            double phiVlepJEC_tmp = phiVlepJEC ; phiVlepJEC = phiVlepJEC_new ; phiVlepJEC_new = phiVlepJEC_tmp;
            double massVlepJEC_tmp = massVlepJEC ; massVlepJEC = massVlepJEC_new ; massVlepJEC_new = massVlepJEC_tmp;
            double mtVlepJEC_tmp = mtVlepJEC ; mtVlepJEC = mtVlepJEC_new ; mtVlepJEC_new = mtVlepJEC_tmp;
            double MET_et_tmp = MET_et ; MET_et = MET_et_new ; MET_et_new = MET_et_tmp;
            double MET_phi_tmp = MET_phi ; MET_phi = MET_phi_new ; MET_phi_new = MET_phi_tmp;
            double m_jlv_tmp = m_jlv ; m_jlv = m_jlv_new ; m_jlv_new = m_jlv_tmp;
            double candMasspuppiJEC_tmp = candMasspuppiJEC ; candMasspuppiJEC = candMasspuppiJEC_new ; candMasspuppiJEC_new = candMasspuppiJEC_tmp;

            lvw[0] = gleptonicV_new;
            lvw[1] = ghadronicVpuppi_new;
            lvw[2] = ghadronicVpuppi_2_new;
	    }
            delPhijetmet = deltaPhi(jetAK8puppi_phi, MET_phi);
            delPhijetlep = deltaPhi(jetAK8puppi_phi, phiVlepJEC);
            delPhijetmet_2 = deltaPhi(jetAK8puppi_phi_2, MET_phi);
            delPhijetlep_2 = deltaPhi(jetAK8puppi_phi_2, phiVlepJEC);
        
            delPhilepmet = deltaPhi(philep1, MET_phi);
            mtVlepJEC       =   sqrt(2*ptlep1*MET_et*(1.0-cos(philep1-MET_phi))); //WLeptonic.mt();

            if (!RunOnMC_){ 
            lvw[0] = gleptonicV;
            lvw[1] = ghadronicVpuppi;
            lvw[2] = ghadronicVpuppi_2;}
            Double_t Wpt[3];
            Wpt[0]=ptVlepJEC;
            Wpt[1]=jetAK8puppi_ptJEC;
            Wpt[2]=jetAK8puppi_ptJEC_2;
            Int_t *indexx=new Int_t[3];
            TMath::Sort(3,Wpt,indexx,1);
            //cout<<Wpt[indexx[0]]<<"   "<<Wpt[indexx[1]]<<"   "<<Wpt[indexx[2]]<<"   "<<endl;
            massww[0] = (lvw[indexx[0]]+lvw[indexx[1]]).Mag();
            massww[1] = (lvw[indexx[0]]+lvw[indexx[2]]).Mag();
            massww[2] = (lvw[indexx[1]]+lvw[indexx[2]]).Mag();

            masslvj1 = (lvw[0]+lvw[1]).Mag();
            masslvj2 = (lvw[0]+lvw[2]).Mag();
            massj1j2 = (lvw[1]+lvw[2]).Mag();
           

        }//1
        if((nLooseEle==1||nLooseMu==1)&&jetAK8puppi_ptJEC>200&&IDLoose==1&& fabs(jetAK8puppi_eta)<2.4 && IDLoose>0&&((jetAK8puppi_ptJEC_2>100)?(fabs(jetAK8puppi_eta_2)<2.4 && IDLoose_2>0):1)&&((jetAK8puppi_ptJEC_3>100)?(fabs(jetAK8puppi_eta_3)<2.4 && IDLoose_3>0):1)) outTree_->Fill();
    outTreew_->Fill();
    //outTree_->Fill();
	}
    
    else {
        outTreew_->Fill();
    //outTree_->Fill();
    }
EDBR_Debug();
//cout<< "test end3" <<endl;
}
//-------------------------------------------------------------------------------------------------------------------------------------//






//define this as a plug-in
DEFINE_FWK_MODULE(EDBRTreeMaker);
