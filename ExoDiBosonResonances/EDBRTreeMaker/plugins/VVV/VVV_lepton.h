#ifndef _EDBRLeptonInfo_
#define _EDBRLeptonInfo_

void
VVVTreeMaker::Lepton_num_Store( edm::Event const & iEvent ){

    edm::Handle<edm::View<pat::Muon>> loosemus;
    iEvent.getByToken(loosemuonToken_,loosemus);
    edm::Handle<edm::View<pat::Electron>> looseels;
    iEvent.getByToken(looseelectronToken_, looseels);
       nLooseMu = loosemus->size();
       nLooseEle = looseels->size();

}

void
VVVTreeMaker::Lepton_fromLeptonicW_Store( edm::Handle<edm::View<reco::Candidate> > leptonicVs , edm::Handle<reco::VertexCollection> vertices){

       const reco::Candidate& leptonicV = leptonicVs->at(0);


        ptlep1       = leptonicV.daughter(0)->pt();
        ptlep2       = leptonicV.daughter(1)->pt();
        etalep1      = leptonicV.daughter(0)->eta();
        etalep2      = leptonicV.daughter(1)->eta();
        philep1      = leptonicV.daughter(0)->phi();
        philep2      = leptonicV.daughter(1)->phi();
        lep          = std::max(abs(leptonicV.daughter(0)->pdgId()), abs(leptonicV.daughter(1)->pdgId()));
        energylep1     = leptonicV.daughter(0)->energy();

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


}

#endif
