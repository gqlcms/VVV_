#ifndef _LeptonicWinfo_
#define _LeptonicWinfo_

void VVVTreeMaker::LeptonicW_Store( edm::Event const & iEvent,   edm::Handle<edm::View<reco::Candidate> > leptonicVs,     edm::Handle<edm::View<reco::Candidate> > metHandle ){

        const reco::Candidate& leptonicV = leptonicVs->at(0);        
        const reco::Candidate& lepton = (*leptonicV.daughter(0));

        TLorentzVector  glepton;

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

}

#endif
