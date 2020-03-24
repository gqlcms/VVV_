#ifndef _METinfo_
#define _METinfo_

// ==================================================== METInfo ===============================================

//-------------------------------------------------------------------------------------------------------------------------------------//
//
// member functions
//
void VVVTreeMaker::addTypeICorr( edm::Event const & event ){
    TypeICorrMap_.clear();
    event.getByToken(jetToken_      , jets_    );
    event.getByToken(rhoToken_      , rho_     );
    //edm::Handle<double> rho_;
    //event.getByLabel("fixedGridRhoFastjetAll",rho_);
    //edm::Handle<reco::VertexCollection> vertices_;
    //event.getByLabel("offlineSlimmedPrimaryVertices", vertices_);
    //event.getByToken(vtxToken_, vertices_);
    edm::Handle<reco::VertexCollection> vertices_;
    event.getByToken(vtxToken_, vertices_);



    bool skipEM_                    = true;
    double skipEMfractionThreshold_ = 0.9;
    bool skipMuons_                 = true;
    
    std::string skipMuonSelection_string = "isGlobalMuon | isStandAloneMuon";
    StringCutObjectSelector<reco::Candidate>* skipMuonSelection_ = new StringCutObjectSelector<reco::Candidate>(skipMuonSelection_string,true);

    double jetCorrEtaMax_           = 9.9;
    double type1JetPtThreshold_     = 15.0; //10.0;

    double corrEx    = 0;
    double corrEy    = 0;
    double corrSumEt = 0;

    for (const pat::Jet &jet : *jets_) {

        double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
        if ( skipEM_ && emEnergyFraction > skipEMfractionThreshold_ ) continue;

        reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0);
        double corr = getJEC(rawJetP4, jet, jetCorrEtaMax_, jetCorrLabel_);

        if ( skipMuons_ ) {
            const std::vector<reco::CandidatePtr> & cands = jet.daughterPtrVector();
            for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();cand != cands.end(); ++cand ) {
                const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
                const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
                if ( mu != 0 && (*skipMuonSelection_)(*mu) ) {
                    reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
                    rawJetP4 -= muonP4;
                }
            }
        }

        reco::Candidate::LorentzVector corrJetP4 = corr*rawJetP4;

        if ( corrJetP4.pt() > type1JetPtThreshold_ ) {
            reco::Candidate::LorentzVector tmpP4 = jet.correctedP4(0);
            corr = getJECOffset(tmpP4, jet, jetCorrEtaMax_, offsetCorrLabel_);
            reco::Candidate::LorentzVector rawJetP4offsetCorr = corr*rawJetP4;

            corrEx    -= (corrJetP4.px() - rawJetP4offsetCorr.px());
            corrEy    -= (corrJetP4.py() - rawJetP4offsetCorr.py());
            corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());
        }
    }
    TypeICorrMap_["corrEx"]    = corrEx;
    TypeICorrMap_["corrEy"]    = corrEy;
    TypeICorrMap_["corrSumEt"] = corrSumEt;
}
void VVVTreeMaker::addTypeICorr_user(edm::Event const& event) {
    TypeICorrMap_user_.clear();
    edm::Handle<pat::JetCollection> jets_;
    event.getByToken(t1jetSrc_userak4_, jets_);
    double corrEx_JEC         = 0;
    double corrEy_JEC         = 0;
    double corrSumEt_JEC      = 0;
    double corrEx_JEC_up      = 0;
    double corrEy_JEC_up      = 0;
    double corrSumEt_JEC_up   = 0;
    double corrEx_JEC_down    = 0;
    double corrEy_JEC_down    = 0;
    double corrSumEt_JEC_down = 0;
    
    double corrEx_JER         = 0;
    double corrEy_JER         = 0;
    double corrSumEt_JER      = 0;
    double corrEx_JER_up      = 0;
    double corrEy_JER_up      = 0;
    double corrSumEt_JER_up   = 0;
    double corrEx_JER_down    = 0;
    double corrEy_JER_down    = 0;
    double corrSumEt_JER_down = 0;
    for (const pat::Jet& jet : *jets_) {
        corrEx_JEC += jet.userFloat("corrEx_MET_JEC");
        corrEy_JEC += jet.userFloat("corrEy_MET_JEC");
        corrSumEt_JEC += jet.userFloat("corrSumEt_MET_JEC");
        corrEx_JEC_up += jet.userFloat("corrEx_MET_JEC_up");
        corrEy_JEC_up += jet.userFloat("corrEy_MET_JEC_up");
        corrSumEt_JEC_up += jet.userFloat("corrSumEt_MET_JEC_up");
        corrEx_JEC_down += jet.userFloat("corrEx_MET_JEC_down");
        corrEy_JEC_down += jet.userFloat("corrEy_MET_JEC_down");
        corrSumEt_JEC_down += jet.userFloat("corrSumEt_MET_JEC_down");
        corrEx_JER += jet.userFloat("corrEx_MET_JER");
        corrEy_JER += jet.userFloat("corrEy_MET_JER");
        corrSumEt_JER += jet.userFloat("corrSumEt_MET_JER");
        corrEx_JER_up += jet.userFloat("corrEx_MET_JER_up");
        corrEy_JER_up += jet.userFloat("corrEy_MET_JER_up");
        corrSumEt_JER_up += jet.userFloat("corrSumEt_MET_JER_up");
        corrEx_JER_down += jet.userFloat("corrEx_MET_JER_down");
        corrEy_JER_down += jet.userFloat("corrEy_MET_JER_down");
        corrSumEt_JER_down += jet.userFloat("corrSumEt_MET_JER_down");
    }
    TypeICorrMap_user_["corrEx_JEC"]         = corrEx_JEC;
    TypeICorrMap_user_["corrEy_JEC"]         = corrEy_JEC;
    TypeICorrMap_user_["corrSumEt_JEC"]      = corrSumEt_JEC;
    TypeICorrMap_user_["corrEx_JEC_up"]      = corrEx_JEC_up;
    TypeICorrMap_user_["corrEy_JEC_up"]      = corrEy_JEC_up;
    TypeICorrMap_user_["corrSumEt_JEC_up"]   = corrSumEt_JEC_up;
    TypeICorrMap_user_["corrEx_JEC_down"]    = corrEx_JEC_down;
    TypeICorrMap_user_["corrEy_JEC_down"]    = corrEy_JEC_down;
    TypeICorrMap_user_["corrSumEt_JEC_down"] = corrSumEt_JEC_down;
    
    TypeICorrMap_user_["corrEx_JER"]         = corrEx_JER;
    TypeICorrMap_user_["corrEy_JER"]         = corrEy_JER;
    TypeICorrMap_user_["corrSumEt_JER"]      = corrSumEt_JER;
    TypeICorrMap_user_["corrEx_JER_up"]      = corrEx_JER_up;
    TypeICorrMap_user_["corrEy_JER_up"]      = corrEy_JER_up;
    TypeICorrMap_user_["corrSumEt_JER_up"]   = corrSumEt_JER_up;
    TypeICorrMap_user_["corrEx_JER_down"]    = corrEx_JER_down;
    TypeICorrMap_user_["corrEy_JER_down"]    = corrEy_JER_down;
    TypeICorrMap_user_["corrSumEt_JER_down"] = corrSumEt_JER_down;
}


//-------------------------------------------------------------------------------------------------------------------------------------//
math::XYZTLorentzVector
VVVTreeMaker::getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType){
    double leppt = lep.Pt();
    double lepphi = lep.Phi();
    double lepeta = lep.Eta();
    double lepenergy = lep.Energy();
    
    double metpt = MetPt;
    double metphi = MetPhi;
    
    double  px = metpt*cos(metphi);
    double  py = metpt*sin(metphi);
    double  pz = 0;
    double  pxl= leppt*cos(lepphi);
    double  pyl= leppt*sin(lepphi);
    double  pzl= leppt*sinh(lepeta);
    double  El = lepenergy;
    double  a = pow(MW_,2) + pow(px+pxl,2) + pow(py+pyl,2) - px*px - py*py - El*El + pzl*pzl;
    double  b = 2.*pzl;
    double  A = b*b -4.*El*El;
    double  B = 2.*a*b;
    double  C = a*a-4.*(px*px+py*py)*El*El;
    
    ///////////////////////////pz for fnal
    double M_mu =  0;
    
    //if(lepType==1)M_mu=0.105658367;//mu
    //if(lepType==0)M_mu=0.00051099891;//electron
    
    int type=2; // use the small abs real root
    
    a = MW_*MW_ - M_mu*M_mu + 2.0*pxl*px + 2.0*pyl*py;
    A = 4.0*(El*El - pzl*pzl);
    B = -4.0*a*pzl;
    C = 4.0*El*El*(px*px + py*py) - a*a;
    
    double tmproot = B*B - 4.0*A*C;
    
    if (tmproot<0) {
        //std::cout << "Complex root detected, taking real part..." << std::endl;
        pz = - B/(2*A); // take real part of complex roots
    }
    else {
        double tmpsol1 = (-B + sqrt(tmproot))/(2.0*A);
        double tmpsol2 = (-B - sqrt(tmproot))/(2.0*A);
        //std::cout << " Neutrino Solutions: " << tmpsol1 << ", " << tmpsol2 << std::endl;
        
        if (type == 0 ) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else { pz = tmpsol1; }
            // if pz is > 300 pick the most central root
            if ( abs(pz) > 300. ) {
                if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                else { pz = tmpsol2; }
            }
        }
        if (type == 1 ) {
            // two real roots, pick the one closest to pz of muon
            if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
            else {pz = tmpsol1; }
        }
        if (type == 2 ) {
            // pick the most central root.
            if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
            else { pz = tmpsol2; }
        }
        /*if (type == 3 ) {
         // pick the largest value of the cosine
         TVector3 p3w, p3mu;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol1);
         p3mu.SetXYZ(pxl, pyl, pzl );
         
         double sinthcm1 = 2.*(p3mu.Perp(p3w))/MW_;
         p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol2);
         double sinthcm2 = 2.*(p3mu.Perp(p3w))/MW_;
         
         double costhcm1 = sqrt(1. - sinthcm1*sinthcm1);
         double costhcm2 = sqrt(1. - sinthcm2*sinthcm2);
         
         if ( costhcm1 > costhcm2 ) { pz = tmpsol1; otherSol_ = tmpsol2; }
         else { pz = tmpsol2;otherSol_ = tmpsol1; }
         
         }*///end of type3
        
    }//endl of if real root
    
    //dont correct pt neutrino
    math::XYZTLorentzVector outP4(px,py,pz,sqrt(px*px+py*py+pz*pz));
    return outP4;
    
}//end neutrinoP4

void VVVTreeMaker::MET_Store( edm::Event const & iEvent, edm::Handle<edm::View<reco::Candidate> > metHandle ){
        // ************************* MET ********************** //

       const reco::Candidate& metCand = metHandle->at(0);
        met          = metCand.pt();
        metPhi         = metCand.phi();



        iEvent.getByToken(metInputToken_ , METs_ );



        addTypeICorr(iEvent);
	if (RunOnMC_) addTypeICorr_user(iEvent);
        for (const pat::MET &met : *METs_) {            

            const float rawPt = met.uncorPt();
            const float rawPhi = met.uncorPhi();
            const float rawSumEt = met.uncorSumEt();
            TVector2 rawMET_;
            rawMET_.SetMagPhi (rawPt, rawPhi );
            Double_t rawPx = rawMET_.Px();
            Double_t rawPy = rawMET_.Py();
            Double_t rawEt = std::hypot(rawPx,rawPy);
            METraw_et = rawEt;
            METraw_phi = rawPhi;
            METraw_sumEt = rawSumEt;
            
            double pxcorr = rawPx+TypeICorrMap_["corrEx"];
            double pycorr = rawPy+TypeICorrMap_["corrEy"];
            double et     = std::hypot(pxcorr,pycorr);
            double sumEtcorr = rawSumEt+TypeICorrMap_["corrSumEt"];
            
            TLorentzVector corrmet;

            corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
            MET_et = et;
            MET_phi = corrmet.Phi();
            MET_sumEt = sumEtcorr;




            Double_t rawPtc       = met.corPt();
            Double_t rawPhic   = met.corPhi();
            Double_t rawSumEtc = met.corSumEt();
            TVector2 rawMET_c;
            rawMET_c.SetMagPhi (rawPtc, rawPhic );
            Double_t rawPxc = rawMET_c.Px();
            Double_t rawPyc = rawMET_c.Py();
            Double_t rawEtc = std::hypot(rawPxc,rawPyc);
            MET_et_m = rawEtc;
            MET_phi_m = rawPhic;
            MET_sumEt_m = rawSumEtc;


            if (RunOnMC_){ 
            double pxcorr_newo= rawPx+TypeICorrMap_user_["corrEx_JEC"];
            double pycorr_newo= rawPy+TypeICorrMap_user_["corrEy_JEC"];
            double et_newo     = std::hypot(pxcorr_newo,pycorr_newo);
	    MET_et_old=et_newo;

            // Marked for debug
            //------------------central value, correction from JetuserDataak4---------------------
            double pxcorr_new= rawPx+TypeICorrMap_user_["corrEx_JEC"]+TypeICorrMap_user_["corrEx_JER"];
            double pycorr_new= rawPy+TypeICorrMap_user_["corrEy_JEC"]+TypeICorrMap_user_["corrEy_JER"];
            double et_new     = std::hypot(pxcorr_new,pycorr_new);
            double sumEtcorr_new = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC"]+TypeICorrMap_user_["corrSumEt_JER"];
            //----for JEC uncertainty study
            double pxcorr_JEC_up = rawPx+TypeICorrMap_user_["corrEx_JEC_up"]+TypeICorrMap_user_["corrEx_JER"];
            double pycorr_JEC_up = rawPy+TypeICorrMap_user_["corrEy_JEC_up"]+TypeICorrMap_user_["corrEy_JER"];
            double et_JEC_up     = std::hypot(pxcorr_JEC_up, pycorr_JEC_up);
            double sumEtcorr_JEC_up = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC_up"]+TypeICorrMap_user_["corrSumEt_JER"];
            double pxcorr_JEC_down = rawPx+TypeICorrMap_user_["corrEx_JEC_down"]+TypeICorrMap_user_["corrEx_JER"];
            double pycorr_JEC_down = rawPy+TypeICorrMap_user_["corrEy_JEC_down"]+TypeICorrMap_user_["corrEy_JER"];
            double et_JEC_down     = std::hypot(pxcorr_JEC_down, pycorr_JEC_down);
            double sumEtcorr_JEC_down = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC_down"]+TypeICorrMap_user_["corrSumEt_JER"];
            //----for JER uncertainty study
            double pxcorr_JER_up = rawPx+TypeICorrMap_user_["corrEx_JEC"]+TypeICorrMap_user_["corrEx_JER_up"];
            double pycorr_JER_up = rawPy+TypeICorrMap_user_["corrEy_JEC"]+TypeICorrMap_user_["corrEy_JER_up"];
            double et_JER_up     = std::hypot(pxcorr_JER_up, pycorr_JER_up);
            double sumEtcorr_JER_up = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC"]+TypeICorrMap_user_["corrSumEt_JER_up"];
            double pxcorr_JER_down = rawPx+TypeICorrMap_user_["corrEx_JEC"]+TypeICorrMap_user_["corrEx_JER_down"];
            double pycorr_JER_down = rawPy+TypeICorrMap_user_["corrEy_JEC"]+TypeICorrMap_user_["corrEy_JER_down"];
            double et_JER_down     = std::hypot(pxcorr_JER_down,pycorr_JER_down);
            double sumEtcorr_JER_down = rawSumEt+TypeICorrMap_user_["corrSumEt_JEC"]+TypeICorrMap_user_["corrSumEt_JER_down"];
            //------------ 
            // Marked for debug
            MET_et_new= et_new;
            MET_et_JEC_up = et_JEC_up;
            MET_et_JEC_down = et_JEC_down;
            MET_et_JER_up = et_JER_up;
            MET_et_JER_down = et_JER_down;
            
            corrmet.SetPxPyPzE(pxcorr_new,pycorr_new,0.,et_new);
            MET_phi_new = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JEC_up,pycorr_JEC_up,0.,et_JEC_up);
            MET_phi_JEC_up = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JEC_down,pycorr_JEC_down,0.,et_JEC_down);
            MET_phi_JEC_down = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JER_up,pycorr_JER_up,0.,et_JER_up);
            MET_phi_JER_up = corrmet.Phi();
            corrmet.SetPxPyPzE(pxcorr_JER_down,pycorr_JER_down,0.,et_JER_down);
            MET_phi_JER_down = corrmet.Phi();
            
            MET_sumEt_new = sumEtcorr_new;
            MET_sumEt_JEC_up = sumEtcorr_JEC_up;
            MET_sumEt_JEC_down = sumEtcorr_JEC_down;
            MET_sumEt_JER_up = sumEtcorr_JER_up;
            MET_sumEt_JER_down = sumEtcorr_JER_down;
            }// Marked for debug
            
        }
        // ***************************************************************** //
}


#endif
