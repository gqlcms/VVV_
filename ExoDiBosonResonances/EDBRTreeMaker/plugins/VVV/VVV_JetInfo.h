#ifndef _EDBRJetID_
#define _EDBRJetID_

// ========================================================== JetID ===============================================================

bool
VVVTreeMaker::looseJetID( const pat::Jet& j ) {
    // refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    double NHF = j.neutralHadronEnergyFraction();
    double NEMF = j.neutralEmEnergyFraction();
    double CHF = j.chargedHadronEnergyFraction();
    //double MUF = j.muonEnergyFraction();
    double CEMF = j.chargedEmEnergyFraction();
    int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    int NumNeutralParticle =j.neutralMultiplicity();
    int CHM = j.chargedMultiplicity();
    double eta = j.eta();
	return ( (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ) || (NEMF<0.90 && NumNeutralParticle>10 && abs(eta)>3.0 )  ;
}

bool
VVVTreeMaker::tightJetID( const pat::Jet& j ) {
    // refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    if(j.pt()>0.){
    double NHF = j.neutralHadronEnergyFraction();
    double NEMF = j.neutralEmEnergyFraction();
    double CHF = j.chargedHadronEnergyFraction();
    //double MUF = j.muonEnergyFraction();
    //double CEMF = j.chargedEmEnergyFraction();
    int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    int NumNeutralParticle =j.neutralMultiplicity();
    int CHM = j.chargedMultiplicity();
    double eta = j.eta();
    return ((  (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 ) || abs(eta)>2.4) && abs(eta)<=2.7 ) || (NHF<0.99 && NEMF>0.02 && NumNeutralParticle>2 && abs(eta)>2.7 && abs(eta)<=3.0 ) || (NEMF<0.90 && NHF>0.02 &&NumNeutralParticle>10 && abs(eta)>3.0) ) ;
}
else{
return (0);
    }
}

bool
VVVTreeMaker::tightJetIDpuppi( const pat::Jet& j ) {
    // refer to https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
    if(j.pt()>0.){
    double NHF = j.neutralHadronEnergyFraction();
    double NEMF = j.neutralEmEnergyFraction();
    double CHF = j.chargedHadronEnergyFraction();
    //double MUF = j.muonEnergyFraction();
    int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    int NumNeutralParticle =j.neutralMultiplicity();
    int CHM = j.chargedMultiplicity();
    double eta = j.eta();
    return ((  (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 ) || (abs(eta)>2.4 && abs(eta)<=2.7) )) || (NHF<0.99 && abs(eta)>2.7 && abs(eta)<=3.0 ) || (NEMF<0.90 && NHF>0.02 &&NumNeutralParticle>2 && NumNeutralParticle<15 && abs(eta)>3.0) ) ;
}
else{
return (0);
    }
}

// ==================================================================== Jec ==========================================================================

void VVVTreeMaker::initJetCorrFactors( void ){

    std::vector<JetCorrectorParameters> vPar;
    vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8puppiPayloadNames_.begin(), payloadEnd = jecAK8puppiPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecAK8puppi_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

    

    vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4PayloadNames_.begin(), payloadEnd = jecAK4PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    // Make the FactorizedJetCorrector
    jecAK4_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

    vPar.clear();
    for ( std::vector<std::string>::const_iterator payloadBegin = offsetCorrLabel_.begin(), payloadEnd = offsetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }
    jecOffset_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

}


double VVVTreeMaker::getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){

    double jetCorrFactor = 1.;
    if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
        jecAK4_->setJetEta( rawJetP4.eta() );
        jecAK4_->setJetPt ( rawJetP4.pt() );
        jecAK4_->setJetE  ( rawJetP4.energy() );
        jecAK4_->setJetPhi( rawJetP4.phi()    );
        jecAK4_->setJetA  ( jet.jetArea() );
        jecAK4_->setRho   ( *(rho_.product()) );
        jecAK4_->setNPV   ( nVtx );
        jetCorrFactor = jecAK4_->getCorrection();
    }
    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;
    return jetCorrFactor;
}

double VVVTreeMaker::getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){

    double jetCorrFactor = 1.;
    if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){
        jecOffset_->setJetEta( rawJetP4.eta()     );
        jecOffset_->setJetPt ( rawJetP4.pt()      );
        jecOffset_->setJetE  ( rawJetP4.energy()  );
        jecOffset_->setJetPhi( rawJetP4.phi()     );
        jecOffset_->setJetA  ( jet.jetArea()      );
        jecOffset_->setRho   ( *(rho_.product())  );
        jecOffset_->setNPV   ( nVtx  );
        jetCorrFactor = jecOffset_->getCorrection();
    }

    reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
    corrJetP4 *= jetCorrFactor;

    return jetCorrFactor;
}


#endif
