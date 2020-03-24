#ifndef _VVVPuppi_
#define _VVVPuppi_

void VVVTreeMaker::Weight_Store( edm::Event const & iEvent ){
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
        for(std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            nBX = PVI->getBunchCrossing();
            if(nBX == 0) { // "0" is the in-time crossing, negative values are the early crossings, positive are late
                npT = PVI->getTrueNumInteractions();
                npIT = PVI->getPU_NumInteractions();
            }
        }


}


void VVVTreeMaker::Filter_Store( edm::Event const & iEvent ){        
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

}

#endif
