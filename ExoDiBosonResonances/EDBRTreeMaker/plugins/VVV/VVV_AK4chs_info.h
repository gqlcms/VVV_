#ifndef _VVVAK4chsInfo_
#define _VVVAK4chsInfo_

void
VVVTreeMaker::VVV_AK4chs_Store( edm::Event const & iEvent ){

    edm::Handle<edm::View<pat::Jet> > ak4jets;
    iEvent.getByToken(ak4jetsSrc_, ak4jets);

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
                    jecAK4_->setNPV ( nVtx );
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


}



#endif
