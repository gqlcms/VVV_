#ifndef _EDBRGENInfo_
#define _EDBRGENInfo_

void
VVVTreeMaker::GENStore( edm::Event const & iEvent ){

    edm::Handle<edm::View<reco::GenParticle> > genParticles;//define genParticle
    iEvent.getByToken(genSrc_, genParticles);

    edm::Handle<edm::View<pat::Muon>> mus;
    iEvent.getByToken(MuSrc_, mus);
    edm::Handle<edm::View<pat::Electron>> eles;
    iEvent.getByToken(EleSrc_, eles);


    // ************************* Gen Level Information******************//
    if(RunOnMC_)
    {//MC Info
        for(size_t ik=0; ik<genParticles->size();ik++)
        {// loop on gen
            const reco::Candidate* ptop0 = &(*genParticles)[ik];
            const reco::Candidate* ptop=findLasttau(ptop0,6);
                if(ptop0->pdgId()== 6 && gentop_pt==-99) {
                    gentop_pt = ptop->pt();
                    gentop_eta = ptop->eta();
                    gentop_phi = ptop->phi();
                    gentop_mass = ptop->mass();
                    for(int i=0;ptop->daughter(i)!=NULL;i++){
                        //if(abs(ptop->daughter(0)->pdgId())!=24&&abs(ptop->daughter(1)->pdgId())!=5) cout<<"no bW  "<<i<<"   "<<ptop->daughter(i)->pdgId()<<"   "<<ptop->daughter(i)->status()<<endl;
                        if(abs(ptop->daughter(i)->pdgId())==24){
                            gent_w_pt=ptop->daughter(i)->pt();
                            gent_w_eta=ptop->daughter(i)->eta();
                            gent_w_phi=ptop->daughter(i)->phi();
                            gent_w_mass=ptop->daughter(i)->mass();
                            const reco::Candidate* ptw0 = ptop->daughter(i);
                            const reco::Candidate* ptw= findLastW(ptw0,24);
                            if(ptw->daughter(0)!=NULL)
                            {
                                //if(abs(ptw->daughter(0)->pdgId())>=5) cout<<"no W-qq   "<<ptw->daughter(0)->pdgId()<<"   "<<ptw->daughter(1)->pdgId()<<endl;
                                if( abs(ptw->daughter(0)->pdgId())<=6 ){
                                    gent_w_tag=4;
                                    gent_w_q1_pt=ptw->daughter(0)->pt();
                                    gent_w_q1_eta=ptw->daughter(0)->eta();
                                    gent_w_q1_phi=ptw->daughter(0)->phi();
                                    gent_w_q1_e=ptw->daughter(0)->energy();
                                    gent_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    gent_w_q2_pt=ptw->daughter(1)->pt();
                                    gent_w_q2_eta=ptw->daughter(1)->eta();
                                    gent_w_q2_phi=ptw->daughter(1)->phi();
                                    gent_w_q2_e=ptw->daughter(1)->energy();
                                    gent_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12 ) gent_w_tag=1;
                                if( abs(ptw->daughter(0)->pdgId())==14 ||abs(ptw->daughter(0)->pdgId())==13 ) gent_w_tag=2;
                                if( abs(ptw->daughter(0)->pdgId())==16 ||abs(ptw->daughter(0)->pdgId())==15 ) gent_w_tag=3;
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12||abs(ptw->daughter(0)->pdgId())==13 ||abs(ptw->daughter(0)->pdgId())==14||abs(ptw->daughter(0)->pdgId())==15 ||abs(ptw->daughter(0)->pdgId())==16)
                                {
                                    gent_w_q1_pt=ptw->daughter(0)->pt();
                                    gent_w_q1_eta=ptw->daughter(0)->eta();
                                    gent_w_q1_phi=ptw->daughter(0)->phi();
                                    gent_w_q1_e=ptw->daughter(0)->energy();
                                    gent_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    gent_w_q2_pt=ptw->daughter(1)->pt();
                                    gent_w_q2_eta=ptw->daughter(1)->eta();
                                    gent_w_q2_phi=ptw->daughter(1)->phi();
                                    gent_w_q2_e=ptw->daughter(1)->energy();
                                    gent_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                        }
                        }
                        if(abs(ptop->daughter(i)->pdgId())==5){
                            gent_b_pt=ptop->daughter(i)->pt();
                            gent_b_eta=ptop->daughter(i)->eta();
                            gent_b_phi=ptop->daughter(i)->phi();
                            gent_b_mass=ptop->daughter(i)->mass();
                        }
                }
                }
                if(ptop0->pdgId()== -6 && genantitop_pt==-99) {
                    genantitop_pt = ptop->pt();
                    genantitop_eta = ptop->eta();
                    genantitop_phi = ptop->phi();
                    genantitop_mass = ptop->mass();
                    for(int i=0;ptop->daughter(i)!=NULL;i++){
                        //cout<<i<<"   "<<ptop->daughter(i)->pdgId()<<"   "<<ptop->daughter(i)->status()<<endl;
                        if(abs(ptop->daughter(i)->pdgId())==24){
                            genantit_w_pt=ptop->daughter(i)->pt();
                            genantit_w_eta=ptop->daughter(i)->eta();
                            genantit_w_phi=ptop->daughter(i)->phi();
                            genantit_w_mass=ptop->daughter(i)->mass();
                            const reco::Candidate* ptw0 = ptop->daughter(i);
                            const reco::Candidate* ptw= findLastW(ptw0,24);
                            if(ptw->daughter(0)!=NULL)
                            {
                                if( abs(ptw->daughter(0)->pdgId())<=6 ){
                                    genantit_w_tag=4;
                                    genantit_w_q1_pt=ptw->daughter(0)->pt();
                                    genantit_w_q1_eta=ptw->daughter(0)->eta();
                                    genantit_w_q1_phi=ptw->daughter(0)->phi();
                                    genantit_w_q1_e=ptw->daughter(0)->energy();
                                    genantit_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    genantit_w_q2_pt=ptw->daughter(1)->pt();
                                    genantit_w_q2_eta=ptw->daughter(1)->eta();
                                    genantit_w_q2_phi=ptw->daughter(1)->phi();
                                    genantit_w_q2_e=ptw->daughter(1)->energy();
                                    genantit_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12 ) genantit_w_tag=1;
                                if( abs(ptw->daughter(0)->pdgId())==14 ||abs(ptw->daughter(0)->pdgId())==13 ) genantit_w_tag=2;
                                if( abs(ptw->daughter(0)->pdgId())==16 ||abs(ptw->daughter(0)->pdgId())==15 ) genantit_w_tag=3;
                                if( abs(ptw->daughter(0)->pdgId())==11 ||abs(ptw->daughter(0)->pdgId())==12||abs(ptw->daughter(0)->pdgId())==13 ||abs(ptw->daughter(0)->pdgId())==14||abs(ptw->daughter(0)->pdgId())==15 ||abs(ptw->daughter(0)->pdgId())==16)
                                {
                                    genantit_w_q1_pt=ptw->daughter(0)->pt();
                                    genantit_w_q1_eta=ptw->daughter(0)->eta();
                                    genantit_w_q1_phi=ptw->daughter(0)->phi();
                                    genantit_w_q1_e=ptw->daughter(0)->energy();
                                    genantit_w_q1_pdg=ptw->daughter(0)->pdgId();
                                    genantit_w_q2_pt=ptw->daughter(1)->pt();
                                    genantit_w_q2_eta=ptw->daughter(1)->eta();
                                    genantit_w_q2_phi=ptw->daughter(1)->phi();
                                    genantit_w_q2_e=ptw->daughter(1)->energy();
                                    genantit_w_q2_pdg=ptw->daughter(1)->pdgId();
                                }
                            }
                        }
                        if(abs(ptop->daughter(i)->pdgId())==5){
                            genantit_b_pt=ptop->daughter(i)->pt();
                            genantit_b_eta=ptop->daughter(i)->eta();
                            genantit_b_phi=ptop->daughter(i)->phi();
                            genantit_b_mass=ptop->daughter(i)->mass();
                        }
                    }
                }


    }//end of loop on gen

    if(gen_mu_pt>0. && mus->size()>0 ){
        double drmumatch=10000.;  size_t mk=0;
        for(size_t ik=0; ik<mus->size();ik++)
        {
            double drtemp=deltaR(gen_mu_eta,gen_mu_phi,(*mus)[ik].eta(),(*mus)[ik].phi());
            if (drtemp<drmumatch) {drmumatch=drtemp; mk=ik;}
        }
        genmatch_mu_pt=(*mus)[mk].pt();
        genmatch_mu_eta=(*mus)[mk].eta();
        genmatch_mu_phi=(*mus)[mk].phi();
        genmatch_mu_e=(*mus)[mk].energy();
        genmatch_mu_dr=drmumatch;
    }
    if(gen_ele_pt>0. && eles->size()>0)
    {
        double drelematch=10000.;  size_t mk=0;
        for(size_t ik=0; ik<eles->size();ik++)
        {
            double drtemp=deltaR(gen_ele_eta,gen_ele_phi,(*eles)[ik].eta(),(*eles)[ik].phi());
            if (drtemp<drelematch) {drelematch=drtemp; mk=ik;}
        }
        genmatch_ele_pt=(*eles)[mk].pt();
        genmatch_ele_eta=(*eles)[mk].eta();
        genmatch_ele_phi=(*eles)[mk].phi();
        genmatch_ele_e=(*eles)[mk].energy();
        genmatch_ele_dr=drelematch;
        }
        
    //w and top info

        int igenw=0;
        int sizew=5;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            if(abs((*genParticles)[ik].pdgId())==24)
                {
                    const reco::Candidate* pwtmp1 = &(*genParticles)[ik];
                    const reco::Candidate* pwtmp=findLastW(pwtmp1,24);

                    int woverlap=0;
                    for (int ia=0;ia<igenw;ia++){
                        if(pwtmp->pt()==ptgenwl[ia]) woverlap=1;
                    }
                    if(pwtmp->pt()>50&&igenw<sizew&&woverlap==0){
                    ptgenwl[igenw] = pwtmp->pt();
                    etagenwl[igenw] = pwtmp->eta();
                    phigenwl[igenw] = pwtmp->phi();
                    massgenwl[igenw] = pwtmp->mass();
                    const reco::Candidate* pwtmp2=findFirstW(pwtmp1,24);
                    ptgenwf[igenw] = pwtmp2->pt();
                    etagenwf[igenw] = pwtmp2->eta();
                    phigenwf[igenw] = pwtmp2->phi();
                    massgenwf[igenw] = pwtmp2->mass();
                        taggenwmother[igenw]=pwtmp2->mother(0)->pdgId();
                        
                    if(pwtmp->daughter(0)!=NULL)//loop on w daughter
                    {
                        const reco::Candidate* pltmp = pwtmp->daughter(0);
                        //std::cout<< "pl pdgId" << pl->pdgId() << std::endl;
                         if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ){
                                taggenwl[igenw]=1;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ){
                            taggenwl[igenw]=2;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16) ){
                            taggenwl[igenw]=3;                    }//end of w daugter loop
                        if(abs(pltmp->pdgId())<6 ) {
                            taggenwl[igenw]=4;
                            genw_q1_pt[igenw]=pwtmp->daughter(0)->pt();
                            genw_q1_eta[igenw]=pwtmp->daughter(0)->eta();
                            genw_q1_phi[igenw]=pwtmp->daughter(0)->phi();
                            genw_q1_e[igenw]=pwtmp->daughter(0)->energy();
                            genw_q1_pdg[igenw]=pwtmp->daughter(0)->pdgId();
                            genw_q2_pt[igenw]=pwtmp->daughter(1)->pt();
                            genw_q2_eta[igenw]=pwtmp->daughter(1)->eta();
                            genw_q2_phi[igenw]=pwtmp->daughter(1)->phi();
                            genw_q2_e[igenw]=pwtmp->daughter(1)->energy();
                            genw_q2_pdg[igenw]=pwtmp->daughter(1)->pdgId();
                            //cout<<"w_q  "<<igenw<<"  "<<pwtmp->daughter(0)->pt()<<"   "<<pwtmp->daughter(1)->pt()<<endl;
                            }
                        if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ||(abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ||(abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16)){
                            genw_q1_pt[igenw]=pwtmp->daughter(0)->pt();
                            genw_q1_eta[igenw]=pwtmp->daughter(0)->eta();
                            genw_q1_phi[igenw]=pwtmp->daughter(0)->phi();
                            genw_q1_e[igenw]=pwtmp->daughter(0)->energy();
                            genw_q1_pdg[igenw]=pwtmp->daughter(0)->pdgId();
                            genw_q2_pt[igenw]=pwtmp->daughter(1)->pt();
                            genw_q2_eta[igenw]=pwtmp->daughter(1)->eta();
                            genw_q2_phi[igenw]=pwtmp->daughter(1)->phi();
                            genw_q2_e[igenw]=pwtmp->daughter(1)->energy();
                            genw_q2_pdg[igenw]=pwtmp->daughter(1)->pdgId();

                        }
                    }
                    cout<<"taggenwl[igenw]"<<taggenwl[igenw]<<endl;
                    igenw+=1;
                    }
                
                }//end of if w

        }

        int igenz=0;
        for(size_t ik=0; ik<genParticles->size();ik++)
        {
            if(abs((*genParticles)[ik].pdgId())==23)
            {
                const reco::Candidate* pztmp1 = &(*genParticles)[ik];
                const reco::Candidate* pztmp=findLasttau(pztmp1,23);

                int zoverlap=0;
                for (int ia=0;ia<igenz;ia++){
                    if(pztmp->pt()==ptgenzl[ia]) zoverlap=1;}
                if(pztmp->pt()>50&&igenz<sizew&&zoverlap==0){
                    ptgenzl[igenz] = pztmp->pt();
                    etagenzl[igenz] = pztmp->eta();
                    phigenzl[igenz] = pztmp->phi();
                    massgenzl[igenz] = pztmp->mass();
                    const reco::Candidate* pztmp2=findFirstW(pztmp1,23);
                    ptgenzf[igenz] = pztmp2->pt();
                    etagenzf[igenz] = pztmp2->eta();
                    phigenzf[igenz] = pztmp2->phi();
                    massgenzf[igenz] = pztmp2->mass();
                    //for(int i=0;pz->daughter(i)!=NULL;i++)//loop on w daughter
                    if(pztmp->daughter(0)!=NULL)//loop on w daughter
                    {
                        const reco::Candidate* pltmp = pztmp->daughter(0);
                        //std::cout<< "pl pdgId" << pl->pdgId() << std::endl;
                        if( (abs(pltmp->pdgId())==11) || (abs(pltmp->pdgId())==12) ){
                            taggenzl[igenz]=1;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==13) || (abs(pltmp->pdgId())==14) ){
                            taggenzl[igenz]=2;                    }//end of w daugter loop
                        if( (abs(pltmp->pdgId())==15) || (abs(pltmp->pdgId())==16) ){
                            taggenzl[igenz]=3;                    }//end of w daugter loop
                        if(abs(pltmp->pdgId())<6 ) {
                            taggenzl[igenz]=4;}
                    }
                    igenz+=1;
                }
            }//end of if w
        }



    }//end of MC Info
    // *************************End of Gen Level Information******************//


}

const reco::Candidate*  VVVTreeMaker::findLastW(const reco::Candidate *particle,int IDpdg){
    int iw=0;
    int pidw=0;
    const reco::Candidate* pw=particle;
    //cout<<"check 1 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<endl;
    for(int ii=0;particle->daughter(ii)!=NULL;ii++){
        if(abs(particle->daughter(ii)->pdgId())>pidw) {
            iw=ii;
            pidw=abs(particle->daughter(ii)->pdgId());
            //cout<<"check 2 "<<iw<<"    "<<pidw<<"   "<<endl;
        }
    }
    if( abs(pidw) == IDpdg ){
        pw = particle->daughter(iw);
        //cout<<"check 5 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<endl;
        return (findLastW(pw,IDpdg));
    }
    //cout<<"check 3 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<pw->daughter(0)->pdgId()<<"    "<<endl;
    return pw;
}

const reco::Candidate*  VVVTreeMaker::findLasttau(const reco::Candidate *particle,int IDpdg){
    int iw=0;
    int pidw=0;
    const reco::Candidate* pw=particle;
    //cout<<"check 1 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<endl;
    for(int ii=0;particle->daughter(ii)!=NULL;ii++){
        if(abs(particle->daughter(ii)->pdgId())== IDpdg) {
            iw=ii;
            pidw=abs(particle->daughter(ii)->pdgId());
            //cout<<"check 2 "<<iw<<"    "<<pidw<<"   "<<endl;
        }
    }
    if( abs(pidw) == IDpdg ){
        pw = particle->daughter(iw);
        //cout<<"check 5 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<endl;
        
        return (findLasttau(pw,IDpdg));
    }
    //cout<<"check 3 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<pw->daughter(0)->pdgId()<<"    "<<endl;
    return pw;
}


const reco::Candidate*  VVVTreeMaker::findFirstW(const reco::Candidate *particle,int IDpdg){
    if (particle->mother(0)!=NULL){
        if(abs(particle->mother(0)->pdgId()) == IDpdg )
        return (findFirstW(particle->mother(0),IDpdg));
    }
    //cout<<"check 3 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<pw->daughter(0)->pdgId()<<"    "<<endl;
    return particle;
}



#endif
