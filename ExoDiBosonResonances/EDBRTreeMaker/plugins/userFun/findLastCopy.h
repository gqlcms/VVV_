#ifndef _userFun_
#define _userFun_

//find last copy

const reco::Candidate*  EDBRTreeMaker::findLastW(const reco::Candidate *particle,int IDpdg){
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

const reco::Candidate*  EDBRTreeMaker::findLasttau(const reco::Candidate *particle,int IDpdg){
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


const reco::Candidate*  EDBRTreeMaker::findFirstW(const reco::Candidate *particle,int IDpdg){
    if (particle->mother(0)!=NULL){
        if(abs(particle->mother(0)->pdgId()) == IDpdg )
        return (findFirstW(particle->mother(0),IDpdg));
    }
    //cout<<"check 3 "<<pw->pdgId()<<"    "<<pw->status()<<"   "<<pw->daughter(0)->pdgId()<<"    "<<endl;
    return particle;
}



#endif
