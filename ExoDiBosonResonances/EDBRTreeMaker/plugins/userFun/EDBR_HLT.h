#ifndef _EDBR_HLT_
#define _EDBR_HLT_

// HLT


void
EDBRTreeMaker::HLTStore(edm::Handle<edm::TriggerResults> trigRes){

    int xtemp1=0;
    for (size_t i=0; i<elPaths1.size();i++) {
        xtemp1 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths1[i]));
        if(HLT_Ele1<xtemp1) HLT_Ele1=xtemp1;
    }
    int xtemp2=0;
    for (size_t i=0; i<elPaths2.size();i++) {
        xtemp2 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths2[i]));
        if(HLT_Ele2<xtemp2) HLT_Ele2=xtemp2;
    }
    int xtemp3=0;
    for (size_t i=0; i<elPaths3.size();i++) {
        xtemp3 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths3[i]));
        if(HLT_Ele3<xtemp3) HLT_Ele3=xtemp3;
    }
    int xtemp4=0;
    for (size_t i=0; i<elPaths4.size();i++) {
        xtemp4 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths4[i]));
        if(HLT_Ele4<xtemp4) HLT_Ele4=xtemp4;
    }
    int xtemp5=0;
    for (size_t i=0; i<elPaths5.size();i++) {
        xtemp5 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths5[i]));
        if(HLT_Ele5<xtemp5) HLT_Ele5=xtemp5;
    }
    int xtemp6=0;
    for (size_t i=0; i<elPaths6.size();i++) {
        xtemp6 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths6[i]));
        if(HLT_Ele6<xtemp6) HLT_Ele6=xtemp6;
    }
    int xtemp7=0;
    for (size_t i=0; i<elPaths7.size();i++) {
        xtemp7 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths7[i]));
        if(HLT_Ele7<xtemp7) HLT_Ele7=xtemp7;
    }
    int xtemp8=0;
    for (size_t i=0; i<elPaths8.size();i++) {
        xtemp8 = (int)trigRes->accept(hltConfig.triggerIndex(elPaths8[i]));
        if(HLT_Ele8<xtemp8) HLT_Ele8=xtemp8;
    }

    int mtemp1=0;
    for (size_t i=0; i<muPaths1.size();i++) {
        mtemp1 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths1[i]));
        if(HLT_Mu1<mtemp1) HLT_Mu1=mtemp1;
    }
    int mtemp2=0;
    for (size_t i=0; i<muPaths2.size();i++) {
        mtemp2 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths2[i]));
        if(HLT_Mu2<mtemp2) HLT_Mu2=mtemp2;
    }
    int mtemp3=0;
    for (size_t i=0; i<muPaths3.size();i++) {
        mtemp3 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths3[i]));
        if(HLT_Mu3<mtemp3) HLT_Mu3=mtemp3;
    }
    int mtemp4=0;
    for (size_t i=0; i<muPaths4.size();i++) {
        mtemp4 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths4[i]));
        if(HLT_Mu4<mtemp4) HLT_Mu4=mtemp4;
    }
    int mtemp5=0;
    for (size_t i=0; i<muPaths5.size();i++) {
        mtemp5 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths5[i]));
        if(HLT_Mu5<mtemp5) HLT_Mu5=mtemp5;
    }
    int mtemp6=0;
    for (size_t i=0; i<muPaths6.size();i++) {
        mtemp6 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths6[i]));
        if(HLT_Mu6<mtemp6) HLT_Mu6=mtemp6;
    }
    int mtemp7=0;
    for (size_t i=0; i<muPaths7.size();i++) {
        mtemp7 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths7[i]));
        if(HLT_Mu7<mtemp7) HLT_Mu7=mtemp7;
    }
    int mtemp8=0;
    for (size_t i=0; i<muPaths8.size();i++) {
        mtemp8 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths8[i]));
        if(HLT_Mu8<mtemp8) HLT_Mu8=mtemp8;
    }
    int mtemp9=0;
    for (size_t i=0; i<muPaths9.size();i++) {
        mtemp9 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths9[i]));
        if(HLT_Mu9<mtemp9) HLT_Mu9=mtemp9;
    }
    int mtemp10=0;
    for (size_t i=0; i<muPaths10.size();i++) {
        mtemp10 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths10[i]));
        if(HLT_Mu10<mtemp10) HLT_Mu10=mtemp10;
    }
    int mtemp11=0;
    for (size_t i=0; i<muPaths11.size();i++) {
        mtemp11 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths11[i]));
        if(HLT_Mu11<mtemp11) HLT_Mu11=mtemp11;
    }
    int mtemp12=0;
    for (size_t i=0; i<muPaths12.size();i++) {
        mtemp12 = (int)trigRes->accept(hltConfig.triggerIndex(muPaths12[i]));
        if(HLT_Mu12<mtemp12) HLT_Mu12=mtemp12;
    }

}


// ------------ method called once each job just before starting event loop  ------------
void 
EDBRTreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void EDBRTreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

    elPaths1.clear();
    elPaths2.clear();
    elPaths3.clear();
    elPaths4.clear();
    elPaths5.clear();
    elPaths6.clear();
    elPaths7.clear();
    elPaths8.clear();
    muPaths1.clear();
    muPaths2.clear();
    muPaths3.clear();
    muPaths4.clear();
    muPaths5.clear();
    muPaths6.clear();
    muPaths7.clear();
    muPaths8.clear();
    muPaths9.clear();
    muPaths10.clear();
    muPaths11.clear();
    muPaths12.clear();

    std::cout<<"-----begin-----"<<std::endl;
    bool changed;
    if ( !hltConfig.init(iRun, iSetup, "HLT", changed) ) {
        edm::LogError("HltAnalysis") << "Initialization of HLTConfigProvider failed!!";
        return;
    }
    for (size_t i = 0; i < elPaths1_.size(); i++) {
        std::vector<std::string> foundPaths1 = hltConfig.matched( hltConfig.triggerNames(), elPaths1_[i] );
        while ( !foundPaths1.empty() ){
            elPaths1.push_back( foundPaths1.back() );
            foundPaths1.pop_back(); }
                                                }
    for (size_t i = 0; i < muPaths1_.size(); i++) {
        std::vector<std::string> foundPaths1 = hltConfig.matched( hltConfig.triggerNames(), muPaths1_[i] );
        while ( !foundPaths1.empty() ){
            muPaths1.push_back( foundPaths1.back() );
            foundPaths1.pop_back();
        }
    }
    std::cout<<"\n************** HLT-1 Information **************\n";
    for (size_t i=0; i < elPaths1.size(); i++) std::cout << "\n Electron paths-1:    " << i<<"  "<<elPaths1[i].c_str() <<"\t"<< std::endl;
    for (size_t i=0; i < muPaths1.size(); i++) std::cout << "\n Muon paths-1:   " << i<<"  "<<muPaths1[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < elPaths2_.size(); i++) {
        std::vector<std::string> foundPaths2 = hltConfig.matched( hltConfig.triggerNames(), elPaths2_[i] );
        while ( !foundPaths2.empty() ){
            elPaths2.push_back( foundPaths2.back() );
            foundPaths2.pop_back();
        }
    }
    for (size_t i = 0; i < muPaths2_.size(); i++) {
        std::vector<std::string> foundPaths2 = hltConfig.matched( hltConfig.triggerNames(), muPaths2_[i] );
        while ( !foundPaths2.empty() ){
            muPaths2.push_back( foundPaths2.back() );
            foundPaths2.pop_back();
        }
    }

    std::cout<<"\n************** HLT-2 Information **************\n";
    for (size_t i=0; i < elPaths2.size(); i++) std::cout << "\n Electron paths-2:    " << i<<"  "<<elPaths2[i].c_str() <<"\t"<< std::endl;
    for (size_t i=0; i < muPaths2.size(); i++) std::cout << "\n Muon paths-2:   " << i<<"  "<<muPaths2[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < elPaths3_.size(); i++) {
        std::vector<std::string> foundPaths3 = hltConfig.matched( hltConfig.triggerNames(), elPaths3_[i] );
        while ( !foundPaths3.empty() ){
            elPaths3.push_back( foundPaths3.back() );
            foundPaths3.pop_back();
        }
    }

    for (size_t i = 0; i < muPaths3_.size(); i++) {
        std::vector<std::string> foundPaths3 = hltConfig.matched( hltConfig.triggerNames(), muPaths3_[i] );
        while ( !foundPaths3.empty() ){
            muPaths3.push_back( foundPaths3.back() );
            foundPaths3.pop_back();
        }
    }

    std::cout<<"\n************** HLT-3 Information **************\n";
    for (size_t i=0; i < elPaths3.size(); i++) std::cout << "\n Electron paths-3:    " << i<<"  "<<elPaths3[i].c_str() <<"\t"<< std::endl;
    for (size_t i=0; i < muPaths3.size(); i++) std::cout << "\n Muon paths-3:   " << i<<"  "<<muPaths3[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < elPaths4_.size(); i++) {
        std::vector<std::string> foundPaths4 = hltConfig.matched( hltConfig.triggerNames(), elPaths4_[i] );
        while ( !foundPaths4.empty() ){
            elPaths4.push_back( foundPaths4.back() );
            foundPaths4.pop_back();
        }
    }

    for (size_t i = 0; i < muPaths4_.size(); i++) {
        std::vector<std::string> foundPaths4 = hltConfig.matched( hltConfig.triggerNames(), muPaths4_[i] );
        while ( !foundPaths4.empty() ){
            muPaths4.push_back( foundPaths4.back() );
            foundPaths4.pop_back();
        }
    }

    std::cout<<"\n************** HLT-4 Information **************\n";
    for (size_t i=0; i < elPaths4.size(); i++) std::cout << "\n Electron paths-4:    " << i<<"  "<<elPaths4[i].c_str() <<"\t"<< std::endl;
    for (size_t i=0; i < muPaths4.size(); i++) std::cout << "\n Muon paths-4:   " << i<<"  "<<muPaths4[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < elPaths5_.size(); i++) {
        std::vector<std::string> foundPaths5 = hltConfig.matched( hltConfig.triggerNames(), elPaths5_[i] );
        while ( !foundPaths5.empty() ){
            elPaths5.push_back( foundPaths5.back() );
            foundPaths5.pop_back();
        }
    }

    for (size_t i = 0; i < muPaths5_.size(); i++) {
        std::vector<std::string> foundPaths5 = hltConfig.matched( hltConfig.triggerNames(), muPaths5_[i] );
        while ( !foundPaths5.empty() ){
            muPaths5.push_back( foundPaths5.back() );
            foundPaths5.pop_back();
        }
    }
    std::cout<<"\n************** HLT-5 Information **************\n";
    for (size_t i=0; i < elPaths5.size(); i++) std::cout << "\n Electron paths-5:    " << i<<"  "<<elPaths5[i].c_str() <<"\t"<< std::endl;
    for (size_t i=0; i < muPaths5.size(); i++) std::cout << "\n Muon paths-5:   " << i<<"  "<<muPaths5[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < elPaths6_.size(); i++) {
        std::vector<std::string> foundPaths6 = hltConfig.matched( hltConfig.triggerNames(), elPaths6_[i] );
        while ( !foundPaths6.empty() ){
            elPaths6.push_back( foundPaths6.back() );
            foundPaths6.pop_back();
        }
    }

    for (size_t i = 0; i < muPaths6_.size(); i++) {
        std::vector<std::string> foundPaths6 = hltConfig.matched( hltConfig.triggerNames(), muPaths6_[i] );
        while ( !foundPaths6.empty() ){
            muPaths6.push_back( foundPaths6.back() );
            foundPaths6.pop_back();
        }
    }

    std::cout<<"\n************** HLT-6 Information **************\n";
    for (size_t i=0; i < elPaths6.size(); i++) std::cout << "\n Electron paths-6:    " << i<<"  "<<elPaths6[i].c_str() <<"\t"<< std::endl;
    for (size_t i=0; i < muPaths6.size(); i++) std::cout << "\n Muon paths-6:   " << i<<"  "<<muPaths6[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < elPaths7_.size(); i++) {
        std::vector<std::string> foundPaths7 = hltConfig.matched( hltConfig.triggerNames(), elPaths7_[i] );
        while ( !foundPaths7.empty() ){
            elPaths7.push_back( foundPaths7.back() );
            foundPaths7.pop_back();
        }
    }

    for (size_t i = 0; i < muPaths7_.size(); i++) {
        std::vector<std::string> foundPaths7 = hltConfig.matched( hltConfig.triggerNames(), muPaths7_[i] );
        while ( !foundPaths7.empty() ){
            muPaths7.push_back( foundPaths7.back() );
            foundPaths7.pop_back();
        }
    }

    std::cout<<"\n************** HLT-7 Information **************\n";
    for (size_t i=0; i < elPaths7.size(); i++) std::cout << "\n Electron paths-7:    " << i<<"  "<<elPaths7[i].c_str() <<"\t"<< std::endl;
    for (size_t i=0; i < muPaths7.size(); i++) std::cout << "\n Muon paths-7:   " << i<<"  "<<muPaths7[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < elPaths8_.size(); i++) {
        std::vector<std::string> foundPaths8 = hltConfig.matched( hltConfig.triggerNames(), elPaths8_[i] );
        while ( !foundPaths8.empty() ){
            elPaths8.push_back( foundPaths8.back() );
            foundPaths8.pop_back();
        }
    }

    for (size_t i = 0; i < muPaths8_.size(); i++) {
        std::vector<std::string> foundPaths8 = hltConfig.matched( hltConfig.triggerNames(), muPaths8_[i] );
        while ( !foundPaths8.empty() ){
            muPaths8.push_back( foundPaths8.back() );
            foundPaths8.pop_back();
        }
    }

    std::cout<<"\n************** HLT-8 Information **************\n";
    for (size_t i=0; i < elPaths8.size(); i++) std::cout << "\n Electron paths-8:    " << i<<"  "<<elPaths8[i].c_str() <<"\t"<< std::endl;
    for (size_t i=0; i < muPaths8.size(); i++) std::cout << "\n Muon paths-8:   " << i<<"  "<<muPaths8[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < muPaths9_.size(); i++) {
        std::vector<std::string> foundPaths9 = hltConfig.matched( hltConfig.triggerNames(), muPaths9_[i] );
        while ( !foundPaths9.empty() ){
            muPaths9.push_back( foundPaths9.back() );
            foundPaths9.pop_back();
        }
    }

    std::cout<<"\n************** HLT-9 Information **************\n";
    for (size_t i=0; i < muPaths9.size(); i++) std::cout << "\n Muon paths-9:   " << i<<"  "<<muPaths9[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < muPaths10_.size(); i++) {
        std::vector<std::string> foundPaths10 = hltConfig.matched( hltConfig.triggerNames(), muPaths10_[i] );
        while ( !foundPaths10.empty() ){
            muPaths10.push_back( foundPaths10.back() );
            foundPaths10.pop_back();
        }
    }

    std::cout<<"\n************** HLT-10 Information **************\n";
    for (size_t i=0; i < muPaths10.size(); i++) std::cout << "\n Muon paths-10:   " << i<<"  "<<muPaths10[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < muPaths11_.size(); i++) {
        std::vector<std::string> foundPaths11 = hltConfig.matched( hltConfig.triggerNames(), muPaths11_[i] );
        while ( !foundPaths11.empty() ){
            muPaths11.push_back( foundPaths11.back() );
            foundPaths11.pop_back();
        }
    }

    std::cout<<"\n************** HLT-11 Information **************\n";
    for (size_t i=0; i < muPaths11.size(); i++) std::cout << "\n Muon paths-11:   " << i<<"  "<<muPaths11[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

    for (size_t i = 0; i < muPaths12_.size(); i++) {
        std::vector<std::string> foundPaths12 = hltConfig.matched( hltConfig.triggerNames(), muPaths12_[i] );
        while ( !foundPaths12.empty() ){
            muPaths12.push_back( foundPaths12.back() );
            foundPaths12.pop_back();
        }
    }

    std::cout<<"\n************** HLT-12 Information **************\n";
    for (size_t i=0; i < muPaths12.size(); i++) std::cout << "\n Muon paths-12:   " << i<<"  "<<muPaths12[i].c_str() <<"\t"<< std::endl;
    std::cout<<"\n*********************************************\n\n";

}

#endif
