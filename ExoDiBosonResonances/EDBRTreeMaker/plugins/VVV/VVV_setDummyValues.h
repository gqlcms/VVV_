#ifndef _Setdummy_
#define _Setdummy_

// setDummyValues

void VVVTreeMaker::setDummyValues() {
    npT=-1.;
    npIT=-1.;
    nBX=-1;
    nLooseEle      =-99;
    nLooseMu       =-99;

    nVtx           = -99;


    ptVlep         = -99;

    L1prefiring = -99;
    L1prefiringup = -99;
    L1prefiringdown = -99;
    
    jetAK8puppi_ptJEC         = -99;
    jetAK8puppi_eta         = -99;
    jetAK8puppi_phi         = -99;
    jetAK8puppi_tau1         = -99;
    jetAK8puppi_tau2         = -99;
    jetAK8puppi_tau3         = -99;
    jetAK8puppi_tau21         = -99;
    jetAK8puppi_tau4         = -99;
    jetAK8puppi_tau42         = -99;
    // DeepAK8
    jetAK8puppi_dnnTop        = -99;
    jetAK8puppi_dnnW          = -99;
    jetAK8puppi_dnnH4q        = -99;
    jetAK8puppi_dnnTop_2        = -99;
    jetAK8puppi_dnnW_2          = -99;
    jetAK8puppi_dnnH4q_2        = -99;
    jetAK8puppi_dnnTop_3        = -99;
    jetAK8puppi_dnnW_3          = -99;
    jetAK8puppi_dnnH4q_3        = -99;
    jetAK8puppi_dnnZ        = -99;
    jetAK8puppi_dnnZbb        = -99;
    jetAK8puppi_dnnHbb        = -99;
    jetAK8puppi_dnnZ_2        = -99;
    jetAK8puppi_dnnZbb_2        = -99;
    jetAK8puppi_dnnHbb_2        = -99;
    jetAK8puppi_dnnZ_3        = -99;
    jetAK8puppi_dnnZbb_3        = -99;
    jetAK8puppi_dnnHbb_3        = -99;
    
    jetAK8puppi_dnnqcd        = -99;
    jetAK8puppi_dnntop        = -99;
    jetAK8puppi_dnnw          = -99;
    jetAK8puppi_dnnz          = -99;
    jetAK8puppi_dnnzbb        = -99;
    jetAK8puppi_dnnhbb        = -99;
    jetAK8puppi_dnnh4q        = -99;
    jetAK8puppi_dnnqcd_2        = -99;
    jetAK8puppi_dnntop_2        = -99;
    jetAK8puppi_dnnw_2          = -99;
    jetAK8puppi_dnnz_2          = -99;
    jetAK8puppi_dnnzbb_2        = -99;
    jetAK8puppi_dnnhbb_2        = -99;
    jetAK8puppi_dnnh4q_2        = -99;
    jetAK8puppi_dnnqcd_3        = -99;
    jetAK8puppi_dnntop_3        = -99;
    jetAK8puppi_dnnw_3          = -99;
    jetAK8puppi_dnnz_3          = -99;
    jetAK8puppi_dnnzbb_3        = -99;
    jetAK8puppi_dnnhbb_3        = -99;
    jetAK8puppi_dnnh4q_3        = -99;

    // Decorrelated DeepAK8
    jetAK8puppi_dnnDecorrTop        = -99;
    jetAK8puppi_dnnDecorrW          = -99;
    jetAK8puppi_dnnDecorrH4q        = -99;
    jetAK8puppi_dnnDecorrTop_2        = -99;
    jetAK8puppi_dnnDecorrW_2          = -99;
    jetAK8puppi_dnnDecorrH4q_2        = -99;
    jetAK8puppi_dnnDecorrTop_3        = -99;
    jetAK8puppi_dnnDecorrW_3          = -99;
    jetAK8puppi_dnnDecorrH4q_3        = -99;
    jetAK8puppi_dnnDecorrZ        = -99;
    jetAK8puppi_dnnDecorrZbb        = -99;
    jetAK8puppi_dnnDecorrHbb        = -99;
    jetAK8puppi_dnnDecorrZ_2        = -99;
    jetAK8puppi_dnnDecorrZbb_2        = -99;
    jetAK8puppi_dnnDecorrHbb_2        = -99;
    jetAK8puppi_dnnDecorrZ_3        = -99;
    jetAK8puppi_dnnDecorrZbb_3        = -99;
    jetAK8puppi_dnnDecorrHbb_3        = -99;
    jetAK8puppi_dnnDecorrbb        = -99;
    jetAK8puppi_dnnDecorrcc        = -99;
    jetAK8puppi_dnnDecorrbbnog        = -99;
    jetAK8puppi_dnnDecorrccnog        = -99;
    jetAK8puppi_dnnDecorrbb_2        = -99;
    jetAK8puppi_dnnDecorrcc_2        = -99;
    jetAK8puppi_dnnDecorrbbnog_2        = -99;
    jetAK8puppi_dnnDecorrccnog_2        = -99;
    jetAK8puppi_dnnDecorrbb_3        = -99;
    jetAK8puppi_dnnDecorrcc_3        = -99;
    jetAK8puppi_dnnDecorrbbnog_3        = -99;
    jetAK8puppi_dnnDecorrccnog_3        = -99;
    
    jetAK8puppi_dnnDecorrqcd        = -99;
    jetAK8puppi_dnnDecorrtop        = -99;
    jetAK8puppi_dnnDecorrw          = -99;
    jetAK8puppi_dnnDecorrz          = -99;
    jetAK8puppi_dnnDecorrzbb        = -99;
    jetAK8puppi_dnnDecorrhbb        = -99;
    jetAK8puppi_dnnDecorrh4q        = -99;
    jetAK8puppi_dnnDecorrqcd_2        = -99;
    jetAK8puppi_dnnDecorrtop_2        = -99;
    jetAK8puppi_dnnDecorrw_2          = -99;
    jetAK8puppi_dnnDecorrz_2          = -99;
    jetAK8puppi_dnnDecorrzbb_2        = -99;
    jetAK8puppi_dnnDecorrhbb_2        = -99;
    jetAK8puppi_dnnDecorrh4q_2        = -99;
    jetAK8puppi_dnnDecorrqcd_3        = -99;
    jetAK8puppi_dnnDecorrtop_3        = -99;
    jetAK8puppi_dnnDecorrw_3          = -99;
    jetAK8puppi_dnnDecorrz_3          = -99;
    jetAK8puppi_dnnDecorrzbb_3        = -99;
    jetAK8puppi_dnnDecorrhbb_3        = -99;
    jetAK8puppi_dnnDecorrh4q_3        = -99;


    jetAK8puppi_sd         = -99;

    
    jetAK8puppi_ptJEC_2         = -99;
    jetAK8puppi_eta_2         = -99;
    jetAK8puppi_phi_2         = -99;
    jetAK8puppi_tau1_2         = -99;
    jetAK8puppi_tau2_2         = -99;
    jetAK8puppi_tau3_2         = -99;
    jetAK8puppi_tau21_2         = -99;
    jetAK8puppi_tau4_2         = -99;
    jetAK8puppi_tau42_2         = -99;
    jetAK8puppi_sd_2         = -99;

    
    jetAK8puppi_ptJEC_3         = -99;
    jetAK8puppi_eta_3         = -99;
    jetAK8puppi_phi_3         = -99;
    jetAK8puppi_tau1_3         = -99;
    jetAK8puppi_tau2_3         = -99;
    jetAK8puppi_tau3_3         = -99;
    jetAK8puppi_tau21_3         = -99;
    jetAK8puppi_tau4_3         = -99;
    jetAK8puppi_tau42_3         = -99;
    jetAK8puppi_sd_3         = -99;

    jetAK8puppi_ptJEC_new         = -99;
    jetAK8puppi_ptJEC_newnew         = -99;
    jetAK8puppi_ptJEC_m         = -99;
    jetAK8puppi_ptJEC_JEC_up         = -99;
    jetAK8puppi_ptJEC_JEC_down         = -99;
    jetAK8puppi_ptJEC_JER_up         = -99;
    jetAK8puppi_ptJEC_JER_down         = -99;
    jetAK8puppi_ptJEC_2_new         = -99;
    jetAK8puppi_ptJEC_2_JEC_up         = -99;
    jetAK8puppi_ptJEC_2_JEC_down         = -99;
    jetAK8puppi_ptJEC_2_JER_up         = -99;
    jetAK8puppi_ptJEC_2_JER_down         = -99;
    jetAK8puppi_ptJEC_3_new         = -99;
    jetAK8puppi_ptJEC_3_JEC_up         = -99;
    jetAK8puppi_ptJEC_3_JEC_down         = -99;
    jetAK8puppi_ptJEC_3_JER_up         = -99;
    jetAK8puppi_ptJEC_3_JER_down         = -99;
    
    jetAK8puppi_e         = -99;
    jetAK8puppi_e_new         = -99;
    jetAK8puppi_e_JEC_up         = -99;
    jetAK8puppi_e_JEC_down         = -99;
    jetAK8puppi_e_JER_up         = -99;
    jetAK8puppi_e_JER_down         = -99;
    jetAK8puppi_e_2         = -99;
    jetAK8puppi_e_2_new         = -99;
    jetAK8puppi_e_2_JEC_up         = -99;
    jetAK8puppi_e_2_JEC_down         = -99;
    jetAK8puppi_e_2_JER_up         = -99;
    jetAK8puppi_e_2_JER_down         = -99;
    jetAK8puppi_e_3         = -99;
    jetAK8puppi_e_3_new         = -99;
    jetAK8puppi_e_3_JEC_up         = -99;
    jetAK8puppi_e_3_JEC_down         = -99;
    jetAK8puppi_e_3_JER_up         = -99;
    jetAK8puppi_e_3_JER_down         = -99;
    

    vbfeta=-10.;
    vbfmjj=-10.;
    vbftag=0;
    nj1=-1;
    nj2=-1;

    yVlep          = -99;
    phiVlep        = -99;
    massVlep       = -99;
    mtVlep         = -99;
    ptlep1         = -99;
    ptlep2         = -99;
    etalep1        = -99;
    etalep2        = -99;
    philep1        = -99;
    philep2        = -99;
    met            = -99;
    metPhi         = -99;
    deltaRlepjet   = -99;
    delPhilepmet   = -99;
    delPhijetmet =  -99;
    delPhijetlep =  -99;
    
    deltaRlepjet_2   = -99;
    
    delPhijetmet_2 =  -99;
    delPhijetlep_2 =  -99;
    
    lep            = -99;

    gen_ele_pt     = -99;
    gen_ele_eta    = -99;
    gen_ele_phi    = -99;
    gen_ele_e      = -99;
    gen_mu_pt     = -99;
    gen_mu_eta    = -99;
    gen_mu_phi    = -99;
    gen_mu_e      = -99;
    genmatch_ele_pt     = -99;
    genmatch_ele_eta    = -99;
    genmatch_ele_phi    = -99;
    genmatch_ele_e      = -99;
    genmatch_ele_dr     =  99;
    genmatch_mu_pt     = -99;
    genmatch_mu_eta    = -99;
    genmatch_mu_phi    = -99;
    genmatch_mu_e      = -99;
    genmatch_mu_dr      = -99;
    gen_ele_pt_2     = -99;
    gen_ele_eta_2    = -99;
    gen_ele_phi_2    = -99;
    gen_ele_e_2      = -99;
    gen_mu_pt_2     = -99;
    gen_mu_eta_2    = -99;
    gen_mu_phi_2    = -99;
    gen_mu_e_2      = -99;
    gen_ele_pt_3     = -99;
    gen_ele_eta_3    = -99;
    gen_ele_phi_3    = -99;
    gen_ele_e_3      = -99;
    gen_mu_pt_3     = -99;
    gen_mu_eta_3    = -99;
    gen_mu_phi_3    = -99;
    gen_mu_e_3      = -99;
    




    

    gentop_pt  = -99;
    gentop_eta  = -99;
    gentop_phi  = -99;
    gentop_mass  = -99;
    genantitop_pt  = -99;
    genantitop_eta  = -99;
    genantitop_phi  = -99;
    genantitop_mass  = -99;



  


    /*for(int j=0; j<882; j++){
        pweight[j]=0.0;
    }*/


    for(Int_t ii=0;ii<8;ii++){
        ak4jet_hf[ii] = -99;
        ak4jet_pf[ii] = -99;
        ak4jet_pt[ii] = -99;
        ak4jet_pt_uncorr[ii] = -99;
        ak4jet_eta[ii] = -99;
        ak4jet_phi[ii] = -99;
        ak4jet_e[ii] = -99;
        ak4jet_dr[ii] = -99;
        ak4jet_csv[ii] = -99;
        ak4jet_icsv[ii] = -99;
        ak4jet_deepcsvudsg[ii] = -99;
        ak4jet_deepcsvb[ii] = -99;
        ak4jet_deepcsvc[ii] = -99;
        ak4jet_deepcsvbb[ii] = -99;
        ak4jet_deepcsvcc[ii] = -99;
        ak4jet_IDLoose[ii] = -99;
        ak4jet_IDTight[ii] = -99;
    }
    
    ak8sj11.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj12.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj13.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj14.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj15.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj21.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj22.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj23.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj24.SetPtEtaPhiM(0,-99,-99,-99);
    ak8sj25.SetPtEtaPhiM(0,-99,-99,-99);


    gleptonicV.SetPtEtaPhiM(0,-99,-99,-99);
    gleptonicV_new.SetPtEtaPhiM(0,-99,-99,-99);
    gleptonicV_JEC_up.SetPtEtaPhiM(0,-99,-99,-99);
    gleptonicV_JEC_down.SetPtEtaPhiM(0,-99,-99,-99);
    gleptonicV_JER_up.SetPtEtaPhiM(0,-99,-99,-99);
    gleptonicV_JER_down.SetPtEtaPhiM(0,-99,-99,-99);


    for(int i=0;i<4;i++){
        jetAK8puppi_mass1[i] = -99;
        jetAK8puppi_pt1[i] = -99;
        jetAK8puppi_eta1[i] = -99;
        jetAK8puppi_pt1_new[i] = -99;
        jetAK8puppi_pt1_m[i] = -99;
        jetAK8puppi_pt1_newnew[i] = -99;
        jetAK8puppi_pt1_JEC_up[i] = -99;
        jetAK8puppi_pt1_JEC_down[i] = -99;
        jetAK8puppi_pt1_JER_up[i] = -99;
        jetAK8puppi_pt1_JER_down[i] = -99;
        jetAK8puppi_e1_new[i] = -99;
        jetAK8puppi_e1_JEC_up[i] = -99;
        jetAK8puppi_e1_JEC_down[i] = -99;
        jetAK8puppi_e1_JER_up[i] = -99;
        jetAK8puppi_e1_JER_down[i] = -99;

        pttau[i] = -99 ;etatau[i] = -99 ;phitau[i] = -99 ;etau[i] = -99 ;pdgidtau[i] = -99 ;
        pttau_2[i] = -99 ;etatau_2[i] = -99 ;phitau_2[i] = -99 ;etau_2[i] = -99 ;pdgidtau_2[i] = -99 ;
        pttau_3[i] = -99 ;etatau_3[i] = -99 ;phitau_3[i] = -99 ;etau_3[i] = -99 ;pdgidtau_3[i] = -99 ;
    }
 
    for(int i=0;i<3;i++){
        ptq[i] = -99 ;etaq[i] = -99 ;phiq[i] = -99 ;eq[i] = -99 ;pdgidq[i] = -99 ;
        ptq_2[i] = -99 ;etaq_2[i] = -99 ;phiq_2[i] = -99 ;eq_2[i] = -99 ;pdgidq_2[i] = -99 ;
        ptq_3[i] = -99 ;etaq_3[i] = -99 ;phiq_3[i] = -99 ;eq_3[i] = -99 ;pdgidq_3[i] = -99 ;
    }

    fastJetRho = -99.0;

    IDLoose = false;
    IDTight = false;
    IDLoose_2 = false;
    IDTight_2 = false;
    IDLoose_3 = false;
    IDTight_3 = false;

    isHighPt = false;
    isHEEP = false;
    //rho = -99;
    iso = -99;
    isoCut = -99;
    et = -99;
    trackIso = -99;
    muchaiso=-99.;
    muneuiso=-99.;
    muphoiso=-99.;
    muPU=-99.;
    muisolation=-99.;

    METraw_et = -99;
    METraw_phi = -99;
    METraw_sumEt = -99;
    MET_et = -99;
    MET_phi = -99;
    MET_et_m = -99;
    MET_et_old = -99;
    MET_phi_m = -99;
    MET_sumEt = -99;

    MET_et_new=-99;
    MET_et_JEC_up=-99;
    MET_et_JEC_down=-99;
    MET_et_JER_up=-99;
    MET_et_JER_down=-99;
    MET_phi_new=-99;
    MET_phi_JEC_up=-99;
    MET_phi_JEC_down=-99;
    MET_phi_JER_up=-99;
    MET_phi_JER_down=-99;
    MET_sumEt_new=-99;
    MET_sumEt_JEC_up=-99;
    MET_sumEt_JEC_down=-99;
    MET_sumEt_JER_up=-99;
    MET_sumEt_JER_down=-99;

    candMasspuppiJEC     =  -99;
    m_jlv     =  -99;
    candMasspuppiJEC_new     =  -99;
    m_jlv_new     =  -99;
    
    candMasspuppiJEC_JEC_up     =  -99;
    m_jlv_JEC_up     =  -99;
    candMasspuppiJEC_JEC_down     =  -99;
    m_jlv_JEC_down     =  -99;
    candMasspuppiJEC_JER_up     =  -99;
    m_jlv_JER_up     =  -99;
    candMasspuppiJEC_JER_down     =  -99;
    m_jlv_JER_down     =  -99;

    ptVlepJEC       =  -99;
    yVlepJEC        =  -99;
    phiVlepJEC      =  -99;
    massVlepJEC     =  -99;
    mtVlepJEC       =  -99;
    ptVlepJEC_new       =  -99;
    yVlepJEC_new        =  -99;
    phiVlepJEC_new      =  -99;
    massVlepJEC_new     =  -99;
    mtVlepJEC_new       =  -99;
    ptVlepJEC_JEC_up       =  -99;
    yVlepJEC_JEC_up        =  -99;
    phiVlepJEC_JEC_up      =  -99;
    massVlepJEC_JEC_up     =  -99;
    mtVlepJEC_JEC_up       =  -99;
    ptVlepJEC_JEC_down       =  -99;
    yVlepJEC_JEC_down        =  -99;
    phiVlepJEC_JEC_down      =  -99;
    massVlepJEC_JEC_down     =  -99;
    mtVlepJEC_JEC_down       =  -99;
    ptVlepJEC_JER_down       =  -99;
    yVlepJEC_JER_down        =  -99;
    phiVlepJEC_JER_down      =  -99;
    massVlepJEC_JER_down     =  -99;
    mtVlepJEC_JER_down       =  -99;
    ptVlepJEC_JER_up       =  -99;
    yVlepJEC_JER_up        =  -99;
    phiVlepJEC_JER_up      =  -99;
    massVlepJEC_JER_up     =  -99;
    mtVlepJEC_JER_up       =  -99;

    massww[0] = -99;
    massww[1] = -99;
    massww[2] = -99;
    masslvj1 = -99;
    masslvj2 = -99;
    massj1j2 = -99;

    HLT_Ele1=-99;
    HLT_Ele2=-99;
    HLT_Ele3=-99;
    HLT_Ele4=-99;
    HLT_Ele5=-99;
    HLT_Ele6=-99;
    HLT_Ele7=-99;
    HLT_Ele8=-99;
    HLT_Mu1=-99;
    HLT_Mu2=-99;
    HLT_Mu3=-99;
    HLT_Mu4=-99;
    HLT_Mu5=-99;
    HLT_Mu6=-99;
    HLT_Mu7=-99;
    HLT_Mu8=-99;
    HLT_Mu9=-99;
    HLT_Mu10=-99;
    HLT_Mu11=-99;
    HLT_Mu12=-99;

    theWeight = -99;
    //nump = 0;
    //numm = 0;
    passFilter_HBHE_                  = false;
    passFilter_HBHEIso_               = false;
    passFilter_GlobalHalo_            = false;
    passFilter_ECALDeadCell_          = false;
    passFilter_GoodVtx_               = false;
    passFilter_EEBadSc_               = false;
    passFilter_badMuon_               = false;
    passFilter_badChargedHadron_      = false;
    passecalBadCalibFilterUpdate_     = false;
    for(int i=0;i<5;i++){
        ptgenwl[i]=-99;etagenwl[i]=-99;phigenwl[i]=-99;massgenwl[i]=-99;taggenwl[i]=-99;taggenwmother[i]=-99;
        genw_q1_pt[i]=-99;genw_q1_phi[i]=-99;genw_q1_eta[i]=-99;genw_q1_e[i]=-99;genw_q1_pdg[i]=-99;
        genw_q2_pt[i]=-99;genw_q2_phi[i]=-99;genw_q2_eta[i]=-99;genw_q2_e[i]=-99;genw_q2_pdg[i]=-99;
        ptgenzl[i]=-99;etagenzl[i]=-99;phigenzl[i]=-99;massgenzl[i]=-99;taggenzl[i]=-99;
        ptgenwf[i]=-99;etagenwf[i]=-99;phigenwf[i]=-99;massgenwf[i]=-99;
        ptgenzf[i]=-99;etagenzf[i]=-99;phigenzf[i]=-99;massgenzf[i]=-99;
    }
    for(int i=0;i<10;i++){
        ptgengl[i]=-99;etagengl[i]=-99;phigengl[i]=-99;egengl[i]=-99;
        ptgengf[i]=-99;etagengf[i]=-99;phigengf[i]=-99;egengf[i]=-99;
        mothergengf[i]=-99;mmothergengf[i]=-99;
    }
    for(int i=0;i<5;i++){
        ptgenq1l[i]=-99;etagenq1l[i]=-99;phigenq1l[i]=-99;egenq1l[i]=-99;
        ptgenq1f[i]=-99;etagenq1f[i]=-99;phigenq1f[i]=-99;egenq1f[i]=-99;
        ptgenq2l[i]=-99;etagenq2l[i]=-99;phigenq2l[i]=-99;egenq2l[i]=-99;
        ptgenq2f[i]=-99;etagenq2f[i]=-99;phigenq2f[i]=-99;egenq2f[i]=-99;
        ptgenq3l[i]=-99;etagenq3l[i]=-99;phigenq3l[i]=-99;egenq3l[i]=-99;
        ptgenq3f[i]=-99;etagenq3f[i]=-99;phigenq3f[i]=-99;egenq3f[i]=-99;
        ptgenq4l[i]=-99;etagenq4l[i]=-99;phigenq4l[i]=-99;egenq4l[i]=-99;
        ptgenq4f[i]=-99;etagenq4f[i]=-99;phigenq4f[i]=-99;egenq4f[i]=-99;
        ptgenq5l[i]=-99;etagenq5l[i]=-99;phigenq5l[i]=-99;egenq5l[i]=-99;
        ptgenq5f[i]=-99;etagenq5f[i]=-99;phigenq5f[i]=-99;egenq5f[i]=-99;
        mmothergenq1f[i]=-99;mmothergenq2f[i]=-99;mmothergenq3f[i]=-99;mmothergenq4f[i]=-99;mmothergenq5f[i]=-99;

    }
    gent_b_pt=-99;gent_b_phi=-99;gent_b_eta=-99;gent_b_mass=-99;
    genantit_b_pt=-99;genantit_b_phi=-99;genantit_b_eta=-99;genantit_b_mass=-99;
    gent_w_pt=-99;gent_w_phi=-99;gent_w_eta=-99;gent_w_mass=-99;
    genantit_w_pt=-99;genantit_w_phi=-99;genantit_w_eta=-99;genantit_w_mass=-99;
    gent_w_q1_pt=-99;gent_w_q1_phi=-99;gent_w_q1_eta=-99;gent_w_q1_e=-99;gent_w_q1_pdg=-99;
    genantit_w_q1_pt=-99;genantit_w_q1_phi=-99;genantit_w_q1_eta=-99;genantit_w_q1_e=-99;genantit_w_q1_pdg=-99;
    gent_w_q2_pt=-99;gent_w_q2_phi=-99;gent_w_q2_eta=-99;gent_w_q2_e=-99;gent_w_q2_pdg=-99;
    genantit_w_q2_pt=-99;genantit_w_q2_phi=-99;genantit_w_q2_eta=-99;genantit_w_q2_e=-99;genantit_w_q2_pdg=-99;
    gent_w_tag=-99;genantit_w_tag=-99;

}



#endif
