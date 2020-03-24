#ifndef _EDBRTreeMakerHead_
#define _EDBRTreeMakerHead_

using namespace edm;

class VVVTreeMaker : public edm::EDAnalyzer {
public:
    explicit VVVTreeMaker(const edm::ParameterSet&);
    ~VVVTreeMaker();
  
private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
    virtual void endRun(const edm::Run&, const edm::EventSetup&) override;

    virtual bool looseJetID( const pat::Jet& j);
    virtual const reco::Candidate* findLastW(const reco::Candidate *particle,int IDpdg);
    virtual const reco::Candidate* findLasttau(const reco::Candidate *particle,int IDpdg);
    virtual const reco::Candidate* findFirstW(const reco::Candidate *particle,int IDpdg);
    virtual bool tightJetID( const pat::Jet& j);
    virtual bool tightJetIDpuppi( const pat::Jet& j);
    virtual float dEtaInSeed( const pat::Electron* ele) ;
    virtual void initJetCorrFactors( void );
    virtual void addTypeICorr( edm::Event const & event );
    virtual void   addTypeICorr_user(edm::Event const& event);  //---for MET,
    virtual double getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
    virtual double getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
    math::XYZTLorentzVector getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType);

    void HLTStore( edm::Event const & );
    void GENStore( edm::Event const & );
    void PuppiAK8_Store( edm::Event const & );
    void MET_Store( edm::Event const & , edm::Handle<edm::View<reco::Candidate> > );
    void EDBR_Debug();
    void Weight_Store( edm::Event const & );
    void LeptonicW_Store( edm::Event const &, edm::Handle<edm::View<reco::Candidate> >, edm::Handle<edm::View<reco::Candidate> > );
    void Filter_Store( edm::Event const & );
    void Lepton_num_Store( edm::Event const & );
    void Lepton_fromLeptonicW_Store( edm::Handle<edm::View<reco::Candidate> > , edm::Handle<reco::VertexCollection>);
    void VVV_AK4chs_Store( edm::Event const & );
    void VVV_inv_mass_Store( void );



    
    boost::shared_ptr<FactorizedJetCorrector>   jecAK8_            ;

    boost::shared_ptr<FactorizedJetCorrector>   jecAK8Groomed_            ;
    boost::shared_ptr<FactorizedJetCorrector>   jecAK8GroomedSD_            ;

    std::vector<std::string>                    jecAK8puppiPayloadNames_;
    boost::shared_ptr<FactorizedJetCorrector>   jecAK8puppi_            ;


    
    std::vector<std::string>                    jecAK4PayloadNames_;
    boost::shared_ptr<FactorizedJetCorrector>   jecAK4_            ;

    
    boost::shared_ptr<FactorizedJetCorrector> jecOffset_;
    edm::Handle< double >  rho_;
    edm::InputTag  METsRawLabel_;
    edm::Handle<pat::METCollection>  METs_;
    edm::Handle<pat::JetCollection> jets_;
    edm::Handle<reco::VertexCollection> vertices_;

    edm::Handle<pat::METCollection>  reclusteredMETs_;
    edm::Handle<edm::View<reco::PFMET> >     pfMET_ ;
    edm::EDGetTokenT<pat::JetCollection> prunedjetInputToken_;
    edm::EDGetTokenT<pat::JetCollection> softdropjetInputToken_;
    edm::EDGetTokenT<pat::JetCollection> fatjetInputToken_;
    edm::EDGetTokenT<pat::JetCollection> puppijetInputToken_;

    // add 3 up
    edm::EDGetTokenT<pat::METCollection>  metInputToken_;
    edm::EDGetTokenT<pat::METCollection>  reclusteredmetInputToken_;
    std::vector<edm::EDGetTokenT<pat::METCollection>> mettokens;
    edm::Handle<pat::JetCollection> puppijets_;

    // add 2 up

    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<pat::METCollection> reclusteredmetToken_;
    edm::EDGetTokenT<pat::JetCollection> jetToken_;



    edm::EDGetTokenT<pat::JetCollection> puppijetToken_;

    // add 4 up
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    std::vector<std::string> jetCorrLabel_;
    std::vector<std::string> jecAK4Labels;
    std::vector<std::string> jecAK8Labels;
    std::vector<std::string> offsetCorrLabel_;
    bool doCorrOnTheFly_;
    // Filter
    edm::EDGetTokenT<edm::TriggerResults> 		     noiseFilterToken_;
    edm::Handle< edm::TriggerResults> 			     noiseFilterBits_;
    std::string HBHENoiseFilter_Selector_;
    std::string HBHENoiseIsoFilter_Selector_;
    std::string GlobalHaloNoiseFilter_Selector_;
    std::string ECALDeadCellNoiseFilter_Selector_;
    std::string GoodVtxNoiseFilter_Selector_;
    std::string EEBadScNoiseFilter_Selector_;
    edm::EDGetTokenT<bool>  badMuon_Selector_;
    edm::EDGetTokenT<bool>  badChargedHadron_Selector_;
    edm::EDGetTokenT<bool>  ecalBadCalibFilterUpdate_token ;

  
    // ----------member data ---------------------------
    TTree* outTree_;
    TTree* outTreew_;

    double MW_;
    int nmetmatch, nmetno;
    int nevent, run, ls;
    int nVtx;

    int nLooseEle, nLooseMu;//Synch
    double ptVlep, yVlep, phiVlep, massVlep;
    double met, metPhi, mtVlep;
    double jetAK8puppi_ptJEC, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_tau1,  jetAK8puppi_tau2, jetAK8puppi_tau3, jetAK8puppi_tau4,jetAK8puppi_tau21, jetAK8puppi_tau42, jetAK8puppi_sd;
    double jetAK8puppi_ptJEC_2, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_tau1_2,  jetAK8puppi_tau2_2, jetAK8puppi_tau3_2,jetAK8puppi_tau4_2, jetAK8puppi_tau21_2,jetAK8puppi_tau42_2,  jetAK8puppi_sd_2;
    double jetAK8puppi_ptJEC_3, jetAK8puppi_eta_3, jetAK8puppi_phi_3, jetAK8puppi_tau1_3,  jetAK8puppi_tau2_3, jetAK8puppi_tau3_3,jetAK8puppi_tau4_3, jetAK8puppi_tau21_3,jetAK8puppi_tau42_3,  jetAK8puppi_sd_3;
    double jetAK8puppi_dnnTop, jetAK8puppi_dnnW,jetAK8puppi_dnnH4q,jetAK8puppi_dnnTop_2, jetAK8puppi_dnnW_2,jetAK8puppi_dnnH4q_2,jetAK8puppi_dnnTop_3, jetAK8puppi_dnnW_3,jetAK8puppi_dnnH4q_3; //DeepAK8
    double jetAK8puppi_dnnqcd,jetAK8puppi_dnntop,jetAK8puppi_dnnw,jetAK8puppi_dnnz,jetAK8puppi_dnnzbb,jetAK8puppi_dnnhbb,jetAK8puppi_dnnh4q,jetAK8puppi_dnnqcd_2,jetAK8puppi_dnntop_2,jetAK8puppi_dnnw_2,jetAK8puppi_dnnz_2,jetAK8puppi_dnnzbb_2,jetAK8puppi_dnnhbb_2,jetAK8puppi_dnnh4q_2,jetAK8puppi_dnnqcd_3,jetAK8puppi_dnntop_3,jetAK8puppi_dnnw_3,jetAK8puppi_dnnz_3,jetAK8puppi_dnnzbb_3,jetAK8puppi_dnnhbb_3,jetAK8puppi_dnnh4q_3;
    double jetAK8puppi_dnnZ,jetAK8puppi_dnnZbb,jetAK8puppi_dnnHbb,jetAK8puppi_dnnZ_2,jetAK8puppi_dnnZbb_2,jetAK8puppi_dnnHbb_2,jetAK8puppi_dnnZ_3,jetAK8puppi_dnnZbb_3,jetAK8puppi_dnnHbb_3;
    double jetAK8puppi_dnnDecorrTop, jetAK8puppi_dnnDecorrW,jetAK8puppi_dnnDecorrH4q,jetAK8puppi_dnnDecorrTop_2, jetAK8puppi_dnnDecorrW_2, jetAK8puppi_dnnDecorrH4q_2,jetAK8puppi_dnnDecorrTop_3, jetAK8puppi_dnnDecorrW_3, jetAK8puppi_dnnDecorrH4q_3; //Decorrelated DeepAK8
    double jetAK8puppi_dnnDecorrZ,jetAK8puppi_dnnDecorrZbb,jetAK8puppi_dnnDecorrHbb,jetAK8puppi_dnnDecorrZ_2,jetAK8puppi_dnnDecorrZbb_2,jetAK8puppi_dnnDecorrHbb_2,jetAK8puppi_dnnDecorrZ_3,jetAK8puppi_dnnDecorrZbb_3,jetAK8puppi_dnnDecorrHbb_3;
    double jetAK8puppi_dnnDecorrbb,jetAK8puppi_dnnDecorrcc,jetAK8puppi_dnnDecorrbbnog,jetAK8puppi_dnnDecorrccnog,jetAK8puppi_dnnDecorrbb_2,jetAK8puppi_dnnDecorrcc_2,jetAK8puppi_dnnDecorrbbnog_2,jetAK8puppi_dnnDecorrccnog_2,jetAK8puppi_dnnDecorrbb_3,jetAK8puppi_dnnDecorrcc_3,jetAK8puppi_dnnDecorrbbnog_3,jetAK8puppi_dnnDecorrccnog_3;
    double jetAK8puppi_dnnDecorrqcd,jetAK8puppi_dnnDecorrtop,jetAK8puppi_dnnDecorrw,jetAK8puppi_dnnDecorrz,jetAK8puppi_dnnDecorrzbb,jetAK8puppi_dnnDecorrhbb,jetAK8puppi_dnnDecorrh4q,jetAK8puppi_dnnDecorrqcd_2,jetAK8puppi_dnnDecorrtop_2,jetAK8puppi_dnnDecorrw_2,jetAK8puppi_dnnDecorrz_2,jetAK8puppi_dnnDecorrzbb_2,jetAK8puppi_dnnDecorrhbb_2,jetAK8puppi_dnnDecorrh4q_2,jetAK8puppi_dnnDecorrqcd_3,jetAK8puppi_dnnDecorrtop_3,jetAK8puppi_dnnDecorrw_3,jetAK8puppi_dnnDecorrz_3,jetAK8puppi_dnnDecorrzbb_3,jetAK8puppi_dnnDecorrhbb_3,jetAK8puppi_dnnDecorrh4q_3;
    double jetAK8puppi_ptJEC_new,jetAK8puppi_ptJEC_JEC_up,jetAK8puppi_ptJEC_JEC_down,jetAK8puppi_ptJEC_JER_up,jetAK8puppi_ptJEC_JER_down,jetAK8puppi_ptJEC_newnew,jetAK8puppi_ptJEC_m;
    double jetAK8puppi_ptJEC_2_new,jetAK8puppi_ptJEC_2_JEC_up,jetAK8puppi_ptJEC_2_JEC_down,jetAK8puppi_ptJEC_2_JER_up,jetAK8puppi_ptJEC_2_JER_down;
    double jetAK8puppi_ptJEC_3_new,jetAK8puppi_ptJEC_3_JEC_up,jetAK8puppi_ptJEC_3_JEC_down,jetAK8puppi_ptJEC_3_JER_up,jetAK8puppi_ptJEC_3_JER_down;
    
    double jetAK8puppi_e,jetAK8puppi_e_new,jetAK8puppi_e_JEC_up,jetAK8puppi_e_JEC_down,jetAK8puppi_e_JER_up,jetAK8puppi_e_JER_down;
    double jetAK8puppi_e_2,jetAK8puppi_e_2_new,jetAK8puppi_e_2_JEC_up,jetAK8puppi_e_2_JEC_down,jetAK8puppi_e_2_JER_up,jetAK8puppi_e_2_JER_down;
    double jetAK8puppi_e_3,jetAK8puppi_e_3_new,jetAK8puppi_e_3_JEC_up,jetAK8puppi_e_3_JEC_down,jetAK8puppi_e_3_JER_up,jetAK8puppi_e_3_JER_down;

    double fastJetRho;

    double ptgenwl[5],etagenwl[5],phigenwl[5],massgenwl[5],taggenwl[5],taggenwmother[5];
    double genw_q1_pt[5],genw_q1_eta[5],genw_q1_phi[5],genw_q1_e[5],genw_q1_pdg[5];
    double genw_q2_pt[5],genw_q2_eta[5],genw_q2_phi[5],genw_q2_e[5],genw_q2_pdg[5];
    double ptgenzl[5],etagenzl[5],phigenzl[5],massgenzl[5],taggenzl[5];
    double ptgengl[10],etagengl[10],phigengl[10],egengl[10];
    double ptgenwf[5],etagenwf[5],phigenwf[5],massgenwf[5];
    double ptgenzf[5],etagenzf[5],phigenzf[5],massgenzf[5];
    double ptgengf[10],etagengf[10],phigengf[10],egengf[10];
    double gent_b_pt,gent_b_phi,gent_b_eta,gent_b_mass;
    double genantit_b_pt,genantit_b_phi,genantit_b_eta,genantit_b_mass;
    double gent_w_pt,gent_w_phi,gent_w_eta,gent_w_mass;
    double genantit_w_pt,genantit_w_phi,genantit_w_eta,genantit_w_mass;
    double gent_w_q1_pt,gent_w_q1_phi,gent_w_q1_eta,gent_w_q1_e,gent_w_q1_pdg;
    double genantit_w_q1_pt,genantit_w_q1_phi,genantit_w_q1_eta,genantit_w_q1_e,genantit_w_q1_pdg;
    double gent_w_q2_pt,gent_w_q2_phi,gent_w_q2_eta,gent_w_q2_e,gent_w_q2_pdg;
    double genantit_w_q2_pt,genantit_w_q2_phi,genantit_w_q2_eta,genantit_w_q2_e,genantit_w_q2_pdg;
    double ptgenq1l[5],etagenq1l[5],phigenq1l[5],egenq1l[5];
    double ptgenq1f[5],etagenq1f[5],phigenq1f[5],egenq1f[5];
    double ptgenq2l[5],etagenq2l[5],phigenq2l[5],egenq2l[5];
    double ptgenq2f[5],etagenq2f[5],phigenq2f[5],egenq2f[5];
    double ptgenq3l[5],etagenq3l[5],phigenq3l[5],egenq3l[5];
    double ptgenq3f[5],etagenq3f[5],phigenq3f[5],egenq3f[5];
    double ptgenq4l[5],etagenq4l[5],phigenq4l[5],egenq4l[5];
    double ptgenq4f[5],etagenq4f[5],phigenq4f[5],egenq4f[5];
    double ptgenq5l[5],etagenq5l[5],phigenq5l[5],egenq5l[5];
    double ptgenq5f[5],etagenq5f[5],phigenq5f[5],egenq5f[5];
    double mothergenq1f[5],mothergenq2f[5],mothergenq3f[5],mothergenq4f[5],mothergenq5f[5];
    
    double gent_w_tag,genantit_w_tag,mothergengf[10],mmothergengf[10],mmothergenq1f[5],mmothergenq2f[5],mmothergenq3f[5],mmothergenq4f[5],mmothergenq5f[5];
    
    double vbfeta, vbfmjj;
    int      vbftag;
    int nj1, nj2;
    int numq,numq_2,numq_3;
    double ptlep1, ptlep2;
    double etalep1, etalep2 ;
    double philep1, philep2 ;
    double energylep1;

    int channel, lep;
    double deltaRlepjet, delPhilepmet, delPhijetmet, delPhijetlep;
    double deltaRlepjet_2,  delPhijetmet_2, delPhijetlep_2;

    double pt_graviton,pt_graviton1;
    double ptVlepJEC, yVlepJEC, phiVlepJEC;
    double ptVlepJEC_new, yVlepJEC_new, phiVlepJEC_new,massVlepJEC_new, mtVlepJEC_new;
    double ptVlepJEC_JEC_up, yVlepJEC_JEC_up, phiVlepJEC_JEC_up,massVlepJEC_JEC_up, mtVlepJEC_JEC_up;
    double ptVlepJEC_JEC_down, yVlepJEC_JEC_down, phiVlepJEC_JEC_down,massVlepJEC_JEC_down, mtVlepJEC_JEC_down;
    double ptVlepJEC_JER_up, yVlepJEC_JER_up, phiVlepJEC_JER_up,massVlepJEC_JER_up, mtVlepJEC_JER_up;
    double ptVlepJEC_JER_down, yVlepJEC_JER_down, phiVlepJEC_JER_down,massVlepJEC_JER_down, mtVlepJEC_JER_down;

    double candMasspuppiJEC,m_jlv;
    double candMasspuppiJEC_new,m_jlv_new,candMasspuppiJEC_JEC_up,m_jlv_JEC_up,candMasspuppiJEC_JEC_down,m_jlv_JEC_down,candMasspuppiJEC_JER_up,m_jlv_JER_up,candMasspuppiJEC_JER_down,m_jlv_JER_down;

    TLorentzVector  gleptonicV,gleptonicV_new,gleptonicV_JEC_up,gleptonicV_JEC_down,gleptonicV_JER_up,gleptonicV_JER_down;

    double massww[3],masslvj1,masslvj2,massj1j2;
    double massVlepJEC, mtVlepJEC;

    double theWeight;
    double  nump=0;
    double  numm=0;
    //double pweight[882];
    double  npT, npIT;
    int     nBX;
    //Gen Level

    double gen_ele_pt, gen_ele_eta, gen_ele_phi, gen_ele_e;


    double pttau[4],etatau[4],phitau[4],etau[4],pdgidtau[4];
    double pttau_2[4],etatau_2[4],phitau_2[4],etau_2[4],pdgidtau_2[4];
    double pttau_3[4],etatau_3[4],phitau_3[4],etau_3[4],pdgidtau_3[4];
   
    double ptq[3],etaq[3],phiq[3],eq[3],pdgidq[3];
    double ptq_2[3],etaq_2[3],phiq_2[3],eq_2[3],pdgidq_2[3];
    double ptq_3[3],etaq_3[3],phiq_3[3],eq_3[3],pdgidq_3[3];





    double gen_mu_pt, gen_mu_eta, gen_mu_phi, gen_mu_e;
    double genmatch_ele_pt, genmatch_ele_eta, genmatch_ele_phi, genmatch_ele_e, genmatch_ele_dr;
    double genmatch_mu_pt, genmatch_mu_eta, genmatch_mu_phi, genmatch_mu_e, genmatch_mu_dr;
    double gen_ele_pt_2, gen_ele_eta_2, gen_ele_phi_2, gen_ele_e_2;
    double gen_mu_pt_2, gen_mu_eta_2, gen_mu_phi_2, gen_mu_e_2;
    double gen_ele_pt_3, gen_ele_eta_3, gen_ele_phi_3, gen_ele_e_3;
    double gen_mu_pt_3, gen_mu_eta_3, gen_mu_phi_3, gen_mu_e_3;
    double gentop_pt, gentop_eta, gentop_phi, gentop_mass;
    double genantitop_pt, genantitop_eta, genantitop_phi, genantitop_mass;



    
    bool IDLoose, IDTight,IDLoose_2, IDTight_2,IDLoose_3, IDTight_3, isHighPt, isHEEP;
    double muchaiso, muneuiso, muphoiso, muPU, muisolation;
    double iso, isoCut, et, trackIso;
    //  double rho,fastJetRho;
    double useless;
    //  JEC
    double corr_AK8puppi[4];
    double jetAK8puppi_pt1[4], jetAK8puppi_mass1[4], jetAK8puppi_eta1[4], jetAK8puppi_jec1[4], jetAK8puppiSD_jec1[4];
    double jetAK8puppi_pt1_new[4],jetAK8puppi_pt1_JEC_up[4],jetAK8puppi_pt1_JEC_down[4],jetAK8puppi_pt1_JER_up[4],jetAK8puppi_pt1_JER_down[4],jetAK8puppi_pt1_newnew[4],jetAK8puppi_pt1_m[4];
    double jetAK8puppi_e1_new[4],jetAK8puppi_e1_JEC_up[4],jetAK8puppi_e1_JEC_down[4],jetAK8puppi_e1_JER_up[4],jetAK8puppi_e1_JER_down[4];

    double corr;
    double METraw_et, METraw_phi, METraw_sumEt;
    double MET_et, MET_phi, MET_sumEt;
    double MET_et_new, MET_phi_new, MET_sumEt_new;
    double MET_et_m,MET_et_old, MET_phi_m, MET_sumEt_m;
    // Marked for debug
    //-------------- Met uncertainty ----------------//
    double MET_et_JEC_up, MET_et_JEC_down, MET_et_JER_up, MET_et_JER_down;
    double MET_phi_JEC_up, MET_phi_JEC_down, MET_phi_JER_up, MET_phi_JER_down;
    double MET_sumEt_JEC_up, MET_sumEt_JEC_down, MET_sumEt_JER_up, MET_sumEt_JER_down;
    // AK4 Jets
    int ak4jet_hf[8],ak4jet_pf[8],ak4jet_hf_2[8],ak4jet_pf_2[8];
    double ak4jet_pt[8],ak4jet_pt_uncorr[8],ak4jet_eta[8],ak4jet_phi[8],ak4jet_e[8], ak4jet_dr[8];
    double ak4jet_csv[8],ak4jet_icsv[8], ak4jet_IDLoose[8], ak4jet_IDTight[8],ak4jet_deepcsvudsg[8],ak4jet_deepcsvb[8],ak4jet_deepcsvc[8],ak4jet_deepcsvbb[8],ak4jet_deepcsvcc[8];
    TLorentzVector ak8sj11,ak8sj12,ak8sj13,ak8sj14,ak8sj15,puppi_softdropj1;
    TLorentzVector ak8sj21,ak8sj22,ak8sj23,ak8sj24,ak8sj25,puppi_softdropj2;

    void setDummyValues();
    
    //// L1 prefiring
    edm::EDGetTokenT< double > prefweight_token;
    edm::EDGetTokenT< double > prefweightup_token;
    edm::EDGetTokenT< double > prefweightdown_token;
    
    /// Parameters to steer the treeDumper

    std::string EDBRChannel_;
    bool isGen_;
    bool isJEC_;
    bool RunOnSig_,RunOnMC_;
    bool CoutOrNot_;

    std::vector<JetCorrectorParameters> vPar;
    std::map<std::string,double>  TypeICorrMap_;
    std::map<std::string, double> TypeICorrMap_user_;

    edm::InputTag mets_;

    //High Level Trigger
    HLTConfigProvider hltConfig;
    edm::EDGetTokenT<edm::TriggerResults> hltToken_;
    std::vector<std::string> elPaths1_, elPaths2_, elPaths3_, elPaths4_, elPaths5_, elPaths6_, elPaths7_, elPaths8_;
    std::vector<std::string> muPaths1_, muPaths2_, muPaths3_, muPaths4_, muPaths5_, muPaths6_, muPaths7_, muPaths8_, muPaths9_, muPaths10_, muPaths11_, muPaths12_;
    std::vector<std::string> elPaths1, elPaths2, elPaths3, elPaths4, elPaths5, elPaths6, elPaths7, elPaths8;
    std::vector<std::string> muPaths1, muPaths2, muPaths3, muPaths4, muPaths5, muPaths6, muPaths7, muPaths8, muPaths9, muPaths10, muPaths11, muPaths12;
    int  HLT_Ele1, HLT_Ele2, HLT_Ele3, HLT_Ele4, HLT_Ele5, HLT_Ele6, HLT_Ele7, HLT_Ele8;
    int  HLT_Mu1, HLT_Mu2, HLT_Mu3, HLT_Mu4, HLT_Mu5, HLT_Mu6, HLT_Mu7, HLT_Mu8, HLT_Mu9, HLT_Mu10, HLT_Mu11,  HLT_Mu12;

    //L1 prefiring
    double L1prefiring,L1prefiringup,L1prefiringdown;
    
    // filter
    bool passFilter_HBHE_                   ;
    bool passFilter_HBHEIso_                ;
    bool passFilter_GlobalHalo_             ;
    bool passFilter_ECALDeadCell_           ;
    bool passFilter_GoodVtx_                ;
    bool passFilter_EEBadSc_                ;
    bool passFilter_badMuon_                ;
    bool passFilter_badChargedHadron_       ;
    bool passecalBadCalibFilterUpdate_      ;

    edm::EDGetTokenT<edm::View<reco::Candidate>> leptonicVSrc_;

    edm::EDGetTokenT<edm::View<pat::Jet>> ak4jetsSrc_;

    edm::EDGetTokenT<LHEEventProduct> LheToken_;
    edm::EDGetTokenT<LHERunInfoProduct> LhestrToken_;

    edm::EDGetTokenT<edm::View<pat::Electron> > looseelectronToken_ ;
    edm::EDGetTokenT<edm::View<pat::Muon>> loosemuonToken_;
    edm::EDGetTokenT<edm::View<pat::Muon>> goodMuSrc_;
    edm::EDGetTokenT<edm::View<pat::Muon>> MuSrc_;
    edm::EDGetTokenT<edm::View<pat::Electron> > EleSrc_;

    edm::EDGetTokenT<pat::JetCollection>             t1jetSrc_userak4_;
    edm::EDGetTokenT<edm::View<reco::Candidate>> metSrc_;
    edm::EDGetTokenT<GenEventInfoProduct> GenToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle>> genSrc_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUToken_;
};


    //
    // constructors and destructor


VVVTreeMaker::VVVTreeMaker(const edm::ParameterSet& iConfig):
    hltToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("hltToken"))),
    elPaths1_(iConfig.getParameter<std::vector<std::string>>("elPaths1")),
    elPaths2_(iConfig.getParameter<std::vector<std::string>>("elPaths2")),
    elPaths3_(iConfig.getParameter<std::vector<std::string>>("elPaths3")),
    elPaths4_(iConfig.getParameter<std::vector<std::string>>("elPaths4")),
    elPaths5_(iConfig.getParameter<std::vector<std::string>>("elPaths6")),
    elPaths6_(iConfig.getParameter<std::vector<std::string>>("elPaths5")),
    elPaths7_(iConfig.getParameter<std::vector<std::string>>("elPaths7")),
    elPaths8_(iConfig.getParameter<std::vector<std::string>>("elPaths8")),
    muPaths1_(iConfig.getParameter<std::vector<std::string>>("muPaths1")),
    muPaths2_(iConfig.getParameter<std::vector<std::string>>("muPaths2")),
    muPaths3_(iConfig.getParameter<std::vector<std::string>>("muPaths3")),
    muPaths4_(iConfig.getParameter<std::vector<std::string>>("muPaths4")),
    muPaths5_(iConfig.getParameter<std::vector<std::string>>("muPaths5")),
    muPaths6_(iConfig.getParameter<std::vector<std::string>>("muPaths6")),
    muPaths7_(iConfig.getParameter<std::vector<std::string>>("muPaths7")),
    muPaths8_(iConfig.getParameter<std::vector<std::string>>("muPaths8")),
    muPaths9_(iConfig.getParameter<std::vector<std::string>>("muPaths9")),
    muPaths10_(iConfig.getParameter<std::vector<std::string>>("muPaths10")),
    muPaths11_(iConfig.getParameter<std::vector<std::string>>("muPaths11")),
    muPaths12_(iConfig.getParameter<std::vector<std::string>>("muPaths12"))
    //  noiseFilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter")))
{
    LheToken_=consumes<LHEEventProduct> (iConfig.getParameter<edm::InputTag>( "lhe") ) ;
    LhestrToken_=consumes<LHERunInfoProduct,edm::InRun> (iConfig.getParameter<edm::InputTag>( "lhe") ) ;

    EDBRChannel_     = iConfig.getParameter<std::string>("EDBRChannel");
    isGen_           = iConfig.getParameter<bool>("isGen");
    isJEC_           = iConfig.getParameter<bool>("isJEC");
    RunOnSig_        = iConfig.getParameter<bool>("RunOnSig");
    RunOnMC_        = iConfig.getParameter<bool>("RunOnMC");
    CoutOrNot_        = iConfig.getParameter<bool>("CoutOrNot");
    // Sources
    //  leptonicVSrc_ = iConfig.getParameter<std::string>("leptonicVSrc");
    leptonicVSrc_=consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>( "leptonicVSrc") ) ;
    looseelectronToken_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("looseElectronSrc"))) ;
    loosemuonToken_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("looseMuonSrc")));

    goodMuSrc_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("goodMuSrc")));
    MuSrc_    = (consumes<edm::View<pat::Muon> > (iConfig.getParameter<edm::InputTag>("MuSrc")));
    EleSrc_    = (consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>("EleSrc")));


    //  goodMuSrc_    = iConfig.getParameter<std::string>("goodMuSrc");
    //  looseMuonSrc_    = iConfig.getParameter<std::string>("looseMuonSrc");
    //  looseElectronsSrc_= iConfig.getParameter<std::string>("looseElectronsSrc");
    muonToken_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
    //  ak4jetsSrc_      = iConfig.getParameter<std::string>("ak4jetsSrc");
    ak4jetsSrc_      = consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>( "ak4jetsSrc") ) ;



    jetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
    puppijetToken_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("puppijets"));



    // add 4 up
    rhoToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
    vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
    GenToken_=consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "generator") ) ;
    genSrc_      = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>( "genSrc") ) ;
    PUToken_=consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileup") ) ;
    t1jetSrc_userak4_   = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("t1jetSrc_userak4"));

    metSrc_      = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>( "metSrc") ) ;
    metToken_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metSrc"));

    //  L1 prefiring
    prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
    prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
    prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
    
    // filter
    noiseFilterToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter"));
    HBHENoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseFilter");
    HBHENoiseIsoFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseIsoFilter");
    GlobalHaloNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_GlobalTightHaloFilter");
    ECALDeadCellNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter");
    GoodVtxNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_goodVertices");
    EEBadScNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_eeBadScFilter");
    badMuon_Selector_ =  consumes<bool>(iConfig.getParameter<edm::InputTag> ("noiseFilterSelection_badMuon"));
    badChargedHadron_Selector_ =  consumes<bool>(iConfig.getParameter<edm::InputTag> ("noiseFilterSelection_badChargedHadron"));
    ecalBadCalibFilterUpdate_token= consumes< bool >(edm::InputTag("ecalBadCalibReducedMINIAODFilter"));

    std::string jecpath = iConfig.getParameter<std::string>("jecpath");
    std::string tmpString;

    std::vector<std::string> tmpVec;


    std::vector<std::string> jecAK8Labelspuppi;
    tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8puppiPayloadNames");
    for( unsigned int v = 0; v < tmpVec.size(); ++v ){
        tmpString = jecpath + tmpVec[v];
        jecAK8Labelspuppi.push_back(tmpString);
    }



    std::vector<std::string> jecAK4Labels;
    tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK4chsPayloadNames");
    for( unsigned int v = 0; v < tmpVec.size(); ++v ){
        tmpString = jecpath + tmpVec[v];
        jecAK4Labels.push_back(tmpString);
    }

    /*=======================================================================================*/
    MW_=80.385;
    nmetmatch = 0;
    nmetno = 0;
    mettokens.push_back( metToken_ );
    mettokens.push_back( reclusteredmetToken_ );

    
    // add 3 up
    metInputToken_ = mettokens[0];
    reclusteredmetInputToken_ = mettokens[1];

    jetCorrLabel_ = jecAK4Labels;
    jetCorrLabel_ = jecAK4Labels;
    offsetCorrLabel_.push_back(jetCorrLabel_[0]);
 
    doCorrOnTheFly_ = false;
    if( jecAK4Labels.size() != 0 && jecAK8Labelspuppi.size() != 0 ){

        jecAK4PayloadNames_ = jecAK4Labels;
        jecAK8puppiPayloadNames_ = jecAK8Labelspuppi;

        puppijetInputToken_ = puppijetToken_;

        initJetCorrFactors();
        doCorrOnTheFly_ = true;
    }

    if(EDBRChannel_ == "VZ_CHANNEL")
        channel=VZ_CHANNEL;
    else if(EDBRChannel_ == "VW_CHANNEL")
        channel=VW_CHANNEL;
    else if(EDBRChannel_ == "VH_CHANNEL")
        channel=VH_CHANNEL;
    else {
        cms::Exception ex("InvalidConfiguration");
        ex << "Unknown channel " << EDBRChannel_<< ". Please check EDBRTreeMaker.cc for allowed values.";
    throw ex;
    }
  

    //now do what ever initialization is needed
    edm::Service<TFileService> fs;

    outTree_ = fs->make<TTree>("EDBRCandidates","EDBR Candidates");
    outTreew_ = fs->make<TTree>("EDBRCandidatesw","EDBR Candidates");
    outTree_->Branch("L1prefiring"           ,&L1prefiring         ,"L1prefiring/D"          );
    outTree_->Branch("L1prefiringup"           ,&L1prefiringup         ,"L1prefiringup/D"          );
    outTree_->Branch("L1prefiringdown"           ,&L1prefiringdown         ,"L1prefiringdown/D"          );
    /// Basic event quantities
    if (RunOnMC_){
        //outTree_->Branch("pweight"           ,pweight         ,"pweight[882]/D"          );
        outTree_->Branch("ptgenwl"           ,ptgenwl         ,"ptgenwl[5]/D"          );
        outTree_->Branch("etagenwl"           ,etagenwl         ,"etagenwl[5]/D"          );
        outTree_->Branch("phigenwl"           ,phigenwl       ,"phigenwl[5]/D"          );
        outTree_->Branch("massgenwl"           ,massgenwl         ,"massgenwl[5]/D"          );
        outTree_->Branch("taggenwl"           ,taggenwl         ,"taggenwl[5]/D"          );
        outTree_->Branch("taggenwmother"           ,taggenwmother         ,"taggenwmother[5]/D"          );
        outTree_->Branch("genw_q1_pt"           ,genw_q1_pt         ,"genw_q1_pt[5]/D"          );
        outTree_->Branch("genw_q1_phi"           ,genw_q1_phi         ,"genw_q1_phi[5]/D"          );
        outTree_->Branch("genw_q1_eta"           ,genw_q1_eta         ,"genw_q1_eta[5]/D"          );
        outTree_->Branch("genw_q1_e"           ,genw_q1_e         ,"genw_q1_e[5]/D"          );
        outTree_->Branch("genw_q1_pdg"           ,genw_q1_pdg         ,"genw_q1_pdg[5]/D"          );
        outTree_->Branch("genw_q2_pt"           ,genw_q2_pt         ,"genw_q2_pt[5]/D"          );
        outTree_->Branch("genw_q2_phi"           ,genw_q2_phi         ,"genw_q2_phi[5]/D"          );
        outTree_->Branch("genw_q2_eta"           ,genw_q2_eta         ,"genw_q2_eta[5]/D"          );
        outTree_->Branch("genw_q2_e"           ,genw_q2_e         ,"genw_q2_e[5]/D"          );
        outTree_->Branch("genw_q2_pdg"           ,genw_q2_pdg         ,"genw_q2_pdg[5]/D"          );

        outTree_->Branch("ptgenzl"           ,ptgenzl         ,"ptgenzl[5]/D"          );
        outTree_->Branch("etagenzl"           ,etagenzl         ,"etagenzl[5]/D"          );
        outTree_->Branch("phigenzl"           ,phigenzl       ,"phigenzl[5]/D"          );
        outTree_->Branch("massgenzl"           ,massgenzl         ,"massgenzl[5]/D"          );
        outTree_->Branch("taggenzl"           ,taggenzl         ,"taggenzl[5]/D"          );
        outTree_->Branch("ptgengl"           ,ptgengl         ,"ptgengl[10]/D"          );
        outTree_->Branch("etagengl"           ,etagengl         ,"etagengl[10]/D"          );
        outTree_->Branch("phigengl"           ,phigengl       ,"phigengl[10]/D"          );
        outTree_->Branch("egengl"           ,egengl         ,"egengl[10]/D"          );
        outTree_->Branch("ptgenwf"           ,ptgenwf         ,"ptgenwf[5]/D"          );
        outTree_->Branch("etagenwf"           ,etagenwf         ,"etagenwf[5]/D"          );
        outTree_->Branch("phigenwf"           ,phigenwf       ,"phigenwf[5]/D"          );
        outTree_->Branch("massgenwf"           ,massgenwf         ,"massgenwf[5]/D"          );
        outTree_->Branch("ptgenzf"           ,ptgenzf         ,"ptgenzf[5]/D"          );
        outTree_->Branch("etagenzf"           ,etagenzf         ,"etagenzf[5]/D"          );
        outTree_->Branch("phigenzf"           ,phigenzf       ,"phigenzf[5]/D"          );
        outTree_->Branch("massgenzf"           ,massgenzf         ,"massgenzf[5]/D"          );
        outTree_->Branch("ptgengf"           ,ptgengf         ,"ptgengf[10]/D"          );
        outTree_->Branch("etagengf"           ,etagengf         ,"etagengf[10]/D"          );
        outTree_->Branch("phigengf"           ,phigengf       ,"phigengf[10]/D"          );
        outTree_->Branch("egengf"           ,egengf         ,"egengf[10]/D"          );
        
        outTree_->Branch("gent_b_pt"           ,&gent_b_pt         ,"gent_b_pt/D"          );
        outTree_->Branch("gent_b_eta"           ,&gent_b_eta         ,"gent_b_eta/D"          );
        outTree_->Branch("gent_b_phi"           ,&gent_b_phi         ,"gent_b_phi/D"          );
        outTree_->Branch("gent_b_mass"           ,&gent_b_mass         ,"gent_b_mass/D"          );
        outTree_->Branch("genantit_b_pt"           ,&genantit_b_pt         ,"genantit_b_pt/D"          );
        outTree_->Branch("genantit_b_eta"           ,&genantit_b_eta         ,"genantit_b_eta/D"          );
        outTree_->Branch("genantit_b_phi"           ,&genantit_b_phi         ,"genantit_b_phi/D"          );
        outTree_->Branch("genantit_b_mass"           ,&genantit_b_mass         ,"genantit_b_mass/D"          );
        outTree_->Branch("gent_w_pt"           ,&gent_w_pt         ,"gent_w_pt/D"          );
        outTree_->Branch("gent_w_eta"           ,&gent_w_eta         ,"gent_w_eta/D"          );
        outTree_->Branch("gent_w_phi"           ,&gent_w_phi         ,"gent_w_phi/D"          );
        outTree_->Branch("gent_w_mass"           ,&gent_w_mass         ,"gent_w_mass/D"          );
        outTree_->Branch("genantit_w_pt"           ,&genantit_w_pt         ,"genantit_w_pt/D"          );
        outTree_->Branch("genantit_w_eta"           ,&genantit_w_eta         ,"genantit_w_eta/D"          );
        outTree_->Branch("genantit_w_phi"           ,&genantit_w_phi         ,"genantit_w_phi/D"          );
        outTree_->Branch("genantit_w_mass"           ,&genantit_w_mass         ,"genantit_w_mass/D"          );
        outTree_->Branch("gent_w_tag"           ,&gent_w_tag         ,"gent_w_tag/D"          );
        outTree_->Branch("gent_w_q1_pt"           ,&gent_w_q1_pt         ,"gent_w_q1_pt/D"          );
        outTree_->Branch("gent_w_q1_eta"           ,&gent_w_q1_eta         ,"gent_w_q1_eta/D"          );
        outTree_->Branch("gent_w_q1_phi"           ,&gent_w_q1_phi         ,"gent_w_q1_phi/D"          );
        outTree_->Branch("gent_w_q1_e"           ,&gent_w_q1_e         ,"gent_w_q1_e/D"          );
        outTree_->Branch("gent_w_q1_pdg"           ,&gent_w_q1_pdg         ,"gent_w_q1_pdg/D"          );
        outTree_->Branch("gent_w_q2_pt"           ,&gent_w_q2_pt         ,"gent_w_q2_pt/D"          );
        outTree_->Branch("gent_w_q2_eta"           ,&gent_w_q2_eta         ,"gent_w_q2_eta/D"          );
        outTree_->Branch("gent_w_q2_phi"           ,&gent_w_q2_phi         ,"gent_w_q2_phi/D"          );
        outTree_->Branch("gent_w_q2_e"           ,&gent_w_q2_e         ,"gent_w_q2_e/D"          );
        outTree_->Branch("gent_w_q2_pdg"           ,&gent_w_q2_pdg         ,"gent_w_q2_pdg/D"          );
        outTree_->Branch("genantit_w_tag"           ,&genantit_w_tag         ,"genantit_w_tag/D"          );
        outTree_->Branch("genantit_w_q1_pt"           ,&genantit_w_q1_pt         ,"genantit_w_q1_pt/D"          );
        outTree_->Branch("genantit_w_q1_eta"           ,&genantit_w_q1_eta         ,"genantit_w_q1_eta/D"          );
        outTree_->Branch("genantit_w_q1_phi"           ,&genantit_w_q1_phi         ,"genantit_w_q1_phi/D"          );
        outTree_->Branch("genantit_w_q1_e"           ,&genantit_w_q1_e         ,"genantit_w_q1_e/D"          );
        outTree_->Branch("genantit_w_q1_pdg"           ,&genantit_w_q1_pdg         ,"genantit_w_q1_pdg/D"          );
        outTree_->Branch("genantit_w_q2_pt"           ,&genantit_w_q2_pt         ,"genantit_w_q2_pt/D"          );
        outTree_->Branch("genantit_w_q2_eta"           ,&genantit_w_q2_eta         ,"genantit_w_q2_eta/D"          );
        outTree_->Branch("genantit_w_q2_phi"           ,&genantit_w_q2_phi         ,"genantit_w_q2_phi/D"          );
        outTree_->Branch("genantit_w_q2_e"           ,&genantit_w_q2_e         ,"gent_w_q2_e/D"          );
        outTree_->Branch("genantit_w_q2_pdg"           ,&genantit_w_q2_pdg         ,"genantit_w_q2_pdg/D"          );

        outTree_->Branch("ptgenq1l"           ,ptgenq1l         ,"ptgenq1l[5]/D"          );
        outTree_->Branch("etagenq1l"           ,etagenq1l         ,"etagenq1l[5]/D"          );
        outTree_->Branch("phigenq1l"           ,phigenq1l       ,"phigenq1l[5]/D"          );
        outTree_->Branch("egenq1l"           ,egenq1l         ,"egenq1l[5]/D"          );
        outTree_->Branch("ptgenq1f"           ,ptgenq1f         ,"ptgenq1f[5]/D"          );
        outTree_->Branch("etagenq1f"           ,etagenq1f         ,"etagenq1f[5]/D"          );
        outTree_->Branch("phigenq1f"           ,phigenq1f       ,"phigenq1f[5]/D"          );
        outTree_->Branch("egenq1f"           ,egenq1f         ,"egenq1f[5]/D"          );
        outTree_->Branch("ptgenq2l"           ,ptgenq2l         ,"ptgenq2l[5]/D"          );
        outTree_->Branch("etagenq2l"           ,etagenq2l         ,"etagenq2l[5]/D"          );
        outTree_->Branch("phigenq2l"           ,phigenq2l       ,"phigenq2l[5]/D"          );
        outTree_->Branch("egenq2l"           ,egenq2l         ,"egenq2l[5]/D"          );
        outTree_->Branch("ptgenq2f"           ,ptgenq2f         ,"ptgenq2f[5]/D"          );
        outTree_->Branch("etagenq2f"           ,etagenq2f         ,"etagenq2f[5]/D"          );
        outTree_->Branch("phigenq2f"           ,phigenq2f       ,"phigenq2f[5]/D"          );
        outTree_->Branch("egenq2f"           ,egenq2f         ,"egenq2f[5]/D"          );
        outTree_->Branch("ptgenq3l"           ,ptgenq3l         ,"ptgenq3l[5]/D"          );
        outTree_->Branch("etagenq3l"           ,etagenq3l         ,"etagenq3l[5]/D"          );
        outTree_->Branch("phigenq3l"           ,phigenq3l       ,"phigenq3l[5]/D"          );
        outTree_->Branch("egenq3l"           ,egenq3l         ,"egenq3l[5]/D"          );
        outTree_->Branch("ptgenq3f"           ,ptgenq3f         ,"ptgenq3f[5]/D"          );
        outTree_->Branch("etagenq3f"           ,etagenq3f         ,"etagenq3f[5]/D"          );
        outTree_->Branch("phigenq3f"           ,phigenq3f       ,"phigenq3f[5]/D"          );
        outTree_->Branch("egenq3f"           ,egenq3f         ,"egenq3f[5]/D"          );
        outTree_->Branch("ptgenq4l"           ,ptgenq4l         ,"ptgenq4l[5]/D"          );
        outTree_->Branch("etagenq4l"           ,etagenq4l         ,"etagenq4l[5]/D"          );
        outTree_->Branch("phigenq4l"           ,phigenq4l       ,"phigenq4l[5]/D"          );
        outTree_->Branch("egenq4l"           ,egenq4l         ,"egenq4l[5]/D"          );
        outTree_->Branch("ptgenq4f"           ,ptgenq4f         ,"ptgenq4f[5]/D"          );
        outTree_->Branch("etagenq4f"           ,etagenq4f         ,"etagenq4f[5]/D"          );
        outTree_->Branch("phigenq4f"           ,phigenq4f       ,"phigenq4f[5]/D"          );
        outTree_->Branch("egenq4f"           ,egenq4f         ,"egenq4f[5]/D"          );
        outTree_->Branch("ptgenq5l"           ,ptgenq5l         ,"ptgenq5l[5]/D"          );
        outTree_->Branch("etagenq5l"           ,etagenq5l         ,"etagenq5l[5]/D"          );
        outTree_->Branch("phigenq5l"           ,phigenq5l       ,"phigenq5l[5]/D"          );
        outTree_->Branch("egenq5l"           ,egenq5l         ,"egenq5l[5]/D"          );
        outTree_->Branch("ptgenq5f"           ,ptgenq5f         ,"ptgenq5f[5]/D"          );
        outTree_->Branch("etagenq5f"           ,etagenq5f         ,"etagenq5f[5]/D"          );
        outTree_->Branch("phigenq5f"           ,phigenq5f       ,"phigenq5f[5]/D"          );
        outTree_->Branch("egenq5f"           ,egenq5f         ,"egenq5f[5]/D"          );
        outTree_->Branch("mothergenq1f"           ,mothergenq1f         ,"mothergenq1f[5]/D"          );
        outTree_->Branch("mothergenq2f"           ,mothergenq2f         ,"mothergenq2f[5]/D"          );
        outTree_->Branch("mothergenq3f"           ,mothergenq3f         ,"mothergenq3f[5]/D"          );
        outTree_->Branch("mothergenq4f"           ,mothergenq4f         ,"mothergenq4f[5]/D"          );
        outTree_->Branch("mothergenq5f"           ,mothergenq5f         ,"mothergenq5f[5]/D"          );
        
        outTree_->Branch("mothergengf"           ,mothergengf         ,"mothergengf[10]/D"          );
        outTree_->Branch("mmothergengf"           ,mmothergengf         ,"mmothergengf[10]/D"          );

        outTree_->Branch("mmothergenq1f"           ,mmothergenq1f         ,"mmothergenq1f[5]/D"          );
        outTree_->Branch("mmothergenq2f"           ,mmothergenq2f         ,"mmothergenq2f[5]/D"          );
        outTree_->Branch("mmothergenq3f"           ,mmothergenq3f         ,"mmothergenq3f[5]/D"          );
        outTree_->Branch("mmothergenq4f"           ,mmothergenq4f         ,"mmothergenq4f[5]/D"          );
        outTree_->Branch("mmothergenq5f"           ,mmothergenq5f         ,"mmothergenq5f[5]/D"          );

    }
    outTree_->Branch("run"             ,&run            ,"run/I");//
    outTree_->Branch("ls"              ,&ls             ,"ls/I"             );//Synch
    outTree_->Branch("nLooseEle"       ,&nLooseEle      ,"nLooseEle/I");//
    outTree_->Branch("nLooseMu"        ,&nLooseMu       ,"nLooseMu/I");//
    outTree_->Branch("event"           ,&nevent         ,"event/I"          );
    outTree_->Branch("nVtx"            ,&nVtx           ,"nVtx/I"           );

    outTree_->Branch("ptVlep"          ,&ptVlep         ,"ptVlep/D"         );

    outTree_->Branch("jetAK8puppi_ptJEC"          ,&jetAK8puppi_ptJEC         ,"jetAK8puppi_ptJEC/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_new"          ,&jetAK8puppi_ptJEC_new         ,"jetAK8puppi_ptJEC_new/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_m"          ,&jetAK8puppi_ptJEC_m         ,"jetAK8puppi_ptJEC_m/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_newnew"          ,&jetAK8puppi_ptJEC_newnew         ,"jetAK8puppi_ptJEC_newnew/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_JEC_up"          ,&jetAK8puppi_ptJEC_JEC_up         ,"jetAK8puppi_ptJEC_JEC_up/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_JEC_down"          ,&jetAK8puppi_ptJEC_JEC_down         ,"jetAK8puppi_ptJEC_JEC_down/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_JER_down"          ,&jetAK8puppi_ptJEC_JER_down         ,"jetAK8puppi_ptJEC_JER_down/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_JER_up"          ,&jetAK8puppi_ptJEC_JER_up         ,"jetAK8puppi_ptJEC_JER_up/D"         );
    outTree_->Branch("jetAK8puppi_eta"          ,&jetAK8puppi_eta         ,"jetAK8puppi_eta/D"         );
    outTree_->Branch("jetAK8puppi_phi"          ,&jetAK8puppi_phi         ,"jetAK8puppi_phi/D"         );
    outTree_->Branch("jetAK8puppi_tau1"          ,&jetAK8puppi_tau1         ,"jetAK8puppi_tau1/D"         );
    outTree_->Branch("jetAK8puppi_tau2"          ,&jetAK8puppi_tau2         ,"jetAK8puppi_tau2/D"         );
    outTree_->Branch("jetAK8puppi_tau3"          ,&jetAK8puppi_tau3         ,"jetAK8puppi_tau3/D"         );
    outTree_->Branch("jetAK8puppi_tau21"          ,&jetAK8puppi_tau21         ,"jetAK8puppi_tau21/D"         );
    outTree_->Branch("jetAK8puppi_tau4"          ,&jetAK8puppi_tau4         ,"jetAK8puppi_tau4/D"         );
    outTree_->Branch("jetAK8puppi_tau42"          ,&jetAK8puppi_tau42         ,"jetAK8puppi_tau42/D"         );
    outTree_->Branch("jetAK8puppi_sd"          ,&jetAK8puppi_sd         ,"jetAK8puppi_sd/D"         );
    // DeepAK8
    outTree_->Branch("jetAK8puppi_dnnTop"         ,&jetAK8puppi_dnnTop       ,"jetAK8puppi_dnnTop/D"         );
    outTree_->Branch("jetAK8puppi_dnnW"           ,&jetAK8puppi_dnnW         ,"jetAK8puppi_dnnW/D"           );
    outTree_->Branch("jetAK8puppi_dnnH4q"         ,&jetAK8puppi_dnnH4q       ,"jetAK8puppi_dnnH4q/D"         );
    outTree_->Branch("jetAK8puppi_dnnTop_2"         ,&jetAK8puppi_dnnTop_2       ,"jetAK8puppi_dnnTop_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnW_2"           ,&jetAK8puppi_dnnW_2         ,"jetAK8puppi_dnnW_2/D"           );
    outTree_->Branch("jetAK8puppi_dnnH4q_2"         ,&jetAK8puppi_dnnH4q_2       ,"jetAK8puppi_dnnH4q_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnTop_3"         ,&jetAK8puppi_dnnTop_3       ,"jetAK8puppi_dnnTop_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnW_3"           ,&jetAK8puppi_dnnW_3         ,"jetAK8puppi_dnnW_3/D"           );
    outTree_->Branch("jetAK8puppi_dnnH4q_3"         ,&jetAK8puppi_dnnH4q_3       ,"jetAK8puppi_dnnH4q_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnZ"         ,&jetAK8puppi_dnnZ       ,"jetAK8puppi_dnnZ/D"         );
    outTree_->Branch("jetAK8puppi_dnnZbb"         ,&jetAK8puppi_dnnZbb       ,"jetAK8puppi_dnnZbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnHbb"         ,&jetAK8puppi_dnnHbb       ,"jetAK8puppi_dnnHbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnZ_2"         ,&jetAK8puppi_dnnZ_2       ,"jetAK8puppi_dnnZ_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnZbb_2"         ,&jetAK8puppi_dnnZbb_2       ,"jetAK8puppi_dnnZbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnHbb_2"         ,&jetAK8puppi_dnnHbb_2       ,"jetAK8puppi_dnnHbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnZ_3"         ,&jetAK8puppi_dnnZ_3       ,"jetAK8puppi_dnnZ_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnZbb_3"         ,&jetAK8puppi_dnnZbb_3       ,"jetAK8puppi_dnnZbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnHbb_3"         ,&jetAK8puppi_dnnHbb_3       ,"jetAK8puppi_dnnHbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnqcd"         ,&jetAK8puppi_dnnqcd       ,"jetAK8puppi_dnnqcd/D"         );
    outTree_->Branch("jetAK8puppi_dnntop"         ,&jetAK8puppi_dnntop       ,"jetAK8puppi_dnntop/D"         );
    outTree_->Branch("jetAK8puppi_dnnw"         ,&jetAK8puppi_dnnw       ,"jetAK8puppi_dnnw/D"         );
    outTree_->Branch("jetAK8puppi_dnnz"         ,&jetAK8puppi_dnnz       ,"jetAK8puppi_dnnz/D"         );
    outTree_->Branch("jetAK8puppi_dnnzbb"         ,&jetAK8puppi_dnnzbb       ,"jetAK8puppi_dnnzbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnhbb"         ,&jetAK8puppi_dnnhbb       ,"jetAK8puppi_dnnhbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnh4q"         ,&jetAK8puppi_dnnh4q       ,"jetAK8puppi_dnnh4q/D"         );
    outTree_->Branch("jetAK8puppi_dnnqcd_2"         ,&jetAK8puppi_dnnqcd_2       ,"jetAK8puppi_dnnqcd_2/D"         );
    outTree_->Branch("jetAK8puppi_dnntop_2"         ,&jetAK8puppi_dnntop_2       ,"jetAK8puppi_dnntop_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnw_2"         ,&jetAK8puppi_dnnw_2       ,"jetAK8puppi_dnnw_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnz_2"         ,&jetAK8puppi_dnnz_2       ,"jetAK8puppi_dnnz_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnzbb_2"         ,&jetAK8puppi_dnnzbb_2       ,"jetAK8puppi_dnnzbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnhbb_2"         ,&jetAK8puppi_dnnhbb_2       ,"jetAK8puppi_dnnhbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnh4q_2"         ,&jetAK8puppi_dnnh4q_2       ,"jetAK8puppi_dnnh4q_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnqcd_3"         ,&jetAK8puppi_dnnqcd_3       ,"jetAK8puppi_dnnqcd_3/D"         );
    outTree_->Branch("jetAK8puppi_dnntop_3"         ,&jetAK8puppi_dnntop_3       ,"jetAK8puppi_dnntop_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnw_3"         ,&jetAK8puppi_dnnw_3       ,"jetAK8puppi_dnnw_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnz_3"         ,&jetAK8puppi_dnnz_3       ,"jetAK8puppi_dnnz_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnzbb_3"         ,&jetAK8puppi_dnnzbb_3       ,"jetAK8puppi_dnnzbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnhbb_3"         ,&jetAK8puppi_dnnhbb_3       ,"jetAK8puppi_dnnhbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnh4q_3"         ,&jetAK8puppi_dnnh4q_3       ,"jetAK8puppi_dnnh4q_3/D"         );

    //Decorrelated DeepAK8
    outTree_->Branch("jetAK8puppi_dnnDecorrTop"         ,&jetAK8puppi_dnnDecorrTop       ,"jetAK8puppi_dnnDecorrTop/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrW"           ,&jetAK8puppi_dnnDecorrW         ,"jetAK8puppi_dnnDecorrW/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrH4q"         ,&jetAK8puppi_dnnDecorrH4q       ,"jetAK8puppi_dnnDecorrH4q/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrTop_2"         ,&jetAK8puppi_dnnDecorrTop_2       ,"jetAK8puppi_dnnDecorrTop_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrW_2"           ,&jetAK8puppi_dnnDecorrW_2         ,"jetAK8puppi_dnnDecorrW_2/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrH4q_2"         ,&jetAK8puppi_dnnDecorrH4q_2       ,"jetAK8puppi_dnnDecorrH4q_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrTop_3"         ,&jetAK8puppi_dnnDecorrTop_3       ,"jetAK8puppi_dnnDecorrTop_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrW_3"           ,&jetAK8puppi_dnnDecorrW_3         ,"jetAK8puppi_dnnDecorrW_3/D"           );
    outTree_->Branch("jetAK8puppi_dnnDecorrH4q_3"         ,&jetAK8puppi_dnnDecorrH4q_3       ,"jetAK8puppi_dnnDecorrH4q_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ"         ,&jetAK8puppi_dnnDecorrZ       ,"jetAK8puppi_dnnDecorrZ/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZbb"         ,&jetAK8puppi_dnnDecorrZbb       ,"jetAK8puppi_dnnDecorrZbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrHbb"         ,&jetAK8puppi_dnnDecorrHbb       ,"jetAK8puppi_dnnDecorrHbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ_2"         ,&jetAK8puppi_dnnDecorrZ_2       ,"jetAK8puppi_dnnDecorrZ_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZbb_2"         ,&jetAK8puppi_dnnDecorrZbb_2       ,"jetAK8puppi_dnnDecorrZbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrHbb_2"         ,&jetAK8puppi_dnnDecorrHbb_2       ,"jetAK8puppi_dnnDecorrHbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZ_3"         ,&jetAK8puppi_dnnDecorrZ_3       ,"jetAK8puppi_dnnDecorrZ_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrZbb_3"         ,&jetAK8puppi_dnnDecorrZbb_3       ,"jetAK8puppi_dnnDecorrZbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrHbb_3"         ,&jetAK8puppi_dnnDecorrHbb_3       ,"jetAK8puppi_dnnDecorrHbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbb"         ,&jetAK8puppi_dnnDecorrbb       ,"jetAK8puppi_dnnDecorrbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrcc"         ,&jetAK8puppi_dnnDecorrcc       ,"jetAK8puppi_dnnDecorrcc/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbbnog"         ,&jetAK8puppi_dnnDecorrbbnog       ,"jetAK8puppi_dnnDecorrbbnog/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrccnog"         ,&jetAK8puppi_dnnDecorrccnog       ,"jetAK8puppi_dnnDecorrccnog/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbb_2"         ,&jetAK8puppi_dnnDecorrbb_2       ,"jetAK8puppi_dnnDecorrbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrcc_2"         ,&jetAK8puppi_dnnDecorrcc_2       ,"jetAK8puppi_dnnDecorrcc_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbbnog_2"         ,&jetAK8puppi_dnnDecorrbbnog_2       ,"jetAK8puppi_dnnDecorrbbnog_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrccnog_2"         ,&jetAK8puppi_dnnDecorrccnog_2       ,"jetAK8puppi_dnnDecorrccnog_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbb_3"         ,&jetAK8puppi_dnnDecorrbb_3       ,"jetAK8puppi_dnnDecorrbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrcc_3"         ,&jetAK8puppi_dnnDecorrcc_3       ,"jetAK8puppi_dnnDecorrcc_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrbbnog_3"         ,&jetAK8puppi_dnnDecorrbbnog_3       ,"jetAK8puppi_dnnDecorrbbnog_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrccnog_3"         ,&jetAK8puppi_dnnDecorrccnog_3       ,"jetAK8puppi_dnnDecorrccnog_3/D"         );

    outTree_->Branch("jetAK8puppi_dnnDecorrqcd"         ,&jetAK8puppi_dnnDecorrqcd       ,"jetAK8puppi_dnnDecorrqcd/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrtop"         ,&jetAK8puppi_dnnDecorrtop       ,"jetAK8puppi_dnnDecorrtop/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrw"         ,&jetAK8puppi_dnnDecorrw       ,"jetAK8puppi_dnnDecorrw/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrz"         ,&jetAK8puppi_dnnDecorrz       ,"jetAK8puppi_dnnDecorrz/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrzbb"         ,&jetAK8puppi_dnnDecorrzbb       ,"jetAK8puppi_dnnDecorrzbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrhbb"         ,&jetAK8puppi_dnnDecorrhbb       ,"jetAK8puppi_dnnDecorrhbb/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrh4q"         ,&jetAK8puppi_dnnDecorrh4q       ,"jetAK8puppi_dnnDecorrh4q/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrqcd_2"         ,&jetAK8puppi_dnnDecorrqcd_2       ,"jetAK8puppi_dnnDecorrqcd_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrtop_2"         ,&jetAK8puppi_dnnDecorrtop_2       ,"jetAK8puppi_dnnDecorrtop_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrw_2"         ,&jetAK8puppi_dnnDecorrw_2       ,"jetAK8puppi_dnnDecorrw_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrz_2"         ,&jetAK8puppi_dnnDecorrz_2       ,"jetAK8puppi_dnnDecorrz_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrzbb_2"         ,&jetAK8puppi_dnnDecorrzbb_2       ,"jetAK8puppi_dnnDecorrzbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrhbb_2"         ,&jetAK8puppi_dnnDecorrhbb_2       ,"jetAK8puppi_dnnDecorrhbb_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrh4q_2"         ,&jetAK8puppi_dnnDecorrh4q_2       ,"jetAK8puppi_dnnDecorrh4q_2/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrqcd_3"         ,&jetAK8puppi_dnnDecorrqcd_3       ,"jetAK8puppi_dnnDecorrqcd_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrtop_3"         ,&jetAK8puppi_dnnDecorrtop_3       ,"jetAK8puppi_dnnDecorrtop_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrw_3"         ,&jetAK8puppi_dnnDecorrw_3       ,"jetAK8puppi_dnnDecorrw_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrz_3"         ,&jetAK8puppi_dnnDecorrz_3       ,"jetAK8puppi_dnnDecorrz_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrzbb_3"         ,&jetAK8puppi_dnnDecorrzbb_3       ,"jetAK8puppi_dnnDecorrzbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrhbb_3"         ,&jetAK8puppi_dnnDecorrhbb_3       ,"jetAK8puppi_dnnDecorrhbb_3/D"         );
    outTree_->Branch("jetAK8puppi_dnnDecorrh4q_3"         ,&jetAK8puppi_dnnDecorrh4q_3       ,"jetAK8puppi_dnnDecorrh4q_3/D"         );


    
    outTree_->Branch("jetAK8puppi_ptJEC_2"          ,&jetAK8puppi_ptJEC_2         ,"jetAK8puppi_ptJEC_2/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_2_new"          ,&jetAK8puppi_ptJEC_2_new         ,"jetAK8puppi_ptJEC_2_new/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_2_JEC_up"          ,&jetAK8puppi_ptJEC_2_JEC_up         ,"jetAK8puppi_ptJEC_2_JEC_up/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_2_JEC_down"          ,&jetAK8puppi_ptJEC_2_JEC_down         ,"jetAK8puppi_ptJEC_2_JEC_down/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_2_JER_down"          ,&jetAK8puppi_ptJEC_2_JER_down         ,"jetAK8puppi_ptJEC_2_JER_down/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_2_JER_up"          ,&jetAK8puppi_ptJEC_2_JER_up         ,"jetAK8puppi_ptJEC_2_JER_up/D"         );
    outTree_->Branch("jetAK8puppi_eta_2"          ,&jetAK8puppi_eta_2         ,"jetAK8puppi_eta_2/D"         );
    outTree_->Branch("jetAK8puppi_phi_2"          ,&jetAK8puppi_phi_2         ,"jetAK8puppi_phi_2/D"         );
    outTree_->Branch("jetAK8puppi_tau1_2"          ,&jetAK8puppi_tau1_2         ,"jetAK8puppi_tau1_2/D"         );
    outTree_->Branch("jetAK8puppi_tau2_2"          ,&jetAK8puppi_tau2_2         ,"jetAK8puppi_tau2_2/D"         );
    outTree_->Branch("jetAK8puppi_tau3_2"          ,&jetAK8puppi_tau3_2         ,"jetAK8puppi_tau3_2/D"         );
    outTree_->Branch("jetAK8puppi_tau21_2"          ,&jetAK8puppi_tau21_2         ,"jetAK8puppi_tau21_2/D"         );
    outTree_->Branch("jetAK8puppi_tau4_2"          ,&jetAK8puppi_tau4_2         ,"jetAK8puppi_tau4_2/D"         );
    outTree_->Branch("jetAK8puppi_tau42_2"          ,&jetAK8puppi_tau42_2         ,"jetAK8puppi_tau42_2/D"         );
    outTree_->Branch("jetAK8puppi_sd_2"          ,&jetAK8puppi_sd_2         ,"jetAK8puppi_sd_2/D"         );


    outTree_->Branch("jetAK8puppi_ptJEC_3"          ,&jetAK8puppi_ptJEC_3         ,"jetAK8puppi_ptJEC_3/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_3_new"          ,&jetAK8puppi_ptJEC_3_new         ,"jetAK8puppi_ptJEC_3_new/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_3_JEC_up"          ,&jetAK8puppi_ptJEC_3_JEC_up         ,"jetAK8puppi_ptJEC_3_JEC_up/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_3_JEC_down"          ,&jetAK8puppi_ptJEC_3_JEC_down         ,"jetAK8puppi_ptJEC_3_JEC_down/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_3_JER_down"          ,&jetAK8puppi_ptJEC_3_JER_down         ,"jetAK8puppi_ptJEC_3_JER_down/D"         );
    outTree_->Branch("jetAK8puppi_ptJEC_3_JER_up"          ,&jetAK8puppi_ptJEC_3_JER_up         ,"jetAK8puppi_ptJEC_3_JER_up/D"         );
    
    outTree_->Branch("jetAK8puppi_e"          ,&jetAK8puppi_e         ,"jetAK8puppi_e/D"         );
    outTree_->Branch("jetAK8puppi_e_new"          ,&jetAK8puppi_e_new         ,"jetAK8puppi_e_new/D"         );
    outTree_->Branch("jetAK8puppi_e_JEC_up"          ,&jetAK8puppi_e_JEC_up         ,"jetAK8puppi_e_JEC_up/D"         );
    outTree_->Branch("jetAK8puppi_e_JEC_down"          ,&jetAK8puppi_e_JEC_down         ,"jetAK8puppi_e_JEC_down/D"         );
    outTree_->Branch("jetAK8puppi_e_JER_down"          ,&jetAK8puppi_e_JER_down         ,"jetAK8puppi_e_JER_down/D"         );
    outTree_->Branch("jetAK8puppi_e_JER_up"          ,&jetAK8puppi_e_JER_up         ,"jetAK8puppi_e_JER_up/D"         );
    outTree_->Branch("jetAK8puppi_e_2"          ,&jetAK8puppi_e_2         ,"jetAK8puppi_e_2/D"         );
    outTree_->Branch("jetAK8puppi_e_2_new"          ,&jetAK8puppi_e_2_new         ,"jetAK8puppi_e_2_new/D"         );
    outTree_->Branch("jetAK8puppi_e_2_JEC_up"          ,&jetAK8puppi_e_2_JEC_up         ,"jetAK8puppi_e_2_JEC_up/D"         );
    outTree_->Branch("jetAK8puppi_e_2_JEC_down"          ,&jetAK8puppi_e_2_JEC_down         ,"jetAK8puppi_e_2_JEC_down/D"         );
    outTree_->Branch("jetAK8puppi_e_2_JER_down"          ,&jetAK8puppi_e_2_JER_down         ,"jetAK8puppi_e_2_JER_down/D"         );
    outTree_->Branch("jetAK8puppi_e_2_JER_up"          ,&jetAK8puppi_e_2_JER_up         ,"jetAK8puppi_e_2_JER_up/D"         );
    outTree_->Branch("jetAK8puppi_e_3"          ,&jetAK8puppi_e_3         ,"jetAK8puppi_e_3/D"         );
    outTree_->Branch("jetAK8puppi_e_3_new"          ,&jetAK8puppi_e_3_new         ,"jetAK8puppi_e_3_new/D"         );
    outTree_->Branch("jetAK8puppi_e_3_JEC_up"          ,&jetAK8puppi_e_3_JEC_up         ,"jetAK8puppi_e_3_JEC_up/D"         );
    outTree_->Branch("jetAK8puppi_e_3_JEC_down"          ,&jetAK8puppi_e_3_JEC_down         ,"jetAK8puppi_e_3_JEC_down/D"         );
    outTree_->Branch("jetAK8puppi_e_3_JER_down"          ,&jetAK8puppi_e_3_JER_down         ,"jetAK8puppi_e_3_JER_down/D"         );
    outTree_->Branch("jetAK8puppi_e_3_JER_up"          ,&jetAK8puppi_e_3_JER_up         ,"jetAK8puppi_e_3_JER_up/D"         );

    outTree_->Branch("jetAK8puppi_eta_3"          ,&jetAK8puppi_eta_3         ,"jetAK8puppi_eta_3/D"         );
    outTree_->Branch("jetAK8puppi_phi_3"          ,&jetAK8puppi_phi_3         ,"jetAK8puppi_phi_3/D"         );
    outTree_->Branch("jetAK8puppi_tau1_3"          ,&jetAK8puppi_tau1_3         ,"jetAK8puppi_tau1_3/D"         );
    outTree_->Branch("jetAK8puppi_tau2_3"          ,&jetAK8puppi_tau2_3         ,"jetAK8puppi_tau2_3/D"         );
    outTree_->Branch("jetAK8puppi_tau3_3"          ,&jetAK8puppi_tau3_3         ,"jetAK8puppi_tau3_3/D"         );
    outTree_->Branch("jetAK8puppi_tau21_3"          ,&jetAK8puppi_tau21_3         ,"jetAK8puppi_tau21_3/D"         );
    outTree_->Branch("jetAK8puppi_tau4_3"          ,&jetAK8puppi_tau4_3         ,"jetAK8puppi_tau4_3/D"         );
    outTree_->Branch("jetAK8puppi_tau42_3"          ,&jetAK8puppi_tau42_3         ,"jetAK8puppi_tau42_3/D"         );
    outTree_->Branch("jetAK8puppi_sd_3"          ,&jetAK8puppi_sd_3         ,"jetAK8puppi_sd_3/D"         );


    outTree_->Branch("vbfeta"    ,&vbfeta   ,"vbfeta/D"   );
    outTree_->Branch("vbfmjj"    ,&vbfmjj   ,"vbfmjj/D"   );

    outTree_->Branch("vbftag"    ,&vbftag   ,"vbftag/I"   );
    outTree_->Branch("nj1"    ,&nj1   ,"nj1/I"   );
    outTree_->Branch("nj2"    ,&nj2   ,"nj2/I"   );
    //outTree_->Branch("ak8sj11"    ,&ak8sj11 ) ;
    //outTree_->Branch("ak8sj12"    ,&ak8sj12 ) ;
    //outTree_->Branch("ak8sj13"    ,&ak8sj13 ) ;
    //outTree_->Branch("ak8sj14"    ,&ak8sj14 ) ;
    //outTree_->Branch("ak8sj15"    ,&ak8sj15  );
    //outTree_->Branch("ak8sj21"    ,&ak8sj21 ) ;
    //outTree_->Branch("ak8sj22"    ,&ak8sj22 ) ;
    //outTree_->Branch("ak8sj23"    ,&ak8sj23 ) ;
    //outTree_->Branch("ak8sj24"    ,&ak8sj24 ) ;
    //outTree_->Branch("ak8sj25"    ,&ak8sj25  );
    //outTree_->Branch("puppi_softdropj1"    ,&puppi_softdropj1 ) ;
    //outTree_->Branch("puppi_softdropj2"    ,&puppi_softdropj2 ) ;

    outTree_->Branch("yVlep"           ,&yVlep          ,"yVlep/D"          );
    outTree_->Branch("phiVlep"         ,&phiVlep        ,"phiVlep/D"        );
    outTree_->Branch("massVlep"        ,&massVlep       ,"massVlep/D"       );
    outTree_->Branch("mtVlep"          ,&mtVlep         ,"mtVlep/D"         );
    outTree_->Branch("lep"             ,&lep            ,"lep/I"            );
    outTree_->Branch("channel"         ,&channel        ,"channel/I"        );


 
    /// Generic kinematic quantities
    outTree_->Branch("ptlep1"          ,&ptlep1         ,"ptlep1/D"         );
    outTree_->Branch("ptlep2"          ,&ptlep2         ,"ptlep2/D"         );
    outTree_->Branch("etalep1"         ,&etalep1        ,"etalep1/D"        );
    outTree_->Branch("etalep2"         ,&etalep2        ,"etalep2/D"        );
    outTree_->Branch("philep1"         ,&philep1        ,"philep1/D"        );
    outTree_->Branch("philep2"         ,&philep2        ,"philep2/D"        );
    outTree_->Branch("met"             ,&met            ,"met/D"            );
    outTree_->Branch("metPhi"          ,&metPhi         ,"metPhi/D"         );

    /// Other quantities
    outTree_->Branch("theWeight", &theWeight, "theWeight/D");
    outTreew_->Branch("theWeight", &theWeight, "theWeight/D");
    outTree_->Branch("nump", &nump, "nump/D");
    outTree_->Branch("numm", &numm, "numm/D");
    outTree_->Branch("npT"           ,&npT         ,"npT/D"          );
    outTree_->Branch("npIT"           ,&npIT         ,"npIT/D"          );
    outTree_->Branch("nBX"           ,&nBX         ,"nBX/I"          );

    outTree_->Branch("delPhilepmet"    ,&delPhilepmet   ,"delPhilepmet/D"   );
    outTree_->Branch("deltaRlepjet"    ,&deltaRlepjet   ,"deltaRlepjet/D"   );
    outTree_->Branch("delPhijetmet"    ,&delPhijetmet   ,"delPhijetmet/D"   );
    outTree_->Branch("delPhijetlep"    ,&delPhijetlep   ,"delPhijetlep/D"   );
  
    outTree_->Branch("deltaRlepjet_2"    ,&deltaRlepjet_2   ,"deltaRlepjet_2/D"   );
    outTree_->Branch("delPhijetmet_2"    ,&delPhijetmet_2   ,"delPhijetmet_2/D"   );
    outTree_->Branch("delPhijetlep_2"    ,&delPhijetlep_2   ,"delPhijetlep_2/D"   );
 
    outTree_->Branch("IDLoose", &IDLoose, "IDLoose/O");
    outTree_->Branch("IDTight", &IDTight, "IDTight/O");
    outTree_->Branch("IDLoose_2", &IDLoose_2, "IDLoose_2/O");
    outTree_->Branch("IDTight_2", &IDTight_2, "IDTight_2/O");
    outTree_->Branch("IDLoose_3", &IDLoose_3, "IDLoose_3/O");
    outTree_->Branch("IDTight_3", &IDTight_3, "IDTight_3/O");
    outTree_->Branch("isHighPt",&isHighPt, "isHighPt/O");
    outTree_->Branch("isHEEP",&isHEEP, "isHEEP/O");
    outTree_->Branch("trackIso",&trackIso,"trackIso/D");
    outTree_->Branch("muchaiso",&muchaiso,"muchaiso/D");
    outTree_->Branch("muneuiso",&muneuiso,"muneuiso/D");
    outTree_->Branch("muphoiso",&muphoiso,"muphoiso/D");
    outTree_->Branch("muPU",&muPU,"muPU/D");
    outTree_->Branch("muisolation",&muisolation,"muisolation/D");
    //after JEC varible
    outTree_->Branch("METraw_et",&METraw_et,"METraw_et/D");
    outTree_->Branch("METraw_phi",&METraw_phi,"METraw_phi/D");
    outTree_->Branch("METraw_sumEt",&METraw_sumEt,"METraw_sumEt/D");
    outTree_->Branch("MET_et",&MET_et,"MET_et/D");
    outTree_->Branch("MET_phi",&MET_phi,"MET_phi/D");
    outTree_->Branch("MET_sumEt",&MET_sumEt,"MET_sumEt/D");
    
    outTree_->Branch("MET_et_new", &MET_et_new, "MET_et_new/D");
    // Marked for debug
    outTree_->Branch("MET_et_JEC_up", &MET_et_JEC_up, "MET_et_JEC_up/D");
    outTree_->Branch("MET_et_JEC_down", &MET_et_JEC_down, "MET_et_JEC_down/D");
    outTree_->Branch("MET_et_JER_up", &MET_et_JER_up, "MET_et_JER_up/D");
    outTree_->Branch("MET_et_JER_down", &MET_et_JER_down, "MET_et_JER_down/D");
    // Marked for debug
    outTree_->Branch("MET_et_m", &MET_et_m, "MET_et_m/D");
    outTree_->Branch("MET_et_old", &MET_et_old, "MET_et_old/D");
    outTree_->Branch("MET_phi_m", &MET_phi_m, "MET_phi_m/D");
    // Marked for debug
    outTree_->Branch("MET_phi_new", &MET_phi_new, "MET_phi_new/D");
    outTree_->Branch("MET_phi_JEC_up", &MET_phi_JEC_up, "MET_phi_JEC_up/D");
    outTree_->Branch("MET_phi_JEC_down", &MET_phi_JEC_down, "MET_phi_JEC_down/D");
    outTree_->Branch("MET_phi_JER_up", &MET_phi_JER_up, "MET_phi_JER_up/D");
    outTree_->Branch("MET_phi_JER_down", &MET_phi_JER_down, "MET_phi_JER_down/D");
    // Marked for debug

    outTree_->Branch("jetAK8puppi_pt1",&jetAK8puppi_pt1,"jetAK8puppi_pt1[4]/D");
    outTree_->Branch("jetAK8puppi_eta1",&jetAK8puppi_eta1,"jetAK8puppi_eta1[4]/D");
    outTree_->Branch("jetAK8puppi_mass1",&jetAK8puppi_mass1,"jetAK8puppi_mass1[4]/D");
    outTree_->Branch("candMasspuppiJEC",&candMasspuppiJEC,"candMasspuppiJEC/D");
    outTree_->Branch("m_jlv",&m_jlv,"m_jlv/D");
    outTree_->Branch("candMasspuppiJEC_new",&candMasspuppiJEC_new,"candMasspuppiJEC_new/D");
    outTree_->Branch("m_jlv_new",&m_jlv_new,"m_jlv_new/D");
    
    outTree_->Branch("candMasspuppiJEC_JEC_up",&candMasspuppiJEC_JEC_up,"candMasspuppiJEC_JEC_up/D");
    outTree_->Branch("m_jlv_JEC_up",&m_jlv_JEC_up,"m_jlv_JEC_up/D");
    outTree_->Branch("candMasspuppiJEC_JEC_down",&candMasspuppiJEC_JEC_down,"candMasspuppiJEC_JEC_down/D");
    outTree_->Branch("m_jlv_JEC_down",&m_jlv_JEC_down,"m_jlv_JEC_down/D");
    outTree_->Branch("candMasspuppiJEC_JER_down",&candMasspuppiJEC_JER_down,"candMasspuppiJEC_JER_down/D");
    outTree_->Branch("m_jlv_JER_down",&m_jlv_JER_down,"m_jlv_JER_down/D");
    outTree_->Branch("candMasspuppiJEC_JER_up",&candMasspuppiJEC_JER_up,"candMasspuppiJEC_JER_up/D");
    outTree_->Branch("m_jlv_JER_up",&m_jlv_JER_up,"m_jlv_JER_up/D");

    outTree_->Branch("massww",&massww,"massww[3]/D");
    outTree_->Branch("masslvj1",&masslvj1,"masslvj1/D");
    outTree_->Branch("masslvj2",&masslvj2,"masslvj2/D");
    outTree_->Branch("massj1j2",&massj1j2,"massj1j2/D");

    outTree_->Branch("ptVlepJEC",&ptVlepJEC,"ptVlepJEC/D");
    outTree_->Branch("yVlepJEC",&yVlepJEC,"yVlepJEC/D");
    outTree_->Branch("phiVlepJEC",&phiVlepJEC,"phiVlepJEC/D");
    outTree_->Branch("massVlepJEC",&massVlepJEC,"massVlepJEC/D");
    outTree_->Branch("mtVlepJEC",&mtVlepJEC,"mtVlepJEC/D");
    outTree_->Branch("ptVlepJEC_new",&ptVlepJEC_new,"ptVlepJEC_new/D");
    outTree_->Branch("yVlepJEC_new",&yVlepJEC_new,"yVlepJEC_new/D");
    outTree_->Branch("phiVlepJEC_new",&phiVlepJEC_new,"phiVlepJEC_new/D");
    outTree_->Branch("massVlepJEC_new",&massVlepJEC_new,"massVlepJEC_new/D");
    outTree_->Branch("mtVlepJEC_new",&mtVlepJEC_new,"mtVlepJEC_new/D");
    
    outTree_->Branch("ptVlepJEC_JEC_up",&ptVlepJEC_JEC_up,"ptVlepJEC_JEC_up/D");
    outTree_->Branch("yVlepJEC_JEC_up",&yVlepJEC_JEC_up,"yVlepJEC_JEC_up/D");
    outTree_->Branch("phiVlepJEC_JEC_up",&phiVlepJEC_JEC_up,"phiVlepJEC_JEC_up/D");
    outTree_->Branch("massVlepJEC_JEC_up",&massVlepJEC_JEC_up,"massVlepJEC_JEC_up/D");
    outTree_->Branch("mtVlepJEC_JEC_up",&mtVlepJEC_JEC_up,"mtVlepJEC_JEC_up/D");
    outTree_->Branch("ptVlepJEC_JEC_down",&ptVlepJEC_JEC_down,"ptVlepJEC_JEC_down/D");
    outTree_->Branch("yVlepJEC_JEC_down",&yVlepJEC_JEC_down,"yVlepJEC_JEC_down/D");
    outTree_->Branch("phiVlepJEC_JEC_down",&phiVlepJEC_JEC_down,"phiVlepJEC_JEC_down/D");
    outTree_->Branch("massVlepJEC_JEC_down",&massVlepJEC_JEC_down,"massVlepJEC_JEC_down/D");
    outTree_->Branch("mtVlepJEC_JEC_down",&mtVlepJEC_JEC_down,"mtVlepJEC_JEC_down/D");
    
    outTree_->Branch("ptVlepJEC_JER_up",&ptVlepJEC_JER_up,"ptVlepJEC_JER_up/D");
    outTree_->Branch("yVlepJEC_JER_up",&yVlepJEC_JER_up,"yVlepJEC_JER_up/D");
    outTree_->Branch("phiVlepJEC_JER_up",&phiVlepJEC_JER_up,"phiVlepJEC_JER_up/D");
    outTree_->Branch("massVlepJEC_JER_up",&massVlepJEC_JER_up,"massVlepJEC_JER_up/D");
    outTree_->Branch("mtVlepJEC_JER_up",&mtVlepJEC_JER_up,"mtVlepJEC_JER_up/D");
    outTree_->Branch("ptVlepJEC_JER_down",&ptVlepJEC_JER_down,"ptVlepJEC_JER_down/D");
    outTree_->Branch("yVlepJEC_JER_down",&yVlepJEC_JER_down,"yVlepJEC_JER_down/D");
    outTree_->Branch("phiVlepJEC_JER_down",&phiVlepJEC_JER_down,"phiVlepJEC_JER_down/D");
    outTree_->Branch("massVlepJEC_JER_down",&massVlepJEC_JER_down,"massVlepJEC_JER_down/D");
    outTree_->Branch("mtVlepJEC_JER_down",&mtVlepJEC_JER_down,"mtVlepJEC_JER_down/D");

    ///HLT bits
    outTree_->Branch("HLT_Ele1"  ,&HLT_Ele1 ,"HLT_Ele1/I" );
    outTree_->Branch("HLT_Ele2"  ,&HLT_Ele2 ,"HLT_Ele2/I" );
    outTree_->Branch("HLT_Ele3"  ,&HLT_Ele3 ,"HLT_Ele3/I" );
    outTree_->Branch("HLT_Ele4"  ,&HLT_Ele4 ,"HLT_Ele4/I" );
    outTree_->Branch("HLT_Ele5"  ,&HLT_Ele5 ,"HLT_Ele5/I" );
    outTree_->Branch("HLT_Ele6"  ,&HLT_Ele6 ,"HLT_Ele6/I" );
    outTree_->Branch("HLT_Ele7"  ,&HLT_Ele7 ,"HLT_Ele7/I" );
    outTree_->Branch("HLT_Ele8"  ,&HLT_Ele8 ,"HLT_Ele8/I" );
    outTree_->Branch("HLT_Mu1"   ,&HLT_Mu1  ,"HLT_Mu1/I"  );
    outTree_->Branch("HLT_Mu2"   ,&HLT_Mu2  ,"HLT_Mu2/I"  );
    outTree_->Branch("HLT_Mu3"   ,&HLT_Mu3  ,"HLT_Mu3/I"  );
    outTree_->Branch("HLT_Mu4"   ,&HLT_Mu4  ,"HLT_Mu4/I"  );
    outTree_->Branch("HLT_Mu5"   ,&HLT_Mu5  ,"HLT_Mu5/I"  );
    outTree_->Branch("HLT_Mu6"   ,&HLT_Mu6  ,"HLT_Mu6/I"  );
    outTree_->Branch("HLT_Mu7"   ,&HLT_Mu7  ,"HLT_Mu7/I"  );
    outTree_->Branch("HLT_Mu8"   ,&HLT_Mu8  ,"HLT_Mu8/I"  );
    outTree_->Branch("HLT_Mu9"   ,&HLT_Mu9  ,"HLT_Mu9/I"  );
    outTree_->Branch("HLT_Mu10"   ,&HLT_Mu10  ,"HLT_Mu10/I"  );
    outTree_->Branch("HLT_Mu11"   ,&HLT_Mu11  ,"HLT_Mu11/I"  );
    outTree_->Branch("HLT_Mu12"   ,&HLT_Mu12  ,"HLT_Mu12/I"  );
    // filter
    outTree_->Branch("passFilter_HBHE"                 ,&passFilter_HBHE_                ,"passFilter_HBHE_/O");
    outTree_->Branch("passFilter_HBHEIso"                 ,&passFilter_HBHEIso_                ,"passFilter_HBHEIso_/O");
    outTree_->Branch("passFilter_GlobalHalo"              ,&passFilter_GlobalHalo_             ,"passFilter_GlobalHalo_/O");
    outTree_->Branch("passFilter_ECALDeadCell"         ,&passFilter_ECALDeadCell_        ,"passFilter_ECALDeadCell_/O");
    outTree_->Branch("passFilter_GoodVtx"              ,&passFilter_GoodVtx_             ,"passFilter_GoodVtx_/O");
    outTree_->Branch("passFilter_EEBadSc"              ,&passFilter_EEBadSc_             ,"passFilter_EEBadSc_/O");
    outTree_->Branch("passFilter_badMuon"                 ,&passFilter_badMuon_                ,"passFilter_badMuon_/O");
    outTree_->Branch("passFilter_badChargedHadron"                 ,&passFilter_badChargedHadron_                ,"passFilter_badChargedHadron_/O");
    outTree_->Branch("passecalBadCalibFilterUpdate"                 ,&passecalBadCalibFilterUpdate_                ,"passecalBadCalibFilterUpdate_/O");

    /// AK4 Jets Info
    outTree_->Branch("ak4jet_hf"        , ak4jet_hf       ,"ak4jet_hf[8]/I"       );
    outTree_->Branch("ak4jet_pf"        , ak4jet_pf       ,"ak4jet_pf[8]/I"       );
    outTree_->Branch("ak4jet_pt"        , ak4jet_pt       ,"ak4jet_pt[8]/D"       );
    outTree_->Branch("ak4jet_pt_uncorr"        , ak4jet_pt_uncorr       ,"ak4jet_pt_uncorr[8]/D"       );
    outTree_->Branch("ak4jet_eta"        , ak4jet_eta       ,"ak4jet_eta[8]/D"       );
    outTree_->Branch("ak4jet_phi"        , ak4jet_phi       ,"ak4jet_phi[8]/D"       );
    outTree_->Branch("ak4jet_e"        , ak4jet_e       ,"ak4jet_e[8]/D"       );
    outTree_->Branch("ak4jet_dr"        , ak4jet_dr       ,"ak4jet_dr[8]/D"       );
    outTree_->Branch("ak4jet_csv"        , ak4jet_csv       ,"ak4jet_csv[8]/D"       );
    outTree_->Branch("ak4jet_icsv"        , ak4jet_icsv       ,"ak4jet_icsv[8]/D"       );
    outTree_->Branch("ak4jet_deepcsvudsg"        , ak4jet_deepcsvudsg       ,"ak4jet_deepcsvudsg[8]/D"       );
    outTree_->Branch("ak4jet_deepcsvb"        , ak4jet_deepcsvb       ,"ak4jet_deepcsvb[8]/D"       );
    outTree_->Branch("ak4jet_deepcsvc"        , ak4jet_deepcsvc       ,"ak4jet_deepcsvc[8]/D"       );
    outTree_->Branch("ak4jet_deepcsvbb"        , ak4jet_deepcsvbb       ,"ak4jet_deepcsvbb[8]/D"       );
    outTree_->Branch("ak4jet_deepcsvcc"        , ak4jet_deepcsvcc       ,"ak4jet_deepcsvcc[8]/D"       );
    outTree_->Branch("ak4jet_IDLoose"        , ak4jet_IDLoose       ,"ak4jet_IDLoose[8]/D"       );
    outTree_->Branch("ak4jet_IDTight"        , ak4jet_IDTight       ,"ak4jet_IDTight[8]/D"       );
    
    /// Gen Level quantities
    if(RunOnMC_){
    outTree_->Branch("pttau"        , pttau       ,"pttau[4]/D"       );
    outTree_->Branch("etatau"        , etatau       ,"etatau[4]/D"       );
    outTree_->Branch("phitau"        , phitau       ,"phitau[4]/D"       );
    outTree_->Branch("etau"        , etau       ,"etau[4]/D"       );
    outTree_->Branch("pdgidtau"        , pdgidtau       ,"pdgidtau[4]/D"       );
    
    outTree_->Branch("pttau_2"        , pttau_2       ,"pttau_2[4]/D"       );
    outTree_->Branch("etatau_2"        , etatau_2       ,"etatau_2[4]/D"       );
    outTree_->Branch("phitau_2"        , phitau_2       ,"phitau_2[4]/D"       );
    outTree_->Branch("etau_2"        , etau_2       ,"etau_2[4]/D"       );
    outTree_->Branch("pdgidtau_2"        , pdgidtau_2       ,"pdgidtau_2[4]/D"       );
    
    outTree_->Branch("pttau_3"        , pttau_3       ,"pttau_3[4]/D"       );
    outTree_->Branch("etatau_3"        , etatau_3       ,"etatau_3[4]/D"       );
    outTree_->Branch("phitau_3"        , phitau_3       ,"phitau_3[4]/D"       );
    outTree_->Branch("etau_3"        , etau_3       ,"etau_3[4]/D"       );
    outTree_->Branch("pdgidtau_3"        , pdgidtau_3       ,"pdgidtau_3[4]/D"       );
    
    outTree_->Branch("ptq"        , ptq       ,"ptq[3]/D"       );
    outTree_->Branch("etaq"        , etaq       ,"etaq[3]/D"       );
    outTree_->Branch("phiq"        , phiq       ,"phiq[3]/D"       );
    outTree_->Branch("eq"        , eq       ,"eq[3]/D"       );
    outTree_->Branch("pdgidq"        , pdgidq       ,"pdgidq[3]/D"       );
    
    outTree_->Branch("ptq_2"        , ptq_2       ,"ptq_2[3]/D"       );
    outTree_->Branch("etaq_2"        , etaq_2       ,"etaq_2[3]/D"       );
    outTree_->Branch("phiq_2"        , phiq_2       ,"phiq_2[3]/D"       );
    outTree_->Branch("eq_2"        , eq_2       ,"eq_2[3]/D"       );
    outTree_->Branch("pdgidq_2"        , pdgidq_2       ,"pdgidq_2[3]/D"       );
    
    outTree_->Branch("ptq_3"        , ptq_3       ,"ptq_3[3]/D"       );
    outTree_->Branch("etaq_3"        , etaq_3       ,"etaq_3[3]/D"       );
    outTree_->Branch("phiq_3"        , phiq_3       ,"phiq_3[3]/D"       );
    outTree_->Branch("eq_3"        , eq_3       ,"eq_3[4]/D"       );
    outTree_->Branch("pdgidq_3"        , pdgidq_3       ,"pdgidq_3[3]/D"       );
    

    outTree_->Branch("gen_ele_pt"        ,&gen_ele_pt       ,"gen_ele_pt/D"       );
    outTree_->Branch("gen_ele_eta"        ,&gen_ele_eta       ,"gen_ele_eta/D"       );
    outTree_->Branch("gen_ele_phi"        ,&gen_ele_phi       ,"gen_ele_phi/D"       );
    outTree_->Branch("gen_ele_e"        ,&gen_ele_e       ,"gen_ele_e/D"       );
    outTree_->Branch("gen_mu_pt"        ,&gen_mu_pt       ,"gen_mu_pt/D"       );
    outTree_->Branch("gen_mu_eta"        ,&gen_mu_eta       ,"gen_mu_eta/D"       );
    outTree_->Branch("gen_mu_phi"        ,&gen_mu_phi       ,"gen_mu_phi/D"       );
    outTree_->Branch("gen_mu_e"        ,&gen_mu_e       ,"gen_mu_e/D"       );
    outTree_->Branch("gen_ele_pt_2"        ,&gen_ele_pt_2       ,"gen_ele_pt_2/D"       );
    outTree_->Branch("gen_ele_eta_2"        ,&gen_ele_eta_2       ,"gen_ele_eta_2/D"       );
    outTree_->Branch("gen_ele_phi_2"        ,&gen_ele_phi_2       ,"gen_ele_phi_2/D"       );
    outTree_->Branch("gen_ele_e_2"        ,&gen_ele_e_2       ,"gen_ele_e_2/D"       );
    outTree_->Branch("gen_mu_pt_2"        ,&gen_mu_pt_2       ,"gen_mu_pt_2/D"       );
    outTree_->Branch("gen_mu_eta_2"        ,&gen_mu_eta_2       ,"gen_mu_eta_2/D"       );
    outTree_->Branch("gen_mu_phi_2"        ,&gen_mu_phi_2       ,"gen_mu_phi_2/D"       );
    outTree_->Branch("gen_mu_e_2"        ,&gen_mu_e_2       ,"gen_mu_e_2/D"       );
    outTree_->Branch("gen_ele_pt_3"        ,&gen_ele_pt_3       ,"gen_ele_pt_3/D"       );
    outTree_->Branch("gen_ele_eta_3"        ,&gen_ele_eta_3       ,"gen_ele_eta_3/D"       );
    outTree_->Branch("gen_ele_phi_3"        ,&gen_ele_phi_3       ,"gen_ele_phi_3/D"       );
    outTree_->Branch("gen_ele_e_3"        ,&gen_ele_e_3       ,"gen_ele_e_3/D"       );
    outTree_->Branch("gen_mu_pt_3"        ,&gen_mu_pt_3       ,"gen_mu_pt_3/D"       );
    outTree_->Branch("gen_mu_eta_3"        ,&gen_mu_eta_3       ,"gen_mu_eta_3/D"       );
    outTree_->Branch("gen_mu_phi_3"        ,&gen_mu_phi_3       ,"gen_mu_phi_3/D"       );
    outTree_->Branch("gen_mu_e_3"        ,&gen_mu_e_3       ,"gen_mu_e_3/D"       );







    
    outTree_->Branch("genmatch_ele_pt"        ,&genmatch_ele_pt       ,"genmatch_ele_pt/D"       );
    outTree_->Branch("genmatch_ele_eta"        ,&genmatch_ele_eta       ,"genmatch_ele_eta/D"       );
    outTree_->Branch("genmatch_ele_phi"        ,&genmatch_ele_phi       ,"genmatch_ele_phi/D"       );
    outTree_->Branch("genmatch_ele_e"        ,&genmatch_ele_e       ,"genmatch_ele_e/D"       );
    outTree_->Branch("genmatch_ele_dr"        ,&genmatch_ele_dr       ,"genmatch_ele_dr/D"       );
    outTree_->Branch("genmatch_mu_pt"        ,&genmatch_mu_pt       ,"genmatch_mu_pt/D"       );
    outTree_->Branch("genmatch_mu_eta"        ,&genmatch_mu_eta       ,"genmatch_mu_eta/D"       );
    outTree_->Branch("genmatch_mu_phi"        ,&genmatch_mu_phi       ,"genmatch_mu_phi/D"       );
    outTree_->Branch("genmatch_mu_e"        ,&genmatch_mu_e       ,"genmatch_mu_e/D"       );
    outTree_->Branch("genmatch_mu_dr"        ,&genmatch_mu_dr       ,"genmatch_mu_dr/D"       );
    outTree_->Branch("gentop_pt"        ,&gentop_pt       ,"gentop_pt/D"       );
    outTree_->Branch("gentop_eta"        ,&gentop_eta       ,"gentop_eta/D"       );
    outTree_->Branch("gentop_phi"        ,&gentop_phi       ,"gentop_phi/D"       );
    outTree_->Branch("gentop_mass"        ,&gentop_mass       ,"gentop_mass/D"       );
    outTree_->Branch("genantitop_pt"        ,&genantitop_pt       ,"genantitop_pt/D"       );
    outTree_->Branch("genantitop_eta"        ,&genantitop_eta       ,"genantitop_eta/D"       );
    outTree_->Branch("genantitop_phi"        ,&genantitop_phi       ,"genantitop_phi/D"       );
    outTree_->Branch("genantitop_mass"        ,&genantitop_mass       ,"genantitop_mass/D"       );





    }
    //outTree_->Branch("");
}


VVVTreeMaker::~VVVTreeMaker()
{
 
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}

void VVVTreeMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
	/*
     edm::Handle<LHERunInfoProduct> runw;
        typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
        
        iRun.getByToken(LhestrToken_,runw);
        LHERunInfoProduct myLHERunInfoProduct = *(runw.product());
        
        for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
            std::cout << iter->tag() << std::endl;
            std::vector<std::string> lines = iter->lines();
            for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
                std::cout << lines.at(iLine);
            }
        }*/
    std::cout << "VVVTreeMaker endJob()... endRun" << std::endl;
}


void
VVVTreeMaker::endJob() {
    std::cout << "VVVTreeMaker endJob()..." << std::endl;
}

void
VVVTreeMaker::EDBR_Debug() {

    std::cout << "========================================================" << std::endl;
    std::cout << "================                   =====================" << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "================                   =====================" << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "========================================================" << std::endl;

    std::cout << "start :  EDBRTreeMaker MET variable..." << std::endl;
    std::cout << "to check MetHandle : ptVlepJEC_new ...: " << ptVlepJEC_new << std::endl;
    std::cout << "end :  EDBRTreeMaker MET variable..." << std::endl;

    std::cout << "start :  EDBRTreeMaker Gen variable ..." << std::endl;
    std::cout << "to check Gen : massgenwf ...: " << massgenwf[0] << std::endl;
    std::cout << "to check Gen : ptgenwl ...: " << ptgenwl[0] << std::endl;
    std::cout << "to check Gen : taggenwl ...: " << taggenwl[0] << std::endl;
    std::cout << "to check Gen : genw_q1_pt ...: " << genw_q1_pt[0] << std::endl;
    std::cout << "to check Gen : genw_q2_eta ...: " << genw_q2_eta[0] << std::endl;
    std::cout << "end :  EDBRTreeMaker Gen variable ..." << std::endl;

    std::cout << "start :  EDBRTreeMaker filter variable ..." << std::endl;
    std::cout << "to check hlt : HLT_Ele2 ...: " << HLT_Ele2 << std::endl;
    std::cout << "end :  EDBRTreeMaker filter variable ..." << std::endl;

    std::cout << "start :  EDBRTreeMaker jet variable ..." << std::endl;
    std::cout << "to check jet : jetAK8puppi_ptJEC ...: " << jetAK8puppi_ptJEC << std::endl;
    std::cout << "to check jet : jetAK8puppi_ptJEC_new ...: " << jetAK8puppi_ptJEC_new << std::endl;
    std::cout << "to check jet : jetAK8puppi_ptJEC_newnew ...: " << jetAK8puppi_ptJEC_newnew << std::endl;
    std::cout << "to check jet : jetAK8puppi_tau1 ...: " << jetAK8puppi_tau1 << std::endl;
    std::cout << "to check jet : jetAK8puppi_sd ...: " << jetAK8puppi_sd << std::endl;
    std::cout << "to check jet : jetAK8puppi_ptJEC_JEC_down ...: " << jetAK8puppi_ptJEC_JEC_down << std::endl;
    std::cout << "to check jet : jetAK8puppi_tau3 ...: " << jetAK8puppi_tau3 << std::endl;
    std::cout << "to check jet : jetAK8puppi_dnnW ...: " << jetAK8puppi_dnnW << std::endl;
    std::cout << "to check jet : jetAK8puppi_ptJEC_2 ...: " << jetAK8puppi_ptJEC_2 << std::endl;
    std::cout << "to check jet : jetAK8puppi_ptJEC_2_JEC_down ...: " << jetAK8puppi_ptJEC_2_JEC_down << std::endl;
    std::cout << "to check jet : jetAK8puppi_sd_2 ...: " << jetAK8puppi_sd_2 << std::endl;
    std::cout << "end :  EDBRTreeMaker jet variable ..." << std::endl;

    std::cout << "start :  EDBRTreeMaker met variable ..." << std::endl;
    std::cout << "to check met : MET_et ...: " << MET_et << std::endl;
    std::cout << "end :  EDBRTreeMaker met variable ..." << std::endl;


    std::cout << "start :  EDBRTreeMaker weight variable ..." << std::endl;
    std::cout << "to check  : L1prefiring ...: " << L1prefiring << std::endl;
    std::cout << "end :  EDBRTreeMaker weight variable ..." << std::endl;

    std::cout << "start :  EDBRTreeMaker lepton variable ..." << std::endl;
    std::cout << "to check  : nLooseEle ...: " << nLooseEle << std::endl;
    std::cout << "to check  : nVtx ...: " << nVtx << std::endl;
    std::cout << "to check  : nLooseMu ...: " << nLooseMu << std::endl;
    std::cout << "to check  : ptVlep ...: " << ptVlep << std::endl;
    std::cout << "to check  : nLooseEle ...: " << nLooseEle << std::endl;
    std::cout << "end :  EDBRTreeMaker lepton variable ..." << std::endl;


    std::cout << "start :  EDBRTreeMaker Generic kinematic quantities ..." << std::endl;
    std::cout << "to check  : ptlep1 ...: " << ptlep1 << std::endl;
    std::cout << "to check  : yVlep ...: " << ptlep1 << std::endl;
    std::cout << "to check  : etalep2 ...: " << etalep2 << std::endl;
    std::cout << "to check  : philep2 ...: " << philep2 << std::endl;
    std::cout << "to check  : lep ...: " << ptlep1 << std::endl;
    std::cout << "to check  : metPhi ...: " << ptlep1 << std::endl;
    std::cout << "end :  EDBRTreeMaker Generic kinematic quantities ..." << std::endl;

    std::cout << "start :  EDBRTreeMaker Other quantities ..." << std::endl;
    std::cout << "to check  : isHighPt ...: " << isHighPt << std::endl;
    std::cout << "to check  : muisolation ...: " << muisolation << std::endl;
    std::cout << "to check  : METraw_et ...: " << METraw_et << std::endl;
    std::cout << "to check  : MET_phi ...: " << MET_phi << std::endl;
    std::cout << "to check  : MET_et_JEC_up ...: " << MET_et_JEC_up << std::endl;
    std::cout << "to check  : MET_phi_m ...: " << MET_phi_m << std::endl;
    std::cout << "to check  : candMasspuppiJEC ...: " << candMasspuppiJEC << std::endl;
    std::cout << "to check  : candMasspuppiJEC_new ...: " << candMasspuppiJEC_new << std::endl;
    std::cout << "to check  : m_jlv_JEC_down ...: " << m_jlv_JEC_down << std::endl;
    std::cout << "to check  : massww ...: " << massww << std::endl;
    std::cout << "to check  : massj1j2 ...: " << massj1j2 << std::endl;
    std::cout << "to check  : ptVlepJEC ...: " << ptVlepJEC << std::endl;
    std::cout << "to check  : massVlepJEC ...: " << massVlepJEC << std::endl;
    std::cout << "to check  : ptVlepJEC_new ...: " << ptVlepJEC_new << std::endl;
    std::cout << "to check  : massVlepJEC_JEC_up ...: " << massVlepJEC_JEC_up << std::endl;
    std::cout << "to check  : mtVlepJEC_JEC_up ...: " << mtVlepJEC_JEC_up << std::endl;
    std::cout << "to check  : ak4jet_hf ...: " << ak4jet_hf[0] << std::endl;
    std::cout << "to check  : ak4jet_pf ...: " << ak4jet_pf[0] << std::endl;
    std::cout << "to check  : ak4jet_pt ...: " << ak4jet_pt[0] << std::endl;
    std::cout << "to check  : ak4jet_phi ...: " << ak4jet_phi[0] << std::endl;
    std::cout << "to check  : ak4jet_deepcsvbb ...: " << ak4jet_deepcsvbb[0] << std::endl;
    std::cout << "to check  : ak4jet_IDLoose ...: " << ak4jet_IDLoose[0] << std::endl;
    std::cout << "to check  : ak4jet_csv ...: " << ak4jet_csv[0] << std::endl;
    std::cout << "to check  : ak4jet_deepcsvcc ...: " << ak4jet_deepcsvcc[0] << std::endl;
    std::cout << "end :  EDBRTreeMaker Other quantities ..." << std::endl;


    std::cout << "========================================================" << std::endl;
    std::cout << "================                   =====================" << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "================                   =====================" << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "========================================================" << std::endl;


}


#endif

