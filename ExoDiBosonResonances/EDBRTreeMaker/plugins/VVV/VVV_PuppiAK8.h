#ifndef _EDBRPuppiAK8_
#define _EDBRPuppiAK8_

void VVVTreeMaker::PuppiAK8_Store( edm::Event const & iEvent ){

    edm::Handle<double> rho;
    iEvent.getByToken(rhoToken_      , rho     );
    double fastJetRho = *(rho.product());
    useless = fastJetRho;

            for(size_t ij=0; ij<puppijets_->size()&&ij<4;ij++){
                corr_AK8puppi[ij] = 1;
                const pat::Jet& hadronicVa = puppijets_->at(ij);
                reco::Candidate::LorentzVector uncorrJet;
                if(not isJEC_) doCorrOnTheFly_ = false;
                if( doCorrOnTheFly_ ){
                    uncorrJet = hadronicVa.correctedP4(0);
                    jecAK8puppi_->setJetEta( uncorrJet.eta()          );
                    jecAK8puppi_->setJetPt ( uncorrJet.pt()           );
                    jecAK8puppi_->setJetE  ( uncorrJet.energy()       );
                    jecAK8puppi_->setRho   (fastJetRho);
                    jecAK8puppi_->setNPV   (nVtx);
                    jecAK8puppi_->setJetA  (hadronicVa.jetArea());
                    corr_AK8puppi[ij] = jecAK8puppi_->getCorrection();
                }
                else{uncorrJet = hadronicVa.p4();}

                if(ij<4){
                    jetAK8puppi_pt1[ij] = corr_AK8puppi[ij]*uncorrJet.pt();
                    jetAK8puppi_mass1[ij] = corr_AK8puppi[ij]*uncorrJet.mass();
                    jetAK8puppi_eta1[ij] = uncorrJet.eta();
                    jetAK8puppi_jec1[ij] = corr_AK8puppi[ij];
                }
		jetAK8puppi_pt1_m[ij]=hadronicVa.p4().pt();
		//cout<<ij<<"   "<<jetAK8puppi_pt1[ij]<<" PF "<<(*puppijets_)[ij].isPFJet()<<endl;
        	if (RunOnMC_){ 
                jetAK8puppi_pt1_newnew[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JEC_central");
                jetAK8puppi_pt1_new[ij]=(*puppijets_)[ij].userFloat("SmearedPt");
                jetAK8puppi_pt1_JEC_up[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JEC_up");
                jetAK8puppi_pt1_JEC_down[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JEC_down");
                jetAK8puppi_pt1_JER_up[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JER_up");
                jetAK8puppi_pt1_JER_down[ij]=(*puppijets_)[ij].userFloat("SmearedPt_JER_down");
                jetAK8puppi_e1_new[ij]=(*puppijets_)[ij].userFloat("SmearedE");
                jetAK8puppi_e1_JEC_up[ij]=(*puppijets_)[ij].userFloat("SmearedE_JEC_up");
                jetAK8puppi_e1_JEC_down[ij]=(*puppijets_)[ij].userFloat("SmearedE_JEC_down");
                jetAK8puppi_e1_JER_up[ij]=(*puppijets_)[ij].userFloat("SmearedE_JER_up");
                jetAK8puppi_e1_JER_down[ij]=(*puppijets_)[ij].userFloat("SmearedE_JER_down");
                }

            }




            int usenumber3 = -1; double pt_larger=0;
            int numvhad = puppijets_->size();
            for( int inum = 0; inum< numvhad; inum++){
                const pat::Jet& Vpuppi = puppijets_->at(inum);
                if(tightJetIDpuppi(Vpuppi)<1) continue;
                if(jetAK8puppi_pt1[inum] > pt_larger && fabs(jetAK8puppi_eta1[inum])<2.4 && inum<4) {pt_larger = jetAK8puppi_pt1[inum]; usenumber3 = inum; continue;}
            }
            //cout<<"usenumber3"<<usenumber3<<endl;
            if (usenumber3>-1) {//2
                const pat::Jet& hadronicVpuppi = puppijets_->at(usenumber3);
                // DeepAK8
                jetAK8puppi_dnnTop       = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnW         = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnH4q       = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnZ         = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnZbb       = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnHbb       = hadronicVpuppi.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnqcd       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDbb")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDcc")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDb")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDc")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnntop       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probTbcq")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnw         = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probWcq")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnz         = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probZcc")+hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnzbb       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnhbb       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnh4q       = hadronicVpuppi.bDiscriminator("pfDeepBoostedJetTags:probHqqqq"); 
                // Decorrelated DeepAK8
                jetAK8puppi_dnnDecorrTop       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnDecorrW         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnDecorrH4q       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnDecorrZ         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnDecorrZbb       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnDecorrHbb       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnDecorrqcd       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnnDecorrbb        = (hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb"))/(hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd); 
                jetAK8puppi_dnnDecorrcc        = (hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc"))/(hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd); 
                jetAK8puppi_dnnDecorrbbnog     = (hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"))/(hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd); 
                jetAK8puppi_dnnDecorrccnog     = (hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc"))/(hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd); 
                jetAK8puppi_dnnDecorrtop       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnDecorrw         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnDecorrz         = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnDecorrzbb       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnDecorrhbb       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnDecorrh4q       = hadronicVpuppi.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq"); 

                 // ------
                jetAK8puppi_ptJEC       = jetAK8puppi_pt1[usenumber3]; // unpruned corrected jet pt
                jetAK8puppi_ptJEC_m       = jetAK8puppi_pt1_m[usenumber3];
                cout << "jetAK8puppi_ptJEC: " << jetAK8puppi_ptJEC << endl;
                cout << "jetAK8puppi_ptJEC_m: " << jetAK8puppi_ptJEC_m << endl;
        	if (RunOnMC_){ 
                jetAK8puppi_ptJEC_new       = jetAK8puppi_pt1_new[usenumber3];
                jetAK8puppi_ptJEC_newnew       = jetAK8puppi_pt1_newnew[usenumber3];
                jetAK8puppi_ptJEC_JEC_up       = jetAK8puppi_pt1_JEC_up[usenumber3];
                jetAK8puppi_ptJEC_JEC_down       = jetAK8puppi_pt1_JEC_down[usenumber3];
                jetAK8puppi_ptJEC_JER_up       = jetAK8puppi_pt1_JER_up[usenumber3];
                jetAK8puppi_ptJEC_JER_down       = jetAK8puppi_pt1_JER_down[usenumber3];
                jetAK8puppi_e_new       = jetAK8puppi_e1_new[usenumber3];
                jetAK8puppi_e_JEC_up       = jetAK8puppi_e1_JEC_up[usenumber3];
                jetAK8puppi_e_JEC_down       = jetAK8puppi_e1_JEC_down[usenumber3];
                jetAK8puppi_e_JER_up       = jetAK8puppi_e1_JER_up[usenumber3];
                jetAK8puppi_e_JER_down       = jetAK8puppi_e1_JER_down[usenumber3];
		}
                jetAK8puppi_eta     = jetAK8puppi_eta1[usenumber3]; // unpruned (w/o jec) jet eta
                jetAK8puppi_phi      = hadronicVpuppi.phi(); // unpruned (w/o jec) jet phi
                jetAK8puppi_tau1         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau1");
                jetAK8puppi_tau2         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau2");
                jetAK8puppi_tau3         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau21        = jetAK8puppi_tau2/jetAK8puppi_tau1;
                jetAK8puppi_tau4         = hadronicVpuppi.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau42        = jetAK8puppi_tau4/jetAK8puppi_tau2;
                jetAK8puppi_sd       =  hadronicVpuppi.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass

                IDLoose = tightJetIDpuppi(hadronicVpuppi);
                IDTight = tightJetIDpuppi(hadronicVpuppi);
            }
            
            int usenumber2 = -1; double pt_larger2=0;
            for( int inum = 0; inum< numvhad; inum++){
                const pat::Jet& Vpuppi = puppijets_->at(inum);
                if(tightJetIDpuppi(Vpuppi)<1) continue;
                if(jetAK8puppi_pt1[inum] > pt_larger2 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum<4) {pt_larger2 = jetAK8puppi_pt1[inum]; usenumber2 = inum; continue;}
            }
            //cout<<"usenumber2"<<usenumber2<<endl;
            if(usenumber2>-1)  {
                const pat::Jet& hadronicVpuppi_2 = puppijets_->at(usenumber2);
                // DeepAK8
                jetAK8puppi_dnnTop_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnW_2         = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnH4q_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnZ_2         = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnZbb_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnHbb_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnqcd_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDbb")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDcc")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDb")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDc")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnntop_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probTbcq")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnw_2         = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probWcq")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnz_2         = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probZcc")+hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnzbb_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnhbb_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnh4q_2       = hadronicVpuppi_2.bDiscriminator("pfDeepBoostedJetTags:probHqqqq"); 
                // Decorrelated DeepAK8
                jetAK8puppi_dnnDecorrTop_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnDecorrW_2         = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnDecorrH4q_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnDecorrZ_2         = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnDecorrZbb_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnDecorrHbb_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnDecorrqcd_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnnDecorrbb_2        = (hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb"))/(hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_2); 
                jetAK8puppi_dnnDecorrcc_2        = (hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc"))/(hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_2); 
                jetAK8puppi_dnnDecorrbbnog_2     = (hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"))/(hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_2); 
                jetAK8puppi_dnnDecorrccnog_2     = (hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc"))/(hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_2); 
                jetAK8puppi_dnnDecorrtop_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnDecorrw_2         = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnDecorrz_2         = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnDecorrzbb_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnDecorrhbb_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnDecorrh4q_2       = hadronicVpuppi_2.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq"); 

                 // ------
                jetAK8puppi_ptJEC_2       = jetAK8puppi_pt1[usenumber2]; // unpruned corrected jet pt
        	if (RunOnMC_){ 
                jetAK8puppi_ptJEC_2_new       = jetAK8puppi_pt1_new[usenumber2];
                jetAK8puppi_ptJEC_2_JEC_up       = jetAK8puppi_pt1_JEC_up[usenumber2];
                jetAK8puppi_ptJEC_2_JEC_down       = jetAK8puppi_pt1_JEC_down[usenumber2];
                jetAK8puppi_ptJEC_2_JER_up       = jetAK8puppi_pt1_JER_up[usenumber2];
                jetAK8puppi_ptJEC_2_JER_down       = jetAK8puppi_pt1_JER_down[usenumber2];
                jetAK8puppi_e_2_new       = jetAK8puppi_e1_new[usenumber2];
                jetAK8puppi_e_2_JEC_up       = jetAK8puppi_e1_JEC_up[usenumber2];
                jetAK8puppi_e_2_JEC_down       = jetAK8puppi_e1_JEC_down[usenumber2];
                jetAK8puppi_e_2_JER_up       = jetAK8puppi_e1_JER_up[usenumber2];
                jetAK8puppi_e_2_JER_down       = jetAK8puppi_e1_JER_down[usenumber2];
		}
                jetAK8puppi_eta_2     = jetAK8puppi_eta1[usenumber2]; // unpruned (w/o jec) jet eta
                jetAK8puppi_phi_2      = hadronicVpuppi_2.phi(); // unpruned (w/o jec) jet phi
                jetAK8puppi_tau1_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau1");
                jetAK8puppi_tau2_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau2");
                jetAK8puppi_tau3_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau21_2        = jetAK8puppi_tau2_2/jetAK8puppi_tau1_2;
                jetAK8puppi_tau4_2         = hadronicVpuppi_2.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau42_2        = jetAK8puppi_tau4_2/jetAK8puppi_tau2_2;
                jetAK8puppi_sd_2       =  hadronicVpuppi_2.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass

                
                IDLoose_2 = tightJetIDpuppi(hadronicVpuppi_2);
                IDTight_2 = tightJetIDpuppi(hadronicVpuppi_2);
            }

            int usenumber1 = -1; double pt_larger1=0;
            for( int inum = 0; inum< numvhad; inum++){
                const pat::Jet& Vpuppi = puppijets_->at(inum);
                if(tightJetIDpuppi(Vpuppi)<1) continue;
                if(jetAK8puppi_pt1[inum] > pt_larger1 && fabs(jetAK8puppi_eta1[inum])<2.4 && inum != usenumber3 && inum != usenumber2 && inum<4) {pt_larger1 = jetAK8puppi_pt1[inum]; usenumber1 = inum; continue;}
            }
            //cout<<"usenumber1"<<usenumber1<<endl;
            if(usenumber1>-1)  {
                const pat::Jet& hadronicVpuppi_3 = puppijets_->at(usenumber1);
                // DeepAK8
                jetAK8puppi_dnnTop_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnW_3         = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnH4q_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnZ_3         = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnZbb_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnHbb_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnqcd_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDbb")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDcc")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDb")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDc")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnntop_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probTbcq")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnw_3         = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probWcq")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnz_3         = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probZcc")+hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnzbb_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnhbb_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnh4q_3       = hadronicVpuppi_3.bDiscriminator("pfDeepBoostedJetTags:probHqqqq"); 
                // Decorrelated DeepAK8
                jetAK8puppi_dnnDecorrTop_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"); 
                jetAK8puppi_dnnDecorrW_3         = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"); 
                jetAK8puppi_dnnDecorrH4q_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD"); 
                jetAK8puppi_dnnDecorrZ_3         = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD"); 
                jetAK8puppi_dnnDecorrZbb_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"); 
                jetAK8puppi_dnnDecorrHbb_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD"); 
                jetAK8puppi_dnnDecorrqcd_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDothers"); 
                jetAK8puppi_dnnDecorrbb_3        = (hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDbb"))/(hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_3); 
                jetAK8puppi_dnnDecorrcc_3        = (hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probQCDcc"))/(hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_3); 
                jetAK8puppi_dnnDecorrbbnog_3     = (hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"))/(hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_3); 
                jetAK8puppi_dnnDecorrccnog_3     = (hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc"))/(hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+jetAK8puppi_dnnDecorrqcd_3); 
                jetAK8puppi_dnnDecorrtop_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbcq")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probTbqq"); 
                jetAK8puppi_dnnDecorrw_3         = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWcq")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probWqq"); 
                jetAK8puppi_dnnDecorrz_3         = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZcc")+hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZqq"); 
                jetAK8puppi_dnnDecorrzbb_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probZbb"); 
                jetAK8puppi_dnnDecorrhbb_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHbb"); 
                jetAK8puppi_dnnDecorrh4q_3       = hadronicVpuppi_3.bDiscriminator("pfMassDecorrelatedDeepBoostedJetTags:probHqqqq"); 

                 // ------
                jetAK8puppi_ptJEC_3       = jetAK8puppi_pt1[usenumber1]; // unpruned corrected jet pt
        	if (RunOnMC_){ 
                jetAK8puppi_ptJEC_3_new       = jetAK8puppi_pt1_new[usenumber1];
                jetAK8puppi_ptJEC_3_JEC_up       = jetAK8puppi_pt1_JEC_up[usenumber1];
                jetAK8puppi_ptJEC_3_JEC_down       = jetAK8puppi_pt1_JEC_down[usenumber1];
                jetAK8puppi_ptJEC_3_JER_up       = jetAK8puppi_pt1_JER_up[usenumber1];
                jetAK8puppi_ptJEC_3_JER_down       = jetAK8puppi_pt1_JER_down[usenumber1];
                jetAK8puppi_e_3_new       = jetAK8puppi_e1_new[usenumber1];
                jetAK8puppi_e_3_JEC_up       = jetAK8puppi_e1_JEC_up[usenumber1];
                jetAK8puppi_e_3_JEC_down       = jetAK8puppi_e1_JEC_down[usenumber1];
                jetAK8puppi_e_3_JER_up       = jetAK8puppi_e1_JER_up[usenumber1];
                jetAK8puppi_e_3_JER_down       = jetAK8puppi_e1_JER_down[usenumber1];
		}
                jetAK8puppi_eta_3     = jetAK8puppi_eta1[usenumber1]; // unpruned (w/o jec) jet eta
                jetAK8puppi_phi_3      = hadronicVpuppi_3.phi(); // unpruned (w/o jec) jet phi
                jetAK8puppi_tau1_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau1");
                jetAK8puppi_tau2_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau2");
                jetAK8puppi_tau3_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau21_3        = jetAK8puppi_tau2_3/jetAK8puppi_tau1_3;
                jetAK8puppi_tau4_3         = hadronicVpuppi_3.userFloat("NjettinessAK8Puppi:tau3");
                jetAK8puppi_tau42_3        = jetAK8puppi_tau4_3/jetAK8puppi_tau2_3;
                jetAK8puppi_sd_3       =  hadronicVpuppi_3.userFloat("ak8PFJetsPuppiSoftDropMass"); // uncorrected pruned mass


                IDLoose_3 = tightJetIDpuppi(hadronicVpuppi_3);
                IDTight_3 = tightJetIDpuppi(hadronicVpuppi_3);
            }

}
#endif
