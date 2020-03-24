#ifndef _inv_massInfo_
#define _inv_massInfo_

void
VVVTreeMaker::VVV_inv_mass_Store( void ){
            
            deltaRlepjet = deltaR(etalep1,philep1,jetAK8puppi_eta,jetAK8puppi_phi);
            deltaRlepjet_2 = deltaR(etalep1,philep1,jetAK8puppi_eta_2,jetAK8puppi_phi_2);
            TLorentzVector ghadronicVpuppi, gravitonpuppiJEC,ghadronicVpuppi_2, gravitonpuppiJEC_2;
            ghadronicVpuppi.SetPtEtaPhiM(jetAK8puppi_ptJEC, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2.SetPtEtaPhiM(jetAK8puppi_ptJEC_2, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC = gleptonicV + ghadronicVpuppi+ ghadronicVpuppi_2;
            candMasspuppiJEC     = gravitonpuppiJEC.Mag();
            m_jlv     = (gleptonicV + ghadronicVpuppi).Mag();

            TLorentzVector lvw[3];
            if (RunOnMC_){ 
            TLorentzVector ghadronicVpuppi_new, gravitonpuppiJEC_new,ghadronicVpuppi_2_new;
            ghadronicVpuppi_new.SetPtEtaPhiM(jetAK8puppi_ptJEC_new, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_new.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_new, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_new = gleptonicV_new + ghadronicVpuppi_new+ ghadronicVpuppi_2_new;
            candMasspuppiJEC_new     = gravitonpuppiJEC_new.Mag();
            m_jlv_new     = (gleptonicV_new + ghadronicVpuppi_new).Mag();
	    //cout<<WLeptonic.pt()<<" old "<<WLeptonic.eta()<<"   "<<WLeptonic.phi()<<"   "<<WLeptonic.mass()<<endl;
	    //cout<<WLeptonic_new.pt()<<" new "<<WLeptonic_new.eta()<<"   "<<WLeptonic_new.phi()<<"   "<<WLeptonic_new.mass()<<endl;
            //cout<<jetAK8puppi_sd<<"  jetAK8puppi_sd  "<<ghadronicVpuppi.Mag()<<"   "<<ghadronicVpuppi_new.Mag()<<"   "<<ghadronicVpuppi.E()<<"   "<<ghadronicVpuppi_new.E()<<endl;    
            TLorentzVector ghadronicVpuppi_JEC_up, gravitonpuppiJEC_JEC_up,ghadronicVpuppi_2_JEC_up;
            ghadronicVpuppi_JEC_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_JEC_up, jetAK8puppi_eta, jetAK8puppi_phi,jetAK8puppi_sd);
            ghadronicVpuppi_2_JEC_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JEC_up, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_JEC_up = gleptonicV_JEC_up + ghadronicVpuppi_JEC_up+ ghadronicVpuppi_2_JEC_up;
            candMasspuppiJEC_JEC_up     = gravitonpuppiJEC_JEC_up.Mag();
            m_jlv_JEC_up     = (gleptonicV_JEC_up + ghadronicVpuppi_JEC_up).Mag();
	    //cout<<ghadronicVpuppi_2_JEC_up.Pt()<<"  "<<candMasspuppiJEC_JEC_up<<endl;
            //cout<<jetAK8puppi_ptJEC_JEC_up<<"  "<<gleptonicV_JEC_up.Pt()<<"  "<<m_jlv_JEC_up<<endl;
            
            TLorentzVector ghadronicVpuppi_JEC_down, gravitonpuppiJEC_JEC_down,ghadronicVpuppi_2_JEC_down;
            ghadronicVpuppi_JEC_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_JEC_down, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_JEC_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JEC_down, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_JEC_down = gleptonicV_JEC_down + ghadronicVpuppi_JEC_down+ ghadronicVpuppi_2_JEC_down;
            candMasspuppiJEC_JEC_down     = gravitonpuppiJEC_JEC_down.Mag();
            m_jlv_JEC_down     = (gleptonicV_JEC_down + ghadronicVpuppi_JEC_down).Mag();
            
            TLorentzVector ghadronicVpuppi_JER_up, gravitonpuppiJEC_JER_up,ghadronicVpuppi_2_JER_up;
            ghadronicVpuppi_JER_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_JER_up, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_JER_up.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JER_up, jetAK8puppi_eta_2, jetAK8puppi_phi_2, jetAK8puppi_sd_2);
            gravitonpuppiJEC_JER_up = gleptonicV_JER_up + ghadronicVpuppi_JER_up+ ghadronicVpuppi_2_JER_up;
            candMasspuppiJEC_JER_up     = gravitonpuppiJEC_JER_up.Mag();
            m_jlv_JER_up     = (gleptonicV_JER_up + ghadronicVpuppi_JER_up).Mag();
            
            TLorentzVector ghadronicVpuppi_JER_down, gravitonpuppiJEC_JER_down,ghadronicVpuppi_2_JER_down;
            ghadronicVpuppi_JER_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_JER_down, jetAK8puppi_eta, jetAK8puppi_phi, jetAK8puppi_sd);
            ghadronicVpuppi_2_JER_down.SetPtEtaPhiM(jetAK8puppi_ptJEC_2_JER_down, jetAK8puppi_eta_2, jetAK8puppi_phi_2,jetAK8puppi_sd_2);
            gravitonpuppiJEC_JER_down = gleptonicV_JER_down + ghadronicVpuppi_JER_down+ ghadronicVpuppi_2_JER_down;
            candMasspuppiJEC_JER_down     = gravitonpuppiJEC_JER_down.Mag();
            m_jlv_JER_down     = (gleptonicV_JER_down + ghadronicVpuppi_JER_down).Mag();

	    //swap var and var_new
            double jetAK8puppi_ptJEC_tmp = jetAK8puppi_ptJEC ; jetAK8puppi_ptJEC = jetAK8puppi_ptJEC_new ; jetAK8puppi_ptJEC_new = jetAK8puppi_ptJEC_tmp;
            double jetAK8puppi_ptJEC_2_tmp = jetAK8puppi_ptJEC_2 ; jetAK8puppi_ptJEC_2 = jetAK8puppi_ptJEC_2_new ; jetAK8puppi_ptJEC_2_new = jetAK8puppi_ptJEC_2_tmp;
            double jetAK8puppi_ptJEC_3_tmp = jetAK8puppi_ptJEC_3 ; jetAK8puppi_ptJEC_3 = jetAK8puppi_ptJEC_3_new ; jetAK8puppi_ptJEC_3_new = jetAK8puppi_ptJEC_3_tmp;
            double jetAK8puppi_e_tmp = jetAK8puppi_e ; jetAK8puppi_e = jetAK8puppi_e_new ; jetAK8puppi_e_new = jetAK8puppi_e_tmp;
            double jetAK8puppi_e_2_tmp = jetAK8puppi_e_2 ; jetAK8puppi_e_2 = jetAK8puppi_e_2_new ; jetAK8puppi_e_2_new = jetAK8puppi_e_2_tmp;
            double jetAK8puppi_e_3_tmp = jetAK8puppi_e_3 ; jetAK8puppi_e_3 = jetAK8puppi_e_3_new ; jetAK8puppi_e_3_new = jetAK8puppi_e_3_tmp;
            double ptVlepJEC_tmp = ptVlepJEC ; ptVlepJEC = ptVlepJEC_new ; ptVlepJEC_new = ptVlepJEC_tmp;
            double yVlepJEC_tmp = yVlepJEC ; yVlepJEC = yVlepJEC_new ; yVlepJEC_new = yVlepJEC_tmp;
            double phiVlepJEC_tmp = phiVlepJEC ; phiVlepJEC = phiVlepJEC_new ; phiVlepJEC_new = phiVlepJEC_tmp;
            double massVlepJEC_tmp = massVlepJEC ; massVlepJEC = massVlepJEC_new ; massVlepJEC_new = massVlepJEC_tmp;
            double mtVlepJEC_tmp = mtVlepJEC ; mtVlepJEC = mtVlepJEC_new ; mtVlepJEC_new = mtVlepJEC_tmp;
            double MET_et_tmp = MET_et ; MET_et = MET_et_new ; MET_et_new = MET_et_tmp;
            double MET_phi_tmp = MET_phi ; MET_phi = MET_phi_new ; MET_phi_new = MET_phi_tmp;
            double m_jlv_tmp = m_jlv ; m_jlv = m_jlv_new ; m_jlv_new = m_jlv_tmp;
            double candMasspuppiJEC_tmp = candMasspuppiJEC ; candMasspuppiJEC = candMasspuppiJEC_new ; candMasspuppiJEC_new = candMasspuppiJEC_tmp;

            lvw[0] = gleptonicV_new;
            lvw[1] = ghadronicVpuppi_new;
            lvw[2] = ghadronicVpuppi_2_new;
	    }
            delPhijetmet = deltaPhi(jetAK8puppi_phi, MET_phi);
            delPhijetlep = deltaPhi(jetAK8puppi_phi, phiVlepJEC);
            delPhijetmet_2 = deltaPhi(jetAK8puppi_phi_2, MET_phi);
            delPhijetlep_2 = deltaPhi(jetAK8puppi_phi_2, phiVlepJEC);
        
            delPhilepmet = deltaPhi(philep1, MET_phi);
            mtVlepJEC       =   sqrt(2*ptlep1*MET_et*(1.0-cos(philep1-MET_phi))); //WLeptonic.mt();

            if (!RunOnMC_){ 
            lvw[0] = gleptonicV;
            lvw[1] = ghadronicVpuppi;
            lvw[2] = ghadronicVpuppi_2;
            }
            Double_t Wpt[3];
            Wpt[0]=ptVlepJEC;
            Wpt[1]=jetAK8puppi_ptJEC;
            Wpt[2]=jetAK8puppi_ptJEC_2;
            Int_t *indexx=new Int_t[3];
            TMath::Sort(3,Wpt,indexx,1);
            //cout<<Wpt[indexx[0]]<<"   "<<Wpt[indexx[1]]<<"   "<<Wpt[indexx[2]]<<"   "<<endl;
            massww[0] = (lvw[indexx[0]]+lvw[indexx[1]]).Mag();
            massww[1] = (lvw[indexx[0]]+lvw[indexx[2]]).Mag();
            massww[2] = (lvw[indexx[1]]+lvw[indexx[2]]).Mag();

            masslvj1 = (lvw[0]+lvw[1]).Mag();
            masslvj2 = (lvw[0]+lvw[2]).Mag();
            massj1j2 = (lvw[1]+lvw[2]).Mag();
}



#endif
