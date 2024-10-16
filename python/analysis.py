import ROOT
import os
import numpy as np

# Enable multi-threading
ROOT.ROOT.EnableImplicitMT(4)

def analyse(fileList):
    
    # Create dataframe(s) from NanoAOD files
    df = {}
    df_tt1 = {}
    df_tt2 = {}
    df_tt3 = {}

    df["2023B"] = ROOT.RDataFrame("Events", fileList)

    h000 = {}
    h001 = {}
    h002 = {}
    h003 = {}
    h004 = {}
    h005 = {}
    h006 = {}
    h007 = {}
    h008 = {}
    h009 = {}
    h010 = {}
    h011 = {}
    h012 = {}
    h013 = {}
    h014 = {}
    h015 = {}
    h016 = {}
    h017 = {}

    # Compile a function to compute the invariant mass of the first 2 particles
    ROOT.gInterpreter.Declare(
    """
    using namespace ROOT;
    using RVecF = ROOT::RVec<float>;
    float InvariantMassOf2(RVecF pt, RVecF eta, RVecF phi, RVecF m) {
        ROOT::Math::PtEtaPhiMVector p0(pt[0], eta[0], phi[0], m[0]);
        ROOT::Math::PtEtaPhiMVector p1(pt[1], eta[1], phi[1], m[1]);
        return (p0 + p1).mass();
    }
    """)

    # Compile a function to compute the invariant mass of 2 different objects
    """
    ROOT.gInterpreter.Declare(

    using namespace ROOT;
    float InvariantMassOfD(float pt0, float eta0, float phi0, float m0, float pt1, float eta1, float phi1, float m1) {
        ROOT::Math::PtEtaPhiMVector p0(pt0, eta0, phi0, m0);
        ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, m1);
        return (p0 + p1).mass();
    }

    """

    ROOT.gInterpreter.Declare(
    """
    using namespace ROOT;
    using RVecF = ROOT::VecOps::RVec<float>;

    // Define a function for the invariant mass calculation
    float InvariantMassOfD(float pt0, float eta0, float phi0, float m0, 
                        float pt1, float eta1, float phi1, float m1) {
        ROOT::Math::PtEtaPhiMVector p0(pt0, eta0, phi0, m0);
        ROOT::Math::PtEtaPhiMVector p1(pt1, eta1, phi1, m1);
        return (p0 + p1).M();  // Use .M() to get the invariant mass
    }

    // Dilepton mass calculation function
    float Dileptonmass(int nTightElectrons, int nTightMuons, 
                    RVecF Electron_pt, 
                    RVecF Electron_eta, 
                    RVecF Electron_phi, 
                    RVecF Electron_mass, 
                    RVecF Muon_pt, 
                    RVecF Muon_eta, 
                    RVecF Muon_phi, 
                    RVecF Muon_mass) {

        // Dielectron mass calculation
        if (nTightElectrons == 2) {
            return InvariantMassOfD(Electron_pt[0], Electron_eta[0], Electron_phi[0], Electron_mass[0], 
                                    Electron_pt[1], Electron_eta[1], Electron_phi[1], Electron_mass[1]);
        }
        // Dimuon mass calculation
        else if (nTightMuons == 2) {
            return InvariantMassOfD(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], 
                                    Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1]);
        }
        // Electron-muon mass calculation
        else if (nTightElectrons == 1 && nTightMuons == 1) {
            return InvariantMassOfD(Electron_pt[0], Electron_eta[0], Electron_phi[0], Electron_mass[0], 
                                    Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0]);
        }

        // Return a default value if no valid combination is found
    //  return -1.0;
    }
    """
    )


    # Compile a function to sum up the first 2 particles
    ROOT.gInterpreter.Declare(
    """
    using namespace ROOT;
    using RVecI = ROOT::RVec<int>;
    int SumOf2(RVecI any) {
        return any[0]+any[1];
    }
    """)

    ROOT.gInterpreter.Declare(
    """
    using namespace ROOT;
    using RVecF = ROOT::RVec<float>;
    float NthOf(RVecF any, int N) {
        return any[N];
    }
    """)

    ROOT.gInterpreter.Declare(
    """
    using namespace ROOT;
    float minp(float pte1) {
        if (pte1 > 44) {
            pte1 = 44;
        }
        return pte1;
    }
    """)

    ROOT.gInterpreter.Declare("""
    float DeltaR(float eta1,float phi1,float eta2,float phi2)
    {
        float deta = eta1-eta2;
        float dphi = phi1-phi2;
        if (dphi > ROOT::Math::Pi()) dphi -= 2*ROOT::Math::Pi();
        if (dphi < -ROOT::Math::Pi()) dphi += 2*ROOT::Math::Pi();
        return sqrt(deta*deta+dphi*dphi);
    }
    """)


    # Declare the cleaningMask function

    ROOT.gInterpreter.Declare("""
    using RVecI = ROOT::RVec<int>;
    using RVecB = ROOT::RVec<bool>;

    RVecB cleaningMask(RVecI indices, int size) {
    RVecB mask(size, true);
    for (int idx : indices) {
        if (idx < 0 ) continue; 
        mask[idx] = false; 
    }
    return mask;
    }
    """)

    ROOT.gInterpreter.Declare("""
    using namespace ROOT;
    using RVecF = ROOT::RVec<float>;

    float compute_jet_lepton_var(RVecF pt, RVecF eta, RVecF phi, RVecF mass, 
                                const RVecF& mu_pt, const RVecF& mu_eta, const RVecF& mu_phi, const RVecF& mu_mass,
                                const RVecF& el_pt, const RVecF& el_eta, const RVecF& el_phi, const RVecF& el_mass,
                                const float met_pt, const float met_phi, unsigned int var)
    {
        if(mu_pt.size() + el_pt.size() == 0) return -1;
        if(pt.size() < 2) return -1;

        ROOT::Math::PtEtaPhiMVector p1(pt[0], eta[0], phi[0], mass[0]);
        ROOT::Math::PtEtaPhiMVector p2(pt[1], eta[1], phi[1], mass[1]);

        for(unsigned int i = 0; i < pt.size(); i++) {
            for(unsigned int j = i + 1; j < pt.size(); j++) {
                if(pt[i] < pt[j]) {
                    std::swap(pt[i], pt[j]);
                    std::swap(eta[i], eta[j]);
                    std::swap(phi[i], phi[j]);
                    std::swap(mass[i], mass[j]);
                }
            }
        }

        float deltaEtaJJ = fabs(p1.Eta() - p2.Eta());
        float maxZ = 0.0;
        float sumHT = p1.Pt() + p2.Pt() + met_pt;

        std::vector<ROOT::Math::PtEtaPhiMVector> p4mom;

        for(unsigned int i = 0; i < mu_pt.size(); i++) {
            p4mom.push_back(ROOT::Math::PtEtaPhiMVector(mu_pt[i], mu_eta[i], mu_phi[i], mu_mass[i]));
            float zVal = fabs(mu_eta[i] - (p1.Eta() + p2.Eta()) / 2.0) / deltaEtaJJ;
            if(zVal > maxZ) maxZ = zVal;
            sumHT += mu_pt[i];
        }

        for(unsigned int i = 0; i < el_pt.size(); i++) {
            p4mom.push_back(ROOT::Math::PtEtaPhiMVector(el_pt[i], el_eta[i], el_phi[i], el_mass[i]));
            float zVal = fabs(el_eta[i] - (p1.Eta() + p2.Eta()) / 2.0) / deltaEtaJJ;
            if(zVal > maxZ) maxZ = zVal;
            sumHT += el_pt[i];
        }

        ROOT::Math::PtEtaPhiMVector p4momVV(met_pt, 0, met_phi, 0);
        ROOT::Math::PtEtaPhiMVector p4momTot = p4momVV + p1 + p2;

        for(const auto& mom : p4mom) {
            p4momVV += mom;
            p4momTot += mom;
        }

        double theVar = 0;
        if(var == 0) theVar = fabs(p4momVV.Eta() - (p1.Eta() + p2.Eta()) / 2.0) / deltaEtaJJ;
        else if(var == 1) theVar = maxZ;
        else if(var == 2) theVar = sumHT;
        else if(var == 3) theVar = p4momVV.Pt();
        else if(var == 4) theVar = p4momTot.Pt();
        else if(var == 5) theVar = fabs(p4momVV.Eta() - p1.Eta());
        else if(var == 6) theVar = fabs(p4momVV.Eta() - p2.Eta());
        else if(var == 7) theVar = (p4momVV.Pt() - (p1 + p2).Pt()) / (p1 + p2).Pt();
        
        return theVar;
    }
    """)




    ROOT.gInterpreter.Declare("""
    using namespace ROOT;
    using RVecF = ROOT::RVec<float>;

        float compute_jet_var(RVecF pt, RVecF eta, RVecF phi, RVecF mass, int var) {
        if(pt.size() < 2) return -1;
        ROOT::Math::PtEtaPhiMVector p1(pt[0], eta[0], phi[0], mass[0]);
        ROOT::Math::PtEtaPhiMVector p2(pt[1], eta[1], phi[1], mass[1]);
        
        // Sort jets by pt (descending order)
        for(unsigned int i = 0; i < pt.size(); i++) {
            for(unsigned int j = i+1; j < pt.size(); j++) {
                if(pt[i] < pt[j]) {
                    std::swap(pt[i], pt[j]);
                    std::swap(eta[i], eta[j]);
                    std::swap(phi[i], phi[j]);
                    std::swap(mass[i], mass[j]);
                }
            }
        }

        double theVar = 0;
        if      (var == 0) theVar = (p1 + p2).M();
        else if (var == 1) theVar = (p1 + p2).Pt();
        else if (var == 2) theVar = fabs(p1.Eta() - p2.Eta());
        else if (var == 4) theVar = p1.Pt();
        else if (var == 5) theVar = p2.Pt();
        else if (var == 6) theVar = fabs(p1.Eta());
        else if (var == 7) theVar = fabs(p2.Eta());
        return theVar;
    }
    """)


    ROOT.gInterpreter.Declare(
    """
    using namespace ROOT;
    using RVecF = ROOT::VecOps::RVec<float>;



    int leadingsub(int electrons, int muons, RVecF fake_Electron_pt, RVecF fake_Muon_pt) {
        float pte1, pte2;

        if (electrons == 2) {
            pte1 = NthOf(fake_Electron_pt, 0);
            pte2 = NthOf(fake_Electron_pt, 1);
            if (pte1 > 25 && pte2 > 20) {
                return 1;
            } else {
                return 0;
            }
        } 
        else if (muons == 2) {
            pte1 = NthOf(fake_Muon_pt, 0);
            pte2 = NthOf(fake_Muon_pt, 1);
            if (pte1 > 25 && pte2 > 20) {
                return 1;
            } else {
                return 0;
            }
        } 
        else if (electrons == 1 && muons == 1) {
            pte1 = NthOf(fake_Electron_pt, 0);
            pte2 = NthOf(fake_Muon_pt, 0);
            if (pte2 > pte1) {
                std::swap(pte1, pte2);
            }
            if (pte1 > 25 && pte2 > 20) {
                return 1;
            } else {
                return 0;
            }
        }

        return 0;  // Default return if none of the above conditions are met
    }
    """
    )

    ROOT.gInterpreter.Declare(
    """
    using namespace ROOT;
    using RVecI = ROOT::VecOps::RVec<int>;  // Assuming genPartFlav is an integer type

    int genflav(int electrons, int muons, RVecI Electron_genPartFlav, RVecI Muon_genPartFlav) {
        int gen1, gen2;

        if (electrons == 2) {
            gen1 = Electron_genPartFlav[0];
            gen2 = Electron_genPartFlav[1];
            if ((gen1 == 1 || gen1 == 15) && (gen2 == 1 || gen2 == 15)) {
                return 1;
            } else {
                return 0;
            }
        } 
        else if (muons == 2) {
            gen1 = Muon_genPartFlav[0];
            gen2 = Muon_genPartFlav[1];
            if ((gen1 == 1 || gen1 == 15) && (gen2 == 1 || gen2 == 15)) {
                return 1;
            } else {
                return 0;
            }
        } 
        else if (electrons == 1 && muons == 1) {
            gen1 = Electron_genPartFlav[0];
            gen2 = Muon_genPartFlav[0];
            if ((gen1 == 1 || gen1 == 15) && (gen2 == 1 || gen2 == 15)) {
                return 1;
            } else {
                return 0;
            }
        }

        return 0;  // Default return if none of the above conditions are met
    }
    """
    )

    for p in ["2023B"]:

        df[p] = df[p].Define("Counts", "1")
        h000[p] = df[p].Histo1D(("Counts_00", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")



    # Triggers

        df[p] = df[p].Define("HLT_mue", "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
        df[p] = df[p].Define("HLT_mumu", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8||HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8")
        df[p] = df[p].Define("HLT_ee", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL||HLT_DoubleEle25_CaloIdL_MW||HLT_DoublePhoton70")
        df[p] = df[p].Define("HLT_mu", "HLT_IsoMu24||HLT_IsoMu27||HLT_Mu50")
        df[p] = df[p].Define("HLT_e", "HLT_Ele30_WPTight_Gsf||HLT_Ele32_WPTight_Gsf||HLT_Ele32_WPTight_Gsf_L1DoubleEG||HLT_Ele35_WPTight_Gsf||HLT_Ele115_CaloIdVT_GsfTrkIdT")


    # Define exclusive mumu, ee and mue triggers
        df[p] = df[p].Define("mumu_trigger", "HLT_mumu || HLT_mu")
        df[p] = df[p].Define("ee_trigger", "HLT_ee || HLT_e")
        df[p] = df[p].Define("mue_trigger", "HLT_mue || (HLT_mu && HLT_e)")
        df[p] = df[p].Define("trigger", "mumu_trigger || ee_trigger || mue_trigger")

        
        df[p] = df[p].Filter("trigger", "events passing Triggers")
        h001[p] = df[p].Histo1D(("Counts_trigger", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")


    # Define loose leptons
        
        df[p] = df[p].Define("loose_electron", "abs(Electron_eta) < 2.5 && Electron_pt > 10 && Electron_cutBased >= 1")
        df[p] = df[p].Define("loose_muon", "abs(Muon_eta) < 2.4 && Muon_pt > 10 && Muon_looseId == true")
        df[p] = df[p].Define("nloosee", "Sum(loose_electron)")
        df[p] = df[p].Define("nloosemu", "Sum(loose_muon)")



    # Defining fake leptons

        df[p] = df[p].Define("fake_muon", "abs(Muon_dxy) < 0.1 && abs(Muon_dz) < 0.2 && abs(Muon_eta) < 2.4 && Muon_pt > 10 && Muon_looseId == true && Muon_mediumId == true && Muon_pfIsoId >= 1 && Muon_jetRelIso < 1.0")
        df[p] = df[p].Define("fake_Muon_pt"              ,"Muon_pt[fake_muon]")
        df[p] = df[p].Define("fake_Muon_eta"             ,"Muon_eta[fake_muon]")
        df[p] = df[p].Define("fake_Muon_phi"             ,"Muon_phi[fake_muon]")
        df[p] = df[p].Define("fake_Muon_dxy"             ,"Muon_dxy[fake_muon]")
        df[p] = df[p].Define("fake_Muon_dz"              ,"Muon_dz[fake_muon]")
        df[p] = df[p].Define("fake_Muon_sip3d"           ,"Muon_sip3d[fake_muon]")
        df[p] = df[p].Define("fake_Muon_mass"            ,"Muon_mass[fake_muon]")
        df[p] = df[p].Define("fake_Muon_charge"          ,"Muon_charge[fake_muon]")
        df[p] = df[p].Define("fake_Muon_jetRelIso"       ,"Muon_jetRelIso[fake_muon]")
        df[p] = df[p].Define("fake_Muon_looseId"         ,"Muon_looseId[fake_muon]")
        df[p] = df[p].Define("fake_Muon_mediumId"        ,"Muon_mediumId[fake_muon]")
        df[p] = df[p].Define("fake_Muon_mediumPromptId"  ,"Muon_mediumPromptId[fake_muon]")
        df[p] = df[p].Define("fake_Muon_tightId"         ,"Muon_tightId[fake_muon]")
        df[p] = df[p].Define("fake_Muon_pfIsoId"         ,"Muon_pfIsoId[fake_muon]")
        df[p] = df[p].Define("fake_Muon_mvaMuID"         ,"Muon_mvaMuID[fake_muon]")
        df[p] = df[p].Define("fake_Muon_miniIsoId"       ,"Muon_miniIsoId[fake_muon]")
        df[p] = df[p].Define("fake_Muon_mvaTTH"          ,"Muon_mvaTTH[fake_muon]")
        df[p] = df[p].Define("fake_Muon_mvaLowPt"        ,"Muon_mvaLowPt[fake_muon]")
        df[p] = df[p].Define("fake_Muon_pfRelIso03_all"  ,"Muon_pfRelIso03_all[fake_muon]")
        df[p] = df[p].Define("fake_Muon_pfRelIso04_all"  ,"Muon_pfRelIso04_all[fake_muon]")
        df[p] = df[p].Define("fake_Muon_miniPFRelIso_all","Muon_miniPFRelIso_all[fake_muon]")
        df[p] = df[p].Define("fake_Muon_nStations"       ,"Muon_nStations[fake_muon]")
        df[p] = df[p].Define("fake_Muon_nTrackerLayers"  ,"Muon_nTrackerLayers[fake_muon]")
        df[p] = df[p].Define("fake_Muon_pfRelIso03_chg"  ,"Muon_pfRelIso03_chg[fake_muon]")
    #   df[p] = df[p].Define("fake_Muon_genpartFlav"     ,"Muon_genPartFlav[fake_muon]")
    #    df[p] = df[p].Define("fake_Muon_p"               ,"computeMomentum(fake_Muon_pt,fake_Muon_eta,fake_Muon_phi,fake_Muon_mass)")
        
        

    #DEFINATION OF FAKE ELECTRON
    
        df[p] = df[p].Define("fake_electron", "(abs(Electron_dxy) < 0.1 && abs(Electron_dz) < 0.2 && abs(Electron_eta) < 2.5 && Electron_pt > 10 && Electron_mvaIso > 0 && Electron_jetRelIso < 1.0)") 
        df[p] = df[p].Define("fake_Electron_pt"                  ,"Electron_pt[fake_electron]")
        df[p] = df[p].Define("fake_Electron_eta"                 ,"Electron_eta[fake_electron]")
        df[p] = df[p].Define("fake_Electron_phi"                 ,"Electron_phi[fake_electron]")
        df[p] = df[p].Define("fake_Electron_mass"                ,"Electron_mass[fake_electron]")
        df[p] = df[p].Define("fake_Electron_charge"              ,"Electron_charge[fake_electron]")
        df[p] = df[p].Define("fake_Electron_jetRelIso"           ,"Electron_jetRelIso[fake_electron]")
        df[p] = df[p].Define("fake_Electron_dxy"                 ,"Electron_dxy[fake_electron]")
        df[p] = df[p].Define("fake_Electron_dz"                  ,"Electron_dz[fake_electron]")
        df[p] = df[p].Define("fake_Electron_sip3d"               ,"Electron_sip3d[fake_electron]")
        df[p] = df[p].Define("fake_Electron_cutBased"            ,"Electron_cutBased[fake_electron]")
        df[p] = df[p].Define("fake_Electron_mvaNoIso_WP80"       ,"Electron_mvaNoIso_WP80[fake_electron]")
        df[p] = df[p].Define("fake_Electron_mvaIso_WP80"         ,"Electron_mvaIso_WP80[fake_electron]")
        df[p] = df[p].Define("fake_Electron_mvaIso_WP90"         ,"Electron_mvaIso_WP90[fake_electron]")
        df[p] = df[p].Define("fake_Electron_tightCharge"         ,"Electron_tightCharge[fake_electron]")
        df[p] = df[p].Define("fake_Electron_mvaTTH"              ,"Electron_mvaTTH[fake_electron]")
        df[p] = df[p].Define("fake_Electron_mvaIso"              ,"Electron_mvaIso[fake_electron]")
        df[p] = df[p].Define("fake_Electron_mvaNoIso"            ,"Electron_mvaNoIso[fake_electron]")
        df[p] = df[p].Define("fake_Electron_pfRelIso03_all"      ,"Electron_pfRelIso03_all[fake_electron]")
        df[p] = df[p].Define("fake_Electron_miniPFRelIso_all"    ,"Electron_miniPFRelIso_all[fake_electron]")
        df[p] = df[p].Define("fake_Electron_hoe"                 ,"Electron_hoe[fake_electron]")
        df[p] = df[p].Define("fake_Electron_r9"                  ,"Electron_r9[fake_electron]")
        df[p] = df[p].Define("fake_Electron_pfRelIso03_chg"      ,"Electron_pfRelIso03_chg[fake_electron]")
    #   df[p] = df[p].Define("fake_Electron_genpartFlav"         ,"Electron_genPartFlav[fake_electron]")
        df[p] = df[p].Define("fake_Electron_seedGain"            ,"Electron_seedGain[fake_electron]")
        df[p] = df[p].Define("nLoose","Sum(loose_muon)+Sum(loose_electron)")
        df[p] = df[p].Define("nFake","Sum(fake_muon)+Sum(fake_electron)")

            
    #Defination of veto 
        df[p] = df[p].Define("veto_electron", "abs(fake_Electron_eta) < 2.5 && fake_Electron_pt > 10 && fake_Electron_cutBased >= 1")
        df[p] = df[p].Define("veto_muon", "abs(fake_Muon_eta) < 2.4 && fake_Muon_pt > 10 && fake_Muon_looseId == true")
        df[p] = df[p].Define("nveto","Sum(veto_muon)+Sum(veto_electron)")      
    
        
    #    Define tight leptons
        df[p] = df[p].Define("tight_muon", "abs(fake_Muon_dxy) < 0.1 && abs(fake_Muon_dz) < 0.2 && abs(fake_Muon_eta) < 2.4 && fake_Muon_pt > 10 && fake_Muon_looseId == true && fake_Muon_mediumPromptId == true && fake_Muon_mvaTTH > 0.5" )   
        df[p] = df[p].Define("tight_electron", " abs(fake_Electron_dxy) < 0.1 && abs(fake_Electron_dz) < 0.2 && abs(fake_Electron_eta) < 2.5 && fake_Electron_pt > 10 && fake_Electron_cutBased >= 2 && fake_Electron_mvaTTH > 0.5 && fake_Electron_tightCharge == 2" ) 
        df[p] = df[p].Define("nTight","Sum(tight_muon)+Sum(tight_electron)") 
        df[p] = df[p].Define("electrons", "Sum(tight_electron)")
        df[p] = df[p].Define("muons", "Sum(tight_muon)")
        


    # Definitions of same-sign ee samples
        df[p] = df[p].Filter("nLoose == 2", "Loose leptons")\
            .Filter("nFake == 2", "Fake leptons")\
            .Filter("nveto == 2")\
            .Filter("Sum(fake_Muon_charge)+Sum(fake_Electron_charge) != 0", "Same-sign leptons")
        h002[p] = df[p].Histo1D(("Counts_fake_lepton", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")

    # Tight lepton filter
        df[p] = df[p].Filter("nTight == 2", "Tight leptons")         
        h003[p] = df[p].Histo1D(("Counts_tight_lepton", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")


    # Defining tau    
        df[p] =df[p].Define("good_tau", "abs(Tau_eta) < 2.5 && Tau_pt > 20 && Tau_idDeepTau2018v2p5VSjet >= 6 && Tau_idDeepTau2018v2p5VSe >= 6 && Tau_idDeepTau2018v2p5VSmu >= 4")\
                .Filter("Sum(good_tau) == 0","No selected hadronic taus")\
                .Define("good_Tau_pt", "Tau_pt[good_tau]")\
                .Define("good_Tau_eta", "Tau_eta[good_tau]")\
                .Define("good_Tau_decayMode", "Tau_decayMode[good_tau]")
        h004[p] = df[p].Histo1D(("Counts_good_taus", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")          
                

    # Leading and subleading leptons
        df[p] = df[p].Define("leading", "leadingsub(electrons, muons, fake_Electron_pt, fake_Muon_pt)")
        df[p] = df[p].Filter("leading == 1", "Filtering on pt1 > 25 and pt2 > 20" )
        h005[p] = df[p].Histo1D(("Counts_leading_subleading", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")

    
    # Dilepton mass
        df[p] = df[p].Define("Dilepton_mass", "Dileptonmass(electrons, muons, fake_Electron_pt, fake_Electron_eta, fake_Electron_phi, fake_Electron_mass, fake_Muon_pt, fake_Muon_eta, fake_Muon_phi, fake_Muon_mass)")
        df[p] = df[p].Filter("Dilepton_mass > 20" , "Dilepton mass > 20")
        h006[p] = df[p].Histo1D(("Counts_dilepton_mass", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")
    #    df[p] = df[p].Filter("nveto == 2")        # I have a doubt with this filter of veto lepton  


    # MC Weight    
    #   df[p] = df[p].Define("genpart", " genflav(electrons, muons, fake_Electron_genpartFlav, fake_Muon_genpartFlav) ")
    #   df[p] = df[p].Filter("genpart == 1", "MC weight")
        
            
            
    # Z veto
        df[p] = df[p].Filter("muons >= 1 || abs(Dilepton_mass - 91.19) > 15" , "Z veto")
        df[p] = df[p].Filter("nJet >= 2") 
        h008[p] = df[p].Histo1D(("Counts_zveto", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")


    # Defining clean jets         
        df[p] = df[p].Define("jet_mask1", "cleaningMask(Muon_jetIdx[loose_muon], nJet)")
        df[p] = df[p].Define("jet_mask2", "cleaningMask(Electron_jetIdx[loose_electron], nJet)")
        df[p] = df[p].Define("clean_jet", "Jet_pt > 10 && jet_mask1 && jet_mask2 && Jet_jetId > 0")
        df[p] = df[p].Define("clean_Jet_pt", "Jet_pt[clean_jet]")
        df[p] = df[p].Define("clean_Jet_eta", "Jet_eta[clean_jet]")
        df[p] = df[p].Define("clean_Jet_phi", "Jet_phi[clean_jet]")
        df[p] = df[p].Define("clean_Jet_mass", "Jet_mass[clean_jet]")
        df[p] = df[p].Define("nclean", "Sum(clean_jet)")


    # Defining VBS jets      
        df[p] = df[p].Define("vbs_Jet", "clean_Jet_pt > 50 && abs(clean_Jet_eta) < 4.7")
        df[p] = df[p].Define("vbs_Jet_pt", "clean_Jet_pt[vbs_Jet]")
        df[p] = df[p].Define("vbs_Jet_eta", "clean_Jet_eta[vbs_Jet]")
        df[p] = df[p].Define("vbs_Jet_phi", "clean_Jet_phi[vbs_Jet]")
        df[p] = df[p].Define("vbs_Jet_mass", "clean_Jet_mass[vbs_Jet]") 
        df[p] = df[p].Define("vbs_jet", "Sum(vbs_Jet)")
        

    # Defining Btag jets
        df[p] = df[p].Define("btagdeepflav", "Jet_btagRobustParTAK4B[clean_jet]")
    #   df[p] = df[p].Define("btagdeepflav", "Jet_btagDeepFlavB[clean_jet]")
        df[p] = df[p].Define("btag_Jet", "clean_Jet_pt > 20 && abs(clean_Jet_eta) < 2.5 && btagdeepflav > 0.0681") 
    #   df[p] = df[p].Define("nbtag", "Sum(btag_Jet)")
        
    

    # Define Zeppenfeld variable
        df[p] =df[p].Define("vbs_zepvv", "compute_jet_lepton_var(vbs_Jet_pt, vbs_Jet_eta, vbs_Jet_phi, vbs_Jet_mass, fake_Muon_pt, fake_Muon_eta, fake_Muon_phi, fake_Muon_mass, fake_Electron_pt, fake_Electron_eta, fake_Electron_phi, fake_Electron_mass, MET_pt, MET_phi, 0)")


    # Filter jets
        df[p] = df[p].Filter("Sum(btag_Jet) == 0", "btag_jet = 0")  
        h007[p] = df[p].Histo1D(("Counts_btag", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")  
        df[p] = df[p].Filter("vbs_jet >= 2", "vbs_jets >= 2")
        h009[p] = df[p].Histo1D(("Counts_vbsjet", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")
        df[p] = df[p].Define("vbs_detajj", "compute_jet_var(vbs_Jet_pt, vbs_Jet_eta, vbs_Jet_phi, vbs_Jet_mass, 2)")    
        df[p] = df[p].Filter("InvariantMassOf2(vbs_Jet_pt, vbs_Jet_eta, vbs_Jet_phi, vbs_Jet_mass) > 500 && vbs_zepvv < 1.0 && vbs_detajj > 2.5", "VBS selection" ) 
        h010[p] = df[p].Histo1D(("Counts_vbs_selection", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")
        df[p] = df[p].Filter("PuppiMET_pt > 30 ", "puppiMET > 30")       
        h011[p] = df[p].Histo1D(("Counts_puppimet", "Counts;1;N_{Events}", 1, 0.5, 1.5), "Counts")



        # Separation of leptons

        df_tt1[p] = df[p].Filter("electrons == 2", "Two electrons")\
                .Define("mom1", "NthOf(fake_Electron_pt,0)")\
                .Define("mom2", "NthOf(fake_Electron_pt,1)")\
                .Define("etae1", "NthOf(fake_Electron_eta,0)")\
                .Define("etae2", "NthOf(fake_Electron_eta,1)")\
                .Define("pte1min", "minp(mom1)")\
                .Define("pte2min", "minp(mom2)")\
                .Define("Dilepton__mass", "Dileptonmass(2, 0, fake_Electron_pt, fake_Electron_eta, fake_Electron_phi, fake_Electron_mass, fake_Muon_pt, fake_Muon_eta, fake_Muon_phi, fake_Muon_mass)")\
                .Define("Dijet_mass", "InvariantMassOf2(vbs_Jet_pt, vbs_Jet_eta, vbs_Jet_phi, vbs_Jet_mass)")

        df_tt2[p] = df[p].Filter("muons == 2", "Two muons")\
                .Define("mom1", "NthOf(fake_Muon_pt,0)")\
                .Define("mom2", "NthOf(fake_Muon_pt,1)")\
                .Define("etae1", "NthOf(fake_Muon_eta,0)")\
                .Define("etae2", "NthOf(fake_Muon_eta,1)")\
                .Define("pte1min", "minp(mom1)")\
                .Define("pte2min", "minp(mom2)")\
                .Define("Dilepton__mass", "Dileptonmass(0, 2, fake_Electron_pt, fake_Electron_eta, fake_Electron_phi, fake_Electron_mass, fake_Muon_pt, fake_Muon_eta, fake_Muon_phi, fake_Muon_mass)")\
                .Define("Dijet_mass", "InvariantMassOf2(vbs_Jet_pt, vbs_Jet_eta, vbs_Jet_phi, vbs_Jet_mass)")
                

        df_tt3[p] = df[p].Filter("electrons == 1 && muons == 1", "One electron + One muon")\
                .Define("mom1", "NthOf(fake_Electron_pt,0)")\
                .Define("mom2", "NthOf(fake_Muon_pt,0)")\
                .Define("etae1", "NthOf(fake_Electron_eta,0)")\
                .Define("etae2", "NthOf(fake_Muon_eta,0)")\
                .Define("pte1min", "minp(mom1)")\
                .Define("pte2min", "minp(mom2)")\
                .Define("Dilepton__mass", "Dileptonmass(1, 1, fake_Electron_pt, fake_Electron_eta, fake_Electron_phi, fake_Electron_mass, fake_Muon_pt, fake_Muon_eta, fake_Muon_phi, fake_Muon_mass)")\
                .Define("Dijet_mass", "InvariantMassOf2(vbs_Jet_pt, vbs_Jet_eta, vbs_Jet_phi, vbs_Jet_mass)")
        

    # Make histograms of dilepton mass and dijet mass
        h012[p] = df_tt1[p].Histo1D(("Dielectron_mass_ee", "Dielectron mass;m_{ee} (GeV);N_{Events}", 100, 0, 1000), "Dilepton__mass")
        h013[p] = df_tt2[p].Histo1D(("Dielectron_mass_mumu", "Dielectron mass;m_{ee} (GeV);N_{Events}", 100, 0, 1000), "Dilepton__mass")
        h014[p] = df_tt3[p].Histo1D(("Dielectron_mass_mue", "Dielectron mass;m_{ee} (GeV);N_{Events}", 100, 0, 1000), "Dilepton__mass")

    
        
        h015[p] = df_tt1[p].Histo1D(("Dijet_mass_ee", "Dijet mass (ee);m_{jj} (GeV);N_{Events}", 100, 0, 2000), "Dijet_mass")
        h016[p] = df_tt2[p].Histo1D(("Dijet_mass_mumu", "Dijet mass (ee);m_{jj} (GeV);N_{Events}", 100, 0, 2000), "Dijet_mass")
        h017[p] = df_tt3[p].Histo1D(("Dijet_mass_mue", "Dijet mass (ee);m_{jj} (GeV);N_{Events}", 100, 0, 2000), "Dijet_mass")



    # Printing reports
    #    report1 = df[p].Report()
    #    report1.Print() 
    #    report2 = df_tt1[p].Report()
    #    report2.Print()
    #    report3 = df_tt2[p].Report()
    #    report3.Print()
    #    report4 = df_tt3[p].Report()
    #    report4.Print()


    # Storing different values
    df_tt1[p].Snapshot("selectionee", "vbsttee.root", ["pte1min", "etae1", "pte2min", "etae2", "Dilepton__mass", "Dijet_mass"])
    df_tt2[p].Snapshot("selectionmumu", "vbsttmumu.root", ["pte1min", "etae1", "pte2min", "etae2", "Dilepton__mass", "Dijet_mass"])
    df_tt3[p].Snapshot("selectionmue", "vbsttmue.root", ["pte1min", "etae1", "pte2min", "etae2", "Dilepton__mass", "Dijet_mass"])

        
    outFile = ROOT.TFile("vbstt.root","recreate")
    outFile.cd()
    for p in ["2023B"]: 
        h000[p].Write()
        h001[p].Write()
        h002[p].Write()
        h003[p].Write()
        h004[p].Write()
        h005[p].Write()
        h006[p].Write()
        h007[p].Write()
        h008[p].Write()
        h009[p].Write()
        h010[p].Write()
        h011[p].Write()
        h012[p].Write()
        h013[p].Write()
        h014[p].Write()
        h015[p].Write()
        h016[p].Write()
        h017[p].Write()
    

    outFile.Close()
######################################################    

