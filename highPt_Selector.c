//
//  highPt_Selector.c
//
//
//  Created by Kalpanie Liyanage on 5/13/20.
//

//To select events for jet studies from SimpleNtupler


#include <stdio.h>
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TText.h"
#include <TString.h>
#include "TEfficiency.h"


void highPt_Selector(){
    
    
    const int mc = 27;
    const int data = 1;
    const int data_f = 3;
    Int_t prev_event = -88;
    
    
    //2017 UL
    TString samp[data] =  {
        "SingleMuonRun2017B-09Aug2019_UL2017-v1",
        //"SingleMuonRun2017C-09Aug2019_UL2017-v1",
       // "SingleMuonRun2017D-09Aug2019_UL2017-v1",
       // "SingleMuonRun2017E-09Aug2019_UL2017-v1",
        //"SingleMuonRun2017F-09Aug2019_UL2017-v1",
        //"SingleMuonRun2017F-09Aug2019_UL2017_EcalRecovery-v1",
     };
     
     
     TString DATA_samples[data_f] =  {
         "zp2mu_histos_1.root",
         "zp2mu_histos_2.root",
         "zp2mu_histos_3.root",
         //"Run2017C_1",
         //"Run2017C_2",
         //"Run2017C_3",
         //"Run2017C_4",
         //"Run2017C_5",
         //"Run2017C_6",
         //"Run2017C_7",
        // "Run2017C_8",
     };
     
 /*    double events_data[data] = {
     1.32161E+8,
     1.579768E+8,
     6.953418E+7,
     1.534541E+8,
     2.370755E+8,
     };
  */
     
     ///////////////////////////////////////MC//////////////////////////////////
     
     TString MC_samples[mc] =  {
         "dy50to120",    //0
         "dy120to200",   //1
         "dy200to400",   //2
         "dy400to800",   //3
         "dy800to1400",  //4
         "dy1400to2300", //5
         "dy2300to3500", //6
         "dy3500to4500", //7
         "dy4500to6000", //8
         "dy6000toInf",  //9
         "ttbar_lep_inclusive",        //10
         "tW",           //11
         "Wantitop",     //12
         "WW",           //13
         "WZ",           //14
         "ZZ",           //15
         "Zjets_filtered", //16 - tautau
         "Zprime_bjets_M-400To1000_8_nLL1",       //17
         "Zprime_bjets_M-400To1000_8_nLL-1",       //18
         "Zprime_bjets_M-400To1000_8_nLR1",       //19
         "Zprime_bjets_M-400To1000_8_nLR-1",       //20
         "Zprime_bjets_M-1000ToInf_8_nLL1",       //21
         "Zprime_bjets_M-1000ToInf_8_nLL-1",       //22
         "Zprime_bjets_M-1000ToInf_8_nLR1",       //23
         "Zprime_bjets_M-1000ToInf_8_nLR-1",       //24
     };
     
   /*
    double events_mc[mc] = {
         2961000,    //0    dy50to120
         100000,     //1    dy120to200
         100000,     //2    dy200to400
         100000,     //3    dy400to800
         100000,     //4    dy800to1400
         100000,     //5    dy1400to2300
         100000,     //6    dy2300to3500
         100000,     //7    dy3500to4500
         100000,     //8    dy4500to6000
         100000,     //9    dy6000toInf
         960752,   //10  ttbar
         4974435,    //11    tW
         5635539,    //12    Wantitop
         7765828,    //13    WW
         3928630,    //14    WZ
         1949768,    //15    ZZ
         30008250,   //16    Wjets
         15296830,  //17 Zjets_filtered
         100000,     //18    Zprime
     };
    
   */
    //2018
    
 /*   TString samp[data] =  {
        "Run2018MuonsOnly_SingleMuonRun2018A_v2-17Sep2018",
        "Run2018MuonsOnly_SingleMuonRun2018B_v1-17Sep2018",
        "Run2018MuonsOnly_SingleMuonRun2018C_v1-17Sep2018",
        "Run2018MuonsOnly_SingleMuonRun2018D_v2-22Jan2019",
    };
    
    
    TString DATA_samples[data] =  {
        "Run2018A",
        "Run2018B",
        "Run2018C",
        "Run2018D",
    };
    
    double events_data[data] = {
        2.389772E+8,
        1.118807E+8,
        1.08034E+8,
        5.118138E+8,
    };
    
    ///////////////////////////////////////MC//////////////////////////////////
    
    TString MC_samples[mc] =  {
        "qcd15to30",        //0
        "qcd30to50",        //1
        "qcd50to80",        //2
        "qcd80to120",       //3
        "qcd120to170",      //4
        "qcd170to300",      //5
        "qcd300to470",      //6
        "qcd470to600",      //7
        "qcd600to800",      //8
        "qcd800to1000",     //9
        "qcd1000to1400",    //10
        "qcd1400to1800",    //11
        "qcd1800to2400",    //12
        "qcd2400to3200",    //13
        "qcd3200toInf",     //14
        "Wjets",            //15
        "Zjets_filtered",   //16
        "ttbar",            //17
        "WW",               //18
        "WZ",               //19
        "ZZ",               //20
        "tW",               //21
        "Wantitop",         //22
        "dy50to120",    //23
        "dy120to200",   //24
        "dy200to400",   //25
        "dy400to800",   //26
        "dy800to1400",  //27
        "dy1400to2300", //28
        "dy2300to3500", //29
        "dy3500to4500", //30
        "dy4500to6000", //31
        "dy6000toInf",  //32
        "Zprime", //33
    };
    
    double events_mc[mc] = {
        19451000,   //0 qcd15to30
        18872000,   //1 qcd30to50
        12909000,   //2 qcd50to80
        29535000,   //3 qcd80to120
        25255000,   //4 qcd120to170
        29710000,   //5 qcd170to300
        41744000,   //6 qcd300to470
        17712000,   //7 qcd470to600
        64061000,   //8 qcd600to800
        37598000,   //9 qcd800to1000
        18485000,   //10    qcd1000to1400
        2160000,    //11    qcd1400to1800
        1445800,    //12    qcd1800to2400
        1440000,    //13    qcd2400to3200
        800000,     //14    qcd3200toInf  151000
        71026861,   //15    Wjets
        31470710,  //16    Zjets
        64310000,   //17    ttbar
        7850000,    //18    WW
        3885000,    //19    WZ
        1979000,    //20    ZZ
        9598000,    //21    tW
        7623000,    //22    Wantitop
        2982000,    //23    dy50to120
        100000,     //24    dy120to200
        100000,     //25    dy200to400
        100000,     //26    dy400to800
        100000,     //27    dy800to1400
        100000,     //28    dy1400to2300
        100000,     //29    dy2300to3500
        100000,     //30    dy3500to4500
        100000,     //31    dy4500to6000
        100000,     //32    dy6000toInf
        100000,     //Zprime
    };
    
*/
    // Declaration of leaf types
    unsigned run;
    unsigned lumi;
    //unsigned long  event;
    ULong64_t  event;
    //int passDileptonSelection;
    UInt_t passDileptonSelection;

    float genWeight;
    float genWeightSign;

    float beamspot_x;
    float beamspot_x_err;
    float beamspot_y;
    float beamspot_y_err;
    float beamspot_z;
    float beamspot_z_err;

    int nvertices;

    float dil_mass;
    float dil_pt;
    float dil_rap;
    float dil_eta;
    float dil_phi;
    float dil_dR;
    float dil_dPhi;
    float dil_lep_pt[2];
    float cos_angle;
    float vertex_chi2;
    float cos_cs;
    float chi_dilepton;
    float phi_cs;

    float vertex_m;
    float vertex_m_err;
    float vertex_x;
    float vertex_x_err;
    float vertex_y;
    float vertex_y_err;
    float vertex_z;
    float vertex_z_err;

    int lep_id[2];
    float lep_p[2];
    float lep_pt[2];
    float lep_et[2];
    float lep_SC_Et[2];
    float lep_pt_err[2];
    float lep_px[2];
    float lep_py[2];
    float lep_pz[2];
    float lep_E[2];
    float lep_eta[2];
    float lep_etaSC[2];
    float lep_phi[2];
    float lep_qOverPt[2];

    //int lep_cutFor[2];
    //int lep_cutFor2018[2];
    Float_t lep_cutFor[2];
    Float_t lep_cutFor2018[2];

    float lep_tk_p[2];
    float lep_tk_pt[2];
    float lep_tk_pt_err[2];
    float lep_tk_px[2];
    float lep_tk_py[2];
    float lep_tk_pz[2];
    float lep_tk_eta[2];
    float lep_tk_phi[2];
    float lep_tk_dz[2];
    float lep_tk_vz[2];
    float lep_tk_chi2[2];
    float lep_tk_ndf[2];
    float lep_tk_qOverPt[2];
    float lep_glb_p[2];
    float lep_glb_pt[2];
    float lep_glb_pt_err[2];
    float lep_glb_px[2];
    float lep_glb_py[2];
    float lep_glb_pz[2];
    float lep_glb_eta[2];
    float lep_glb_phi[2];
    float lep_glb_chi2[2];
    float lep_glb_ndf[2];
    float lep_glb_qOverPt[2];
    float lep_tuneP_p[2];
    float lep_tuneP_pt[2];
    float lep_tuneP_pt_err[2];
    float lep_tuneP_px[2];
    float lep_tuneP_py[2];
    float lep_tuneP_pz[2];
    float lep_tuneP_eta[2];
    float lep_tuneP_phi[2];
    float lep_tuneP_dz[2];
    float lep_tuneP_vz[2];
    float lep_tuneP_chi2[2];
    float lep_tuneP_ndf[2];
    float lep_tuneP_qOverPt[2];

    float lep_Mu50_triggerMatchPt[2];
    float lep_Mu50_triggerMatchEta[2];
    float lep_Mu50_triggerMatchPhi[2];
    float lep_OldMu100_triggerMatchPt[2];
    float lep_OldMu100_triggerMatchEta[2];
    float lep_OldMu100_triggerMatchPhi[2];
    float lep_TkMu100_triggerMatchPt[2];
    float lep_TkMu100_triggerMatchEta[2];
    float lep_TkMu100_triggerMatchPhi[2];
    float lep_Mu17Mu8_Mu8_triggerMatchPt[2];
    float lep_Mu17Mu8_Mu8_triggerMatchEta[2];
    float lep_Mu17Mu8_Mu8_triggerMatchPhi[2];
    float lep_Mu17Mu8_Mu17_triggerMatchPt[2];
    float lep_Mu17Mu8_Mu17_triggerMatchEta[2];
    float lep_Mu17Mu8_Mu17_triggerMatchPhi[2];
    float lep_Mu17Mu8_TrkIsoVVL_triggerMatchPt[2];
    float lep_Mu17Mu8_TrkIsoVVL_triggerMatchEta[2];
    float lep_Mu17Mu8_TrkIsoVVL_triggerMatchPhi[2];
    float lep_Mu17Mu8_DZ_triggerMatchPt[2];
    float lep_Mu17Mu8_DZ_triggerMatchEta[2];
    float lep_Mu17Mu8_DZ_triggerMatchPhi[2];
    float lep_Mu17Mu8_Mass3p8_triggerMatchPt[2];
    float lep_Mu17Mu8_Mass3p8_triggerMatchEta[2];
    float lep_Mu17Mu8_Mass3p8_triggerMatchPhi[2];
    float lep_Mu17Mu8_Mass8_triggerMatchPt[2];
    float lep_Mu17Mu8_Mass8_triggerMatchEta[2];
    float lep_Mu17Mu8_Mass8_triggerMatchPhi[2];
    float lep_Mu37TkMu27_Mu37_triggerMatchPt[2];
    float lep_Mu37TkMu27_Mu37_triggerMatchEta[2];
    float lep_Mu37TkMu27_Mu37_triggerMatchPhi[2];
    float lep_Mu37TkMu27_TkMu27_triggerMatchPt[2];
    float lep_Mu37TkMu27_TkMu27_triggerMatchEta[2];
    float lep_Mu37TkMu27_TkMu27_triggerMatchPhi[2];
    float lep_Mu37TkMu27_DZ_triggerMatchPt[2];
    float lep_Mu37TkMu27_DZ_triggerMatchEta[2];
    float lep_Mu37TkMu27_DZ_triggerMatchPhi[2];

    // for electron
    float lep_L1_triggerMatchEt[2];
    float lep_L1_triggerMatchPt[2];
    float lep_L1_triggerMatchEta[2];
    float lep_L1_triggerMatchPhi[2];
    float lep_HLT_triggerMatchEt[2];
    float lep_HLT_triggerMatchPt[2];
    float lep_HLT_triggerMatchEta[2];
    float lep_HLT_triggerMatchPhi[2];
    float lep_lowestUnprescaledL1[2];

    float lep_chi2dof[2];
    float lep_dB[2];
    float lep_sumPt[2];
    float lep_emEt[2];
    float lep_hadEt[2];
    float lep_hoEt[2];
    float lep_pfIso[2];
    float lep_pfIsoDB[2];
    int lep_timeNdof[2];
    float lep_timeInOut[2];
    float lep_timeOutIn[2];
    float lep_timeInOutErr[2];
    float lep_timeOutInErr[2];
    int lep_heep_id[2];
    float lep_gen_match[2];
    float lep_min_muon_dR[2];
    short lep_tk_numberOfValidTrackerHits[2];
    short lep_tk_numberOfValidTrackerLayers[2];
    short lep_tk_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidTrackerHits[2];
    short lep_glb_numberOfValidTrackerLayers[2];
    short lep_glb_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidMuonHits[2];
    short lep_glb_numberOfValidMuonDTHits[2];
    short lep_glb_numberOfValidMuonCSCHits[2];
    short lep_glb_numberOfValidMuonRPCHits[2];
    short lep_glb_muonStationsWithValidHits[2];
    short lep_glb_dtStationsWithValidHits[2];
    short lep_glb_cscStationsWithValidHits[2];
    short lep_glb_rpcStationsWithValidHits[2];
    short lep_glb_innermostMuonStationWithValidHits[2];
    short lep_glb_outermostMuonStationWithValidHits[2];
    short lep_tuneP_numberOfValidMuonHits[2];
    short lep_tuneP_numberOfValidMuonDTHits[2];
    short lep_tuneP_numberOfValidMuonCSCHits[2];
    short lep_tuneP_numberOfValidMuonRPCHits[2];
    short lep_tuneP_muonStationsWithValidHits[2];
    short lep_tuneP_dtStationsWithValidHits[2];
    short lep_tuneP_cscStationsWithValidHits[2];
    short lep_tuneP_rpcStationsWithValidHits[2];
    short lep_tuneP_innermostMuonStationWithValidHits[2];
    short lep_tuneP_outermostMuonStationWithValidHits[2];
    short lep_numberOfMatches[2];
    short lep_numberOfMatchedStations[2];
    short lep_expectedNnumberOfMatchedStations[2];
    short lep_numberOfMatchedRPCLayers[2];
    //unsigned int lep_stationMask[2];
    Int_t lep_stationMask[2];
    int lep_numberOfChambers[2];
    int lep_numberOfChambersNoRPC[2];
    //unsigned int lep_stationGapMaskDistance[2];
    Int_t lep_stationGapMaskDistance[2];
    //unsigned int lep_stationGapMaskPull[2];
    Int_t lep_stationGapMaskPull[2];
    bool lep_isGlobalMuon[2];
    bool lep_isTrackerMuon[2];

    bool GoodDataRan;
    bool GoodVtx;
    bool METFilter;

    std::vector<std::string> *met_filters;

    float gen_res_mass;
    float gen_res_pt;
    float gen_res_rap;
    float gen_res_eta;
    float gen_res_phi;
    float gen_dil_mass;
    float gen_dil_pt;
    float gen_dil_rap;
    float gen_dil_eta;
    float gen_dil_phi;
    float gen_dil_dR;
    float gen_dil_dPhi;
    float gen_dil_noib_mass;
    float gen_dil_noib_pt;
    float gen_dil_noib_rap;
    float gen_dil_noib_eta;
    float gen_dil_noib_phi;
    float gen_dil_noib_dR;
    float gen_dil_noib_dPhi;
    float gen_lep_p[2];
    float gen_lep_pt[2];
    float gen_lep_px[2];
    float gen_lep_py[2];
    float gen_lep_pz[2];
    float gen_lep_E[2];
    float gen_lep_eta[2];
    float gen_lep_phi[2];
    float gen_lep_qOverPt[2];
    float gen_lep_noib_p[2];
    float gen_lep_noib_pt[2];
    float gen_lep_noib_px[2];
    float gen_lep_noib_py[2];
    float gen_lep_noib_pz[2];
    float gen_lep_noib_E[2];
    float gen_lep_noib_eta[2];
    float gen_lep_noib_phi[2];
    float gen_lep_noib_qOverPt[2];

    float met_pt;
    float met_phi;

    int nJets;
    std::vector<float> *jet_pt;
    std::vector<float> *jet_pt_Uncorrected;
    std::vector<float> *jet_pt_L1FastJet;
    std::vector<float> *jet_pt_L2Relative;
    std::vector<float> *jet_pt_L3Absolute;
    std::vector<float> *jet_pt_L2L3Residual;
    std::vector<float> *jet_eta;
    std::vector<float> *jet_phi;
    std::vector<float> *jet_partonFlavour;
    std::vector<float> *jet_hadronFlavour;
    std::vector<float> *jet_NHF;
    std::vector<float> *jet_NEMF;
    std::vector<float> *jet_CHF;
    std::vector<float> *jet_MUF;
    std::vector<float> *jet_CEMF;
    std::vector<int> *jet_NumConst;
    std::vector<int> *jet_NumNeutralParticles;
    std::vector<int> *jet_CHM;
    std::vector<int> *jet_pileupJetId_fullId;
    std::vector<float> *jet_pileupJetId_fullDiscriminant;
    std::vector<float> *jet_pfDeepCSVJetTags_probb;
    std::vector<float> *jet_pfDeepCSVJetTags_probc;
    std::vector<float> *jet_pfDeepCSVJetTags_probudsg;
    std::vector<float> *jet_pfDeepCSVJetTags_probbb;
    std::vector<float> *jet_pfDeepCSVJetTags_probcc;
    std::vector<float> *jet_pfDeepFlavourJetTags_probb;
    std::vector<float> *jet_pfDeepFlavourJetTags_probbb;
    std::vector<float> *jet_pfDeepFlavourJetTags_problepb;
    std::vector<float> *jet_pfDeepFlavourJetTags_probc;
    std::vector<float> *jet_pfDeepFlavourJetTags_probuds;
    std::vector<float> *jet_pfDeepFlavourJetTags_probg;
    std::vector<float> *jet_genp_pt;
    std::vector<float> *jet_genp_eta;
    std::vector<float> *jet_genp_phi;
    std::vector<float> *jet_genp_energy;
    std::vector<int> *jet_genp_charge;
    std::vector<int> *jet_genp_pdgId;
    std::vector<int> *jet_genp_status;
    std::vector<float> *jet_genj_pt;
    std::vector<float> *jet_genj_eta;
    std::vector<float> *jet_genj_phi;
    std::vector<float> *jet_genj_energy;

    int nGenp;
    std::vector<float> *genp_pt;
    std::vector<float> *genp_eta;
    std::vector<float> *genp_phi;
    std::vector<float> *genp_energy;
    std::vector<int> *genp_charge;
    std::vector<int> *genp_pdgId;
    std::vector<int> *genp_status;
    std::vector<int> *genp_isPrompt;
    std::vector<int> *genp_isTauDecayProduct;
    std::vector<int> *genp_isPromptTauDecayProduct;
    std::vector<int> *genp_isDecayedLeptonHadron;
    std::vector<int> *genp_isPromptFinalState;
    std::vector<int> *genp_isDirectPromptTauDecayProductFinalState;
    std::vector<int> *genp_isHardProcess;
    std::vector<int> *genp_isLastCopy;
    std::vector<int> *genp_isLastCopyBeforeFSR;
    std::vector<int> *genp_isPromptDecayed;
    std::vector<int> *genp_fromHardProcessBeforeFSR;
    std::vector<int> *genp_fromHardProcessDecayed;
    std::vector<int> *genp_fromHardProcessFinalState;
    
    
 
    
    
    
    //strat of looping over Data samples
    for(int j=0; j<data; j++){
        
        for(int r=0; r<data_f; r++){
        
        std::cout<<"opening.. "<<DATA_samples[r]<<".root"<<std::endl;
        TChain *treeDATA = new TChain("SimpleNtupler/t");
        //treeDATA->Add("/eos/user/k/kaliyana/2017_data/Zprime_bbll/UL/"+ samp[j] +"/"+ DATA_samples[r] +".root");
        treeDATA->Add(samp[j] +"/"+ DATA_samples[r] +".root");
            
        
        // Set branch addresses and branch pointers
        
        
        treeDATA->SetBranchAddress("run", &run);
        treeDATA->SetBranchAddress("lumi", &lumi);
        treeDATA->SetBranchAddress("event", &event);
        treeDATA->SetBranchAddress("passDileptonSelection", &passDileptonSelection);
        treeDATA->SetBranchAddress("beamspot_x", &beamspot_x);
        treeDATA->SetBranchAddress("beamspot_x_err", &beamspot_x_err);
        treeDATA->SetBranchAddress("beamspot_y", &beamspot_y);
        treeDATA->SetBranchAddress("beamspot_y_err", &beamspot_y_err);
        treeDATA->SetBranchAddress("beamspot_z", &beamspot_z);
        treeDATA->SetBranchAddress("beamspot_z_err", &beamspot_z_err);
        treeDATA->SetBranchAddress("nvertices", &nvertices);
        treeDATA->SetBranchAddress("dil_mass", &dil_mass);
        treeDATA->SetBranchAddress("dil_pt", &dil_pt);
        treeDATA->SetBranchAddress("dil_rap", &dil_rap);
        treeDATA->SetBranchAddress("dil_eta", &dil_eta);
        treeDATA->SetBranchAddress("dil_phi", &dil_phi);
        treeDATA->SetBranchAddress("dil_dR", &dil_dR);
        treeDATA->SetBranchAddress("dil_dPhi", &dil_dPhi);
        treeDATA->SetBranchAddress("dil_lep_pt", dil_lep_pt);
        treeDATA->SetBranchAddress("cos_angle", &cos_angle);
        treeDATA->SetBranchAddress("vertex_chi2", &vertex_chi2);
        treeDATA->SetBranchAddress("cos_cs", &cos_cs);
        treeDATA->SetBranchAddress("chi_dilepton", &chi_dilepton);
        treeDATA->SetBranchAddress("phi_cs", &phi_cs);
        treeDATA->SetBranchAddress("vertex_m", &vertex_m);
        treeDATA->SetBranchAddress("vertex_m_err", &vertex_m_err);
        treeDATA->SetBranchAddress("vertex_x", &vertex_x);
        treeDATA->SetBranchAddress("vertex_x_err", &vertex_x_err);
        treeDATA->SetBranchAddress("vertex_y", &vertex_y);
        treeDATA->SetBranchAddress("vertex_y_err", &vertex_y_err);
        treeDATA->SetBranchAddress("vertex_z", &vertex_z);
        treeDATA->SetBranchAddress("vertex_z_err", &vertex_z_err);
        treeDATA->SetBranchAddress("lep_id", lep_id);
        treeDATA->SetBranchAddress("lep_heep_id", lep_heep_id);
        treeDATA->SetBranchAddress("lep_gen_match", lep_gen_match);
        treeDATA->SetBranchAddress("lep_p", lep_p);
        treeDATA->SetBranchAddress("lep_pt", lep_pt);
        treeDATA->SetBranchAddress("lep_pt_err", lep_pt_err);
        treeDATA->SetBranchAddress("lep_px", lep_px);
        treeDATA->SetBranchAddress("lep_py", lep_py);
        treeDATA->SetBranchAddress("lep_pz", lep_pz);
        treeDATA->SetBranchAddress("lep_E", lep_E);
        treeDATA->SetBranchAddress("lep_eta", lep_eta);
        treeDATA->SetBranchAddress("lep_etaSC", lep_etaSC);
        treeDATA->SetBranchAddress("lep_et", lep_et);
        treeDATA->SetBranchAddress("lep_SC_Et", lep_SC_Et);
        treeDATA->SetBranchAddress("lep_phi", lep_phi);
        treeDATA->SetBranchAddress("lep_qOverPt", lep_qOverPt);
        treeDATA->SetBranchAddress("lep_cutFor", lep_cutFor);
        treeDATA->SetBranchAddress("lep_cutFor2018", lep_cutFor2018);
        treeDATA->SetBranchAddress("lep_tk_p", lep_tk_p);
        treeDATA->SetBranchAddress("lep_tk_pt", lep_tk_pt);
        treeDATA->SetBranchAddress("lep_tk_pt_err", lep_tk_pt_err);
        treeDATA->SetBranchAddress("lep_tk_px", lep_tk_px);
        treeDATA->SetBranchAddress("lep_tk_py", lep_tk_py);
        treeDATA->SetBranchAddress("lep_tk_pz", lep_tk_pz);
        treeDATA->SetBranchAddress("lep_tk_eta", lep_tk_eta);
        treeDATA->SetBranchAddress("lep_tk_phi", lep_tk_phi);
        treeDATA->SetBranchAddress("lep_tk_dz", lep_tk_dz);
        treeDATA->SetBranchAddress("lep_tk_vz", lep_tk_vz);
        treeDATA->SetBranchAddress("lep_tk_chi2", lep_tk_chi2);
        treeDATA->SetBranchAddress("lep_tk_ndf", lep_tk_ndf);
        treeDATA->SetBranchAddress("lep_tk_qOverPt", lep_tk_qOverPt);
        treeDATA->SetBranchAddress("lep_glb_p", lep_glb_p);
        treeDATA->SetBranchAddress("lep_glb_pt", lep_glb_pt);
        treeDATA->SetBranchAddress("lep_glb_pt_err", lep_glb_pt_err);
        treeDATA->SetBranchAddress("lep_glb_px", lep_glb_px);
        treeDATA->SetBranchAddress("lep_glb_py", lep_glb_py);
        treeDATA->SetBranchAddress("lep_glb_pz", lep_glb_pz);
        treeDATA->SetBranchAddress("lep_glb_eta", lep_glb_eta);
        treeDATA->SetBranchAddress("lep_glb_phi", lep_glb_phi);
        treeDATA->SetBranchAddress("lep_glb_chi2", lep_glb_chi2);
        treeDATA->SetBranchAddress("lep_glb_ndf", lep_glb_ndf);
        treeDATA->SetBranchAddress("lep_glb_qOverPt", lep_glb_qOverPt);
        treeDATA->SetBranchAddress("lep_tuneP_p", lep_tuneP_p);
        treeDATA->SetBranchAddress("lep_tuneP_pt", lep_tuneP_pt);
        treeDATA->SetBranchAddress("lep_tuneP_pt_err", lep_tuneP_pt_err);
        treeDATA->SetBranchAddress("lep_tuneP_px", lep_tuneP_px);
        treeDATA->SetBranchAddress("lep_tuneP_py", lep_tuneP_py);
        treeDATA->SetBranchAddress("lep_tuneP_pz", lep_tuneP_pz);
        treeDATA->SetBranchAddress("lep_tuneP_eta", lep_tuneP_eta);
        treeDATA->SetBranchAddress("lep_tuneP_phi", lep_tuneP_phi);
        treeDATA->SetBranchAddress("lep_tuneP_chi2", lep_tuneP_chi2);
        treeDATA->SetBranchAddress("lep_tuneP_ndf", lep_tuneP_ndf);
        treeDATA->SetBranchAddress("lep_tuneP_qOverPt", lep_tuneP_qOverPt);
        treeDATA->SetBranchAddress("lep_Mu50_triggerMatchPt", lep_Mu50_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu50_triggerMatchEta", lep_Mu50_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu50_triggerMatchPhi", lep_Mu50_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_OldMu100_triggerMatchPt", lep_OldMu100_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_OldMu100_triggerMatchEta", lep_OldMu100_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_OldMu100_triggerMatchPhi", lep_OldMu100_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_TkMu100_triggerMatchPt", lep_TkMu100_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_TkMu100_triggerMatchEta", lep_TkMu100_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_TkMu100_triggerMatchPhi", lep_TkMu100_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mu8_triggerMatchPt", lep_Mu17Mu8_Mu8_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mu8_triggerMatchEta", lep_Mu17Mu8_Mu8_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mu8_triggerMatchPhi", lep_Mu17Mu8_Mu8_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mu17_triggerMatchPt", lep_Mu17Mu8_Mu17_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mu17_triggerMatchEta", lep_Mu17Mu8_Mu17_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mu17_triggerMatchPhi", lep_Mu17Mu8_Mu17_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_TrkIsoVVL_triggerMatchPt", lep_Mu17Mu8_TrkIsoVVL_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_TrkIsoVVL_triggerMatchEta", lep_Mu17Mu8_TrkIsoVVL_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_TrkIsoVVL_triggerMatchPhi", lep_Mu17Mu8_TrkIsoVVL_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_DZ_triggerMatchPt", lep_Mu17Mu8_DZ_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_DZ_triggerMatchEta", lep_Mu17Mu8_DZ_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_DZ_triggerMatchPhi", lep_Mu17Mu8_DZ_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mass3p8_triggerMatchPt", lep_Mu17Mu8_Mass3p8_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mass3p8_triggerMatchEta", lep_Mu17Mu8_Mass3p8_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mass3p8_triggerMatchPhi", lep_Mu17Mu8_Mass3p8_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mass8_triggerMatchPt", lep_Mu17Mu8_Mass8_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mass8_triggerMatchEta", lep_Mu17Mu8_Mass8_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu17Mu8_Mass8_triggerMatchPhi", lep_Mu17Mu8_Mass8_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_Mu37TkMu27_Mu37_triggerMatchPt", lep_Mu37TkMu27_Mu37_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu37TkMu27_Mu37_triggerMatchEta", lep_Mu37TkMu27_Mu37_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu37TkMu27_Mu37_triggerMatchPhi", lep_Mu37TkMu27_Mu37_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_Mu37TkMu27_TkMu27_triggerMatchPt", lep_Mu37TkMu27_TkMu27_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu37TkMu27_TkMu27_triggerMatchEta", lep_Mu37TkMu27_TkMu27_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu37TkMu27_TkMu27_triggerMatchPhi", lep_Mu37TkMu27_TkMu27_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_Mu37TkMu27_DZ_triggerMatchPt", lep_Mu37TkMu27_DZ_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_Mu37TkMu27_DZ_triggerMatchEta", lep_Mu37TkMu27_DZ_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_Mu37TkMu27_DZ_triggerMatchPhi", lep_Mu37TkMu27_DZ_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_L1_triggerMatchEt", lep_L1_triggerMatchEt);
        treeDATA->SetBranchAddress("lep_L1_triggerMatchPt", lep_L1_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_L1_triggerMatchEta", lep_L1_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_L1_triggerMatchPhi", lep_L1_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_HLT_triggerMatchEt", lep_HLT_triggerMatchEt);
        treeDATA->SetBranchAddress("lep_HLT_triggerMatchPt", lep_HLT_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_HLT_triggerMatchEta", lep_HLT_triggerMatchEta);
        treeDATA->SetBranchAddress("lep_HLT_triggerMatchPhi", lep_HLT_triggerMatchPhi);
        treeDATA->SetBranchAddress("lep_lowestUnprescaledL1", lep_lowestUnprescaledL1);
        treeDATA->SetBranchAddress("lep_chi2dof", lep_chi2dof);
        treeDATA->SetBranchAddress("lep_dB", lep_dB);
        treeDATA->SetBranchAddress("lep_sumPt", lep_sumPt);
        treeDATA->SetBranchAddress("lep_emEt", lep_emEt);
        treeDATA->SetBranchAddress("lep_hadEt", lep_hadEt);
        treeDATA->SetBranchAddress("lep_hoEt", lep_hoEt);
        treeDATA->SetBranchAddress("lep_pfIso", lep_pfIso);
        treeDATA->SetBranchAddress("lep_pfIsoDB", lep_pfIsoDB);
        treeDATA->SetBranchAddress("lep_timeNdof", lep_timeNdof);
        treeDATA->SetBranchAddress("lep_timeInOut", lep_timeInOut);
        treeDATA->SetBranchAddress("lep_timeOutIn", lep_timeOutIn);
        treeDATA->SetBranchAddress("lep_timeInOutErr", lep_timeInOutErr);
        treeDATA->SetBranchAddress("lep_timeOutInErr", lep_timeOutInErr);
        treeDATA->SetBranchAddress("lep_min_muon_dR", lep_min_muon_dR);
        treeDATA->SetBranchAddress("lep_tk_numberOfValidTrackerHits", lep_tk_numberOfValidTrackerHits);
        treeDATA->SetBranchAddress("lep_tk_numberOfValidTrackerLayers", lep_tk_numberOfValidTrackerLayers);
        treeDATA->SetBranchAddress("lep_tk_numberOfValidPixelHits", lep_tk_numberOfValidPixelHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidTrackerHits", lep_glb_numberOfValidTrackerHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidTrackerLayers", lep_glb_numberOfValidTrackerLayers);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidPixelHits", lep_glb_numberOfValidPixelHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidMuonHits", lep_glb_numberOfValidMuonHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidMuonDTHits", lep_glb_numberOfValidMuonDTHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidMuonCSCHits", lep_glb_numberOfValidMuonCSCHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidMuonRPCHits", lep_glb_numberOfValidMuonRPCHits);
        treeDATA->SetBranchAddress("lep_glb_muonStationsWithValidHits", lep_glb_muonStationsWithValidHits);
        treeDATA->SetBranchAddress("lep_glb_dtStationsWithValidHits", lep_glb_dtStationsWithValidHits);
        treeDATA->SetBranchAddress("lep_glb_cscStationsWithValidHits", lep_glb_cscStationsWithValidHits);
        treeDATA->SetBranchAddress("lep_glb_rpcStationsWithValidHits", lep_glb_rpcStationsWithValidHits);
        treeDATA->SetBranchAddress("lep_glb_innermostMuonStationWithValidHits", lep_glb_innermostMuonStationWithValidHits);
        treeDATA->SetBranchAddress("lep_glb_outermostMuonStationWithValidHits", lep_glb_outermostMuonStationWithValidHits);
        treeDATA->SetBranchAddress("lep_tuneP_numberOfValidMuonHits", lep_tuneP_numberOfValidMuonHits);
        treeDATA->SetBranchAddress("lep_tuneP_numberOfValidMuonDTHits", lep_tuneP_numberOfValidMuonDTHits);
        treeDATA->SetBranchAddress("lep_tuneP_numberOfValidMuonCSCHits", lep_tuneP_numberOfValidMuonCSCHits);
        treeDATA->SetBranchAddress("lep_tuneP_numberOfValidMuonRPCHits", lep_tuneP_numberOfValidMuonRPCHits);
        treeDATA->SetBranchAddress("lep_tuneP_muonStationsWithValidHits", lep_tuneP_muonStationsWithValidHits);
        treeDATA->SetBranchAddress("lep_tuneP_dtStationsWithValidHits", lep_tuneP_dtStationsWithValidHits);
        treeDATA->SetBranchAddress("lep_tuneP_cscStationsWithValidHits", lep_tuneP_cscStationsWithValidHits);
        treeDATA->SetBranchAddress("lep_tuneP_rpcStationsWithValidHits", lep_tuneP_rpcStationsWithValidHits);
        treeDATA->SetBranchAddress("lep_tuneP_innermostMuonStationWithValidHits", lep_tuneP_innermostMuonStationWithValidHits);
        treeDATA->SetBranchAddress("lep_tuneP_outermostMuonStationWithValidHits", lep_tuneP_outermostMuonStationWithValidHits);
        treeDATA->SetBranchAddress("lep_numberOfMatches", lep_numberOfMatches);
        treeDATA->SetBranchAddress("lep_numberOfMatchedStations", lep_numberOfMatchedStations);
        treeDATA->SetBranchAddress("lep_expectedNnumberOfMatchedStations",lep_expectedNnumberOfMatchedStations);
        treeDATA->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
        treeDATA->SetBranchAddress("lep_stationMask", lep_stationMask);
        treeDATA->SetBranchAddress("lep_numberOfChambers", lep_numberOfChambers);
        treeDATA->SetBranchAddress("lep_numberOfChambersNoRPC", lep_numberOfChambersNoRPC);
        treeDATA->SetBranchAddress("lep_stationGapMaskDistance", lep_stationGapMaskDistance);
        treeDATA->SetBranchAddress("lep_stationGapMaskPull", lep_stationGapMaskPull);
        treeDATA->SetBranchAddress("lep_isGlobalMuon", lep_isGlobalMuon);
        treeDATA->SetBranchAddress("lep_isTrackerMuon", lep_isTrackerMuon);
        treeDATA->SetBranchAddress("GoodDataRan", &GoodDataRan);
        treeDATA->SetBranchAddress("GoodVtx", &GoodVtx);
        treeDATA->SetBranchAddress("METFilter", &METFilter);
        treeDATA->SetBranchAddress("met_filters", &met_filters);
        treeDATA->SetBranchAddress("met_pt", &met_pt);
        treeDATA->SetBranchAddress("met_phi", &met_phi);
        treeDATA->SetBranchAddress("nJets", &nJets);
        treeDATA->SetBranchAddress("jet_pt", &jet_pt);
        treeDATA->SetBranchAddress("jet_pt_Uncorrected", &jet_pt_Uncorrected);
        treeDATA->SetBranchAddress("jet_pt_L1FastJet", &jet_pt_L1FastJet);
        treeDATA->SetBranchAddress("jet_pt_L2Relative", &jet_pt_L2Relative);
        treeDATA->SetBranchAddress("jet_pt_L3Absolute", &jet_pt_L3Absolute);
        treeDATA->SetBranchAddress("jet_pt_L2L3Residual", &jet_pt_L2L3Residual);
        treeDATA->SetBranchAddress("jet_eta", &jet_eta);
        treeDATA->SetBranchAddress("jet_phi", &jet_phi);
        treeDATA->SetBranchAddress("jet_partonFlavour", &jet_partonFlavour);
        treeDATA->SetBranchAddress("jet_hadronFlavour", &jet_hadronFlavour);
        treeDATA->SetBranchAddress("jet_NHF", &jet_NHF);
        treeDATA->SetBranchAddress("jet_NEMF", &jet_NEMF);
        treeDATA->SetBranchAddress("jet_CHF", &jet_CHF);
        treeDATA->SetBranchAddress("jet_MUF", &jet_MUF);
        treeDATA->SetBranchAddress("jet_CEMF", &jet_CEMF);
        treeDATA->SetBranchAddress("jet_NumConst", &jet_NumConst);
        treeDATA->SetBranchAddress("jet_NumNeutralParticles", &jet_NumNeutralParticles);
        treeDATA->SetBranchAddress("jet_CHM", &jet_CHM);
        treeDATA->SetBranchAddress("jet_pileupJetId_fullId", &jet_pileupJetId_fullId);
        treeDATA->SetBranchAddress("jet_pileupJetId_fullDiscriminant", &jet_pileupJetId_fullDiscriminant);
        treeDATA->SetBranchAddress("jet_pfDeepCSVJetTags_probb", &jet_pfDeepCSVJetTags_probb);
        treeDATA->SetBranchAddress("jet_pfDeepCSVJetTags_probc", &jet_pfDeepCSVJetTags_probc);
        treeDATA->SetBranchAddress("jet_pfDeepCSVJetTags_probudsg", &jet_pfDeepCSVJetTags_probudsg);
        treeDATA->SetBranchAddress("jet_pfDeepCSVJetTags_probbb", &jet_pfDeepCSVJetTags_probbb);
        treeDATA->SetBranchAddress("jet_pfDeepCSVJetTags_probcc", &jet_pfDeepCSVJetTags_probcc);
        treeDATA->SetBranchAddress("jet_pfDeepFlavourJetTags_probb", &jet_pfDeepFlavourJetTags_probb);
        treeDATA->SetBranchAddress("jet_pfDeepFlavourJetTags_probbb", &jet_pfDeepFlavourJetTags_probbb);
        treeDATA->SetBranchAddress("jet_pfDeepFlavourJetTags_problepb", &jet_pfDeepFlavourJetTags_problepb);
        treeDATA->SetBranchAddress("jet_pfDeepFlavourJetTags_probc", &jet_pfDeepFlavourJetTags_probc);
        treeDATA->SetBranchAddress("jet_pfDeepFlavourJetTags_probuds", &jet_pfDeepFlavourJetTags_probuds);
        treeDATA->SetBranchAddress("jet_pfDeepFlavourJetTags_probg", &jet_pfDeepFlavourJetTags_probg);
        treeDATA->SetBranchAddress("jet_genp_pt", &jet_genp_pt);
        treeDATA->SetBranchAddress("jet_genp_eta", &jet_genp_eta);
        treeDATA->SetBranchAddress("jet_genp_phi", &jet_genp_phi);
        treeDATA->SetBranchAddress("jet_genp_energy", &jet_genp_energy);
        treeDATA->SetBranchAddress("jet_genp_charge", &jet_genp_charge);
        treeDATA->SetBranchAddress("jet_genp_pdgId", &jet_genp_pdgId);
        treeDATA->SetBranchAddress("jet_genp_status", &jet_genp_status);
        treeDATA->SetBranchAddress("jet_genj_pt", &jet_genj_pt);
        treeDATA->SetBranchAddress("jet_genj_eta", &jet_genj_eta);
        treeDATA->SetBranchAddress("jet_genj_phi", &jet_genj_phi);
        treeDATA->SetBranchAddress("jet_genj_energy", &jet_genj_energy);

        
        cout << " Creating new root-file ..."<< std::endl;
        TFile *newFile = new TFile(DATA_samples[r]+"_Zprime.root","recreate");
        
        cout << " Creating new tree ..."<< endl;
        //TTree *tree = (TTree*)treeDATA->GetTree()->CloneTree(0);
        TChain *newchain = (TChain*)treeDATA->CloneTree(0);
        TTree *tree = newchain->GetTree();
        
        if (treeDATA == 0 && tree == 0) return;
        Long64_t nentries = Long64_t(treeDATA->GetEntries());
        cout << "Number of entries : " << nentries << endl;
        
        
        
        Int_t selected = 0;
        //looping over entries
        for (Long64_t jentry=0; jentry<nentries; jentry++) {
            //cout << jentry << endl;
            //Int_t ientry = LoadTree(jentry); //in case of a TChain,
            treeDATA->GetEntry(jentry);
            cout << jentry << endl;
            if (jentry % 10000 == 0) cout << jentry << endl;
            
           
            
            //hightpT selection cuts or muons
            if(
               GoodVtx &&
               fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 &&
               lep_pt[0]>53. && lep_pt[1]>53. &&  //offline reconstruction pt threshold
               lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 &&
               lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 &&
               fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 &&
               (lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2))) &&
               (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2))) &&
               (lep_glb_numberOfValidMuonHits[0]>0 || lep_tuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_tuneP_numberOfValidMuonHits[1]>0) &&
               lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 &&
               lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 &&
               lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
               lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 &&
               cos_angle>-0.9998 &&
               lep_id[0]*lep_id[1]<0 &&
               //trigger pt threshold
               ((lep_Mu50_triggerMatchPt[0]>50. || lep_OldMu100_triggerMatchPt[0]>50. || lep_TkMu100_triggerMatchPt[0]>50.)
                                 || (lep_Mu50_triggerMatchPt[1]>50. || lep_OldMu100_triggerMatchPt[1]>50. || lep_TkMu100_triggerMatchPt[1]>50.))
               //fabs(lep_tk_dz[0]) < 1.0 && fabs(lep_tk_dz[1]) < 1.0
                
              
               
               )
          {
                
                if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
                prev_event = event;
                
                selected++;
                tree->Fill();
                
         }
            
            
        }// end loop over entries
        
        std::cout<< DATA_samples[r]<<".root"<<std::endl;
        std::cout<<"Number of entires in the raw data sample : "<<nentries<< " Number of entires after selection : " <<selected<< " Percentage : " << (double) selected/nentries <<std::endl;
        
        newFile->cd();
        tree->Write();
        newFile->Close();
        
        
      }//end looping over combined files
        
    } //end looping over Data samples
    
    
    
    
    
    
    
}


