#ifndef OpenDataTreeProducer_h
#define OpenDataTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;

class OpenDataTreeProducer : public edm::EDAnalyzer 
{
  public:

    explicit OpenDataTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void beginRun(edm::Run const &, edm::EventSetup const& iSetup);
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endRun(edm::Run const &, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~OpenDataTreeProducer();


  private:  

    // Function to help sort the jet wrt. pT
    static bool cmp_patjets(const pat::Jet &pj1, const pat::Jet &pj2) {
        return pj1.pt() > pj2.pt();
    }

    //---- Configurable parameters --------  
    bool            mIsMCarlo;
    bool            mUseGenInfo;
    bool            mPrintTriggerMenu;
    int             mMinNPFJets;
    double          mMinPFPtJets, mMinGenPt, mMaxY, mMinJJMass;
    int             mGoodVtxNdof;
    double          mGoodVtxZ; 
    edm::InputTag   mPFak5JetsName;
    edm::InputTag   mPFak7JetsName;
    
    // ---- PF Jet input tags ----- //
    edm::InputTag   mGenJetsName;
    edm::InputTag   mSrcPFRho;
    edm::InputTag   mPFMET; 
    edm::InputTag   mOfflineVertices;

    // ---- BTag discriminator input tag
    std::string     mBTagDiscriminator;

    // ---- Muon and Electron input tags
    edm::InputTag   mMuonName;
    edm::InputTag   mElectronName;

    // ---- Electorn selection criteria
    double          mMinPtElectrons;
    double          mMaxEtaElectrons;
    double          mElectronTIP; // Electron transverse impact parameter
    std::string     mElectronID;
    double          mElectronDeltaR;
    double          mMaxREI;

    // ---- Muon selection criteria
    double          mMinPtMuons;
    double          mMaxEtaMuons;
    bool            mGlobalMuon;
    bool            mTrackerMuon;
    std::string     mMuonID;
    int             mNumValidHitsMuon;
    double          mChi2OverNdof;
    double          mMuonTIP; // Muon transverse impact parameter
    double          mMaxRMI;
    
    //---- Trigger----------------------
    std::string                 processName_;
    std::vector<std::string>    triggerNames_;
    std::vector<unsigned int>   triggerIndex_;
    edm::InputTag               triggerResultsTag_;
    HLTConfigProvider           hltConfig_;
    
    // Output variables
    edm::Service<TFileService>  fs;
    TTree                       *mTree;

    
    //---- TTree variables --------
    
    static const UInt_t kMaxNjet = 128;
    static const UInt_t kMaxNtrg = 32;
    static const UInt_t kMaxNmu = 128;
    static const UInt_t kMaxNele = 128;

    // PF jets
    UInt_t njet;
    Float_t jet_pt[kMaxNjet];
    Float_t jet_eta[kMaxNjet];
    Float_t jet_phi[kMaxNjet];
    Float_t jet_E[kMaxNjet];
    Bool_t jet_tightID[kMaxNjet];
    Float_t jet_area[kMaxNjet];
    Float_t jet_jes[kMaxNjet];
    Int_t jet_igen[kMaxNjet];

    Float_t jet_btag[kMaxNjet];

    // PF jets
    UInt_t njet_ak7;
    Float_t jet_pt_ak7[kMaxNjet];
    Float_t jet_eta_ak7[kMaxNjet];
    Float_t jet_phi_ak7[kMaxNjet];
    Float_t jet_E_ak7[kMaxNjet];
    Float_t jet_area_ak7[kMaxNjet];
    Float_t jet_jes_ak7[kMaxNjet];
    Int_t ak7_to_ak5[kMaxNjet];

    // Jet composition
    Float_t chf[kMaxNjet];
   	Float_t nhf[kMaxNjet];
   	Float_t phf[kMaxNjet];
   	Float_t elf[kMaxNjet];
   	Float_t muf[kMaxNjet];
   	Float_t hf_hf[kMaxNjet];
   	Float_t hf_phf[kMaxNjet];
   	Int_t hf_hm[kMaxNjet];
   	Int_t hf_phm[kMaxNjet];
   	Int_t chm[kMaxNjet];
   	Int_t nhm[kMaxNjet];
   	Int_t phm[kMaxNjet];
   	Int_t elm[kMaxNjet];
   	Int_t mum[kMaxNjet];   
    Float_t beta[kMaxNjet];   
    Float_t bstar[kMaxNjet];

    // Generated jets
    UInt_t ngen;
    Float_t gen_pt[kMaxNjet];
    Float_t gen_eta[kMaxNjet];
    Float_t gen_phi[kMaxNjet];
    Float_t gen_E[kMaxNjet];

    // Event identification
    UInt_t run;
    UInt_t lumi;
    ULong64_t event;

    // Triggers
    UInt_t ntrg;
    Bool_t triggers[kMaxNtrg];
    std::vector<std::string> triggernames;
    UInt_t prescales[kMaxNtrg];

    // MET, SuMET, rho, eventfilter
    Float_t met;
    Float_t sumet;
    Float_t rho;

    // MC variables
    Float_t pthat;
    Float_t mcweight;
    
    // Detailed MET
    Float_t met_et;
    Float_t met_phi;

    // Muons
    UInt_t nmu;
    Float_t muon_pt[kMaxNmu];
    Float_t muon_eta[kMaxNmu];
    Float_t muon_phi[kMaxNmu];
    Float_t muon_E[kMaxNmu];
    Int_t muon_charge[kMaxNmu];
    Bool_t muon_ID[kMaxNmu];
    Float_t muon_TIP[kMaxNmu];

    // Electrons
    UInt_t nele;
    Float_t electron_pt[kMaxNele];
    Float_t electron_eta[kMaxNele];
    Float_t electron_phi[kMaxNele];
    Float_t electron_E[kMaxNele];
    Int_t electron_charge[kMaxNele];
    Float_t electron_ID[kMaxNele];
    Float_t electron_TIP[kMaxNmu];

    // Variables before cut
    UInt_t b_njet;
    Float_t b_jet_pt[kMaxNjet];
    Float_t b_jet_eta[kMaxNjet];
    Float_t b_jet_phi[kMaxNjet];
    Float_t b_jet_E[kMaxNjet];
    Bool_t b_jet_tightID[kMaxNjet];
    Float_t b_jet_area[kMaxNjet];
    Float_t b_jet_jes[kMaxNjet];
    Int_t b_jet_igen[kMaxNjet];
    Float_t b_jet_btag[kMaxNjet];

    UInt_t b_nmu;
    Float_t b_muon_pt[kMaxNmu];
    Float_t b_muon_eta[kMaxNmu];
    Float_t b_muon_phi[kMaxNmu];
    Float_t b_muon_E[kMaxNmu];
    Int_t b_muon_charge[kMaxNmu];
    Bool_t b_muon_ID[kMaxNmu];
    Float_t b_muon_TIP[kMaxNmu];

    UInt_t b_nele;
    Float_t b_electron_pt[kMaxNele];
    Float_t b_electron_eta[kMaxNele];
    Float_t b_electron_phi[kMaxNele];
    Float_t b_electron_E[kMaxNele];
    Int_t b_electron_charge[kMaxNele];
    Float_t b_electron_ID[kMaxNele];
    Float_t b_electron_TIP[kMaxNmu];

};

#endif
