#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <stdio.h>
#include <string>

const float bTagDiscrimLevel = 0.75;

void condenseNTuple(const char* fileName, const char* treeName="ak5ak7/OpenDataTree")
{
    TFile *file = new TFile(fileName);
    TTree *tree = (TTree*) file->Get(treeName);

    UInt_t njet, nmu, nele;
    const int maxElements = 64;

    Float_t jet_pt[maxElements];
    Float_t jet_eta[maxElements];
    Float_t jet_phi[maxElements];
    Float_t jet_E[maxElements];
    Float_t jet_bTag[maxElements];

    Float_t muon_pt[maxElements];
    Float_t muon_eta[maxElements];
    Float_t muon_phi[maxElements];
    Float_t muon_E[maxElements];
    Int_t muon_charge[maxElements];

    Float_t electron_pt[maxElements];
    Float_t electron_eta[maxElements];
    Float_t electron_phi[maxElements];
    Float_t electron_E[maxElements];
    Int_t electron_charge[maxElements];

    Float_t met_pt;
    //Float_t met_eta;
    Float_t met_phi;

    tree->SetBranchAddress("njet",&njet);
    tree->SetBranchAddress("nmu",&nmu);
    tree->SetBranchAddress("nele",&nele);

    tree->SetBranchAddress("jet_pt",jet_pt);
    tree->SetBranchAddress("jet_eta",jet_eta);
    tree->SetBranchAddress("jet_phi",jet_phi);
    tree->SetBranchAddress("jet_E",jet_E);
    tree->SetBranchAddress("jet_btag",jet_bTag);

    tree->SetBranchAddress("muon_pt",muon_pt);
    tree->SetBranchAddress("muon_eta",muon_eta);
    tree->SetBranchAddress("muon_phi",muon_phi);
    tree->SetBranchAddress("muon_E",muon_E);
    tree->SetBranchAddress("muon_charge",muon_charge);

    tree->SetBranchAddress("electron_pt",electron_pt);
    tree->SetBranchAddress("electron_eta",electron_eta);
    tree->SetBranchAddress("electron_phi",electron_phi);
    tree->SetBranchAddress("electron_E",electron_E);
    tree->SetBranchAddress("electron_charge",electron_charge);

    tree->SetBranchAddress("met_pt",&met_pt);
    //tree->SetBranchAddress("met_eta",&met_eta);
    tree->SetBranchAddress("met_phi",&met_phi);

    TFile outFile((string("condensed_"+string(fileName)).c_str()),"RECREATE");
    TTree outTree("condensedNTuple","Condensed NTuple");

    Float_t jet1_pt;
    Float_t jet1_eta;
    Float_t jet1_phi;
    Float_t jet1_E;
    Float_t jet1_bTag;

    Float_t jet2_pt;
    Float_t jet2_eta;
    Float_t jet2_phi;
    Float_t jet2_E;
    Float_t jet2_bTag;

    Float_t lepton1_pt;
    Float_t lepton1_eta;
    Float_t lepton1_phi;
    Float_t lepton1_E;
    Int_t lepton1_charge;
    Bool_t lepton1_isMuon;

    Float_t lepton2_pt;
    Float_t lepton2_eta;
    Float_t lepton2_phi;
    Float_t lepton2_E;
    Int_t lepton2_charge;
    Bool_t lepton2_isMuon;

    outTree.Branch("jet1_pt",&jet1_pt,"jet1_pt/F");
    outTree.Branch("jet1_eta",&jet1_eta,"jet1_eta/F");
    outTree.Branch("jet1_phi",&jet1_phi,"jet1_phi/F");
    outTree.Branch("jet1_E",&jet1_E,"jet1_E/F");
    outTree.Branch("jet1_bTag",&jet1_bTag,"jet1_bTag/F");

    outTree.Branch("jet2_pt",&jet2_pt,"jet2_pt/F");
    outTree.Branch("jet2_eta",&jet2_eta,"jet2_eta/F");
    outTree.Branch("jet2_phi",&jet2_phi,"jet2_phi/F");
    outTree.Branch("jet2_E",&jet2_E,"jet2_E/F");
    outTree.Branch("jet2_bTag",&jet2_bTag,"jet2_bTag/F");

    outTree.Branch("lepton1_pt",&lepton1_pt,"lepton1_pt/F");
    outTree.Branch("lepton1_eta",&lepton1_eta,"lepton1_eta/F");
    outTree.Branch("lepton1_phi",&lepton1_phi,"lepton1_phi/F");
    outTree.Branch("lepton1_E",&lepton1_E,"lepton1_E/F");
    outTree.Branch("lepton1_charge",&lepton1_charge,"lepton1_charge/I");
    outTree.Branch("lepton1_isMuon",&lepton1_isMuon,"lepton1_isMuon/O");

    outTree.Branch("lepton2_pt",&lepton2_pt,"lepton2_pt/F");
    outTree.Branch("lepton2_eta",&lepton2_eta,"lepton2_eta/F");
    outTree.Branch("lepton2_phi",&lepton2_phi,"lepton2_phi/F");
    outTree.Branch("lepton2_E",&lepton2_E,"lepton2_E/F");
    outTree.Branch("lepton2_charge",&lepton2_charge,"lepton2_charge/I");
    outTree.Branch("lepton2_isMuon",&lepton2_isMuon,"lepton2_isMuon/O");

    outTree.Branch("met_pt",&met_pt,"met_pt/F");
    //outTree.Branch("met_eta",&met_eta,"met_eta/F");
    outTree.Branch("met_phi",&met_phi,"met_phi/F");

    Long64_t nentries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nentries; entry++)
    {
        tree->GetEntry(entry);

        if (njet < 2 || nele + nmu < 2) continue;

        bool jet1filled = false;
        bool jet2filled = false;
        bool lepton1filled = false;
        bool lepton2filled = false;

        for (int j = 0; j < njet; j++)
        {
            if (jet_bTag[j] < bTagDiscrimLevel) continue;
            if (!jet1filled)
            {
                jet1_pt = jet_pt[j];
                jet1_eta = jet_eta[j];
                jet1_phi = jet_phi[j];
                jet1_E = jet_E[j];
                jet1_bTag = jet_bTag[j];

                jet1filled = true;
            }
            else
            {
                if (jet_pt[j] > jet1_pt)
                {
                    jet2_pt = jet1_pt;
                    jet2_eta = jet1_eta;
                    jet2_phi = jet1_phi;
                    jet2_E = jet1_E;
                    jet2_bTag = jet1_bTag;

                    jet2filled = true;

                    jet1_pt = jet_pt[j];
                    jet1_eta = jet_eta[j];
                    jet1_phi = jet_phi[j];
                    jet1_E = jet_E[j];
                    jet1_bTag = jet_bTag[j];

                    jet1filled = true;
                }
                else if (!jet2filled)
                {
                    jet2_pt = jet_pt[j];
                    jet2_eta = jet_eta[j];
                    jet2_phi = jet_phi[j];
                    jet2_E = jet_E[j];
                    jet2_bTag = jet_bTag[j];

                    jet2filled = true;
                }
                else if (jet_pt[j] > jet2_pt)
                {
                    jet2_pt = jet_pt[j];
                    jet2_eta = jet_eta[j];
                    jet2_phi = jet_phi[j];
                    jet2_E = jet_E[j];
                    jet2_bTag = jet_bTag[j];

                    jet2filled = true;
                }
            }
        }

        if (!(jet1filled && jet2filled)) for (int j = 0; j < njet; j++)
        {
            if (!jet1filled)
            {
                jet1_pt = jet_pt[j];
                jet1_eta = jet_eta[j];
                jet1_phi = jet_phi[j];
                jet1_E = jet_E[j];
                jet1_bTag = jet_bTag[j];

                jet1filled = true;
            }
            else
            {
                if (jet_pt[j] > jet1_pt)
                {
                    jet2_pt = jet1_pt;
                    jet2_eta = jet1_eta;
                    jet2_phi = jet1_phi;
                    jet2_E = jet1_E;
                    jet2_bTag = jet1_bTag;

                    jet2filled = true;

                    jet1_pt = jet_pt[j];
                    jet1_eta = jet_eta[j];
                    jet1_phi = jet_phi[j];
                    jet1_E = jet_E[j];
                    jet1_bTag = jet_bTag[j];

                    jet1filled = true;
                }
                else if (!jet2filled)
                {
                    jet2_pt = jet_pt[j];
                    jet2_eta = jet_eta[j];
                    jet2_phi = jet_phi[j];
                    jet2_E = jet_E[j];
                    jet2_bTag = jet_bTag[j];

                    jet2filled = true;
                }
                else if (jet_pt[j] > jet2_pt)
                {
                    jet2_pt = jet_pt[j];
                    jet2_eta = jet_eta[j];
                    jet2_phi = jet_phi[j];
                    jet2_E = jet_E[j];
                    jet2_bTag = jet_bTag[j];

                    jet2filled = true;
                }
            }
        }

        if (jet1_bTag < 0) jet1_bTag = 0;
        if (jet2_bTag < 0) jet2_bTag = 0;

        UInt_t p_count = 0;
        UInt_t m_count = 0;

        for (int i = 0; i < nmu; i++)
        {
            if (muon_charge[i] > 0)p_count++;
            else m_count++;
        }
        
        for (int i = 0; i < nele; i++)
        {
            if (electron_charge[i] > 0) p_count++;
            else m_count++;
        }

        if (p_count == 0 || m_count == 0) continue;
        
        const UInt_t leptonSize = 16;

        Float_t leptonP_pt[leptonSize];
        Float_t leptonP_eta[leptonSize];
        Float_t leptonP_phi[leptonSize];
        Float_t leptonP_E[leptonSize];
        Bool_t  leptonP_isMuon[leptonSize];

        Float_t leptonM_pt[leptonSize];
        Float_t leptonM_eta[leptonSize];
        Float_t leptonM_phi[leptonSize];
        Float_t leptonM_E[leptonSize];
        Bool_t  leptonM_isMuon[leptonSize];

        p_count = 0;
        m_count = 0;

        for (int i = 0; i < nmu; i++)
        {
            if (muon_charge[i] > 0)
            {
                leptonP_pt[p_count]  = muon_pt[i];
                leptonP_eta[p_count] = muon_eta[i];
                leptonP_phi[p_count] = muon_phi[i];
                leptonP_E[p_count]   = muon_E[i];
                leptonP_isMuon[p_count] = true;
                p_count++;
            }
            else
            {
                leptonM_pt[m_count]  = muon_pt[i];
                leptonM_eta[m_count] = muon_eta[i];
                leptonM_phi[m_count] = muon_phi[i];
                leptonM_E[m_count]   = muon_E[i];
                leptonM_isMuon[p_count] = true;
                m_count++;
            }
        }

        for (int i = 0; i < nele; i++)
        {
            if (electron_charge[i] > 0)
            {
                leptonP_pt[p_count]  = electron_pt[i];
                leptonP_eta[p_count] = electron_eta[i];
                leptonP_phi[p_count] = electron_phi[i];
                leptonP_E[p_count]   = electron_E[i];
                leptonP_isMuon[p_count] = false;
                p_count++;
            }
            else
            {
                leptonM_pt[m_count]  = electron_pt[i];
                leptonM_eta[m_count] = electron_eta[i];
                leptonM_phi[m_count] = electron_phi[i];
                leptonM_E[m_count]   = electron_E[i];
                leptonM_isMuon[p_count] = false;
                m_count++;
            }
        }

        UInt_t bestlepP;
        UInt_t bestlepM;
        Float_t bestPt = -1;

        TLorentzVector leptonP, leptonM;

        for (int i = 0; i < p_count; i++) for (int j = 0; j < m_count; j++)
        {
            leptonP.SetPtEtaPhiE(leptonP_pt[i], leptonP_eta[i], leptonP_phi[i], leptonP_E[i]);
            leptonM.SetPtEtaPhiE(leptonM_pt[j], leptonM_eta[j], leptonM_phi[j], leptonM_E[j]);

            if ((leptonP + leptonM).Pt() > bestPt)
            {
                bestlepP = i;
                bestlepM = j;
                bestPt   = (leptonP + leptonM).Pt();
            }
        }

        if (bestPt > 0)
        {
            if (leptonP_pt[bestlepP] > leptonM_pt[bestlepM])
            {
                lepton1_pt  = leptonP_pt[bestlepP];
                lepton1_eta = leptonP_eta[bestlepP];
                lepton1_phi = leptonP_phi[bestlepP];
                lepton1_E   = leptonP_E[bestlepP];
                lepton1_isMuon = leptonP_isMuon[bestlepP];
                lepton1_charge = 1;

                lepton2_pt  = leptonM_pt[bestlepM];
                lepton2_eta = leptonM_eta[bestlepM];
                lepton2_phi = leptonM_phi[bestlepM];
                lepton2_E   = leptonM_E[bestlepM];
                lepton2_isMuon = leptonM_isMuon[bestlepM];
                lepton2_charge = -1;
            }
            else
            {
                lepton1_pt  = leptonM_pt[bestlepM];
                lepton1_eta = leptonM_eta[bestlepM];
                lepton1_phi = leptonM_phi[bestlepM];
                lepton1_E   = leptonM_E[bestlepM];
                lepton1_isMuon = leptonM_isMuon[bestlepP];
                lepton1_charge = -1;

                lepton2_pt  = leptonP_pt[bestlepP];
                lepton2_eta = leptonP_eta[bestlepP];
                lepton2_phi = leptonP_phi[bestlepP];
                lepton2_E   = leptonP_E[bestlepP];
                lepton2_isMuon = leptonP_isMuon[bestlepM];
                lepton2_charge = 1;
            }
        }

        outTree.Fill();
    }

    outFile.cd();
    outTree.Write();
    //outFile.Close();
    //file->Close();
}
