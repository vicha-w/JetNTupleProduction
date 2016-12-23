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

    TFile outFile(strcat("condensed_",fileName),"RECREATE");
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

    Float_t lepton2_pt;
    Float_t lepton2_eta;
    Float_t lepton2_phi;
    Float_t lepton2_E;
    Int_t lepton2_charge;

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

    outTree.Branch("lepton2_pt",&lepton2_pt,"lepton2_pt/F");
    outTree.Branch("lepton2_eta",&lepton2_eta,"lepton2_eta/F");
    outTree.Branch("lepton2_phi",&lepton2_phi,"lepton2_phi/F");
    outTree.Branch("lepton2_E",&lepton2_E,"lepton2_E/F");
    outTree.Branch("lepton2_charge",&lepton2_charge,"lepton2_charge/I");

    outTree.Branch("met_pt",&met_pt,"met_pt/F");
    //outTree.Branch("met_eta",&met_eta,"met_eta/F");
    outTree.Branch("met_phi",&met_phi,"met_phi/F");

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

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

        for (int j = 0; j < nmu; j++)
        {
            if (!lepton1filled)
            {
                lepton1_pt = muon_pt[j];
                lepton1_eta = muon_eta[j];
                lepton1_phi = muon_phi[j];
                lepton1_E = muon_E[j];
                lepton1_charge = muon_charge[j];

                lepton1filled = true;
            }
            else
            {
                if (jet_pt[j] > lepton1_pt)
                {
                    lepton2_pt = lepton1_pt;
                    lepton2_eta = lepton1_eta;
                    lepton2_phi = lepton1_phi;
                    lepton2_E = lepton1_E;
                    lepton2_charge = lepton1_charge;

                    lepton2filled = true;

                    lepton1_pt = muon_pt[j];
                    lepton1_eta = muon_eta[j];
                    lepton1_phi = muon_phi[j];
                    lepton1_E = muon_E[j];
                    lepton1_charge = muon_charge[j];

                    lepton1filled = true;
                }
                else if (!lepton2filled)
                {
                    lepton2_pt = muon_pt[j];
                    lepton2_eta = muon_eta[j];
                    lepton2_phi = muon_phi[j];
                    lepton2_E = muon_E[j];
                    lepton2_charge = muon_charge[j];

                    lepton2filled = true;
                }
                else if (muon_pt[j] > lepton2_pt)
                {
                    lepton2_pt = muon_pt[j];
                    lepton2_eta = muon_eta[j];
                    lepton2_phi = muon_phi[j];
                    lepton2_E = muon_E[j];
                    lepton2_charge = muon_charge[j];

                    lepton2filled = true;
                }
            }
        }

        for (int j = 0; j < nele; j++)
        {
            if (!lepton1filled)
            {
                lepton1_pt = electron_pt[j];
                lepton1_eta = electron_eta[j];
                lepton1_phi = electron_phi[j];
                lepton1_E = electron_E[j];
                lepton1_charge = electron_charge[j];

                lepton1filled = true;
            }
            else
            {
                if (jet_pt[j] > lepton1_pt)
                {
                    lepton2_pt = lepton1_pt;
                    lepton2_eta = lepton1_eta;
                    lepton2_phi = lepton1_phi;
                    lepton2_E = lepton1_E;
                    lepton2_charge = lepton1_charge;

                    lepton2filled = true;

                    lepton1_pt = electron_pt[j];
                    lepton1_eta = electron_eta[j];
                    lepton1_phi = electron_phi[j];
                    lepton1_E = electron_E[j];
                    lepton1_charge = electron_charge[j];

                    lepton1filled = true;
                }
                else if (!lepton2filled)
                {
                    lepton2_pt = electron_pt[j];
                    lepton2_eta = electron_eta[j];
                    lepton2_phi = electron_phi[j];
                    lepton2_E = electron_E[j];
                    lepton2_charge = electron_charge[j];

                    lepton2filled = true;
                }
                else if (electron_pt[j] > lepton2_pt)
                {
                    lepton2_pt = electron_pt[j];
                    lepton2_eta = electron_eta[j];
                    lepton2_phi = electron_phi[j];
                    lepton2_E = electron_E[j];
                    lepton2_charge = electron_charge[j];

                    lepton2filled = true;
                }
            }
        }

        outTree.Fill();
    }

    outFile.cd();
    outTree.Write();
    //outFile.Close();
    //file->Close();
}
