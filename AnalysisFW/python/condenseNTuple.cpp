#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <stdio.h>
#include <string>
#include <Math/Polynomial.h>

const float bTagDiscrimLevel = 0.244; // CombinedSecondaryVertexLoose
const double massW = 80.385; // PDG 2016
const bool offShell = true; // Assume top quarks are off-shell. (Don't satisfy E^2 = m^2 + p^2)

void condenseNTuple(const char* fileName, const char* treeName="ak5ak7/OpenDataTree")
{
    TFile *file = new TFile(fileName);
    TTree *tree = (TTree*) file->Get(treeName);

    UInt_t njet, nmu, nele;
    const int maxElements = 128;

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

    tree->SetBranchAddress("met_et",&met_pt);
    //tree->SetBranchAddress("met_eta",&met_eta);
    tree->SetBranchAddress("met_phi",&met_phi);

    TFile outFile((string("condensed_"+string(fileName)).c_str()),"RECREATE");
    TTree outTree("condensedNTuple","Condensed NTuple");

    // Jets
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

    // Leptons
    Float_t lepton1_pt;
    Float_t lepton1_eta;
    Float_t lepton1_phi;
    Float_t lepton1_E;
    Float_t lepton1_charge;
    Bool_t lepton1_isMuon;

    Float_t lepton2_pt;
    Float_t lepton2_eta;
    Float_t lepton2_phi;
    Float_t lepton2_E;
    Float_t lepton2_charge;
    Bool_t lepton2_isMuon;

    // Top candidates
    /*
    Float_t top1_pt;
    Float_t top1_eta;
    Float_t top1_phi;
    Float_t top1_E;
    Float_t tbar1_pt;
    Float_t tbar1_eta;
    Float_t tbar1_phi;
    Float_t tbar1_E;

    Float_t top2_pt;
    Float_t top2_eta;
    Float_t top2_phi;
    Float_t top2_E;
    Float_t tbar2_pt;
    Float_t tbar2_eta;
    Float_t tbar2_phi;
    Float_t tbar2_E;

    Float_t top1_mass;
    Float_t top2_mass;
    */

    Float_t top_pt;
    Float_t top_eta;
    Float_t top_phi;
    Float_t top_E;
    Float_t tbar_pt;
    Float_t tbar_eta;
    Float_t tbar_phi;
    Float_t tbar_E;
    
    Float_t top_mass;

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
    outTree.Branch("lepton1_charge",&lepton1_charge,"lepton1_charge/F");
    outTree.Branch("lepton1_isMuon",&lepton1_isMuon,"lepton1_isMuon/O");

    outTree.Branch("lepton2_pt",&lepton2_pt,"lepton2_pt/F");
    outTree.Branch("lepton2_eta",&lepton2_eta,"lepton2_eta/F");
    outTree.Branch("lepton2_phi",&lepton2_phi,"lepton2_phi/F");
    outTree.Branch("lepton2_E",&lepton2_E,"lepton2_E/F");
    outTree.Branch("lepton2_charge",&lepton2_charge,"lepton2_charge/F");
    outTree.Branch("lepton2_isMuon",&lepton2_isMuon,"lepton2_isMuon/O");

    outTree.Branch("met_pt",&met_pt,"met_pt/F");
    //outTree.Branch("met_eta",&met_eta,"met_eta/F");
    outTree.Branch("met_phi",&met_phi,"met_phi/F");

    /*
    outTree.Branch("top1_pt",&top1_pt,"top1_pt/F");
    outTree.Branch("top1_eta",&top1_eta,"top1_eta/F");
    outTree.Branch("top1_phi",&top1_phi,"top1_phi/F");
    outTree.Branch("top1_E",&top1_E,"top1_E/F");
    outTree.Branch("tbar1_pt",&tbar1_pt,"tbar1_pt/F");
    outTree.Branch("tbar1_eta",&tbar1_eta,"tbar1_eta/F");
    outTree.Branch("tbar1_phi",&tbar1_phi,"tbar1_phi/F");
    outTree.Branch("tbar1_E",&tbar1_E,"tbar1_E/F");

    outTree.Branch("top2_pt",&top2_pt,"top2_pt/F");
    outTree.Branch("top2_eta",&top2_eta,"top2_eta/F");
    outTree.Branch("top2_phi",&top2_phi,"top2_phi/F");
    outTree.Branch("top2_E",&top2_E,"top2_E/F");
    outTree.Branch("tbar2_pt",&tbar2_pt,"tbar2_pt/F");
    outTree.Branch("tbar2_eta",&tbar2_eta,"tbar2_eta/F");
    outTree.Branch("tbar2_phi",&tbar2_phi,"tbar2_phi/F");
    outTree.Branch("tbar2_E",&tbar2_E,"tbar2_E/F");

    outTree.Branch("top1_mass",&top1_mass,"top1_mass/F");
    outTree.Branch("top2_mass",&top2_mass,"top2_mass/F");
    */

    outTree.Branch("top_pt",&top_pt,"top_pt/F");
    outTree.Branch("top_eta",&top_eta,"top_eta/F");
    outTree.Branch("top_phi",&top_phi,"top_phi/F");
    outTree.Branch("top_E",&top_E,"top_E/F");
    outTree.Branch("tbar_pt",&tbar_pt,"tbar_pt/F");
    outTree.Branch("tbar_eta",&tbar_eta,"tbar_eta/F");
    outTree.Branch("tbar_phi",&tbar_phi,"tbar_phi/F");
    outTree.Branch("tbar_E",&tbar_E,"tbar_E/F");

    outTree.Branch("top_mass",&top_mass,"top_mass/F");

    Long64_t nentries = tree->GetEntries();
    int validLEntries = 0;
    for (Long64_t entry = 0; entry < nentries; entry++)
    {
        tree->GetEntry(entry);

        if (njet < 2 || nele + nmu < 2) continue;

        bool jet1filled = false;
        bool jet2filled = false;
        bool lepton1filled = false;
        bool lepton2filled = false;

        int chooseJetInd[2] = {-1, -1};

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
                chooseJetInd[0] = j;
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
                    chooseJetInd[1] = chooseJetInd[0];

                    jet1_pt = jet_pt[j];
                    jet1_eta = jet_eta[j];
                    jet1_phi = jet_phi[j];
                    jet1_E = jet_E[j];
                    jet1_bTag = jet_bTag[j];

                    jet1filled = true;
                    chooseJetInd[0] = j;
                }
                else if (!jet2filled)
                {
                    jet2_pt = jet_pt[j];
                    jet2_eta = jet_eta[j];
                    jet2_phi = jet_phi[j];
                    jet2_E = jet_E[j];
                    jet2_bTag = jet_bTag[j];

                    jet2filled = true;
                    chooseJetInd[1] = j;
                }
                else if (jet_pt[j] > jet2_pt)
                {
                    jet2_pt = jet_pt[j];
                    jet2_eta = jet_eta[j];
                    jet2_phi = jet_phi[j];
                    jet2_E = jet_E[j];
                    jet2_bTag = jet_bTag[j];

                    jet2filled = true;
                    chooseJetInd[1] = j;
                }
            }
        }

        if (!(jet1filled && jet2filled)) for (int j = 0; j < njet; j++)
        {
            if (j == chooseJetInd[0] || j == chooseJetInd[1]) continue;

            if (!jet1filled)
            {
                jet1_pt = jet_pt[j];
                jet1_eta = jet_eta[j];
                jet1_phi = jet_phi[j];
                jet1_E = jet_E[j];
                jet1_bTag = jet_bTag[j];

                jet1filled = true;
                chooseJetInd[0] = j;
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
                    chooseJetInd[1] = chooseJetInd[0];

                    jet1_pt = jet_pt[j];
                    jet1_eta = jet_eta[j];
                    jet1_phi = jet_phi[j];
                    jet1_E = jet_E[j];
                    jet1_bTag = jet_bTag[j];

                    jet1filled = true;
                    chooseJetInd[0] = j;
                }
                else if (!jet2filled)
                {
                    jet2_pt = jet_pt[j];
                    jet2_eta = jet_eta[j];
                    jet2_phi = jet_phi[j];
                    jet2_E = jet_E[j];
                    jet2_bTag = jet_bTag[j];

                    jet2filled = true;
                    chooseJetInd[1] = j;
                }
                else if (jet_pt[j] > jet2_pt)
                {
                    jet2_pt = jet_pt[j];
                    jet2_eta = jet_eta[j];
                    jet2_phi = jet_phi[j];
                    jet2_E = jet_E[j];
                    jet2_bTag = jet_bTag[j];

                    jet2filled = true;
                    chooseJetInd[1] = j;
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
        else validLEntries++;

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
                leptonM_isMuon[m_count] = true;
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
                leptonM_isMuon[m_count] = false;
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
                lepton1_charge = 1.;

                lepton2_pt  = leptonM_pt[bestlepM];
                lepton2_eta = leptonM_eta[bestlepM];
                lepton2_phi = leptonM_phi[bestlepM];
                lepton2_E   = leptonM_E[bestlepM];
                lepton2_isMuon = leptonM_isMuon[bestlepM];
                lepton2_charge = -1.;
            }
            else
            {
                lepton1_pt  = leptonM_pt[bestlepM];
                lepton1_eta = leptonM_eta[bestlepM];
                lepton1_phi = leptonM_phi[bestlepM];
                lepton1_E   = leptonM_E[bestlepM];
                lepton1_isMuon = leptonM_isMuon[bestlepP];
                lepton1_charge = -1.;

                lepton2_pt  = leptonP_pt[bestlepP];
                lepton2_eta = leptonP_eta[bestlepP];
                lepton2_phi = leptonP_phi[bestlepP];
                lepton2_E   = leptonP_E[bestlepP];
                lepton2_isMuon = leptonP_isMuon[bestlepM];
                lepton2_charge = 1.;
            }
        }

        Float_t *jetPtPt[2] = {&jet1_pt, &jet2_pt};
        Float_t *jetEtaPt[2] = {&jet1_eta, &jet2_eta};
        Float_t *jetPhiPt[2] = {&jet1_phi, &jet2_phi};
        Float_t *jetEPt[2] = {&jet1_E, &jet2_E};

        Float_t *lepPtPt[2] = {&lepton1_pt, &lepton2_pt};
        Float_t *lepEtaPt[2] = {&lepton1_eta, &lepton2_eta};
        Float_t *lepPhiPt[2] = {&lepton1_phi, &lepton2_phi};
        Float_t *lepEPt[2] = {&lepton1_E, &lepton2_E};
        Float_t *lepChargePt[2] = {&lepton1_charge, &lepton2_charge};

        /*
        Float_t *topPtPt[2] = {&top1_pt, &top2_pt};
        Float_t *topEtaPt[2] = {&top1_eta, &top2_eta};
        Float_t *topPhiPt[2] = {&top1_phi, &top2_phi};
        Float_t *topEPt[2] = {&top1_E, &top2_E};

        Float_t *tbarPtPt[2] = {&tbar1_pt, &tbar2_pt};
        Float_t *tbarEtaPt[2] = {&tbar1_eta, &tbar2_eta};
        Float_t *tbarPhiPt[2] = {&tbar1_phi, &tbar2_phi};
        Float_t *tbarEPt[2] = {&tbar1_E, &tbar2_E};

        Float_t *topMassPt[2] = {&top1_mass, &top2_mass};
        */

        unsigned short useJetInd[2] = {0,1};

        double bestDelta = -1.;
        double bestDeltaTop = 0.;
        double bestMassTop = -1.;
        double deltaRBestMassTop1, deltaRBestMassTop2;

        for (int ind = 0; ind < 2; ind++)
        {
            TLorentzVector jetb, jetbbar, leptonPVect, leptonMVect, METVect;
            jetb.SetPtEtaPhiE(*jetPtPt[useJetInd[ind]], *jetEtaPt[useJetInd[ind]], *jetPhiPt[useJetInd[ind]], *jetEPt[useJetInd[ind]]);
            jetbbar.SetPtEtaPhiE(*jetPtPt[useJetInd[ind^1]], *jetEtaPt[useJetInd[ind^1]], *jetPhiPt[useJetInd[ind^1]], *jetEPt[useJetInd[ind^1]]);

            if (lepton1_charge > 0)
            {
                leptonPVect.SetPtEtaPhiE(*lepPtPt[0], *lepEtaPt[0], *lepPhiPt[0], *lepEPt[0]);
                leptonMVect.SetPtEtaPhiE(*lepPtPt[1], *lepEtaPt[1], *lepPhiPt[1], *lepEPt[1]);
            }
            else
            {
                leptonPVect.SetPtEtaPhiE(*lepPtPt[1], *lepEtaPt[1], *lepPhiPt[1], *lepEPt[1]);
                leptonMVect.SetPtEtaPhiE(*lepPtPt[0], *lepEtaPt[0], *lepPhiPt[0], *lepEPt[0]);
            }

            METVect.SetPtEtaPhiE(met_pt, 0, met_phi, met_pt);

            /*
            double bestDelta = -1.;
            double bestDeltaTop = 0.;
            double bestMassTop = -1.;
            double deltaRBestMassTop1, deltaRBestMassTop2;
            */

            double a[4], b[4], c[3][3], dp[3][3], d[3][3], h[5];
			double Eb, Ebbar, Elp, Elm;
			double Emetx, Emety;

			double pbx, pby, pbz;
			double pbbarx, pbbary, pbbarz;
			double plpx, plpy, plpz;
			double plmx, plmy, plmz;

			double massb, massbbar, masslp, masslm;
			double massWp, massWm, massnu, massnubar;

            Eb    = jetb.E();
            Ebbar = jetbbar.E();
            Elp   = leptonPVect.E();
            Elm   = leptonMVect.E();

            Emetx = METVect.X();
            Emety = METVect.Y();

            pbx = jetb.Px();
            pby = jetb.Py();
            pbz = jetb.Pz();
            pbbarx = jetbbar.Px();
            pbbary = jetbbar.Py();
            pbbarz = jetbbar.Pz();
            plpx = leptonPVect.Px();
            plpy = leptonPVect.Py();
            plpz = leptonPVect.Pz();
            plmx = leptonMVect.Px();
            plmy = leptonMVect.Py();
            plmz = leptonMVect.Pz();

            massb = jetb.M();
            massbbar = jetbbar.M();
            masslp = leptonPVect.M();
            masslm = leptonMVect.M();

            massWp = massW;
            massWm = massW;
            massnu = 0;
            massnubar = 0;

            TLorentzVector topCand;
            TLorentzVector tbarCand;

            for (double massTop = 100.; massTop <= 400.; massTop += 0.1)
            {
                double pnux, pnuy, pnuz;
				double pnubarx, pnubary, pnubarz;
				TLorentzVector pnu, pnubar;

                a[0] = (Eb + Elp)*(massW*massW - masslp*masslp - massnu*massnu)
						- Elp*(massTop*massTop - massb*massb - masslp*masslp - massnu*massnu)
						+ 2*Eb*Elp*Elp - 2*Elp*(pbx*plpx + pby*plpy + pbz*plpz);
                a[1] = 2*(Eb*plpx - Elp*pbx);
				a[2] = 2*(Eb*plpy - Elp*pby);
				a[3] = 2*(Eb*plpz - Elp*pbz);
				b[0] = (Ebbar + Elm)*(massW*massW - masslm*masslm - massnubar*massnubar)
                        - Elm*(massTop*massTop - massbbar*massbbar - masslm*masslm - massnubar*massnubar)
                        + 2*Ebbar*Elm*Elm - 2*Elm*(pbbarx*plmx + pbbary*plmy + pbbarz*plmz);
				b[1] = 2*(Ebbar*plmx - Elm*pbbarx);
				b[2] = 2*(Ebbar*plmy - Elm*pbbary);
				b[3] = 2*(Ebbar*plmz - Elm*pbbarz);
				c[2][2] = (massWp*massWp - masslp*masslp - massnu*massnu)*(massWp*massWp - masslp*masslp - massnu*massnu)
                        - 4*(Elp*Elp - plpz*plpz)*(a[0]*a[0]/a[3]/a[3])
                        - 4*(massWp*massWp - masslp*masslp - massnu*massnu)*plpz*a[0]/a[3];
				c[2][1] = 4*(massWp*massWp - masslp*masslp - massnu*massnu)*(plpx - plpz*a[1]/a[3])
                        -8*(Elp*Elp - plpz*plpz)*a[0]*a[1]/a[3]/a[3] - 8*(plpx*plpz*a[0]/a[3]);
				c[2][0] = -4*(Elp*Elp - plpx*plpx) - 4*(Elp*Elp - plpz*plpz)*a[1]*a[1]/a[3]/a[3]
                        -8*plpx*plpz*a[1]/a[3];
				c[1][1] = 4*(massWp*massWp - masslp*masslp - massnu*massnu)*(plpy - plpz*a[2]/a[3])
                        -8*(Elp*Elp - plpz*plpz)*a[0]*a[2]/a[3]/a[3] - 8*plpy*plpz*a[0]/a[3];
				c[1][0] = -8*(Elp*Elp - plpz*plpz)*a[1]*a[2]/a[3]/a[3] + 8*plpz*plpy
                        -8*plpx*plpz*a[2]/a[3] - 8*plpy*plpz*a[1]/a[3];
				c[0][0] = -4*(Elp*Elp - plpy*plpy) - 4*(Elp*Elp - plpz*plpz)*a[2]*a[2]/a[3]/a[3]
                        -8*plpy*plpz*a[2]/a[3];
				dp[2][2] = (massWm*massWm - masslm*masslm - massnubar*massnubar)*(massWm*massWm - masslm*masslm - massnubar*massnubar)
                        - 4*(Elm*Elm - plmz*plmz)*b[0]*b[0]/b[3]/b[3]
                        - 4*(massWm*massWm - masslm*masslm - massnubar*massnubar)*plmz*b[0]/b[3];
				dp[2][1] = 4*(massW*massW - masslm*masslm - massnubar*massnubar)*(plmx - plmz*b[1]/b[3])
                        - 8*(Elm*Elm - plmz*plmz)*b[0]*b[1]/b[3]/b[3] - 8*plmx*plmz*b[0]/b[3];
				dp[2][0] = -4*(Elm*Elm - plmx*plmx) - 4*(Elm*Elm - plmz*plmz)*b[1]*b[1]/b[3]/b[3]
                        - 8*plmx*plmz*b[1]/b[3];
				dp[1][1] = 4*(massWm*massWm - masslm*masslm - massnubar*massnubar)*(plmy - plmz*b[2]/b[3])
                        - 8*(Elm*Elm - plmz*plmz)*b[0]*b[2]/b[3]/b[3] - 8*plmy*plmz*b[0]/b[3];
				dp[1][0] = -8*(Elm*Elm - plmz*plmz)*b[1]*b[2]/b[3]/b[3] + 8*plmx*plmy
                        - 8*plmx*plmz*b[2]/b[3] - 8*plmy*plmz*b[1]/b[3];
				dp[0][0] = -4*(Elm*Elm - plmy*plmy) - 4*(Elm*Elm - plmz*plmz)*b[2]*b[2]/b[3]/b[3]
                        - 8*plmy*plmz*b[2]/b[3];
				d[2][2] = dp[2][2] + Emetx*Emetx*dp[2][0] + Emety*Emety*dp[0][0] +Emetx*Emety*dp[1][0]
                        + Emetx*dp[2][1] + Emety*dp[1][1];
				d[2][1] = -dp[2][1] - 2*Emetx*dp[2][0] - Emety*dp[1][0];
				d[2][0] = dp[2][0];
				d[1][1] = -dp[1][1] - 2*Emety*dp[0][0] - Emetx*dp[1][0];
				d[1][0] = dp[1][0];
				d[0][0] = dp[0][0];

                h[4] = c[0][0]*c[0][0]*d[2][2]*d[2][2]
                        + c[1][1]*d[2][2]*(c[1][1]*d[0][0] - c[0][0]*d[1][1])
                        + c[0][0]*c[2][2]*(d[1][1]*d[1][1] - 2*d[0][0]*d[2][2])
                        + c[2][2]*d[0][0]*(c[2][2]*d[0][0] - c[1][1]*d[1][1]);
				h[3] = c[0][0]*d[2][1]*(2*c[0][0]*d[2][2] - c[1][1]*d[1][1])
                        + c[0][0]*d[1][1]*(2*c[2][2]*d[1][0] + c[2][1]*d[1][1])
                        + c[2][2]*d[0][0]*(2*c[2][1]*d[0][0] - c[1][1]*d[1][0])
                        - c[0][0]*d[2][2]*(c[1][1]*d[1][0] + c[1][0]*d[1][1])
                        - 2*c[0][0]*d[0][0]*(c[2][2]*d[2][1] + c[2][1]*d[2][2])
                        - d[0][0]*d[1][1]*(c[1][1]*c[2][1] + c[1][0]*c[2][2])
                        + c[1][1]*d[0][0]*(c[1][1]*d[2][1] + 2*c[1][0]*d[2][2]);
				h[2] = c[0][0]*c[0][0]*(2*d[2][2]*d[2][0] + d[2][1]*d[2][1])
                        - c[0][0]*d[2][1]*(c[1][1]*d[1][0] + c[1][0]*d[1][1])
                        + c[1][1]*d[2][0]*(c[1][1]*d[0][0] - c[0][0]*d[1][1])
                        + c[0][0]*d[1][0]*(c[2][2]*d[1][0] - c[1][0]*d[2][2])
                        + c[0][0]*d[1][1]*(2*c[2][1]*d[1][0] + c[2][0]*d[1][1])
                        + (2*c[2][2]*c[2][0] + c[2][1]*c[2][1])*d[0][0]*d[0][0]
                        - 2*c[0][0]*d[0][0]*(c[2][2]*d[2][0] + c[2][1]*d[2][1] + c[2][0]*d[2][2])
                        + c[1][0]*d[0][0]*(2*c[1][1]*d[2][1] + c[1][0]*d[2][2])
                        - d[0][0]*d[1][0]*(c[1][1]*c[2][1] + c[1][0]*c[2][2])
                        - d[0][0]*d[1][1]*(c[1][1]*c[2][0] + c[1][0]*c[2][1]);
				h[1] = c[0][0]*d[2][1]*(2*c[0][0]*d[2][0] - c[1][0]*d[1][0])
                        - c[0][0]*d[2][0]*(c[1][1]*d[1][0] + c[1][0]*d[1][1])
                        + c[0][0]*d[1][0]*(c[2][1]*d[1][0] + 2*c[2][0]*d[1][1])
                        - 2*c[0][0]*d[0][0]*(c[2][1]*d[2][0] + c[2][0]*d[2][1])
                        + c[1][0]*d[0][0]*(2*c[1][1]*d[2][0] + c[1][0]*d[2][1])
                        - c[2][0]*d[0][0]*(2*c[2][1]*d[0][0] - c[1][0]*d[1][1])
                        - d[0][0]*d[1][0]*(c[1][1]*c[2][0] + c[1][0]*c[2][1]);
				h[0] = c[0][0]*c[0][0]*d[2][0]*d[2][0]
                        + c[1][0]*d[2][0]*(c[1][0]*d[0][0] - c[0][0]*d[1][0])
                        + c[2][0]*d[1][0]*(c[0][0]*d[1][0] - c[1][0]*d[0][0])
                        + c[2][0]*d[0][0]*(c[2][0]*d[0][0] - 2*c[0][0]*d[2][0]);

				vector<double> pnuxSol;

				//! Get the solution by using ROOT::Math::Polynomial object
				//! and ROOT::Math::Polynomial.FindRealRoots() here.

				ROOT::Math::Polynomial pnuxPoly(h[0],h[1],h[2],h[3],h[4]);
				pnuxSol = pnuxPoly.FindRealRoots();
				//printf("pnuxSol has %d solutions\n",(int)pnuxSol.size());

				//for (int i=0; i < (int) pnuxSol.size(); i++) printf("%lf\t",(double) pnuxSol.at(i));
				//printf("\n");

				//solutionCount += (int) pnuxSol.size();
				//histpnuxsolCount->Fill((int)pnuxSol.size());

                for (int solCount = 0; solCount < (int) pnuxSol.size(); solCount++)
				{
					bool uniqueSol = true;
					for (int i=0;i<solCount;i++) if ((double)pnuxSol.at(i) == (double)pnuxSol.at(solCount)) uniqueSol = false;
					if (!uniqueSol) continue;

					pnux    = (double) pnuxSol.at(solCount);
					pnubarx = Emetx - pnux;

					double c0, c1, c2;
					double d0, d1, d2;

					c0 = c[0][0];
					c1 = c[1][0]*pnux + c[1][1];
					c2 = c[2][0]*pnux*pnux + c[2][1]*pnux + c[2][2];

					d0 = d[0][0];
					d1 = d[1][0]*pnux +d[1][1];
					d2 = d[2][0]*pnux*pnux + d[2][1]*pnux + d[2][2];

					pnuy    = (c0*d2 - c2*d0)/(c1*d0 - c0*d1);
					pnubary = Emety - pnuy;

					pnuz    = -(a[0] + a[1]*pnux + a[2]*pnuy)/a[3];
					pnubarz = -(b[0] + b[1]*pnubarx + b[2]*pnubary)/b[3];

					pnu.SetPxPyPzE(pnux,pnuy,pnuz,TMath::Sqrt(pnux*pnux + pnuy*pnuy + pnuz*pnuz));
					pnubar.SetPxPyPzE(pnubarx,pnubary,pnubarz,TMath::Sqrt(pnubarx*pnubarx + pnubary*pnubary + pnubarz*pnubarz));

					double Delta = 0;
                    double DeltaTop = 0;
					/*
					for (int i=0;i<4;i++) Delta += (pnu[i] - pnubar[i])*(pnu[i] - pnubar[i]);
					Delta = TMath::Sqrt(Delta);
					histDeltaDistrib->Fill(Delta);
					*/
					Delta = (pnu.Vect() + pnubar.Vect() - METVect.Vect()).Mag();
					//printf("Delta = %lf\n",Delta);
					double WLeptonRatio1 = ((pnu + leptonPVect).M())/massW;
					double WLeptonRatio2 = ((pnubar + leptonMVect).M())/massW;

                    TLorentzVector topVect, tbarVect;
                    topVect = leptonPVect + jetb + pnu;
                    tbarVect = leptonMVect + jetbbar + pnubar;

                    DeltaTop = TMath::Abs(massTop - topVect.M());
                    DeltaTop = TMath::Abs(massTop - tbarVect.M()) > DeltaTop ? TMath::Abs(massTop - tbarVect.M()) : DeltaTop;

                    bool DeltaIsLower;
                    if (!offShell) DeltaIsLower = Delta*DeltaTop < bestDelta*bestDeltaTop;
                    else DeltaIsLower = Delta < bestDelta;

					if (bestDelta < 0. || DeltaIsLower)
					{
						bestDelta = Delta;
                        if (!offShell) bestDeltaTop = DeltaTop;
						bestMassTop = massTop;
						deltaRBestMassTop1 = leptonPVect.DeltaR(pnu);
						deltaRBestMassTop2 = leptonMVect.DeltaR(pnubar);
                        topCand = leptonPVect + jetb + pnu;
                        tbarCand = leptonMVect + jetbbar + pnubar;
					}
				}
            }

            /*
            *topMassPt[ind] = bestMassTop;

            if (bestMassTop > 0.)
            {
                *topPtPt[ind]  = topCand.Pt();
                *topEtaPt[ind] = topCand.Eta();
                *topPhiPt[ind] = topCand.Phi();
                *topEPt[ind]   = topCand.E();

                *tbarPtPt[ind]  = tbarCand.Pt();
                *tbarEtaPt[ind] = tbarCand.Eta();
                *tbarPhiPt[ind] = tbarCand.Phi();
                *tbarEPt[ind]   = tbarCand.E();
            }
            */
            
            top_mass = bestMassTop;
            if (bestMassTop > 0.)
            {
                top_pt  = topCand.P();
                top_eta = topCand.Eta();
                top_phi = topCand.Phi();
                top_e   = topCand.E();
                tbar_pt  = tbarCand.P();
                tbar_eta = tbarCand.Eta();
                tbar_phi = tbarCand.Phi();
                tbar_e   = tbarCand.E();
            }
        }

        //if (!(top1_mass > 0 && top2_mass > 0)) continue;
        if (top_mass < 0) continue;

        outTree.Fill();
    }

    outFile.cd();
    outTree.Write();
    printf("All entries in file: %lld\n",nentries);
    printf("Valid entries with opposite charge leptons: %d\n", validLEntries);
    printf("Output tree contains %lld events.\n",outTree.GetEntries());
    //outFile.Close();
    //file->Close();
}
