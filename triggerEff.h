//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 16 14:26:09 2024 by ROOT version 6.26/06
// from TTree towerntup/Towers
// found on file: outputfiles/G4sPHENIX_merged.root
//////////////////////////////////////////////////////////

#ifndef triggerEff_h
#define triggerEff_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "set"

// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;

class triggerEff {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   Int_t           fseed;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         vz;
   vector<float>   *emcal_energy;
   vector<int>     *emcal_etabin;
   vector<int>     *emcal_phibin;
   vector<float>   *emcal_etacorr;
   Int_t           truthpar_n;
   vector<float>   *truth_energy;
   vector<float>   *truth_pt;
   vector<float>   *truth_eta;
   vector<float>   *truth_phi;
   vector<int>     *truth_id;
   Int_t           npart;
   Int_t           ncoll;
   Float_t         bimp;
   Float_t         Cent_impactparam;
   Int_t           Centrality;

   // List of branches
   TBranch        *b_vz;   //!
   TBranch        *b_emcal_energy;   //!
   TBranch        *b_emcal_etabin;   //!
   TBranch        *b_emcal_phibin;   //!
   TBranch        *b_emcal_etacorr;   //!
   TBranch        *b_truthpar_n;   //!
   TBranch        *b_truth_energy;   //!
   TBranch        *b_truth_pt;   //!
   TBranch        *b_truth_eta;   //!
   TBranch        *b_truth_phi;   //!
   TBranch        *b_truth_id;   //!
   TBranch        *b_npart;   //!
   TBranch        *b_ncoll;   //!
   TBranch        *b_bimp;   //!
   TBranch        *b_Cent_impactparam;   //!
   TBranch        *b_Centrality;   //!

   triggerEff(TTree *tree=0, int seed=1);
   virtual ~triggerEff();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree,int seed);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Double_t getEnergy(int eta, int phi, vector<int> *eta_bin, vector<int> *phi_bin, vector<float> *energy);

   static const int GRID_ETA = 96;
   static const int GRID_PHI = 256;
   static const int SIZE = GRID_ETA * GRID_PHI;
   static const int nBins = 50;

   static const int NCASES = 3;
   string histname[NCASES] = {"4ol", "4nol", "2nol"};
   std::pair<int,int> cl_comb[NCASES] = {{4,1},{4,0},{2,0}};

   struct TriggerResult {
     int numClustersAboveThreshold;
     double maxEnergyInTriggeredCluster;
     int maxEnergyEtaBin;
     int maxEnergyPhiBin;
   };

   virtual TriggerResult AnalyzeEvent(vector<int> *eta_bin, vector<int> *phi_bin, vector<float> *energy, double threshold, int clusterSize, bool allowOverlap);


};

#endif

#ifdef triggerEff_cxx
triggerEff::triggerEff(TTree *tree, int seed) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("../outputfiles/G4sPHENIX_seed%d.root",seed));
      if (!f || !f->IsOpen()) {
         f = new TFile(Form("../outputfiles/G4sPHENIX_seed%d.root",seed));
      }
      f->GetObject("towerntup",tree);

   }
   Init(tree,seed);
}

triggerEff::~triggerEff()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t triggerEff::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t triggerEff::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void triggerEff::Init(TTree *tree, int seed)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   emcal_energy = 0;
   emcal_etabin = 0;
   emcal_phibin = 0;
   emcal_etacorr = 0;
   truth_energy = 0;
   truth_pt = 0;
   truth_eta = 0;
   truth_phi = 0;
   truth_id = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fseed  = seed;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("emcal_energy", &emcal_energy, &b_emcal_energy);
   fChain->SetBranchAddress("emcal_etabin", &emcal_etabin, &b_emcal_etabin);
   fChain->SetBranchAddress("emcal_phibin", &emcal_phibin, &b_emcal_phibin);
   fChain->SetBranchAddress("emcal_etacorr", &emcal_etacorr, &b_emcal_etacorr);
   fChain->SetBranchAddress("truthpar_n", &truthpar_n, &b_truthpar_n);
   fChain->SetBranchAddress("truth_energy", &truth_energy, &b_truth_energy);
   fChain->SetBranchAddress("truth_pt", &truth_pt, &b_truth_pt);
   fChain->SetBranchAddress("truth_eta", &truth_eta, &b_truth_eta);
   fChain->SetBranchAddress("truth_phi", &truth_phi, &b_truth_phi);
   fChain->SetBranchAddress("truth_id", &truth_id, &b_truth_id);
   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("ncoll", &ncoll, &b_ncoll);
   fChain->SetBranchAddress("bimp", &bimp, &b_bimp);
   fChain->SetBranchAddress("Cent_impactparam", &Cent_impactparam, &b_Cent_impactparam);
   fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
   Notify();
}

Bool_t triggerEff::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void triggerEff::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t triggerEff::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Double_t triggerEff::getEnergy(int eta, int phi, vector<int> *eta_bin, vector<int> *phi_bin, vector<float> *energy) {
  for (int i = 0; i < SIZE; ++i) {
    if ((*eta_bin)[i] == eta && (*phi_bin)[i] == phi) {
      return (*energy)[i];
    }
  }
  return 0;
}

#endif // #ifdef triggerEff_cxx
