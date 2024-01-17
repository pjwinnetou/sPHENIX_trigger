#define triggerEff_cxx
#include "triggerEff.h"
#include "TEfficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include "set"

using namespace std;

void triggerEff::Loop()
{
//   In a ROOT session, you can do:
//      root> .L triggerEff.C
//      root> triggerEff t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  ROOT::EnableImplicitMT();
  TFile *f = new TFile(Form("OutEff_seed%d.root",fseed),"recreate");
  TH1D* hden = new TH1D("hden",";Truth Energy;All",nBins,0,15);
  
  map<TString,TH1D*> hnum_th5;
  map<TString,TH1D*> hnum_th4;
  map<TString,TH1D*> hnum_th3;
  map<TString,TriggerResult> result_th5;
  map<TString,TriggerResult> result_th4;
  map<TString,TriggerResult> result_th3;
  for(int i=0; i<NCASES; i++){
    hnum_th5[i] = new TH1D(Form("hnum_th5_%s",histname[i].c_str()),";Truth Energy;Pass Counts",nBins,0,15);
    hnum_th4[i] = new TH1D(Form("hnum_th4_%s",histname[i].c_str()),";Truth Energy;Pass Counts",nBins,0,15);
    hnum_th3[i] = new TH1D(Form("hnum_th3_%s",histname[i].c_str()),";Truth Energy;Pass Counts",nBins,0,15);
  }
  std::map<int, std::pair<int, int>> binEfficiencies;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    if(jentry%10==0) std::cout << "event : " << jentry << " / " << nentries << " (" << (double)jentry/nentries*100. << "%)" << std::endl;
    fChain->GetEntry(jentry);
    
    if(abs((*truth_eta)[0])>1.1) continue;
    hden->Fill((*truth_energy)[0]);

    for(int ir=0; ir<NCASES;ir++){
      result_th5[ir] = AnalyzeEvent(emcal_etabin, emcal_phibin, emcal_energy,5,cl_comb[ir].first,cl_comb[ir].second);
      result_th4[ir] = AnalyzeEvent(emcal_etabin, emcal_phibin, emcal_energy,4,cl_comb[ir].first,cl_comb[ir].second);
      result_th3[ir] = AnalyzeEvent(emcal_etabin, emcal_phibin, emcal_energy,3,cl_comb[ir].first,cl_comb[ir].second);

      if(result_th5[ir].numClustersAboveThreshold > 0) hnum_th5[ir]->Fill((*truth_energy)[0]);
      if(result_th4[ir].numClustersAboveThreshold > 0) hnum_th4[ir]->Fill((*truth_energy)[0]);
      if(result_th3[ir].numClustersAboveThreshold > 0) hnum_th3[ir]->Fill((*truth_energy)[0]);
    }

    //if(result.numClustersAboveThreshold > 0){eff->Fill(true, (*truth_energy)[0]);}
    //else{eff->Fill(false, (*truth_energy)[0]);}
  }
  f->cd();
  hden->Write();
  for(int i=0; i<NCASES; i++){
    hnum_th5[i]->Write();
    hnum_th4[i]->Write();
    hnum_th3[i]->Write();
  }
  f->Close();
}

triggerEff::TriggerResult triggerEff::AnalyzeEvent(vector<int> *eta_bin,vector<int> *phi_bin, vector<float> *energy, double threshold, int clusterSize, bool allowOverlap) {
  TriggerResult result = {0, 0.0, -1, -1};
  std::set<std::pair<int, int>> usedTowers;

  int nClusters = 23529;
  if(allowOverlap){
    if(clusterSize==2) nClusters = 24225;
    else if(clusterSize==4) nClusters = 23529;
    else if(clusterSize==8) nClusters = 22161;
    else{ std::cerr << "ERROR incorrect cluster size"; exit(1);}
  }
  else if(!allowOverlap){
    if(clusterSize==2) nClusters = 6144;
    else if(clusterSize==4) nClusters = 1536;
    else if(clusterSize==8) nClusters = 384;
    else{ std::cerr << "ERROR incorrect cluster size"; exit(1);}
  }

  int stepSize = allowOverlap ? 1 : clusterSize;

  for (int startEta = 0; startEta <= GRID_ETA - clusterSize; startEta += stepSize) {
    for (int startPhi = 0; startPhi <= GRID_PHI - clusterSize; startPhi += stepSize) {
      double clusterEnergy = 0.0;
      double maxEnergyInCluster = 0.0;
      int maxEnergyEta = -1;
      int maxEnergyPhi = -1;

      for (int dEta = 0; dEta < clusterSize; ++dEta) {
        for (int dPhi = 0; dPhi < clusterSize; ++dPhi) {
          int eta = startEta + dEta;
          int phi = startPhi + dPhi;

          double energyAtBin = getEnergy(eta, phi, eta_bin, phi_bin, energy);
          clusterEnergy += energyAtBin;

          if (energyAtBin > maxEnergyInCluster) {
            maxEnergyInCluster = energyAtBin;
            maxEnergyEta = eta;
            maxEnergyPhi = phi;
          }
        }
      }

      if (clusterEnergy > threshold) {
        result.numClustersAboveThreshold++;
        if (maxEnergyInCluster > result.maxEnergyInTriggeredCluster) {
          result.maxEnergyInTriggeredCluster = maxEnergyInCluster;
          result.maxEnergyEtaBin = maxEnergyEta;
          result.maxEnergyPhiBin = maxEnergyPhi;
        }
      }
    }
  }

  return result;
}

