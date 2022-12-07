#include <iostream>
#include <array>
#include "Riostream.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include <fstream>
#include "TGraph.h"
#include "TGraphErrors.h"



//this is just a place to test our cpp code without the clutter of other stuff going on

void code_environment(){


  //in the original file numdets = 18, numE=10, numcrystals=72
  //changed numcrystals to 4 from 64 for testing purposes
  static const int numdets = 16;
  static const int numE = 12;
  static const int numcrystals = 64;

  //Energies
  int Energy[numE];  

  Energy[0] = 43;  //in keV
  Energy[1] = 87; //in keV
  Energy[2] = 105; //in keV
  Energy[3] = 123; //in keV
  //Energy[4] = 176; //in keV
  Energy[4] = 247; //in keV
  //Energy[6] = 428; //in keV
  //Energy[7] = 463; //in keV
  Energy[5] = 591; //in keV
  //Energy[9] = 601; //in keV
  //Energy[10] = 636; //in keV
  Energy[6] = 723;  //in keV
  Energy[7] = 873; //in keV
  Energy[8] = 996; //in keV
  Energy[9] = 1004; //in keV
  Energy[10] = 1274; //in keV
  Energy[11] = 1596; //in keV



  //Number of Gammas in Simulation
  //changing this nGamma value for testing purposes on Jan 26 2021. 
  //was originally 10000000
  double nGamma = 10000000;

  // Efficiencies [Energy][Det Number]
  double Peak_Efficiency[numE][numdets] = {};
  double Peak_Efficiency_Err[numE][numdets] = {};

  double Peak_Efficiency_crystal[numE][numcrystals] = {};
  double Peak_Efficiency_Err_crystal[numE][numcrystals] = {};
  
  double Total_Efficiency[numE][numdets] = {};
  double Total_Efficiency_Err[numE][numdets] = {};

  double Total_Efficiency_crystal[numE][numcrystals] = {};
  double Total_Efficiency_Err_crystal[numE][numcrystals] = {};

  //FIles
  TFile *f1t[numE];
  
  for(int i=0; i<numE; i++) {
    //f1t[i] = new TFile(Form("/data/e17011/simulations/analysis/e17011_Eff_%dkeV_4p8in_0.root",Energy[i]));
    //f1t[i] = new TFile(Form("/data/e17011/simulations/analysis/e17011_Eff_%dkeV_5p5in_0.root",Energy[i]));

    //changing this file path from /data/e17011/simulation/analysis/e17011_Eff_%dkeV_test_0.root on Jan 26 2021
    f1t[i] = new TFile(Form("/home/dcs411/e17011/Efficiency/sim_analysis/simrootfiles/%dkeVtest_0.root",Energy[i]));
  }

  //Fits & Histos
  TF1 *fClover[numE][numdets];
  TH1D *hClover[numE][numdets];
  TH1D *hClover_crystal[numE][numcrystals];
  TF1 *fClover_crystal[numE][numcrystals];
  
  for(int i=0; i<numE; i++) {
    for(int j=0; j<numdets; j++) {
      hClover[i][j] = (TH1D*) f1t[i]->Get(Form("hEnergyDepositClover_addback_%d",j));

      if(Energy[i]<87) {
  	fClover[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus",0.80*Energy[i],0.95*Energy[i]);
      }
      else if(Energy[i]<100) {
  	fClover[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus",0.95*Energy[i],1.10*Energy[i]);
      }
      else if(Energy[i]<200) {
  	fClover[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus",0.95*Energy[i],1.05*Energy[i]);
      }
      else {
  	fClover[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus",0.98*Energy[i],1.02*Energy[i]);
      }
    }
  }

  for(int i=0; i<numE; i++) {
    for(int j=0; j<numcrystals; j++) {
      hClover_crystal[i][j] = (TH1D*) f1t[i]->Get(Form("hEnergyDepositClover_%d",j));
      
      if(Energy[i]<100) {
	fClover_crystal[i][j] = new TF1(Form("fc_%dkev_%d",Energy[i],j),"gaus",0.90*Energy[i],1.10*Energy[i]);
      }
      else if(Energy[i]<200) {
	fClover_crystal[i][j] = new TF1(Form("fc%dkev_%d",Energy[i],j),"gaus",0.95*Energy[i],1.05*Energy[i]);
      }
      else {
	fClover_crystal[i][j] = new TF1(Form("fc%dkev_%d",Energy[i],j),"gaus",0.98*Energy[i],1.02*Energy[i]);
      }
    }
  }
  
  
  //TCanvas
  // TCanvas *c1t[numE];
  
  // for(int i=0; i<numE; i++) {
  //   c1t[i] = new TCanvas(Form("c_%d_keV",Energy[i]),Form("c_%d_keV",Energy[i]),1024,768);
  //   c1t[i]->Divide(4,5);
  // }

  //output root file
  //changing file path from /data/e17011/simulations/analysis/e17011_Eff_out.root Jan 26 2021
  TFile *fout = new TFile("/home/dcs411/e17011/Efficiency/sim_analysis/simrootfiles/test_0.root","RECREATE");

  for(int i=0; i<numE; i++) {
    for(int j=0; j<numdets; j++) {

      // c1t[i]->cd(j+1);
      hClover[i][j]->GetXaxis()->SetRangeUser(0.9*Energy[i],1.10*Energy[i]);
      hClover[i][j]->Draw();
    
      fClover[i][j]->SetParLimits(0,1,1000000);
      fClover[i][j]->SetParameter(0,100);
      fClover[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      fClover[i][j]->SetParameter(1,Energy[i]);   
      fClover[i][j]->SetParLimits(2,0,5);
      fClover[i][j]->SetParameter(2,5);
      
      hClover[i][j]->Fit(fClover[i][j],"RQ");
      
      fClover[i][j]->SetParLimits(0,1,1000000);
      fClover[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      fClover[i][j]->SetParLimits(2,0,5);
      
      hClover[i][j]->Fit(fClover[i][j],"RQ");
      
      fClover[i][j]->SetParLimits(0,1,1000000);
      fClover[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      fClover[i][j]->SetParLimits(2,0,5);
      
      hClover[i][j]->Fit(fClover[i][j],"RQ");
      
      fClover[i][j]->SetParLimits(0,1,1000000);
      fClover[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      fClover[i][j]->SetParLimits(2,0,5);
      
      hClover[i][j]->Fit(fClover[i][j],"RQ+");

      double area = 2.506628*(fClover[i][j]->GetParameter(0))*(fClover[i][j]->GetParameter(2));
      double area_err = area*pow((pow((fClover[i][j]->GetParError(0)/fClover[i][j]->GetParameter(0)),2)+pow((fClover[i][j]->GetParError(2)/fClover[i][j]->GetParameter(2)),2)),0.5);
      double eff = 100.0*area/(1.0*nGamma);
      double eff_err = 100.0*area_err/(1.0*nGamma);

      double total_eff = 0;
      for(int k=2; k<Energy[i]*1.10; k++) {
  	total_eff += hClover[i][j]->GetBinContent(k);
      }
      double total_eff_err = pow(total_eff,0.5)/(1.0*nGamma);
      total_eff = total_eff/(1.0*nGamma);
      

      //cout<<i<<"  "<<j<<"  "<<eff<<" +/- "<<eff_err<<"  "<<100*total_eff<<" +/- "<<100*total_eff_err<<endl;

  //     Efficiencies [Energy][Det Number]
      // if(j>=8 && (Energy[i]==42 || Energy[i] == 86 || Energy[i] ==105)) {
      if(j<0){
  	cout<<"Special Condition"<<endl;
  	Peak_Efficiency[i][j] = 0;
  	Peak_Efficiency_Err[i][j] = 0;
	
  	Total_Efficiency[i][j] = 0;
  	Total_Efficiency_Err[i][j] = 0;
	
      }
      else {
  	Peak_Efficiency[i][j] = eff;
  	Peak_Efficiency_Err[i][j] = eff_err;
	
  	Total_Efficiency[i][j] = total_eff;
  	Total_Efficiency_Err[i][j] = total_eff_err;
	
      }
    }
  }

  for(int i=0; i<numE; i++) {
    for(int j=0; j<numcrystals; j++) {

      //c1t[i]->cd(j+1);
      hClover_crystal[i][j]->GetXaxis()->SetRangeUser(0.9*Energy[i],1.10*Energy[i]);
      //hClover_crystal[i][j]->Draw();
    
      fClover_crystal[i][j]->SetParLimits(0,1,1000000);
      fClover_crystal[i][j]->SetParameter(0,100);
      fClover_crystal[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      fClover_crystal[i][j]->SetParameter(1,Energy[i]);   
      fClover_crystal[i][j]->SetParLimits(2,0,5);
      fClover_crystal[i][j]->SetParameter(2,5);
      
      hClover_crystal[i][j]->Fit(fClover_crystal[i][j],"RQ");
      
      fClover_crystal[i][j]->SetParLimits(0,1,1000000);
      fClover_crystal[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      fClover_crystal[i][j]->SetParLimits(2,0,5);
      
      hClover_crystal[i][j]->Fit(fClover_crystal[i][j],"RQ");
      
      fClover_crystal[i][j]->SetParLimits(0,1,1000000);
      fClover_crystal[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      fClover_crystal[i][j]->SetParLimits(2,0,5);
      
      hClover_crystal[i][j]->Fit(fClover_crystal[i][j],"RQ");
      
      fClover_crystal[i][j]->SetParLimits(0,1,1000000);
      fClover_crystal[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      fClover_crystal[i][j]->SetParLimits(2,0,5);
      
      hClover_crystal[i][j]->Fit(fClover_crystal[i][j],"RQ+");

      double area = 2.506628*(fClover_crystal[i][j]->GetParameter(0))*(fClover_crystal[i][j]->GetParameter(2));
      double area_err = area*pow((pow((fClover_crystal[i][j]->GetParError(0)/fClover_crystal[i][j]->GetParameter(0)),2)+pow((fClover_crystal[i][j]->GetParError(2)/fClover_crystal[i][j]->GetParameter(2)),2)),0.5);
      double eff = 100.0*area/(1.0*nGamma);
      double eff_err = 100.0*area_err/(1.0*nGamma);

      double total_eff = 0;
      for(int k=2; k<Energy[i]*1.10; k++) {
	total_eff += hClover_crystal[i][j]->GetBinContent(k);
      }
      double total_eff_err = pow(total_eff,0.5)/(1.0*nGamma);
      total_eff = total_eff/(1.0*nGamma);
      

      //cout<<i<<"  "<<j<<"  "<<eff<<" +/- "<<eff_err<<"  "<<100*total_eff<<" +/- "<<100*total_eff_err<<endl;

      // Efficiencies [Energy][Det Number]
      //need to update this for clover. This is prob for SeGa
      // if(j>=8 && (Energy[i]==42 || Energy[i] == 86 || Energy[i] ==105)) {
      if (j<0){
	cout<<"Special Condition"<<endl;
	Peak_Efficiency_crystal[i][j] = 0;
	Peak_Efficiency_Err_crystal[i][j] = 0;
	
	Total_Efficiency_crystal[i][j] = 0;
	Total_Efficiency_Err_crystal[i][j] = 0;
	
      }
      else {
	Peak_Efficiency_crystal[i][j] = eff;
	Peak_Efficiency_Err_crystal[i][j] = eff_err;
	
	Total_Efficiency_crystal[i][j] = total_eff;
	Total_Efficiency_Err_crystal[i][j] = total_eff_err;
      }
    }
  }

  // for(int i = 0; i < numE; i++){
  //   for(int j = 0; j < numcrystals; j++){
  //     if(Peak_Efficiency_Err_crystal[i][j] > 0.003){
  // 	cout << i << " " << j << " " << Peak_Efficiency_Err_crystal[i][j] << endl;}
  //   }
  // }
  
  ofstream data;
  data.open("Sim_Eff.txt");
  
  for(int j=0; j<numdets; j++) {
    data<<"Clover"<<j<<"\n";
    for(int i=0; i<numE; i++) {
      data<<Peak_Efficiency[i][j]<<"  "<< Peak_Efficiency_Err[i][j]<<"  "<<Total_Efficiency[i][j]<<"  "<<Total_Efficiency_Err[i][j]<<"\n";	  
    }
    data<<"\n";	
  }
  for(int j=0; j<numcrystals;j++){
    data<<"Crystal"<<j<<"\n";
    for(int i=0; i<numE; i++) {
      data<<Peak_Efficiency_crystal[i][j]<<"  "<< Peak_Efficiency_Err_crystal[i][j]<<"  "<<Total_Efficiency_crystal[i][j]<<"  "<<Total_Efficiency_Err_crystal[i][j]<<"\n";	  
    }
    data<<"\n";	
  }

  //make graph of efficiency curve
  TCanvas *c1graph1 = new TCanvas("addback","addback",800,600);
  TCanvas *c1graph2 = new TCanvas("total","total",800,600);
  TCanvas *c1graph3 = new TCanvas("upstream","upstream",800,600);
  TCanvas *c1graph4 = new TCanvas("downstream","downstream",800,600);
  TCanvas *c1graph5 = new TCanvas("above","above",800,600);
  TCanvas *c1graph6 = new TCanvas("below","below",800,600);
  TCanvas *c1graph7 = new TCanvas("left","left",800,600);
  TCanvas *c1graph8 = new TCanvas("right","right",800,600);
  Double_t sum_eff[numE] = {};
  Double_t sum_eff_err[numE]= {};
  Double_t sum_eff_crystal[numE] = {};
  Double_t sum_eff_crystal_err[numE] = {};
  Double_t err_temp1[numE] = {};
  Double_t err_temp2[numE] = {};


  for(int j = 0; j < numdets; j++){
    for(int i = 0; i < numE; i++){
      sum_eff[i] += Peak_Efficiency[i][j];
      err_temp1[i] += pow(Peak_Efficiency_Err[i][j],2);
    }
  }


  for(int j = 0; j < numcrystals; j++){
    for(int i = 0; i < numE; i++){
      sum_eff_crystal[i] += Peak_Efficiency_crystal[i][j];
      err_temp2[i] += pow(Peak_Efficiency_Err_crystal[i][j],2);
      //cout << err_temp2[i] << endl;
    }
  }


  for(int i = 0; i < numE; i++){
    if(i<=1){
      sum_eff_err[i] = 0;
      sum_eff_crystal_err[i] = 0;
    }
    else{
    sum_eff_err[i] = err_temp1[i]*sqrt(err_temp1[i]);
    sum_eff_crystal_err[i] = err_temp2[i]*sqrt(err_temp2[i]);
    }
    //cout << sum_eff_crystal_err[i] << endl;
  }

  //declaring crystal number vectors for regions of interest in the simulation
  vector<int> simupcrystals = {8,9,10,11,16,17,18,19,32,33,34,35,36,37,38,39,4,5,6,7,48,49,50,51,28,29,30,31,52,53,54,55,16,17,18,19,56,57,58,59,40,41,42,43,60,61,62,63};
  vector<int> simdowncrystals = {0,1,2,3,24,25,26,27,28,29,30,31,48,49,50,51};
  vector<int> simleftcrystals = {0,1,2,3,4,5,6,7,8,9,10,11,48,49,50,51,60,61,62,63};
  vector<int> simrightcrystals = {12,13,14,15,16,17,18,19,20,21,22,23,52,53,54,55,56,57,58,59};
  vector<int> simabovecrystals = {24,25,26,27,28,29,30,31,32,33,34,35,48,49,50,51,52,53,54,55};
  vector<int> simbelowcrystals = {36,37,38,39,40,41,42,43,44,45,46,47,56,57,58,59,60,61,62,63};

  double simupeff[numE] = {};
  double simupefferr[numE] = {};
  double simdowneff[numE] = {};
  double simdownefferr[numE] = {};
  double simaboveeff[numE] = {};
  double simaboveefferr[numE] = {};
  double simbeloweff[numE] = {};
  double simbelowefferr[numE] = {};
  double simlefteff[numE] = {};
  double simleftefferr[numE] = {};
  double simrighteff[numE] = {};
  double simrightefferr[numE] = {};

  //populating region vectors
    //looping through the crystals
  for(int i = 0; i < numcrystals; i++){
    //looping through energies
    for(int j = 0; j < numE; j++){
      //checking to see if the crystal is downstream
      if(std::find(simdowncrystals.begin(), simdowncrystals.end(), i) != simdowncrystals.end()) {
	//populating downstream efficiency and error vectors
	simdowneff[j] += Peak_Efficiency_crystal[j][i];
	simdownefferr[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is upstream
      if(std::find(simupcrystals.begin(), simupcrystals.end(), i) != simupcrystals.end()) {
	//populating upstream efficiency and error vectors
        simupeff[j] += Peak_Efficiency_crystal[j][i];
        simupefferr[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is on the right
      if(std::find(simrightcrystals.begin(), simrightcrystals.end(), i) != simrightcrystals.end()) {
	//populating right efficiency and error vectors
	simrighteff[j] += Peak_Efficiency_crystal[j][i];
	simrightefferr[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is on the left
      if(std::find(simleftcrystals.begin(), simleftcrystals.end(), i) != simleftcrystals.end()) {
	//populating left efficiency and error vectors
	simlefteff[j] += Peak_Efficiency_crystal[j][i];
	simleftefferr[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is above
      if(std::find(simabovecrystals.begin(), simabovecrystals.end(), i) != simabovecrystals.end()) {
	//populating above efficiency and error vectors
	simaboveeff[j] += Peak_Efficiency_crystal[j][i];
	simaboveefferr[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is below
      if(std::find(simbelowcrystals.begin(), simbelowcrystals.end(), i) != simbelowcrystals.end()) {
	//populating below efficiency and error vectors
	simbeloweff[j] += Peak_Efficiency_crystal[j][i];
	simbelowefferr[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }

    }

  }

  //calculating errors
  //setting 43keV energy err to 0 right now because not enough statistics and throws off the graph
    for(int i = 0; i<numE; i++){
      if(i<1){
	simbelowefferr[i] = 0;
	simaboveefferr[i] = 0;
	simdownefferr[i] = 0;
	simupefferr[i] = 0;
	simleftefferr[i] = 0;
	simrightefferr[i] = 0;
      }
      else{
	simbelowefferr[i] = sqrt(simbelowefferr[i]);
	simaboveefferr[i] = sqrt(simaboveefferr[i]);
	simdownefferr[i] = sqrt(simdownefferr[i]);
	simupefferr[i] = sqrt(simupefferr[i]);
	simleftefferr[i] = sqrt(simleftefferr[i]);
	simrightefferr[i] = sqrt(simrightefferr[i]);
      }
  }




  Double_t ex[numE] = {};
  Double_t Edouble[numE];


  for(int i = 0; i < numE; i++){
    Edouble[i] = (Double_t)Energy[i];
  }


  TGraphErrors *gre1 = new TGraphErrors(numE,Edouble,sum_eff,ex,sum_eff_err);
  TGraphErrors *gre2 = new TGraphErrors(numE,Edouble,sum_eff_crystal,ex,sum_eff_crystal_err);
  TGraphErrors *gre3 = new TGraphErrors(numE,Edouble,simupeff,ex,simupefferr);
  TGraphErrors *gre4 = new TGraphErrors(numE,Edouble,simdowneff,ex,simdownefferr);
  TGraphErrors *gre5 = new TGraphErrors(numE,Edouble,simaboveeff,ex,simaboveefferr);
  TGraphErrors *gre6 = new TGraphErrors(numE,Edouble,simbeloweff,ex,simbelowefferr);
  TGraphErrors *gre7 = new TGraphErrors(numE,Edouble,simlefteff,ex,simleftefferr);
  TGraphErrors *gre8 = new TGraphErrors(numE,Edouble,simrighteff,ex,simrightefferr);

  c1graph1->cd();
  gre1->SetTitle("Clover Array Add-back Efficiency");
  gre1->GetXaxis()->SetTitle("Energy (keV)");
  gre1->GetYaxis()->SetTitle("Efficiency (percent)");
  gre1->SetMarkerColor(4);
  gre1->SetMarkerStyle(21);
  gre1->Draw("AP");

  c1graph2->cd();
  gre2->SetTitle("Total Crystal Efficiency");
  gre2->GetXaxis()->SetTitle("Energy (keV)");
  gre2->GetYaxis()->SetTitle("Efficiency (percent)");
  gre2->SetMarkerColor(4);
  gre2->SetMarkerStyle(21);
  gre2->Draw("AP");

  c1graph3->cd();
  gre3->SetTitle("Upstream Crystal Efficiency");
  gre3->GetXaxis()->SetTitle("Energy (keV)");
  gre3->GetYaxis()->SetTitle("Efficiency (percent)");
  gre3->SetMarkerColor(4);
  gre3->SetMarkerStyle(21);
  gre3->Draw("AP");

  c1graph4->cd();
  gre4->SetTitle("Downstream Crystal Efficiency");
  gre4->GetXaxis()->SetTitle("Energy (keV)");
  gre4->GetYaxis()->SetTitle("Efficiency (percent)");
  gre4->SetMarkerColor(4);
  gre4->SetMarkerStyle(21);
  gre4->Draw("AP");

  c1graph5->cd();
  gre5->SetTitle("Above Crystal Efficiency");
  gre5->GetXaxis()->SetTitle("Energy (keV)");
  gre5->GetYaxis()->SetTitle("Efficiency (percent)");
  gre5->SetMarkerColor(4);
  gre5->SetMarkerStyle(21);
  gre5->Draw("AP");

  c1graph6->cd();
  gre6->SetTitle("Below Crystal Efficiency");
  gre6->GetXaxis()->SetTitle("Energy (keV)");
  gre6->GetYaxis()->SetTitle("Efficiency (percent)");
  gre6->SetMarkerColor(4);
  gre6->SetMarkerStyle(21);
  gre6->Draw("AP");

  c1graph7->cd();
  gre7->SetTitle("Left Crystal Efficiency");
  gre7->GetXaxis()->SetTitle("Energy (keV)");
  gre7->GetYaxis()->SetTitle("Efficiency (percent)");
  gre7->SetMarkerColor(4);
  gre7->SetMarkerStyle(21);
  gre7->Draw("AP");

  c1graph8->cd();
  gre8->SetTitle("Right Crystal Efficiency");
  gre8->GetXaxis()->SetTitle("Energy (keV)");
  gre8->GetYaxis()->SetTitle("Efficiency (percent)");
  gre8->SetMarkerColor(4);
  gre8->SetMarkerStyle(21);
  gre8->Draw("AP");

  gre1->Write();
  gre2->Write();
  gre3->Write();
  gre4->Write();
  gre5->Write();
  gre6->Write();
  gre7->Write();
  gre8->Write();
  
  cout<<"Done!"<<endl;
  
  cout<<"Writing output root file" << endl;
  fout->Write();
  
}

