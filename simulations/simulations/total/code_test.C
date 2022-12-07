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



//main function

void code_test(){

  //expected energies = {42.8,86.5,105.3,123.1,247,591,723,873,996,1004,1274,1596}
  static const int numE = 12;
  static const int numcrystals = 64;
  static const int numdets = 16;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //----------Beginning of experimental analysis-------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  
  Double_t energy, chan, dchan, cts, dcts, eff, deff;
  Double_t dum1, dum2, dum3, dum4;
  ifstream infile;



  //declaring 2D vectors to store efficiency and efficiency error information

  Double_t Efficiency[12][64];
  Double_t Efficiency_Err[12][64];

  //setting initial values to zero since a few crystals are missing
  for(int i = 0; i < 12; i++){
    for(int j = 0; j < 64; j++){
      Efficiency[i][j] = 0;
      Efficiency_Err[i][j] = 0;
    }
  }


  //looping through results files and crystal numbers

  for(int i = 0; i < numcrystals; i++){

    infile.open(Form("Clover_%i_Results.txt",i));

    if(!infile.is_open()) {
      cout << "File " << Form("Clover_%i_Results.txt",i) << " could not be opened, moving on....\n" ;
    } else {
      // Skip first two lines of calibration pars
      int counter = 0;
      infile.ignore(2000,'\n');
      infile.ignore(2000,'\n');
     
	// Reset
	energy = chan = dchan = cts = dcts = eff = deff = 0;

	// Read in data to create master array of efficiencies and efficiency errors

	while(!infile.eof()) {

	  //cout << energy << " " << chan << " "  << dchan << " " << cts << " " << dcts << " " << eff << " " << deff << endl;

	  infile >> energy >> chan >> dchan >> cts >> dcts >> eff >> deff;

	  if(energy == 42.8){
	    Efficiency[0][i] = eff;
	    Efficiency_Err[0][i] = deff;
	  }
	  if(energy == 86.5){
	    Efficiency[1][i] = eff;
	    Efficiency_Err[1][i] = deff;
	  }
	  if(energy == 105.3){
	    Efficiency[2][i] = eff;
	    Efficiency_Err[2][i] = deff;
	  }
	  if(energy == 123.1){
	    Efficiency[3][i] = eff;
	    Efficiency_Err[3][i] = deff;
	  }
	  if(energy == 247.7){
	    Efficiency[4][i] = eff;;
	    Efficiency_Err[4][i] = deff;
	  }
	  if(energy == 591.8){
	    Efficiency[5][i] = eff;
	    Efficiency_Err[5][i] = deff;
	  }
	  if(energy == 723.3){
	    Efficiency[6][i] = eff;
	    Efficiency_Err[6][i] = deff;
	  }
	  if(energy == 873.2){
	    Efficiency[7][i] = eff;
	    Efficiency_Err[7][i] = deff;
	  }
	  if(energy == 996.3){
	    Efficiency[8][i] = eff;
	    Efficiency_Err[8][i] = deff;
	  }
	  if(energy == 1004.7){
	    Efficiency[9][i] = eff;
	    Efficiency_Err[9][i] = deff;
	  }
	  if(energy == 1274.5){
	    Efficiency[10][i] = eff;
	    Efficiency_Err[10][i] = deff;
	  }
	  if(energy == 1596.4){
	    Efficiency[11][i] = eff;
	    Efficiency_Err[11][i] = deff;
	  }
	  if(infile.eof()){
	    infile.close();
	  }
	}
    }
  }  


  //defining crystal regions for experiment (simulation crystal numbers are different)
  vector<int> downcrystals = {48,49,50,51,52,53,54,55,56,57,58,60,61,62,63};
  vector<int> upcrystals = {0,1,2,3,4,6,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,39,40,41,42,43,44,45,46,47};
  vector<int> abovecrystals = {4,6,24,25,26,27,28,29,30,31,32,33,34,35,52,53,54,55};
  vector<int> belowcrystals = {12,13,14,15,16,17,18,19,40,41,42,43,44,45,46,47,60,61,62,63};
  vector<int> rightcrystals = {0,1,2,3,16,17,18,19,20,21,22,23,24,25,26,27,48,49,50,51};
  vector<int> leftcrystals = {8,10,11,28,29,30,31,32,33,34,35,36,37,39,56,57,58};

  //initializing arrays for regions of interest
  double toteff[numE] = {};
  double totefferr[numE] = {};
  double downeff[numE] = {};
  double downefferr[numE] = {};
  double upeff[numE] = {};
  double upefferr[numE] = {};
  double righteff[numE] = {};
  double rightefferr[numE] = {};
  double lefteff[numE] = {};
  double leftefferr[numE] = {};
  double aboveeff[numE] = {};
  double aboveefferr[numE] = {};
  double beloweff[numE] = {};
  double belowefferr[numE] = {};
  double energies[12] = {42.8,86.5,105.3,123.1,247.7,591.8,723.3,873.2,996.3,1004.7,1274.5,1596.4};
  double energieserr[12] = {};
  
  //looping through the crystals
  for(int i = 0; i < numcrystals; i++){
    //looping through energies
    for(int j = 0; j < numE; j++){
      //cout << "Efficiency[" << j << "][" << i << "]: " << Efficiency[j][i] << endl;
      //checking to see if the crystal is downstream
      if(std::find(downcrystals.begin(), downcrystals.end(), i) != downcrystals.end()) {
	//populating downstream efficiency and error vectors
	downeff[j] += Efficiency[j][i];
	downefferr[j] += Efficiency_Err[j][i]*Efficiency_Err[j][i];
      }
      //checking to see if the crystal is upstream
      if(std::find(upcrystals.begin(), upcrystals.end(), i) != upcrystals.end()) {
	//populating upstream efficiency and error vectors
        upeff[j] += Efficiency[j][i];
        upefferr[j] += Efficiency_Err[j][i]*Efficiency_Err[j][i];
      }
      //checking to see if the crystal is on the right
      if(std::find(rightcrystals.begin(), rightcrystals.end(), i) != rightcrystals.end()) {
	//populating right efficiency and error vectors
	righteff[j] += Efficiency[j][i];
	rightefferr[j] += Efficiency_Err[j][i]*Efficiency_Err[j][i];
      }
      //checking to see if the crystal is on the left
      if(std::find(leftcrystals.begin(), leftcrystals.end(), i) != leftcrystals.end()) {
	//populating left efficiency and error vectors
	lefteff[j] += Efficiency[j][i];
	leftefferr[j] += Efficiency_Err[j][i]*Efficiency_Err[j][i];
      }
      //checking to see if the crystal is above
      if(std::find(abovecrystals.begin(), abovecrystals.end(), i) != abovecrystals.end()) {
	//populating above efficiency and error vectors
	aboveeff[j] += Efficiency[j][i];
	aboveefferr[j] += Efficiency_Err[j][i]*Efficiency_Err[j][i];
      }
      //checking to see if the crystal is below
      if(std::find(belowcrystals.begin(), belowcrystals.end(), i) != belowcrystals.end()) {
	//populating below efficiency and error vectors
	beloweff[j] += Efficiency[j][i];
	belowefferr[j] += Efficiency_Err[j][i]*Efficiency_Err[j][i];
      }
      //populating efficiency and error vectors for all the crystals
      toteff[j] += Efficiency[j][i];
      totefferr[j] += Efficiency_Err[j][i]*Efficiency_Err[j][i];
    }

  }



  //calculating errors
  for(int i = 0; i<numE; i++){
    totefferr[i] = sqrt(totefferr[i]);
    belowefferr[i] = sqrt(belowefferr[i]);
    aboveefferr[i] = sqrt(aboveefferr[i]);
    downefferr[i] = sqrt(downefferr[i]);
    upefferr[i] = sqrt(upefferr[i]);
    leftefferr[i] = sqrt(leftefferr[i]);
    rightefferr[i] = sqrt(rightefferr[i]);
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //-------------end of experiment analysis------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------




  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //------------beginning of sim analysis--------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


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
      if(i<=1){
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


 

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //----------end of sim analysis---------------------------------------------
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
 



  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //----------Plotting stuff---------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


  //creating root file to place graphs in
  TFile *f = new TFile("efficiencies.root","recreate");

    //sim plots
    //sim plots

  // TCanvas *c1graph1 = new TCanvas("addback","addback",800,600);
  // TCanvas *c1graph2 = new TCanvas("simTotal","simTotal",800,600);
  // TCanvas *c1graph3 = new TCanvas("simUpstream","simUpstream",800,600);
  // TCanvas *c1graph4 = new TCanvas("simDownstream","simDownstream",800,600);
  // TCanvas *c1graph5 = new TCanvas("simAbove","simAbove",800,600);
  // TCanvas *c1graph6 = new TCanvas("simBelow","simBelow",800,600);
  // TCanvas *c1graph7 = new TCanvas("simLeft","simLeft",800,600);
  // TCanvas *c1graph8 = new TCanvas("simRight","simRight",800,600);

  //canvases for comparing experimental and simulated crystal efficiencies
  TCanvas *Total = new TCanvas("Total", "Total Comparison", 800, 600);
  TCanvas *Downstream = new TCanvas("Downstream", "Downstream Comparison", 800, 600);
  TCanvas *Upstream = new TCanvas("Upstream", "Upstream Comparison", 800, 600);
  TCanvas *Right = new TCanvas("Right", "Right Comparison", 800, 600);
  TCanvas *Left = new TCanvas("Left", "Left Comparison", 800, 600);
  TCanvas *Above = new TCanvas("Above", "Above Comparison", 800, 600);
  TCanvas *Below = new TCanvas("Below", "Below Comparison", 800, 600);

  //canvas for displaying all the graphs at the same time
  TCanvas *All = new TCanvas("All", "All", 1000, 1200);
  All->DivideSquare(7,0.005,0.005);

  //multigraphs for showing both TGraphs at once for each region
  TMultiGraph *totmg = new TMultiGraph();
  TMultiGraph *dowmg = new TMultiGraph();
  TMultiGraph *upmg = new TMultiGraph();
  TMultiGraph *rigmg = new TMultiGraph();
  TMultiGraph *lefmg = new TMultiGraph();
  TMultiGraph *abomg = new TMultiGraph();
  TMultiGraph *belmg = new TMultiGraph();

    
  Double_t ex[numE] = {};
  Double_t Edouble[numE];


  for(int i = 0; i < numE; i++){
    Edouble[i] = (Double_t)Energy[i];
  }

  //simulation TGraphs
  auto gre1 = new TGraphErrors(numE,Edouble,sum_eff,ex,sum_eff_err);
  auto gre2 = new TGraphErrors(numE,Edouble,sum_eff_crystal,ex,sum_eff_crystal_err);
  auto gre3 = new TGraphErrors(numE,Edouble,simupeff,ex,simupefferr);
  auto gre4 = new TGraphErrors(numE,Edouble,simdowneff,ex,simdownefferr);
  auto gre5 = new TGraphErrors(numE,Edouble,simaboveeff,ex,simaboveefferr);
  auto gre6 = new TGraphErrors(numE,Edouble,simbeloweff,ex,simbelowefferr);
  auto gre7 = new TGraphErrors(numE,Edouble,simlefteff,ex,simleftefferr);
  auto gre8 = new TGraphErrors(numE,Edouble,simrighteff,ex,simrightefferr);
  //experiment TGraphs
  auto grtot = new TGraphErrors(12,energies,toteff,energieserr,totefferr);
  auto grdow = new TGraphErrors(12,energies,downeff,energieserr,downefferr);
  auto grup = new TGraphErrors(12,energies,upeff,energieserr,upefferr);
  auto grrig = new TGraphErrors(12,energies,righteff,energieserr,rightefferr);
  auto grlef = new TGraphErrors(12,energies,lefteff,energieserr,leftefferr);
  auto grabo = new TGraphErrors(12,energies,aboveeff,energieserr,aboveefferr);
  auto grbel = new TGraphErrors(12,energies,beloweff,energieserr,belowefferr);

  //setting simulation points in blue
  gre2->SetMarkerColor(kBlue);
  gre3->SetMarkerColor(kBlue);
  gre4->SetMarkerColor(kBlue);
  gre5->SetMarkerColor(kBlue);
  gre6->SetMarkerColor(kBlue);
  gre7->SetMarkerColor(kBlue);
  gre8->SetMarkerColor(kBlue);

  //setting experiment points in red
  grtot->SetMarkerColor(kRed);
  grdow->SetMarkerColor(kRed);
  grup->SetMarkerColor(kRed);
  grrig->SetMarkerColor(kRed);
  grlef->SetMarkerColor(kRed);
  grabo->SetMarkerColor(kRed);
  grbel->SetMarkerColor(kRed);

  //setting simulation marker styles
  gre2->SetMarkerStyle(21);
  gre3->SetMarkerStyle(21);
  gre4->SetMarkerStyle(21);
  gre5->SetMarkerStyle(21);
  gre6->SetMarkerStyle(21);
  gre7->SetMarkerStyle(21);
  gre8->SetMarkerStyle(21);

  //setting experiment marker styles
  grtot->SetMarkerStyle(20);
  grdow->SetMarkerStyle(20);
  grup->SetMarkerStyle(20);
  grrig->SetMarkerStyle(20);
  grlef->SetMarkerStyle(20);
  grabo->SetMarkerStyle(20);
  grbel->SetMarkerStyle(20);


  //adding individual graphs to their multigraph and drawing those on their canvas

  Total->cd();
  totmg->Add(gre2);
  totmg->Add(grtot);
  totmg->GetHistogram()->SetTitle("Total");
  totmg->GetXaxis()->SetTitle("Energy (keV)");
  totmg->GetYaxis()->SetTitle("Efficiency (percent)");
  totmg->Draw("ap");
  totmg->Write();

  Upstream->cd();
  upmg->Add(gre3);
  upmg->Add(grup);
  upmg->GetHistogram()->SetTitle("Upstream");
  upmg->GetXaxis()->SetTitle("Energy (keV)");
  upmg->GetYaxis()->SetTitle("Efficiency (percent)");
  upmg->Draw("ap");
  upmg->Write();

  Downstream->cd();
  dowmg->Add(gre4);
  dowmg->Add(grdow);
  dowmg->GetHistogram()->SetTitle("Downstream");
  dowmg->GetXaxis()->SetTitle("Energy (keV)");
  dowmg->GetYaxis()->SetTitle("Efficiency (percent)");
  dowmg->Draw("ap");
  dowmg->Write();

  Above->cd();
  abomg->Add(gre5);
  abomg->Add(grabo);
  abomg->GetHistogram()->SetTitle("Above");
  abomg->GetXaxis()->SetTitle("Energy (keV)");
  abomg->GetYaxis()->SetTitle("Efficiency (percent)");
  abomg->Draw("ap");
  abomg->Write();

  Below->cd();
  belmg->Add(gre6);
  belmg->Add(grbel);
  belmg->GetHistogram()->SetTitle("Below");
  belmg->GetXaxis()->SetTitle("Energy (keV)");
  belmg->GetYaxis()->SetTitle("Efficiency (percent)");
  belmg->Draw("ap");
  belmg->Write();

  Left->cd();
  lefmg->Add(gre7);
  lefmg->Add(grlef);
  lefmg->GetHistogram()->SetTitle("Left");
  lefmg->GetXaxis()->SetTitle("Energy (keV)");
  lefmg->GetYaxis()->SetTitle("Efficiency (percent)");
  lefmg->Draw("ap");
  lefmg->Write();

  Right->cd();
  rigmg->Add(gre8);
  rigmg->Add(grrig);
  rigmg->GetHistogram()->SetTitle("Right");
  rigmg->GetXaxis()->SetTitle("Energy (keV)");
  rigmg->GetYaxis()->SetTitle("Efficiency (percent)");
  rigmg->Draw("ap");
  rigmg->Write();

  All->cd(1);
  totmg->Draw("ap");
  All->cd(2);
  upmg->Draw("ap");
  All->cd(3);
  dowmg->Draw("ap");
  All->cd(4);
  rigmg->Draw("ap");
  All->cd(5);
  lefmg->Draw("ap");
  All->cd(6);
  abomg->Draw("ap");
  All->cd(7);
  belmg->Draw("ap");
  All->Write();


  // c1graph1->cd();
  // gre1->SetTitle("Sim Clover Array Add-back Efficiency");
  // gre1->GetXaxis()->SetTitle("Energy (keV)");
  // gre1->GetYaxis()->SetTitle("Efficiency (percent)");
  // gre1->SetMarkerColor(4);
  // gre1->SetMarkerStyle(21);
  // gre1->Draw("AP");

  // c1graph2->cd();
  // gre2->SetTitle("Sim Total Crystal Efficiency");
  // gre2->GetXaxis()->SetTitle("Energy (keV)");
  // gre2->GetYaxis()->SetTitle("Efficiency (percent)");
  // gre2->SetMarkerColor(4);
  // gre2->SetMarkerStyle(21);
  // gre2->Draw("AP");

  // c1graph3->cd();
  // gre3->SetTitle("Sim Upstream Crystal Efficiency");
  // gre3->GetXaxis()->SetTitle("Energy (keV)");
  // gre3->GetYaxis()->SetTitle("Efficiency (percent)");
  // gre3->SetMarkerColor(4);
  // gre3->SetMarkerStyle(21);
  // gre3->Draw("AP");

  // c1graph4->cd();
  // gre4->SetTitle("Sim Downstream Crystal Efficiency");
  // gre4->GetXaxis()->SetTitle("Energy (keV)");
  // gre4->GetYaxis()->SetTitle("Efficiency (percent)");
  // gre4->SetMarkerColor(4);
  // gre4->SetMarkerStyle(21);
  // gre4->Draw("AP");

  // c1graph5->cd();
  // gre5->SetTitle("Above Crystal Efficiency");
  // gre5->GetXaxis()->SetTitle("Energy (keV)");
  // gre5->GetYaxis()->SetTitle("Efficiency (percent)");
  // gre5->SetMarkerColor(4);
  // gre5->SetMarkerStyle(21);
  // gre5->Draw("AP");

  // c1graph6->cd();
  // gre6->SetTitle("Sim Below Crystal Efficiency");
  // gre6->GetXaxis()->SetTitle("Energy (keV)");
  // gre6->GetYaxis()->SetTitle("Efficiency (percent)");
  // gre6->SetMarkerColor(4);
  // gre6->SetMarkerStyle(21);
  // gre6->Draw("AP");

  // c1graph7->cd();
  // gre7->SetTitle("Sim Left Crystal Efficiency");
  // gre7->GetXaxis()->SetTitle("Energy (keV)");
  // gre7->GetYaxis()->SetTitle("Efficiency (percent)");
  // gre7->SetMarkerColor(4);
  // gre7->SetMarkerStyle(21);
  // gre7->Draw("AP");

  // c1graph8->cd();
  // gre8->SetTitle("Sim Right Crystal Efficiency");
  // gre8->GetXaxis()->SetTitle("Energy (keV)");
  // gre8->GetYaxis()->SetTitle("Efficiency (percent)");
  // gre8->SetMarkerColor(4);
  // gre8->SetMarkerStyle(21);
  // gre8->Draw("AP");

  // gre1->Write();
  // gre2->Write();
  // gre3->Write();
  // gre4->Write();
  // gre5->Write();
  // gre6->Write();
  // gre7->Write();
  // gre8->Write();
  


    //experiment plots
    //experiment plots

// // Total Efficiency vs energy plot
//   auto grtot = new TGraphErrors(12,energies,toteff,energieserr,totefferr);
//   grtot->SetMarkerStyle(21);
//   grtot->SetTitle("Experiment Total Crystal Efficiency");
//   grtot->GetXaxis()->SetTitle("Energy [keV]");
//   grtot->GetYaxis()->SetTitle("Efficiency [%]");
//   auto Total = new TCanvas("Total","Total",800,600);
//   grtot->Draw("ap");
//   //ci->Draw("e3 same");
//   Total->Write();


// // Downstream Efficiency vs energy plot
//   auto grdow = new TGraphErrors(12,energies,downeff,energieserr,downefferr);
//   grdow->SetMarkerStyle(21);
//   grdow->SetTitle("Experiment Downstream Crystal Efficiency");
//   grdow->GetXaxis()->SetTitle("Energy [keV]");
//   grdow->GetYaxis()->SetTitle("Efficiency [%]");
//   auto downstream = new TCanvas("downstream","downstream",800,600);
//   grdow->Draw("ap");
//   //ci->Draw("e3 same");
//   downstream->Write();


// // Upstream Efficiency vs energy plot
//   auto grup = new TGraphErrors(12,energies,upeff,energieserr,upefferr);
//   grup->SetMarkerStyle(21);
//   grup->SetTitle("Experiment Upstream Crystal Efficiency");
//   grup->GetXaxis()->SetTitle("Energy [keV]");
//   grup->GetYaxis()->SetTitle("Efficiency [%]");
//   auto upstream = new TCanvas("upstream","upstream",800,600);
//   grup->Draw("ap");
//   //ci->Draw("e3 same");
//   upstream->Write();


// // Right Efficiency vs energy plot
//   auto grrig = new TGraphErrors(12,energies,righteff,energieserr,rightefferr);
//   grrig->SetMarkerStyle(21);
//   grrig->SetTitle("Experiment Right Crystal Efficiency");
//   grrig->GetXaxis()->SetTitle("Energy [keV]");
//   grrig->GetYaxis()->SetTitle("Efficiency [%]");
//   auto right = new TCanvas("right","right",800,600);
//   grrig->Draw("ap");
//   //ci->Draw("e3 same");
//   right->Write();


// // Left Efficiency vs energy plot
//   auto grlef = new TGraphErrors(12,energies,lefteff,energieserr,leftefferr);
//   grlef->SetMarkerStyle(21);
//   grlef->SetTitle("Experiment Left Crystal Efficiency");
//   grlef->GetXaxis()->SetTitle("Energy [keV]");
//   grlef->GetYaxis()->SetTitle("Efficiency [%]");
//   auto left = new TCanvas("Left","Left",800,600);
//   grlef->Draw("ap");
//   //ci->Draw("e3 same");
//   left->Write();


// // Above Efficiency vs energy plot
//   auto grabo = new TGraphErrors(12,energies,aboveeff,energieserr,aboveefferr);
//   grabo->SetMarkerStyle(21);
//   grabo->SetTitle("Experiment Above Crystal Efficiency");
//   grabo->GetXaxis()->SetTitle("Energy [keV]");
//   grabo->GetYaxis()->SetTitle("Efficiency [%]");
//   auto above = new TCanvas("Above","Above",800,600);
//   grabo->Draw("ap");
//   //ci->Draw("e3 same");
//   above->Write();


// // Below Efficiency vs energy plot
//   auto grbel = new TGraphErrors(12,energies,beloweff,energieserr,belowefferr);
//   grbel->SetMarkerStyle(21);
//   grbel->SetTitle("Experiment Below Crystal Efficiency");
//   grbel->GetXaxis()->SetTitle("Energy [keV]");
//   grbel->GetYaxis()->SetTitle("Efficiency [%]");
//   auto below = new TCanvas("Below","Below",800,600);
//   grbel->Draw("ap");
//   //ci->Draw("e3 same");
//   below->Write();











    
}
