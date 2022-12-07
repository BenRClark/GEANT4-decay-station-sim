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
#include <math.h>


//function for discriminating between experimental and simulated energy lists. Simulation needs more energies analyzed for use in summing corrections, but only want to compare with energies that are analyzed in experimental results. To change which ones are compared, add of remove indices as necessary. In this function, the "i" index refers to the simulated energies array Energy[i], so we are ignoring Energy[4], Energy[6], and so on.

//Credit to Lee Ward for this function idea 05/13/2021-DCS
bool is_exp_energy(int i){
  if(i == 4 || i == 6 || i == 7 || i == 9 || i == 10 || i == 12 || i == 13 || i == 14 || i == 15 || i == 16 || i == 17 || i == 19 || i == 21 || i == 22 || i == 25 || i == 26){
    return false;
  }
  else{
    return true;
  }
}




//main function

void analysis_comparison(){

  //expected energies = {42.8,86.5,105.3,123.1,247,591,723,873,996,1004,1274,1596}
  static const int numE = 12;
  static const int simnumE = 29;
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
  //array of energies analyzed in the experimental results
  double energies[12] = {42.8,86.5,105.3,123.1,247.7,591.8,723.3,873.2,996.3,1004.7,1274.5,1596.4};
  double energieserr[12] = {};
  
  //looping through the crystals
  //this is probably unecessary since the summing corrections were added, but it might be convenient to have on hand if you want to see what the summing corrections did -DCS
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

  //declaring array of simulated energies
  int Energy[simnumE];  

  Energy[0] = 43;  //in keV
  Energy[1] = 87; //in keV
  Energy[2] = 105; //in keV
  Energy[3] = 123; //in keV
  Energy[4] = 176; //in keV. not in e17011 experimental results 
  Energy[5] = 247; //in keV
  Energy[6] = 428; //in keV. not in e17011 experimental results 
  Energy[7] = 444; //in keV. not in e17011 experimental results
  Energy[8] = 463; //in keV. not in e17011 experimental results 
  Energy[9] = 478; //in keV. not in e17011 experimental results
  Energy[10] = 582; //in keV. not in e17011 experimental results
  Energy[11] = 591; //in keV
  Energy[12] = 601; //in keV. not in e17011 experimental results
  Energy[13] = 612; //in keV. not in e17011 experimental results
  Energy[14] = 625; //in keV. not in e17011 experimental results
  Energy[15] = 636; //in keV. not in e17011 experimental results 
  Energy[16] = 677; //in keV. not in e17011 experimental results
  Energy[17] = 692; //in keV. not in e17011 experimental results
  Energy[18] = 723;  //in keV
  Energy[19] = 757; //in keV. not in e17011 experimental results
  Energy[20] = 873; //in keV
  Energy[21] = 893; //in keV. not in e17011 experimental results
  Energy[22] = 904; //in keV. not in e17011 experimental results
  Energy[23] = 996; //in keV
  Energy[24] = 1004; //in keV
  Energy[25] = 1119; //in keV. not in e17011 experimental results
  Energy[26] = 1246; //in keV. not in e17011 experimental results
  Energy[27] = 1274; //in keV
  Energy[28] = 1596; //in keV



  //Number of Gammas in Simulation

  //was originally 10000000
  double nGamma = 10000000;

  // Efficiencies [Energy][Det Number]
  double Peak_Efficiency[simnumE][numdets] = {};
  double Peak_Efficiency_Err[simnumE][numdets] = {};

  double Peak_Efficiency_crystal[simnumE][numcrystals] = {};
  double Peak_Efficiency_Err_crystal[simnumE][numcrystals] = {};
  
  double Total_Efficiency[simnumE][numdets] = {};
  double Total_Efficiency_Err[simnumE][numdets] = {};

  double Total_Efficiency_crystal[simnumE][numcrystals] = {};
  double Total_Efficiency_Err_crystal[simnumE][numcrystals] = {};

  //FIles
  TFile *f1t[simnumE];
  
  for(int i=0; i<simnumE; i++) {
    //change this file path to data drive once simulated root files are moved over there. -DCS 05/20/2021
    f1t[i] = new TFile(Form("/home/dcs411/e17011/Efficiency/sim_analysis/simrootfiles/%dkeVtest_0.root",Energy[i]));
  }

  //Fits & Histos
  TF1 *fClover[simnumE][numdets] = {};
  TH1D *hClover[simnumE][numdets];
  TH1D *hClover_crystal[simnumE][numcrystals];
  TF1 *fClover_crystal[simnumE][numcrystals] = {};

  
  for(int i=0; i<simnumE; i++) {
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

  for(int i=0; i<simnumE; i++) {
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
  


  //output root file
  //move this to data drive whenever given permission DCS 05/20/2021
  TFile *fout = new TFile("/home/dcs411/e17011/Efficiency/sim_analysis/simrootfiles/test_0.root","RECREATE");

  for(int i=0; i<simnumE; i++) {
    for(int j=0; j<numdets; j++) {

   
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
      cout << "Special Condition" << endl;
      // if(j == 0 || j == 4){
      // 	cout<<"Clover Missing"<<endl;
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

  for(int i=0; i<simnumE; i++) {
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
      // double area;
      // double area_err;
      // double eff;
      // double eff_err;
      // double total_eff;
      // double total_eff_err;
      // if(j == 0 || j == 1 || j == 2 || j == 3 || j == 16 || j == 17 || j == 18 || j == 19){
      // 	area = 0;
      // 	area_err = 0;
      // 	eff = 0;
      // 	eff_err = 0;
      // 	total_eff = 0;
      // 	total_eff_err = 0;
      // }
      // else{
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
      //}
      

      //cout<<i<<"  "<<j<<"  "<<eff<<" +/- "<<eff_err<<"  "<<100*total_eff<<" +/- "<<100*total_eff_err<<endl;
      cout << i << "  " << j << "  " << area << "  " << total_eff*1.0*nGamma << endl;

      // Efficiencies [Energy][Det Number]

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


  //writing simulated data to output text file
  //need to path this to data drive later
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

  //master efficiency arrays are populated. now we can manipulate them

  Double_t sum_eff_Full[simnumE] = {};
  Double_t sum_eff_err_Full[simnumE]= {};
  Double_t sum_eff_crystal_Full[simnumE] = {};
  Double_t sum_eff_crystal_err_Full[simnumE] = {};
  Double_t err_temp1[simnumE] = {};
  Double_t err_temp2[simnumE] = {};

  //summing sim efficiencies for each energy and storing in array so we can plot later

  for(int j = 0; j < numdets; j++){
    for(int i = 0; i < simnumE; i++){
	sum_eff_Full[i] += Peak_Efficiency[i][j];
	err_temp1[i] += pow(Peak_Efficiency_Err[i][j],2);
    }
  }


  for(int j = 0; j < numcrystals; j++){
    for(int i = 0; i < simnumE; i++){
	sum_eff_crystal_Full[i] += Peak_Efficiency_crystal[i][j];
	err_temp2[i] += pow(Peak_Efficiency_Err_crystal[i][j],2);
    }
  }


  for(int i = 0; i < numE; i++){
    //seting 43keV error to 0 for now
    if(i<=1){
      sum_eff_err_Full[i] = 0;
      sum_eff_crystal_err_Full[i] = 0;
    }
    else{
      sum_eff_err_Full[i] = err_temp1[i]*sqrt(err_temp1[i]);
      sum_eff_crystal_err_Full[i] = err_temp2[i]*sqrt(err_temp2[i]);
    }
  }
  //cout << sum_eff_crystal_err[i] << endl;

  //declaring crystal number vectors for regions of interest in the simulation
  vector<int> simupcrystals = {8,9,10,11,16,17,18,19,32,33,34,35,36,37,38,39,4,5,6,7,48,49,50,51,28,29,30,31,52,53,54,55,16,17,18,19,56,57,58,59,40,41,42,43,60,61,62,63};
  vector<int> simdowncrystals = {0,1,2,3,24,25,26,27,28,29,30,31,48,49,50,51};
  vector<int> simleftcrystals = {0,1,2,3,4,5,6,7,8,9,10,11,48,49,50,51,60,61,62,63};
  vector<int> simrightcrystals = {12,13,14,15,16,17,18,19,20,21,22,23,52,53,54,55,56,57,58,59};
  vector<int> simabovecrystals = {24,25,26,27,28,29,30,31,32,33,34,35,48,49,50,51,52,53,54,55};
  vector<int> simbelowcrystals = {36,37,38,39,40,41,42,43,44,45,46,47,56,57,58,59,60,61,62,63};

  double simupeff_Full[simnumE] = {};
  double simupefferr_Full[simnumE] = {};
  double simdowneff_Full[simnumE] = {};
  double simdownefferr_Full[simnumE] = {};
  double simaboveeff_Full[simnumE] = {};
  double simaboveefferr_Full[simnumE] = {};
  double simbeloweff_Full[simnumE] = {};
  double simbelowefferr_Full[simnumE] = {};
  double simlefteff_Full[simnumE] = {};
  double simleftefferr_Full[simnumE] = {};
  double simrighteff_Full[simnumE] = {};
  double simrightefferr_Full[simnumE] = {};

  //populating region vectors
  //looping through the crystals
  for(int i = 0; i < numcrystals; i++){
    //looping through energies
    for(int j = 0; j < simnumE; j++){
      //checking to see if the crystal is downstream
      if(std::find(simdowncrystals.begin(), simdowncrystals.end(), i) != simdowncrystals.end()) {
	//populating downstream efficiency and error vectors
	simdowneff_Full[j] += Peak_Efficiency_crystal[j][i];
	simdownefferr_Full[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is upstream
      if(std::find(simupcrystals.begin(), simupcrystals.end(), i) != simupcrystals.end()) {
	//populating upstream efficiency and error vectors
	simupeff_Full[j] += Peak_Efficiency_crystal[j][i];
	simupefferr_Full[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is on the right
      if(std::find(simrightcrystals.begin(), simrightcrystals.end(), i) != simrightcrystals.end()) {
	//populating right efficiency and error vectors
	simrighteff_Full[j] += Peak_Efficiency_crystal[j][i];
	simrightefferr_Full[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is on the left
      if(std::find(simleftcrystals.begin(), simleftcrystals.end(), i) != simleftcrystals.end()) {
	//populating left efficiency and error vectors
	simlefteff_Full[j] += Peak_Efficiency_crystal[j][i];
	simleftefferr_Full[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is above
      if(std::find(simabovecrystals.begin(), simabovecrystals.end(), i) != simabovecrystals.end()) {
	//populating above efficiency and error vectors
	simaboveeff_Full[j] += Peak_Efficiency_crystal[j][i];
	simaboveefferr_Full[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
      //checking to see if the crystal is below
      if(std::find(simbelowcrystals.begin(), simbelowcrystals.end(), i) != simbelowcrystals.end()) {
	//populating below efficiency and error vectors
	simbeloweff_Full[j] += Peak_Efficiency_crystal[j][i];
	simbelowefferr_Full[j] += pow(Peak_Efficiency_Err_crystal[j][i],2);
      }
    }
   

  }

  //calculating errors
  //setting 43keV energy err to 0 right now because not enough statistics and throws off the graph
  for(int i = 0; i<simnumE; i++){
    if(i<=1){
      simbelowefferr_Full[i] = 0;
      simaboveefferr_Full[i] = 0;
      simdownefferr_Full[i] = 0;
      simupefferr_Full[i] = 0;
      simleftefferr_Full[i] = 0;
      simrightefferr_Full[i] = 0;
    }
    else{
      simbelowefferr_Full[i] = sqrt(simbelowefferr_Full[i]);
      simaboveefferr_Full[i] = sqrt(simaboveefferr_Full[i]);
      simdownefferr_Full[i] = sqrt(simdownefferr_Full[i]);
      simupefferr_Full[i] = sqrt(simupefferr_Full[i]);
      simleftefferr_Full[i] = sqrt(simleftefferr_Full[i]);
      simrightefferr_Full[i] = sqrt(simrightefferr_Full[i]);
    }
  }


  //need to remove energies that aren't included in experiment for comparison
  //can do that by removing the elements of the Full arrays associated with those energies we aren't interested in. 

  Double_t sum_eff_trunc[numE] = {};
  Double_t sum_eff_err_trunc[numE]= {};
  Double_t sum_eff_crystal_trunc[numE] = {};
  Double_t sum_eff_crystal_err_trunc[numE] = {};
  double simupeff_trunc[simnumE] = {};
  double simupefferr_trunc[simnumE] = {};
  double simdowneff_trunc[simnumE] = {};
  double simdownefferr_trunc[simnumE] = {};
  double simaboveeff_trunc[simnumE] = {};
  double simaboveefferr_trunc[simnumE] = {};
  double simbeloweff_trunc[simnumE] = {};
  double simbelowefferr_trunc[simnumE] = {};
  double simlefteff_trunc[simnumE] = {};
  double simleftefferr_trunc[simnumE] = {};
  double simrighteff_trunc[simnumE] = {};
  double simrightefferr_trunc[simnumE] = {};

  int counter = 0;
  

  //truncating efficiency arrays to include only energies we care about for comparison purposes. In sim there are 29 energies. We only want to compare with the 12 we have for experiment. The is_exp_energy function defined at the beginning of this code takes care of this
  for(int j = 0; j < simnumE; j++){
    if(is_exp_energy(j)){
      sum_eff_trunc[counter] = sum_eff_Full[j];
      sum_eff_err_trunc[counter] = sum_eff_err_Full[j];
      sum_eff_crystal_trunc[counter] = sum_eff_crystal_Full[j];
      sum_eff_crystal_err_trunc[counter] = sum_eff_crystal_err_Full[j];
      simupeff_trunc[counter] = simupeff_Full[j];
      simupefferr_trunc[counter] = simupefferr_Full[j];
      simdowneff_trunc[counter] = simdowneff_Full[j];
      simdownefferr_trunc[counter] = simdownefferr_Full[j];
      simaboveeff_trunc[counter] = simaboveeff_Full[j];
      simaboveefferr_trunc[counter] = simaboveefferr_Full[j];
      simbeloweff_trunc[counter] = simbeloweff_Full[j];
      simbelowefferr_trunc[counter] = simbelowefferr_Full[j];
      simlefteff_trunc[counter] = simlefteff_Full[j];
      simleftefferr_trunc[counter] = simleftefferr_Full[j];
      simrighteff_trunc[counter] = simrighteff_Full[j];
      simrightefferr_trunc[counter] = simrightefferr_Full[j];
      counter++;
      }
  }

 

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //----------end of sim analysis---------------------------------------------
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------



  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //----------summing corrections---------------------------------------------
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------

  //declare new summing corrected arrays. Names are cumbersome, but alleviate any anbiguity
  
  Double_t Summing_Corrections[numE][numcrystals] = {};
  Double_t Efficiency_corrected[numE][numcrystals] = {};
  Double_t Efficiency_Err_corrected[numE][numcrystals] = {};
  Double_t Sum_Efficiency_corrected[numE] = {};
  Double_t Sum_Efficiency_Err_corrected[numE] = {};
  Double_t Sum_Eff_Up_corr[numE] = {};
  Double_t Sum_Eff_Down_corr[numE] = {};
  Double_t Sum_Eff_Left_corr[numE] = {};
  Double_t Sum_Eff_Right_corr[numE] = {};
  Double_t Sum_Eff_Above_corr[numE] = {};
  Double_t Sum_Eff_Below_corr[numE] = {};
  Double_t Sum_Err_temp1[numE] = {};
  Double_t Sum_Err_temp2[numE] = {};
  Double_t Sum_Err_temp3[numE] = {};
  Double_t Sum_Err_temp4[numE] = {};
  Double_t Sum_Err_temp5[numE] = {};
  Double_t Sum_Err_temp6[numE] = {};
  Double_t Sum_Err_temp7[numE] = {};
  Double_t Sum_Eff_Err_Tot_corr[numE] = {};
  Double_t Sum_Eff_Err_Down_corr[numE] = {};
  Double_t Sum_Eff_Err_Up_corr[numE] = {};
  Double_t Sum_Eff_Err_Right_corr[numE] = {};
  Double_t Sum_Eff_Err_Left_corr[numE] = {};
  Double_t Sum_Eff_Err_Above_corr[numE] = {};
  Double_t Sum_Eff_Err_Below_corr[numE] = {};

  //checking to see how good the peak fits are on the 1596keV
  for(int i = 0; i < numcrystals; i++){
    cout << "Peak_Efficiency_crystal[11][" << i << "]: " << Peak_Efficiency_crystal[11][i] << endl;
  }


  
  //computationally cumbersome. can't think of a better way at the moment
  //looping through crystals and applying summing corrections
  for(int i = 0; i < numcrystals; i++){

    //43keV no correction needed
    Summing_Corrections[0][i] = 1.0;
    Efficiency_corrected[0][i] = Efficiency[0][i]/Summing_Corrections[0][i];
    Efficiency_Err_corrected[0][i] = Efficiency_Err[0][i];

    //87keV no correction needed
    Summing_Corrections[1][i] = 1.0;
    Efficiency_corrected[1][i] = Efficiency[1][i]/Summing_Corrections[1][i];
    Efficiency_Err_corrected[1][i] = Efficiency_Err[1][i];

    //105keV no correction needed
    Summing_Corrections[2][i] = 1.0;
    Efficiency_corrected[2][i] = Efficiency[2][i]/Summing_Corrections[2][i];
    Efficiency_Err_corrected[2][i] = Efficiency_Err[2][i];

    //i have separated terms inside the "sqrt()" with a "   +   " so they're easily identifiable

    //123keV------------------------------------------------------------------
    Summing_Corrections[3][i] = (1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]);

    Efficiency_corrected[3][i] = Efficiency[3][i]/Summing_Corrections[3][i];  

  //yikes
    Efficiency_Err_corrected[3][i] = sqrt( (pow((1.0/(1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]))*Efficiency_Err[3][i],2))   +   (pow(((0.055*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[5][i],2))   +   (pow(((0.072*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[11][i],2))   +   (pow(((0.019*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[17][i],2))   +   (pow(((0.120*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[18][i],2))   +   (pow(((0.049*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[19][i],2))   +   (pow(((0.130*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[20][i],2))   +   (pow(((0.201*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[24][i],2))   +   (pow(((0.010*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[26][i],2))    +   (pow(((0.401*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[27][i],2))   +   (pow(((0.021*Efficiency[3][i])/pow((1.0 - 0.072*Total_Efficiency_crystal[5][i] - 0.055*Total_Efficiency_crystal[8][i] - 0.019*Total_Efficiency_crystal[17][i] - 0.120*Total_Efficiency_crystal[18][i] - 0.049*Total_Efficiency_crystal[19][i] - 0.130*Total_Efficiency_crystal[20][i] - 0.201*Total_Efficiency_crystal[24][i] - 0.010*Total_Efficiency_crystal[26][i] - 0.401*Total_Efficiency_crystal[27][i] - 0.021*Total_Efficiency_crystal[28][i]),2))*Total_Efficiency_Err_crystal[28][i],2)) );


    //247keV----------------------------------------------------------------------------
    Summing_Corrections[4][i] = (1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]);

    Efficiency_corrected[4][i] = Efficiency[4][i]/Summing_Corrections[4][i];

    //bigger yikes
    Efficiency_Err_corrected[4][i] = sqrt( pow((1.0/(1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]))*Efficiency_Err[4][i],2)   +   pow(((0.287*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[0][i],2)   +   pow(((0.455*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[3][i],2)   +   pow(((0.072*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[7][i],2)   +   pow(((0.022*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[10][i],2)   +   pow(((0.134*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[11][i],2)   +   pow(((0.015*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[13][i],2)   +   pow(((0.043*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[14][i],2)   +   pow(((0.022*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[16][i],2)   +   pow(((0.039*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[18][i],2)   +   pow(((0.613*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[19][i],2)   +   pow(((0.059*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[21][i],2)   +   pow(((0.022*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[22][i],2)   +   pow(((0.130*Efficiency[4][i])/pow((1.0 - 0.287*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.072*Total_Efficiency_crystal[7][i] - 0.022*Total_Efficiency_crystal[10][i] - 0.134*Total_Efficiency_crystal[11][i] - 0.015*Total_Efficiency_crystal[13][i] - 0.043*Total_Efficiency_crystal[14][i] - 0.022*Total_Efficiency_crystal[16][i] - 0.039*Total_Efficiency_crystal[18][i] - 0.613*Total_Efficiency_crystal[19][i] - 0.059*Total_Efficiency_crystal[21][i] - 0.022*Total_Efficiency_crystal[22][i] - 0.130*Total_Efficiency_crystal[26][i]),2))*Total_Efficiency_Err_crystal[26][i],2) );

    //591keV------------------------------------------------------------------------------------------------------
    Summing_Corrections[5][i] = (1.0 - 0.297*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.178*Total_Efficiency_crystal[5][i] - 0.196*Total_Efficiency_crystal[19][i] - 0.800*Total_Efficiency_crystal[24][i]);

    Efficiency_corrected[5][i] = Efficiency[5][i]/Summing_Corrections[5][i];

    Efficiency_Err_corrected[5][i] = sqrt( pow((1.0/(1.0 - 0.297*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.178*Total_Efficiency_crystal[5][i] - 0.196*Total_Efficiency_crystal[19][i] - 0.800*Total_Efficiency_crystal[24][i])*Efficiency_Err[5][i]),2)   +   pow(((0.297*Efficiency[5][i]/pow((1.0 - 0.297*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.178*Total_Efficiency_crystal[5][i] - 0.196*Total_Efficiency_crystal[19][i] - 0.800*Total_Efficiency_crystal[24][i]),2)))*Total_Efficiency_Err_crystal[0][i],2)   +   pow(((0.455*Efficiency[5][i]/pow((1.0 - 0.297*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.178*Total_Efficiency_crystal[5][i] - 0.196*Total_Efficiency_crystal[19][i] - 0.800*Total_Efficiency_crystal[24][i]),2)))*Total_Efficiency_Err_crystal[3][i],2)   +   pow(((0.178*Efficiency[5][i]/pow((1.0 - 0.297*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.178*Total_Efficiency_crystal[5][i] - 0.196*Total_Efficiency_crystal[19][i] - 0.800*Total_Efficiency_crystal[24][i]),2)))*Total_Efficiency_Err_crystal[5][i],2)   +   pow(((0.196*Efficiency[5][i]/pow((1.0 - 0.297*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.178*Total_Efficiency_crystal[5][i] - 0.196*Total_Efficiency_crystal[19][i] - 0.800*Total_Efficiency_crystal[24][i]),2)))*Total_Efficiency_Err_crystal[19][i],2)   +   pow(((0.800*Efficiency[5][i]/pow((1.0 - 0.297*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.178*Total_Efficiency_crystal[5][i] - 0.196*Total_Efficiency_crystal[19][i] - 0.800*Total_Efficiency_crystal[24][i]),2)))*Total_Efficiency_Err_crystal[24][i],2) );

    //723keV-----------------------------------------------------------------------------------------------------
    Summing_Corrections[6][i] = (1.0 - 0.154*Total_Efficiency_crystal[0][i] - 0.243*Total_Efficiency_crystal[3][i] - 0.013*Total_Efficiency_crystal[5][i] - 0.014*Total_Efficiency_crystal[14][i] - 0.518*Total_Efficiency_crystal[20][i] - 0.465*Total_Efficiency_crystal[23][i]);

    Efficiency_corrected[6][i] = Efficiency[6][i]/Summing_Corrections[6][i];

    Efficiency_Err_corrected[6][i] = sqrt( pow(((1.0/(1.0 - 0.154*Total_Efficiency_crystal[0][i] - 0.243*Total_Efficiency_crystal[3][i] - 0.013*Total_Efficiency_crystal[5][i] - 0.014*Total_Efficiency_crystal[14][i] - 0.518*Total_Efficiency_crystal[20][i] - 0.465*Total_Efficiency_crystal[23][i]))*Efficiency_Err[6][i]),2)   +   pow((((0.154*Efficiency[6][i])/pow((1.0 - 0.154*Total_Efficiency_crystal[0][i] - 0.243*Total_Efficiency_crystal[3][i] - 0.013*Total_Efficiency_crystal[5][i] - 0.014*Total_Efficiency_crystal[14][i] - 0.518*Total_Efficiency_crystal[20][i] - 0.465*Total_Efficiency_crystal[23][i]),2))*Total_Efficiency_Err_crystal[0][i]),2)   +   pow((((0.243*Efficiency[6][i])/pow((1.0 - 0.154*Total_Efficiency_crystal[0][i] - 0.243*Total_Efficiency_crystal[3][i] - 0.013*Total_Efficiency_crystal[5][i] - 0.014*Total_Efficiency_crystal[14][i] - 0.518*Total_Efficiency_crystal[20][i] - 0.465*Total_Efficiency_crystal[23][i]),2))*Total_Efficiency_Err_crystal[3][i]),2)   +   pow((((0.013*Efficiency[6][i])/pow((1.0 - 0.154*Total_Efficiency_crystal[0][i] - 0.243*Total_Efficiency_crystal[3][i] - 0.013*Total_Efficiency_crystal[5][i] - 0.014*Total_Efficiency_crystal[14][i] - 0.518*Total_Efficiency_crystal[20][i] - 0.465*Total_Efficiency_crystal[23][i]),2))*Total_Efficiency_Err_crystal[5][i]),2)   +   pow((((0.014*Efficiency[6][i])/pow((1.0 - 0.154*Total_Efficiency_crystal[0][i] - 0.243*Total_Efficiency_crystal[3][i] - 0.013*Total_Efficiency_crystal[5][i] - 0.014*Total_Efficiency_crystal[14][i] - 0.518*Total_Efficiency_crystal[20][i] - 0.465*Total_Efficiency_crystal[23][i]),2))*Total_Efficiency_Err_crystal[14][i]),2)   +   pow((((0.518*Efficiency[6][i])/pow((1.0 - 0.154*Total_Efficiency_crystal[0][i] - 0.243*Total_Efficiency_crystal[3][i] - 0.013*Total_Efficiency_crystal[5][i] - 0.014*Total_Efficiency_crystal[14][i] - 0.518*Total_Efficiency_crystal[20][i] - 0.465*Total_Efficiency_crystal[23][i]),2))*Total_Efficiency_Err_crystal[20][i]),2)   +   pow((((0.465*Efficiency[6][i])/pow((1.0 - 0.154*Total_Efficiency_crystal[0][i] - 0.243*Total_Efficiency_crystal[3][i] - 0.013*Total_Efficiency_crystal[5][i] - 0.014*Total_Efficiency_crystal[14][i] - 0.518*Total_Efficiency_crystal[20][i] - 0.465*Total_Efficiency_crystal[23][i]),2))*Total_Efficiency_Err_crystal[23][i]),2) );

    //873keV-------------------------------------------------------------------------------------------------------
    Summing_Corrections[7][i] = ((1.0 + 0.024*((Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[14][i])/Peak_Efficiency_crystal[20][i]))*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.894*Total_Efficiency_crystal[18][i]));

    Efficiency_corrected[7][i] = Efficiency[7][i]/Summing_Corrections[7][i];

    Efficiency_Err_corrected[7][i] = sqrt( pow(((1.0/((1.0 + 0.024*(Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[14][i]/Peak_Efficiency_crystal[20][i]))*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.894*Total_Efficiency_crystal[18][i]))))*Efficiency_Err[7][i],2)   +   pow((0.024*Efficiency[7][i]*Peak_Efficiency_crystal[14][i]/(pow((1.0 + 0.024*Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[14][i]/Peak_Efficiency_crystal[20][i]),2)*Peak_Efficiency_crystal[20][i]*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.894*Total_Efficiency_crystal[18][i])))*Peak_Efficiency_Err_crystal[5][i],2)   +   pow((0.024*Efficiency[7][i]*Peak_Efficiency_crystal[5][i]/(pow((1.0 + 0.024*Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[14][i]/Peak_Efficiency_crystal[20][i]),2)*Peak_Efficiency_crystal[20][i]*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.894*Total_Efficiency_crystal[18][i])))*Peak_Efficiency_Err_crystal[14][i],2)   +   pow((0.024*Efficiency[7][i]*Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[14][i]/(pow((1.0 + 0.024*Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[14][i]/Peak_Efficiency_crystal[20][i]),2)*pow(Peak_Efficiency_crystal[20][i],2)*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.894*Total_Efficiency_crystal[18][i])))*Peak_Efficiency_Err_crystal[20][i],2)   +   pow(((0.282*Efficiency[7][i]/((1.0 + 0.024*(Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[14][i]/Peak_Efficiency_crystal[20][i]))*pow((1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.894*Total_Efficiency_crystal[18][i]),2))))*Total_Efficiency_Err_crystal[0][i],2)   +   pow(((0.455*Efficiency[7][i]/((1.0 + 0.024*(Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[14][i]/Peak_Efficiency_crystal[20][i]))*pow((1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.894*Total_Efficiency_crystal[18][i]),2))))*Total_Efficiency_Err_crystal[3][i],2)   +   pow(((0.894*Efficiency[7][i]/((1.0 + 0.024*(Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[14][i]/Peak_Efficiency_crystal[20][i]))*pow((1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.894*Total_Efficiency_crystal[18][i]),2))))*Total_Efficiency_Err_crystal[18][i],2) );

    //996keV---------------------------------------------------------------------------------------------------------
    Summing_Corrections[8][i] = ((1.0 + 0.507*((Peak_Efficiency_crystal[3][i]*Peak_Efficiency_crystal[20][i])/Peak_Efficiency_crystal[23][i]))*(1.0 - 0.894*Total_Efficiency_crystal[18][i]));

    Efficiency_corrected[8][i] = Efficiency[8][i]/Summing_Corrections[8][i];

    Efficiency_Err_corrected[8][i] = sqrt( pow((1.0/((1.0 + 0.507*((Peak_Efficiency_crystal[3][i]*Peak_Efficiency_crystal[20][i])/Peak_Efficiency_crystal[23][i]))*(1.0 - 0.894*Total_Efficiency_crystal[18][i])))*Efficiency_Err[8][i],2)   +   pow((0.507*Efficiency[8][i]*Peak_Efficiency_crystal[20][i]/(pow((1.0 + 0.507*((Peak_Efficiency_crystal[3][i]*Peak_Efficiency_crystal[20][i])/Peak_Efficiency_crystal[23][i])),2)*Peak_Efficiency_crystal[23][i]*(1.0 - 0.894*Total_Efficiency_crystal[18][i])))*Peak_Efficiency_Err_crystal[3][i],2)   +   pow((0.507*Efficiency[8][i]*Peak_Efficiency_crystal[3][i]/(pow((1.0 + 0.507*((Peak_Efficiency_crystal[3][i]*Peak_Efficiency_crystal[20][i])/Peak_Efficiency_crystal[23][i])),2)*Peak_Efficiency_crystal[23][i]*(1.0 - 0.894*Total_Efficiency_crystal[18][i])))*Peak_Efficiency_Err_crystal[20][i],2)   +   pow((0.507*Efficiency[8][i]*Peak_Efficiency_crystal[3][i]*Peak_Efficiency_crystal[20][i]/(pow((1.0 + 0.507*((Peak_Efficiency_crystal[3][i]*Peak_Efficiency_crystal[20][i])/Peak_Efficiency_crystal[23][i])),2)*pow(Peak_Efficiency_crystal[23][i],2)*(1.0 - 0.894*Total_Efficiency_crystal[18][i])))*Peak_Efficiency_Err_crystal[23][i],2)   +   pow((0.894*Efficiency[8][i]/((1.0 + 0.507*((Peak_Efficiency_crystal[3][i]*Peak_Efficiency_crystal[20][i])/Peak_Efficiency_crystal[23][i]))*pow((1.0 - 0.894*Total_Efficiency_crystal[18][i]),2)))*Efficiency_Err[8][i],2) );

    //1004keV----------------------------------------------------------------------------------------------------------
    Summing_Corrections[9][i] = ((1.0 + 0.221*((Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[19][i])/Peak_Efficiency_crystal[24][i]))*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.217*Total_Efficiency_crystal[11][i]));

    Efficiency_corrected[9][i] = Efficiency[9][i]/Summing_Corrections[9][i];

    Efficiency_Err_corrected[9][i] = sqrt( pow((1.0/((1.0 + 0.221*((Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[19][i])/Peak_Efficiency_crystal[24][i]))*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.217*Total_Efficiency_crystal[11][i])))*Efficiency_Err[9][i],2)   +   pow((0.221*Efficiency[9][i]*Peak_Efficiency_crystal[19][i]/(pow((1.0 + 0.221*((Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[19][i])/Peak_Efficiency_crystal[24][i])),2)*Peak_Efficiency_crystal[24][i]*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.217*Total_Efficiency_crystal[11][i])))*Peak_Efficiency_Err_crystal[5][i],2)   +   pow((0.221*Efficiency[9][i]*Peak_Efficiency_crystal[5][i]/(pow((1.0 + 0.221*((Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[19][i])/Peak_Efficiency_crystal[24][i])),2)*Peak_Efficiency_crystal[24][i]*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.217*Total_Efficiency_crystal[11][i])))*Peak_Efficiency_Err_crystal[19][i],2)   +   pow((0.221*Efficiency[9][i]*Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[19][i]/(pow((1.0 + 0.221*((Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[19][i])/Peak_Efficiency_crystal[24][i])),2)*pow(Peak_Efficiency_crystal[24][i],2)*(1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.217*Total_Efficiency_crystal[11][i])))*Peak_Efficiency_Err_crystal[24][i],2)   +   pow((0.282*Efficiency[9][i]/((1.0 + 0.221*((Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[19][i])/Peak_Efficiency_crystal[24][i]))*pow((1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.217*Total_Efficiency_crystal[11][i]),2)))*Total_Efficiency_Err_crystal[0][i],2)   +   pow((0.455*Efficiency[9][i]/((1.0 + 0.221*((Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[19][i])/Peak_Efficiency_crystal[24][i]))*pow((1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.217*Total_Efficiency_crystal[11][i]),2)))*Total_Efficiency_Err_crystal[3][i],2)   +   pow((0.217*Efficiency[9][i]/((1.0 + 0.221*((Peak_Efficiency_crystal[5][i]*Peak_Efficiency_crystal[19][i])/Peak_Efficiency_crystal[24][i]))*pow((1.0 - 0.282*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i] - 0.217*Total_Efficiency_crystal[11][i]),2)))*Total_Efficiency_Err_crystal[11][i],2) );

    //1274keV------------------------------------------------------------------------------------------------------------
    Summing_Corrections[10][i] = ((1.0 + 0.014*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[10][i])/Peak_Efficiency_crystal[27][i]))*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i]));

    Efficiency_corrected[10][i] = Efficiency[10][i]/Summing_Corrections[10][i];

    Efficiency_Err_corrected[10][i] = sqrt( pow((1.0/((1.0 + 0.014*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[10][i])/Peak_Efficiency_crystal[27][i]))*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Efficiency_Err[10][i],2)   +   pow((0.014*Efficiency[10][i]*Peak_Efficiency_crystal[10][i]/(pow((1.0 + 0.014*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[10][i])/Peak_Efficiency_crystal[27][i])),2)*Peak_Efficiency_crystal[27][i]*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[17][i],2)   +   pow((0.014*Efficiency[10][i]*Peak_Efficiency_crystal[17][i]/(pow((1.0 + 0.014*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[10][i])/Peak_Efficiency_crystal[27][i])),2)*Peak_Efficiency_crystal[27][i]*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[10][i],2)   +   pow((0.014*Efficiency[10][i]*Peak_Efficiency_crystal[10][i]*Peak_Efficiency_crystal[17][i]/(pow((1.0 + 0.014*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[10][i])/Peak_Efficiency_crystal[27][i])),2)*pow(Peak_Efficiency_crystal[27][i],2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[27][i],2)   +   pow((0.281*Efficiency[10][i]/((1.0 + 0.014*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[10][i])/Peak_Efficiency_crystal[27][i]))*pow((1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i]),2)))*Total_Efficiency_Err_crystal[0][i],2)   +   pow((0.455*Efficiency[10][i]/((1.0 + 0.014*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[10][i])/Peak_Efficiency_crystal[27][i]))*pow((1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i]),2)))*Total_Efficiency_Err_crystal[3][i],2) );

    //1596keV----------------------------------------------------------------------------------------------------------
    Summing_Corrections[11][i] = ((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i]))*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i]));

    Efficiency_corrected[11][i] = Efficiency[11][i]/Summing_Corrections[11][i];

    //good luck troubleshooting
    Efficiency_Err_corrected[11][i] = sqrt( pow((1.0/((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i]))*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Efficiency_Err[11][i],2)   +   pow((0.275*Efficiency[11][i]*Peak_Efficiency_crystal[22][i]/(Peak_Efficiency_crystal[28][i]*pow((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i])),2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[17][i],2)   +   pow((0.275*Efficiency[11][i]*Peak_Efficiency_crystal[17][i]/(Peak_Efficiency_crystal[28][i]*pow((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i])),2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[22][i],2)   +   pow((5.568*Efficiency[11][i]*Peak_Efficiency_crystal[18][i]/(Peak_Efficiency_crystal[28][i]*pow((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i])),2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[20][i],2)   +   pow((5.568*Efficiency[11][i]*Peak_Efficiency_crystal[20][i]/(Peak_Efficiency_crystal[28][i]*pow((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i])),2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[18][i],2)   +   pow((2.094*Efficiency[11][i]*Peak_Efficiency_crystal[11][i]/(Peak_Efficiency_crystal[28][i]*pow((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i])),2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[24][i],2)   +   pow((2.094*Efficiency[11][i]*Peak_Efficiency_crystal[24][i]/(Peak_Efficiency_crystal[28][i]*pow((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i])),2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[11][i],2)   +   pow((0.052*Efficiency[11][i]*Peak_Efficiency_crystal[9][i]/(Peak_Efficiency_crystal[28][i]*pow((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i])),2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[25][i],2)   +   pow((0.052*Efficiency[11][i]*Peak_Efficiency_crystal[25][i]/(Peak_Efficiency_crystal[28][i]*pow((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i])),2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[9][i],2)   +   pow((0.281*Efficiency[11][i]/((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i]))*pow((1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i]),2)))*Total_Efficiency_Err_crystal[0][i],2)   +   pow((0.455*Efficiency[11][i]/((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i]))*pow((1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i]),2)))*Total_Efficiency_Err_crystal[3][i],2)   +   pow((Efficiency[11][i]*(0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i]))/(pow((1.0 + 0.275*((Peak_Efficiency_crystal[17][i]*Peak_Efficiency_crystal[22][i])/Peak_Efficiency_crystal[28][i]) + 5.568*((Peak_Efficiency_crystal[20][i]*Peak_Efficiency_crystal[18][i])/Peak_Efficiency_crystal[28][i]) + 2.094*((Peak_Efficiency_crystal[24][i]*Peak_Efficiency_crystal[11][i])/Peak_Efficiency_crystal[28][i]) + 0.052*((Peak_Efficiency_crystal[25][i]*Peak_Efficiency_crystal[9][i])/Peak_Efficiency_crystal[28][i])),2)*(1.0 - 0.281*Total_Efficiency_crystal[0][i] - 0.455*Total_Efficiency_crystal[3][i])))*Peak_Efficiency_Err_crystal[28][i],2) );

    cout << "Efficiency[11][" << i << "]: " << Efficiency[11][i] << " +- " << Efficiency_Err[11][i] << "  " << "Efficiency_corrected[11][" << i << "]: " << Efficiency_corrected[11][i] << " +- " << Efficiency_Err_corrected[11][i] << "  " << "Summing_Corrections[11][" << i << "]: " << Summing_Corrections[11][i] << endl;

  } 

  cout << "------------------------------------" << endl;


  //adding up corrected efficiencies and temporary error terms for each energy for each region of interest
  //looping through crystals
  for(int i = 0; i < numcrystals; i++){
    //looping through energies
    for(int j = 0; j < numE; j++){
      //total efficiencies
      Sum_Efficiency_corrected[j] += Efficiency_corrected[j][i];
      Sum_Err_temp1[j] += pow(Efficiency_Err_corrected[j][i],2);
      //checking to see if the crystal is downstream
      if(std::find(downcrystals.begin(), downcrystals.end(), i) != downcrystals.end()) {
	//populating downstream efficiency and error vectors
	Sum_Eff_Down_corr[j] += Efficiency_corrected[j][i];
	Sum_Err_temp2[j] += pow(Efficiency_Err_corrected[j][i],2);
      }
      //checking to see if the crystal is upstream
      if(std::find(upcrystals.begin(), upcrystals.end(), i) != upcrystals.end()) {
	//populating upstream efficiency and error vectors
        Sum_Eff_Up_corr[j] += Efficiency_corrected[j][i];
	Sum_Err_temp3[j] += pow(Efficiency_Err_corrected[j][i],2);
      }
      //checking to see if the crystal is on the right
      if(std::find(rightcrystals.begin(), rightcrystals.end(), i) != rightcrystals.end()) {
	//populating right efficiency and error vectors
        Sum_Eff_Right_corr[j] += Efficiency_corrected[j][i];
	Sum_Err_temp4[j] += pow(Efficiency_Err_corrected[j][i],2);
      }
      //checking to see if the crystal is on the left
      if(std::find(leftcrystals.begin(), leftcrystals.end(), i) != leftcrystals.end()) {
	//populating left efficiency and error vectors
        Sum_Eff_Left_corr[j] += Efficiency_corrected[j][i];
	Sum_Err_temp5[j] += pow(Efficiency_Err_corrected[j][i],2);
      }
      //checking to see if the crystal is above
      if(std::find(abovecrystals.begin(), abovecrystals.end(), i) != abovecrystals.end()) {
	//populating above efficiency and error vectors
        Sum_Eff_Above_corr[j] += Efficiency_corrected[j][i];
	Sum_Err_temp6[j] += pow(Efficiency_Err_corrected[j][i],2);
      }
      //checking to see if the crystal is below
      if(std::find(belowcrystals.begin(), belowcrystals.end(), i) != belowcrystals.end()) {
	//populating below efficiency and error vectors
        Sum_Eff_Below_corr[j] += Efficiency_corrected[j][i];
	Sum_Err_temp7[j] += pow(Efficiency_Err_corrected[j][i],2);
      }

    }

  }
  
  //calculating summing corrected errors
  for(int i = 0; i < numE; i++){

    Sum_Eff_Err_Tot_corr[i] = sqrt(Sum_Err_temp1[i]);
    Sum_Eff_Err_Down_corr[i] = sqrt(Sum_Err_temp2[i]);
    Sum_Eff_Err_Up_corr[i] = sqrt(Sum_Err_temp3[i]);
    Sum_Eff_Err_Right_corr[i] = sqrt(Sum_Err_temp4[i]);
    Sum_Eff_Err_Left_corr[i] = sqrt(Sum_Err_temp5[i]);
    Sum_Eff_Err_Above_corr[i] = sqrt(Sum_Err_temp6[i]);
    Sum_Eff_Err_Below_corr[i] = sqrt(Sum_Err_temp7[i]);

  }

  
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //----------end of summing corrections---------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


  //printing some quantitative values for comparison purposes

  //comparing exp to sim in each region. Ideally, these ratios should be close to 1
  for(int i = 0; i < numE; i++){
    cout << "Ratio of exp total/sim total: " << energies[i] << "keV: " << Sum_Efficiency_corrected[i]/sum_eff_crystal_trunc[i] << endl;
  }
  cout << "-----------------------" << endl;
  for(int i = 0; i < numE; i++){
    cout << "Ratio of sim downstream/exp downstream: " << energies[i] << "keV: " << simdowneff_trunc[i]/Sum_Eff_Down_corr[i] << endl;
  }
  cout << "-----------------------" << endl;
  for(int i = 0; i < numE; i++){
    cout << "Ratio of sim upstream/exp upstream: " << energies[i] << "keV: " << simupeff_trunc[i]/Sum_Eff_Up_corr[i] << endl;
  }
  cout << "-----------------------" << endl;
  for(int i = 0; i < numE; i++){
    cout << "Ratio of sim left/exp left: " << energies[i] << "keV: " << simlefteff_trunc[i]/Sum_Eff_Left_corr[i] << endl;
  }
  cout << "-----------------------" << endl;
  for(int i = 0; i < numE; i++){
    cout << "Ratio of sim right/exp right: " << energies[i] << "keV: " << simrighteff_trunc[i]/Sum_Eff_Right_corr[i] << endl;
  }
  cout << "-----------------------" << endl;
  for(int i = 0; i < numE; i++){
    cout << "Ratio of sim above/exp above: " << energies[i] << "keV: " << simaboveeff_trunc[i]/Sum_Eff_Above_corr[i] << endl;
  }
  cout << "-----------------------" << endl;
  for(int i = 0; i < numE; i++){
    cout << "Ratio of sim below/exp below: " << energies[i] << "keV: " << simbeloweff_trunc[i]/Sum_Eff_Below_corr[i] << endl;
  }
  cout << "-----------------------" << endl;
 
  //ignoring differences in individual crystals, right/left and above/below should be
  //symmetric. That is, their ratios should be pretty close to 1 for each energy.
  //upstream and downstream comparison not included because those are not symmetric
  for(int i = 0; i < numE; i++){
    cout << "Ratio of sim right/sim left: " << energies[i] << "keV: " << simrighteff_trunc[i]/simlefteff_trunc[i] << endl;
  }
  cout << "-----------------------" << endl;
  for(int i = 0; i < numE; i++){
    cout << "Ratio of sim above/sim below: " << energies[i] << "keV: " << simaboveeff_trunc[i]/simbeloweff_trunc[i] << endl;
  }
  cout << "-----------------------" << endl;



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
  //DivideSquare function allows for uneven number of graphs to be displayed
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
  // Double_t Edouble_Full[simnumE];
  // Double_t Edouble_trunc[numE];


  // for(int i = 0; i < simnumE; i++){
  //   Edouble_Full[i] = (Double_t)Energy[i];
  // }







 

  //simulation TGraphs
  auto gre1 = new TGraphErrors(numE,energies,sum_eff_trunc,ex,sum_eff_err_trunc);
  auto gre2 = new TGraphErrors(numE,energies,sum_eff_crystal_trunc,ex,sum_eff_crystal_err_trunc);
  auto gre3 = new TGraphErrors(numE,energies,simupeff_trunc,ex,simupefferr_trunc);
  auto gre4 = new TGraphErrors(numE,energies,simdowneff_trunc,ex,simdownefferr_trunc);
  auto gre5 = new TGraphErrors(numE,energies,simaboveeff_trunc,ex,simaboveefferr_trunc);
  auto gre6 = new TGraphErrors(numE,energies,simbeloweff_trunc,ex,simbelowefferr_trunc);
  auto gre7 = new TGraphErrors(numE,energies,simlefteff_trunc,ex,simleftefferr_trunc);
  auto gre8 = new TGraphErrors(numE,energies,simrighteff_trunc,ex,simrightefferr_trunc);
  //experiment TGraphs
  //auto grtot = new TGraphErrors(12,energies,toteff,energieserr,totefferr);
  auto grtot = new TGraphErrors(12,energies,Sum_Efficiency_corrected,energieserr,Sum_Eff_Err_Tot_corr);
  //auto grdow = new TGraphErrors(12,energies,downeff,energieserr,downefferr);
  auto grdow = new TGraphErrors(12,energies,Sum_Eff_Down_corr,energieserr,Sum_Eff_Err_Down_corr);
  //auto grup = new TGraphErrors(12,energies,upeff,energieserr,upefferr);
  auto grup = new TGraphErrors(12,energies,Sum_Eff_Up_corr,energieserr,Sum_Eff_Err_Up_corr);
  //auto grrig = new TGraphErrors(12,energies,righteff,energieserr,rightefferr);
  auto grrig = new TGraphErrors(12,energies,Sum_Eff_Right_corr,energieserr,Sum_Eff_Err_Right_corr);
  //auto grlef = new TGraphErrors(12,energies,lefteff,energieserr,leftefferr);
  auto grlef = new TGraphErrors(12,energies,Sum_Eff_Left_corr,energieserr,Sum_Eff_Err_Left_corr);
  //auto grabo = new TGraphErrors(12,energies,aboveeff,energieserr,aboveefferr);
  auto grabo = new TGraphErrors(12,energies,Sum_Eff_Above_corr,energieserr,Sum_Eff_Err_Above_corr);
  //auto grbel = new TGraphErrors(12,energies,beloweff,energieserr,belowefferr);
  auto grbel = new TGraphErrors(12,energies,Sum_Eff_Below_corr,energieserr,Sum_Eff_Err_Below_corr);

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









    
}
