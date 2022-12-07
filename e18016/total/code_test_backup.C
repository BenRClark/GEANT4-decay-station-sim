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

  //creating root file to place graphs in
  TFile *f = new TFile("efficiencies.root","recreate");

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

  //sim currently has more energies that we looked for in experiment

  



  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //----------Plotting stuff---------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------





// Total Efficiency vs energy plot
  auto grtot = new TGraphErrors(12,energies,toteff,energieserr,totefferr);
  grtot->SetMarkerStyle(21);
  grtot->GetXaxis()->SetTitle("Energy [keV]");
  grtot->GetYaxis()->SetTitle("Efficiency [%]");
  auto Total = new TCanvas("Total","Total",800,600);
  grtot->Draw("ap");
  //ci->Draw("e3 same");
  Total->Write();
  //Total->SaveAs("efficiencyTOTAL.png");

// Downstream Efficiency vs energy plot
  auto grdow = new TGraphErrors(12,energies,downeff,energieserr,downefferr);
  grdow->SetMarkerStyle(21);
  grdow->GetXaxis()->SetTitle("Energy [keV]");
  grdow->GetYaxis()->SetTitle("Efficiency [%]");
  auto downstream = new TCanvas("downstream","downstream",800,600);
  grdow->Draw("ap");
  //ci->Draw("e3 same");
  downstream->Write();
  //Total->SaveAs("efficiencyTOTAL.png");

// Upstream Efficiency vs energy plot
  auto grup = new TGraphErrors(12,energies,upeff,energieserr,upefferr);
  grup->SetMarkerStyle(21);
  grup->GetXaxis()->SetTitle("Energy [keV]");
  grup->GetYaxis()->SetTitle("Efficiency [%]");
  auto upstream = new TCanvas("upstream","upstream",800,600);
  grup->Draw("ap");
  //ci->Draw("e3 same");
  upstream->Write();
  //Total->SaveAs("efficiencyTOTAL.png");

// Right Efficiency vs energy plot
  auto grrig = new TGraphErrors(12,energies,righteff,energieserr,rightefferr);
  grrig->SetMarkerStyle(21);
  grrig->GetXaxis()->SetTitle("Energy [keV]");
  grrig->GetYaxis()->SetTitle("Efficiency [%]");
  auto right = new TCanvas("right","right",800,600);
  grrig->Draw("ap");
  //ci->Draw("e3 same");
  right->Write();
  //Total->SaveAs("efficiencyTOTAL.png");

// Left Efficiency vs energy plot
  auto grlef = new TGraphErrors(12,energies,lefteff,energieserr,leftefferr);
  grlef->SetMarkerStyle(21);
  grlef->GetXaxis()->SetTitle("Energy [keV]");
  grlef->GetYaxis()->SetTitle("Efficiency [%]");
  auto left = new TCanvas("left","left",800,600);
  grlef->Draw("ap");
  //ci->Draw("e3 same");
  left->Write();
  //Total->SaveAs("efficiencyTOTAL.png");

// Above Efficiency vs energy plot
  auto grabo = new TGraphErrors(12,energies,aboveeff,energieserr,aboveefferr);
  grabo->SetMarkerStyle(21);
  grabo->GetXaxis()->SetTitle("Energy [keV]");
  grabo->GetYaxis()->SetTitle("Efficiency [%]");
  auto above = new TCanvas("above","above",800,600);
  grabo->Draw("ap");
  //ci->Draw("e3 same");
  above->Write();
  //Total->SaveAs("efficiencyTOTAL.png");

// Below Efficiency vs energy plot
  auto grbel = new TGraphErrors(12,energies,beloweff,energieserr,belowefferr);
  grbel->SetMarkerStyle(21);
  grbel->GetXaxis()->SetTitle("Energy [keV]");
  grbel->GetYaxis()->SetTitle("Efficiency [%]");
  auto below = new TCanvas("below","below",800,600);
  grbel->Draw("ap");
  //ci->Draw("e3 same");
  below->Write();
  //Total->SaveAs("efficiencyTOTAL.png");

 









    
}
