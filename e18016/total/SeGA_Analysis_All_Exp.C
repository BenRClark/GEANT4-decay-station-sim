#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TAxis.h"
#include "TROOT.h"

#include <iostream>
#include <fstream>

void SeGA_Analysis_All_Exp() {

  //Energies (Naming Convention)
  int Energy[12];  
  Energy[0] = 43; //in keV
  Energy[1] = 87; //in keV
  Energy[2] = 105; //in keV
  Energy[3] = 123; //in keV
  Energy[4] = 248; //in keV
  Energy[5] = 592; //in keV
  Energy[6] = 723; //in keV
  Energy[7] = 873; //in keV
  Energy[8] = 996; //in keV
  Energy[9] = 1005; //in keV
  Energy[10] = 1275; //in keV
  Energy[11] = 1596; //in keV
  //Energy[12] = 2152; //in keV

  //Actual Energy Values
  double Peak_Location[12];
  Peak_Location[0] = 42.8;  
  Peak_Location[1] = 86.5;  
  Peak_Location[2] = 105.3;  
  Peak_Location[3] = 123.1;  
  Peak_Location[4] = 247.7;  
  Peak_Location[5] = 591.8;  
  Peak_Location[6] = 723.3;  
  Peak_Location[7] = 873.2;  
  Peak_Location[8] = 996.2;  
  Peak_Location[9] = 1004.7;  
  Peak_Location[10] = 1274.5;  
  Peak_Location[11] = 1596.4;
  //Peak_Location[11] = 2151.5;

  int fit_high;
  int fit_low;
  
  //Number of Gammas in Simulation
  double nGamma = 1000000;
  //double nGamma = 500000;

  // Efficiencies [Energy][Det Number]
  double Peak_Efficiency[12][16];
  double Peak_Efficiency_Err[12][16];
  
  double Total_Efficiency[12][16];
  double Total_Efficiency_Err[12][16];

  //FIles
  TFile *f1t[12];
  for(int i=0; i<12; i++) {
    f1t[i] = new TFile(Form("/mnt/analysis/e16032/tho/simulation_rootfiles/geant4_output_rootfiles_calibrated_2/e16032_final_may_2020_outer_DL_1_micron/e16032_Eff_%dkeV_0.root",Energy[i]));

    //f1t[i] = new TFile(Form("/mnt/analysis/e16032/tho/simulation_rootfiles/geant4_output_rootfiles_calibrated_2/testing_DL_thickness_with_tentative_final_sim_conditions_0p575_active_layer_not_reduced_by_increment/e16032_Eff_%dkeV_0.root",Energy[i]));
  }

  //Fits & Histos
  TF1 *fSeGA1[12][16];
  TF1 *fSeGA[12][16];
  TH1D *hSeGA[12][16];

  TF1 *fgaussian[12][16];
  TF1 *fbackground[12][16];
  
  for(int i=0; i<12; i++) {
    for(int j=0; j<16; j++) {
      hSeGA[i][j] = (TH1D*) f1t[i]->Get(Form("hEnergyDepositSega_%d",j));

      //Fit Regions
      if(Energy[i]==43) {
	fit_high= 2;
	fit_low= 20;
      }      
      if(Energy[i]==87) {
	fit_high= 3;
	fit_low= 20;
      }     
      if(Energy[i]==105) {
	fit_high= 5;
	fit_low= 20;
      }
      if(Energy[i]==123) {
	fit_high= 7;
	fit_low= 20;
      }
      if(Energy[i]==248) {
	fit_high= 7;
	fit_low= 20;
      }
      if(Energy[i]==592) {
	fit_high= 9;
	fit_low = 20;
      }
      if(Energy[i]==723) {
	fit_high= 10;
	fit_low= 30;
      }
      if(Energy[i]==873) {
	fit_high= 10;
	fit_low= 30;
      }
      if(Energy[i]==996) {
	fit_high= 10;
	fit_low= 30;
      }
      if(Energy[i]==1005) {
	fit_high= 10;
	fit_low= 30;
      }
      if(Energy[i]==1275) {
	fit_high= 16;
	fit_low= 30;
      }
      if(Energy[i]==1596) {
	fit_high= 15;
	fit_low= 20;
      }
      //  if(Energy[i]==2152) {
      // 	fit_high= 15;
      // 	fit_low= 20;
      // }


      //fSeGA1[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus(0)+pol1(3)",Energy[i]-fit_low,Energy[i]+fit_high);
      fSeGA1[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus(0)",Energy[i]-fit_low,Energy[i]+fit_high);

      if(Energy[i]==1596||Energy[i]==2152){
	fSeGA[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus(0)",Energy[i]-fit_low,Energy[i]+fit_high);	
      }
      else{      
	fSeGA[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus(0)+pol1(3)",Energy[i]-fit_low,Energy[i]+fit_high);
      }

    fgaussian[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus",0.96*Energy[i],1.04*Energy[i]);
    fbackground[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"pol1",0.85*Energy[i],1.10*Energy[i]);

      // if(Energy[i]<100) {
      // 	fSeGA[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus(0)+pol1(3)",0.83*Energy[i],1.10*Energy[i]);
      // }
      // else if(Energy[i]<300) {
      // 	fSeGA[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus(0)+pol1(3)",0.83*Energy[i],1.06*Energy[i]);
      // }
      // else {
      // 	fSeGA[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus(0)+pol1(3)",0.84*Energy[i],1.02*Energy[i]);
      // }
    // if(Energy[i]<100) {
    //   fgaussian[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus",0.95*Energy[i],1.04*Energy[i]);
    // }
    // else if(Energy[i]<300) {
    //   fgaussian[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus",0.96*Energy[i],1.04*Energy[i]);
    // }
    // else {
    //   fgaussian[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"gaus",0.99*Energy[i],1.01*Energy[i]);
    // }

    // fbackground[i][j] = new TF1(Form("f%dkev_%d",Energy[i],j),"pol1",0.85*Energy[i],1.10*Energy[i]);
   
    }
  }
  
  TCanvas *c1t[12];  
  for(int i=0; i<12; i++) {
    c1t[i] = new TCanvas(Form("c_%d_keV",Energy[i]),Form("c_%d_keV",Energy[i]),1024,768);
    c1t[i]->Divide(4,4);

  }
  for(int i=0; i<12; i++) {
    for(int j=0; j<16; j++) {
      c1t[i]->cd(j+1);
      //gPad->SetLogy(i+1);
      hSeGA[i][j]->GetXaxis()->SetRangeUser(0.8*Energy[i],1.10*Energy[i]);
      hSeGA[i][j]->Draw();
    
      fSeGA1[i][j]->SetParLimits(0,1,1000000);
      fSeGA1[i][j]->SetParameter(0,100);
      
      //fSeGA1[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      if(Energy[i]== 43){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[0]);
      }
      if(Energy[i]== 87){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[1]);
      }
      if(Energy[i]== 105){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[2]);
      }
      if(Energy[i]== 123){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[3]);
      }
      if(Energy[i]== 248){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[4]);
      }
      if(Energy[i]== 592){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[5]);
      }
      if(Energy[i]== 723){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[6]);
      }
      if(Energy[i]== 873){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[7]);
      }
      if(Energy[i]== 996){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[8]);
      }
      if(Energy[i]== 1005){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[9]);
      }
      if(Energy[i]== 1275){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[10]);
      }
      if(Energy[i]== 1596){
	fSeGA1[i][j]->FixParameter(1,Peak_Location[11]);
      }
      
      fSeGA1[i][j]->SetParLimits(2,0,5);
      fSeGA1[i][j]->SetParameter(2,5);

      // if(Energy[i]==1596){
      // 	fSeGA1[i][j]->SetParLimits(2,0,5);
      // 	fSeGA1[i][j]->SetParameter(2,3);
      // }
      
      hSeGA[i][j]->Fit(fSeGA[i][j],"RMNQ");
      
      // fSeGA[i][j]->SetParLimits(0,1,1000000);
      // fSeGA[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      // fSeGA[i][j]->SetParLimits(2,0,5);
      
      //hSeGA[i][j]->Fit(fSeGA[i][j],"RQ");
      
      // fSeGA[i][j]->SetParLimits(0,1,1000000);
      // fSeGA[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      // fSeGA[i][j]->SetParLimits(2,0,5);
      
      // hSeGA[i][j]->Fit(fSeGA[i][j],"RQ");
      
      // fSeGA[i][j]->SetParLimits(0,1,1000000);
      // fSeGA[i][j]->SetParLimits(1,0.98*Energy[i],1.02*Energy[i]);
      // fSeGA[i][j]->SetParLimits(2,0,5);
      
      // hSeGA[i][j]->Fit(fSeGA[i][j],"RQ+");

      fSeGA[i][j]->SetParameter(0,fSeGA1[i][j]->GetParameter(0));
      fSeGA[i][j]->SetParameter(1,fSeGA1[i][j]->GetParameter(1));
      fSeGA[i][j]->SetParameter(2,fSeGA1[i][j]->GetParameter(2));
      //fSeGA[i][j]->SetParameter(3,fSeGA1[i][j]->GetParameter(3));
      //fSeGA[i][j]->SetParameter(4,fSeGA1[i][j]->GetParameter(4));

      hSeGA[i][j]->Fit(fSeGA[i][j],"RMQL+");
      fSeGA[i][j]->SetLineColor(kViolet);
      fSeGA[i][j]->Draw("SAME");
      
      gPad->Modified();
      gPad->Update();

      fgaussian[i][j]->SetParameter(0,fSeGA[i][j]->GetParameter(0));
      fgaussian[i][j]->SetParameter(1,fSeGA[i][j]->GetParameter(1));
      fgaussian[i][j]->SetParameter(2,fSeGA[i][j]->GetParameter(2));
      fgaussian[i][j]->SetLineColor(kViolet);
      //fgaussian[i][j]->Draw("SAME");

      gPad->Modified();
      gPad->Update();
  
      //fbackground[i][j]->SetParameter(0,fSeGA[i][j]->GetParameter(3));
      //fbackground[i][j]->SetParameter(1,fSeGA[i][j]->GetParameter(4));
      //fbackground[i][j]->SetLineColor(kGreen);
      //fbackground[i][j]->Draw("SAME");


      // double area = 2.506628*(fSeGA[i][j]->GetParameter(0))*(fSeGA[i][j]->GetParameter(2));
      // double area_err = area*pow((pow((fSeGA[i][j]->GetParError(0)/fSeGA[i][j]->GetParameter(0)),2)+pow((fSeGA[i][j]->GetParError(2)/fSeGA[i][j]->GetParameter(2)),2)),0.5);

      double area = 2.506628*(fSeGA[i][j]->GetParameter(0))*fabs((fSeGA[i][j]->GetParameter(2)));
      double area_err = area*pow((pow((fSeGA[i][j]->GetParError(0)/fSeGA[i][j]->GetParameter(0)),2)+pow((fabs(fSeGA[i][j]->GetParError(2))/fSeGA[i][j]->GetParameter(2)),2)),0.5);
      
      double eff = 100.0*area/(1.0*nGamma);
      double eff_err = 100.0*area_err/(1.0*nGamma);

      double total_eff = 0;
      for(int k=2; k<Energy[i]*1.10; k++) {
	total_eff += hSeGA[i][j]->GetBinContent(k);
      }
      double total_eff_err = pow(total_eff,0.5)/(1.0*nGamma);
      total_eff = total_eff/(1.0*nGamma);
      

      cout<<i<<"  "<<j<<"  "<<eff<<" +/- "<<eff_err<<"  "<<100*total_eff<<" +/- "<<100*total_eff_err<<endl;

      // Efficiencies [Energy][Det Number]
      if(j<=7 && (Energy[i]==43 || Energy[i] == 87 || Energy[i] == 105)) {
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
  
  ofstream data;
  data.open("/mnt/analysis/e16032/tho/simulation_rootfiles/geant4_output_rootfiles_calibrated_2/e16032_final_may_2020_outer_DL_1_micron/Sim_Eff.txt");
  //data.open("/mnt/analysis/e16032/tho/simulation_rootfiles/geant4_output_rootfiles_calibrated_2/testing_DL_thickness_with_tentative_final_sim_conditions_0p575_active_layer_not_reduced_by_increment/Sim_Eff_2.txt");
  
  for(int j=0; j<16; j++) {
    data<<"SeGA"<<j<<"\n";
    for(int i=0; i<12; i++) {
      data<<Peak_Efficiency[i][j]<<"  "<< Peak_Efficiency_Err[i][j]<<"  "<<Total_Efficiency[i][j]<<"  "<<Total_Efficiency_Err[i][j]<<"\n";	  
    }
    data<<"\n"<<"\n";	
  }
  
  cout<<"Done!"<<endl;
}

// double Sigma[12][16];
// double Sigma_Error[12][16];
// double Mean[12][16];
// double Mean_Error[12][16];
// double Height[12][16];
// double Height_Error[12][16];




//       //Get the Parameters of the fit
//       Sigma[i][j] = fabs(fgaussian[i][j]->GetParameter(2));
//       Mean[i][j] = fgaussian[i][j]->GetParameter(1);
//       Height[i][j] = fgaussian[i][j]->GetParameter(0);
      
//       //Get the associated Errors
//       Sigma_Error[i][j] = fgaussian[i][j]->GetParError(2);
//       Mean_Error[i][j] = fgaussian[i][j]->GetParError(1);
//       Height_Error[i][j] = fgaussian[i][j]->GetParError(0);

      
//     }
//  }

  //     for(int k=5; k<Energy[i]+50; k++) {	
  // 	total_eff+= SeGA[i][j]->GetBinContent(k);
  //     }

  //   // Missing detector 6 (5 counting from 0 in data)
  //   if(i==5 || i==12 || (Energy==42 && i== 2)) {
  //     cout<<"Special Condition"<<endl;
  //   }
  //   else {
  //     cout<<"SeGA :"<<i<<" Area: "<<area<<"  Eff: "<<eff<<endl;
  //     nGammaPhotoPeak+= area;
  //     totaleff+= eff; 
  //     total_eff = total_eff/(0.01*nGamma);		    
  //     cout<<"SeGA :"<<i<<" Total_eff: "<<total_eff<<"  %"<<endl;
  //     total_array_eff += total_eff;
  //   }
    
  //   if(i==7) {
  //     cout<<"Downstream Ring: "<<totaleff<<endl;
  //     totaleff = 0;
  //     cout<<"Total Downstream Ring: "<<total_array_eff<<" %"<<endl;	
  //     total_arrray_eff = 0;	
  //   }

  //   total_eff = 0;
  // }

  // cout<<"Upstream Ring: "<<totaleff<<endl;
  // cout<<"Total Upstream Ring: "<<total_array_eff<<endl;
  // cout<<"There were "<<nGammaPhotoPeak<<" Gammas detected in the photopeak"<<endl;
  // cout<<"Efficiency: "<<nGammaPhotoPeak*100.0/nGamma<< "%"<<endl;

//}
