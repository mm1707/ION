#include <iostream>
#include <vector>

#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TCutG.h"

using namespace std;
auto background(double *x, double *par) {
    double func=par[0]*x[0]*x[0]*x[0]*x[0]+par[1]*x[0]*x[0]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]+par[4];
    if(func>0)
      return func;
    else
        return 0.0;
    }

auto gauss(double *x, double *par){
        return par[0]*TMath::Gaus(x[0],par[1],par[2]);
    }

auto signal(double *x, double *par) {
      return gauss(x, par);
    }

auto combined(double* x, double* par){
        return background(x,par) + signal(x,&par[5]);
    }


void analysis(){
    TFile* inFile = new TFile("results/step0/AnalysisResults.root"); //getting the root file
    vector <TH1F*> histograms;
    vector <float> integrals, mi;
    vector <TF1*> masses, backgrounds, signals;

    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt1"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt2"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt3"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt4"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt5"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt6"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt7"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt8"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt9"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt10"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt11"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt12"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt13"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt14"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt15"));
    histograms.push_back(inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambdaPt16"));

    TH1F* hMass = inFile->Get<TH1F>("strangeness_tutorial/Lambda/hMassLambda");
    TCanvas* c1=new TCanvas("ptLambda","Histogram of pT Lambda", 2000, 1000);
    TCanvas* c2=new TCanvas("Lambda from Pt","Histogram of Lambda particles in given Pt", 2000, 1000);
    c1->Divide(4, 4);


    for (size_t i = 0; i < histograms.size(); i++ ){
        histograms[i]->GetXaxis()->SetRangeUser(1.05, 1.2);
        TF1* fittedMass=new TF1("signal and background fitted", combined, 1.05, 1.2, 8);
        fittedMass->SetParameters(0, 0, 0, 0, 0, 3000, 1.115, -2.22915e-03);
        histograms[i]->Fit(fittedMass, "MER0");
        double* fparsMass = fittedMass -> GetParameters();
        TF1* backgroundFromFitMass=new TF1("signal and background fitted", background, 1.05, 1.2, 5); //separating function of background and signal after fitting
        backgroundFromFitMass->SetParameters(fparsMass[0], fparsMass[1], fparsMass[2], fparsMass[3], fparsMass[4]);
        TF1* signalFromFitMass=new TF1("signal and background fitted", signal, 1.05, 1.2, 3);
        signalFromFitMass->SetParameters(fparsMass[5], fparsMass[6], fparsMass[7]);
        backgroundFromFitMass->SetLineColor(kGreen);
        fittedMass->SetLineColor(kBlue);
        masses.push_back(fittedMass);
        signals.push_back(signalFromFitMass);
        backgrounds.push_back(backgroundFromFitMass);
        
        float integral=0;
        float minus=0;
        for (int j=0; j<histograms[i]->GetNbinsX(); j++){
            
            
            if(signalFromFitMass->Eval(histograms[i]->GetBinCenter(j))>1)
            {
               integral+=signalFromFitMass->Eval(histograms[i]->GetBinCenter(j));
                minus+=backgroundFromFitMass->Eval(histograms[i]->GetBinCenter(j));
            
            }
      
            
        }
        integrals.push_back(integral-minus);
     
        mi.push_back(minus);
    }
    float pt=0.5;
    TGraphErrors* hLambdaFromPt =new TGraphErrors();
    TH1F* hbackground =new TH1F("hbackground", "Background Particles; Pt ; number of lambda particles", 30, 0, 3);
    hbackground->SetLineColor(kRed);
    for(int i =0; i<integrals.size(); i++){
        hLambdaFromPt->AddPoint(pt, integrals[i]);
        hLambdaFromPt->SetPointError(i, 0.125/2, TMath::Sqrt(integrals[i]));
        char name[60];
        gStyle->SetTitleFontSize(0.07);
        sprintf(name, "Histogram of Minv of lambda particles with pT in range (%.3f, %.3f) GeV/c", pt, pt+0.125);
        hbackground->Fill(pt, mi[i]);
        
        histograms[i]->SetTitle(name);
        pt+=0.125;
        c1->cd(i+1);
        histograms[i]->Draw();
        backgrounds[i]->Draw("same");
        signals[i]->Draw("same");
        masses[i]->Draw("same");
        
    }
    for (int j=0; j<hbackground->GetNbinsX(); j++){
        
        hbackground->SetBinError(j, 0);
    }

    TCanvas* c11=new TCanvas("MLambda1","Histograms of invariant mass of Lambda", 2000, 1000);
    c11->Divide(2,2);
    TCanvas* c22=new TCanvas("MLambda2","Histograms of invariant mass of Lambda", 2000, 1000);
    c22->Divide(2,2);
    TCanvas* c3=new TCanvas("MLambda3","Histograms of invariant mass of Lambda", 2000, 1000);
    c3->Divide(2,2);
    TCanvas* c4=new TCanvas("MLambda4","Histograms of invariant mass of Lambda", 2000, 1000);
    c4->Divide(2,2);
   

    hLambdaFromPt->SetMarkerStyle(21);
    hLambdaFromPt->SetTitle("Number of lambda particles for given pT;pT (GeV/c); number of lambda particles");
    c2->cd();
    hLambdaFromPt->Draw("AP");

    for(int i=0; i<4; i++){
        c11->cd(i+1);
        histograms[i]->Draw();
        backgrounds[i]->Draw("same");
        signals[i]->Draw("same");
        masses[i]->Draw("same");
        c22->cd(i+1);
        histograms[i+4]->Draw();
        backgrounds[i+4]->Draw("same");
        signals[i+4]->Draw("same");
        masses[i+4]->Draw("same");
        c3->cd(i+1);
        histograms[i+8]->Draw();
        backgrounds[i+8]->Draw("same");
        signals[i+8]->Draw("same");
        masses[i+8]->Draw("same");
        c4->cd(i+1);
        histograms[i+12]->Draw();
        backgrounds[i+12]->Draw("same");
        signals[i+12]->Draw("same");
        masses[i+12]->Draw("same");
        
    }
    
}
