#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include <iostream>
#include "TF2.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "RooUniform.h"
#include "RooClassFactory.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "TGraphErrors.h"
#include "RooDataHist.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "RooNLLVar.h"
#include "RooMCStudy.h"
#include "RooChi2Var.h"
#include "RooAddPdf.h"
#include "RooMinimizer.h"
#include <cmath>
#include "string.h"
#include <iostream>

using namespace RooFit ;

double myPDF(double *x, double *par){ //to define TF2 obj
  return par[3]*(3./(4.*M_PI))*(0.5*(1.-par[0]) + (0.5)*(3.*par[0]-1)*cos(x[0])*cos(x[0]) - par[1]*sin(x[0])*sin(x[0])*cos(2.*x[1])- sqrt(2.)*par[2]*sin(2.*x[0])*cos(x[1]));
}


int main(){
    
    //-------------------------------
    //----Initialization-------------
    //-------------------------------

    gStyle->SetOptStat(0);
    //Defining Real variables according to RooFit
    //Angles physical defined in [0-pi], [0, 2*pi]
    RooRealVar  t("t","t",0,M_PI) ;
    RooRealVar  p("p","p", 0, 2*M_PI) ;
    
    //Defining Parameters according to RooFit.
    //Initial value is the 3rd argument, 4th and 5th arg
    //are the range of existance of the parameter.
    //RooRealVar alpha("alpha", "alpha", 0.65, -1,1);
    RooRealVar alpha("alpha", "alpha", 0.65, 0.62, 0.67);

    //RooRealVar beta("beta", "beta", 0.06, -1, 1);
    RooRealVar beta("beta", "beta", 0.06, 0.05, 0.075);
    RooRealVar gamma("gamma", "gamma", -0.18, -.19, -0.16);
    
    //ClassFactory generates a pdf instance and normalize it.
    RooAbsPdf* genpdf = RooClassFactory::makePdfInstance("GenPdf","(3./(4.*M_PI))*(0.5*(1.-alpha) + (0.5)*(3.*alpha-1)*cos(t)*cos(t) - beta*sin(t)*sin(t)*cos(2.*p)- sqrt(2.)*gamma*sin(2.*t)*cos(p))",RooArgSet(t,p,alpha, beta, gamma)) ;
    
    //plotting the pdf in a binned frame
    TH2F* h2_pdf = new TH2F("F(#theta, #phi) distribution", "F(#theta, #phi) distribution", 100, 0, M_PI, 100, 0, 2*M_PI);
    genpdf->fillHistogram(h2_pdf, RooArgList(t,p));
    h2_pdf->GetXaxis()->SetTitle("#theta [rad]");
    h2_pdf->GetYaxis()->SetTitle("#phi [rad]");
    h2_pdf->GetXaxis()->SetTitleOffset(2);
    h2_pdf->GetYaxis()->SetTitleOffset(2);
    h2_pdf->GetZaxis()->SetTitle("F(#theta, #phi)");
    h2_pdf->GetZaxis()->SetTitleOffset(1.5);
    TCanvas* c_pdf = new TCanvas("c_pdf", "c_pdf", 1000,1000,1000,800);
    c_pdf->SetLeftMargin(2);
    h2_pdf->Draw("surf");
    c_pdf->Draw();
    c_pdf->SaveAs("./graphs/real_pdf.png", "png");

    // Generate a toy MC dataset from the interpreted p.d.f
    //data2 will contain unbinned data.
    RooDataSet* data = genpdf->generate(RooArgSet(t,p),50000) ;
    
    //visualizing the events generated through a ROOT TH2F.
    TH2F* dh2 = new TH2F("h2", "h2", 100, 0, M_PI, 100, 0, 2*M_PI);
    data->fillHistogram(dh2, RooArgList(t,p));
    TCanvas*c = new TCanvas("c", "c", 1000,1000,1000,800);
    dh2->SetTitle("Binned #theta, #phi distribution");
    dh2->GetXaxis()->SetTitle("#theta [rad]");
    dh2->GetYaxis()->SetTitle("#phi [rad]");
    dh2->Draw("colz");
    c->Draw();
    c->SaveAs("./graphs/generated_distribution.png", "png");

    //-------------------------------
    //----Unbinned ML fit-------------
    //-------------------------------

    //Unbinned maximum likelihood fit:
    RooFitResult* r_ML = genpdf->fitTo(*data, Save()) ;

    //Retireve statistics:
    auto alpha_est = alpha.getVal();
    auto beta_est = beta.getVal();
    auto gamma_est = gamma.getVal();

    TH2F* h2_pdf_ml = new TH2F("F(#theta, #phi) distribution ml", "F(#theta, #phi) distribution", 100, 0, M_PI, 100, 0, 2*M_PI);
    h2_pdf_ml->SetLineColor(kRed);
    genpdf->fillHistogram(h2_pdf_ml, RooArgList(t,p));

    //defining a TF2 object with estimated max LL parameters
    const int npar = 4;
    double dist_params[npar] = {alpha_est, beta_est, gamma_est, 5.};
    TF2* myPdf = new TF2("myPdf", myPDF, 0, M_PI, 0, 2*M_PI, npar);
    myPdf->SetParameters(dist_params);
    myPdf->SetParNames("#alpha", "#beta", "#gamma", "#scaling");

    gStyle->SetOptStat(1111);
    TCanvas* from_ML = new TCanvas("from_ML", "from_ML", 1000,1000,1000,800);
    from_ML->SetTopMargin(0.15);
    dh2->SetTitle(Form("#splitline{Unbinned ML fit: #alpha = %.3f#pm%.3f, #beta = %.3f#pm%.3f, #gamma = %.3f#pm%.3f}{minimized FCN value: %.5f, EDM: %5f}",alpha_est, alpha.getError(), beta_est, beta.getError(), gamma_est, gamma.getError(),r_ML->minNll() ,r_ML->edm() ));
    dh2->GetXaxis()->SetTitle("#theta [rad]");
    dh2->GetYaxis()->SetTitle("#phi [rad]");
    dh2->Draw("lego");
    //myPdf->Draw("same cont1");
    h2_pdf_ml->Draw("surf3 same");
    from_ML->Draw();
    from_ML->SaveAs("./graphs/result_ML_unbinned.pdf", "pdf");

    TCanvas* from_ML_2 = new TCanvas("from_ML_2", "from_ML_2", 1000,1000,1000,800);
    from_ML_2->SetTopMargin(0.15);
    dh2->SetTitle(Form("#splitline{Unbinned ML fit: #alpha = %.3f#pm%.3f, #beta = %.3f#pm%.3f, #gamma = %.3f#pm%.3f}{minimized FCN value: %.5f, EDM: %5f}",alpha_est, alpha.getError(), beta_est, beta.getError(), gamma_est, gamma.getError(),r_ML->minNll() ,r_ML->edm() ));
    dh2->GetXaxis()->SetTitle("#theta [rad]");
    dh2->GetXaxis()->SetTitleOffset(2);
    dh2->GetYaxis()->SetTitle("#phi [rad]");
    dh2->GetYaxis()->SetTitleOffset(2);
    dh2->Draw("lego");
    myPdf->Draw("same cont1");
    from_ML_2->Draw();
    from_ML_2->SaveAs("./graphs/result_ML_unbinned_2.pdf", "pdf");


 /* 
    //-------------------------------
    //----Likelihood plots-------------
    //-------------------------------

    //obtaining Likelihhod from model.
    RooAbsReal* nll = genpdf->createNLL(*data, NumCPU(4));

    //instantiating minimizer for likelihood
    RooMinimizer m(*nll) ;
    m.migrad();

    //plotting likelihood for parameter alpha
    RooPlot* frame1 = alpha.frame(Title("Likelihood #alpha")) ;
    nll->plotOn(frame1,ShiftToZero()) ;

    //plotting likelihood for parameter beta
    RooPlot* frame2 = beta.frame(Title("Likelihood #beta")) ;
    nll->plotOn(frame2,ShiftToZero()) ;

    //plotting likelihood for parameter gamma
    RooPlot* frame3= gamma.frame(Title("Likelihood #gamma")) ;
    nll->plotOn(frame3,ShiftToZero()) ;
    
    TCanvas*c1 = new TCanvas("c1", "c1", 1000, 1000, 1000, 800);
    frame1->Draw();
    c1->Draw();
    c1->SaveAs("./graphs/likelihood_alpha.png", "png");

    TCanvas*c2 = new TCanvas("c2", "c2", 1000, 1000, 1000, 800);
    frame2->Draw();
    c2->Draw();
    c2->SaveAs("./graphs/likelihood_beta.png", "png");

    TCanvas*c3 = new TCanvas("c3", "c3", 1000, 1000, 1000, 800);
    frame3->Draw();
    c3->Draw();
    c3->SaveAs("./graphs/likelihood_gamma.png", "png");


    //-------------------------------
    //----Minuit Cont plot-----------
    //-------------------------------

    //plotting Minuit contour plot.
    RooPlot* frame_cont = m.contour(alpha, beta, 1, 2, 3);
    
    TCanvas* ct = new TCanvas("ct", "ct", 1000, 1000, 1000, 800);
    gPad->SetLeftMargin(0.15);
    frame_cont->GetYaxis()->SetTitleOffset(2);
    frame_cont->GetYaxis()->SetTitle("#gamma");
    frame_cont->GetXaxis()->SetTitle("#beta");
    frame_cont->Draw();
    ct->Draw();
    ct->SaveAs("./graphs/cont_gb.png", "png");


    //Creating many montecarlo events with this setup
    RooMCStudy mcs(*genpdf, RooArgSet(t,p));
    mcs.generateAndFit(1,50000) ;  // generate & fit 100 samples of 500 events each

    TCanvas* mc_can = new TCanvas("mc_can", "mc_can", 1000,1000,1000,800);
    RooPlot* mc_frame1 =  mcs.plotParam(alpha) ; // Draw distribution of parameter
    mc_frame1->Draw() ;
    mc_can->Draw();
    mc_can->SaveAs("./graphs/alpha_dist.png", "png");

    RooPlot* mc_frame2 =  mcs.plotError(alpha) ; // Draw distribution of parameter error
    mc_frame2->Draw() ;
    mc_can->Draw();
    mc_can->SaveAs("./graphs/alpha_error_dist.png", "png");

    RooPlot* mc_frame3 =  mcs.plotPull(alpha) ; // Draw distribution of parameter pull
    mc_frame3->Draw() ;
    mc_can->Draw();
    mc_can->SaveAs("./graphs/alpha_pull.png", "png");
*/


//----------------------
//---Many MC analysis---
//----------------------

/* 
double num = 50000.;
std::vector<double> num_vec;
std::vector<double> err_num_vec;
std::vector<double> al;
std::vector<double> bel;
std::vector<double> gam;
std::vector<double> al_er;
std::vector<double> bel_er;
std::vector<double> gam_er;

for(int i = 0; i < 5; i++){

  RooMCStudy mcs(*genpdf, RooArgSet(t,p));
  mcs.generateAndFit(200,num, kTRUE) ;
  RooDataSet  fitParData = mcs.fitParDataSet();
  RooDataSet all_alpha_est("all_alpha_est", "all_alpha_est",RooArgSet(alpha));
  RooDataSet all_beta_est("all_beta_est", "all_beta_est",RooArgSet(beta));
  RooDataSet all_gamma_est("all_gamma_est", "all_gamma_est",RooArgSet(gamma));

  for(int i = 0; i < 100; i++){
      const RooArgSet* f = fitParData.get(i);
      alpha = f->getRealValue("alpha");
      beta = f->getRealValue("beta");
      gamma = f->getRealValue("gamma");
      all_alpha_est.add(RooArgSet(alpha));
      all_beta_est.add(RooArgSet(beta));
      all_gamma_est.add(RooArgSet(gamma));
    }

  RooRealVar meana("mean","mean",0.65,0.5,1) ;
  RooRealVar sigmaa("sigma","sigma",0.001,0.0001,10) ;
  RooGaussian gaussa("gauss","gauss",alpha, meana, sigmaa) ;

  gaussa.fitTo(all_alpha_est);

  RooRealVar meanb("mean","mean",0.06,0.,1) ;
  RooRealVar sigmab("sigma","sigma",0.001,0.0001,10) ;
  RooGaussian gaussb("gauss","gauss",beta, meanb, sigmab) ;

  gaussb.fitTo(all_beta_est);
  
  RooRealVar meang("mean","mean",-0.18,-0.3,0.) ;
  RooRealVar sigmag("sigma","sigma",0.001,0.0001,10) ;
  RooGaussian gaussg("gauss","gauss",gamma, meang, sigmag) ;

  gaussg.fitTo(all_gamma_est);

  al.push_back(abs(meana.getVal()-0.65));
  bel.push_back(abs(meanb.getVal()-0.06));
  gam.push_back(abs(meang.getVal()+0.18));

  al_er.push_back(sigmaa.getVal()/sqrt(100));
  bel_er.push_back(sigmab.getVal()/sqrt(100));
  gam_er.push_back(sigmag.getVal()/sqrt(100));

  num_vec.push_back(num);
  err_num_vec.push_back(0.);

  num += 20000.;

}

TCanvas* c_eval = new TCanvas("c_eval", "c_eval", 10000, 1000, 1000, 800);

TGraphErrors* g_a = new TGraphErrors(al_er.size(), &num_vec[0], &al[0], &err_num_vec[0], &al_er[0]);
TGraphErrors* g_b = new TGraphErrors(bel_er.size(), &num_vec[0], &bel[0], &err_num_vec[0], &bel_er[0]);
TGraphErrors* g_g = new TGraphErrors(gam_er.size(), &num_vec[0], &gam[0], &err_num_vec[0], &gam_er[0]);
g_a->SetTitle("Estimation of #alpha varying MC event number");
g_a->GetXaxis()->SetTitle("MC event number");
g_a->GetYaxis()->SetTitle("#alpha est.");
g_a->SetMarkerColor(4);
g_a->SetMarkerStyle(21);
g_a->Draw("AP");
c_eval->Draw();
c_eval->SaveAs("./graphs/prova_many_a.png", "png");

g_b->SetTitle("Estimation of #beta varying MC event number");
g_b->GetXaxis()->SetTitle("MC event number");
g_b->GetYaxis()->SetTitle("#beta est.");
g_b->SetMarkerColor(4);
g_b->SetMarkerStyle(21);
g_b->Draw("AP");
c_eval->Draw();
c_eval->SaveAs("./graphs/prova_many_b.png", "png");

g_g->SetTitle("Estimation of #gamma varying MC event number");
g_g->GetXaxis()->SetTitle("MC event number");
g_g->GetYaxis()->SetTitle("#gamma est.");
g_g->SetMarkerColor(4);
g_g->SetMarkerStyle(21);
g_g->Draw("AP");
c_eval->Draw();
c_eval->SaveAs("./graphs/prova_many_g.png", "png");

*/

/* 
    std::cout<< "----->MCS<--------" << std::endl;
    
    //defiying a RooMCStudy object
    RooMCStudy mcs(*genpdf, RooArgSet(t,p));
    //generate & fit 1000 samples of 50000 events each
    mcs.generateAndFit(10,50000, kTRUE) ;  

    //defininying gauss to fit
    RooDataSet  fitParData = mcs.fitParDataSet();
    RooDataSet all_alpha_est("all_alpha_est", "all_alpha_est",RooArgSet(alpha));
    RooDataSet all_beta_est("all_beta_est", "all_beta_est",RooArgSet(beta));
    RooDataSet all_gamma_est("all_gamma_est", "all_gamma_est",RooArgSet(gamma));

    for(int i = 0; i < 100; i++){
      const RooArgSet* f = fitParData.get(i);
      alpha = f->getRealValue("alpha");
      beta = f->getRealValue("beta");
      gamma = f->getRealValue("gamma");
      all_alpha_est.add(RooArgSet(alpha));
      all_beta_est.add(RooArgSet(beta));
      all_gamma_est.add(RooArgSet(gamma));
      //f->Print("v"); 
    }

    RooRealVar meana("mean","mean",0.65,0.5,1) ;
    RooRealVar sigmaa("sigma","sigma",0.001,0.0001,10) ;
    RooGaussian gaussa("gauss","gauss",alpha, meana, sigmaa) ;

    gaussa.fitTo(all_alpha_est);
    //-------
    //--Alpha
    //-------

    TCanvas* mc_can = new TCanvas("mc_can", "mc_can", 1000,1000,1000,800);
    RooPlot* mc_frame1 =  mcs.plotParam(alpha) ; 
    mc_frame1->SetTitle("#alpha MC estimations");
    gaussa.plotOn(mc_frame1);
    gaussa.paramOn(mc_frame1, Layout(0.55,0.99,0.87));
    //gauss->fitTo()
    mc_frame1->Draw() ;
    mc_can->Draw();
    mc_can->SaveAs("./graphs/alpha_dist.png", "png");

    

    RooPlot* mc_frame2 =  mcs.plotError(alpha) ; // Draw distribution of parameter error
    mc_frame2->SetTitle("#alpha errors MC estimations");
    mc_frame2->Draw() ;
    mc_can->Draw();
    mc_can->SaveAs("./graphs/alpha_error_dist.png", "png");

    RooPlot* mc_frame3 =  mcs.plotPull(alpha) ; // Draw distribution of parameter pull
    mc_frame3->SetTitle("#alpha pull MC estimations");
    mc_frame3->Draw() ;
    mc_can->Draw();
    mc_can->SaveAs("./graphs/alpha_pull.png", "png");

    //-------
    //--beta
    //-------

    RooRealVar meanb("mean","mean",0.06,0.,1) ;
    RooRealVar sigmab("sigma","sigma",0.001,0.0001,10) ;
    RooGaussian gaussb("gauss","gauss",beta, meanb, sigmab) ;

    gaussb.fitTo(all_beta_est);

    TCanvas* mc_can1 = new TCanvas("mc_can1", "mc_can1", 1000,1000,1000,800);
    RooPlot* mc1_frame1 =  mcs.plotParam(beta) ; // Draw distribution of parameter
    mc1_frame1->SetTitle("#beta MC estimations");
    gaussb.plotOn(mc1_frame1);
    gaussb.paramOn(mc1_frame1, Layout(0.55,0.99,0.87));
    mc1_frame1->Draw() ;
    mc_can1->Draw();
    mc_can1->SaveAs("./graphs/beta_dist.png", "png");

    RooPlot* mc1_frame2 =  mcs.plotError(beta) ; // Draw distribution of parameter error
    mc1_frame2->SetTitle("#beta errors MC estimations");
    mc1_frame2->Draw() ;
    mc_can1->Draw();
    mc_can1->SaveAs("./graphs/beta_error_dist.png", "png");

    RooPlot* mc1_frame3 =  mcs.plotPull(beta) ; // Draw distribution of parameter pull
    mc1_frame3->SetTitle("#beta pulls MC estimations");
    mc1_frame3->Draw() ;
    mc_can1->Draw();
    mc_can1->SaveAs("./graphs/beta_pull.png", "png");

    //-------
    //--gamma
    //-------

    RooRealVar meang("mean","mean",-0.18,-0.3,0.) ;
    RooRealVar sigmag("sigma","sigma",0.001,0.0001,10) ;
    RooGaussian gaussg("gauss","gauss",gamma, meang, sigmag) ;

    gaussg.fitTo(all_gamma_est);

    TCanvas* mc_can2 = new TCanvas("mc_can2", "mc_can2", 1000,1000,1000,800);
    RooPlot* mc2_frame1 =  mcs.plotParam(gamma) ; // Draw distribution of parameter
    mc2_frame1->SetTitle("#gamma MC estimations");
    gaussg.plotOn(mc2_frame1);
    gaussg.paramOn(mc2_frame1, Layout(0.55,0.99,0.87));
    mc2_frame1->Draw() ;
    mc_can2->Draw();
    mc_can2->SaveAs("./graphs/gamma_dist.png", "png");

    RooPlot* mc2_frame2 =  mcs.plotError(gamma) ; // Draw distribution of parameter error
    mc2_frame2->SetTitle("#gamma errors MC estimations");
    mc2_frame2->Draw() ;
    mc_can2->Draw();
    mc_can2->SaveAs("./graphs/gamma_error_dist.png", "png");

    RooPlot* mc2_frame3 =  mcs.plotPull(gamma) ; // Draw distribution of parameter pull
    mc2_frame3->SetTitle("#gamma pulls MC estimations");
    mc2_frame3->Draw() ;
    mc_can2->Draw();
    mc_can2->SaveAs("./graphs/gamma_pull.png", "png");

  */

    //-------------------------------
    //----Binning and chi2-----------
    //-------------------------------

    //RooDataHist data_hist("data_hist","binned version of data",RooArgSet(t,p),*data) ;

    TH2F* dh2_chi = new TH2F("", "",100, 0, M_PI, 100, 0, 2*M_PI);
    data->fillHistogram(dh2_chi, RooArgList(t,p));

/* 
    std::cout << "Imposing errors" << "\n";
    int empty_beans = 0;
    for(int i = 1; i < dh2_chi->GetNbinsX()+1; i++){
      for(int j = 1; j <  dh2_chi->GetNbinsY(); j++){
          int bin = dh2_chi->GetBin(i,j);
          if(dh2_chi->GetBinContent(bin) == 0){
            empty_beans += 1;
            dh2_chi->SetBinError(bin, 0.01);          
            }
      }
    }

    std::cout << empty_beans << "\n";
    RooDataHist data_hist("data_hist","binned version of data",RooArgList(t,p),dh2_chi) ;
    std::cout<< " \n";
    std::cout<< " \n";
    std::cout<< " \n";
    std::cout<< " \n";
    std::cout<< " \n";
    std::cout<< " \n";
    //RooChi2Var chi2("chi2","chi2",*genpdf,data_hist,DataError(RooAbsData::SumW2)) ;
    //RooMinimizer n(chi2) ;
    //n.migrad();
  */
  RooDataHist data_hist("data_hist","binned version of data",RooArgList(t,p),dh2_chi) ;
  

    RooFitResult* r_chi2 = genpdf->chi2FitTo(data_hist, Save());

    auto alpha_est_chi = alpha.getVal();
    auto beta_est_chi = beta.getVal();
    auto gamma_est_chi = gamma.getVal();
    
    //defining a TF2 object with estimated chi2 parameters
    const int npar_chi = 4;
    double dist_params_chi[npar_chi] = {alpha_est_chi, beta_est_chi, gamma_est_chi, 25.};
    TF2* myPdf_chi = new TF2("myPdf_chi", myPDF, 0, M_PI, 0, 2*M_PI, npar_chi);
    myPdf_chi->SetParameters(dist_params_chi);
    myPdf_chi->SetParNames("#alpha", "#beta", "#gamma", "#scaling");

    //plotting chi2 results
    gStyle->SetOptStat(0);

    TCanvas* from_chi2 = new TCanvas("from_chi", "from_chi", 1000,1000,1000,800);
    from_chi2->SetTopMargin(0.15);
    dh2_chi->SetTitle(Form("#splitline{Binned #chi^{2} fit: #alpha = %.3f#pm%.3f, #beta = %.3f#pm%.3f, #gamma = %.3f#pm%.3f}{minimized FCN value: %.5f, EDM: %5f}",alpha_est_chi, alpha.getError(), beta_est_chi, beta.getError(), gamma_est_chi, gamma.getError(),r_chi2->minNll() ,r_chi2->edm() ));
    dh2_chi->SetTitleOffset(5.);
    dh2_chi->GetXaxis()->SetTitle("#theta [rad]");
    dh2_chi->GetXaxis()->SetTitleOffset(2);
    dh2_chi->GetYaxis()->SetTitle("#phi [rad]");
    dh2_chi->GetYaxis()->SetTitleOffset(2);
    dh2_chi->Draw("lego");
    myPdf->Draw("same cont1");
    from_chi2->Draw();
    from_chi2->SaveAs("./graphs/result_CHI2_binned.pdf", "pdf");

    TCanvas* from_chi2_2 = new TCanvas("from_chi2", "from_chi2", 1000,1000,1000,800);
    from_chi2_2->SetTopMargin(0.15);
    dh2_chi->GetXaxis()->SetTitleOffset(0);
    dh2_chi->GetYaxis()->SetTitleOffset(0);
    dh2_chi->Draw("colz");
    myPdf->Draw("same cont1");
    from_chi2_2->Draw();
    from_chi2_2->SaveAs("./graphs/result_CHI2_binned_2.pdf", "pdf");

    //Printing results
    std::cout << "==> Chi2 Fit results" << std::endl ;
    r_chi2->Print() ;
    std::cout << "==> ML Fit results" << std::endl ;
    r_ML->Print() ;

    std::cout << r_ML->edm() << std::endl;
    std::cout << "-log(L) at minimum = " << r_ML->minNll() << std::endl ;
    std::cout << r_chi2->edm() << std::endl;
    std::cout << "-log(L) at minimum = " << r_chi2->minNll() << std::endl ;

    double chi = genpdf->createChi2(data_hist)->getVal();
    std::cout << chi << std::endl;
/* 
    //----------------------------------
    //----Generating uniform pdf--------
    //----------------------------------

    //Defining new variables
    RooRealVar t_uni("t_uni", "t_uni", 0, 2*M_PI);
    RooRealVar p_uni("p_uni", "p_uni", 0, 2*M_PI);
    RooRealVar scale_uni("scale_uni", "scale_uni", -10, 10);
    RooUniform uni_t("uniformt", "uniformt", t_uni);
    RooUniform uni_p("uniformp", "uniformp", p_uni);
    
    //UNIFORM PDF
    RooAddPdf model("model", "uni_t+uni_p", RooArgList(uni_t, uni_p), RooArgList(scale_uni));

    //Plotting data generated from uniform and the pdf itself
    RooDataSet* data2 = model.generate(RooArgSet(t_uni,p_uni),50000) ;
    TH2F* dh3 = new TH2F("h3", "h3", 100, 0, 2*M_PI, 100, 0, 2*M_PI);
    data2->fillHistogram(dh3, RooArgList(t_uni,p_uni));
    dh3->SetTitle("#theta #phi uniform distribution");
    dh3->GetXaxis()->SetTitle("#theta [rad]");
    dh3->GetYaxis()->SetTitle("#phi [rad]");
    
    TH1* hh_pdf = model.createHistogram("t_uni,p_uni", 100,100);
    hh_pdf->SetLineColor(kRed) ;

    TCanvas* c_uni = new TCanvas("c_uni", "c_uni", 1000,1000,1000,800);
    hh_pdf->Draw("surf");
    dh3->Draw("same lego");
    c_uni->Draw();
    c_uni->SaveAs("p.png");

    //----------------------------------
    //----Fitting to previous data------
    //----------------------------------
    RooFitResult* r_unif = model.fitTo(*data2, Save());

    //printing fit result of uniform distribution.
    std::cout << "==> Uniform Fit results" << std::endl ;
    r_unif->Print() ;

    RooFormulaVar llratio_func("llratio","log10(@0)-log10(@1)",RooArgList(*genpdf,model));
    

//Plotting contours
    
    RooPlot* frame_cont3 = n.contour(gamma, beta, 1, 2, 3);
    RooPlot* frame_cont2 = m.contour(alpha, beta, 1, 2, 3);

    frame_cont->SetTitle("Minuit contour plot #alpha #gamma");
    frame_cont2->SetTitle("Minuit contour plot #alpha #beta");
    frame_cont3->SetTitle("Minuit contour plot #gamma #beta");

    TCanvas* ct = new TCanvas("ct", "ct", 1000, 1000, 1000, 800);
    gPad->SetLeftMargin(0.15);
    frame_cont->GetYaxis()->SetTitleOffset(2);
    frame_cont->GetYaxis()->SetTitle("#gamma");
    frame_cont->GetXaxis()->SetTitle("#alpha");
    frame_cont->Draw();
    ct->Draw();
    ct->SaveAs("./cont_ag.png", "png");

    TCanvas* ct2 = new TCanvas("ct2", "ct2", 1000, 1000, 1000, 800);
    gPad->SetLeftMargin(0.15);
    frame_cont2->GetYaxis()->SetTitleOffset(2);
    frame_cont2->GetYaxis()->SetTitle("#beta");
    frame_cont2->GetXaxis()->SetTitle("#alpha");
    frame_cont2->Draw();
    ct2->Draw();
    ct2->SaveAs("./cont_ab.png", "png");

    TCanvas* ct3 = new TCanvas("ct3", "ct3", 1000, 1000, 1000, 800);
    gPad->SetLeftMargin(0.15);
    frame_cont3->GetYaxis()->SetTitleOffset(2);
    frame_cont3->GetYaxis()->SetTitle("#beta");
    frame_cont3->GetXaxis()->SetTitle("#alpha");
    frame_cont3->Draw();
    ct3->Draw();
    ct3->SaveAs("./cont_bg.png", "png");
    
    */
    std::cout<<"Over" << "\n";

    
    return 0;

}
