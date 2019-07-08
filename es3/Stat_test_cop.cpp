#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/HybridCalculatorOriginal.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooRealVar.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TLatex.h"
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
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
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
using namespace RooStats;

double myPDF(double *x, double *par){ //to define TF2 obj
  return par[3]*(3./(4.*M_PI))*(0.5*(1.-par[0]) + (0.5)*(3.*par[0]-1)*cos(x[0])*cos(x[0]) - par[1]*sin(x[0])*sin(x[0])*cos(2.*x[1])- sqrt(2.)*par[2]*sin(2.*x[0])*cos(x[1]));
}

int main(){
    
    //defining Real variables according to RooFit
    //Angles physical defined in [0-2*pi]
    RooRealVar  t("t","t",0,2*M_PI) ;
    RooRealVar  p("p","p", 0, M_PI) ;
    
    //Defining Parameters according to RooFit.
    //Initial value is the 3rd argument, 4th and 5th arg
    //are the range of existance of the parameter.
    RooRealVar alpha("alpha", "alpha", 0.65, 0., 0.7);

    RooRealVar beta("beta", "beta", 0.06, 0., 0.1);
    RooRealVar gamma("gamma", "gamma", -0.18, -0.5, 0.);
    RooRealVar scale("scale", "scale", 1., 0, 10);
    RooRealVar scale_uni("scale_uni", "scale_uni", 1.);
    RooRealVar for_comp("for_comp", "for_comp", 0.1, 0, 1);

    RooUniform uni_t("uniformt", "uniformt", t);
    RooUniform uni_p("uniformp", "uniformp", p);
    
    //UNIFORM PDF
    RooAddPdf unif("unif", "uni_t+uni_p", RooArgList(uni_t, uni_p), RooArgList(scale_uni));
    
    //ClassFactory generates a pdf instance and normalize it.
    RooRealVar n("n","number of events",10000,0,20000) ;
    RooAbsPdf* genpdf = RooClassFactory::makePdfInstance("GenPdf","(3./(4.*M_PI))*(0.5*(1.-alpha) + (0.5)*(3.*alpha-1)*cos(t)*cos(t) - beta*sin(t)*sin(t)*cos(2.*p)- sqrt(2.)*gamma*sin(2.*t)*cos(p))",RooArgSet(t,p,alpha, beta, gamma)) ;
    //RooExtendPdf egenpdf("egenpdf","extended gen PDF",*genpdf,n);
    RooDataSet* data = genpdf->generate(RooArgSet(t,p), 50000);


    std::cout<< "--------------------------"<<std::endl;
    std::cout<< "Statistical_Test_3"<<std::endl;
    std::cout<< "--------------------------"<<std::endl;

    RooWorkspace*w = new RooWorkspace("w");
    w->import(t);
    w->import(p);
    w->import(*genpdf);
    w->import(*data);

    w->defineSet("obs",RooArgSet(t,p));
    w->defineSet("arg",RooArgSet(alpha, beta, gamma));
    w->defineSet("nuisParams","");
    ModelConfig zero_model("zero_model", w);
    zero_model.SetObservables(*w->set("obs"));
    zero_model.SetPdf(*w->pdf("GenPdf"));
    zero_model.SetParametersOfInterest(*w->set("arg"));
    w->var("alpha")->setVal(1./3.);  // important only uniform
    w->var("beta")->setVal(0.0);
    w->var("gamma")->setVal(0.0);
    zero_model.SetSnapshot(*w->set("arg"));

    // create the alternate meaning uno boson.
    ModelConfig uno_model("uno_model", w);
    uno_model.SetPdf(*w->pdf("GenPdf"));
    uno_model.SetObservables(*w->set("obs"));
    uno_model.SetParametersOfInterest(*w->set("arg"));
    w->var("alpha")->setVal(0.65);  // important for distribution genpdf
    w->var("beta")->setVal(0.06);
    w->var("gamma")->setVal(-0.18);
    uno_model.SetSnapshot(*w->set("arg"));

    w->Print();

//-----------------------------------
//----Computing Likelihood ration----
//-----------------------------------

double prod_H1; //first sum variable
double prod_H0 ; //second sum variable

long double ln_r = 0.; //final statistic
double unif_prob = 1./(M_PI*2*M_PI); //a uniform distribution
RooArgSet* pdfObs = genpdf->getObservables(*data) ; //for retrieving data

//lambda = (prod(f(x|H1))/prod(f(x|H0)))
//log(lambda) = sum(log(f(x|H1))) - sum(log(f(x|H0)))

    //cycle on all data
    for(int i = 0; i < 50000; i++){

      //retireving the row
      auto t_p_ref = *data->get(i);
      double t_val = t_p_ref.getRealValue("t");

      //retrieving the value of the variables
      double p_val = t_p_ref.getRealValue("p");

      //the probability is uniform for each entry, does not depend on 
      //theta and phi.
      prod_H0 = log(unif_prob);

      //Retrieving height of the pdf at the dataset point
      *pdfObs = *data->get(i);
      prod_H1 = log(genpdf->getVal());

      //subracting logarithm and summing them to the statistics
      ln_r += prod_H0-prod_H1;
      std::cout << prod_H0 << std::endl;
      std::cout << prod_H1 << std::endl;
      std::cout << ln_r << std::endl;

    }

    std::cout << "------------------" << std::endl;
    std::cout << "Likelihood Results:" << std::endl;
    std::cout << "------------------" << std::endl;

    
    std::cout << "Ratio: " << ln_r << std::endl;




    double confLevel = 0.95;
    int nScanPoints = 50;
    bool plotAsTF1 = false;
    double poiMinPlot = 1;
    double poiMaxPlot = 0;
    bool doHypoTest = false;
    double nullValue = 0;
    
    ProfileLikelihoodCalculator pl(*data, uno_model);
    pl.SetConfidenceLevel(confLevel); // 95% interval
    LikelihoodInterval *interval = pl.GetInterval();
    RooRealVar *firstPOI = (RooRealVar *)uno_model.GetParametersOfInterest()->first();
    std::cout << "\n>>>> RESULT : " << confLevel * 100 << "% interval on " << firstPOI->GetName() << " is : ["
        << interval->LowerLimit(*firstPOI) << ", " << interval->UpperLimit(*firstPOI) << "]\n " << std::endl;
    
    

return 0;

}