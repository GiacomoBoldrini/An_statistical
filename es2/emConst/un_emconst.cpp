#include <iostream>
#include "sum_alg.h"
#include "gamma_const.h"
#include <cmath>
#include <numeric>
#include <iomanip> //setPrecision
#include <vector>
#include <gsl/gsl_math.h> //M_EULER
#include <limits>
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TApplication.h"


int main(int argc, char** argv){

    //Object TApplication istance
    TApplication* Grafica = new TApplication("",NULL,NULL);

    //Object gamma_const istance (class implemented)
    gamma_const gammaEM;
    

    std::cout << "Euler constant: " << std::setprecision(std::numeric_limits<long double>::digits10 +1) << M_EULER << std::endl;

    int size = 9; //Iterations order of magnitu
    double* x_axis = new double[size];
    double* emF = new double[size];
    double* emLD = new double[size];
    double* stdemF = new double[size];    
    double* stdemLD = new double[size];


    double* KemF = new double[size];
    double* KemLD = new double[size];
    double* KstdemF = new double[size];    
    double* KstdemLD = new double[size];

    double* PWemF = new double[size];
    double* PWemLD = new double[size];
    double* PWstdemF = new double[size];    
    double* PWstdemLD = new double[size];


    for(int i = 1; i < size+1; i++){
        
        double x = pow(10,(double)i);
        std::cout << "Iterazione: " << "10E+" << i <<"  Constant computed: " <<
        std::setprecision(std::numeric_limits<long double>::digits10 +1) << gammaEM.PWgamma_constLDI(x) << std::endl;
        x_axis[i-1]= x;

        KemF[i-1]=gammaEM.Kgamma_constFI(x);
        emF[i-1]=gammaEM.gamma_constFI(x);
        KemLD[i-1]=gammaEM.Kgamma_constLDI(x);
        emLD[i-1]=gammaEM.gamma_constLDI(x);
        KstdemF[i-1]=abs(gammaEM.Kgamma_constFI(x) - M_EULER);
        stdemF[i-1]=abs(gammaEM.gamma_constFI(x) - M_EULER)/M_EULER;
        KstdemLD[i-1]=abs(gammaEM.Kgamma_constLDI(x) - M_EULER);
        stdemLD[i-1]=abs(gammaEM.gamma_constLDI(x) - M_EULER)/M_EULER;
        PWemF[i-1]=gammaEM.PWgamma_constFI(x);
        PWstdemF[i-1]=(gammaEM.PWgamma_constFI(x)-M_EULER);        
        PWemLD[i-1]=gammaEM.PWgamma_constLDI(x);
        PWstdemLD[i-1]=(gammaEM.PWgamma_constLDI(x)-M_EULER);
    }


    TMultiGraph* stdmg = new TMultiGraph();
    TLegend * legstd = new TLegend(0.1,0.1,0.3,0.3);  
    
    TGraph * stdgF = new TGraph(size,x_axis,stdemF);    
    TGraph * stdgLD = new TGraph(size,x_axis,stdemLD);
    
    TGraph * KstdgF = new TGraph(size,x_axis,KstdemF);    
    TGraph * KstdgLD = new TGraph(size,x_axis,KstdemLD);
    
    TGraph * PWstdgF = new TGraph(size,x_axis,PWstdemF);    
    TGraph * PWstdgLD = new TGraph(size,x_axis,PWstdemLD);    


    stdgF->SetLineWidth(2);
    stdgLD->SetLineWidth(2);
    stdgF->SetLineColor(kRed);
    stdgLD->SetLineColor(kBlue);

    KstdgF->SetLineWidth(2);
    KstdgLD->SetLineWidth(2);
    KstdgF->SetLineColor(kOrange);
    KstdgLD->SetLineColor(kGreen);
    
    PWstdgF->SetLineWidth(2);
    PWstdgLD->SetLineWidth(2);
    PWstdgF->SetLineColor(kTeal);
    PWstdgLD->SetLineColor(kMagenta-4);
    
    
    stdmg->Add(stdgF);   
    stdmg->Add(stdgLD);

    stdmg->Add(KstdgF);   
    stdmg->Add(KstdgLD);
   
    stdmg->Add(PWstdgF);
    stdmg->Add(PWstdgLD);

    TCanvas* c1 = new TCanvas("Euler-Mascheroni_std","c1",800,800);  
        
    legstd->AddEntry(stdgF,"Direct-Float");
    legstd->AddEntry(stdgLD,"Direct-Long double");

    legstd->AddEntry(KstdgF,"Kahan-Float");
    legstd->AddEntry(KstdgLD,"Kahan-Long double");

    legstd->AddEntry(PWstdgF,"Pairwise-Float");
    legstd->AddEntry(PWstdgLD,"Pairwise-Long double");
    stdmg->SetTitle("Relative error vs n_{iter}");
    stdmg->GetXaxis()->SetTitle("n_{iter}");
    stdmg->GetYaxis()->SetTitle("#gamma - #gamma_{true}/#gamma");
    stdmg->GetYaxis()->SetTitleOffset(1.5);
    stdmg->Draw("APC*");
    legstd->Draw("SAME");
    c1->SetLogx();
    c1->SetLogy();
    c1->SaveAs("./results/all_euler_masc_std.pdf","pdf");


    Grafica->Run();
    return 0;

}
