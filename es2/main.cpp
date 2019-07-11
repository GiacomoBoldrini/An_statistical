#include <iostream>
#include "sum_alg.h"
#include "gamma_const.h"
#include <cmath>
#include <numeric>
#include <iomanip>
#include <vector>
#include <gsl/gsl_math.h> //M_EULER
#include <limits>
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLegend.h"


int main(int argc, char** argv){

    gamma_const gammaEM;
    long long int N = 12;
    std::cout << "Euler constant: " << std::setprecision(std::numeric_limits<long double>::digits10 +1) << M_EULER << std::endl;

    long int size = 12;
    double* emF = new double[size];
    double* emD = new double[size];
    long double* emLD = new long double[size];
    double* 
    double* stdemF = new double[size];
    double* stdemD = new double[size];
    double* stdemLD = new double[size];
    double* x_axis = new double[size];

    for(int i = 1; i < size+1; i++){
        std::cout << i << std::endl;
        long double x = pow(10,(long double)i);
        std::cout << gammaEM.gamma_constLDLLI(x) << std::endl;
        x_axis[i-1]=x;

        emF[i-1]=gammaEM.gamma_constFI(x);
        emD[i-1]=gammaEM.gamma_constDI(x);
        emLD[i-1]=gammaEM.gamma_constLDI(x);
        stdemF[i-1]=abs(gammaEM.gamma_constFI(x) - M_EULER);
        stdemD[i-1]=abs(gammaEM.gamma_constDI(x) - M_EULER);
        stdemLD[i-1]=abs(gammaEM.gamma_constLDI(x) - M_EULER);

    }


    TCanvas* c1 = new TCanvas("Euler-Mascheroni","c1",800,800);
    TMultiGraph* mg = new TMultiGraph();
    TMultiGraph* stdmg = new TMultiGraph();
    TLegend * leg = new TLegend();
    TLegend * legstd = new TLegend(0.1,0.7,0.3,0.9);
    

    TGraph * gF = new TGraph(size,x_axis,emF);
    TGraph * stdgF = new TGraph(size,x_axis,stdemF);
    TGraph * gD = new TGraph(size,x_axis,emD);
    TGraph * stdgD = new TGraph(size,x_axis,stdemD);
    TGraph * gLD = new TGraph(size,x_axis,emLD);
    TGraph * stdgLD = new TGraph(size,x_axis,stdemLD);
    
    mg->Add(gF);
    mg->Add(gLD);
    mg->Add(gD);

    stdmg->Add(stdgF);
    stdmg->Add(stdgD);    
    stdmg->Add(stdgLD);

    leg->AddEntry(gF, "Float");
    leg->AddEntry(gD, "Double");
    leg->AddEntry(gLD, "Long double");
    

    mg->SetTitle("#gamma vs n_{iter}");
    mg->Draw("AC* PMC PLC");
    leg->Draw("same");
    c1->SetLogx();

    c1->SaveAs("euler_masc.pdf","pdf");
    c1->Close();
    
    TCanvas* c2 = new TCanvas("Euler-Mascheroni_std","c2",800,800);  
        
    legstd->AddEntry(stdgF,"Float");
    legstd->AddEntry(stdgD,"Double");
    legstd->AddEntry(stdgLD,"Long double");
    stdmg->SetMaximum(5);
    stdmg->SetMinimum(-1.);
    stdmg->Draw("AC* PMC PLC");
    legstd->Draw("Same");
    c2->SetLogx();

    c2->SaveAs("euler_masc_std.pdf","pdf");


    









    return 0;

}
