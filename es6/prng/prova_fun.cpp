#include <iostream>
#include "simulation.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <vector>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>

int main(int argc, char** argv){

    //simulation exp;
    int N = 3;
    //int N = 50000;
    double diameter = 0.75e-5;
    double lambda = 500e-9;
    double c = 1;


    double i = 0;
    int k = 0;
    double* x_axis = new double[30];
    double* y_axis = new double[30];

    simulation exp(N,c,lambda,diameter);

    while(i<N){

        y_axis[k]=exp.diffraction_fun(i);
        x_axis[k]=i;
        i +=0.1;
        k+=1;              
    }

    TCanvas* c1 = new TCanvas("c","c",700,700);
    TGraph * g = new TGraph(30,x_axis,y_axis);
    g->Draw("APC");
    c1->Draw();
    c1->Print("./try.pdf","pdf");




    return 0;
}