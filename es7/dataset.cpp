#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "prng.h"
#include "xorshiro.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TLegend.h"
#include "TEllipse.h"


//funzione gaussiana 2d
double gaussian2d(double x, double y, double mux,double muy, double sx, double sy,double rho){


    double t_rho = 1 - rho*rho;
    double norm = 1./(2*M_PI * sx*sy*sqrt(t_rho));
    double t_mux = (x - mux)/sx;
    double t_muy = (y - muy)/sy;
    double exponent = (-1./(2*t_rho))*(t_mux*t_mux + t_muy*t_muy - 2*rho*(t_mux*t_muy));

    return norm*exp(exponent);

    }

//genero numeri casuali con algoritmo xorshiro128p
xorshiro random_number; //Nota istanza da fare fuori dalla funzione!

double rand_uni(double min = 0., double max = 1.){
    
    return random_number.rand(min,max);

}
//struct punti x e y
struct xy{
    double x;
    double y;
};

//genero coppie di numeri che seguono distribuzione gauss2d con try and catch
//Nota ritorno una struct. 
xy gauss_tc(double mux,double muy, double sx, double sy,double rho,
                double xmin, double xmax, double ymin, double ymax){
    
    xy rndxy; 
    //double rndx;
    //double rndy;
    
    int i = 0;
    while(1){

        rndxy.x = rand_uni(xmin,xmax);
        rndxy.y = rand_uni(ymin,ymax);

        double rndz = rand_uni();
        //tc
        double rndGauss = gaussian2d(rndxy.x,rndxy.y,mux,muy,sx,sy,rho);

        if(rndGauss > rndz){ //try and catch

            //restituisce struct punti x,y
            return rndxy; 
        }        
        
    }
    
}
/* 
struct ellipse{

    double x0;
    double y0;
    double theta;
    double r1;
    double r2;
};

ellipse quadratic_discriminator(double mux_s,double muy_s, double sx_s, double sy_s,double rho_s, double mux_b,double muy_b, double sx_b, double sy_b,double rho_b){

    //determinanti matrici covarianza
    //segnale e fondo
    double det_s = pow(sx_s,2)*pow(sy_s,2)-pow(rho_s,2)*pow(sx_s,2)*pow(sy_s,2);
    double det_b = pow(sx_b,2)*pow(sy_b,2)-pow(rho_b,2)*pow(sx_b,2)*pow(sy_b,2);


    
    // x^2 coeff
    double a = pow(sy_b,2)/det_b -pow(sy_s, 2)/det_s; 
    // y^2 coeff
    double c = pow(sx_b,2)/det_b -pow(sx_s, 2)/det_s; 
    // 2*xy coeff
    double b = 2*(-(rho_b*sx_b*sy_b)/det_b + (rho_s*sx_s*sy_s)/det_s); 
    // 2*x coeff
    double d = ((2/det_s)*(mux_s*sy_s*sy_s-rho_s*sx_s*sy_s*muy_s)) - ((2/det_b)*(mux_b*sy_b*sy_b-rho_b*sx_b*sy_b*muy_b));
    // 2*y coeff
    double f = ((2/det_s)*(muy_s*sy_s*sy_s-rho_s*sx_s*sy_s*mux_s)) - ((2/det_b)*(muy_b*sy_b*sy_b-rho_b*sx_b*sy_b*mux_b)); 
    // costante
    double g = (1/det_b)*mux_b*mux_b*sy_b*sy_b - (1/det_s)*mux_s*mux_s*sy_s*sy_s - (2/det_b)*(rho_b*sx_b*sy_b)*muy_b*mux_b + (2/det_s)*(rho_s*sx_s*sy_s)*muy_s*mux_s +(1/det_b)*sx_b*sx_b*muy_b*muy_b - (1/det_s)*sx_s*sx_s*muy_s*muy_s +log(det_b/det_s);

    // troviamo i parametri veri
    b /= 2;
    d /= 2;
    f /= 2;

    double x0 = (c*d-b*f)/(b*b-a*c);
    double y0 = (a*f-b*d)/(b*b-a*c);

    double r_one = sqrt((2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g))/((b*b-a*c)*(sqrt(pow(a-c,2)+4*b*b)-(a+c))));
    double r_two = sqrt((2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g))/((b*b-a*c)*(-sqrt(pow(a-c,2)+4*b*b)-(a+c))));

    double theta;

    if (a == c){
        theta = 45.;
    }

    ellipse result;
    result.x0 = x0;
    result.y0 = y0;
    result.theta = theta;
    result.r1 = r_one;
    result.r2 = r_two;

    return result;
}

*/

struct ellipse{

    double a;
    double b;
    double c;
    double d;
    double f;
    double g;
};

ellipse quadratic_discriminator(double mux_s,double muy_s, double sx_s, double sy_s,double rho_s, double mux_b,double muy_b, double sx_b, double sy_b,double rho_b){

    //determinanti matrici covarianza
    //segnale e fondo
    double det_s = pow(sx_s,2)*pow(sy_s,2)-pow(rho_s,2)*pow(sx_s,2)*pow(sy_s,2);
    double det_b = pow(sx_b,2)*pow(sy_b,2)-pow(rho_b,2)*pow(sx_b,2)*pow(sy_b,2);


    
    // x^2 coeff
    double a = pow(sy_b,2)/det_b -pow(sy_s, 2)/det_s; 
    // y^2 coeff
    double c = pow(sx_b,2)/det_b -pow(sx_s, 2)/det_s; 
    // 2*xy coeff
    double b = 2*(-(rho_b*sx_b*sy_b)/det_b + (rho_s*sx_s*sy_s)/det_s); 
    // 2*x coeff
    double d = ((2/det_s)*(mux_s*sy_s*sy_s-rho_s*sx_s*sy_s*muy_s)) - ((2/det_b)*(mux_b*sy_b*sy_b-rho_b*sx_b*sy_b*muy_b));
    // 2*y coeff
    double f = ((2/det_s)*(muy_s*sy_s*sy_s-rho_s*sx_s*sy_s*mux_s)) - ((2/det_b)*(muy_b*sy_b*sy_b-rho_b*sx_b*sy_b*mux_b)); 
    // costante
    double g = (1/det_b)*mux_b*mux_b*sy_b*sy_b - (1/det_s)*mux_s*mux_s*sy_s*sy_s - (2/det_b)*(rho_b*sx_b*sy_b)*muy_b*mux_b + (2/det_s)*(rho_s*sx_s*sy_s)*muy_s*mux_s +(1/det_b)*sx_b*sx_b*muy_b*muy_b - (1/det_s)*sx_s*sx_s*muy_s*muy_s +log(det_b/det_s);

    // troviamo i parametri veri
    b /= 2;
    d /= 2;
    f /= 2;

    double x0 = (c*d-b*f)/(b*b-a*c);
    double y0 = (a*f-b*d)/(b*b-a*c);

    double r_one = sqrt((2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g))/((b*b-a*c)*(sqrt(pow(a-c,2)+4*b*b)-(a+c))));
    double r_two = sqrt((2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g))/((b*b-a*c)*(-sqrt(pow(a-c,2)+4*b*b)-(a+c))));

    double theta;

    if (a == c){
        theta = 45.;
    }

    ellipse result;
    result.a = a;
    result.b = b;
    result.c = c;
    result.d = d;
    result.f = f;
    result.g = g;

    return result;
}

int disc_evaluator(ellipse params, double x, double y){

    if(params.a*x*x + 2*params.b*x*y +params.c*y*y+ 2*params.d*x+2*params.f*y +params.g >= 0){ 
        return 0;
    }
    else{
        return 1;
    }
}





int main(int argc, char** argv){

    //Signal gaussian params
    double mux_sig = 0.;
    double muy_sig = 0.;
    double sx_sig = 0.3;
    double sy_sig = 0.3;
    double rho_sig = 0.5;

    double xmin = -6;
    double xmax = 10;
    double ymin = -4;
    double ymax = 10;

    int N_sig = 10000;

    //Background gaussian params
    double mux_bkg = 4.;
    double muy_bkg = 4.;
    double sx_bkg = 1;
    double sy_bkg = 1;
    double rho_bkg = 0.4;

    int N_bkg = 10000;  
    
    TCanvas* c1 = new TCanvas("c1","c1",1000,1000, 1000, 800);
    TH2F* gauss_sig = new TH2F("gauss_sig","gauss_sig", 100, xmin, xmax,100,ymin,ymax);
    gauss_sig->SetFillColor(kRed);
    TH2F* gauss_bkg = new TH2F("gauss_bkg","gauss_bkg",100,xmin,xmax,100,ymin,ymax);
    gauss_bkg->SetFillColor(kBlue);
    TLegend* legend = new TLegend(.1, .8, .4, .9);
    legend->AddEntry(gauss_sig, "Signal", "f");
    legend->AddEntry(gauss_bkg, "Background", "f");
    
    xy xy_sig; //Struct per il segnale
    xy xy_bkg; //Struct per il fondo


    double xentry_sig;
    double yentry_sig;
    double xentry_bkg;
    double yentry_bkg;

    //TMVA Legge dati dai tree, ne creo uno con dentro le variabili x e y per segnale e fondo.
    TTree* tree_sig = new TTree("Signal","Signal");
    tree_sig->Branch("x", &xentry_sig);
    tree_sig->Branch("y", &yentry_sig);

    TTree* tree_bkg = new TTree("Background","Background");
    tree_bkg->Branch("x",&xentry_bkg);
    tree_bkg->Branch("y",&yentry_bkg);

    TFile* outfile = new TFile("dataset.root", "RECREATE");

    ellipse disc = quadratic_discriminator(mux_sig, muy_sig, sx_sig, sy_sig, rho_sig, mux_bkg, muy_bkg, sx_bkg, sy_bkg, rho_bkg);
    
    int signal_true = 0;
    int signal_neg = 0;
    int bkg_true = 0;
    int bkg_neg = 0;

    for(int i = 0; i < N_sig; i++){

        xy_sig = gauss_tc(mux_sig,muy_sig,sx_sig,sy_sig,rho_sig,xmin,xmax,ymin,ymax);

        xentry_sig = xy_sig.x;
        yentry_sig = xy_sig.y;

        if(disc_evaluator(disc, xy_sig.x, xy_sig.y) == 0){
            signal_true += 1;
        }
        else{
            signal_neg += 1;
        }

        gauss_sig->Fill(xy_sig.x,xy_sig.y);
        tree_sig->Fill();
    }

    std::cout<< "True positive: " << signal_true << " " << "Accuracy sig: " << signal_true/N_sig << std::endl;
    std::cout<< "False negative: " << signal_neg << " " << "Test Error Rate sig: " << signal_neg/N_sig << std::endl;
    

    for(int i = 0; i < N_bkg; i++){

        xy_bkg = gauss_tc(mux_bkg,muy_bkg,sx_bkg,sy_bkg,rho_bkg,xmin,xmax,ymin,ymax);

        xentry_bkg = xy_bkg.x;
        yentry_bkg = xy_bkg.y;

        if(disc_evaluator(disc, xy_bkg.x, xy_bkg.y) == 1){
            bkg_true += 1;
        }
        else{
            bkg_neg += 1;
        }

        gauss_bkg->Fill(xy_bkg.x,xy_bkg.y);
        tree_bkg->Fill();
    }
    std::cout<< "True negative: " << signal_true << " " << "Accuracy bkg: " << signal_true/N_sig << std::endl;
    std::cout<< "False positive: " << signal_neg << " " << "Test Error Rate bkg : " << signal_neg/N_sig << std::endl;

    std::cout << " Total accuracy in classifying sig+bkg: " << (signal_true+bkg_true)/(N_sig+N_bkg) << std::endl;
    std::cout << " Total error rate in classifying sig+bkg: " << 1 - (signal_true+bkg_true)/(N_sig+N_bkg) << std::endl;



    //TEllipse* el = new TEllipse(-0.426991, -0.426991, 2.12586, 1.21307, 0, 360, 45);
    //TEllipse* el = new TEllipse(disc.x0, disc.y0, disc.r2, disc.r1, 0, 360, disc.theta);
    //el->SetFillStyle(0);
    //legend->AddEntry(el, "Quadratic classifier", "l");
    
    gauss_sig->SetMarkerColor(kRed);
    gauss_bkg->SetMarkerColor(kBlue);
    gauss_sig->GetXaxis()->SetTitle("X [a.u.]");
    gauss_sig->GetYaxis()->SetTitle("Y [a.u.]");
    gauss_sig->SetTitle("Distributions of Signal: N(#mu = #bar{0}, #sigma = #bar{0.3}, #rho = 0.5), Background: N(#mu = #bar{4}, #sigma = #bar{1}, #rho = 4 )");
    gauss_sig->SetStats(0000);
    gauss_sig->Draw();
    gauss_bkg->Draw("SAME");
    //el->Draw("same");
    legend->Draw();
    c1->Print("./analysis/dataset.pdf","pdf");
    
    tree_sig->Write();
    tree_bkg->Write();
    outfile->Close();

    

 return 0;   
}
