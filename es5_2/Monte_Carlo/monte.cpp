#include "monte.h"
#include "prng.h"
#include "xorshiro.h"
#include "function.h"
#include "one.h"
#include <iostream>
#include <vector>
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include <cmath>

void monte_carlo::run_simulation(){

    //inizializzo seed e vettore finale.
    srand(time(NULL));
    std::vector<double> sim_res;
    std::vector<double> sim_var;

    xorshiro generator;

    //ciclo sul numero di eventi per generare la simulazione.
    for(int i=0; i<events; i++){
        
        std::vector<double> myarr;
        for(int j = 0; j < measures; j++){
            my_arr.push_back(generator.rand());
        }

        //definisco oggetto funzione passandogli le x random
        func expo(myarr);

        //computo valore delle y, verrà salvato nella classe
        std::vector<double> y = expo.compute(myarr);

        //computo integrale e varianza e appendo alla lista finale
        
        double S_n = 0; //integrale
        double Error = 0; //varianza
        for(int j = 0; j < measures; j++){
            S_n += y[j];
            Error += pow(y[j],2);

        }
        S_n = S_n/measures;
        double variance = Error/measures -pow(S_n, 2);

        sim_res.push_back(S_n);
        sim_var.push_back(variance);
    }
    mc_est = sim_res;
    mc_var = sim_var;


    return;
}



void monte_carlo::run_stratified_simulation(double k){

   //inizializzo seed e vettore finale.
    srand(time(NULL));
    std::vector<double> sim_res;
    std::vector<double> sim_var;

    //Pmf uniforme
    double p = 1./k;
    //Numero di elementi per ogni intervallo
    int elem = int(p*measures);

    xorshiro generator;

    //ciclo sul numero di eventi per generare la simulazione.
    for(int i=0; i<events; i++){

        std::vector<double> var;
        std::vector<double> means;

        //ciclo sugli intervalli
        for(int r = 1; r < k+1; r++){
            
            std::vector<double> myarr;
            for(int h = 0; h < elem; h++){
                my_arr.push_back(generator.rand());
            }

            //numeri casuali nell'intervallo
            std::vector<double>interval_sampling;
            for(int h = 0; h < elem; h++){
                interval_sampling.push_back((r-1)/k + myarr[h]/k);
            }
            func expo(interval_sampling);
            std::vector<double> y = expo.compute(interval_sampling);

            //media nell'intervallo
            double m = expo.mc_integral();
            means.push_back(m);

            //varianza
            double v = 0;
            for(int h = 0; h < elem; h++){
                v += pow(y[h]-m, 2);
            }
            var.push_back(v/elem);
        }

        //calcolo medie e varianze globali su tutti gli intervalli
        double tot_means = 0;
        for(int h = 0; h < means.size(); h++){
            tot_means += means[h];
        }
        tot_means*=p;

        double tot_var = 0;
        for(int h = 0; h < var.size(); h++){
            tot_var += var[h];
        }
        tot_var*=pow(p,2);

        sim_res.push_back(tot_means);
        sim_var.push_back(tot_var);

    }

    mc_est = sim_res;
    mc_var = sim_var;

    return;
};



void monte_carlo::run_importance_expo_simu(){

    srand(time(NULL));
    std::vector<double> sim_res;
    std::vector<double> sim_var;

    xorshiro generator;

    for(int i=0; i<events; i++){
        int s = rand();
        ran evento(measures, rand());

        //routine per generare dati secondo distribuzione esponenziale.
        evento.generate_expo();

        std::vector<double> myarr = evento.get_arr();

        func expo(myarr);

        std::vector<double> y = expo.compute(myarr);

        double integral = expo.mc_integral();
        sim_res.push_back(integral);
        //double var = expo.mc_variance();
        
        sim_res.push_back(integral);
        //sim_var.push_back(var);

    }
    mc_est = sim_res;
    mc_var = sim_var;

    return; 
}


void monte_carlo::run_importance_par_simu(){

    srand(time(NULL));
    std::vector<double> sim_res;
    std::vector<double> sim_var;
    

    for(int i=0; i<events; i++){
        
        ran evento(measures, rand());
        evento.generate_par();

        std::vector<double> para = evento.get_arr();

        func expo(para);
        func parabola(para);

        std::vector<double> y1 = expo.compute(para);
        std::vector<double> y2 = parabola.compute_par(para);

        //computing inline integral.
        double integral = 0;

        for(int j = 0; j < measures; j++){
            integral += y1[j]/y2[j];
        }

        //pushing in vector.
        integral /= measures;
        
        sim_res.push_back(integral);

        //computing Integral and variance:
        double S_n = 0;
        double Err = 0;
        for(int j = 0; j < measures; j++){
            S_n+= y1[j]/y2[j];
            Err+= pow(y1[j]/y2[j], 2);
        }

        S_n /= measures;
        double var = Err/measures - pow(S_n, 2);
        sim_var.push_back(var);

    }
    mc_est = sim_res;
    mc_var = sim_var;

    return;
}

void monte_carlo::run_antithetic_simu(){

    srand(time(NULL));
    std::vector<double> sim_res;
    std::vector<double> sim_var;

    xorshiro generator;

    for(int i = 0; i < events; i++){
        
        int anti_meas = int(measures/2);

        std::vector<double> myarr;
        for(int j = 0; j < anti_meas; j++){
            my_arr.push_back(generator.rand());
        }
        
        std::vector<double> anti = myarr;
        std::for_each(anti.begin(), anti.end(), [](double& d) { d=1.0-d;});

        func expo1(myarr);
        func expo2(anti);
        std::vector<double> y1 = expo1.compute(myarr);
        std::vector<double> y2 = expo2.compute(anti);

        //computo integrale e varianza e appendo alla lista finale
        
        double S_n = 0; //integrale
        double Error = 0; //varianza
        for(int j = 0; j < anti_meas; j++){
            S_n += (y1[j]+y2[j])/2;
            Error += pow((y1[j]+y2[j])/2,2);

        }
        S_n = S_n/anti_meas;
        double variance = Error/anti_meas -pow(S_n, 2);

        sim_res.push_back(S_n);
        sim_var.push_back(variance);

    }

    mc_est = sim_res;
    mc_var = sim_var;

    return;
}


void monte_carlo::plot_estimators(const char* title, const char* output){

    int bins = int(sqrt(events));
    auto min = *std::min_element(mc_est.begin(), mc_est.end());
    auto max = *std::max_element(mc_est.begin(), mc_est.end());

    TH1F* h = new TH1F("Mc integral estimator", "Mc integral estimator", bins, min, max);
    h->SetTitle(title);
    h->GetXaxis()->SetTitle("MC Estimator");
    h->GetYaxis()->SetTitle("Events");
    h->SetLineWidth(2);
    TCanvas*c = new TCanvas("c", "c", 1000,1000,1000,800);

    for(int i=0; i< events; i++){
        h->Fill(mc_est[i]);
    }


    TF1* fit = new TF1("fit", "gaus", min, max);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    h->Fit(fit);
    h->Draw("hist");
    fit->Draw("same");
    c->Draw();
    c->SaveAs(output, "png");
    std::cout << "Results--------------" << std::endl;
    std::cout << "MEan est: " <<fit->GetParameter(1) << std::endl;
    std::cout << "Error est: "<< fit->GetParameter(2)/sqrt(10000) << std::endl;

    return;

}

void monte_carlo::plot_variance(const char* title, const char* output){

    int bins = int(sqrt(events));
    auto min = *std::min_element(mc_var.begin(), mc_var.end());
    auto max = *std::max_element(mc_var.begin(), mc_var.end());

    TH1F* h = new TH1F("Variance distribution", "Variance distribution", bins, min, max);
    h->SetTitle(title);
    h->GetXaxis()->SetTitle("Variance");
    h->GetYaxis()->SetTitle("Events");
    h->SetLineWidth(2);
    TCanvas*c = new TCanvas("c", "c", 1000,1000,1000,800);

    for(int i=0; i< events; i++){
        h->Fill(mc_var[i]);
    }

    TF1* fit = new TF1("fit", "gaus", min, max);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    h->Fit(fit);
    h->Draw("hist");
    fit->Draw("same");
    c->Draw();
    c->SaveAs(output, "png");

    std::cout << "Results--------------" << std::endl;
    std::cout << std::scientific << std::setprecision(15) << "MEan variance: "<< fit->GetParameter(1) << std::endl;
    std::cout << std::scientific << std::setprecision(15) << "Error variance: "<< fit->GetParameter(2)/sqrt(10000) << std::endl;

    return;

}

void monte_carlo::run_importance_par_simu_2(){

    srand(time(NULL));
    std::vector<double> sim_res;
    std::vector<double> sim_var;

    xorshiro generator;
    

    for(int i=0; i<events; i++){

        std::vector<double> para;

        for(int j =0; j < measures; j++){
            double y = generator.rand(0, 1.86274); //Numero casuale in [0,1.8627].
            para.push_back(pow(2.5*y +1, 1./2.5)-1.); //Metodo funzione inversa secondo 2.5*x^(1.5).

        }
        

        func expo(para);
        func parabola(para);

        std::vector<double> y1 = expo.compute(para);
        std::vector<double> y2 = parabola.compute_par_2(para);

        //computing inline integral.
        double integral = 0;

        for(int j = 0; j < measures; j++){
            
            integral += y1[j]/y2[j];
            
        }

        integral /= (measures);

        //computing variance
        double var = 0;
        for(int j = 0 ; j < measures; j++){
            var += pow(y1[j]/y2[j]-integral, 2);
        }
        var /= measures;
        integral*=integral_importance;
        
        sim_res.push_back(integral);
        sim_var.push_back(var);

    }
    mc_var = sim_var;
    mc_est = sim_res;

    return;
}
