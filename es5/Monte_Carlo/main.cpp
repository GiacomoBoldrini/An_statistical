#include "monte.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include <cmath>
#include "TGraphErrors.h"
#include "TLine.h"

int main(){

    
    //Setto numero eventi, numero punti random e vettore estimatori finali.
    int events = 5000;
    int measures = 100000;

    monte_carlo mc(events, measures); //Istanza mc con 100 eventi ognuno da 10000 misure.

    //Eseguo simulazioni dei dati e ne traggo l'integrale con vari metodi
    //di riduzione della varianza:

    //crude_mc
    mc.run_simulation();
    mc.plot_estimators("Crude MC estimators distribution", "./mc_est_crude.png");
    mc.plot_variance("Crude MC variance distribution", "./mc_var_crude.png");


    //stratified
    mc.run_stratified_simulation(5); //Stratified sampling.
    mc.plot_estimators("Stratified MC estimators distribution", "./mc_est_strat.png");
    mc.plot_variance("Stratified MC variance distribution", "./mc_var_strat.png");

    //importance 2
    mc.run_importance_par_simu_2(); //Importance sampling.
    mc.plot_estimators("Importance MC estimators distribution", "./mc_est_imp.png");
    mc.plot_variance("Importance MC variance distribution", "./mc_var_imp.png");

    //antithetic
    mc.run_antithetic_simu(); //Metodo antithetic variates.
    mc.plot_estimators("Antithetic MC estimators distribution", "./mc_est_anti.png");
    mc.plot_variance("Antithetic MC variance distribution", "./mc_var_anti.png");
    
    


    //Analisi della stratificazione

    int k = 1;

    std::vector<double> all_int;
    std::vector<double> all_var;
    std::vector<double> all_errors;
    std::vector<double> all_k;

    for(int m = 0; m < 20; m++){
        std::cout << "Running analysis" << std::endl;
        monte_carlo mc_strat(1, 10000);
        mc_strat.run_stratified_simulation(k);
        std::vector<double> inte = mc_strat.get_est();
        std::vector<double> var = mc_strat.get_var();

        all_int.push_back(inte[0]);
        all_errors.push_back(var[0]/sqrt(10000));
        all_var.push_back(var[0]);
        all_k.push_back(k);

        k++;

    }

    TGraphErrors* g = new TGraphErrors(all_k.size(), &all_k[0], &all_int[0], 0, &all_errors[0]);
    g->SetMarkerStyle(21);
    g->SetMarkerColor(4);
    g->GetXaxis()->SetTitle("Stratifications [k]");
    g->GetYaxis()->SetTitleOffset(1.5);
    g->GetYaxis()->SetTitle("#Theta estimate");
    g->SetTitle("Trend of #Theta stratified estimator varying k");

    TGraph* g_var = new TGraph(all_k.size(), &all_k[0], &all_var[0]);
    g_var->SetMarkerStyle(21);
    g_var->SetMarkerColor(4);
    g_var->SetTitle("Trend of the variance of #Theta stratified estimator");

    TLine* line = new TLine(0, exp(1)-1, 20, exp(1)-1);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);

    TCanvas* c = new TCanvas("c", "c", 1000,1000,1000,800);
    g->Draw("AP");
    line->Draw("same");
    c->Draw();
    c->SaveAs("./Stratified_comp.png");

    g_var->Draw("AP");
    g_var->GetHistogram()->GetXaxis()->SetTitle("Stratifications [k]");
    g_var->GetHistogram()->GetYaxis()->SetTitle("Var[#Theta]");
    g_var->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
    c->SetLogy();
    c->Update();
    c->Draw();
    c->SaveAs("./Stratified_variance.png");

    return 0;
    
}
