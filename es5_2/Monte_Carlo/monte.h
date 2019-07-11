#ifndef MONTE_H
#define MONTE_H
#include <vector>

class monte_carlo{
    int events;
    int measures;
    int seed;
    std::vector<double> mc_est;
    std::vector<double> mc_var;
    public:
        monte_carlo(int n, int m){ events = n ; measures = m; seed = 0; };
        ~monte_carlo(){};
        void set_seed(int s){seed = s ;}
        void run_simulation();

        void run_stratified_simulation(double p);
        void run_importance_expo_simu();
        void run_importance_par_simu();
        void run_importance_par_simu_2();
        void run_antithetic_simu();
        std::vector<double> get_var(){return mc_var; } ;
        std::vector<double> get_est(){return mc_est; } ;
        void plot_estimators(const char* title, const char* output);
        void plot_variance(const char* title, const char* output);
        

};

#endif


