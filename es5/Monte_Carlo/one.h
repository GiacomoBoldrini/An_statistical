#ifndef ONE_H
#define ONE_H

#include <vector>

class ran{
    int size;
    std::vector<double> arr;
    int seed;
    int interval_elements;
    std::vector<double> strat_arr;

    public:
    ran(int s, int se){ size = s ; seed = se ; };
    ~ran(){};
    //methods
    
    int get_size(){ return size ; };
    int get_strat_size(){ return interval_elements; };
    std::vector<double> get_arr(){ return arr ; };
    std::vector<double> get_stratified_arr(){ return strat_arr ; };
    void generate();
    void generate_interval(std::vector<double> intervals);
    void generate_expo();
    void generate_par();
    void generate_par_2();

};

#endif