#ifndef FUNCTION_H
#define FUNCTION_H

#include <vector>

class func{

    int size;
    std::vector<double> x;
    std::vector<double> y;
    public:
        func(std::vector<double> x_gen){ x = x_gen; size = x_gen.size() ; };
        ~func(){};
        std::vector<double> compute(std::vector<double> arr);
        std::vector<double> compute_par(std::vector<double> arr);
        std::vector<double> compute_par_2(std::vector<double> arr);
        double mc_integral();

};

#endif