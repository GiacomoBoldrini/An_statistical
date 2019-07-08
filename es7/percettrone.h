#ifndef PERCETTRONE_H
#define PERCETTRONE_H
#include <vector>

struct line_bound{

    double m;
    double q;
};

class p_clf{
    private:

    int dimension;
    int size;
    std::vector<double> weights;

    public:

    p_clf(int d, int s){ size = s; dimension = d ; } ;
    ~p_clf(){};
    std::vector<double> get_w(){ return weights ; };
    void set_weights();
    int predict(double x, double y);
    void train(double x , double y, int label, double eta);
    line_bound get_boundary();
    
};

#endif // PERCETTRONE_H