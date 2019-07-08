#ifndef SVM_H
#define SVM_H
#include <vector>


class svm_clf{
    private:
        double reg;
        int size;
        int dimension;
        int f_space_dim = 6;
        std::vector<double> weight;
        
    public:
        svm_clf(int s, int d, int r){ size = s; dimension = d; reg = r; } ;
        ~svm_clf(){};
        int predict(double x, double y);
        double predict_train(double x, double y);
        double kernel(double x, double y);
        void set_weights();
        void train(double x, double y, int label, double eta);
        std::vector<double> get_weight(){ return weight ;} ;
        std::vector<double> mapping(double x, double y);
        void momentum();

};

#endif //SVM_H