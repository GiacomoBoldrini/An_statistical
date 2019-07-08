#include "TRandom3.h"
#include "svm.h"
#include <vector>
#include <iostream>
#include <numeric>

void svm_clf::set_weights(){

    //inizializza random il vettore dei pesi
    TRandom3* gen = new TRandom3();
    std::vector<double> tr_w {27.7986, -5.714, -5.714, -13.62, +13.85, -13.62};
    for(int i = 0; i < f_space_dim; i++){
        //meglio settare pesi randomici 
        //grandi e iterare più volte su tutto
        //il dataset.
        weight.push_back(gen->Uniform(-20,20));
        //weight.push_back(tr_w[i]);
    }

    return;

}

std::vector<double> svm_clf::mapping(double x, double y){

    std::vector<double> pol_map;
    pol_map.push_back(1);
    pol_map.push_back(sqrt(2)*x);
    pol_map.push_back(sqrt(2)*y);
    pol_map.push_back(x*x);
    pol_map.push_back(sqrt(2)*x*y);
    pol_map.push_back(y*y);

    return pol_map;
}

double svm_clf::predict_train(double x, double y){

    double r1 = 0;

    std::vector<double> phi_x = mapping(x, y);
    for(int i = 0; i < f_space_dim; i++){
        r1+= weight[i]*phi_x[i];
    }
    return r1;

}

int svm_clf::predict(double x, double y){

    std::vector<double> phi_x = mapping(x, y);
    double r1 = std::inner_product(weight.begin(), weight.end(), phi_x.begin(), 0);
    if(r1 >= 0) return 1;
    else return -1;

}

void svm_clf::train(double x, double y, int label, double eta){

    double p = predict_train(x,y);
    
    std::vector<double> phi = mapping(x,y);
    if(label*p < 1 ){
        for(int i = 0; i < f_space_dim; i++){
            //weight[i] = weight[i] - eta*((2*weight[i])/(reg*size) - label*phi[i]);
            weight[i] = weight[i] - eta*((2*weight[i])/(reg*size*100) - label*phi[i]);
            
        }
    }
    else{
        for(int i = 0; i < f_space_dim; i++){
            weight[i] = weight[i]-eta*(2/(reg*size))*weight[i];
        }
    }

}

TRandom3* gen = new TRandom3();

void svm_clf::momentum(){

    for(int i = 0; i < f_space_dim; i++){
        weight[i] += gen->Uniform(-20, 20);
    }
    return;
}

