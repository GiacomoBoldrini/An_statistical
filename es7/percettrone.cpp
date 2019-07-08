#include <iostream>
#include "TRandom3.h"
#include "percettrone.h"

void p_clf::set_weights(){

    //inizializza random il vettore dei pesi
    TRandom3* gen = new TRandom3();
    for(int i = 0; i < dimension+1; i++){
        weights.push_back(gen->Uniform(0,1));
    }

    return;

}

int p_clf::predict(double x, double y){

    double activation = weights[0];
    double sum = x*weights[1] + y*weights[2] + activation;

    if(sum >= 0){
        return 1;
    }
    else{
        return 0;
    }

}

void p_clf::train(double x , double y, int label, double eta){

    int pred = predict(x,y);
    double error = label-pred;
    weights[0] = weights[0] + eta*error;
    weights[1] = weights[1] + eta*error*x;
    weights[2] = weights[2] + eta*error*y;

    return;
}

line_bound p_clf::get_boundary(){

    line_bound results;
    results.m = -weights[1]/weights[2];
    results.q = -weights[0]/weights[2];

    return results;
}





