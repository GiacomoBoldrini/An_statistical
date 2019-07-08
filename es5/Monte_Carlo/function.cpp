#include "function.h"
#include <iostream>
#include "stdlib.h"
#include <cmath>

std::vector<double> func::compute(std::vector<double> arr){

    int size = arr.size(); //Numero di elementi.
    std::vector<double> final_arr; //Istanzio array finale.
    for(int i = 0; i< size; i++){ //Ciclo su tutti gli elementi.
        final_arr.push_back(exp(arr[i])); //Push nel vettore dell'esponenziale valutata nell'elemento.
    }

    y = final_arr; //Salvo nell'attributo della classe.

    return final_arr;
}

std::vector<double> func::compute_par(std::vector<double> arr){

    int size = arr.size();
    std::vector<double> final_arr;
    for(int i = 0; i < size; i++){
        final_arr.push_back(2.5 * pow(arr[i], 1.5)); //Push nel vettore della funzione valutata nell'elemento.
    }

    y=final_arr; //Salvo nell'attributo della classe;
    return final_arr;
}

std::vector<double> func::compute_par_2(std::vector<double> arr){

    int size = arr.size();
    std::vector<double> final_arr;
    for(int i = 0; i < size; i++){
        //Push nel vettore della funzione valutata nell'elemento. 
        final_arr.push_back(pow(arr[i]+1, 1.5)); 
    }
    y=final_arr; //Salvo nell'attributo della classe;
    return final_arr;
}


double func::mc_integral(){

    double integral = 0; //Istanzio variabile integrale a zero.
    for(int i=0; i<size; i++){ //Ciclo su elementi.
        integral += y[i]; //Aggiungo valutazione funzione nell'elemento a integrale.
    }
    integral /= size; //Media aritmetica.

    return integral;
}
