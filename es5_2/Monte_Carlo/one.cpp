#include "one.h"
#include <iostream>
#include "xorshiro.h"
#include "prng.h"
#include <ctime>
#include <cmath>

void ran::generate(){

    xorshiro gen; //Creo generatore TRandom3. Seed da attributo classe.
    for(int i = 0; i < size; i++){ //Ciclo su size, attributo della classe.
        gen.rand();
        //arr.push_back(gen->rand(0, 1)); //Puso nell'array attributo della classe il numero casuale uniforme.
    }
    
    return; //Void.
}

void ran::generate_interval(std::vector<double> intervals){

    int el_interval = intervals.size() - 1; //Ottengo il numero di intervalli.
    int size_first_int = int(size*intervals[1]);
    int size_second_int = int(size*(intervals[2]-intervals[1])); //Numeri random equidistribuiti.

    std::vector<int> sizes;
    sizes.push_back(size_first_int);
    sizes.push_back(size_second_int);

    xorshiro gen;
    std::vector<double> stratified; //Vettore numeri random stratificati.

    for(int i = 0; i < el_interval; i++){ //Ciclo sul numero di intervalli.
        for(int j = 0; j < sizes[i]; j++){ //Ciclo sul numero di elementi per intervallo da generare.
            stratified.push_back(gen.rand(intervals[i], intervals[i+1])); //Push nell'array del numero generato uniformemente nell'intervallo specificato.
        }
    }

    //Salvo nella classe il vettore e la dimensione equispaziata 
    //degli intervalli. Servirà per valutare l'integrali a pezzi.

    strat_arr = stratified;
    

    return; //Void.
}

void ran::generate_expo(){
   
    xorshiro gen;
    
    for(int i=0; i < size; i++){
        double y = gen.rand(0,exp(1)-1); //Genero in modo da avere valori in [0,e-1].
        arr.push_back(log(y+1)); //Metodo funzione inversa da valori in [0,1].
    }

    return; //Void.
}

void ran::generate_par(){

    xorshiro gen;;

    for(int i =0; i < size; i++){
        double y = gen.rand(0,1); //Numero casuale in [0,1].
        arr.push_back(pow(y, 1./2.5)); //Metodo funzione inversa secondo 2.5*x^(1.5).
    }

    return; //Void.
}  

void ran::generate_par_2(){

    xorshiro gen;;

    for(int i =0; i < size; i++){
        double y = gen.rand(0, 1.86274); //Numero casuale in [0,1.8627].
        arr.push_back(pow(2.5*y +1, 1./2.5)-1.); //Metodo funzione inversa secondo 2.5*x^(1.5).

    }

    return; //Void.
}  