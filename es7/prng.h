#ifndef PRNG_H
#define PRNG_H

//definisco la variabile randtype. 
//unsigned (solo 0 o var positive), long long -> 64bit di memoria.
typedef unsigned long long int randtype;

class prng{

    private:
        randtype a;
        randtype c;
        randtype x;

    public:
        prng(); //constructor
        ~prng(); 
        
        randtype randInt(); //numero casuale
        double rand(double min=0., double max =1.); //numero casuale tra 0 e 1
};

#endif
