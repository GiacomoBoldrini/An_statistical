#include "monte_carlo.h"
#include "TRandom3.h"
#include <cmath>

generation generate_mc::simulation(){

    generation result;

    TRandom3* gen_p = new TRandom3(seed); //generatore causale momento
    TRandom3* gen_cos = new TRandom3(seed + rand()); //generatore causale angoli

    double p = gen_p->Uniform(0,10);
    double coseno = cos( gen_cos->Uniform(0,2*M_PI) );

    result.pl = abs(p*coseno);
    result.pt = p*(sqrt(1-coseno*coseno));

    return result;

}