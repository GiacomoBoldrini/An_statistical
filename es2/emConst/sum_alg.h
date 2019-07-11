 #ifndef SUM_
 #define SUM_

#include <iostream>
#include <cmath>
#include <numeric>
#include <limits>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_math.h>

template <typename float_num, typename int_num>
class sum{

    private: 
        float_num sum_in;
        float_num c;

    public:
        sum(){sum_in = 0;};
        ~sum(){};
        float_num norm_sum(int_num N = 100){
            
            sum_in = 0;
            for(int i = 1; i < N; i++ ){
                
                sum_in += 1./(float_num)i;                                
            } 
            return sum_in;

        };

        float_num kahan_sum(int_num N = 100){

            sum_in = 0;
            c = 0;
            for(int i = 1; i < N; i++){

                float_num y = 1./((float_num)i) - c;
                float_num t = sum_in + y;
                c = (t-sum_in) - y;
                sum_in = t;
            }
            return sum_in;

        };

        float_num pairwise_sum(int start, int stop){
            if ((stop - start ) < 100){ //Solo se numero è abbastanza grande
                 sum_in = 0.;
                    for (int i = start; i<=stop; i++){
                        sum_in += 1 / (double)i;
                    }
                return sum_in;
             }
             else{
                int m = floor((stop - start) /2 );
                return pairwise_sum(start, start+m) + pairwise_sum(start+m+1, stop);
             }
        
        
        }
};



 #endif //SUM_
