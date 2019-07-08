#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

struct generation{
    double pl;
    double pt;
};

class generate_mc{

    int seed;

    public:
        generate_mc(int s){ seed = s; };
        ~generate_mc(){};
        generation simulation();

};
#endif //MONTE_CARLO_H