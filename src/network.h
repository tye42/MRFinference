#ifndef NETWORK_H
#define NETWORK_H

#include "factor.h"
#include <map>
#include <set>
#include <vector>
#include <random>

class network { // class representing a Bayesian or Markov network
public:
    network();
    explicit network(unsigned int sd);
    typedef std::map<int, std::vector<int> > nbmap;
    // first is label of variable, second is index of factor
    typedef std::pair<int, int> edge;
    // f needs to be a factor that includes v in its scope
    // (and represents the conditional distribution of v)
    void addfactor(const factor &f);
    void readfromfile(const char *filename);

    unsigned int getseed() const { return seed; };
    std::vector<factor> getfactors() const { return factors; };
    std::map<int, int> getvars() const { return vars; };
    int getnumvars() const { return vars.size(); };

    class bethegraph {
    public:
        nbmap nbv; // neighbors of variable node
        nbmap nbf; // neighbors of factor node
        std::vector<edge> edges;
    };

    // mean field returned results
    class mfresults {
    public:
        double time;
        int converged;
        int iteration;
        double energy;
        std::map<int, factor> beliefs;
    };

    // loop belief returned results
    class lbpresults {
    public:
        double time;
        int converged;
        int iteration;
        double energy;
        std::map<int, factor> betav; // belief for each variable
        std::vector<factor> betaf; // belief for each factor
    };

    // run mean field algorithm on network
    mfresults mf(bool uniform=true, int maxiter=1000, double tol=1e-9) const;

    // run loopy belief propagation on bethe graph
    lbpresults lbp(int maxiter=1000, double tol=1e-9) const;

private:
    std::vector<factor> factors;
    std::map<int, int> vars;
    unsigned int seed;

    factor mfupdatevar(int i, std::vector<int>& nbi, std::map<int, factor>& beliefs) const;
    double mfenergy(std::map<int, factor>& beliefs) const;

    bethegraph createbethegraph() const;
    factor initmessage(int i) const;
    factor lbpupdatemsg(int vi, int fi, bethegraph& g, std::map<edge, factor>& messages) const;
    factor lbpbetav(int vi, std::map<edge, factor>& messages, bethegraph& g) const;
    factor lbpbetaf(int fi, std::map<edge, factor>& messages, bethegraph& g) const;
    double lbpenergy(bethegraph& g, lbpresults& pt) const;
};

#endif
