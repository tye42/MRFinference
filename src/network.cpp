#include "network.h"
#include <exception>
#include <string>
#include <sstream>
#include <algorithm>
#include <limits>
#include <fstream>
#include <ctime>
#include <random>

using namespace std;

network::network() {
    std::random_device rd;
    seed = rd();
}

network::network(unsigned int sd) {
    seed = sd;
}

void network::addfactor(const factor &f) {
    factor::scope s = f.getscope();
    // add variables in f to vars
    for(factor::scope::const_iterator it=s.begin(); it!=s.end(); it++){
        if(vars.find(it->first) == vars.end()){
            vars[it->first] = it->second;
        }
    }
    // add f to factors
    factors.push_back(f);
}

void network::readfromfile(const char *filename) {
    string line;
    int numFactors;

    ifstream file;
    file.open(filename);
    if(file.is_open()){
        while(file.peek() == '#'){
            getline(file, line);
        }
        file >> numFactors;
        getline(file, line);
        for(int i=0; i<numFactors; i++){
            int numVars;
            file >> numVars;
            vector<int> labels;
            vector<int> dims;
            for(int j=0; j<numVars; j++){
                int label;
                file >> label;
                labels.push_back(label);
            }
            for(int j=0; j<numVars; j++){
                int dim;
                file >> dim;
                dims.push_back(dim);
            }
            // add variables to vars
            factor::scope s;
            vector<int> stride;
            int stepsize = 1;
            for(int j=0; j<numVars; j++){
                s[labels[j]] = dims[j];
                stride.push_back(stepsize);
                stepsize *= dims[j];
                if(vars.find(labels[j]) == vars.end()){
                    vars[labels[j]] = dims[j];
                }
            }
            // add factor to factors
            factor f(s);
            int numEntries;
            file >> numEntries;
            for(int j=0; j<numEntries; j++){
                int idx;
                double val;
                file >> idx;
                file >> val;
                factor::assign a;
                for(int k=0; k<numVars; k++){
                    a[labels[k]] = (idx/stride[k])%dims[k];
                }
                f(a) = val;
            }
            factors.push_back(f);
        }
    }
    else{
        cout << "Cannot open file!" << endl;
    }
}

network::bethegraph network::createbethegraph() const {
    nbmap nbv, nbf;
    vector<edge> edges;
    for(int i = 0; i < factors.size(); i++){
        factor::scope s = factors[i].getscope();
        for(factor::scope::const_iterator j=s.begin(); j!=s.end(); j++){
            nbv[j->first].push_back(i);
            nbf[i].push_back(j->first);
            edges.push_back(make_pair(j->first, i));
        }
    }
    bethegraph ret;
    ret.nbv = nbv;
    ret.nbf = nbf;
    ret.edges = edges;
    return ret;
}

factor network::initmessage(int i) const {
    factor::scope s;
    s[i] = vars.find(i)->second;
    return factor(s, 1);
}

factor network::lbpupdatemsg(int vi, int fi, bethegraph& g, map<edge, factor>& messages) const {
    if(factors[fi].getscope().size() == 1){
        return factors[fi];
    }else{
        const factor::scope sempty;
        factor ret = factors[fi];
        vector<int>& nbfi = g.nbf.find(fi)->second;
        for(vector<int>::const_iterator it=nbfi.begin(); it!=nbfi.end(); it++){
            int vj = (*it);
            if(vj != vi){
                factor vjfi(sempty, 1);
                vector<int>& nbvj = g.nbv.find(vj)->second;
                for(vector<int>::const_iterator jit=nbvj.begin(); jit!=nbvj.end(); jit++){
                    if((*jit) != fi){
                        vjfi = vjfi * messages.find(make_pair(vj, *jit))->second;
                    }
                }
                ret = ret * vjfi;
            }
        }
        factor::scope s = ret.getscope();
        s.erase(vi);
        ret = ret.marginalize(s);
        return ret.normalize(); // need to normalzie the message to prevent overflow
    }
}

factor network::lbpbetav(int vi, map<edge, factor>& messages, bethegraph& g) const {
    const factor::scope sempty;
    factor ret(sempty, 1);
    vector<int>& nbvi = g.nbv.find(vi)->second;
    for(vector<int>::const_iterator it=nbvi.begin(); it!=nbvi.end(); it++){
        ret = ret * messages.find(make_pair(vi, *it))->second;
    }
    return ret.normalize();
}

factor network::lbpbetaf(int fi, std::map<edge, factor>& messages, bethegraph& g) const {
    const factor::scope sempty;
    factor ret = factors[fi];
    vector<int>& nbfi = g.nbf.find(fi)->second;
    for(vector<int>::const_iterator it=nbfi.begin(); it!=nbfi.end(); it++){
        int vj = (*it);
        factor vjfi(sempty, 1);
        vector<int>& nbvj = g.nbv.find(vj)->second;
        for(vector<int>::const_iterator jit=nbvj.begin(); jit!=nbvj.end(); jit++){
            if((*jit) != fi){
                vjfi = vjfi * messages.find(make_pair(vj, *jit))->second;
            }
        }
        ret = ret * vjfi;
    }
    return ret.normalize();
}

double network::lbpenergy(bethegraph &g, lbpresults &pt) const {
    double f = 0;
    for(map<int, factor>::const_iterator it=pt.betav.begin(); it!=pt.betav.end(); it++){
        f += (1.0 - g.nbv[it->first].size()) * (it->second).entropy();
    }
    for(int i=0; i<pt.betaf.size(); i++){
        f += pt.betaf[i].entropy();
        f += (pt.betaf[i] * factors[i].log()).sum();
    }
    return f;
}

network::lbpresults network::lbp(int maxiter, double tol) const {
    std::mt19937 eng(seed);
    lbpresults ret;
    int _iter = 0;
    ret.converged = 0;
    clock_t rt = clock();
    bethegraph g = createbethegraph();
    map<edge, factor> messages;
    vector<int> schedule; // index of edges to update
    for(int i=0; i<g.edges.size(); i++){
        // initial message
        messages.insert(make_pair(g.edges[i], initmessage(g.edges[i].first)));
        // initial schedule
        schedule.push_back(i);
    }
    // init betav
    for(map<int, int>::const_iterator it=vars.begin(); it!=vars.end(); it++){
        factor::scope s;
        s[it->first] = it->second;
        ret.betav.insert(make_pair(it->first, factor(s, 1)));
    }
    // init betaf
    for(int fi=0; fi<factors.size(); fi++){
        factor::scope s = factors[fi].getscope();
        ret.betaf.push_back(factor(s, 1));
    }

    double maxdiff = numeric_limits<double>::infinity();
    for(; _iter<maxiter && maxdiff>tol; _iter++){
        shuffle(schedule.begin(), schedule.end(), eng);
        maxdiff = - numeric_limits<double>::infinity();
        for(vector<int>::const_iterator it=schedule.begin(); it!=schedule.end(); it++){
            edge e = g.edges[*it];
            factor msg = lbpupdatemsg(e.first, e.second, g, messages);
            messages.find(e)->second = msg;
        }
        // update betav
        for(map<int, factor>::iterator it=ret.betav.begin(); it!=ret.betav.end(); it++){
            factor _b = lbpbetav(it->first, messages, g);
            maxdiff = max(maxdiff, _b.dist(it->second));
            it->second = _b;
        }
        // update betaf
        for(int fi=0; fi<factors.size(); fi++){
            factor _b = lbpbetaf(fi, messages, g);
            maxdiff = max(maxdiff, _b.dist(ret.betaf[fi]));
            ret.betaf[fi] = _b;
        }
    }
    ret.energy = lbpenergy(g, ret);
    rt = clock() - rt;
    ret.time = ((double)rt)/CLOCKS_PER_SEC;
    if(maxdiff <= tol){
        ret.converged = 1;
    }
    ret.iteration = _iter;
    return ret;
}

factor network::mfupdatevar(int i, vector<int>& nbi, map<int, factor>& beliefs) const {
    const factor::scope sempty;
    factor ret(sempty, 1);
    for(vector<int>::const_iterator it=nbi.begin(); it!=nbi.end(); it++){
        factor::scope s = factors[*it].getscope();
        factor uminusi(sempty, 1);
        for(factor::scope::const_iterator sit=s.begin(); sit!=s.end(); sit++){
            if(sit->first != i){
                uminusi = uminusi * (beliefs.find(sit->first)->second);
            }
        }
        s.erase(i);
        factor msg = (uminusi * factors[*it].log()).marginalize(s);
        ret = ret * msg.exp();
    }
    return ret.normalize();
}

double network::mfenergy(std::map<int, factor> &beliefs) const {
    double f = 0;
    for(map<int, factor>::const_iterator it=beliefs.begin(); it!=beliefs.end(); it++){
        f += (it->second).entropy();
    }
    const factor::scope sempty;
    for(int i=0; i<factors.size(); i++){
        factor::scope s = factors[i].getscope();
        factor q(sempty, 1);
        for(factor::scope::const_iterator it=s.begin(); it!=s.end(); it++){
            q = q * (beliefs.find(it->first)->second);
        }
        f += (q.normalize() * factors[i].log()).sum();
    }
    return f;
}

network::mfresults network::mf(bool uniform, int maxiter, double tol) const {
    mfresults ret;
    int _iter = 0;
    ret.converged = 0;
    clock_t rt;
    rt = clock();
    // find all the neighbors of each variable
    nbmap nb;
    for(int i = 0; i < factors.size(); i++){
        factor::scope s = factors[i].getscope();
        for(factor::scope::const_iterator j=s.begin(); j!=s.end(); j++){
            nb[j->first].push_back(i);
        }
    }
    vector<int> schedule;
    // initialize beliefs
    std::mt19937 eng(seed);
    uniform_real_distribution<double> dist(0.0, 1.0);
    for(map<int, int>::const_iterator it=vars.begin(); it!=vars.end(); it++){
        factor::scope s;
        s[it->first] = it->second;
        if(uniform) {
            ret.beliefs.insert(make_pair(it->first, factor(s, 1)));
        }
        else{
            // randomly init belief
            factor initb(s, 1);
            for(int i=0; i<initb.getsize(); i++){
                double val = dist(eng);
                while(val==0){
                    val = dist(eng);
                }
                initb.set(i, val);
            }
            ret.beliefs.insert(make_pair(it->first, initb));
        }
        schedule.push_back(it->first);
    }
    double maxdiff = numeric_limits<double>::infinity();
    for(; _iter < maxiter && maxdiff > tol; _iter++){
        shuffle(schedule.begin(), schedule.end(), eng);
        maxdiff = - numeric_limits<double>::infinity();
        for(vector<int>::const_iterator it=schedule.begin(); it!=schedule.end(); it++){
            factor _bf = mfupdatevar(*it, nb[*it], ret.beliefs);
            maxdiff = max(maxdiff, _bf.dist(ret.beliefs.find(*it)->second));
            ret.beliefs.find(*it)->second = _bf;
        }
    }
    ret.energy = mfenergy(ret.beliefs);
    rt = clock() - rt;
    ret.time = ((double)rt)/CLOCKS_PER_SEC;
    if(maxdiff <= tol){
        ret.converged = 1;
    }
    ret.iteration = _iter;
    return ret;
}
