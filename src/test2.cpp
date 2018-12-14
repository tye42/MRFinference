#include "network.h"
#include <cstdarg>
#include <algorithm>
#include <fstream>

using namespace std;

void printmfresults(network::mfresults& ret, ofstream& out){
    out << "[MF]" << endl;
    out << "TIME:" << ret.time << endl;
    out << "NITER:" << ret.iteration << endl;
    out << "CONVERGED:" << ret.converged << endl;
    out << "ENERGY:" << ret.energy << endl;
}

void printlbpresults(network::lbpresults& ret, ofstream& out){
    out << "[LBP]" << endl;
    out << "TIME:" << ret.time << endl;
    out << "NITER:" << ret.iteration << endl;
    out << "CONVERGED:" << ret.converged << endl;
    out << "ENERGY:" << ret.energy << endl;
}


int main(int argc, char *argv[]){
    network net;
    net.readfromfile(argv[1]);
    network::mfresults mfret = net.mf();
    ofstream out(argv[2]);
    out << "GRAPH:" << argv[1] << endl;
    out << "SEED:" << net.getseed() << endl;
    printmfresults(mfret, out);
    cout << "mf time:" << mfret.time << endl;
    cout << "mf niter:" << mfret.iteration << endl;
    cout << "mf converged:" << mfret.converged << endl;
    cout << "mf energy:" << mfret.energy << endl;
//    for(map<int, factor>::const_iterator it=mfret.beliefs.begin(); it!=mfret.beliefs.end(); it++){
//        (it->second).print(cout);
//    }
    network::lbpresults lbpret = net.lbp();
    printlbpresults(lbpret, out);
    cout << "lbp time:" << lbpret.time << endl;
    cout << "lbp niter:" << lbpret.iteration << endl;
    cout << "lbp converged:" << lbpret.converged << endl;
    cout << "lbp energy:" << lbpret.energy << endl;
//    for(map<int, factor>::const_iterator it=lbpret.betav.begin(); it!=lbpret.betav.end(); it++){
//        (it->second).print(cout);
//    }
    out.close();
    return 0;
}