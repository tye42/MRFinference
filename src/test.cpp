#include "network.h"
#include <cstdarg>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[]){
    network net;
    net.readfromfile(argv[1]);
    network::mfresults mfret = net.mf();

    cout << "MF running time: " << mfret.time << endl;
    cout << "MF number of iteration: " << mfret.iteration << endl;
    cout << "MF converged: " << mfret.converged << endl;
    cout << "MF energy: " << mfret.energy << endl;
    cout << endl;

    network::lbpresults lbpret = net.lbp();
    cout << "LBP running time: " << lbpret.time << endl;
    cout << "LBP number of iteration: " << lbpret.iteration << endl;
    cout << "LBP converged: " << lbpret.converged << endl;
    cout << "LBP energy: " << lbpret.energy << endl;
    cout << endl;
    return 0;
}