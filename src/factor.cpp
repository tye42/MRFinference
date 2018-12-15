#include "factor.h"
#include <assert.h>
#include <string.h>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <numeric>

using namespace std;

factor::factor(const scope &sc, double d) {
	scp = sc;
	size = getsize(sc);
	values.assign(size, d);
	int step_size = 1;
	for(scope::const_iterator it=sc.begin(); it!=sc.end(); it++){
		stride[it->first] = step_size;
		step_size *= it->second;
	}
}

double factor::operator()(const factor::assign &i) const {
	return values[getindex(i)];
}

double &factor::operator()(const factor::assign &i) {
	return values[getindex(i)];
}

factor &factor::fill(double val){
	std::fill(values.begin(), values.end(), val);
	return *this;
}

void factor::print(ostream &os) const {
	for(scope::const_iterator it=scp.begin(); it!=scp.end(); it++){
		os << 'X' << it->first << string(3, ' ');
	}
	os << string(3, ' ') << "value" << endl;
	os << string(5*scp.size()+16, '-') << endl;
	for(int i=0; i<size; i++){
		for(scope::const_iterator it=scp.begin(); it!=scp.end(); it++){
			os << string(3, ' ') << (i/stride.find(it->first)->second)%it->second << ' ';
		}
		os << string(4, ' ') << values[i] << endl;
	}
}

factor factor::reduce(const assign &a) const {
	scope scr = scopeminus(a);
	factor fr(scr);
	assign ar = initassign(scr, 0);
	int index = 0;
	// fix the variables assignment to a
	for(assign::const_iterator it=a.begin(); it!=a.end(); it++){
		index += it->second * stride.find(it->first)->second;
	}
	// set each row of the reduced factor
	for(int i=0; i<fr.size; i++){
		fr.values[i] = values[index];
		// step to the next assignment
		nextassign(index, scr, ar);
	}
	return fr;
}

factor factor::marginalize(const scope &tosumout) const {
	scope scm = scopeminus(tosumout);
	factor fm(scm);
	assign am = initassign(scp, 0);
	int index = 0;
	int excludesize = getsize(tosumout);
	// linear in inputsize = fm.size * excludesize
	for(int i=0; i<fm.size; i++){
		double margin = 0;
		for(int j=0; j<excludesize; j++){
			margin += values[index];
			nextassign(index, tosumout, am);
		}
		fm.values[i] = margin;
		nextassign(index, scm, am);
	}
	return fm;
}

factor factor::operator*(const factor &f2) const {
	scope scu = scopeunion(f2.scp);
	factor prod(scu);
	assign a = initassign(scu, 0);
	int j=0, k=0;
	for(int i=0; i<prod.size; i++){
		prod.values[i] = values[j] * f2.values[k];
		nextassign(j, k, scu, a, f2);
	}
	return prod;
}

factor factor::operator/(const factor &f2) const {
	scope scu = scopeunion(f2.scp);
	factor div(scu);
	assign a = initassign(scu, 0);
	int j=0, k=0;
	for(int i=0; i<div.size; i++){
		div.values[i] = values[j] / f2.values[k];
		nextassign(j, k, scu, a, f2);
	}
	return div;
}

factor factor::operator+(const factor &f2) const {
	scope scu = scopeunion(f2.scp);
	factor sum(scu);
	assign a = initassign(scu, 0);
	int j=0, k=0;
	for(int i=0; i<sum.size; i++){
		sum.values[i] = values[j] + f2.values[k];
		nextassign(j, k, scu, a, f2);
	}
	return sum;
}

factor factor::operator-(const factor &f2) const {
	scope scu = scopeunion(f2.scp);
	factor sub(scu);
	assign a = initassign(scu, 0);
	int j=0, k=0;
	for(int i=0; i<sub.size; i++){
		sub.values[i] = values[j] - f2.values[k];
		nextassign(j, k, scu, a, f2);
	}
	return sub;
}

factor factor::log() const {
	factor ret(scp);
	std::transform(values.begin(), values.end(), ret.values.begin(), (double(*)(double))std::log);
	return ret;
}

factor factor::exp() const {
	factor ret(scp);
	std::transform(values.begin(), values.end(), ret.values.begin(), (double(*)(double))std::exp);
	return ret;
}

double factor::entropy() const {
	double h = 0;
	for(int i=0; i<values.size(); i++){
		h += values[i] * std::log(values[i]);
	}
	return -h;
}

double factor::sum() const {
	double s = 0;
	for(int i=0; i<values.size(); i++){
		s += values[i];
	}
	return s;
}

double factor::dist(const factor &f2) const {
	double diff = 0;
	for(int i=0; i<values.size(); i++){
		diff = max(diff, abs(values[i] - f2.values[i]));
	}
	return diff;
}

factor factor::normalize() const {
	double z = 0;
	for(vector<double>::const_iterator it = values.begin(); it!= values.end(); it++){
		z += (*it);
	}
	factor ret = *this;
	ret /= z;
	return ret;
}

factor &factor::operator*=(double d) {
	for(int i=0;i<size;i++) values[i] *= d;
	return *this;
}
factor &factor::operator/=(double d) {
	for(int i=0;i<size;i++) values[i] /= d;
	return *this;
}
factor &factor::operator+=(double d) {
	for(int i=0;i<size;i++) values[i] += d;
	return *this;
}
factor &factor::operator-=(double d) {
	for(int i=0;i<size;i++) values[i] -= d;
	return *this;
}

// helper functions

// get the size of the scope, which is # rows in the table
int factor::getsize(const scope &s) const {
	int sz = 1;
	for(scope::const_iterator it=s.begin(); it!=s.end(); it++){
		sz *= it->second;
	}
	// if the scope is empty, size = 1
	return sz;
}

// given an assignment, find the index of value,
// extra variables are ignored
int factor::getindex(const assign &i) const {
	int index = 0;
	for(scope::const_iterator it=scp.begin(); it!=scp.end(); it++){
		index += i.find(it->first)->second * stride.find(it->first)->second;
	}
	// if scope is empty, index = 0
	return index;
}

// return the union with scope s2
factor::scope factor::scopeunion(const scope &s2) const {
	scope scu;
	for(scope::const_iterator it=scp.begin(); it!=scp.end(); it++){
		scu[it->first] = it->second;
	}
	for(scope::const_iterator it=s2.begin(); it!=s2.end(); it++){
		if(scu.count(it->first)==0){
			scu[it->first] = it->second;
		}
	}
	return scu;
}

// exclude the variables in assignment or scope
factor::scope factor::scopeminus(const map<int,int> &a) const {
	scope remain;
	for(scope::const_iterator it=scp.begin(); it!=scp.end(); it++){
		if(a.count(it->first)==0){
			remain[it->first] = it->second;
		}
	}
	return remain;
}

// initialize a assign of scope s with assignment val
factor::assign factor::initassign(const scope &s, int val) const {
	assign a;
	for(scope::const_iterator it=s.begin(); it!=s.end(); it++){
		a[it->first] = val;
	}
	return a;
}

// step to the next assignment, change the variables in scope, update current index
void factor::nextassign(int &index, const scope &s, assign &a) const {
	for(scope::const_iterator it=s.begin(); it!=s.end(); it++){
		a[it->first] += 1;
		if(a[it->first]==it->second){
			a[it->first] = 0;
			index -= (it->second-1)*stride.find(it->first)->second;
		}else{
			index += stride.find(it->first)->second;
			break;
		}
	}
}

void factor::nextassign(int &j, int &k, const scope &s, assign &a, const factor &f2) const {
	for(scope::const_iterator it=s.begin(); it!=s.end(); it++){
		a[it->first] += 1;
		if(a[it->first] == it->second){
			a[it->first] = 0;
			if(scp.count(it->first)>0){
				j -= (it->second-1)*stride.find(it->first)->second;
			}
			if(f2.scp.count(it->first)>0){
				k -= (it->second-1)*f2.stride.find(it->first)->second;
			}
		}else{
			if(scp.count(it->first)>0){
				j += stride.find(it->first)->second;
			}
			if(f2.scp.count(it->first)>0){
				k += f2.stride.find(it->first)->second;
			}
			break;
		}
	}
}
