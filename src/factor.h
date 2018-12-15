#ifndef FACTOR_H
#define FACTOR_H

#include <iostream>
#include <map>
#include <cstdarg>
#include <vector>

class factor {
public:
	typedef std::map<int,int> scope;
	typedef std::map<int,int> assign;

	factor(const scope &sc, double d=0.0);

	double operator()(const assign &i) const;
	double &operator()(const assign &i);

	// get and set value of the i-th entry
	double get(size_t i) const { return  values[i]; }
	void set(size_t i, double val) { values[i] = val; }
	int getsize() const { return values.size(); }

	// fill the entries with same value
	factor &fill(double val);

	// reduce the factor to only those assignments consistent with a
	factor reduce(const assign &a) const;
	// sum out those variables mentioned in the scope tosumout
	factor marginalize(const scope &tosumout) const;

	factor operator*(const factor &f2) const;
	factor operator/(const factor &f2) const;
	factor operator+(const factor &f2) const;
	factor operator-(const factor &f2) const;

	// take the log of each entry
	factor log() const;
	// take the exp of each entry
	factor exp() const;
	// calculate entropy - sum {p * log(p)}
	double entropy() const;
	double sum() const;

    // calculate distance between two factors
    double dist(const factor &f2) const;

	// normalize the factor
	factor normalize() const;

	factor &operator*=(double d);
	factor &operator/=(double d);
	factor &operator+=(double d);
	factor &operator-=(double d);

	// returns the scope of this factor
	scope getscope() const { return scp; }

	void print(std::ostream &os) const;
private:
	scope scp;
	std::vector<double> values;
	std::map<int,int> stride;
	int size;

	int getsize(const scope &s) const;
	int getindex(const assign &i) const;
	scope scopeunion(const scope &s2) const;
	scope scopeminus(const std::map<int,int> &a) const;
	assign initassign(const scope &s, int val) const;
	void nextassign(int &index, const scope &s, assign &a) const;
	void nextassign(int &j, int &k, const scope &s, assign &a, const factor &f2) const;
};

#endif
