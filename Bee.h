#pragma once

#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>

extern double a;
extern double b;
extern unsigned int D;
extern unsigned int CS;
extern double L;

using namespace std;

class Bee {
private:
	vector<double>* position;
	unsigned int trialCounter;

public:
	Bee();
	~Bee();

	vector<double>* getPos();
	void setPos(vector<double>* newPos);

	void genRandomPos();

	int getTrialCounter();
	void setTrialCounter(int newTrialCounter);
	void incTrialCounter();
};