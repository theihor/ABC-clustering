#include "Bee.h"

double a = -5.0;
double b = 5.0;
unsigned int D = 8;
unsigned int CS = 50;
double L = CS * D / 2.0;

Bee::Bee() {
	setTrialCounter(0);
	this->position = new vector<double>(D);
}

Bee::~Bee() {
	delete this->position;
}

vector<double>* Bee::getPos() { return position; }

void Bee::setPos(vector<double>* newPos) { position = newPos;}

int Bee::getTrialCounter() { return trialCounter; }

void Bee::setTrialCounter(int newTrialCounter) { trialCounter = newTrialCounter; }

void Bee::incTrialCounter() { trialCounter++; }

void Bee::genRandomPos() {
	getPos()->resize(D);

	for (unsigned int i = 0; i < getPos()->size(); i++)
		getPos()->at(i) = a + (b - a) * (double)rand() / RAND_MAX;
}