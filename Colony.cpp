#include "Colony.h"

using std::cout;

void showVector(vector<double> *x) {
	for (unsigned int i = 0; i < x->size(); i++) {
		cout << x->at(i) << " ";
	}
	cout << endl;
}

Colony::Colony(double (*func)(vector<double>*)) {
	//cout << "Colony constructor\n";
	f = func;

	X.resize(CS);
	
	for (int i = 0; i < CS; i++) {
		X.at(i) = new Bee();
		X.at(i)->genRandomPos();
	}
	
	fit.resize(CS / 2);
	genFitnessVector();

	p.resize(CS / 2);
	genProbabilities();

	best.resize(D);
	best = *getCurBestSolution();

	minimums.resize(D);
	maximums.resize(D);
	updateMinsAndMaxs();
}

Colony::Colony(PointSet *ps) {
	//cout << "Colony constructor\n";
	//f = func;
	this->ps = ps;

	X.resize(CS);
	
	for (int i = 0; i < CS; i++) {
		X.at(i) = new Bee();
		vector<double>* pos = new vector<double>(D);
		for (int j = 0; j < D / 2; j++) {
			point* p = ps->points->at(rand() % ps->size);
			pos->at(j*2) = p->x;
			pos->at(j*2+1) = p->y;
		}
		X.at(i)->setPos(pos);
	}
	
	fit.resize(CS / 2);
	genFitnessVector();

	p.resize(CS / 2);
	genProbabilities();

	abstractBee = new vector<double>(D);
	*abstractBee = *X.at(0)->getPos();
	best.resize(D);
	best = *getCurBestSolution();

	minimums.resize(D);
	maximums.resize(D);
	updateMinsAndMaxs();

	*abstractBee = best;
}


Colony::~Colony() {
	//cout << "Colony destructor\n";
	for (int i = 0; i < X.size(); i++)
		delete X.at(i);
}


double Colony::getFit(Bee *bee) {
	return getFit(bee->getPos());
}

double Colony::sse(vector<double> *x) {
	// SSE
	double sum = 0;

	// for each point from input PointSet
	for (int k = 0; k < ps->size; k++) {

		// looking for nearest potential centroid
		double px = ps->points->at(k)->x;
		double py = ps->points->at(k)->y;
		double cx = x->at(0);
		double cy = x->at(1);
		double d = E_DIST(px, py, cx, cy);
		for (int i = 1; i < D / 2; i++) {
			double new_d = E_DIST(px, py, x->at(i*2), x->at(i*2+1));
			if (new_d < d) {
				cx = x->at(i*2);
				cy = x->at(i*2 + 1);
				d = new_d;
			}
		}	

		sum += SQR(d);
	}

	return sum;
}

double Colony::getFit(vector<double> *x) {
	double sum = sse(x);
	return 1.0 / (1.0 + sum);
}

/*
double Colony::getFit(vector<double> *x) {
	//cout << "Colony::getFit(vector<double> *x)\n";
	double curF = f(x);
	if(curF >= 0)
		return 1.0 / (1 + curF);
	else
		return 1.0 + abs(curF);
}
*/

void Colony::genFitnessVector() {
	//cout << "Colony::genFitnessVector()\n";
	for (int i = 0; i < fit.size(); i++)
		fit.at(i) = getFit(X.at(i));
}

void Colony::genProbabilities() {
	//cout << "Colony::genProbabilities()\n";
	// generating p-vector
	double fitSum = 0;
	for (int i = 0; i < CS / 2; i++)
		fitSum += fit.at(i);
			
	for (int i = 0; i < CS / 2; i++) {
		p.at(i) = fit.at(i) / fitSum;
	}
}

unsigned int Colony::chooseRandomBee() {
	//cout << "Colony::chooseRandomBee()\n";
	double P = (double)rand() / RAND_MAX;
	unsigned int j = 0; double w = p.at(j);
	while (P >= w) {
		w += p.at(j);
		j++;
	}
	//j--; // random selected bee = X.at(j)

	return j;
}

void Colony::updateBee(int i, int j, double value) {
	vector<double>* pos = X.at(i)->getPos();
	
	// processing abstract bee
	int k = rand() % (D / 2);
	double x = abstractBee->at(k*2);
	double y = abstractBee->at(k*2+1);

	double fitness1 = getFit(abstractBee);
	abstractBee->at(k*2) = pos->at(k*2);
	abstractBee->at(k*2+1) = pos->at(k*2+1);
	double fitness2 = getFit(abstractBee);
	if (fitness1 > fitness2) {
		abstractBee->at(k*2) = x;
		abstractBee->at(k*2+1) = y;
	}

	// modifying bee
	X.at(i)->getPos()->at(j) = value;
	this->updateMinsAndMaxs(value, j);
	X.at(i)->setTrialCounter(0);
}


void Colony::resetBee(int i) {
	vector<double>* newPos = new vector<double>(D);

	for (int j = 0; j < D; j++) {
		newPos->at(j) = minimums.at(j) + 
			(1.0 * rand() / RAND_MAX) * (maximums.at(j) - minimums.at(j));
	}
	X.at(i)->setPos(newPos);
	X.at(i)->setTrialCounter(0);
}

/*
void Colony::resetBee(int i) {
	vector<double>* newPos = new vector<double>(D);

	vector<double>* pos = new vector<double>(D);
	for (int j = 0; j < D / 2; j++) {
		point* p = ps->points->at(rand() % ps->size);
		pos->at(j*2) = p->x;
		pos->at(j*2+1) = p->y;
	}

	X.at(i)->setPos(pos);
	X.at(i)->setTrialCounter(0);
}
*/
void Colony::run() {
	//cout << "Colony::run()\n";
	vector<double> V;
	V.resize(D);
	
	unsigned int j, k;
	double R;

	// employee phase
	for (int i = 0; i < CS / 2; i++) {
		V = *X.at(i)->getPos();
		
		// V = Xij + R*(Xij-Xkj)
		j = rand() % D;
		k = rand() % ((CS / 2) - 1);
		if (k >= i) {
			k += 1;
		}
		R = -1.0 + 2.0 * rand() / RAND_MAX;
			
		V.at(j) += R * (V.at(j) - X.at(k)->getPos()->at(j));
	
		if (getFit(&V) > getFit(X.at(i))) {
			// improve
			updateBee(i, j, V.at(j));
		} else X.at(i)->incTrialCounter();
	}
	
	// onlookers' phase
	genFitnessVector();
	genProbabilities();

	for (int i = CS / 2; i < CS; i++) {
		V = *X.at(i)->getPos();
	
		// V = Xij + R*(Xij-Xkj)
		j = rand() % D;
		k = chooseRandomBee();
		R = -1.0 + 2.0 * rand() / RAND_MAX;
			
		V.at(j) += R * (V.at(j) - X.at(k)->getPos()->at(j));

		if (getFit(&V) > getFit(X.at(i))) {
			// improve
			updateBee(i, j, V.at(j));
		} else X.at(i)->incTrialCounter();

	}

	for (int i = 0; i < CS / 2; i++) {
		if (X.at(i)->getTrialCounter() > L) {
			resetBee(i);
		}
	}

	genFitnessVector();
	genProbabilities();

	vector<double>* curBest = getCurBestSolution();
	if (getFit(&best) < getFit(curBest)) {
		best = *curBest;
	}
}

vector<double>* Colony::getCurBestSolution() {
	//cout << "Colony::getCurBestSolution()\n";
	unsigned int i, index = -1;
	double max = getFit(abstractBee);

	for (i = 0; i < CS; i++) {
		double fitness = getFit(X.at(i)->getPos());
		if (fitness > max) {
			max = fitness;
			index = i;
		}
	}
	
	if (index == -1) return abstractBee;
	else return X.at(index)->getPos();
}

vector<double>* Colony::getEverBestSolution() { return &best; }

void Colony::updateMinsAndMaxs() {
	for (int j = 0; j < D; j++) {
		minimums.at(j) = X.at(0)->getPos()->at(j);
		maximums.at(j) = X.at(0)->getPos()->at(j);
	}

	for (int i = 1; i < CS; i++) {
		updateMinsAndMaxs(X.at(i)->getPos());
	}
}

void Colony::updateMinsAndMaxs(vector<double>* x) {
	for (int j = 0; j < D; j++) {
		updateMinsAndMaxs(x->at(j), j);
	}
}

void Colony::updateMinsAndMaxs(double value, int j) {
	if (minimums.at(j) > value) {
		minimums.at(j) = value;
	}
	if (maximums.at(j) < value) {
		maximums.at(j) = value;
	}
}