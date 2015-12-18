#include "Bee.h"
#include "PointSet.h"

void showVector(vector<double> *x);

class Colony {
private:
	vector<Bee*> X;
	vector<double> p;
	vector<double> fit;
	vector<double> best;
	vector<double> minimums;
	vector<double> maximums;

	vector<double> *abstractBee;

	double (*f)(vector<double>*);
	
	double getFit(Bee *bee);
	
	void genFitnessVector();	
	void genProbabilities();
	unsigned int chooseRandomBee();

	void updateMinsAndMaxs();
	void updateMinsAndMaxs(vector<double>* x);
	void updateMinsAndMaxs(double value, int j);

	void updateBee(int i, int j, double value);
	void resetBee(int i);

public:
	Colony(double (*func)(vector<double>*));
	Colony(PointSet *ps);
	~Colony();
	
	PointSet *ps;

	void run();
	
	vector<double>* getCurBestSolution();
	vector<double>* getEverBestSolution();
	double getFit(vector<double>* x);
	double sse(vector<double> *x);
};