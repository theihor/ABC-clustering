#include "Colony.h"
#include <time.h>
#include "PointSet.h"

using namespace std;

double hyperSphere(vector<double> *x) {
	double f = 0;

	for (unsigned int i = 0; i < x->size(); i++)
		f += x->at(i) * x->at(i);

	return f;
}

double aphef(vector<double> *x)
{
	double f = 0;

	for (int i = 0; i < x->size(); i++)
		f += (i + 1) * x->at(i) * x->at(i);

	return f;
}

double rhef(vector<double> *x)
{
	double f = 0, sum = 0;

	for (int i = 0; i < x->size(); i++)
	{
		sum = 0;
		for (int j = 0; j <= i; j++) sum += x->at(j);
		f += sum * sum;
	}

	return f;
}
double sqr (double x) {
	return x * x;
}


double rosenbrock (vector<double> *x)
{
	double f = 0;

	for (int i = 0; i < x->size() - 1; i++) {
		f += 100 * sqr(sqr(x->at(i)) - x->at(i + 1)) +
			sqr(1 - x->at(i));
	}

	return f;
}

// Branin


double branin (vector<double> *x) {
	return sqr(x->at(1) - 5.1 * x->at(0) / (4 * sqr(PI)) +
		5 * x->at(0) / PI - 6.0) + 
		10.0 * (1 - 1 / (8 * PI)) * cos(x->at(0)) + 10.0;
}

double zero (vector<double> *x) {
	return 0;
}

/*
int main() {
	system("color a");

	srand(time(NULL));

	double (*f)(vector<double>*);
	f = hyperSphere;//branin;

	Colony c(f);
	double min = f(c.getEverBestSolution());
	cout << "i = -1" << " min = " << min << endl;
	for (int i = 0; i < 300; i++) {
		c.run();
		//cout << "vector: ";
		if (min > f(c.getEverBestSolution())) {
			min = f(c.getEverBestSolution());
			cout << "i = " << i << " min = " << min << endl;
			showVector(c.getEverBestSolution());
		}
	}
	return 0;
}
*/

int main() {
	PointSet *ps = new PointSet();
	string fname = "Skewdistribution_3.dat";
	ps->load(fname);


	a = ps->points->at(0)->x;
	if (a > ps->points->at(0)->y) { a = ps->points->at(0)->x; }
	for (int k = 1; k < ps->size; k++) {
		point *p = ps->points->at(k);
		if(a > p->x) { a = p->x; }
		if(a > p->y) { a = p->y; }
	}

	b = ps->points->at(0)->x;
	if (b < ps->points->at(0)->y) { b = ps->points->at(0)->x; }
	for (int k = 1; k < ps->size; k++) {
		point *p = ps->points->at(k);
		if(b < p->x) { b = p->x; }
		if(b < p->y) { b = p->y; }
	}

	D = 3 * 2;
	CS = 20;
	L = CS * D / 2.0;

	Colony c(ps);

	double max = c.getFit(c.getEverBestSolution());
	vector<double>* best = c.getEverBestSolution();
	cout << "i = -1" << " best fitness = " << max << endl;
	for (int i = 0; i < 16; i++) {
		c.run();
		//cout << "i = " << i << endl;
		if (max < c.getFit(c.getEverBestSolution())) {
			max = c.getFit(c.getEverBestSolution());
			cout << "i = " << i << " best fitness = " << max << endl;
			best = c.getEverBestSolution();
			showVector(c.getEverBestSolution());
		}
	}
	
	ps->paintWith(best);
	cout << "SSE: " << c.sse(best) << endl;
	cout << "Silhouette index: " << ps->silhouetteIndex(D / 2) << endl;
	cout << "Index Xie-Beni: " << ps->xieBeni(best) << endl;

	vector<point*>* points = new vector<point*>(ps->size + D / 2);
	for (int i = 0; i < ps->size; i++) {
		points->at(i) = ps->points->at(i);
	}

	for (int k = 0; k < D / 2; k++) {
		point *p = new point(best->at(k*2), best->at(k*2+1));
		Argb argb;
		argb.alpha = 255;
		argb.b = 255;
		argb.g = 255;
		argb.r = 255;
		p->iscentroid = 1;
		p->argb = argb;
		points->at(ps->size + k) = p;
	}

	ps->points = points;
	ps->size = points->size();

	ps->dump2file(fname + ".out.dat");
	return 0;
}