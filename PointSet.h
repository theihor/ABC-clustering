#pragma once

#include <string>
#include <vector>

#define PI 3.14159265
#define _USE_MATH_DEFINES
#include <math.h> //M_PI

using std::string;
using std::vector;

#define SQR(x) (x) * (x) 
#define E_DIST(x1, y1, x2, y2) sqrt(SQR(x2 - x1) + SQR(y2 - y1))

struct Argb{
	int alpha,// always 255
 	r, //red
	g, //green
	b; //blue
};

class point{
public:
	point(double x, double y);

	double x, y;
	int iscentroid;//centroid mark
	Argb argb;//color scheme
	int cluster; 
};

class PointSet
{
public:
	PointSet(void);
	~PointSet(void);

	vector<point*>* points;
	int size;

	void load( string fname  );
	void dump2file( string fname );

	void paintWith(vector<double> *x);

	double silhouetteIndex(int K);
	double xieBeni(vector<double>* cnt);
};
