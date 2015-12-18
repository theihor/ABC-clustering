#include "PointSet.h"
#include <fstream>
#include <iostream>

using namespace std;

point::point(double x, double y) {
	this->x = x;
	this->y = y;
}

PointSet::PointSet(void)
{
	size = 0;
}

PointSet::~PointSet(void)
{
}

void PointSet::load( string fname  ) {
	ifstream file;
	file.open(fname.c_str(), ios::in|ios::binary);

	if (file.fail()) {
		cerr << "unable to open file input.dat for reading" << std::endl;
		exit(1);
    }

	file.read((char*)&size, sizeof(int));
	cout << "Number of points: " << size << endl;

	points = new vector<point*>(size);

	for(int i = 0; i < size; i++){

		double x, y;
		file.read((char*)&x, sizeof(double));
		file.read((char*)&y, sizeof(double));

		int alpha, r, g, b, cnt;
		file.read((char*)&alpha, sizeof(int));
		file.read((char*)&r, sizeof(int));
		file.read((char*)&g, sizeof(int));
		file.read((char*)&b, sizeof(int));
		file.read((char*)&cnt, sizeof(int));

		point *p = new point(x, y);
		p->argb.alpha = alpha;
		p->argb.r = r;
		p->argb.g = g;
		p->argb.b = b;
		p->iscentroid = cnt;
		points->at(i) = p;	
	}

	file.close();
}

void PointSet::dump2file( string fname ){
	ofstream file;
	file.open(fname.c_str(), ios::out|ios::binary);

	file.write((char*)&size, sizeof(int));

	for(int i = 0; i < size; i++){
		point *p = points->at(i);

		file.write((char*)&p->x, sizeof(double));
		file.write((char*)&p->y, sizeof(double));
		file.write((char*)&p->argb.alpha, sizeof(int));
		file.write((char*)&p->argb.r, sizeof(int));
		file.write((char*)&p->argb.g, sizeof(int));
		file.write((char*)&p->argb.b, sizeof(int));
		file.write((char*)&p->iscentroid, sizeof(int));
	}
	file.close();
}

void PointSet::paintWith(vector<double> *x) {
	int K = x->size() / 2;
	vector<Argb> colors(K);
	for(int i = 0; i < K; i++) {
		Argb argb;
		argb.alpha = 255;
		argb.r = i * 255 / K;
		argb.g = 255 - i * 255 / K;
		argb.b = i * 255 / K;
		colors.at(i) = argb;
	}

	// for each point from input PointSet
	for (int k = 0; k < this->size; k++) {
		point *p = points->at(k);
		// looking for nearest potential centroid
		double px = points->at(k)->x;
		double py = points->at(k)->y;
		double cx = x->at(0);
		double cy = x->at(1);
		double d = E_DIST(px, py, cx, cy);

		double color_index = 0;
		for (int i = 1; i < K; i++) {
			double new_d = E_DIST(px, py, x->at(i*2), x->at(i*2+1));
			if (new_d < d) {
				cx = x->at(i * 2);
				cy = x->at(i * 2 + 1);
				d = new_d;
				color_index = i;
			}
		}
		p->argb = colors.at(color_index);
		p->cluster = color_index;
	}
}

double point_dist(point *p1, point *p2) {
	return E_DIST(p1->x, p1->y, p2->x, p2->y);
}

double PointSet::silhouetteIndex(int K) {
	vector<double> res(K, 0);

	vector<double> cluster_sizes(K, 0);

	for (int i = 0; i < this->size; i++) {
		cluster_sizes.at(points->at(i)->cluster) += 1;
	}

	for (int k = 0; k < K; k++) {
		for (int i = 0; i < this->size; i++) {
			double ai = 0;
			vector<double> bis(K, 0);
			point *pi = points->at(i);
			for (int j = 0; j < this->size; j++) {
				if (j == i) continue;
				point *pj = points->at(j);
				if (pj->cluster == pi->cluster) {
					ai += point_dist(pi, pj);
				} else {
					bis.at(pj->cluster) += point_dist(pi, pj);
				}
			}
			ai /= cluster_sizes.at(pi->cluster);

			double bi;
			if (pi->cluster == 0) bi = bis.at(1) / cluster_sizes.at(1);
			else bi = bis.at(0) / cluster_sizes.at(0);
			for (int j = 1; j < K; j++) {
				if (j == pi->cluster) continue;
				bis.at(j) /= cluster_sizes.at(j);
				if (bis.at(j) < bi) bi = bis.at(j);
				//cout << "bis/at(" << j << ") = " << bis.at(j) << endl;
			}

			//cout << "ai = " << ai << endl;
			//cout << "bi = " << bi << endl;
			double si = bi - ai;
			if (bi > ai) si /= bi;
			else si /= ai;
			
			res.at(k) += si;
		}
	}

	double sum = 0;
	for (int i = 0; i < res.size(); i++) sum += res.at(i);

	return sum;
}

double PointSet::xieBeni(vector<double>* cnt) {
	int K = cnt->size() / 2;

	double msdcc;
	msdcc = SQR(E_DIST(cnt->at(0), cnt->at(1), cnt->at(2), cnt->at(3)));
	for (int i = 0; i < K; i++) {
		for (int j = i + 2; j < K; j++) {
			double sdcc = SQR(E_DIST(cnt->at(i*2), cnt->at(i*2+1), cnt->at(j*2), cnt->at(j*2+1)));
			if (sdcc < msdcc) msdcc = sdcc;
		}
	}

	double sum = 0;
	for (int i = 0; i < this->size; i++) {
		point *p = points->at(i);
		double cx = cnt->at(p->cluster*2);
		double cy = cnt->at(p->cluster*2+1);
		sum += SQR(E_DIST(p->x, p->y, cx, cy));
	}
	double msdoc = sum / this->size;

	double xb = msdoc / (this->size * msdcc);

	return xb;
}