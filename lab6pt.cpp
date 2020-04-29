#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
double * xt = new double[200];
double * yt = new double[200];

double f(double x, int num) {
	double x1;
	double x2;
	double x3;
	if (num == 0) {
		x1 = x;
		x2 = 7;
	}
	if (num == 1) {
		x1 = 6;
		x2 = x;
	}
	return pow(x1, 4) + 20 * pow(x1, 3) + 2 * pow(x1, 2) * pow(x2, 2) + 36 * pow(x1, 2) * x2 +
		312 * pow(x1, 2) + 20 * x1 * pow(x2, 2) + 360 * x1 * x2 + 2121 * x1 + pow(x2, 4) + 36 * pow(x2, 3) + 537 * pow(x2, 2) + 3834 * x2 + 11308.;
}

double g(double * x) {
	return pow(x[0], 4) + 20 * pow(x[0], 3) + 2 * pow(x[0], 2) * pow(x[1], 2) + 36 * pow(x[0], 2) * x[1] +
		312 * pow(x[0], 2) + 20 * x[0] * pow(x[1], 2) + 360 * x[0] * x[1] + 2121 * x[0] + pow(x[1], 4) + 36 * pow(x[1], 3) + 537 * pow(x[1], 2) + 3834 * x[1] + 11308.;
}

double paul(double x1, double Dx, int n);

int main() {
	ofstream fout;
	fout.open("pt.dat");
	double * x = new double[3];
	x[0] = 6;
	x[1] = 7;
	fout << x[0] << " " << x[1] << endl;
	double alpha = 0.62;
	double beta = 10000;
	for (int i = 0; i < 2; i++) {
		if (i == 0)  x[i] = (paul(x[i], 0.01, i) - alpha);
		else x[i] = paul(x[i], 0.01, i);
	}
	cout << "Method Pauella: x1min = " << x[0] << " x2min =  " << x[1] << " fmin = " << g(x) << endl << "  Iterations " << xt[199] << endl;
	for (int i = 0; i < 16; i++) {
		if (i == 0 && fabs(xt[0])<beta && fabs(yt[0])<beta) {
			fout << xt[0]<< " " << yt[0] << endl;
		}
		if (i < 16 && i >= 1 && fabs(xt[i]) < beta && fabs(yt[i]) < beta) {
			fout << xt[i]  << " " << yt[i] << endl;
		}
	}
	fout.close();
	system("pt.py");
	return 0;
}

double paul(double x1, double Dx, int n) {
	double * pt = new double[200];
	pt[199] = 0;
	int iter = 0;
	double e1 = 0.1;
	double e2 = 0.1;
	double x, x2, x3, xmin;
	double fx1, fx2, fx3, fmin;
	double num, denum;
	x2 = x1 + Dx;
	fx1 = f(x1, n);
	fx2 = f(x2, n);
	if (fx1 > fx2) {
		x3 = x1 + 2 * Dx;
	}
	else {
		x3 = x1 - Dx;
	}
	//ШАГ 4
	for (;; ) {
		fx3 = f(x3, n);
		if (fx1 < fx2 && fx1 < fx3) {
			fmin = fx1;
			xmin = x1 + 0.62;
		}
		if (fx2 < fx1 && fx2 < fx3) {
			fmin = fx2;
			xmin = x2;
		}
		if (fx3 < fx1 && fx3 < fx2) {
			fmin = fx3;
			xmin = x3;
		}
		//ШАГ 5
		num = (pow(x2, 2) - pow(x3, 2))*fx1 + (pow(x3, 2) - pow(x1, 2))*fx2 + (pow(x1, 2) - pow(x2, 2))*fx3;
		denum = (x2 - x3)*fx1 + (x3 - x1)*fx2 + (x1 - x2)*fx3;
		x = 0.5*(num / denum);
		//ШАГ 6
		pt[iter] = x;
		if (n == 0) {
			xt = pt;
		}
		if (n == 1) {
			yt = pt;
		}
		if (pt[199] < iter) {
			pt[199] = iter;
		}
		if (fabs(xmin - x) < e1 && fabs(fmin - f(x, n)) < e2) {
			return x;
		}
		else {
			iter++;
			if (x < xmin) {
				xmin = x;
				x1 = xmin - 0.62;
				x2 = x1 + 0.62 + Dx;
				fx1 = f(x1, n);
				fx2 = f(x2, n);
				if (fx1 > fx2) {
					x3 = x1 + 2 * Dx;
				}
				else {
					x3 = x1 - Dx;
				}
			}
			else {
				x1 = xmin - 0.62;
				x2 = x1 + 0.62 + Dx;
				fx1 = f(x1, n);
				fx2 = f(x2, n);
				if (fx1 > fx2) {
					x3 = x1 + 2 * Dx;
				}
				else {
					x3 = x1 - Dx;
				}
			}
		}
	}
	return x;
}


