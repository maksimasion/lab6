#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>
#include <limits>
#include "vector.h"
using namespace std;
#include <cmath>

Vector::Vector() : x1(0), x2(0) {};

Vector::Vector(double a, double b) : x1(a), x2(b) {};
double Vector::norm() const {
	return hypot(x1, x2);
}

Vector Vector::operator+(const Vector& a) const {
	return Vector(x1 + a.x1, x2 + a.x2);
}

Vector Vector::operator*(double a) const {
	return Vector(a * x1, a * x2);
}

Vector Vector::operator-() const {
	return Vector{ -x1, -x2 };
}

Vector Vector::operator-(const Vector& a) const {
	return *this + (-a);
}

Vector operator*(double a, const Vector& v) {
	return v * a;
}

Vector Vector::operator/(double a) const {
	return *this * (1 / a);
}

std::ostream& operator<<(std::ostream& out, const Vector& v) {
	out << v.x1 << ' ' << v.x2;
	return out;
}


const double epsilon = 1e-4;

double f(Vector x) {
	return pow(x.x1, 4) + 20 * pow(x.x1, 3) + 2 * pow(x.x1, 2) * pow(x.x2, 2) + 36 * pow(x.x1, 2) * x.x2 +
		312 * pow(x.x1, 2) + 20 * x.x1 * pow(x.x2, 2) + 360 * x.x1 * x.x2 + 2121 * x.x1 + pow(x.x2, 4) + 36 * pow(x.x2, 3) + 537 * pow(x.x2, 2) + 3834 * x.x2 + 11308.;
}

int main() {
	// random number generator initialization
	random_device dev;
	mt19937 rng(dev());
	uniform_real_distribution<> dist(-1., 1.);

	const double alpha = 1.618;
	const double beta = 0.618;
	const int M = 6;
	const int N = 5000;
	double t = 1;
	int k = 0;
	int j = 1;
	Vector x{ 6, 7 };
	Vector x0 = x;

	ofstream fout;
	fout << x << '\n';
	fout.open("direction.dat");
	do {
		bool success = false;
		Vector d, z;
		d.x1 = dist(rng);
		d.x2 = dist(rng);
		Vector y = x + t * d / d.norm();
		if (f(y) < f(x)) {
			z = x + alpha * (y - x);
			if (f(z) < f(x)) {
				success = true;
				x = z;
				fout << x << '\n';
				t = alpha * t;
				k++;
			}
		}
		if (f(y) >= f(x) || !success) {
			if (j < M) {
				j++;
			}
			else {
				t = beta * t;
				j = 1;
			}
		}
	} while (k < N && t > epsilon);
	fout.close();

	cout <<"Min:"<< x << ' ' <<"Iter:"<< k << endl;

	vector<double> x1;
	vector<double> x2;
	fout.open("x1.dat");

	for (double x = 0; x < 500; x += 0.1) {
		x1.push_back(x);
		fout << x << '\n';
	}

	fout.close();
	fout.open("x2.dat");

	for (double x = 0; x < 70; x += 0.1) {
		x2.push_back(x);
		fout << x << '\n';
	}

	fout.close();
	fout.open("f.dat");

	for (size_t j = 0; j < x2.size(); j++) {
		for (size_t i = 0; i < x1.size(); i++) {
			fout << f(Vector(x1[i], x2[j])) << ' ';
		}
		fout << '\n';
	}
	fout.close();

	system("linurrs.py");
}
