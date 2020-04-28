#include <iostream>  
#include <cmath>
#include <fstream>

using namespace std;

double F(double x1, double x2)
{

	return pow(x1, 4) + 20 * pow(x1, 3) + 2 * pow(x1, 2) * pow(x2, 2) + 36 * pow(x1, 2) * x2 +
		312 * pow(x1, 2) + 20 * x1 * pow(x2, 2) + 360 * x1 * x2 + 2121 * x1 + pow(x2, 4) + 36 * pow(x2, 3) + 537 * pow(x2, 2) + 3834 * x2 + 11308.;
}


double grad1(double x1, double x2)
{
	return 4 * pow(x1, 3) + 60 * pow(x1, 2) + 4 * x1 * pow(x2, 2) + 72 * x1 * x2 + 624 * x1 + 20 * pow(x2, 2) + 360 * x2 + 2121.;
}
double grad2(double x1, double x2)
{
	return 4 * pow(x1, 2) * x2 + 36 * pow(x1, 2) + 40 * x1 * x2 + 360 * x1 + 4 * pow(x2, 3) + 108 * pow(x2, 2) + 1074 * x2 + 3834.;
}



class Sopr
{
private:
	double  h = 0.001, xmin;
	int k = 0;
	double x0[2]{ 6.,7. };
	double x1[2]{ 0,0 };
	double p0[2]{ 0,0 }, p1[2]{ 0,0 };
	ofstream fout;

	double f(double alfa) {
		return F(x0[0] + alfa * p0[0], x0[1] + alfa * p0[1]);
	}

	double interpMnogochlen(double& x1, double& x2, double& x3)
	{
		if (((x2 - x3) * f(x1) + (x3 - x1) * f(x2) + (x1 - x2) * f(x3)) == 0) {//знаменатель обращается в нуль
			x1 = xmin;
			x2 = x1 + h;
			if (f(x1) > f(x2)) x3 = x1 + 2 * h;
			else if (f(x1) <= f(x2)) x3 = x1 - h;
			xmin = argmin(x1, x2, x3);
		}
		return (1. / 2.) * ((x2 * x2 - x3 * x3) * f(x1) + (x3 * x3 - x1 * x1) * f(x2) + (x1 * x1 - x2 * x2) * f(x3)) / ((x2 - x3) * f(x1) + (x3 - x1) * f(x2) + (x1 - x2) * f(x3));
	}

	double argmin(double x1, double x2, double x3)
	{
		if (f(x1) <= f(x2) && f(x1) <= f(x3)) return x1;
		else if (f(x2) <= f(x1) && f(x2) <= f(x3)) return x2;
		else { return x3; }
	}

	double argmin(double x1, double x2)
	{
		if (f(x1) <= f(x2)) return x1;
		else { return x2; }
	}

public:
	int GetIter()
	{
		return k;
	}

	double* soprGrad()
	{
		double  beta, alfa, a = 0, b = 0;
		fout.open("tr.txt");
		fout << x0[0] << " " << x0[1] << endl;
		double e = 0.001;
		int n = 2;
		double l = sqrt(pow(grad1(x0[0], x0[1]), 2) + pow(grad2(x0[0], x0[1]), 2));
		do {

			if (l <= e)
			{
				return x1;
			}
			else {
				if (k == 0) {
					p0[0] = -grad1(x0[0], x0[1]);
					p0[1] = -grad2(x0[0], x0[1]);
					p0[0] = p0[0] / sqrt(pow(p0[0], 2) + pow(p0[1], 2));
					p0[1] = p0[1] / sqrt(pow(p0[0], 2) + pow(p0[1], 2));
				}
				else
				{
					p1[0] = -grad1(x1[0], x1[1]) + beta * p0[0];
					p1[1] = -grad2(x1[0], x1[1]) + beta * p0[1];

					p1[0] = p1[0] / sqrt(pow(p1[0], 2) + pow(p1[1], 2));
					p1[1] = p1[1] / sqrt(pow(p1[0], 2) + pow(p1[1], 2));



					x0[0] = x1[0]; x0[1] = x1[1];
					x1[0] = 0; x1[1] = 0;
					p0[0] = p1[0]; p0[1] = p1[1];
					p1[0] = 0; p1[1] = 0;
				}

				if (k % n != 0)
				{

					beta = (grad1(x1[0], x1[1]) * (grad1(x1[0], x1[1]) - grad1(x0[0], x0[1])) + grad2(x1[0], x1[1]) * (grad1(x1[0], x1[1]) - grad1(x0[0], x0[1]))) /
						(grad1(x1[0], x1[1]) * grad1(x1[0], x1[1]) + grad2(x1[0], x1[1]) * grad2(x1[0], x1[1]));

				}
				else { beta = 0; }

				unimodal(7., a, b);
				alfa = KvadrInterp(a, b);
				x1[0] = x0[0] + alfa * p0[0];
				x1[1] = x0[1] + alfa * p0[1];
				fout << x1[0] << " " << x1[1] << endl;

				k = k + 1;

				l = sqrt(pow(grad1(x1[0], x1[1]), 2) + pow(grad2(x1[0], x1[1]), 2));
			}
		} while (true);
		fout.close();
	}


	double KvadrInterp(double a, double b)
	{
		double  e = 0.0001;
		double x1 = (a + b) / 2, x2, x3, x_;
		x2 = x1 + h;
		if (f(x1) > f(x2))
		{
			x3 = x1 + 2 * h;
		}
		else if (f(x1) <= f(x2))
		{
			x3 = x1 - h;
		}

		do {
			xmin = argmin(x1, x2, x3);
			x_ = interpMnogochlen(x1, x2, x3);

			if (abs(f(xmin) - f(x_)) < e && abs(xmin - x_) < e) return x_;
			else {
				if (x_ >= x1 && x_ <= x3) {
					x1 = argmin(xmin, x_);
					x2 = x1 + h;
					if (f(x1) > f(x2)) x3 = x1 + 2 * h;
					else if (f(x1) <= f(x2)) x3 = x1 - h;
				}
				else {
					x1 = x_;
					x2 = x1 + h;
					if (f(x1) > f(x2)) x3 = x1 + 2 * h;
					else if (f(x1) <= f(x2)) x3 = x1 - h;
				}
			}
		} while (true);
	}

	void unimodal(double x0, double& a, double& b)
	{
		bool r;
		do {
			r = f(x0) < f(x0 + h) && f(x0) < f(x0 + h);
			if (r) {
				a = x0 - h;
				b = x0 + h;
			}
			else if (f(x0 - h) > f(x0 + h)) {
				x0 = x0 + h / 2;
			}
			else {
				x0 = x0 - h / 2;
			}
		} while (!r);
	}
};


int main()
{
	Sopr s;
	double* arr = s.soprGrad();
	cout << "Minimum: " << arr[0] << " " << arr[1] << endl;
	cout << "Iteracii: " << s.GetIter() << endl;
	system("linur.py");

}
