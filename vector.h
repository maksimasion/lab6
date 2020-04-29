#pragma once
#include <iostream>

class Vector {
public:
    double x1;
    double x2;
    Vector();
    Vector(double a, double b);

    double norm() const;
    Vector operator+(const Vector& a) const;
    Vector operator*(double a) const;
    Vector operator/(double a) const;
    Vector operator-() const;
    Vector operator-(const Vector& a) const;
    double& operator[](std::size_t i);
    const double& operator[](std::size_t i) const;
    friend std::ostream& operator<<(std::ostream& out, const Vector& v);
};

Vector operator*(double a, const Vector& v);
