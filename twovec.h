#ifndef TWOVEC_H_INCLUDED
#define TWOVEC_H_INCLUDED

#include <cmath>

struct TwoVec{
    TwoVec() : x(0), y(0), z(0) {}
    TwoVec(double xin, double yin) : x(xin), y(yin), z(0) {}

    double x;
    double y;
    double z;

    TwoVec operator+(TwoVec v) {
        return TwoVec(x + v.x, y + v.y);
    }

    TwoVec& operator+=(const TwoVec& v) {
        this->x += v.x;
        this->y += v.y;
        this->z += v.z;
        return *this;
    }

    TwoVec operator-(TwoVec v) {
        return TwoVec(x - v.x, y - v.y);
    }

    TwoVec& operator-=(const TwoVec& v) {
        this->x -= v.x;
        this->y -= v.y;
        this->z -= v.z;
        return *this;
    }

    TwoVec operator*(double d) {
        return TwoVec(d*x, d*y);
    }

    TwoVec operator/(double d) {
        return TwoVec(x/d, y/d);
    }

    double operator*(TwoVec myvec) {
        return (x*myvec.x + y*myvec.y);
    }

    void str();
};

static TwoVec operator*(double d, TwoVec v) {
    return TwoVec(d*v.x, d*v.y);
}

static TwoVec cross(TwoVec u, TwoVec v) {
    TwoVec ans;
    ans.x = (u.y * v.z) - (u.z * v.y);
    ans.y = (u.z * v.x) - (u.x * v.z);
    ans.z = (u.x * v.y) - (u.y * v.x);
    return ans;
}
static double determinant(TwoVec a, TwoVec b, TwoVec c) {
    double x = b.y * c.z - c.y * b.z;
    double y = a.y * c.z - c.y * a.z;
    double z = a.y * b.z - b.y * a.z;
    return a.x*x - b.x*y + c.x*z;
}
static double norm(TwoVec a) {
    return sqrt(a*a);
}
double angleBetweenVectors(TwoVec u, TwoVec v);

double mySin(double input);

#endif // TWOVEC_H_INCLUDED
