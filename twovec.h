#ifndef TWOVEC_H_INCLUDED
#define TWOVEC_H_INCLUDED

struct TwoVec{
    TwoVec();
    TwoVec(double xin, double yin);

    double x;
    double y;
    double z;

    TwoVec operator+(TwoVec myvec);
    TwoVec& operator+=(const TwoVec& v);
    TwoVec operator-(TwoVec myvec);
    TwoVec operator*(double d);
    double operator*(TwoVec myvec);

    void str();
};

TwoVec cross(TwoVec u, TwoVec v);

double mySin(double input);

#endif // TWOVEC_H_INCLUDED
