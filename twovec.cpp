#include <iostream>
#include <cmath>
#include "twovec.h"

///Empty constructor
TwoVec::TwoVec(){
    x = 0;
    y = 0;
    z = 0;
}

///Non-empty constructor
TwoVec::TwoVec(double xin, double yin){
    x = xin;
    y = yin;
    z = 0;
}

///Vector addition
TwoVec TwoVec::operator+(TwoVec v){
    return TwoVec(x + v.x, y + v.y);
}
TwoVec& TwoVec::operator+=(const TwoVec& v){
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
}

///Vector subtraction
TwoVec TwoVec::operator-(TwoVec v){
    return TwoVec(x - v.x, y - v.y);
}
TwoVec& TwoVec::operator-=(const TwoVec& v){
    this->x -= v.x;
    this->y -= v.y;
    this->z -= v.z;
    return *this;
}

///Scalar multiplication
TwoVec TwoVec::operator*(double d){
    return TwoVec(d*x, d*y);
}
TwoVec operator*(double d, TwoVec v){
    return TwoVec(d*v.x, d*v.y);
}
TwoVec TwoVec::operator/(double d){
    return TwoVec(x/d, y/d);
}

///Dot product
double TwoVec::operator*(TwoVec myvec){
    return (x*myvec.x + y*myvec.y);
}

///Print method
void TwoVec::str(){
    std::cout << "<" << x << "," << y << ">" << std::endl;
}

//Functions involving TwoVecs
///Cross product
TwoVec cross(TwoVec u, TwoVec v){
    TwoVec ans;
    ans.x = (u.y * v.z) - (u.z * v.y);
    ans.y = (u.z * v.x) - (u.x * v.z);
    ans.z = (u.x * v.y) - (u.y * v.x);
    return ans;
}

///3x3 determinant
double determinant(TwoVec a, TwoVec b, TwoVec c){
    double x = b.y * c.z - c.y * b.z;
    double y = a.y * c.z - c.y * a.z;
    double z = a.y * b.z - b.y * a.z;
    return a.x*x - b.x*y + c.x*z;
}

///Returns length of a TwoVec
double norm(TwoVec a){
    return sqrt(a*a);
}

double angleBetweenVectors(TwoVec u, TwoVec v){
    double alpha = atan2(u.y, u.x);
    double beta = atan2(v.y, v.x);
    if(alpha > beta){
        return alpha - beta;
    }
    return beta - alpha;
}

//Not related to anything, but often needed so placed high in the dependency tree
///Sine that cuts off
double mySin(double input){
    if(fabs(input) < 0.0001 || fabs(input-3.1415) < 0.0001){
        return 0;
    }
    return sin(input);
}
