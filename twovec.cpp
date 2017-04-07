#include <iostream>
#include <cmath>
#include "twovec.h"

//Empty constructor
TwoVec::TwoVec(){
    x = 0;
    y = 0;
    z = 0;
}

//Non-empty constructor
TwoVec::TwoVec(double xin, double yin){
    x = xin;
    y = yin;
    z = 0;
}

//Vector addition
TwoVec TwoVec::operator+(TwoVec v){
    return TwoVec(x + v.x, y + v.y);
}
TwoVec& TwoVec::operator+=(const TwoVec& v){
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
}

//Vector subtraction
TwoVec TwoVec::operator-(TwoVec v){
    return TwoVec(x - v.x, y - v.y);
}

//Scalar multiplication
TwoVec TwoVec::operator*(double d){
    return TwoVec(d*x, d*y);
}

//Dot product
double TwoVec::operator*(TwoVec myvec){
    return (x*myvec.x + y*myvec.y);
}

void TwoVec::str(){
    std::cout << "<" << x << "," << y << ">" << std::endl;
}

TwoVec cross(TwoVec u, TwoVec v){
    TwoVec ans;
    ans.x = (u.y * v.z) - (u.z * v.y);
    ans.y = (u.z * v.x) - (u.x * v.z);
    ans.z = (u.x * v.y) - (u.y * v.x);
    return ans;
}

double mySin(double input){
    if(fabs(input) < 0.0001 || fabs(input-3.1415) < 0.0001){
        return 0;
    }
    return sin(input);
}
