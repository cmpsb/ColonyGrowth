#include <iostream>
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
