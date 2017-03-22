#include "twovec.h"
#include "coordinate.h"
#include <cmath>
#include <iostream>

//Empty constructor
Coordinate::Coordinate(){
    x = 9001;
    y = 0;
}

//Non-empty constructor
Coordinate::Coordinate(double xin, double yin){
    x = xin;
    y = yin;
}

//Get difference vector between two points
TwoVec Coordinate::operator-(Coordinate mycoord){
    return TwoVec(x - mycoord.x, y - mycoord.y);
}

//Shift point by a TwoVec
Coordinate Coordinate::operator+(TwoVec v){
    return Coordinate(x + v.x, y + v.y);
}

Coordinate Coordinate::operator-(TwoVec v){
    return Coordinate(x + v.x, y + v.y);
}

void Coordinate::str(){
    std::cout << "<" << x << "," << y << ">" << std::endl;
}

///Functions taking Coordinate
//Returns Euclidean distance between 2 coordinates
double dist(Coordinate s1, Coordinate s2){
    return sqrt((s1.y - s2.y)*(s1.y - s2.y) + (s1.x - s2.x)*(s1.x - s2.x));
}

//Returns angle between 2 coordinates
double ang(Coordinate s1, Coordinate s2){
    return atan2(s1.y - s2.y, s1.x - s2.x);
}
