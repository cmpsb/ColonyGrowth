#include "global.h"
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

//Finds the internal angle between 3 points
double internalAngle(Coordinate a, Coordinate b, Coordinate c){
    double alpha = atan2(b.y - a.y, b.x - a.x);
    double beta = atan2(b.y - c.y, b.x - c.x);
    if(alpha > beta){
        double theta = alpha - beta;
        if(theta > pi){
            return 2*pi - theta;
        }
        return theta;
    }
    double theta = beta - alpha;
    if(theta > pi){
        return 2*pi - theta;
    }
    return theta;

}

//Rotate point around center counterclockwise by angle radians
Coordinate spinAround(Coordinate point, Coordinate center, double angle){
    double cx = center.x;
    double cy = center.y;

    double s = sin(angle);
    double c = cos(angle);

    // translate point back to origin:
    point.x -= cx;
    point.y -= cy;

    // rotate point
    float xnew = point.x * c - point.y * s;
    float ynew = point.x * s + point.y * c;

    // translate point back:
    point.x = xnew + cx;
    point.y = ynew + cy;
  return point;
}
