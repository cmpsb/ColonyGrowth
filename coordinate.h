#ifndef COORDINATE_H_INCLUDED
#define COORDINATE_H_INCLUDED

#include <vector>

struct Coordinate{
    Coordinate();
    Coordinate(double xin, double yin);
    double x;
    double y;
    Coordinate operator+(TwoVec v);
    Coordinate operator-(TwoVec v);
    TwoVec operator-(Coordinate mycoord);
    void str();
};

double dist(Coordinate s1, Coordinate s2);
double ang(Coordinate s1, Coordinate s2);
double internalAngle(Coordinate s0, Coordinate s1, Coordinate s2);
double getTotalLength(std::vector<Coordinate> &myPoints);
Coordinate rotateAroundPoint(Coordinate p, double cx, double cy, double angle);

#endif // COORDINATE_H_INCLUDED
