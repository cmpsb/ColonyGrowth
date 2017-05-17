#ifndef COORDINATE_H_INCLUDED
#define COORDINATE_H_INCLUDED

#include <vector>
#include <cmath>

struct Coordinate{
    double x;
    double y;

    Coordinate() : x(9001), y(0) {}
    Coordinate(double xin, double yin) : x(xin), y(yin) {}

    Coordinate operator+(TwoVec v) {
        return Coordinate(x + v.x, y + v.y);
    }
    Coordinate operator-(TwoVec v) {
        return Coordinate(x + v.x, y + v.y);
    }

    TwoVec operator-(Coordinate mycoord) {
        return TwoVec(x - mycoord.x, y - mycoord.y);
    }


    void str();
};

static double dist(Coordinate s1, Coordinate s2) {
    return sqrt((s1.y - s2.y)*(s1.y - s2.y) + (s1.x - s2.x)*(s1.x - s2.x));
}
static double ang(Coordinate s1, Coordinate s2) {
    return atan2(s1.y - s2.y, s1.x - s2.x);
}
double internalAngle(Coordinate s0, Coordinate s1, Coordinate s2);
double getTotalLength(std::vector<Coordinate> &myPoints);
Coordinate rotateAroundPoint(Coordinate p, double cx, double cy, double angle);

#endif // COORDINATE_H_INCLUDED
