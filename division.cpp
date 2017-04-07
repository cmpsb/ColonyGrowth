#include "twovec.h" //Coordinate needs to know what a TwoVec is
#include "coordinate.h" //Particle needs to know what a Coordinate is
#include "global.h"
#include <iostream>
#include <list>
#include <cmath>
#include <vector>
#include <array>
#include <deque>
#include "particle.h"
#include "division.h"

double getTotalLength(std::vector<Coordinate> &myPoints){
    double l = 0;
    for(int i = 0; i < myPoints.size()-1; ++i){
        l+=dist(myPoints[i], myPoints[i + 1]);
    }
    l /= (myPoints.size()+1);
}

///Finds Coordinate on a line that is a distance l removed from point circleCenter, returning bool to show if it succeeded
std::pair<Coordinate, bool> equidistantPointOnLine(Coordinate p1, Coordinate p2, Coordinate circleCenter, double radius){
    TwoVec d = p2 - p1;
    TwoVec f = p1 - circleCenter;
    double a = d*d;
    double b = f*d*2;
    double c = f*f - radius*radius;

    double D = b*b-4*a*c;
    if(D >= 0){
        D = sqrt(D);
        float t1 = (-b - D)/(2*a);
        float t2 = (-b + D)/(2*a);
        if(t1 >= 0 && t1 <= 1){
            return std::pair<Coordinate, bool>{p1 + (d*t1), true};
        }
        if(t2 >= 0 && t2 <= 1){
            return std::pair<Coordinate, bool>{p1 + (d*t2), true};
        }
    }
    return std::pair<Coordinate, bool>{Coordinate(0,0), false}; // If no equidistant point is found, a coordinate that is over nine thousand will be returned
}

///Makes it such that all points besides the last one are equidistant and therefore relaxed
void relax(std::vector<Coordinate> &myArray){
    std::deque<Coordinate> myDeq;
    for(Coordinate i : myArray){
        myDeq.push_back(i);
    }
    std::vector<Coordinate> fixed;
    fixed.push_back(myDeq[0]);
    myDeq.pop_front();

    double l = getTotalLength(myArray);
    bool flag = false;
    std::pair<Coordinate, bool> myData;

    for(int i = 0; i < (npivot+1)/2+1; ++i){
        while(!flag){
            myData = equidistantPointOnLine(fixed[i], myDeq[0], fixed[i], l);
            flag = myData.second;
            if(!flag){
                myDeq.pop_front();
            }
        }
        fixed.push_back(myData.first);
        flag = false;
    }
    fixed.push_back(myArray[(npivot+1)/2]);
    myArray = fixed;
}

void divide(Particle &pOld, Particle &pNew){

    // Find middle point of particle
    int split = (npivot+1)/2;
    /*
    Got to add angular noise to the model, but not sure if on particle or spring scale
    */
    pNew.mu = pOld.mu; //Growth noise instead of length noise
    pNew.Lmax = pOld.Lmax;
    pNew.D = pOld.D;

    // Insert points of old particle for relaxation
    std::vector<Coordinate> newPositions;
    newPositions.push_back(pOld.positions[0]);
    for(int i = 0; i < (npivot-1)/2; i++){
        newPositions.push_back(pOld.positions[i+1]);
    }
    double d = dist(pOld.positions[split], pOld.positions[split-1]);
    double phi = ang(pOld.positions[split], pOld.positions[split-1]);
    TwoVec e{cos(phi), mySin(phi)};
    Coordinate shifted = pOld.positions[split-1] + e*(d-pOld.D/2);
    newPositions.push_back(shifted);

    // Remove tension from all springs but one
    relax(newPositions);
    std::array<Coordinate, npivot+2> newPositionArray;
    for(int i = 0; i < npivot+2; ++i){
        newPositionArray[i] = newPositions[i];
    }
    pNew.positions = newPositionArray;

    // Set rest length equal to total particle length divided by number of springs
    double totalLength;
    for(int i = 0; i < npivot + 1; i++){
        totalLength += dist(pNew.positions[i], pNew.positions[i+1]);
    }
    pNew.len = totalLength;
    pNew.L = totalLength/(npivot+1);

    // ---Right particle pOld---

    // ---Update old particle from the right half of the particle positions---
    newPositions.clear();
    newPositions.push_back(pOld.positions[npivot+1]);
    for(int i = 0; i < (npivot-1)/2; i++){
        newPositions.push_back(pOld.positions[i+1+split]);
    }
    d = dist(pOld.positions[split], pOld.positions[split+1]);
    phi = ang(pOld.positions[split], pOld.positions[split+1]);
    e = TwoVec(cos(phi), mySin(phi));
    shifted = pOld.positions[split+1] + e*(d-pOld.D/2);
    newPositions.push_back(shifted);

    for(auto foo : newPositions) foo.str();

    // Remove tension from all springs but one
    relax(newPositions);
    for(int i = 0; i < npivot+2; ++i){
        newPositionArray[i] = newPositions[i];
    }
    pOld.positions = newPositionArray;

    // Set rest length equal to total particle length divided by number of springs
    totalLength = 0;
    for(int i = 0; i < npivot + 1; i++){
        totalLength += dist(pOld.positions[i], pOld.positions[i+1]);
    }
    pOld.len = totalLength;
    pOld.L = totalLength/(npivot+1);

}

