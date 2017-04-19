#include "twovec.h" //Coordinate needs to know what a TwoVec is
#include "coordinate.h" //Particle needs to know what a Coordinate is
#include "global.h"
#include <iostream>
#include <list>
#include <cmath>
#include <vector>
#include <array>
#include <deque>
#include <tuple>
#include <boost/random.hpp>
#include "particle.h"
#include "division.h"

	//Random number generator
boost::mt19937 generator(time(0)); //number generator from random list
boost::normal_distribution<> //setup distributions
normalDistGrowth(0.0, growthRateDev), //growth distribution added to daughter cells, sigma = 0.277
normalDistAngle(0.0, 0), //Noise in orientation for daughter cells
normalDistLength(maxLength, maxLengthDev); //distribution of maximum length that particles can reach, 4.54, sigma = 0.46
boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
randomMu(generator, normalDistGrowth), //generates deviation in growth rate
randomTheta(generator, normalDistAngle), //generates noise in orientation
randomLmax(generator, normalDistLength); //generates division lengths

double getTotalLength(std::vector<Coordinate> &myPoints){
    double l = 0;
    for(int i = 0; i < myPoints.size()-1; ++i){
        l+=dist(myPoints[i], myPoints[i + 1]);
    }
    l /= (npivot+1);
	return l;
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

void correctHead(std::vector<Coordinate> &pos, double D){
    std::vector<Coordinate> myVec;
    for(Coordinate p : pos){
        myVec.push_back(p);
    }
    Coordinate endPoint = myVec.back();
    int s = 1;
    while(s != 0){
        s = myVec.size();
        std::pair<Coordinate, bool> ans = equidistantPointOnLine(myVec[s - 1], myVec[s - 2], endPoint, D/2);
        if(ans.second){
            myVec.pop_back();
            myVec.push_back(ans.first);
            break;
        }
        else myVec.pop_back();
    }
    pos.clear();
    for(Coordinate foo : myVec) pos.push_back(foo);
}

///Makes it such that all points besides the last one are equidistant and therefore relaxed
void relax(std::vector<Coordinate> &myArray){
    std::deque<Coordinate> myDeq;
    for(Coordinate i : myArray){
        myDeq.push_back(i); // Fill the double-ended queue (sort of like a two-headed vector) up with all coordinates
    }
    std::vector<Coordinate> fixed; // Contains the relaxed points
    fixed.push_back(myDeq[0]);
    myDeq.pop_front(); // Move the zeroth element, that is relaxed by definition, to the fixed vector
	
    double l = getTotalLength(myArray); // Relaxation length
    bool flag = false;
    std::pair<Coordinate, bool> myData; // Return type of equidistantPointOnLine

    for(int i = 0; i < (npivot+1)/2+1; ++i){
        while(!flag){
            myData = equidistantPointOnLine(fixed[i], myDeq[0], fixed[i], l);
            flag = myData.second;
            if(!flag){
                myDeq.pop_front(); // Pop the line if the point is passed
            }
        }
        fixed.push_back(myData.first);
        flag = false;
    }
    if(fixed.size() != npivot+2) fixed.push_back(myArray.back()); // To prevent off-by-one errors
	myArray.clear();
	for(Coordinate foo : fixed) myArray.push_back(foo);
}

void divide(Particle &pOld, Particle &pNew){

    // Find middle point of particle
    int split = (npivot+1)/2;
    pNew.mu = pOld.mu; //Growth noise instead of length noise
    pNew.Lmax = pOld.Lmax;
    pNew.D = pOld.D;
	
    // Insert points of old particle for relaxation
    std::vector<Coordinate> newPositions;
    newPositions.push_back(pOld.positions[0]);
    for(int i = 0; i < (npivot-1)/2; i++){
        newPositions.push_back(pOld.positions[i+1]);
    }
    newPositions.push_back(pOld.positions[split]);
    correctHead(newPositions, pOld.D);
    // Remove tension from all springs but one
    relax(newPositions);

    std::array<Coordinate, npivot+2> newPositionArray;
    for(int i = 0; i < npivot+2; ++i){
        newPositionArray[i] = newPositions[i];
    }
    pNew.positions = newPositionArray;
    pNew.positions[0] = pNew.positions[0] - TwoVec(0, 0.001); //Angle noise!

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
    newPositions.push_back(pOld.positions[split]);
    correctHead(newPositions, pOld.D);

    // Remove tension from all springs but one
    relax(newPositions);
    for(int i = 0; i < npivot+2; ++i){
        newPositionArray[i] = newPositions[i];
    }
    pOld.positions = newPositionArray;
    pOld.positions[0] = pOld.positions[0] + TwoVec(0, 0.001);

    // Set rest length equal to total particle length divided by number of springs
    totalLength = 0;
    for(int i = 0; i < npivot + 1; i++){
        totalLength += dist(pOld.positions[i], pOld.positions[i+1]);
    }
    pOld.len = totalLength;
    pOld.L = totalLength/(npivot+1);

}
