#include <iostream>
#include <cmath>
#include "twovec.h"

///Print method
void TwoVec::str(){
    std::cout << "<" << x << "," << y << ">" << std::endl;
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
