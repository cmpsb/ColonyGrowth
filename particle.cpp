#include <vector>
#include <array>
#include <iostream>
#include <boost/random.hpp>
#include "twovec.h"
#include "coordinate.h"
#include "global.h"
#include "particle.h"


/// Constructor of the Particle structure
Particle::Particle(double xstart, double ystart, double angle, double restlength, double diameter, double growth){

    ID = 0;
    D = diameter;
    L = restlength;
    mu = growth;
    Lmax = maxLength;
    len = L*(npivot+1);

    theta = angle;
    positions[0].x = xstart;
    positions[0].y = ystart;
    positions[npivot + 1].x = xstart + L*(npivot + 1)*cos(theta);
    positions[npivot + 1].y = ystart + L*(npivot + 1)*sin(theta);

    for (int i = 1; i <= npivot; i++) {
        positions[i].x = xstart + (positions[npivot + 1].x - xstart)/(npivot + 1)*i;
        positions[i].y = ystart + (positions[npivot + 1].y - ystart)/(npivot + 1)*i;
    }

    for (int i = 0; i <= npivot + 1; i++){
        forces[i].x = 0;
        forces[i].y = 0;
    }
}

///Grows every spring with distance mu
void Particle::grow(){
    L += mu*dt;
}

///Output a couple of interesting facts about a particle
void Particle::str(){
    std::cout << "The ID is: " << ID << '\n';
    std::cout << "The diameter is: " << D << '\n';
    std::cout << "The growth rate is: " << mu << '\n';
    std::cout << "The head-to-head length is: " << len << '\n';
    std::cout << "The rest length each spring needs to have before dividing is: " << Lmax << '\n';
    std::cout << "Each spring has rest length: " << L << '\n';
    for(int i = 0; i <= npivot+1; i++ ){
        std::cout << "The coordinates of pivot " << i << " is: (" << positions[i].x << "," << positions[i].y << ")" << '\n';
    }
    for(int i = 0; i <= npivot+1; i++ ){
    std::cout << "The force on pivot " << i << " is: (" << forces[i].x << "," << forces[i].y << ")" << '\n';
    }
}

///Calculates the internal spring force
void Particle::forceInternal(){
    double dx;
    double dy;
    double len;
    double fx;
    double fy;
    for(int i = 0; i < npivot+1; i++){
        dx = positions[i].x - positions[i+1].x;
        dy = positions[i].y - positions[i+1].y;
        len = dist(positions[i], positions[i+1]);
        fx = -ki * (len - L) * (dx / len);
        fy = -ki * (len - L) * (dy / len);
        forces[i] = TwoVec(forces[i].x + fx, forces[i].y + fy);
        forces[i+1] = TwoVec(forces[i+1].x - fx, forces[i+1].y - fy);
    }
}

///Moves the particle
void Particle::move(){
    double xvel;
    double yvel;
    for(int i = 0; i < npivot + 2; i++){
        xvel = 2*forces[i].x / (dzeta * D);
        yvel = 2*forces[i].y / (dzeta * D);
        positions[i].x = positions[i].x + xvel*dt;
        positions[i].y = positions[i].y + yvel*dt;
    }
}

///Clears the forces after moving
void Particle::clear(){
    TwoVec naught = TwoVec(0, 0);
    for(TwoVec &f : forces){
        f = naught;
    }
}
