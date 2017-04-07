#include <vector>
#include <cmath>
#include <array>
#include <iostream>
#include <boost/random.hpp>
#include "twovec.h"
#include "coordinate.h"
#include "global.h"
#include "shortestdistance.h"
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
        torques[i] = TwoVec(0,0);
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

///Prevents the particle from crossing through itself
void Particle::removeSelfOverlap(){
    double d;
    double Ftot;
    TwoVec Fp;
    TwoVec Fq;
    std::array<double, 5> st;
    for(int i = 0; i < npivot + 1; ++i){ //Loop over segments of p1
        for (int j = i+2; j < npivot + 1; ++j){ //Loop over segments of other particles
            st = dist3D_Segment_to_Segment(positions[i], positions[i+1], positions[j], positions[j+1]);
            d = st[0];
            if (d < D){ /* Is D a right value to use, and should it not depend on such things as maximum length and number of pivots */
                Ftot = ko*(d - D);
                Fp = TwoVec(st[1], st[2])*Ftot; //Force of a segment of particle p1 on p2
                Fq = Fp*-1;
                forces[i] += Fq*(1 - st[3]);
                forces[i + 1] += Fq*st[3];
                forces[j] += Fp*(1 - st[4]);
                forces[j + 1] += Fp*st[4];
            }
        }
    }
}

///Torsion spring
void Particle::torsionForce(){
    double gamma;
    TwoVec myTorque;
    TwoVec F;
    double newAngle;
    for(int i = 1; i < npivot + 1; ++i){
        gamma = internalAngle(positions[i - 1], positions[i], positions[i + 1]);
        myTorque = TwoVec(0.0,0.0);
        myTorque.z = -kappa*(gamma - restAngle); //Torque is divided evenly over both ends of the bar
        //To the "left"
        F = cross(myTorque, (positions[i-1] - positions[i]));

        newAngle = internalAngle(positions[i-1] + F, positions[i], positions[i + 1]);
        //std::cout << i  << "{ " << newAngle << " , " << gamma << " }" << std::endl;
        if(fabs(restAngle - newAngle) < fabs(restAngle - gamma) && i-1){
            if(i-1 != npivot-1){
                forces[i-1] += F;
            }

        }
        else{
            if(i-1 != npivot-1){
                    forces[i-1] += F*(-1);
            }
        }
        //To the "right"
        F = cross(myTorque, positions[i+1] - positions[i]);
        newAngle = internalAngle(positions[i-1], positions[i], positions[i + 1] + F);
        if(fabs(restAngle - newAngle) < fabs(restAngle - gamma)){
            forces[i+1] += F;
        }
        else{
            forces[i+1] += F*(-1);
        }
    }

}

///Moves the particle
void Particle::move(){
    TwoVec vel;
    for(int i = 0; i < npivot + 2; i++){
        vel = forces[i]*(2/dzeta*D);
        positions[i] = positions[i] + vel*dt;
    }
}

///Clears the forces after moving
void Particle::clear(){
    TwoVec naught = TwoVec(0, 0);
    for(TwoVec &f : forces){
        f = naught;
    }
    for(TwoVec &myTau : torques){
        myTau = TwoVec(0,0);
    }
}
