// romano@romano-HP-ZBook-15 ~ $ ssh rvangenderen@hpc05.tudelft.net
// Copy to cluster  scp helloWorld.cpp hpc05.tudelft.net:mydir/mydir


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <boost/random.hpp>
#include <ctime>
#include <vector>
#include <string>
#include "print.h"
#include <ctime> // time_t
#include <cstdio>
#include <array>

#include <fstream>

const double pi = 4 * atan(1);

	//Characteristics
const double diameter = 1.0; //diameter should equal 0.75 +- 0.0375 micrometer in real life
const double dzeta = 1.0; //viscous drag coefficient
const double dt = 1.0; //timestep
const double tau = diameter * dzeta / dt;

	//Parameters
const double maxLoverD = 5.0; //gives relation between maxL and D
const double ki = 0.25 * tau; //internal spring constant, must be smaller dan dzeta*D*dt/2 otherwise system will explode, good value is 0.25
const double ko = 0.1 * tau; //overlap spring constant
const int relaxTime = 5; //set time to let system relax after growth step
const double growthRate = 0.00005 * tau * relaxTime; //number gives growth rate per timestep. exp: 1.23 per hour, small compared to ki for relaxation
const double growthRateDev = 0.05 * growthRate; //sets deviation in growth rate
const double maxLengthDev = 0.02 * diameter * maxLoverD; //sets deviation in max length
const double orientNoise = 0.01; //sets value for noise in orientation of daughter cells to prevent growing in one line
const int npivot = 3; //number of pivots

	//Constants
const double maxLength = diameter * maxLoverD; //sets mean maximum length of particles
const double startLength = maxLength / 2; //starting length of first particle
const int Nmax = 200; //maximum amount of particles

	//Random number generator
boost::mt19937 generator(time(0)); //number generator from random list
boost::normal_distribution<> //setup distributions
normalDistGrowth(0.0, growthRateDev), //growth distribution added to daughter cells, sigma = 0.277
normalDistAngle(0.0, orientNoise), //Noise in orientation for daughter cells
normalDistLength(maxLength, maxLengthDev); //distribution of maximum length that particles can reach, 4.54, sigma = 0.46
boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
randomMu(generator, normalDistGrowth), //generates deviation in growth rate
randomTheta(generator, normalDistAngle), //generates noise in orientation
randomLmax(generator, normalDistLength); //generates division length

struct Coordinate{
  double x;
  double y;
};

struct Particle{
  /* Constructor of the Particle structure*/
  Particle(double xstart, double ystart, double angle, double length, double diameter, double growth, double numberofpivots){

    ID = 0;
    D = diameter;
    L = length;
    mu = growth;
    Lmax = randomLmax();

    theta = angle;
    positions[0].x = xstart;
    positions[0].y = ystart;
    positions[numberofpivots + 1].x = xstart + L*cos(theta);
    positions[numberofpivots + 1].y = ystart + L*sin(theta);

    for (int i = 1; i <= numberofpivots; i++) {
      positions[i].x = xstart + (positions[numberofpivots + 1].x - xstart)/(numberofpivots + 1)*i;
      positions[i].y = ystart + (positions[numberofpivots + 1].y - ystart)/(numberofpivots + 1)*i;
    }

    for (int i = 0; i <= numberofpivots + 1; i++){
      velocities[i].x = 0;
      velocities[i].y = 0;
      forces[i].x = 0;
      forces[i].y = 0;
    }
  }

  std::vector<int> neighbours; //List of ID's of neighbouring particles
  int ID; //ID number
  double D, mu, Lmax; //diameter, growth rate, maximum length of a spring
  double theta; //angular noise in initial orientation
  double L; //rest length of springs

  std::array<Coordinate, npivot + 2> positions;
  std::array<Coordinate, npivot + 2> velocities;
  std::array<Coordinate, npivot + 2> forces;

  int colour = 240;

	void grow(){
    for (int i = 0; i <= npivot; i++){
      L += mu*dt;
    }
  }

  void str(){
    std::cout << "The ID is: " << ID << '\n';
    std::cout << "The diameter is: " << D << '\n';
    std::cout << "The growth rate is: " << mu << '\n';
    std::cout << "The maximum length before dividing is: " << mu << '\n';
    std::cout << "Spring has rest length: " << L << '\n';
    for(int i = 0; i <= npivot+1; i++ ){
      std::cout << "The coordinates of pivot " << i << " is: (" << positions[i].x << "," << positions[i].y << ")" << '\n';
    }
  }
	// void forceInternal();
	// void move();
	// void clear();
};

double dist(Coordinate s1, Coordinate s2){
  return sqrt((s1.y - s2.y)*(s1.y - s2.y) + (s1.x - s2.x)*(s1.x - s2.x));
}

void divide(Particle &pOld, Particle &pNew){
  //divide new points such that all L's are equal
  int split = npivot/2 + 1;
  pNew.theta += randomTheta();
  // Make new particle from the left half of the particle positions
  pNew.positions[0] = pOld.positions[0];
  // Correct for distance D/2
  double d = dist(pOld.positions[split], pOld.positions[split-1]);
  double phi = atan2((pOld.positions[split].y - pOld.positions[split-1].y), (pOld.positions[split].x - pOld.positions[split-1].x));
  pNew.positions[npivot+1].x = pOld.positions[split-1].x + (d - pOld.D/2)*cos(phi);
  pNew.positions[npivot+1].y = pOld.positions[split-1].y + (d - pOld.D/2)*sin(phi);
  //make space for an interpolated particle
  for(int i = 0; i < (npivot-1)/2; i++){
    pNew.positions[2*i+2] = pOld.positions[i+1];
  }
  //interpolation
  for(int i = 1; i <= npivot; i+=2){
    pNew.positions[i].x = (pNew.positions[i-1].x + pNew.positions[i+1].x)/2;
    pNew.positions[i].y = (pNew.positions[i-1].y + pNew.positions[i+1].y)/2;
  }
  // ---Update old particle from the right half of the particle positions---
  d = dist(pOld.positions[split+1], pOld.positions[split]);
  phi = atan2((pOld.positions[split].y - pOld.positions[split+1].y), (pOld.positions[split].x - pOld.positions[split+1].x));
  std::cout << "D " << d << '\n';
  pOld.positions[0].x = pOld.positions[split+1].x + (d - pOld.D/2)*cos(phi);
  pOld.positions[0].y = pOld.positions[split+1].y + (d - pOld.D/2)*sin(phi);
  //make space for an interpolated particle
  for(int i = 0; i < (npivot-1)/2; i++){
    pOld.positions[2*i+2] = pOld.positions[i+1+split];
  }
  //interpolation
  for(int i = 1; i <= npivot; i+=2){
    pOld.positions[i].x = (pOld.positions[i-1].x + pOld.positions[i+1].x)/2;
    pOld.positions[i].y = (pOld.positions[i-1].y + pOld.positions[i+1].y)/2;
  }

}


int main() {
  std::vector<Particle> p;
  p.reserve((int)(1.1*Nmax));
  Particle test = Particle(0, 0, 3.14159265, 8, 1, 5, npivot);
  test.positions[1].x = -2;
  test.positions[1].y = -2.5;
  test.positions[2].x = -4;
  test.positions[2].y = -5;
  test.positions[3].x = -6;
  test.positions[3].y = -2.5;
  p.push_back(test);
  p.push_back(Particle(0, 0, 3.14, 5, 1, 5, npivot));
  p.push_back(test);
  divide(p[0], p[1]);
  p[0].str();
  //p[1].str();
  // -- WRITING TEMP -- //
  ofstream outStream;
  outStream.open("/home/romano/PycharmProjects/BEP/Results.txt");
  for(int k = 0; k < 3; ++k)
  {
    outStream << p[k].D << std::endl;
    for(int i = 0; i < npivot+2; ++i){
      outStream << p[k].positions[i].x << "," << p[k].positions[i].y << std::endl;
    }
  }
  outStream.close();
  return 0;
}

//
// Particle::Particle()
//
// class Master{
//   Master(Particle)
//   void append();
//   void growAll();
//   void moveAll();
// };
