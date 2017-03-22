#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED

const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;

	//Characteristics
const double diameter = 1.0; //diameter should equal 0.75 +- 0.0375 micrometer in real life
const double dzeta = 1.0; //viscous drag coefficient
const double dt = 1; //timestep
const double tau = diameter * dzeta / dt;

	//Parameters
const double maxLoverD = 2.0; //gives relation between maxL and D
const double ki = 0.25 * tau; //internal spring constant, must be smaller than zeta*D*dt/2 otherwise system will explode, good value is 0.25
const double ko = 0.1 * tau; //overlap spring constant
const int relaxTime = 5; //set time to let system relax after growth step
const double growthRate = 0.00005 * tau * relaxTime; //number gives growth rate per time step. exp: 1.23 per hour, small compared to ki for relaxation
const double growthRateDev = 0.05 * growthRate; //sets deviation in growth rate
const double maxLengthDev = 0.02 * diameter * maxLoverD; //sets deviation in max length
const double orientNoise = 0.01; //sets value for noise in orientation of daughter cells to prevent growing in one line
const int npivot = 3; //number of pivots

	//Constants
const double maxLength = diameter * maxLoverD; //sets mean maximum length of particles
const double startLength = maxLength / 2; //starting length of first particle
const int Nmax = 200; //maximum amount of particles

#endif // GLOBAL_H_INCLUDED
