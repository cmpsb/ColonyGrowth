/* TODO:
- Aren't you looping over neighbours twice?
- p array vs p vector speed
- Local vs global p
- why cstdio?
*/

/* SETUP:
- set parameters
- set max number of particles
- change folder name
- choose number of loops (number of data files)
- optional: set lag time to save, ...
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <boost/random.hpp>
#include <ctime>
#include <vector>
#include <string>
#include "print.h"
//#include <cstdio>

using namespace std;
const double pi = 4 * atan(1);

//Characteristics
const double diameter = 1.0; //diameter should equal 0.75 +- 0.0375 micrometer in real life
const double dzeta = 1.0; //viscous drag coefficient
const double dt = 1.0; //timestep
const double tau = diameter * dzeta / dt;

//Parameters
const double maxLoverD = 5; //gives relation between maxL and D, for rod-shaped particles: LD > 1
const double ki = 0.1 * tau; //internal spring constant, must be smaller dan dzeta*D*dt/2 otherwise system will explode, good value is 0.25
const double ko = 0.2 * tau; //overlap spring constant
const int relaxTime = 5; //set time to let system relax after growth step
const double growthRate = 0.00005 * (0.25 * (maxLoverD + diameter)) * tau * relaxTime; //number gives growth rate per timestep. experiments: 1.23 per hour, must be small compared to ki for relaxation
const double growthRateDev = 0.1 * growthRate; //sets deviation in growth rate
const double maxLengthDev = 0.1 * diameter * maxLoverD; //sets deviation in max length
const double orientNoise = 0.1; //sets value for noise in orientation of daughter cells to prevent growing in one line

								//Constants
const double maxLength = diameter * maxLoverD; //sets mean maximum length of particles
const double startLength = maxLength / 2; //starting length of first particle

const int Nmax = 2000; //maximum amount of particles

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

struct Particle
	//Structure of particles with its properties
{
	Particle(double, double, double, double, double, double); //constructor

	void grow();
	void forceInternal();
	void move();
	void clear();

	vector<int> neighbours;

	int ID;
	double D, L, len, mu, Lmax; //size properties: diameter, restlength, actual length, growth rate, max length
	double xcpos, ycpos, theta; //spatial variables
	double x1pos, y1pos, x2pos, y2pos, xvel1, yvel1, xvel2, yvel2, Fx1, Fy1, Fx2, Fy2; //more spatial variables
	double oop; int ram; //sets orientational order parameter and defect as property of particle
	int colour;
};

Particle::Particle(double xstart, double ystart, double angle, double length, double diameter, double growth)
{
	ID = 0;
	D = diameter;
	L = length; //rest length of particle (internal spring is at rest)
	len = L;
	mu = growth; //initial growth rate of particle
	Lmax = randomLmax();

	theta = angle;
	x1pos = xstart; y1pos = ystart;
	x2pos = x1pos + L * cos(theta); y2pos = y1pos + L * sin(theta);
	xcpos = (x2pos + x1pos) / 2; ycpos = (y2pos + y1pos) / 2;

	xvel1 = 0.0; yvel1 = 0.0; xvel2 = 0.0; yvel2 = 0.0; //velocities of both ends of particle
	Fx1 = 0.0; Fy1 = 0.0; Fx2 = 0.0; Fy2 = 0.0;

	oop = 1.0; ram = 0;
	colour = 180; //used for visualisation
}

void Particle::grow()
{
	//double dx = x1pos - x2pos; double dy = y1pos - y2pos;
	//double len = sqrt(dx*dx + dy*dy); //compute distance between discs to determine spring force

	L += (mu * dt);

	//x1pos -= ((mu * dt) / 2) * cos(theta); y1pos -= ((mu * dt) / 2) * sin(theta);
	//x2pos += ((mu * dt) / 2) * cos(theta); y2pos += ((mu * dt) / 2) * sin(theta);
}

void Particle::forceInternal()
{
	double f_x = 0.0; double f_y = 0.0;
	double dx = x1pos - x2pos; double dy = y1pos - y2pos;
	len = sqrt(dx*dx + dy*dy); //compute distance between discs to determine spring force

	f_x = -ki * (len - L) * (dx / len);
	f_y = -ki * (len - L) * (dy / len);

	Fx1 += f_x; Fy1 += f_y;
	Fx2 -= f_x; Fy2 -= f_y;
}

void Particle::move()
{
	xvel1 = 2 * Fx1 / (dzeta * D); yvel1 = 2 * Fy1 / (dzeta * D); //calculate velocities
	xvel2 = 2 * Fx2 / (dzeta * D); yvel2 = 2 * Fy2 / (dzeta * D);

	x1pos += xvel1 * dt; y1pos += yvel1 * dt; //update positions with velocities
	x2pos += xvel2 * dt; y2pos += yvel2 * dt;

	xcpos = (x2pos + x1pos) / 2; ycpos = (y2pos + y1pos) / 2; //update central positions

	theta = atan((y2pos - y1pos) / (x2pos - x1pos)); //update orientation
}

void Particle::clear()
{
	xvel1 = 0.0; yvel1 = 0.0; xvel2 = 0.0; yvel2 = 0.0; //velocities of both ends of particle
	Fx1 = 0.0; Fy1 = 0.0; Fx2 = 0.0; Fy2 = 0.0;
}

void divide(Particle &pOld, Particle &pNew)
{
	double oldL = pOld.L;
	double oldDx = pOld.x2pos - pOld.x1pos; double oldDy = pOld.y2pos - pOld.y1pos;

	pNew.L = (oldL - pNew.D) / 2; //update new positions daughter cell
	pNew.x2pos = pOld.x2pos; pNew.y2pos = pOld.y2pos; //needed to make sure positions are chosen in right direction
	pNew.x1pos = pNew.x2pos - (oldDx * pNew.L / oldL); pNew.y1pos = pNew.y2pos - (oldDy * pNew.L / oldL);

	double dtheta = randomTheta(); //give daughter cell new orientation (noise)
	pNew.x2pos += dtheta * sin(pNew.theta); pNew.y2pos += dtheta * cos(pNew.theta);

	pNew.mu += randomMu(); //change growthrate of dauhter cell

	pOld.L = (oldL - pOld.D) / 2; //update new positions mother cell
	pOld.x2pos = pOld.x1pos + (oldDx * pOld.L / oldL); pOld.y2pos = pOld.y1pos + (oldDy * pOld.L / oldL);
	pOld.mu += randomMu();
}

double inRange(double d0)
{
	double rp = d0;
	if (rp > 1.0) //check if relative point is in range [0,1]
		rp = 1.0; //outer right point of line segment, p.x2pos
	else if (rp < 0.0)
		rp = 0.0; //outer left point of line segment,  p.x1pos
	return rp;
}

void forceOverlap(Particle &p1, Particle &p2, double Rx, double Ry, double dist, double rp1, double rp2)
{
	double f_x = 0.0; double f_y = 0.0;
	double rlen = (p1.D + p2.D) / 2; //compute restlength (corresponds to D) of overlap spring

	f_x = -ko * (dist - rlen) * Rx;
	f_y = -ko * (dist - rlen) * Ry;

	p1.Fx1 += (1 - rp1) * f_x; p1.Fx2 += rp1 * f_x; //Force is distributed over ends of line segments (xy1pos, xy2pos)
	p1.Fy1 += (1 - rp1) * f_y; p1.Fy2 += rp1 * f_y;
	p2.Fx1 -= (1 - rp2) * f_x; p2.Fx2 -= rp2 * f_x; //Newtons second law f1 = -f2
	p2.Fy1 -= (1 - rp2) * f_y; p2.Fy2 -= rp2 * f_y;
}

void distForce(Particle &p1, Particle &p2)
{
	double d1 = 0.0, d2; //relative point [0,1] on line segments of particles to determine closest distance, d1 starts as 0 when lines are parallel

						 //calculate dot products, s# corresponds to line segment of particle #, r is vector between first point of s1, s2
	double s1s1 = p1.len*p1.len; double s2s2 = p2.len*p2.len;
	double s1s2 = (p1.x2pos - p1.x1pos) * (p2.x2pos - p2.x1pos) + (p1.y2pos - p1.y1pos) * (p2.y2pos - p2.y1pos);
	double s1r = (p1.x2pos - p1.x1pos) * (p1.x1pos - p2.x1pos) + (p1.y2pos - p1.y1pos) * (p1.y1pos - p2.y1pos);
	double s2r = (p1.x1pos - p2.x1pos) * (p2.x2pos - p2.x1pos) + (p1.y1pos - p2.y1pos) * (p2.y2pos - p2.y1pos);

	double deler = s1s1*s2s2 - s1s2*s1s2; //denominator to determine d1, d2. d1 is set 0 if lines are parallel
	if (deler != 0.0)
	{
		d1 = inRange((s1s2*s2r - s2s2*s1r) / deler); //d1 = (s1s2*s2r - s2s2*s1r) / deler, check if in range
	}

	d2 = (s1s2*d1 + s2r) / s2s2;
	if (d2 > 1.0) //recalculate d1 if d2 does not lie on second line segment (not in range [0,1]), fill in for d2.
	{
		d2 = 1.0; //outer right point of second line segment, p2.x2pos
		d1 = inRange((s1s2 - s1r) / s1s1); //d1 = (s1s2 - s1r) / s1s1, check if in range
	}
	else if (d2 < 0.0)
	{
		d2 = 0.0; //outer left point of second line segment,  p2.x1pos
		d1 = inRange(-s1r / s1s1); //d1 = -s1r / s1s1, check if in range
	}
	//determine points on line segments that give smallest distance between line segments
	double D1x = p1.x1pos + d1*(p1.x2pos - p1.x1pos); double D1y = p1.y1pos + d1*(p1.y2pos - p1.y1pos);
	double D2x = p2.x1pos + d2*(p2.x2pos - p2.x1pos); double D2y = p2.y1pos + d2*(p2.y2pos - p2.y1pos);
	double sqd = (D1x - D2x)*(D1x - D2x) + (D1y - D2y)*(D1y - D2y); //square of distance between segments

	double touch = (p1.D + p2.D) / 2; //distance for which particles can touch each other
	if (sqd < (touch * touch)) //check if the particles overlap
	{
		double Rx = 0.0;
		double Ry = 0.0;
		double dist = sqrt(sqd);
		if (dist != 0.0)
		{
			Rx = (D1x - D2x) / dist; //gives the relative difference in length for x
			Ry = (D1y - D2y) / dist;
		}
		forceOverlap(p1, p2, Rx, Ry, dist, d1, d2); //calculate forces between n and nb
	}

	if (sqd < 2 * p1.D*p2.D) //Check if particles are neighbours, if so, append to neighbour vector
	{
		p1.neighbours.push_back(p2.ID); //only this line is needed for ram, otherwise ram is performed 3 times more
		p2.neighbours.push_back(p1.ID); //use only for oop
	}
}

void growAll(vector<Particle> &p, int nop)
{
	for (int g = 0;g < nop;g++) //first loop to enable growth
	{
		p[g].neighbours.clear(); //clear neighbours vector
		p[g].grow();
		if (p[g].len > p[g].Lmax) //division, choose restlength (L) or actual Length (len)
		{
			p.push_back(Particle(p[g].x1pos, p[g].y1pos, p[g].theta, p[g].L, p[g].D, p[g].mu)); //new particle is made with same properties as mother
			divide(p[g], p[p.size() - 1]); //function that sets new properties to daughter particles
			p[p.size() - 1].ID = p.size() - 1; //set ID of new particle
		}
	}
}

void printAll(vector<Particle> &p, Print &myPrint, int nop)
{
	for (int w = 0; w < nop; w++)
	{
		myPrint.print_data_js(p[w].x1pos, p[w].y1pos, p[w].x2pos, p[w].y2pos, p[w].D / 2, p[w].colour, p[w].ID);
	}
}

void moveAll(vector<Particle> &p, int nop)
{
	double distx, disty, range;
	for (int m = 0; m < nop; m++) //second loop for movement of particles
	{
		for (int nb = m + 1; nb < nop; nb++)
		{
			distx = (p[m].xcpos - p[nb].xcpos) * (p[m].xcpos - p[nb].xcpos);
			disty = (p[m].ycpos - p[nb].ycpos) * (p[m].ycpos - p[nb].ycpos);
			range = 4 * p[m].Lmax * p[m].Lmax;
			if ((distx < range) && (disty < range))
				distForce(p[m], p[nb]); //calculate distances and overlapping forces between neighbours
		}
		p[m].forceInternal(); //calculate internal spring forces
		p[m].move(); //move the particles according to the net force
		p[m].clear(); //clear values for velocities and forces for next step
	}
}

void analyzeAll(vector<Particle> &p, int nop, Print &myPrint, int save, int ts)
{
	int amonb; //amount of neighbours
	double sumsop; //sum of sop over triangles
				   //int notr; //number of triangles
	for (int d = 0; d < nop; d++) //third loop to find defects
	{
		amonb = p[d].neighbours.size();
		sumsop = 0.0;
		//notr = 0;

		if (amonb > 0) //if using triangles than amonb > 1, if using neighbours than amonb > 0
		{
			for (int a = 0; a < amonb; a++) //loop over neighbours
			{
				//use to determine oop or sop
				int nb1 = p[d].neighbours[a]; //ID of neighbour 1

				double dtheta = p[d].theta - p[nb1].theta; //use for sop nb
				if (dtheta > pi / 2)
					dtheta = pi - dtheta;
				else if (dtheta < -pi / 2)
					dtheta += pi;

				sumsop += cos(2 * dtheta);

				int nonb1 = p[nb1].neighbours.size(); //set number of neighbours of nb1

				for (int b = 0; b < amonb; b++) //use loop for finding defects with rotational angle method
				{
					if (b > a)
					{
						int nb2 = p[d].neighbours[b]; //ID of neighbour 2
						int checknn = 0; //check if neighbours are also neighbours
						for (int c = 0;c < nonb1;c++)
						{
							int nb12 = p[nb1].neighbours[c]; //neighbours of first neighbour nb1
							if (nb2 == nb12)
								checknn = 1;
						}

						if (checknn == 1) //neighbours 1 and 2 are also neighbours of each other
						{
							double theta12 = (p[d].theta - p[nb1].theta) * (p[d].theta - p[nb1].theta);
							double theta23 = (p[nb1].theta - p[nb2].theta) * (p[nb1].theta - p[nb2].theta);
							double theta31 = (p[nb2].theta - p[d].theta) * (p[nb2].theta - p[d].theta);
							double ta = (pi / 2) * (pi / 2);
							//if only one angle of the three is larger than pi/2 than there is a defect
							if ((theta12 < ta && theta23 < ta && theta31 > ta) || (theta23 < ta && theta31 < ta && theta12 > ta) || (theta12 < ta && theta31 < ta && theta23 > ta))
							{
								p[d].ram = 1; //defect
							}
							//notr++;
							//double sumcos = cos(2 * p[d].theta) + cos(2 * p[nb1].theta) + cos(2 * p[nb2].theta); //use for sop triangle
							//double sumsin = sin(2 * p[d].theta) + sin(2 * p[nb1].theta) + sin(2 * p[nb2].theta);
							//sumsop += 2 * sqrt((sumcos / 6) * (sumcos / 6) + (sumsin / 6) * (sumsin / 6));
						}
					}
				}
			}
			p[d].oop = sumsop / amonb; //use for sop nb
									   //if (notr > 0)
									   //	p[d].oop = sumsop / notr; //gives sop from 0 to 1 for triangles
			if (p[d].ram == 1)
				p[d].colour = 10;
			else
				p[d].colour = (int)((1 - p[d].oop) * 80 + 180);
		}

		if (ts % save == 0)
			myPrint.print_data(ts, p[d].ID, p[d].xcpos, p[d].ycpos, p[d].theta, p[d].L, p[d].oop, p[d].ram); //save densities, in cluster 8 inputs are needed

		p[d].ram = 0; //reset defect after printing
		p[d].oop = 1.0; //reset oop after printing
	}
}

int main()
{
	int start = time(0);
	string dataNumber = "0123456789abcdefghijklmnopqrstuvwxyz"; //able to create up to 36 datafiles during one simulation
	int makefolder = 1;
	int nop; //number of particles
	int ts; //timestep
	int save = 300; //determines how many times data is saved, for nop = 3000 time runs till 650000. 400 gives nice set of data

	string folderName = "Test"; //name of data folder
	int startRound = 0;
	int noRounds = 1; //number of rounds
	
	for (int round = startRound; round < (startRound + noRounds); round++) //multiple rounds to gather data, set number after round < ..., to set number of data folders
	{
		cout << "Type of sim: loop nr = " << round << ", nop = " << Nmax << endl;
		cout << "Mu = " << growthRate / relaxTime << ", muDev = " << growthRateDev << ", L/D = " << maxLoverD << ", Ldev = " << maxLengthDev << ", k_i = " << ki / tau << ", k_o = " << ko / tau << endl;
		cout << "DF = RAM+SOPnb" << ", ON = " << orientNoise << endl;

		nop = 1;
		ts = 0;

		//setup printing for visualisation and data writing
		string dataName(1, dataNumber[round]);
		if (round != startRound)
			makefolder = 0; //only make run directory for first loop
		Print print;
		print.init(folderName, dataName, nop, makefolder);
		print.N = nop;

		//initialise particle
		vector<Particle> p(1, Particle(0.0, 0.0, 0.0, startLength, diameter, growthRate));

		while (nop < Nmax + 1) //timestep
		{
			if (ts % relaxTime == 0) //grow and relax the system
				growAll(p, nop); // t < 1 ns

			nop = p.size(); //determine new number of particles after growth
			if (ts % save == 0)
				print.N = nop; //update nop after growth and division

			moveAll(p, nop); //find forces, velocities and movement of particles

			if (ts % save == 0)
				analyzeAll(p, nop, print, save, ts); //analyze system to find order parameter and defects

			if ((round == startRound) && (ts % (save * 8) == 0))
				printAll(p, print, nop); //printing for visualisation

			if (ts % (save * 20) == 0)
				cout << "At time " << ts << " nop is " << nop << endl;

			ts++;
		}
		int finish = time(0);
		cout << "nop is " << nop << ", end time is " << ts << ", elapsed time is: " << finish - start << endl;
	}
	return 0;
}
