#include <iostream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <random>

using namespace std;
//height
#define WORLDSIZEX	100
//Width
#define WORLDSIZEY	100
#define NUMLANDMARKS 4
#define NUMPARTICLE 1000

Nw=sqrt(NUMPARTICLE*WORLDSIZEY/WORLDSIZEX);
Nh=sqrt(NUMPARTICLE*WORLDSIZEX/WORLDSIZEY);

//les balises
class Landmarks {
public:
	double x;
	double y;
	Landmarks();
	void set(double xx, double yy);
};
Landmarks::Landmarks()
{
	x = 0;
	y = 0;
}
void Landmarks::set(double xx, double yy)
{
	x = xx;
	y = yy;
	cout << "New landmarks (" << x << "," << y << ") setup!" << endl;
}

//le téléphone et les particules
class Robot {
public:
	double x;
	double y;
	double orientation;
	double forwardnoise;
	double turnnoise;
	double sensenoise;
	

	double sensors[NUMLANDMARKS];
	double weight;
	
	Robot();
	void set(double xx,double yy,double o);
	void setnoise(double fnoise, double tnoise, double snoise);
	void move(double turn, double forward);	
	void sense(Landmarks *landmark);
	double Gaussian(double mu,double sigma,double xh);
	double Euc();
	void w(double* sensors, Landmarks *landmark);
};

Robot::Robot(void)
{
	x = (double)rand()/RAND_MAX*WORLDSIZEX;
	y = (double)rand()/RAND_MAX*WORLDSIZEY;
	orientation = (double)rand()/RAND_MAX*2*M_PI;
	forwardnoise = 0;
	turnnoise = 0;
	sensenoise = 0;
	weight = 0;

//	cout << "x:" << x << endl;
//	cout << "y:" << y << endl;
//	cout << "Orientation:" << orientation << endl;
}
void Robot::set(double xx,double yy,double o)
{
	if(xx<0||xx>WORLDSIZEX)
	{
		cout << "X coordinate out of bound!" << endl;
	}
	if(yy<0||yy>WORLDSIZEY)
	{
		cout << "Y coordinate out of bound!" << endl;
	}
	if(o<0||o>2*M_PI)
	{
		cout << "Orientation must be in [0..2pi]!";
	}
	x = xx;
	y = yy;
	orientation = o;

	cout << "x:" << x << " y:" << y << " Orientation:" << orientation << endl;
}
void Robot::setnoise(double fnoise, double tnoise, double snoise)
{
	forwardnoise = fnoise;
	turnnoise = tnoise;
	sensenoise = snoise;
}
void Robot::move(double turn, double forward)
{	
	std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,turnnoise);
	double dist;
	if(forward < 0)
	{
		cout << "Robot cannot move backward!" << endl;
	}
	
	orientation += turn;
	orientation += distribution(generator);
	orientation = fmod(orientation, 2*M_PI);
	
	std::normal_distribution<double> distributionf(0.0,forwardnoise);
	dist = forward;
	dist += distributionf(generator);
	x += cos(orientation)*dist;
	y += sin(orientation)*dist;
	x = fmod(x,WORLDSIZEX);
	y = fmod(y,WORLDSIZEY);
	
		

	cout << "x:" << x << " y:" << y << " Orientation:" << orientation << endl;
}
void Robot::sense(Landmarks *landmark)
{
	std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0,sensenoise);
	for(int i=0;i<NUMLANDMARKS;i++)
	{
		sensors[i] = sqrt((x-landmark[i].x)*(x-landmark[i].x) + (y-landmark[i].y)*(y-landmark[i].y));
		sensors[i] += distribution(generator);
		cout << "sensors " << i << ":" << sensors[i]<< endl;
	}
}
double Robot::Gaussian(double mu,double sigma,double xh)
{
	return exp(-((mu - xh) * (mu - xh) / (sigma * sigma) / 2.0) / sqrt(2.0 * M_PI * (sigma * sigma)));
}
void Robot::w(double* sensors, Landmarks *landmark)
{	
	double dist = 0;
	weight = 1;
	for(int i=0;i<NUMLANDMARKS;i++)
	{
		dist = sqrt((x-landmark[i].x)*(x-landmark[i].x) + (y-landmark[i].y)*(y-landmark[i].y));
		weight *= Gaussian(dist,sensenoise, sensors[i]); 
	}
	//cout << "weight " << ":" << weight << endl;
}

/*
double Robot::Euc()
{
	double min = distolandmarks[0];
	for(int i=1;i<NUMLANDMARKS;i++)
	{
		if(min < distolandmarks[i])
		{
			min = distolandmarks[i];
		}
	}
	return min;
}
void Robot::w(Landmarks *landmark)
{
	sense(landmark);
	weight = 1 / (double)NUMPARTICLE;
	weight /= Euc();
	cout << "weight: " << weight << endl;
}
int test()
{
	Robot robot;
	robot.set(30.0, 50.0, M_PI/2);
	cout << robot.x << endl;
	cout << robot.y << endl;
	cout << robot.orientation << endl;
	robot.move(-M_PI/2,15.0);
	cout << robot.x << endl;
	cout << robot.y << endl;
	cout << robot.orientation << endl;
	double landmarks[NUMLANDMARKS][2] = {{20,20},{80,80},{20,80},{80,20}};
	Landmarks landmark[NUMLANDMARKS];
	for(int i=0;i<NUMLANDMARKS;i++)
	{
		landmark[i].set(landmarks[i][0],landmarks[i][1]);
	}
	robot.sense(landmark);
}
*/


void evaluation(Robot r, Robot* p)
{
	double sum = 0,dx = 0,dy = 0,err = 0;
	for(int i=0;i<NUMPARTICLE;i++) {
		dx = (p[i].x - r.x);
		dy = (p[i].y - r.y);
		err = sqrt(dx * dx + dy * dy);
        sum += err;
	}
	cout << "error:" << sum / (double)NUMPARTICLE << endl;
}

void ParticleFilter(Robot robot, Landmarks* landmark)
{
	//creat particles__Initialisation
	Robot p[NUMPARTICLE];
	for(int i=0;i<NUMPARTICLE;i++)
	{
		p[i].setnoise(0.05,0.05,5.0);
		for(int w=0;w<Nw-1;w++){
			for(int h=0;h<Nh-1;h++){
				p[i].set(w*(WORLDSIZEX/Nh),h/(WORLDSIZEY/Nw));
			}
		}

	}



	//use neff and not all repeat
	int T = 10; //times of repeates
	for(int i=0; i<T; i++) 
	{
		//Prédiction

		for(int i=0;i<NUMPARTICLE;i++) {
			p[i].x=
			p[i].y=
		}
		//calculate the weights;
		double sumweight = 0;
		for(int i=0;i<NUMPARTICLE;i++)
		{	
			//cout << i << endl;
			p[i].w(robot.sensors, landmark);
			sumweight += p[i].weight;
		}
		//normalize and cumulate the weights;
		double w[NUMPARTICLE];
		double cumulatesum = 0;
		for(int i=0;i<NUMPARTICLE;i++)
		{	
			//cout << i << endl;
			cumulatesum += p[i].weight / sumweight;
			w[i] = cumulatesum;
			//cout << w[i] << endl;	
		}	
		//create random table;
		double sampletable[NUMPARTICLE];
		srand (time(NULL));
		for(int i=0;i<NUMPARTICLE;i++)
		{		
			sampletable[i] = (double)rand()/RAND_MAX;
			//cout << i << ":" << sampletable[i] << endl;
		}
		//resampling
		Robot psample[NUMPARTICLE];
		for(int i=0;i<NUMPARTICLE;i++) {
			int pointer = (int) NUMPARTICLE / 2;
			int delta = (int) pointer / 2;
			while(1) {
				if(delta==0) {
					delta = 1;
				}
				if(sampletable[i] > w[pointer]) {
					if(pointer==NUMPARTICLE-1) {
						psample[i] = p[NUMPARTICLE-1];
						break;
					}
					if(sampletable[i] < w[pointer+1]) {
						psample[i] = p[pointer+1];
						break;
					} else {
						pointer += delta;
						delta = (int) delta / 2;
					}
				} else if(sampletable[i] < w[pointer]) {
					if(pointer==0) {
						psample[i] = p[0];
						break;
					}
					if(sampletable[i] > w[pointer-1]) {
						psample[i] = p[pointer];
						break;	
					} else {					
						pointer -= delta;
						delta = (int) delta / 2;
					}
				} else {
					psample[i] = p[pointer];
					break;
				}
			}
		}
		for(int i=0;i<NUMPARTICLE;i++) {
			p[i] = psample[i];
			//cout << p[i].x << "," << p[i].y << endl;
		}

		evaluation(robot, p);
	}
}

int main()
{
	//creat landmarks
	double landmarksNum[NUMLANDMARKS][2] = {{20,20},{80,80},{20,80},{80,20}};
	Landmarks mylandmark[NUMLANDMARKS];
	for(int i=0;i<NUMLANDMARKS;i++)
	{
		mylandmark[i].set(landmarksNum[i][0],landmarksNum[i][1]);
	}
	//creat robot
	Robot robot;
	robot.set(33.0, 54.0, M_PI/2);
	robot.setnoise(0.05,0.05,5.0);
	robot.sense(mylandmark);
	//particle filter
	ParticleFilter(robot, mylandmark);
}

//W
//Neff
//Estimation
