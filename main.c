#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PhaseField_InitialCondition.h"

int main(int argc,char *argv[])
{
	int nx = 300;
	int ny = 300; 
	int nt = 20000;
	int PhaseField_nt = 50000; // Increase this.
	int Circle_x = 150;
	int Circle_y = 150;
	int ObstacleRadius = 10;
	int PulseStartTime = 0 ; // The count starts from zero not starttime. 
	int PulseInterval = 30000;
	int savingInterval = 500;
	int startTime = 0; //atoi(argv[2]);

	double dx = 0.1; //1./6.;
	double dy = 0.1; //1./6.;
	double dt = 0.0001; //1./400.;
	double PhaseField_dt = 0.0001; //1./400.;
	double zi = 0.4;

	/*****************************************************************************************************
	 *
	 * The parameters are taken from the paper "Forced parallel drift of spiral waves in the 
	 * Belousov-Zhabotinsky reaction" by Bernd Schmidt and Stefan C. Mu ̈ller, PRE (1997)
	 *
	*****************************************************************************************************/

	double q = 0.002; 
	double f = 1.4;
	double epsilon = 0.01;
	double epsilonDash = 0.0001;
	double E = 0.0; //atof(argv[1]); 
// 	double Ey = 0.0;
	double Mu = 0.0; // Setting the value of the drift term of u to zero as advised by the referee of JPC.  
	double Mv = -2.0; 
	double Mw = -2.0; // CHANGE THIS VALUE ACCORDING TO THE PAPER 
//    	double T = atoi(argv[2]);

	/* Modify the code so that not to run the first two functions below, if the relevent files already exists */

	PhaseField_InitialCondition(nx, ny, Circle_x, Circle_y, ObstacleRadius);
	PhaseFieldTimeEvolution(nx, ny, nt, dx, dy, dt, zi, Circle_x, Circle_y, ObstacleRadius);
	OregonatorTimeEvolution(nx, ny, dx, dy, dt, 
			nt, f, q, epsilon, epsilonDash, Circle_x, 
			Circle_y, ObstacleRadius, E, 
			PulseInterval,
			savingInterval, startTime, Mu, Mv, Mw, PulseStartTime);
	
	return(0);

}
