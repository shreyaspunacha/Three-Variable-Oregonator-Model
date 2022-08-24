#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PhaseField_InitialCondition.h"

int main(int argc,char *argv[])
{
	int nx = 300; 				// Grid points along x direction.
	int ny = 300; 				// Grid points along x direction.
	int nt = 50000; 			// Total number of time steps.
	int PhaseField_nt = 100000; 		// Total number of time steps for generation of phase field file. Increase this while doing the simulations.
	int Circle_x = 150;			// x coordinate of obstacle center. 
	int Circle_y = 150;			// y coordinate of obstacle center. 
	int ObstacleRadius = 10;		// Radius of the obstacle.
	int PulseStartTime = 0; 		// The time at which the electric field must be started. The count starts from zero, not starttime. 
	int PulseInterval = 100000;     	// The frequency at which pulses are delivered. 
	int savingInterval = 100;       	// The figure plotting interval.
	int startTime = 100000;			// The starttime for the code. if 0, it will run the initial conditions, else it will pickup the files from Data folder 
	// int startTime = atoi(argv[2]);       // The command in case a shell script is used to run the program.

	double dx = 0.1; 			// Grid spacing along x direction.
	double dy = 0.1;			// Grid spacing along y direction. 
	double dt = 0.0001;			// time step according to the Von Neumann stability condition. Calculate it yourself! 
	double PhaseField_dt = 0.0001;		// time step for the phase field generation, according to the Von Neumann stability condition. 
	double zi = 0.4;			// The value of the phase field method parameter.

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

	double E = 1.1;				// The strength of the electric field. 
	// double E = atof(argv[1]); 		// The command in case a shell script is used to run the program.
	// double Ey = 0.0; 			// In case the electric field is applied along the y direction. In that case change the previous line to E=Ex.
	double Mu = 0; 				// Setting the value of the drift term of u to zero as advised by the referee of Journal of Physical Chemistry. 
						// Look at correspondense e-mails for details.  
	double Mv = -2.0;			// Drift term for v variable. 
	double Mw = 1.0; 			// Drift term for w variable. 
	//double T = atoi(argv[2]);

	/****************************************************************************************************************************
	 *
	 * YOU PROBABLY NEED NOT CHANGE ANYTHING BELOW THIS.
	 *
	 * The function of the if-else condition below is as follows: 
	 * if the file "Phi_Equilibrium already exists from the previous run, it skips running the phase field generation code.
	 * In case it does not exist, it goes to the else section and generates the Phi_Equilibrium file. 
	 * Note that, Phi_Equilibrium need to be generated only once.
	 *
	****************************************************************************************************************************/

	FILE *file = fopen("Phi_Equilibrium.txt", "r");
	if (file)
	{
		printf("Phase Field Equilibrium File Already Exists\n");
		printf("Integrating Oregonator Equations\n");

		OregonatorTimeEvolution(nx, ny, dx, dy, dt, 
			nt, f, q, epsilon, epsilonDash, Circle_x, 
			Circle_y, ObstacleRadius, E, 
			PulseInterval,
			savingInterval, startTime, Mu, Mv, Mw, PulseStartTime);

		fclose(file);

	}
	else
	{
		printf("Phase Field Equilibrium does not exist\n");
		printf("Computing Phase Filed Equilibrium\n");

		PhaseField_InitialCondition(nx, ny, Circle_x, Circle_y, ObstacleRadius);
		PhaseFieldTimeEvolution(nx, ny, nt, dx, dy, dt, zi, Circle_x, Circle_y, ObstacleRadius);


		printf("Integrating Oregonator Equations\n");
		OregonatorTimeEvolution(nx, ny, dx, dy, dt, 
			nt, f, q, epsilon, epsilonDash, Circle_x, 
			Circle_y, ObstacleRadius, E, 
			PulseInterval,
			savingInterval, startTime, Mu, Mv, Mw, PulseStartTime);
			
	}
		
		

	
	return(0);

}
