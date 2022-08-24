#include <stdio.h>
#include "PhaseField_InitialCondition.h"

double PhaseField_InitialCondition(int nx, int ny, int Circle_x, int Circle_y, int ObstacleRadius)
{
	double Phi_0[nx+2][ny+2];
	FILE *f1;

	/* Replacing all the garbage values by 1.0 in the array */

	for(int i=0;i<=nx+1;i++)
	{
		for(int j=0;j<=ny+1;j++)
		{
			Phi_0[i][j] = 1.0;
		}
	}

	/* Defining a circle with center (Circle_x,Circle_y) and radius "ObstacleRadius"
	 * and setting the value inside the circle to very small. This will create 
	 * an initial condition where the value of Phi is almost zero inside and 
	 * 1. outside with a discontinuity at the boundary of the obstacle.  
	 * */

	for(int i=(Circle_x-ObstacleRadius);i<=(Circle_x+ObstacleRadius);i++)
	{

		for(int j=(Circle_y-ObstacleRadius);j<=(Circle_y+ObstacleRadius);j++)
		{
			if ( ( ( (i-Circle_x)*(i-Circle_x) ) + ( (j-Circle_y)*(j-Circle_y) ) ) <= (ObstacleRadius*ObstacleRadius) )
			{
				Phi_0[i][j] = 0.00001;
			}
		}
	}

	/* Save the phase filed initial conditions into a file */

	f1 = fopen("Phi_0.txt", "w");
	for(int i=1;i<=nx;i++)
	{
		for(int j=1;j<=ny;j++)
		{
			fprintf(f1, "%lf\n", Phi_0[i][j]);
		}
	}
	fclose(f1);

	return(0);
}
