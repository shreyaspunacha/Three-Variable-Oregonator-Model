#include <stdio.h>
#include "PhaseField_InitialCondition.h"


double PhaseFieldTimeEvolution(int nx, int ny, int nt, double dx, double dy, double dt, double zi,int Circle_x,int Circle_y,int ObstacleRadius) 
{
	/* Define */
	double c1 = (dt*(zi*zi))/(dx*dx);
	double c2 = (dt*(zi*zi))/(dy*dy);

	/* Define arrays */
	double Phi_0[nx+2][ny+2];
	double Phi[nx+2][ny+2];
	double PhiBackup[nx+2][ny+2];

	/* Load the initial conditions from Phi_0.txt file into Phi_0 and Phi arrays. */

	FILE *f1, *f2;
	double buffer;

	f1 = fopen("Phi_0.txt", "r");

	for(int i=1;i<=nx;i++)
	{
		for(int j=1;j<=ny;j++)
		{
			fscanf(f1, "%lf", &buffer);
			Phi_0[i][j] = buffer;
			PhiBackup[i][j] = buffer;
		}
	}
	fclose(f1);

	/* Time Loop */


	for(int t=0;t<=nt;t++)
	{

		/* Array Exchange */
		
		
		for(int i=1;i<=nx;i++)
                {
                        for(int j=1;j<=ny;j++)
                                {
                                        Phi[i][j] = PhiBackup[i][j];
                                }
                }

		/* Neumann Boundary Conditions */

		for(int i=1;i<=nx;i++)
		{
		
			Phi[i][0] = Phi[i][2];
			Phi[i][ny+1] = Phi[i][ny-1];
		}

		for(int j=1;j<=ny;j++)
		{
			Phi[0][j] = Phi[2][j];
			Phi[nx+1][j] = Phi[nx-1][j];
		}
		
		
		for(int i=1;i<=nx;i++)
		{
				for(int j=1;j<=ny;j++)
				{
					PhiBackup[i][j] = Phi[i][j] + ( dt *( Phi_0[i][j]- Phi[i][j] )) + ( (c1)*(Phi[i+1][j] - 2.*Phi[i][j] + Phi[i-1][j]) ) + ( (c2)*(Phi[i][j+1] - 2.*Phi[i][j] + Phi[i][j-1]) );
				}
		}
	}

	/* Saving the last value of Phi, i.e Phi_Equilibrium into a file */

	f2 = fopen("Phi_Equilibrium.txt", "w");

	for(int i=1;i<=nx;i++)
	{
		for(int j=1;j<=ny;j++)
		{
			fprintf(f2, "%lf\n", PhiBackup[i][j]);
		}
	}
	fclose(f2);

	return(0);
}
