#include <stdio.h>
#include <math.h>
#include "PhaseField_InitialCondition.h"

double OregonatorTimeEvolution(int nx, int ny, double dx, double dy, double dt, int nt, double f, double q, double epsilon, double epsilonDash, int Circle_x, int Circle_y, int ObstacleRadius, double E, int PulseInterval, int savingInterval, int startTime, double Mu, double Mv, double Mw, int PulseStartTime)
{
	/* Define the vaiables to be used later */
	double uReactionTerm;
	double vReactionTerm;
	double wReactionTerm;
	double c1 = (dt)/(dx*dx);
	double c2 = (dt)/(dy*dy);
	double c3 = (c1/4.0);
	double c4 = (c2/4.0);
	double c5 = (dt/(2.0*dx));
	double c6 = (dt/(2.0*dy));
    

    	double pi = 3.14159265359;

	/***************************************************************************************************/	

	/* Memory Allocation */
	double u[nx+2][ny+2];
	double v[nx+2][ny+2];
	double w[nx+2][ny+2];

	double ubackup[nx+2][ny+2];
	double vbackup[nx+2][ny+2];
	double wbackup[nx+2][ny+2];

    	double PhiEquilibrium[nx+2][ny+2];
	double Ln_PhiEquilibrium[nx+2][ny+2];
	double Du[nx+2][ny+2];
	double Dv[nx+2][ny+2];
	double Dw[nx+2][ny+2];


	/***************************************************************************************************/	

	/* Load PhiEquilibrium values */

	/***************************************************************************************************/	

	FILE *f1;
	double buffer;
	f1 = fopen("Phi_Equilibrium.txt", "r");

	for(int i=1;i<=nx;i++)
	{
		for(int j=1;j<=ny;j++)
		{
			fscanf(f1, "%lf", &buffer);
			PhiEquilibrium[i][j] = buffer;
			//printf("%lf\t",PhiEquilibrium[i][j]);
		}
	}
	fclose(f1);

	/***************************************************************************************************/	

	/* Calculate the log of PhiEquilibrium array */
	// In C programming, function log is  natural logarithm, i.e its equivalent of ln.
	
	/***************************************************************************************************/	

	for(int i=1;i<=nx;i++)
	{
		for(int j=1;j<=ny;j++)
		{
			Ln_PhiEquilibrium[i][j] = log(PhiEquilibrium[i][j]); 
		//	printf("%lf\t",Ln_PhiEquilibrium[i][j]);
		}
	}

	/***************************************************************************************************/	

	/* Set the values of diffusion coefficient to 1 for all grid points (i,j) */

	/***************************************************************************************************/	

	
	for(int i=0;i<=nx+1;i++)
	{
		for(int j=0;j<=ny+1;j++)
		{
			Du[i][j] = 1.0;
		}
	}

	/***************************************************************************************************/	
	
	/* Setting the diffusion coefficient value very low inside the circle. 
	 * This creates an obstacle in excitable media. */

	/***************************************************************************************************/	

	for(int i=(Circle_x-ObstacleRadius);i<=(Circle_x+ObstacleRadius);i++)
	{

		for(int j=(Circle_y-ObstacleRadius);j<=(Circle_y+ObstacleRadius);j++)
		{
			if ( ( ( (i-Circle_x)*(i-Circle_x) ) + ( (j-Circle_y)*(j-Circle_y) ) ) <= (ObstacleRadius*ObstacleRadius) )
			{
				Du[i][j] = 0.0001;
			}
		}
	}

	/***************************************************************************************************/	

	/* Defining the obstacle array Dv here */ 

	/***************************************************************************************************/	

	for(int i=0;i<=nx+1;i++)
	{
		for(int j=0;j<=ny+1;j++)
		{
			Dv[i][j] = 0.6;
		}
	}


	/* Setting the diffusion coefficient value, Dv very low inside the obstacle. 
	 * THIS CAUSES ABNORMAL WAVE EMISSION FROM THE OBSTACLE. 
	 * This section is retained for future reference. Do not make this mistake. 
	 * */

	 // for(int i=(Circle_x-ObstacleRadius);i<=(Circle_x+ObstacleRadius);i++)
	 // {

	 // 	for(int j=(Circle_y-ObstacleRadius);j<=(Circle_y+ObstacleRadius);j++)
	 // 	{
	 // 		if ( ( ( (i-Circle_x)*(i-Circle_x) ) + ( (j-Circle_y)*(j-Circle_y) ) ) <= (ObstacleRadius*ObstacleRadius) )
	 // 		{

	 // 			Dv[i][j] = 0.6;
	 // 		}
	 // 	}
	 // }

	/***************************************************************************************************/	

	/* Defining the obstacle array Dw here */ 
	 
	/***************************************************************************************************/	

	for(int i=0;i<=nx+1;i++)
	{
		for(int j=0;j<=ny+1;j++)
		{
			Dw[i][j] = 1.12;
		}
	}
	
	/***************************************************************************************************/	

	/* Files for saving u and v data */

	/***************************************************************************************************/	

	FILE *fu;
	FILE *fv;
	FILE *fw;


	/***************************************************************************************************/	

	/* Target and Spiral initial conditions. 
	 * Note that this condition may not always give spiral in the medium. 
	 * In that case, chop-off the portion of the wave to generate the spiral. 
	 * */

	/* Replacing all the garbage values of the arrays with zeros.  */

	/***************************************************************************************************/	

	for(int i=0;i<=nx+1;i++)
	{
		for(int j=0;j<=ny+1;j++)
		{
			ubackup[i][j] = 0.0;
			vbackup[i][j] = 0.0;
			wbackup[i][j] = 0.0;
		}
	}


	if (startTime == 0)
	{
		 printf("Start time is zero\n");
	
		/* Load Initial Conditions */

		/* Target Wave Initial Conditions */

		// for(int i=200;i<=220;i++)
		// {
		// 	// printf("i = %d\n",i);
		// 	for(int j=200;j<=220;j++)
		// 	{
		// 		ubackup[i][j] = 0.9;
		// 	}
		// }

		/* Spiral Initial Conditons */

		for(int i=200;i<=nx;i++)
		{
			// printf("i = %d\n",i);
			for(int j=0;j<=ny;j++)
			{
				ubackup[i][j] = 0.9;
			}
		}

		for(int i=0;i<=300;i++)
		{
			for(int j=100;j<=ny;j++)
			{
				vbackup[i][j] = 0.05;
			}
		}
		
		/***************************************************************************************************/	

		/* Note that we are not giving any initial values for v and w variables. Again, a lesson! */

		/***************************************************************************************************/	

		// for(int i=0;i<=((nx/2)-1);i++)
		// {
		// 	for(int j=0;j<=ny;j++)
		// 	{
		// 		vbackup[i][j] = 0.1;
		// 	}
		// }


		// for(int i=0;i<=((nx/2)-1);i++) /* May have to change the values of i's and j's here */
		// {
		// 	for(int j=0;j<=ny;j++)
		// 	{
		// 		wbackup[i][j] = 0.1; 
		// 	}
		// }


		for(int i=(Circle_x-ObstacleRadius);i<=(Circle_x+ObstacleRadius);i++)
	 	{

	 		for(int j=(Circle_y-ObstacleRadius);j<=(Circle_y+ObstacleRadius);j++)
	 		{
	 			if ( ( ( (i-Circle_x)*(i-Circle_x) ) + ( (j-Circle_y)*(j-Circle_y) ) ) <= (ObstacleRadius*ObstacleRadius) )
	 			{

	 				ubackup[i][j] = 0;
					vbackup[i][j]=0;
					wbackup[i][j]=0;
	 			}
	 		}
	 	}
	}
	
	if (startTime != 0)
	{
		printf("Start time is not zero\n");
		FILE *f2, *f3, *f4;
		double ubuffer;
		double vbuffer;
		double wbuffer;

		char Ufile[32];
		char Vfile[32];
		char Wfile[32];
		sprintf(Ufile, "Data/u_%.7d.txt",startTime);
		sprintf(Vfile, "Data/v_%.7d.txt",startTime);
		sprintf(Wfile, "Data/w_%.7d.txt",startTime);
		
		f2 = fopen(Ufile, "r");
		f3 = fopen(Vfile, "r");
		f4 = fopen(Wfile, "r");
		

		// printf("Entering the loop\n");

		for(int i=1;i<=nx;i++)
		{
			for(int j=1;j<=ny;j++)
			{
				//printf("i=%d j=%d\n",i,j);
				fscanf(f2, "%lf", &ubuffer);
				ubackup[i][j] = ubuffer;
				fscanf(f3, "%lf", &vbuffer);
				vbackup[i][j] = vbuffer;
				fscanf(f4, "%lf", &wbuffer);
				wbackup[i][j] = wbuffer;
			}
		}
		fclose(f2);
		fclose(f3);
		fclose(f4);
	}




	/***************************************************************************************************/	

	/* Time Loop */

	/***************************************************************************************************/	

	for(int t=0;t<=nt;t++)
	{
		/* Array Exchange */

		for(int i=1;i<=nx;i++)
		{
			for(int j=1;j<=ny;j++)
			{
				u[i][j] = ubackup[i][j];
				v[i][j] = vbackup[i][j];
				w[i][j] = wbackup[i][j];
			}
		}
	

		/* Neumann Boundary Conditions */

		for(int i=1;i<=nx;i++)
		{
			u[i][0] = u[i][2];
			u[i][ny+1] = u[i][ny-1];
			v[i][0] = v[i][2];
			v[i][ny+1] = v[i][ny-1];
			w[i][0] = w[i][2];
			w[i][ny+1] = w[i][ny-1];
		}

		for(int j=1;j<=ny;j++)
		{
			u[0][j] = u[2][j];
			u[nx+1][j] = u[nx-1][j];
			v[0][j] = v[2][j];
			v[nx+1][j] = v[nx-1][j];
			w[0][j] = w[2][j];
			w[nx+1][j] = w[nx-1][j];
		}

		if ( ( t >= PulseStartTime) && (t <= PulseStartTime+PulseInterval) ) /* Apply E field */
		{
		
			printf("t=%d, applying E=%0.2f field\n",t,E);

			for(int i=1;i<=nx;i++)
			{
				for(int j=1;j<=ny;j++)
        			{
                      
						uReactionTerm = (q*w[i][j] - u[i][j]*w[i][j] + u[i][j] - (u[i][j]*u[i][j]))/epsilon;
						vReactionTerm = (u[i][j] - v[i][j]);
						wReactionTerm = (-q*w[i][j] - u[i][j]*w[i][j] + f*v[i][j])/epsilonDash; 
						ubackup[i][j] =   u[i][j] + (dt * uReactionTerm) + 
								( (Du[i][j]  * c1) * (u[i+1][j] - 2.0*u[i][j] + u[i-1][j]) ) + 
								( (Du[i][j]  * c2) * (u[i][j+1] - 2.0*u[i][j] + u[i][j-1]) ) + 
								( (Du[i][j]  * c3) * (Ln_PhiEquilibrium[i+1][j] - Ln_PhiEquilibrium[i-1][j]) * (u[i+1][j] - u[i-1][j]) ) + 
								( (Du[i][j]  * c4) * (Ln_PhiEquilibrium[i][j+1] - Ln_PhiEquilibrium[i][j-1]) * (u[i][j+1] - u[i][j-1]) ) + 
								( (Mu*E*c5) * (u[i+1][j] - u[i-1][j]));// -((Mu*E*c6) * (u[i][j+1] - u[i][j-1]));
						vbackup[i][j] = v[i][j] +  (dt * vReactionTerm) + 
							        ( (Dv[i][j] * c1)*(v[i+1][j] - 2.*v[i][j] + v[i-1][j]) ) + 
								( (Dv[i][j] * c2)*(v[i][j+1] - 2.*v[i][j] + v[i][j-1]) ) + 
								( (Dv[i][j]  * c3) * (Ln_PhiEquilibrium[i+1][j] - Ln_PhiEquilibrium[i-1][j]) * (v[i+1][j] - v[i-1][j]) ) + 
								( (Dv[i][j]  * c4) * (Ln_PhiEquilibrium[i][j+1] - Ln_PhiEquilibrium[i][j-1]) * (v[i][j+1] - v[i][j-1]) ) +
								( (Mv*E*c5)* (v[i+1][j] - v[i-1][j]));// -( (Mv*E*c6) * (v[i][j+1] - v[i][j-1]));  

						wbackup[i][j] = w[i][j] +  (dt * wReactionTerm) + 
							        ( (Dw[i][j] * c1)*(w[i+1][j] - 2.*w[i][j] + w[i-1][j]) ) + 
								( (Dw[i][j] * c2)*(w[i][j+1] - 2.*w[i][j] + w[i][j-1]) ) + 
								( (Dw[i][j]  * c3) * (Ln_PhiEquilibrium[i+1][j] - Ln_PhiEquilibrium[i-1][j]) * (w[i+1][j] - w[i-1][j]) ) + 
								( (Dw[i][j]  * c4) * (Ln_PhiEquilibrium[i][j+1] - Ln_PhiEquilibrium[i][j-1]) * (w[i][j+1] - w[i][j-1]) ) +
								( (Mw*E*c5)* (w[i+1][j] - w[i-1][j]));// -( (Mw*E*c6) * (w[i][j+1] - w[i][j-1]));  
				}
			}

		}

		else /* No electric field */
		{
			printf("t=%d\n",t);

			for(int i=1;i<=nx;i++)
			{
				for(int j=1;j<=ny;j++)
				{		 

						uReactionTerm = (q*w[i][j] - u[i][j]*w[i][j] + u[i][j] - (u[i][j]*u[i][j]))/epsilon;
						vReactionTerm = (u[i][j] - v[i][j]);
						wReactionTerm = (-q*w[i][j] - u[i][j]*w[i][j] + f*v[i][j])/epsilonDash; 
					
						ubackup[i][j] =   u[i][j] + (dt * uReactionTerm) + 
								( (Du[i][j]  * c1) * (u[i+1][j] - 2.0*u[i][j] + u[i-1][j]) ) + 
								( (Du[i][j]  * c2) * (u[i][j+1] - 2.0*u[i][j] + u[i][j-1]) ) + 
								( (Du[i][j]  * c3) * (Ln_PhiEquilibrium[i+1][j] - Ln_PhiEquilibrium[i-1][j]) * (u[i+1][j] - u[i-1][j]) ) + 
								( (Du[i][j]  * c4) * (Ln_PhiEquilibrium[i][j+1] - Ln_PhiEquilibrium[i][j-1]) * (u[i][j+1] - u[i][j-1]) );
						  
						vbackup[i][j] = v[i][j] + (dt* vReactionTerm) + 
							        ( (Dv[i][j] * c1)*(v[i+1][j] - 2.*v[i][j] + v[i-1][j]) ) + 
								( (Dv[i][j] * c2)*(v[i][j+1] - 2.*v[i][j] + v[i][j-1]) ) +
								( (Dv[i][j]  * c3) * (Ln_PhiEquilibrium[i+1][j] - Ln_PhiEquilibrium[i-1][j]) * (v[i+1][j] - v[i-1][j]) ) + 
								( (Dv[i][j]  * c4) * (Ln_PhiEquilibrium[i][j+1] - Ln_PhiEquilibrium[i][j-1]) * (v[i][j+1] - v[i][j-1]) ) ;

						wbackup[i][j] = w[i][j] + (dt* wReactionTerm) + 
							        ( (Dw[i][j] * c1)*(w[i+1][j] - 2.*w[i][j] + w[i-1][j]) ) + 
								( (Dw[i][j] * c2)*(w[i][j+1] - 2.*w[i][j] + w[i][j-1]) ) +
								( (Dw[i][j] * c3) * (Ln_PhiEquilibrium[i+1][j] - Ln_PhiEquilibrium[i-1][j]) * (w[i+1][j] - w[i-1][j]) ) + 
								( (Dw[i][j] * c4) * (Ln_PhiEquilibrium[i][j+1] - Ln_PhiEquilibrium[i][j-1]) * (w[i][j+1] - w[i][j-1]) ) ;
				}
			}

		}


		/* Saving Data into files */

		if((t%savingInterval) == 0)
		{
			char ufile[32];
			char vfile[32];
			char wfile[32];

			sprintf(ufile, "u_%.7d.txt",t+startTime);
			sprintf(vfile, "v_%.7d.txt",t+startTime);
			sprintf(wfile, "w_%.7d.txt",t+startTime);

			fu = fopen(ufile, "w");
			fv = fopen(vfile, "w");
			fw = fopen(wfile, "w");

			for(int i=1;i<=nx;i++)
			{
				for(int j=1;j<=ny;j++)
				{
					fprintf(fu, "%lf\n", u[i][j]);
					fprintf(fv, "%lf\n", v[i][j]);
					fprintf(fw, "%lf\n", w[i][j]);
				}
			}
			fclose(fu);
			fclose(fv);
			fclose(fw);
		}
	
	}
	return(0);
	}
