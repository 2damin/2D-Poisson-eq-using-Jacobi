#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#pragma warning(disable:4996)

#define pi 3.14159265358979323846264338327950288
#define np 50                   /*node number*/
#define del_x 2.0/(np-1.0)
#define del_y 2.0/(np-1.0)


double exact(int i, int j)                                             /*exact solution*/
{
	return cos(pi*i*del_x)*sin(pi + j*del_y);
}

double f(int i, int j)                                                    /*Poisson eq*/
{
	return cos(pi*i*del_x)*sin(j*del_y) + pi*pi*cos(pi*i*del_x)*sin(j*del_y);
}


/*---------------------------------------------------------------------------------------------------------------------*/

int main(void) {

	int i = 0;
	int j = 0;
	double L2error = 0.0;
	double maxerror = 0.0;
	double sum = 0.0;

	double *X = (double*)malloc(sizeof(double)*np*np);
	double *Y = (double*)malloc(sizeof(double)*np*np);
	double *u_differ = (double*)malloc(sizeof(double)*np*np);
	double *u = (double*)malloc(sizeof(double)*np*np);
	double *u_k = (double*)malloc(sizeof(double)*np*np);
	


	////////initial value////////
	for (j = 0; j<np; j++)
	{
		for (i = 0; i<np; i++)
		{
			u[j * np + i] = 0.0;
			u_k[j * np + i] = 0.0;
			u_differ[j * np + i] = 0.0;
			

			X[j * np + i] = del_x*i;
			Y[j * np + i] = del_y*j;
		}
	}

	/*Loop start*/
	do
	{
		//dirichlet boundary condition//
		for (j = 0; j < np; j++)
		{
			u[j*np] = exact(0, j);
			u_k[j*np] = exact(0, j);
			u[j*np + np - 1] = exact(np-1, j);
			u_k[j*np + np - 1] = exact(np-1, j);
		}

		for (i = 0; i < np; i++)
		{
			u[i] = exact(i, 0);
			u_k[i] = exact(i, 0);
			u[(np - 1)*np + i] = exact(i, np-1);
			u_k[(np - 1)*np + i] = exact(i, np-1);
		}

		//Jacobi method//
		for (j = 1; j<np-1; j++)
		{
			for (i = 1; i<np-1; i++)
			{
				u_k[j*np + i] = -(f(i, j) - (u[j*np + i - 1] + u[j*np + i + 1])/pow(del_x, 2) - (u[(j - 1)*np + i] + u[(j + 1)*np + i])/pow(del_y, 2))/(2 * (1/pow(del_x, 2) + 1/pow(del_y, 2)));
			}
		}

		
		//Maxerror//
		maxerror = 0.0;
		for (j = 1; j<np-1; j++)
		{
			for (i = 1; i<np-1; i++)
			{
				u_differ[j * np + i] = fabs(u_k[j * np + i] - u[j * np + i]);
				if (u_differ[j*np + i] > maxerror)
				{
					maxerror = u_differ[j*np + i];
				}
			}
		}

		printf("error : %f\n", maxerror);

		//Initialize//
		for (j = 0; j<np; j++)
		{
			for (i = 0; i<np; i++)
			{
				u[j * np + i] = u_k[j * np + i];
				u_differ[j * np + i] = 0.0;
			}
		}

	} while (maxerror >= 0.0001);


	//L2 Error//
	sum = 0.0;
	for (j = 0; j < np; j++)
	{
		for (i = 0; i < np; i++)
		{
			sum += pow((u[j*np + i] - exact(i, j)), 2);
		}
	}
	L2error = sqrt(sum)/ pow(np, 2);


	printf("end");


	////////////////////////////////////print out////////////////////////////////////////
	FILE*fp;
	char name[50];
	sprintf(name, "jacobi_serial_np%d.dat", np);
	fp = fopen(name, "w");
	fprintf(fp, "variables= X Y u  \n");
	fprintf(fp, "zone i=%d j=%d  \n",np,np);

	fprintf(fp, "L2 error: %f \n", L2error);
	for (j = 0; j<np; j++)
	{
		for (i = 0; i<np; i++)

		{
			fprintf(fp, " %f   %f   %f \n", X[j * np + i], Y[j * np + i], u[j*np+i]);
		}
	}
	fclose(fp);

	free(X);
	free(Y);
	free(u_differ);
	free(u);
	free(u_k);

	system("pause");
	return 0;
}
