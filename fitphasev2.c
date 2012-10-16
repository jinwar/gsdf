/* This program is used to fit the velocity using least square method. 
 * it also fit a inception as well as the slope(velocity)
 * The input txt file is a list of distance and travel time
 * written by Ge Jin, Columbia Univ.
 * jinwar@gmail.com
 */

#include <stdio.h>
#include <gsl/gsl_fit.h>

#define N 100000
     
int main(int argc,char *argv[]){

	int i, n;
	char stemp[100];
	double x[N], y[N];
	FILE *fp;

	double c1, c0, cov00, cov01, cov11, chisq;

	if (argc == 1)
	{
		printf("Usage: fitphasev2 distance_time_table \n");
		return(0);
	}

	// open in the data file
	fp=fopen(argv[1],"r");
	if (fp==NULL)
	{
		printf("Can't open the distance_time list file!\n");
		return 1;
	}

	n=0;
	fgets(stemp,100,fp);
	while(!feof(fp)){
		sscanf(stemp,"%lf %lf\n", &x[n],&y[n]);
		n++;
		fgets(stemp,100,fp);
	}
	fclose(fp);
     
	gsl_fit_linear(x, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

	fp=fopen("fitphasev2.out","w");
	fprintf(fp,"%e %e %e %e %e %e\n", c0, c1, cov00, cov01, cov11, chisq);
	fclose(fp);

	return 0;
}

