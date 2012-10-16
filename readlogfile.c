/* This program is written to read the log files which are the output of 
 * gsdfmain. 
 * The output of this program is in text file "logsumfile"
 */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#include "parameter.h"

int main(int argc,char *argv[]){

	int i,ifile;
	int period;
	char stemp[100];
	char stemp1[100], stemp2[100];
	FILE *logfp, *outfp;
	float evla,evlo;
	float phv;
	int pcsnum;

	if (argc == 1)
	{
		printf("Usage: readlogfile logfiles \n");
		return(0);
	}

	outfp=fopen("logsumfile","w");

	for (ifile=0;ifile<argc-1;ifile++)
	{
		logfp=fopen(argv[ifile+1],"r");
		if (logfp==NULL)
		{
			printf("Can't open log file %s, skip. \n", argv[ifile+1]);
			continue;
		}
		if(!feof(logfp))
		{
			fgets(stemp,100,logfp);
			sscanf(stemp,"%s %s %f %f\n", stemp1, stemp2, \
					&evla, &evlo);
			fprintf(outfp, "%f %f ", evla, evlo);
			for (i=0;i<5;i++)
				fgets(stemp,100,logfp);
			for (i=0;i<PNUM;i++)
			{
				fgets(stemp,100,logfp);
				sscanf(stemp,"Period: %d, New PhaseV: %f, Goodnum: %d \n", 
						&period, &phv, &pcsnum);
				fprintf(outfp, "%f %d ", phv, pcsnum);
			}
			fprintf(outfp, "\n");
		}
		fclose(logfp);
	}
	fclose(outfp);
}

