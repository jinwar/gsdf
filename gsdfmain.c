/*
 * This problem is used to meas[0]ure the relative travel time of surface wave to different 
 * stations by calculating and fitting the cross-correlation of the real data.
 * It can automatically solve the "circle skipping" problem by grading each meas[0]urement
 * in GSDF programs.
 * The basic idea is using 3 grading system:
	 * 1. Measurement comparing to the reference model
	 * 2. Measurement comparing to neighbor frequency
	 * 3. Measurement comparing to neighbor stations.
 */ 

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
//#include<sacio.h>
//#include<sac.h>

//files to set parameters
#include "parameter.h"
// Define the functions to calculate the absolute time to the 2000.1.1
#include "abtime.h"
// define data structure
#include "datastruct.h" 
// define functions to do station operations
#include "stafunction.h"
// Create the data base from the sacfiles
#include "createdatabase.h"
// Define the function to calculate the path grade
#include "pathgrading.h"

int main(int argc,char *argv[]){
	STA *stahead,*stap;
	CS *cshead,*csp;
	EV *evp, *evhead=NULL,*evq=NULL;
	int i;
	int skipn,skipnbefore;
	int csnum=0;
	int stationnum;
	int eventn;
	char stemp[100];
	double r,max,min;
	double grade[3];
	FILE *evfp, *opfp;

	if (argc == 1)
	{
		printf("Usage: gsdfmain event_path_list \n");
		return(0);
	}

	evfp=fopen(argv[1],"r");

	if (evfp==NULL)
	{
		printf("Can't open the event path list file!\n");
		return 1;
	}

	/* Read in the model files */
	ReadMod();
	
/* Begin the loop of events	*/
	eventn=0;
	fgets(stemp,100,evfp);
while(!feof(evfp)){

	/* build event chain */
	eventn++;
	evp=(EV *)malloc(EVLEN);	

	sscanf(stemp,"%s %f %f %f %f\n",evp->path,&evp->grv0,&evp->t0, &evp->grv1, &evp->t1);

	if (eventn==1)
	{
		evhead=evp;
		evq=evp;
		evp->n=NULL;
	}
	else
	{
		evq->n=evp;
		evp->n=NULL;
		evq=evp;
	}

	/* Build station chain for this event and fill the event information */
	stahead=CreateStationChain(evp);
	evp->stachain=stahead;

	// Count the number of the staitons for this event
	stap=stahead;
	stationnum=0;
	while (stap!=NULL)
	{
		//printf("%s %f %f %s\n",stap->filename, stap->la, stap->lo, stap->staname);
		stap=stap->n;
		stationnum++;
	}
	printf("Station Number is %d \n", stationnum);
	evp->stationnum=stationnum;

	// Build the cross-correlation meas[0]urement chain 
	cshead=CreateMainCSChain(stahead);
	evp->cschain=cshead;
	csp=cshead;
	skipnbefore=0;
	while (csp!=NULL)
	{
		csp=csp->n;
	}

	// Average the both way measurement
	//AverageBothWay(cshead);


	// Grading the path using common friends and calculate the total grade
	printf("Grading path.....\n");
	pathgrading(cshead);

	ErrLevel(cshead);

	/*
	printf("Grading model.....\n");
	RefModGrading(cshead);
	*/

	/*
	printf("Grading station.....\n");
	StationGrading(evp);
	*/

	printf("Writing output.....\n");
	PrintOutputEV(evp);

	fgets(stemp,100,evfp);
	}  /* End of loop event */

	
/*
	printf("gathering station information.....\n");
	StaInformGather(evhead);
*/

	// Print output
}
