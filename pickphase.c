/* This function is used to get the average group velocity for each event, 
 * should be run before the gsdfmain.
 * The output file is the input for gsdfmain.
 * written by Ge Jin, jinwar@gmail.com
 * Jul. 2010
 */
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include "parameter.h"

struct station{
	float la, lo;
	float dist;
	char filename[50];
};

double degtorad(double deg) {
  return (deg * PI / 180);
}

double distance(double lat1, double lon1, double lat2, double lon2) {
  double theta, dist;
  theta = lon1 - lon2;
  dist = sin(degtorad(lat1)) * sin(degtorad(lat2)) + cos(degtorad(lat1))\
		 * cos(degtorad(lat2)) * cos(degtorad(theta));
  dist = acos(dist);
  dist = dist * 6371;
  return (dist);
}

int copysta(struct station *sta1, struct station *sta2){
	sta1->la=sta2->la;
	sta1->lo=sta2->lo;
	sta1->dist=sta2->dist;
	strcpy(sta1->filename, sta2->filename);
	return 1;
}

int swapsta(struct station *sta1, struct station *sta2){
	struct station temp;
	copysta(&temp,sta1);
	copysta(sta1,sta2);
	copysta(sta2,&temp);
}

int main(int argc,char *argv[]){

	char stemp[100];
	char path[40];
	char sacfilename[100];
	char key,ctemp;
	int i,j;
	int nlen, nerr, max = MAXDP;		// for reading the sac
	int Isgood;
	int stanum=0;
	FILE *evfp, *ofp, *sacmacrofp, *saclistfp;
	float neardist, fardist;
	float epidist;
	float neart0,fart0,neart1,fart1;
	float grv0,grv1;
	float grvmax,grvmin;
	float t0,t1,t3;
	float tt0,tt1;
	float data[MAXDP];
	float beg, del;
	float stla,stlo;
	float evla,evlo;
	struct station sta[1000];
	struct station tempsta;

	if (argc !=3)
	{
		printf("Usage: pickphase eventlist outputfile \n");
		exit(0);
	}
	
	evfp=fopen(argv[1],"r");
	ofp=fopen(argv[2],"w");
	fgets(stemp,40,evfp);
	while(!feof(evfp)){  // loop over the events
		sscanf(stemp,"%s\n", path);
		sprintf(stemp,"ls %s/*%3s*sac > saclist.txt\n",path,COMP);
		system(stemp);

		// Loop over the station list and sort the stations by epi-dist
		saclistfp=fopen("saclist.txt","r");
		fgets(stemp,100,saclistfp);
		stanum=0;
		while(!feof(saclistfp))
		{
			sscanf(stemp,"%s\n",sacfilename);
			// Check whether this station is within the interested area
			rsac1( sacfilename, data, &nlen, &beg, &del, &max,\
					&nerr, strlen(sacfilename) ) ;
			getfhv("stla", &stla, &nerr, 4);
			getfhv("stlo", &stlo, &nerr, 4);
			getfhv("evla", &evla, &nerr, 4);
			getfhv("evlo", &evlo, &nerr, 4);
			if (ISSELECTSTA)
				if (stla < MINSTALA || stla > MAXSTALA \
						|| stlo<MINSTALO || stlo>MAXSTALO)
				{
					printf("%s is not in the area!\n",sacfilename);
					fgets(stemp,100,saclistfp);
					continue;
				}
			sta[stanum].la=stla;
			sta[stanum].lo=stlo;
			sta[stanum].dist=distance(stla,stlo,evla,evlo);
			strcpy(sta[stanum].filename,sacfilename);
			stanum++;
			fgets(stemp,100,saclistfp);
		}
		// sort the stations by distance
		for (i=0;i<stanum;i++)
			for(j=0;j<stanum-1-i;j++)
				if (sta[j].dist>sta[j+1].dist)
					swapsta(&sta[j],&sta[j+1]);

		// pick the nearest station
		for (i=0;i<stanum;i++)
		{
			strcpy(sacfilename, sta[i].filename);

			sacmacrofp=fopen("sacmacrotemp.csh","w");
			fprintf(sacmacrofp,"sac <<!\n");
			fprintf(sacmacrofp,"r %s\n", sacfilename);
			fprintf(sacmacrofp,"w over\n");
			fprintf(sacmacrofp,"q\n");
			fprintf(sacmacrofp,"!\n");
			fclose(sacmacrofp);
			system("csh sacmacrotemp.csh");

			//getfhv("dist", &epidist, &nerr, 4);
			epidist=sta[i].dist;
			tt0=epidist/GRVMAXU;
			tt1=epidist/GRVMINL;

			// Change sac head
			sacmacrofp=fopen("sacmacrotemp.csh","w");
			fprintf(sacmacrofp,"sac <<!\n");
			fprintf(sacmacrofp,"r %s\n", sacfilename);
			fprintf(sacmacrofp,"ch t0 %7.1f\n",tt0);
			fprintf(sacmacrofp,"ch kt0 RB\n");
			fprintf(sacmacrofp,"ch t1 %7.1f\n",tt1);
			fprintf(sacmacrofp,"ch kt1 RE\n");
			fprintf(sacmacrofp,"ch t3 0\n");
			fprintf(sacmacrofp,"ch kt3 Quit\n");
			fprintf(sacmacrofp,"wh\n");
			fprintf(sacmacrofp,"q\n");
			fprintf(sacmacrofp,"!\n");
			fclose(sacmacrofp);
			system("csh sacmacrotemp.csh");

			// Band filter and make temp files
			sacmacrofp=fopen("sacmacrotemp.csh","w");
			fprintf(sacmacrofp,"sac <<!\n");
			for (j=0;j<PNUM;j++)
			{
				fprintf(sacmacrofp,"r %s\n", sacfilename);
				fprintf(sacmacrofp, "bp n 4 co %f %f passes 2 \n", \
						1.0/periods[j]*0.9, 1.0/periods[j]*1.1);
				//fprintf(sacmacrofp,"envelope\n");
				fprintf(sacmacrofp,"w temp_%1d.sac\n", j);
			}
			fprintf(sacmacrofp,"q\n");
			fprintf(sacmacrofp,"!\n");
			fclose(sacmacrofp);
			system("csh sacmacrotemp.csh");
			
			// Manually select t0 and t1
			sacmacrofp=fopen("sacmacrotemp.csh","w");
			fprintf(sacmacrofp,"sac <<!\n");
			fprintf(sacmacrofp,"r %s temp_?.sac\n", sacfilename);
			fprintf(sacmacrofp, "bp n 4 co 0.002 0.1 passes 1 \n");
			fprintf(sacmacrofp,"title &1,dist\n");
			fprintf(sacmacrofp,"qdp off\n");
			fprintf(sacmacrofp,"xlim t0 -500 t1 +3200\n");
			fprintf(sacmacrofp,"window x 0.05 0.9\n");
			fprintf(sacmacrofp,"beginwindow\n");
			fprintf(sacmacrofp,"ppk\n");
			fprintf(sacmacrofp,"wh\n");
			fprintf(sacmacrofp,"sc rm temp_?.sac\n");
			fprintf(sacmacrofp,"q\n");
			fprintf(sacmacrofp,"!\n");
			fclose(sacmacrofp);
			system("csh sacmacrotemp.csh");


			rsac1( sacfilename, data, &nlen, &beg, &del, &max,\
					&nerr, strlen(sacfilename) ) ;
			getfhv("t0", &t0, &nerr, 2);
			getfhv("t1", &t1, &nerr, 2);
			getfhv("t3", &t3, &nerr, 2);
			if (fabs(t0-tt0) > 1 || fabs(t1-tt1) > 1)
			{
				neart0=t0;
				neart1=t1;
				neardist=epidist;
				grvmax=epidist/t0;
				grvmin=epidist/t1;
				printf("grvmax=%f grvmin=%f \n",grvmax,grvmin);
				printf("neardist=%f \n",epidist);
				printf("Is it ok? y/n, q for quit\n");
				key=getchar();
				while ((ctemp=getchar())!='\n');
				if (key=='q')
				{
					fclose(evfp);
					fclose(ofp);
					exit(0);
				}

				if (key=='y')
				{
					Isgood=1;
					break;
				}
			}
			if (t3 > 1)
			{
				Isgood=0;
				break;
			}
			fgets(stemp,100,saclistfp);
		}

		// pick the farthest station
		for (i=stanum-1;i>-1;i--)
		{
			strcpy(sacfilename, sta[i].filename);

			sacmacrofp=fopen("sacmacrotemp.csh","w");
			fprintf(sacmacrofp,"sac <<!\n");
			fprintf(sacmacrofp,"r %s\n", sacfilename);
			fprintf(sacmacrofp,"w over\n");
			fprintf(sacmacrofp,"q\n");
			fprintf(sacmacrofp,"!\n");
			fclose(sacmacrofp);
			system("csh sacmacrotemp.csh");

			//getfhv("dist", &epidist, &nerr, 4);
			epidist=sta[i].dist;
			tt0=epidist/GRVMAXU;
			tt1=epidist/GRVMINL;

			// Change sac head
			sacmacrofp=fopen("sacmacrotemp.csh","w");
			fprintf(sacmacrofp,"sac <<!\n");
			fprintf(sacmacrofp,"r %s\n", sacfilename);
			fprintf(sacmacrofp,"ch t0 %7.1f\n",tt0);
			fprintf(sacmacrofp,"ch kt0 RB\n");
			fprintf(sacmacrofp,"ch t1 %7.1f\n",tt1);
			fprintf(sacmacrofp,"ch kt1 RE\n");
			fprintf(sacmacrofp,"ch t3 0\n");
			fprintf(sacmacrofp,"ch kt3 Quit\n");
			fprintf(sacmacrofp,"wh\n");
			fprintf(sacmacrofp,"q\n");
			fprintf(sacmacrofp,"!\n");
			fclose(sacmacrofp);
			system("csh sacmacrotemp.csh");

			// Band filter and make temp files
			sacmacrofp=fopen("sacmacrotemp.csh","w");
			fprintf(sacmacrofp,"sac <<!\n");
			for (j=0;j<PNUM;j++)
			{
				fprintf(sacmacrofp,"r %s\n", sacfilename);
				fprintf(sacmacrofp, "bp n 4 co %f %f passes 2 \n", \
						1.0/periods[j]*0.9, 1.0/periods[j]*1.1);
				//fprintf(sacmacrofp,"envelope\n");
				fprintf(sacmacrofp,"w temp_%1d.sac\n", j);
			}
			fprintf(sacmacrofp,"q\n");
			fprintf(sacmacrofp,"!\n");
			fclose(sacmacrofp);
			system("csh sacmacrotemp.csh");
			
			// Manually select t0 and t1
			sacmacrofp=fopen("sacmacrotemp.csh","w");
			fprintf(sacmacrofp,"sac <<!\n");
			fprintf(sacmacrofp,"r %s temp_?.sac\n", sacfilename);
			fprintf(sacmacrofp, "bp n 4 co 0.002 0.1 passes 1 \n");
			fprintf(sacmacrofp,"title &1,dist\n");
			fprintf(sacmacrofp,"qdp off\n");
			fprintf(sacmacrofp,"xlim t0 -500 t1 +3200\n");
			fprintf(sacmacrofp,"window x 0.05 0.9\n");
			fprintf(sacmacrofp,"beginwindow\n");
			fprintf(sacmacrofp,"ppk\n");
			fprintf(sacmacrofp,"wh\n");
			fprintf(sacmacrofp,"sc rm temp_?.sac\n");
			fprintf(sacmacrofp,"q\n");
			fprintf(sacmacrofp,"!\n");
			fclose(sacmacrofp);
			system("csh sacmacrotemp.csh");


			rsac1( sacfilename, data, &nlen, &beg, &del, &max,\
					&nerr, strlen(sacfilename) ) ;
			getfhv("t0", &t0, &nerr, 2);
			getfhv("t1", &t1, &nerr, 2);
			getfhv("t3", &t3, &nerr, 2);
			if (fabs(t0-tt0) > 1 || fabs(t1-tt1) > 1)
			{
				fart0=t0;
				fart1=t1;
				fardist=epidist;
				grv0=(fardist-neardist)/(fart0-neart0);
				grv1=(fardist-neardist)/(fart1-neart1);
				t0=neart0-neardist/grv0;
				t1=neart1-neardist/grv1;
				printf("fardist=%f \n",epidist);
				printf("grv0=%f t0=%f \n",grv0,t0);
				printf("grv1=%f t1=%f \n",grv1,t1);
				printf("Is it ok? y/n, q for quit\n");
				key=getchar();
				while ((ctemp=getchar())!='\n');
				if (key=='q')
				{
					fclose(evfp);
					fclose(ofp);
					exit(0);
				}

				if (key=='y')
				{
					Isgood=1;
					break;
				}
			}
			if (t3 > 1)
			{
				Isgood=0;
				break;
			}
			fgets(stemp,100,saclistfp);
		}

	if (Isgood==1)
		fprintf(ofp,"%s %f %f %f %f\n", path, grv0, t0,grv1, t1);

	fgets(stemp,40,evfp);

	}
	
	fclose(evfp);
	fclose(ofp);
}
