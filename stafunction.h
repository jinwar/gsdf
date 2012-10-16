
#define N 30		// Constant control the number of circle skipping

// pre-declare the functions
double degtorad(double deg);
double distance(double lat1, double lon1, double lat2, double lon2);
double StationDistance(STA *sta1, STA *sta2);
double GaussianRand(int n);
STA *GetOtherStaion (STA *sta, CS *csp);
int IsGoodData(CS *csp, int pid);
float GetMinStaSNR(CS *csp,int i);

/* 
 * Function to find the stations within a certain distance from one station
 * It will create a station chain connected by tempn and return the head.
*/
STA *FindStaByDist( STA *sta, STA *stahead, double dist, double mindist){
	int targetstanum;  // number of stations within this dist
	double stadist,maxdist;
	STA *p,*q;
	STA *head,*tail;
	STA *maxdiststap;

	p=stahead;
	targetstanum=0;
	head=NULL;
	while (p!=NULL)
	{
		if (p==sta)
			stadist=dist+1;
		else
			stadist=StationDistance(sta,p);
		//printf("%f\n",stadist);
		if (stadist<dist && stadist > mindist)
		{
			if (targetstanum==0) 
			{
				head=p;
				tail=p;
				p->tempn=NULL;
			}
			else
			{
				tail->tempn=p;
				tail=p;
				p->tempn=NULL;
			}
			targetstanum++;
		}
		p=p->n;
	}
	while (targetstanum > MAXCSNUM)
	{
		maxdist=0;
		p=head;
		while (p!=NULL)
		{
			stadist=StationDistance(sta,p);
			if (stadist > maxdist)
			{
				maxdist=stadist;
				maxdiststap=p;
			}
			p=p->tempn;
		}
		if (head==maxdiststap)
			head=maxdiststap->tempn;
		else
		{
			p=head;
			while (p->tempn!=maxdiststap && p->tempn!=NULL){
				p=p->tempn;
			}
			if (p->tempn!=NULL)
				p->tempn=p->tempn->tempn;
		}
		targetstanum=targetstanum-1;
	}
	//printf("%d\n",targetstanum);
	return(head);
}

/* Function to put \0 into the end of station name */
int EndStaname(STA *sta){
	int i;
	
	for (i=0;i<8;i++)
		if (sta->staname[i] == ' ')
		{
			sta->staname[i] = '\0';
			return 1;
		}

	return 0;

}

/* Function to find the cross co-relation meas[0]urement between two stations */
CS *FindStaCo( STA *sta1, STA *sta2){
	CSCH *p,*q;

	if (sta1==sta2)
		return NULL;
	p=sta1->csch;
	while (p!=NULL)
	{
		q=sta2->csch;
		while (q!=NULL)
		{
			if (p->cs==q->cs)
				return(p->cs);
			q=q->n;
		}
		p=p->n;
	}
	return(NULL);
}


/* Function to fill costa structure with example data */
/*
void FillCostaDataTest(CS *csp){
	double rightdt,refdt,contdt;
	double v=3.0;
	double a,b,c;
	int i;
	
	rightdt=(csp->sta2->la-csp->sta1->la)*111.0/v+GaussianRand(2*N);
	refdt=rightdt+GaussianRand(N);
	contdt=rightdt+GaussianRand(N);

	csp->dist=StationDistance(csp->sta1,csp->sta2);
	csp->CurveFitErr=1;

	//Set the time for each circle
	for (i=0;i<3;i++)
		csp->meas[0][i].dt=rightdt+(i-1)*10;
	//Set the Grade for reference model
	for (i=0;i<3;i++)
	{
		a=fabs(csp->meas[0][i].dt-refdt);
		b=fabs(csp->meas[0][(i+1)%3].dt-refdt);
		c=fabs(csp->meas[0][(i+2)%3].dt-refdt);
		csp->meas[0][i].refG=100/(1+a/b+a/c);
	}
	//Set the Grade for continuity
	for (i=0;i<3;i++)
	{
		a=fabs(csp->meas[0][i].dt-contdt);
		b=fabs(csp->meas[0][(i+1)%3].dt-contdt);
		c=fabs(csp->meas[0][(i+2)%3].dt-contdt);
		csp->meas[0][i].contG=100/(1+a/b+a/c);
	}
	// Initialize the path grade
	for (i=0;i<3;i++)
	{
		csp->meas[0][i].pathG=0;
		csp->meas[0][i].localG=csp->meas[0][i].refG+csp->meas[0][i].contG;
	}
}
*/


/* Function to find common friends between two stations.
 * Return a station chain's head connected by tempn*/
STA *FindCommonFriends (CS *csp){
	CSCH *cschpsta1, *cschpsta2;
	STA *sta1,*sta2;
	STA *friendsta1, *friendsta2;
	STA *cfhead=NULL, *cftail, *cfp;
	int cfnumber;

	sta1=csp->sta1;
	sta2=csp->sta2;

	cfnumber=0;
	cschpsta1=csp->sta1->csch;
	while (cschpsta1!=NULL)
	{
		friendsta1=GetOtherStaion(sta1,cschpsta1->cs);
		cschpsta2=csp->sta2->csch;
		while (cschpsta2!=NULL)
		{
			friendsta2=GetOtherStaion(sta2,cschpsta2->cs);
			if (friendsta1==friendsta2)
			{
				if (cfnumber==0)
				{
					cfhead=friendsta1;
					cftail=friendsta1;
				}
				else
				{
					cftail->tempn=friendsta1;
					cftail=friendsta1;
				}
				friendsta1->tempn=NULL;
				cfnumber++;
			}
			cschpsta2=cschpsta2->n;
		}
		cschpsta1=cschpsta1->n;
	}
	return (cfhead);
}
	
/* Get the other station from a costa */
STA *GetOtherStaion (STA *sta, CS *csp) {
	if (csp->sta1 == sta )
		return csp->sta2;
	else
		return csp->sta1;
}

/* Get the time difference FROM sta1 TO sta2 */
double Getdt(STA *sta1, STA *sta2, int i){
	CS *csp;
	csp=FindStaCo(sta1,sta2);
	if (csp->sta1==sta1)
		return csp->meas[0][i].dt;
	else
		return -csp->meas[0][i].dt;
}

/* Get the local best dt FROM sta1 TO sta2 */
double Getbestdt(STA *sta1, STA *sta2, int pid){
	CS *csp;
	int i, besti=1;
	double maxgrade=0;
	csp=FindStaCo(sta1,sta2);
	for (i=0;i<CIRCNUM;i++)
		if (maxgrade<csp->meas[pid][i].localG)
		{
			maxgrade=csp->meas[pid][i].localG;
			besti=i;
		}
	if (csp->sta1==sta1)
		return csp->meas[pid][besti].dt;
	else
		return -csp->meas[pid][besti].dt;
}

/* Get the local best dt of CS for a certain period  */
double GetBestdtCS(CS *csp, int pid){
	int i, besti=1;
	double maxgrade=-1;
	for (i=0;i<CIRCNUM;i++)
		if (maxgrade<csp->meas[pid][i].localG)
		{
			maxgrade=csp->meas[pid][i].localG;
			besti=i;
		}
	return csp->meas[pid][besti].dt;
}

/* Get the local best index of CS for a certain period  */
int GetBestCircCS(CS *csp, int pid){
	int i, besti=1;
	double maxgrade=0;
	for (i=0;i<CIRCNUM;i++)
		if (maxgrade<csp->meas[pid][i].localG)
		{
			maxgrade=csp->meas[pid][i].localG;
			besti=i;
		}
	return besti;
}

/* Get the total best dt of CS for a certain period  */
double GetTotalBestdtCS(CS *csp, int pid){
	int i, besti=1;
	double maxgrade=0;
	for (i=0;i<CIRCNUM;i++)
		if (maxgrade<csp->meas[pid][i].totalG)
		{
			maxgrade=csp->meas[pid][i].totalG;
			besti=i;
		}
	return csp->meas[pid][besti].dt;
}

/* Get the total best index of CS for a certain period  */
int GetTotalBestCircCS(CS *csp, int pid){
	int i, besti=1;
	double maxgrade=0;
	for (i=0;i<CIRCNUM;i++)
		if (maxgrade<csp->meas[pid][i].totalG)
		{
			maxgrade=csp->meas[pid][i].totalG;
			besti=i;
		}
	return besti;
}

/* Get the total best dt FROM sta1 TO sta2 */
double GetTotalBestdtSTA(STA *sta1, STA *sta2, int pid){
	CS *csp;
	int i, besti=1;
	double maxgrade=0;
	csp=FindStaCo(sta1,sta2);
	for (i=0;i<CIRCNUM;i++)
		if (maxgrade<csp->meas[pid][i].totalG)
		{
			maxgrade=csp->meas[pid][i].totalG;
			besti=i;
		}
	if (csp->sta1==sta1)
		return csp->meas[pid][besti].dt;
	else
		return -csp->meas[pid][besti].dt;
}

/* Function to calculate the distance between stations */
double StationDistance(STA *sta1, STA *sta2) {
	return (distance(sta1->la,sta1->lo,sta2->la,sta2->lo));
}

/* Function to calculate the distance between the epicenter distance */
double epidiffCS(CS *csp) {
	float epidist1, epidist2;
	epidist1=distance(csp->sta1->la, csp->sta1->lo, csp->ev->la, csp->ev->lo);
	epidist2=distance(csp->sta2->la, csp->sta2->lo, csp->ev->la, csp->ev->lo);
	csp->sta1->epidist=epidist1;
	csp->sta2->epidist=epidist2;

	return epidist1-epidist2;
}

/* Function to calculate the distance between two points based on their 
 * latitude and longtitude. */
double distance(double lat1, double lon1, double lat2, double lon2) {
  double theta, dist;
  theta = lon1 - lon2;
  dist = sin(degtorad(lat1)) * sin(degtorad(lat2)) + cos(degtorad(lat1))\
		 * cos(degtorad(lat2)) * cos(degtorad(theta));
  dist = acos(dist);
  dist = dist * 6371;
  return (dist);
}

double degtorad(double deg) {
  return (deg * PI / 180);
}

/* Function to generate a approximate zero mean Gaussian distributed random 
 * number. */
double GaussianRand(int n){
	int i;
	double sum;

	sum=0;
	for (i=0;i<n;i++)
		sum+=rand()%100;
	sum=sum*1.0;
	sum=sum/n;
	sum=sum-50;
	return sum;
}

/* Calculate Gaussian distributed grade */
double GaussianGrade(double x, double x0, double sigma){
	double g;
	g = 100*exp(-(x-x0)*(x-x0)/sigma/sigma);
	return g;
}

/* Function to check whether a station within station tolerant distance (STA_TOL_DIST) exist
 */
int IsExistSta(STA *stahead, STA *sta){
	STA *stap;
	float dist=STA_TOL_DIST;

	stap=stahead;
	while (stap!=NULL) 
	{
		if (StationDistance(stap, sta) < dist)
			return 1;
		stap=stap->n;
	}

	return 0;
}

/* Function to check whether a station exist in the station information chain */
int IsExistStaStaInform(STAINF *stainfhead, STA *sta){
	STAINF *stainfp;
	float dist=STA_TOL_DIST;

	stainfp=stainfhead;
	while (stainfp!=NULL) 
	{
		if (distance(stainfp->la,stainfp->lo, sta->la, sta->lo) < dist)
			return 1;
		stainfp=stainfp->n;
	}

	return 0;
}
		
/* Function to check whether a station within station tolerant distance (STA_TOL_DIST) exist
   in the temperate station chain*/
int IsExistStaTemp(STA *stahead, STA *sta){
	STA *stap;
	float dist=STA_TOL_DIST;

	stap=stahead;
	while (stap!=NULL) 
	{
		if (StationDistance(stap, sta) < dist)
			return 1;
		stap=stap->tempn;
	}

	return 0;
}

/* Get phase velocity for a certern period */
float GetPhV(int T){
	int i;
	
	for (i=0;i<MODNUM;i++)
		if (T - phvmodel[i][0] < 0.1)
			return phvmodel[i][1];

	if (T < phvmodel[0][0])
		return phvmodel[0][1];
	else
		return phvmodel[MODNUM-1][1];
}

/* Get group velocity for a certern period */
float GetGrV(int T){
	int i;
	
	for (i=0;i<MODNUM;i++)
		if (T - grvmodel[i][0] < 0.1)
			return grvmodel[i][1];

	if (T < grvmodel[0][0])
		return grvmodel[0][1];
	else
		return grvmodel[MODNUM-1][1];
}

/* windowing the data using Hanning pater  */
int HanWindow(char *filename, float center, float halfwidth)
{
	char newfilename[100];
	int nlen, nerr, max = MAXDP;	// for sac file reading
	int i;
	int n;
	float beg, del;
	float b,e;
	float data[MAXDP];
	float begin, end;

	begin=center-halfwidth;
	end=center+halfwidth;

	//Read in the sac file
	rsac1( filename, data, &nlen, &beg, &del, &max, &nerr, strlen(filename) ) ;
	if ( nerr != 0 ) {
		fprintf(stderr, "Error reading in SAC file: %s\n", filename);
		exit ( nerr ) ;
	}

	getfhv("b", &b, &nerr, 1);
	n=floor(2*halfwidth/del);
	
	// Do the windowing
	for (i=0;i<nlen;i++)
	{
		if ( i*del+b < begin || i*del+b > end )
			data[i]=0;
		else
		{
			data[i]=data[i]*(0.5-0.5*cos(2*PI*(floor((i*del+b-begin)/del))/(n-1)));
		}
	}
	sprintf(newfilename, "%s.w", filename);

	wsac1( newfilename, data, &nlen, &beg, &del, &nerr, strlen(newfilename) ) ;

}

/* windowing the data using Gaussian taper  */
int GausWindow(char *filename, float center, float halfwidth)
{
	char newfilename[100];
	int nlen, nerr, max = MAXDP;	// for sac file reading
	int i;
	int n;
	float beg, del;
	float b,e;
	float data[MAXDP];
	float begin, end;

	begin=center-halfwidth;
	end=center+halfwidth;

	//Read in the sac file
	rsac1( filename, data, &nlen, &beg, &del, &max, &nerr, strlen(filename) ) ;
	if ( nerr != 0 ) {
		fprintf(stderr, "Error reading in SAC file: %s\n", filename);
		exit ( nerr ) ;
	}

	getfhv("b", &b, &nerr, 1);
	n=floor(2*halfwidth/del);
	
	// Do the windowing
	for (i=0;i<nlen;i++)
	{
		if ( i*del+b < begin || i*del+b > end )
			data[i]=0;
		else
		{
			data[i]=data[i];
		}
	}
	sprintf(newfilename, "%s.w", filename);

	wsac1( newfilename, data, &nlen, &beg, &del, &nerr, strlen(newfilename) ) ;

}

/* Find the peak in order to window the xcor waveform */
float FindXcorPeak(char *filename, float center, float halfwidth)
{
	char newfilename[100];
	int nlen, nerr, max = MAXDP;	// for sac file reading
	int i,maxi;
	int n;
	float beg, del;
	float b,e;
	float data[MAXDP];
	float begin, end;
	float maxdata;

	begin=center-halfwidth;
	end=center+halfwidth;

	//Read in the sac file
	rsac1( filename, data, &nlen, &beg, &del, &max, &nerr, strlen(filename) ) ;
	if ( nerr != 0 ) {
		fprintf(stderr, "Error reading in SAC file: %s\n", filename);
		exit ( nerr ) ;
	}

	getfhv("b", &b, &nerr, 1);
	n=floor(2*halfwidth/del);
	maxdata=0;maxi=0;
	
	// Do the windowing
	for (i=0;i<nlen;i++)
	{
		if ( i*del+b > begin && i*del+b < end )
			if ( maxdata < data[i] )
			{
				maxdata = data[i];
				maxi= i;
			}
	}
	return maxi*del+b;
}

/* Function to calculate total energy contained by this sacfile */
float AvgEnergy(char *filename){
	int nlen, nerr, max = MAXDP;	// for sac file reading
	int i;
	int n;
	float beg, del;
	float b,e;
	float data[MAXDP];
	float begin, end;
	float sum=0;

	rsac1( filename, data, &nlen, &beg, &del, &max, &nerr, strlen(filename) ) ;
	if (nerr!=0)
	{
		printf("Error in reading %s \n",filename);
		return 0;
	}

	for (i=0;i<nlen;i++)
		sum=sum+data[i]*data[i];

	if (filename[strlen(filename)-1]=='w')
		return sum/HALFWINDOW/2*del;
	else
		return sum/nlen;
}

/* windowing the data using Hanning pater  */
int HanWindowFlat(char *filename, float center, float halfwidth,float taperwidth){
	char newfilename[100];
	int nlen, nerr, max = MAXDP;	// for sac file reading
	int i;
	int n;
	float beg, del;
	float b,e;
	float data[MAXDP];
	float begin, end;

	begin=center-halfwidth;
	end=center+halfwidth;

	//Read in the sac file
	rsac1( filename, data, &nlen, &beg, &del, &max, &nerr, strlen(filename) ) ;
	if ( nerr != 0 ) {
		fprintf(stderr, "Error reading in SAC file: %s\n", filename);
		exit ( nerr ) ;
	}

	getfhv("b", &b, &nerr, 1);
	n=floor(2*taperwidth/del);
	
	// Do the windowing
	for (i=0;i<nlen;i++)
	{
		if ( i*del+b < begin || i*del+b > end )
			data[i]=0;
		else if (i*del+b < begin + taperwidth)
			data[i]=data[i]*(0.5-0.5*cos(2*PI*(floor((i*del+b-begin)/del))/(n-1)));
		else if (i*del+b > end - taperwidth)
			data[i]=data[i]*(0.5-0.5*cos(2*PI*(floor((i*del+b-(end-taperwidth)+taperwidth)/del))/(n-1)));
		else
			data[i]=data[i];
	}
	sprintf(newfilename, "%s.w", filename);

	wsac1( newfilename, data, &nlen, &beg, &del, &nerr, strlen(newfilename) ) ;

}

/* Find the nearest circle to the theoretical phase difference */
double FindNearestCircle( double dt, double tdt, double t){
	double ndt;
	double period;
	
	
	period=fabs(t);
	ndt=dt;
	while ( fabs(ndt-tdt) > period/2 )
	{
		if (ndt > tdt)
			ndt = ndt - period;
		else
			ndt = ndt + period;
	}
	return ndt;
}

/* This function is used to remove the effect of windowing and filtering of cross-correlation functions. The corrected group and phase differences are stored in variances taug and taup
 */
int CSRecover(CS *csp){
	int i;
	int n;
	char stemp[100];
	FILE *sacfp, *fitfp;
	float xi,wf;
	float maxt;
	float epidist1, epidist2;

	/*
	xi=csp->sta2->sc/csp->sta2->swc;
	for(i=0;i<PNUM;i++)
	{
		wf=csp->fit[i].wo;
		csp->taug[i]=(csp->fit[i].tg - csp->sta2->fit[i].tg \
				-(1-xi*xi)*csp->tc)/xi/xi;
		csp->taup[i]=csp->fit[i].tp \
			-FindNearestCircle(csp->sta2->fit[i].tp, 0, 2*PI/csp->sta2->fit[i].wo)\
			-(1-xi*xi)*(wf-csp->sta2->wc)/wf*(csp->tc - csp->taug[i]);
		csp->tprecover[i] = (1-xi*xi)*(wf-csp->sta2->wc)/wf*(csp->tc - csp->taug[i])\
			+FindNearestCircle(csp->sta2->fit[i].tp, 0, 2*PI/csp->sta2->fit[i].wo);
	}
	*/

	for(i=0;i<PNUM;i++)
		csp->taug[i] = csp->fit[i].tg - csp->sta2->fit[i].tg;

	for(i=0;i<PNUM;i++)
		for(n=0;n<ITN[i];n++)
		{
			sprintf(stemp,"%s/%s_%s.xcor",csp->ev->path, csp->sta1->staname, csp->sta2->staname);
			maxt = csp->taug[i];
			if (csp->fit[i].ierr!=1)
				break;
			
			epidist1=distance(csp->sta1->la, csp->sta1->lo, csp->sta1->ev->la, csp->sta1->ev->lo);
			epidist2=distance(csp->sta2->la, csp->sta2->lo, csp->sta2->ev->la, csp->sta2->ev->lo);
			if (maxt<(epidist1-epidist2)/grvavg-HALFWINDOW/2 || \
					maxt>(epidist1-epidist2)/grvavg+HALFWINDOW/2)
				break;
			
			if (ISMAKEFILES)
				HanWindow(stemp, maxt, HALFWINDOWXCOR);

			// Narrow band filtering
			sacfp=fopen("sacmacrotemp.csh","w");
			if (sacfp==NULL)
				printf("Error in reading sacmacrotemp.csh\n");
			fprintf(sacfp, "sac<<! \n");
			fprintf(sacfp, "r %s/%s_%s.xcor.w \n", csp->ev->path, \
					csp->sta1->staname, csp->sta2->staname);
			fprintf(sacfp, "bp n 4 co %f %f passes 2 \n", 1.0/periods[i]*0.9, 1.0/periods[i]*1.1);
			fprintf(sacfp, "w %s/%s_%s_%1d.%1d.xcor \n",csp->ev->path, \
					csp->sta1->staname, csp->sta2->staname, i, n);
			fprintf(sacfp, "r %s/%s_%s_%1d.%1d.xcor \n",csp->ev->path, \
					csp->sta1->staname, csp->sta2->staname, i, n);
			fprintf(sacfp, "cut %d %d \n", -6*periods[PNUM-1], 6*periods[PNUM-1] );
			fprintf(sacfp, "r \n");
			fprintf(sacfp, "w over \n");
			fprintf(sacfp, "cut off \n");
			fprintf(sacfp, "q \n");
			fprintf(sacfp, "! \n");
			fprintf(sacfp, "\n");
			fclose(sacfp);
			if (ISMAKEFILES)
				system("csh sacmacrotemp.csh");

			// Fit the waveform
			sprintf(stemp, "gsdf_fit %s/%s_%s_%1d.%1d.xcor %d %d %s/%s_%s_%1d.%1d.xcor.fit\n",\
					csp->ev->path, csp->sta1->staname, csp->sta2->staname, i, n, periods[i], NFIT,\
					csp->ev->path, csp->sta1->staname, csp->sta2->staname, i, n);
			if (ISMAKEFILES)
			{
				printf("%s\n",stemp);
				system(stemp);
			}
			sprintf(stemp, "%s/%s_%s_%1d.%1d.xcor.fit",\
					csp->ev->path, csp->sta1->staname, csp->sta2->staname, i,n);
			fitfp=fopen(stemp,"r");
			if (fitfp==NULL)
				printf("Error in reading %s\n",stemp);
			fscanf(fitfp, "%lf %lf %lf %lf %lf %lf %d", \
					&csp->fit[i].ao,&csp->fit[i].so,&csp->fit[i].wo, \
					&csp->fit[i].tp,&csp->fit[i].tg,&csp->fit[i].chi, &csp->fit[i].ierr);
			fclose(fitfp);

			csp->taug[i] = csp->fit[i].tg - csp->sta2->fit[i].tg;

		}
	
	for(i=0;i<PNUM;i++)
		csp->taup[i]=csp->fit[i].tp;
	// Test result shows that it's better not to make this correction.
	//- FindNearestCircle(csp->sta2->fit[i].tp, 0, 2*PI/csp->sta2->fit[i].wo);
	if (csp->pair != NULL)
		for(i=0;i<PNUM;i++)
		{
			csp->taup[i]=(csp->taup[i]- FindNearestCircle(csp->pair->taup[i],\
						-csp->taup[i], 2*PI/csp->pair->fit[i].wo))/2;
		}
}

/* This function is used to fill the measurement into different periods and circles.
 * and grade the refG and contG and the same time.
 */
int GradeLocal(CS *csp){
	char xcorfname[100],wxcorfname[100];
	double epidist1, epidist2;
	double tdt,dt;
	double a,b,c;
	double contdt;
	int i;
	int j;

	// Calculate the theoretical time difference
	epidist2=distance(csp->sta2->la, csp->sta2->lo, csp->sta2->ev->la, csp->sta2->ev->lo);
	epidist1=distance(csp->sta1->la, csp->sta1->lo, csp->sta1->ev->la, csp->sta1->ev->lo);
	tdt=(epidist1-epidist2)/GetPhV(periods[CENTP]);

	// Set the center period
	dt = FindNearestCircle(csp->taup[CENTP], tdt, 2*PI/csp->fit[CENTP].wo );
	csp->meas[CENTP][0].dt = dt - 2*PI/csp->fit[CENTP].wo;
	csp->meas[CENTP][1].dt = dt;
	csp->meas[CENTP][2].dt = dt + 2*PI/csp->fit[CENTP].wo;
	// Grade center period
	for (i=0;i<CIRCNUM;i++)
	{
		csp->meas[CENTP][i].refG=GaussianGrade(csp->meas[CENTP][i].dt, tdt, periods[CENTP]);
		csp->meas[CENTP][i].contG=0; 	// For center period, do not consider continuity.
		csp->meas[CENTP][i].localG = csp->meas[CENTP][i].refG;
		csp->meas[CENTP][i].pathG = 0;
		csp->meas[CENTP][i].totalG = csp->meas[CENTP][i].localG;
	}

	// Loop over the periods before center period
	for (i=CENTP-1;i>=0;i--)
	{
		tdt=(epidist1-epidist2)/GetPhV(periods[i]);
		contdt= GetBestdtCS(csp,i+1);
		//decide the interested circles by either tdt or contdt, decided by weight
		if (localGweight[i][0] > localGweight[i][1])
			dt = FindNearestCircle(csp->taup[i], tdt, 2*PI/csp->fit[i].wo );
		else
			dt = FindNearestCircle(csp->taup[i], contdt, 2*PI/csp->fit[i].wo );
		csp->meas[i][0].dt = dt - 2*PI/csp->fit[i].wo;
		csp->meas[i][1].dt = dt;
		csp->meas[i][2].dt = dt + 2*PI/csp->fit[i].wo;
		// Calculate the refG
		for (j=0;j<CIRCNUM;j++)
		{
			csp->meas[i][j].refG=GaussianGrade(csp->meas[i][j].dt, tdt, periods[i]);
		}
		// Calculate the contG
		for (j=0;j<CIRCNUM;j++)
		{
			csp->meas[i][j].contG=GaussianGrade(csp->meas[i][j].dt, contdt, periods[i]);
		}
		// Calculate the localG
		for (j=0;j<CIRCNUM;j++)
		{
			csp->meas[i][j].localG \
				= (localGweight[i][0]*csp->meas[i][j].refG \
				+ localGweight[i][1]*csp->meas[i][j].contG)\
				/(localGweight[i][0]+localGweight[i][1]);
			csp->meas[i][j].pathG = 0;
			csp->meas[i][j].totalG = csp->meas[i][j].localG;
		}
	}

	// Loop over the periods after center period
	for (i=CENTP+1;i<PNUM;i++)
	{
		tdt=(epidist1-epidist2)/GetPhV(periods[i]);
		contdt= GetBestdtCS(csp,i-1);
		//decide the interested circles by either tdt or contdt, decided by weight
		if (localGweight[i][0] > localGweight[i][1])
			dt = FindNearestCircle(csp->taup[i], tdt, 2*PI/csp->fit[i].wo );
		else
			dt = FindNearestCircle(csp->taup[i], contdt, 2*PI/csp->fit[i].wo );
		// remove the effect of windowing
		//dt = dt - FindNearestCircle(csp->sta2->fit[i].tp, 0, 2*PI/csp->sta2->fit[i].wo);
		csp->meas[i][0].dt = dt - 2*PI/csp->fit[i].wo;
		csp->meas[i][1].dt = dt;
		csp->meas[i][2].dt = dt + 2*PI/csp->fit[i].wo;
		// Calculate the refG
		for (j=0;j<CIRCNUM;j++)
		{
			csp->meas[i][j].refG=GaussianGrade(csp->meas[i][j].dt, tdt, periods[i]);
		}
		// Calculate the contG
		for (j=0;j<CIRCNUM;j++)
		{
			csp->meas[i][j].contG=GaussianGrade(csp->meas[i][j].dt, contdt, periods[i]);
		}
		// Calculate the localG
		for (j=0;j<CIRCNUM;j++)
		{
			csp->meas[i][j].localG \
				= (localGweight[i][0]*csp->meas[i][j].refG \
				+ localGweight[i][1]*csp->meas[i][j].contG)\
				/(localGweight[i][0]+localGweight[i][1]);
			csp->meas[i][j].totalG = csp->meas[i][j].localG;
		}
	}

	// Calculate the rate of energy that windowed cross-correlation to the all xcor
	sprintf(xcorfname,"%s/%s_%s.xcor",csp->ev->path,csp->sta1->staname,csp->sta2->staname);
	sprintf(wxcorfname,"%s/%s_%s.xcor.w",csp->ev->path,csp->sta1->staname,csp->sta2->staname);
	csp->snr=AvgEnergy(wxcorfname)/AvgEnergy(xcorfname);

	// Calculate the cohere of each periods
	for (i=0;i<PNUM;i++)
		csp->cohere[i] = csp->fit[i].ao*csp->fit[i].ao \
					     /csp->sta1->fit[i].ao/csp->sta2->fit[i].ao;

}

/* The function to decide whether a cross-correlation measurement is good or not */
int IsGoodData(CS *csp, int pid){
	double epidist1,epidist2;
	float tt;
	float ftemp;
	int j;

	if (csp->fit[pid].ierr != 1)
		return -1;

	if (csp->sta2->fit[pid].ierr != 1)
		return -2;

	if (fabs(2*PI/csp->fit[pid].wo - periods[pid]) > periods[pid]*0.1 )
		return -3;

	if (fabs(2*PI/csp->sta2->fit[pid].wo - periods[pid]) > periods[pid]*0.1 )
		return -4;

	if (csp->fit[pid].chi > MAXCHI)
		return -5;
	if (csp->sta2->fit[pid].chi > MAXCHI)
		return -6;

	if (fabs(csp->sta2->fit[pid].tp) > periods[pid]*0.1)
		return -7;
	
	if (csp->cohere[pid] < MINCOHERE)
		return -10;

	if (csp->snr < MINENERGYRATE)
		return -8;
	
	if (GetMinStaSNR(csp, pid) < 2)
		return -9;

	/*
	j=GetTotalBestCircCS(csp,pid);
	if (fabs(csp->epidiff/GetPhV(periods[pid]) - csp->meas[pid][j].dt) > MAXDTMISFIT )
		return -11;
		*/

	if (fabs(csp->dterracpt[pid])==0)
		return -12;

/*	if (csp->sta2->cfffit.ierr != 1)
		return -13;

	if (csp->sta2->wcfffit.ierr != 1)
		return -14;

	if (csp->sta2->sc/csp->sta2->swc > 1 || csp->sta2->sc/csp->sta2->swc < 0.6)
		return -15;

	ftemp=csp->sta2->cfffit.wo/csp->sta2->wcfffit.wo;
	if (ftemp > 1.1 || ftemp < 0.9)
		return -16;
*/

	/* Check other period using recursion method */
	if (pid == CENTP)
		return 1;
	else if (pid < CENTP)
		if (IsGoodData(csp, pid+1) == 1)
			return 1;
		else
			return -99;
	else
		return 1;
		/*
		if (IsGoodData(csp, pid-1) == 1)
			return 1;
		else
			return -99;
			*/

	return 1;
}

/* The function to make auto-correlation sac files for each station */
int MakeAcorFiles(STA *stap){
	char stemp[100];
	float epidist;
	float tt1,tt2;
	FILE *sacfp,*fitfp;
	int i;

	epidist=distance(stap->la, stap->lo, stap->ev->la, stap->ev->lo);
	tt1=stap->ev->t0 + epidist/stap->ev->grv0;
	tt2=stap->ev->t1 + epidist/stap->ev->grv1;
	//tt2=stap->ev->t2 + (epidist-stap->ev->epidist)/grvmin;
	//tt1=stap->ev->t1 + (epidist-stap->ev->epidist)/grvmax;

	if (ISMAKEFILES){
	// Window the waveform of station2 to get the isolate filter
	HanWindowFlat(stap->filename, (tt1-CUTBEFORE+tt2+CUTAFTER)/2, (tt2+CUTAFTER-tt1+CUTBEFORE)/2, TAPERWIDTH);

	sacfp=fopen("sacmacrotemp.csh","w");
	if (sacfp==NULL)
		printf("Error in reading sacmacrotemp.csh\n");
	// Read in waveforms and do the pre-filter and auto-correlation
	fprintf(sacfp, "sac<<! \n");
	fprintf(sacfp, "r %s %s.w \n", stap->filename,stap->filename);
	if (ISPREFILTER)
		fprintf(sacfp, "bp n 4 co %f %f passes 1 \n", 1/prefilter[1], 1/prefilter[0]);
	fprintf(sacfp, "correlate master 2 \n");
	fprintf(sacfp, "w %s.acor %s.facor \n", stap->staname, stap->staname);
	fprintf(sacfp, "q \n");
	fprintf(sacfp, "! \n");
	fprintf(sacfp, "\n");
	fclose(sacfp);
	system("csh sacmacrotemp.csh");

	sprintf(stemp,"%s.acor",stap->staname);
	HanWindow(stemp, 0, HALFWINDOWXCOR);

	sacfp=fopen("sacmacrotemp.csh","w");
	if (sacfp==NULL)
		printf("Error in reading sacmacrotemp.csh\n");
	fprintf(sacfp, "sac<<! \n");
	// Narrow filter the auto-correlation waveforms
	for (i=0;i<PNUM;i++)
	{
		fprintf(sacfp, "r %s.acor.w %s.facor \n", stap->staname, stap->staname);

		fprintf(sacfp, "bp n 4 co %f %f passes 2 \n", 1.0/periods[i]*0.9, 1.0/periods[i]*1.1);
		fprintf(sacfp, "w %s_%1d.acor %s_%1d.facor \n",\
				stap->staname, i,stap->staname,i);
		fprintf(sacfp, "r %s_%1d.acor %s_%1d.facor \n",\
				stap->staname, i, stap->staname, i);
		fprintf(sacfp, "cut %d %d \n", -6*periods[PNUM-1], 6*periods[PNUM-1] );
		fprintf(sacfp, "r \n");
		fprintf(sacfp, "w over \n");
		fprintf(sacfp, "cut off \n");
	}

	// End of sac script file
	fprintf(sacfp, "mv *.acor %s/ \n", stap->ev->path);
	fprintf(sacfp, "mv *.facor %s/ \n", stap->ev->path);
	fprintf(sacfp, "mv *.w %s/ \n", stap->ev->path);
	fprintf(sacfp, "q \n");
	fprintf(sacfp, "! \n");
	fprintf(sacfp, "\n");
	fclose(sacfp);
	system("csh sacmacrotemp.csh");

	}
	// Fit the auto-correlation curve and get the phase delay 
	for (i=0;i<PNUM;i++)	// Loop over the frequency
	{
		sprintf(stemp, "gsdf_fit %s/%s_%1d.acor %d %d %s/%s_%1d.acor.fit\n",\
				stap->ev->path, stap->staname, i, periods[i], NFIT,\
				stap->ev->path, stap->staname, i);
		if (ISMAKEFILES)
		{
			printf("%s\n",stemp);
			system(stemp);
		}

		sprintf(stemp, "%s/%s_%1d.acor.fit",\
				stap->ev->path, stap->staname, i);
		fitfp=fopen(stemp,"r");
		if (fitfp==NULL)
			printf("Error in reading %s\n",stemp);
		fscanf(fitfp, "%lf %lf %lf %lf %lf %lf %d", \
				&stap->fit[i].ao,&stap->fit[i].so,&stap->fit[i].wo, \
				&stap->fit[i].tp,&stap->fit[i].tg,&stap->fit[i].chi, &stap->fit[i].ierr);
		fclose(fitfp);
	}
	/* old version, kick it out
	// Fit cff to get \sigma_c, \sigma_wc, \omega_c, then use them to exclude the windowing and filter effect
	sprintf(stemp, "gsdf_fit %s/%s.facor %d %d %s/%s.facor.fit\n",\
			stap->ev->path, stap->staname, 40, NFIT,\
			stap->ev->path, stap->staname);
	if (ISMAKEFILES)
	{
		printf("%s\n",stemp);
		system(stemp);
	}

	sprintf(stemp, "%s/%s.facor.fit",\
			stap->ev->path, stap->staname);
	fitfp=fopen(stemp,"r");
	fscanf(fitfp, "%lf %lf %lf %lf %lf %lf %d", \
			&stap->cfffit.ao,&stap->cfffit.so,&stap->cfffit.wo, \
			&stap->cfffit.tp,&stap->cfffit.tg,&stap->cfffit.chi, &stap->cfffit.ierr);
	fclose(fitfp);
	stap->sc=stap->cfffit.so;
	stap->wc=stap->cfffit.wo;

	// Window the cff to get the swc
	sprintf(stemp, "%s/%s.facor",\
			stap->ev->path, stap->staname);
	HanWindow(stemp, 0, HALFWINDOWXCOR);
	sprintf(stemp, "gsdf_fit %s/%s.facor.w %d %d %s/%s.facor.w.fit\n",\
			stap->ev->path, stap->staname, 40, NFIT,\
			stap->ev->path, stap->staname);
	if (ISMAKEFILES)
	{
		printf("%s\n",stemp);
		system(stemp);
	}

	sprintf(stemp, "%s/%s.facor.w.fit",\
			stap->ev->path, stap->staname);
	fitfp=fopen(stemp,"r");
	fscanf(fitfp, "%lf %lf %lf %lf %lf %lf %d", \
			&stap->wcfffit.ao,&stap->wcfffit.so,&stap->wcfffit.wo, \
			&stap->wcfffit.tp,&stap->wcfffit.tg,&stap->wcfffit.chi, &stap->wcfffit.ierr);
	fclose(fitfp);
	stap->swc=stap->wcfffit.so;
	*/

}

/* Function to fill costa structure with real data
 * This is one of the most important function in this program. 
 * It reads in the original waveform from sac files,
 * window the two waveforms in a large range, 
 * cut the interested waveform part of the second station data as isolated filter
 * doing cross-correlation of the filter and the waveform from both stations
 * and remove the effect of windowing. 
 * Windowed the cross-correlaiton waveform,
 * narrow band filter the cross-correlation waveform
 * fit the curve
 * grading the measurement
 */
void FillCostaData(CS *csp){
	char commandline[100];
	char stemp[100];
	double epidist2, epidist1;
	int nlen, nerr, max = MAXDP;	// for sac file reading
	int i,j;
	float beg, del;
	float b,e,o;
	float data[MAXDP];
	float dist;
	float tt1,tt2;
	float maxt;
	FILE *sacfp, *fitfp, *tempfp;

	// Build up csh shell to window, filter and cross-correlation using sac
	epidist2=distance(csp->sta2->la, csp->sta2->lo, csp->sta2->ev->la, csp->sta2->ev->lo);
	epidist1=distance(csp->sta1->la, csp->sta1->lo, csp->sta1->ev->la, csp->sta1->ev->lo);

	if (ISMAKEFILES)
	{

	printf("Making %s_%s.xcor\n",csp->sta1->staname, csp->sta2->staname);
	// Create a temporary file to run sac
	sacfp=fopen("sacmacrotemp.csh","w");
	if (sacfp==NULL)
		printf("Error in reading sacmacrotemp.csh\n");

	// Read in all waveforms and do the pre-filter and cross-correlation
	fprintf(sacfp, "sac<<! \n");
	fprintf(sacfp, "r %s %s.w \n", \
				csp->sta1->filename,csp->sta2->filename);
	if (ISPREFILTER)
		fprintf(sacfp, "bp n 4 co %f %f passes 2 \n", 1/prefilter[1], 1/prefilter[0]);
	fprintf(sacfp, "correlate master 2 \n");
	fprintf(sacfp, "w %s_%s.xcor %s.facor \n",\
				csp->sta1->staname, csp->sta2->staname, csp->sta2->staname);

	fprintf(sacfp, "q \n");
	fprintf(sacfp, "! \n");
	fprintf(sacfp, "\n");
	fclose(sacfp);

	system("csh sacmacrotemp.csh");
	}

	if (ISMAKEFILES)
		sprintf(stemp,"%s_%s.xcor",csp->sta1->staname, csp->sta2->staname);
	else
		sprintf(stemp,"%s/%s_%s.xcor",csp->ev->path, csp->sta1->staname, csp->sta2->staname);

	if (ISWINDOWMAX)
		maxt = FindXcorPeak(stemp, (epidist1-epidist2)/grvavg, HALFWINDOW);
	else 
		maxt = (epidist1-epidist2)/grvavg;
	csp->tc=maxt;

	if (ISMAKEFILES)
	{
	HanWindow(stemp, maxt, HALFWINDOWXCOR);

	// Narrow band filtering
	sacfp=fopen("sacmacrotemp.csh","w");
	if (sacfp==NULL)
		printf("Error in reading sacmacrotemp.csh\n");

	// Window the waveform of station2 to get the isolate filter
	fprintf(sacfp, "sac<<! \n");
	// Narrow filter the correlation waveforms
	for (i=0;i<PNUM;i++)
	{
		fprintf(sacfp, "r %s_%s.xcor.w \n",\
				csp->sta1->staname, csp->sta2->staname);
		fprintf(sacfp, "bp n 4 co %f %f passes 2 \n", 1.0/periods[i]*0.9, 1.0/periods[i]*1.1);
		fprintf(sacfp, "w %s_%s_%1d.xcor \n",\
				csp->sta1->staname, csp->sta2->staname, i);
		fprintf(sacfp, "r %s_%s_%1d.xcor \n",\
				csp->sta1->staname, csp->sta2->staname, i);
		fprintf(sacfp, "cut %d %d \n", -6*periods[PNUM-1], 6*periods[PNUM-1] );
		fprintf(sacfp, "r \n");
		fprintf(sacfp, "w over \n");
		fprintf(sacfp, "cut off \n");
	}

	// End of sac script file
	fprintf(sacfp, "mv *.xcor %s/ \n", csp->ev->path);
	fprintf(sacfp, "mv *.facor %s/ \n", csp->ev->path);
	fprintf(sacfp, "mv *.w %s/ \n", csp->ev->path);
	fprintf(sacfp, "q \n");
	fprintf(sacfp, "! \n");
	fprintf(sacfp, "\n");
	fclose(sacfp);

	system("csh sacmacrotemp.csh");
	}


	// Now fit the cross-correlation curve and get the phase delay 
	for (i=0;i<PNUM;i++)	// Loop over the frequency
	{
		sprintf(stemp, "gsdf_fit %s/%s_%s_%1d.xcor %d %d %s/%s_%s_%1d.xcor.fit\n",\
				csp->ev->path, csp->sta1->staname, csp->sta2->staname, i, periods[i], NFIT,\
				csp->ev->path, csp->sta1->staname, csp->sta2->staname, i);
		if (ISMAKEFILES)
		{
			printf("%s\n",stemp);
			system(stemp);
		}
		sprintf(stemp, "%s/%s_%s_%1d.xcor.fit",\
				csp->ev->path, csp->sta1->staname, csp->sta2->staname, i);
		fitfp=fopen(stemp,"r");
		if (fitfp==NULL)
			printf("Error in reading %s\n",stemp);
		fscanf(fitfp, "%lf %lf %lf %lf %lf %lf %d", \
				&csp->fit[i].ao,&csp->fit[i].so,&csp->fit[i].wo, \
				&csp->fit[i].tp,&csp->fit[i].tg,&csp->fit[i].chi, &csp->fit[i].ierr);
		fclose(fitfp);
	}


	/* Print out and check
	for (i=0;i<PNUM;i++)	// Loop over the frequency
	{
		printf("%e %lf %lf %lf %lf %lf %d\n", \
				csp->meas[i][0].ao,csp->meas[i][0].so,csp->meas[i][0].wo, \
				csp->meas[i][0].tp,csp->meas[i][0].tg,csp->meas[i][0].chi, csp->meas[i][0].ierr);
	}
	*/


	// Output the result to check
/*	
	sprintf(stemp,"%s.out",csp->ev->path);
	tempfp=fopen("tempoutput.txt","a+");
	fprintf(tempfp,"\n");
	//fprintf(tempfp, "%s_%s\n",csp->sta1->staname,csp->sta2->staname);
	for (i=0;i<PNUM;i++)
	{
		if ( IsGoodData(csp,i)!=1 )
		{
			//fprintf(tempfp,"Period %d is skipped because of %d \n",periods[i], IsGoodData(csp,i));
		}
		else {
			tt1=(epidist1-epidist2)/GetPhV(periods[i]);
			fprintf(tempfp," %d %f ",periods[i], tt1);
			fprintf(tempfp,"%5.1f %3.1f ", \
					GetBestdtCS(csp,i),csp->meas[i][GetBestCircCS(csp,i)].localG);
			fprintf(tempfp,"\n");
		}
			
	}
	fclose(tempfp);
	*/
}

/* Function to evaluate the error level based on the agreement of common friends measurement */
int ErrLevel(CS *cshead){

	CS *csp;
	int errlevel,i;
	double err;

	csp=cshead;
	while(csp!=NULL)
	{
		for (i=0;i<PNUM;i++)
		{
			if (csp->cfnum[i] == 0)
				csp->errlevel[i]=4;
			else
			{
				err=fabs(GetTotalBestdtCS(csp,i)-csp->cfdt[i]);
				if (err > 0.3)
					csp->errlevel[i]=4;
				else if (err > 0.2)
					csp->errlevel[i]=3;
				else if (err > 0.1)
					csp->errlevel[i]=2;
				else 
					csp->errlevel[i]=1;
			}
		}
		csp=csp->n;
	}
}

/* Function to find the best average phase velocity using least square method. Another program
 * called fitphasev is used.
 */
float FitPhaseVCS(CS *cshead, int pid){
	FILE *fp;
	CS *csp;
	char stemp[100];
	int csnum=0;
	int j;
	float phasev, cov, sumsq;

	// Ready the input file for program "fitphasev"
	fp=fopen("fitphasev.in","w");
	if (fp==NULL)
		printf("Error in reading fitphasev.in\n");
	csp=cshead;
	while(csp!=NULL)
	{
		if (IsGoodData(csp,pid) > 0)
		{
			j=GetTotalBestCircCS(csp,pid);
			fprintf(fp,"%f %f\n",csp->epidiff,csp->meas[pid][j].dt);
			csnum++;
		}
		csp=csp->n;
	}
	fclose(fp);

	// Run the program
	sprintf(stemp,"fitphasev fitphasev.in\n");
	system(stemp);

	// Read in the output file
	fp=fopen("fitphasev.out","r");
	if (fp==NULL)
		printf("Error in reading fitphasev.out\n");
	fgets(stemp,100,fp);
	sscanf(stemp,"%e %e %e\n", &phasev, &cov, &sumsq);
	phasev=1/phasev;
	fclose(fp);

	return phasev;

}

/* Function to find the best average phase velocity using least square method. Another program
 * called fitphasev is used.
 */
float FitPhaseVCSCH(CSCH *cschhead, int pid){
	FILE *fp;
	CSCH *cschp;
	CS *csp;
	char stemp[100];
	int csnum=0;
	int j;
	float phasev, cov, sumsq;

	// Ready the input file for program "fitphasev"
	fp=fopen("fitphasev.in","w");
	if (fp==NULL)
		printf("Error in reading fitphasev.in\n");
	cschp=cschhead;
	while(cschp!=NULL)
	{
		csp=cschp->cs;
		if (IsGoodData(csp,pid) > 0)
		{
			j=GetTotalBestCircCS(csp,pid);
			fprintf(fp,"%f %f\n",csp->epidiff,csp->meas[pid][j].dt);
			csnum++;
		}
		cschp=cschp->n;
	}
	fclose(fp);

	// Run the program
	sprintf(stemp,"fitphasev fitphasev.in\n");
	system(stemp);

	// Read in the output file
	fp=fopen("fitphasev.out","r");
	if (fp==NULL)
		printf("Error in reading fitphasev.out\n");
	fgets(stemp,100,fp);
	sscanf(stemp,"%e %e %e\n", &phasev, &cov, &sumsq);
	phasev=1/phasev;
	fclose(fp);

	return phasev;

}

/* Function to find the best average phase velocity based on the measurement of 
 * one station. 
 * program called fitphasev2 is used.
 * This program calculate the slope as well as the offset, which is considered
 * as one of the standard to see whether the station is good or not.
 */
int FitPhaseVSTA(STA *stap, int pid){
	FILE *fp;
	CSCH *cschp;
	CS *csp;
	char stemp[100],stemp2[100];
	int csnum=0;
	int j;
	float phasev, cov, sumsq;
	float r;

	// Ready the input file for program "fitphasev"
	sprintf(stemp,"%s/%s_%1d.phv",stap->ev->path,stap->staname,pid);
	fp=fopen(stemp,"w");
	if (fp==NULL)
		printf("Error in reading %s\n",stemp);
	cschp=stap->csch;
	while(cschp!=NULL)
	{
		csp=cschp->cs;
		if (IsGoodData(csp,pid) > 0)
		{
			j=GetTotalBestCircCS(csp,pid);
			fprintf(fp,"%f %f\n",csp->epidiff,csp->meas[pid][j].dt);
			csnum++;
		}
		cschp=cschp->n;
	}
	fclose(fp);

	// Run the program
	sprintf(stemp2,"fitphasev2 %s\n",stemp);
	system(stemp2);

	// Read in the output file
	fp=fopen("fitphasev2.out","r");
	if (fp==NULL)
		printf("Error in reading fitphasev2.out");
	fgets(stemp,100,fp);
	sscanf(stemp,"%e %e %e %e %e %e\n", \
			&stap->gsl[pid].c0,&stap->gsl[pid].c1,\
			&stap->gsl[pid].cov00, &stap->gsl[pid].cov01,&stap->gsl[pid].cov11,\
			&stap->gsl[pid].sumsq);
	fclose(fp);
	stap->avgphv[pid]=1/stap->gsl[pid].c1;
	stap->offset[pid]=stap->gsl[pid].c0;
	// Calculate the uncertainty of the slope
	r = sqrt(stap->gsl[pid].cov11)/stap->gsl[pid].c1;
	stap->avgphvdv[pid]=stap->avgphv[pid]*r;
	return 1;
}

/* Function to calculate the total err of a phase velocity */
float GetPhvErr(float phv, CS *cshead, int pid){
	CS *csp;
	float err=0;
	int csnum=0;
	int j;

	csp=cshead;
	while(csp!=NULL)
	{
		if (IsGoodData(csp,pid) > 0)
		{
			j=GetTotalBestCircCS(csp,pid);
			err+=fabs(csp->meas[pid][j].dt - csp->epidiff/phv);
			csnum++;
		}
		csp=csp->n;
	}
	
	if (csnum!=0)
		err=err/csnum;
	else
		err=0;
	return err;
}

/* Function to calculate a best fit phase velocity and exclude the outliers. */
int RefModGrading(CS *cshead){
	CS *csp;
	float uv,dv,mv;
	float apv;
	float uerr,derr;
	float dt,dtl,dth;
	int i,j;

	// Find the best fitted phase velocity for each frequency.
	for (i=0;i<PNUM;i++)
	{
		mv=GetPhV(periods[i]);
		apv=FitPhaseVCS(cshead,i);
		if (apv < mv*1.1 && apv > mv*0.9)
			cshead->ev->avgphv[i]=apv;
		else
			cshead->ev->avgphv[i]=mv;
	}
	/*   This is the old version, forget it!
	for (i=0;i<PNUM;i++)
	{
		mv=GetPhV(periods[i]);
		uv=mv*1.1;
		dv=mv*0.9;
		uerr=GetPhvErr(uv,cshead,i);
		if (uerr==0)
		{
			cshead->ev->avgphv[i]=0;
			continue;
		}
		while (uv-dv > 1e-2)
		{
			mv=(uv+dv)/2;
			uerr=GetPhvErr(uv,cshead,i);
			derr=GetPhvErr(dv,cshead,i);
			if (uerr > derr)
				uv=mv;
			else
				dv=mv;
		}
		cshead->ev->avgphv[i]=mv;
	}
	*/

	// Check whether this cs is a outlier, depended on +-10% velocity + 5s
	for(i=0;i<PNUM;i++)
	{
		csp=cshead;
		while(csp!=NULL)
		{
			if (IsGoodData(csp,i) > 0)
			{
				j=GetTotalBestCircCS(csp,i);
				if (csp->epidiff>0)
				{
					dtl=csp->epidiff/csp->ev->avgphv[i]/1.1-5;
					dth=csp->epidiff/csp->ev->avgphv[i]/0.9+5;
				}
				else
				{
					dtl=csp->epidiff/csp->ev->avgphv[i]/0.9-5;
					dth=csp->epidiff/csp->ev->avgphv[i]/1.1+5;
				}
				dt=csp->meas[i][j].dt;
				if (dt < dtl || dt > dth)
					csp->dterracpt[i]=0;
			}
			csp=csp->n;
		}
	}

}

/* Function to grade each station based on their data quality. */
int StationGrading(EV *evp){
	STA *stap;
	CS *csp;
	CSCH *cschp;
	int csnum=0, goodnum[PNUM];
	int i;

	stap=evp->stachain;
	
	// Grade by percent of good cs measurement
	while(stap!=NULL)
	{
		for (i=0;i<PNUM;i++)
			goodnum[i]=0;
		csnum=0;

		cschp=stap->csch;
		while (cschp!=NULL)
		{
			csnum++;
			for (i=0;i<PNUM;i++)
				if (IsGoodData(cschp->cs,i)>0 )
					goodnum[i]=goodnum[i]+1;
			cschp=cschp->n;
		}

		for (i=0;i<PNUM;i++)
			if (csnum!=0)
				stap->stagrade[i]=goodnum[i]*100.0/csnum;
			else
				stap->stagrade[i]=0;

		// Calculate the station average phase velocity
		for (i=0;i<PNUM;i++)
		{
			stap->goodnum[i]=goodnum[i];
			if ( goodnum[i]>FITPHASEVMINNUM )
			{
				FitPhaseVSTA(stap,i);
				if (stap->avgphv[i]<GetPhV(periods[i])*1.1 && \
						stap->avgphv[i]>GetPhV(periods[i])*0.9)
					stap->gsl[i].Isgood=1;
				else
					stap->gsl[i].Isgood=0;

			}
			else
			{
				stap->avgphv[i]=GetPhV(periods[i]);
				stap->gsl[i].Isgood=0;
			}
		}
		stap=stap->n;
	}

}

/* Function to calculate station's snr on each period */
int SetStationSNR(STA *stap){
	char filename[100];
	int nlen, nerr, max = MAXDP;		// for reading the sac
	FILE *sacmacrofp;
	float beg, del;
	float epidist;
	float data[MAXDP];
	float t0, t1;
	double sum;
	float signalavg, noiseavg;
	float tt1,tt2;
	float b,e;
	int i,j,k;
	int num;

	// Narrow band the original data and generate sacfiles
	
	sacmacrofp=fopen("sacmacrotemp.csh","w");
	if (sacmacrofp==NULL)
		printf("Error in reading %s\n","sacmacrotemp.csh");
	fprintf(sacmacrofp, "sac <<!\n");
	for (i=0;i<PNUM;i++){
		fprintf(sacmacrofp, "r %s\n",stap->filename);
		fprintf(sacmacrofp, "bp n 4 co %f %f passes 1 \n", 1.0/periods[i]*0.9, 1.0/periods[i]*1.1);
		fprintf(sacmacrofp, "w %s_%d.bp\n",stap->filename,i);
	}
	fprintf(sacmacrofp, "q\n");
	fprintf(sacmacrofp, "!\n");
	fclose(sacmacrofp);

	if (ISMAKEFILES)
		system("csh sacmacrotemp.csh\n");

	// Calculate the snr for each period
	epidist=distance(stap->la, stap->lo, stap->ev->la, stap->ev->lo);
	tt1=stap->ev->t0 + epidist/stap->ev->grv0;
	tt2=stap->ev->t1 + epidist/stap->ev->grv1;

	for (i=0;i<PNUM;i++){
		sprintf(filename,"%s_%d.bp",stap->filename,i);
		// Read in sac file
		rsac1( filename, data, &nlen, &beg, &del, &max, &nerr, strlen(filename) ) ;
		if ( nerr != 0 ) {
			fprintf(stderr, "Error reading in SAC file: %s\n", filename);
			exit ( nerr ) ;
		}

		getfhv("b", &b, &nerr, 1);
		// Calculate the signal average energy
		sum=0; num=0;
		for (j=floor((tt1-b)/del);j<floor((tt2-b)/del);j++){
			sum+=data[j]*data[j];
			num++;
		}
		signalavg=sum/num;
		//Calculate the noise average energy
		sum=0; num=0;
		for (j=floor((tt2-b)/del);j<nlen;j++){
			sum+=data[j]*data[j];
			num++;
		}
		noiseavg=sum/num;
		stap->snr[i]=signalavg/noiseavg;
	}
	return 1;
}

/* Function to get the minimum station snr */
float GetMinStaSNR(CS *csp,int i){
	if (csp->sta1->snr[i] < csp->sta2->snr[i])
		return csp->sta1->snr[i];
	else
		return csp->sta2->snr[i];
}

/* Function to get all the stations' relative phase difference to one station */
int GetStaDtNet(STA *sta0, STA *stahead, int pid){
	STA *statemptail,*stap, *staq;
	CS *csp;
	int addnum=1;
	int netnum=0;

	// Clean the mark
	stap=stahead;
	while(stap!=NULL){
		stap->isinchain=0;
		stap->dt=0;
		stap->amp=1;
		stap=stap->n;
	}

	sta0->dt=0;
	sta0->isinchain=1;
	sta0->tempn=NULL;
	statemptail=sta0;

	// Loop
	while (addnum!=0)
	{
		addnum=0;
		stap=stahead;
		while (stap!=NULL)
		{
			if (stap->isinchain==0)
			{
				staq=sta0;
				while (staq!=NULL)
				{
					csp=FindStaCo(stap,staq);
					if ( csp!=NULL && IsGoodData(csp,pid)>0 && csp->errlevel[pid]<4 )
					{
						statemptail->tempn=stap;
						statemptail=stap;
						stap->isinchain=1;
						stap->tempn=NULL;
						stap->dt = staq->dt + GetTotalBestdtSTA(stap,staq,pid);
						addnum++;
						netnum++;
					}
					staq=staq->tempn;
				}
			}
			stap=stap->n;
		}
	}
	return netnum;
}

/* Function to summerize each stations information from all events */
int StaInformGather(EV *evhead){
	STAINF *stainfhead, *stainftail;
	STAINF *stainfp;
	STA *stap, *staq;
	STA *statemptail;
	EV *evp,*evq;
	int i;
	int stanum;
	FILE *fp;
	char stemp[100];

	// initial everything
	evp=evhead;
	stainfhead=NULL;
	stainftail=NULL;
	stanum=0;

	// Loop of all events to find new stations and put them in the chain.
	while(evp!=NULL)
	{
		// stap is for find first station to be added into the station information chain
		stap=evp->stachain;
		while(stap!=NULL)
		{
			if (IsExistStaStaInform(stainfhead, stap)!=1)
			{
				stainfp=(STAINF *)malloc(sizeof(STAINF));
				// Fill in information
				stainfp->la=stap->la;
				stainfp->lo=stap->lo;
				strcpy(stainfp->staname, stap->staname);
				// Put the station information in chain
				stainfp->n=NULL;
				if (stainfhead==NULL)
				{
					stainfhead=stainfp;
					stainftail=stainfp;
				}
				else
				{
					stainftail->n=stainfp;
					stainftail=stainfp;
				}
				// Find the same station in other events and put them together using tempn
				stainfp->stachain=stap;
				stap->tempn=NULL;
				statemptail=stap;
				evq=evp->n;  // evq is for all the events after evp
				while(evq!=NULL)
				{
					staq=evq->stachain;  //staq is used to find same station
					while(staq!=NULL)
					{
						if(StationDistance(stainfp->stachain,staq) < STA_TOL_DIST)
						{
							statemptail->tempn=staq;
							statemptail=staq;
							staq->tempn=NULL;
						}
						staq=staq->n;
					}
					evq=evq->n;
				}
			}
			stap=stap->n;
		}	// End of station loop
		evp=evp->n;
	}	// End of event loop

	// Output everything into the ./stainform direction
	system("mkdir stainform\n");
	stainfp=stainfhead;
	while(stainfp!=NULL)
	{
		for(i=0;i<PNUM;i++)
		{
			sprintf(stemp,"stainform/%s_%1d.phv",stainfp->staname,i);
			fp=fopen(stemp,"w");
			if (fp==NULL)
				printf("Error in reading %s\n",stemp);
			if (fp==NULL)
			{
				printf("Cannot open %s\n",stemp);
				exit(0);
			}
			stap=stainfp->stachain;
			while(stap!=NULL)
			{
				fprintf(fp,"%f %f %f %f ", \
						stap->ev->la, stap->ev->lo, \
						stap->la, stap->lo);
				fprintf(fp,"%f %f ",stap->avgphv[i], stap->avgphvdv[i]);
				fprintf(fp,"%f ", stap->offset[i]);
				fprintf(fp,"%f %d ",stap->gsl[i].sumsq, stap->goodnum[i]);
				fprintf(fp,"%d ",stap->gsl[i].Isgood);
				fprintf(fp,"\n");
				stap=stap->tempn;
			}
			fclose(fp);
		}
		stainfp=stainfp->n;
	}
	
}

/* Function to format print csp data into txt files */
int FormatPrintCS(FILE *fp, CS *csp, int i){
	double epidist1,epidist2;
	float tt;
	int j;

	epidist2=distance(csp->sta2->la, csp->sta2->lo, csp->sta2->ev->la, csp->sta2->ev->lo);
	epidist1=distance(csp->sta1->la, csp->sta1->lo, csp->sta1->ev->la, csp->sta1->ev->lo);
	tt=(epidist1-epidist2)/GetPhV(periods[i]);
	//tt=csp->epidiff/csp->ev->avgphv[i];
	j=GetTotalBestCircCS(csp,i);
	fprintf(fp,"%s %s ",csp->sta1->staname, csp->sta2->staname);
	fprintf(fp,"%d %d ",csp->sta1->id, csp->sta2->id);
	fprintf(fp,"%8.5f %10.5f %8.5f %10.5f ",csp->sta1->la,csp->sta1->lo,csp->sta2->la,csp->sta2->lo);
	fprintf(fp,"%6.1f ",epidist1-epidist2);
	fprintf(fp,"%3d %f %3.1f %3.1f %3.1f %3.1f %1f ",periods[i], \
				csp->meas[i][j].dt,csp->meas[i][j].refG, \
				csp->meas[i][j].contG,csp->meas[i][j].pathG, \
				csp->snr,csp->errlevel[i]);
	fprintf(fp,"%f ",GetMinStaSNR(csp,i));
	fprintf(fp,"%f ",csp->tprecover[i]);
	fprintf(fp,"%f ",csp->cohere[i]);
	fprintf(fp,"%d ",IsGoodData(csp,i));
	fprintf(fp,"\n");

}

/* Function to write the output files */
int PrintOutputEV(EV *evp){
	char stemp[100];
	FILE *fp;
	CS *csp;
	CSCH *cschp;
	STA *beststap[PNUM];
	int maxnetnum[PNUM],netnum;
	int skipn,csnum,correctn,badn;
	int pcsnum[PNUM];
	int i,j,k;
	int stanum=0;
	int goodnum;
	double epidist1, epidist2;
	double grade[3];
	double tt;
	float mv,apv;
	float avgdt;
	float **dtmat[PNUM];
	STA *stp,*stq;

	if (evp!=NULL)
	{

		// Output station list
		sprintf(stemp,"%s.sta",evp->path);
		fp=fopen(stemp,"w");
		if (fp==NULL)
			printf("Error in reading %s\n",stemp);
		stp=evp->stachain;
		while(stp!=NULL){
			fprintf(fp,"%s %d %f %f ", stp->staname, stp->id, stp->la, stp->lo);
			for (i=0;i<PNUM;i++)
				fprintf(fp,"%e ",stp->fit[i].ao/stp->fit[i].so);
			fprintf(fp,"\n");
			stanum++;
			stp=stp->n;
		}
		fclose(fp);

		sprintf(stemp,"%s.stainv",evp->path);
		fp=fopen(stemp,"w");
		if (fp==NULL)
			printf("Error in reading %s\n",stemp);
		stp=evp->stachain;
		while(stp!=NULL){
			fprintf(fp,"%d %f %f ", stp->id, stp->la, stp->lo);
			for (i=0;i<PNUM;i++)
				fprintf(fp,"%e ",stp->fit[i].ao/stp->fit[i].so);
			fprintf(fp,"\n");
			stp=stp->n;
		}
		fclose(fp);

		// Output station average phase velocity
		stp=evp->stachain;
		while(stp!=NULL)
		{
			sprintf(stemp,"%s/%s.phv",evp->path, stp->staname,i);
			fp=fopen(stemp,"w");
			if (fp==NULL)
				printf("Error in reading %s\n",stemp);
			for (i=0;i<PNUM;i++)
			{
				fprintf(fp,"%f %f %f ",stp->avgphv[i],stp->avgphvdv[i],stp->offset[i]);
				fprintf(fp,"%f %d ",stp->gsl[i].sumsq,stp->goodnum[i]);
				fprintf(fp,"%d\n",stp->gsl[i].Isgood);
			}
			fclose(fp);
			stp=stp->n;
		}

		/*
	// Output station phase net
		
		// Build up huge matrix to record the station dt information
		for (i=0;i<PNUM;i++)
		{
			dtmat[i]=(float**)malloc(stanum*sizeof(float*));
			for (j=0;j<stanum;j++)
				dtmat[i][j]=(float*)malloc(stanum*sizeof(float));
		}
		//init the matrix
		for (i=0;i<PNUM;i++)
			for (j=0;j<PNUM;j++)
				for (k=0;k<PNUM;k++)
					dtmat[i][j][k]=0;

		printf("Calculating the station nets\n");
		stp=evp->stachain;
		for (i=0;i<PNUM;i++)
			maxnetnum[i]=0;
		while(stp!=NULL){
			for (i=0;i<PNUM;i++)
			{
				sprintf(stemp,"%s/%s_%1d.stadt",evp->path, stp->staname,i);
				fp=fopen(stemp,"w");
				if (fp==NULL)
					printf("Error in reading %s\n",stemp);
				netnum=GetStaDtNet(stp, evp->stachain, i);
				if (maxnetnum[i]<=netnum)
				{
					maxnetnum[i]=netnum;
					beststap[i]=stp;
				}
				stq=stp;
				while (stq!=NULL)
				{
					fprintf(fp,"%f %f %f %f\n", stq->la, stq->lo, stq->dt, stq->epidist);
					dtmat[i][stp->id][stq->id]=stq->dt;
					stq=stq->tempn;
				}
				fclose(fp);
			}
			stp=stp->n;
		}

		// output the best station net
		printf("Output the best station net......");
		for (i=0;i<PNUM;i++)
		{
			sprintf(stemp,"%s/best_%1d.stadt",evp->path,i);
			fp=fopen(stemp,"w");
			if (fp==NULL)
				printf("Error in reading %s\n",stemp);
			//netnum=GetStaDtNet(beststap[i], evp->stachain, i);
			stq=beststap[i];
			stp=evp->stachain;
			// average the nets from different stations

			while (stp!=NULL)
			{
				goodnum=0;
				avgdt=0;
				for (j=0;j<stanum;j++)
					if (dtmat[i][j][stp->id]!=0 || j==stp->id)
					if (dtmat[i][j][stq->id]!=0 || j==stq->id)
					{
						avgdt=avgdt+dtmat[i][j][stp->id]-dtmat[i][j][stq->id];
						goodnum++;
					}
				if (goodnum>0)
				{
					avgdt=avgdt/goodnum;
					fprintf(fp,"%f %f %f %f\n", stp->la, stp->lo, avgdt, stp->epidist);
				}
				stp=stp->n;
			}

			while (stq!=NULL)
			{
				fprintf(fp,"%f %f %f %f\n", stq->la, stq->lo, stq->dt, stq->epidist);
				stq=stq->tempn;
			}
			fclose(fp);
		}
		*/

		// Generate CS file for the inverse method.
		printf("Generate csinv files\n");
		for (i=0;i<PNUM;i++)
		{
			// Clear the isinchain mark
			csp=evp->cschain;
			while(csp!=NULL)
			{
				csp->isinchain=0;
				csp=csp->n;
			}
			// Get the best station net and output
			sprintf(stemp,"%s_%1d.csinv",evp->path,i);
			fp=fopen(stemp,"w");

		if (fp==NULL)
			printf("Error in reading %s\n",stemp);

			csp=evp->cschain;
			while(csp!=NULL)
			{
				j=GetTotalBestCircCS(csp,i);
				if (IsGoodData(csp, i)>0 )
					fprintf(fp,"%i %i %f\n", csp->sta1->id, csp->sta2->id, \
							csp->meas[i][j].dt);
				csp=csp->n;
			}

			/* Old version, now have matlab do the job.
			stp=beststap[i];
			netnum=GetStaDtNet(stp, evp->stachain, i);
			stq=stp;
			while (stq!=NULL)
			{
				cschp=stq->csch;
				while (cschp!=NULL)
				{
					j=GetTotalBestCircCS(cschp->cs,i);
					if (cschp->cs->isinchain == 0 && IsGoodData(cschp->cs, i)>0 )
						fprintf(fp,"%i %i %f\n", cschp->cs->sta1->id, cschp->cs->sta2->id, \
								cschp->cs->meas[i][j].dt);
					cschp->cs->isinchain=1;
					cschp=cschp->n;
				}
				stq=stq->tempn;
			} */
			fclose(fp);
		}

		printf("done!\n");

		//output log file
		sprintf(stemp,"%s.log",evp->path);
		fp=fopen(stemp,"w");
		if (fp==NULL)
			printf("Error in reading %s\n",stemp);
		fprintf(fp,"Event: %s %f %f \n", evp->path, evp->la, evp->lo);
		fprintf(fp,"Station Number: %d \n", evp->stationnum);
		// Count the number of circles that is corrected
		csp=evp->cschain;
		skipn=0;
		csnum=0;
		correctn=0;
		badn=0;
		for (i=0;i<PNUM;i++)
			pcsnum[i]=0;
		while (csp!=NULL)
		{
			for (i=0;i<PNUM;i++)
			{
				if (IsGoodData(csp,i)>0)
				{
					csnum++;
					pcsnum[i]=pcsnum[i]+1;
					for (j=0;j<3;j++)
						grade[j] = csp->meas[i][j].pathG;
					k = GetBestCircCS(csp,i);
					if (grade[k] < grade[(k+1)%3] || grade[k] < grade[(k+2)%3])
						skipn++;
					if (GetBestCircCS(csp,i)!=GetTotalBestCircCS(csp,i))
						correctn++;
				}
				else
					badn++;
			}
			csp=csp->n;
		}
		fprintf(fp,"Total Good measurement number is %d \n",csnum);
		fprintf(fp,"Bad measurement number is %d \n",badn);
		fprintf(fp,"circle skipping number based on pathG is %d.\n",skipn);
		fprintf(fp,"circle skipping number based on totalG is %d.\n",correctn);
		// Re-calculate the average phase velocity of each event
		for (i=0;i<PNUM;i++)
		{
			mv=GetPhV(periods[i]);
			apv=FitPhaseVCS(evp->cschain,i);
			if (apv < mv*1.1 && apv > mv*0.9)
				evp->avgphv[i]=apv;
			else
				evp->avgphv[i]=mv;
		}

		// Output the average phase velocity of each event
		for (i=0;i<PNUM;i++)
			fprintf(fp, "Period: %d, New PhaseV: %f, Goodnum: %d \n",\
					periods[i],evp->avgphv[i],pcsnum[i]); 
		fclose(fp);

		//output circle skip log file
		sprintf(stemp,"%s.circ",evp->path);
		fp=fopen(stemp,"w");
		if (fp==NULL)
			printf("Error in reading %s\n",stemp);
		csp=evp->cschain;
		while (csp!=NULL)
		{
			for (i=0;i<PNUM;i++)
			{
				if ((IsGoodData(csp,i)>0))
				{
					for (j=0;j<3;j++)
						grade[j] = csp->meas[i][j].pathG;
					k = GetBestCircCS(csp,i);
					if (grade[k] < grade[(k+1)%3] || grade[k] < grade[(k+2)%3])
					{
						fprintf(fp,"%s_%s_%1d.xcor ",csp->sta1->staname,csp->sta2->staname, i);
						fprintf(fp,"localbest: %1d path: %4.1f %4.1f %4.1f \n", \
								k, csp->meas[i][0].pathG,csp->meas[i][1].pathG,\
								csp->meas[i][2].pathG);
					}
					if (GetBestCircCS(csp,i)!=GetTotalBestCircCS(csp,i))
						fprintf(fp," corrected\n ");
				}
			}
			csp=csp->n;
		}
		fclose(fp);

		//output circle skip log file
		sprintf(stemp,"%s.bad",evp->path);
		fp=fopen(stemp,"w");
		if (fp==NULL)
		{
			printf("Error in reading %s\n",stemp);
		}
		csp=evp->cschain;
		while (csp!=NULL)
		{
			for (i=0;i<PNUM;i++)
				if (IsGoodData(csp,i)<0)
				{
					FormatPrintCS(fp, csp, i);
				}
			csp=csp->n;
		}
		fclose(fp);

		// Output CS data
		sprintf(stemp,"%s.cs",evp->path);
		fp=fopen(stemp,"w");
		if (fp==NULL)
		{
			printf("Error in reading %s\n",stemp);
		}
		csp=evp->cschain;
		while (csp!=NULL)
		{
			for (i=0;i<PNUM;i++)
				if (IsGoodData(csp,i)>0)
					FormatPrintCS(fp, csp, i);
			csp=csp->n;
		}
		fclose(fp);
	}

}	

/* Function to average the measurement from sta1 to sta2 and from sta2 to sta1*/

