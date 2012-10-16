/* This is the head file contain the functions to construct a data base
 * structure from sac file system for GSDF automatic searching program gsdfmain.c
 */
STA *CreateStationChain(EV *evp) {
	STA *p,*q,*head,*tail;
	int stationnum;
	int i,j,k;
	int columnnum=20, rownum=20;
	int nlen, nerr, max = MAXDP;		// for reading the sac
	int nzsec, nzmsec;
	char stemp[100];
	char sacfilename[100];
	float beg, del;
	float epidist;
	float data[MAXDP];
	float tt1,tt2;
	float b, e;
	FILE *saclistfp;
	FILE *sacmacrofp;


	// Create the sac file list
	printf("Creating station chain......\n");
	sprintf(stemp,"ls %s/*%3s*sac > saclist.txt\n",evp->path,COMP);
	printf("%s\n",stemp);
	system(stemp);

	// read the sac file list
	saclistfp=fopen("saclist.txt","r");
	if (saclistfp==NULL)
		printf("Error in reading saclist.txt\n");

	if (saclistfp==NULL)
	{
		printf("Can't open the sac list file!\n");
		return NULL;
	}
			
	// loop over the sac files to build the station chain
	// and fill the event information
	stationnum=0;
	fgets(stemp,100,saclistfp);
	while(!feof(saclistfp))
	{
		sscanf(stemp,"%s\n",sacfilename);
		// Read in sac file
		rsac1( sacfilename, data, &nlen, &beg, &del, &max, &nerr, strlen(sacfilename) ) ;
		if ( nerr != 0 ) {
			fprintf(stderr, "Error reading in SAC file: %s\n", sacfilename);
			//exit ( nerr ) ;
			fgets(stemp,100,saclistfp);
			continue;
		}

		// Check whether the sac file contain the data we need.
		getfhv("dist", &epidist, &nerr, 4);
		getfhv("b", &b, &nerr, 1);
		e = b + nlen*del;
		tt1=evp->t0 + epidist/evp->grv0;
		tt2=evp->t1 + epidist/evp->grv1;
		
		if (b > tt1-CUTBEFORE || e < tt2 + 100)
		{
			printf(" The sac files do not contain this phase! \n");
			printf(" in the file %s \n", sacfilename);
			fgets(stemp,100,saclistfp);
			continue;
		}

		// Create station structure
		p= (STA *)malloc(STALEN);

		// Fill station info 
		getfhv("stla", &p->la, &nerr, 4);
		getfhv("stlo", &p->lo, &nerr, 4);

		// Check whether the station is in the range
		if (ISSELECTSTA)
			if (p->la < MINSTALA || p->la > MAXSTALA || \
				p->lo < MINSTALO || p->lo > MAXSTALO)
		{
			printf(" Station is not in the range \n");
			printf(" in the file %s \n", sacfilename);
			free(p);
			fgets(stemp,100,saclistfp);
			continue;
		}

		getkhv("kstnm", p->staname, &nerr, 5,8);
		EndStaname(p);			// put \0 to the end of station name
		strcpy(p->filename,sacfilename);

		// Link the station structure to the event
		p->ev=evp;
		/*printf("%s %f %f %s\n",p->filename, p->la, p->lo, p->staname);*/

		// Put this station into the station chain
		if (stationnum==0)
		{
			head=p;
			tail=p;
			// Fill the event information 
			getfhv("evla", &evp->la, &nerr, 4);
			getfhv("evlo", &evp->lo, &nerr, 4);
			getfhv("evdp", &evp->depth, &nerr, 4);
			getnhv("nzyear", &evp->year, &nerr, 6);
			getnhv("nzjday", &evp->jday, &nerr, 6);
			getnhv("nzhour", &evp->hour, &nerr, 6);
			getnhv("nzmin", &evp->min, &nerr, 5);
			getnhv("nzsec", &nzsec, &nerr, 5);
			getnhv("nzmsec", &nzmsec, &nerr, 6);
			evp->sec=nzsec+nzmsec/1000.0;
			evp->abt=n2abtimej(evp->year, evp->jday, evp->hour, evp->min, evp->sec);
		}
		else
		{
			// Check whether
			getnhv("nzsec", &nzsec, &nerr, 5);
			if ( evp->sec - nzsec > 1 )
			{
				printf(" The original time for sac files are not same! \n");
				printf(" in the file %s \n", sacfilename);
				free(p);
				fgets(stemp,100,saclistfp);
				continue;
			}

			// Check whether there's station too close 
			if (!IsExistSta(head,p))
			{
				tail->n=p;
				tail=p;
			}
			else
			{
				free(p);
				fgets(stemp,100,saclistfp);
				continue;
			}
		}

		p->n=NULL;
		p->tempn=NULL;
		p->csch=NULL;
		p->id=stationnum;
		p->epidist=distance(p->la,p->lo,evp->la,evp->lo);

		//station index
		stationnum++;
		fgets(stemp,100,saclistfp);
	}

	// Calculate the SNR of each frequency.
	p=head;
	while (p!=NULL)
	{
		SetStationSNR(p);
		MakeAcorFiles(p);
		p=p->n;
	}
	return(head);
}

CS *CreateMainCSChain(STA *stahead){
	int csnum,stationnum;
	int i;
	double dist,mindist;
	CS *csp, *csq;
	CS *cshead,*cstail;
	CSCH *cschp, *cschq;
	STA *neighborstahead;
	STA *stap,*staq;

	printf("Creating CS chain......\n");

	dist=DIST;
	mindist=MINDIST;

	csnum=0;
	stap=stahead;
	while (stap!=NULL)  //looping of all stations
	{
		neighborstahead=FindStaByDist(stap,stahead,dist,mindist);
		staq=neighborstahead;
		while(staq!=NULL)  //looping of nearby stations
		{
			csp=FindStaCo(stap,staq);

			// For test
			//csp=NULL;

			if (csp==NULL)  
			{
				//create new costa
				csp=(CS *)malloc(CSLEN);
				if (csnum==0)
				{
					cshead=csp;
					cstail=csp;
				}
				else
				{
					cstail->n=csp;
					cstail=csp;
				}
				csp->n=NULL;
				csnum++;
				//link cs to stations
				if (stap->epidist > staq->epidist)
				{
					csp->sta1=stap;
					csp->sta2=staq;
				}
				else
				{
					csp->sta2=stap;
					csp->sta1=staq;
				}
				//link cs to event
				csp->ev=stap->ev;
				csp->epidiff=epidiffCS(csp);
				for (i=0;i<PNUM;i++)
					csp->dterracpt[i]=1;


				// put cs into sta1's cs chain.
				cschp=(CSCH *)malloc(CSCHLEN);
				cschp->cs=csp;
				cschp->n=NULL;
				if (stap->csch==NULL)
					stap->csch=cschp;
				else
				{
					cschq=stap->csch;
					while (cschq->n!=NULL)
						cschq=cschq->n;
					cschq->n=cschp;
				}

				// put cs into sta2's cs chain.
				cschp=(CSCH *)malloc(CSCHLEN);
				cschp->cs=csp;
				cschp->n=NULL;
				if (staq->csch==NULL)
					staq->csch=cschp;
				else
				{
					cschq=staq->csch;
					while (cschq->n!=NULL)
						cschq=cschq->n;
					cschq->n=cschp;
				}

				// fill the measurement data of these cs
				// Core part of this program
				FillCostaData(csp);

				// Do it again from sta2 to sta1, and average the result.
				csq=(CS *)malloc(CSLEN);
				csq->sta1=csp->sta2;
				csq->sta2=csp->sta1;
				csq->ev=csp->ev;
				csq->epidiff=epidiffCS(csq);
				for (i=0;i<PNUM;i++)
					csq->dterracpt[i]=1;

				csp->pair=csq;
				csq->pair=NULL;

				FillCostaData(csq);

				// Remove the effect of windowing and filter
				CSRecover(csq);
				CSRecover(csp);

				// Select the best local circle 
				GradeLocal(csp);

				}
				staq=staq->tempn;
			}
			stap=stap->n;
		}
		return(cshead);
	}

	int ReadMod(){
		FILE *pfp, *gfp;
		int i;
		char stemp[100];

		pfp=fopen(PHMODFILE,"r");
		if (pfp==NULL)
		{
			printf("Can't open the phase model file!\n");
			exit(0);
		}

		for (i=0;i<MODNUM;i++)
		{
			if(!feof(pfp))
			{
				fgets(stemp,100,pfp);
				sscanf(stemp,"%f %f\n",&phvmodel[i][0],&phvmodel[i][1]);
			}
		}
		fclose(pfp);

		gfp=fopen(GRMODFILE,"r");
		if (gfp==NULL)
		{
			printf("Can't open the group model file!\n");
			exit(0);
		}
		for (i=0;i<MODNUM;i++)
		{
			if(!feof(gfp))
			{
				fgets(stemp,100,gfp);
				sscanf(stemp,"%f %f\n",&grvmodel[i][0],&grvmodel[i][1]);
			}
		}
		fclose(gfp);

		grvmax=0;
		for (i=0;i<MODNUM;i++)
			if (grvmodel[i][0]>=periods[0] && grvmodel[i][0]<=periods[PNUM-1])
				if (grvmax < grvmodel[i][1])
					grvmax = grvmodel[i][1];

		grvmin=999;
		for (i=0;i<MODNUM;i++)
			if (grvmodel[i][0]>=periods[0] && grvmodel[i][0]<=periods[PNUM-1])
				if (grvmin > grvmodel[i][1])
					grvmin = grvmodel[i][1];

		grvavg=(grvmax+grvmin)/2;

		return 1;
	}
