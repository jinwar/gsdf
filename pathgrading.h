/* This is the main function of this program to remove the circle skip for 
 * GSDF method. */

void pathgrading (CS *cshead) {
	CS *csp,*cfp1,*cfp2;
	STA *sta1,*sta2,*cfsta;
	int i,j,k;
	int bestpath[3];
	int commonfriendnum[PNUM];
	double maxgrade;
	double dterr;
	double dterrG,totalG;
	double sumpathG;
	double averagedt[PNUM],nthroot;
	double a,b,c;

	// Find the best solution for each path and add points
	csp=cshead;
	while (csp!=NULL)
	{
		sta1=csp->sta1;
		sta2=csp->sta2;
		cfsta=FindCommonFriends(csp);
		//printf("%f\n",csp->meas[1].dt);
		if (cfsta==NULL)
		{
			for (j=0;j<PNUM;j++)
			{
				for (i=0;i<3;i++)
				{
					csp->meas[j][i].pathG=0;
					csp->meas[j][i].totalG = csp->meas[j][i].localG;
				}
				csp->cfnum[j]=0;
			}
			csp=csp->n;
			continue;
		}
		for (i=0;i<PNUM;i++)
		{
			averagedt[i]=0;
			commonfriendnum[i]=0;
		}
		cfsta=FindCommonFriends(csp);
		while (cfsta!=NULL)
		{
			cfp1=FindStaCo(sta1,cfsta);
			cfp2=FindStaCo(sta2,cfsta);
			for (i=0;i<PNUM;i++)
			{
				if (IsGoodData(cfp1,i)>0 && IsGoodData(cfp2,i)>0)
				{
					nthroot = Getbestdt(sta1,cfsta,i) + Getbestdt(cfsta,sta2,i);
					//nthroot = nthroot/fabs(nthroot)*pow(fabs(nthroot), 1/commonfriendnum);
					averagedt[i] = averagedt[i] + nthroot;
					commonfriendnum[i]=commonfriendnum[i]+1;
				}
			}
			cfsta=cfsta->tempn;
		}
		//Calculate the average travel time from sta1 to sta2 through their friends by using
		//  nth root method
		for (i=0;i<PNUM;i++)
		{
			csp->cfnum[i]=commonfriendnum[i];
			if (commonfriendnum[i]!=0)
				averagedt[i] = averagedt[i]/commonfriendnum[i];
		}
		//averagedt = averagedt/fabs(averagedt)*pow(fabs(averagedt),commonfriendnum);
		//printf("%f %f\n", averagedt, csp->meas[1].dt);
		//Grading the path
		for (i=0;i<PNUM;i++)
			for (j=0;j<3;j++)
			{
				csp->cfdt[i]=averagedt[i];
				if (commonfriendnum[i]!=0)
				{
					a=fabs(csp->meas[i][j].dt-averagedt[i]);
					b=fabs(csp->meas[i][(j+1)%3].dt-averagedt[i]);
					c=fabs(csp->meas[i][(j+2)%3].dt-averagedt[i]);
					csp->meas[i][j].pathG=100/(1+a/b+a/c);
				} 
				else
					csp->meas[i][j].pathG=0;
			}
		// Calculate the total grade
		for (i=0;i<PNUM;i++)
			for (j=0;j<3;j++)
				csp->meas[i][j].totalG=(csp->meas[i][j].localG*totalGweight[0]\
							+ csp->meas[i][j].pathG*totalGweight[1])\
							/ (totalGweight[0]+totalGweight[1]);
		csp=csp->n;
	}

}
	
