/* This is the main function of this program to remove the circle skip of
 * GSDF method. */
void pathgrading (CS *cshead) {
	CS *csp,*cfp1,*cfp2;
	STA *sta1,*sta2,*cfsta;
	int i,j,k;
	int bestpath[3];
	double maxgrade;
	double dterr;
	double dterrG,totalG;
	double sumpathG;

	// Find the best solution for each path and add points
	csp=cshead;
	while (csp!=NULL)
	{
		sta1=csp->sta1;
		sta2=csp->sta2;
		cfsta=FindCommonFriends(csp);
		while (cfsta!=NULL)
		{
			cfp1=FindStaCo(sta1,cfsta);
			cfp2=FindStaCo(sta2,cfsta);
			maxgrade=0;
			for (i=0;i<CIRCNUM;i++)
				for (j=0;j<CIRCNUM;j++)
					for (k=0;k<CIRCNUM;k++)
					{
						dterr=Getdt(sta1,cfsta,j)+Getdt(cfsta,sta2,k) \
							  +Getdt(sta2,sta1,i);
						dterr=fabs(dterr);
						if (dterr > csp->centerF)
							dterrG=0;
						else
							dterrG=100-dterr/csp->centerF*100;
						totalG=csp->meas[i].localG+cfp1->meas[j].localG\
							   +cfp2->meas[k].localG;
						totalG=totalG/3;
						totalG+=dterrG;
						if (maxgrade<totalG)
						{
							maxgrade=totalG;
							bestpath[0]=i;
							bestpath[1]=j;
							bestpath[2]=k;
						}
					}
			if (bestpath[0]!=1 || bestpath[1]!=1 || bestpath[2]!=1)
			{
			printf("%d %d %d \n",bestpath[0],bestpath[1],bestpath[2]);
			printf("%f %f %f \n",csp->meas[bestpath[0]].localG, \
					cfp1->meas[bestpath[1]].localG, \
					cfp2->meas[bestpath[2]].localG);
			printf("%f %f %f \n",csp->meas[1].localG, \
					cfp1->meas[1].localG, \
					cfp2->meas[1].localG);
			}
			csp->meas[bestpath[0]].pathG++;
			cfp1->meas[bestpath[1]].pathG++;
			cfp2->meas[bestpath[2]].pathG++;
			cfsta=cfsta->tempn;
		}
		csp=csp->n;
	}

	// Normalize the grade into percent
	csp=cshead;
	while (csp!=NULL)
	{
		sumpathG=0;
		for (i=0;i<CIRCNUM;i++)
			sumpathG+=csp->meas[i].pathG;
		for (i=0;i<CIRCNUM;i++)
			csp->meas[i].pathG=csp->meas[i].pathG*100/sumpathG;
		csp=csp->n;
	}
}

	
