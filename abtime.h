//Declare the functions


int Isleapyear(int year)
{
	if ((year % 4) != 0)
		return 0;
	else if ((year % 100) != 0)
		return 1;
	else if ((year % 400) != 0)
		return 0;
	else
		return 1;
}
//This function is used to calculate absolute seconds to 20000101:00:00:00
double n2abtime(int year, int month, int day, int hour, int min, double sec)
{
	int i,sumday;
	int monthday[12]={31,28,31,30,31,30,31,31,30,31,30,31};
	double abst;

	sumday=0;
	for(i=2000;i<year;i++)
	{
		if (Isleapyear(i))
			sumday+=366;
		else
			sumday+=365;
	}
	if (Isleapyear(year))
		monthday[1]=29;
	for(i=0;i<month-1;i++)
		sumday+=monthday[i];
	sumday+=day-1;

	abst=sumday*24*3600.0;
	abst=abst+hour*3600+min*60+sec;
	return abst;

}

double n2abtimej(int year, int day, int hour, int min, double sec)
{
	int i,sumday;
	double abst;

	sumday=0;
	for(i=2000;i<year;i++)
	{
		if (Isleapyear(i))
			sumday+=366;
		else
			sumday+=365;
	}
	sumday+=day-1;

	abst=sumday*24*3600.0;
	abst=abst+hour*3600+min*60+sec;
	return abst;
}

double ab2ntimej(int* year, int* day, int* hour, int* min, double* sec,double abtime)
{

	double absec=abtime;

	*year=1999;
	while (absec > 0)
	{
		*year=*year+1;
		if (Isleapyear(*year))
			absec=absec-366*24*3600.0;
		else
			absec=absec-365*24*3600.0;
	}
	if (Isleapyear(*year))
		absec=absec+366*24*3600.0;
	else absec=absec+365*24*3600.0;

	*day=0;
	while (absec>0)
	{
		*day=*day+1;
		absec=absec-24*3600;
	}
	absec+=24*3600;
	//printf("%f",absec);

	*hour=-1;
	while (absec>0)
	{
		*hour=*hour+1;
		absec=absec-3600;
	}
	absec+=3600;

	*min=-1;
	while (absec>0)
	{
		*min=*min+1;
		absec=absec-60;
	}
	absec+=60;

	*sec=absec;
	
	return 1;

}

double ab2ntime(int* year,int* month, int* day, int* hour, int* min, double* sec,double abtime)
{
	double absec=abtime;
	int monthday[12]={31,28,31,30,31,30,31,31,30,31,30,31};

	*year=1999;
	while (absec > 0)
	{
		*year=*year+1;
		if (Isleapyear(*year))
			absec=absec-366*24*3600.0;
		else
			absec=absec-365*24*3600.0;
	}
	if (Isleapyear(*year))
		absec=absec+366*24*3600.0;
	else absec=absec+365*24*3600.0;

	*month=0;
	if (Isleapyear(*year))
		monthday[1]=29;
	while (absec>0)
	{
		absec=absec-24*3600*monthday[*month];
		*month=*month+1;
	}
	absec+=24*3600*monthday[*month-1];
	
	*day=0;
	while (absec>0)
	{
		*day=*day+1;
		absec=absec-24*3600;
	}
	absec+=24*3600;
	//printf("%f",absec);

	*hour=-1;
	while (absec>0)
	{
		*hour=*hour+1;
		absec=absec-3600;
	}
	absec+=3600;

	*min=-1;
	while (absec>0)
	{
		*min=*min+1;
		absec=absec-60;
	}
	absec+=60;

	*sec=absec;
	
	return 1;

}
