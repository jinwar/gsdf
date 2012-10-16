#define STA struct station
#define STAINF struct stationinform
#define CS struct costa
#define CSCH struct cschain
#define EV struct event
#define EVLEN sizeof(struct event)
#define STALEN sizeof(struct station)
#define CSLEN sizeof(struct costa)
#define CSCHLEN sizeof(struct cschain)

struct event{
	float la;
	float lo;
	float depth;
	float epidist;
	int year;
	int jday;
	int month;
	int day;
	int hour;
	int min;
	int stationnum;
	double sec;
	double abt;
	float grvmax,grvmin;
	float grv0,grv1,t0,t1;
	float avgphv[PNUM];
	char path[40];
	struct event *n;
	STA *stachain;
	CS *cschain;
};

struct fitting{
	double ao, so, wo, tp, tg, chi;
	int ierr;
};

struct gslfit{
	float c0, c1, cov00, cov01, cov11, sumsq;
	int datanum;
	int Isgood;
};

struct station{
	float la;
	float lo;
	float stagrade[PNUM];
	float snr[PNUM];
	float avgphv[PNUM];
	float avgphvdv[PNUM];
	float offset[PNUM];
	float sc, swc, wc;
	int goodnum[PNUM];
	int id;
	float dt;
	float amp;
	float epidist;
	int isinchain;
	char staname[8];
	char filename[100];
	struct gslfit gsl[PNUM];
	struct cschain *csch;
	struct station *n;
	struct station *tempn;
	struct event *ev;
	struct fitting fit[PNUM];
	struct fitting cfffit;
	struct fitting wcfffit;
};

struct stationinform{
	char staname[8];
	float la;
	float lo;
	STA *stachain;
	struct stationinform *n;
};

struct cschain {
	struct costa *cs;
	struct cschain *n;
};


struct measurement{
	double dt;
	float refG;
	float contG;
	float localG;
	float pathG;
	float totalG;
};

struct costa {
	float snr;
	float tc;
	int dterracpt[PNUM];
	float epidiff;
	float tprecover[PNUM];
	double errlevel[PNUM];
	float cohere[PNUM];
	double dist;
	int cfnum[PNUM];
	int isinchain;
	double cfdt[PNUM];
	double taup[PNUM];
	double taug[PNUM];
	struct fitting fit[PNUM];
	struct measurement meas[PNUM][CIRCNUM];
	struct costa *n;
	struct costa *pair;
	struct station *sta1, *sta2;
	struct event *ev;
};


