// Define how many circles concerned, has to be 3
#define CIRCNUM 3 

// Define the distance within which the measurements are made
#define DIST 250 
#define MINDIST 10 

// Opinion whether to make xcor files. set to 0 on the second run.
#define ISMAKEFILES 1

// Define pre-filter
#define ISPREFILTER 1
float prefilter[2]={ 15,300 };

// When windowing the xcor, whether window it center at peak or at theoritical time based on average group velocity
#define ISWINDOWMAX 1

// Define the maximum points in sac files
#define MAXDP 100000

// Define the smallest station distance
#define STA_TOL_DIST 1

#define PI 3.1415926

// reference phase velocity file
#define PHMODFILE "tna_phv.mod"

// reference group velocity file
#define GRMODFILE "tna_grv.mod"

// Define the number of periods that model files include
#define MODNUM 281

// Define the number of period being measured. 
#define PNUM 3

// Define interested periods and narrow band filters
int periods[PNUM]={25, 50, 100};

// Define the ID of center period, relate to period[ID]
#define CENTP 1

// Define the model
float phvmodel[MODNUM][2];
float grvmodel[MODNUM][2];
float grvmax,grvmin,grvavg;

#define ISLIMITGRV 0
#define GRVMAXU 4.3
#define GRVMAXL 3.6
#define GRVMINU 3.3
#define GRVMINL 2.7


// Define the half width of windowing process
#define HALFWINDOW 70

// Define the cut window before and after the theoritical arrival time
#define CUTBEFORE 30
#define CUTAFTER 60

// Define the half width of taper windowing process
#define TAPERWIDTH 30

// Define the number of circle of cross-correlation curve to be fit
#define NFIT 1

// define the weight of each local grade for each period: refG contG
float localGweight[PNUM][2] = { {20, 80},{100, 0}, {70, 30}};

// define the weight of the local grade and path grade
float totalGweight[2] = {70, 30};

// Define the maximum error that are allowed in the fiting function
#define MAXCHI 0.1

// Define the minimum rate of energy contained by the windowed xcor waveform to the whole xcor waveform that can be output
#define MINENERGYRATE 5

// Define the max misfit allowed from the model prediction
//#define MAXDTMISFIT 10

// Define the minimum number of measure to fit the average phase velocity
#define FITPHASEVMINNUM 5

// Define whether turn on the station location selection
#define ISSELECTSTA 1
#define MINSTALA 27
#define MAXSTALA 51
#define MINSTALO -114
#define MAXSTALO -95
