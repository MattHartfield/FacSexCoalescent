/* FacSexCoalescent.c 
The facultative sex coalescent program...in C!!!

< Add further preamble here once the program is near release - e.g. runtime instructions, etc. >

Simulation uses routines found with the GNU Scientific Library (GSL)
(http://www.gnu.org/software/gsl/)
Since GSL is distributed under the GNU General Public License 
(http://www.gnu.org/copyleft/gpl.html), you must download it 
separately from this file.

*/

/* Preprocessor statements */
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Function prototypes */
/* TO BE ADDED */

/* Global variable declaration */
unsigned int N = 0;			/* Population Size */
double rec = 0;				/* Per-site recombination rate */
unsigned int nsites = 0;	/* Number of sites (for recombination) */
double g = 0;				/* Broad-scale gene conversion */
double theta = 0;			/* Scaled mutation rate, 4Nmu */
unsigned int pSTIN = 0;		/* Determine type of heterogeneity (0=fluctuating sex, 1=stepwise change, 2 = constant sex) */
double pLH = 0;				/* Prob of low-sex to high-sex transition, OR time of transition if stepwise change */
double pHL = 0;				/* Prob of high-sex to low-sex transition */
double mig = 0;				/* Migration rate between demes */
unsigned int d = 0;			/* Number of demes */
unsigned int Itot = 0;		/* Number of initial samples */
unsigned int Iindv = 0;		/* Number of initial individuals */
unsigned int Nreps = 0;		/* Number of simulation replicates */

/* Main program */
int main(int argc, char *argv[]){
	unsigned int x;		/* Assignment counter */

	/* GSL random number definitions */
	const gsl_rng_type * T; 
	gsl_rng * r;
	
	/* Reading in data from command line. */
	/*
	if(argc != 8){
		fprintf(stderr,"Invalid number of input values.\n");
		exit(1);
	}
	*/
	N = strtod(argv[1],NULL);
	rec = strtod(argv[2],NULL);
	nsites = strtod(argv[3],NULL);
	g = strtod(argv[4],NULL);
	theta = strtod(argv[5],NULL);
	pSTIN = strtod(argv[6],NULL);
	pLH = strtod(argv[7],NULL);
	pHL = strtod(argv[8],NULL);
	mig = strtod(argv[9],NULL);
	d = strtod(argv[10],NULL);
	if(d == 1){
		mig = 0;	/* Set migration to zero if only one deme, as a precaution */
	}
	if(rec == 0){
		nsites = 1; /* Set no sites to 1 if no recombination, as a precaution */
	}
	if(N%d != 0){
		fprintf(stderr,"Population size must be a multiple of deme number.\n");
		exit(1);
	}
	N = N/(d*1.0);	/* Scaling NT to a demetic size, for more consistent use in calculations */
	
	/* Initial Error checking */
	if(N <= 0){
		fprintf(stderr,"Total Population size N is zero or negative, not allowed.\n");
		exit(1);
	}
	if(g < 0 || g > 1){
		fprintf(stderr,"Rate of gene conversion has to lie between 0 and 1.\n");
		exit(1);
	}
	if(g > (10.0/(1.0*N))){
		printf("WARNING: Analytical transitions assume gene conversion is weak.\n");
		printf("Current input is much larger than O(1/NT) - transitions may be inaccurate.\n");
	}
	if(theta < 0){
		fprintf(stderr,"Mutation rate must be a positive (or zero) value.\n");
		exit(1);
	}
	if(rec < 0){
		fprintf(stderr,"Recombination rate must be a positive (or zero) value.\n");
		exit(1);
	}
	if(nsites <= 0){
		fprintf(stderr,"Number of sites must be a positive integer.\n");
		exit(1);
	}
	if(mig < 0){
		fprintf(stderr,"Migration rate must be a positive (or zero) value.\n");
		exit(1);
	}
	if(mig > (10.0/(1.0*N))){
		printf("WARNING: Analytical transitions assume migration is weak.\n");
		printf("Current input is much larger than O(1/NT) - transitions may be inaccurate.\n");
	}
	if(d <= 0){
		fprintf(stderr,"Number of demes has to be a positive integer.\n");
		exit(1);
	}
	if((d > 1) && (mig == 0)){
		fprintf(stderr,"Migration rate cannot be zero with multiple demes.\n");
		exit(1);
	}
	if(pSTIN != 1){
		if(pHL < 0 || pHL > 1 || pLH < 0 || pLH > 1){
			fprintf(stderr,"Sex transition probabilities have to lie between 0 and 1 (if there is no stepwise change in sex).\n");
			exit(1);
		}
	}
	if( (pHL == 0 || pLH == 0 ) && pSTIN == 0){
		fprintf(stderr,"Sex transition probabilities have to lie between 0 and 1.\n");
		exit(1);
	}
	if(pSTIN != 0 && pSTIN != 1 && pSTIN != 2){
		fprintf(stderr,"pSTIN has to equal 0, 1, or 2.\n");
		exit(1);
	}
	
	/*
	Initial state of two samples;
	Iwith = No. of within host
	Ibet = No. of between host
	So total number of initial samples in a deme = 2*Iwith + Ibet
	THEN sexL[i], sexH[i] = low-sex, high-sex rate in each deme
	*/
	
	double *Iwith = calloc(d,sizeof(unsigned int));		/* Within-individual samples */
	double *Ibet = calloc(d,sizeof(unsigned int));		/* Between-individual samples */
	double *sexL = calloc(d,sizeof(double));		/* Low-sex rates */	
	double *sexH = calloc(d,sizeof(double));		/* High-sex rates */
	
	for(x = 0; x < d; x++){
		*(Iwith + x) = strtod(argv[11 + (4*x + 0)],NULL);
		*(Ibet + x) = strtod(argv[11 + (4*x + 1)],NULL);
		*(sexL + x) = strtod(argv[11 + (4*x + 2)],NULL);
		*(sexH + x) = strtod(argv[11 + (4*x + 3)],NULL);
		Itot += 2*(*(Iwith + x)) + (*(Ibet + x));
		Iindv += (*(Iwith + x)) + (*(Ibet + x));
		
		/* More error checking as inputting */
		if( *(sexL + x) < 0 || *(sexL + x) > 1 || *(sexH + x) < 0 || *(sexH + x) > 1){
			fprintf(stderr,"Rate of sexual reproduction has to lie between 0 and 1.\n");
			exit(1);
		}
		if((*(sexH + x) == 0) && (pSTIN != 2)){
			fprintf(stderr,"Rate of sexual reproduction has to lie between 0 and 1.\n");
			exit(1);
		}
		if((*(sexL + x) > *(sexH + x)) && pSTIN == 0){
			fprintf(stderr,"All low-sex values must be less than or equal to high-sex values. Please re-check.\n");
			exit(1);
		}
		if( (*(Iwith + x) < 0) || (*(Ibet + x) < 0) ){
			fprintf(stderr,"Number of within- and between-host samples has to be positive.\n");
			exit(1);
		}
	}
	
	/* Number of samples/reps to take */
	Nreps = strtod(argv[4*d + 11],NULL);
	
	/* Final Error Checking */
	if(Itot <= 1){
		fprintf(stderr,"More than one sample must be initially present to execute the simulation.\n");
		exit(1);		
	}
	if(Nreps <= 0){
		fprintf(stderr,"Must set positive number of repetitions.\n");
		exit(1);
	}
	
	/* Arrays definition and memory assignment */
	  
	/* create a generator chosen by the 
    environment variable GSL_RNG_TYPE */
     
	gsl_rng_env_setup();
	if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	printf("%lu\n",gsl_rng_default_seed);
	
	/* Freeing memory and wrapping up */
 	gsl_rng_free(r);
 	free(sexH);
 	free(sexL);
 	free(Ibet);
 	free(Iwith); 	 	 	
 	
	return 0;
	
}	/* End of main program */

/* End of File */