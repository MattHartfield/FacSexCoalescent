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
unsigned int trig(unsigned int x);
double P23(unsigned int y, unsigned int k, unsigned int Na);
double P4(unsigned int x, unsigned int k, unsigned int Na);
double P56(unsigned int y, unsigned int Na);
double P7(unsigned int x, unsigned int k, unsigned int Na);
double P8(unsigned int x, unsigned int y, unsigned int k, unsigned int Na);
double P9(unsigned int x, unsigned int k, double gee);
double P10(unsigned int x, unsigned int y, unsigned int k, double mee);
double P11(unsigned int y, unsigned int k, double sexC, double ree, unsigned int lrec, unsigned int nlrec, unsigned int nlrec2);
void probset2(unsigned int N, double g, double *sexC, double rec, unsigned int lrec, unsigned int *nlrec, unsigned int *nlrec2, double mig, unsigned int *Nwith, unsigned int *Nbet, unsigned int *kin, unsigned int sw, double **pr);
void rate_change(unsigned int pST,double pLH, double pHL, double *sexH, double *sexL, unsigned int Na, unsigned int d, unsigned int switch1, double *sexCN, double *tts, unsigned int *npST,const gsl_rng *r);
unsigned int isanyUI(unsigned int *vin, unsigned int size_t, unsigned int match);
unsigned int isanyD(double *vin, unsigned int size_t, double match);
unsigned int isallUI(unsigned int *vin, unsigned int size_t, unsigned int match);
unsigned int isallD(double *vin, unsigned int size_t, double match);

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

/* Function to replicate the 'any' func in R (unsigned int) */
unsigned int isanyUI(unsigned int *vin, unsigned int size_t, unsigned int match){
	unsigned int i;
	unsigned int res = 0;
	for(i = 0; i < size_t; i++){
		if(*(vin + 1) == match){
			res = 1;
		}
	}
	return res;
}

/* Function to replicate the 'any' func in R (double) */
unsigned int isanyD(double *vin, unsigned int size_t, double match){
	unsigned int i;
	unsigned int res = 0;
	for(i = 0; i < size_t; i++){
		if(*(vin + 1) == match){
			res = 1;
		}
	}
	return res;
}

/* Function to replicate the 'all' func in R (unsigned int) */
unsigned int isallUI(unsigned int *vin, unsigned int size_t, unsigned int match){
	unsigned int i;
	unsigned int res = 1;
	for(i = 0; i < size_t; i++){
		if(*(vin + 1) != match){
			res = 0;
		}
	}
	return res;
}

/* Function to replicate the 'all' func in R (double) */
unsigned int isallD(double *vin, unsigned int size_t, double match){
	unsigned int i;
	unsigned int res = 1;
	for(i = 0; i < size_t; i++){
		if(*(vin + 1) != match){
			res = 0;
		}
	}
	return res;
}

/* 'Triangle function' calculation */
unsigned int trig(unsigned int x){
	return (x*(x-1))/2.0;
}

/* Functions of each transition probability calculation */
double P23(unsigned int y, unsigned int k, unsigned int Na){
	/* One of the paired samples is recreated: (x,y) -> (x-k+1, y + 2(k-1)). 
	OR A unique samples coaleses with another (either pre-existing or new): (x,y) -> (x-k, y + 2k - 1) */
	return (k*y)/(1.0*Na) + trig(2.0*k)/(2.0*Na);
}
double P4(unsigned int x, unsigned int k, unsigned int Na){
	/* A paired sample coaleses with new unique sample: (x,y) -> (x-k,y + 2k -1) */
	return ((x-k)*2.0*k)/(1.0*Na);
}
double P56(unsigned int y, unsigned int Na){
	/* Two pre-existing unique samples re-create a paired sample: (x,y) -> (x - k + 1, y + 2(k-1)). 
	OR Two pre-existing paired samples coalesce: (x,y) -> (x - k, y + 2k-1) */
	return trig(y)/(2.0*Na);
}
double P7(unsigned int x, unsigned int k, unsigned int Na){
	/* Two remaining paired samples doubly coalesce asexually: (x,y) -> (x-k-1,y+2k) */
	return trig(x-k)/(1.0*Na);
}
double P8(unsigned int x, unsigned int y, unsigned int k, unsigned int Na){
	/* One of the x - k remaining paired samples can coalesce with a unique sample: (x,y) -> (x-k,y+2k-1) */
	return ((x-k)*y)/(1.0*Na);
}
double P9(unsigned int x, unsigned int k, double gee){
	/* Paired sample coaleses via gene conversion: (x,y) -> (x-k-1,y+2k+1) */
	/* NEEDS TO BE UPDATED TO ACCOUNT FOR MULTIPLE SITES */
	return gee*(x-k);
}
double P10(unsigned int x, unsigned int y, unsigned int k, double mee){
	/* A sample migrates to another deme */
	return mee*(x + y + k);
}
double P11(unsigned int y, unsigned int k, double sexC, double ree, unsigned int lrec, unsigned int nlrec, unsigned int nlrec2){
	/* One of the single samples splits by recombination, creates two new single samples: (x,y) -> (x-k,y+2k+1) */
	return (sexC*rec*((lrec - 1)*(y + 2*k) - nlrec) + rec*((lrec - 1)*(2*k) - nlrec2));
}

/* Calculate probability change vectors each time OVER EACH DEME */
void probset2(unsigned int N, double g, double *sexC, double rec, unsigned int lrec, unsigned int *nlrec, unsigned int *nlrec2, double mig, unsigned int *Nwith, unsigned int *Nbet, unsigned int *kin, unsigned int sw, double **pr){
	unsigned int x;				/* Deme counter */
	unsigned int ksum = 0;		/* Total number of segregating events */
	
	/* First calculate number of splits... */
	for(x = 0; x < d; x++){
		ksum += *(kin + x);
	}
	
	/* Calculate each transition probability, per deme */
	for(x = 0; x < d; x++){
		*((*(pr + 5)) + x) = P56(*(Nbet + x),N);
		*((*(pr + 6)) + x) = P56(*(Nbet + x),N);
		*((*(pr + 7)) + x) = P7(*(Nwith + x),*(kin + x),N);
		*((*(pr + 8)) + x) = P8(*(Nwith + x),*(Nbet + x),*(kin + x),N);
		*((*(pr + 9)) + x) = P9(*(Nwith + x),*(kin + x),g);
		*((*(pr + 10)) + x) = P10(*(Nwith + x),*(Nbet + x),*(kin + x),mig);		
		*((*(pr + 11)) + x) = P11(*(Nbet + x),*(kin + x),*(sexC + x),rec,lrec,*(nlrec + x),*(nlrec2 + x));
		
		/* Only activate the first three events if need to consider segregation via sex 
		(fourth is 'split pairs remain split') */
		if(sw == 1){
			if(ksum != 1){
				*((*(pr + 2)) + x) = P23(*(Nbet + x),*(kin + x),N);
			}else if(ksum == 1){
				*((*(pr + 2)) + x) = 0;
			}
			*((*(pr + 3)) + x) = P23(*(Nbet + x),*(kin + x),N);
			*((*(pr + 4)) + x) = P4(*(Nwith + x),*(kin + x),N);
		}else if(sw == 0){
			*((*(pr + 2)) + x) = 0;
			*((*(pr + 3)) + x) = 0;
			*((*(pr + 4)) + x) = 0;
		}
	}
}	/* End of 'probset2' function */

/* Function to change rates of sex given a state change */
void rate_change(unsigned int pST,double pLH, double pHL, double *sexH, double *sexL, unsigned int Na, unsigned int d, unsigned int switch1, double *sexCN, double *tts, unsigned int *npST, const gsl_rng *r){
	unsigned int x = 0;			/* Deme counter */
	
	/* Setting up transition time (tts, or 'time to switch')	*/
	if(pST == 0){
		for(x = 0; x < d; x++){
			*(sexCN + x) = *(sexL + x);
		}
		*tts = gsl_ran_geometric(r,pLH)/(2.0*Na*d);
		*npST = 1;
	}else if(pST == 1){
		for(x = 0; x < d; x++){
			*(sexCN + x) = *(sexH + x);
		}
		*tts = gsl_ran_geometric(r,pHL)/(2.0*Na*d);
		*npST = 0;
	}else if(pST == 2){		/* If stepwise change, alter depending on whether already switched or not */
		*npST = pST;
		if(switch1 == 0){
			for(x = 0; x < d; x++){
				*(sexCN + x) = *(sexL + x);
			}
			*tts = pLH;
		}else if(switch1 == 1){
			for(x = 0; x < d; x++){
				*(sexCN + x) = *(sexH + x);
			}
			*tts = (0.0/0.0);
		}
	}else if(pST == 3){		/* If constant, no time to sex switch */
		for(x = 0; x < d; x++){
				*(sexCN + x) = *(sexL + x);
			}
		*tts = (0.0/0.0);
		*npST = pST;
	}
	
	printf("%lf %d \n",*tts,*npST);
	
}	/* End of 'rate_change' function */

/* Main program */
int main(int argc, char *argv[]){
	unsigned int x, i, j;		/* Assignment counter, rep counter, indv counter */
	unsigned int pST, npST = 0;		/* State of reproduction heterogeneity */	
	double Ttot = 0;			/* Time in past, initiate at zero */
	unsigned int Ntot = 0;		/* Total number of samples at time */
	unsigned int Nindv = 0;		/* Total number of individuals */
	unsigned int lrec = 0;		/* Length of recombinable genome */
	unsigned int IwithT = 0;	/* Total within individual samples */
	unsigned int IbetT = 0;		/* Total between individual samples */	
	unsigned int IwithC = 0;	/* Cumulative Iwith sum */
	unsigned int IbetC = 0;		/* Cumulative Ibet sum */
	double tls = 0;				/* 'Time since Last Switch' or tls */
	double tts = 0;				/* 'Time to switch' */
	unsigned int done = 0;		/* Is simulation complete? */
	unsigned int nbreaks = 0;	/* Number of non-rec tracts */

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
	
	unsigned int *Iwith = calloc(d,sizeof(unsigned int));		/* Within-individual samples */
	unsigned int *Ibet = calloc(d,sizeof(unsigned int));		/* Between-individual samples */
	unsigned int *Nwith = calloc(d,sizeof(unsigned int));		/* To be used in sim */
	unsigned int *Nbet = calloc(d,sizeof(unsigned int));		/* To be used in sim */
	double *sexL = calloc(d,sizeof(double));		/* Low-sex rates */	
	double *sexH = calloc(d,sizeof(double));		/* High-sex rates */
	
	for(x = 0; x < d; x++){
		*(Iwith + x) = strtod(argv[11 + (4*x + 0)],NULL);
		*(Ibet + x) = strtod(argv[11 + (4*x + 1)],NULL);
		*(sexL + x) = strtod(argv[11 + (4*x + 2)],NULL);
		*(sexH + x) = strtod(argv[11 + (4*x + 3)],NULL);
		IwithT += (*(Iwith + x));
		IbetT += (*(Ibet + x));
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
		/*
		if( (*(Iwith + x) < 0) || (*(Ibet + x) < 0) ){
			fprintf(stderr,"Number of within- and between-host samples has to be positive.\n");
			exit(1);
		}
		*/
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
	unsigned int *nlrec = calloc(d,sizeof(unsigned int));			/* Non-recombinable samples 1 */
	unsigned int *nlrec2 = calloc(d,sizeof(unsigned int));			/* Non-recombinable samples 2 */	
	double *sexC = calloc(d,sizeof(double));						/* Current rates of sex per deme */	
	double **pr = calloc(11,sizeof(double *));						/* Probability matrix per deme */
	for (x = 0; x < 11; x++){										/* Assigning space for each population within each deme */
		pr[x] = calloc(d,sizeof(double));
	}
	  
	/* create a generator chosen by the 
    environment variable GSL_RNG_TYPE */
     
	gsl_rng_env_setup();
	if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	printf("%lu\n",gsl_rng_default_seed);
	
	/* 
	In R code I defined lists here for:
	(1) NEWICK tree of each run; and
	(2) List of mutation matrices.
	I'll leave these for now: my plan now is to form each of these after each run, 
	and print out to individual files instead.
	(If I need to create such an object, probably need to create a 'struct' - need
	to relearn how to do such things!)
	*/
	
	/* Running the simulation Nreps times */
	for(i = 0; i < Nreps; i++){

		/* Setting up type of sex heterogeneity */
		if(pSTIN == 0){
			pST = gsl_ran_bernoulli(r,pLH/(pLH+pHL));	/* Randomly assigning initial low-sex or high-sex state */
		}else if(pSTIN == 1){
			pST = 2;
		}else if(pSTIN == 2){
			pST = 3;
		}
		Ttot = 0;
		Ntot = Itot;
		Nindv = Iindv;
		lrec = nsites;
	
		for(x = 0; x < d; x++){
			*(Nwith + x) = *(Iwith + x);	/* Resetting number of within-host samples */
			*(Nbet + x) = *(Ibet + x);		/* Resetting number of between-host samples */
			*(nlrec + x) = 0;				/* Genome not affected by recombination */
		}
		
		/* Setting up temporal heterogeneity */
		rate_change(pST,pLH,pHL,sexH,sexL,N,d,0,sexC,&tts,&npST,r);
		pST = npST;
		tls = 0;
		
		/* Setting up summary table of individual samples */
		/* ASSIGNING MEMORY FROM SCRATCH HERE, SINCE TABLES WILL BE MODIFIED FOR EACH SIM */
		
		unsigned int **indvs = calloc(Itot,sizeof(unsigned int *));		/* Table of individual samples */
		unsigned int **GType = calloc(Itot,sizeof(unsigned int *));		/* Table of sample genotypes */
		unsigned int **CTms = calloc(Itot,sizeof(unsigned int *));		/* Coalescent times per sample */
		unsigned int **TAnc = calloc(Itot,sizeof(unsigned int *));		/* Table of ancestors for each sample */
		unsigned int **breaks = calloc(2,sizeof(unsigned int *));		/* Table of breakpoints created in the simulation */
		for (x = 0; x < 11; x++){										/* Assigning space for each population within each deme */
			pr[x] = calloc(d,sizeof(double));
		}
		for(j = 0; j < Itot; j++){										/* Assigning space for each genome sample */
			indvs[j] = calloc(4,sizeof(unsigned int));
			GType[j] = calloc(2,sizeof(unsigned int));
			CTms[j] = calloc(2,sizeof(unsigned int));
			TAnc[j] = calloc(2,sizeof(unsigned int));
		}
		breaks[0] = calloc(1,sizeof(unsigned int));
		breaks[1] = calloc(1,sizeof(unsigned int));
		unsigned int *bcoal = calloc(1,sizeof(unsigned int));
		nbreaks = 1;
		
		IwithC = 0;
		IbetC = 0;
		for(j = 0; j < Itot; j++){
			*((*(indvs + j)) + 0) = j;
			*((*(GType + j)) + 0) = j;
			*((*(GType + j)) + 1) = j;
			*((*(CTms + j)) + 0) = j;
			*((*(TAnc + j)) + 0) = j;
		}
		if(IwithT != 0){
			for(j = 0; j < IwithT; j++){
				*((*(indvs + 2*j)) + 1) = j;
				*((*(indvs + (2*j+1))) + 1) = j;
				*((*(indvs + 2*j)) + 2) = 0;
				*((*(indvs + (2*j+1))) + 2) = 0;
			}
			for(x = 0; x < d; x++){
				for(j = IwithC; j < (IwithC + *(Iwith + x)); j++){
					*((*(indvs + 2*j)) + 3) = x;
					*((*(indvs + (2*j+1))) + 3) = x;
				}
				IwithC += *(Iwith + x);
			}
		}
		if(IbetT != 0){
			for(j = 0; j < IbetT; j++){
				*((*(indvs + (2*IwithT + j))) + 1) = IwithT + j;
				*((*(indvs + (2*IwithT + j))) + 2) = 1;
			}
			for(x = 0; x < d; x++){
				for(j = IbetC; j < (IbetC + *(Ibet + x)); j++){
					*((*(indvs + (2*IwithT + j))) + 3) = x;
				}
				IbetC += *(Ibet + x);
			}
		}
		
		for(j = 0; j < Itot; j++){
			for(x = 0; x < 4; x++){
				printf("%d ",*((*(indvs + j)) + x));
			}
			printf("\n");
		}
		
		done = 0;
		while(done != 1){
			/* Simulation code here... */
			
			/* Testing if all sites coalesced or not */
			for(x = 0; x < nbreaks; x++){
				*(bcoal + x) = *((*(breaks + 1)) + x);
			}
			done = isallUI(bcoal,nbreaks,1);
		}
		
		/* ...then freeing the memory */
		free(bcoal);
		free(breaks[1]);
		free(breaks[0]);
		for(j = 0; j < Itot; j++){
			free(TAnc[j]);
			free(CTms[j]);
			free(GType[j]);
			free(indvs[j]);
		}
		free(breaks);
		free(TAnc);
		free(CTms);
		free(GType);		
		free(indvs);
	}
	
	/* Freeing memory and wrapping up */
 	gsl_rng_free(r);
 	for(x = 0; x < 11; x++){
		free(pr[x]);
	}
	free(pr);
	free(sexC);
	free(nlrec2);
	free(nlrec);
 	free(sexH);
 	free(sexL);
 	free(Nbet);
 	free(Nwith);
 	free(Ibet);
 	free(Iwith); 	 	 	
 	
	return 0;
	
}	/* End of main program */

/* End of File */