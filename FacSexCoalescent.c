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
/* TO BE ADDED */
/* Probably the input variables... */

/* Main program */
int main(int argc, char *argv[]){

	/* GSL random number definitions */
	const gsl_rng_type * T; 
	gsl_rng * r;
	
	/* This reads in data from command line. Random seed is piped in
	before program is read. */
	if(argc != 8){
		fprintf(stderr,"Invalid number of input values.\n");
		exit(1);
	}
	N = strtod(argv[1],NULL);
	s = strtod(argv[2],NULL);
	rec = strtod(argv[3],NULL);
	sex = strtod(argv[4],NULL);
	self = strtod(argv[5],NULL);
	gc = strtod(argv[6],NULL);
	reps = strtod(argv[7],NULL);
	
	/* Arrays definition and memory assignment */
	  
	/* create a generator chosen by the 
    environment variable GSL_RNG_TYPE */
     
	gsl_rng_env_setup();
	if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	printf("%d\n",gsl_rng_default_seed);
	
	/* Freeing memory and wrapping up */
 	gsl_rng_free(r);
	return 0;
}	/* End of main program */

/* End of File */