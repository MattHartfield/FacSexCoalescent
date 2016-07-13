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
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_lambert.h>
#include <sys/stat.h>
#include <sys/types.h>

#define INITBR 10
#define HUGEVAL 1000000
#define WARNQ 0.05

/* Function prototypes */
unsigned int isanyUI(unsigned int *vin, unsigned int size_t, unsigned int match);
unsigned int isanyD(double *vin, unsigned int size_t, double match, unsigned int offset);
unsigned int isallUI(unsigned int *vin, unsigned int size_t, unsigned int match, unsigned int offset);
unsigned int isallI(int *vin, unsigned int size_t, int match, unsigned int offset);
unsigned int isallD(double *vin, unsigned int size_t, double match);
double powDUI(double *Va, unsigned int *Vb, unsigned int size_t);
unsigned int sumUI(unsigned int *Vin, unsigned int size_t);
double sumT_D(double **Tin, unsigned int nrow, unsigned int ncol);
unsigned int isanylessD_2D(double **Tin, unsigned int nrow, unsigned int ncol, double match);
void vcopyUI(unsigned int *Vout, unsigned int *Vin, unsigned int size_t);
void vcopyI(int *Vout, int *Vin, unsigned int size_t);
void vcopyUI_I(int *Vout, unsigned int *Vin, unsigned int size_t);
void vcopyD(double *Vout, double *Vin, unsigned int size_t);
void rowsumD(double **Tin, unsigned int nrow, unsigned int ncol, double *Vout);
unsigned int matchUI(unsigned int *Vin, unsigned int size_t, unsigned int match);
void smultI_UI(int *Vout, unsigned int *Vin, unsigned int size_t, int scale);
void vsum_UI_I(unsigned int *Vbase, int *Vadd, unsigned int size_t);
void sselect_UI(unsigned int **Tin, unsigned int *Vout, unsigned int nrow, unsigned int matchcol, unsigned int datcol, unsigned int match, unsigned int dcol, unsigned int deme);
void sselect_UIV(unsigned int **Tin, unsigned int *Vout, unsigned int nrow, unsigned int matchcol, unsigned int datcol, unsigned int *match, unsigned int mlength, unsigned int mlength2, unsigned int dcol, unsigned int deme);
double sumrep(unsigned int n);
double sumrepsq(unsigned int n);
unsigned int maxUI(unsigned int *vin, unsigned int size_t);
unsigned int minUI(unsigned int *vin, unsigned int size_t, unsigned int offset);
int minI(int *vin, unsigned int size_t, unsigned int offset);
unsigned int first_neI(int *vin, unsigned int size_t, int target, unsigned int offset);
unsigned int first_neUI(unsigned int *vin, unsigned int size_t, unsigned int target, unsigned int offset);
unsigned int last_neI(int *vin, unsigned int size_t, int target, unsigned int offset);
unsigned int last_neUI(unsigned int *vin, unsigned int size_t, unsigned int target, unsigned int offset);

unsigned int trig(unsigned int x);
double KQ(double Q);
double singGC(unsigned int y, unsigned int k, double geemi, double geeme, double Qmi, double Qme, double sexC);
double pairGC(unsigned int x, unsigned int k, double geemi, double Qmi);
double P23(unsigned int y, unsigned int k, unsigned int Na);
double P4(unsigned int x, unsigned int k, unsigned int Na);
double P56(unsigned int y, unsigned int Na);
double P7(unsigned int x, unsigned int k, unsigned int Na);
double P8(unsigned int x, unsigned int y, unsigned int k, unsigned int Na);
double P9(unsigned int x, unsigned int y, unsigned int k, double geemi, double geeme, double Qmi, double Qme, double sexC);
double P10(unsigned int x, unsigned int y, unsigned int k, double mee);
double P11(unsigned int y, unsigned int k, double sexC, double ree, unsigned int lrec, unsigned int nlrec, unsigned int nlrec2);
double P12(unsigned int x, unsigned int k, double geemi, double Qmi);
double invs1(double Si, double Qin);
double invt1(double Ti, double Qin);
double invs2(double Si, double Qin);
void probset2(unsigned int N, double gmi, double gme, double *sexC, double rec, double Qmi, double Qme, unsigned int lrec, unsigned int *nlrec, unsigned int *nlrec2, double mig, unsigned int *Nwith, unsigned int *Nbet, unsigned int *kin, unsigned int sw, double **pr);
void rate_change(unsigned int N, unsigned int pST,double pLH, double pHL, double *sexH, double *sexL, unsigned int switch1, double *sexCN, double *sexCNInv, double *tts, unsigned int *npST,const gsl_rng *r);
void stchange2(unsigned int ev, unsigned int deme, unsigned int *kin, int *WCH, int *BCH);
void sexconv(unsigned int **Tin, unsigned int *rsex, unsigned int nsum, unsigned int Ntot, unsigned int Nid, unsigned int ex);
unsigned int coalesce(unsigned int **indvs, int **GType, double **CTms , int **TAnc, unsigned int **nlri, double Ttot, unsigned int *Nwith, unsigned int *Nbet, unsigned int deme, unsigned int *rsex, unsigned int *nsex, unsigned int ex, unsigned int drec, unsigned int e2, unsigned int **breaks, unsigned int nsites, unsigned int *nbreaks, unsigned int NMax, unsigned int Itot, double lambdami, double lambdame, double gmi, double gme, double Qmi, double Qme, unsigned int *gcalt, double *sexC, int *WCHex, int *BCHex, const gsl_rng *r, unsigned int ei);
void cchange(unsigned int **indvs, int **GType, double **CTms, int **TAnc, unsigned int **breaks, unsigned int *csamp, unsigned int *par, unsigned int lsamp, unsigned int Ntot, unsigned int cst, unsigned int *cend, double Ttot, unsigned int isall);
unsigned int ccheck(unsigned int **indvs, int **GType, unsigned int **breaks, unsigned int nsites, unsigned int Ntot, unsigned int nbreaks);
void excoal(unsigned int **indvs, int **GType, unsigned int *par, unsigned int nbreaks, unsigned int npar, unsigned int Ntot, int *WCHex, int *BCHex, unsigned int deme);
unsigned int coalcalc(unsigned int **breaks, unsigned int nsites, unsigned int nbreaks, unsigned int start);
void sexsamp(unsigned int **indvs, unsigned int *rsex, unsigned int *nsex, unsigned int *Nwith, unsigned int Ntot, const gsl_rng *r);
void indv_sort(unsigned int **indvs, unsigned int nrow);
void indv_sortD(double **Tin, unsigned int nrow, unsigned int ncol, unsigned int tcol);
void Wait();
void TestTabs(unsigned int **indvs, int **GType, double **CTms, int **TAnc, unsigned int **breaks, unsigned int **nlri, unsigned int NMax, unsigned int Itot, unsigned int Nbet, unsigned int nbreaks);
char * treemaker(double **TFin, double thetain, unsigned int mind2, unsigned int maxd2, double mind, double maxd, unsigned int Itot, unsigned int run, double gmi, double gme, unsigned int ismsp, unsigned int *nmutT, unsigned int prtrees, unsigned int ismut, double pburst, unsigned int mburst, double bdist, const gsl_rng *r);
void reccal(unsigned int **indvs, int **GType, unsigned int **breaks, unsigned int **nlri, unsigned int *Nbet, unsigned int *Nwith, unsigned int *rsex, unsigned int esex, unsigned int *lnrec, unsigned int nbreaks, unsigned int NMax, unsigned int sw, unsigned int run);
void proberr(unsigned int est, double **pr);
void proberr2(double **pr);
void proberr3(double **pr);
void manyr();
void usage();

/* Global variable declaration */
double rec = 0;				/* Per-site recombination rate */
unsigned int nsites = 1;	/* Number of sites (for recombination) */
unsigned int d = 1;			/* Number of demes */

/* Function to replicate the 'any' func in R (unsigned int) */
unsigned int isanyUI(unsigned int *vin, unsigned int size_t, unsigned int match){
	unsigned int i;
	unsigned int res = 0;
	for(i = 0; i < size_t; i++){
		if(*(vin + i) == match){
			res = 1;
		}
	}
	return res;
}

/* Function to replicate the 'any' func in R (double) */
unsigned int isanyD(double *vin, unsigned int size_t, double match, unsigned int offset){
	unsigned int i;
	unsigned int res = 0;
	for(i = offset; i < size_t; i++){
		if(*(vin + i) == match){
			res = 1;
		}
	}
	return res;
}

/* Function to replicate the 'all' func in R (unsigned int) */
unsigned int isallUI(unsigned int *vin, unsigned int size_t, unsigned int match, unsigned int offset){
	unsigned int i;
	unsigned int res = 1;
	for(i = offset; i < size_t; i++){
		if(*(vin + i) != match){
			res = 0;
		}
	}
	return res;
}

/* Function to replicate the 'all' func in R (normal int) */
unsigned int isallI(int *vin, unsigned int size_t, int match, unsigned int offset){
	unsigned int i;
	unsigned int res = 1;
	for(i = offset; i < size_t; i++){
	if(*(vin + i) != match){
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
		if(*(vin + i) != match){
			res = 0;
		}
	}
	return res;
}

/* Dot product of two vectors (double X UI) */
double powDUI(double *Va, unsigned int *Vb, unsigned int size_t){
	unsigned int i;
	double res = 1;
	for(i = 0; i < size_t; i++){
		res *= pow((*(Va + i)),(*(Vb + i)));
	}
	return res;
}

/* Summing vector (UI) */
unsigned int sumUI(unsigned int *Vin, unsigned int size_t){
	unsigned int i;
	unsigned int res = 0;
	for(i = 0; i < size_t; i++){
		res += *(Vin + i);
	}
	return res;
}

/* Summing entire table (double) */
double sumT_D(double **Tin, unsigned int nrow, unsigned int ncol){
	unsigned int i, j;
	double res = 0;
	for(i = 0; i < nrow; i++){
		for(j = 0; j < ncol; j++){
			res += (*((*(Tin + i)) + j));
		}
	}
	return res;
}

/* Is any entry of table less than input? */
unsigned int isanylessD_2D(double **Tin, unsigned int nrow, unsigned int ncol, double match){
	unsigned int i, j;
	double res = 0;
	for(i = 0; i < nrow; i++){
		for(j = 0; j < ncol; j++){
			if(*((*(Tin + i)) + j) < match){
				res = 1;
			}
		}
	}
	return res;
}

/* Copying vectors (UI) */
void vcopyUI(unsigned int *Vout, unsigned int *Vin, unsigned int size_t){
	unsigned int x;
	for(x = 0; x < d; x++){
		*(Vout + x) = *(Vin + x);
	}
}

/* Copying vectors (Int) */
void vcopyI(int *Vout, int *Vin, unsigned int size_t){
	unsigned int x;
	for(x = 0; x < d; x++){
		*(Vout + x) = *(Vin + x);
	}
}

/* Copying vectors (UInt -> Int) */
void vcopyUI_I(int *Vout, unsigned int *Vin, unsigned int size_t){
	unsigned int x;
	for(x = 0; x < d; x++){
		*(Vout + x) = *(Vin + x);
	}
}

/* Copying vectors (double) */
void vcopyD(double *Vout, double *Vin, unsigned int size_t){
	unsigned int x;
	for(x = 0; x < d; x++){
		*(Vout + x) = *(Vin + x);
	}
}

/* Calculating rowsums (double) */
void rowsumD(double **Tin, unsigned int nrow, unsigned int ncol, double *Vout){
	unsigned int i, j;
	for(i = 0; i < nrow; i++){
		*(Vout + i) = 0;
		for(j = 0; j < ncol; j++){
			*(Vout + i) += *((*(Tin + i)) + j);
		}
	}
}

/* Replicating 'match' R function (UI) */
unsigned int matchUI(unsigned int *Vin, unsigned int size_t, unsigned int match){
	unsigned int i;
	unsigned int res = 0;
	for(i = 0; i < size_t; i++){
		if(*(Vin + i) == match){
			res = i;
		}
	}
	return res;
}

/* Multiplying vector by a scalar (Int) */
void smultI_UI(int *Vout, unsigned int *Vin, unsigned int size_t, int scale){
	unsigned int i;
	for(i = 0; i < size_t; i++){
		*(Vout + i) = (scale)*(*(Vin + i));
	}
}

/* Summing two vector (UI + Int) */
void vsum_UI_I(unsigned int *Vbase, int *Vadd, unsigned int size_t){
	unsigned int i;
	for(i = 0; i < size_t; i++){
		*(Vbase + i) += *(Vadd + i);
	}
}

/* Choosing UNIQUE elements from array that match certain pattern (UI) */
void sselect_UI(unsigned int **Tin, unsigned int *Vout, unsigned int nrow, unsigned int matchcol, unsigned int datcol, unsigned int match, unsigned int dcol, unsigned int deme){
	unsigned int j, x = 0;
	unsigned int count = 0;
	unsigned int ny = 1;	/* 'Not yet' - is element already present? */
	for(j = 0; j < nrow; j++){
		if((*((*(Tin + j)) + matchcol) == match) && (*((*(Tin + j)) + dcol) == deme) ){
			/* Only add if not yet present */
			ny = 1;
			for(x = 0; x < count; x++){
				if( *(Vout + x) == *((*(Tin + j)) + datcol)){
					ny = 0;
				}
			}
			if(ny == 1){
				*(Vout + count) = *((*(Tin + j)) + datcol);
				count++;
			}
		}
	}
}

/* Choosing elements from array that match certain vector (UI) */
void sselect_UIV(unsigned int **Tin, unsigned int *Vout, unsigned int nrow, unsigned int matchcol, unsigned int datcol, unsigned int *match, unsigned int mlength, unsigned int mlength2, unsigned int dcol, unsigned int deme){
	unsigned int j = 0;
	unsigned int count = 0;
	unsigned int count2 = 0;	
	while(count < mlength && count2 < mlength2){
		for(j = 0; j < nrow; j++){
			if((*((*(Tin + j)) + matchcol) == (*(match + count2))) && (*((*(Tin + j)) + dcol) == deme) ){
				*(Vout + 2*count) = *((*(Tin + j)) + datcol);
				*(Vout + 2*count+1) = *((*(Tin + j + 1)) + datcol);
				count++;
				count2++;
				break;
			}else if((*((*(Tin + j)) + matchcol) == *(match + count2)) && (*((*(Tin + j)) + dcol) != deme) ){
				count2++;
				break;
			}
		}
	}
}

/* Sum of 1/i */
double sumrep(unsigned int n){
	unsigned int count = 0;
	unsigned int i;
	for(i = 1; i < n; i++){
		count += 1/(1.0*i);
	}
	return(count);
}

/* Sum of 1/i^2 */
double sumrepsq(unsigned int n){
	unsigned int count = 0;
	unsigned int i;
	for(i = 1; i < n; i++){
		count += 1/(1.0*i*i);
	}
	return(count);
}

/* Finding max of samples (UI) */
unsigned int maxUI(unsigned int *vin, unsigned int size_t){
	unsigned int ret = 0;
	unsigned int i = 0;
	for(i = 0; i < size_t; i++){
		if(*(vin + i) > ret){
			ret = *(vin + i);
		}
	}
	return ret;
}

/* Finding min of samples (UI) */
unsigned int minUI(unsigned int *vin, unsigned int size_t, unsigned int offset){
	unsigned int ret = 0;
	unsigned int i = 0;
	for(i = 0; i < size_t; i++){
		if(*(vin + i + offset) < ret){
			ret = *(vin + i + offset);
		}
	}
	return ret;
}

/* Finding min of samples (I) */
int minI(int *vin, unsigned int size_t, unsigned int offset){
	int ret = 0;
	unsigned int i = 0;
	for(i = 0; i < size_t; i++){
		if(*(vin + i + offset) < ret){
			ret = *(vin + i + offset);
		}
	}
	return ret;
}

/* First element not of type (I) */
unsigned int first_neI(int *vin, unsigned int size_t, int target, unsigned int offset){
	int ret = 0;
	unsigned int i = 0;
	for(i = offset; i < size_t; i++){
		if(*(vin + i) != target){
			ret = i;
			break;
		}
	}
	return ret;
}

/* First element not of type (UI) */
unsigned int first_neUI(unsigned int *vin, unsigned int size_t, unsigned int target, unsigned int offset){
	unsigned int ret = 0;
	unsigned int i = 0;
	for(i = offset; i < size_t; i++){
		if(*(vin + i) != target){
			ret = i;
			break;
		}
	}
	return ret;
}

/* Last element not of type (I) */
unsigned int last_neI(int *vin, unsigned int size_t, int target, unsigned int offset){
	int ret = 0;
	unsigned int i = 0;
	for(i = offset; i < size_t; i++){
		if(*(vin + i) != target){
			ret = i;
		}
	}
	return ret;
}

/* Last element not of type (UI) */
unsigned int last_neUI(unsigned int *vin, unsigned int size_t, unsigned int target, unsigned int offset){
	unsigned int ret = 0;
	unsigned int i = 0;
	for(i = offset; i < size_t; i++){
		if(*(vin + i) != target){
			ret = i;
		}
	}
	return ret;
}

/* 'Triangle function' calculation */
unsigned int trig(unsigned int x){
	return (x*(x-1))/2.0;
}

/* 'K' function, for GC calculations */
double KQ(double Q){
	return 1.0 - ((1.0 - exp(-Q))/(1.0*Q));
}

/* Partial GC prob calculations. For calculating overall GC prob, and subprobs of different events */
double singGC(unsigned int y, unsigned int k, double geemi, double geeme, double Qmi, double Qme, double sexC){
	return y*(geemi*(2.0-KQ(Qmi)) + sexC*geeme*(2.0-KQ(Qme))) + 2*k*(geemi*(2.0-KQ(Qmi)) + geeme*(2.0-KQ(Qme)));
}

double pairGC(unsigned int x, unsigned int k, double geemi, double Qmi){
	return geemi*(2.0-KQ(Qmi))*(x-k);
}

/* Functions of each transition probability calculation */
double P23(unsigned int y, unsigned int k, unsigned int Na){
	/* One of the paired samples is recreated: (x,y) -> (x-k+1, y + 2(k-1)). 
	OR A unique samples coalesces with another (either pre-existing or new): (x,y) -> (x-k, y + 2k - 1) */
	return (k*y)/(1.0*Na) + trig(2.0*k)/(2.0*Na);
}
double P4(unsigned int x, unsigned int k, unsigned int Na){
	/* A paired sample coaleses with new unique sample: (x,y) -> (x-k,y + 2k -1) */
	return ((x-k)*2.0*k)/(1.0*Na);
}
double P56(unsigned int y, unsigned int Na){
	/* Two pre-existing unique samples re-create a paired sample: (x,y) -> (x - k + 1, y + 2(k-1)). 
	OR Two pre-existing unique samples coalesce: (x,y) -> (x - k, y + 2k-1) */
	return trig(y)/(2.0*Na);
}
double P7(unsigned int x, unsigned int k, unsigned int Na){
	/* Two remaining paired samples doubly coalesce asexually: (x,y) -> (x-k-1,y+2k) */
	unsigned int diff = 0;
	diff = x-k;
	return trig(diff)/(1.0*Na);
}
double P8(unsigned int x, unsigned int y, unsigned int k, unsigned int Na){
	/* One of the x - k remaining paired samples can coalesce with a unique sample: (x,y) -> (x-k,y+2k-1) */
	unsigned int diff = 0;
	diff = x-k;
	return (diff*y)/(1.0*Na);
}
double P9(unsigned int x, unsigned int y, unsigned int k, double geemi, double geeme, double Qmi, double Qme, double sexC){
	/* 'Disruptive' gene conversion event (I.e. where at least one bp is in sampled material) */
	return pairGC(x, k, geemi, Qmi) + singGC(y, k, geemi, geeme, Qmi, Qme, sexC);
}
double P10(unsigned int x, unsigned int y, unsigned int k, double mee){
	/* A sample migrates to another deme */
	return mee*(x + y + k);
}
double P11(unsigned int y, unsigned int k, double sexC, double ree, unsigned int lrec, unsigned int nlrec, unsigned int nlrec2){
	/* One of the single samples splits by recombination, creates two new single samples: (x,y) -> (x-k,y+2k+1) */
	return (sexC*rec*((lrec - 1)*(y) - nlrec) + rec*((lrec - 1)*(2*k) - nlrec2));
}
double P12(unsigned int x, unsigned int k, double geemi, double Qmi){
	/* Complete gene conversion, coalesces paired sample into single sample */
	return (geemi*(x-k)*(exp(-Qmi)/(1.0*Qmi)));
}

double invs1(double Si, double Qin){
	/* Inversion of start point, if starting but not ending in tract */
	return log(1 + (exp(Qin)-1)*Si)/(1.0*Qin);
}

double invt1(double Ti, double Qin){
	/* Inversion of end point, if starting outside but ending in tract */
	return (Qin - log(exp(Qin) + Ti - exp(Qin)*Ti))/(1.0*Qin);
}


double invs2(double Si, double Qin){
	/* Inversion of start point, double GC event (2 bps) */
	return exp(-Qin)*(-1.0 + Si + exp(Qin)*(Si*(Qin - 1.0) - gsl_sf_lambert_W0((-1.0)*exp(-Qin + Si*(Qin - 1.0) + exp(-Qin)*(Si-1.0)))))/(1.0*Qin);
}


/* Calculate probability change vectors each time OVER EACH DEME */
void probset2(unsigned int N, double gmi, double gme, double *sexC, double rec, double Qmi, double Qme, unsigned int lrec, unsigned int *nlrec, unsigned int *nlrec2, double mig, unsigned int *Nwith, unsigned int *Nbet, unsigned int *kin, unsigned int sw, double **pr){
	unsigned int x;				/* Deme counter */
	unsigned int ksum = 0;		/* Total number of segregating events */
	
	/* First calculate number of splits... */
	for(x = 0; x < d; x++){
		ksum += *(kin + x);
	}
	
	/* Calculate each transition probability, per deme */
	for(x = 0; x < d; x++){
		*((*(pr + 4)) + x) = P56(*(Nbet + x),N);
		*((*(pr + 5)) + x) = P56(*(Nbet + x),N);
		*((*(pr + 6)) + x) = P7(*(Nwith + x),*(kin + x),N);
		*((*(pr + 7)) + x) = P8(*(Nwith + x),*(Nbet + x),*(kin + x),N);
		*((*(pr + 8)) + x) = P9(*(Nwith + x),*(Nbet + x),*(kin + x),gmi,gme,Qmi,Qme,*(sexC + x));
		*((*(pr + 9)) + x) = P10(*(Nwith + x),*(Nbet + x),*(kin + x),mig);
		*((*(pr + 10)) + x) = P11(*(Nbet + x),*(kin + x),*(sexC + x),rec,lrec,*(nlrec + x),*(nlrec2 + x));
		*((*(pr + 11)) + x) = P12(*(Nwith + x),*(kin + x),gmi,Qmi);
		
		/* Only activate the first three events if need to consider segregation via sex 
		(fourth is 'split pairs remain split') */
		if(sw == 1){
			if(ksum != 1){
				*((*(pr + 1)) + x) = P23(*(Nbet + x),*(kin + x),N);
			}else if(ksum == 1){
				*((*(pr + 1)) + x) = 0;
			}
			*((*(pr + 2)) + x) = P23(*(Nbet + x),*(kin + x),N);
			*((*(pr + 3)) + x) = P4(*(Nwith + x),*(kin + x),N);

			/* Last entry is simply 1-(sum all other probs) */
			if(x == 0){
				*((*(pr + 0)) + x) = (1-sumT_D(pr,12,d));
			}
			
		}else if(sw == 0){
			*((*(pr + 1)) + x) = 0;
			*((*(pr + 2)) + x) = 0;
			*((*(pr + 3)) + x) = 0;
			*((*(pr + 0)) + x) = 0;
		}
	}
}	/* End of 'probset2' function */

/* Function to change rates of sex given a state change */
void rate_change(unsigned int N, unsigned int pST,double pLH, double pHL, double *sexH, double *sexL, unsigned int switch1, double *sexCN, double *sexCNInv, double *tts, unsigned int *npST, const gsl_rng *r){
	unsigned int x = 0;			/* Deme counter */
	
	/* Setting up transition time (tts, or 'time to switch') */
	if(pST == 0){
		for(x = 0; x < d; x++){
			*(sexCN + x) = *(sexL + x);
			*(sexCNInv + x) = 1.0 - (*(sexL + x));
		}
		*tts = gsl_ran_geometric(r,pLH)/(2.0*N*d);
		*npST = 1;
	}else if(pST == 1){
		for(x = 0; x < d; x++){
			*(sexCN + x) = *(sexH + x);
			*(sexCNInv + x) = 1.0 - (*(sexH + x));
		}
		*tts = gsl_ran_geometric(r,pHL)/(2.0*N*d);
		*npST = 0;
	}else if(pST == 2){		/* If stepwise change, alter depending on whether already switched or not */
		*npST = pST;
		if(switch1 == 0){
			for(x = 0; x < d; x++){
				*(sexCN + x) = *(sexL + x);
				*(sexCNInv + x) = 1.0 - (*(sexL + x));
			}
			*tts = pLH;
		}else if(switch1 == 1){
			for(x = 0; x < d; x++){
				*(sexCN + x) = *(sexH + x);
				*(sexCNInv + x) = 1.0 - (*(sexH + x));
			}
			*tts = (1.0/0.0);
		}
	}else if(pST == 3){		/* If constant, no time to sex switch */
		for(x = 0; x < d; x++){
				*(sexCN + x) = *(sexL + x);
				*(sexCNInv + x) = 1.0 - (*(sexL + x));
			}
		*tts = (1.0/0.0);
		*npST = pST;
	}
	
}	/* End of 'rate_change' function */

/* Function to determine how to change state numbers following an event,
taking into account events over all demes*/
void stchange2(unsigned int ev, unsigned int deme, unsigned int *kin, int *WCH, int *BCH){

	int *oo3 = calloc(2,sizeof(int));		/* Extra change in pop due to event */
	int *negk = calloc(d,sizeof(int));		/* Negative of k */
	int *dblek = calloc(d,sizeof(int));		/* Double k */
	
	/* Rescaling k */
	smultI_UI(negk, kin, d, (-1));
	smultI_UI(dblek, kin, d, 2);
	
	/* Baseline sex events */
	vcopyI(WCH, negk, d);
	vcopyI(BCH, dblek, d);
	
	/* Now deciding extra events depending on deme and event */
	switch(ev)
	{
		case 0:
			*(oo3 + 0) = 0;
			*(oo3 + 1) = 0;
			break;
		case 1:
			*(oo3 + 0) = 1;
			*(oo3 + 1) = -2;
			break;
		case 2:
			*(oo3 + 0) = 0;
			*(oo3 + 1) = -1;
			break;
		case 3:
			*(oo3 + 0) = 0;
			*(oo3 + 1) = -1;
			break;
		case 4:
			*(oo3 + 0) = 1;
			*(oo3 + 1) = -2;
			break;
		case 5:
			*(oo3 + 0) = 0;
			*(oo3 + 1) = -1;
			break;
		case 6:
			*(oo3 + 0) = -1;
			*(oo3 + 1) = 0;
			break;
		case 7:
			*(oo3 + 0) = 0;
			*(oo3 + 1) = -1;
			break;
		case 8:
			*(oo3 + 0) = 0;
			*(oo3 + 1) = 0;
			break;
		case 9:
			*(oo3 + 0) = 0;
			*(oo3 + 1) = 0;
			break;
		case 10:
			*(oo3 + 0) = 1;
			*(oo3 + 1) = -1;
			break;
		case 11:
			*(oo3 + 0) = -1;
			*(oo3 + 1) = 1;
			break;			
		default:	/* If none of these cases chosen, exit with error message */
			fprintf(stderr,"Error: Non-standard coalescent case selected ('stchange2').\n");
			exit(1);
            break;
	}
	
	*(WCH + deme) += *(oo3 + 0);
	*(BCH + deme) += *(oo3 + 1);	
	
	free(dblek);
	free(negk);
	free(oo3);
	
}	/* End of 'stchange' function */

/* For converting WH to BH */
void sexconv(unsigned int **Tin, unsigned int *rsex, unsigned int nsum, unsigned int Ntot, unsigned int Nid, unsigned int ex){
	unsigned int j;
	unsigned int count = 0;
	unsigned int npc = 0;	
	while(count < nsum){
		for(j = 0; j < Ntot; j++){
			if( *((*(Tin + j)) + 1) == *(rsex + count) ){
				*((*(Tin + j)) + 2) = 1;
				*((*(Tin + j + 1)) + 2) = 1;				/* Since other paired sample also split */
				if(ex == 8 || ex == 10){
					*((*(Tin + j + 1)) + 1) = Nid + npc;		/* Placing sample in new individual */
					npc++;
				}
				count++;
				break;
			}
		}
	}
}

/* Function to change status of samples following event change */
unsigned int coalesce(unsigned int **indvs, int **GType, double **CTms , int **TAnc, unsigned int **nlri, double Ttot, unsigned int *Nwith, unsigned int *Nbet, unsigned int deme, unsigned int *rsex, unsigned int *nsex, unsigned int ex, unsigned int drec, unsigned int e2, unsigned int **breaks, unsigned int nsites, unsigned int *nbreaks, unsigned int NMax, unsigned int Itot, double lambdami, double lambdame, double gmi, double gme, double Qmi, double Qme, unsigned int *gcalt, double *sexC, int *WCHex, int *BCHex, const gsl_rng *r, unsigned int ei){
	
	unsigned int NWtot = 2*sumUI(Nwith,d);
	unsigned int NBtot = sumUI(Nbet,d);
	unsigned int Ntot = NWtot + NBtot;
	unsigned int Nindv = sumUI(Nwith,d) + sumUI(Nbet,d);
	unsigned int nsum = sumUI(nsex,d);
	unsigned int j = 0;
	unsigned int a = 0;
	int x = 0;
	unsigned int count = 0;		/* For converting WH to BH samples */
	unsigned int done = 0;		/* For sampling right individual */
	unsigned int rands = 0;		/* Sample split by rec (event 10) */
	unsigned int rands2 = 0;	/* Sample that does not split fully (event 1; also used in event 8) */
	unsigned int nos = 0;		/* Sub-sample that does not split fully (event 1) */
	unsigned int bhc = 0;		/* Single sample that repairs (ev 1) or migrates (ev 9); */
	unsigned int csamp = 0;		/* Sample that coalesces */
	unsigned int par = 0;		/* Parental sample in coalescence */
	unsigned int par2 = 0;		/* Parental sample of paired coalescence (ev 2) */
	unsigned int WHsel = 0;		/* WH sample involved in event (ev 7) */
	unsigned int isWH = 0;		/* Is the sample from WH? (ev 7) */
	unsigned int parNo = 0;		/* Parent where coalescent occurs (event 6) */
	unsigned int yesrec = 0;	/* Has a suitable recombination site been chosen? */
	unsigned int rsite = 0;		/* Position of recombination breakpoint (event 10) */
	unsigned int isyetbp = 0;	/* Is breakpoint already present? (event 10) */
	unsigned int isyetbp2 = 0;	/* Is 2nd breakpoint already present? (event 8) */	
/*	unsigned int isbpend = 0;	 Is bp at end of table? (event 10) */
	unsigned int maxtr = 0;		/* Max site in bp table before breakpoint (event 10) */
	unsigned int mintr = 0;		/* Start of GC event (event 8) */	
	unsigned int gt = 0;		/* GC acting on single or paired sample? (event 8) */
	unsigned int gcst = 0;		/* GC start point (event 8) */
	int gcend = 0; 				/* GC end point (event 8) */
	unsigned int gcsamp = 0;	/* Index of GC'ed sample (event 8) */
	unsigned int gcsamp2 = 0;	/* Index of GC'ed sample if paired sample involved (event 8) */	
	unsigned int gcln = 0;		/* Length of GC event (event 8) */
	double NWd = 0; 			/* Weighted WH sample chosen (event 8) */
	double NTd = 0; 			/* Total GC prob (ev 8) */	
	double lambda = 0;			/* Assignment of lambda after choosing GC type (event 8) */
	double gcMI = 0;			/* Probability that single samp GC event is mitotic (event 8) */
	double Qin = 0;				/* Q value used in subsequent calcs (event 8) */
	double gcst2 = 0;			/* Initial start site for GC (event 8) */
	double gcend2 = 0;			/* Initial end site for GC (event 8) */
	double p1bp = 0;			/* Prob 1 breakpoint (event 8) */
	unsigned int gcst3 = 0;		/* */
	unsigned int gcS = 0;		/* Type of GC evening on unpaired samples (event 8) */
	unsigned int mindr = 0;		/* Done choosing min tract point? (event 8) */
	unsigned int proceed = 0;	/* Proceed with gene conversion? (event 8) */
	unsigned int maxck = 0;		/* Check if end of bp correctly chosen (event 8) */
	unsigned int iscoal = 0;	/* Check if paired GC event leads to coalescence (event 8) */
	unsigned int rpar = 0;		/* Indv number of recombined sample (event 10) */
	unsigned int gcbp = 0;		/* Type of GC event (1 or 2 breakpoints; event 8) */
	unsigned int achange = 0;	/* Note whether to run extra coal check */
	unsigned int length = 0;	/* Length of sampled chrom (event 10) */
	unsigned int start = 0;		/* Start of sampled chrom  (event 10) */
	
	/* Then further actions based on other event */
	switch(ex)
	{
		case 0:		/* Event 0: 2k new samples created from paired samples. Nothing else to be done */
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);
			break;
		case 1:		/* Event 1: One of the paired samples is recreated, no coalescence */
			/* First choose sample that does not split fully */
			while(done == 0){
				gsl_ran_choose(r,&rands2,1,rsex,nsum,sizeof(unsigned int));
				/* Then checking it is in the same deme as the action */
				for(j = 0; j < Ntot; j++){
					if( *((*(indvs + j)) + 1) == rands2 ){
						if(*((*(indvs + j)) + 3) == deme){
							done = 1;
						}
					}
				}
			}
			
			/* Then setting BH samples */
			nos = gsl_ran_bernoulli(r,0.5);	/* Only sample that splits fully */
			while(count < nsum){
				for(j = 0; j < Ntot; j++){
					if( *((*(indvs + j)) + 1) == *(rsex + count) ){
						if(*(rsex + count) != rands2){
							*((*(indvs + j)) + 2) = 1;
							*((*(indvs + j + 1)) + 2) = 1;
							*((*(indvs + j + 1)) + 1) = (Nindv + count);
							count++;
						}else if(*(rsex + count) == rands2){
							*((*(indvs + j + nos)) + 2) = 1;
							*((*(indvs + j + nos)) + 1) = (Nindv + count);
							count++;						
						}
						break;
					}
				}
			}
			
			/* Now; out of all single samples, choose one to rebind with the single sample */
			unsigned int *singsamps = calloc((*(Nbet + deme) + 2*(*(nsex + deme)) - 1),sizeof(unsigned int));			/* For storing BH samples */
			sselect_UI(indvs, singsamps, Ntot, 2, 0, 1, 3, deme);
			gsl_ran_choose(r,&bhc,1,singsamps,(*(Nbet + deme) + 2*(*(nsex + deme)) - 1),sizeof(unsigned int));
			for(j = 0; j < Ntot; j++){
				if( *((*(indvs + j)) + 0) == bhc ){
					*((*(indvs + j)) + 2) = 0;
					*((*(indvs + j)) + 1) = rands2;	/* Ensuring paired samples have same parents*/
					break;
				}
			}
			
			free(singsamps);
			break;
		case 2:		/* Event 2: One of the unique samples coaleses with another unique one (either pre-existing or new) */
			done = 0;
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);
			
			unsigned int *singsamps2 = calloc((*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(unsigned int));					/* For storing BH samples */
			sselect_UI(indvs, singsamps2, Ntot, 2, 0, 1, 3, deme);
			
			while(done == 0){
				gsl_ran_choose(r,&csamp,1,singsamps2,(*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(unsigned int));			/* One sample involved in coalescence (csamp) */
				par = csamp;
				while(par == csamp){	/* Ensuring par != csamp */
					gsl_ran_choose(r,&par,1,singsamps2,(*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(unsigned int));			/* Other sample involved in coalescence (par) */
				}
		
				/* Now checking that at least one of the two is from sample split by sex */
				for(j = 0; j < Ntot; j++){
					if( (*((*(indvs + j)) + 0) == csamp) ||  (*((*(indvs + j)) + 0) == par) ){
						for(a = 0; a < nsum; a++){
							if( (*((*(indvs + j)) + 1)) == *(rsex + a)){
								done = 1;
								break;
							}
						}
					}
				}
			}
			
			/* Now updating coalescent times */
			cchange(indvs, GType, CTms, TAnc, breaks, &csamp, &par, 1, Ntot, 0, nbreaks, Ttot, 1);
			
			/* Check if tracts have coalesced */
			achange = ccheck(indvs,GType,breaks,nsites,Ntot,*nbreaks);
			
			if(achange == 1){
				excoal(indvs, GType, &par, *nbreaks, 1, Ntot, WCHex, BCHex, deme);
			}
			
			free(singsamps2);
			break;
		case 3:		/* Event 3: A paired sample coaleses with new unique sample */
			
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);
			unsigned int *singsamps3N = calloc(2*(*(nsex + deme)),sizeof(unsigned int));			/* For storing new BH samples */
			unsigned int *singsamps3B = calloc( 2*((*(Nwith + deme) - (*(nsex + deme)))) ,sizeof(unsigned int));			/* For storing WH samples */
			unsigned int *twosamps = calloc(2,sizeof(unsigned int));		/* Two samples in coalescence */
				
			/* Creating vector of new unique samples created */
			sselect_UIV(indvs, singsamps3N, Ntot, 1, 0, rsex, 2*(*(nsex + deme)), nsum, 3, deme);
			/* Creating vector of existing paired samples */
			sselect_UI(indvs, singsamps3B, Ntot, 2, 0, 0, 3, deme);
			
			gsl_ran_choose(r,&twosamps[0],1,singsamps3N,(2*(*(nsex + deme))),sizeof(unsigned int));			/* Unique sample involved in coalescence */
			gsl_ran_choose(r,&twosamps[1],1,singsamps3B,(2*((*(Nwith + deme) - (*(nsex + deme))))),sizeof(unsigned int));			/* Paired sample involved in coalescence */
			gsl_ran_choose(r,&csamp,1,twosamps,2,sizeof(unsigned int));			/* One sample involved in coalescence (csamp) */
			par = csamp;
			while(par == csamp){	/* Ensuring par != csamp */
				gsl_ran_choose(r,&par,1,twosamps,2,sizeof(unsigned int));			/* Other sample involved in coalescence (par) */
			}
		
			/* Correction if parential sample is unique sample */
			unsigned int pcorr = 0;
			for(j = 0; j < Ntot; j++){
				if( (*((*(indvs + j)) + 0) == par) && (*((*(indvs + j)) + 2) == 1) ){
					pcorr = 1;
					par2 = j;
					(*((*(indvs + j)) + 2) = 0);
					break;
				}
			}
			if(pcorr == 1){
				for(j = 0; j < Ntot; j++){
					if(*((*(indvs + j)) + 0) == csamp){
						*((*(indvs + par2)) + 1) = *((*(indvs + j)) + 1);
						break;
					}
				}
			}
		
			/* Now updating coalescent times */
			cchange(indvs, GType, CTms, TAnc, breaks, &csamp, &par, 1, Ntot, 0, nbreaks, Ttot, 1);
			
			/* Check if tracts have coalesced */
			achange = ccheck(indvs,GType,breaks,nsites,Ntot,*nbreaks);
			
			if(achange == 1){
				excoal(indvs, GType, &par, *nbreaks, 1, Ntot, WCHex, BCHex, deme);
			}
			
			free(twosamps);
			free(singsamps3B);
			free(singsamps3N);

			break;
		case 4:			/* Event 4: Two pre-existing unique samples re-create paired sample. */
			done = 0;
			unsigned int *singsamps4 = calloc(*(Nbet + deme),sizeof(unsigned int));			/* For storing BH samples */
			unsigned int *twosing = calloc(2,sizeof(unsigned int));	

			sselect_UI(indvs, singsamps4, Ntot, 2, 0, 1, 3, deme);

			/* Two sample involved in pairing */
			gsl_ran_choose(r,twosing,2,singsamps4,*(Nbet + deme),sizeof(unsigned int));
			gsl_ran_shuffle(r, twosing, 2, sizeof(unsigned int));
			
			/* Now going through (twice!) and updating ancestry */
			for(j = 0; j < Ntot; j++){
				if(*((*(indvs + j)) + 0) == *(twosing + 0)){
					par2 = j;
					*((*(indvs + j)) + 2) = 0;
					break;
				}
			}
			
			for(j = 0; j < Ntot; j++){
				if(*((*(indvs + j)) + 0) == *(twosing + 1)){
					*((*(indvs + j)) + 2) = 0;
					*((*(indvs + j)) + 1) = *((*(indvs + par2)) + 1);
					break;
				}
			}
			
			/* THEN convert WH to BH samples */
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);
			
			free(twosing);
			free(singsamps4);
			break;
		case 5:			/* Event 5: Two pre-existing unique samples coalesce */	
			csamp = 0;
			unsigned int *singsamps5 = calloc(*(Nbet + deme),sizeof(unsigned int));			/* For storing BH samples */
			sselect_UI(indvs, singsamps5, Ntot, 2, 0, 1, 3, deme);
			
			gsl_ran_choose(r,&csamp,1,singsamps5,(*(Nbet + deme)),sizeof(unsigned int));			/* One sample involved in coalescence (csamp) */
			par = csamp;
			while(par == csamp){	/* Ensuring par != csamp */
				gsl_ran_choose(r,&par,1,singsamps5,(*(Nbet + deme)),sizeof(unsigned int));			/* Other sample involved in coalescence (par) */
			}
/*			printf("Csamp, par is %d %d\n",csamp,par);*/
		
			/* Now updating coalescent times */
			cchange(indvs, GType, CTms, TAnc, breaks, &csamp, &par, 1, Ntot, 0, nbreaks, Ttot, 1);
			
			/* Check if tracts have coalesced */
			achange = ccheck(indvs,GType,breaks,nsites,Ntot,*nbreaks);
			
			if(achange == 1){
				excoal(indvs, GType, &par, *nbreaks, 1, Ntot, WCHex, BCHex, deme);
			}
		
			/* THEN convert WH to BH samples */
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);
			
			free(singsamps5);
			break;
		case 6:			/* Event 6: Two remaining paired samples doubly coalesce asexually. */
			/* Converting WH to BH samples */
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);
			
			unsigned int *parsamps6 = calloc((*(Nwith + deme) - *(nsex + deme)),sizeof(unsigned int));			/* For storing WH indvs */
			unsigned int *twopars6 = calloc(2,sizeof(unsigned int));
			unsigned int *lhs = calloc(2,sizeof(unsigned int));
			unsigned int *rhs = calloc(2,sizeof(unsigned int));			
			unsigned int *csamp2 = calloc(2,sizeof(unsigned int));
			unsigned int *parT = calloc(2,sizeof(unsigned int));
				
			sselect_UI(indvs, parsamps6, Ntot, 2, 1, 0, 3, deme);
			
			/* Two parents involved in coalescence */
			gsl_ran_choose(r,twopars6,2,parsamps6,(*(Nwith + deme) - *(nsex + deme)),sizeof(unsigned int));
			gsl_ran_shuffle(r,twopars6, 2, sizeof(unsigned int));
			
			/* Finding parents and partitioning samples (twice) */
			for(j = 0; j < Ntot; j++){
				if(*((*(indvs + j)) + 1) == *(twopars6 + 0)){
					*(lhs + 0) = *((*(indvs + j)) + 0);
					*(rhs + 0) = *((*(indvs + j + 1)) + 0);
					break;
				}
			}
			for(j = 0; j < Ntot; j++){
				if(*((*(indvs + j)) + 1) == *(twopars6 + 1)){
					*(lhs + 1) = *((*(indvs + j)) + 0);
					*(rhs + 1) = *((*(indvs + j + 1)) + 0);
					break;
				}
			}
			
			/* Now assigning coalesced and parental samples respectively */
			gsl_ran_choose(r,&csamp2[0],1,lhs,2,sizeof(unsigned int));
			parT[0] = csamp2[0];
			while(parT[0] == csamp2[0]){
				gsl_ran_choose(r,&parT[0],1,lhs,2,sizeof(unsigned int));
			}
			
			gsl_ran_choose(r,&csamp2[1],1,rhs,2,sizeof(unsigned int));
			parT[1] = csamp2[1];
			while(parT[1] == csamp2[1]){
				gsl_ran_choose(r,&parT[1],1,rhs,2,sizeof(unsigned int));
			}
			
			/* Now updating coalescent times */
			cchange(indvs, GType, CTms, TAnc, breaks, csamp2, parT, 2, Ntot, 0, nbreaks, Ttot, 1);
			
			/* Check if tracts have coalesced */
			achange = ccheck(indvs,GType,breaks,nsites,Ntot,*nbreaks);
			
			/* Making sure parent samples are in same individual */
			for(j = 0; j < Ntot; j++){
				if(*((*(indvs + j)) + 0) == parT[0]){
					parNo = *((*(indvs + j)) + 1);
					break;
				}
			}
			for(j = 0; j < Ntot; j++){
				if(*((*(indvs + j)) + 0) == parT[1]){
					*((*(indvs + j)) + 1) = parNo;
					break;
				}
			}
			
			if(achange == 1){
				excoal(indvs, GType, parT, *nbreaks, 2, Ntot, WCHex, BCHex, deme);
			}
			
			free(parT);
			free(csamp2);
			free(rhs);
			free(lhs);
			free(twopars6);
			free(parsamps6);
			break;
		case 7: 	/* Event 7: One of the x - k remaining paired samples coalesces with a unique sample */
		
			/* Converting WH to BH samples 
			(idea being that I later check if chosen BH originated from WH 
			- discard if so) */
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);

			/* For storing WH indvs */
			unsigned int *parsamps7 = calloc((*(Nwith + deme) - *(nsex + deme)),sizeof(unsigned int));
			unsigned int *singsamps7 = calloc((*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(unsigned int));			/* For storing BH samples */
			unsigned int *twosamps7 = calloc(2,sizeof(unsigned int));
			
			sselect_UI(indvs, parsamps7, Ntot, 2, 1, 0, 3, deme);
			sselect_UI(indvs, singsamps7, Ntot, 2, 0, 1, 3, deme);
						
			/* A paired sample involved in coalescence */
			gsl_ran_choose(r,&WHsel,1,parsamps7,(*(Nwith + deme) - *(nsex + deme)),sizeof(unsigned int));
			nos = gsl_ran_bernoulli(r,0.5);		/* Which side involved in event */
			/* Finding sample */
			for(j = 0; j < Ntot; j++){
				if(*((*(indvs + j)) + 1) == WHsel){
					*(twosamps7 + 0) = *((*(indvs + j + nos)) + 0);
					break;
				}
			}
/*			printf("Nsum, WHSel, nos, samp1 are %d %d %d %d\n",nsum,WHsel,nos,*(twosamps7 + 0));*/
			
			/* Choosing BH sample */
			while(done == 0){	
				gsl_ran_choose(r,&twosamps7[1],1,singsamps7,(*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(unsigned int));
				/* Checking that it didn't originate from WH */
				for(j = 0; j < Ntot; j++){
					if(*((*(indvs + j)) + 0) == *(twosamps7 + 1)){
						/* ACTIONS */
						isWH = 1;
						for(a = 0; a < nsum; a++){
							if( *((*(indvs + j)) + 1) == *(rsex + a)){
								isWH = 0;
/*								printf("Rechange\n");*/
							}
						}
						if(isWH == 1){
							done = 1;
						}
						break;
					}
				}
			}
/*			printf("Samp2 is %d\n",*(twosamps7 + 1));*/
			
			/* Choosing csamp, par */
			gsl_ran_choose(r,&csamp,1,twosamps7,2,sizeof(unsigned int));			/* One sample involved in coalescence (csamp) */
			par = csamp;
			while(par == csamp){	/* Ensuring par != csamp */
				gsl_ran_choose(r,&par,1,twosamps7,2,sizeof(unsigned int));			/* Other sample involved in coalescence (par) */
			}
/*			printf("Csamp, par are %d %d\n",csamp,par);*/
			
			/* Correction if parental sample is BH */
			if(par == *(twosamps7 + 1)){
				for(j = 0; j < Ntot; j++){
					if(*((*(indvs + j)) + 0) == par){
						*((*(indvs + j)) + 2) = 0;
						*((*(indvs + j)) + 1) = WHsel;		/* Same parent as WH sample */
						break;
					}
				}
			}
			
			/* Now updating coalescent times */
			cchange(indvs, GType, CTms, TAnc, breaks, &csamp, &par, 1, Ntot, 0, nbreaks, Ttot, 1);
			
			/* Check if tracts have coalesced */
			achange = ccheck(indvs,GType,breaks,nsites,Ntot,*nbreaks);
			
			if(achange == 1){
				excoal(indvs, GType, &par, *nbreaks, 1, Ntot, WCHex, BCHex, deme);
			}
		
			free(twosamps7);
			free(singsamps7);
			free(parsamps7);
			break;
		case 8:		/* Event 8: (Partial) Gene conversion */

			/* Converting WH to BH samples */
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);
			
			/* First, is it a paired or single sample that is affected? */
			NWd = pairGC(*(Nwith + deme), *(nsex + deme), gmi, Qmi);
			NTd = NWd + singGC(*(Nbet + deme), *(nsex + deme), gmi, gme, Qmi, Qme, *(sexC + deme));
			gt = gsl_ran_bernoulli(r,(NWd/(1.0*NTd)));
			
			if(gt == 0){	/* Acts on single sample */
				unsigned int *singsamps8 = calloc((*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(unsigned int));
				
				/* Obtaining list of samples to choose from */
				sselect_UI(indvs, singsamps8, Ntot, 2, 0, 1, 3, deme);
				
				/* Sample that undergoes GC */
				gsl_ran_choose(r,&rands,1,singsamps8,(*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(unsigned int));
				for(j = 0; j < NMax; j++){
					if( *((*(GType + j)) + 0) == rands){
						gcsamp = j;
						break;
					}
				}
				
				/* Then choosing type of event */
				gcMI = ((gmi*(2.0-KQ(Qmi))*((*(Nbet + deme)) + 2*(*(nsex + deme))))/(1.0*singGC(*(Nbet + deme), *(nsex + deme), gmi, gme, Qmi, Qme, *(sexC + deme))));
				gcS = gsl_ran_bernoulli(r,gcMI);
				if(gcS == 0){
					lambda = lambdame;
					Qin = Qme;
				}else if(gcS == 1){
					lambda = lambdami;
					Qin = Qmi;					
				}
				free(singsamps8);
			}else if(gt == 1){
				lambda = lambdami;
				Qin = Qmi;
				unsigned int *parsamps8 = calloc((*(Nwith + deme) - *(nsex + deme)),sizeof(unsigned int));
				
				/* Obtaining list of samples to choose from */
				sselect_UI(indvs, parsamps8, Ntot, 2, 1, 0, 3, deme);
				
				/* Sample that undergoes GC */
				gsl_ran_choose(r,&rands,1,parsamps8,(*(Nwith + deme) - *(nsex + deme)),sizeof(unsigned int));
				rands2 = gsl_ran_bernoulli(r,0.5);
				for(j = 0; j < Ntot; j++){
					if( *((*(indvs + j)) + 1) == rands){
						gcsamp = *((*(indvs + j + rands2)) + 0);
						gcsamp2 = *((*(indvs + j + (rands2 + 1)%2)) + 0);
						break;
					}
				}
				free(parsamps8);
			}
			
			/* Now determining type of GC event. */
			/* Does there exist one or two breakpoints? */
			if(nsites > 1){
				p1bp = 2.0*((1-KQ(Qin))/(1.0*(2.0-KQ(Qin))));
				gcbp = gsl_ran_bernoulli(r,p1bp);
			
				if(gcbp == 0){			/* Two breakpoints */
					if(nsites > 2){
						gcst = 0;
						gcend = nsites;
						while(gcst == 0 || gcst >= (nsites-1)){
							gcst2 = gsl_ran_flat(r, 0, 1);
							gcst = (unsigned int)(invs2(gcst2,Qin)*nsites);
						}
						while(gcend >= nsites){
							gcln = 0;
							while(gcln == 0){
								gcln = (gsl_ran_geometric(r,(1.0/lambda)));
							}
							gcend = gcst + gcln;
						}
					}else if(nsites == 2){
						gcst = 0;
						gcend = nsites;
					}
				}else if(gcbp == 1){ 	/* One breakpoint */
					gcst3 = gsl_ran_bernoulli(r,0.5);
					if(gcst3 == 0){		/* Starts outside, ends in */
						gcst = 0;
						gcend = 0;
						while(gcend == 0 || gcend == nsites){
							gcend2 = gsl_ran_flat(r, 0, 1);
							gcend = (unsigned int)(invt1(gcend2,Qin)*nsites);
						}
					}else if(gcst3 == 1){		/* Starts inside, ends out */
						gcend = nsites;
						gcst = 0;
						while(gcst == 0 || gcst == nsites){
							gcst2 = gsl_ran_flat(r, 0, 1);
							gcst = (unsigned int)(invs1(gcst2,Qin)*nsites);
						}
					}
				}
			}else if(nsites == 1){
				gcst = 0;
				gcend = nsites;				
			}
			
			/* Finding block number associated with start point */
			mindr = 0;
			for(x = 0; x < *nbreaks; x++){
				if( *((*(breaks + 0)) + x) == gcst){
					isyetbp = 1;
				}
				if( *((*(breaks + 0)) + x) > gcst){
					mintr = (x-1);
					mindr = 1;
					break;
				}
			}
			
			if(mindr == 0){		/* GC start past end of current breakpoints, so need to force it */
				mintr = ((*nbreaks)-1);
			}
			
			/* Inserting new, start BP if needed */
			if((isyetbp != 1) && (*((*(GType + gcsamp)) + mintr + 1) != (-1))){
				(*nbreaks)++;
				for(x = *nbreaks-2; x >= (int)(mintr+1); x--){
					*((*(breaks + 0)) + x + 1) = *((*(breaks + 0)) + x);
					*((*(breaks + 1)) + x + 1) = *((*(breaks + 1)) + x);						
				}
				*((*(breaks + 0)) + mintr + 1) = gcst;
				*((*(breaks + 1)) + mintr + 1) = 0;
				/* Adding new site to genotype; coalescent time; ancestry table */
				for(j = 0; j < NMax; j++){
					for(x = (*nbreaks-1); x > (int)(mintr); x--){
						*((*(GType + j)) + x + 1) = *((*(GType + j)) + x);
						if(j < Itot){
							*((*(CTms + j)) + x + 1) = *((*(CTms + j)) + x);
							*((*(TAnc + j)) + x + 1) = *((*(TAnc + j)) + x);
						}
					}
				}
				mintr++;
			}
			
			/* Finding block number associated with end point */
			maxck = 0;
			for(x = 0; x < *nbreaks; x++){
				if( *((*(breaks + 0)) + x) == gcend){
					isyetbp2 = 1;
				}
				if( *((*(breaks + 0)) + x) > gcend){
					maxtr = (x-1);
					maxck = 1;
					break;
				}
			}

			if(maxck == 0 && (gcend != nsites) && (gcend != 0)){ 		/* If endpoint not found amongst existing bps, must be at end */
				maxtr = ((*nbreaks)-1);
			}else if(maxck == 0 && (gcend == nsites)){
				maxtr = (*nbreaks);
				isyetbp2 = 1;
			}
			
			/* Inserting new, end BP if needed */
			if((isyetbp2 != 1) && (*((*(GType + gcsamp)) + maxtr + 1) != (-1))){
				(*nbreaks)++;
				for(x = *nbreaks-2; x >= (int)(maxtr + 1); x--){
					*((*(breaks + 0)) + x + 1) = *((*(breaks + 0)) + x);
					*((*(breaks + 1)) + x + 1) = *((*(breaks + 1)) + x);						
				}
				*((*(breaks + 0)) + maxtr + 1) = gcend;
				*((*(breaks + 1)) + maxtr + 1) = 0;
				/* Adding new site to genotype; coalescent time; ancestry table */
				for(j = 0; j < NMax; j++){
					for(x = (*nbreaks-1); x > (int)(maxtr); x--){
						*((*(GType + j)) + x + 1) = *((*(GType + j)) + x);
						if(j < Itot){
							*((*(CTms + j)) + x + 1) = *((*(CTms + j)) + x);
							*((*(TAnc + j)) + x + 1) = *((*(TAnc + j)) + x);
						}
					}
				}
				maxtr++;
			}

			/* ONLY PROCEED (in single-sample case) IF NOT ALL SITES EMPTY */
			if(gt == 0 && (isallI((*(GType + gcsamp)), (maxtr + 1), (-1), (mintr + 1)) != 1) && ((isallI((*(GType + gcsamp)), (mintr + 1), (-1), 1) != 1) || (isallI((*(GType + gcsamp)), (*nbreaks + 1), (-1), (maxtr + 1)) != 1))){
				proceed = 1;
			}
		
			if(proceed == 1){
			
				*gcalt = 1;
				
				/* First, checking parent number of existing sample, setting as paired */
				for(j = 0; j < Ntot; j++){
					if( *((*(indvs + j)) + 0) == gcsamp){
						rpar = *((*(indvs + j)) + 1);
						*((*(indvs + j)) + 2) = 0;
						break;
					}
				}
			
				/* Adding new sample to indv table */
				*((*(indvs + NMax)) + 0) = NMax;
				*((*(indvs + NMax)) + 1) = rpar;
				*((*(indvs + NMax)) + 2) = 0;
				*((*(indvs + NMax)) + 3) = deme;
		
				/* Now creating the new sample genotype; updating all other tables */
				for(x = (maxtr); x > (int)(mintr); x--){
					*((*(GType + NMax)) + x) = *((*(GType + gcsamp)) + x);
					*((*(GType + gcsamp)) + x) = (-1);
				}
		
				for(x = (*nbreaks); x > (int)(maxtr); x--){
					*((*(GType + NMax)) + x) = (-1);
				}
				for(x = (mintr); x > (int)0; x--){
					*((*(GType + NMax)) + x) = (-1);
				}
			
				*((*(GType + NMax)) + 0) = NMax;				
			}
			
			/* With paired GC, we have part coalescence instead */
			if(gt == 1){
			
				csamp = gcsamp;
				par = gcsamp2;

				if(
				( (isallI((*(GType + gcsamp)), (mintr+1), (-1), 1) == 1) )
				&& ( (isallI((*(GType + gcsamp)), (*nbreaks+1), (-1), (maxtr+1)) == 1) )
				){
					iscoal = 1;
				}
				
				if(iscoal == 1){
					*gcalt = 2;
					/* Making sure remaining sample becomes BH */
					for(j = 0; j < Ntot; j++){
						if( *((*(indvs + j)) + 0) == gcsamp2){
							*((*(indvs + j)) + 2) = 1;
							break;
						}
					}
				}
				
				/* Now updating coalescent times */
				cchange(indvs, GType, CTms, TAnc, breaks, &csamp, &par, 1, Ntot, mintr, &maxtr, Ttot, iscoal);
			
				/* Check if tracts have coalesced */
				achange = ccheck(indvs,GType,breaks,nsites,Ntot,*nbreaks);
				
				if(achange == 1){
					excoal(indvs, GType, &par, *nbreaks, 1, Ntot, WCHex, BCHex, deme);
				}
				
			}
			
			break;
		case 9:		/* Event 9: Migration of a sample */
		
			/* Converting WH to BH samples  */
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);
			
			if(e2 == 0){		/* WH sample migrates */
				/* For storing WH indvs */
				unsigned int *parsamps9 = calloc((*(Nwith + deme) + 1),sizeof(unsigned int));
				
				/* Obtaining list of samples to choose from */
				sselect_UI(indvs, parsamps9, Ntot, 2, 1, 0, 3, deme);
				
				/* Sample that migrates */
				gsl_ran_choose(r,&bhc,1,parsamps9,(*(Nwith + deme) + 1),sizeof(unsigned int));
				
				for(j = 0; j < Ntot; j++){
					if(*((*(indvs + j)) + 1) == bhc){
						*((*(indvs + j)) + 3) = drec;
						*((*(indvs + j + 1)) + 3) = drec;		/* Altering deme accordingly */
						break;
					}
				}
				
				free(parsamps9);
				
			}else if(e2 == 1){	/* BH sample migrates */
				/* For storing BH samples */
				unsigned int *singsamps9 = calloc((*(Nbet + deme) + 1),sizeof(unsigned int));

				/* Obtaining list of samples to choose from */
				sselect_UI(indvs, singsamps9, Ntot, 2, 0, 1, 3, deme);
				
				/* Sample that migrates */
				gsl_ran_choose(r,&bhc,1,singsamps9,(*(Nbet + deme) + 1),sizeof(unsigned int));

				for(j = 0; j < Ntot; j++){
					if(*((*(indvs + j)) + 0) == bhc){
						*((*(indvs + j)) + 3) = drec;			/* Altering deme accordingly */
						break;
					}
				}
				
				free(singsamps9);
				
			}
			break;
		case 10:	/* Event 10: Recombination - splitting a sample to create a new one */
			
			/* Converting WH to BH samples  */
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);
			
			/* For storing BH weights */
			double *weights10 = calloc((*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(double));
			double *singsamps10 = calloc((*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(double));
			double *startp = calloc((*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(double));			
			
			/* Obtaining weights of each sample */
			count = 0;
			while(count < (*(Nbet + deme) + 2*(*(nsex + deme)))){
				for(j = 0; j < (NBtot + 2*nsum); j++){
/*					printf("%d %d %d\n",*((*(nlri + j)) + 0),*((*(nlri + j)) + 1),*((*(nlri + j)) + 2),*((*(nlri + j)) + 3));	*/
					if( *((*(nlri + j)) + 3) == deme){
						*(weights10 + count) = *((*(nlri + j)) + 2);
						*(startp + count) = *((*(nlri + j)) + 1);						
						*(singsamps10 + count) = *((*(nlri + j)) + 0);
						count++;
					}
				}
			}
			
			/* For storing BH samples */
			unsigned int *singsamps10a = calloc((*(Nbet + deme) + 2*(*(nsex + deme))),sizeof(unsigned int));
				
			/* Choosing single sample to be split in recombination */
			gsl_ran_multinomial(r,(*(Nbet + deme) + 2*(*(nsex + deme))),1,weights10,singsamps10a);
				
			/* Matching choice to sample */
			for(j = 0; j < (*(Nbet + deme) + 2*(*(nsex + deme))); j++){
				if(*(singsamps10a + j) == 1){
					rands = *(singsamps10 + j);
					length = *(weights10 + j);
					start = *(startp + j);
					break;
				}
			}
					
			free(singsamps10a);
			
			rsite = gsl_rng_uniform_int(r,length) + (start + 1);
			/*
			if(rsite == 1 || rsite == (nsites-1)){
				printf("rsite is %d\n",rsite);
			}
			*/
			if(*nbreaks == 1){
				yesrec = 1;
				isyetbp = 0;
/*				isbpend = 1;	*/
				maxtr = 1;
			}else if(*nbreaks > 1){
				/* First, check if (1) rsite is already in existing table of breakpoints;
				and (2) if so, breakpoint included at end of existing table */
				isyetbp = 0;
/*				isbpend = 0;	*/
				maxtr = 0;
				for(x = 0; x < *nbreaks; x++){
					if( *((*(breaks + 0)) + x) == rsite){
						isyetbp = 1;
					}
					if( *((*(breaks + 0)) + x) <= rsite){
						maxtr = (x+1);
					}
				}
			}
			
/*			printf("Samp, length, Start, RS, maxtr, iybp are %d %d %d %d %d %d\n",rands,length,start,rsite,maxtr,isyetbp);	*/
				
			/*	
			yesrec = 0;
			while(yesrec != 1){
				
				 Initial choice for rec site */
				/* rsite = (unsigned int)(gsl_ran_flat(r, 0.51, (nsites-0.51)));	
				gsl_rng_uniform_int				if(rsite == 1 || rsite == (nsites-1)){
					printf("Rsite is %d\n",rsite);
				}
				*/
				
					/*
					if(maxtr == *nbreaks){
						isbpend = 1;
					}
					*/
					
					/* Next, checking if valid breakpoint depending on location */
					/*
					if( (isbpend == 1) && (isyetbp == 0) ){
						for(j = 0; j < NMax; j++){
							if( *((*(GType + j)) + 0) == rands){
								if( (*((*(GType + j)) + *nbreaks) != (-1)) ){
									yesrec = 1;
								}
								break;
							}
						}
					}
					else if( (isbpend == 1) && (isyetbp == 1) ){
						for(j = 0; j < NMax; j++){
							if( *((*(GType + j)) + 0) == rands){
								if((*((*(GType + j)) + *nbreaks) != (-1)) && (isallI((*(GType + j)), *nbreaks, -1, 1) != 1)){
									yesrec = 1;
								}
								break;
							}
						}
					}
					else 
					if( isyetbp == 0 ){			
						for(j = 0; j < NMax; j++){
							if( *((*(GType + j)) + 0) == rands){
								if( (isallI((*(GType + j)), (maxtr+1), (-1), 1) != 1) && (isallI((*(GType + j)), (*nbreaks+1), (-1), maxtr) != 1)){
									yesrec = 1;
								}
								break;
							}
						}
					}
					else if( isyetbp == 1 ){		
						for(j = 0; j < NMax; j++){
							if( *((*(GType + j)) + 0) == rands){
								if( (isallI((*(GType + j)), maxtr, (-1), 1) != 1) && (isallI((*(GType + j)), (*nbreaks+1), (-1), maxtr) != 1) ){
									yesrec = 1;
								}
								break;
							}
						}
					}					
				}	
			}		 End of BP verification step 
			*/
			/*					printf("Rands, rsite, iybp are %d %d %d\n",rands,rsite,isyetbp);	*/
			
			/* printf("rsite is %d\n",rsite);	*/
					
			/* Adding new site and re-ordering tracts in other tables */
			if((isyetbp != 1) && (*((*(GType + rands)) + maxtr) != (-1)) ){
				(*nbreaks)++;
				for(x = *nbreaks-2; x >= (int)(maxtr-1); x--){
					*((*(breaks + 0)) + x + 1) = *((*(breaks + 0)) + x);
					*((*(breaks + 1)) + x + 1) = *((*(breaks + 1)) + x);						
				}
				*((*(breaks + 0)) + maxtr) = rsite;
				*((*(breaks + 1)) + maxtr) = *((*(breaks + 1)) + maxtr - 1);
				/* Adding new site to genotype; coalescent time; ancestry table */
				for(j = 0; j < NMax; j++){
					for(x = (*nbreaks-1); x >= (int)(maxtr); x--){
						*((*(GType + j)) + x + 1) = *((*(GType + j)) + x);
						if(j < Itot){
							*((*(CTms + j)) + x + 1) = *((*(CTms + j)) + x);
							*((*(TAnc + j)) + x + 1) = *((*(TAnc + j)) + x);
						}
					}
				}
			}else if((isyetbp == 1) || (*((*(GType + rands)) + maxtr) == (-1) )){
				maxtr--;
			}
			
			/* First, checking parent number of existing sample, setting at paired */
			for(j = 0; j < Ntot; j++){
				if( *((*(indvs + j)) + 0) == rands){
					rpar = *((*(indvs + j)) + 1);
					*((*(indvs + j)) + 2) = 0;
					break;
				}
			}
			
			/* Adding new sample to indv table */
			*((*(indvs + NMax)) + 0) = NMax;
			*((*(indvs + NMax)) + 1) = rpar;
			*((*(indvs + NMax)) + 2) = 0;
			*((*(indvs + NMax)) + 3) = deme;
			
			/* Now creating the new sample genotype */
			for(j = 0; j < NMax; j++){
				if( *((*(GType + j)) + 0) == rands ){
					for(x = (*nbreaks); x > (int)maxtr; x--){
						*((*(GType + NMax)) + x) = *((*(GType + j)) + x);
						*((*(GType + j)) + x) = (-1);
					}
					break;
				}
			}
			
			for(x = maxtr; x > (int)0; x--){
				*((*(GType + NMax)) + x) = (-1);
			}
			*((*(GType + NMax)) + 0) = NMax;
			
			free(startp);
			free(singsamps10);
			free(weights10);
			break;
		case 11:
			/* Converting WH to BH samples */
			sexconv(indvs, rsex, nsum, Ntot, Nindv, ex);

			/* For storing WH indvs */
			unsigned int *parsamps11 = calloc((*(Nwith + deme) - *(nsex + deme)),sizeof(unsigned int));
			sselect_UI(indvs, parsamps11, Ntot, 2, 1, 0, 3, deme);
						
			/* A paired sample involved in coalescence */
			gsl_ran_choose(r,&WHsel,1,parsamps11,(*(Nwith + deme) - *(nsex + deme)),sizeof(unsigned int));
			nos = gsl_ran_bernoulli(r,0.5);		/* Which side involved in event */
			/* Finding sample */
			for(j = 0; j < Ntot; j++){
				if(*((*(indvs + j)) + 1) == WHsel){
					csamp = *((*(indvs + j + nos)) + 0);
					par = *((*(indvs + j + (nos+1)%2)) + 0);
					*((*(indvs + j + (nos+1)%2)) + 2) = 1;	/* Setting 'par' as BH */
					break;
				}
			}
			
			/* Now updating coalescent times */
			cchange(indvs, GType, CTms, TAnc, breaks, &csamp, &par, 1, Ntot, 0, nbreaks, Ttot, 1);
			
			/* Check if tracts have coalesced */
			ccheck(indvs,GType,breaks,nsites,Ntot,*nbreaks);
			
			if(achange == 1){
				excoal(indvs, GType, &par, *nbreaks, 1, Ntot, WCHex, BCHex, deme);
			}	
		
			free(parsamps11);
		
			break;	
		default:	/* If none of these cases chosen, exit with error message */
			fprintf(stderr,"Error: Non-standard coalescent case selected ('coalesce').\n");
			exit(1);
            break;
	}
	
	return(achange);
	
}	/* End of coalescent routine */

/* Updating information following coalescent event */
void cchange(unsigned int **indvs, int **GType, double **CTms, int **TAnc, unsigned int **breaks, unsigned int *csamp, unsigned int *par, unsigned int lsamp, unsigned int Ntot, unsigned int cst, unsigned int *cend, double Ttot, unsigned int isall){
	unsigned int j, i, x;
	unsigned int crow = 0;
	/*
	unsigned int prow = 0;
	*/
	int tempval = 0;
	
	/* Code works on assumption that columns 'align' 
	(i.e. column with 'csamp' in indvs is same column in other tables) */
	
	for(i = 0; i < lsamp; i++){
		/* Finding crow, prow */
		for(j = 0; j < Ntot; j++){
			if(*((*(indvs + j)) + 0) == *(csamp + i)){
				crow = j;
				break;
			}
		}
		/*
		for(j = 0; j < Ntot; j++){
			if(*((*(indvs + j)) + 0) == *(par + i)){
				prow = j;
				break;
			}
		}
		*/
		
		/* 'csamp' coalesces if ALL parts coalesce (so it's removed from the simulation) */
		if(isall == 1){
 			*((*(indvs + crow)) + 1) = HUGEVAL;
			*((*(indvs + crow)) + 2) = 2;
		}
		
		for(x = cst; x < (*cend); x++){
			
			/* Updating coalescent time and parental sample */
			tempval = -1;
			if( (*((*(GType + (*(csamp + i)))) + (x+1))) != (-1) && (*((*(GType + (*(par + i)))) + (x+1))) != (-1) ){
				tempval = *((*(GType + (*(csamp + i)))) + (x+1));
				*((*(CTms + (tempval))) + (x+1)) = Ttot;
				*((*(TAnc + (tempval))) + (x+1)) = (*((*(GType + (*(par + i)))) + (x+1)));
			}
			
			/* Updating genotype table */
			if( (*((*(GType + (*(csamp + i)))) + (x+1))) != (-1) && (*((*(GType + (*(par + i)))) + (x+1))) == (-1)){
				(*((*(GType + (*(par + i)))) + (x+1))) = (*((*(GType + (*(csamp + i)))) + (x+1)));
			}
			
			/* If coal event was via GC, then making coal'ed sample non-ancestral */
			if(isall == 0 && (*((*(GType + (*(csamp + i)))) + (x+1))) != (-1)){
				(*((*(GType + (*(csamp + i)))) + (x+1))) = (-1);
			}
			
		}
	}
	
}	/* End of 'cchange' function */

/* After coalescence, check if tracts have coalesced and make non-ancestral if need be */
unsigned int ccheck(unsigned int **indvs, int **GType, unsigned int **breaks, unsigned int nsites, unsigned int Ntot, unsigned int nbreaks){
	unsigned int j, x;			/* counters */
	unsigned int achange = 0;	/* Has there been 'a change'? */
	unsigned int gcount = 0;	/* count of extant tracts present */
	unsigned int ridx = 0;		/* Index of sample */
	int isfst = -1;				/* Is first indv that carries coal'ed material? */
	
	/* Checking if tract has completely coalesced */
	
	for(x = 0; x < nbreaks; x++){
		isfst = (-1);
		if( *((*(breaks + 1)) + x) != 1 ){		/* If not yet coalesced... */
			gcount = 0;
			for(j = 0; j < Ntot; j++){
				if( (*((*(indvs + j)) + 2) != 2 ) ){
					ridx = *((*(indvs + j)) + 0);
					if( (*((*(GType + ridx)) + (x + 1)) != (-1))){
						gcount++;
						if( isfst == (-1) ){
							isfst = ridx;
						}
					}
				}
			}
			
			/* If only one individual exists which carries that tract, then it has coalesced */
			/* Setting last ancestral piece to non-ancestral to reflect this */
			if(gcount == 1){
				*((*(breaks + 1)) + x) = 1;
				*((*(GType + isfst)) + (x + 1)) = (-1);
				achange = 1;
			}
		}
	}
	
	return(achange);
	
}	/* End of 'ccheck' function */

/* Check if extra samples are removed, as all extant material coalesces */
void excoal(unsigned int **indvs, int **GType, unsigned int *par, unsigned int nbreaks, unsigned int npar, unsigned int Ntot, int *WCHex, int *BCHex, unsigned int deme){
	
	unsigned int i, j;
	unsigned int prow = 0;
	unsigned int prent = 0;
	
	for(i = 0; i < npar; i++){
		if( isallI(*(GType + (*(par + i))), nbreaks + 1, (-1), 1) == 1 ) {
			for(j = 0; j < Ntot; j++){
				if( *((*(indvs + j)) + 0) == *(par + i) ){
					prow = j;
					break;
				}
			}
			
			if(*((*(indvs + prow)) + 2) == 1){
/*				printf("excoal 1\n");	*/
				*(BCHex + deme) -= 1;
			}else if(*((*(indvs + prow)) + 2) == 0){
/*				printf("excoal 2\n");			*/
				prent = *((*(indvs + prow)) + 1);
				for(j = 0; j < Ntot; j++){
					if( (*((*(indvs + j)) + 1) == prent) && (*((*(indvs + j)) + 0) != *(par + i) ) ){
						*((*(indvs + j)) + 2) = 1;
						break;
					}
				}		
				*(WCHex + deme) -= 1;
				*(BCHex + deme) += 1;
			}

			*((*(indvs + prow)) + 1) = HUGEVAL;
			*((*(indvs + prow)) + 2) = 2;
			
		}
	}
	
}

/* Routine to calculate length of coalesced tracts */
unsigned int coalcalc(unsigned int **breaks, unsigned int nsites, unsigned int nbreaks, unsigned int start){
	
	/* Calculating length of coalesced samples: 
	deducting 1 to account for edge effects; 
	adding on breakpoints lying between two coalesced samples */
	
	unsigned int x;
	unsigned int val1, val2;
	unsigned int lcoal = 0;
	
	for(x = start; x < (nbreaks - 1); x++){
		val1 = 0;
		val2 = 0;		
		if(*((*(breaks + 1)) + x) == 1){
			val1 = *((*(breaks + 0)) + x);
			val2 = *((*(breaks + 0)) + x + 1);
			lcoal += (val2 - val1 - 1);
			if(*((*(breaks + 1)) + x + 1) == 1){
				lcoal++;
			}
		}
	}
	
	val1 = 0;
	val2 = 0;		
	if(*((*(breaks + 1)) + (nbreaks - 1)) == 1){
		val1 = *((*(breaks + 0)) + (nbreaks - 1));
		val2 = nsites;
		lcoal += (val2 - val1 - 1);
	}

	return(lcoal);

}	/* End of coalcalc function */

/* Function to choose which samples split by sex */
void sexsamp(unsigned int **indvs, unsigned int *rsex, unsigned int *nsex, unsigned int *Nwith, unsigned int Ntot, const gsl_rng *r){
	unsigned int count = 0;
	unsigned int x, a;
	
	for(x = 0; x < d; x++){
		if(*(Nwith + x) != 0 && *(nsex + x) != 0){		/* Avoiding errors if no WH in deme */
			unsigned int *WHsamp = calloc(*(Nwith + x),sizeof(unsigned int));
			unsigned int *WHchosen = calloc(*(nsex + x),sizeof(unsigned int));
			/* Obtaining WH in that deme */
			sselect_UI(indvs, WHsamp, Ntot, 2, 1, 0, 3, x);
			/* Sampling those to be split */
			gsl_ran_choose(r,WHchosen,(*(nsex + x)),WHsamp,(*(Nwith + x)),sizeof(unsigned int));
			gsl_ran_shuffle(r, WHchosen, (*(nsex + x)), sizeof(unsigned int));
			
			/* Assigning samples to 'rsex' */
			for(a = 0; a < *(nsex + x); a++){
				*(rsex + (count + a)) = *(WHchosen + a);
			}
			count += *(nsex + x);
			
			free(WHchosen);
			free(WHsamp);
		}
	}
}

/* Insertion-sort rows by individual no., to ensure paired samples are together after an action */
void indv_sort(unsigned int **indvs, unsigned int nrow){
	unsigned int j, i;		/* Sorting indices */
	unsigned int temp0, temp1, temp2, temp3;	/* For swapping */
	unsigned int count = 0;
	unsigned int icount = 0;
	
	for(j = 1; j < nrow; j++){
		i = j;
		while( (i > 0) &&  ( *((*(indvs + (i - 1) )) + 1) > *((*(indvs + i)) + 1) )){
			/* Swapping entries */
			temp0 = *((*(indvs + (i - 1) )) + 0);
			temp1 = *((*(indvs + (i - 1) )) + 1);
			temp2 = *((*(indvs + (i - 1) )) + 2);
			temp3 = *((*(indvs + (i - 1) )) + 3);
				
			*((*(indvs + (i-1) )) + 0) =  *((*(indvs + (i) )) + 0);
			*((*(indvs + (i-1) )) + 1) =  *((*(indvs + (i) )) + 1);
			*((*(indvs + (i-1) )) + 2) =  *((*(indvs + (i) )) + 2);
			*((*(indvs + (i-1) )) + 3) =  *((*(indvs + (i) )) + 3);
				
			*((*(indvs + (i) )) + 0) = temp0;
			*((*(indvs + (i) )) + 1) = temp1;
			*((*(indvs + (i) )) + 2) = temp2;
			*((*(indvs + (i) )) + 3) = temp3;
			
			i--;
		}
	}
	
	/* Then renumbering indvs */
	while(count < nrow){
		if(*((*(indvs + count)) + 2) == 1){
			*((*(indvs + count)) + 1) = icount;
			icount++;
			count++;
		}else if(*((*(indvs + count)) + 2) == 0){
			*((*(indvs + count)) + 1) = icount;
			*((*(indvs + count + 1)) + 1) = icount;
			count++;
			count++;
			icount++;
		}else if(*((*(indvs + count)) + 2) == 2){
			count++;
		}
	}
	
}

/* Insertion-sort for double-type tables */
void indv_sortD(double **Tin, unsigned int nrow, unsigned int ncol, unsigned int tcol){
	unsigned int j, i, k;		/* Sorting indices */
	
	double *tempcol = calloc(ncol,sizeof(double));		/* temp entries for swapping */
	
	for(j = 1; j < nrow; j++){
		i = j;
		while( (i > 0) &&  ( *((*(Tin + (i - 1) )) + tcol) > *((*(Tin + i)) + tcol) )){
			/* Swapping entries */
			for(k = 0; k < ncol; k++){
				*(tempcol + k) = *((*(Tin + (i - 1) )) + k);
				*((*(Tin + (i - 1) )) + k) =  *((*(Tin + (i) )) + k);
				*((*(Tin + (i) )) + k) = *(tempcol + k);
			}			
			i--;
		}
	}
	free(tempcol);
}

void Wait(){
	printf("Press Enter to Continue");
	while( getchar() != '\n' );
	printf("\n");	
}

void TestTabs(unsigned int **indvs, int **GType, double **CTms, int **TAnc, unsigned int **breaks, unsigned int **nlri, unsigned int NMax, unsigned int Itot, unsigned int Nbet, unsigned int nbreaks){

	unsigned int j, x;
	
	printf("INDV TABLE\n");
	for(j = 0; j < NMax; j++){
		for(x = 0; x < 4; x++){
			printf("%d ",*((*(indvs + j)) + x));
		}
		printf("\n");
	}
	printf("\n");

	printf("GTYPE TABLE\n");
	for(j = 0; j < NMax; j++){
		for(x = 0; x <= nbreaks; x++){
			printf("%d ",*((*(GType + j)) + x));
		}
		printf("\n");
	}
	printf("\n");
/*
	printf("CTMS TABLE\n");
	for(j = 0; j < Itot; j++){
		for(x = 0; x <= nbreaks; x++){
			printf("%lf ",*((*(CTms + j)) + x));
		}
		printf("\n");
	}
	printf("\n");
	
	printf("TANC TABLE\n");
	for(j = 0; j < Itot; j++){
		for(x = 0; x <= nbreaks; x++){
			printf("%d ",*((*(TAnc + j)) + x));
		}
		printf("\n");
	}
	printf("\n");

	printf("NLRI TABLE\n");
	for(j = 0; j < Nbet; j++){
		for(x = 0; x < 4; x++){
			printf("%d ",*((*(nlri + j)) + x));
		}
		printf("\n");
	}
	printf("\n");
*/
	printf("BREAKS TABLE\n");
	for(j = 0; j < 2; j++){
		for(x = 0; x < nbreaks; x++){
			printf("%d ",*((*(breaks + j)) + x));
		}
		printf("\n");
	}
	printf("\n");	

/*	Wait();		*/
/*	exit(1);	*/

}

/* Function to reconstruct genealogy and to add mutation to branches */
char * treemaker(double **TFin, double thetain, unsigned int mind2, unsigned int maxd2, double mind, double maxd, unsigned int Itot, unsigned int run, double gmi, double gme, unsigned int ismsp, unsigned int *nmutT, unsigned int prtrees, unsigned int ismut, double pburst, unsigned int mburst, double bdist, const gsl_rng *r){
	unsigned int i, j, k, a, x;
	unsigned int lct = Itot;
	unsigned int lct2 = lct-1;
	unsigned int nc = 0;			/* {N}umber of {c}lades in current reconstruction */
	double birthtime = 0;			/* Coalescent time */
	double cbst = 0;				/* Centre of burst event */
	unsigned int child1 = 0;		/* Coalesced sample */
	unsigned int parent1 = 0;		/* Parental sample */	
	unsigned int csum = 0;			/* How many of each have been sampled, to decide action*/
	unsigned int ischild = 0;		/* Is it a child sample? */
	unsigned int rmut1 = 0;			/* Mutations along first branch */
	unsigned int rmut2 = 0;			/* Mutations along second branch */	
	unsigned int cc = 0;			/* Child clade */
	unsigned int pc = 0;			/* Parental clade */
	unsigned int exsamps = 0;		/* Extra samples */
	unsigned int cs = 0;			/* 'Clade' samps */
	unsigned int minc = 0;			/* Min, max clade (for resorting) */
	unsigned int maxc = 0;
	unsigned int ccM = 0;
	unsigned int brk = 0;			/* Breakpoint where tract starts */
	unsigned int clen = 4096;		/* Space allocated to clades array */
	unsigned int nmut = 0;			/* Number of mutants added */
	unsigned int n;				 	/* sprintf counter */
	unsigned int newr = 0;			/* New rows in table if need to expand */
	unsigned int isbst = 0;			/* Are next mutants clustered ? */
	unsigned int bsze = 0;			/* Burst counter */
	unsigned int isrb = 0;			/* Checking if burst distance lies in actual distance */
	
	static const char lbr[] = "(";
	static const char rbr[] = ")";
	static const char com[] = ",";
	static const char cln[] = ":";
	static const char scln[] = ";";
	static const char lsq[] = "[";
	static const char rsq[] = "]";
	static const char spa[] = " ";	
	char p1char[10];
	char c1char[10];
	char btchar1[16];
	char btchar2[16];	
	char brkchar[16];
	char Mout[32];				 /* String to hold filename in (Mutations) */
	FILE *ofp_mut = NULL;		 /* Pointer for data output */
		
	/* Defining necessary tables */
	double *Cheight = calloc(lct,sizeof(double));					/* Current 'height' (time) of each clade */
	char **clades = calloc(lct, sizeof(char *));					/*	Vector of possible clades */
	unsigned int **samps = calloc(lct,sizeof(unsigned int *));		/* Table of samples present in each clade (row = each clade) */
	for(j = 0; j < lct; j++){
		clades[j] = calloc((clen+1),sizeof(char));
		samps[j] = calloc(lct,sizeof(unsigned int));
		for(k = 0; k < lct; k++){
			*((*(samps + j)) + k) = Itot;
		}
	}
	
	/* Allocating space for mutation table */
	unsigned int MTRows = 30;
	double **MTab = calloc(MTRows,sizeof(double *));			/* Mutation table */
	for(j = 0; j < MTRows; j++){
		MTab[j] = calloc((Itot+1),sizeof(double));
	}


	for(i = 0; i < lct2; i++){
	
		birthtime = *((*(TFin + i)) + 1);
	    child1 = *((*(TFin + i)) + 0);
    	parent1 = *((*(TFin + i)) + 2);
    	ischild = 0;
    	csum = 0;
    	
    	if(i == 0){
    	
	    	*((*(samps + 0)) + 0) = parent1;
    		*((*(samps + 0)) + 1) = child1;
    		*(Cheight + 0) = birthtime;
    		
    		/* Converting values to strings */
    		sprintf(p1char, "%d", parent1);
	    	sprintf(c1char, "%d", child1);
	    	sprintf(btchar1,"%0.5lf",birthtime);    	

	    	strcpy(*(clades + 0),lbr);
	    	strcat(*(clades + 0),p1char);
	    	strcat(*(clades + 0),cln);
	    	strcat(*(clades + 0),btchar1);
	    	strcat(*(clades + 0),com);
	    	strcat(*(clades + 0),c1char);
	    	strcat(*(clades + 0),cln);
	    	strcat(*(clades + 0),btchar1);
	    	strcat(*(clades + 0),rbr);
    		
			/* Assigning mutations */
			rmut1 = gsl_ran_poisson(r,(0.5*thetain*birthtime));
			rmut2 = gsl_ran_poisson(r,(0.5*thetain*birthtime));
			if(rmut1 != 0){
				isbst = 0;
				bsze = 0;
				if((nmut + rmut1) > MTRows){
					newr = MTRows;
					while(newr < (nmut + rmut1)){
						newr = 2*newr;
					}
					MTab = (double **)realloc(MTab, newr*sizeof(double *));
					for(j = 0; j < MTRows; j++){
						MTab[j] = (double *)realloc( *(MTab + j) ,(Itot+1)*sizeof(double));
					}
					for(j = MTRows; j < newr; j++){
						MTab[j] = (double *)calloc((Itot + 1),sizeof(double));
					}
					MTRows = newr;
				}
				for(a = 0; a < rmut1; a++){
					/* Indicating location of mutants */
					if(bsze == 0){
						*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, mind, maxd);
						isbst = gsl_ran_bernoulli(r,pburst);
						if(isbst == 1){
							bsze = mburst;
							cbst = *((*(MTab + (nmut+a))) + 0);
						}
					}else if(bsze > 0){
						isrb = 0;
						while(isrb == 0){
							*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, cbst-bdist, cbst+bdist);
							if( (mind < *((*(MTab + (nmut+a))) + 0)) && (*((*(MTab + (nmut+a))) + 0) < maxd) ){
								isrb = 1;
							}
						}
						bsze--;
					}
					*((*(MTab + (nmut+a))) + (parent1+1)) = 1;
				}
				nmut += rmut1;
			}
			
			if(rmut2 != 0){
				isbst = 0;
				bsze = 0;			
				if((nmut + rmut2) > MTRows){
					newr = MTRows;
					while(newr < (nmut + rmut2)){
						newr = 2*newr;
					}
					MTab = (double **)realloc(MTab, newr*sizeof(double *));
					for(j = 0; j < MTRows; j++){
						MTab[j] = (double *)realloc( *(MTab + j) ,(Itot+1)*sizeof(double));
					}
					for(j = MTRows; j < newr; j++){
						MTab[j] = (double *)calloc((Itot + 1),sizeof(double));
					}
					MTRows = newr;
				}
				for(a = 0; a < rmut2; a++){
					/* Indicating location of mutants */
					if(bsze == 0){
						*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, mind, maxd);
						isbst = gsl_ran_bernoulli(r,pburst);
						if(isbst == 1){
							bsze = mburst;
							cbst = *((*(MTab + (nmut+a))) + 0);
						}
					}else if(bsze > 0){
						isrb = 0;
						while(isrb == 0){
							*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, cbst-bdist, cbst+bdist);
							if( (mind < *((*(MTab + (nmut+a))) + 0)) && (*((*(MTab + (nmut+a))) + 0) < maxd) ){
								isrb = 1;
							}
						}
						bsze--;
					}
					*((*(MTab + (nmut+a))) + (child1+1)) = 1;
				}
				nmut += rmut2;
			}
			
    	}else if(i > 0){
   			
    		/* There can be three cases: Merge clades if child already listed; 
    		Add to clade if child new but parent already listed; 
    		Create new clade otherwise.
    		
    	  	Testing how many of the pair have already been sampled, to decide tree reconstruction */
    	  	cc = Itot;
    	  	pc = Itot;
    	  	csum = 0;
	    	for(j = 0; j < lct; j++){
	    		if( *((*(samps + j)) + 0) == Itot ){
					break;
    			}
    			for(k = 0; k < lct; k++){
    				if( *((*(samps + j)) + k) == child1 ){
    					cc = j;
    					csum++;
    				}
    				if( *((*(samps + j)) + k) == parent1 ){
    					pc = j;
    					csum++;
    				}
    				if( *((*(samps + j)) + k) == Itot ){
						break;
    				}
    			}
	    	}
	    	  	
	    	if(csum==0){	/* Create a new clade */
    			nc++;
    			
	   			*((*(samps + nc)) + 0) = parent1;
    			*((*(samps + nc)) + 1) = child1;
    			*(Cheight + nc) = birthtime;

    			/* Converting values to strings */
	    		sprintf(p1char, "%d", parent1);
		    	sprintf(c1char, "%d", child1);
	    		sprintf(btchar1,"%0.5lf",birthtime);

	    		strcpy(*(clades + nc),lbr);
		    	strcat(*(clades + nc),p1char);
		    	strcat(*(clades + nc),cln);
	    		strcat(*(clades + nc),btchar1);
	    		strcat(*(clades + nc),com);
		    	strcat(*(clades + nc),c1char);
		    	strcat(*(clades + nc),cln);
	    		strcat(*(clades + nc),btchar1);
		    	strcat(*(clades + nc),rbr);
	    
	    		/* Assigning mutations */
				rmut1 = gsl_ran_poisson(r,(0.5*thetain*birthtime));
				rmut2 = gsl_ran_poisson(r,(0.5*thetain*birthtime));
				if(rmut1 != 0){
					isbst = 0;
					bsze = 0;
					if((nmut + rmut1) > MTRows){
						newr = MTRows;
						while(newr < (nmut + rmut1)){
							newr = 2*newr;
						}
						MTab = (double **)realloc(MTab, newr*sizeof(double *));
						for(j = 0; j < MTRows; j++){
							MTab[j] = (double *)realloc( *(MTab + j) ,(Itot+1)*sizeof(double));
						}
						for(j = MTRows; j < newr; j++){
							MTab[j] = (double *)calloc((Itot + 1),sizeof(double));
						}
						MTRows = newr;
					}
					for(a = 0; a < rmut1; a++){
						/* Indicating location of mutants */
						if(bsze == 0){
							*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, mind, maxd);
							isbst = gsl_ran_bernoulli(r,pburst);
							if(isbst == 1){
								bsze = mburst;
								cbst = *((*(MTab + (nmut+a))) + 0);
							}
						}else if(bsze > 0){
							isrb = 0;
							while(isrb == 0){
								*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, cbst-bdist, cbst+bdist);
								if( (mind < *((*(MTab + (nmut+a))) + 0)) && (*((*(MTab + (nmut+a))) + 0) < maxd) ){
									isrb = 1;
								}
							}
							bsze--;
						}
						*((*(MTab + (nmut+a))) + (parent1+1)) = 1;
					}
					nmut += rmut1;
				}
		
				if(rmut2 != 0){
					isbst = 0;
					bsze = 0;				
					if((nmut + rmut2) > MTRows){
						newr = MTRows;
						while(newr < (nmut + rmut2)){
							newr = 2*newr;
						}
						MTab = (double **)realloc(MTab, newr*sizeof(double *));
						for(j = 0; j < MTRows; j++){
							MTab[j] = (double *)realloc( *(MTab + j) ,(Itot+1)*sizeof(double));
						}
						for(j = MTRows; j < newr; j++){
							MTab[j] = (double *)calloc((Itot + 1),sizeof(double));
						}
						MTRows = newr;
					}
					for(a = 0; a < rmut2; a++){
						/* Indicating location of mutants */
						if(bsze == 0){
							*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, mind, maxd);
							isbst = gsl_ran_bernoulli(r,pburst);
							if(isbst == 1){
								bsze = mburst;
								cbst = *((*(MTab + (nmut+a))) + 0);
							}
						}else if(bsze > 0){
							isrb = 0;
							while(isrb == 0){
								*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, cbst-bdist, cbst+bdist);
								if( (mind < *((*(MTab + (nmut+a))) + 0)) && (*((*(MTab + (nmut+a))) + 0) < maxd) ){
									isrb = 1;
								}
							}
							bsze--;
						}
						*((*(MTab + (nmut+a))) + (child1+1)) = 1;
					}
					nmut += rmut2;
				}		
				
    		}else if(csum == 1){	/* Add to existing clade */
				/* Choosing row (and therefore clade) containing existing parent (or child) */
				if(pc==Itot){
   					ischild = 1;
   					pc = cc;
   				}
   				
   				/* Converting values to strings */
				sprintf(c1char, "%d", child1);
				sprintf(p1char, "%d", parent1);				
				sprintf(btchar1,"%0.5lf",birthtime);
				sprintf(btchar2,"%0.5lf",(birthtime - (*(Cheight+pc))));
				char *tc = malloc( (strlen((*(clades + pc))) + 40) * sizeof(char) );
				if( tc == NULL ) {
    				fprintf(stderr, "Error - unable to allocate required memory\n");
				}
    			
    			if(ischild == 0){
    			
					/* Concatenating new clade */
					strcpy(tc,lbr);
					strcat(tc,c1char);
					strcat(tc,cln);
					strcat(tc,btchar1);
					strcat(tc,com);
					strcat(tc,(*(clades + pc)));
					strcat(tc,cln);
					strcat(tc,btchar2);
					strcat(tc,rbr);
					
					for(k = 0; k < lct; k++){
    					if( *((*(samps + pc)) + k) == Itot ){
    						*((*(samps + pc)) + k) = child1;
    						break;
	    				}
    				}
    				exsamps = child1;
    			}else if(ischild==1){
    			
					strcpy(tc,lbr);
					strcat(tc,p1char);
					strcat(tc,cln);
					strcat(tc,btchar1);
					strcat(tc,com);
					strcat(tc,(*(clades + pc)));
					strcat(tc,cln);
					strcat(tc,btchar2);
					strcat(tc,rbr);
					
					for(k = 0; k < lct; k++){
    					if( *((*(samps + pc)) + k) == Itot ){
    						*((*(samps + pc)) + k) = parent1;
    						break;
	    				}
    				}
    				exsamps = parent1;
    			}
   				
			   	/* Assigning mutations */
				rmut1 = gsl_ran_poisson(r,(0.5*thetain*(birthtime - (*(Cheight+pc)))));
				rmut2 = gsl_ran_poisson(r,(0.5*thetain*birthtime));
				if(rmut1 != 0){
					isbst = 0;
					bsze = 0;
					if((nmut + rmut1) > MTRows){
						newr = MTRows;
						while(newr < (nmut + rmut1)){
							newr = 2*newr;
						}
						MTab = (double **)realloc(MTab, newr*sizeof(double *));
						for(j = 0; j < MTRows; j++){
							MTab[j] = (double *)realloc( *(MTab + j) ,(Itot+1)*sizeof(double));
						}
						for(j = MTRows; j < newr; j++){
							MTab[j] = (double *)calloc((Itot + 1),sizeof(double));
						}
						MTRows = newr;
					}
					for(a = 0; a < rmut1; a++){
						/* Indicating location of mutants */
						if(bsze == 0){
							*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, mind, maxd);
							isbst = gsl_ran_bernoulli(r,pburst);
							if(isbst == 1){
								bsze = mburst;
								cbst = *((*(MTab + (nmut+a))) + 0);
							}
						}else if(bsze > 0){
							isrb = 0;
							while(isrb == 0){
								*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, cbst-bdist, cbst+bdist);
								if( (mind < *((*(MTab + (nmut+a))) + 0)) && (*((*(MTab + (nmut+a))) + 0) < maxd) ){
									isrb = 1;
								}
							}
							bsze--;
						}
						
						for(k = 0; k < lct; k++){
							if( (*((*(samps + pc)) + k) != Itot) && (*((*(samps + pc)) + k) != exsamps) ){
								cs = *((*(samps + pc)) + k);
								*((*(MTab + (nmut+a))) + (cs + 1)) = 1;
							}
						}
					}
					nmut += rmut1;
				}
		
				if(rmut2 != 0){
					isbst = 0;
					bsze = 0;
					if((nmut + rmut2) > MTRows){
						newr = MTRows;
						while(newr < (nmut + rmut2)){
							newr = 2*newr;
						}
						MTab = (double **)realloc(MTab, newr*sizeof(double *));
						for(j = 0; j < MTRows; j++){
							MTab[j] = (double *)realloc( *(MTab + j) ,(Itot+1)*sizeof(double));
						}
						for(j = MTRows; j < newr; j++){
							MTab[j] = (double *)calloc((Itot + 1),sizeof(double));
						}
						MTRows = newr;
					}
					for(a = 0; a < rmut2; a++){
						/* Indicating location of mutants */
						if(bsze == 0){
							*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, mind, maxd);
							isbst = gsl_ran_bernoulli(r,pburst);
							if(isbst == 1){
								bsze = mburst;
								cbst = *((*(MTab + (nmut+a))) + 0);
							}
						}else if(bsze > 0){
							isrb = 0;
							while(isrb == 0){
								*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, cbst-bdist, cbst+bdist);
								if( (mind < *((*(MTab + (nmut+a))) + 0)) && (*((*(MTab + (nmut+a))) + 0) < maxd) ){
									isrb = 1;
								}
							}
							bsze--;
						}
						
						*((*(MTab + (nmut+a))) + (exsamps + 1)) = 1;
					}
					nmut += rmut2;
				}
				
				memset((*(clades + pc)),'\0',(clen+1));
				if(strlen(tc) > clen){
					while(strlen(tc) > clen){
						clen = 2*clen;
					}						
					for(j = 0; j < lct; j++){
						clades[j] = realloc((*(clades + j)), (clen+1) * sizeof(char) );
					}
				}				
				strcpy((*(clades + pc)),tc);
   				*(Cheight + pc) = birthtime;
				free(tc);   				
   				
    		}else if((csum == 2) && (nc > 0)){		/* Add to existing clade */
    		
    			/* Converting values to strings */
				sprintf(c1char, "%d", child1);
				sprintf(p1char, "%d", parent1);				
				sprintf(btchar1,"%0.5lf",(birthtime - (*(Cheight + pc))));
				sprintf(btchar2,"%0.5lf",(birthtime - (*(Cheight + cc))));
				/* memset(tc,'\0',sizeof(tc)); */
				char *tc2 = malloc((strlen((*(clades + pc))) + strlen((*(clades + cc))) + 40) * sizeof(char));
				if( tc2 == NULL ) {
    				fprintf(stderr, "Error - unable to allocate required memory\n");
				}
				
				strcpy(tc2,lbr);
				strcat(tc2,(*(clades + pc)));
				strcat(tc2,cln);
				strcat(tc2,btchar1);
				strcat(tc2,com);
				strcat(tc2,(*(clades + cc)));
				strcat(tc2,cln);
				strcat(tc2,btchar2);
				strcat(tc2,rbr);
				
				/* Assigning mutations */
				rmut1 = gsl_ran_poisson(r,(0.5*thetain*(birthtime - (*(Cheight+pc)))));
				rmut2 = gsl_ran_poisson(r,(0.5*thetain*(birthtime - (*(Cheight+cc)))));
				if(rmut1 != 0){
					isbst = 0;
					bsze = 0;
					if((nmut + rmut1) > MTRows){
						newr = MTRows;
						while(newr < (nmut + rmut1)){
							newr = 2*newr;
						}
						MTab = (double **)realloc(MTab, newr*sizeof(double *));
						for(j = 0; j < MTRows; j++){
							MTab[j] = (double *)realloc( *(MTab + j) ,(Itot+1)*sizeof(double));
						}
						for(j = MTRows; j < newr; j++){
							MTab[j] = (double *)calloc((Itot + 1),sizeof(double));
						}
						MTRows = newr;
					}
					for(a = 0; a < rmut1; a++){
						/* Indicating location of mutants */
						if(bsze == 0){
							*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, mind, maxd);
							isbst = gsl_ran_bernoulli(r,pburst);
							if(isbst == 1){
								bsze = mburst;
								cbst = *((*(MTab + (nmut+a))) + 0);
							}
						}else if(bsze > 0){
							isrb = 0;
							while(isrb == 0){
								*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, cbst-bdist, cbst+bdist);
								if( (mind < *((*(MTab + (nmut+a))) + 0)) && (*((*(MTab + (nmut+a))) + 0) < maxd) ){
									isrb = 1;
								}
							}
							bsze--;
						}

						for(k = 0; k < lct; k++){
							if( *((*(samps + pc)) + k) != Itot ){
								cs = *((*(samps + pc)) + k);
								*((*(MTab + (nmut+a))) + (cs + 1)) = 1;
							}
						}
					}
					nmut += rmut1;
				}
				
				if(rmut2 != 0){
					isbst = 0;
					bsze = 0;				
					if((nmut + rmut2) > MTRows){
						newr = MTRows;
						while(newr < (nmut + rmut2)){
							newr = 2*newr;
						}
						MTab = (double **)realloc(MTab, newr*sizeof(double *));
						for(j = 0; j < MTRows; j++){
							MTab[j] = (double *)realloc( *(MTab + j) ,(Itot+1)*sizeof(double));
						}
						for(j = MTRows; j < newr; j++){
							MTab[j] = (double *)calloc((Itot + 1),sizeof(double));
						}
						MTRows = newr;
					}
					for(a = 0; a < rmut2; a++){
						/* Indicating location of mutants */
						if(bsze == 0){
							*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, mind, maxd);
							isbst = gsl_ran_bernoulli(r,pburst);
							if(isbst == 1){
								bsze = mburst;
								cbst = *((*(MTab + (nmut+a))) + 0);
							}
						}else if(bsze > 0){
							isrb = 0;
							while(isrb == 0){
								*((*(MTab + (nmut+a))) + 0) = gsl_ran_flat(r, cbst-bdist, cbst+bdist);
								if( (mind < *((*(MTab + (nmut+a))) + 0)) && (*((*(MTab + (nmut+a))) + 0) < maxd) ){
									isrb = 1;
								}
							}
							bsze--;
						}
						
						for(k = 0; k < lct; k++){
							if( *((*(samps + cc)) + k) != Itot ){
								cs = *((*(samps + cc)) + k);
								*((*(MTab + (nmut+a))) + (cs + 1)) = 1.0;
							}
						}
					}
					nmut += rmut2;
				}
				
				/* Deleting old clade data */
				if(cc < pc){
					minc = cc;
					maxc = pc;
				}else if (cc > pc){
					minc = pc;
					maxc = cc;
				}
				
				if(strlen(tc2) > clen){
					while(strlen(tc2) > clen){
						clen = 2*clen;
					}	
					for(j = 0; j < lct; j++){
						clades[j] = realloc((*(clades + j)), (clen+1) * sizeof(char) );
					}
				}
				memset((*(clades + minc)),'\0',(clen+1));
				memset((*(clades + maxc)),'\0',(clen+1));				
				strcpy((*(clades + minc)),tc2);
   				*(Cheight + minc) = birthtime;
   				*(Cheight + maxc) = 0;
   				free(tc2);
   				
   				for(k = 0; k < lct; k++){
   					if(*((*(samps + minc)) + k) == Itot){
   						ccM = k;
   						break;
   					}
   				}
   				for(k = 0; k < lct; k++){
   					if(*((*(samps + maxc)) + k) == Itot){
   						break;
   					}
   					*((*(samps + minc)) + (k + ccM)) = *((*(samps + maxc)) + k);
/*					*((*(samps + maxc)) + k) = Itot; */
   				}
   				
   				/* Now to reorder (to prevent overflow/going out of array length) */  			
    			/* Then re-writing over original arrays */  			
    			for(j = maxc; j < nc; j++){
					for(k = 0; k < lct; k++){					
						*((*(samps + j)) + k) = *((*(samps + j + 1)) + k);
   					}
					memset((*(clades + j)),'\0',(clen+1));
					strcpy((*(clades + j)),(*(clades + j + 1)));
					*(Cheight + j) = *(Cheight + j + 1);
    			}
    			
    			for(k = 0; k < lct; k++){
					*((*(samps + nc)) + k) = Itot;
				}
				memset((*(clades + nc)),'\0',(clen+1));
				*(Cheight + nc) = 0;
				
    			nc--;
    			
    		}
   			
    	}
	}
	
	/* Printing out Mutations to file */
	indv_sortD(MTab,nmut,(Itot+1),0);
	if(thetain > 0){
		n = sprintf(Mout,"Mutations/Muts_%d.dat",run);
		if(mind != 0){
			ofp_mut = fopen(Mout,"a");
		}else if(mind == 0){
			ofp_mut = fopen(Mout,"w");
		}
		for(j = 0; j < nmut; j++){
			fprintf(ofp_mut,"%lf ",*((*(MTab + j)) + 0));
			for(a = 1; a < (Itot + 1); a++){
				fprintf(ofp_mut,"%d ",(unsigned int)(*((*(MTab + j)) + a)));
			}
			fprintf(ofp_mut,"\n");
		}
		fclose(ofp_mut);
	}
	
	*nmutT += nmut;
	
	strcat((*(clades + 0)),scln);
    char *str = malloc(clen * sizeof(char));
    if(str == NULL) {
    	fprintf(stderr, "Error - unable to allocate required memory (str1) \n");
	}

    char *str2 = malloc(clen * sizeof(char));
    if(str2 == NULL) {
    	fprintf(stderr, "Error - unable to allocate required memory (str2) \n");
	}
    	
    if( (rec == 0 && gmi == 0 && gme == 0) || nsites == 1){
    	strcpy(str,(*(clades + 0)));
    	strcpy(str2,(*(clades + 0)));    	
    }else if((rec > 0 || gmi > 0 || gme > 0) && nsites > 1){
    	brk = (maxd2 - mind2);
    	sprintf(brkchar, "%d", brk);
    	strcpy(str,lsq);
    	strcpy(str2,lsq);    	
    	strcat(str,brkchar);
    	strcat(str2,brkchar);    	
    	strcat(str,rsq);
    	strcat(str2,rsq);    	
    	strcat(str,spa);
    	strcat(str,(*(clades + 0)));
	    strcat(str2,(*(clades + 0)));	
    }
    
    if(ismsp == 1){
    	if(prtrees == 1){
	    	printf("%s\n",str2);
    	}
    	/* If last tree, then print out mutation table */
    	if(maxd == 1 && ismut == 1){
    		/* If last tree and using MS-style printout, read in data and print to screen */
    		double **MTabF = calloc(*nmutT,sizeof(double *));			/* Mutation table */
			for(j = 0; j < *nmutT; j++){
				MTabF[j] = calloc((Itot+1),sizeof(double));
			}
			
			/*	Reading in table */
			n = sprintf(Mout,"Mutations/Muts_%d.dat",run);
			ofp_mut = fopen(Mout,"r");
			x = 1;
		    for(j = 0; j < *nmutT; j++){
		    	for(a = 0; a < (Itot + 1); a++){
					x = fscanf(ofp_mut,"%lf", &(MTabF[j][a]) );
				}
		    }
			fclose(ofp_mut);
			
			/* Print off sites and positions */
    		printf("segsites: %d\n",*nmutT);
    		printf("positions: ");
    		for(j = 0; j < *nmutT; j++){
    			printf("%.4lf ",*((*(MTabF + j)) + 0));
    		}
    		printf("\n");
    		
    		/* Now printing off mutation tables */
	    	for(a = 1; a < (Itot + 1); a++){
	    		for(j = 0; j < *nmutT; j++){
		    		printf("%d",(unsigned int)(*((*(MTabF + j)) + a)));
				}
			    printf("\n");
		    }

			/* Then free tables */
			for(j = 0; j < *nmutT; j++){
				free(MTabF[j]);
			}
			free(MTabF);
    	}
    }
   	free(str2);
    
    for(j = 0; j < MTRows; j++){
		free(MTab[j]);
	}
	for(j = 0; j < lct; j++){										
		free(samps[j]);
		free(clades[j]);
	}
	free(MTab);
	free(samps);
	free(clades);
	free(Cheight);
	return(str);
	free(str);

}	/* End of treemaker routine */

/* Reccal: calculating effective recombination rate over potential sites */
void reccal(unsigned int **indvs, int **GType, unsigned int **breaks, unsigned int **nlri, unsigned int *Nbet, unsigned int *Nwith, unsigned int *rsex, unsigned int esex, unsigned int *lnrec, unsigned int nbreaks, unsigned int NMax, unsigned int sw, unsigned int run){
	unsigned int j, i;
	unsigned int count = 0;
	unsigned int count2 = 0;	
	unsigned int vl = 0;
	unsigned int mintr = 0;
	unsigned int maxtr = 0;
	unsigned int Ntot = sumUI(Nbet,d) + 2*sumUI(Nwith,d);
	unsigned int is0l = 0;
	unsigned int is0r = 0;
	unsigned int ridx = 0;
	unsigned int brec = 0; 		/* Accounted breakpoints */
	unsigned int intot = 0;
	unsigned int count3 = 0;
	if(sw == 1){
		count3 = sumUI(Nbet,d);
	}
	
	if(sw == 0){
		vl = sumUI(Nbet,d);
	}else if(sw == 1){
		vl = 2*esex;
	}
	unsigned int *BHi = calloc(vl,sizeof(unsigned int));
	unsigned int *BHid = calloc(vl,sizeof(unsigned int));
	
	/* Vector of samples */
	count = 0;
	count2 = 0;
	if(sw == 0){
		while(count < vl){
			for(j = 0; j < Ntot; j++){
				if( *((*(indvs + j)) + 2) == 1){
					*(BHi + count) = *((*(indvs + j)) + 0);
					*(BHid + count) = *((*(indvs + j)) + 3);
					count++;
				}
			}
		}
	}else if(sw == 1){
		while(count < vl){
			for(j = 0; j < Ntot; j++){
				if( *((*(indvs + j)) + 1) == *(rsex + count2)){
					*(BHi + count) = *((*(indvs + j)) + 0);
					*(BHid + count) = *((*(indvs + j)) + 3);
					*(BHi + count + 1) = *((*(indvs + j + 1)) + 0);
					*(BHid + count + 1) = *((*(indvs + j + 1)) + 3);			
					count++;
					count++;
					count2++;					
					break;
				}
			}
		}
	}
	
	for(i = 0; i < d; i++){
		*(lnrec + i) = 0;
	}
	
	if(vl > 0){
		for(i = 0; i < vl; i++){
			mintr = 0;
			maxtr = 0;
			is0l = 0;
			is0r = 0;			
			ridx = 0;
			brec = 0;
			intot = 0;
			*((*(nlri + count3 + i)) + 0) = (*(BHi + i));
			*((*(nlri + count3 + i)) + 1) = 0;
			/* Determining case to run */
			for(j = 0; j < NMax; j++){
				if( *((*(GType + j)) + 0) == *(BHi + i) ){
					ridx = j;
					if( *((*(GType + j)) + 1) == -1 ){
						is0l = 1;
					}
					if( *((*(GType + j)) + nbreaks) == -1 ){
						is0r = 1;
					}
					break;
				}
			}
			if( is0l == 1 ){
				mintr = first_neI(*(GType + ridx), nbreaks + 1, (-1), 1);
				mintr--;	/* So concordant with 'breaks' table */
				*((*(nlri + count3 + i)) + 1) = *((*(breaks + 0)) + mintr);
				brec = *((*(breaks + 0)) + mintr);
				*(lnrec + (*(BHid + i))) += brec;
				intot = brec;
			}
			
			brec = 0;
			if( is0r == 1 ){
				maxtr = last_neI(*(GType + ridx), nbreaks+1, (-1), 1);
				maxtr--;	/* So concordant with 'breaks' table */
				brec = nsites - *((*(breaks + 0)) + maxtr + 1);
				*(lnrec + (*(BHid + i))) += brec;
				intot += brec;
			}
			*((*(nlri + count3 + i)) + 2) = (nsites - 1 - intot);
			*((*(nlri + count3 + i)) + 3) = (*(BHid + i));
		}
	}

	free(BHid);
	free(BHi);
}

void proberr(unsigned int est, double **pr){

	unsigned int j, x;
	
	fprintf(stderr,"A negative probability exists, you need to double-check your algebra (or probability inputs) - esex %d.\n",est);
	fprintf(stderr,"This likely arises due to having excessively high recombination or gene conversion rates - please check.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"For error reporting, the final probability calculations per deme are:\n");
	
	for(j = 0; j < 12; j++){
		fprintf(stderr,"Event %d: ",j);
		for(x = 0; x < d; x++){			
			fprintf(stderr,"%0.10lf ",(*((*(pr + j)) + x)));
		}
		fprintf(stderr,"\n");
	}
	
	fprintf(stderr,"\n");
	exit(1);				

}

void proberr2(double **pr){

	unsigned int j, x;
	
	fprintf(stderr,"Summed probabilities exceed one, you need to double-check your algebra (or probability inputs).\n");
	fprintf(stderr,"This likely arises due to having excessively high recombination or gene conversion rates - please check.\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"For error reporting, the final probability calculations per deme are:\n");
	
	for(j = 0; j < 12; j++){
		fprintf(stderr,"Event %d: ",j);
		for(x = 0; x < d; x++){			
			fprintf(stderr,"%0.10lf ",(*((*(pr + j)) + x)));
		}
		fprintf(stderr,"\n");
	}
	
	fprintf(stderr,"\n");
	exit(1);				

}

void proberr3(double **pr){

	unsigned int j, x;
	
	fprintf(stderr,"Summed probabilites are zero or negative, you need to double-check your algebra (or probability inputs).\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"For error reporting, the final probability calculations per deme are:\n");
	
	for(j = 0; j < 12; j++){
		fprintf(stderr,"Event %d: ",j);
		for(x = 0; x < d; x++){			
			fprintf(stderr,"%0.10lf ",(*((*(pr + j)) + x)));
		}
		fprintf(stderr,"\n");
	}
	
	fprintf(stderr,"\n");
	exit(1);				

}

void manyr(){
	fprintf(stderr,"Too many recombinants (exceeds HUGEVAL), exiting program.\n");
	fprintf(stderr,"This is likely due to having excessively large recombination\n");
	fprintf(stderr,"or gene conversion rates.\n");
	exit(1);
}

void usage(){
fprintf(stderr,"\n");
fprintf(stderr,"Command: FacSexCoalescent <Population Size> <Paired Samples> <Single Samples> <Rate of Sex> <Reps> <Other parameters>\n");
fprintf(stderr,"\n");
fprintf(stderr,"'Other parameters' include:\n");
fprintf(stderr," -t: [4Nmu] defines neutral mutation rate\n");
fprintf(stderr," -T: prints out individuals trees to file\n");
fprintf(stderr," (Note that one of -t or -T MUST be used)\n");
fprintf(stderr," -P: prints out data to screen using MS-style format\n");
fprintf(stderr," -r: [2Nr nsites] to specify recombination\n");
fprintf(stderr," -c: [2N(g_me) lambda_me] specifies rate and mean length of MEIOTIC gene conversion\n");
fprintf(stderr," -m: [2N(g_mi) lambda_mi] specifies rate and mean length of MITOTIC gene conversion\n");
fprintf(stderr," -b: [p_burst max_burst dist] specifies parameters for some polymorphisms to cluster ('burst')\n");
fprintf(stderr," (For -c, -m, -b: need to first use -r to specify number of sites, even if 2Nr = 0)\n");
/*fprintf(stderr," -G: [alpha] specifies population growth (scaled by 2N), N(t) = N0 exp(-alpha*t) \n");*/
fprintf(stderr," -I: d [Paired_j Single_j]...[2Nm] for defining island model.\n");
fprintf(stderr,"     d is number of demes: [Paired_j Single_j] are d pairs of samples per deme;\n");
fprintf(stderr,"     2Nm is net migration rate between demes.\n");
fprintf(stderr," -H: [Type of heterogeneity] [pLH pHL] [SexL_j SexH_j] for heterogeneity in rates of sex.\n");
fprintf(stderr,"     Type of heterogeneity is 0 for constant switching; 1 for stepwise change; 2 for no temporal change (only spatial change).\n");
fprintf(stderr,"     If type = 0; pLH pHL is probability of switching from low-sex to high-sex state (and vice versa).\n");
fprintf(stderr,"     If type = 1; pLH is time (units of 2N generations) when switch occurs.\n");
fprintf(stderr,"     If type = 2; these do not need to be defined.\n");
fprintf(stderr,"     [SexL_j SexH_j] are d pairs of low-sex, high-sex rates in deme j (if type = 0).\n");
fprintf(stderr,"     OR they are d pairs of initial-rate, changed-rate of sex (if type = 1).\n");
fprintf(stderr,"     Only SexL_j is defined if just spatial heterogeneity is present (if type = 2).\n");
fprintf(stderr,"     Note that with -H, one must first specify -I to specify subdivision (even if only one population).\n");
fprintf(stderr,"\n");
fprintf(stderr,"See documentation for further details.\n");
fprintf(stderr,"\n");
exit(1);
}



/* Main program */
int main(int argc, char *argv[]){
	unsigned int x, i, j;		/* Assignment counter, rep counter, indv counter */
	unsigned int pST = 0;		/* State of reproduction heterogeneity */	
	unsigned int npST = 0;		/* State of reproduction heterogeneity */	
	unsigned int Ntot = 0;		/* Total number of samples at time */
	unsigned int Nindv = 0;		/* Total number of individuals */
	unsigned int IwithT = 0;	/* Total within individual samples */
	unsigned int IbetT = 0;		/* Total between individual samples */
	unsigned int IwithT2 = 0;	/* Sum for checking structured population case */
	unsigned int IbetT2 = 0;	/* Sum for checking structured population case */
	unsigned int IwithC = 0;	/* Cumulative Iwith sum */
	unsigned int IbetC = 0;		/* Cumulative Ibet sum */
	unsigned int esex = 0;		/* Total sex events (segregation of paired samples) */
	unsigned int CsexS = 0;		/* Switch once sex samples chosen */
	unsigned int done = 0;		/* Is simulation complete? */
	unsigned int nbreaks = 0;	/* Number of non-rec tracts */
	unsigned int event = 0;		/* What event happens? */
	unsigned int deme = 0;		/* Which deme does event happen in? */
	unsigned int drec = 0;		/* Receiving deme for migration event */
	unsigned int e2 = 0;		/* Outcome of mig sampling, type of deme that migrates */
	unsigned int count = 0;		/* For creating ancestry table */	
	unsigned int exr = INITBR;	/* Extra rows */
	unsigned int exc = INITBR;	/* Extra columns */
	unsigned int NMax = 0;		/* Max samples present (for correct table searching!) */
	unsigned int Iindv = 0;		/* Number of initial individuals */
	unsigned int Itot = 0;		/* Number of initial samples */	
	unsigned int Nreps = 0;		/* Number of simulation replicates */	
	unsigned int pSTIN = 2;		/* Determine type of heterogeneity (0 = fluctuating sex, 1 = stepwise change, 2 = constant sex) */
	unsigned int N = 0;			/* Population Size */
	unsigned int gcalt = 0;		/* How to alter number samples following GC event */
	unsigned int argx = 6;		/* Extra command line inputs, for parameter definitions */
	unsigned int prtrees = 0;	/* Print out trees to file */
	unsigned int isrec = 0;		/* Has rec been defined? */
	unsigned int ismig = 0;		/* Has population structure been defined? */	
	unsigned int isburst = 0;	/* Has a burst of mutations been defined? */
	unsigned int ismsp = 0;		/* Use MS-style printout? */
	unsigned int nmutT = 0;		/* Total mutants (for printout) */
	unsigned int mind2 = 0;
	unsigned int maxd2 = 0;	
	unsigned int ismut = 0;		/* Mutation defined? */
	unsigned int mburst = 0;	/* Max burst size */
	unsigned int isexp = 0;		/* Assume exponential growth/decay? */
	unsigned int achange = 0;	/* Has there been a coal event? If so then extra checks needed */
	double bsex = 0;			/* Baseline rate of sex (for initial inputting) */
	double pLH = 0;				/* Prob of low-sex to high-sex transition, OR time of transition if stepwise change */
	double pHL = 0;				/* Prob of high-sex to low-sex transition */
	double Ttot = 0;			/* Time in past, initiate at zero */
	double NextT = 0;			/* Next time, after drawing event */	
	double tls = 0;				/* 'Time since Last Switch' or tls */
	double tts = 0;				/* 'Time to switch' */
	double nosex = 0;			/* Probability of no sexual reproduction over all demes */
	double psum = 0;			/* Sum of transition probs (first go) */
	double tjump = 0;			/* Time until next event */
	double mind = 0;			/* Min tract dist */
	double maxd = 0;			/* Max tract dist */
	double gmi = 0;				/* MITOTIC gene conversion probability */
	double gme = 0;				/* MEIOTIC gene conversion probability */	
	double theta = 0;			/* Scaled mutation rate, 4Nmu */	
	double mig = 0;				/* Migration rate between demes */	
	double lambdami = 0;		/* Average length of MITOTIC GC event */
	double lambdame = 0;		/* Average length of MEIOTIC GC event */	
	double bigQmi = 1/(1.0*0);	/* Relative length of lambdami (for GC prob calculations) */	
	double bigQme = 1/(1.0*0);	/* Relative length of lambdame (for GC prob calculations) */		
	double pburst = 0;			/* Probability that subsequent mutations will cluster */
	double bdist = 0;			/* Max distance if mutations cluster */
	double alpha = 0;			/* Population growth/shrink parameter */
	char Tout[32];				/* String to hold filename in (Trees) */
	FILE *ofp_tr = NULL;		/* Pointer for tree output */
	FILE *ofp_sd;				/* Pointer for seed output */
		
	/* GSL random number definitions */
	const gsl_rng_type * T;
	gsl_rng * r;
	
	/* Reading in data from command line */
	if(argc < 6){
		fprintf(stderr,"At least 5 inputs are needed.\n");
		usage();
	}
	
	N = atoi(argv[1]);
	if(N <= 0){
		fprintf(stderr,"Total Population size N is zero or negative, not allowed.\n");
		usage();
	}
	
	IwithT = atoi(argv[2]);
	IbetT = atoi(argv[3]);
	Itot = 2*IwithT + IbetT;
	Iindv = IwithT + IbetT;
	if(Itot <= 1){
		fprintf(stderr,"More than one sample must be initially present to execute the simulation.\n");
		usage();
	}
	
	bsex = strtod(argv[4],NULL);
	if( bsex <= 0 || bsex > 1){
		fprintf(stderr,"Baseline rate of sexual reproduction has to lie between 0 and 1.\n");
		usage();
	}

	/* Number of samples/reps to take */
	Nreps = atoi(argv[5]);
	if(Nreps <= 0){
		fprintf(stderr,"Must set positive number of repetitions.\n");
		usage();
	}
	
	/*
	Initial state of two samples;
	Iwith = No. of within host
	Ibet = No. of between host
	So total number of initial samples in a deme = 2*Iwith + Ibet
	THEN sexL[i], sexH[i] = low-sex, high-sex rate in each deme
	*/
	
	unsigned int *Iwith = calloc(1,sizeof(unsigned int));		/* Within-individual samples */
	unsigned int *Ibet = calloc(1,sizeof(unsigned int));		/* Between-individual samples */
	unsigned int *Nwith = calloc(1,sizeof(unsigned int));		/* To be used in individual rep */
	unsigned int *Nbet = calloc(1,sizeof(unsigned int));		/* To be used in individual rep */
	unsigned int *zeros = calloc(1,sizeof(unsigned int));		/* Placeholder array of zeros */
	unsigned int *demes = calloc(1,sizeof(unsigned int));		/* Indices of demes (for sampling) */
	double *sexL = calloc(1,sizeof(double));					/* Low-sex rates */	
	double *sexH = calloc(1,sizeof(double));					/* High-sex rates */
	
	*(Iwith + 0) = IwithT;
	*(Ibet + 0) = IbetT;
	*(sexL + 0) = bsex;
	*(sexH + 0) = 0;
	*(zeros + 0) = 0;
	*(demes + 0) = 0;
	
	/* Now moving onto case-by-case parameter setting */
	if(argc > 6){
		while(argx < argc){
			if(argv[argx][0] != '-'){
				fprintf(stderr,"Further arguments should be of form -C, for C a defined character.\n");
				usage();
			}
			switch ( argv[argx][1] ){
				case 't':		/* Defining Mutation Rate (theta) */
					ismut = 1;
					argx++;
					theta = strtod(argv[argx],NULL);
					if(theta <= 0){
						fprintf(stderr,"Mutation rate must be a positive value.\n");
						usage();
					}
					argx++;
					break;
				case 'r':		/* Defining recombination (cross over) */
					isrec = 1;
					argx++;
					rec = strtod(argv[argx],NULL);
					if(rec < 0){
						fprintf(stderr,"Must define a positive recombination rate.\n");
						usage();
					}
					argx++;
					nsites = atoi(argv[argx]);
					if(nsites < 2 && rec > 0){
						fprintf(stderr,"Must define at least two sites with recombination.\n");
						usage();
					}
					argx++;
					break;
				case 'T': 		/* Print trees? */
					argx++;
					prtrees = 1;
					break;
				case 'c':		/* MEIOTIC gene conversion */
					if(isrec == 0){
						fprintf(stderr,"Have to define number of sites first (using -r) before defining conversion rates.\n");
						usage();
					}
					argx++;
					gme = strtod(argv[argx],NULL);
					if(gme < 0){
						fprintf(stderr,"Must define a positive meiotic gene conversion rate.\n");
						usage();
					}
					if(nsites < 2 && gme > 0){
						fprintf(stderr,"Must define at least two sites with meiotic gene conversion.\n");
						usage();
					}
					argx++;
					lambdame = strtod(argv[argx],NULL);
					if(lambdame < 1){
						fprintf(stderr,"With meiotic gene conversion, average length (lambda) has to be at least 1.\n");
						usage();
					}
					bigQme = nsites/(1.0*lambdame);
					if(bigQme < WARNQ){
						fprintf(stderr,"WARNING: Simulation assumes that lambda ~ nsites.\n");
						fprintf(stderr,"Too large a lambda might produce inaccurate coalescent times.\n");
					}
					argx++;
					break;
				case 'm':		/* MITOTIC gene conversion */
					if(isrec == 0){
						fprintf(stderr,"Have to define number of sites first (using -r) before defining conversion rates.\n");
						usage();
					}
					argx++;
					gmi = strtod(argv[argx],NULL);
					if(gmi < 0){
						fprintf(stderr,"Must define a positive mitotic gene conversion rate.\n");
						usage();
					}
					argx++;
					lambdami = strtod(argv[argx],NULL);
					if(lambdami < 1){
						fprintf(stderr,"With mitotic gene conversion, average length (lambda) has to be at least 1.\n");
						usage();
					}
					bigQmi = nsites/(1.0*lambdami);
					if(bigQmi < WARNQ){
						fprintf(stderr,"WARNING: Simulation assumes that lambda ~ nsites.\n");
						fprintf(stderr,"Too large a lambda might produce inaccurate coalescent times.\n");
					}
					argx++;
					break;
				case 'I':		/* Population subdivision (Island model) */
					ismig = 1;
					argx++;
					d = atoi(argv[argx]);
					if(d <= 0){
						fprintf(stderr,"Number of demes has to be a positive integer.\n");
						usage();
					}
					if(N%d != 0){
						fprintf(stderr,"Population size must be a multiple of deme number.\n");
						usage();
					}
					N = (unsigned int)(N/(d*1.0));	/* Scaling NT to a demetic size, for more consistent use in calculations */
					if(d > 1){		/*	Assigning per-deme samples */
					
						Iwith = (unsigned int *)realloc(Iwith,d*sizeof(unsigned int));
						Ibet = (unsigned int *)realloc(Ibet,d*sizeof(unsigned int));
						Nwith = (unsigned int *)realloc(Nwith,d*sizeof(unsigned int));
						Nbet = (unsigned int *)realloc(Nbet,d*sizeof(unsigned int));
						zeros = (unsigned int *)realloc(zeros,d*sizeof(unsigned int));
						demes = (unsigned int *)realloc(demes,d*sizeof(unsigned int));
						sexL = (double *)realloc(sexL,d*sizeof(double));
						sexH = (double *)realloc(sexH,d*sizeof(double));
	
						for(x = 0; x < d; x++){
							argx++;
							*(Iwith + x) = atoi(argv[argx]);
							argx++;
							*(Ibet + x) = atoi(argv[argx]);
							*(sexL + x) = bsex;
							*(sexH + x) = 0;
							*(zeros + x) = 0;
							*(demes + x) = x;
							IwithT2 += (*(Iwith + x));
							IbetT2 += (*(Ibet + x));
						}
						if( IwithT2 != IwithT ){
							fprintf(stderr,"Sum of paired per-deme samples should equal total samples.\n");
							usage();
						}
						if( IbetT2 != IbetT ){
							fprintf(stderr,"Sum of unpaired per-deme samples should equal total samples.\n");
							usage();
						}
					}else if(d == 1){
						argx += 2;
					}
					argx++;
					mig = strtod(argv[argx],NULL);
					mig = mig/(2.0*N*d);
					if(mig < 0){
						fprintf(stderr,"Migration rate must be a positive (or zero) value.\n");
						usage();
					}
					if((d > 1) && (mig == 0)){
						fprintf(stderr,"Migration rate cannot be zero with multiple demes.\n");
						usage();
					}
					argx++;
					break;	
				case 'H':
					if(ismig == 0){
						fprintf(stderr,"Have to define number of demes first (using -I) before defining heterogeneity.\n");
						usage();
					}
					
					/* Defining type of heterogeneity */
					argx++;
					pSTIN = atoi(argv[argx]);
					argx++;
					pLH = strtod(argv[argx],NULL);
					argx++;
					pHL = strtod(argv[argx],NULL);
					
					if(pSTIN != 0 && pSTIN != 1 && pSTIN != 2){
						fprintf(stderr,"pSTIN has to equal 0, 1 or 2.\n");
						usage();
					}
					if(pSTIN != 1){
						if(pHL < 0 || pHL > 1 || pLH < 0 || pLH > 1){
							fprintf(stderr,"Sex transition probabilities have to lie between 0 and 1 (if there is no stepwise change in sex).\n");
							usage();
						}
					}
					if( (pHL == 0 || pLH == 0 ) && pSTIN == 0){
						fprintf(stderr,"Sex transition probabilities have to lie between 0 and 1 with fluctuating sex.\n");
						usage();
					}
					
					for(x = 0; x < d; x++){
						argx++;
						*(sexL + x) = strtod(argv[argx],NULL);
						argx++;						
						*(sexH + x) = strtod(argv[argx],NULL);

						if( *(sexL + x) < 0 || *(sexL + x) > 1 || *(sexH + x) < 0 || *(sexH + x) > 1){
							fprintf(stderr,"Rate of sexual reproduction has to lie between 0 and 1.\n");
							usage();
						}
						if( *(sexH + x) == 0 && pSTIN == 0){
							fprintf(stderr,"With temporally heterogeneous sex, high rate cannot be zero.\n");
							usage();
						}
						if((*(sexL + x) > *(sexH + x)) && pSTIN == 0){
							fprintf(stderr,"All low-sex values must be less than or equal to high-sex values. Please re-check.\n");
							usage();
						}
					}
					argx++;
					break;
				case 'P':
					ismsp = 1;
					argx++;
					break;
				case 'b':
					if(isrec == 0){
						fprintf(stderr,"Have to define number of breakpoints first (using -r) before defining burst of mutations.\n");
						usage();
					}
					if(ismut == 0){
						fprintf(stderr,"Need to define a mutation rate first.\n");
						usage();
					}					
				
					isburst = 1;
					argx++;
					pburst = strtod(argv[argx],NULL);
					argx++;
					mburst = atoi(argv[argx]);
					argx++;
					bdist = strtod(argv[argx],NULL);
					
					if( (pburst < 0) || (pburst > 1) ){
						fprintf(stderr,"Burst probability has to lie between 0 and 1.\n");
						usage();
					}
					if(mburst == 0){
						fprintf(stderr,"Max burst size has to be greater than zero.\n");
						usage();
					}
					if(bdist == 0){
						fprintf(stderr,"Burst distance has to be greater than zero.\n");
						usage();
					}
					argx++;
					break;
				case 'G':
					isexp = 1;
					argx++;
					alpha = strtod(argv[argx],NULL);
					if(alpha <= 0){
						fprintf(stderr,"Growth rate has to be a positive, non-zero value.\n");
						usage();
					}
					argx++;					
					break;
				default:	/* If none of these cases chosen, exit with error message */
					fprintf(stderr,"Error: Non-standard input given.\n");
					usage();
					break;
			}
		}
	}
	
	if(ismut == 0 && prtrees == 0){
		fprintf(stderr,"Error: Have to define a mutation rate or print out trees.\n");
		usage();
	}
	
	if(isrec == 1){
		rec = rec/(2.0*(nsites-1)*N*d);
		gme = ((gme*nsites)/(2.0*N*d*(nsites + 1)));
		gmi = ((gmi*nsites)/(2.0*N*d*(nsites + 1)));
	}
	
	if(isburst == 1){
		bdist = bdist/(1.0*nsites);
	}
	
	if(d == 1){
		mig = 0;	/* Set migration to zero if only one deme, as a precaution */
	}
	if(rec == 0 && gmi == 0 && gme == 0 && isburst == 0){
		nsites = 1; /* Set number of sites to 1 if no recombination OR gc, as a precaution */
	}
	if(nsites == 1){
		rec = 0;
	}
	
	/* Arrays definition and memory assignment */
	unsigned int *nlrec = calloc(d,sizeof(unsigned int));			/* Non-recombinable samples 1 */
	unsigned int *nlrec2 = calloc(d,sizeof(unsigned int));			/* Non-recombinable samples 2 */
	unsigned int *evsex = calloc(d,sizeof(unsigned int));			/* Number of sex events per deme */
	unsigned int *csex = calloc(2,sizeof(unsigned int));			/* Does sex occur or not? */
	unsigned int *draw = calloc(12,sizeof(unsigned int));			/* Event that happens */
	unsigned int *draw2 = calloc(d,sizeof(unsigned int));			/* Deme in which event happens */
	unsigned int *draw3 = calloc(2,sizeof(unsigned int));			/* Which type of sample migrates */	
	double *Nsamps = calloc(2,sizeof(double));						/* Within and between-indv samples in deme */
	int *WCH = calloc(d,sizeof(int));								/* How within-indv samples change */
	int *BCH = calloc(d,sizeof(int));								/* How between-indv samples change */
	int *WCHex = calloc(d,sizeof(int));								/* Extra calc if null values formed */
	int *BCHex = calloc(d,sizeof(int));								/* Extra calc if null values formed */
	double *sexC = calloc(d,sizeof(double));						/* Current rates of sex per deme */	
	double *sexCInv = calloc(d,sizeof(double));						/* Inverse of current rates of sex (1-sexC) */
	double *psex = calloc(2,sizeof(double));						/* Individual probabilities if individuals undergo sex or not */
	double *pr_rsums = calloc(12,sizeof(double));					/* Rowsums of probs (for event choosing) */
	double **pr = calloc(12,sizeof(double *));						/* Probability matrix per deme */
	for (j = 0; j < 12; j++){										/* Assigning space for each population within each deme */
		pr[j] = calloc(d,sizeof(double));
	}
	  
	/* create a generator chosen by the 
    environment variable GSL_RNG_TYPE */
     
	gsl_rng_env_setup();
	if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	ofp_sd = fopen("Seed.dat","w");
	fprintf(ofp_sd,"%lu\n",gsl_rng_default_seed);
	fclose(ofp_sd);
	if(ismsp == 1){
		printf("%lu\n",gsl_rng_default_seed);
	}
	
	/* Creating necessary directories */
	if((rec > 0 || gmi > 0 || gme > 0) && (nsites != 1) && (prtrees == 1)){
		mkdir("Trees/", 0777);
	}
	if(ismut == 1){
		mkdir("Mutations/", 0777);
	}
	
	/* Running the simulation Nreps times */
	for(i = 0; i < Nreps; i++){
		
/*		printf("Starting Run %d\n",i);	*/
		nmutT = 0;

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
		NMax = Itot;		
		Nindv = Iindv;
		exr = INITBR;
	    exc = INITBR;
	
		for(x = 0; x < d; x++){
			*(Nwith + x) = *(Iwith + x);	/* Resetting number of within-host samples */
			*(Nbet + x) = *(Ibet + x);		/* Resetting number of between-host samples */
			*(nlrec + x) = 0;				/* Genome not affected by recombination */
			*(nlrec2 + x) = 0;				/* Genome not affected by recombination */
		}
		
		/* Setting up temporal heterogeneity */
		rate_change(N,pST,pLH,pHL,sexH,sexL,0,sexC,sexCInv,&tts,&npST,r);
		pST = npST;
		tls = 0;
		
		/* Setting up summary table of individual samples */
		/* ASSIGNING MEMORY FROM SCRATCH HERE, SINCE TABLES WILL BE MODIFIED FOR EACH SIM */
		
		unsigned int **indvs = calloc(Itot+exr,sizeof(unsigned int *));		/* Table of individual samples */
		int **GType = calloc(Itot+exr,sizeof(int *));						/* Table of sample genotypes */
		double **CTms = calloc(Itot+exr,sizeof(double *));					/* Coalescent times per sample */
		int **TAnc = calloc(Itot+exr,sizeof(int *));						/* Table of ancestors for each sample */
		unsigned int **breaks = calloc(2,sizeof(unsigned int *));			/* Table of breakpoints created in the simulation */
		double **TFin = calloc((Itot-1),sizeof(double *));					/* Final ancestry table, for tree reconstruction */
		unsigned int **nlri = calloc(Itot+exr,sizeof(unsigned int *));		/* Per-sample valid bp table (for weighted sampling) */
		for(j = 0; j < (Itot+exr); j++){									/* Assigning space for each genome sample */
			indvs[j] = calloc(4,sizeof(unsigned int));
			GType[j] = calloc(exc+1,sizeof(int));
			nlri[j] = calloc(4,sizeof(unsigned int));			
			if(j < Itot){
				CTms[j] = calloc(exc+1,sizeof(double));
				TAnc[j] = calloc(exc+1,sizeof(int));
				if(j < (Itot - 1)){
					TFin[j] = calloc(3,sizeof(double));
				}
			}
		}
		breaks[0] = calloc(exc,sizeof(unsigned int));
		breaks[1] = calloc(exc,sizeof(unsigned int));
		nbreaks = 1;
				
		IwithC = 0;
		IbetC = 0;
		for(j = 0; j < Itot; j++){
			*((*(indvs + j)) + 0) = j;
			*((*(GType + j)) + 0) = j;
			*((*(GType + j)) + 1) = j;
			*((*(CTms + j)) + 0) = j;
			*((*(CTms + j)) + 1) = -1;			
			*((*(TAnc + j)) + 0) = j;
			/* Entry of "-1" equivalent to NA in old code, 
			i.e. it is not ancestral, reflects extant tracts */
			*((*(TAnc + j)) + 1) = (-1);
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
				*((*(nlri + j)) + 0) = j;
				*((*(nlri + j)) + 1) = 0;
				*((*(nlri + j)) + 2) = (nsites-1);
			}
			for(x = 0; x < d; x++){
				for(j = IbetC; j < (IbetC + *(Ibet + x)); j++){
					*((*(indvs + (2*IwithT + j))) + 3) = x;
					*((*(nlri + j)) + 3) = x;
				}
				IbetC += *(Ibet + x);
			}
		}
		done = 0;
		while(done != 1){
					
			/* Setting up vector of state-change probabilities *without sex* */
			probset2(N, gmi, gme, sexC, rec, bigQmi, bigQme, nsites, nlrec, zeros, mig, Nwith, Nbet, zeros, 0, pr);
			nosex = powDUI(sexCInv,Nwith,d);				/* Probability of no segregation via sex, accounting for within-deme variation */
			psum = (1-nosex) + nosex*(sumT_D(pr,12,d));		/* Sum of all event probabilities, for drawing random time */
							
			/* Intermediate error checking */
			if(psum > 1){
				proberr2(pr);
			}
			if((psum <= 0) && (isallD(sexC,d,0) != 1)){
				proberr3(pr);
			}
			if(isanylessD_2D(pr,12,d,0) == 1){
				proberr(0, pr);
			}
			
			/* Drawing time to next event, SCALED TO 2NT GENERATIONS */
			if(psum == 1){
				tjump = 0.0;
			}else if(psum == 0){
				tjump = (1.0/0.0);
			}else{
				tjump = (gsl_ran_geometric(r,psum))/(2.0*N*d);
			}
			NextT = (Ttot + tjump);
			
			/* Outcomes depends on what's next: an event or change in rates of sex	*/
			if(NextT > (tls + tts)){ 	/* If next event happens after a switch, change rates of sex */
				tls = (tls + tts);		/* 'Time since Last Switch' or tls	*/
				Ttot = tls;
				rate_change(N,pST,pLH,pHL,sexH,sexL,1,sexC,sexCInv,&tts,&npST,r);
				pST = npST;
			}else if (NextT <= (tls + tts)){	/* If next event happens before a switch, draw an action	*/
				Ttot = NextT;

				/* Determines if sex occurs; if so, which samples are chosen */
				/* (deme-independent Binomial draws) */
				esex = 0;
				CsexS = 0;
				vcopyUI(evsex,zeros,d);
				*(psex + 0) = nosex*(sumT_D(pr,12,d));
				*(psex + 1) = 1-nosex;
				gsl_ran_multinomial(r,2,1,psex,csex);
				if(*(csex + 1) == 1){				/* Working out number of sex events IF it does occur */
					while(CsexS == 0){
						for(x = 0; x < d; x++){
							*(evsex + x) = gsl_ran_binomial(r,*(sexC + x),*(Nwith + x));
							esex += *(evsex + x);
						}
						if(esex > 0){
							CsexS = 1;
						}
					}
				}
				
				unsigned int *rsex = calloc(esex,sizeof(unsigned int));
				/* Now redrawing probabilities with changed configuration */
				if(esex >= 1){
					/* New in ARG simulation: 
					already determining which samples have split 
					(so can calculate recombination prob accurately) */
					
					/* First, choosing samples to split by sex */
					sexsamp(indvs, rsex, evsex, Nwith, Ntot, r);
					
					/* Then calculating relevant breakpoints in each new sample */					
					reccal(indvs, GType, breaks, nlri, Nbet, Nwith, rsex, esex, nlrec2, nbreaks, NMax, 1, i);
					
					/* Then recalculating probability of events */
					probset2(N, gmi, gme, sexC, rec, bigQmi, bigQme, nsites, nlrec, nlrec2, mig, Nwith, Nbet, evsex, 1, pr);
					if(isanylessD_2D(pr,12,d,0) == 1){
						proberr(1, pr);
					}
				}
				
				/* Given event happens, what is that event? 
				Weighted average based on above probabilities. 
				Then drawing deme of event. */
				rowsumD(pr,12,d,pr_rsums);
				gsl_ran_multinomial(r,12,1,pr_rsums,draw);			
				event = matchUI(draw,12,1);
				gsl_ran_multinomial(r,d,1,(*(pr + event)),draw2);
				deme = matchUI(draw2,d,1);
				
/*				printf("Event is %d\n",event);	*/

				if(event == 9){		/* Choosing demes to swap NOW if there is a migration */
					stchange2(event,deme,evsex,WCH,BCH);
					vsum_UI_I(Nwith, WCH, d);
					vsum_UI_I(Nbet, BCH, d);
					Ntot = 2*(sumUI(Nwith,d)) + sumUI(Nbet,d);
					*(Nsamps + 0) = *(Nwith + deme);
					*(Nsamps + 1) = *(Nbet + deme);
					drec = deme;
					while(drec == deme){
						gsl_ran_choose(r,&drec,1,demes,d,sizeof(unsigned int));
					}
					gsl_ran_multinomial(r,2,1,Nsamps,draw3);
					e2 = matchUI(draw3,2,1);
					if(e2 == 0){	/* Paired sample migrates */
						(*(Nwith + deme))--;
						(*(Nwith + drec))++;
					}else if(e2 == 1){	/* Single sample migrates */
						(*(Nbet + deme))--;
						(*(Nbet + drec))++;
					}
				}
				
				/* Changing ancestry accordingly */
				/* MOVED CODE HERE TO AVOID ERRORS WHEN NUMBER OF SAMPLES MISMATCH */
				gcalt = 0;
				achange = coalesce(indvs, GType, CTms, TAnc, nlri, Ttot, Nwith, Nbet, deme, rsex, evsex, event, drec, e2, breaks, nsites, &nbreaks, NMax, Itot, lambdami, lambdame, gmi, gme, bigQmi, bigQme, &gcalt, sexC, WCHex, BCHex, r, i);
				
				/* Based on outcome, altering (non-mig) states accordingly */
				if(event != 9){		/* Since already done for state = 9 above... */
					stchange2(event,deme,evsex,WCH,BCH);
					vsum_UI_I(Nwith, WCH, d);
					vsum_UI_I(Nbet, BCH, d);
					Ntot = 2*(sumUI(Nwith,d)) + sumUI(Nbet,d);
					*(Nsamps + 0) = *(Nwith + deme);
					*(Nsamps + 1) = *(Nbet + deme);
					if(gcalt != 0){
						if(gcalt == 1){		/* GC led to new BH sample being produced */
							(*(Nwith + deme))++;
							(*(Nsamps + 0))++;
							(*(Nbet + deme))--;
							(*(Nsamps + 1))--;				
							Ntot++;
							NMax++;
							if(NMax > HUGEVAL){
								manyr();			
							}
						}
						if(gcalt == 2){		/* GC led to WH sample coalescing */
							(*(Nbet + deme))++;
							(*(Nsamps + 1))++;				
							(*(Nwith + deme))--;
							(*(Nsamps + 0))--;
							Ntot--;
						}
					}

					if(achange == 1){
						vsum_UI_I(Nwith, WCHex, d);
						vsum_UI_I(Nbet, BCHex, d);
						vcopyUI_I(WCHex, zeros, d);
						vcopyUI_I(BCHex, zeros, d);
						Ntot = 2*(sumUI(Nwith,d)) + sumUI(Nbet,d);
						*(Nsamps + 0) = *(Nwith + deme);
						*(Nsamps + 1) = *(Nbet + deme);
						achange = 0;
					}
					
				}
				if(event == 10){
					NMax++;
					if(NMax > HUGEVAL){
						manyr();		
					}
				}
				
				/* Sorting table afterwards to ensure paired samples are together */
				indv_sort(indvs, NMax);

				/* Updating baseline recombinable material depending on number single samples */
				if(isallUI(*(breaks+1),nbreaks,1,0) == 0){			
					reccal(indvs, GType, breaks, nlri, Nbet, Nwith, rsex, esex, nlrec, nbreaks, NMax, 0, i);
					for(x = 0; x < d; x++){
						*(nlrec2 + x) = 0;
					}
				}
				free(rsex);		/* Can be discarded once used to change ancestry */

				/*
				printf("Lrec is %d, Nbet is %d\n",lrec,sumUI(Nbet,d));
				for(x = 0; x < d; x++){
					printf("nlrec in deme %d is %d\n",x,*(nlrec + x));
				}
				printf("\n");
				if(i == 3 && event == 10){
					TestTabs(indvs, GType, CTms, TAnc, breaks, nlri, NMax, Itot, sumUI(Nbet,d), nbreaks);
				}						
				
				if(i == 0){
					TestTabs(indvs, GType, CTms, TAnc, breaks, nlri, NMax, Itot, sumUI(Nbet,d), nbreaks);
				}

				*/
																							
				/* Checking if need to expand tables */
				if(NMax == (exr+Itot-1)){
					exr += INITBR;
					indvs = (unsigned int **)realloc(indvs,(Itot+exr)*sizeof(unsigned int *));
					GType = (int **)realloc(GType, (Itot+exr)*sizeof(int *));
					nlri = (unsigned int **)realloc(nlri,(Itot+exr)*sizeof(unsigned int *));
					for(j = 0; j < (Itot+exr-INITBR); j++){
						indvs[j] = (unsigned int *)realloc(*(indvs + j),4*sizeof(unsigned int));
						GType[j] = (int *)realloc( *(GType + j) ,(exc + 1)*sizeof(int));
						nlri[j] = (unsigned int *)realloc(*(nlri + j),4*sizeof(unsigned int));
					}
					for(j = (Itot+exr-INITBR); j < (Itot+exr); j++){
						indvs[j] = (unsigned int *)calloc(4,sizeof(unsigned int));
						GType[j] = (int *)calloc((exc + 1),sizeof(int));
						nlri[j] = (unsigned int *)calloc(4,sizeof(unsigned int));						
					}
				}
				
				if(nbreaks >= exc-1){
					exc += INITBR;
					for(j = 0; j < (Itot+exr); j++){
						GType[j] = (int *)realloc( *(GType + j) ,(exc + 1)*sizeof(int));
						if(j < Itot){
							CTms[j] = (double *)realloc( *(CTms + j) ,(exc + 1)*sizeof(double));
							TAnc[j] = (int *)realloc( *(TAnc + j) ,(exc + 1)*sizeof(int));
						}
					}
					breaks[0] = (unsigned int *)realloc(*(breaks + 0),exc*sizeof(unsigned int));
					breaks[1] = (unsigned int *)realloc(*(breaks + 1),exc*sizeof(unsigned int));
				}
				
				/* Testing if all sites coalesced or not */
				done = isallUI(*(breaks + 1),nbreaks,1,0);
			}
		}
		
		if(ismsp == 1){
	    	printf("\n");
			printf("// \n");
		}		

		for(x = 1; x <= nbreaks; x++){
			
			/* Creating ancestry table */
			count = 0;
			for(j = 0; j < Itot; j++){
				if((*((*(CTms + j)) + x)) != (-1.0)){
					*((*(TFin + count)) + 0) = *((*(CTms + j)) + 0);
					if(isexp == 0){
						*((*(TFin + count)) + 1) = *((*(CTms + j)) + x);
					}else if(isexp == 1){
						*((*(TFin + count)) + 1) = (1/(1.0*alpha))*log(1 + alpha*(*((*(CTms + j)) + x)));
					}
					*((*(TFin + count)) + 2) = *((*(TAnc + j)) + x);
					count++;
				}
			}
			indv_sortD(TFin,(Itot-1),3,1);

			/* Using ancestry table to build tree and mutation table */
			if(x < nbreaks){
				maxd2 = (*((*(breaks + 0)) + x));
				maxd = maxd2/(1.0*nsites);
				mind2 = (*((*(breaks + 0)) + (x-1)));
				mind = mind2/(1.0*nsites);
			}else if(x == nbreaks){
				maxd2 = nsites;
				maxd = 1;
				mind2 = *((*(breaks + 0)) + (x-1));
				mind = (mind2)/(1.0*nsites);
			}
			char *ret_tree = treemaker(TFin, theta*(maxd-mind), mind2, maxd2 ,mind, maxd, Itot, i, gmi, gme, ismsp, &nmutT, prtrees, ismut, pburst, mburst, bdist, r);
			if(prtrees == 1){
				if((rec == 0 && gmi == 0 && gme == 0) || (nsites == 1) ){
					if(i == 0){
						ofp_tr = fopen("Trees.dat","w");
					}else if(i > 0){
						ofp_tr = fopen("Trees.dat","a+");
					}
					fprintf(ofp_tr,"%s\n",ret_tree);
				}else if((rec > 0 || gmi > 0 || gme > 0) && (nsites != 1)){
					sprintf(Tout,"Trees/Trees_%d.dat",i);
					if(x == 1){
						ofp_tr = fopen(Tout,"w");
					}else if(x > 1){
						ofp_tr = fopen(Tout,"a+");
					}
					fprintf(ofp_tr,"%s\n",ret_tree);
				}
				fclose(ofp_tr);
			}
			free(ret_tree);
		}
		
		/* Freeing memory at end of particular run */		
		free(breaks[1]);
		free(breaks[0]);
		for(j = 0; j < (Itot + exr); j++){
			if(j < Itot){
				free(TAnc[j]);
				free(CTms[j]);
				if(j < (Itot - 1)){
					free(TFin[j]);
				}
			}
			free(nlri[j]);			
			free(GType[j]);
			free(indvs[j]);
		}
		free(nlri);
		free(TFin);
		free(breaks);
		free(TAnc);
		free(CTms);
		free(GType);
		free(indvs);
	}
	
	/* Freeing memory and wrapping up */
 	gsl_rng_free(r);
 	for(x = 0; x < 12; x++){
		free(pr[x]);
	}
	free(pr);
	free(pr_rsums);
	free(psex);
	free(sexCInv);	
	free(sexC);
	free(nlrec2);
	free(nlrec);
 	free(sexH);
 	free(sexL);
 	free(BCHex);
 	free(WCHex);
 	free(BCH);
 	free(WCH);
 	free(Nsamps);
 	free(draw3);
 	free(draw2);
 	free(draw);
 	free(csex);
 	free(evsex);
  	free(demes);	
 	free(zeros);
 	free(Nbet);
 	free(Nwith);
 	free(Ibet);
 	free(Iwith); 	 	 	
 	
	return 0;
	
}	/* End of main program */

/* End of File */