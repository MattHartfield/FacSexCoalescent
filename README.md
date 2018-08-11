README FOR FAC SEX COALESCENT

FacSexCoalescent is a simulation program, written in C, to simulate genealogies from individuals that alternate between sexual and asexual (parthenogenetic) reproduction. It is the source code for the Hartfield, Wright and Agrawal manuscript "Coalescence and linkage disequilibrium in facultatively sexual diploids". This updated version of the program can simulate genealogies at multiple sites, and is much faster since it is written in C (the 2015 program was written in R and only simulated genealogies at a single site).

COMPILATION:
This is a program written in C and needs to be compiled before execution. Compile using a program like gcc using the following command:

gcc -o FacSexCoalescent -lm -lgsl -lgslcblas -I/usr/local/include -L/usr/local/lib FacSexCoalescent.c

Simulation uses routines found with the GNU Scientific Library (GSL) (http://www.gnu.org/software/gsl/) Since GSL is distributed under the GNU General Public License (http://www.gnu.org/copyleft/gpl.html), you must download it separately from this file.

EXECUTION:
For basic usage, run the program without any input parameters to obtain the basic command line.
Otherwise, more information is available in the user documentation.

Comments to matthew.hartfield@birc.au.dk.
