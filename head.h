#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <dlfcn.h>

#define PI 3.141593

// structures

typedef struct result{
	double survival;
	double Va;
	double Vg;
} result;

typedef struct model{
	int id;
	double mod[3][3];
	double modg[4][4];
} model;

// input globals
int USEED, UNID, UNMOD, UMAXG, UGTHRESH, UNFREQS1, UNFREQS2, UNBESTMODS, UMAXRUNS, *UNEXTMODS, UNPERT, URANGE, USTARTVA;
double *UFREQS1, *UFREQS2, *UPERT, UOVERRIDE[3][3], UINIT[3][3];
char UFILENAME[50];
result ***URES;

// function declarations
void inputsummary(char *inputfile);
double **get_space_f(int m, int n);
void release_space_f(double **X);
double ***make_cube_d(int x, int y, int z);
result ***make_cube_r(int x, int y, int z);
void free_cube_d(double ***array);
void free_cube_r(result ***array);
void rnorm(double *x, int n);
void inputsummary(char *inputfile);
void convert_modg(double mod[3][3], double modg[4][4]);
void gfreq(double A, double B, double Rsq, double C[4]);
double calc_vg(double mod[3][3], double Cf[4]);
double calc_va(double mod[3][3], double Cf[4]);
void nextgeneration(double modg[4][4], double *Cf, double *Cn);
int checkfix(double Cf[4], double n);
void override_pattern(double mod[3][3], double override[3][3]);
void init_pattern(double mod[3][3], double modg[4][4]);
void mutate_pattern(double modo[3][3], double mod[3][3], double modg[4][4]);
void iteration(model *X);
void optimise_ga(int *bestmods, double *max);
void print_bestmods(int *bestmods, model *X, int iter, double *max);
void mutate(model *X, int *bestmods);
void perform_ga(int maxruns);
void printmodels(model *X);