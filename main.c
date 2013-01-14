/* Gib Hemani
 * Genetic algorithm to generate epistatic interactions that survive for a long time whilst making large additive effects
 * 1/12/2010
 */

/* input:
UNID
UMAXG
UMAXRUNS
UGTHRESH
USTARTVA
URMIN URMAX
UNFREQS
UNFREQ[0] UNFREQ[1] ...
UNMOD
UNBESTMODS
UNEXTMODS[0] UNEXTMODS[1] ...
UNPERT
UPERT[0] UPERT[1] ...
USEED
UFILENAME
*/


/*
Genetic algorithm:

Maximise additive variance and survival time
- Create UNMOD random patterns with values in the range URMIN-URMAX
- frequencies - UFREQS1 UFREQS2
- UNFREQS1 x UNFREQS2 runs, if > UGTHRESH survive for at least UMAXG generations then eligable to continue
- If none survive then choose a new random set
- Additive variance is summed from all the runs of eligable patterns starting from generation USTARTVA, top UNBESTMOD proceed
- Each of the two has UNEXTMODS[i]-1 mutant versions (plus original version), then make UNEXTMODS[UNBESTMODS-1] new random patterns
- To mutate, each element is changed by randomly sampling the array UPERT
- All patterns are scaled to be within range URMIN-URMAX


record-
pattern
survival time for each frequency
total variance for each frequency
*/

#include "head.h"


double **get_space_f(int m, int n)
{
	int i;
	double *p, **a;
	p = malloc(m * n * sizeof(double));
	a = malloc(m * sizeof(double *));
	for(i = 0; i < m; i++) a[i] = p + (i * n);
	return a;
}

void release_space_f(double **X)
{
	double * p;
	p = (double *) X[0];
	free(p);
	free(X);
}

double ***make_cube_d(int x, int y, int z)
{
	int i,j;
	double *p, **p1, ***a;
	p = malloc(x * y * z * sizeof(double));
	p1 = malloc(x * y * sizeof(double *));
	a = malloc(x * sizeof(double **));
	
	for(i = 0; i < x; i++)
	{
		a[i] = p1 + i*y;
		for(j = 0; j < y; j++)
		{
			a[i][j] = p + i*y*z + j*y;
		}
	}
	return(a);
}

void free_cube_d(double ***array)
{
	double **p1, *p2;
	
	p1 = (double **) array[0];
	p2 = (double *) array[0][0];
	free(p2);
	free(p1);
	free(array);
}

result ***make_cube_r(int x, int y, int z)
{
	int i,j;
	result *p, **p1, ***a;
	p = malloc(x * y * z * sizeof(result));
	p1 = malloc(x * y * sizeof(result *));
	a = malloc(x * sizeof(result **));
	
	for(i = 0; i < x; i++)
	{
		a[i] = p1 + i*y;
		for(j = 0; j < y; j++)
		{
			a[i][j] = p + i*y*z + j*y;
		}
	}
	return(a);
}

void free_cube_r(result ***array)
{
	result **p1, *p2;
	
	p1 = (result **) array[0];
	p2 = (result *) array[0][0];
	free(p2);
	free(p1);
	free(array);
}

// Box-Muller method approximation. mean = 0; var = 1
// Marsaglia polar method may be faster
void rnorm(double *x, int n)
{
	int
		i;
	double
		tmp1f, tmp2f,
		x1,x2;

	for(i = 0; i < n; i++)
	{
		x1 = (double)rand() / RAND_MAX;
		x2 = (double)rand() / RAND_MAX;
		tmp1f = sqrtf(-2 * logf(x1));
		tmp2f = cosf(2*PI*x2);
		x[i] = tmp1f * tmp2f;
	}
}


void inputsummary(char *inputfile)
{
	int i,j;
	FILE *in;

	// READ INPUT

	in = fopen(inputfile, "r");

	fscanf(in,"%d\n",&UNID);
	fscanf(in,"%d\n",&UMAXG);
	fscanf(in,"%d\n",&UMAXRUNS);
	fscanf(in,"%d\n",&UGTHRESH);
	fscanf(in,"%d\n",&USTARTVA);
	fscanf(in,"%d\n",&URANGE);
	fscanf(in,"%d\n",&UNFREQS1);
	UFREQS1 = (double *)malloc(sizeof(double)*UNFREQS1);
	for(i = 0; i < UNFREQS1; i++)
	{
		fscanf(in,"%lf ",&UFREQS1[i]);
	}
	fscanf(in,"%d\n",&UNFREQS2);
	UFREQS2 = (double *)malloc(sizeof(double)*UNFREQS2);
	for(i = 0; i < UNFREQS2; i++)
	{
		fscanf(in,"%lf ",&UFREQS2[i]);
	}
	fscanf(in,"%d\n",&UNMOD);
	fscanf(in,"%d\n",&UNBESTMODS);
	UNEXTMODS = (int *)malloc(sizeof(int)*(UNBESTMODS+1));
	for(i = 0; i < (UNBESTMODS+1); i++)
	{
		fscanf(in,"%d ",&UNEXTMODS[i]);
	}
	fscanf(in,"%d\n",&UNPERT);
	UPERT = (double *)malloc(sizeof(double)*UNPERT);
	for(i = 0; i < UNPERT; i++)
	{
		fscanf(in,"%lf ",&UPERT[i]);
	}
	fscanf(in,"%d\n",&USEED);
	USEED ? srand(USEED) : srand(time(NULL));
	fscanf(in,"%s\n",UFILENAME);
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			fscanf(in,"%lf ",&UOVERRIDE[i][j]);
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			fscanf(in,"%lf ",&UINIT[i][j]);

	fclose(in);
	//URES = make_cube_r(UNMOD,UNFREQS1,UNFREQS2);
		URES = make_cube_r(UNMOD,10,10);
	printf("\n\n");
	printf("UNID = %d\n",UNID);
	printf("UNMOD = %d\n",UNMOD);
	printf("UMAXG = %d\n",UMAXG);
	printf("UMAXRUNS = %d\n",UMAXRUNS);
	printf("UGTHRESH = %d\n",UGTHRESH);
	printf("URANGE = %d\n",URANGE);
	printf("UNFREQS1 = %d\n",UNFREQS1);
	for(i = 0; i < UNFREQS1; i++)
	{
		printf("%f ",UFREQS1[i]);
	}
	printf("\nUNFREQS2 = %d\n",UNFREQS2);
	for(i = 0; i < UNFREQS2; i++)
	{
		printf("%f ",UFREQS2[i]);
	}
	printf("\nUBESTMODS = %d\n",UNBESTMODS);
	for(i = 0; i < (UNBESTMODS+1); i++)
	{
		printf("%d ",UNEXTMODS[i]);
	}
	printf("\nUNPERT = %d\n",UNPERT);
	for(i = 0; i < UNPERT; i++)
	{
		printf("%f ",UPERT[i]);
	}
	printf("\nUSEED = %d\n",USEED);
	printf("%s\n",UFILENAME);
	printf("override:\n");
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			printf("%0.2lf ",UOVERRIDE[i][j]);
		}
		printf("\n");
	}
	printf("initial:\n");
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			printf("%0.2lf ",UINIT[i][j]);
		}
		printf("\n");
	}
		
	printf("\n\n");
}

void convert_modg(double mod[3][3], double modg[4][4])
{
	int i,j;
	
	modg[0][0] = mod[0][0];
	modg[1][0] = mod[0][1];
	modg[1][1] = mod[0][2];
	modg[2][0] = mod[1][0];
	modg[2][1] = mod[1][1];
	modg[2][2] = mod[2][0];
	modg[3][0] = mod[1][1];
	modg[3][1] = mod[1][2];
	modg[3][2] = mod[2][1];
	modg[3][3] = mod[2][2];
	
	for(i = 0; i < 4; i++)
		for(j = 0; j < i; j++)
			modg[j][i] = modg[i][j];

}

void gfreq(double A, double B, double Rsq, double C[4])
{
	double prod;
	prod = sqrt(A*B*(1-A)*(1-B));
	C[0] = Rsq * prod + A*B;
	C[1] = A*(1-B) - Rsq * prod;
	C[2] = (1-A)*B - Rsq * prod;
	C[3] = Rsq * prod + (1-A)*(1-B);
}

double calc_vg(double mod[3][3], double Cf[4])
{
	int i,j;
	double fA[3], fB[3], Vg, meanW;
	
	fA[0] = (Cf[0]+Cf[1])*(Cf[0]+Cf[1]);
	fA[1] = 2*(Cf[0]+Cf[1])*(Cf[2]+Cf[3]);
	fA[2] = (Cf[2]+Cf[3])*(Cf[2]+Cf[3]);
	
	fB[0] = (Cf[0]+Cf[2])*(Cf[0]+Cf[2]);
	fB[1] = 2*(Cf[0]+Cf[2])*(Cf[1]+Cf[3]);
	fB[2] = (Cf[1]+Cf[3])*(Cf[1]+Cf[3]);

	Vg = 0;
	meanW = 0;
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
		{
			Vg += mod[i][j] * mod[i][j] * fA[i] * fB[j];
			meanW += mod[i][j] * fA[i] * fB[j];
		}
	Vg -= meanW * meanW;
	return Vg;
}

double calc_va(double mod[3][3], double Cf[4])
{
	int i;
	double U[3], V[3], Ga, Gb, D, X, Y, Ha, Hb, x, y, Va;
	
	x = Cf[0]+Cf[1];
	y = Cf[0]+Cf[2];
	for(i = 0; i < 3; i++)
	{
		U[i] = y*y*mod[i][0] + 2*y*(1-y)*mod[i][1] + (1-y)*(1-y)*mod[i][2];
		V[i] = x*x*mod[0][i] + 2*x*(1-x)*mod[1][i] + (1-x)*(1-x)*mod[2][i];
	}

	Ga = U[0]*(Cf[0]+Cf[1]) + U[1]*(1-2*Cf[0]-2*Cf[1]) - U[2]*(Cf[2]+Cf[3]);
	Gb = V[0]*(Cf[0]+Cf[2]) + V[1]*(1-2*Cf[0]-2*Cf[2]) - V[2]*(Cf[1]+Cf[3]);
	
	D = Cf[0] - (Cf[0]+Cf[1])*(Cf[0]+Cf[2]);
	X = (Cf[0]+Cf[1])*(Cf[2]+Cf[3]);
	Y = (Cf[0]+Cf[2])*(Cf[1]+Cf[3]);
	Ha = (Ga - D*Gb/X) / (1 - D*D/(X*Y));
	Hb = (Gb - D*Ga/Y) / (1 - D*D/(X*Y));
	
	Va = 2*(X*Ha*Ha + 2*Ha*Hb*D + Y*Hb*Hb);
	return Va;
}

void nextgeneration(double modg[4][4], double *Cf, double *Cn)
{
	int i,j;
	double meanW = 0, Wmarg[4], nu[4] = {1,-1,-1,1};
	
	
	for(i = 0; i < 4; i++)
	{
		Wmarg[i] = 0;
		for(j = 0; j < 4; j++)
		{
			Wmarg[i] += modg[i][j] * Cf[j];
		}
		meanW += Wmarg[i] * Cf[i];
	}

	for(i = 0; i < 4; i++)
	{
		Cn[i] = (Cf[i]*Wmarg[i] + nu[i]*0.5*modg[0][3]*(Cf[1]*Cf[2]-Cf[0]*Cf[3]))/meanW;
	}
}

int checkfix(double Cf[4], double n)
{
	double A, B, uthr, lthr, fixed = 0;
	
	//uthr = 1 - (1/(2*n));
	//lthr = 1/(2*n);
	uthr = 0.95;
	lthr = 0.05;
	A = Cf[0]+Cf[1];
	B = Cf[0]+Cf[2];
	if(A > uthr || A < lthr)
		fixed = 1;
	if(B > uthr || B < lthr)
		fixed = 1;

	return fixed;
}

void override_pattern(double mod[3][3], double override[3][3])
{
	int i,j;
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			mod[i][j] = override[i][j] >= 0 ? override[i][j] : mod[i][j];
		}
	}
}

void init_pattern(double mod[3][3], double modg[4][4])
{
	int i,j;
	double min, max;
	min = URANGE;
	max = 0;
	//double pat[3][3] = {{0,2,0},{1,1,1},{2,0,2}};
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			mod[i][j] = (double)(rand() % (URANGE+1))/URANGE;
			min = mod[i][j] < min ? mod[i][j] : min;
			max = mod[i][j] > max ? mod[i][j] : max;
		}
	}
	
	// scale values to stay within URMIN/URMAX and maximise range
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			mod[i][j] = (mod[i][j] - min)/(max-min);
	override_pattern(mod,UOVERRIDE);
	convert_modg(mod,modg);
}

void mutate_pattern(double modo[3][3], double mod[3][3], double modg[4][4])
{
	int i,j,k,temp;
	double min=URANGE, max=0, pert[9];
//	rnorm(pert,9);
	k = 0;
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			temp = rand() % UNPERT;
			mod[i][j] = modo[i][j] + UPERT[temp];
			//mod[i][j] = modo[i][j]*URMAX + pert[k++]*0.1;
		}
	}
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			min = mod[i][j] < min ? mod[i][j] : min;
			max = mod[i][j] > max ? mod[i][j] : max;
		}
	}

	// scale values to stay within URMIN/URMAX and maximise range
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			mod[i][j] = (mod[i][j] - min)/(max-min);

	override_pattern(mod,UOVERRIDE);
	convert_modg(mod,modg);
}

void iteration(model *X)
{
	int i,j,k,l,m,fixed;
	double Cf[4],Cn[4],sumvg,sumva;
	
	for(i = 0; i < UNMOD; i++)
	{
		for(j = 0; j < UNFREQS1; j++)
		{
			for(k = 0; k < UNFREQS2; k++)
			{
				gfreq(UFREQS1[j],UFREQS2[k],0,Cf);
				sumvg = 0;
				sumva = 0;
				fixed = 0;
				for(l = 0; l < UMAXG; l++)
				{
					if(l >= USTARTVA)
					{
						sumva += calc_va(X[i].mod,Cf);
						sumvg += calc_vg(X[i].mod,Cf);
					}
//					printf("%f %f\n",Cf[0]+Cf[1],Cf[0]+Cf[2]);
					nextgeneration(X[i].modg,Cf,Cn);
					fixed += checkfix(Cn,UNID);
					if(fixed) break;
					for(m = 0; m < 4; m++) Cf[m] = Cn[m];
				}
//				printf("\n");
				URES[i][j][k].survival = (double)l/UMAXG;
				URES[i][j][k].Va = sumva;
				URES[i][j][k].Vg = sumvg;
			}
		}
	}
}

void optimise_ga(int *bestmods, double *max)
{
	int i,j,k,*count;
	double *sumva;
	count = (int *)malloc(sizeof(int) * UNMOD);
	sumva = (double *)malloc(sizeof(double) * UNMOD);
	
	for(i = 0; i < UNMOD; i++)
	{
		count[i] = 0;
		sumva[i] = 0;
		for(j = 0; j < UNFREQS1; j++)
		{
			for(k = 0; k < UNFREQS2; k++)
			{
				if(URES[i][j][k].survival > 0.99)
				{
					(count[i])++;
				}
				//sumva[i] += (URES[i][j][k].Va/URES[i][j][k].Vg);
				//sumva[i] += URES[i][j][k].Va;
				sumva[i] += URES[i][j][k].Vg;
			}
		}
	//	printf("%d\t%f\n",count[i],sumva[i]);
	}
	
	for(i = 0; i < UNBESTMODS; i++){ bestmods[i] = -1; max[i] = 0; }
	
	for(i = 0; i < UNMOD; i++)
	{
		if(count[i] >= UGTHRESH)
		{
			for(j = 0; j < UNBESTMODS; j++)
			{
				if(sumva[i] >= max[j])
				{
					bestmods[j] = i;
					max[j] = sumva[i];
					break;
				}
			}
		}
	}
	free(count);
	free(sumva);
}

void print_bestmods(int *bestmods, model *X, int iter, double *max)
{
	int i,j,k;
	char filename2[50];
	sprintf(filename2,"%spats",UFILENAME);
	FILE *out, *outpats;
	out = fopen(UFILENAME,"a");
	outpats = fopen(filename2,"a");
	
	for(i = 0; i < UNBESTMODS; i++)
	{
		if(bestmods[i] > -1)
		{
			printf("%d: %d\t%f\n",iter,X[bestmods[i]].id, max[i]);
			fprintf(outpats,"%d\t%f\n",X[bestmods[i]].id, max[i]);
			for(j = 0; j < 3; j++)
			{
				for(k = 0; k < 3; k++)
				{
					fprintf(outpats,"%0.4f ",X[bestmods[i]].mod[j][k]);
				}
				fprintf(outpats,"\n");
			}
			fprintf(outpats,"\n");
	
			for(j = 0; j < UNFREQS1; j++)
			{
				for(k = 0; k < UNFREQS2; k++)
				{
					// iteration pattern freq1 freq2 survival Va Vg
					fprintf(out,"%d\t%d\t%0.2f\t%0.2f\t%f\t%f\t%f\n",
						iter,X[bestmods[i]].id,
						UFREQS1[j],UFREQS2[k],
						URES[bestmods[i]][j][k].survival,
						URES[bestmods[i]][j][k].Va,
						URES[bestmods[i]][j][k].Vg);
				}
			}
		}
	}
	printf("\n");
	fclose(out);
	fclose(outpats);
}

void printmodels(model *X)
{
	int i,j,k;
	for(i = 0; i < UNMOD; i++)
	{
		printf("model %d:\n",X[i].id);
		for(j = 0; j < 3; j++)
		{
			for(k = 0; k < 3; k++)
			{
				printf("%0.2f ",X[i].mod[j][k]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void mutate(model *X, int *bestmods)
{
	int i,j,k,l,*nextmods, nmods, lastid;
	model *tempmod;
	//tempmod = make_cube_d(UNBESTMODS,3,3);
	tempmod = malloc(sizeof(model) * UNBESTMODS);
	nextmods = malloc(sizeof(int)*UNBESTMODS);
	
	
	lastid = X[UNMOD-1].id;
	nmods = 0;
	for(i = 0; i < UNBESTMODS; i++)
	{
		if(bestmods[i] > -1)
		{
			nextmods[nmods++] = bestmods[i];
			for(j = 0; j < 3; j++)
			{
				for(k = 0; k < 3; k++)
				{
					tempmod[i].mod[j][k] = X[bestmods[i]].mod[j][k];
				}
			}
		}
	}

	l = 0;
	for(i = 0; i < nmods; i++)
	{
		for(j = 0; j < 3; j++)
		{
			for(k = 0; k < 3; k++)
			{
				X[l].mod[j][k] = tempmod[i].mod[j][k];
			}
		}
		X[l].id = X[nextmods[i]].id;
		convert_modg(X[l].mod,X[l].modg);
		l++;
		for(j = 1; j < UNEXTMODS[i]; j++)
		{
			mutate_pattern(tempmod[i].mod,X[l].mod,X[l].modg);
			X[l].id = ++lastid;
			l++;
		}
	}
	for(i = l; i < UNMOD; i++)
	{
		init_pattern(X[i].mod, X[i].modg);
		X[i].id = ++lastid;
	}
	
	free(tempmod);
	free(nextmods);
	
}

void perform_ga(int maxruns)
{
	int i,j,k,*bestmods;
	double *max;
	model *X;
	X = malloc(sizeof(model) * UNMOD);
	bestmods = malloc(sizeof(int) * UNBESTMODS);
	max = malloc(sizeof(double) * UNBESTMODS);

	for(i = 0; i < UNMOD; i++)
	{
		init_pattern(X[i].mod, X[i].modg);
		X[i].id = i;
	}

	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			X[0].mod[i][j] = UINIT[i][j] >= 0 ? UINIT[i][j] : X[0].mod[i][j];
	convert_modg(X[0].mod,X[0].modg);
	for(i = 0; i < maxruns; i++)
	{
		iteration(X);
		optimise_ga(bestmods,max);
		k = 0;
		for(j = 0; j < UNBESTMODS; j++) k += bestmods[j];
		if(k > (UNBESTMODS*-1))
		{
			print_bestmods(bestmods,X,i,max);
		} else {
			printf("%d: scratch\n\n",i);
		}
		mutate(X,bestmods);
	}
	
}

int main(int argc, char **argv)
{
	if(argc != 2)
	{
		printf("\n%s input.dat\n\n",argv[0]);
		exit(1);
	}
	inputsummary(argv[1]);

	perform_ga(UMAXRUNS);
	
	return 0;
}


