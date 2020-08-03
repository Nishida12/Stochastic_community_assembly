// N_Species_M_Resources_Figure2A.cpp
// 20-Sep-2019
// All resources have some microbial species specialized for each resource.
// Each microorganism takes up resources and secretes metabolic by-products.

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <random>
#include <sstream>
#include <string>
#include <ctime>

using namespace std;

#define TIMESTEP 40000000	// Simulation time
#define dt 0.001			// Clock tick
#define Th 0.00001			// Threshold of biomass
#define SYM_NUM 10			// Number of essential-symbionts for each microorganism
#define SPENUM_MAX 500		// Initial species
#define RESNUM_MAX 10		// Initial number of resources
#define TRI_PARAMS 100		// Number of repetitions by re-choosing parameters
#define TRI 10				// Number of repetitions by changing the initial biomasses
#define FileHead "Result/"	// Storage location

void DistSet(double C[SPENUM_MAX][RESNUM_MAX], double D[SPENUM_MAX][RESNUM_MAX][RESNUM_MAX]);
void SymTarget(int Sym[SPENUM_MAX][SPENUM_MAX]);
void VarReset(double spe[SPENUM_MAX], double res[RESNUM_MAX], double K);
void Calc(double M, double K, double C[SPENUM_MAX][RESNUM_MAX], double D[SPENUM_MAX][RESNUM_MAX][RESNUM_MAX], int Sym[SPENUM_MAX][SPENUM_MAX], double spe[SPENUM_MAX], double res[RESNUM_MAX], int n, int cnt);
void PrintParam(double M, double K, double C[SPENUM_MAX][RESNUM_MAX], double D[SPENUM_MAX][RESNUM_MAX][RESNUM_MAX], int Sym[SPENUM_MAX][SPENUM_MAX], const char* FileName);
mt19937 mt((int)time(0));
uniform_real_distribution<double> rnd1(0.0, 1.0);
normal_distribution<double> rnd4(0.4, 0.01);
uniform_real_distribution<double> spe_rnd(0.01, 1.0);

int main()
{
	double M;							// Death rate
	double K;							// Resource supply
	double C[SPENUM_MAX][RESNUM_MAX];	// Species-specific resource consuming rate
	double D[SPENUM_MAX][RESNUM_MAX][RESNUM_MAX];	// Stoichiometric matrix
	double spe[SPENUM_MAX];				// Biomass
	double res[RESNUM_MAX];				// Resource abundance
	int spenum, resnum;					// Numbers of microbial species and resources
	double t;							// Time
	int Sym[SPENUM_MAX][SPENUM_MAX];	// If species j is essential for species i, Sym[i][j]=1
	int cnt, i, n;

	ofstream f_sum_bio, f_sum_res;
	stringstream FileName;
	
	// Parameter setting
	M = 0.05;
	K = 100;
	
	// Repeating by re-choosing parameters
	for(n=0; n<TRI_PARAMS; n++){
		FileName << FileHead << "summary_biomass_" << n << ".txt";
		f_sum_bio.open(FileName.str());
		FileName.str("");						// Buffer clear
		FileName.clear(stringstream::goodbit);	// Stream state clear
		f_sum_bio << "N";
		for(i=0; i<SPENUM_MAX; i++){
			f_sum_bio << "\t" << "N" << i+1;
		}
		f_sum_bio << endl;

		FileName << FileHead << "summary_resource_" << n << ".txt";
		f_sum_res.open(FileName.str());
		FileName.str("");						// Buffer clear
		FileName.clear(stringstream::goodbit);	// Stream state clear
		f_sum_res << "R";
		for(i=0; i<RESNUM_MAX; i++){
			f_sum_res << "\t" << "R" << i+1;
		}
		f_sum_res << endl;
		
		// Parameter setting
		DistSet(C, D);
		SymTarget(Sym);

		// Repeating by changing the initial biomasses
		for(cnt=0; cnt<TRI; cnt++){
			// Initialization
			VarReset(spe, res, K);

			// Simulation
			Calc(M, K, C, D, Sym, spe, res, n, cnt);

			// Record
			f_sum_bio << (cnt+1);
			for(i=0; i<SPENUM_MAX; i++){
				f_sum_bio << "\t" << spe[i];
			}
			f_sum_bio << endl;
			f_sum_res << (cnt+1);
			for(i=0; i<RESNUM_MAX; i++){
				f_sum_res << "\t" << res[i];
			}
			f_sum_res << endl;
		}

		FileName << FileHead << "TimeLapse/" << n << "/params.txt";
		PrintParam(M, K, C, D, Sym, FileName.str().c_str());
		FileName.str("");						// buffer clear
		FileName.clear(stringstream::goodbit);	// stream state clear

		f_sum_bio.close();
		f_sum_res.close();
	}

	return 0;
}

// Parameter setting of C (species-specific resource consuming rate) and D (Stoichiometric matrix)
void DistSet(double C[SPENUM_MAX][RESNUM_MAX], double D[SPENUM_MAX][RESNUM_MAX][RESNUM_MAX])
{
	int i, j, k;
	double sum;
	double tmp[RESNUM_MAX];

	for(i=0; i<SPENUM_MAX; i++){
		k = i / (SPENUM_MAX / RESNUM_MAX);
		C[i][k] = rnd4(mt);
		sum = 0.0;
		for(j=0; j<RESNUM_MAX; j++){
			if(j != k){
				tmp[j] = rnd1(mt);
				sum += tmp[j];
			}
		}
		for(j=0; j<RESNUM_MAX; j++){
			if(j != k){
				C[i][j] = (1.0-C[i][k])*tmp[j]/sum;
			}
		}
	}

	for(k=0; k<SPENUM_MAX; k++){
		for(i=0; i<RESNUM_MAX; i++){
			for(j=0; j<RESNUM_MAX; j++){
				D[k][i][j] = rnd1(mt)/RESNUM_MAX;
			}
		}
	}
}

// Parameter setting of Sym (If species j is essential for species i, Sym[i][j]=1)
void SymTarget(int Sym[SPENUM_MAX][SPENUM_MAX]){
	int i, j, k;
	int target[SYM_NUM], tmp;

	uniform_int_distribution<> rnd_int(0,SPENUM_MAX-1);

	for(i=0; i<SPENUM_MAX; i++){
		for(j=0; j<SYM_NUM; j++){
			k = i / (SPENUM_MAX / RESNUM_MAX);
			target[j] = rnd_int(mt);
			tmp = target[j]  / (SPENUM_MAX / RESNUM_MAX);
			while(k == tmp){
				target[j] = rnd_int(mt);
				tmp = target[j] / (SPENUM_MAX / RESNUM_MAX);
			}
		}
		for(j=0; j<SPENUM_MAX; j++){
			Sym[i][j] = 0;
		}
		for(j=0; j<SYM_NUM; j++){
			Sym[i][target[j]] = 1;
		}
	}
}

// Initialization of biomasses and resource abundances
void VarReset(double spe[SPENUM_MAX], double res[RESNUM_MAX], double K)
{
	int i;

	for(i=0; i<SPENUM_MAX; i++){
		spe[i] = spe_rnd(mt);
	}

	for(i=0; i<RESNUM_MAX; i++){
		res[i] = 0.0;
	}
	res[0] = K;
}

// Simulation
void Calc(double M, double K, double C[SPENUM_MAX][RESNUM_MAX], double D[SPENUM_MAX][RESNUM_MAX][RESNUM_MAX], int Sym[SPENUM_MAX][SPENUM_MAX], double spe[SPENUM_MAX], double res[RESNUM_MAX], int n, int cnt)
{
	int RecCnt;						// Record frequency
	double copy_spe[SPENUM_MAX];
	double copy_res[RESNUM_MAX];
	int i, j, k, t;
	double tmp1, tmp2;
	double max_spe[SPENUM_MAX];		// Total amount of essential microorganisms
	double sym_K = 0.01;

	ofstream f_bio, f_res;
	stringstream FileName;

	// Initialization
	RecCnt = 5000;

	FileName << FileHead << "TimeLapse/" << n << "/t_biomass_" << cnt+1 << ".txt";
	f_bio.open(FileName.str());
	FileName.str("");	// Buffer clear
	FileName.clear(stringstream::goodbit);	// Stream state clear
	f_bio << "t";
	for(int i=0; i<SPENUM_MAX; i++){
		f_bio << "\t" << "N" << i+1;
	}
	f_bio << endl;

	FileName << FileHead << "TimeLapse/" << n << "/t_resource_" << cnt+1 << ".txt";
	f_res.open(FileName.str());
	FileName.str("");	// Buffer clear
	FileName.clear(stringstream::goodbit);	// Stream state clear
	f_res << "t";
	for(int i=0; i<RESNUM_MAX; i++){
		f_res << "\t" << "R" << i+1;
	}
	f_res << endl;

	// Timelapse
	for(t=0; t<TIMESTEP; t++){
		// Record
		if(t%RecCnt == 0){
			f_bio << t;
			for(i=0; i<SPENUM_MAX; i++){
				f_bio << '\t' << spe[i];
			}
			f_bio << endl;

			f_res << t;
			for(i=0; i<RESNUM_MAX; i++){
				f_res << '\t' << res[i];
			}
			f_res << endl;
		}
		
		// Update
		for(i=0; i<SPENUM_MAX; i++){
			if(spe[i] > Th){
				max_spe[i] = 0;
				for(j=0; j<SPENUM_MAX; j++){
					if(Sym[i][j] == 1){
						max_spe[i] += spe[j];
					}
				}
			}
		}
		
		for(i=0; i<SPENUM_MAX; i++){
			if(spe[i] > Th){
				tmp1 = 0.0;
				for(j=0; j<RESNUM_MAX; j++){
					tmp1 += C[i][j]*res[j];
				}
				copy_spe[i] = spe[i]*(max_spe[i]/(sym_K+max_spe[i]))*tmp1 - spe[i]*M;
			}
		}
				
		for(i=0; i<RESNUM_MAX; i++){
			tmp1 = 0.0;
			for(j=0; j<SPENUM_MAX; j++){
				tmp1 += C[j][i]*res[i]*spe[j]*(max_spe[j]/(sym_K+max_spe[j]));
			}
			
			tmp2 = 0.0;
			for(j=0; j<SPENUM_MAX; j++){
				for(k=0; k<RESNUM_MAX; k++){
					tmp2 += D[j][i][k]*C[j][k]*res[k]*spe[j]*(max_spe[j]/(sym_K+max_spe[j]));
				}
			}

			if(i == 0){
				copy_res[i] = K - res[i] - tmp1 + tmp2;
			}else{
				copy_res[i] = tmp2 - tmp1;
			}
		}

		for(i=0; i<SPENUM_MAX; i++){
			if(spe[i] > Th){
				spe[i] += copy_spe[i]*dt;
				if(spe[i] <= Th) spe[i] = 0.0;
			}
		}
		for(i=0; i<RESNUM_MAX; i++){
			res[i] += copy_res[i]*dt;
			if(res[i] < 0.0) res[i] = 0.0;
		}
	}

	// Record
	f_bio << t;
	for(i=0; i<SPENUM_MAX; i++){
		f_bio << '\t' << spe[i];
	}
	f_bio << endl;
	f_res << t;
	for(i=0; i<RESNUM_MAX; i++){
		f_res << '\t' << res[i];
	}
	f_res << endl;

	f_bio.close();
	f_res.close();
}

// Parameters record
void PrintParam(double M, double K, double C[SPENUM_MAX][RESNUM_MAX], double D[SPENUM_MAX][RESNUM_MAX][RESNUM_MAX], int Sym[SPENUM_MAX][SPENUM_MAX], const char* FileName)
{
	int i, j, k;
	ofstream f_param;

	f_param.open(FileName);

	f_param << "M: " << M << endl;
	f_param << "K: " << K << endl;
	f_param << "TIMESTEP: " << TIMESTEP << endl;
	f_param << "dt: " << dt << endl;
	f_param << "SPENUM_MAX: " << SPENUM_MAX << endl;
	f_param << "RESNUM_MAX: " << RESNUM_MAX << endl;
	f_param << "TRI_PARAMS: " << TRI_PARAMS << endl;
	f_param << "TRI: " << TRI << endl;

	f_param << "SymbiosisTarget" << endl;
	for(i=0; i<SPENUM_MAX; i++){
		for(j=0; j<(SPENUM_MAX-1); j++){
			f_param << Sym[i][j] << '\t';
		}
		f_param << Sym[i][j] << endl;
	}
	f_param << endl;

	f_param << "C" << endl;
	for(i=0; i<SPENUM_MAX; i++){
		for(j=0; j<(RESNUM_MAX-1); j++){
			f_param << setprecision(4) << C[i][j] << '\t';
		}
		f_param << setprecision(4) << C[i][j] << endl;
	}
	f_param << endl;

	for(k=0; k<SPENUM_MAX; k++){
		f_param << "D" << k+1 << endl;
		for(i=0; i<RESNUM_MAX; i++){
			for(j=0; j<(RESNUM_MAX-1); j++){
				f_param << setprecision(4) << D[k][i][j] << '\t';
			}
			f_param << setprecision(4) << D[k][i][j] << endl;
		}
	}

	f_param.close();
}