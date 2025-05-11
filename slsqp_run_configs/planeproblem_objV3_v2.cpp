// planeproblem.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<set>
#include<cmath>
#include<nlopt.hpp>
#include<omp.h>
#include <bits/stdc++.h>
#include <ctime>

#include <algorithm>
#include <numeric>

#include <cstdlib>
#include <cfloat>
#include <cfenv>

using namespace std;

#pragma STDC FENV_ACCESS ON
void show_fe_exceptions(void)
{
/*
    if(fetestexcept(FE_DIVBYZERO)) cout << "exceptions raised: FE_DIVBYZERO" << endl;
    if(fetestexcept(FE_INEXACT))   cout << "exceptions raised: FE_INEXACT" << endl;
    if(fetestexcept(FE_INVALID))   cout << "exceptions raised: FE_INVALID" << endl;
    if(fetestexcept(FE_OVERFLOW))  cout << "exceptions raised: FE_OVERFLOW" << endl;
    if(fetestexcept(FE_UNDERFLOW)) cout << "exceptions raised: FE_UNDERFLOW" << endl;
    feclearexcept(FE_ALL_EXCEPT);
*/
}

double obj_exec_time = 0;
int order;
double TOL = 1e-10;
vector<double> weights;
vector<double> values;
struct BB {
	double x = 0;
	double y = 0;
	double detJ;
	vector<vector<double> > Jmat;
	vector<vector<double> > Bd;
	vector<vector<double> > Bsmat;
};

struct RetBC {
	vector<vector<double>> Kmat;
	vector<vector<double>> Fmat;
};

struct node {
	int id;
	double x_cord;
	double y_cord;
};

typedef struct boundry_condition {
	int node1, node2, ele, dc_x1, dc_x2, dc_y1, dc_y2;
	double tractvec[4]; 
	double dispvec[4];
}bc;

struct element {
	int node1, node2, node3, node4, ele, mattype, nummatpar;
	vector<double> matpar;
};

vector<bc> bc_arr;
vector<node> nodes;
vector<element> elements;

struct constdata {
	vector<vector<double>> Amat, intmat;
	vector<double> Fvec;
	int rows;
};
struct objdata {
	vector<element> elemdat;
	vector<node> nodecord;
	double dispscalefac;
	int numnode;
	int numele;
	int lnpos;
};

constdata constraint_data;
objdata obj_data;

struct Errret {
	double errv;
	vector<double> grdval;
};

struct Cnstrelret {
	vector<double> delval;
	vector<vector<double>> dfdstrain;
	vector<vector<double>> dfdstress;
};


double myfunc(unsigned, const double*, double*, void*);
void myconstraint(unsigned, double*, unsigned, const double*, double*, void* );
// double myconstraint(unsigned, const double*, double*, void* );
struct GQVandW {
	vector<double> value;
	vector<double> weight;
};

RetBC boundry_condition_fun(bc, GQVandW, int);
GQVandW getWeightAndValues(int order)
{
	GQVandW result;
	vector<double> weights;
	vector<double> values;
	if (order == 2) {
		weights.push_back(1.0000000000000000);
		weights.push_back(1.0000000000000000);
		values.push_back(-0.5773502691896257);
		values.push_back(0.5773502691896257);
	}

	if (order == 3) {
		weights.push_back(0.8888888888888888);
		weights.push_back(0.5555555555555556);
		weights.push_back(0.5555555555555556);
		values.push_back(0.0000000000000000);
		values.push_back(-0.7745966692414834);
		values.push_back(0.7745966692414834);
	}

	if (order == 4) {
		weights.push_back(0.6521451548625461);
		weights.push_back(0.6521451548625461);
		weights.push_back(0.3478548451374538);
		weights.push_back(0.3478548451374538);

		values.push_back(-0.3399810435848563);
		values.push_back(0.3399810435848563);
		values.push_back(-0.8611363115940526);
		values.push_back(0.8611363115940526);

	}

	if (order == 5) {
		weights.push_back(0.5688888888888889);
		weights.push_back(0.4786286704993665);
		weights.push_back(0.4786286704993665);
		weights.push_back(0.2369268850561891);
		weights.push_back(0.2369268850561891);

		values.push_back(0.0000000000000000);
		values.push_back(-0.5384693101056831);
		values.push_back(0.5384693101056831);
		values.push_back(-0.9061798459386640);
		values.push_back(0.9061798459386640);

	}

	if (order == 6) {
		weights.push_back(0.3607615730481386);
		weights.push_back(0.3607615730481386);
		weights.push_back(0.4679139345726910);
		weights.push_back(0.4679139345726910);
		weights.push_back(0.1713244923791704);
		weights.push_back(0.1713244923791704);

		values.push_back(-0.6612093864662645);
		values.push_back(0.6612093864662645);
		values.push_back(-0.2386191860831969);
		values.push_back(0.2386191860831969);
		values.push_back(-0.9324695142031521);
		values.push_back(0.9324695142031521);
	}

	if (order == 7) {
		weights.push_back(0.4179591836734694);
		weights.push_back(0.3818300505051189);
		weights.push_back(0.3818300505051189);
		weights.push_back(0.2797053914892766);
		weights.push_back(0.2797053914892766);
		weights.push_back(0.1294849661688697);
		weights.push_back(0.1294849661688697);

		values.push_back(0.0000000000000000);
		values.push_back(-0.4058451513773972);
		values.push_back(0.4058451513773972);
		values.push_back(-0.7415311855993945);
		values.push_back(0.7415311855993945);
		values.push_back(-0.9491079123427585);
		values.push_back(0.9491079123427585);

	}

	if (order == 8) {
		weights.push_back(0.3626837833783620);
		weights.push_back(0.3626837833783620);
		weights.push_back(0.3137066458778873);
		weights.push_back(0.3137066458778873);
		weights.push_back(0.2223810344533745);
		weights.push_back(0.2223810344533745);
		weights.push_back(0.1012285362903763);
		weights.push_back(0.1012285362903763);

		values.push_back(-0.1834346424956498);
		values.push_back(0.1834346424956498);
		values.push_back(-0.5255324099163290);
		values.push_back(0.5255324099163290);
		values.push_back(-0.7966664774136267);
		values.push_back(0.7966664774136267);
		values.push_back(-0.9602898564975363);
		values.push_back(0.9602898564975363);
	}

	if (order == 9) {
		weights.push_back(0.3302393550012598);
		weights.push_back(0.1806481606948574);
		weights.push_back(0.1806481606948574);
		weights.push_back(0.0812743883615744);
		weights.push_back(0.0812743883615744);
		weights.push_back(0.3123470770400029);
		weights.push_back(0.3123470770400029);
		weights.push_back(0.2606106964029354);
		weights.push_back(0.2606106964029354);

		values.push_back(0.0000000000000000);
		values.push_back(-0.8360311073266358);
		values.push_back(0.8360311073266358);
		values.push_back(-0.9681602395076261);
		values.push_back(0.9681602395076261);
		values.push_back(-0.3242534234038089);
		values.push_back(0.3242534234038089);
		values.push_back(-0.6133714327005904);
		values.push_back(0.6133714327005904);
	}

	if (order == 10) {
		weights.push_back(0.2955242247147529);
		weights.push_back(0.2955242247147529);
		weights.push_back(0.2692667193099963);
		weights.push_back(0.2692667193099963);
		weights.push_back(0.2190863625159820);
		weights.push_back(0.2190863625159820);
		weights.push_back(0.1494513491505806);
		weights.push_back(0.1494513491505806);
		weights.push_back(0.0666713443086881);
		weights.push_back(0.0666713443086881);

		values.push_back(-0.1488743389816312);
		values.push_back(0.1488743389816312);
		values.push_back(-0.4333953941292472);
		values.push_back(0.4333953941292472);
		values.push_back(-0.6794095682990244);
		values.push_back(0.6794095682990244);
		values.push_back(-0.8650633666889845);
		values.push_back(0.8650633666889845);
		values.push_back(-0.9739065285171717);
		values.push_back(0.9739065285171717);
	}
	result.weight = weights;
	result.value = values;
	return result;
}

BB Bsigma_Bd(double *xcoord, double  *ycoord, double s, double t) {

	vector<double>Nis(4);
	Nis[0] = (t - 1) / 4;
	Nis[1] = (1 - t) / 4;
	Nis[2] = (t + 1) / 4;
	Nis[3] = -(t + 1) / 4;

	vector<double>Nit(4);
	Nit[0] = (s - 1) / 4;
	Nit[1] = -(s + 1) / 4;
	Nit[2] = (s + 1) / 4;
	Nit[3] = (1 - s) / 4;

	vector<double>N(4);
	N[0] = (1 - s) * (1 - t) / 4;
	N[1] = (1 + s) * (1 - t) / 4;
	N[2] = (1 + s) * (1 + t) / 4;
	N[3] = (1 - s) * (1 + t) / 4;

	double matX = 0, matY = 0;
	for (int i = 0; i < 4; i++) {
		matX += N[i] * xcoord[i];
		matY += N[i] * ycoord[i];
	}

	vector<vector<double> > Jmat;
	vector<double> Jmatrow(2);
	for (int i = 0; i < 2; i++) {
		Jmat.push_back(Jmatrow);
		for (int j = 0; j < 2; j++)
			Jmat[i][j] = 0;
	}

	for (int i = 0; i < 4; i++) {
		Jmat[0][0] += Nis[i] * xcoord[i];
		Jmat[0][1] += Nit[i] * xcoord[i];
		Jmat[1][0] += Nis[i] * ycoord[i];
		Jmat[1][1] += Nit[i] * ycoord[i];
	}

	vector<vector<double> > invJ(2, vector<double>(2, 0));
	double detJ = (Jmat[0][0] * Jmat[1][1]) - (Jmat[0][1] * Jmat[1][0]);
	invJ[0][0] = Jmat[1][1] / detJ;
	invJ[0][1] = -(Jmat[0][1]) / detJ;
	invJ[1][0] = -(Jmat[1][0]) / detJ;
	invJ[1][1] = Jmat[0][0] / detJ;

	vector<vector<double> > Bd;
	vector<double>Bdrow(8);
	vector<vector<double> > Bs;
	vector<double> Bsrow(24);
	for (int i = 0; i < 3; i++) {
		Bd.push_back(Bdrow);
		Bs.push_back(Bsrow);
		for (int j = 0; j < 8; j++)
			Bd[i][j] = 0;
		for (int j = 0; j < 24; j++)
			Bs[i][j] = 0;
	}

	vector<vector<double> > Tmat;
	vector<double>Tmatrow(24);
	for (int i = 0; i < 24; i++) {
		Tmat.push_back(Tmatrow);
		for (int j = 0; j < 24; j++)
			Tmat[i][j] = 0;
	}
	vector<vector<double> > T1mat;
	vector<double>T1matrow(6);
	for (int i = 0; i < 6; i++) {
		T1mat.push_back(T1matrow);
		for (int j = 0; j < 6; j++)
			T1mat[i][j] = 0;
	}

	Bd[0][0] = ((Nis[0] * invJ[0][0]) + (Nit[0] * invJ[1][0]));
	Bd[0][2] = ((Nis[1] * invJ[0][0]) + (Nit[1] * invJ[1][0]));
	Bd[0][4] = ((Nis[2] * invJ[0][0]) + (Nit[2] * invJ[1][0]));
	Bd[0][6] = ((Nis[3] * invJ[0][0]) + (Nit[3] * invJ[1][0]));

	Bd[1][1] = ((Nis[0] * invJ[0][1]) + (Nit[0] * invJ[1][1]));
	Bd[1][3] = ((Nis[1] * invJ[0][1]) + (Nit[1] * invJ[1][1]));
	Bd[1][5] = ((Nis[2] * invJ[0][1]) + (Nit[2] * invJ[1][1]));
	Bd[1][7] = ((Nis[3] * invJ[0][1]) + (Nit[3] * invJ[1][1]));

	Bd[2][0] = ((Nis[0] * invJ[0][1]) + (Nit[0] * invJ[1][1]));
	Bd[2][2] = ((Nis[1] * invJ[0][1]) + (Nit[1] * invJ[1][1]));
	Bd[2][4] = ((Nis[2] * invJ[0][1]) + (Nit[2] * invJ[1][1]));
	Bd[2][6] = ((Nis[3] * invJ[0][1]) + (Nit[3] * invJ[1][1]));

	Bd[2][1] = ((Nis[0] * invJ[0][0]) + (Nit[0] * invJ[1][0]));
	Bd[2][3] = ((Nis[1] * invJ[0][0]) + (Nit[1] * invJ[1][0]));
	Bd[2][5] = ((Nis[2] * invJ[0][0]) + (Nit[2] * invJ[1][0]));
	Bd[2][7] = ((Nis[3] * invJ[0][0]) + (Nit[3] * invJ[1][0]));

	vector<vector<double> > Dmat(3, vector<double>(3, 0));
	Dmat[0][0] = (invJ[1][1] * invJ[1][1]);
	Dmat[0][1] = (invJ[0][1] * invJ[0][1]);
	Dmat[0][2] = -2 * (invJ[1][1] * invJ[0][1]);
	Dmat[1][0] = (invJ[1][0] * invJ[1][0]);
	Dmat[1][1] = (invJ[0][0] * invJ[0][0]);
	Dmat[1][2] = -2 * (invJ[0][0] * invJ[1][0]);
	Dmat[2][0] = -(invJ[1][0] * invJ[1][1]);
	Dmat[2][1] = -(invJ[0][0] * invJ[0][1]);
	Dmat[2][2] = ((invJ[0][0] * invJ[1][1]) + (invJ[0][1] * invJ[1][0]));

	T1mat[0][0] = 1;
	T1mat[1][1] = Jmat[0][0];
	T1mat[1][2] = Jmat[1][0];
	T1mat[2][1] = Jmat[0][1];
	T1mat[2][2] = Jmat[1][1];
	T1mat[3][3] = (Jmat[0][0] * Jmat[0][0]);
	T1mat[3][4] = (Jmat[1][0] * Jmat[1][0]);
	T1mat[3][5] = 2 * (Jmat[0][0] * Jmat[1][0]);
	T1mat[4][3] = (Jmat[0][1] * Jmat[0][1]);
	T1mat[4][4] = (Jmat[1][1] * Jmat[1][1]);
	T1mat[4][5] = 2 * (Jmat[0][1] * Jmat[1][1]);
	T1mat[5][3] = (Jmat[0][0] * Jmat[0][1]);
	T1mat[5][4] = (Jmat[1][0] * Jmat[1][1]);
	T1mat[5][5] = ((Jmat[0][0] * Jmat[1][1]) + (Jmat[0][1] * Jmat[1][0]));

	Tmat[0][0] = Tmat[6][6] = Tmat[12][12] = Tmat[18][18] = 1;
	Tmat[1][1] = Tmat[7][7] = Tmat[13][13] = Tmat[19][19] = T1mat[1][1];
	Tmat[1][2] = Tmat[7][8] = Tmat[13][14] = Tmat[19][20] = T1mat[1][2];
	Tmat[2][1] = Tmat[8][7] = Tmat[14][13] = Tmat[20][19] = T1mat[2][1];
	Tmat[2][2] = Tmat[8][8] = Tmat[14][14] = Tmat[20][20] = T1mat[2][2];
	Tmat[3][3] = Tmat[9][9] = Tmat[15][15] = Tmat[21][21] = T1mat[3][3];
	Tmat[3][4] = Tmat[9][10] = Tmat[15][16] = Tmat[21][22] = T1mat[3][4];
	Tmat[3][5] = Tmat[9][11] = Tmat[15][17] = Tmat[21][23] = T1mat[3][5];
	Tmat[4][3] = Tmat[10][9] = Tmat[16][15] = Tmat[22][21] = T1mat[4][3];
	Tmat[4][4] = Tmat[10][10] = Tmat[16][16] = Tmat[22][22] = T1mat[4][4];
	Tmat[4][5] = Tmat[10][11] = Tmat[16][17] = Tmat[22][23] = T1mat[4][5];
	Tmat[5][3] = Tmat[11][9] = Tmat[17][15] = Tmat[23][21] = T1mat[5][3];
	Tmat[5][4] = Tmat[11][10] = Tmat[17][16] = Tmat[23][22] = T1mat[5][4];
	Tmat[5][5] = Tmat[11][11] = Tmat[17][17] = Tmat[23][23] = T1mat[5][5];

	Bs[0][0] = (3 * t * (s - 1) * (s * s + s + 5 * t * t - 5)) / 8;
	Bs[0][1] = (3 * t * (s - 1) * (s - 1) * (s + 1)) / 8;
	Bs[0][2] = -((s - 1) * (-3 * s * s * t + s * s - 3 * s * t + s - 15 * t * t * t + 3 * t * t + 15 * t - 3)) / 8;
	Bs[0][4] = ((s - 1) * (t - 1) * (5 * t * t + 2 * t - 1)) / 8;
	Bs[0][5] = ((3 * t - 1) * (s - 1) * (s - 1) * (s + 1)) / 8;
	Bs[0][6] = (3 * t * (s + 1) * (-s * s + s - 5 * t * t + 5)) / 8;
	Bs[0][7] = (3 * t * (s - 1) * (s + 1) * (s + 1)) / 8;
	Bs[0][8] = -((s + 1) * (3 * s * s * t - s * s - 3 * s * t + s + 15 * t * t * t - 3 * t * t - 15 * t + 3)) / 8;
	Bs[0][10] = -((s + 1) * (t - 1) * (5 * t * t + 2 * t - 1)) / 8;
	Bs[0][11] = ((3 * t - 1) * (s - 1) * (s + 1) * (s + 1)) / 8;
	Bs[0][12] = -(3 * t * (s + 1) * (-s * s + s - 5 * t * t + 5)) / 8;
	Bs[0][13] = -(3 * t * (s - 1) * (s + 1) * (s + 1)) / 8;
	Bs[0][14] = ((s + 1) * (-3 * s * s * t - s * s + 3 * s * t + s - 15 * t * t * t - 3 * t * t + 15 * t + 3)) / 8;
	Bs[0][16] = ((s + 1) * (t + 1) * (5 * t * t - 2 * t - 1)) / 8;
	Bs[0][17] = ((3 * t + 1) * (s - 1) * (s + 1) * (s + 1)) / 8;
	Bs[0][18] = -(3 * t * (s - 1) * (s * s + s + 5 * t * t - 5)) / 8;
	Bs[0][19] = -(3 * t * (s - 1) * (s - 1) * (s + 1)) / 8;
	Bs[0][20] = ((s - 1) * (3 * s * s * t + s * s + 3 * s * t + s + 15 * t * t * t + 3 * t * t - 15 * t - 3)) / 8;
	Bs[0][22] = -((s - 1) * (t + 1) * (5 * t * t - 2 * t - 1)) / 8;
	Bs[0][23] = ((3 * t + 1) * (s - 1) * (s - 1) * (s + 1)) / 8;

	Bs[1][0] = (3 * s * (t - 1) * (5 * s * s + t * t + t - 5)) / 8;
	Bs[1][1] = -((t - 1) * (-15 * s * s * s + 3 * s * s - 3 * s * t * t - 3 * s * t + 15 * s + t * t + t - 3)) / 8;
	Bs[1][2] = (3 * s * (t - 1) * (t - 1) * (t + 1)) / 8;
	Bs[1][3] = ((s - 1) * (t - 1) * (5 * s * s + 2 * s - 1)) / 8;
	Bs[1][5] = ((3 * s - 1) * (t - 1) * (t - 1) * (t + 1)) / 8;
	Bs[1][6] = -(3 * s * (t - 1) * (5 * s * s + t * t + t - 5)) / 8;
	Bs[1][7] = ((t - 1) * (15 * s * s * s + 3 * s * s + 3 * s * t * t + 3 * s * t - 15 * s + t * t + t - 3)) / 8;
	Bs[1][8] = -(3 * s * (t - 1) * (t - 1) * (t + 1)) / 8;
	Bs[1][9] = ((s + 1) * (t - 1) * (-5 * s * s + 2 * s + 1)) / 8;
	Bs[1][11] = ((3 * s + 1) * (t - 1) * (t - 1) * (t + 1)) / 8;
	Bs[1][12] = -(3 * s * (t + 1) * (-5 * s * s - t * t + t + 5)) / 8;
	Bs[1][13] = ((t + 1) * (-15 * s * s * s - 3 * s * s - 3 * s * t * t + 3 * s * t + 15 * s - t * t + t + 3)) / 8;
	Bs[1][14] = -(3 * s * (t - 1) * (t + 1) * (t + 1)) / 8;
	Bs[1][15] = -((s + 1) * (t + 1) * (-5 * s * s + 2 * s + 1)) / 8;
	Bs[1][17] = ((3 * s + 1) * (t - 1) * (t + 1) * (t + 1)) / 8;
	Bs[1][18] = (3 * s * (t + 1) * (-5 * s * s - t * t + t + 5)) / 8;
	Bs[1][19] = ((t + 1) * (-15 * s * s * s + 3 * s * s - 3 * s * t * t + 3 * s * t + 15 * s + t * t - t - 3)) / 8;
	Bs[1][20] = (3 * s * (t - 1) * (t + 1) * (t + 1)) / 8;
	Bs[1][21] = -((s - 1) * (t + 1) * (5 * s * s + 2 * s - 1)) / 8;
	Bs[1][23] = ((3 * s - 1) * (t - 1) * (t + 1) * (t + 1)) / 8;

	Bs[2][0] = -(15 * s * s * s * s + 18 * s * s * t * t - 36 * s * s + 15 * t * t * t * t - 36 * t * t + 24) / 32;
	Bs[2][1] = -((3 * s + 1) * (s - 1) * (5 * s * s + 2 * s + 6 * t * t - 9)) / 32;
	Bs[2][2] = -((3 * t + 1) * (t - 1) * (6 * s * s + 5 * t * t + 2 * t - 9)) / 32;
	Bs[2][3] = -((5 * s + 1) * (s - 1) * (s - 1) * (s + 1)) / 32;
	Bs[2][4] = -((5 * t + 1) * (t - 1) * (t - 1) * (t + 1)) / 32;
	Bs[2][5] = -((3 * s + 1) * (3 * t + 1) * (s - 1) * (t - 1)) / 16;
	Bs[2][6] = (15 * s * s * s * s + 18 * s * s * t * t - 36 * s * s + 15 * t * t * t * t - 36 * t * t + 24) / 32;
	Bs[2][7] = ((3 * s - 1) * (s + 1) * (-5 * s * s + 2 * s - 6 * t * t + 9)) / 32;
	Bs[2][8] = ((3 * t + 1) * (t - 1) * (6 * s * s + 5 * t * t + 2 * t - 9)) / 32;
	Bs[2][9] = ((5 * s - 1) * (s - 1) * (s + 1) * (s + 1)) / 32;
	Bs[2][10] = ((5 * t + 1) * (t - 1) * (t - 1) * (t + 1)) / 32;
	Bs[2][11] = -((3 * s - 1) * (3 * t + 1) * (s + 1) * (t - 1)) / 16;
	Bs[2][12] = -(15 * s * s * s * s + 18 * s * s * t * t - 36 * s * s + 15 * t * t * t * t - 36 * t * t + 24) / 32;
	Bs[2][13] = -((3 * s - 1) * (s + 1) * (-5 * s * s + 2 * s - 6 * t * t + 9)) / 32;
	Bs[2][14] = -((3 * t - 1) * (t + 1) * (-6 * s * s - 5 * t * t + 2 * t + 9)) / 32;
	Bs[2][15] = -((5 * s - 1) * (s - 1) * (s + 1) * (s + 1)) / 32;
	Bs[2][16] = -((5 * t - 1) * (t - 1) * (t + 1) * (t + 1)) / 32;
	Bs[2][17] = -((3 * s - 1) * (3 * t - 1) * (s + 1) * (t + 1)) / 16;
	Bs[2][18] = (15 * s * s * s * s + 18 * s * s * t * t - 36 * s * s + 15 * t * t * t * t - 36 * t * t + 24) / 32;
	Bs[2][19] = ((3 * s + 1) * (s - 1) * (5 * s * s + 2 * s + 6 * t * t - 9)) / 32;
	Bs[2][20] = ((3 * t - 1) * (t + 1) * (-6 * s * s - 5 * t * t + 2 * t + 9)) / 32;
	Bs[2][21] = ((5 * s + 1) * (s - 1) * (s - 1) * (s + 1)) / 32;
	Bs[2][22] = ((5 * t - 1) * (t - 1) * (t + 1) * (t + 1)) / 32;
	Bs[2][23] = -((3 * s + 1) * (3 * t - 1) * (s - 1) * (t + 1)) / 16;

	vector<vector<double> > Bs1mat(3, vector<double>(24, 0));  // Bs1mat = Dmat*Bs
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 24; j++) {
			Bs1mat[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				Bs1mat[i][j] += Dmat[i][k] * Bs[k][j];
			}
		}
	}


	vector<vector<double> > Bsmat(3, vector<double>(24, 0)); // Bsmat = Dmat*Bs*Tmat = Bs1mat*Tmat 
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 24; j++) {
			Bsmat[i][j] = 0;
			for (int k = 0; k < 24; k++) {
				Bsmat[i][j] += Bs1mat[i][k] * Tmat[k][j];
			}
		}
	}

	BB BBresult;
	BBresult.x = matX;
	BBresult.y = matY;
	BBresult.Jmat = Jmat;
	BBresult.detJ = detJ;
	BBresult.Bd = Bd;
	BBresult.Bsmat = Bsmat;
	return BBresult;
}

vector<BB> get_BB_results_mat(int k, int order, double *xcoord, double *ycoord, vector<double> values) {
	vector<BB> BB_results_mat(order);
	
	if (k == 2)
	{
		for (int i = 0; i < order; i++)
		{
			BB_results_mat[i] = Bsigma_Bd(xcoord, ycoord, 1, values[i]);
		}
	}
	else if (k == 4)
	{
		for (int i = 0; i < order; i++)
		{
			BB_results_mat[i] = Bsigma_Bd(xcoord, ycoord, -1, values[i]);
		}
	}
	else if (k == 3)
	{
		for (int i = 0; i < order; i++)
		{
			BB_results_mat[i] = Bsigma_Bd(xcoord, ycoord, values[i], 1);
		}
	}
	else if (k == 1)
	{
		for (int i = 0; i < order; i++)
		{
			BB_results_mat[i] = Bsigma_Bd(xcoord, ycoord, values[i], -1);
		}
	}
	return BB_results_mat;
}

bool cmp_for_desc(int x, int y)
{
    if (x > y)
        return true;
    else
        return false;
}
int nn = 0, nl = 0;
const double EPS = 1E-13;
std::vector<int> rset;
int compute_rank(vector<vector<double>> A) {
    int n = A.size();
    int m = A[0].size();

    int rank = 0;
    vector<bool> row_selected(n, false);
    for (int i = 0; i < m; ++i) {
        int j;
        for (j = 0; j < n; ++j) {
            if (!row_selected[j] && abs(A[j][i]) > EPS)
                break;
        }

        if (j != n) {
            ++rank;
	    rset.push_back(j);
            row_selected[j] = true;
            for (int p = i + 1; p < m; ++p)
                A[j][p] /= A[j][i];
            for (int k = 0; k < n; ++k) {
                if (k != j && abs(A[k][i]) > EPS) {
                    for (int p = i + 1; p < m; ++p)
                        A[k][p] -= A[j][p] * A[k][i];
                }
            }
        }
    }
    return rank;
}
int Amat_rows = 0, Amat_cols = 0;
int main(int argc, char* argv[]) {
	// reading codrinates of nodes
	vector<vector<double>> Amat, intmat;
	int colm;
	int row;
	int count_dc = 0;
	int thereflg;	
	double tempmatpar;
	vector<int> rest_dispdof;
	int dof, numnode;
	numnode = 0;
	ifstream infile_node_cord(argv[1]);
	//ifstream infile_node_cord("../Data_640_dof/NodeData.txt");
	int i;
	double ix, iy;
	while (infile_node_cord >>i >>ix >>iy) {
		//int i;
		//double x, y;
		node n;
		//infile_node_cord >> i >> x >> y;
		n.id = i;
		std::cout<<n.id<<" ";
		n.x_cord = ix;
		n.y_cord = iy;
		nodes.push_back(n);
		numnode++;
	}
	infile_node_cord.close();
	//numnode = numnode-1;
	cout<<"Number of nodes = " << nodes.size() << endl;
	/*for (int a=0;a<numnode;a++)
	{
		cout << "id : " << nodes[a].id << " x value : " << nodes[a].x_cord << " y value : " << nodes[a].y_cord << endl;
	}*/

	//reading node connectivity or elements
	int numele = 0;
	ifstream infile_ele(argv[2]);
	//ifstream infile_ele("../Data_640_dof/NodeConnect.txt ");
	while (!infile_ele.eof()) {
		element tempe;
		vector<double> matpar_ele;
		infile_ele >> tempe.ele >> tempe.node1 >> tempe.node2 >> tempe.node3 >> tempe.node4 >> tempe.mattype >> tempe.nummatpar;
		for (int i=0; i < tempe.nummatpar; i++)
		{
			infile_ele >> tempmatpar;
			matpar_ele.push_back(tempmatpar);
		}
		tempe.matpar = matpar_ele;
		elements.push_back(tempe);
		numele++;
	}
	infile_ele.close();
	//numele = numele-1;
	cout << "Finished reading element data "<<numele<<endl;
	
	cout<<"Starting boundary conditions"<<endl;
	// reading boundary condition of elements
	int bndsdes = 0;
	ifstream infile_boundry_con(argv[3]);
	//ifstream infile_boundry_con("Data_640_dof/TractionBC.txt");
	while (!infile_boundry_con.eof()) {
		bc temp;
		double dispvecval[4];
		double tracvecval[4];
		infile_boundry_con >> temp.ele >> temp.node1 >> temp.node2; 
		temp.ele = temp.ele-1;
		temp.node1 = temp.node1-1;
		temp.node2 = temp.node2-1;
		infile_boundry_con >> temp.tractvec[0]>>temp.tractvec[1]>>temp.tractvec[2]>>temp.tractvec[3];
		infile_boundry_con >> temp.dc_x1 >> temp.dc_y1 >> temp.dc_x2 >> temp.dc_y2;
		infile_boundry_con >> temp.dispvec[0]>>temp.dispvec[1]>>temp.dispvec[2]>>temp.dispvec[3];
		bc_arr.push_back(temp);
		bndsdes++;
	}
	infile_boundry_con.close();
	//bndsdes = bndsdes-1;
	bc temp;
	vector<double> knwndisp;
	for (int a=0;a<bndsdes;a++)
	{
		temp = bc_arr[a];
		if (temp.dc_x1 == 1){
			dof = 2*temp.node1;
			cout<<"dof at x1 "<<dof<<endl;
			thereflg = 0;
			for (int i=0; i<rest_dispdof.size(); i++)
			{
				if (rest_dispdof[i] == dof){
					thereflg = 1;
					if (knwndisp[i]!=temp.dispvec[0])
						cout<<"Inconsistent displacement specification 1"<<endl;
				}
			}
			if (thereflg==0){
				rest_dispdof.push_back(dof);
				knwndisp.push_back(temp.dispvec[0]);
				count_dc++;
			}
		}
		if (temp.dc_y1 == 1){
			dof = 2*temp.node1+1;
			cout<<"dof at y1 "<<dof<<endl;			
			thereflg = 0;
			for (int i=0; i<rest_dispdof.size(); i++)
			{
				if (rest_dispdof[i] == dof){
					thereflg = 1;
					if (knwndisp[i]!=temp.dispvec[1])
						cout<<"Inconsistent displacement specification 2"<<endl;
				}
			}
			if (thereflg==0){
				rest_dispdof.push_back(dof);
				knwndisp.push_back(temp.dispvec[1]);
				count_dc++;
			}
		}
		if (temp.dc_x2 == 1){
			dof = 2*temp.node2;
			cout<<"dof at x2 "<<dof<<endl;			
			thereflg = 0;
			for (int i=0; i<rest_dispdof.size(); i++)
			{
				if (rest_dispdof[i] == dof){
					thereflg = 1;
					if (knwndisp[i]!=temp.dispvec[2])
						cout<<"Inconsistent displacement specification 3"<<endl;
				}
			}
			if (thereflg==0){
				rest_dispdof.push_back(dof);
				knwndisp.push_back(temp.dispvec[2]);				
				count_dc++;
			}
		}
		if (temp.dc_y2 == 1){
			dof = 2*temp.node2+1;
			cout<<"dof at y2 "<<dof<<endl;			
			thereflg = 0;
			for (int i=0; i<rest_dispdof.size(); i++)
			{
				if (rest_dispdof[i] == dof){
					thereflg = 1;
					if (knwndisp[i]!=temp.dispvec[3])
						cout<<"Inconsistent displacement specification 4"<<endl;
				}
			}
			if (thereflg==0){
				rest_dispdof.push_back(dof);
				knwndisp.push_back(temp.dispvec[3]);
				count_dc++;
			}
		}
	}

	// Sorting
	for(int i=0;i<count_dc;i++)
	{	
		int tempdof;
		double	tempdisp;
		for(int j=i+1;j<count_dc;j++)
		{
			if(rest_dispdof[i]>rest_dispdof[j])
			{
				tempdof=rest_dispdof[i];
				tempdisp = knwndisp[i];
				rest_dispdof[i]=rest_dispdof[j];
				knwndisp[i] = knwndisp[j];
				rest_dispdof[j]=tempdof;
				knwndisp[j] = tempdisp;
			}
		}
	}
	cout << "printing resting disp" << endl;
	for(int i = 0 ; i < count_dc ; i++){
		cout << rest_dispdof[i] << " , ";
	}
	cout << endl;

	cout <<"Finished reading boundary conditions"<<endl;
	std::cout<<"count_dc: "<<count_dc<<'\n';
	/*cout<<"Number of boundary sides "<<bndsdes<<endl;
	for (int i = 0; i < bndsdes; i++)
	{
		cout << "element " << bc_arr[i].ele << " first_node " << bc_arr[i].node1 << " second_node " << bc_arr[i].node2 << bc_arr[i].tractvec[0] << bc_arr[i].tractvec[1] << bc_arr[i].tractvec[2]<<bc_arr[i].tractvec[3]<< endl;
	}
	for(int a=0; a<count_dc; a++)
	{
		cout<<"Restrained disp dof"<<rest_dispdof[a]<<endl;
	}
	cout << " count_dc " << count_dc <<endl;
	*/

	order = 6;
	GQVandW GQdata;
	
	GQdata = getWeightAndValues(order);
	std::cout<<"numnode: "<<numnode<<" count_dc: "<<count_dc<<'\n';
	colm = 8*numnode + count_dc;                   
	row = 7*bndsdes + count_dc;

	cout << "column " << colm << " rows " << row << endl;
	vector<double> Fvector(row);

	vector<double> temp_Amatrow(colm);
	for (int i = 0; i < row; i++) {
		Amat.push_back(temp_Amatrow);
		Fvector.push_back(0);
		for (int j = 0; j < colm; j++)
			Amat[i][j] = 0;
	}

	int rowindx = 0;
	int tmpindx;
	double tractionvec[4];
	double displvec[4];
	vector<vector<double>> tempkmat;
	vector<vector<double>> tempFmat;
	vector<vector<double>> temptransmat;
	cout<<"Going to assemble Amat"<<endl;
	for (int ab=0; ab<bndsdes; ab++)
	{
		vector<int> colindx_sts;
		vector<int> colindx_disp;
		vector<int> colindx_trct;
		vector<int> ddof;
		temp = bc_arr[ab];
		RetBC res = boundry_condition_fun(temp, GQdata, order);
		tempkmat = res.Kmat;
		tempFmat = res.Fmat;
		/*for (int a=0; a<7; a++) {
			for (int b=0; b<12; b++) {
				cout << tempkmat[a][b] << "\t";
			}
			cout<<endl;
		}
		for (int a=0; a<7; a++) {
			for (int b=0; b<4; b++) {
				cout << tempFmat[a][b] << "\t";
			}
			cout<<endl;
		}*/
		tractionvec[0] = temp.tractvec[0]; 
		tractionvec[1] = temp.tractvec[1];
		tractionvec[2] = temp.tractvec[2];
		tractionvec[3] = temp.tractvec[3];
		//cout << tractionvec[0] << "\t"<< tractionvec[1] << "\t"<< tractionvec[2] << "\t"<< tractionvec[3] << endl;
		displvec[0] = temp.dispvec[0];
		displvec[1] = temp.dispvec[1];
		displvec[2] = temp.dispvec[2];
		displvec[3] = temp.dispvec[3];
		//cout<<"On side " << ab <<endl;
		tmpindx = 2*numnode + 6*temp.node1;
		for (int b=0; b<6; b++)
		{
			colindx_sts.push_back(tmpindx + b);
		}
		tmpindx = 2*numnode + 6*temp.node2;
		for (int b=0; b<6; b++)
		{
			colindx_sts.push_back(tmpindx + b);
		}
		if (temp.dc_x1 == 0 && temp.dc_x2 == 0 && temp.dc_y1 == 0 && temp.dc_y2 == 0)
		{
			cout<<"Going to assign AMat pure tract bc"<<endl;
			for (int a=0; a<7; a++)
			{
				for(int b=0;b<12;b++)
					{
						Amat[rowindx+a][colindx_sts[b]] = tempkmat[a][b];
					}
				for(int c=0;c<4;c++)
					{
						Fvector[rowindx+a] = Fvector[rowindx+a]+tempFmat[a][c]*tractionvec[c];
					}
			}
			rowindx = rowindx+7;			
		}
		else
		{
			if (temp.dc_x1 == 1 && temp.dc_x2 == 1 && temp.dc_y1 == 1 && temp.dc_y2 == 1)
			{
				cout<<"Going to assign AMat pure disp bc "<<endl;
				ddof.push_back(2*temp.node1);
				ddof.push_back(2*temp.node1+1);
				ddof.push_back(2*temp.node2);
				ddof.push_back(2*temp.node2+1);
				for (int e=0;e<4;e++)
				{
					for (int d=0;d<count_dc;d++)
					{
						if (ddof[e]==rest_dispdof[d])
							colindx_trct.push_back(d);
					}
				}
				for (int a=0; a<7; a++)
				{
					for(int b=0;b<12;b++)
						{
							Amat[rowindx+a][colindx_sts[b]] = tempkmat[a][b];
						}
					for(int c=0;c<4;c++)
						{
							Amat[rowindx+a][8*numnode+colindx_trct[c]] = -tempFmat[a][c];
							Fvector[rowindx+a] = 0;
						}
				}
			rowindx = rowindx+7;
			}
			else
			{
				cout<<"Going to assign AMat mixed bc "<<endl;
				int cntd = 0;
				int cntu = 0;
				vector<int> relind;
				vector<int> othind;
				vector<int> ddof;
				vector<int> udof; 
				if (temp.dc_x1==1){
					ddof.push_back(2*temp.node1);
					relind.push_back(0);
					cntd++;}
				else
					{udof.push_back(2*temp.node1);othind.push_back(0);cntu++;}
				if (temp.dc_y1==1)
					{ddof.push_back(2*temp.node1+1);relind.push_back(1); cntd++;}
				else
					{udof.push_back(2*temp.node1+1);othind.push_back(1);cntu++;}
				if (temp.dc_x2==1)
					{ddof.push_back(2*temp.node2);relind.push_back(2);cntd++;}
				else
					{udof.push_back(2*temp.node2);othind.push_back(2);cntu++;}
				if (temp.dc_y2==1)
					{ddof.push_back(2*temp.node2+1);relind.push_back(3);cntd++;}
				else
					{udof.push_back(2*temp.node2+1);othind.push_back(3);cntu++;}
				for (int e=0;e<cntd;e++)
				{
					for (int d=0;d<count_dc;d++)
					{
						if (ddof[e]==rest_dispdof[d])
							colindx_trct.push_back(d);
					}
				}
				for (int a=0; a<7; a++)
				{
					for(int b=0;b<12;b++)
						{
							Amat[rowindx+a][colindx_sts[b]] = tempkmat[a][b];
						}
					for(int c=0;c<cntd;c++)
						{
							cout<<"IN col " <<8*numnode+colindx_trct[c]<<"\t"<<relind[c]<<"\t"<<-tempFmat[a][relind[c]]<<endl;
							Amat[rowindx+a][8*numnode+colindx_trct[c]] = -tempFmat[a][relind[c]];
						
						}
					for(int c=0;c<cntu;c++)
						{
							Fvector[rowindx+a] = Fvector[rowindx+a] + tempFmat[a][othind[c]]*tractionvec[othind[c]];
						}
				}
			rowindx = rowindx+7;
			}
		}
	cout<<"Finished side - " << ab << endl;
	}	
	cout << "Assignment of Amat is completed" << endl;	
	for (int a=0;a<count_dc;a++)
	{
		Amat[rowindx+a][rest_dispdof[a]] = 1;
		Fvector[rowindx+a] = knwndisp[a];
	}

	int rmcount = 0;
	std::set<int> rmrws;
	int counter1 = 0, counter2 = 0, counter3 = 0;
	cout << "Starting inconstistency elimination logic" << endl;
	cout << "Amat.size : " << Amat.size() << endl;
	cout << "Amat[0].size : " << Amat[0].size() << endl;

	//std::map<int, int> deps;
	// Updated Rank Consistency logic
	constexpr double zeroset = 1e-13;
	for(int a = 0; a < 7 * bndsdes; a++) {

		//std::vector<double> bsevec(Amat[a].begin(), Amat[a].begin() + ((8 * numnode) - count_dc));
		std::vector<double> bsevec(Amat[a].begin(), Amat[a].begin() + (8 * numnode));
		
		std::vector<int> bseind (bsevec.size(), 0);
		// bseind = abs(bsevec)>0;
		for (std::size_t i = 0; i < bsevec.size(); i++){
			if (abs(bsevec[i]) > 0){
				bseind[i] = 1;
			}
		}
		
		for(int b = a + 1; b < 7 * bndsdes; b++) {
			//reqvec = kknew(b,1:(8*numnode));
			std::vector<double> reqvec(Amat[b].begin(), Amat[b].begin() + (8 * numnode));
			//std::vector<double> reqvec(Amat[b].begin(), Amat[b].begin() + ((8 * numnode) - count_dc));
			
			double diff_value1 = 0;
			//double diff_value2 = 0;
			
			for(int diff_index = 0 ; diff_index < bsevec.size() ; diff_index++) {
				diff_value1 += abs(abs(reqvec[diff_index]) - abs(bsevec[diff_index]));
				//diff_value2 += abs(reqvec[diff_index] + bsevec[diff_index]);
			}
			
			if (diff_value1 < zeroset) {
				rmrws.insert(b);
				rmcount++;
				counter1++;
				//deps[b] +=1;
				//cout<<a<<','<<b<<'\n';
			}
/*
			else if (diff_value2 < zeroset) {
				rmrws.insert(b);
				rmcount++;
				counter2++;
				//deps[b] += 1;
				cout<<a<<','<<b<<'\n';
			}
*/
			else {
				double sum, sum1 = 0, sum2 = 0;
				bool dep = true;
				std::vector<int> reqind (reqvec.size(), 0);
				// reqind = abs(reqvec)>0;
				for (std::size_t i = 0; i < reqvec.size(); i++){
					if (abs(reqvec[i]) > 0){
						reqind[i] = 1;
					}
				}

				// if length(find(reqind==1)) == length(find(bseind==1))
				//if (std::accumulate(reqind.begin(), reqind.end(), 0) == std::accumulate(bseind.begin(), bseind.end(), 0)){
				
				for (std::size_t i = 0; i < reqind.size(); i++){
					if (reqind[i] == bseind[i]){
						continue;
					}
					else{
						dep = false;
						break;
					}
				}
					
				//cout<<"Dependent vector index: "<<b<<'\n';
				if (dep == true){
					
					cout<<"Dependent vector index: "<<b<<'\n';
					//double min = 0.0, max = 0.0;
					std::vector<double> req_n, req_d, reqrat;
					auto bseind_it1 = bseind.begin(), bseind_it2 = bseind.begin();

					// reqvec(bseind)
					std::copy_if(reqvec.begin(), reqvec.end(), std::back_inserter(req_n), [&](int){return *bseind_it1++;});
					// bsevec(bseind)
					std::copy_if(bsevec.begin(), bsevec.end(), std::back_inserter(req_d), [&](int){return *bseind_it2++;});
					// reqrat = reqvec(bseind)./bsevec(bseind);
					std::transform(req_n.begin(), req_n.end(), req_d.begin(), std::back_inserter(reqrat), std::divides<double>());

					// nzerind = abs(reqrat)>0;
					//std::vector<int> nzerind (reqrat.size(), 0);
					/*
					for (std::size_t i = 0; i < reqrat.size(); i++){
						if (abs(reqrat[i]) > 0){
							nzerind[i] = 1;
						}
					}
					*/
					// sum(reqrat)
					sum = std::accumulate(reqrat.begin(), reqrat.end(), 0.0);

					//std::vector<double> f_reqrat;
					//auto nzerind_it = nzerind.begin();

					// reqrat(nzerind)
					//std::copy_if(reqrat.begin(), reqrat.end(), std::back_inserter(f_reqrat), [&](int){return *nzerind_it++;});

					//if (f_reqrat.size() > 0){
						// min(reqrat(nzerind))
					//min = *std::min_element(f_reqrat.begin(), f_reqrat.end());
					double min = *std::min_element(reqrat.begin(), reqrat.end());
					
						//max(reqrat(nzerind))
					//double max = *std::max_element(f_reqrat.begin(), f_reqrat.end());
					double max = *std::max_element(reqrat.begin(), reqrat.end());
					//}

					for (const auto& elem: reqrat){
						sum1 += abs(elem - max);
					}

					for (const auto& elem: reqrat){
						sum2 += abs(elem - min);
					}

					if (sum > 0 && sum1 < zeroset){
						rmrws.insert(b);
						rmcount++;
						counter3++;	
					}
					else if (sum < 0 && sum2 < zeroset){
						rmrws.insert(b);
						rmcount++;
						counter3++;
					}

				}

			}
		}
	}

	cout << "Rank before removing rows: "<<compute_rank(Amat)<<'\n';

	cout<<"counter1 "<<counter1<<'\n';
	cout<<"counter2 "<<counter2<<'\n';
	cout<<"counter3 "<<counter3<<'\n';
	cout<<"counter sum values "<<counter1+counter2+counter3<<'\n';

/*
	ofstream outfile_Am;
	outfile_Am.open("Amat.csv");
	for (int i=0;i<Amat.size();i++) {
		for (int j =0;j<Amat[i].size() - 1;j++) {
			outfile_Am << Amat[i][j] << ",";
		}
		outfile_Am << Amat[i][colm-1] << "\n";
	}
	outfile_Am.close();
	cout<<"Wrote AMat "<<endl;
*/	
	//for (const auto& [key, value] : deps)
        //	std::cout << '[' << key << "] = " << value << "; ";
	//cout<<endl;
	//cout<<"size of rset "<<rset.size()<<'\n';
	std::sort(rset.begin(), rset.end());
	//for (const auto &e: rset){
	//	cout<<e<<" ";
	//}
	//cout<<'\n';
	
	//rmrws.erase(486);
	//cout << "Size of set " << rmrws.size()<<'\n';
	//for (const auto &e: rmrws){
	//	cout<<e<<" ";
	//}
	//cout<<'\n';
/*
	std::vector<int> diff;

	std::set_intersection(rset.begin(), rset.end(), rmrws.begin(), rmrws.end(), std::back_inserter(diff));
	cout<<"Diff ";
	for (const auto &e: diff){
		cout<<e<<" ";
	}
	cout<<'\n';	
*/

	std::set<int>::reverse_iterator it = rmrws.rbegin();

	std::vector<int> all(7*bndsdes+count_dc);
	std::iota(all.begin(), all.end(), 0);

	std::vector<int> diff;

	std::set_difference(all.begin(), all.end(), rset.begin(), rset.end(), std::back_inserter(diff));

	std::sort(diff.begin(), diff.end());

	auto it2 = diff.rbegin();	

	cout << "Removing rows\n";
	while (it2 != diff.rend()) {
		Amat.erase(Amat.begin() + *it2);
		Fvector.erase(Fvector.begin() + *it2);
		it2++;
	}
	cout << "Just completed inconsistent rank code\n";
	cout << "Amat.size : " << Amat.size() << endl;
	cout << "Amat[0].size : " << Amat[0].size() << endl;
	row = Amat.size();
	cout << "Rank after removing rows: "<<compute_rank(Amat)<<'\n';
	
	//std::vector<int> all(7*bndsdes+count_dc);
	//std::iota(all.begin(), all.end(), 0);

	//std::vector<int> diff;

	//std::set_difference(all.begin(), all.end(), rmrws.begin(), rmrws.end(), std::back_inserter(diff));

	//cout<<"Length rqrws " <<diff.size()<<'\n';

	if((Amat.size() - compute_rank(Amat)) > 0) {
		cout << "Inconsistent boundary condition" << endl;
		cout << "Amat.size : " << Amat.size() << endl;
		cout << "compute_rank(Amat) : " << compute_rank(Amat) << endl;

		//exit(0);
	}

	ofstream outfile_A;
	outfile_A.open("Amat.csv");
	for (int i=0;i<Amat.size();i++) {
		for (int j =0;j<Amat[i].size() - 1;j++) {
			outfile_A << Amat[i][j] << ",";
		}
		outfile_A << Amat[i][colm-1] << "\n";
	}
	outfile_A.close();
	cout<<"Wrote AMat "<<endl;

	ofstream outfile_B;
	outfile_B.open("Fvec.csv");
	for (int i=0;i<row;i++) {
		outfile_B << Fvector[i] << "\n";
	}
	outfile_B.close();
	cout<<"Wrote Fvec "<<endl;
	//int keycmd;
	cout<<"Going to optimize. Do you want to continue?"<<endl;
	//cin>>keycmd;
	Amat_rows = Amat.size();
	Amat_cols = Amat[0].size();
	obj_data.elemdat = elements;
	obj_data.nodecord = nodes;
	obj_data.dispscalefac = 1;
	obj_data.numnode = numnode;
	obj_data.numele = numele;
	obj_data.lnpos = count_dc;
	
	nn = numnode;
	nl = count_dc;

	const unsigned n = 8*numnode+count_dc;

	vector<double> temp_intmat(n, 0);
	for (int i = 0; i < n; i++) {
		intmat.push_back(temp_intmat);
		for (int j = 0; j < n; j++) {
			for(int k=0; k<Amat.size(); k++){
				intmat[i][j] = intmat[i][j] + Amat[k][i]*Amat[k][j];
			}
		}
	}

	constraint_data.Amat = Amat;
	constraint_data.Fvec = Fvector;
	constraint_data.rows = row;
	constraint_data.intmat = intmat;
	const double toler = 1e-14;
	unsigned m = row;
	//double tol[m];
	double *tol = new double[m];
	for(int i=0;i<m;i++){
		tol[i] = toler;
	}
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LD_SLSQP, n);
	nlopt_set_min_objective(opt, myfunc, &obj_data);
	nlopt_add_equality_mconstraint(opt, m, myconstraint, &constraint_data, tol);
	// nlopt_add_equality_constraint(opt, myconstraint, &constraint_data, toler);
	nlopt_set_xtol_rel(opt, toler);
	nlopt_set_ftol_rel(opt, toler);
	nlopt_set_ftol_abs(opt, toler/10000);
	double x[n];
	for (int i=0; i<2*numnode; i++){
		x[i] = 0;
	}
	for(int i=2*numnode; i<(8*numnode+count_dc); i++){
		x[i] = 20;
	}	
	double minf;
	double *xreq;
	xreq = x;
	double start = omp_get_wtime();
	int answer= nlopt_optimize(opt, xreq, &minf);
	double end = omp_get_wtime();

	cout << "Time taken : " << end - start << endl;
	if (answer<0){
		cout<<"Optimization did not converge. Exit code: "<<answer<<endl;
	}
	else {
		cout<<"Optimization converged"<<endl<<"minf = "<<minf<<"exit flag = " << answer<<endl;
	}
	cout << "Overall time taken by objective function : " << obj_exec_time << endl; 
	ofstream outoptfile;

        time_t tt; 
   
        // localtime() 
        struct tm* ti; 
  
        time(&tt); 

        ti = localtime(&tt); 
  	char filename[80];  // Choose an appropriate size
    	strftime(filename, sizeof(filename), "file_%Y%m%d_%H%M%S.txt", ti);
		std::cout << ti->tm_hour << ti->tm_min << ti->tm_sec << std::endl;
	std::string filenameStr(filename);
         
	outoptfile.open(filenameStr);
	for (int i=0; i<(8*numnode+count_dc); i++){
		outoptfile<< x[i]<< endl;
	}
	outoptfile.close();
	nlopt_destroy(opt);

	return 0;
}


RetBC boundry_condition_fun(bc boundry, GQVandW GQdat, int order) {
	cout << "In boundary condition" << endl;
	RetBC result;
	double bc_x1, bc_y1, bc_x2, bc_y2;
	int c1, c2, nd1, nd2;
	vector<double> values;
	vector<double> weights;
	values = GQdat.value;
	weights = GQdat.weight;
	
	nd1 = boundry.node1;
	nd2 = boundry.node2;
	bc_x1 = nodes[nd1].x_cord;
	bc_y1 = nodes[nd1].y_cord;
	bc_x2 = nodes[nd2].x_cord;
	bc_y2 = nodes[nd2].y_cord;
	
	vector<vector<double> > kmat;
	vector<double> kmatrow(12);
	// inti with Zero
	for (int i = 0; i < 7; i++) {
		kmat.push_back(kmatrow);
		for (int j = 0; j < 12; j++)
			kmat[i][j] = 0;
	}

	vector<vector<double> > fmat;
	vector<double>fmatrow(4);
	//init with zero
	for (int i = 0; i < 7; i++) {
		fmat.push_back(fmatrow);
		for (int j = 0; j < 4; j++)
			fmat[i][j] = 0;
	}

	double len, nx, ny, f1, f2, m1, m2, m3, m4, var, value, value2, value3;
	len = sqrt(pow((bc_x2 - bc_x1), 2) + pow((bc_y2 - bc_y1), 2));
/*
	cout<<"nd1: "<<nd1<<'\n';
	cout<<"nd2: "<<nd2<<'\n';
	cout<<"bc_x1: "<<bc_x1<<'\n';
	cout<<"bc_x2: "<<bc_x2<<'\n';
	
	cout<<"bc_y1: "<<bc_y1<<'\n';
	cout<<"bc_y2: "<<bc_y2<<'\n';	

	cout<<"len : "<<len<<'\n';
*/
	nx = (bc_y2 - bc_y1) / len;
	ny = -(bc_x2 - bc_x1) / len;

	if (std::isnan(nx)){
		cout<<"nx is Nan\n";
	}

	if (std::isnan(ny)){
		cout<<"ny is Nan\n";
	}

	int condition = 0;
	int ele = boundry.ele;

	if (boundry.node1 == (elements[ele].node1-1) && boundry.node2 == (elements[ele].node2-1)) // t = -1
	{
		condition = 1;
		c1 = 0, c2 = 6;
	}
	else if (boundry.node1 == (elements[ele].node2-1) && boundry.node2 == (elements[ele].node3-1)) // s = 1
	{
		condition = 2;
		c1 = 6, c2 = 12;
	}
	else if (boundry.node1 == (elements[ele].node3-1) && boundry.node2 == (elements[ele].node4-1)) // t = 1
	{
		condition = 3;
		c1 = 12, c2 = 18;
	}
	else if (boundry.node1 == (elements[ele].node4-1) && boundry.node2 == (elements[ele].node1-1)) // s = -1
	{
		condition = 4;
		c1 = 18, c2 = 0;
	}
	cout << "condition " << condition << endl;

	vector<BB> BB_results_mat;
	
	double xcoord[4];
	xcoord[0] = nodes[elements[ele].node1-1].x_cord;
	xcoord[1] = nodes[elements[ele].node2-1].x_cord;
	xcoord[2] = nodes[elements[ele].node3-1].x_cord;
	xcoord[3] = nodes[elements[ele].node4-1].x_cord;

	double ycoord[4];
	ycoord[0] = nodes[elements[ele].node1-1].y_cord;
	ycoord[1] = nodes[elements[ele].node2-1].y_cord;
	ycoord[2] = nodes[elements[ele].node3-1].y_cord;
	ycoord[3] = nodes[elements[ele].node4-1].y_cord;
	
	BB_results_mat = get_BB_results_mat(condition, order, xcoord, ycoord, values);

	if (condition == 2 || condition == 4) {
	value = 0.0;
	value2 = 0.0;
	value3 = 0.0;
	for (int i = 0; i < order; i++) {
		if (condition==2)
			var = sqrt(pow(BB_results_mat[i].Jmat[0][1], 2) + pow(BB_results_mat[i].Jmat[1][1], 2));
		else
		 	var = -sqrt(pow(BB_results_mat[i].Jmat[0][1], 2) + pow(BB_results_mat[i].Jmat[1][1], 2));
		value += weights[i] * var;
		value2 += weights[i] * var * BB_results_mat[i].x;
		value3 += weights[i] * var * BB_results_mat[i].y;
	}
	f1 = value; // h1
	m1 = value2; // x*h1
	m3 = value3; // y*h1
	value = 0.0;
	value2 = 0.0;
	value3 = 0.0;
	for (int i = 0; i < order; i++) {
		if (condition==2)
			var = sqrt((BB_results_mat[i].x - bc_x1) * (BB_results_mat[i].x - bc_x1) + (BB_results_mat[i].y - bc_y1) * (BB_results_mat[i].y - bc_y1)) * sqrt(pow(BB_results_mat[i].Jmat[0][1], 2) + pow(BB_results_mat[i].Jmat[1][1], 2))/ len;
		else
			var = -sqrt((BB_results_mat[i].x - bc_x1) * (BB_results_mat[i].x - bc_x1) + (BB_results_mat[i].y - bc_y1) * (BB_results_mat[i].y - bc_y1)) * sqrt(pow(BB_results_mat[i].Jmat[0][1], 2) + pow(BB_results_mat[i].Jmat[1][1], 2))/ len;
		value += weights[i] * var;
		value2 += weights[i] * var * BB_results_mat[i].x;
		value3 += weights[i] * var * BB_results_mat[i].y;
	}
	f2 = value; // h2 
	m2 = value2;  // x*h2
	m4 = value3; // y*h2
	fmat[0][0] = 1;
	fmat[1][2] = 1;
	fmat[2][1] = 1;
	fmat[3][3] = 1;

	}
	
	else if (condition == 1 || condition == 3) {
	value = 0.0;
	value2 = 0.0;
	value3 = 0.0;
	for (int i = 0; i < order; i++) {
		if (condition==1)
			var = sqrt(pow(BB_results_mat[i].Jmat[1][0], 2) + pow(BB_results_mat[i].Jmat[0][0], 2));
		else
			var = -sqrt(pow(BB_results_mat[i].Jmat[1][0], 2) + pow(BB_results_mat[i].Jmat[0][0], 2));
		value += weights[i] * var;
		value2 += weights[i] * var * BB_results_mat[i].x;
		value3 += weights[i] * var * BB_results_mat[i].y;
	}
	f1 = value; // h1
	m1 = value2; // x*h1
	m3 = value3; // y*h1
	value = 0.0;
	value2 = 0.0;
	value3 = 0.0;
	for (int i = 0; i < order; i++) {
		if (condition==1)
			var = sqrt((BB_results_mat[i].x - bc_x1) * (BB_results_mat[i].x - bc_x1) + (BB_results_mat[i].y - bc_y1) * (BB_results_mat[i].y - bc_y1)) * sqrt(pow(BB_results_mat[i].Jmat[1][0], 2) + pow(BB_results_mat[i].Jmat[0][0], 2))/ len;
		else
			var = -sqrt((BB_results_mat[i].x - bc_x1) * (BB_results_mat[i].x - bc_x1) + (BB_results_mat[i].y - bc_y1) * (BB_results_mat[i].y - bc_y1)) * sqrt(pow(BB_results_mat[i].Jmat[1][0], 2) + pow(BB_results_mat[i].Jmat[0][0], 2))/ len;
		value += weights[i] * var;
		value2 += weights[i] * var * BB_results_mat[i].x;
		value3 += weights[i] * var * BB_results_mat[i].y;
	}
	f2 = value; // h2 
	m2 = value2;  // x*h2
	m4 = value3; // y*h2
	fmat[0][1] = 1;
	fmat[1][3] = 1;
	fmat[2][0] = 1;
	fmat[3][2] = 1;
	}

	fmat[4][0] = f1 - f2;
	fmat[4][2] = f2;
	fmat[5][1] = f1 - f2;
	fmat[5][3] = f2;
	fmat[6][0] = m4 - m3;
	fmat[6][1] = m1 - m2;
	fmat[6][2] = -m4;
	fmat[6][3] = m2;
	
	if (condition == 2 || condition == 4) {
	kmat[0][4] = nx;
	kmat[0][5] = -ny;
	kmat[2][3] = ny;
	kmat[2][5] = -nx;
	kmat[1][10] = nx;
	kmat[1][11] = -ny;
	kmat[3][9] = ny;
	kmat[3][11] = -nx;
	
	for (int i = 0; i < 6; i++) {
		value = 0;
		value2 = 0;
		value3 = 0;
		for (int j = 0; j < order; j++)
		{
			if (condition==2)
				var = sqrt(pow(BB_results_mat[j].Jmat[0][1], 2) + pow(BB_results_mat[j].Jmat[1][1], 2));
			else
				var = -sqrt(pow(BB_results_mat[j].Jmat[0][1], 2) + pow(BB_results_mat[j].Jmat[1][1], 2));
			value += weights[j] * (BB_results_mat[j].Bsmat[0][i + c1] * nx + BB_results_mat[j].Bsmat[2][i + c1] * ny) * var;
			value2 += weights[j] * (BB_results_mat[j].Bsmat[2][i + c1] * nx + BB_results_mat[j].Bsmat[1][i + c1] * ny) * var;
			value3 += weights[j] * (BB_results_mat[j].x * BB_results_mat[j].Bsmat[1][i + c1] * ny - BB_results_mat[j].y * BB_results_mat[j].Bsmat[0][i + c1] * nx + BB_results_mat[j].Bsmat[2][i + c1] * (BB_results_mat[j].x * nx - BB_results_mat[j].y * ny)) * var;
		}
		kmat[4][i] = value; // (Bsmat[0][i]*nx*h1) + (Bsmat[2][i]*ny*h1))
		kmat[5][i] = value2; // (Bsmat[2][i]*nx*h1) + (Bsmat[1][i]*ny*h1))
		kmat[6][i] = value3;	// (x*Bsmat[1][i]*ny*h1) - (y*Bsmat[0][i]*nx*h1) + Bsmat[2][i]*(x*nx - y*ny)*h1)
		value = 0;
		value2 = 0;
		value3 = 0;
		for (int j = 0; j < order; j++)
		{
			if (condition==2)
				var = sqrt(pow(BB_results_mat[j].Jmat[0][1], 2) + pow(BB_results_mat[j].Jmat[1][1], 2));
			else
				var = -sqrt(pow(BB_results_mat[j].Jmat[0][1], 2) + pow(BB_results_mat[j].Jmat[1][1], 2));
			value += weights[j] * (BB_results_mat[j].Bsmat[0][i + c2] * nx + BB_results_mat[j].Bsmat[2][i + c2] * ny) * var;
			value2 += weights[j] * (BB_results_mat[j].Bsmat[2][i + c2] * nx + BB_results_mat[j].Bsmat[1][i + c2] * ny) * var;
			value3 += weights[j] * (BB_results_mat[j].x * BB_results_mat[j].Bsmat[1][i + c2] * ny - BB_results_mat[j].y * BB_results_mat[j].Bsmat[0][i + c2] * nx + BB_results_mat[j].Bsmat[2][i + c2] * (BB_results_mat[j].x * nx - BB_results_mat[j].y * ny)) * var;
		}
		kmat[4][i + 6] = value; // (Bsmat[0][i]*nx*h1) + (Bsmat[2][i]*ny*h1))
		kmat[5][i + 6] = value2; // (Bsmat[2][i]*nx*h1) + (Bsmat[1][i]*ny*h1))
		kmat[6][i + 6] = value3;	//(x*Bsmat[1][i]*ny*h1) - (y*Bsmat[0][i]*nx*h1) + Bsmat[2][i]*(x*nx - y*ny)*h1)
	}
	}

	else if (condition == 1 || condition == 3) {
	kmat[0][5] = -nx;
	kmat[0][3] = ny;
	kmat[2][5] = -ny;
	kmat[2][4] = nx;
	kmat[1][11] = -nx;
	kmat[1][9] = ny;
	kmat[3][11] = -ny;
	kmat[3][10] = nx;

		
		for (int i = 0; i < 6; i++) {
		value = 0;
		value2 = 0;
		value3 = 0;
		for (int j = 0; j < order; j++)
		{
			if (condition==1)
				var = sqrt(pow(BB_results_mat[j].Jmat[1][0], 2) + pow(BB_results_mat[j].Jmat[0][0], 2));
			else
				var = -sqrt(pow(BB_results_mat[j].Jmat[1][0], 2) + pow(BB_results_mat[j].Jmat[0][0], 2));
			value += weights[j] * (BB_results_mat[j].Bsmat[0][i + c1] * nx + BB_results_mat[j].Bsmat[2][i + c1] * ny) * var;
			value2 += weights[j] * (BB_results_mat[j].Bsmat[2][i + c1] * nx + BB_results_mat[j].Bsmat[1][i + c1] * ny) * var;
			value3 += weights[j] * (BB_results_mat[j].x * BB_results_mat[j].Bsmat[1][i + c1] * ny - BB_results_mat[j].y * BB_results_mat[j].Bsmat[0][i + c1] * nx + BB_results_mat[j].Bsmat[2][i + c1] * (BB_results_mat[j].x * nx - BB_results_mat[j].y * ny)) * var;
		}
		kmat[4][i] = value; // (Bsmat[0][i]*nx*h1) + (Bsmat[2][i]*ny*h1))
		kmat[5][i] = value2; // (Bsmat[2][i]*nx*h1) + (Bsmat[1][i]*ny*h1))
		kmat[6][i] = value3;	// (x*Bsmat[1][i]*ny*h1) - (y*Bsmat[0][i]*nx*h1) + Bsmat[2][i]*(x*nx - y*ny)*h1)
		value = 0;
		value2 = 0;
		value3 = 0;
		for (int j = 0; j < order; j++)
		{
			if (condition==1)
				var = sqrt(pow(BB_results_mat[j].Jmat[1][0], 2) + pow(BB_results_mat[j].Jmat[0][0], 2));
			else
				var = -sqrt(pow(BB_results_mat[j].Jmat[1][0], 2) + pow(BB_results_mat[j].Jmat[0][0], 2));
			value += weights[j] * (BB_results_mat[j].Bsmat[0][i + c2] * nx + BB_results_mat[j].Bsmat[2][i + c2] * ny) * var;
			value2 += weights[j] * (BB_results_mat[j].Bsmat[2][i + c2] * nx + BB_results_mat[j].Bsmat[1][i + c2] * ny) * var;
			value3 += weights[j] * (BB_results_mat[j].x * BB_results_mat[j].Bsmat[1][i + c2] * ny - BB_results_mat[j].y * BB_results_mat[j].Bsmat[0][i + c2] * nx + BB_results_mat[j].Bsmat[2][i + c2] * (BB_results_mat[j].x * nx - BB_results_mat[j].y * ny)) * var;
		}
		kmat[4][i + 6] = value; // (Bsmat[0][i]*nx*h1) + (Bsmat[2][i]*ny*h1))
		kmat[5][i + 6] = value2; // (Bsmat[2][i]*nx*h1) + (Bsmat[1][i]*ny*h1))
		kmat[6][i + 6] = value3;	// (x*Bsmat[1][i]*ny*h1) - (y*Bsmat[0][i]*nx*h1) + Bsmat[2][i]*(x*nx - y*ny)*h1)
	}
	}
	for (int a=0; a<7; a++) {
		for (int b=0; b<12; b++) {
			if (kmat[a][b] <= TOL && kmat[a][b] >= -TOL)
				kmat[a][b] = 0;
		}
	}
	for (int a=0; a<7; a++) {
		for (int b=0; b<4; b++) {
			if (fmat[a][b] <= TOL && fmat[a][b] >= -TOL)
				fmat[a][b] = 0;
		}
	}
	cout<<"Finished Kmat and Fmat"<<endl;
	result.Fmat = fmat;
	result.Kmat = kmat;
	return result;
}


Cnstrelret LinElast(double *stress, double *strain, vector<double> matpar) {
	vector<vector<double>> Cmat;
	vector<vector<double>> Emat;
	vector<double> delvalue;
	vector<double> tempmat(3);
	//cout<<"In constitutive relation : s1 = "<<stress[0]<<" s2 = "<<stress[1]<<"s3= "<<stress[2]<<endl;
	for (int i = 0; i < 3; i++) {
		delvalue.push_back(0);
		Cmat.push_back(tempmat);
		Emat.push_back(tempmat);
		for (int j = 0; j < 3; j++){
			Cmat[i][j] = 0;
			Emat[i][j] = 0;
		}
	}
	Cmat[0][0] = Cmat[1][1] = -1/matpar[0];
	Cmat[0][1] = Cmat[1][0] = matpar[1]/matpar[0];
	Cmat[2][2] = -2*(1 + matpar[1])/matpar[0];
	Emat[0][0] = Emat[1][1] = Emat[2][2] = 1;
	
	for (int i = 0; i < 3; i++) {
		delvalue[i] = strain[i];
		for (int j = 0; j < 3; j++) {
			delvalue[i] = delvalue[i] + Cmat[i][j]*stress[j];
		}
	}
	Cnstrelret LEresult;
	LEresult.delval = delvalue;
	LEresult.dfdstrain = Emat;
	LEresult.dfdstress = Cmat;
	return LEresult;
}

Errret ErrGQ(double *uval, double *phival, double *xcord, double *ycord, int mattype, vector<double> mtpar, double dispscalefac, int order) {
	GQVandW GQresult;
	Cnstrelret Cnstreldat;
	vector<double> values;
	vector<double> weights;
	GQresult = getWeightAndValues(order);
	values = GQresult.value;
	weights = GQresult.weight;
	double sts[3];	
	double strn[3];
	vector<vector<double>> dfdstrn;
	vector<vector<double>> dfdsts;
	double errval = 0;
	double reqBd[8][3];
	double reqBs[24][3];	
	double ki[32];
	vector<double> Ivals;
	vector<double> delvalue;
	double objint = 0;
	double detJ=0;
	Errret ErrAGradval;
	for (int k=0;k<32;k++){
		Ivals.push_back(0);
	}
	for(int i=0; i<order; i++) {
		for (int j = 0; j<order; j++) {			
			BB BBresult = Bsigma_Bd(xcord,ycord,values[i],values[j]);
			detJ = BBresult.detJ;
			for (int k = 0; k < 3; k++) {
				sts[k] = 0;
				for (int l = 0; l < 24; l++) {
					sts[k] = sts[k] + BBresult.Bsmat[k][l]*phival[l];
				}
			}
			for (int k = 0; k < 3; k++) {
				strn[k] = 0;
				for (int l = 0; l < 8; l++) {
					strn[k] = strn[k] + BBresult.Bd[k][l]*uval[l];
				}
			}
			switch (mattype){
				case 1: Cnstreldat = LinElast(sts,strn,mtpar);
					break;
				//case 2: Cnstreldat = concmod(sts,strn,mtpar);
				//	break;
					
			}
			delvalue = Cnstreldat.delval;
			dfdstrn = Cnstreldat.dfdstrain;
			dfdsts = Cnstreldat.dfdstress;
			//cout<<"Printing delvalue : "<<delvalue[0] << "\t"<<delvalue[1] << "\t"<<delvalue[2] << "\t"<<endl;
			//int brkvar;
			//cin>>brkvar;
			objint = 0;
			for(int k=0; k<3; k++){
				objint = objint + delvalue[k]*delvalue[k];
			}
			//cout<<objint<<endl;
			errval = errval + detJ*weights[i]*weights[j]*sqrt(objint);
			//cout<<"OBjint = "<<sqrt(objint)<<endl;
			//cout<<"detJ = "<<detJ<<endl;			
			for (int k = 0; k < 8; k++) {		
				//cout<<endl;		
				for (int l = 0; l <3; l++) {
					reqBd[k][l] = 0;
					for(int m=0; m<3; m++){
						reqBd[k][l] = reqBd[k][l] + BBresult.Bd[m][k]*dfdstrn[l][m];
					}
					//cout<<reqBd[k][l]<<" ";
				}
			}
			//cout<<"Printing ki values: "<< endl;
			
			for (int k = 0; k < 8; k++) {
				ki[k] = 0;
				for (int l = 0; l < 3; l++) {
					ki[k] = ki[k] + reqBd[k][l]*delvalue[l]/dispscalefac/sqrt(objint);
				}
				Ivals[k] = Ivals[k]+detJ*weights[i]*weights[j]*ki[k];
				//cout<<ki[k]<<"\t";
			}
			for (int k = 0; k < 24; k++) {	
				//cout<<endl;			
				for (int l = 0; l <3; l++) {
					reqBs[k][l] = 0;
					for(int m=0; m<3; m++){
						reqBs[k][l] = reqBs[k][l] + BBresult.Bsmat[m][k]*dfdsts[l][m];
					}
					//cout<<reqBs[k][l]<<",\t";
				}
			}
			for (int k = 0; k < 24; k++) {
				ki[k+8] = 0;
				for (int l = 0; l < 3; l++) {
					ki[k+8] = ki[k+8] + reqBs[k][l]*delvalue[l]/sqrt(objint);
				}
				Ivals[k+8] = Ivals[k+8] + detJ*weights[i]*weights[j]*ki[k+8];
				//cout<<ki[k]<<"\t";
			}
			//cout<<"Finished printing ki for i = " <<i<<" and j = " <<j<<endl;
		}
		//cout<<weights[i]<<endl;
	}
	
	/*cout << "Error from myfunc : " << errval << endl;
	for(int i = 0 ; i < 32 ; i++) {
		cout << Ivals[i] << ", ";
	}
	cout << endl;*/
	ErrAGradval.errv = errval;
	ErrAGradval.grdval = Ivals;
	return ErrAGradval;
}

//vector<element> global_ele;
//int ele_ctr = 0;
double *grad2 = (double*) malloc((8*nn+nl) * sizeof(double));
double myfunc(unsigned n, const double *x, double *grad, void *f_data)
{
	//cout<<"Objective function value: ";
	double s1 = omp_get_wtime();

	objdata *d = (objdata *)f_data;
	vector<node> nodecords = d->nodecord;

	vector<element> ele = d->elemdat;
	Errret ErrAGrad;
	int numnodes = d->numnode;
	int numelem = d->numele;
	int lenpos = d->lnpos;
	double dispscalefac = d->dispscalefac;
	double xcoord[4];
	double ycoord[4];
	int u_pos[8];
	int phi_pos[24];
	vector<double> gradval;
	double u_val[8];
	double phi_val[24];
	int elenodes[4];
	int numpar;
	int GQorder = 10;
	double errval=0;
	double ScaleFac = 1e3;
	if (grad == nullptr){
		//grad = (double*) malloc((8*numnodes+lenpos) * sizeof(double));
		free(grad2);
		grad2 = (double*) malloc((8*numnodes+lenpos) * sizeof(double));
		
		grad = grad2;	
	}
	for (int i=0; i<(8*numnodes+lenpos); i++){
		grad[i] = 0;
	}
	for (int e=0; e<numelem; e++){
		vector<double> matparv;
		elenodes[0] = ele[e].node1-1;
		elenodes[1] = ele[e].node2-1;
		elenodes[2] = ele[e].node3-1;
		elenodes[3] = ele[e].node4-1;		
		xcoord[0] = nodes[elenodes[0]].x_cord;
		xcoord[1] = nodes[elenodes[1]].x_cord;
		xcoord[2] = nodes[elenodes[2]].x_cord;
		xcoord[3] = nodes[elenodes[3]].x_cord;
		ycoord[0] = nodes[elenodes[0]].y_cord;
		ycoord[1] = nodes[elenodes[1]].y_cord;
		ycoord[2] = nodes[elenodes[2]].y_cord;
		ycoord[3] = nodes[elenodes[3]].y_cord;
		numpar = ele[e].nummatpar;
		for(int i=0;i<numpar;i++){
			matparv.push_back(ele[e].matpar[i]);
		}
		for(int i=0; i<4; i++){
			u_pos[2*i] = 2*elenodes[i];
			u_pos[2*i+1] = 2*elenodes[i]+1;
			u_val[2*i] = x[u_pos[2*i]]/dispscalefac;
			u_val[2*i+1] = x[u_pos[2*i+1]]/dispscalefac;			
			for (int c=0; c<6;c++){
				phi_pos[6*i+c] = numnodes*2+6*elenodes[i]+c;
				phi_val[6*i+c] = x[phi_pos[6*i+c]];
			}
		}
		ErrAGrad = ErrGQ(u_val,phi_val,xcoord,ycoord,ele[e].mattype,matparv,dispscalefac,GQorder);
		gradval = ErrAGrad.grdval;
		for (int i = 0; i < 8; i++) {
			grad[u_pos[i]] = grad[u_pos[i]] + gradval[i]*ScaleFac;
		}
		for (int i=0; i<24; i++){
			grad[phi_pos[i]] = grad[phi_pos[i]] + gradval[i+8]*ScaleFac;
		}
		errval = errval + ErrAGrad.errv;			
	}
	errval = errval*ScaleFac;

	//cout<<errval<<endl;
	double e1 = omp_get_wtime();
	obj_exec_time += e1 - s1;
	return errval;

}

int ctr = 0;
double* grad1 = (double*) malloc(Amat_rows * Amat_cols * sizeof(double));
void myconstraint(unsigned m, double *result, unsigned n, const double *x, double* grad, void *cdat)
{
	constdata *d = (constdata *) cdat;
	vector<vector<double>> Amatrix = d->Amat;
	vector<double> frcvec = d->Fvec;
	double ScaleFac = 1;
	if (grad == nullptr){
		//ctr++;
		//std::cout<<"Grad nullptr: "<<ctr<<'\n';
		free(grad1);
		//grad1 = (double*) malloc(Amat_rows * Amat_cols * sizeof(double));
		grad1 = (double*) malloc(Amatrix.size() * Amatrix[0].size() * sizeof(double));
		for (int i = 0; i < Amatrix.size()*Amatrix[0].size(); i++) grad1[i] = 0;
		
		grad = grad1;	
	}
	//else{
	//	ctr++;
	//	std::cout<<"Grad not nullptr: "<<ctr<<'\n';

	//}
	for (int i = 0; i < Amatrix.size(); i++) {
		result[i] = -frcvec[i];
		for (int j = 0; j < n; j++) {
			result[i] = result[i] + Amatrix[i][j]*x[j];
			grad[i*n+j] = Amatrix[i][j];
		}
		result[i] = result[i]*ScaleFac;
	}
	//cout << "completed constrained funtion from planeproblem\n";
	//free(grad);
	return;

}

/*double myconstraint(unsigned n, const double *x, double *grad, void *cdat)
{
	constdata *d = (constdata *) cdat;
	vector<vector<double>> Amatrix = d->Amat;
	vector<vector<double>> intmat = d->intmat;
	vector<double> frcvec = d->Fvec;
	int m = d->rows;
	double ScaleFac = 1;
	cout<<"Constraint result: ";
	vector<double> intvec(n, 0);
	vector<double> intxvec(n, 0);
	double result = 0;
	if (grad == nullptr)
		grad = (double*) malloc(n * sizeof(double));
	for(int i=0; i<n;i++){
		for (int k=0; k<n; k++){
			intxvec[i] = intxvec[i] + intmat[i][k]*x[k];
		}
		for(int k=0; k<m;k++){
			intvec[i] = intvec[i] + Amatrix[k][i]*frcvec[k];
		}
		grad[i] = 2*intxvec[i] -2*intvec[i];
		result = result - 2*x[i]*intvec[i] + x[i]*intxvec[i];
	}
	for(int k=0; k<m; k++){
		result = result  + frcvec[k]*frcvec[k];
	}
	cout<<result<<endl;
	return result;
}*/
