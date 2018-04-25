#include <iomanip>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>
#include <Windows.h>
#include <cmath>
#include <fstream>
#include "Diploma.h"


using namespace std;

const int N = 4;

long double k1 = 0.51;
long double k2 = 0.07;

long double **y;
long double **y_;
long double H;
long double B_;
int numP;


void show(long double** A, int N){
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			cout << setprecision(4) << setw(8) << A[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

string LDToStr(long double one){
   ostringstream ss;
    ss << one;

    return ss.str();
}

void show(ofstream &fout, long double* a, int N){
	for (int j = 0; j < N; j++)
	{
		fout << a[j] << " ";
		cout << a[j] << " ";
	}
	fout << " \n";
	cout << endl;

}

void show(long double* a, int N){
	for (int j = 0; j < N; j++)
	{
		cout << a[j] << " ";
	}
	cout << endl;

}

long double* Diploma::f(long double* y, int N){
	long double* result = new long double[N];
	result[0] = - k1 * y[0] - 2 * k2 * y[0] * y[0];
	result[1] = k1 * y[0] + k2 * y[0] * y[0];
	result[2] = k1 * y[0];
	result[3] = 2 * k2 * y[0] * y[0];

	return result;
}

long double* Diploma::Add(long double* x, long double* y, int N){
	long double* result = new long double[N];
	for (int i = 0; i < N; i++)
		{
		    result[i] = x[i] + y[i];
		}

	return result;
}

long double* Diploma::Mult(long double a, long double* y, int N){
	long double* result = new long double[N];
	for (int i = 0; i < N; i++)
		result[i] = a * y[i];
	return result;
}

long double** Diploma::Add(long double** A, long double** B, int N){
	long double** Result = new long double*[N];
	for (int i = 0; i < N; i++)
		Result[i] = new long double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Result[i][j] = A[i][j] + B[i][j];
	return Result;
}

long double** Diploma::Mult(long double** A, long double** B, int N){
	long double** Result = new long double*[N];
	for (int i = 0; i < N; i++)
		Result[i] = new long double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			Result[i][j] = 0;
			for (int k = 0; k < N; k++)
				Result[i][j] += A[i][k] * B[k][j];
		}
	return Result;
}

long double** Diploma::Mult(long double a, long double** B, int N){
	long double** Result = new long double*[N];
	for (int i = 0; i < N; i++)
		Result[i] = new long double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Result[i][j] = a * B[i][j];

	return Result;
}

long double* Diploma::Mult(long double** A, long double* v, int N){
	long double* Result = new long double[N];
	long double temp = 0;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			temp += A[i][j] * v[j];
		Result[i] = temp;
		temp = 0;
	}
	return Result;
}

long double** Diploma::CGInversion(long double** a_src, int N){
	int	*	mask;
	long double	fmaxval;
	int		maxind;
	int		tmpi;
	long double	tmp;
	//long double	a[N][N];

	long double	**a;
	long double	**am;

	mask = new int[N];
	a = new long double*[N];
    am = new long double*[N];

	for (int i = 0; i < N; i++)
	{
		a[i] = new long double[N];
		for (int j = 0; j < N; j++)
		{
			a[i][j] = a_src[i][j];
		}
	}

	for (int i = 0; i < N; i++)
		am[i] = new long double[N];

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == j)
			{
				am[i][j] = 1.0;
			}
			else {
				am[i][j] = 0.0;
			}
		}
	}
	for (int i = 0; i < N; i++)
	{
		mask[i] = i;
	}
	for (int i = 0; i < N; i++)
	{
		maxind = i;
		fmaxval = fabs(a[i][i]);
		for (int ni = i + 1; ni < N; ni++)
		{
			if (fabs(fmaxval) <= fabs(a[ni][i]))
			{
				fmaxval = fabs(a[ni][i]);
				maxind = ni;
			}
		}
		fmaxval = a[maxind][i];
		if (fmaxval == 0)
		{
			return a_src;
		}
		if (i != maxind)
		{
			for (int nj = 0; nj < N; nj++)
			{
				tmp = a[i][nj];
				a[i][nj] = a[maxind][nj];
				a[maxind][nj] = tmp;

				tmp = am[i][nj];
				am[i][nj] = am[maxind][nj];
				am[maxind][nj] = tmp;
			}
			tmpi = mask[i];
			mask[i] = mask[maxind];
			mask[maxind] = tmpi;
		}
		long double aii = a[i][i];
		for (int j = 0; j < N; j++)
		{
			a[i][j] = a[i][j] / aii;
			am[i][j] = am[i][j] / aii;
		}
		for (int ni = 0; ni < N; ni++)
		{
			if (ni != i)
			{
				long double fconst = a[ni][i];
				for (int nj = 0; nj < N; nj++)
				{
					a[ni][nj] = a[ni][nj] - fconst *  a[i][nj];
					am[ni][nj] = am[ni][nj] - fconst * am[i][nj];
				}
			}
		}
	}
	/**/
	for (int i = 0; i < N; i++)
	{
		if (mask[i] != i)
		{
			for (int j = 0; j < N; j++)
			{
				tmp				= a[i][j];
				a[i][j]			= a[mask[i]][j];
				a[mask[i]][j]	= tmp;
			}
		}
	}
	/**/
	for (int i = 0; i < N; i++)
	{
		delete[] a[i];
	}
	delete[] a;
	delete[] mask;
	return am;
}

long double** Diploma::GInversion(long double **Ai){
	long double** B = new long double*[N];
	for (int i = 0; i < N; i++)
		B[i] = new long double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			B[i][j] = 0.0;
			if (i == j)
				B[i][j] = 1.0;
		}

	long double** A = new long double*[N];
	for (int i = 0; i < N; i++)
		A[i] = new long double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			A[i][j] = Ai[i][j];
		}

	long double a, b;
	for (int i = 0; i<N; i++)
	{
		a = A[i][i];
		for (int j = i + 1; j<N; j++)
		{
			b = A[j][i];
			for (int k = 0; k<N; k++)
			{
				A[j][k] = A[i][k] * b - A[j][k] * a;
				B[j][k] = B[i][k] * b - B[j][k] * a;
			}
		}
	}

	long double sum;
	for (int i = 0; i<N; i++)
	{
		for (int j = N - 1; j >= 0; j--)
		{
			sum = 0;
			for (int k = N - 1; k > j; k--)
			{
				sum += A[j][k] * B[k][i];
			}

			B[j][i] = B[j][i] / A[j][j] - sum / A[j][j];
		}
	}

	return B;
}

long double** Diploma::LUInversion(long double **Ai){

	long double sum;

	long double **mObr = new long double *[N];
	for (int i = 0; i < N; i++)
		mObr[i] = new long double[N];

	long double **LU = new long double *[N];//создаём массив под матрицу LU
	for (int i = 0; i < N; i++)
		LU[i] = new long double[N];

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j< N; j++)
		{
			sum = 0;
			if (i <= j)
			{
				for (int k = 0; k< i; k++)
					sum += LU[i][k] * LU[k][j];
				LU[i][j] = Ai[i][j] - sum;//вычисляем элементы верхней треугольной матрицы
			}
			else
			{
				for (int k = 0; k< j; k++)
					sum += LU[i][k] * LU[k][j];
				if (LU[j][j] == 0)
					return 0;
				LU[i][j] = (Ai[i][j] - sum) / LU[j][j];//вычисляем элементы нижней треугольной матрицы
			}
		}
	}


	for (int i = N - 1; i >= 0; i--)//нахождение обратной матрицы
	{
		for (int j = N - 1; j >= 0; j--)
		{
			sum = 0;
			if (i == j)
			{
				for (int p = j + 1; p< N; p++)
					sum += LU[j][p] * mObr[p][j];
				mObr[j][j] = (1 - sum) / LU[j][j];
			}
			else if (i< j)
			{
				for (int p = i + 1; p< N; p++)
					sum += LU[i][p] * mObr[p][j];
				mObr[i][j] = -sum / LU[i][i];
			}
			else
			{
				for (int p = j + 1; p< N; p++)
					sum += mObr[i][p] * LU[p][j];
				mObr[i][j] = -sum;
			}
		}
	}

	return mObr;
}

long double** Diploma::Jacobian(long double* y, int N){
	long double** J = new long double*[N];
	for (int i = 0; i < N; i++)
		J[i] = new long double[N];

	J[0][0] = -k1-4*k2*y[0];
	J[0][1] = 0;
	J[0][2] = 0;
	J[0][3] = 0;

	J[1][0] = k1+2*k2*y[0];
	J[1][1] = 0;
	J[1][2] = 0;
	J[1][3] = 0;

	J[2][0] = k1;
	J[2][1] = 0;
	J[2][2] = 0;
	J[2][3] = 0;

	J[3][0] = 4*k2*y[0];
	J[3][1] = 0;
	J[3][2] = 0;
	J[3][3] = 0;

	return J;
}

long double Diploma::norm(long double* x, long double* y, int N){
	long double temp = 0;

	for (int i = 0; i < N; i++)
		temp += pow(x[i] - y[i],2);
	return sqrt(temp);
}

long double** Y(long double t, int N){
    long double** Ym = new long double*[N];
	for (int i = 0; i < N; i++)
		Ym[i] = new long double[N];

    Ym[0][0] = exp(-k1*t);
    Ym[0][1] = 0;
    Ym[0][2] = 0;
    Ym[0][3]  = 0;

    Ym[1][0] = 1 - exp(-k1*t);
    Ym[1][1] = 1;
    Ym[1][2] = 0;
    Ym[1][3] = 0;

    Ym[2][0] = 1 - exp(-k1*t);
    Ym[2][1] = 0;
    Ym[2][2] = 1;
    Ym[2][3] = 0;

    Ym[3][0] = 0;
    Ym[3][1] = 0;
    Ym[3][2] = 0;
    Ym[3][3] = 1;

    return Ym;
}

long double* P(long double* x, int N){
    long double * Pv = new long double[N];

    Pv[0] = -2 * k2 * x[0] * x[0];
    Pv[1] = k2 * x[0] * x[0];
    Pv[2] = 0;
    Pv[3] = 2 * k2 * x[0] * x[0];

    return Pv;
}

long double* Diploma::LpTransformation(){


    //cout << numP;
    //show(y, numP);
    //show(y_, numP-1);
    long double *y0 = new long double[N];
    for (int i = 0; i < N; i++)
        y0[i] = y[0][i];
    long double t = 0;

    for( int i = 0; i < numP - 1; i++)
    {
        long double **TempYa = Y(-t, N);
        long double **TempYab = Y(-(t + H/2), N);
        long double **TempYb = Y(-(t + H), N);

        long double *TempPa = P(y[i], N);
        long double *TempPab = P(y_[i], N);
        long double *TempPb = P(y[i+1], N);

        long double *dY = Mult(
                               H/6.0,
                               Add(
                                   Add(
                                       Mult(TempYa, TempPa, N),
                                       Mult(4.0,
                                            Mult(TempYab, TempPab, N),
                                            N),
                                       N),
                                   Mult(TempYb, TempPb,N),
                                   N),
                               N);
        y0 = Add(y0, dY, N);

        t += H;
    }

    return y0;

}

void Diploma::CalcCoefs(){
    long double *A = new long double[2];
    A[0] = 1.08*pow(10,16);
    A[1] = 3.16*pow(10,16);

    long double *E = new long double[2];
    E[0] = 2.5*pow(10,5);
    E[1] = 2.7*pow(10,5);

    long double R = 8.31;

    k1 = A[0] * exp(-E[0]/(R*Temperature));
    k2 = A[1] * exp(-E[1]/(R*Temperature));

    cout << endl << k1 << " " << k2 << endl;

    double test;
    cin >> test;
}

void Diploma::mkMethod(long double Eps, long double h, long double B){

    numP = B/h + 1;
    H = h;
    B_ = B;
    y = new long double*[numP];
    y_ = new long double*[numP];
    for (int i = 0; i < numP; i++)
    {
        y[i] = new long double[N];
        y_[i] = new long double[N];
    }

    h = h/2;

    string sT = LDToStr(Temperature);
    string path = "resultMK(" + sT + ").txt";
    const char* cPath = path.c_str();
    ofstream fout(cPath);

	long double m, a, b, c, d;

	long double *yo = new long double[N];
	long double *v = new long double[N];
	long double *p = new long double[N];

	long double **k = new long double*[4];
	for (int j = 0; j < 4; j++)
		k[j] = new long double[N];

	long double** Dn = new long double*[N];
	for (int i = 0; i < N; i++)
		Dn[i] = new long double[N];

	long double** E = new long double*[N];
	for (int i = 0; i < N; i++)
		E[i] = new long double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;

	m = 0.57281606248213;
	a = 1.00900469029922;
	b = -0.25900469029921;
	c = -0.49552206416578;
	d = -1.28777648233922;

	p[0] = 1.27836939012447;
	p[1] = -1.00738680980438;
	p[2] = 0.92655391093950;
	p[3] = -0.33396131834691;

	yo[0] = 1; //.086;
	yo[1] = 0; //.0;
	yo[2] = 0; //.903;
	yo[3] = 0; //.011;

	for (int i = 0; i < N; i++)
        y[0][i] = yo[i];

	cout << "t = 0" << endl;
    show(fout, yo, N);

	long double time = h;
	int counter = 1;
	int l = 1;

	do
	{

		Dn = Add(E, Mult(-m*h, Jacobian(f(yo, N), N), N), N);
		Dn = LUInversion(Dn);

		k[0] = Mult(Dn, Mult(h, f(yo, N),N), N);
		k[1] = Mult(Dn, k[0], N);
		v = f(Add(yo, Add(Mult(a, k[0], N), Mult(b, k[1], N), N), N), N);
		k[2] = Mult(Dn, Add(Mult(h, v, N), Mult(c, k[1], N), N), N);
		k[3] = Mult(Dn, Add(k[2], Mult(d, k[1], N), N), N);


		for (int i = 0; i < 4; i++)
		{
			yo = Add(yo, Mult(p[i], k[i], N), N);
		}
/*
		double mMass = yo[0];

        for (int i = 1; i < N; i++)
            mMass += yo[i];

        for (int i = 0; i < N; i++)
            yo[i] /= mMass;
*/
		if(counter % 2)
            for (int i = 0; i < N; i++)
                y_[l-1][i] = yo[i];
        else
            {
                for (int i = 0; i < N; i++)
                y[l][i] = yo[i];
                l++;
                show(fout, yo, N);
            }

		cout << "t = " << time << endl;
		time+=h;
		counter++;

	} while (time <= B);
    cout << endl<< endl<< endl;
    show(LpTransformation(), N);
	fout.close();
	system("pause");
}

void Diploma::EuMethod(long double Eps, long double h, long double B){
    long double *yo = LpTransformation();
/*
    yo[0] = 0.78;
    yo[1] = 0.16;
    yo[2] = 0.1;
    yo[3] = 0.12;*/

    long double** A = new long double*[N];
	for (int i = 0; i < N; i++)
		A[i] = new long double[N];
    A[0][0] = -k1;
	A[0][1] = 0;
	A[0][2] = 0;
	A[0][3] = 0;

	A[1][0] = k1;
	A[1][1] = 0;
	A[1][2] = 0;
	A[1][3] = 0;

	A[2][0] = k1;
	A[2][1] = 0;
	A[2][2] = 0;
	A[2][3] = 0;

	A[3][0] = 0;
	A[3][1] = 0;
	A[3][2] = 0;
	A[3][3] = 0;

	string sT = LDToStr(Temperature);
    string path = "resultEU(" + sT + ").txt";
    const char* cPath = path.c_str();
    ofstream fout(cPath);

	cout << "t = 0" << endl;
    show(fout, yo, N);

	long double time = h;

    do
	{
        yo = Add(yo, Mult(A, yo, N), N);

		cout << "t = " << time << endl;
		time+=h;

		show(fout, yo, N);
	} while (time <= B);

    cout << endl<< endl<< endl;

	fout.close();
	system("pause");
}

