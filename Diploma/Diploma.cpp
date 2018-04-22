#include <iomanip>
#include <cstdio>
#include <iostream>
#include <Windows.h>
#include <cmath>
#include <fstream>
#include "Diploma.h"

using namespace std;

const int N = 4;

const long double k1 = 0.51;
const long double k2 = 0.07;

void show(long double** A, int N){
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			cout << setprecision(4) << setw(8) << A[i][j] << " ";
		cout << endl;
	}
	cout << endl;
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


void Diploma::Method(long double Eps)
{
    ofstream fout("result.txt");

	long double m, a, b, c, d;
	long double h = 0.01;

	long double *yo = new long double[N];
	long double *yn = new long double[N];
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

	yo[0] = .086;
	yo[1] = .0;
	yo[2] = .903;
	yo[3] = .011;
    show(fout, yo, N);
	p[0] = 1.27836939012447;
	p[1] = -1.00738680980438;
	p[2] = 0.92655391093950;
	p[3] = -0.33396131834691;
	int I = 0;

	cout << "CAAAMON" << endl;
	double g;
	cin >> g;
	do
	{
		yn = yo;

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
		cout << "t = " << I*h << endl;
		I++;
		show(fout, yo, N);
	} while (norm(yo, yn, N) > Eps);
	fout.close();
	system("pause");
}
