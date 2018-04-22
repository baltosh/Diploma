#pragma once

class Diploma
{
public:
    virtual long double* f(long double*, int); // правая часть системы дифференциальных уравнений
    virtual long double* Add(long double*, long double*, int); // сложение векторов
    virtual long double* Mult(long double, long double*, int); // умножение скаляра на вектор
    virtual long double** Add(long double**, long double**, int); // сложение матриц
    virtual long double** Mult(long double**, long double**, int); // умножение матрицы на матрицу
    virtual long double** Mult(long double, long double**, int); // умножение скаляра на матрицу
    virtual long double* Mult(long double**, long double*, int); // умножение матрицы на вектор
    virtual long double** CGInversion(long double**, int); // вычисление обратной матрицы
    virtual long double** GInversion(long double**); // вычисление обратной матрицы методом Гаусса
    virtual long double** LUInversion(long double**); // вычисление обратной матрицы с помощью LU разложения
    virtual long double** Jacobian(long double*, int); // расчет матрицы якоби
    virtual long double norm(long double*, long double*, int); // подсчет расстояния между векторами
    virtual void Method(long double, long double, long double);
    virtual long double* LpTransformation();

    ~Diploma();
protected:


};
