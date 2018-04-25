#pragma once

class Diploma
{
public:
    long double Temperature;

    virtual long double* f(long double*, int); // ������ ����� ������� ���������������� ���������
    virtual long double* Add(long double*, long double*, int); // �������� ��������
    virtual long double* Mult(long double, long double*, int); // ��������� ������� �� ������
    virtual long double** Add(long double**, long double**, int); // �������� ������
    virtual long double** Mult(long double**, long double**, int); // ��������� ������� �� �������
    virtual long double** Mult(long double, long double**, int); // ��������� ������� �� �������
    virtual long double* Mult(long double**, long double*, int); // ��������� ������� �� ������
    virtual long double** CGInversion(long double**, int); // ���������� �������� �������
    virtual long double** GInversion(long double**); // ���������� �������� ������� ������� ������
    virtual long double** LUInversion(long double**); // ���������� �������� ������� � ������� LU ����������
    virtual long double** Jacobian(long double*, int); // ������ ������� �����
    virtual long double norm(long double*, long double*, int); // ������� ���������� ����� ���������
    virtual void mkMethod(long double, long double, long double);
    virtual void EuMethod(long double, long double, long double);
    virtual long double* LpTransformation();
    virtual void CalcCoefs();
    ~Diploma();
protected:


};
