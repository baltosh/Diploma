#include "Diploma.h"

using namespace std;

const int N = 4;

const long double Eps = .000001;


int main()
{
    Diploma* n = new Diploma;
    n->Temperature = 800;
    n->CalcCoefs();
	n->mkMethod(Eps, .01, 10);
    n->EuMethod(Eps, .01, 10);
	return 0;
}
