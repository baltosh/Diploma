#include "Diploma.h"

using namespace std;

const int N = 4;

const long double Eps = .000001;


int main()
{
    Diploma* n = new Diploma;
	n->mkMethod(Eps, .001, 10);
    n->EuMethod(Eps, .001, 10);
	return 0;
}
