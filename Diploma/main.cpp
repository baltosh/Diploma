#include "Diploma.h"

using namespace std;

const int N = 4;

const long double Eps = .000001;


int main()
{
    Diploma* n = new Diploma;
	n->Method(Eps, 1.0, 5);

	return 0;
}
