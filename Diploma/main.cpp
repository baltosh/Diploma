#include "Diploma.h"

using namespace std;

const int N = 4;

const long double Eps = .00000001;


int main()
{
    Diploma* n = new Diploma;
	n->Method(Eps);

	return 0;
}
