#include <iostream>
#include "intel.cpp"
using namespace std;
int main(int argc, char *argv[])
{
    double mat[] {1, 2, 3,
                    2, 5, 5,
                    3, 5, 2};
    double *temp{mat};
    double *up = (double *) calloc(9,sizeof(double));
    dagqrf(temp,3,3,up);
    return 0;
    
}