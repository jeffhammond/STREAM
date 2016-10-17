#include <stdio.h>
#include <stdlib.h>
static int cmpdouble(const void * ptr1, const void * ptr2)
{
    const double * dp1 = (const double*) ptr1;
    const double * dp2 = (const double*) ptr2;
    if      (*dp1 < *dp2) return -1;
    else if (*dp1 > *dp2) return  1;
    else                  return  0;
}

int main(int argc, char* argv[])
{
    int n = (argc>1) ? atoi(argv[1]) : 10;
    double * x = malloc(n*sizeof(double));

    srand(n);
    for (int i=0; i<n; i++) {
        x[i] = (double)rand()/(double)RAND_MAX;
    }

    printf("BEFORE\n");
    for (int i=0; i<n; i++) {
        printf("x[%d]=%lf\n", i, x[i]);
    }

    qsort(x, n, sizeof(double), &cmpdouble);

    printf("AFTER\n");
    for (int i=0; i<n; i++) {
        printf("x[%d]=%lf\n", i, x[i]);
    }

    free(x);

    return 0;
}
