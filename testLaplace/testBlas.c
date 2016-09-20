#include <stdio.h>
#include <blas.h>

double m[] = {
3, 1, 3,
1, 5, 9,
2, 6, 5
};

double b[] = {
3, 1, 3,
1, 5, 9,
2, 6, 5
};

double e[] = {
0, 0, 0,
0, 0, 0,
0, 0, 0
};

int
main()
{
int i, j;

for (i=0; i<3; ++i) {
for (j=0; j<3; ++j) printf("%5.1f", m[i*3+j]);
putchar('\n');
}

cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, m, 3, b, 3, 0.0, e, 3);
putchar('\n');
for (i=0; i<3; ++i) {
for (j=0; j<3; ++j) printf("%5.1f", e[i*3+j]);
putchar('\n');
}

return 0;
}
