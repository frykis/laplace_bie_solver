#include <stdio.h>
#include <cblas.h>

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

double f[] = {
1, 2, 3,
4, 5, 6,
7, 8, 9
};

double g[] = {
0, 0, 0,
};

double identity[] = {
1, 2, 3,
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
for (i=0; i<3; ++i) {
for (j=0; j<3; ++j) printf("%5.1f", m[i*3+j]);
putchar('\n');
}
//C←αAB + βC
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1.0, m, 3, identity, 1, 0.0, g, 1);
putchar('\n');






for (j=0; j<3; ++j){ printf("%5.1f", g[j]);
putchar('\n');
}


return 0;
}
