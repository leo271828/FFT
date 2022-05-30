#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>


complex F[11][11];
void FFT(complex*, complex*, int);

int main(){
    time_t T;
    int k, j, N = 1200*20; // N = 2^p * 3^q * 5^r
    double theta;
    complex *z, *u;
    z = malloc( N*sizeof(complex) ); 
    u = malloc( N*sizeof(complex) );

    theta = 2.0 * M_PI / N; // Omega_n = cos(theta) + i sin(theta)
    for(k=0;k<N;k++){
        u[k] = (double)k;
        z[k] = 0;
    }

    // Compute F matrix
    F[2][0] = 1, F[2][1] = -1;
    for(int i=0;i<3;i++) F[3][i] = cexp(-I*2*i*M_PI/3);
    for(int i=0;i<5;i++) F[5][i] = cexp(-I*2*i*M_PI/5);

    T = clock();
    printf("Start FFT with N : %d ...\n\n", N);

    FFT(z, u, N);
    printf("You can check last 3 values are \n");
    for(int i=N-3;i<N;i++){
        printf("%d : %.2f + %.2f \n", i, z[i], z[i]);
    }
    printf("\nTotal time : %d\n", clock() - T);

    free(z); free(u); 
    return 0;
}

void FFT(complex* y, complex *x, int N){
    // Part 1
    int base;
    if( N % 2 == 0 ) base = 2;
    else if ( N % 3 == 0 ) base = 3;
    else if ( N % 5 == 0 ) base = 5;
    else {
        base = N;
        for(int i=0;i<5;i++) F[5][i] = cexp(-I*2*i*M_PI/5);
    }

    if( N == base ){
        for(int i=0;i<base;i++){
            y[i] = 0;
            for(int j=0;j<base;j++)
                y[i] += x[j] * F[base][(i*j)%base];
        }
        return ;
    }
    complex *Nx = malloc( N * sizeof(complex) );
    complex *Ny = malloc( N * sizeof(complex) );
    if( Nx == NULL || Ny == NULL ) printf("None");
    for(int k=0;k<N/base;k++){
        for(int j=0;j<base;j++)
            Nx[k + j*N/base] = x[k*base + j];
    }

    // Part 2
    for(int k=0;k<base;k++)
        FFT(Ny + k*N/base, Nx + k*N/base, N/base);

    // Part 3
    complex wn = 1;
    complex w = cexp(-I*2*M_PI/N);
    for(int k=0;k<N/base;k++){
        for(int j=0;j<base;j++)
            Ny[k + j*N/base] *= cpow(wn, j);
        for(int j=0;j<base;j++){
            y[k+j*N/base] = 0;
            for(int i=0;i<base;i++)
                y[k+j*N/base] += F[base][i*j%base]*Ny[k+i*N/base];
        }
        wn *= w;
    }
    free(Ny); free(Nx);
}
