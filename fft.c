#include "fft.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define PI 3.1415926535897932384626434

int power(int down, int up) {
    int ret = 1; 
    for(int i = 0; i < up; i++)
    {
       ret *= down; 
    }
    return ret; 
}

int inverse_index(int n, int L) {
    int ret = 0;    
    for(int i = 0; i < L; i++)
    {
        ret <<= 1;
        ret |= (n & 1);
        n >>= 1;
    }
    return ret;
}

void FFT_butterfly(complex* loc1, complex* loc2, int N, int k) {
    complex tmp = *loc1;       
    *loc2 = multiply(*loc2, conv_from_polar(1, -2*PI*k/N));
    *loc1 = add(*loc1, *loc2); // calculate loc1
    loc2->re = -(loc2->re);   // mul -1
    loc2->im = -(loc2->im);
    *loc2 = add(tmp, *loc2);  // calculate loc2
}

/* Implements radix 2 decimation of time FFT algorithm.
 * 
 * @input: pointer to a complex array  
 * @L: number of stages in data flow graph of FFT. Length of the array should equal to 2^N
 */
void FFT_radix2_decimation_in_time(complex* array, int L, int inverse) {
    // change the order of input array
    int inverse_multiplier = inverse?-1:1;
    int N = power(2, L);
    complex* tmp = (complex *)malloc(sizeof(struct complex_t) * N);
    memcpy(tmp, array, sizeof(struct complex_t) * N);
    for(int i = 0; i < N; i++)
        array[i] = tmp[inverse_index(i, L)];
    free(tmp);

    // FFT
    for(int stage = 0; stage < L; stage++)
    {
        int sep = power(2, stage+1);
        for(int j = 0; j < N/sep; j++)
        {
            for(int k = 0; k < sep/2; k++)
                FFT_butterfly(array+j*sep+k, array+j*sep+sep/2+k, N, inverse_multiplier*k*N/sep);
        }
    }

    if (inverse) 
    {
        for(int i = 0; i < N; i++)
        {
            array[i].re /= N;
            array[i].im /= N;
        }
    }
}


/*
 *  Only support square matrix (rows = cols)
 * 
 */
void FFT2D_radix2_decimation_in_time(complex* array[], int L, int inverse) {
    int N = power(2, L);
    for (int i = 0; i < N; i++)
    {
        FFT_radix2_decimation_in_time(array[i], L, inverse);
    }

    complex *tmp = (complex*) malloc(sizeof(struct complex_t) * N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            tmp[j] = array[j][i];
        FFT_radix2_decimation_in_time(tmp, L, inverse);
        for (int j = 0; j < N; j++)
            array[j][i] = tmp[j]; 
    }
    free(tmp);
}