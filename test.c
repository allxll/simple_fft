#include "fft.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define L 18
#define n 512
#define N (n*n) 

int main(void) {

    FILE* fin = fopen("./noise.txt", "r");
    FILE* fout = fopen("./result/noise_freq.txt", "w");
    FILE* fout2 = fopen("./result/noise_filtered.txt", "w");
    if (fin == NULL || fout == NULL) {
        printf("null\n");
        return 1;
    }
    
    // Load matrix
    complex *array[n]; 
    for(int i = 0; i < n; i++)  
    {
        array[i] = (complex*) malloc(sizeof(struct complex_t) * n);
        memset(array[i], 0, sizeof(struct complex_t) * n);
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            fscanf(fin, "%lf", &(array[i][j].re));
        }
    }

    // FFT
    FFT2D_radix2_decimation_in_time(array, L/2, 0);

    // Save FFT result
    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fprintf(fout, "%d ", (int)sqrt(array[i][j].re*array[i][j].re+array[i][j].im*array[i][j].im));
        }
        fprintf(fout, "\n");
    }

    // Gaussian Filter
    
    int col, row;
    for(int i = 0; i < n; i++)
    {
        row = n/2-abs(i - n/2);
        // printf("*** row: %d\n", row);
        for(int j = 0; j < n; j++)
        {
            col = n/2-abs(j - n/2);
            // printf("col: %d\n", col);
            array[i][j].re *= exp(-(col*col+row*row)/(2*30*30));
            array[i][j].im *= exp(-(col*col+row*row)/(2*30*30));
        }
    }
     
    

    // IDFT
    FFT2D_radix2_decimation_in_time(array, L/2, 1);

    // Save filtered matrix
    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fprintf(fout2, "%d ", (int)sqrt(array[i][j].re*array[i][j].re+array[i][j].im*array[i][j].im));
        }
        fprintf(fout2, "\n");
    }
    
    for(int i = 0; i < n; i++)
        free(array[i]); 
    
    fclose(fin);
    fclose(fout);
    fclose(fout2);
    return 0;
}

