/*
 ============================================================================
 Name        : testneon.c
 Author      : Tzvi
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */




#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "test.h"

FILE *file;
FILE *output_file;

#ifdef __ARM

#include <NE10.h>

#define NE10_ENABLE_DSP
#define NE10_ENABLE_MATH



int real_fft(ne10_float32_t src[], int length, ne10_fft_cpx_float32_t dst [])
{

	if(length != ROW_LEN)
	{
		return 1;
	}
    ne10_fft_r2c_cfg_float32_t cfg;                     // An FFT "configuration structure"

    // Initialise Ne10, using hardware auto-detection to set library function pointers
    if (ne10_init() != NE10_OK)
    {
        fprintf(stderr, "Failed to initialise Ne10.\n");
        return 1;
    }

    // Prepare the real-to-complex single precision floating point FFT configuration
    // structure for inputs of length `SAMPLES`. (You need only generate this once for a
    // particular input size.)
    cfg = ne10_fft_alloc_r2c_float32(length);

    // Perform the FFT
    ne10_fft_r2c_1d_float32(dst, src, cfg);

    // Display the results
    for (int i = 0; i < length; i++)
    {
        //printf( "IN[%2d]: %10.4f\t", i, src[i]);
        if (i <= length / 2)
            printf("OUT[%2d]: %10.4f + %10.4fi", i, dst[i].r, dst[i].i);
        printf("\n");
    }

    // Free the allocated configuration structure
    ne10_fft_destroy_r2c_float32(cfg);

    return 0;
}

int fir(ne10_float32_t src [], int length, ne10_float32_t coeffs [], ne10_float32_t gain)
{
/*
	if(length != ROW_LEN)
	{
		return 1;
	}
*/
    ne10_float32_t dst[ROW_LEN]; // A destination array for the transformed data
    ne10_float32_t st[NUMTAPS + BLOCKSIZE - 1]; // A "state" buffer for use within the FIR
    ne10_fir_instance_f32_t cfg;    // An FIR "instance structure"

    // Initialize Ne10, using hardware auto-detection to set library function pointers
//    if (ne10_init() != NE10_OK)
//    {
//        fprintf(stderr, "Failed to initialise Ne10.\n");
//        return 1;
//    }

    // Prepare the FIR instance structure, storing `NUMTAPS`, `coeffs`, and `st` within
    // it, and clearing the state buffer. (For constant parameters, this process can
    // instead be performed manually.)
    if (ne10_fir_init_float(&cfg, NUMTAPS, coeffs, st, BLOCKSIZE) != NE10_OK)
    {
        fprintf(stderr, "Failed to initialize FIR instance structure.\n");
        return 1;
    }

    // Perform the FIR filtering of the input buffer in `NUMBLOCKS` blocks of `BLOCKSIZE`
    // elements using the parameters set up in the FIR instance structure `cfg`.
    for (int b = 0; b < NUMBLOCKS; b++)
    {
        ne10_fir_float(&cfg, src + (b * BLOCKSIZE), dst + (b * BLOCKSIZE), BLOCKSIZE);
    }

    // multiplying the result values by the gain constant
    for(int i = 0; i < ACTUAL_ROW_LEN; i++)
    {
    	dst[i] = gain * dst[i];
    }

    // Display the results (dst[i] = b[0] * src[i] + b[1] * src[i - 1] + b[2] * src[i - 2]
    //                               + ... + b[NUMTAPS - 1] * src[i - (NUMTAPS - 1)])
    printf("Coefficients:\n");
    for (int i = NUMTAPS - 1; i >= 0; i--)
    {
        printf("\tb[%d] = %5.4f\n", NUMTAPS - (i + 1), coeffs[i]);
    }
    for (int i = 0; i < ACTUAL_ROW_LEN; i++)
    {
        //printf( "IN[%2d]: %9.4f\t", i, src[i]);
        printf("OUT[%2d]: %9.4f\n", i, dst[i]);
    }

    return 0;
}

// ne10_float32_t fft_input [ROW_LEN] = {};
ne10_float32_t fft_input [ROW_LEN]={};
ne10_fft_cpx_float32_t fft_output[ROW_LEN] = {};
ne10_float32_t fir_coeffs2 [NUMTAPS] = {0.000408,-0.041049,0.0086493,0.22937,-0.056319,-0.52033,0.16501,0.89238,-0.067669,-0.9717,0.15692,0.71792,-0.1289,-0.34998,0.091564,0.090303,-0.022514,-0.013823,0.003896};
ne10_float32_t fir_coeffs[NUMTAPS];

#endif

int main888(void)
{


	char FilePath[MAX_FILE_PATH];
	char OutFilePath[MAX_FILE_PATH];

#ifdef __ARM
	int i;
	for ( i = 0 ; i < NUMTAPS; i++)
	{
		fir_coeffs[i] = fir_coeffs2[NUMTAPS-i-1];
	}
#endif

//	covar_opencv();
//	fir(fft_input,ROW_LEN,fir_coeffs,FIR_GAIN);
//covariance_matrix_main();
//	test_fft_r2c_1d_float32_performance();

/*****  Working implementations *****/
#ifdef __ARM
	// from libNE10
	test_fir();
#endif
	// from opencv library - covariance


//	test_fft();//me add



	opencv_tests();
	/* from gsl library */
	sprintf (FilePath,"%s/%s",INPUTS_BASE_DIR,INPUT_MEDIAN1_FILE);
	sprintf (OutFilePath,"%s/%s",OUTPUTS_BASE_DIR,OUTPUT_MEDIAN1_FILE);
	median_main(FilePath,OutFilePath);
//	sprintf (FilePath,"%s/%s",INPUTS_BASE_DIR,INPUT_MEDIAN2_FILE);
//	sprintf (OutFilePath,"%s/%s",OUTPUTS_BASE_DIR,OUTPUT_MEDIAN2_FILE);
//	median_main(FilePath,OutFilePath);
	// from FFTW library
	test_fft();

	return 0;
}
//

