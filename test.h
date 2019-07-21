/*
 * test.h
 *
 *  Created on: Apr 3, 2019
 *      Author: echocare
 */

#ifndef TEST_H_
#define TEST_H_


//FIR-definitions
#define PERFORMANCE_TEST
#define ENABLE_NE10_FIR_FLOAT_NEON
// #define ENABLE_NE10_FIR_DECIMATE_FLOAT_NEON
//#define ENABLE_NE10_FIR_INTERPOLATE_FLOAT_NEON
//#define ENABLE_NE10_FIR_LATTICE_FLOAT_NEON
//#define ENABLE_NE10_FIR_SPARSE_FLOAT_NEON
#define NUMTAPS   19
#define NUMBLOCKS  1
#define BLOCKSIZE 288
#define FIR_GAIN  0.59926


//General-definitions - FFT
#define ACTUAL_ROW_LEN 288
#define ROW_COUNT 800
// with padding
#define ROW_LEN 512



#define TEST_LENGTH_SAMPLES (ROW_LEN + 32)
#define TEST_COUNT ROW_COUNT
#define NUM_PERF_TESTS ROW_COUNT
#define MAX_NUMTAPS   NUMTAPS
#define MAX_BLOCKSIZE BLOCKSIZE


#ifdef __ARM
#include <NE10.h>
extern ne10_float32_t fir_coeffs [NUMTAPS];
#endif

#ifdef __ARM
#define INPUTS_BASE_DIR "/home/debian/inputs/"
#define OUTPUTS_BASE_DIR "/home/debian/outputs/"
#else
/* x86 */
#define INPUTS_BASE_DIR "/home/echocare/Desktop/winshare/"
#define OUTPUTS_BASE_DIR "/home/echocare/Desktop/winshare/"
#endif

#define MAX_FILE_PATH 100

#define INPUT_FFT_FILE "input_FIR_FFT.csv"
#define INPUT_FIR_FILE "input_FIR_FFT.csv"
#define INPUT_MEDIAN1_FILE "input_MEDIAN_1.csv"
#define INPUT_MEDIAN2_FILE "input_MEDIAN_2.csv"
#define INPUT_COV_FILE "input_COV.csv"

#define OUTPUT_FFT_FILE "output_FFT.csv"
#define OUTPUT_FIR_FILE "output_FIR.csv"
#define OUTPUT_MEDIAN1_FILE "output_MEDIAN_1.csv"
#define OUTPUT_MEDIAN2_FILE "output_MEDIAN_2.csv"
#define OUTPUT_COV_FILE "output_COV.csv"



//covar_opencv();

void test_fft_r2c_1d_float32_performance();

int covariance_matrix_main();
int median_main();

void test_fir();
void test_fir_decimate();
void test_fir_interpolate();
void test_fir_lattice();
void test_fir_sparse();

void my_fir_test_setup (void);
void my_test_setup();
void ReadLineFromFile();

void opencv_tests();
int test_fft();

extern FILE *file;
extern FILE *output_file;


#endif /* TEST_H_ */
