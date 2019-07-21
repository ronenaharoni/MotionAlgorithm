/*
 *  Copyright 2012-16 ARM Limited and Contributors.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of ARM Limited nor the
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY ARM LIMITED AND CONTRIBUTORS "AS IS" AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL ARM LIMITED AND CONTRIBUTORS BE LIABLE FOR ANY
 *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * NE10 Library : test_suite_fir.c
 */
#ifdef __ARM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "test.h"

#include "NE10_dsp.h"
#include "seatest.h"




/* ----------------------------------------------------------------------
** Global defines
** ------------------------------------------------------------------- */

/* Max data Length and block size, numTaps */
/*
#define TEST_LENGTH_SAMPLES 320
#define MAX_BLOCKSIZE 320
#define MAX_NUMTAPS 100
*/
/*
#define TEST_LENGTH_SAMPLES 288
#define MAX_NUMTAPS   20
#define NUMBLOCKS 2
#define MAX_BLOCKSIZE 8
*/


// #define TEST_COUNT 800

//input and output
static ne10_float32_t * guarded_in_c = NULL;
static ne10_float32_t * guarded_in_neon = NULL;
static ne10_float32_t * in_c = NULL;
static ne10_float32_t * in_neon = NULL;

static ne10_float32_t * guarded_out_c = NULL;
static ne10_float32_t * guarded_out_neon = NULL;
static ne10_float32_t * out_c = NULL;
static ne10_float32_t * out_neon = NULL;

static ne10_float32_t * guarded_fir_state_c = NULL;
static ne10_float32_t * guarded_fir_state_neon = NULL;
static ne10_float32_t * fir_state_c = NULL;
static ne10_float32_t * fir_state_neon = NULL;

#if defined (SMOKE_TEST)||(REGRESSION_TEST)
static ne10_float32_t snr = 0.0f;
#endif
#ifdef PERFORMANCE_TEST
static ne10_int64_t time_c = 0;
static ne10_int64_t time_neon = 0;
static ne10_float32_t time_speedup = 0.0f;
static ne10_float32_t time_savings = 0.0f;
#endif


void FirReadLineFromFile();

static ne10_float32_t testInput_f32[TEST_LENGTH_SAMPLES];

/* ----------------------------------------------------------------------
** Defines each of the tests performed
** ------------------------------------------------------------------- */
typedef struct
{
    ne10_uint32_t blockSize;
    ne10_uint32_t numTaps;
    ne10_uint32_t numFrames;
    ne10_float32_t *coeffsF32;
    ne10_float32_t *inputF32;
} test_config;

/* Test configurationsfor conformance test, 100% Code Coverage */
#if defined (SMOKE_TEST)||(REGRESSION_TEST)
static test_config CONFIG[] =
{
    {64, 32, 5, &testCoeffs32_f32[0], &testInput_f32[0]},
    {64, 3, 5, &testCoeffs3_f32[0], &testInput_f32[0]},
    {64, 7, 5, &testCoeffs7_f32[0], &testInput_f32[0]},
    {64, 1, 5, &testCoeffs1_f32, &testInput_f32[0]},
    {5, 3, 64, &testCoeffs3_f32[0], &testInput_f32[0]},
    {2, 7, 160, &testCoeffs7_f32[0], &testInput_f32[0]},
    {4, 1, 80, &testCoeffs1_f32, &testInput_f32[0]},
    {32, 32, 10, &testCoeffs32_f32[0], &testInput_f32[0]},
    {7, 4, 1, &testCoeffs4_f32[0], &testInput_f32[0]},
    {8, 4, 1, &testCoeffs4_f32[0], &testInput_f32[0]},
    {9, 4, 1, &testCoeffs4_f32[0], &testInput_f32[0]},
};
#define NUM_TESTS (sizeof(CONFIG) / sizeof(CONFIG[0]) )
#endif
/* Test configurations for performance test */
#ifdef PERFORMANCE_TEST
#if 0
static test_config CONFIG_PERF[] =
{
    {64, 32, 5, &testCoeffs32_f32[0], &testInput_f32[0]},
    {64, 3, 5, &testCoeffs3_f32[0], &testInput_f32[0]},
    {64, 7, 5, &testCoeffs7_f32[0], &testInput_f32[0]},
};
#define NUM_PERF_TESTS (sizeof(CONFIG_PERF) / sizeof(CONFIG_PERF[0]) )
#endif
static test_config CONFIG_PERF[] = { BLOCKSIZE, NUMTAPS, NUMBLOCKS, &fir_coeffs[0],&testInput_f32[0] };
#endif
unsigned int SaveFirIndex = 0;

void test_fir_case0()
{

    ne10_fir_instance_f32_t SC, SN;

    ne10_uint16_t loop = 0;
    ne10_uint16_t block = 0;
    ne10_uint16_t i = 0;

    test_config *config;

    char FilePath[MAX_FILE_PATH];

    fprintf (stdout, "----------%30s start\n", __FUNCTION__);

    /* init input memory */
    NE10_SRC_ALLOC (in_c, guarded_in_c, TEST_LENGTH_SAMPLES ); // 16 extra bytes at the begining and 16 extra bytes at the end
    NE10_SRC_ALLOC (in_neon, guarded_in_neon, TEST_LENGTH_SAMPLES); // 16 extra bytes at the begining and 16 extra bytes at the end

    /* init dst memory */
    NE10_DST_ALLOC (out_c, guarded_out_c, TEST_LENGTH_SAMPLES);
    NE10_DST_ALLOC (out_neon, guarded_out_neon, TEST_LENGTH_SAMPLES);

    /* init state memory */
    NE10_DST_ALLOC (fir_state_c, guarded_fir_state_c, MAX_NUMTAPS + MAX_BLOCKSIZE);
    NE10_DST_ALLOC (fir_state_neon, guarded_fir_state_neon, MAX_NUMTAPS + MAX_BLOCKSIZE);

    /*
     * @TODO Separate neon and c version test. Mixing of these 2 tests makes
     * it difficult to disable and enable one of them. Current macro
     * ENABLE_NE10_FIR_FLOAT_NEON disables both of them, which should only
     * disable the neon version ideally.
     */
#ifdef ENABLE_NE10_FIR_FLOAT_NEON
#if defined (SMOKE_TEST)||(REGRESSION_TEST)
    ne10_uint16_t pos = 0;
    for (loop = 0; loop < NUM_TESTS; loop++)
    {
        config = &CONFIG[loop];

        /* Initialize the CFFT/CIFFT module */
        ne10_fir_init_float (&SC, config->numTaps, config->coeffsF32, fir_state_c, config->blockSize);
        ne10_fir_init_float (&SN, config->numTaps, config->coeffsF32, fir_state_neon, config->blockSize);

        /* copy input to input buffer */
        for (i = 0; i < TEST_LENGTH_SAMPLES; i++)
        {
            in_c[i] = testInput_f32[i];
            in_neon[i] = testInput_f32[i];
        }

        GUARD_ARRAY (out_c, TEST_LENGTH_SAMPLES);
        GUARD_ARRAY (out_neon, TEST_LENGTH_SAMPLES);

        for (block = 0; block < config->numFrames; block++)
        {
            ne10_fir_float_c (&SC, in_c + (block * config->blockSize), out_c + (block * config->blockSize), config->blockSize);
        }

        for (block = 0; block < config->numFrames; block++)
        {
            ne10_fir_float_neon (&SN, in_neon + (block * config->blockSize), out_neon + (block * config->blockSize), config->blockSize);
        }

        assert_true (CHECK_ARRAY_GUARD (out_c, TEST_LENGTH_SAMPLES));
        assert_true (CHECK_ARRAY_GUARD (out_neon, TEST_LENGTH_SAMPLES));

        //conformance test 1: compare snr
        snr = CAL_SNR_FLOAT32 (out_c, out_neon, TEST_LENGTH_SAMPLES);
        assert_false ( (snr < SNR_THRESHOLD));

        //conformance test 2: compare output of C and neon
#if defined (DEBUG_TRACE)
        printf ("--------------------config %d\n", loop);
        printf ("snr %f\n", snr);
#endif
        for (pos = 0; pos < TEST_LENGTH_SAMPLES; pos++)
        {
#if defined (DEBUG_TRACE)
            printf ("pos %d \n", pos);
            printf ("c %f (0x%04X) neon %f (0x%04X)\n", out_c[pos], * (ne10_uint32_t*) &out_c[pos], out_neon[pos], * (ne10_uint32_t*) &out_neon[pos]);
#endif
            assert_float_vec_equal (&out_c[pos], &out_neon[pos], ERROR_MARGIN_SMALL, 1);
        }
    }
#endif
#endif // ENABLE_NE10_FIR_FLOAT_NEON

#ifdef PERFORMANCE_TEST
    ne10_uint16_t k;
    fprintf (stdout, "%25s%20s%20s%20s%20s\n", "FIR Length&Taps", "C Time (micro-s)", "NEON Time (micro-s)", "Time Savings", "Performance Ratio");
    for (loop = 0; loop < NUM_PERF_TESTS; loop++)
    {
       // config = &CONFIG_PERF[loop];
    	config = &CONFIG_PERF[0];
    	FirReadLineFromFile();
    	if (loop == 312)
		{
			printf ("just for debug\n");
		}
        /* Initialize the CFFT/CIFFT module */
        ne10_fir_init_float (&SC, config->numTaps, config->coeffsF32, fir_state_c, config->blockSize);
        ne10_fir_init_float (&SN, config->numTaps, config->coeffsF32, fir_state_neon, config->blockSize);

        /* copy input to input buffer */
        for (i = 0; i < TEST_LENGTH_SAMPLES; i++)
        {
            in_c[i] = testInput_f32[i];
            in_neon[i] = testInput_f32[i];
        }

        GET_TIME
        (
            time_c,
        {
 //           for (k = 0; k < TEST_COUNT; k++)
 //           {
                for (block = 0; block < config->numFrames; block++)
                {
                    ne10_fir_float_c (&SC, in_c + (block * config->blockSize), out_c + (block * config->blockSize), config->blockSize);
                }
 //           }
        }
        );

#ifdef ENABLE_NE10_FIR_FLOAT_NEON
        GET_TIME
        (
            time_neon,
        {
//            for (k = 0; k < TEST_COUNT; k++)
//            {
                for (block = 0; block < config->numFrames; block++)
                {
                    ne10_fir_float_neon (&SN, in_neon + (block * config->blockSize), out_neon + (block * config->blockSize), config->blockSize);
                }
 //           }
        }
        );
#endif // ENABLE_NE10_FIR_FLOAT_NEON
        SaveFirIndex = loop;
        time_speedup = (ne10_float32_t) time_c / time_neon;
        time_savings = ( ( (ne10_float32_t) (time_c - time_neon)) / time_c) * 100;
        ne10_log (__FUNCTION__, "%20d,%4d%20lld%20lld%19.2f%%%18.2f:1\n", config->numTaps, time_c, time_neon, time_savings, time_speedup);

    	for(k = 0 ; k < ACTUAL_ROW_LEN  && !feof(file); k++)
    	{
    		// fprintf(output_file, "%lf\n", gsl_vector_get(output, j));
    		fprintf(output_file, "%lf,",out_neon[k]*FIR_GAIN);
    		// print results
    	    //printf("vector num %d: %lf\n",k,out_neon[k]*FIR_GAIN);

    	}
    	fprintf(output_file, "\n");
    }

#endif

	fclose (output_file);

    free (guarded_in_c);
    free (guarded_in_neon);
    free (guarded_out_c);
    free (guarded_out_neon);
    free (guarded_fir_state_c);
    free (guarded_fir_state_neon);
    fprintf (stdout, "----------%30s end\n", __FUNCTION__);
}

void test_fir()
{
	my_fir_test_setup();
    test_fir_case0();
}

/* ----------------------------------------------------------------------
** end of fir test
** ------------------------------------------------------------------- */
#if 0
static void my_test_setup (void)
{
    ne10_log_buffer_ptr = ne10_log_buffer;
}
#endif


void FirReadLineFromFile()
{
    ne10_int32_t j;

//	for(;j < (ACTUAL_ROW_LEN + NUMTAPS - 2) && !feof(file); j++)
    for(j=0;j < (ACTUAL_ROW_LEN - 1) && !feof(file); j++)
	{
		if (fscanf(file,"%f,",&testInput_f32[j]) == EOF)
		{
			printf ("got EOF\n");
		}
		//fft_input[i] = (ne10_float32_t)rand() / RAND_MAX * 50.0f;
	}
	if (fscanf(file,"%f/n",&testInput_f32[j]) == EOF)
	{
		printf ("got EOF\n");
	}
}


void my_fir_test_setup (void)
{
	ne10_log_buffer_ptr = &ne10_log_buffer[0];

	char FilePath[MAX_FILE_PATH];

	sprintf (FilePath,"%s/%s",INPUTS_BASE_DIR,INPUT_FIR_FILE);
	file = fopen (FilePath, "r");

    // write results to file
	sprintf (FilePath,"%s/%s",OUTPUTS_BASE_DIR,OUTPUT_FIR_FILE);
	output_file = fopen (FilePath, "w+");
}

void test_fixture_fir (void)
{
    test_fixture_start();               // starts a fixture

    fixture_setup (my_test_setup);

    run_test (test_fir);       // run tests

    test_fixture_end();                 // ends a fixture
}
#endif

