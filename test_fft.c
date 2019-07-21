/*
 * test_fft.c
 *
 *  Created on: Apr 17, 2019
 *      Author: echocare
 */



#include <fftw3.h>
#include "test.h"

float fft_in[ACTUAL_ROW_LEN];
void FFTReadLineFromFile()
{
	int j;

    for(j=0;j < (ACTUAL_ROW_LEN - 1) && !feof(file); j++)
	{
		if (fscanf(file,"%f,",&fft_in[j]) == EOF)
		{
			printf ("got EOF\n");
		}
	}
	if (fscanf(file,"%f/n",&fft_in[j]) == EOF)
	{
		printf ("got EOF\n");
	}
}
int test_fft()
{
	int i,k ;
	// int PaddingSize;
	char FilePath[MAX_FILE_PATH];
	char OutputFilePath[MAX_FILE_PATH];
	fftwf_complex out[ACTUAL_ROW_LEN];


	fftw_plan my_plan;
	int flags = 1;
	float fftin_local_ptr_complex[ACTUAL_ROW_LEN];

	memset(fftin_local_ptr_complex,0,sizeof(fftin_local_ptr_complex));

	sprintf (FilePath,"%s%s",INPUTS_BASE_DIR,INPUT_FFT_FILE);
	sprintf (OutputFilePath,"%s%s",OUTPUTS_BASE_DIR,OUTPUT_FFT_FILE);

	file = fopen (FilePath, "r");
	output_file = fopen (OutputFilePath, "w+");
	for(i = 0; i < ROW_COUNT; i++)
	{

		FFTReadLineFromFile();
		memset(out,0,sizeof(out));


		my_plan = fftwf_plan_dft_r2c_2d(ACTUAL_ROW_LEN,1,fftin_local_ptr_complex,out,flags);
		fftwf_execute_dft_r2c(my_plan,fft_in,out);
		for ( k = 0; k < ACTUAL_ROW_LEN; k++)
		{
//			printf("index: %d: %f %fi",k,out[k][0],out[k][1]);
			fprintf(output_file, "%f %fi,",out[k][0],out[k][1]);
		}
//		printf ("\n");
		fprintf(output_file, "\n");
	}
	fclose(file);
	fclose(output_file);
	return 0;

}

