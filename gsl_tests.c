/*
 * gsl_median.c
 *
 *  Created on: Mar 27, 2019
 *      Author: echocare
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_double.h>

#include "test.h"

//General Median definitions
#define MEDIAN_LEN 155
#define WINDOW_SIZE 20

// General Cov definitions
#define COL_AMOUNT 200
#define COL_LEN 288

#define print_matrix(M)	\
	for(i = 0; i < M->size1; i++) \
	{ \
		for(j = 0; j < (M->size2 - 1); j++) \
		{ \
			printf("%lf\t",gsl_matrix_get(M, i, j)); \
		} \
		printf("%lf\n",gsl_matrix_get(M, i, j)); \
	}

double data[MEDIAN_LEN];

int cov_calculate(gsl_matrix *r, gsl_matrix *m)
{
	gsl_vector_view a, b;
	size_t i, j;
	for (i = 0; i < m->size1; i++)
	{
		for (j = 0; j < m->size1; j++)
		{
			double v;
			a = gsl_matrix_column (m, i);
			b = gsl_matrix_column (m, j);
			v = gsl_stats_covariance (a.vector.data, a.vector.stride,
			b.vector.data, b.vector.stride, a.vector.size);
			gsl_matrix_set (r, i, j, v);
		}
	}

	return 0;
}

int covariance_matrix_main()
{
	FILE *file;
	FILE *output_file;
	int i, j;
	double input;
	char FilePath[MAX_FILE_PATH];

	gsl_matrix * input_matrix = gsl_matrix_alloc(COL_AMOUNT, COL_LEN);
	gsl_matrix * output_matrix = gsl_matrix_alloc(COL_AMOUNT, COL_AMOUNT);
	sprintf (FilePath,"%s/%s",INPUTS_BASE_DIR,INPUT_COV_FILE);
	file = fopen (FilePath, "r");
	output_file = fopen(FilePath, "w+");
	for(i = 0; i < input_matrix->size1; i++)
	{
		for(j = 0; j < (input_matrix->size2 - 1) && !feof(file); j++)
		{
			fscanf(file,"%lf,",&input);
			gsl_matrix_set(input_matrix, i, j, input);
			//printf("%lf\t",gsl_matrix_get(input_matrix, i, j));
		}
		fscanf(file,"%lf\n,",&input);
		gsl_matrix_set(input_matrix, i, j, input);
		//printf("%lf\n",gsl_matrix_get(input_matrix, i, j));
	}

	cov_calculate(output_matrix, input_matrix);
	sprintf (FilePath,"%s/%s",OUTPUTS_BASE_DIR,OUTPUT_COV_FILE);

	for(i = 0; i < output_matrix->size1; i++)
	{
		for(j = 0; j < (output_matrix->size2 - 1); j++)
		{
			fprintf(output_file, "%lf,",gsl_matrix_get(output_matrix, i, j));
		}
		fprintf(output_file, "%lf\n",gsl_matrix_get(output_matrix, i, j));
	}

	gsl_matrix_free(input_matrix);
	gsl_matrix_free(output_matrix);

	return 0;
}

int median_main(char *InFile, char *OutFile)
{
	FILE *file = NULL;
	FILE * output_file = NULL;
	double myout[MEDIAN_LEN];
	int j;
//	gsl_vector * input = gsl_vector_alloc(MEDIAN_LEN + WINDOW_SIZE);
//	gsl_vector * output = gsl_vector_alloc(MEDIAN_LEN + WINDOW_SIZE);
	gsl_vector * input = gsl_vector_alloc(MEDIAN_LEN);
	gsl_vector * output = gsl_vector_alloc(MEDIAN_LEN);
	gsl_movstat_workspace * w = gsl_movstat_alloc(WINDOW_SIZE+1);

	// sprintf (FilePath,"%s/%s",INPUTS_BASE_DIR,INPUT_MEDIAN1_FILE);
	file = fopen (InFile, "r");
/*
	for(j = 0; j < WINDOW_SIZE / 2  ; j++)
	{
		gsl_vector_set(input, j, 0);
		printf("%lf\n", gsl_vector_get(input, j));
	}
*/
//	for(; j < (MEDIAN_LEN + (WINDOW_SIZE / 2) ) && !feof(file); j++)
	for(j = 0 ; j < (MEDIAN_LEN) && !feof(file); j++)
	{
		fscanf(file,"%lf\n",&data[j]);
		gsl_vector_set(input, j, data[j]);
		printf("%lf\n", gsl_vector_get(input, j));
	}
/*
	for(; j < (MEDIAN_LEN + WINDOW_SIZE); j++)
	{
		gsl_vector_set(input, j, 0);
		printf("%lf\n", gsl_vector_get(input, j));
	}
*/
	fclose(file);

	printf("Now outputing results:\n");


// 	gsl_movstat_median(GSL_MOVSTAT_END_TRUNCATE, input, output, w);
	gsl_movstat_median(GSL_MOVSTAT_END_PADZERO, input, output, w);
	myout[0] = gsl_vector_get(output, 0) / 2;
	printf("value 0 : %f\n",myout[0]);
	for(j = 1 ; j < (MEDIAN_LEN ) ; j++)
	{
		if ((j&0x1) == 0)
			myout[j] =  gsl_vector_get(output, j);
		else
			myout[j] = (gsl_vector_get(output, j) + gsl_vector_get(output, j-1)) / 2;
		printf("value %d : %f\n",j,myout[j] );
	}



//	sprintf (FilePath,"%s/%s",OUTPUTS_BASE_DIR,OUTPUT_MEDIAN1_FILE);
	output_file = fopen (OutFile, "w+");

	for(j = 0 ; j < (MEDIAN_LEN ) && !feof(file); j++)
	{
		fprintf(output_file, "%lf\n",myout[j]);
	}
	fprintf(output_file, "\n");

	fclose (output_file);
	return 0;
}


