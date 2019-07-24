#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
#define I _Complex_I (0.0f+1.0if)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "MotionHeader.h"
#include <math.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_statistics_double.h>

#include "polyfit.h"
#include <fftw3.h>
#include <NE10.h>
#include "seatest.h"
#include "NE10_dsp.h"
#include <sys/time.h>
#ifndef NULL
#define NULL   ((void *) 0)
#endif

struct timeval tpStart, tpStop;
struct timeval tp;
float f1 = 0;


int MotionPreprocess(float* Mscan[],float* Mscan_PostProcessReal[],_Complex float* Mscan_PostProcessHilbert[],FgParams_Struct* FgParams){

	//PreProcessing of Mscan
	MatchedFilter(Mscan, FgParams);
	SlowProcessing(Mscan, Mscan_PostProcessReal, FgParams); //remove DC
	NotchFilter(Mscan_PostProcessReal, FgParams); //filter the 50&100Hz disturbance that occurred from the fluorescent lamp


	Hilbert(Mscan, Mscan_PostProcessHilbert, FgParams);	//get the IQ signal with hilbert transform
	SlowProcessingHilbert(Mscan_PostProcessHilbert, Mscan_PostProcessHilbert, FgParams); //remove DC
	NotchFilterHilbert(Mscan_PostProcessHilbert, FgParams); //filter the 50&100Hz disturbance that occurred from the fluorescent lamp

	return 0;
}



int MotionAnalyzerTwoMscans(Motion_Struct* MotionStruct0,Motion_Struct* MotionStruct1,Motion_Struct* MotionStruct2,FgParams_Struct* FgParams,Tree_Struct** All_Trees,SVM_Struct* SVM_Model,RF_Struct* RF_Model,int *y_hat_M,float* MotionDistribution){
	int i;
	Motion_Struct UnitedMotionStruct;
	Edge2_Struct Edge2_United;
	Edge2_Struct Edge2_Plus_United;
	Edge2_Struct Edge2_Minus_United;
	Features_Struct Plus, Minus, Both;
	AllFeatures_Struct FeatureSet;

	Event_Struct EventStruct_United;

	FeatureSet.Plus = &Plus;
	FeatureSet.Minus = &Minus;
	FeatureSet.Both = &Both;

	///////Create UnitedMotionStruct for the gap interpolation/////////
	UnitedMotionStruct.Edge2 = &Edge2_United;
	UnitedMotionStruct.Edge2_Plus = &Edge2_Plus_United;
	UnitedMotionStruct.Edge2_Minus = &Edge2_Minus_United;
	UnitedMotionStruct.EventStruct = &EventStruct_United;
	UnitedMotionStruct.EventPlusPassEvent = 0;
	UnitedMotionStruct.EventMinusPassEvent = 0;

	//////////////////////Gap Interpolation/////////////////

	GapInterpolation(MotionStruct0, MotionStruct1, MotionStruct2,
			&UnitedMotionStruct, FgParams);


	//////////////////////Feature Extraction/////////////////
	if (FgParams->Motion.FirstTimeMotion == 1) {
		FgParams->Motion.SpectrogramTimeBins =
				FgParams->Motion.SpectrogramTimeBinsTwoMotions;

		FeatureExtractionBasedCurves(&UnitedMotionStruct, &FeatureSet,
				FgParams);

	} else {
		FgParams->Motion.SpectrogramTimeBins =
				FgParams->Motion.SpectrogramTimeBinsThreeMotions;
		FeatureExtractionBasedCurves(&UnitedMotionStruct, &FeatureSet,
				FgParams);

	}

	/////////Classification///////
	if (UnitedMotionStruct.EventPlusPassEvent
			|| UnitedMotionStruct.EventMinusPassEvent) {//there was some event
		*y_hat_M = RandomForrestClassifier(&FeatureSet, All_Trees,
				MotionDistribution, RF_Model);
		//		printf("y_hat %d\n", y_hat_M);
		ClassifierCorrection(&FeatureSet, MotionDistribution, y_hat_M,
				SVM_Model);

		if ((MotionDistribution[2] < 0.625) && (*y_hat_M == 3)) {
			*y_hat_M = 1;
			printf("Changed to normal motion from fall, not fully sure\n");
			MotionDistribution[0] = 1;
			MotionDistribution[1] = 0;
			MotionDistribution[2] = 0;
			MotionDistribution[3] = 0;
		}

	} else {	//there is no event, set quasistatic automatically
		*y_hat_M = 4;
		MotionDistribution[3] = 1;	//MotionDistribution=[0,0,0,1]
	}



	int SaveResults = 1;
	if (SaveResults == 1 && FgParams->Motion.FirstTimeMotion == 0) {

		SaveToCsv(*y_hat_M, MotionDistribution, &FeatureSet, RF_Model,
				FgParams);
	}


	/////////Save MotionStruct2 into MotionStruc0 for the next loop

	memcpy(MotionStruct0->Edge2->Peak_Filtered,
			MotionStruct2->Edge2->Peak_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));

	memcpy(MotionStruct0->Edge2_Plus->Fmax, MotionStruct2->Edge2_Plus->Fmax,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->FiftyPrecent_Filtered,
			MotionStruct2->Edge2_Plus->FiftyPrecent_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->Peak_Filtered,
			MotionStruct2->Edge2_Plus->Peak_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->SumEnergy_Post,
			MotionStruct2->Edge2_Plus->SumEnergy_Post,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->T1_t, MotionStruct2->Edge2_Plus->T1_t,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));

	memcpy(MotionStruct0->Edge2_Minus->Fmax, MotionStruct2->Edge2_Minus->Fmax,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->FiftyPrecent_Filtered,
			MotionStruct2->Edge2_Minus->FiftyPrecent_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->Peak_Filtered,
			MotionStruct2->Edge2_Minus->Peak_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->SumEnergy_Post,
			MotionStruct2->Edge2_Minus->SumEnergy_Post,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->T1_t, MotionStruct2->Edge2_Minus->T1_t,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));

	//save the last AvgValue-1 bins of FiftyPrecent for the initial state of the AvgFilter OF motionstruct1, the inital is of motionstruct2 that becomes motionstruct0 at the next loop
	for (i = FgParams->Motion.AvgValue - 1; i > 0; i--) {	//it's done only to Plus and Minus
		MotionStruct1->Edge2_Plus->PrevLastFiftyPrecent[(FgParams->Motion.AvgValue - 1)
														- i] =
																MotionStruct2->Edge2_Plus->FiftyPrecent[FgParams->Motion.SpectrogramTimeBinsSingleMotion
																									   - i];
		MotionStruct1->Edge2_Minus->PrevLastFiftyPrecent[(FgParams->Motion.AvgValue
				- 1) - i] =
						MotionStruct2->Edge2_Minus->FiftyPrecent[FgParams->Motion.SpectrogramTimeBinsSingleMotion
																- i];
	}

	return 0;
}


int MotionAnalyzer(float* Mscan1[],float* Mscan2[],float* Mscan_PostProcess1[],float* Mscan_PostProcess2[],FgParams_Struct* FgParams,Tree_Struct** All_Trees,SVM_Struct* SVM_Model,RF_Struct* RF_Model,Motion_Struct* MotionStruct0,BgRadarParams *bgParams,int *y_hat_M,float* MotionDistribution)
{
	int i;
	float* Mscan_abs_FFT[FgParams->Motion.DFTLengthForPSD];
	float* Pxx2_dB[FgParams->Motion.SpectrogramFreqBins],
	*Pxx2[FgParams->Motion.SpectrogramFreqBins];
	float * Pxx2_Hilbert[FgParams->Motion.SpectrogramFreqBinsHilbert];

	Motion_Struct MotionStruct1;
	Motion_Struct MotionStruct2;
	Motion_Struct UnitedMotionStruct;

	Edge2_Struct Edge2_1;
	Edge2_Struct Edge2_Plus_1;
	Edge2_Struct Edge2_Minus_1;

	Edge2_Struct Edge2_2;
	Edge2_Struct Edge2_Plus_2;
	Edge2_Struct Edge2_Minus_2;

	Edge2_Struct Edge2_United;
	Edge2_Struct Edge2_Plus_United;
	Edge2_Struct Edge2_Minus_United;

	Features_Struct Plus, Minus, Both;
	AllFeatures_Struct FeatureSet;

	Event_Struct EventStruct_United;

	FeatureSet.Plus = &Plus;
	FeatureSet.Minus = &Minus;
	FeatureSet.Both = &Both;

	//////////////////////////Calculation of MotionStruct1//////////////////////////

	FgParams->Motion.SpectrogramTimeBins = FgParams->Motion.SpectrogramTimeBinsSingleMotion; //first there is extraction on single motion (=75)

	for (i = 0; i < FgParams->Motion.DFTLengthForPSD; i++) {
		Mscan_abs_FFT[i] = (float *) malloc(FgParams->Motion.Nbins * sizeof(float));
	}
	for (i = 0; i < FgParams->Motion.SpectrogramFreqBins; i++) {
		Pxx2_dB[i] = (float *) malloc(
				FgParams->Motion.SpectrogramTimeBins * sizeof(float));
		Pxx2[i] = (float *) malloc(
				FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	}
	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert; i++) {

		Pxx2_Hilbert[i] = (float *) malloc(
				FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	}


	//Memory Allocation for the curves

	MemoryAllocation(&Edge2_1, &Edge2_2, &Edge2_Plus_1, &Edge2_Minus_1,
			&Edge2_Plus_2, &Edge2_Minus_2, FgParams);

	if (FgParams->Motion.FirstTimeMotion == 1) { //if it's first time there is no history and therefore set all 0 for the initial state of the AvgFilter
		memset(Edge2_1.PrevLastFiftyPrecent, 0,
				(FgParams->Motion.AvgValue - 1) * sizeof(float));
		memset(Edge2_Plus_1.PrevLastFiftyPrecent, 0,
				(FgParams->Motion.AvgValue - 1) * sizeof(float));
		memset(Edge2_Minus_1.PrevLastFiftyPrecent, 0,
				(FgParams->Motion.AvgValue - 1) * sizeof(float));
	} else {
		memcpy(Edge2_Plus_1.PrevLastFiftyPrecent,
				MotionStruct0->Edge2_Plus->PrevLastFiftyPrecent,
				(FgParams->Motion.AvgValue - 1) * sizeof(float));
		memcpy(Edge2_Minus_1.PrevLastFiftyPrecent,
				MotionStruct0->Edge2_Minus->PrevLastFiftyPrecent,
				(FgParams->Motion.AvgValue - 1) * sizeof(float));
		//notice that PrevLastFiftyPrecent is actually the last 50precent of MotionStruct0;
	}
	//		gettimeofday(&tpStop,0); f1 = ( (float)( tpStop.tv_sec-tpStart.tv_sec)+ (float)(tpStop.tv_usec)/1000000 ) -  ((float)(tpStart.tv_usec)/1000000); printf (" %f sec\n", f1 );

	//Extraction of the curves
	MotionCurveExtraction(&MotionStruct1, Mscan1, Mscan_abs_FFT,Mscan_PostProcess1, Pxx2_Hilbert,
			Pxx2, Pxx2_dB, &Edge2_1, &Edge2_Plus_1, &Edge2_Minus_1, FgParams);



	//////////////////////////Calculation of MotionStruct2//////////////////////////

	//save the last AvgValue-1 bins of FiftyPrecent for the initial state of the AvgFilter

	for (i = FgParams->Motion.AvgValue - 1; i > 0; i--) { //it's done only to Plus and Minus
		Edge2_Plus_2.PrevLastFiftyPrecent[(FgParams->Motion.AvgValue - 1) - i] =
				Edge2_Plus_1.FiftyPrecent[FgParams->Motion.SpectrogramTimeBins - i];
		Edge2_Minus_2.PrevLastFiftyPrecent[(FgParams->Motion.AvgValue - 1) - i] =
				Edge2_Minus_1.FiftyPrecent[FgParams->Motion.SpectrogramTimeBins - i];
	}

	//Extraction of the curves
	MotionCurveExtraction(&MotionStruct2, Mscan2, Mscan_abs_FFT,Mscan_PostProcess2, Pxx2_Hilbert,
			Pxx2, Pxx2_dB, &Edge2_2, &Edge2_Plus_2, &Edge2_Minus_2, FgParams);







	///////Create UnitedMotionStruct for the gap interpolation/////////
	UnitedMotionStruct.Edge2 = &Edge2_United;
	UnitedMotionStruct.Edge2_Plus = &Edge2_Plus_United;
	UnitedMotionStruct.Edge2_Minus = &Edge2_Minus_United;
	UnitedMotionStruct.EventStruct = &EventStruct_United;
	UnitedMotionStruct.EventPlusPassEvent = 0;
	UnitedMotionStruct.EventMinusPassEvent = 0;

	//////////////////////Gap Interpolation/////////////////

	GapInterpolation(MotionStruct0, &MotionStruct1, &MotionStruct2,
			&UnitedMotionStruct, FgParams);


	//////////////////////Feature Extraction/////////////////
	if (FgParams->Motion.FirstTimeMotion == 1) {
		FgParams->Motion.SpectrogramTimeBins =
				FgParams->Motion.SpectrogramTimeBinsTwoMotions;

		FeatureExtractionBasedCurves(&UnitedMotionStruct, &FeatureSet,
				FgParams);

	} else {
		FgParams->Motion.SpectrogramTimeBins =
				FgParams->Motion.SpectrogramTimeBinsThreeMotions;
		FeatureExtractionBasedCurves(&UnitedMotionStruct, &FeatureSet,
				FgParams);

	}

	/////////Classification///////
	if (UnitedMotionStruct.EventPlusPassEvent
			|| UnitedMotionStruct.EventMinusPassEvent) {//there was some event
		*y_hat_M = RandomForrestClassifier(&FeatureSet, All_Trees,
				MotionDistribution, RF_Model);
		//		printf("y_hat %d\n", y_hat_M);
		ClassifierCorrection(&FeatureSet, MotionDistribution, y_hat_M,
				SVM_Model);

		if ((MotionDistribution[2] < 0.625) && (*y_hat_M == 3)) {
			*y_hat_M = 1;
			printf("Changed to normal motion from fall, not fully sure\n");
			MotionDistribution[0] = 1;
			MotionDistribution[1] = 0;
			MotionDistribution[2] = 0;
			MotionDistribution[3] = 0;
		}

	} else {	//there is no event, set quasistatic automatically
		*y_hat_M = 4;
		MotionDistribution[3] = 1;	//MotionDistribution=[0,0,0,1]
	}



	int SaveResults = 0;
	if (SaveResults == 1 && FgParams->Motion.FirstTimeMotion == 0) {

		SaveToCsv(*y_hat_M, MotionDistribution, &FeatureSet, RF_Model,
				FgParams);
	}

	/////////Save MotionStruct2 into MotionStruc0 for the next loop
	//	memcpy(MotionStruct0->Edge2->Fmax, MotionStruct2.Edge2->Fmax,
	//			FgParams->Motion.Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	//	memcpy(MotionStruct0->Edge2->FiftyPrecent_Filtered,
	//			MotionStruct2.Edge2->FiftyPrecent_Filtered,
	//			FgParams->Motion.Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2->Peak_Filtered,
			MotionStruct2.Edge2->Peak_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));

	memcpy(MotionStruct0->Edge2_Plus->Fmax, MotionStruct2.Edge2_Plus->Fmax,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->FiftyPrecent_Filtered,
			MotionStruct2.Edge2_Plus->FiftyPrecent_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->Peak_Filtered,
			MotionStruct2.Edge2_Plus->Peak_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->SumEnergy_Post,
			MotionStruct2.Edge2_Plus->SumEnergy_Post,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Plus->T1_t, MotionStruct2.Edge2_Plus->T1_t,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));

	memcpy(MotionStruct0->Edge2_Minus->Fmax, MotionStruct2.Edge2_Minus->Fmax,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->FiftyPrecent_Filtered,
			MotionStruct2.Edge2_Minus->FiftyPrecent_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->Peak_Filtered,
			MotionStruct2.Edge2_Minus->Peak_Filtered,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->SumEnergy_Post,
			MotionStruct2.Edge2_Minus->SumEnergy_Post,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	memcpy(MotionStruct0->Edge2_Minus->T1_t, MotionStruct2.Edge2_Minus->T1_t,
			FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));

	//save the last AvgValue-1 bins of FiftyPrecent for the initial state of the AvgFilter
	for (i = FgParams->Motion.AvgValue - 1; i > 0; i--) {	//it's done only to Plus and Minus
		MotionStruct0->Edge2_Plus->PrevLastFiftyPrecent[(FgParams->Motion.AvgValue - 1)
														- i] =
																MotionStruct2.Edge2_Plus->FiftyPrecent[FgParams->Motion.SpectrogramTimeBinsSingleMotion
																									   - i];
		MotionStruct0->Edge2_Minus->PrevLastFiftyPrecent[(FgParams->Motion.AvgValue
				- 1) - i] =
						MotionStruct2.Edge2_Minus->FiftyPrecent[FgParams->Motion.SpectrogramTimeBinsSingleMotion
																- i];
	}

	//FREE MEMORY

	for (i = 0; i < FgParams->Motion.DFTLengthForPSD; i++) {
		free(Mscan_abs_FFT[i]);
	}
	for (i = 0; i < FgParams->Motion.SpectrogramFreqBins; i++) {
		free(Pxx2_dB[i]);
		free(Pxx2[i]);
	}
	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert; i++) {

		free(Pxx2_Hilbert[i]);
	}

	FreeMemory(&Edge2_1, &Edge2_2, &Edge2_Plus_1, &Edge2_Minus_1, &Edge2_Plus_2,
			&Edge2_Minus_2, &UnitedMotionStruct);

	return 0;
}

int GapInterpolation(Motion_Struct *MotionStruct0,
		Motion_Struct *MotionStruct1, Motion_Struct *MotionStruct2,
		Motion_Struct *UnitedMotionStruct, FgParams_Struct* FgParams) {
	int i, FirstGap;
	if (FgParams->Motion.FirstTimeMotion == 1) {	//there are only 2 motion structs
		FirstGap = 1;
		//Memory allocation for the new curves with length SpectrogramTimeBinsTwoMotions=2*SpectrogramTimeBins + GapLength
		//		UnitedMotionStruct->Edge2->Fmax = (float *) malloc(
		//				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->Fmax = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->Fmax = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));

		UnitedMotionStruct->Edge2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));

		//		UnitedMotionStruct->Edge2->FiftyPrecent_Filtered = (float *) malloc(
		//				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered =
				(float *) malloc(
						(FgParams->Motion.SpectrogramTimeBinsTwoMotions)
						* sizeof(float));
		UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered =
				(float *) malloc(
						(FgParams->Motion.SpectrogramTimeBinsTwoMotions)
						* sizeof(float));

		UnitedMotionStruct->Edge2_Plus->SumEnergy_Post = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->SumEnergy_Post = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));

		UnitedMotionStruct->Edge2_Plus->T1_t = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->T1_t = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsTwoMotions) * sizeof(float));

		//		GapInterpolation_2Curves(MotionStruct1->Edge2->Fmax,
		//				MotionStruct2->Edge2->Fmax, UnitedMotionStruct->Edge2->Fmax,
		//				FgParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Plus->Fmax,
				MotionStruct2->Edge2_Plus->Fmax,
				UnitedMotionStruct->Edge2_Plus->Fmax, FgParams, FirstGap,1);//1 Becasue isHilbert=1
		GapInterpolation_2Curves(MotionStruct1->Edge2_Minus->Fmax,
				MotionStruct2->Edge2_Minus->Fmax,
				UnitedMotionStruct->Edge2_Minus->Fmax, FgParams, FirstGap,1);

		GapInterpolation_2Curves(MotionStruct1->Edge2->Peak_Filtered,
				MotionStruct2->Edge2->Peak_Filtered,
				UnitedMotionStruct->Edge2->Peak_Filtered, FgParams, FirstGap,0);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Plus->Peak_Filtered,
				MotionStruct2->Edge2_Plus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Plus->Peak_Filtered, FgParams,
				FirstGap,1);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Minus->Peak_Filtered,
				MotionStruct2->Edge2_Minus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Minus->Peak_Filtered, FgParams,
				FirstGap,1);

		//		GapInterpolation_2Curves(MotionStruct1->Edge2->FiftyPrecent_Filtered,
		//				MotionStruct2->Edge2->FiftyPrecent_Filtered,
		//				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered, FgParams,
		//				FirstGap);
		GapInterpolation_2Curves(
				MotionStruct1->Edge2_Plus->FiftyPrecent_Filtered,
				MotionStruct2->Edge2_Plus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered,
				FgParams, FirstGap,1);
		GapInterpolation_2Curves(
				MotionStruct1->Edge2_Minus->FiftyPrecent_Filtered,
				MotionStruct2->Edge2_Minus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered,
				FgParams, FirstGap,1);

		GapInterpolation_2Curves(MotionStruct1->Edge2_Plus->SumEnergy_Post,
				MotionStruct2->Edge2_Plus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Plus->SumEnergy_Post, FgParams,
				FirstGap,1);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Minus->SumEnergy_Post,
				MotionStruct2->Edge2_Minus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Minus->SumEnergy_Post, FgParams,
				FirstGap,1);

		GapInterpolation_2Curves(MotionStruct1->Edge2_Plus->T1_t,
				MotionStruct2->Edge2_Plus->T1_t,
				UnitedMotionStruct->Edge2_Plus->T1_t, FgParams, FirstGap,1);
		GapInterpolation_2Curves(MotionStruct1->Edge2_Minus->T1_t,
				MotionStruct2->Edge2_Minus->T1_t,
				UnitedMotionStruct->Edge2_Minus->T1_t, FgParams, FirstGap,1);

		FgParams->Motion.SpectrogramTimeBins =
				FgParams->Motion.SpectrogramTimeBinsTwoMotions;//now the SpectrogramTimebins is longer

		for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {//FiftyPrecent_Filtered corrections, 50precent curve cannot be bigger than the Fmax
			//			if (UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i]
			//																 > UnitedMotionStruct->Edge2->Fmax[i])
			//				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i] =
			//						UnitedMotionStruct->Edge2->Fmax[i];

			if (UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																	  > UnitedMotionStruct->Edge2_Plus->Fmax[i])
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Plus->Fmax[i];

			if (UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]
																	   > UnitedMotionStruct->Edge2_Minus->Fmax[i])
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Minus->Fmax[i];

		}

	} else {	//there are 3 motion structs
		FirstGap = 1;
		//Memory allocation for the new curves with length SpectrogramTimeBinsThreeMotions=3*SpectrogramTimeBins + 2*GapLength
		//		UnitedMotionStruct->Edge2->Fmax = (float *) malloc(
		//				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->Fmax = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->Fmax = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));

		UnitedMotionStruct->Edge2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));

		//		UnitedMotionStruct->Edge2->FiftyPrecent_Filtered = (float *) malloc(
		//				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered =
				(float *) malloc(
						(FgParams->Motion.SpectrogramTimeBinsThreeMotions)
						* sizeof(float));
		UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered =
				(float *) malloc(
						(FgParams->Motion.SpectrogramTimeBinsThreeMotions)
						* sizeof(float));

		UnitedMotionStruct->Edge2_Plus->SumEnergy_Post = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->SumEnergy_Post = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));

		UnitedMotionStruct->Edge2_Plus->T1_t = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));
		UnitedMotionStruct->Edge2_Minus->T1_t = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBinsThreeMotions) * sizeof(float));

		//implement gap interpolation between the first two motion structs (MostionStruct0&MostionStruct1)
		//		GapInterpolation_2Curves(MotionStruct0->Edge2->Fmax,
		//				MotionStruct1->Edge2->Fmax, UnitedMotionStruct->Edge2->Fmax,
		//				FgParams, FirstGap);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Plus->Fmax,
				MotionStruct1->Edge2_Plus->Fmax,
				UnitedMotionStruct->Edge2_Plus->Fmax, FgParams, FirstGap,1);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Minus->Fmax,
				MotionStruct1->Edge2_Minus->Fmax,
				UnitedMotionStruct->Edge2_Minus->Fmax, FgParams, FirstGap,1);

		GapInterpolation_2Curves(MotionStruct0->Edge2->Peak_Filtered,
				MotionStruct1->Edge2->Peak_Filtered,
				UnitedMotionStruct->Edge2->Peak_Filtered, FgParams, FirstGap,0);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Plus->Peak_Filtered,
				MotionStruct1->Edge2_Plus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Plus->Peak_Filtered, FgParams,
				FirstGap,1);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Minus->Peak_Filtered,
				MotionStruct1->Edge2_Minus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Minus->Peak_Filtered, FgParams,
				FirstGap,1);

		//		GapInterpolation_2Curves(MotionStruct0->Edge2->FiftyPrecent_Filtered,
		//				MotionStruct1->Edge2->FiftyPrecent_Filtered,
		//				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered, FgParams,
		//				FirstGap);
		GapInterpolation_2Curves(
				MotionStruct0->Edge2_Plus->FiftyPrecent_Filtered,
				MotionStruct1->Edge2_Plus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered,
				FgParams, FirstGap,1);
		GapInterpolation_2Curves(
				MotionStruct0->Edge2_Minus->FiftyPrecent_Filtered,
				MotionStruct1->Edge2_Minus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered,
				FgParams, FirstGap,1);

		GapInterpolation_2Curves(MotionStruct0->Edge2_Plus->SumEnergy_Post,
				MotionStruct1->Edge2_Plus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Plus->SumEnergy_Post, FgParams,
				FirstGap,1);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Minus->SumEnergy_Post,
				MotionStruct1->Edge2_Minus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Minus->SumEnergy_Post, FgParams,
				FirstGap,1);

		GapInterpolation_2Curves(MotionStruct0->Edge2_Plus->T1_t,
				MotionStruct1->Edge2_Plus->T1_t,
				UnitedMotionStruct->Edge2_Plus->T1_t, FgParams, FirstGap,1);
		GapInterpolation_2Curves(MotionStruct0->Edge2_Minus->T1_t,
				MotionStruct1->Edge2_Minus->T1_t,
				UnitedMotionStruct->Edge2_Minus->T1_t, FgParams, FirstGap,1);

		for (i = 0; i < FgParams->Motion.SpectrogramTimeBinsTwoMotions; i++) {//FiftyPrecent_Filtered corrections, 50precent curve cannot be bigger than the Fmax
			//			if (UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i]
			//																 > UnitedMotionStruct->Edge2->Fmax[i])
			//				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i] =
			//						UnitedMotionStruct->Edge2->Fmax[i];

			if (UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																	  > UnitedMotionStruct->Edge2_Plus->Fmax[i])
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Plus->Fmax[i];

			if (UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]
																	   > UnitedMotionStruct->Edge2_Minus->Fmax[i])
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Minus->Fmax[i];

		}

		//now take the result and implement interpolation with the last MostionStruct2
		FirstGap = 0;
		//		GapInterpolation_2Curves(UnitedMotionStruct->Edge2->Fmax,
		//				MotionStruct2->Edge2->Fmax, UnitedMotionStruct->Edge2->Fmax,
		//				FgParams, FirstGap);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Plus->Fmax,
				MotionStruct2->Edge2_Plus->Fmax,
				UnitedMotionStruct->Edge2_Plus->Fmax, FgParams, FirstGap,1);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Minus->Fmax,
				MotionStruct2->Edge2_Minus->Fmax,
				UnitedMotionStruct->Edge2_Minus->Fmax, FgParams, FirstGap,1);

		GapInterpolation_2Curves(UnitedMotionStruct->Edge2->Peak_Filtered,
				MotionStruct2->Edge2->Peak_Filtered,
				UnitedMotionStruct->Edge2->Peak_Filtered, FgParams, FirstGap,0);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Plus->Peak_Filtered,
				MotionStruct2->Edge2_Plus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Plus->Peak_Filtered, FgParams,
				FirstGap,1);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Minus->Peak_Filtered,
				MotionStruct2->Edge2_Minus->Peak_Filtered,
				UnitedMotionStruct->Edge2_Minus->Peak_Filtered, FgParams,
				FirstGap,1);

		//		GapInterpolation_2Curves(
		//				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered,
		//				MotionStruct2->Edge2->FiftyPrecent_Filtered,
		//				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered, FgParams,
		//				FirstGap);
		GapInterpolation_2Curves(
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered,
				MotionStruct2->Edge2_Plus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered,
				FgParams, FirstGap,1);
		GapInterpolation_2Curves(
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered,
				MotionStruct2->Edge2_Minus->FiftyPrecent_Filtered,
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered,
				FgParams, FirstGap,1);

		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Plus->SumEnergy_Post,
				MotionStruct2->Edge2_Plus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Plus->SumEnergy_Post, FgParams,
				FirstGap,1);
		GapInterpolation_2Curves(
				UnitedMotionStruct->Edge2_Minus->SumEnergy_Post,
				MotionStruct2->Edge2_Minus->SumEnergy_Post,
				UnitedMotionStruct->Edge2_Minus->SumEnergy_Post, FgParams,
				FirstGap,1);

		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Plus->T1_t,
				MotionStruct2->Edge2_Plus->T1_t,
				UnitedMotionStruct->Edge2_Plus->T1_t, FgParams, FirstGap,1);
		GapInterpolation_2Curves(UnitedMotionStruct->Edge2_Minus->T1_t,
				MotionStruct2->Edge2_Minus->T1_t,
				UnitedMotionStruct->Edge2_Minus->T1_t, FgParams, FirstGap,1);

		for (i = 0; i < FgParams->Motion.SpectrogramTimeBinsThreeMotions; i++) {//FiftyPrecent_Filtered corrections, 50precent curve cannot be bigger than the Fmax
			//			if (UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i]
			//																 > UnitedMotionStruct->Edge2->Fmax[i])
			//				UnitedMotionStruct->Edge2->FiftyPrecent_Filtered[i] =
			//						UnitedMotionStruct->Edge2->Fmax[i];

			if (UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																	  > UnitedMotionStruct->Edge2_Plus->Fmax[i])
				UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Plus->Fmax[i];

			if (UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]
																	   > UnitedMotionStruct->Edge2_Minus->Fmax[i])
				UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i] =
						UnitedMotionStruct->Edge2_Minus->Fmax[i];

//printf("%d %f\n",i,UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]);
		}
	}

	return 0;
}

int GapInterpolation_2Curves(float* LeftCurve, float* RightCurve,
		float *IntrpolatedCurve, FgParams_Struct* FgParams, int FirstGap,int isHilbert) {
	// we want to find a polynom of 3'rd order y = a3*x^3 + a2*x^2 + a1*x + a0 that best interpolate 2 curves,
	// There are 4 constrains -> and therefore 4 coeff are estimated:
	// 1. p1: continuity in left  curve: the polynom should get the same value in x1 as in curve1(x1):
	// 2. p2: continuity in right curve: the polynom should get the same value in x2 as in curve2(x2):
	//3. p3: same slope/derivative as in curve1 from left:
	// 4. p4: same slope/derivative as in curve2 from right:
	//p1 = a3*x^3 + a2*x^2 + a1*x + a0 @ x1
	// p2 = a3*x^3 + a2*x^2 + a1*x + a0 @ x2
	// p3 = 3a3*x^2 + 2a2*x + a1   + 0  @ x1
	// p4 = 3a3*x^2 + 2a2*x + a1   + 0  @ x2
	int i, m;
	float SumDiffRight = 0, SumDiffLeft = 0;
	float VectorOfValues[4];//0 - LeftVal, 1 - RightVal, 2 - LeftSlope, 3 -RightSlope
	float PolynomCoeffs[4] = { 0 };
	float GapInterpolation[FgParams->Motion.GapLength];
	if (FirstGap == 1) {//this is the first interpolation between MotionStruct0 and  MotionStruct1
		VectorOfValues[0] = LeftCurve[FgParams->Motion.SpectrogramTimeBinsSingleMotion
									  - 1];
		VectorOfValues[1] = RightCurve[0];

		for (i = 0; i < FgParams->Motion.NumSamplesForDerivativeEstimation; i++) {//Calculate mean of the diff
			if (i < FgParams->Motion.NumSamplesForDerivativeEstimation - 1)//take 4 elements like in matlab
				SumDiffRight += (RightCurve[i + 1] - RightCurve[i]);
			SumDiffLeft += (LeftCurve[FgParams->Motion.SpectrogramTimeBinsSingleMotion
									  - FgParams->Motion.NumSamplesForDerivativeEstimation + i]
									  - LeftCurve[FgParams->Motion.SpectrogramTimeBinsSingleMotion
												  - FgParams->Motion.NumSamplesForDerivativeEstimation + i
												  - 1]);
			//		printf("%d %lf\n",i, LeftCurve[FgParams->Motion.SpectrogramTimeBins-FgParams->Motion.NumSamplesForDerivativeEstimation+i]);
		}
		VectorOfValues[2] = SumDiffLeft
				/ (FgParams->Motion.NumSamplesForDerivativeEstimation);//the diffrenc ein purpose because it's like that in matlab
		VectorOfValues[3] = SumDiffRight
				/ (FgParams->Motion.NumSamplesForDerivativeEstimation - 1);

		for (i = 0; i < 4; i++) {	//P=inv(A)*V;
			PolynomCoeffs[i] = FgParams->Motion.A12_Inverse[i][0] * VectorOfValues[0]
																					+ FgParams->Motion.A12_Inverse[i][1] * VectorOfValues[1]
																																		   + FgParams->Motion.A12_Inverse[i][2] * VectorOfValues[2]
																																																  + FgParams->Motion.A12_Inverse[i][3] * VectorOfValues[3];
		}
		for (m = 0; m < FgParams->Motion.GapLength; m++) {//implement polyval  i.e insert the missed timebins in the polynom.
			GapInterpolation[m] = PolynomCoeffs[0]
												* powf((FgParams->Motion.SpectrogramTimeBinsSingleMotion + m + 1),
														3)
														+ PolynomCoeffs[1]
																		* powf(
																				(FgParams->Motion.SpectrogramTimeBinsSingleMotion
																						+ m + 1), 2)
																						+ PolynomCoeffs[2]
																										* (FgParams->Motion.SpectrogramTimeBinsSingleMotion + m
																												+ 1) + PolynomCoeffs[3];
			if(isHilbert==1){
				if (GapInterpolation[m] < FgParams->Motion.MinFreqPlusMinus)	//max(gap, minFreq);
					GapInterpolation[m] = FgParams->Motion.MinFreqPlusMinus;
			}
			else{//isHilbert==0
				if (GapInterpolation[m] < FgParams->Motion.MinFreqReal)	//max(gap, minFreq);
					GapInterpolation[m] = FgParams->Motion.MinFreqReal;//+0.01 to be sure that the event will not be cropped in the gap

			}
		}

		//Insert LeftCurve,GapInterpolation,RightCurve into IntrpolatedCurve
		memcpy(IntrpolatedCurve, LeftCurve,
				FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
		memcpy(&IntrpolatedCurve[FgParams->Motion.SpectrogramTimeBinsSingleMotion],
				GapInterpolation, FgParams->Motion.GapLength * sizeof(float));
		memcpy(
				&IntrpolatedCurve[(FgParams->Motion.SpectrogramTimeBinsSingleMotion
						+ FgParams->Motion.GapLength)], RightCurve,
						FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	} else {//FirstGap=0, this is the second interpolation between MotionStruct1 and  MotionStruct2
		VectorOfValues[0] = LeftCurve[FgParams->Motion.SpectrogramTimeBinsTwoMotions
									  - 1];
		VectorOfValues[1] = RightCurve[0];

		for (i = 0; i < FgParams->Motion.NumSamplesForDerivativeEstimation; i++) {//Calculate mean of the diff
			if (i < FgParams->Motion.NumSamplesForDerivativeEstimation - 1)//take 4 elements like in matlab
				SumDiffRight += (RightCurve[i + 1] - RightCurve[i]);
			SumDiffLeft += (LeftCurve[FgParams->Motion.SpectrogramTimeBinsTwoMotions
									  - FgParams->Motion.NumSamplesForDerivativeEstimation + i]
									  - LeftCurve[FgParams->Motion.SpectrogramTimeBinsTwoMotions
												  - FgParams->Motion.NumSamplesForDerivativeEstimation + i
												  - 1]);
		}
		VectorOfValues[2] = SumDiffLeft
				/ (FgParams->Motion.NumSamplesForDerivativeEstimation);//the diffrenc ein purpose because it's like that in matlab
		VectorOfValues[3] = SumDiffRight
				/ (FgParams->Motion.NumSamplesForDerivativeEstimation - 1);

		for (i = 0; i < 4; i++) {	//P=inv(A)*V;
			PolynomCoeffs[i] = FgParams->Motion.A12_Inverse[i][0] * VectorOfValues[0]
																					+ FgParams->Motion.A12_Inverse[i][1] * VectorOfValues[1]
																																		   + FgParams->Motion.A12_Inverse[i][2] * VectorOfValues[2]
																																																  + FgParams->Motion.A12_Inverse[i][3] * VectorOfValues[3];
		}
		for (m = 0; m < FgParams->Motion.GapLength; m++) {//implement polyval  i.e insert the missed timebins in the polynom.
			GapInterpolation[m] = PolynomCoeffs[0]
												* powf((FgParams->Motion.SpectrogramTimeBinsSingleMotion + m + 1),
														3)
														+ PolynomCoeffs[1]
																		* powf(
																				(FgParams->Motion.SpectrogramTimeBinsSingleMotion
																						+ m + 1), 2)
																						+ PolynomCoeffs[2]
																										* (FgParams->Motion.SpectrogramTimeBinsSingleMotion + m
																												+ 1) + PolynomCoeffs[3];

			if(isHilbert==1){
				if (GapInterpolation[m] < FgParams->Motion.MinFreqPlusMinus)	//max(gap, minFreq);
					GapInterpolation[m] = FgParams->Motion.MinFreqPlusMinus;
			}
			else{//isHilbert==0
				if (GapInterpolation[m] < FgParams->Motion.MinFreqReal)	//max(gap, minFreq);
					GapInterpolation[m] = FgParams->Motion.MinFreqReal;

			}
		}

		//Insert LeftCurve,GapInterpolation,RightCurve into IntrpolatedCurve
		memcpy(IntrpolatedCurve, LeftCurve,
				FgParams->Motion.SpectrogramTimeBinsTwoMotions * sizeof(float));
		memcpy(&IntrpolatedCurve[FgParams->Motion.SpectrogramTimeBinsTwoMotions],
				GapInterpolation, FgParams->Motion.GapLength * sizeof(float));
		memcpy(
				&IntrpolatedCurve[(FgParams->Motion.SpectrogramTimeBinsTwoMotions
						+ FgParams->Motion.GapLength)], RightCurve,
						FgParams->Motion.SpectrogramTimeBins * sizeof(float));

	}
	return 0;
}

int ClassifierCorrection(AllFeatures_Struct* FeatureSet,
		float * MotionDistribution, int *y_hat_M, SVM_Struct* SVM_Model) {
	//This function look at specific features in order to correct the classification
	int SNR_thresh = 10;
	if (((FeatureSet->Plus->Max_Black_Curve > 10)
			|| (FeatureSet->Minus->Max_Black_Curve > 10)) && (*y_hat_M == 4)) {
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf(
				"Changed to Normal motion from quasi static, there is a movement\n");
	}
	if ((FeatureSet->Both->SNR_Both < SNR_thresh) && (*y_hat_M != 4)) {
		*y_hat_M = 4;
		MotionDistribution[0] = 0;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 1;
		printf("Decision changed to Quasi Static because of low SNR\n");
	}
	if ((FeatureSet->Both->p1 < 0.1) && (*y_hat_M == 3)) {
		*y_hat_M = 1;
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf("Changed to Normal Motion - probably Standing (P1 < 0.2)\n");
	}
	if (((FeatureSet->Plus->FmaxFpeakMultiplication > 16.5) && (*y_hat_M != 3))) {
		*y_hat_M = 3;
		MotionDistribution[0] = 0;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 1;
		MotionDistribution[3] = 0;
		printf(
				"Changed to Fall, high FmaxFpeakMultiplication_Plus ( > 16.5)\n");
	}
	if (((FeatureSet->Plus->yDiff_Fall_Max_FiftyPrecent_Filtered > 6)
			&& (*y_hat_M != 3))) {
		*y_hat_M = 3;
		MotionDistribution[0] = 0;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 1;
		MotionDistribution[3] = 0;
		printf("Changed to Fall, high polynomialFeatures_Plus\n");
	}
	if ((FeatureSet->Plus->Edge2_50Precent_PeakToAvg > 2.2)
			&& (*y_hat_M == 2)) {
		*y_hat_M = 3;
		MotionDistribution[0] = 0;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 1;
		MotionDistribution[3] = 0;
		printf("Changed to Fall from Sitting,from sitting, big peak2avgplus\n");
	}
	if ((FeatureSet->Plus->Edge2_50Precent_Fmax < 30)
			&& (FeatureSet->Minus->Edge2_50Precent_Fmax < 30)
			&& (*y_hat_M == 3)) {
		*y_hat_M = 1;
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf("Changed to Normal motion, max blue freq smaller than 30\n");
	}
	if (*y_hat_M == 3)//determine between Fall or False alarms by Garbage(=fast hand movements and get ups form the floor)
		SVMClassifier(FeatureSet, MotionDistribution, y_hat_M, SVM_Model);

	if (MotionDistribution[2] < 0.625 && *y_hat_M == 3) {
		*y_hat_M = 1;
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf("Changed to Normal motion, not fully sure that it is a fall\n");

	}

	return 0;
}
int SVMClassifier(AllFeatures_Struct* FeatureSet, float * MotionDistribution,
		int *y_hat_M, SVM_Struct *SVM_Model) {
	//This function Normalize the relevant features and calculates the SVM score by
	//score=FeatureSet_normalized*Beta+Bias, if bigger than 0 it's garabge, to avoid misdetect the threshold will be 2
	//The features are: [17	25	26	30	36	41	42]
	float Score = 0;
	int i, NumOfFeatures = 7;
	float NormalizedSVMFeatures[NumOfFeatures];
	NormalizedSVMFeatures[0] =
			(FeatureSet->Plus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq
					- SVM_Model->x_mean[0]) / SVM_Model->x_std[0];
	NormalizedSVMFeatures[1] =
			(FeatureSet->Plus->yDiff_Raise_Max_FiftyPrecent_Filtered
					- SVM_Model->x_mean[1]) / SVM_Model->x_std[1];
	NormalizedSVMFeatures[2] =
			(FeatureSet->Minus->yDiff_Raise_Max_FiftyPrecent_Filtered
					- SVM_Model->x_mean[2]) / SVM_Model->x_std[2];
	NormalizedSVMFeatures[3] =
			(FeatureSet->Minus->yDiff_Raise_Mean_FiftyPrecent_Filtered
					- SVM_Model->x_mean[3]) / SVM_Model->x_std[3];
	NormalizedSVMFeatures[4] = (FeatureSet->Both->p1 - SVM_Model->x_mean[4])
																																					/ SVM_Model->x_std[4];
	NormalizedSVMFeatures[5] = (FeatureSet->Both->maxEventLength
			- SVM_Model->x_mean[5]) / SVM_Model->x_std[5];
	NormalizedSVMFeatures[6] = (FeatureSet->Both->SumEnergy42
			- SVM_Model->x_mean[6]) / SVM_Model->x_std[6];

	for (i = 0; i < NumOfFeatures; i++) {
		Score += SVM_Model->Beta[i] * NormalizedSVMFeatures[i];
	}
	Score += SVM_Model->Bias;	//add the bias

	if (Score > 2 || (Score > 0 && FeatureSet->Both->SumEnergy42 > 0.5)
			|| (Score > 0 && FeatureSet->Plus->Edge2_50Precent_PeakToAvg > 2.83)) {
		*y_hat_M = 1;	//change to normal motion if think that it is garbage
		MotionDistribution[0] = 1;
		MotionDistribution[1] = 0;
		MotionDistribution[2] = 0;
		MotionDistribution[3] = 0;
		printf("Changed to Normal motion due to SVMClassifier\n");
	}
	return 0;
}

int RandomForrestClassifier(AllFeatures_Struct* FeatureSet,
		Tree_Struct** All_Trees, float * MotionDistribution,
		RF_Struct* RF_Model) {
	float x[6] = { 0 };	//FeatureSet->NumOFRelevantFeatures=6
	int NaNvalue = -1;
	int y_hat_M = 1;
	int TotalTrees = All_Trees[0]->TotalTrees;	//=8
	int treeIdx, currentPredictor, foundFlag = 0;
	int nextNode, currentNode;	//starting point
	int ProbIdx, NumOfClasses = 4;
	float currentCutPoint, MaxProbabilty;
	//The relevant features in this version is according to MotionClassifierRF_NF_Model8 [11,12,13,18,21,39]
	//The features in MotionClassifierRF_NF_Model9: [11,12,18,21,37,40,41,42], not relevant for this version
	//The decisions are: 1- Normal Motion, 2-Sitting, 3 - Fall, 4- Quasistatic
	//Normalize the features: (feature-mean)/std
	x[0] = (FeatureSet->Plus->Edge2_50Precent_Fmax - RF_Model->x_mean[0])
																																					/ RF_Model->x_std[0];	//#11
	x[1] = (FeatureSet->Minus->Edge2_50Precent_Fmax - RF_Model->x_mean[1])
																																					/ RF_Model->x_std[1];	//#12
	x[2] = (FeatureSet->Plus->Edge2_50Precent_AvgTopFive - RF_Model->x_mean[2])
																																					/ RF_Model->x_std[2];	//#13
	x[3] = (FeatureSet->Minus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq
			- RF_Model->x_mean[3]) / RF_Model->x_std[3];	//#18
	x[4] = (FeatureSet->Both->HilbertRatio - RF_Model->x_mean[4])
																																					/ RF_Model->x_std[4];	//#21
	x[5] = (FeatureSet->Plus->FmaxFpeakMultiplication - RF_Model->x_mean[5])
																																					/ RF_Model->x_std[5];	//#39

	for (treeIdx = 0; treeIdx < TotalTrees; treeIdx++) {
		foundFlag = 0;
		currentNode = 0;
		while (foundFlag == 0) {
			currentPredictor = All_Trees[treeIdx]->CutPredictor[currentNode];//what number of feature
			//			printf("currentPredictor %d\n", currentPredictor);

			if (currentPredictor == NaNvalue) {	//Equivalent to if isnan(currentPredictor)
				foundFlag = 1;	//we are at the end and got decision
				for (ProbIdx = 0; ProbIdx < NumOfClasses; ProbIdx++) {
					*(MotionDistribution + ProbIdx) +=
							All_Trees[treeIdx]->ClassProb[4 * currentNode
														  + ProbIdx];	//sum the final decision probabilities for each tree
					//					printf("%d %d %f\n", treeIdx, ProbIdx,
					//							*(MotionDistribution + ProbIdx));
				}
				break;
			} else {
				currentCutPoint = All_Trees[treeIdx]->CutPoint[currentNode];
				//				printf("cutpoint %lf\n", currentCutPoint);

				if (x[currentPredictor] < currentCutPoint)	//go left
					nextNode = All_Trees[treeIdx]->Children[2 * currentNode];//childern is 2 column array
				else
					//go right
					nextNode =
							All_Trees[treeIdx]->Children[2 * currentNode + 1];

				currentNode = nextNode;
			}
		}

	}

	for (ProbIdx = 0; ProbIdx < NumOfClasses; ProbIdx++) {//calculate mean of the probabilities over all the trees
		*(MotionDistribution + ProbIdx) = *(MotionDistribution + ProbIdx)
																																						/ TotalTrees;
	}
	MaxProbabilty = *(MotionDistribution);
	for (ProbIdx = 1; ProbIdx < NumOfClasses; ProbIdx++) {//find maximum probabilty and y_hat_M
		//		printf("%lf\n",*(MotionDistribution+ProbIdx));
		if (*(MotionDistribution + ProbIdx) > MaxProbabilty) {
			MaxProbabilty = *(MotionDistribution + ProbIdx);
			y_hat_M = ProbIdx + 1;		//final decision
		}
	}

	return y_hat_M;
}

int cmpfunc_for_signal(const void * a, const void * b) {//Comparator function for the quick sort(qsort)
	return (*(float*) b < *(float*) a) - (*(float*) b > *(float*) a);
}


int MotionCurveExtraction(Motion_Struct* MotionStruct, float* Mscan[],
		float* Mscan_abs_FFT[],float* Mscan_PostProcess[], float* Pxx2_Hilbert[], float* Pxx2[],
		float* Pxx2_dB[], Edge2_Struct* Edge2, Edge2_Struct* Edge2_Plus,
		Edge2_Struct* Edge2_Minus, FgParams_Struct* FgParams) {

	int p1;
	float T2_dB;

	AbsOfFFT(Mscan_PostProcess, Mscan_abs_FFT, FgParams);

	GET_ROI(Mscan_abs_FFT, FgParams, &p1); //calculate the region of the bins with max energy


#if 0
	//old

	Spectrogram2(Pxx2, Pxx2_dB, &T2_dB, Mscan_PostProcess, p1, FgParams);//calcualte real spectrogram

	CalcCurves(Pxx2, Pxx2_dB, &T2_dB, Edge2, FgParams);//CalcRedStarsCurve3 in Matlab

	HilbertSpectrogram2(Pxx2_Hilbert,Pxx2_dB,&T2_dB,Mscan,p1,Edge2_Plus, Edge2_Minus,FgParams);

	//
#endif
	//	gettimeofday(&tpStart,0);

	HilbertSpectrogram4(Pxx2_Hilbert, Pxx2, Pxx2_dB, &T2_dB, Mscan, p1,
			Edge2_Plus, Edge2_Minus, Edge2, FgParams); //calculate spectrogram after hilbert (complex spectrogram)

	MotionStruct->Edge2 = Edge2;
	MotionStruct->Edge2_Plus = Edge2_Plus;
	MotionStruct->Edge2_Minus = Edge2_Minus;


	return 0;
}




int MotionCurveExtractionPerMscan(Motion_Struct* MotionStruct, float* Mscan[],
		float* Mscan_abs_FFT[],float* Mscan_PostProcess[], float* Pxx2_Hilbert[], float* Pxx2[],
		float* Pxx2_dB[], FgParams_Struct* FgParams) {

	int p1;
	float T2_dB;

	AbsOfFFT(Mscan_PostProcess, Mscan_abs_FFT, FgParams);

	GET_ROI(Mscan_abs_FFT, FgParams, &p1); //calculate the region of the bins with max energy

	HilbertSpectrogram4(Pxx2_Hilbert, Pxx2, Pxx2_dB, &T2_dB, Mscan, p1,
			MotionStruct->Edge2_Plus, MotionStruct->Edge2_Minus, MotionStruct->Edge2, FgParams); //calculate spectrogram after hilbert (complex spectrogram)



	return 0;
}





int MotionCurveExtractionPerMscan2(Motion_Struct* MotionStruct,Pxx2_Plus_Struct* Pxx2_Plus
		,	Pxx2_Minus_Struct* Pxx2_Minus,
		float* Mscan_abs_FFT[],float* Mscan_PostProcessReal[],_Complex float* Mscan_PostProcessHilbert[],  float* Pxx2_Hilbert[], float* Pxx2[],
		float* Pxx2_dB[], FgParams_Struct* FgParams) {

	int p1;
	float T2_dB;

	AbsOfFFT(Mscan_PostProcessReal, Mscan_abs_FFT, FgParams);

	GET_ROI(Mscan_abs_FFT, FgParams, &p1); //calculate the region of the bins with max energy

	HilbertSpectrogram5(Pxx2_Plus,Pxx2_Minus, Pxx2_Hilbert, Pxx2, Pxx2_dB, &T2_dB, Mscan_PostProcessHilbert, p1,
			MotionStruct->Edge2_Plus, MotionStruct->Edge2_Minus, MotionStruct->Edge2, FgParams); //calculate spectrogram after hilbert (complex spectrogram)



	return 0;
}








int FeatureExtractionBasedCurves(Motion_Struct* MotionStruct,
		AllFeatures_Struct* FeatureSet, FgParams_Struct* FgParams) {
	int Type;		//1 - ALL ,2 - only Plus, 3- only Minus
	float sum_Edge2_Minus_Fmax = 0, sum_Edge2_Plus_Fmax = 0;
	int i, maxEventLength;

	//Set default values that without event, if there will be event they will be changed
	FeatureSet->Plus->Edge2_50Precent_Fmax = FgParams->Motion.MinFreqPlusMinus;//=2.5 Hz   #11;
	FeatureSet->Minus->Edge2_50Precent_Fmax = FgParams->Motion.MinFreqPlusMinus; // #12;
	FeatureSet->Plus->Edge2_50Precent_AvgTopFive = FgParams->Motion.MinFreqPlusMinus; // #13
	FeatureSet->Plus->Edge2_50Precent_PeakToAvg = 1; //#15
	FeatureSet->Plus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq =
			FgParams->Motion.MinFreqPlusMinus; //#17
	FeatureSet->Minus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq =
			FgParams->Motion.MinFreqPlusMinus; //#18
	FeatureSet->Plus->yDiff_Raise_Max_FiftyPrecent_Filtered = 0; //#25
	FeatureSet->Minus->yDiff_Raise_Max_FiftyPrecent_Filtered = 0; //#26
	FeatureSet->Minus->yDiff_Raise_Mean_FiftyPrecent_Filtered = 0; //#30
	FeatureSet->Plus->Max_Black_Curve = FgParams->Motion.MinFreqPlusMinus; //#37
	FeatureSet->Minus->Max_Black_Curve = FgParams->Motion.MinFreqPlusMinus; //#38
	FeatureSet->Plus->FmaxFpeakMultiplication = 1; //#39 %SHOULD BE 2.5*2.5/100=0.0625 instead 1 in the origin
	//

	CurveLength(MotionStruct, FgParams);//Calculate the event that defined by the peak (black) curve: maxlegnth, start index and end index

	if (MotionStruct->EventStruct->LengthOfMaxEvent
			> FgParams->Motion.MinEventDuration) {//if the maxlength is bigger than minEventduaration=30
		MotionStruct->EventPlusPassEvent = 1;
		MotionStruct->EventMinusPassEvent = 1;
		Type = 1;	//ALL
		ExtractFeatures(MotionStruct, FeatureSet, Type, FgParams);
		maxEventLength = MotionStruct->EventStruct->LengthOfMaxEvent;

	} else {
		//try calculate the event based only on Plus
		CurveLength(MotionStruct, FgParams);//Calculate the event that defined by the peak (black) curve: maxlegnth, start index and end index

		if (MotionStruct->EventStruct->LengthOfMaxEvent
				> FgParams->Motion.MinEventDuration) {//if the maxlength is bigger than minEventduaration=30
			MotionStruct->EventPlusPassEvent = 1;
			Type = 2;	//Only Plus
			ExtractFeatures(MotionStruct, FeatureSet, Type, FgParams);
			maxEventLength = MotionStruct->EventStruct->LengthOfMaxEvent;
		}
		//try calculate the event based only on Minus
		CurveLength(MotionStruct, FgParams);//Calculate the event that defined by the peak (black) curve: maxlegnth, start index and end index

		if (MotionStruct->EventStruct->LengthOfMaxEvent
				> FgParams->Motion.MinEventDuration) {//if the maxlength is bigger than minEventduaration=30
			MotionStruct->EventMinusPassEvent = 1;
			Type = 3;	//Only Minus
			ExtractFeatures(MotionStruct, FeatureSet, Type, FgParams);
			if (MotionStruct->EventStruct->LengthOfMaxEvent > maxEventLength)//for taking the maximum event from plus and minus
				maxEventLength = MotionStruct->EventStruct->LengthOfMaxEvent;
		}

	}

	//#21 HilbertRatio
	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
		sum_Edge2_Plus_Fmax += MotionStruct->Edge2_Plus->Fmax[i];
		sum_Edge2_Minus_Fmax += MotionStruct->Edge2_Minus->Fmax[i];
		//				printf("50 fmax %d, %f\n", i,
		//						MotionStruct->Edge2_Minus->Fmax[i]);
	}
	FeatureSet->Both->HilbertRatio = log(
			sum_Edge2_Plus_Fmax / sum_Edge2_Minus_Fmax);
	//

	//#30 SNR Both AND #36  Probability of positive movement
	//	FeatureSet->Both->p1=(MotionStruct->Edge2_Plus->SNR_Linear)/(MotionStruct->Edge2_Plus->SNR_Linear+MotionStruct->Edge2_Minus->SNR_Linear);
	SNRFeature(MotionStruct->Edge2_Plus, MotionStruct->Edge2_Minus, FeatureSet,
			FgParams);
	//

	//#40 maxEventLength
	FeatureSet->Both->maxEventLength = maxEventLength;
	//
	//#42, Energy at high frequencies
	//		gettimeofday(&tpStart,0);

	Feature42(MotionStruct, FeatureSet,FgParams);

	//	gettimeofday(&tpStop,0);

	//			f1 = ( (float)( tpStop.tv_sec-tpStart.tv_sec)+ (float)(tpStop.tv_usec)/1000000 ) -  ((float)(tpStart.tv_usec)/1000000) ;
	//				printf (" %f sec\n", f1 );

	return 0;
}



int Feature42(Motion_Struct* MotionStruct, AllFeatures_Struct* FeatureSet,FgParams_Struct* FgParams) {
	//This function is taking the filtered 50precent curve, zero-padds to 256,perform DFT,normalize and caclulates the energy from min_bin
	//perform to Plus and Minus and take the maximum
	//fast movements like hands have energy at high frequencies so the outcome of this function will be big

	int i, k, FFTLength = 256;
	gsl_complex FFTResultForAbs;


	float preFFTFiftyPrecentFiltered[FFTLength];
	memset(preFFTFiftyPrecentFiltered, 0, FFTLength * sizeof(float));//initialize with zeros

	int min_bin = 6;//this is the minimum bin that the energy summation will be, equal to 2.7 Hz
	int startIdx = 4;
	float sum_preFFT_FiftyPrecent_Filtered = 0,
			mean_preFFT_FiftyPrecent_Filtered;
	float sum_nonzero_preFFT_FiftyPrecent_Filtered = 0;

	float SumEnergyPlus = 0, SumEnergyMinus = 0, sumFFTFiftyPrecentFiltered;
	float absFFTFiftyPrecentFiltered[FFTLength / 2 + 1];




	ne10_fft_cpx_float32_t FFTResult[(256 / 2) + 1]={};//FFTLength=256
	//	ne10_fft_r2c_cfg_float32_t FFT_cfg;                     // An FFT "configuration structure"
	//	FFT_cfg = ne10_fft_alloc_r2c_float32(FFTLength);


	//Plus
	if (MotionStruct->EventPlusPassEvent == 1) {
		//first zero-pad the signal to Length=256
		if (startIdx < MotionStruct->EventStruct->chosenStart)//max(Event_Plus.chosenStart,5)
			startIdx = MotionStruct->EventStruct->chosenStart;
		for (i = 0; i <= MotionStruct->EventStruct->chosenEnd - startIdx; i++) {
			preFFTFiftyPrecentFiltered[i] =
					MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i + startIdx];

			sum_preFFT_FiftyPrecent_Filtered +=
					MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i + startIdx];

		}
		mean_preFFT_FiftyPrecent_Filtered = sum_preFFT_FiftyPrecent_Filtered
				/ FFTLength;
		for (i = 0; i < FFTLength; i++) {			//remove DC
			//			preFFT_FiftyPrecent_Filtered[i] -= mean_preFFT_FiftyPrecent_Filtered;
			preFFTFiftyPrecentFiltered[i] = preFFTFiftyPrecentFiltered[i]
																	   - mean_preFFT_FiftyPrecent_Filtered;
			if (MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i] != 0)
				sum_nonzero_preFFT_FiftyPrecent_Filtered += 1;//sum of non zero elements
		}
		if (sum_nonzero_preFFT_FiftyPrecent_Filtered != 0) {
			//FFT 256 PTS
			sumFFTFiftyPrecentFiltered = 0;

			ne10_fft_r2c_1d_float32_neon(FFTResult,(ne10_float32_t*)preFFTFiftyPrecentFiltered,FgParams->Motion.FFT_Feature42);

			for (k = 0; k <= FFTLength / 2; k++) {//Length/2 because only one side is wanted
				FFTResultForAbs.dat[0]=FFTResult[k].r;
				FFTResultForAbs.dat[1]=FFTResult[k].i;
				absFFTFiftyPrecentFiltered[k]=gsl_complex_abs2(FFTResultForAbs);
				sumFFTFiftyPrecentFiltered += absFFTFiftyPrecentFiltered[k];
			}
			for (k = min_bin; k <= FFTLength / 2; k++) {
				SumEnergyPlus += absFFTFiftyPrecentFiltered[k]
															/ sumFFTFiftyPrecentFiltered;//normalize and sum from min_bin
			}


		} else
			SumEnergyPlus = 0;

	} else
		SumEnergyPlus = 0;
	//minus
	sum_nonzero_preFFT_FiftyPrecent_Filtered = 0;
	sum_preFFT_FiftyPrecent_Filtered = 0;
	memset(preFFTFiftyPrecentFiltered, 0, FFTLength * sizeof(float));//initialize with zeros
	if (MotionStruct->EventMinusPassEvent == 1) {
		//first zero-pad the signal to Length=256
		if (startIdx < MotionStruct->EventStruct->chosenStart)//max(Event_Plus.chosenStart,5)
			startIdx = MotionStruct->EventStruct->chosenStart;
		for (i = 0; i <= MotionStruct->EventStruct->chosenEnd - startIdx; i++) {
			preFFTFiftyPrecentFiltered[i] =
					MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i
																	 + startIdx];
			sum_preFFT_FiftyPrecent_Filtered +=
					MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i
																	 + startIdx];
		}
		mean_preFFT_FiftyPrecent_Filtered = sum_preFFT_FiftyPrecent_Filtered
				/ FFTLength;
		for (i = 0; i < FFTLength; i++) {			//remove DC
			preFFTFiftyPrecentFiltered[i] -= mean_preFFT_FiftyPrecent_Filtered;

			//			preFFT_FiftyPrecent_Filtered[i] = preFFT_FiftyPrecent_Filtered[i]
			//																		   - mean_preFFT_FiftyPrecent_Filtered;
			if (MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i] != 0)
				sum_nonzero_preFFT_FiftyPrecent_Filtered += 1;//sum of non zero elements
		}
		if (sum_nonzero_preFFT_FiftyPrecent_Filtered != 0) {
			//FFT 256 PTS
			sumFFTFiftyPrecentFiltered = 0;

			ne10_fft_r2c_1d_float32_neon(FFTResult,(ne10_float32_t*)preFFTFiftyPrecentFiltered,FgParams->Motion.FFT_Feature42);


			for (k = 0; k <= FFTLength / 2; k++) {//Length/2 because only one side is wanted
				FFTResultForAbs.dat[0]=FFTResult[k].r;
				FFTResultForAbs.dat[1]=FFTResult[k].i;
				absFFTFiftyPrecentFiltered[k]=gsl_complex_abs2(FFTResultForAbs);
				sumFFTFiftyPrecentFiltered += absFFTFiftyPrecentFiltered[k];
			}
			for (k = min_bin; k <= FFTLength / 2; k++) {
				SumEnergyMinus += absFFTFiftyPrecentFiltered[k]
															 / sumFFTFiftyPrecentFiltered;//normalize and sum from min_bin
			}


		} else
			SumEnergyMinus = 0;
	} else
		SumEnergyMinus = 0;

	if (SumEnergyPlus > SumEnergyMinus)
		FeatureSet->Both->SumEnergy42 = SumEnergyPlus;
	else
		FeatureSet->Both->SumEnergy42 = SumEnergyMinus;

	//	ne10_fft_destroy_r2c_float32(FFT_cfg);
	return 0;
}









int ExtractFeatures(Motion_Struct* MotionStruct,
		AllFeatures_Struct* FeatureSet, int Type, FgParams_Struct *FgParams) {
	int TimeBin_50Precent_Fmax_Plus = MotionStruct->EventStruct->chosenStart, TimeBin_50Precent_Fmax_Minus = MotionStruct->EventStruct->chosenStart;
	float Arr_for_sort_Plus[MotionStruct->EventStruct->LengthOfChosenEvent + 1];
	float Edge2_AvgTopFive_Plus = 0;
	float Edge2_noDC_Favg_Plus = 0;
	int num_of_Edge2_noDC_Favg_Plus = 0;
	int i;

	//	for(i=0;i<FgParams->Motion.SpectrogramTimeBins;i++)
	//		printf("fifty fmax filtered %d %lf\n", i,
	//				MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]);
	for (i = MotionStruct->EventStruct->chosenStart;
			i <= MotionStruct->EventStruct->chosenEnd; i++) {//run on the timebins in the event



		if (Type == 1 || Type == 2) {			//All or only Plus
			Arr_for_sort_Plus[i - MotionStruct->EventStruct->chosenStart] =
					MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i];	//for Edge2_AvgTopFive

			//#11 Edge2_50Precent_Fmax_Plus
			if (MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																> FeatureSet->Plus->Edge2_50Precent_Fmax) {
				FeatureSet->Plus->Edge2_50Precent_Fmax =
						MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i];
				TimeBin_50Precent_Fmax_Plus = i;			//timebin of Fmax
			}

			//#15 - Edge2+ 50% Peak2Avg: Step 1
			if (MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]
																> FgParams->Motion.MinFreqPlusMinus) {//calculate mean for fifty_precent_filtered without Fmin bins
				Edge2_noDC_Favg_Plus +=
						MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i];
				num_of_Edge2_noDC_Favg_Plus += 1;
			}

			//#37 - max black curve plus
			if (MotionStruct->Edge2_Plus->Peak_Filtered[i]
														> FeatureSet->Plus->Max_Black_Curve)
				FeatureSet->Plus->Max_Black_Curve =
						MotionStruct->Edge2_Plus->Peak_Filtered[i];
			//
		}

		if (Type == 1 || Type == 3) {			//All or only Minus

			//#12 Edge2_50Precent_Fmax_Minus
			if (MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]
																 > FeatureSet->Minus->Edge2_50Precent_Fmax) {
				FeatureSet->Minus->Edge2_50Precent_Fmax =
						MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i];
				TimeBin_50Precent_Fmax_Minus = i;			//timebin of Fmax

			}
			//			printf("line 1147 %d  ,   %lf\n",i,MotionStruct->Edge2_Minus->FiftyPrecent_Filtered[i]);

			//

			//#38 - max black curve Minus
			if (MotionStruct->Edge2_Minus->Peak_Filtered[i]
														 > FeatureSet->Minus->Max_Black_Curve)
				FeatureSet->Minus->Max_Black_Curve =
						MotionStruct->Edge2_Minus->Peak_Filtered[i];
			//
		}

	}
	if (Type == 1 || Type == 2) {			//All or only Plus
		//#13 Edge2_50Precent_AvgTopFive_Plus

		qsort(Arr_for_sort_Plus,
				MotionStruct->EventStruct->LengthOfChosenEvent + 1,
				sizeof(float), cmpfunc_for_signal);	//quick sort the 20(=MedianValue) values
		for (i = 0; i < FgParams->Motion.TopMaxFreq; i++)
			Edge2_AvgTopFive_Plus +=
					Arr_for_sort_Plus[MotionStruct->EventStruct->LengthOfChosenEvent
									  - i];
		FeatureSet->Plus->Edge2_50Precent_AvgTopFive = Edge2_AvgTopFive_Plus
				/ (FgParams->Motion.TopMaxFreq);
		//		//
		//#15 - Edge2+ 50% Peak2Avg: Step 2
		if (num_of_Edge2_noDC_Favg_Plus == 0)//equivalent to if isempty(Edge2_noDC) in matlab
			num_of_Edge2_noDC_Favg_Plus = FgParams->Motion.MinFreqPlusMinus;
		Edge2_noDC_Favg_Plus = Edge2_noDC_Favg_Plus
				/ (num_of_Edge2_noDC_Favg_Plus);
		//
		if (Edge2_noDC_Favg_Plus != 0)	//otherwise Edge2_PeakToAvg_Plus stay 0
			FeatureSet->Plus->Edge2_50Precent_PeakToAvg =
					FeatureSet->Plus->Edge2_50Precent_AvgTopFive
					/ Edge2_noDC_Favg_Plus;
		else
			FeatureSet->Plus->Edge2_50Precent_PeakToAvg = 0;//THINK ABOUT IT AGAIN
		//

		//#17 Edge2_maxPeakFreq_fromTimeBinWithMaxFreq_Plus
		FeatureSet->Plus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq =
				MotionStruct->Edge2_Plus->Peak_Filtered[TimeBin_50Precent_Fmax_Plus];
		//

		//#25 yDiff_Raise_Max__Edge2_Plus.FiftyPrecent_Filtered & #27 - yDiff_Fall__Max__Edge2_Plus.FiftyPrecent_Filtered
		PolynomialFeatures2(MotionStruct->Edge2_Plus, MotionStruct,
				FeatureSet->Plus);
		//
		//#39 FmaxFpeakMultiplication_Plus
		FeatureSet->Plus->FmaxFpeakMultiplication =
				FeatureSet->Plus->Edge2_50Precent_Fmax
				* (MotionStruct->Edge2_Plus->Peak_Filtered[TimeBin_50Precent_Fmax_Plus])
				* 0.01;
		//
	}
	if (Type == 1 || Type == 3) {		//All or only Minus
		//#18 Edge2_maxPeakFreq_fromTimeBinWithMaxFreq_Minus
		FeatureSet->Minus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq =
				MotionStruct->Edge2_Minus->Peak_Filtered[TimeBin_50Precent_Fmax_Minus];
		//
		//#26 yDiff_Raise_Max__Edge2_Minus.FiftyPrecent_Filtered, #30 yDiff_Raise_Mean_Edge2_Minus.FiftyPrecent_Filtered
		PolynomialFeatures2(MotionStruct->Edge2_Minus, MotionStruct,
				FeatureSet->Minus);
		//
		//#40 FmaxFpeakMultiplication_Minus
		FeatureSet->Minus->FmaxFpeakMultiplication =
				FeatureSet->Minus->Edge2_50Precent_Fmax
				* (MotionStruct->Edge2_Minus->Peak_Filtered[TimeBin_50Precent_Fmax_Minus])
				* 0.01;
		//
	}

	return 0;
}

int PolynomialFeatures2(Edge2_Struct* Edge2, Motion_Struct *MotionStruct,
		Features_Struct *Features) {
	int MedianValue = 20, peakIdx, NumBinsForPolynom;
	int howManySamplesEachSide = 20;
	int timeBin_start = 0;
	int timeBin_stop = MotionStruct->EventStruct->LengthOfChosenEvent;
	int startIdx, stopIdx, i;
	int polynomialOrder = 3; // p4*x^3 + p3*x^2 + p2*x + p1
	double p[polynomialOrder + 1]; //coefficients
	float Edge2_50Precent_MedianFiltered[MotionStruct->EventStruct->LengthOfChosenEvent
										 + MedianValue - 1];
	float MaxValue;
	float yDiff_Raise = 0, yDiff_Fall = 0, yDiff_Raise_mean = 0,
			num_yDiff_Raise_mean = 0;
	//float  yDiff_Fall=0, yDiff_minus_mean=0,num_yDiff_minus_mean=0;//not used

	double *TimeBinsForPolynom, *ValuesForPolynom, yDiff;
	if (howManySamplesEachSide
			> floor(MotionStruct->EventStruct->LengthOfChosenEvent / 2) - 1) //% if the Event is less than 41 samples
		howManySamplesEachSide = floor(
				MotionStruct->EventStruct->LengthOfChosenEvent / 2) - 1;
	MedianFilterFor50Precent(Edge2, MotionStruct,
			Edge2_50Precent_MedianFiltered, MedianValue); //Filter the 50precentfiltered curve
	MaxOfArr(Edge2_50Precent_MedianFiltered, &peakIdx, &MaxValue,
			MotionStruct->EventStruct->LengthOfChosenEvent); //find max idx of 50Precent_MedianFiltered
	//GetMargin in Matlab:
	if (timeBin_start < (peakIdx - howManySamplesEachSide))
		timeBin_start = (peakIdx - howManySamplesEachSide);
	if (timeBin_stop > peakIdx + howManySamplesEachSide)
		timeBin_stop = peakIdx + howManySamplesEachSide;

	if (timeBin_start == 0) {
		startIdx = 0;	//was 1
		stopIdx = 2 * howManySamplesEachSide;	//was+1
	} else if (timeBin_stop == MotionStruct->EventStruct->LengthOfChosenEvent) {
		startIdx = MotionStruct->EventStruct->LengthOfChosenEvent
				- (2 * howManySamplesEachSide) - 1;
		stopIdx = MotionStruct->EventStruct->LengthOfChosenEvent - 1;
	} else {
		startIdx = timeBin_start;
		stopIdx = timeBin_stop;
	}
	//

	NumBinsForPolynom = stopIdx - startIdx + 1;
	TimeBinsForPolynom = (double*) malloc((NumBinsForPolynom) * sizeof(double));//xdata
	ValuesForPolynom = (double*) malloc((NumBinsForPolynom) * sizeof(double));//ydata

	for (i = 0; i < NumBinsForPolynom; i++) {
		TimeBinsForPolynom[i] = startIdx + i + 1;//1 FOR BE LIKE MATLAB????THINKABOUT IT....
		ValuesForPolynom[i] = Edge2_50Precent_MedianFiltered[startIdx + i];
		//		printf("vforpo %d %f\n",i,ValuesForPolynom[i]);
	}
	polyfit(TimeBinsForPolynom, ValuesForPolynom, NumBinsForPolynom,
			polynomialOrder, p);	//MUST TO BE float

	for (i = 0; i <= NumBinsForPolynom; i++) {
		//yDiff_plus = max(yDiff);yDiff_minus = max(-yDiff); yDiff_plus_mean = mean(yDiff(yDiff>=0));yDiff_minus_mean = -mean(yDiff(yDiff<=0));
		if (i < NumBinsForPolynom)
			yDiff = 3 * p[3] * powf(TimeBinsForPolynom[i], 2)
			+ 2 * p[2] * TimeBinsForPolynom[i] + p[1];
		else
			yDiff = 0;

		if (yDiff > yDiff_Raise)	//max(yDiff)
			yDiff_Raise = yDiff;
		if (-yDiff > yDiff_Fall)	//max(-yDiff)
			yDiff_Fall = -yDiff;

		if (yDiff >= 0) {
			yDiff_Raise_mean += yDiff;
			num_yDiff_Raise_mean += 1;
		}

		//		if(-yDiff>yDiff_Fall)
		//			yDiff_Fall=-yDiff;
		//		if (yDiff[i]<=0){//not used
		//			yDiff_minus_mean-=yDiff[i];
		//			num_yDiff_minus_mean+=1;
		//		}
	}
	//#25 yDiff_Raise_Max__Edge2_Plus.FiftyPrecent_Filtered #26yDiff_Raise_Max__Edge2_Minus.FiftyPrecent_Filtered:
	Features->yDiff_Raise_Max_FiftyPrecent_Filtered = yDiff_Raise;
	//
	//#27 yDiff_Fall_Max__Edge2_Plus.FiftyPrecent_Filtered #28yDiff_Fall_Max__Edge2_Minus.FiftyPrecent_Filtered:
	Features->yDiff_Fall_Max_FiftyPrecent_Filtered = yDiff_Fall;
	//
	//#30 - yDiff_Raise_Mean_Edge2_Minus.FiftyPrecent_Filtered
	Features->yDiff_Raise_Mean_FiftyPrecent_Filtered = yDiff_Raise_mean
			/ num_yDiff_Raise_mean;

	//	//#27 - yDiff_Fall__Max__Edge2_Plus.FiftyPrecent_Filtered
	//	Features->yDiff_Fall_Max_FiftyPrecent_Filtered=yDiff_Fall;
	//
	//	Features->yDiff_Fall_Mean_FiftyPrecent_Filtered=yDiff_minus_mean/num_yDiff_minus_mean; not used
	free(TimeBinsForPolynom);
	free(ValuesForPolynom);

	return 0;
}

int SNRFeature(Edge2_Struct* Edge2_Plus, Edge2_Struct* Edge2_Minus,
		AllFeatures_Struct *FeatureSet, FgParams_Struct *FgParams) {
	int i;
	float SumOf_SumEnergy_Post_Plus = 0, MeanOf_SumEnergy_Post_Plus;
	float SumOf_SumEnergy_Post_Minus = 0, MeanOf_SumEnergy_Post_Minus;
	float sum_T1_t_Plus = 0, sum_T1_t_Minus = 0;
	float T1_t_mean_Plus, T1_t_mean_Minus;

	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
		SumOf_SumEnergy_Post_Plus += Edge2_Plus->SumEnergy_Post[i];	//for mean calculation
		SumOf_SumEnergy_Post_Minus += Edge2_Minus->SumEnergy_Post[i];
		sum_T1_t_Plus += Edge2_Plus->T1_t[i];	//for mean calculation
		sum_T1_t_Minus += Edge2_Minus->T1_t[i];	//for mean calculation

	}
	T1_t_mean_Plus = (sum_T1_t_Plus / FgParams->Motion.SpectrogramTimeBins);
	T1_t_mean_Minus = (sum_T1_t_Minus / FgParams->Motion.SpectrogramTimeBins);

	MeanOf_SumEnergy_Post_Plus = SumOf_SumEnergy_Post_Plus
			/ FgParams->Motion.SpectrogramTimeBins;
	MeanOf_SumEnergy_Post_Minus = SumOf_SumEnergy_Post_Minus
			/ FgParams->Motion.SpectrogramTimeBins;

	Edge2_Plus->SNR_Linear = powf(10,
			(10 * log10(MeanOf_SumEnergy_Post_Plus) - T1_t_mean_Plus) * 0.1);//mean in dB =10*log10(SumOf_SumEnergy_Post_Plus/SpectrogramTimeBins)
	Edge2_Minus->SNR_Linear = powf(10,
			(10 * log10(MeanOf_SumEnergy_Post_Minus) - T1_t_mean_Minus) * 0.1);

	FeatureSet->Both->p1 = (Edge2_Plus->SNR_Linear)
																																					/ (Edge2_Plus->SNR_Linear + Edge2_Minus->SNR_Linear);
	FeatureSet->Both->SNR_Both = (10 * log10(MeanOf_SumEnergy_Post_Plus)
	- T1_t_mean_Plus)
																																					+ (10 * log10(MeanOf_SumEnergy_Post_Minus) - T1_t_mean_Minus);//caluclation for feature SNR_Both

	return 0;
}

int CurveLength(Motion_Struct *MotionStruct, FgParams_Struct *FgParams) {
	//This function searches for events with minimal length of MinEventDuration=30
	//	        | x(n-1)=0    | x(n-1)=1    |
	//	 x(n)=0 | DO NOTHING  | STOP EVENT  |
	//	 x(n)=1 | START EVENT |  DO NOTHING |
	int MaxIdx;
	float FmaxPerEventPlus, FmaxPerEventMinus;
	//	float minValue = 2.512562814070352; //Fmin in Hz
	int i, x_before, x_current;
	int size_list = floor(FgParams->Motion.SpectrogramTimeBins / 2);
	int startList[size_list];
	int endList[size_list];
	int LengthOfMaxEvent = 0, idxChosenEvent;
	int pointer = 0; //number of the event in startlist & endlist
	float FmaxPerEvent, GlobalFmax = 0;
	int oneHotEdge2_Plus_Peak_Filtered[FgParams->Motion.SpectrogramTimeBins + 2]; //+2 for inital state and end state
	memset(oneHotEdge2_Plus_Peak_Filtered, 0,
			(FgParams->Motion.SpectrogramTimeBins + 2) * sizeof(int));

	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
//						printf("peak fil %d,    %f\n", i,
//								MotionStruct->Edge2_Minus->max[i]);
		if (MotionStruct->Edge2->Peak_Filtered[i] > FgParams->Motion.MinFreqReal) { //put 1's in all the bins that pass the minValue otherwise stay 0
			oneHotEdge2_Plus_Peak_Filtered[i + 1] = 1;

		} //first and last stayed 0
	}

	for (i = 1; i < FgParams->Motion.SpectrogramTimeBins + 2; i++) {
		x_before = oneHotEdge2_Plus_Peak_Filtered[i - 1];
		x_current = oneHotEdge2_Plus_Peak_Filtered[i];

		if ((x_before == 0) && (x_current == 1)) { //start of event
			startList[pointer] = i - 1;
		} else if ((x_before == 1) && (x_current == 0)) { //end of event and
			if (((i - 2) - startList[pointer]) > FgParams->Motion.MinEventDuration) { //the event is bigger than minduration
				endList[pointer] = i - 2;
				if ((i - 2) - startList[pointer] > LengthOfMaxEvent) {
					LengthOfMaxEvent = (i - 2) - startList[pointer]; //was i-1
					//					idx_of_LengthOfMaxEvent=pointer;
				}

				//find the maximum fiftyprecent_filtered in the event between plus and minus
				MaxOfArr(
						MotionStruct->Edge2_Plus->FiftyPrecent_Filtered
						+ startList[pointer], &MaxIdx,
						&FmaxPerEventPlus,
						endList[pointer] - startList[pointer] + 1);
				MaxOfArr(
						MotionStruct->Edge2_Minus->FiftyPrecent_Filtered
						+ startList[pointer], &MaxIdx,
						&FmaxPerEventMinus,
						endList[pointer] - startList[pointer] + 1);

				FmaxPerEvent = Max(FmaxPerEventPlus, FmaxPerEventMinus);

				if (FmaxPerEvent > GlobalFmax) {//the event is chosen by the the globalfmax
					GlobalFmax = FmaxPerEvent;
					idxChosenEvent = pointer;
				}
				pointer += 1;
			}

		}
	}
	//	printf("%d\n", LengthOfMaxEvent);

	if (LengthOfMaxEvent > 0)
		MotionStruct->EventStruct->LengthOfMaxEvent = LengthOfMaxEvent;
	else
		//event not found
		MotionStruct->EventStruct->LengthOfMaxEvent = 0;

	MotionStruct->EventStruct->chosenStart = startList[idxChosenEvent];
	MotionStruct->EventStruct->chosenEnd = endList[idxChosenEvent];
	MotionStruct->EventStruct->LengthOfChosenEvent = endList[idxChosenEvent]
															 - startList[idxChosenEvent];

	return 0;
}


int CalcCurvesHilbert(Pxx2_Plus_Struct *Pxx2_Plus,
		Pxx2_Minus_Struct* Pxx2_Minus, Edge2_Struct* Edge2_Plus,
		Edge2_Struct* Edge2_Minus, FgParams_Struct* FgParams) {//CalcRedStarsCurve3_FirstStep2 in Matlab
	int k, i, idx_FreqBins_Spectogram = 0;
	int MedianType, IsHilbert;
	int AllFreqAboveThresh_Plus[FgParams->Motion.SpectrogramFreqBinsHilbert / 2],
	AllFreqAboveThresh_Minus[FgParams->Motion.SpectrogramFreqBinsHilbert / 2];
	int idx_FreqAboveThresh_Plus, idx_FreqAboveThresh_Minus;
	int MaxPxx2_dB_idx_Plus, MaxPxx2_dB_idx_Minus;
	float freq_increment = 1.25;					//=(Fs/LengthOfDFT);
	float MaxPxx2_dB_Plus, MaxPxx2_dB_Minus;
	//	float SumEnergy_Post_Minus,SumEnergy_Post_Plus;
	//	int Freq_50_PrecentIdx_Plus,Freq_50_PrecentIdx_Minus;
	int flag_50_PrecentIdx_Plus, flag_50_PrecentIdx_Minus;
	float Energy_50_Precent_Plus, Energy_50_Precent_Minus;
	float SpectrogramFreqBinsInHz[FgParams->Motion.SpectrogramFreqBinsHilbert / 2];//, Energy_50_Precent;
	for (float i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert / 2; i++) {	//create the frequency vector inHz
		SpectrogramFreqBinsInHz[idx_FreqBins_Spectogram] = i * freq_increment;
		idx_FreqBins_Spectogram += 1;
	}
	//	memset(Edge2_Plus->Peak,0,174*sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	//	memset(Edge2_Minus->Peak,0,174*sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions

	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
		idx_FreqAboveThresh_Plus = 0;
		idx_FreqAboveThresh_Minus = 0;
		MaxPxx2_dB_Plus = Pxx2_Plus->dB[FgParams->Motion.FminBin][i];
		MaxPxx2_dB_Minus = Pxx2_Minus->dB[FgParams->Motion.FminBin][i];
		MaxPxx2_dB_idx_Plus = FgParams->Motion.FminBin;
		MaxPxx2_dB_idx_Minus = FgParams->Motion.FminBin;

		for (k = FgParams->Motion.FminBin;
				k
				< FgParams->Motion.SpectrogramFreqBinsHilbert / 2
				- FgParams->Motion.SpectrogramNoiseFreqBins; k++) {//frequency range [FminBin,SpectrogramFreqBins-SpectrogramNoiseFreqBins]
			if (Pxx2_Plus->dB[k][i] > Edge2_Plus->T1_t[i]) {//if the powfer is passing the threshold
				AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus] = k;
				idx_FreqAboveThresh_Plus += 1;
				//				printf("%d %d\n",idx_FreqAboveThresh_Plus,AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus]);
			}

			if (Pxx2_Plus->dB[k][i] > MaxPxx2_dB_Plus) {//find MaxPxx2_dB_Plus
				MaxPxx2_dB_Plus = Pxx2_Plus->dB[k][i];
				MaxPxx2_dB_idx_Plus = k;
			}

			if (Pxx2_Minus->dB[k][i] > Edge2_Minus->T1_t[i]) {//if the powfer is passing the threshold
				AllFreqAboveThresh_Minus[idx_FreqAboveThresh_Minus] = k;
				idx_FreqAboveThresh_Minus += 1;
			}

			if (Pxx2_Minus->dB[k][i] > MaxPxx2_dB_Minus) {//find MaxPxx2_dB_Plus
				MaxPxx2_dB_Minus = Pxx2_Minus->dB[k][i];
				MaxPxx2_dB_idx_Minus = k;
			}
		}

		if (idx_FreqAboveThresh_Plus > 0) {	//there were frequencies above the threshold
			Edge2_Plus->maxFreqIdxs[i] =
					AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus - 1];
			Edge2_Plus->maxFreqEnergy[i] =
					Pxx2_Plus->dB[AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus
														  - 1]][i];
			Edge2_Plus->PeakEnergy[i] = MaxPxx2_dB_Plus;
			Edge2_Plus->maxPeakIdxs[i] = MaxPxx2_dB_idx_Plus;
			Edge2_Plus->Fmax[i] =
					SpectrogramFreqBinsInHz[AllFreqAboveThresh_Plus[idx_FreqAboveThresh_Plus
																	- 1]];

			if (FgParams->Motion.TruncateHilbert == 0)//in this case there is zero-pad at the beginning
				Edge2_Plus->Peak[i + FgParams->Motion.MedianValue / 2] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx_Plus];//first MedianValue/2(=10) bins are 0 for the inital sates of the median filter
			else {
				//there is no-zero padding
				Edge2_Plus->Peak[i] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx_Plus];
				//				printf("%d %f peak\n",i,Edge2_Plus->Peak[i]);
				//		Edge2_Plus->SumEnergy[i]=SumEnergy_Plus;
			}
		} else {	//take default values
			Edge2_Plus->maxFreqIdxs[i] = FgParams->Motion.FminBin;
			Edge2_Plus->maxFreqEnergy[i] = 0;//should be Pxx2_Plus->dB[FminBin][i]! but in matlab 0..
			Edge2_Plus->PeakEnergy[i] = Pxx2_Plus->dB[FgParams->Motion.FminBin][i];
			Edge2_Plus->maxPeakIdxs[i] = FgParams->Motion.FminBin;
			Edge2_Plus->Fmax[i] = SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			if (FgParams->Motion.TruncateHilbert == 0)//in this case there is zero-pad at the beginning
				Edge2_Plus->Peak[i + FgParams->Motion.MedianValue / 2] =
						SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			else
				//there is no-zero padding
				Edge2_Plus->Peak[i] =
						SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			//		Edge2_Plus->SumEnergy[i]=0;
		}

		if (idx_FreqAboveThresh_Minus > 0) {//there were frequencies above the threshold
			Edge2_Minus->maxFreqIdxs[i] =
					AllFreqAboveThresh_Minus[idx_FreqAboveThresh_Minus - 1];
			Edge2_Minus->maxFreqEnergy[i] =
					Pxx2_Minus->dB[AllFreqAboveThresh_Minus[idx_FreqAboveThresh_Minus
															- 1]][i];
			Edge2_Minus->PeakEnergy[i] = MaxPxx2_dB_Minus;
			Edge2_Minus->maxPeakIdxs[i] = MaxPxx2_dB_idx_Minus;
			Edge2_Minus->Fmax[i] =
					SpectrogramFreqBinsInHz[AllFreqAboveThresh_Minus[idx_FreqAboveThresh_Minus
																	 - 1]];
			if (FgParams->Motion.TruncateHilbert == 0)//in this case there is zero-pad at the beginning
				Edge2_Minus->Peak[i + FgParams->Motion.MedianValue / 2] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx_Minus];//first MedianValue/2(=10) bins are 0 for the inital sates of the median filter
			else
				//there is no-zero padding
				Edge2_Minus->Peak[i] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx_Minus];//		Edge2_Plus->SumEnergy[i]=SumEnergy_Plus;
		} else {
			Edge2_Minus->maxFreqIdxs[i] = FgParams->Motion.FminBin;
			Edge2_Minus->maxFreqEnergy[i] = 0;// should be Pxx2_Minus->dB[FminBin][i];
			Edge2_Minus->PeakEnergy[i] = Pxx2_Minus->dB[FgParams->Motion.FminBin][i];
			Edge2_Minus->maxPeakIdxs[i] = FgParams->Motion.FminBin;
			Edge2_Minus->Fmax[i] = SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			if (FgParams->Motion.TruncateHilbert == 0)//in this case there is zero-pad at the beginning
				Edge2_Minus->Peak[i + FgParams->Motion.MedianValue / 2] =
						SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			else
				//there is no-zero padding
				Edge2_Minus->Peak[i] =
						SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];//		Edge2_Plus->SumEnergy[i]=0;
		}

		if (Edge2_Plus->PeakEnergy[i] > Edge2_Minus->PeakEnergy[i]) {//store the max energy for each timebin
			Edge2_Plus->PeakEnergy_PM[i] = Edge2_Plus->PeakEnergy[i];
			Edge2_Minus->PeakEnergy_PM[i] = Edge2_Plus->PeakEnergy[i];
		} else {
			Edge2_Plus->PeakEnergy_PM[i] = Edge2_Minus->PeakEnergy[i];
			Edge2_Minus->PeakEnergy_PM[i] = Edge2_Minus->PeakEnergy[i];

		}
	}




	MedianType = 1;		//TRUNCATE
	MedianFilter(Edge2_Plus, FgParams, MedianType);//median filter on the peak(black) curve with MedianValue=20 , MedianType=0 ZERPOAD,MedianType=1 TRUNCATE
	MedianFilter(Edge2_Minus, FgParams, MedianType);//median filter on the peak(black) curve with MedianValue=20,  MedianType=0 ZERPOAD,MedianType=1 TRUNCATE

	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
		Energy_50_Precent_Plus = (Edge2_Plus->PeakEnergy_PM[i]
															+ Edge2_Plus->maxFreqEnergy[i]) / 2;
		Energy_50_Precent_Minus = (Edge2_Minus->PeakEnergy_PM[i]
															  + Edge2_Minus->maxFreqEnergy[i]) / 2;
		Edge2_Plus->SumEnergy_Post[i]=0;
		Edge2_Minus->SumEnergy_Post[i]=0;

		flag_50_PrecentIdx_Plus = 0;
		flag_50_PrecentIdx_Minus = 0;

		for (k = Edge2_Plus->maxFreqIdxs[i]; k >= Edge2_Plus->maxPeakIdxs[i];
				k--) {//search backwards the first freq index that passes Energy_50_Precent
			if (Pxx2_Plus->dB[k][i] > Energy_50_Precent_Plus) {
				Edge2_Plus->FiftyPrecent[i] = SpectrogramFreqBinsInHz[k];
				Edge2_Plus->FiftyPrecentIdxs[i] = k;
				flag_50_PrecentIdx_Plus = 1;
				break;
			}
		}

		if (flag_50_PrecentIdx_Plus == 0) {
			Edge2_Plus->FiftyPrecent[i] =
					SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			Edge2_Plus->Peak_Filtered[i] =
					SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			Edge2_Plus->FiftyPrecentIdxs[i] = FgParams->Motion.FminBin;
			Edge2_Plus->maxPeakIdxs[i] = FgParams->Motion.FminBin;

		}

		for (k = Edge2_Minus->maxFreqIdxs[i]; k >= Edge2_Minus->maxPeakIdxs[i];
				k--) {//search backwards the first freq index that passes Energy_50_Precent
			//			printf("%d %d %lf ",k,i,Pxx2_Minus->dB[k][i]);

			if (Pxx2_Minus->dB[k][i] > Energy_50_Precent_Minus) {
				Edge2_Minus->FiftyPrecent[i] = SpectrogramFreqBinsInHz[k];
				Edge2_Minus->FiftyPrecentIdxs[i] = k;
				flag_50_PrecentIdx_Minus = 1;
				break;
			}
		}
		//		if(Freq_50_PrecentIdx_Minus>Edge2_Minus->maxPeakIdxs[i]){//Freq_50_PrecentIdx_Minus was found
		//			Edge2_Minus->FiftyPrecent[i]=SpectrogramFreqBinsInHz[Freq_50_PrecentIdx_Minus];
		//			Edge2_Minus->FiftyPrecentIdxs[i]=Freq_50_PrecentIdx_Minus;
		//		}
		if (flag_50_PrecentIdx_Minus == 0) {
			Edge2_Minus->FiftyPrecent[i] =
					SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			Edge2_Minus->Peak_Filtered[i] =
					SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			Edge2_Minus->FiftyPrecentIdxs[i] = FgParams->Motion.FminBin;
			Edge2_Minus->maxPeakIdxs[i] = FgParams->Motion.FminBin;
		}
		//	}

		if (fabs(Edge2_Plus->Peak_Filtered[i])
				<= fabs(SpectrogramFreqBinsInHz[FgParams->Motion.FminBin])) {// if the peak has no energy, give the 50% also the same energy
			Edge2_Plus->FiftyPrecent[i] =
					SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			Edge2_Plus->FiftyPrecentIdxs[i] = FgParams->Motion.FminBin;
		}
		if (fabs(Edge2_Minus->Peak_Filtered[i])
				<= fabs(SpectrogramFreqBinsInHz[FgParams->Motion.FminBin])) {// if the peak has no energy, give the 50% also the same energy
			Edge2_Minus->FiftyPrecent[i] =
					SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			Edge2_Minus->FiftyPrecentIdxs[i] = FgParams->Motion.FminBin;
		}
		//		SumEnergy_Post_Plus=0;
		//		SumEnergy_Post_Minus=0;
		for (k = FgParams->Motion.FminBin; k <= Edge2_Plus->maxPeakIdxs[i]; k++) {
			Edge2_Plus->SumEnergy_Post[i] += Pxx2_Plus->Linear[k][i];//sum the energy until black star curve, in linear
		}
		for (k = FgParams->Motion.FminBin; k <= Edge2_Minus->maxPeakIdxs[i]; k++) {
			Edge2_Minus->SumEnergy_Post[i] += Pxx2_Minus->Linear[k][i];	//sum the energy until black star curve, in linear
		}
	}

//			for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
//				printf("##peak idxs %d %d\n", i+1, Edge2_Plus->maxPeakIdxs[i]+1);
//			}

	IsHilbert = 1;
	AvgFilter(Edge2_Plus, FgParams, IsHilbert);//Average filter of the FiftyPrecent(blue) curve with AvgValue=5
	AvgFilter(Edge2_Minus, FgParams, IsHilbert);//Average filter of the FiftyPrecent(blue) curve with AvgValue=5
	///PRINT HERE FIFTY////////////////
	//	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
	//		printf("###fifty FILTERED %d %lf\n", i,
	//				Edge2_Minus->FiftyPrecent_Filtered[i]);
	//	}
	return 0;
}




int HilbertSpectrogram4(float* Pxx2_Hilbert[], float* Pxx2[], float* Pxx2_dB[],
		float *T2_dB, float* Mscan[], int p1, Edge2_Struct *Edge2_Plus,
		Edge2_Struct * Edge2_Minus, Edge2_Struct * Edge2,
		FgParams_Struct* FgParams) {

	double VecForMeanPlus[4],VecForMeanMinus[4],VecForMeanReal[4],MeanNoisePlus,MeanNoiseMinus,MeanNoiseRealPerTimeBin[FgParams->Motion.SpectrogramTimeBins];

	int i, j, m;
	int k;
	int p2 = p1 + FgParams->Motion.SpectrogramRangeWidth - 1;

	int NoiseThreshHilbert = (5 + 6);
	int NoiseThresh = 5;
	float MeanFactor= 1.0/ (FgParams->Motion.SpectrogramRangeWidth);

	_Complex float *MscanIQ[FgParams->Motion.Nscans];

	float SpectrogramPerTimeBin[FgParams->Motion.SpectrogramFreqBinsHilbert];

	Pxx2_Plus_Struct Pxx2_Plus;
	Pxx2_Minus_Struct Pxx2_Minus;

	gsl_complex FFTResultForAbs;
	memset(FgParams->Motion.SignalForFFT, 0,
			sizeof(fftwf_complex) * FgParams->Motion.SpectrogramFreqBinsHilbert);	//was in original


	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 + 1; i++) {//=101
		Pxx2_Plus.Linear[i] = (float *) malloc(
				FgParams->Motion.SpectrogramTimeBins * sizeof(float));
		Pxx2_Plus.dB[i] = (float *) malloc(
				FgParams->Motion.SpectrogramTimeBins * sizeof(float));


		for(int m=0;m < FgParams->Motion.SpectrogramTimeBins; m++ )
			Pxx2_Plus.Linear[i][m]=0;
	}

	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert / 2;i++){
		Pxx2_Minus.Linear[i] = (float *) malloc(
				FgParams->Motion.SpectrogramTimeBins * sizeof(float));
		Pxx2_Minus.dB[i] = (float *) malloc(
				FgParams->Motion.SpectrogramTimeBins * sizeof(float));

		for(int m=0;m < FgParams->Motion.SpectrogramTimeBins; m++ ){
			Pxx2_Minus.Linear[i][m]=0;
			Pxx2[i][m]=0;
		}
	}
	//	gettimeofday(&tpStop,0);


	for (i = 0; i < FgParams->Motion.Nscans; i++) {
		MscanIQ[i] = (_Complex float *) malloc(
				FgParams->Motion.Nbins * sizeof(_Complex float));
	}

	Hilbert(Mscan, MscanIQ, FgParams);	//get the IQ signal with hilbert transform


	//PreProcess


	SlowProcessingHilbert(MscanIQ, MscanIQ, FgParams); //remove DC

	NotchFilterHilbert(MscanIQ, FgParams); //filter the 50&100Hz disturbance that occurred from the fluorescent lamp

	//	gettimeofday(&tpStart,0);



	for (j = p1; j <= p2; j++) {	//relevant bins [p1,p2]

		for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {//take the next 30 points of slow for the FFT

			for (m = 0; m < FgParams->Motion.SpectrogramWinLength; m++) {//windowing with hamming

				FgParams->Motion.SignalForFFT[m] = FgParams->Motion.Hamming[m]
																			  * MscanIQ[m + i * (FgParams->Motion.SpectrogramTimeShift)][j];
			}

			fftwf_execute(FgParams->Motion.FFT_HilbertSpectrogram);	//short time fft on 30 time bins

			//
			for (k = 1; k < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 + 1; k++){//fill the Plus part of the spectrogram
				FFTResultForAbs.dat[0]=FgParams->Motion.FFTResult[k];//real part
				FFTResultForAbs.dat[1]=cimag(FgParams->Motion.FFTResult[k]);
				SpectrogramPerTimeBin[k]=gsl_complex_abs2(FFTResultForAbs);
				Pxx2_Plus.Linear[k][i] += (SpectrogramPerTimeBin[k]);//* MeanFactor;
				Pxx2[k-1][i] +=  (SpectrogramPerTimeBin[k]);
			}

			FFTResultForAbs.dat[0]=FgParams->Motion.FFTResult[0];
			FFTResultForAbs.dat[1]=cimag(FgParams->Motion.FFTResult[0]);
			SpectrogramPerTimeBin[0]=gsl_complex_abs2(FFTResultForAbs);
			Pxx2_Plus.Linear[0][i] += (SpectrogramPerTimeBin[0]);//* MeanFactor;

			for(k = FgParams->Motion.SpectrogramFreqBinsHilbert / 2+1; k < FgParams->Motion.SpectrogramFreqBinsHilbert; k++) {//fill the Minus part of the spectrogram

				//				SpectrogramPerTimeBin[k] = powf(cabsf(FFTResult[k]), 2);
				FFTResultForAbs.dat[0]=FgParams->Motion.FFTResult[k];
				FFTResultForAbs.dat[1]=cimag(FgParams->Motion.FFTResult[k]);
				SpectrogramPerTimeBin[k]=gsl_complex_abs2(FFTResultForAbs);
				//				printf("%d %d %f\n",i,k,SpectrogramPerTimeBin[k]);


				Pxx2_Minus.Linear[FgParams->Motion.SpectrogramFreqBinsHilbert
								  - k][i] += (SpectrogramPerTimeBin[k]);
				Pxx2[FgParams->Motion.SpectrogramFreqBinsHilbert - k][i] +=
						(SpectrogramPerTimeBin[k]);
			}

			Pxx2_Minus.Linear[0][i] = Pxx2_Plus.Linear[0][i];
			Pxx2[0][i] +=  (SpectrogramPerTimeBin[0]);


		}
	}




	//	gettimeofday(&tpStop,0); f1 = ( (float)( tpStop.tv_sec-tpStart.tv_sec)+ (float)(tpStop.tv_usec)/1000000 ) -  ((float)(tpStart.tv_usec)/1000000) ;
	//	printf (" %f sec\n", f1 );

	//caclucaltuation of the mean of the last 4 bins for the noise floor estimation
	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {

		for (k = 0; k < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 + 1; k++) {
			Pxx2_Plus.Linear[k][i]=Pxx2_Plus.Linear[k][i]*MeanFactor;//the mean factor is equal o 1/N=1/120 for the averagring
			Pxx2_Plus.dB[k][i] = 10 * log10(Pxx2_Plus.Linear[k][i]);
		}
		for (k = 0; k < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 ; k++) {
			Pxx2_Minus.Linear[k][i]=Pxx2_Minus.Linear[k][i]*MeanFactor;
			Pxx2_Minus.dB[k][i] = 10 * log10(Pxx2_Minus.Linear[k][i]);
			Pxx2[k][i]=0.5*Pxx2[k][i]*MeanFactor;
			Pxx2_dB[k][i] = 10 * log10(Pxx2[k][i]);
			//		printf("%d %d %f\n",i,k,Pxx2_dB[k][i]);
		}
		//Noise estimation based on the 4 last bins in the spectrogram
		VecForMeanPlus[0]=Pxx2_Plus.Linear[97][i];
		VecForMeanPlus[1]=Pxx2_Plus.Linear[98][i];
		VecForMeanPlus[2]=Pxx2_Plus.Linear[99][i];
		VecForMeanPlus[3]=Pxx2_Plus.Linear[100][i];
		VecForMeanMinus[0]=Pxx2_Minus.Linear[96][i];
		VecForMeanMinus[1]=Pxx2_Minus.Linear[97][i];
		VecForMeanMinus[2]=Pxx2_Minus.Linear[98][i];
		VecForMeanMinus[3]=Pxx2_Minus.Linear[99][i];
		VecForMeanReal[0]=Pxx2[96][i];
		VecForMeanReal[1]=Pxx2[97][i];
		VecForMeanReal[2]=Pxx2[98][i];
		VecForMeanReal[3]=Pxx2[99][i];

		MeanNoisePlus=gsl_stats_mean(VecForMeanPlus,1,4);
		MeanNoiseMinus=gsl_stats_mean(VecForMeanMinus,1,4);
		MeanNoiseRealPerTimeBin[i]=gsl_stats_mean(VecForMeanReal,1,4);
		Edge2_Plus->T1_t[i]=10*log10(MeanNoisePlus)+ NoiseThreshHilbert;
		Edge2_Minus->T1_t[i]=10*log10(MeanNoiseMinus)+ NoiseThreshHilbert;


	}
	//	gettimeofday(&tpStop, 0);

	//	f1 = ((float) (tpStop.tv_sec - tpStart.tv_sec)
	//			+ (float) (tpStop.tv_usec) / 1000000)
	//																					- ((float) (tpStart.tv_usec) / 1000000);
	*T2_dB = 10 * log10(gsl_stats_mean(MeanNoiseRealPerTimeBin,1,FgParams->Motion.SpectrogramTimeBins))+ NoiseThresh;

	CalcCurves(Pxx2, Pxx2_dB, T2_dB, Edge2, FgParams);//CalcRedStarsCurve3 in Matlab

	CalcCurvesHilbert(&Pxx2_Plus, &Pxx2_Minus, Edge2_Plus, Edge2_Minus,
			FgParams);				//CalcRedStarsCurve3_FirstStep2 in Matlab

	//FREE MEMORY
	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 + 1; i++) {//=101
		free(Pxx2_Plus.Linear[i]);
		free(Pxx2_Plus.dB[i]);
	}
	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 ; i++) {
		free(Pxx2_Minus.Linear[i]);
		free(Pxx2_Minus.dB[i]);
	}

	for (i = 0; i < FgParams->Motion.Nscans; i++) {
		free(MscanIQ[i]);
	}


	return 0;
}





int HilbertSpectrogram5(Pxx2_Plus_Struct* Pxx2_Plus
,	Pxx2_Minus_Struct* Pxx2_Minus, float* Pxx2_Hilbert[], float* Pxx2[], float* Pxx2_dB[],
		float *T2_dB, _Complex float* Mscan_PostProcessHilbert[], int p1, Edge2_Struct *Edge2_Plus,
		Edge2_Struct * Edge2_Minus, Edge2_Struct * Edge2,
		FgParams_Struct* FgParams) {

	double VecForMeanPlus[4],VecForMeanMinus[4],VecForMeanReal[4],MeanNoisePlus,MeanNoiseMinus,MeanNoiseRealPerTimeBin[FgParams->Motion.SpectrogramTimeBins];

	int i, j, m;
	int k;
	int p2 = p1 + FgParams->Motion.SpectrogramRangeWidth - 1;

	int NoiseThreshHilbert = (5 + 6);
	int NoiseThresh = 5;
	float MeanFactor= 1.0/ (FgParams->Motion.SpectrogramRangeWidth);

//	_Complex float *MscanIQ[FgParams->Motion.Nscans];

	float SpectrogramPerTimeBin[FgParams->Motion.SpectrogramFreqBinsHilbert];

//	Pxx2_Plus_Struct Pxx2_Plus;
//	Pxx2_Minus_Struct Pxx2_Minus;

	gsl_complex FFTResultForAbs;
	memset(FgParams->Motion.SignalForFFT, 0,
			sizeof(fftwf_complex) * FgParams->Motion.SpectrogramFreqBinsHilbert);	//was in original


	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 + 1; i++) {//=101
		for(int m=0;m < FgParams->Motion.SpectrogramTimeBins; m++ )
			Pxx2_Plus->Linear[i][m]=0;
	}
	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert / 2;i++){

		for(int m=0;m < FgParams->Motion.SpectrogramTimeBins; m++ ){
			Pxx2_Minus->Linear[i][m]=0;
			Pxx2[i][m]=0;
		}
	}



	for (j = p1; j <= p2; j++) {	//relevant bins [p1,p2]

		for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {//take the next 30 points of slow for the FFT

			for (m = 0; m < FgParams->Motion.SpectrogramWinLength; m++) {//windowing with hamming

				FgParams->Motion.SignalForFFT[m] = FgParams->Motion.Hamming[m]
																			  * Mscan_PostProcessHilbert[m + i * (FgParams->Motion.SpectrogramTimeShift)][j];
			}

			fftwf_execute(FgParams->Motion.FFT_HilbertSpectrogram);	//short time fft on 30 time bins

			//
			for (k = 1; k < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 + 1; k++){//fill the Plus part of the spectrogram
				FFTResultForAbs.dat[0]=FgParams->Motion.FFTResult[k];//real part
				FFTResultForAbs.dat[1]=cimag(FgParams->Motion.FFTResult[k]);
				SpectrogramPerTimeBin[k]=gsl_complex_abs2(FFTResultForAbs);
				Pxx2_Plus->Linear[k][i] += (SpectrogramPerTimeBin[k]);//* MeanFactor;
				Pxx2[k-1][i] +=  (SpectrogramPerTimeBin[k]);
			}

			FFTResultForAbs.dat[0]=FgParams->Motion.FFTResult[0];
			FFTResultForAbs.dat[1]=cimag(FgParams->Motion.FFTResult[0]);
			SpectrogramPerTimeBin[0]=gsl_complex_abs2(FFTResultForAbs);
			Pxx2_Plus->Linear[0][i] += (SpectrogramPerTimeBin[0]);//* MeanFactor;

			for(k = FgParams->Motion.SpectrogramFreqBinsHilbert / 2+1; k < FgParams->Motion.SpectrogramFreqBinsHilbert; k++) {//fill the Minus part of the spectrogram

				//				SpectrogramPerTimeBin[k] = powf(cabsf(FFTResult[k]), 2);
				FFTResultForAbs.dat[0]=FgParams->Motion.FFTResult[k];
				FFTResultForAbs.dat[1]=cimag(FgParams->Motion.FFTResult[k]);
				SpectrogramPerTimeBin[k]=gsl_complex_abs2(FFTResultForAbs);
				//				printf("%d %d %f\n",i,k,SpectrogramPerTimeBin[k]);


				Pxx2_Minus->Linear[FgParams->Motion.SpectrogramFreqBinsHilbert
								  - k][i] += (SpectrogramPerTimeBin[k]);
				Pxx2[FgParams->Motion.SpectrogramFreqBinsHilbert - k][i] +=
						(SpectrogramPerTimeBin[k]);
			}

			Pxx2_Minus->Linear[0][i] = Pxx2_Plus->Linear[0][i];
			Pxx2[0][i] +=  (SpectrogramPerTimeBin[0]);


		}
	}




	//	gettimeofday(&tpStop,0); f1 = ( (float)( tpStop.tv_sec-tpStart.tv_sec)+ (float)(tpStop.tv_usec)/1000000 ) -  ((float)(tpStart.tv_usec)/1000000) ;
	//	printf (" %f sec\n", f1 );

	//caclucaltuation of the mean of the last 4 bins for the noise floor estimation
	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {

		for (k = 0; k < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 + 1; k++) {
			Pxx2_Plus->Linear[k][i]=Pxx2_Plus->Linear[k][i]*MeanFactor;//the mean factor is equal o 1/N=1/120 for the averagring
			Pxx2_Plus->dB[k][i] = 10 * log10(Pxx2_Plus->Linear[k][i]);
		}
		for (k = 0; k < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 ; k++) {
			Pxx2_Minus->Linear[k][i]=Pxx2_Minus->Linear[k][i]*MeanFactor;
			Pxx2_Minus->dB[k][i] = 10 * log10(Pxx2_Minus->Linear[k][i]);
			Pxx2[k][i]=0.5*Pxx2[k][i]*MeanFactor;
			Pxx2_dB[k][i] = 10 * log10(Pxx2[k][i]);
			//		printf("%d %d %f\n",i,k,Pxx2_dB[k][i]);
		}
		//Noise estimation based on the 4 last bins in the spectrogram
		VecForMeanPlus[0]=Pxx2_Plus->Linear[97][i];
		VecForMeanPlus[1]=Pxx2_Plus->Linear[98][i];
		VecForMeanPlus[2]=Pxx2_Plus->Linear[99][i];
		VecForMeanPlus[3]=Pxx2_Plus->Linear[100][i];
		VecForMeanMinus[0]=Pxx2_Minus->Linear[96][i];
		VecForMeanMinus[1]=Pxx2_Minus->Linear[97][i];
		VecForMeanMinus[2]=Pxx2_Minus->Linear[98][i];
		VecForMeanMinus[3]=Pxx2_Minus->Linear[99][i];
		VecForMeanReal[0]=Pxx2[96][i];
		VecForMeanReal[1]=Pxx2[97][i];
		VecForMeanReal[2]=Pxx2[98][i];
		VecForMeanReal[3]=Pxx2[99][i];

		MeanNoisePlus=gsl_stats_mean(VecForMeanPlus,1,4);
		MeanNoiseMinus=gsl_stats_mean(VecForMeanMinus,1,4);
		MeanNoiseRealPerTimeBin[i]=gsl_stats_mean(VecForMeanReal,1,4);
		Edge2_Plus->T1_t[i]=10*log10(MeanNoisePlus)+ NoiseThreshHilbert;
		Edge2_Minus->T1_t[i]=10*log10(MeanNoiseMinus)+ NoiseThreshHilbert;


	}
	//	gettimeofday(&tpStop, 0);

	//	f1 = ((float) (tpStop.tv_sec - tpStart.tv_sec)
	//			+ (float) (tpStop.tv_usec) / 1000000)
	//																					- ((float) (tpStart.tv_usec) / 1000000);
	*T2_dB = 10 * log10(gsl_stats_mean(MeanNoiseRealPerTimeBin,1,FgParams->Motion.SpectrogramTimeBins))+ NoiseThresh;

	CalcCurves(Pxx2, Pxx2_dB, T2_dB, Edge2, FgParams);//CalcRedStarsCurve3 in Matlab

	CalcCurvesHilbert(Pxx2_Plus, Pxx2_Minus, Edge2_Plus, Edge2_Minus,
			FgParams);				//CalcRedStarsCurve3_FirstStep2 in Matlab

//	//FREE MEMORY
//	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 + 1; i++) {//=101
//		free(Pxx2_Plus->Linear[i]);
//		free(Pxx2_Plus->dB[i]);
//	}
//	for (i = 0; i < FgParams->Motion.SpectrogramFreqBinsHilbert / 2 ; i++) {
//		free(Pxx2_Minus->Linear[i]);
//		free(Pxx2_Minus->dB[i]);
//	}

//	for (i = 0; i < FgParams->Motion.Nscans; i++) {
//		free(MscanIQ[i]);
//	}


	return 0;
}





















int Hilbert(float* Mscan[], _Complex float* MscanIQ[],
		FgParams_Struct* FgParams) {
	//This function is the implemantion of hilbert function in Matlab that based on the article
	//Marple, S. L. �Computing the Discrete-Time Analytic Signal via FFT.�
	int k, n, i;

	float SignalForFFT[FgParams->Motion.Nbins];
	//	fftwf_complex MscanFFT[FgParams->Motion.Nbins];
	//	fftwf_complex MscanIFFT[FgParams->Motion.Nbins];
	float NbinsFactor=1.0/(FgParams->Motion.Nbins);
	for (k = 1; k < FgParams->Motion.Nbins ; k++)
		FgParams->Motion.MscanFFT[k]=0;
	for (i = 0; i < FgParams->Motion.Nscans; i++) {//implement FFT, erase the negative frequencies and implement IFFT.

		memcpy(SignalForFFT, Mscan[i], FgParams->Motion.Nbins * sizeof(Mscan[0][0]));	//prepare the signal for FFT

		fftwf_execute_dft_r2c(FgParams->Motion.FFT_Hilbert, SignalForFFT, FgParams->Motion.MscanFFT);	//FFT

		for (k = 1; k < FgParams->Motion.Nbins / 2; k++) {
			FgParams->Motion.MscanFFT[k] = 2 * FgParams->Motion.MscanFFT[k]; //in range K=[1,Nscans/2] as the algorithm says and normalize by N
		}

		//		for (k = 0; k < FgParams->Motion.Nbins ; k++) {
		//			printf("%d %f+%f\n",k,creal(FgParams->Motion.MscanFFT[k]),cimag(FgParams->Motion.MscanFFT[k]));
		//
		//		}

		fftwf_execute(FgParams->Motion.IFFT_Hilbert); //IFFT

		for (n = 0; n < FgParams->Motion.Nbins; n++) {
			MscanIQ[i][n] = conj(FgParams->Motion.MscanIFFT[n]) * NbinsFactor;// NbinsFactor=1/Nbins;
		}

	}


	return 0;
}




int MedianFilter(Edge2_Struct *Edge2, FgParams_Struct* FgParams,
		int MedianType) { //Filter the black curve (peak curve)
	int i;

	gsl_vector * input = gsl_vector_alloc(FgParams->Motion.SpectrogramTimeBins);
	gsl_vector * output = gsl_vector_alloc(FgParams->Motion.SpectrogramTimeBins);
	gsl_movstat_workspace * w = gsl_movstat_alloc2(FgParams->Motion.MedianValue / 2,
			FgParams->Motion.MedianValue / 2 - 1);

	if (MedianType == 1) { //TRUNCATE , Computes medians of smaller segments as it reaches the signal edges.

		for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
			gsl_vector_set(input, i, Edge2->Peak[i]);
		}

		gsl_movstat_median(GSL_MOVSTAT_END_TRUNCATE, input, output, w);

		for (i = 0; i < (FgParams->Motion.SpectrogramTimeBins); i++) {
			Edge2->Peak_Filtered[i] = gsl_vector_get(output, i);

		}
	}

	else {// MedianType=0,ZEROPAD, Considers the signal to be zero beyond the endpoints.

		for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
			gsl_vector_set(input, i, Edge2->Peak[i]);			//need +20
		}

		gsl_movstat_median(GSL_MOVSTAT_END_PADZERO, input, output, w);

		for (i = 0; i < (FgParams->Motion.SpectrogramTimeBins); i++) {
			Edge2->Peak_Filtered[i] = gsl_vector_get(output, i);

		}
	}
	gsl_movstat_free(w);
	return 0;
}

int MedianFilterFor50Precent(Edge2_Struct *Edge2, Motion_Struct *MotionStruct,
		float *Edge2_50Precent_MedianFiltered, int MedianValue) {//Filter the 50precentfiltered curve
	int i;
	int SignalLength = MotionStruct->EventStruct->LengthOfChosenEvent + 1;
	gsl_vector * input = gsl_vector_alloc(SignalLength);
	gsl_vector * output = gsl_vector_alloc(SignalLength);
	gsl_movstat_workspace * w = gsl_movstat_alloc2(MedianValue / 2,
			MedianValue / 2 - 1);

	//ZEROPAD, Considers the signal to be zero beyond the endpoints.

	for (i = 0; i < SignalLength; i++) {
		gsl_vector_set(input, i,
				Edge2->FiftyPrecent_Filtered[i
											 + MotionStruct->EventStruct->chosenStart]);	//need +20
	}


	gsl_movstat_median(GSL_MOVSTAT_END_PADZERO, input, output, w);

	for (i = 0; i < (SignalLength); i++) {
		Edge2_50Precent_MedianFiltered[i] = gsl_vector_get(output, i);
		//		printf("%d %f\n",i,Edge2_50Precent_MedianFiltered[i]);
	}

	gsl_movstat_free(w);
	gsl_vector_free(input);
	gsl_vector_free(output);

	return 0;
}


int SlowProcessing(float* Mscan[], float* Mscan_PostProcess[],
		FgParams_Struct* FgParams) {
	int i, j, WindowLength = 10;
	float Alpha = 0.1;
	float AdaptiveDC[FgParams->Motion.Nscans + WindowLength];
	float MeanOf10FirstTimeBins;
	for (j = 0; j < FgParams->Motion.Nbins; j++) {
		MeanOf10FirstTimeBins = 0;
		for (int k = 0; k < WindowLength; k++) {//Calculate mean of 10 first bins for the initial conditions
			MeanOf10FirstTimeBins += Mscan[k][j];
		}
		MeanOf10FirstTimeBins = MeanOf10FirstTimeBins / WindowLength;
		AdaptiveDC[0] = (1 - Alpha) * MeanOf10FirstTimeBins
				+ Alpha * Mscan[0][j];

		for (i = 1; i < FgParams->Motion.Nscans + WindowLength; i++) {
			if (i < FgParams->Motion.Nscans - 1) {
				AdaptiveDC[i] = (1 - Alpha) * AdaptiveDC[i - 1]
														 + Alpha * Mscan[i][j];
			} else {		//if in the last 10 time bins don't change the value
				AdaptiveDC[i] = AdaptiveDC[i - 1];
			}
		}
		for (i = 0; i < FgParams->Motion.Nscans; i++) {
			//			Mscan[i][j] -= AdaptiveDC[i + WindowLength];
			Mscan_PostProcess[i][j] = Mscan[i][j]
											   - AdaptiveDC[i + WindowLength];
			//						printf("%d %d %lf\n",i, j,Mscan[i][j]);

		}
	}

	return 0;
}

int SlowProcessingHilbert(_Complex float* Mscan[],
		_Complex float* Mscan_PostProcess[], FgParams_Struct* FgParams) {
	int i, j, WindowLength = 10;
	_Complex float Alpha = 0.1;
	_Complex float AdaptiveDC[FgParams->Motion.Nscans + WindowLength];
	_Complex float MeanOf10FirstTimeBins;
	for (j = 0; j < FgParams->Motion.Nbins; j++) {
		MeanOf10FirstTimeBins = 0;
		for (int k = 0; k < WindowLength; k++) {//Calculate mean of 10 first bins for the initial conditions
			MeanOf10FirstTimeBins += Mscan[k][j];
		}
		MeanOf10FirstTimeBins = MeanOf10FirstTimeBins / WindowLength;
		AdaptiveDC[0] = (1 - Alpha) * MeanOf10FirstTimeBins
				+ Alpha * Mscan[0][j];

		for (i = 1; i < FgParams->Motion.Nscans + WindowLength; i++) {
			if (i < FgParams->Motion.Nscans - 1) {
				AdaptiveDC[i] = (1 - Alpha) * AdaptiveDC[i - 1]
														 + Alpha * Mscan[i][j];
			} else {		//if in the last 10 time bins don't change the value
				AdaptiveDC[i] = AdaptiveDC[i - 1];
			}
		}
		for (i = 0; i < FgParams->Motion.Nscans; i++) {
			//			Mscan[i][j] -= AdaptiveDC[i + WindowLength];
			Mscan_PostProcess[i][j] = Mscan[i][j]
											   - AdaptiveDC[i + WindowLength];
			//						printf("%d %d %lf\n",i, j,Mscan[i][j]);

		}
	}

	return 0;
}

int NotchFilter(float* Mscan[], FgParams_Struct* FgParams) {

	//first num_concatenate_bins the first 50 bins to the beginning of the Mscan for the initial conditions
	//then filter the slows with notch filter at 50 and 100 Hz.
	int i, j, num_concatenate_scans = 50;
	float concatenated_slow[FgParams->Motion.Nscans + num_concatenate_scans],
	filter_result_50[FgParams->Motion.Nscans + num_concatenate_scans],
	filter_result_100[FgParams->Motion.Nscans + num_concatenate_scans];
	float a_50[3], a_100[3], b_50[3], b_100[3];
	//filter coefficients:
	a_50[0] = 1;
	a_50[1] = -0.546618432826014;
	a_50[2] = 0.768894406379384;
	b_50[0] = 0.884447203189692;
	b_50[1] = -0.546618432826014;
	b_50[2] = 0.884447203189692;
	a_100[0] = 1;
	a_100[1] = 1.43106563601571;
	a_100[2] = 0.768894406379384;
	b_100[0] = 0.884447203189692;
	b_100[1] = 1.43106563601571;
	b_100[2] = 0.884447203189692;

	for (j = 0; j < FgParams->Motion.Nbins; j++) {	//run on the columns
		for (i = 0; i < FgParams->Motion.Nscans + num_concatenate_scans; i++) {
			if (i < num_concatenate_scans) {//concatenate the first 50 nscans
				concatenated_slow[i] = Mscan[i][j];
			} else {
				concatenated_slow[i] = Mscan[i - num_concatenate_scans][j];
			}
			//notch filter @50Hz:
			if (i == 0) {	//initial conditions for the first bins
				filter_result_50[0] = b_50[0] * concatenated_slow[0] / a_50[0];
			} else if (i == 1) {
				filter_result_50[1] = (b_50[0] * concatenated_slow[1]
																   + b_50[1] * concatenated_slow[0]
																								 - a_50[1] * filter_result_50[0]) / a_50[0];
			} else {

				filter_result_50[i] = (b_50[0] * concatenated_slow[i]
																   + b_50[1] * concatenated_slow[i - 1]
																								 + b_50[2] * concatenated_slow[i - 2]
																															   - a_50[1] * filter_result_50[i - 1]
																																							- a_50[2] * filter_result_50[i - 2]) / a_50[0];
			}

			//notch filter @100Hz:
			if (i == 0) {	//initial conditions for the first bins
				filter_result_100[0] = b_100[0] * filter_result_50[0]
																   / a_100[0];
			} else if (i == 1) {
				filter_result_100[1] = (b_100[0] * filter_result_50[1]
																	+ b_100[1] * filter_result_50[0]
																								  - a_100[1] * filter_result_100[0]) / a_100[0];
			} else {

				filter_result_100[i] = (b_100[0] * filter_result_50[i]
																	+ b_100[1] * filter_result_50[i - 1]
																								  + b_100[2] * filter_result_50[i - 2]
																																- a_100[1] * filter_result_100[i - 1]
																																							   - a_100[2] * filter_result_100[i - 2]) / a_100[0];
			}
			//			printf("%d %d %lf\n",i, j,filter_result_100[i]);
			if (i >= num_concatenate_scans) {	//crop the the relevant Nscans
				Mscan[i - num_concatenate_scans][j] = filter_result_100[i];
				//								printf("%d %d %lf \n",i-num_concatenate_scans, j,Mscan[i-num_concatenate_scans][j]);
			}
		}
	}

	return 0;
}

int NotchFilterHilbert(_Complex float* Mscan[], FgParams_Struct* FgParams) {

	//first num_concatenate_bins the first 50 bins to the beginning of the Mscan for the initial conditions
	//then filter the slows with notch filter at 50 and 100 Hz.
	int i, j, num_concatenate_scans = 50;
	_Complex float concatenated_slow[FgParams->Motion.Nscans + num_concatenate_scans],
	filter_result_50[FgParams->Motion.Nscans + num_concatenate_scans],
	filter_result_100[FgParams->Motion.Nscans + num_concatenate_scans];
	_Complex float a_50[3], a_100[3], b_50[3], b_100[3];
	//filter coefficients:
	a_50[0] = 1;
	a_50[1] = -0.546618432826014;
	a_50[2] = 0.768894406379384;
	b_50[0] = 0.884447203189692;
	b_50[1] = -0.546618432826014;
	b_50[2] = 0.884447203189692;
	a_100[0] = 1;
	a_100[1] = 1.43106563601571;
	a_100[2] = 0.768894406379384;
	b_100[0] = 0.884447203189692;
	b_100[1] = 1.43106563601571;
	b_100[2] = 0.884447203189692;

	for (j = 0; j < FgParams->Motion.Nbins; j++) {	//run on the columns
		for (i = 0; i < FgParams->Motion.Nscans + num_concatenate_scans; i++) {
			if (i < num_concatenate_scans) {//concatenate the first 50 nscans
				concatenated_slow[i] = Mscan[i][j];
			} else {
				concatenated_slow[i] = Mscan[i - num_concatenate_scans][j];
			}
			//notch filter @50Hz:
			if (i == 0) {	//initial conditions for the first bins
				filter_result_50[0] = b_50[0] * concatenated_slow[0] / a_50[0];
			} else if (i == 1) {
				filter_result_50[1] = (b_50[0] * concatenated_slow[1]
																   + b_50[1] * concatenated_slow[0]
																								 - a_50[1] * filter_result_50[0]) / a_50[0];
			} else {

				filter_result_50[i] = (b_50[0] * concatenated_slow[i]
																   + b_50[1] * concatenated_slow[i - 1]
																								 + b_50[2] * concatenated_slow[i - 2]
																															   - a_50[1] * filter_result_50[i - 1]
																																							- a_50[2] * filter_result_50[i - 2]) / a_50[0];
			}

			//notch filter @100Hz:
			if (i == 0) {	//initial conditions for the first bins
				filter_result_100[0] = b_100[0] * filter_result_50[0]
																   / a_100[0];
			} else if (i == 1) {
				filter_result_100[1] = (b_100[0] * filter_result_50[1]
																	+ b_100[1] * filter_result_50[0]
																								  - a_100[1] * filter_result_100[0]) / a_100[0];
			} else {

				filter_result_100[i] = (b_100[0] * filter_result_50[i]
																	+ b_100[1] * filter_result_50[i - 1]
																								  + b_100[2] * filter_result_50[i - 2]
																																- a_100[1] * filter_result_100[i - 1]
																																							   - a_100[2] * filter_result_100[i - 2]) / a_100[0];
			}
			//			printf("%d %d %lf\n",i, j,filter_result_100[i]);
			if (i >= num_concatenate_scans) {	//crop the the relevant Nscans
				Mscan[i - num_concatenate_scans][j] = filter_result_100[i];
				//								printf("%d %d %lf \n",i-num_concatenate_scans, j,Mscan[i-num_concatenate_scans][j]);
			}
		}
	}

	return 0;
}



int AbsOfFFT(float* Mscan[], float* Mscan_abs_FFT[],
		FgParams_Struct* FgParams) {
	int k, n, j;
	gsl_complex FFTResultForAbs;




	for (j = 0; j < FgParams->Motion.Nbins; j++) {
		for (n = 0; n < FgParams->Motion.Nscans; n++) {
			FgParams->Motion.SignalForFFTABS[n] = Mscan[n][j];
		}
		fftwf_execute(FgParams->Motion.FFT_ABS);
		for (k = 0; k < FgParams->Motion.DFTLengthForPSD; k++) {//Length/2 because only one side is wanted

			FFTResultForAbs.dat[0]=FgParams->Motion.FFTResultABS[k];
			FFTResultForAbs.dat[1]=cimag(FgParams->Motion.FFTResultABS[k]);
			Mscan_abs_FFT[k][j]=gsl_complex_abs2(FFTResultForAbs);
		}

	}

	return 0;
}








int MatchedFilter(float* Mscan[], FgParams_Struct* FgParams) {

	ne10_fir_instance_f32_t SN; // An FIR "instance structure"
	int FilterLength = 19;
	//	int SignalLength = FgParams->Motion.Nbins + FilterLength - 1;
	ne10_float32_t mf_coeffs[FilterLength];

	ne10_float32_t fir_state_neon[FilterLength + FgParams->Motion.Nbins-1];
	ne10_float32_t FilteredSignal[FgParams->Motion.Nbins];
	ne10_float32_t SignalForFilter[FgParams->Motion.Nbins];

	int i, j;


	//the coeffs are after multipling on InverseFilterGain = 0.599264554413850;

	mf_coeffs[0] = 0.002334734703996;
	mf_coeffs[1] = -0.0082837663;
	mf_coeffs[2] = -0.01349205;
	mf_coeffs[3] = 0.0541153871;
	mf_coeffs[4] = 0.054871207;
	mf_coeffs[5] = -0.2097325325;
	mf_coeffs[6] = -0.0772440025;
	mf_coeffs[7] = 0.4302240491;
	mf_coeffs[8] = 0.0940381154;
	mf_coeffs[9] = -0.5823053675;
	mf_coeffs[10] = -0.0405517524;
	mf_coeffs[11] = 0.534771969;
	mf_coeffs[12] = 0.0988852434;
	mf_coeffs[13] = -0.3118127661;
	mf_coeffs[14] = -0.0337499283;
	mf_coeffs[15] = 0.1374527116;
	mf_coeffs[16] = 0.0051831959;
	mf_coeffs[17] = -0.0245994931;
	mf_coeffs[18] = 0.0002444999;

	//
	//	NE10_SRC_ALLOC (SignalForFilter, guarded_SignalForFilter,FgParams->Motion.Nbins); // 16 extra bytes at the begining and 16 extra bytes at the end
	//	GUARD_ARRAY (SignalForFilter, FgParams->Motion.Nbins);
	//	NE10_DST_ALLOC (FilteredSignal, guarded_FilteredSignal, FgParams->Motion.Nbins);
	//	GUARD_ARRAY (FilteredSignal,FgParams->Motion.Nbins);
	//	NE10_DST_ALLOC (fir_state_neon, guarded_fir_state_neon, FilterLength + FgParams->Motion.Nbins);

	for (i = 0; i < FgParams->Motion.Nscans; i++) {
		for (j = 0; j < FgParams->Motion.Nbins; j++)
			SignalForFilter[j] = Mscan[i][j];

		ne10_fir_init_float(&SN, FilterLength, mf_coeffs, fir_state_neon,
				FgParams->Motion.Nbins);
		ne10_fir_float_c(&SN, SignalForFilter, FilteredSignal,
				FgParams->Motion.Nbins);
		for (j = 0; j < FgParams->Motion.Nbins; j++)
			Mscan[i][j] = FilteredSignal[j];

	}



	return 0;
}



int GET_ROI(float *Mscan_abs_FFT[], FgParams_Struct* FgParams, int *p1) {
	int FrequncyRange = ceil(FgParams->Motion.DFTLengthForPSD / 2);//Rangepowfer is calucalted over this range of frequencies
	float Rangepower[FgParams->Motion.Nbins];
	float MaxValue;
	int LegnthOfSignal=FgParams->Motion.Nbins-FgParams->Motion.SpectrogramRangeWidth+1;
	int i, j;


	ne10_fir_instance_f32_t SN; // An FIR "instance structure"
	int FilterLength = FgParams->Motion.SpectrogramRangeWidth;
	ne10_float32_t fir_coeffs[FgParams->Motion.SpectrogramRangeWidth];

	ne10_float32_t fir_state_neon[FgParams->Motion.SpectrogramRangeWidth + FgParams->Motion.Nbins-1];
	ne10_float32_t FilteredSignal[FgParams->Motion.Nbins];

	for(i=0;i<FgParams->Motion.SpectrogramRangeWidth;i++){
		fir_coeffs[i]=1;
	}

	ne10_fir_init_float(&SN, FilterLength, fir_coeffs, fir_state_neon,
			FgParams->Motion.Nbins);



	for (j = 0; j < FgParams->Motion.Nbins; j++) {
		Rangepower[j] = 0;
		for (i = 0; i <= FrequncyRange; i++) {
			Rangepower[j] += Mscan_abs_FFT[i][j];
		}

	}
	//implement with FIR filter a moving sum over SpectrogramRangeWidth and choose the maximum
	ne10_fir_float_c(&SN, (ne10_float32_t*)Rangepower, FilteredSignal,
			FgParams->Motion.Nbins);


	MaxValue = 0;
	for (i = 0; i < LegnthOfSignal; i++) {
		if (FilteredSignal[i+FgParams->Motion.SpectrogramRangeWidth-1] > MaxValue) {
			MaxValue = FilteredSignal[i+FgParams->Motion.SpectrogramRangeWidth-1];
			*p1 = i;
		}
	}


	return 0;
}



int Spectrogram2(float* Pxx2[], float* Pxx2_dB[], float *T2_dB, float* Mscan[],

		int p1, FgParams_Struct* FgParams) {

	fftwf_complex out[FgParams->Motion.SpectrogramFreqBins];

	float SignalForFFT[FgParams->Motion.DFTLengthForSpectrogram];
	memset(SignalForFFT, 0, sizeof(SignalForFFT));	//was in original

	float SpectrogramPerRangeBin[FgParams->Motion.SpectrogramFreqBins][FgParams->Motion.SpectrogramTimeBins];

	fftw_plan my_plan;
	int flags = 1;
	float fftin_local_ptr_complex[FgParams->Motion.DFTLengthForSpectrogram];
	memset(fftin_local_ptr_complex, 0, sizeof(fftin_local_ptr_complex));

	int i, j, m;
	int k, p2;
	//	float Signal_for_DFT[FgParams->Motion.SpectrogramWinLength];

	int SizeForMean = FgParams->Motion.SpectrogramTimeBins
			* (FgParams->Motion.SpectrogramNoiseFreqBins);
	//	float Spectrogram_per_RangeBin[FgParams->Motion.DFTLengthForSpectrogram / 2 + 1][FgParams->Motion.SpectrogramTimeBins];
	float T2_Linear = 0;	//average of all the spectrogram in linear
	//	_Complex float ResultDFT[FgParams->Motion.DFTLengthForSpectrogram];
	//	_Complex float exp_arg = -I * 2 * M_PI / FgParams->Motion.DFTLengthForSpectrogram;

	p2 = p1 + FgParams->Motion.SpectrogramRangeWidth - 1;
	//	my_plan = fftwf_plan_dft_r2c_2d(FgParams->Motion.DFTLengthForSpectrogram,1,fftin_local_ptr_complex,out,flags);
	my_plan = fftwf_plan_dft_r2c_1d(FgParams->Motion.DFTLengthForSpectrogram,
			fftin_local_ptr_complex, out, flags);

	//	gettimeofday(&tpStart,0);
	for (j = p1; j <= p2; j++) {	//relevant bins [p1,p2]
		for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {//take the next 30 points of slow for the FFT

			for (m = 0; m < FgParams->Motion.SpectrogramWinLength; m++) {//windowing with hamming

				SignalForFFT[m] = FgParams->Motion.Hamming[m]
															* Mscan[m + i * (FgParams->Motion.SpectrogramTimeShift)][j];

			}

			//FFT
			//	memset(out,0,sizeof(out));//was in original
			fftwf_execute_dft_r2c(my_plan, SignalForFFT, out);
			for (k = 0; k < FgParams->Motion.SpectrogramFreqBins; k++) {//Perform short-time fourier transform take in range [0,pi]

				SpectrogramPerRangeBin[k][i] = powf(cabs(out[k]), 2);
				if (j == p1) {//if it's the first spectrogram per range bin, put 0 first
					Pxx2[k][i] = 0;
				}
				Pxx2[k][i] += (SpectrogramPerRangeBin[k][i])
																																								/ (FgParams->Motion.SpectrogramRangeWidth);//Calculate the average for the final spectrogram
			}

		}
	}
	//    gettimeofday(&tpStop,0); f1 = ( (float)( tpStop.tv_sec-tpStart.tv_sec)+ (float)(tpStop.tv_usec)/1000000 ) -  ((float)(tpStart.tv_usec)/1000000); printf (" %f sec\n", f1 );

	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
		//		for(k=SpectrogramFreqBins-SpectrogramNoiseFreqBins;k<SpectrogramFreqBins;k++){
		for (k = 0; k < FgParams->Motion.SpectrogramFreqBins; k++) {

			if ((k >= FgParams->Motion.SpectrogramFreqBins - FgParams->Motion.SpectrogramNoiseFreqBins)
					&& (k < FgParams->Motion.SpectrogramFreqBins)) {
				T2_Linear += Pxx2[k][i];
			}
			Pxx2_dB[k][i] = 10 * log10(Pxx2[k][i]);
			//									printf(" %d %d %lf\n",k,i,Pxx2_dB[k][i]);
		}
	}

	T2_Linear = T2_Linear / SizeForMean;

	*T2_dB = 10 * log10(T2_Linear) + FgParams->Motion.NoiseThresh;//average of all the spectrogram in dB + NoiseThresh

	return 0;
}

int AvgFilter(Edge2_Struct *Edge2, FgParams_Struct* FgParams, int IsHilbert) {
	int SignalLength = FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.AvgValue - 1;
	int i;
	ne10_fir_instance_f32_t SN; // An FIR "instance structure"
	ne10_float32_t SignalForAvg[SignalLength];
	ne10_float32_t FilteredSignal[SignalLength];
	ne10_float32_t coeffs[FgParams->Motion.AvgValue];
	; //1/AvgValue = 1/5HARD CODED
	ne10_float32_t fir_state_neon[FgParams->Motion.AvgValue + SignalLength];
	for (i = 0; i < FgParams->Motion.AvgValue; i++)
		coeffs[i] = 1.0 / FgParams->Motion.AvgValue; //=1/5=0.2 HARD CODED

	if (IsHilbert == 1) { //in the spectrogram we contacte the last avgvalues of the previous motionstruct to the current

		for (i = 0; i <= SignalLength; i++) { //contacte the end of the previous curve
			if (i < FgParams->Motion.AvgValue - 1)
				SignalForAvg[i] = Edge2->PrevLastFiftyPrecent[i];

			else
				SignalForAvg[i] = Edge2->FiftyPrecent[i - FgParams->Motion.AvgValue
													  + 1];
		}

		ne10_fir_init_float(&SN, FgParams->Motion.AvgValue, coeffs, fir_state_neon,
				SignalLength);

		ne10_fir_float_neon(&SN, SignalForAvg, FilteredSignal, SignalLength);
		memcpy(Edge2->FiftyPrecent_Filtered,
				FilteredSignal + (FgParams->Motion.AvgValue - 1),
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(ne10_float32_t));
//				for(i=0;i<FgParams->Motion.SpectrogramTimeBins;i++)
//					printf("50 %d %f\n",i,Edge2->FiftyPrecent[i]);
	}

	else {		//IsHilbert==0
		ne10_fir_init_float(&SN, FgParams->Motion.AvgValue, coeffs, fir_state_neon,
				FgParams->Motion.SpectrogramTimeBins);
		ne10_fir_float_neon(&SN, Edge2->FiftyPrecent,
				Edge2->FiftyPrecent_Filtered, FgParams->Motion.SpectrogramTimeBins);
		//		for(i=0;i<FgParams->Motion.SpectrogramTimeBins;i++)
		//			printf("%d %f\n",i,Edge2->FiftyPrecent_Filtered[i]);
	}

	return 0;
}

int CalcCurves(float* Pxx2[], float* Pxx2_dB[], float *T2_dB,
		Edge2_Struct *Edge2, FgParams_Struct* FgParams) {	//CalcRedStarsCurve3 in Matlab
	int k, i, MaxPxx2_dB_idx, idx_FreqAboveThresh, idx_FreqBins_Spectogram = 0,
			AllFreqAboveThresh[FgParams->Motion.SpectrogramFreqBins];
	int IsHilbert;
	int MedianType;
	//	int AvgValue = 5;
	float freq_increment = 1.256281407035176;	//=(Fs/LengthOfDFT);
	float SpectrogramFreqBinsInHz[FgParams->Motion.DFTLengthForSpectrogram],
	MaxPxx2_dB, Energy_50_Precent;
	for (float i = 0; i < FgParams->Motion.SpectrogramFreqBins; i++) {//create the frequency vector inHz
		SpectrogramFreqBinsInHz[idx_FreqBins_Spectogram] = i * freq_increment;
		idx_FreqBins_Spectogram += 1;
	}
	//	memset(Edge2->Peak,0,174*sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions




	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {

		idx_FreqAboveThresh = 0;
		MaxPxx2_dB = Pxx2_dB[FgParams->Motion.FminBin][i];
		MaxPxx2_dB_idx = FgParams->Motion.FminBin;
		//		MaxPxx2_dB=FminBin;
		//		SumEnergy=0;
		//		if(i==28){
		//			int f=0;
		//			for(int m=0;m<100;m++)
		////				printf("%d %f\n",m,Pxx2_dB[m][i]);
		//		}
		for (k = FgParams->Motion.FminBin;
				k < FgParams->Motion.SpectrogramFreqBins - FgParams->Motion.SpectrogramNoiseFreqBins;
				k++) {//do the work in frequency range [FminBin,SpectrogramFreqBins-SpectrogramNoiseFreqBins]
			if (Pxx2_dB[k][i] > *T2_dB) {
				AllFreqAboveThresh[idx_FreqAboveThresh] = k;
				//				SumEnergy+=Pxx2[k][i];
				idx_FreqAboveThresh += 1;
			}

			if (Pxx2_dB[k][i] > MaxPxx2_dB) {		//find MaxPxx2_dB
				MaxPxx2_dB = Pxx2_dB[k][i];
				MaxPxx2_dB_idx = k;
			}
		}

		//		printf("idx above thr %d %d\n",i,idx_FreqAboveThresh);
		if (idx_FreqAboveThresh > 0) {//there were frequencies above the threshold
			//			Edge2->maxFreqIdxs[i] = AllFreqAboveThresh[idx_FreqAboveThresh - 1];
			//			Edge2->maxFreqEnergy[i] =
			//					Pxx2_dB[AllFreqAboveThresh[idx_FreqAboveThresh - 1]][i];
			//			Edge2->PeakEnergy[i] = MaxPxx2_dB;
			Edge2->maxPeakIdxs[i] = MaxPxx2_dB_idx;
			//			Edge2->Fmax[i] =
			//					SpectrogramFreqBinsInHz[AllFreqAboveThresh[idx_FreqAboveThresh
			//															   - 1]];
			if (FgParams->Motion.TruncateReal == 0) {	//in this case there is zero-pad at the beginning
				Edge2->Peak[i + FgParams->Motion.MedianValue / 2] =
						SpectrogramFreqBinsInHz[MaxPxx2_dB_idx];//first MedianValue/2(=10) bins are 0 for the inital sates of the median filter
				//				printf("peak %d %d %f\n",i,MaxPxx2_dB_idx+1,Edge2->Peak[i]);

			} else {
				//there is no-zero padding
				Edge2->Peak[i] = SpectrogramFreqBinsInHz[MaxPxx2_dB_idx];
				//				printf("peak %d %d %f\n",i+1,MaxPxx2_dB_idx+1,Edge2->Peak[i]);

			}
		} else {		//take defaultes
			//			Edge2->maxFreqIdxs[i] = FgParams->Motion.FminBin;
			//			Edge2->maxFreqEnergy[i] = Pxx2_dB[FgParams->Motion.FminBin][i];
			//			Edge2->PeakEnergy[i] = Pxx2_dB[FgParams->Motion.FminBin][i];
			Edge2->maxPeakIdxs[i] = FgParams->Motion.FminBin;
			//			Edge2->Fmax[i] = SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
			if (FgParams->Motion.TruncateReal == 0) {	//in this case there is zero-pad at the beginning
				Edge2->Peak[i + FgParams->Motion.MedianValue / 2] =
						SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
				//				printf("peak %d  %f\n",i,Edge2->Peak[i]);

			} else
				//there is no-zero padding
				Edge2->Peak[i] = SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
		}

		Energy_50_Precent = (Edge2->PeakEnergy[i] + Edge2->maxFreqEnergy[i])
																																						/ 2;

		//		for (k = Edge2->maxFreqIdxs[i]; k > Edge2->maxPeakIdxs[i]; k--) {//search backwards the first freq index that passes Energy_50_Precent
		//			if (Pxx2_dB[k][i] > Energy_50_Precent) {
		//				break;
		//			}
		//		}
		//		Freq_50_PrecentIdx = k;
		//		//****think about the case: if ~isempty(Freq_50_PrecentIdx)*****
		//		Edge2->FiftyPrecent[i] = SpectrogramFreqBinsInHz[Freq_50_PrecentIdx];
		//		Edge2->FiftyPrecentIdxs[i] = Freq_50_PrecentIdx;

	}
	//	MedianFilter2(Edge2, FgParams->Motion.SpectrogramTimeBins, FgParams->Motion.MedianValue,
	//			FgParams->Motion.truncate);//median filter on the peak(black) curve with MedianValue=20




	MedianType = 0;	//ZEROPOAD
	MedianFilter(Edge2, FgParams, MedianType);	//median filter on the peak(black) curve with MedianValue=20, MedianType=0 ZERPOAD,MedianType=1 TRUNCATE
	//	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
	//	printf("peak pre fil %d  %f\n",i,Edge2->Peak[i]);
	//	}

	//	for (i = 0; i < FgParams->Motion.SpectrogramTimeBins; i++) {
	//
	//		if (fabs(Edge2->Peak_Filtered[i])
	//				<= fabs(SpectrogramFreqBinsInHz[FgParams->Motion.FminBin])) {// if the peak has no energy, give the 50% also the same energy
	//			Edge2->FiftyPrecent[i] =
	//					SpectrogramFreqBinsInHz[FgParams->Motion.FminBin];
	//			Edge2->FiftyPrecentIdxs[i] = FgParams->Motion.FminBin;
	//		}
	//
	//	}

	//NOT RELEVANT DELETE
	//	IsHilbert=0;
	//	AvgFilter(Edge2, FgParams,IsHilbert);//Average filter of the FiftyPrecent(blue) curve with AvgValue=5

	return 0;
}

int MaxOfArr(float *Arr, int *MaxIdx, float *MaxValue, int Length) {
	int i;
	*MaxValue = 0;
	for (i = 0; i < Length; i++) {
		if (Arr[i] > *MaxValue) {
			*MaxValue = Arr[i];
			*MaxIdx = i;
		}
	}
	return 0;
}

void MemoryAllocation(Edge2_Struct* Edge2_1, Edge2_Struct* Edge2_2,
		Edge2_Struct* Edge2_Plus_1, Edge2_Struct* Edge2_Minus_1,
		Edge2_Struct* Edge2_Plus_2, Edge2_Struct* Edge2_Minus_2,
		FgParams_Struct *FgParams) {

	Edge2_1->FiftyPrecent = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_1->FiftyPrecent_Filtered = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_1->Fmax = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_1->PeakEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_1->PeakEnergy_PM = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_1->maxFreqEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_1->FreqBins = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_1->maxFreqIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_1->maxPeakIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_1->FiftyPrecentIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_1->PrevLastFiftyPrecent = (float *) malloc(
			(FgParams->Motion.AvgValue - 1) * sizeof(float));

	if (FgParams->Motion.TruncateReal == 1) {	// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and begining
		Edge2_1->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_1->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_1->Peak,0,(FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_1->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		Edge2_1->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		memset(Edge2_1->Peak, 0,
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Plus_1->FiftyPrecent = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->FiftyPrecent_Filtered = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->Fmax = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->PeakEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->PeakEnergy_PM = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->SumEnergy_Post = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->maxFreqEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_1->FreqBins = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_1->maxFreqIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_1->maxPeakIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_1->FiftyPrecentIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_1->PrevLastFiftyPrecent = (float *) malloc(
			(FgParams->Motion.AvgValue - 1) * sizeof(float));
	Edge2_Plus_1->T1_t = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));

	memset(Edge2_Plus_1->SumEnergy_Post, 0,
			(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));

	if (FgParams->Motion.TruncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Plus_1->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_Plus_1->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Plus_1->Peak,0,(FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Plus_1->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		Edge2_Plus_1->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Plus_1->Peak, 0,
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Minus_1->FiftyPrecent = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->FiftyPrecent_Filtered = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->Fmax = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->PeakEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->PeakEnergy_PM = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->SumEnergy_Post = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->maxFreqEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_1->FreqBins = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_1->maxFreqIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_1->maxPeakIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_1->FiftyPrecentIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_1->PrevLastFiftyPrecent = (float *) malloc(
			(FgParams->Motion.AvgValue - 1) * sizeof(float));
	Edge2_Minus_1->T1_t = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));

	memset(Edge2_Minus_1->SumEnergy_Post, 0,
			(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));

	if (FgParams->Motion.TruncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Minus_1->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_Minus_1->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Minus_1->Peak,0,(FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the initial conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Minus_1->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		Edge2_Minus_1->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Minus_1->Peak, 0,
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_2->FiftyPrecent = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_2->FiftyPrecent_Filtered = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_2->Fmax = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_2->PeakEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_2->PeakEnergy_PM = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_2->maxFreqEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_2->FreqBins = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_2->maxFreqIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_2->maxPeakIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_2->FiftyPrecentIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_2->PrevLastFiftyPrecent = (float *) malloc(
			(FgParams->Motion.AvgValue - 1) * sizeof(float));

	if (FgParams->Motion.TruncateReal == 1) {	// Computes medians of smaller segments as it reaches the signal edges-> no zeropadding at the end and begining
		Edge2_2->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_1->Peak,0,(FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	} else {		//Considers the signal to be zero beyond the endpoints->
		Edge2_2->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		Edge2_2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		memset(Edge2_2->Peak, 0,
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Plus_2->FiftyPrecent = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->FiftyPrecent_Filtered = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->Fmax = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->PeakEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->PeakEnergy_PM = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->SumEnergy_Post = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->maxFreqEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_2->FreqBins = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_2->maxFreqIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_2->maxPeakIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_2->FiftyPrecentIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus_2->PrevLastFiftyPrecent = (float *) malloc(
			(FgParams->Motion.AvgValue - 1) * sizeof(float));
	Edge2_Plus_2->T1_t = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));

	memset(Edge2_Plus_2->SumEnergy_Post, 0,
			(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));

	if (FgParams->Motion.TruncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Plus_2->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_Plus_2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Plus_2->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		Edge2_Plus_2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Plus_2->Peak, 0,
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Minus_2->FiftyPrecent = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->FiftyPrecent_Filtered = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->Fmax = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->PeakEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->PeakEnergy_PM = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->SumEnergy_Post = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->maxFreqEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_2->FreqBins = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_2->maxFreqIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_2->maxPeakIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_2->FiftyPrecentIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus_2->PrevLastFiftyPrecent = (float *) malloc(
			(FgParams->Motion.AvgValue - 1) * sizeof(float));
	Edge2_Minus_2->T1_t = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	memset(Edge2_Minus_2->SumEnergy_Post, 0,
			(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));

	if (FgParams->Motion.TruncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Minus_2->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_Minus_2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Minus_1.Peak,0,(FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the initial conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Minus_2->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		Edge2_Minus_2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Minus_2->Peak, 0,
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

}



void FreeMemory(Edge2_Struct* Edge2_1, Edge2_Struct* Edge2_2,
		Edge2_Struct* Edge2_Plus_1, Edge2_Struct* Edge2_Minus_1,
		Edge2_Struct* Edge2_Plus_2, Edge2_Struct* Edge2_Minus_2,
		Motion_Struct *UnitedMotionStruct) {

	free(Edge2_1->FiftyPrecent);
	free(Edge2_1->FiftyPrecent_Filtered);
	free(Edge2_1->Fmax);
	free(Edge2_1->PeakEnergy);
	free(Edge2_1->PeakEnergy_PM);
	free(Edge2_1->maxFreqEnergy);
	free(Edge2_1->FreqBins);
	free(Edge2_1->maxPeakIdxs);
	free(Edge2_1->maxFreqIdxs);
	free(Edge2_1->FiftyPrecentIdxs);
	free(Edge2_1->PrevLastFiftyPrecent);
	free(Edge2_1->Peak);
	free(Edge2_1->Peak_Filtered);

	free(Edge2_Plus_1->FiftyPrecent);
	free(Edge2_Plus_1->FiftyPrecent_Filtered);
	free(Edge2_Plus_1->Fmax);
	free(Edge2_Plus_1->PeakEnergy);
	free(Edge2_Plus_1->PeakEnergy_PM);
	free(Edge2_Plus_1->maxFreqEnergy);
	free(Edge2_Plus_1->FreqBins);
	free(Edge2_Plus_1->maxPeakIdxs);
	free(Edge2_Plus_1->maxFreqIdxs);
	free(Edge2_Plus_1->FiftyPrecentIdxs);
	free(Edge2_Plus_1->PrevLastFiftyPrecent);
	free(Edge2_Plus_1->Peak);
	free(Edge2_Plus_1->Peak_Filtered);

	free(Edge2_Minus_1->FiftyPrecent);
	free(Edge2_Minus_1->FiftyPrecent_Filtered);
	free(Edge2_Minus_1->Fmax);
	free(Edge2_Minus_1->PeakEnergy);
	free(Edge2_Minus_1->PeakEnergy_PM);
	free(Edge2_Minus_1->maxFreqEnergy);
	free(Edge2_Minus_1->FreqBins);
	free(Edge2_Minus_1->maxPeakIdxs);
	free(Edge2_Minus_1->maxFreqIdxs);
	free(Edge2_Minus_1->FiftyPrecentIdxs);
	free(Edge2_Minus_1->PrevLastFiftyPrecent);
	free(Edge2_Minus_1->Peak);
	free(Edge2_Minus_1->Peak_Filtered);

	free(Edge2_2->FiftyPrecent);
	free(Edge2_2->FiftyPrecent_Filtered);
	free(Edge2_2->Fmax);
	free(Edge2_2->PeakEnergy);
	free(Edge2_2->PeakEnergy_PM);
	free(Edge2_2->maxFreqEnergy);
	free(Edge2_2->FreqBins);
	free(Edge2_2->maxPeakIdxs);
	free(Edge2_2->maxFreqIdxs);
	free(Edge2_2->FiftyPrecentIdxs);
	free(Edge2_2->PrevLastFiftyPrecent);
	free(Edge2_2->Peak);
	free(Edge2_2->Peak_Filtered);

	free(Edge2_Plus_2->FiftyPrecent);
	free(Edge2_Plus_2->FiftyPrecent_Filtered);
	free(Edge2_Plus_2->Fmax);
	free(Edge2_Plus_2->PeakEnergy);
	free(Edge2_Plus_2->PeakEnergy_PM);
	free(Edge2_Plus_2->maxFreqEnergy);
	free(Edge2_Plus_2->FreqBins);
	free(Edge2_Plus_2->maxPeakIdxs);
	free(Edge2_Plus_2->maxFreqIdxs);
	free(Edge2_Plus_2->FiftyPrecentIdxs);
	free(Edge2_Plus_2->PrevLastFiftyPrecent);
	free(Edge2_Plus_2->Peak);
	free(Edge2_Plus_2->Peak_Filtered);

	free(Edge2_Minus_2->FiftyPrecent);
	free(Edge2_Minus_2->FiftyPrecent_Filtered);
	free(Edge2_Minus_2->Fmax);
	free(Edge2_Minus_2->PeakEnergy);
	free(Edge2_Minus_2->PeakEnergy_PM);
	free(Edge2_Minus_2->maxFreqEnergy);
	free(Edge2_Minus_2->FreqBins);
	free(Edge2_Minus_2->maxPeakIdxs);
	free(Edge2_Minus_2->maxFreqIdxs);
	free(Edge2_Minus_2->FiftyPrecentIdxs);
	free(Edge2_Minus_2->PrevLastFiftyPrecent);
	free(Edge2_Minus_2->Peak);
	free(Edge2_Minus_2->Peak_Filtered);

	//	free(UnitedMotionStruct->Edge2->Fmax);
	free(UnitedMotionStruct->Edge2_Plus->Fmax);
	free(UnitedMotionStruct->Edge2_Minus->Fmax);

	free(UnitedMotionStruct->Edge2->Peak_Filtered);
	free(UnitedMotionStruct->Edge2_Plus->Peak_Filtered);
	free(UnitedMotionStruct->Edge2_Minus->Peak_Filtered);

	//	free(UnitedMotionStruct->Edge2->FiftyPrecent_Filtered);
	free(UnitedMotionStruct->Edge2_Plus->FiftyPrecent_Filtered);
	free(UnitedMotionStruct->Edge2_Minus->FiftyPrecent_Filtered);

	free(UnitedMotionStruct->Edge2_Plus->SumEnergy_Post);
	free(UnitedMotionStruct->Edge2_Minus->SumEnergy_Post);

	free(UnitedMotionStruct->Edge2_Plus->T1_t);
	free(UnitedMotionStruct->Edge2_Minus->T1_t);
}

float Max(float Val1, float Val2) {
	if (Val1 > Val2) {
		return Val1;
	} else {
		return Val2;
	}
}

int SaveToCsv(int y_hat, float* MotionDistribution,
		AllFeatures_Struct* FeatureSet, RF_Struct* RF_Model,
		FgParams_Struct *FgParams) {
	float x[6] = { 0 };
	int i;
	char Filename[35];
	sprintf(Filename, "/home/debian/Results/Result%d.csv",
			FgParams->RecordNumber);
	FILE *f = fopen(Filename, "w+");

	fprintf(f, "%d\n", y_hat);
	for (i = 0; i < 4; i++) {
		fprintf(f, "%f\n", MotionDistribution[i]);
	}
	x[0] = (FeatureSet->Plus->Edge2_50Precent_Fmax - RF_Model->x_mean[0])
																																					/ RF_Model->x_std[0];	//#11
	x[1] = (FeatureSet->Minus->Edge2_50Precent_Fmax - RF_Model->x_mean[1])
																																					/ RF_Model->x_std[1];	//#12
	x[2] = (FeatureSet->Plus->Edge2_50Precent_AvgTopFive - RF_Model->x_mean[2])
																																					/ RF_Model->x_std[2];	//#13
	x[3] = (FeatureSet->Minus->Edge2_maxPeakFreq_fromTimeBinWithMaxFreq
			- RF_Model->x_mean[3]) / RF_Model->x_std[3];	//#18
	x[4] = (FeatureSet->Both->HilbertRatio - RF_Model->x_mean[4])
																																					/ RF_Model->x_std[4];	//#21
	x[5] = (FeatureSet->Plus->FmaxFpeakMultiplication - RF_Model->x_mean[5])
																																					/ RF_Model->x_std[5];	//#39

	for (i = 0; i < 6; i++) {
		fprintf(f, "%f\n", x[i]);

	}
	fclose(f);

	return 0;
}

int SavePxx2ToCsv(float* Pxx2[],
		FgParams_Struct *FgParams) {
	char Filename[35];
	sprintf(Filename, "/home/debian/Results/Pxx%d.csv",
			FgParams->RecordNumber);
	FILE *f = fopen(Filename, "w+");

	for(int i=0;i<100;i++){
		for(int j=0;j<75;j++){
			if(j<75-1){

				//					printf("%d %d %f\n",i,j,Pxx2[i][j]);
				fprintf(f, "%f,", Pxx2[i][j]);

			}
			else{
				fprintf(f, "%f\n", Pxx2[i][j]);
				//					printf("%d %d %f\n",i,j,Pxx2[i][j]);
			}
		}

	}

	return 0;
}



int SaveMscanIQToCsv(_Complex float* MscanIQ[],
		FgParams_Struct *FgParams) {
	char FilenameReal[35];
	sprintf(FilenameReal, "/home/debian/Results/MscanIQReal%d.csv",
			FgParams->RecordNumber);
	FILE *fReal = fopen(FilenameReal, "w+");

	char FilenameImag[35];
	sprintf(FilenameImag, "/home/debian/Results/MscanIQImag%d.csv",
			FgParams->RecordNumber);
	FILE *fImag = fopen(FilenameImag, "w+");
	for(int i=0;i<400;i++){
		for(int j=0;j<288;j++){
			if(j<288-1){

				//					printf("%d %d %f\n",i,j,Pxx2[i][j]);
				fprintf(fReal, "%f,", creal(MscanIQ[i][j]));
				fprintf(fImag, "%f,", cimag(MscanIQ[i][j]));

			}
			else{
				fprintf(fReal, "%f\n", creal(MscanIQ[i][j]));
				fprintf(fImag, "%f\n", cimag(MscanIQ[i][j]));

				//					printf("%d %d %f\n",i,j,Pxx2[i][j]);
			}
		}

	}

	return 0;
}

//more features
//		if(MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i]>minfreq){//calculate mean for fifty_precent_filtered without Fmin bins
//			Edge2_noDC_Favg_Plus+=MotionStruct->Edge2_Plus->FiftyPrecent_Filtered[i];
//			num_of_Edge2_noDC_Favg_Plus+=1;
//		}

//if (num_of_Edge2_noDC_Favg_Plus==0)//equivalent to if isempty(Edge2_noDC) in matlab
//	num_of_Edge2_noDC_Favg_Plus=minfreq;
//Edge2_noDC_Favg_Plus=Edge2_noDC_Favg_Plus/(num_of_Edge2_noDC_Favg_Plus);

//	if(Edge2_noDC_Favg_Plus!=0)//otherwise Edge2_PeakToAvg_Plus stay 0
//		Features_Plus.Edge2_50Precent_PeakToAvg=Edge2_AvgTopFive_Plus/Edge2_noDC_Favg_Plus;
