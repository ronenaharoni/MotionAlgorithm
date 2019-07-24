#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "MotionHeader.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include "test.h"
#include "seatest.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "BgXmlParser.h"

struct timeval tpStart , tpStop ;
struct timeval tp ;
float f1;
int i;


FgParams_Struct FgParams;
Tree_Struct Tree[8];
Tree_Struct* All_Trees[8]={&Tree[0],&Tree[1],&Tree[2],&Tree[3],&Tree[4],&Tree[5],&Tree[6],&Tree[7]};
RF_Struct RF_Model;
SVM_Struct SVM_Model;
BgRadarParams bgParams;
MotionXML_Params MotionParamsFromXML;
int main( void ){//init Motion

	//Load configuration from XML//
	char  Filename_ECHO_Configurations[]= "/home/debian/XML_Configurations/ECHO_Configurations.xml";
	parseDoc(Filename_ECHO_Configurations ,&MotionParamsFromXML);

	FgParams.Motion.Nbins=MotionParamsFromXML.Nbins;
	FgParams.Motion.Nscans=MotionParamsFromXML.Nscans; // 400 for 1.6 sec record, 800 in the single 3.2 sec record
	FgParams.Motion.Fs=MotionParamsFromXML.Fs;
	FgParams.Motion.RF_NumOfClasses=MotionParamsFromXML.RF_NumOfClasses;
	FgParams.Motion.SpectrogramRangeWidth=MotionParamsFromXML.SpectrogramRangeWidth;
	FgParams.Motion.SpectrogramWinLength=MotionParamsFromXML.SpectrogramWinLength;
	FgParams.Motion.SpectrogramN_Overlap=MotionParamsFromXML.SpectrogramN_Overlap;
	FgParams.Motion.SpectrogramNoiseFreqBins=MotionParamsFromXML.SpectrogramNoiseFreqBins;
	FgParams.Motion.FminBin=MotionParamsFromXML.FminBin;
	FgParams.Motion.Fbins=MotionParamsFromXML.Fbins;
	FgParams.Motion.TopMaxFreq=MotionParamsFromXML.TopMaxFreq;
	FgParams.Motion.NoiseThresh=MotionParamsFromXML.NoiseThresh;
	FgParams.Motion.MinFreqPlusMinus=MotionParamsFromXML.MinFreqPlusMinus;//in Hz
	FgParams.Motion.MinFreqReal = MotionParamsFromXML.MinFreqReal; //Fmin in Hz
	FgParams.Motion.MinEventDuration=MotionParamsFromXML.MinEventDuration;
	FgParams.Motion.DFTLengthForPSD=floor(FgParams.Motion.Nscans/2)+1;
	FgParams.Motion.DFTLengthForSpectrogram=2*FgParams.Motion.Fbins-1;//=199
	FgParams.Motion.GapLength=MotionParamsFromXML.GapLength;
	FgParams.Motion.DFTLengthForSpectrogramHilbert=2*FgParams.Motion.Fbins;//=200
	FgParams.Motion.MedianValue=MotionParamsFromXML.MedianValue;//for MedianFilter
	FgParams.Motion.AvgValue=MotionParamsFromXML.AvgValue;//for AvgFilter
	FgParams.Motion.TruncateReal=MotionParamsFromXML.TruncateReal;// Computes medians of smaller segments as it reaches the signal edges in Pxx2 (regular) Spectrogram
	FgParams.Motion.TruncateHilbert=MotionParamsFromXML.TruncateHilbert;// Computes medians of smaller segments as it reaches the signal edges in Hilbert Spectrogram
	FgParams.Motion.SpectrogramTimeShift=FgParams.Motion.SpectrogramWinLength-FgParams.Motion.SpectrogramN_Overlap;//=5
	FgParams.Motion.NumSamplesForDerivativeEstimation=MotionParamsFromXML.NumSamplesForDerivativeEstimation;

	FgParams.Motion.SpectrogramFreqBins=FgParams.Motion.DFTLengthForSpectrogram/2+1;//=100
	FgParams.Motion.SpectrogramTimeBinsSingleMotion = floor((FgParams.Motion.Nscans-FgParams.Motion.SpectrogramN_Overlap)/(FgParams.Motion.SpectrogramWinLength - FgParams.Motion.SpectrogramN_Overlap)); // total time bins of the spectrogram=155
	FgParams.Motion.SpectrogramTimeBinsTwoMotions=(2*FgParams.Motion.SpectrogramTimeBinsSingleMotion+FgParams.Motion.GapLength);//=2*75+14=164
	FgParams.Motion.SpectrogramTimeBinsThreeMotions=(3*FgParams.Motion.SpectrogramTimeBinsSingleMotion+2*FgParams.Motion.GapLength);//=3*75+2*14=253
	FgParams.Motion.SpectrogramFreqBinsHilbert=2*FgParams.Motion.Fbins;//=200
	FgParams.Motion.Hamming[0]=0.08; FgParams.Motion.Hamming[1]=0.0907545443733602; FgParams.Motion.Hamming[2]=0.122515306951360;
	FgParams.Motion.Hamming[3]=0.173797189775404; FgParams.Motion.Hamming[4]=0.242202309000359; FgParams.Motion.Hamming[5]=0.324532117278097;
	FgParams.Motion.Hamming[6]=0.416936964276558;FgParams.Motion.Hamming[7]=0.515096102050708; FgParams.Motion.Hamming[8]=0.614419718414272;
	FgParams.Motion.Hamming[9]=0.710263551456361; FgParams.Motion.Hamming[10]=0.798146050066696;FgParams.Motion.Hamming[11]=0.873957926284640;
	FgParams.Motion.Hamming[12]=0.934154301037091; FgParams.Motion.Hamming[13]=0.975920458744089; FgParams.Motion.Hamming[14]=0.997303460291006;
	FgParams.Motion.Hamming[15]=FgParams.Motion.Hamming[14]; FgParams.Motion.Hamming[16]=FgParams.Motion.Hamming[13]; FgParams.Motion.Hamming[17]=FgParams.Motion.Hamming[12]; FgParams.Motion.Hamming[18]=FgParams.Motion.Hamming[11];
	FgParams.Motion.Hamming[19]=FgParams.Motion.Hamming[10]; FgParams.Motion.Hamming[20]=FgParams.Motion.Hamming[9]; FgParams.Motion.Hamming[21]=FgParams.Motion.Hamming[8]; FgParams.Motion.Hamming[22]=FgParams.Motion.Hamming[7];
	FgParams.Motion.Hamming[23]=FgParams.Motion.Hamming[6]; FgParams.Motion.Hamming[24]=FgParams.Motion.Hamming[5]; FgParams.Motion.Hamming[25]=FgParams.Motion.Hamming[4]; FgParams.Motion.Hamming[26]=FgParams.Motion.Hamming[3];
	FgParams.Motion.Hamming[27]=FgParams.Motion.Hamming[2]; FgParams.Motion.Hamming[28]=FgParams.Motion.Hamming[1]; FgParams.Motion.Hamming[29]=FgParams.Motion.Hamming[0];
	FgParams.Motion.A12_Inverse[0][0]=5.925925925925995e-04;FgParams.Motion.A12_Inverse[0][1]=-5.925925925925996e-04;FgParams.Motion.A12_Inverse[0][2]=0.004444444444445;FgParams.Motion.A12_Inverse[0][3]=0.004444444444444;
	FgParams.Motion.A12_Inverse[1][0]=-0.146666666666668;FgParams.Motion.A12_Inverse[1][1]=0.146666666666668;FgParams.Motion.A12_Inverse[1][2]=-1.133333333333347;FgParams.Motion.A12_Inverse[1][3]=-1.066666666666678;
	FgParams.Motion.A12_Inverse[2][0]=12.0000000000001;FgParams.Motion.A12_Inverse[2][1]=-12.0000000000001;FgParams.Motion.A12_Inverse[2][2]=96.0000000000011;FgParams.Motion.A12_Inverse[2][3]=85.0000000000010;
	FgParams.Motion.A12_Inverse[3][0]=-324.000000000004;FgParams.Motion.A12_Inverse[3][1]=325.000000000004;FgParams.Motion.A12_Inverse[3][2]=-2.700000000000030e+03;FgParams.Motion.A12_Inverse[3][3]=-2.250000000000026e+03;


	//////Global Parameters/////

	bgParams.Rmin=1;//need to change this from the configuration XML!!!!
	bgParams.Rmax=3.5; //need to change this from the configuration XML!!!!

	FgParams.Motion.RealRingSizeForMotionTracking=2.62572418212891; //the ring width in meters; FOR DEMO ONLY maybe 2.5??? why 2.6

	bgParams.Rstart_m=2.231486865000000;//for demo only!!! it will be the current Rstart

	int y_hat_M;
	float MotionDistribution[4]; //Num Of classes, satic because dynamic did porblems

	//create Rbin_m
	FgParams.Motion.Rbin_m=(float *)malloc(FgParams.Motion.Nbins * sizeof(float));
	float delta_R=(FgParams.Motion.RealRingSizeForMotionTracking)/((float)FgParams.Motion.Nbins-1);//the delta of the bins inside the ring(=Rstop-Rstart)=2.5m

	for(i=0;i<FgParams.Motion.Nbins;i++){
		FgParams.Motion.Rbin_m[i]=bgParams.Rstart_m+delta_R*i;
	}

	///////Classifiers Preparation////////
	int NumOfTrees=8,NumOfFeatures=6;//RandomForrest
	RF_Params_Import(NumOfTrees, NumOfFeatures,All_Trees,&RF_Model);
	SVM_Params_Import(&SVM_Model);

	//////Create FFT Plans//////
	CreateFFTPlans(&FgParams);


	////////////////// change to 0 in order to go back
#if 0
	Edge2_PlusonStruct0 Preparation////////
	Motion_Struct MotionStruct0;

	Edge2_Struct Edge2_0;
	Edge2_Struct Edge2_Plus_0;
	Edge2_Struct Edge2_Minus_0;
	FgParams.Motion.SpectrogramTimeBins=FgParams.Motion.SpectrogramTimeBinsSingleMotion;

	if(FgParams.Motion.TruncateReal==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and begining
		Edge2_0.Peak=(float *)malloc((FgParams.Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_0.Peak_Filtered=(float *)malloc((FgParams.Motion.SpectrogramTimeBins) * sizeof(float));
	}
	else{//Considers the signal to be zero beyond the endpoints.
		Edge2_0.Peak_Filtered=(float *)malloc((FgParams.Motion.SpectrogramTimeBins+FgParams.Motion.MedianValue-1) * sizeof(float));
	}

	Edge2_Plus_0.FiftyPrecent_Filtered=(float *)malloc(FgParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_0.Fmax=(float *)malloc(FgParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_0.PrevLastFiftyPrecent=(float *)malloc((FgParams.Motion.AvgValue-1) * sizeof(float));
	Edge2_Plus_0.SumEnergy_Post=(float *)malloc(FgParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_0.T1_t=(float *)malloc(FgParams.Motion.SpectrogramTimeBins * sizeof(float));

	if(FgParams.Motion.TruncateHilbert==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Plus_0.Peak_Filtered=(float *)malloc((FgParams.Motion.SpectrogramTimeBins) * sizeof(float));
	}
	else{//Considers the signal to be zero beyond the endpoints.
		Edge2_Plus_0.Peak_Filtered=(float *)malloc((FgParams.Motion.SpectrogramTimeBins+FgParams.Motion.MedianValue-1) * sizeof(float));
	}

	Edge2_Minus_0.FiftyPrecent_Filtered=(float *)malloc(FgParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_0.Fmax=(float *)malloc(FgParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_0.PrevLastFiftyPrecent=(float *)malloc((FgParams.Motion.AvgValue-1) * sizeof(float));
	Edge2_Minus_0.SumEnergy_Post=(float *)malloc(FgParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_0.T1_t=(float *)malloc(FgParams.Motion.SpectrogramTimeBins * sizeof(float));

	if(FgParams.Motion.TruncateHilbert==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Minus_0.Peak_Filtered=(float *)malloc((FgParams.Motion.SpectrogramTimeBins) * sizeof(float));
	}
	else{//Considers the signal to be zero beyond the endpoints.
		Edge2_Minus_0.Peak_Filtered=(float *)malloc((FgParams.Motion.SpectrogramTimeBins+FgParams.Motion.MedianValue-1) * sizeof(float));
	}

	MotionStruct0.Edge2=&Edge2_0;
	MotionStruct0.Edge2_Plus=&Edge2_Plus_0;
	MotionStruct0.Edge2_Minus=&Edge2_Minus_0;

#endif

	Motion_Struct MotionStruct0;
	Motion_Struct MotionStruct1;
	Motion_Struct MotionStruct2;
	Edge2_Struct Edge2_0;
	Edge2_Struct Edge2_Plus_0;
	Edge2_Struct Edge2_Minus_0;
	Edge2_Struct Edge2_1;
	Edge2_Struct Edge2_Plus_1;
	Edge2_Struct Edge2_Minus_1;
	Edge2_Struct Edge2_2;
	Edge2_Struct Edge2_Plus_2;
	Edge2_Struct Edge2_Minus_2;

	AllocateMemoryForEdges(&Edge2_0, &Edge2_Plus_0,&Edge2_Minus_0,&FgParams);
	AllocateMemoryForEdges(&Edge2_1, &Edge2_Plus_1,&Edge2_Minus_1,&FgParams);
	AllocateMemoryForEdges(&Edge2_2, &Edge2_Plus_2,&Edge2_Minus_2,&FgParams);
	MotionStruct0.Edge2=&Edge2_0;
	MotionStruct0.Edge2_Plus=&Edge2_Plus_0;
	MotionStruct0.Edge2_Minus=&Edge2_Minus_0;
	MotionStruct1.Edge2=&Edge2_1;
	MotionStruct1.Edge2_Plus=&Edge2_Plus_1;
	MotionStruct1.Edge2_Minus=&Edge2_Minus_1;
	MotionStruct2.Edge2=&Edge2_2;
	MotionStruct2.Edge2_Plus=&Edge2_Plus_2;
	MotionStruct2.Edge2_Minus=&Edge2_Minus_2;



	float* Mscan_abs_FFT[FgParams.Motion.DFTLengthForPSD];
	float* Pxx2_dB[FgParams.Motion.SpectrogramFreqBins],
	*Pxx2[FgParams.Motion.SpectrogramFreqBins];
	float * Pxx2_Hilbert[FgParams.Motion.SpectrogramFreqBinsHilbert];


	for (i = 0; i < FgParams.Motion.DFTLengthForPSD; i++) {
		Mscan_abs_FFT[i] = (float *) malloc(FgParams.Motion.Nbins * sizeof(float));
	}
	for (i = 0; i < FgParams.Motion.SpectrogramFreqBins; i++) {
		Pxx2_dB[i] = (float *) malloc(
				FgParams.Motion.SpectrogramTimeBins * sizeof(float));
		Pxx2[i] = (float *) malloc(
				FgParams.Motion.SpectrogramTimeBins * sizeof(float));
	}
	for (i = 0; i < FgParams.Motion.SpectrogramFreqBinsHilbert; i++) {

		Pxx2_Hilbert[i] = (float *) malloc(
				FgParams.Motion.SpectrogramTimeBins * sizeof(float));
	}

	float* Mscan_PostProcessReal[FgParams.Motion.Nscans];

	float* Mscan_PostProcess0[FgParams.Motion.Nscans];
	float* Mscan_PostProcess1[FgParams.Motion.Nscans];

	for (i = 0; i < FgParams.Motion.Nscans; i++)
		Mscan_PostProcess0[i] = (float *) malloc(
				FgParams.Motion.Nbins * sizeof(float));
	for (i = 0; i < FgParams.Motion.Nscans; i++)
		Mscan_PostProcess1[i] = (float *) malloc(
				FgParams.Motion.Nbins * sizeof(float));
	for (i = 0; i < FgParams.Motion.Nscans; i++)
		Mscan_PostProcessReal[i] = (float *) malloc(
				FgParams.Motion.Nbins * sizeof(float));

	Pxx2_Plus_Struct Pxx2_Plus;
	Pxx2_Minus_Struct Pxx2_Minus;

	_Complex float *Mscan_PostProcessHilbert[FgParams.Motion.Nscans];


	for (i = 0; i < FgParams.Motion.SpectrogramFreqBinsHilbert / 2 + 1; i++) {//=101
		Pxx2_Plus.Linear[i] = (float *) malloc(
				FgParams.Motion.SpectrogramTimeBins * sizeof(float));
		Pxx2_Plus.dB[i] = (float *) malloc(
				FgParams.Motion.SpectrogramTimeBins * sizeof(float));

	}

	for (i = 0; i < FgParams.Motion.SpectrogramFreqBinsHilbert / 2;i++){
		Pxx2_Minus.Linear[i] = (float *) malloc(
				FgParams.Motion.SpectrogramTimeBins * sizeof(float));
		Pxx2_Minus.dB[i] = (float *) malloc(
				FgParams.Motion.SpectrogramTimeBins * sizeof(float));

	}

	for (i = 0; i < FgParams.Motion.Nscans; i++) {
		Mscan_PostProcessHilbert[i] = (_Complex float *) malloc(
				FgParams.Motion.Nbins * sizeof(_Complex float));
	}


	///////////////END OF INIT//////////////////////



	////////GET MSCAN0 AND MSCAN1////////

	int j;


	//
	//	float *Mscan_PostProcess3;
	//	Mscan_PostProcess3=(float*)malloc(FgParams.Motion.Nbins *FgParams.Motion.Nscans* sizeof(float));
	///
	float Mscan_flat[FgParams.Motion.Nscans*FgParams.Motion.Nbins];
	float* Mscan0[FgParams.Motion.Nscans],* Mscan1[FgParams.Motion.Nscans],*Mscan2[FgParams.Motion.Nscans],* Mscan3[FgParams.Motion.Nscans];//,*Mscan_old[FgParams.Motion.Nscans];
	for (i=0; i<FgParams.Motion.Nscans; i++){
		Mscan0[i] = (float *)malloc(FgParams.Motion.Nbins * sizeof(float));
		Mscan1[i] = (float *)malloc(FgParams.Motion.Nbins * sizeof(float));
		Mscan2[i] = (float *)malloc(FgParams.Motion.Nbins * sizeof(float));
		Mscan3[i] = (float *)malloc(FgParams.Motion.Nbins * sizeof(float));

	}


	for(FgParams.RecordNumber=700;FgParams.RecordNumber<701;FgParams.RecordNumber++){
		printf("record number %d\n",FgParams.RecordNumber);
		char  Filename_Mscan1[50] ;
		sprintf(Filename_Mscan1, "/home/debian/Records/Mscan_withoutMF%d_1.csv" , FgParams.RecordNumber);
		FILE *MscanFile1=fopen(Filename_Mscan1, "r");

		char  Filename_Mscan2[50] ;
		sprintf(Filename_Mscan2, "/home/debian/Records/Mscan_withoutMF%d_2.csv" , FgParams.RecordNumber);
		FILE *MscanFile2=fopen(Filename_Mscan2, "r");

		char  Filename_Mscan3[50] ;
		sprintf(Filename_Mscan3, "/home/debian/Records/Mscan_withoutMF%d_3.csv" , FgParams.RecordNumber);
		FILE *MscanFile3=fopen(Filename_Mscan3, "r");

		read_data_from_file(MscanFile1,FgParams.Motion.Nscans,FgParams.Motion.Nbins,Mscan_flat);//get the MSCAN in flat form

		for ( i=0;i<FgParams.Motion.Nscans;i++)
		{//get the  MSCAN from flat
			for ( j=0;j<FgParams.Motion.Nbins;j++)
			{
				Mscan0[i][j]=Mscan_flat[i*FgParams.Motion.Nbins+j];
			}
		}

		read_data_from_file(MscanFile1,FgParams.Motion.Nscans,FgParams.Motion.Nbins,Mscan_flat);//get the MSCAN in flat form
		for ( i=0;i<FgParams.Motion.Nscans;i++)
		{//get the  MSCAN from flat
			for ( j=0;j<FgParams.Motion.Nbins;j++)
			{
				Mscan1[i][j]=Mscan_flat[i*FgParams.Motion.Nbins+j];
			}
		}

		read_data_from_file(MscanFile2,FgParams.Motion.Nscans,FgParams.Motion.Nbins,Mscan_flat);//get the MSCAN in flat form
		for ( i=0;i<FgParams.Motion.Nscans;i++)
		{//get the  MSCAN from flat
			for ( j=0;j<FgParams.Motion.Nbins;j++)
			{
				Mscan2[i][j]=Mscan_flat[i*FgParams.Motion.Nbins+j];
			}
		}


		read_data_from_file(MscanFile3,FgParams.Motion.Nscans,FgParams.Motion.Nbins,Mscan_flat);//get the MSCAN in flat form
		for ( i=0;i<FgParams.Motion.Nscans;i++)
		{//get the  MSCAN from flat
			for ( j=0;j<FgParams.Motion.Nbins;j++)
			{
				Mscan3[i][j]=Mscan_flat[i*FgParams.Motion.Nbins+j];
			}
		}



		gettimeofday(&tpStart,0);


		//for each Mscan: PreProcess->MotionTracking->Curveecraction and then for two mscans Analyze and classify the motion and save Motionstrct2 in Motionstruct0 for the next cycle
		FgParams.Motion.FirstTimeMotion=1;//1 if it's the first 2 frames
		FgParams.Motion.SpectrogramTimeBins=FgParams.Motion.SpectrogramTimeBinsSingleMotion;//=74



		MotionPreprocess(Mscan0,Mscan_PostProcessReal,Mscan_PostProcessHilbert,&FgParams);
		MotionTracking(Mscan_PostProcessReal, &FgParams,&bgParams);//Send Mscan0 to tracking//
		MotionCurveExtractionPerMscan2(&MotionStruct1,&Pxx2_Plus,&Pxx2_Minus, Mscan_abs_FFT,Mscan_PostProcessReal,Mscan_PostProcessHilbert, Pxx2_Hilbert,
				Pxx2, Pxx2_dB, &FgParams);


		//save the last AvgValue-1 bins of FiftyPrecent of motionstruct1 for the initial state of the AvgFilter of mostionstruct2
		for (i = FgParams.Motion.AvgValue - 1; i > 0; i--) { //it's done only to Plus and Minus
			MotionStruct2.Edge2_Plus->PrevLastFiftyPrecent[(FgParams.Motion.AvgValue - 1) - i] =
					Edge2_Plus_1.FiftyPrecent[FgParams.Motion.SpectrogramTimeBins - i];
			MotionStruct2.Edge2_Minus->PrevLastFiftyPrecent[(FgParams.Motion.AvgValue - 1) - i] =
					Edge2_Minus_1.FiftyPrecent[FgParams.Motion.SpectrogramTimeBins - i];
		}


		MotionPreprocess(Mscan1,Mscan_PostProcessReal,Mscan_PostProcessHilbert,&FgParams);
		MotionTracking(Mscan_PostProcessReal, &FgParams,&bgParams);//Send Mscan1 to tracking//
		MotionCurveExtractionPerMscan2(&MotionStruct2,&Pxx2_Plus,&Pxx2_Minus, Mscan_abs_FFT,Mscan_PostProcessReal,Mscan_PostProcessHilbert, Pxx2_Hilbert,
				Pxx2, Pxx2_dB, &FgParams);

		memset(MotionDistribution, 0, FgParams.Motion.RF_NumOfClasses * sizeof(float));	//initialize with zeros

		MotionAnalyzerTwoMscans(&MotionStruct0,&MotionStruct1,&MotionStruct2,&FgParams,All_Trees,&SVM_Model,&RF_Model,&y_hat_M,MotionDistribution);

		gettimeofday(&tpStop,0); f1 = ( (float)( tpStop.tv_sec-tpStart.tv_sec)+ (float)(tpStop.tv_usec)/1000000 ) -  ((float)(tpStart.tv_usec)/1000000); printf (" %f sec\n", f1 );



		//second loop for demo
		FgParams.Motion.FirstTimeMotion=0;
		FgParams.Motion.SpectrogramTimeBins=FgParams.Motion.SpectrogramTimeBinsSingleMotion;//=74



		MotionPreprocess(Mscan2,Mscan_PostProcessReal,Mscan_PostProcessHilbert,&FgParams);
		MotionTracking(Mscan_PostProcessReal, &FgParams,&bgParams);//Send Mscan0 to tracking//
		MotionCurveExtractionPerMscan2(&MotionStruct1,&Pxx2_Plus,&Pxx2_Minus, Mscan_abs_FFT,Mscan_PostProcessReal,Mscan_PostProcessHilbert, Pxx2_Hilbert,
				Pxx2, Pxx2_dB, &FgParams);

		//save the last AvgValue-1 bins of FiftyPrecent of motionstruct1 for the initial state of the AvgFilter of mostionstruct2

		for (i = FgParams.Motion.AvgValue - 1; i > 0; i--) { //it's done only to Plus and Minus
			MotionStruct2.Edge2_Plus->PrevLastFiftyPrecent[(FgParams.Motion.AvgValue - 1) - i] =
					Edge2_Plus_1.FiftyPrecent[FgParams.Motion.SpectrogramTimeBins - i];
			MotionStruct2.Edge2_Minus->PrevLastFiftyPrecent[(FgParams.Motion.AvgValue - 1) - i] =
					Edge2_Minus_1.FiftyPrecent[FgParams.Motion.SpectrogramTimeBins - i];
		}


		MotionPreprocess(Mscan3,Mscan_PostProcessReal,Mscan_PostProcessHilbert,&FgParams);
		MotionTracking(Mscan_PostProcessReal, &FgParams,&bgParams);//Send Mscan1 to tracking//
		MotionCurveExtractionPerMscan2(&MotionStruct2,&Pxx2_Plus,&Pxx2_Minus, Mscan_abs_FFT,Mscan_PostProcessReal,Mscan_PostProcessHilbert, Pxx2_Hilbert,
				Pxx2, Pxx2_dB, &FgParams);

		memset(MotionDistribution, 0, FgParams.Motion.RF_NumOfClasses * sizeof(float));	//initialize with zeros

		MotionAnalyzerTwoMscans(&MotionStruct0,&MotionStruct1,&MotionStruct2,&FgParams,All_Trees,&SVM_Model,&RF_Model,&y_hat_M,MotionDistribution);



		for(i=0;i<4;i++){

			printf("%f\n",MotionDistribution[i]);
		}

		//free memory

	}

	for (i=0; i<FgParams.Motion.Nscans; i++){
		free(Mscan0[i]);
		free(Mscan1[i]);
		free(Mscan2[i]);
		free(Mscan3[i]);
	}
	fftwf_destroy_plan(FgParams.Motion.FFT_Hilbert);
	fftwf_destroy_plan(FgParams.Motion.IFFT_Hilbert);
	fftwf_destroy_plan(FgParams.Motion.FFT_HilbertSpectrogram);
	fftwf_destroy_plan(FgParams.Motion.FFT_ABS);
	ne10_fft_destroy_r2c_float32(FgParams.Motion.FFT_Feature42);


	//	FreeMemoryEdeges(&Edge2_0,
	//			&Edge2_Plus_0,&Edge2_Minus_0);//think about it again

	return 0;
}


int read_data_from_file(FILE *fp_read ,int Nscans,int Nbins,float* RespMscan_flat)
{
	int i;


	for (i=0;i<Nscans*Nbins;i++)
	{
		fscanf(fp_read,"%f",RespMscan_flat+i);
		//				printf("%d %lf\n",i,*(RespMscan_flat+i));
	}
	return 0;
}


int RF_Params_Import(int NumOfTrees,int NumOfFeatures, Tree_Struct** All_Trees,RF_Struct* RF_Model){
	FILE *x_mean_File;
	FILE *x_std_File;
	double Value0;
	int i;
	x_mean_File=fopen("/home/debian/RandomForestClassifier/RF_x_mean.csv", "r");
	x_std_File=fopen("/home/debian/RandomForestClassifier/RF_x_std.csv", "r");

	for(i=0;i<NumOfFeatures;i++){//import the mean and std for the normalization
		fscanf(x_mean_File,"%lf",&Value0);
		RF_Model->x_mean[i]=Value0;
		fscanf(x_std_File,"%lf",&Value0);
		RF_Model->x_std[i]=Value0;
	}

	TreesCreator(NumOfTrees, All_Trees);//create 8 trees of the classifier
	return 0;
}

int TreesCreator(int NumOfTrees, Tree_Struct** All_Trees){

	//notice: some variabels are hard coded
	char  Filename_Children[50] ;
	char Filename_CutPoint[50];
	char Filename_CutPredictor[54];
	char Filename_ClassProb[50];
	FILE *Children_File;
	FILE *CutPoint_File;
	FILE *CutPredictor_File;
	FILE *ClassProb_File;
	FILE *Lengths_File;
	int TreeNum,i;
	int Current_Length;//the length of each tree
	float Value0,Value1,Value2,Value3;
	int ValueCutPredictor;
	int Value0Children,Value1Children;
	Lengths_File=fopen("/home/debian/RandomForestClassifier/Lengths.csv", "r");
	All_Trees[0]->TotalTrees=8;
	for(TreeNum=0;TreeNum<All_Trees[0]->TotalTrees;TreeNum++){
		fscanf(Lengths_File,"%d",&Current_Length);
		memset(Filename_Children , 0 , sizeof(Filename_Children));
		memset(Filename_CutPoint , 0 , sizeof(Filename_CutPoint));
		memset(Filename_CutPredictor , 0 , sizeof(Filename_CutPredictor));
		memset(Filename_ClassProb , 0 , sizeof(Filename_ClassProb));

		sprintf(Filename_CutPredictor, "/home/debian/RandomForestClassifier/CutPredictor%d.csv" , TreeNum);
		CutPredictor_File=fopen(Filename_CutPredictor, "r");


		sprintf(Filename_CutPoint, "/home/debian/RandomForestClassifier/CutPoint%d.csv" , TreeNum);
		CutPoint_File=fopen(Filename_CutPoint, "r");


		sprintf(Filename_Children, "/home/debian/RandomForestClassifier/Children%d.csv" , TreeNum);
		Children_File=fopen(Filename_Children, "r");

		sprintf(Filename_ClassProb, "/home/debian/RandomForestClassifier/ClassProb%d.csv" , TreeNum);
		ClassProb_File=fopen(Filename_ClassProb, "r");


		All_Trees[TreeNum]->CutPoint =(float *)malloc(Current_Length * sizeof(float));
		All_Trees[TreeNum]->CutPredictor =(int *)malloc(Current_Length * sizeof(int));
		All_Trees[TreeNum]->Children=(int *)malloc(2 *Current_Length* sizeof(int));
		All_Trees[TreeNum]->ClassProb=(float *)malloc(4 *Current_Length* sizeof(float));

		for(i=0;i<Current_Length;i++){
			fscanf(CutPoint_File,"%f",&Value0);//Value=-1 is for NaN
			//     		printf("%f cut\n",Value0);
			All_Trees[TreeNum]->CutPoint[i]=Value0;
			fscanf(CutPredictor_File,"%d",&ValueCutPredictor);
			All_Trees[TreeNum]->CutPredictor[i]=ValueCutPredictor;//Value=-1 is for NaN
			//			     		printf("%d cutp\n",ValueCutPredictor);
			fscanf(Children_File,"%d,%d",&Value0Children,&Value1Children);

			All_Trees[TreeNum]->Children[2*i]=Value0Children;//children is 2 columns left & right
			All_Trees[TreeNum]->Children[2*i+1]=Value1Children;

			fscanf(ClassProb_File,"%f,%f,%f,%f",&Value0,&Value1,&Value2,&Value3);

			All_Trees[TreeNum]->ClassProb[4*i]=Value0;//classprob is 4 columns because 4 classes
			All_Trees[TreeNum]->ClassProb[4*i+1]=Value1;
			All_Trees[TreeNum]->ClassProb[4*i+2]=Value2;
			All_Trees[TreeNum]->ClassProb[4*i+3]=Value3;
			//			printf("%f,%f,%f,%f class\n",Value0,Value1,Value2,Value3);


		}
	}
	return 0;
}


int SVM_Params_Import(SVM_Struct* SVM_Model){
	FILE *Bias_File;
	FILE *Beta_File;
	FILE *x_mean_File;
	FILE *x_std_File;
	float Value0;
	int i,NumOfFeatures=6;

	Bias_File=fopen("/home/debian/SVMClassifier/SVM_Bias.csv", "r");
	Beta_File=fopen("/home/debian/SVMClassifier/SVM_Beta.csv", "r");
	x_mean_File=fopen("/home/debian/SVMClassifier/SVM_x_mean.csv", "r");
	x_std_File=fopen("/home/debian/SVMClassifier/SVM_x_std.csv", "r");

	for(i=0;i<NumOfFeatures;i++){
		fscanf(Beta_File,"%f",&Value0);
		SVM_Model->Beta[i]=Value0;
		fscanf(x_mean_File,"%f",&Value0);
		SVM_Model->x_mean[i]=Value0;
		fscanf(x_std_File,"%f",&Value0);
		SVM_Model->x_std[i]=Value0;
	}
	fscanf(Bias_File,"%f",&Value0);
	SVM_Model->Bias=Value0;

	return 0;
}

int CreateFFTPlans(FgParams_Struct *FgParams){


	int flags = 1;
	float fftin_local_ptr_complex[FgParams->Motion.Nbins];
	memset(fftin_local_ptr_complex, 0, sizeof(fftin_local_ptr_complex));
	FgParams->Motion.FFT_Hilbert = fftwf_plan_dft_r2c_1d(FgParams->Motion.Nbins,
			fftin_local_ptr_complex, FgParams->Motion.MscanFFT, flags);

	FgParams->Motion.IFFT_Hilbert = fftwf_plan_dft_1d(FgParams->Motion.Nbins, FgParams->Motion.MscanFFT, FgParams->Motion.MscanIFFT,
			FFTW_BACKWARD, FFTW_MEASURE);	//MAYBE ESTIMATEEEEE!!!!!




	FgParams->Motion.FFT_HilbertSpectrogram = fftwf_plan_dft_r2c_1d(FgParams->Motion.Nbins,
			fftin_local_ptr_complex, FgParams->Motion.MscanFFT, flags);


	FgParams->Motion.SignalForFFT = (fftwf_complex*) fftwf_malloc( sizeof(fftwf_complex) * FgParams->Motion.SpectrogramFreqBinsHilbert);	//pay attention

	FgParams->Motion.SignalForFFTABS = (fftwf_complex*) fftwf_malloc( sizeof(fftwf_complex) * FgParams->Motion.Nscans);	//pay attention



	FgParams->Motion.FFTResult = (fftwf_complex*) fftwf_malloc(
			sizeof(fftwf_complex) * FgParams->Motion.SpectrogramFreqBinsHilbert);


	FgParams->Motion.FFTResultABS = (fftwf_complex*) fftwf_malloc(
			sizeof(fftwf_complex) * (FgParams->Motion.Nscans));

	FgParams->Motion.FFT_HilbertSpectrogram=fftwf_plan_dft_1d(FgParams->Motion.SpectrogramFreqBinsHilbert,
			FgParams->Motion.SignalForFFT, FgParams->Motion.FFTResult, FFTW_FORWARD, FFTW_ESTIMATE);

	FgParams->Motion.FFT_ABS= fftwf_plan_dft_1d(FgParams->Motion.Nscans,
			FgParams->Motion.SignalForFFTABS, FgParams->Motion.FFTResultABS, FFTW_FORWARD, FFTW_ESTIMATE);

	// An FFT "configuration structure"
	FgParams->Motion.FFT_Feature42 = ne10_fft_alloc_r2c_float32(256);

	return 0;
}




int AllocateMemoryForEdges(Edge2_Struct* Edge2, Edge2_Struct* Edge2_Plus,Edge2_Struct* Edge2_Minus,FgParams_Struct *FgParams){


	FgParams->Motion.SpectrogramTimeBins=FgParams->Motion.SpectrogramTimeBinsSingleMotion;

	//	if(FgParams->Motion.TruncateReal==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and begining
	//		Edge2->Peak=(float *)malloc((FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
	//		Edge2->Peak_Filtered=(float *)malloc((FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
	//	}
	//	else{//Considers the signal to be zero beyond the endpoints.
	//		Edge2->Peak_Filtered=(float *)malloc((FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue-1) * sizeof(float));
	//	}
	//
	//	Edge2_Plus->FiftyPrecent_Filtered=(float *)malloc(FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	//	Edge2_Plus->Fmax=(float *)malloc(FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	//	Edge2_Plus->PrevLastFiftyPrecent=(float *)malloc((FgParams->Motion.AvgValue-1) * sizeof(float));
	//	Edge2_Plus->SumEnergy_Post=(float *)malloc(FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	//	Edge2_Plus->T1_t=(float *)malloc(FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	//
	//	if(FgParams->Motion.TruncateHilbert==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
	//		Edge2_Plus->Peak_Filtered=(float *)malloc((FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
	//	}
	//	else{//Considers the signal to be zero beyond the endpoints.
	//		Edge2_Plus->Peak_Filtered=(float *)malloc((FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue-1) * sizeof(float));
	//	}
	//
	//	Edge2_Minus->FiftyPrecent_Filtered=(float *)malloc(FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	//	Edge2_Minus->Fmax=(float *)malloc(FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	//	Edge2_Minus->PrevLastFiftyPrecent=(float *)malloc((FgParams->Motion.AvgValue-1) * sizeof(float));
	//	Edge2_Minus->SumEnergy_Post=(float *)malloc(FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	//	Edge2_Minus->T1_t=(float *)malloc(FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	//
	//	if(FgParams->Motion.TruncateHilbert==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
	//		Edge2_Minus->Peak_Filtered=(float *)malloc((FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
	//	}
	//	else{//Considers the signal to be zero beyond the endpoints.
	//		Edge2_Minus->Peak_Filtered=(float *)malloc((FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue-1) * sizeof(float));
	//	}
	//
	//	//	MotionStruct->Edge2=&Edge2;
	//	//	MotionStruct->Edge2_Plus=&Edge2_Plus;
	//	//	MotionStruct->Edge2_Minus=&Edge2_Minus;
	//
	//

	Edge2->FiftyPrecent = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2->FiftyPrecent_Filtered = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2->Fmax = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2->PeakEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2->PeakEnergy_PM = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2->maxFreqEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2->FreqBins = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2->maxFreqIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2->maxPeakIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2->FiftyPrecentIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2->PrevLastFiftyPrecent = (float *) malloc(
			(FgParams->Motion.AvgValue - 1) * sizeof(float));

	if (FgParams->Motion.TruncateReal == 1) {	// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and begining
		Edge2->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2->Peak,0,(FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		Edge2->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		memset(Edge2->Peak, 0,
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Plus->FiftyPrecent = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus->FiftyPrecent_Filtered = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus->Fmax = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus->PeakEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus->PeakEnergy_PM = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus->SumEnergy_Post = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus->maxFreqEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus->FreqBins = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus->maxFreqIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus->maxPeakIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus->FiftyPrecentIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Plus->PrevLastFiftyPrecent = (float *) malloc(
			(FgParams->Motion.AvgValue - 1) * sizeof(float));
	Edge2_Plus->T1_t = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));

	memset(Edge2_Plus->SumEnergy_Post, 0,
			(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));

	if (FgParams->Motion.TruncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Plus->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_Plus->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Plus->Peak,0,(FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the inital conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Plus->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		Edge2_Plus->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Plus->Peak, 0,
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}

	Edge2_Minus->FiftyPrecent = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus->FiftyPrecent_Filtered = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus->Fmax = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus->PeakEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus->PeakEnergy_PM = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus->SumEnergy_Post = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus->maxFreqEnergy = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus->FreqBins = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus->maxFreqIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus->maxPeakIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus->FiftyPrecentIdxs = (int *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(int));
	Edge2_Minus->PrevLastFiftyPrecent = (float *) malloc(
			(FgParams->Motion.AvgValue - 1) * sizeof(float));
	Edge2_Minus->T1_t = (float *) malloc(
			FgParams->Motion.SpectrogramTimeBins * sizeof(float));

	memset(Edge2_Minus->SumEnergy_Post, 0,
			(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));

	if (FgParams->Motion.TruncateHilbert == 1) {// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Minus->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_Minus->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins) * sizeof(float));
		//		memset(Edge2_Minus->Peak,0,(FgParams->Motion.SpectrogramTimeBins+FgParams->Motion.MedianValue/2)*sizeof(float));//set 0 the first MedianValue/2 for the initial conditions
	} else {		//Considers the signal to be zero beyond the endpoints.
		Edge2_Minus->Peak = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		Edge2_Minus->Peak_Filtered = (float *) malloc(
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));
		memset(Edge2_Minus->Peak, 0,
				(FgParams->Motion.SpectrogramTimeBins + FgParams->Motion.MedianValue - 1)
				* sizeof(float));//set 0 the first MedianValue/2 and last MedianValue/2-1 for the inital conditions
	}


	return 0;
}

int FreeMemoryEdeges(Edge2_Struct* Edge2_0,
		Edge2_Struct* Edge2_Plus_0, Edge2_Struct* Edge2_Minus_0){
	free(Edge2_0->FiftyPrecent);
	free(Edge2_0->FiftyPrecent_Filtered);
	free(Edge2_0->Fmax);
	free(Edge2_0->PeakEnergy);
	free(Edge2_0->PeakEnergy_PM);
	free(Edge2_0->maxFreqEnergy);
	free(Edge2_0->FreqBins);
	free(Edge2_0->maxPeakIdxs);
	free(Edge2_0->maxFreqIdxs);
	free(Edge2_0->FiftyPrecentIdxs);
	free(Edge2_0->PrevLastFiftyPrecent);
	free(Edge2_0->Peak);
	free(Edge2_0->Peak_Filtered);

	free(Edge2_Plus_0->FiftyPrecent);
	free(Edge2_Plus_0->FiftyPrecent_Filtered);
	free(Edge2_Plus_0->Fmax);
	free(Edge2_Plus_0->PeakEnergy);
	free(Edge2_Plus_0->PeakEnergy_PM);
	free(Edge2_Plus_0->maxFreqEnergy);
	free(Edge2_Plus_0->FreqBins);
	free(Edge2_Plus_0->maxPeakIdxs);
	free(Edge2_Plus_0->maxFreqIdxs);
	free(Edge2_Plus_0->FiftyPrecentIdxs);
	free(Edge2_Plus_0->PrevLastFiftyPrecent);
	free(Edge2_Plus_0->Peak);
	free(Edge2_Plus_0->Peak_Filtered);

	free(Edge2_Minus_0->FiftyPrecent);
	free(Edge2_Minus_0->FiftyPrecent_Filtered);
	free(Edge2_Minus_0->Fmax);
	free(Edge2_Minus_0->PeakEnergy);
	free(Edge2_Minus_0->PeakEnergy_PM);
	free(Edge2_Minus_0->maxFreqEnergy);
	free(Edge2_Minus_0->FreqBins);
	free(Edge2_Minus_0->maxPeakIdxs);
	free(Edge2_Minus_0->maxFreqIdxs);
	free(Edge2_Minus_0->FiftyPrecentIdxs);
	free(Edge2_Minus_0->PrevLastFiftyPrecent);
	free(Edge2_Minus_0->Peak);
	free(Edge2_Minus_0->Peak_Filtered);
	return 0;
}

#ifdef __cplusplus
}
#endif
