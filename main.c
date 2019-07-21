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


SysParams_Struct SysParams;
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

	SysParams.Motion.Nbins=MotionParamsFromXML.Nbins;
	SysParams.Motion.Nscans=MotionParamsFromXML.Nscans; // 400 for 1.6 sec record, 800 in the single 3.2 sec record
	SysParams.Motion.Fs=MotionParamsFromXML.Fs;
	SysParams.Motion.RF_NumOfClasses=MotionParamsFromXML.RF_NumOfClasses;
	SysParams.Motion.SpectrogramRangeWidth=MotionParamsFromXML.SpectrogramRangeWidth;
	SysParams.Motion.SpectrogramWinLength=MotionParamsFromXML.SpectrogramWinLength;
	SysParams.Motion.SpectrogramN_Overlap=MotionParamsFromXML.SpectrogramN_Overlap;
	SysParams.Motion.SpectrogramNoiseFreqBins=MotionParamsFromXML.SpectrogramNoiseFreqBins;
	SysParams.Motion.FminBin=MotionParamsFromXML.FminBin;
	SysParams.Motion.Fbins=MotionParamsFromXML.Fbins;
	SysParams.Motion.TopMaxFreq=MotionParamsFromXML.TopMaxFreq;
	SysParams.Motion.NoiseThresh=MotionParamsFromXML.NoiseThresh;
	SysParams.Motion.MinFreqPlusMinus=MotionParamsFromXML.MinFreqPlusMinus;//in Hz
	SysParams.Motion.MinFreqReal = MotionParamsFromXML.MinFreqReal; //Fmin in Hz
	SysParams.Motion.MinEventDuration=MotionParamsFromXML.MinEventDuration;
	SysParams.Motion.DFTLengthForPSD=floor(SysParams.Motion.Nscans/2)+1;
	SysParams.Motion.DFTLengthForSpectrogram=2*SysParams.Motion.Fbins-1;//=199
	SysParams.Motion.GapLength=MotionParamsFromXML.GapLength;
	SysParams.Motion.DFTLengthForSpectrogramHilbert=2*SysParams.Motion.Fbins;//=200
	SysParams.Motion.MedianValue=MotionParamsFromXML.MedianValue;//for MedianFilter
	SysParams.Motion.AvgValue=MotionParamsFromXML.AvgValue;//for AvgFilter
	SysParams.Motion.TruncateReal=MotionParamsFromXML.TruncateReal;// Computes medians of smaller segments as it reaches the signal edges in Pxx2 (regular) Spectrogram
	SysParams.Motion.TruncateHilbert=MotionParamsFromXML.TruncateHilbert;// Computes medians of smaller segments as it reaches the signal edges in Hilbert Spectrogram
	SysParams.Motion.SpectrogramTimeShift=SysParams.Motion.SpectrogramWinLength-SysParams.Motion.SpectrogramN_Overlap;//=5
	SysParams.Motion.NumSamplesForDerivativeEstimation=MotionParamsFromXML.NumSamplesForDerivativeEstimation;

	SysParams.Motion.SpectrogramFreqBins=SysParams.Motion.DFTLengthForSpectrogram/2+1;//=100
	SysParams.Motion.SpectrogramTimeBinsSingleMotion = floor((SysParams.Motion.Nscans-SysParams.Motion.SpectrogramN_Overlap)/(SysParams.Motion.SpectrogramWinLength - SysParams.Motion.SpectrogramN_Overlap)); // total time bins of the spectrogram=155
	SysParams.Motion.SpectrogramTimeBinsTwoMotions=(2*SysParams.Motion.SpectrogramTimeBinsSingleMotion+SysParams.Motion.GapLength);//=2*75+14=164
	SysParams.Motion.SpectrogramTimeBinsThreeMotions=(3*SysParams.Motion.SpectrogramTimeBinsSingleMotion+2*SysParams.Motion.GapLength);//=3*75+2*14=253
	SysParams.Motion.SpectrogramFreqBinsHilbert=2*SysParams.Motion.Fbins;//=200
	SysParams.Motion.Hamming[0]=0.08; SysParams.Motion.Hamming[1]=0.0907545443733602; SysParams.Motion.Hamming[2]=0.122515306951360;
	SysParams.Motion.Hamming[3]=0.173797189775404; SysParams.Motion.Hamming[4]=0.242202309000359; SysParams.Motion.Hamming[5]=0.324532117278097;
	SysParams.Motion.Hamming[6]=0.416936964276558;SysParams.Motion.Hamming[7]=0.515096102050708; SysParams.Motion.Hamming[8]=0.614419718414272;
	SysParams.Motion.Hamming[9]=0.710263551456361; SysParams.Motion.Hamming[10]=0.798146050066696;SysParams.Motion.Hamming[11]=0.873957926284640;
	SysParams.Motion.Hamming[12]=0.934154301037091; SysParams.Motion.Hamming[13]=0.975920458744089; SysParams.Motion.Hamming[14]=0.997303460291006;
	SysParams.Motion.Hamming[15]=SysParams.Motion.Hamming[14]; SysParams.Motion.Hamming[16]=SysParams.Motion.Hamming[13]; SysParams.Motion.Hamming[17]=SysParams.Motion.Hamming[12]; SysParams.Motion.Hamming[18]=SysParams.Motion.Hamming[11];
	SysParams.Motion.Hamming[19]=SysParams.Motion.Hamming[10]; SysParams.Motion.Hamming[20]=SysParams.Motion.Hamming[9]; SysParams.Motion.Hamming[21]=SysParams.Motion.Hamming[8]; SysParams.Motion.Hamming[22]=SysParams.Motion.Hamming[7];
	SysParams.Motion.Hamming[23]=SysParams.Motion.Hamming[6]; SysParams.Motion.Hamming[24]=SysParams.Motion.Hamming[5]; SysParams.Motion.Hamming[25]=SysParams.Motion.Hamming[4]; SysParams.Motion.Hamming[26]=SysParams.Motion.Hamming[3];
	SysParams.Motion.Hamming[27]=SysParams.Motion.Hamming[2]; SysParams.Motion.Hamming[28]=SysParams.Motion.Hamming[1]; SysParams.Motion.Hamming[29]=SysParams.Motion.Hamming[0];
	SysParams.Motion.A12_Inverse[0][0]=5.925925925925995e-04;SysParams.Motion.A12_Inverse[0][1]=-5.925925925925996e-04;SysParams.Motion.A12_Inverse[0][2]=0.004444444444445;SysParams.Motion.A12_Inverse[0][3]=0.004444444444444;
	SysParams.Motion.A12_Inverse[1][0]=-0.146666666666668;SysParams.Motion.A12_Inverse[1][1]=0.146666666666668;SysParams.Motion.A12_Inverse[1][2]=-1.133333333333347;SysParams.Motion.A12_Inverse[1][3]=-1.066666666666678;
	SysParams.Motion.A12_Inverse[2][0]=12.0000000000001;SysParams.Motion.A12_Inverse[2][1]=-12.0000000000001;SysParams.Motion.A12_Inverse[2][2]=96.0000000000011;SysParams.Motion.A12_Inverse[2][3]=85.0000000000010;
	SysParams.Motion.A12_Inverse[3][0]=-324.000000000004;SysParams.Motion.A12_Inverse[3][1]=325.000000000004;SysParams.Motion.A12_Inverse[3][2]=-2.700000000000030e+03;SysParams.Motion.A12_Inverse[3][3]=-2.250000000000026e+03;


	//////Global Parameters/////

	bgParams.Rmin=1;//need to change this from the configuration XML!!!!
	bgParams.Rmax=3.5; //need to change this from the configuration XML!!!!

	SysParams.Motion.RealRingSizeForMotionTracking=2.62572418212891; //the ring width in meters; FOR DEMO ONLY maybe 2.5??? why 2.6

	bgParams.Rstart_m=2.231486865000000;//for demo only!!! it will be the current Rstart

	//create Rbin_m
	SysParams.Motion.Rbin_m=(float *)malloc(SysParams.Motion.Nbins * sizeof(float));
	float delta_R=(SysParams.Motion.RealRingSizeForMotionTracking)/((float)SysParams.Motion.Nbins-1);//the delta of the bins inside the ring(=Rstop-Rstart)=2.5m

	for(i=0;i<SysParams.Motion.Nbins;i++){
		SysParams.Motion.Rbin_m[i]=bgParams.Rstart_m+delta_R*i;
	}

	int y_hat_M;
		float MotionDistribution[SysParams.Motion.RF_NumOfClasses]; //=4

	///////Classifiers Preparation////////
	int NumOfTrees=8,NumOfFeatures=6;//RandomForrest
	RF_Params_Import(NumOfTrees, NumOfFeatures,All_Trees,&RF_Model);
	SVM_Params_Import(&SVM_Model);

	//////Create FFT Plans//////
	CreateFFTPlans(&SysParams);

	////////MotionStruct0 Preparation////////
	Motion_Struct MotionStruct0;

	Edge2_Struct Edge2_0;
	Edge2_Struct Edge2_Plus_0;
	Edge2_Struct Edge2_Minus_0;
	SysParams.Motion.SpectrogramTimeBins=SysParams.Motion.SpectrogramTimeBinsSingleMotion;

	if(SysParams.Motion.TruncateReal==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and begining
		Edge2_0.Peak=(float *)malloc((SysParams.Motion.SpectrogramTimeBins) * sizeof(float));
		Edge2_0.Peak_Filtered=(float *)malloc((SysParams.Motion.SpectrogramTimeBins) * sizeof(float));
	}
	else{//Considers the signal to be zero beyond the endpoints.
		Edge2_0.Peak_Filtered=(float *)malloc((SysParams.Motion.SpectrogramTimeBins+SysParams.Motion.MedianValue-1) * sizeof(float));
	}

	Edge2_Plus_0.FiftyPrecent_Filtered=(float *)malloc(SysParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_0.Fmax=(float *)malloc(SysParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_0.PrevLastFiftyPrecent=(float *)malloc((SysParams.Motion.AvgValue-1) * sizeof(float));
	Edge2_Plus_0.SumEnergy_Post=(float *)malloc(SysParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Plus_0.T1_t=(float *)malloc(SysParams.Motion.SpectrogramTimeBins * sizeof(float));

	if(SysParams.Motion.TruncateHilbert==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Plus_0.Peak_Filtered=(float *)malloc((SysParams.Motion.SpectrogramTimeBins) * sizeof(float));
	}
	else{//Considers the signal to be zero beyond the endpoints.
		Edge2_Plus_0.Peak_Filtered=(float *)malloc((SysParams.Motion.SpectrogramTimeBins+SysParams.Motion.MedianValue-1) * sizeof(float));
	}

	Edge2_Minus_0.FiftyPrecent_Filtered=(float *)malloc(SysParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_0.Fmax=(float *)malloc(SysParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_0.PrevLastFiftyPrecent=(float *)malloc((SysParams.Motion.AvgValue-1) * sizeof(float));
	Edge2_Minus_0.SumEnergy_Post=(float *)malloc(SysParams.Motion.SpectrogramTimeBins * sizeof(float));
	Edge2_Minus_0.T1_t=(float *)malloc(SysParams.Motion.SpectrogramTimeBins * sizeof(float));

	if(SysParams.Motion.TruncateHilbert==1){// Computes medians of smaller segments as it reaches the signal edges. no zeropadding at the end and beginning
		Edge2_Minus_0.Peak_Filtered=(float *)malloc((SysParams.Motion.SpectrogramTimeBins) * sizeof(float));
	}
	else{//Considers the signal to be zero beyond the endpoints.
		Edge2_Minus_0.Peak_Filtered=(float *)malloc((SysParams.Motion.SpectrogramTimeBins+SysParams.Motion.MedianValue-1) * sizeof(float));
	}

	MotionStruct0.Edge2=&Edge2_0;
	MotionStruct0.Edge2_Plus=&Edge2_Plus_0;
	MotionStruct0.Edge2_Minus=&Edge2_Minus_0;



	///////////////END OF INIT//////////////////////



	////////GET MSCAN0 AND MSCAN1////////

	int j;
	float* Mscan_PostProcess0[SysParams.Motion.Nscans];
	float* Mscan_PostProcess1[SysParams.Motion.Nscans];

	for (i = 0; i < SysParams.Motion.Nscans; i++)
		Mscan_PostProcess0[i] = (float *) malloc(
				SysParams.Motion.Nbins * sizeof(float));
	for (i = 0; i < SysParams.Motion.Nscans; i++)
		Mscan_PostProcess1[i] = (float *) malloc(
				SysParams.Motion.Nbins * sizeof(float));

	float Mscan_flat[SysParams.Motion.Nscans*SysParams.Motion.Nbins];
	float* Mscan0[SysParams.Motion.Nscans],* Mscan1[SysParams.Motion.Nscans],*Mscan2[SysParams.Motion.Nscans],* Mscan3[SysParams.Motion.Nscans];//,*Mscan_old[SysParams.Motion.Nscans];
	for (i=0; i<SysParams.Motion.Nscans; i++){
		Mscan0[i] = (float *)malloc(SysParams.Motion.Nbins * sizeof(float));
		Mscan1[i] = (float *)malloc(SysParams.Motion.Nbins * sizeof(float));
		Mscan2[i] = (float *)malloc(SysParams.Motion.Nbins * sizeof(float));
		Mscan3[i] = (float *)malloc(SysParams.Motion.Nbins * sizeof(float));

	}


	for(SysParams.RecordNumber=700;SysParams.RecordNumber<701;SysParams.RecordNumber++){
		printf("record number %d\n",SysParams.RecordNumber);
		char  Filename_Mscan1[50] ;
		sprintf(Filename_Mscan1, "/home/debian/Records/Mscan_withoutMF%d_1.csv" , SysParams.RecordNumber);
		FILE *MscanFile1=fopen(Filename_Mscan1, "r");

		char  Filename_Mscan2[50] ;
		sprintf(Filename_Mscan2, "/home/debian/Records/Mscan_withoutMF%d_2.csv" , SysParams.RecordNumber);
		FILE *MscanFile2=fopen(Filename_Mscan2, "r");

		char  Filename_Mscan3[50] ;
		sprintf(Filename_Mscan3, "/home/debian/Records/Mscan_withoutMF%d_3.csv" , SysParams.RecordNumber);
		FILE *MscanFile3=fopen(Filename_Mscan3, "r");

		read_data_from_file(MscanFile1,SysParams.Motion.Nscans,SysParams.Motion.Nbins,Mscan_flat);//get the MSCAN in flat form

		for ( i=0;i<SysParams.Motion.Nscans;i++)
		{//get the  MSCAN from flat
			for ( j=0;j<SysParams.Motion.Nbins;j++)
			{
				Mscan0[i][j]=Mscan_flat[i*SysParams.Motion.Nbins+j];
			}
		}

		read_data_from_file(MscanFile1,SysParams.Motion.Nscans,SysParams.Motion.Nbins,Mscan_flat);//get the MSCAN in flat form
		for ( i=0;i<SysParams.Motion.Nscans;i++)
		{//get the  MSCAN from flat
			for ( j=0;j<SysParams.Motion.Nbins;j++)
			{
				Mscan1[i][j]=Mscan_flat[i*SysParams.Motion.Nbins+j];
			}
		}

		read_data_from_file(MscanFile2,SysParams.Motion.Nscans,SysParams.Motion.Nbins,Mscan_flat);//get the MSCAN in flat form
		for ( i=0;i<SysParams.Motion.Nscans;i++)
		{//get the  MSCAN from flat
			for ( j=0;j<SysParams.Motion.Nbins;j++)
			{
				Mscan2[i][j]=Mscan_flat[i*SysParams.Motion.Nbins+j];
			}
		}


		read_data_from_file(MscanFile3,SysParams.Motion.Nscans,SysParams.Motion.Nbins,Mscan_flat);//get the MSCAN in flat form
		for ( i=0;i<SysParams.Motion.Nscans;i++)
		{//get the  MSCAN from flat
			for ( j=0;j<SysParams.Motion.Nbins;j++)
			{
				Mscan3[i][j]=Mscan_flat[i*SysParams.Motion.Nbins+j];
			}
		}




		gettimeofday(&tpStart,0);
		//		SysParams.Motion.FirstMscan=1;


		memset(MotionDistribution, 0, SysParams.Motion.RF_NumOfClasses * sizeof(float));	//initialize with zeros


		SysParams.Motion.FirstTimeMotion=1;//1 if it's the first 2 frames

		MotionPreprocess(Mscan0,Mscan_PostProcess0,&SysParams);
		MotionTracking(Mscan_PostProcess0, &SysParams,&bgParams);//Send Mscan0 to tracking//


		MotionPreprocess(Mscan1,Mscan_PostProcess1,&SysParams);
		MotionTracking(Mscan_PostProcess1, &SysParams,&bgParams);//Send Mscan1 to tracking//

		MotionAnalyzer(Mscan0, Mscan1, Mscan_PostProcess0,Mscan_PostProcess1,&SysParams,All_Trees,&SVM_Model,&RF_Model,&MotionStruct0,&bgParams,&y_hat_M,MotionDistribution);


		gettimeofday(&tpStop,0); f1 = ( (float)( tpStop.tv_sec-tpStart.tv_sec)+ (float)(tpStop.tv_usec)/1000000 ) -  ((float)(tpStart.tv_usec)/1000000) ;printf (" %f sec\n", f1 );



		MotionPreprocess(Mscan2,Mscan_PostProcess0,&SysParams);
		MotionTracking(Mscan_PostProcess0, &SysParams,&bgParams);//Send Mscan0 to tracking//


		MotionPreprocess(Mscan3,Mscan_PostProcess1,&SysParams);
		MotionTracking(Mscan_PostProcess1, &SysParams,&bgParams);//Send Mscan1 to tracking//


		SysParams.Motion.FirstTimeMotion=0;////after 1st time

		memset(MotionDistribution, 0, SysParams.Motion.RF_NumOfClasses * sizeof(float));	//initialize with zeros

		MotionAnalyzer(Mscan2, Mscan3, Mscan_PostProcess0,Mscan_PostProcess1,&SysParams,All_Trees,&SVM_Model,&RF_Model,&MotionStruct0,&bgParams,&y_hat_M,MotionDistribution);

		for(i=0;i<4;i++){

			printf("%f\n",MotionDistribution[i]);
		}


		//		MotionAlgorithm(Mscan2,Mscan3,&SysParams,All_Trees,&SVM_Model,&RF_Model,&MotionStruct0,&bgParams);


		//free memory



	}

	for (i=0; i<SysParams.Motion.Nscans; i++){
				free(Mscan0[i]);
				free(Mscan1[i]);
				free(Mscan2[i]);
				free(Mscan3[i]);
			}
	fftwf_destroy_plan(SysParams.Motion.FFT_Hilbert);
	fftwf_destroy_plan(SysParams.Motion.IFFT_Hilbert);
	fftwf_destroy_plan(SysParams.Motion.FFT_HilbertSpectrogram);
	fftwf_destroy_plan(SysParams.Motion.FFT_ABS);
	ne10_fft_destroy_r2c_float32(SysParams.Motion.FFT_Feature42);


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
	int i,NumOfFeatures=7;

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

int CreateFFTPlans(SysParams_Struct *SysParams){


	int flags = 1;
	float fftin_local_ptr_complex[SysParams->Motion.Nbins];
	memset(fftin_local_ptr_complex, 0, sizeof(fftin_local_ptr_complex));
	SysParams->Motion.FFT_Hilbert = fftwf_plan_dft_r2c_1d(SysParams->Motion.Nbins,
			fftin_local_ptr_complex, SysParams->Motion.MscanFFT, flags);

	SysParams->Motion.IFFT_Hilbert = fftwf_plan_dft_1d(SysParams->Motion.Nbins, SysParams->Motion.MscanFFT, SysParams->Motion.MscanIFFT,
			FFTW_BACKWARD, FFTW_MEASURE);	//MAYBE ESTIMATEEEEE!!!!!




	SysParams->Motion.FFT_HilbertSpectrogram = fftwf_plan_dft_r2c_1d(SysParams->Motion.Nbins,
			fftin_local_ptr_complex, SysParams->Motion.MscanFFT, flags);


	SysParams->Motion.SignalForFFT = (fftwf_complex*) fftwf_malloc( sizeof(fftwf_complex) * SysParams->Motion.SpectrogramFreqBinsHilbert);	//pay attention

	SysParams->Motion.SignalForFFTABS = (fftwf_complex*) fftwf_malloc( sizeof(fftwf_complex) * SysParams->Motion.Nscans);	//pay attention



	SysParams->Motion.FFTResult = (fftwf_complex*) fftwf_malloc(
			sizeof(fftwf_complex) * SysParams->Motion.SpectrogramFreqBinsHilbert);


	SysParams->Motion.FFTResultABS = (fftwf_complex*) fftwf_malloc(
			sizeof(fftwf_complex) * (SysParams->Motion.Nscans));

	SysParams->Motion.FFT_HilbertSpectrogram=fftwf_plan_dft_1d(SysParams->Motion.SpectrogramFreqBinsHilbert,
			SysParams->Motion.SignalForFFT, SysParams->Motion.FFTResult, FFTW_FORWARD, FFTW_ESTIMATE);

	SysParams->Motion.FFT_ABS= fftwf_plan_dft_1d(SysParams->Motion.Nscans,
			SysParams->Motion.SignalForFFTABS, SysParams->Motion.FFTResultABS, FFTW_FORWARD, FFTW_ESTIMATE);

	// An FFT "configuration structure"
	SysParams->Motion.FFT_Feature42 = ne10_fft_alloc_r2c_float32(256);

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
