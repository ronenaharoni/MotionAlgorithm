#include <stdio.h>
#include <fftw3.h>
#include <NE10.h>
#include "seatest.h"
#include "NE10_dsp.h"



typedef struct {//need to change to flexiable size
	float* Fmax;
	float* Peak;
	float* Peak_Filtered;
	float* FiftyPrecent;
	float* FiftyPrecent_Filtered;
	int* FiftyPrecentIdxs;
	float* PeakEnergy;
	//	float* SumEnergy;
	float* maxFreqEnergy;
	int* maxFreqIdxs;
	int* maxPeakIdxs;
	int* FreqBins;
	float* PeakEnergy_PM;
	float* SumEnergy_Post;
	float SNR_Linear;
	float SNR_Both;
	float* PrevLastFiftyPrecent;
	float* T1_t;

} Edge2_Struct;


typedef struct
{
	float        Rstart_m ;               // 1 (null)
	float        Rstop_m ;                // 3.5 (null)
	float          Rmin;
	float          Rmax;//ADD THIS AT RONENES XML
}BgRadarParams;


//typedef struct {
//	int Nscans;
//	int Nbins;
//	int Fs;
//	float Rmin;
//	float Rmax;
//	float Rstart_corrected;
//	float Rstop_corrected;
//	float *Rbin_m;
//	float Ring;
//	int NumOfClasses;
//	int RangeWide;
//	int N_Overlap;
//	int winLength;
//	int TimeShift;
//	int Fbins;
//	int noiseFreqBins;
//	int Fmin_bin;
//	float MinFreqPlusMinus;
//	float MinFreqReal;
//	float Hamming[30];
//	int DFTLengthForPSD;
//	int noiseThresh;
//	int DFTLengthForSpectrogram;
//	int SpectrogramFreqBins;
//	int SpectrogramTimeBins;
//	int SpectrogramTimeBinsSingleMotion;
//	int SpectrogramTimeBinsTwoMotions;
//	int SpectrogramTimeBinsThreeMotions;
//	int SpectrogramFreqBinsHilbert;
//	int DFTLengthForSpectrogramHilbert;
//	int MedianValue;
//	int AvgValue;
//	int truncate;
//	int truncateHilbert;
//	int FirstTimeMotion;
//	float A12_Inverse[4][4];
//	int NumSamplesForDerivativeEstimation;
//	int GapLength;
//	int minEventDuration;
//	int topMaxFreq;
//	int RecordNumber;
//	fftwf_plan FFT_Hilbert;
//	fftwf_plan IFFT_Hilbert;
//	fftwf_plan FFT_HilbertSpectrogram;
//	fftwf_plan FFT_ABS;
//
//	fftwf_complex MscanFFT[288];
//	fftwf_complex MscanIFFT[288];
//	fftwf_complex* SignalForFFT;
//	fftwf_complex* SignalForFFTABS;
//	fftwf_complex*  FFTResult;
//	fftwf_complex*  FFTResultABS;
//	ne10_fft_r2c_cfg_float32_t FFT_Feature42;
//
//
//} SysParams_Struct;


typedef struct {
	int FirstTimeMotion;
	int Nscans;
	int Nbins;
	int Fs;
	int RF_NumOfClasses;
	int SpectrogramRangeWidth;
	int SpectrogramWinLength;
	int SpectrogramN_Overlap;
	int SpectrogramNoiseFreqBins;
	int SpectrogramTimeShift;
	int FminBin;
	int Fbins;
	int NoiseThresh;
	int TopMaxFreq;
	int MinFreqPlusMinus;
	float MinFreqReal;
	int MinEventDuration;
	int GapLength;
	int NumSamplesForDerivativeEstimation;
	int MedianValue;
	int AvgValue;
	int TruncateReal;
	int TruncateHilbert;
	int DFTLengthForPSD;
	int DFTLengthForSpectrogram;
	int SpectrogramFreqBins;
	int SpectrogramTimeBins;
	int SpectrogramTimeBinsSingleMotion;
	int SpectrogramTimeBinsTwoMotions;
	int SpectrogramTimeBinsThreeMotions;
	int SpectrogramFreqBinsHilbert;
	int DFTLengthForSpectrogramHilbert;
	float A12_Inverse[4][4];
	float Hamming[30];
	float RealRingSizeForMotionTracking;
	float *Rbin_m;
	fftwf_plan FFT_Hilbert;
	fftwf_plan IFFT_Hilbert;
	fftwf_plan FFT_HilbertSpectrogram;
	fftwf_plan FFT_ABS;
	fftwf_complex MscanFFT[288];
	fftwf_complex MscanIFFT[288];
	fftwf_complex* SignalForFFT;
	fftwf_complex* SignalForFFTABS;
	fftwf_complex*  FFTResult;
	fftwf_complex*  FFTResultABS;
	ne10_fft_r2c_cfg_float32_t FFT_Feature42;
int FirstMscan;

} Motion_Params;



typedef struct {
	int RecordNumber;
	Motion_Params Motion;

} SysParams_Struct;







typedef struct {
	float* Linear[101];
	float* dB[101];
	float T1_t[155];
	float T1_t_mean;
} Pxx2_Plus_Struct;

typedef struct {
	float* Linear[100];
	float* dB[100];
	float T1_t[155];
	float T1_t_mean;

} Pxx2_Minus_Struct;


typedef struct {//need to change to flexiable size
	float Fmax[155];
	float Peak[174];
	float Peak_Filtered[155];
	float FiftyPrecent[155];
	float FiftyPrecent_Filtered[155];
	int FiftyPrecentIdxs[155];
	float PeakEnergy[155];
	float SumEnergy[155];
	float maxFreqEnergy[155];
	int maxFreqIdxs[155];
	int maxPeakIdxs[155];
	int FreqBins[155];
	float PeakEnergy_PM[155];
	//	float SumEnergy_Post[155];
	float SNR_Linear;
	float SNR_Both;
} Edge_Struct;



typedef struct {
	int chosenStart;
	int chosenEnd;
	int LengthOfMaxEvent;
	int LengthOfChosenEvent;
} Event_Struct;





typedef struct {
	Edge2_Struct* Edge2;
	Edge2_Struct* Edge2_Plus;
	Edge2_Struct* Edge2_Minus;
	Event_Struct *EventStruct;
	int EventPlusPassEvent;
	int EventMinusPassEvent;

}
Motion_Struct;



typedef struct {
	Motion_Struct* MotionStruct0;
	Motion_Struct* MotionStruct1;
	Motion_Struct* MotionStruct2;
	Motion_Struct* UnitedMotionStruct;

}
AllMotion_Structs;


typedef struct {
	float Edge2_50Precent_Fmax;
	float Edge2_50Precent_AvgTopFive;
	float Edge2_50Precent_PeakToAvg;
	float FmaxFpeakMultiplication;
	float Edge2_maxPeakFreq_fromTimeBinWithMaxFreq;
	float HilbertRatio;
	float yDiff_Raise_Max_FiftyPrecent_Filtered;
	float yDiff_Fall_Max_FiftyPrecent_Filtered;
	float yDiff_Raise_Mean_FiftyPrecent_Filtered;
	float yDiff_Fall_Mean_FiftyPrecent_Filtered;
	float  Max_Black_Curve;
	int maxEventLength;
	float SumEnergy42;
	float p1;
	float SNR_Both;
} Features_Struct;

typedef struct {
	Features_Struct* Plus;
	Features_Struct* Minus;
	Features_Struct* Both;

}AllFeatures_Struct;



typedef struct {
	int* Children;
	float* CutPoint;
	int* CutPredictor;
	float* ClassProb;
	int TotalTrees;
} Tree_Struct;

typedef struct {
	double x_mean[8];
	double x_std[8];
	//	Tree_Struct *All_Trees;
} RF_Struct;

typedef struct {
	float Bias;
	float Beta[7];
	float x_mean[7];
	float x_std[7];
} SVM_Struct;



int RF_Params_Import(int NumOfTrees,int NumOfFeatures, Tree_Struct** All_Trees,RF_Struct* RF_Model);
int TreesCreator(int NumOfTrees,Tree_Struct** All_Trees);
int SVM_Params_Import(SVM_Struct* SVM_Model);
int CreateFFTPlans(SysParams_Struct *SysParams);
int CreateEdges(SysParams_Struct* SysParams);
int MotionPreprocess(float* Mscan[],float* Mscan_PostProcess[],SysParams_Struct* SysParams);

int MotionAnalyzer(float* Mscan1[],float* Mscan2[],float* Mscan_PostProcess1[],float* Mscan_PostProcess2[],SysParams_Struct* SysParams,Tree_Struct** All_Trees,SVM_Struct* SVM_Model,RF_Struct* RF_Model,Motion_Struct* MotionStruct0,BgRadarParams *bgParams,int *y_hat_M,float* MotionDistribution);

void MemoryAllocation(Edge2_Struct* Edge2_1,Edge2_Struct* Edge2_2,Edge2_Struct* Edge2_Plus_1,Edge2_Struct* Edge2_Minus_1,Edge2_Struct* Edge2_Plus_2,Edge2_Struct* Edge2_Minus_2,SysParams_Struct *SysParams);
void FreeMemory(Edge2_Struct* Edge2_1,Edge2_Struct* Edge2_2,Edge2_Struct* Edge2_Plus_1,Edge2_Struct* Edge2_Minus_1,Edge2_Struct* Edge2_Plus_2,Edge2_Struct* Edge2_Minus_2,Motion_Struct *UnitedMotionStruct);
int MatchedFilter(float* Mscan[], SysParams_Struct* SysParams);
int SlowProcessing2(float* Mscan[],float* Mscan_PostProcess[], SysParams_Struct* SysParams);
int SlowProcessingHilbert(_Complex float* Mscan[],_Complex float* Mscan_PostProcess[], SysParams_Struct* SysParams);

int NotchFilter2(float* Mscan[],SysParams_Struct* SysParams);
int NotchFilterHilbert(_Complex float* Mscan[], SysParams_Struct* SysParams);
int AbsOfFFT(float* Mscan[], float* Mscan_abs_FFT[],SysParams_Struct* SysParams);

int GET_ROI(float *Mscan_abs_FFT[], SysParams_Struct* SysParams, int *p1);

int Spectrogram2(float* Pxx2[],float* Pxx2_dB[],float *T2_dB,float* Mscan[],int p1,SysParams_Struct* SysParams);
int Hilbert(float* Mscan[], _Complex float* MscanIQ[], SysParams_Struct* SysParams);
int CurveLength(Motion_Struct *MotionStruct,SysParams_Struct *SysParams);

int MotionCurveExtraction(Motion_Struct* MotionStruct,float* Mscan[],float* Mscan_abs_FFT[],float* Pxx2_Hilbert[],float* Pxx2[],float* Pxx2_dB[],Edge2_Struct* Edge2,Edge2_Struct* Edge2_Plus,Edge2_Struct* Edge2_Minus,SysParams_Struct* SysParams,BgRadarParams *bgParams);
int MotionCurveExtraction2(Motion_Struct* MotionStruct, float* Mscan[],
		float* Mscan_abs_FFT[],float* Mscan_PostProcess[], float* Pxx2_Hilbert[], float* Pxx2[],
		float* Pxx2_dB[], Edge2_Struct* Edge2, Edge2_Struct* Edge2_Plus,
		Edge2_Struct* Edge2_Minus, SysParams_Struct* SysParams,BgRadarParams *bgParams);

int FeatureExtractionBasedCurves2(Motion_Struct* MotionStruct,AllFeatures_Struct* FeatureSet,SysParams_Struct* SysParams);
int CalcCurves2(float* Pxx2[],float* Pxx2_dB[],float *T2_dB,Edge2_Struct *Edge2,SysParams_Struct* SysParams);


int HilbertSpectrogram4(float* Pxx2_Hilbert[],float* Pxx2[], float* Pxx2_dB[], float *T2_dB,float* Mscan[], int p1, Edge2_Struct *Edge2_Plus,Edge2_Struct * Edge2_Minus,Edge2_Struct * Edge2, SysParams_Struct* SysParams);

int CalcCurvesHilbert2(Pxx2_Plus_Struct *Pxx2_Plus,Pxx2_Minus_Struct* Pxx2_Minus,Edge2_Struct* Edge2_Plus,Edge2_Struct* Edge2_Minus,SysParams_Struct* SysParams);
int MedianFilter(Edge2_Struct *Edge2, SysParams_Struct* SysParams,int MedianType);
int MedianFilter2(Edge2_Struct *Edge2,int SpectrogramTimeBins, int MedianValue, int truncate);//Filter the black curve (peak curve)
int read_data_from_file(FILE *fp_read ,int Nscans,int Nbins,float* RespMscan_flat);
int AvgFilter(Edge2_Struct *Edge2,SysParams_Struct* SysParams, int IsHilbert);

int AvgFilter2(Edge2_Struct *Edge2,int SpectrogramTimeBins,int AvgValue);
int ExtractFeatures2(Motion_Struct* MotionStruct,AllFeatures_Struct* FeatureSet,int Type,SysParams_Struct *SysParams);
int MedianFilterFor50Precent(Edge2_Struct *Edge2, Motion_Struct *MotionStruct, float *Edge2_50Precent_MedianFiltered,int MedianValue);
float MedianFilterFor50Precent2(Edge2_Struct *Edge2, Motion_Struct *MotionStruct, float *Edge2_50Precent_MedianFiltered,int MedianValue);
int MaxOfArr(float *Arr,int *MaxIdx, float *MaxValue, int Length);
int PolynomialFeatures2(Edge2_Struct* Edge2,Motion_Struct *MotionStruct,Features_Struct *Features);
int Feature42(Motion_Struct* MotionStruct, AllFeatures_Struct* FeatureSet,SysParams_Struct* SysParams);
int Feature42_2(Motion_Struct* MotionStruct,AllFeatures_Struct* FeatureSet);
int SNRFeature(Edge2_Struct* Edge2_Plus,Edge2_Struct* Edge2_Minus,AllFeatures_Struct *Featureset,SysParams_Struct *SysParams);
int RandomForrestClassifier(AllFeatures_Struct* FeatureSet,Tree_Struct** All_Trees,float * MotionDistribution,RF_Struct* RF_Model);
int ClassifierCorrection(AllFeatures_Struct* FeatureSet,float * MotionDistribution, int *y_hat_M,SVM_Struct* SVM_Model);
int SVMClassifier(AllFeatures_Struct* FeatureSet,float * MotionDistribution, int *y_hat_M,SVM_Struct *SVM_Model);
int GapInterpolation2( Motion_Struct *MotionStruct0,Motion_Struct *MotionStruct1,Motion_Struct *MotionStruct2,Motion_Struct *UnitedMotionStruct,SysParams_Struct* SysParams);
int GapInterpolation_2Curves(float* LeftCurve, float* RightCurve,
		float *IntrpolatedCurve, SysParams_Struct* SysParams, int FirstGap,int isHilbert);
float Max(float Val1, float Val2);
int MotionTracking(float* Mscan[],SysParams_Struct* SysParams,BgRadarParams *bgParams);
int K_means(float* energyPerTime_dB[],SysParams_Struct* SysParams,float* ChosenCentoird1,float* ChosenCentoird2);
int SaveToCsv(int y_hat,float* MotionDistribution,AllFeatures_Struct* FeatureSet,
		RF_Struct* RF_Model,SysParams_Struct *SysParams);

int SavePxx2ToCsv(float* Pxx2[],
		SysParams_Struct *SysParams);

int hilbert(float *signal, int N, float _Complex *output);
int Hilbert2(float* Mscan[],_Complex float* MscanIQ[],int Nbins,int Nscans);
int FreeMemoryEdeges(Edge2_Struct* Edge2_0,
		Edge2_Struct* Edge2_Plus_0, Edge2_Struct* Edge2_Minus_0);
